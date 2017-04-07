#include "PlaneGEM.h"
#include "GeomPlaneGEM.h"
#include "ChipAPV.h"

#include <TF1.h>
#include <cmath>
#include <stdexcept>

ClassImp(PlaneGEM);

PlaneGEM::PlaneGEM(const char *detname,int nchan, int center, int width, bool pixel)
  : PlaneAPV(detname,nchan,center,width, pixel), fPrintedMissCalib(false) {
  INeedGeom();

#if USE_DATABASE==1
  calib_data.clear();
  calib_time.clear();
#endif
}

void PlaneGEM::Init(TTree* tree) {
  // override the default settings for hittimes from PlaneAPV before
  // calling the init method, as there the histograms are created
  fVt->SetRange(-100., 100.);
  fVt->SetNbins(200);

  PlaneAPV::Init(tree);

  // Geometry and calibrations
  if(fGeom) {
    GeomPlaneGEM* geompl = (GeomPlaneGEM*) fGeom;
    geompl->SetParameters();
#if USE_DATABASE==1
    geompl->SetCalibrations(calib_data, fRunMaps, calib_time);
#else
    geompl->SetCalibrations(fRunMaps);
#endif
  }

  std::string name;

  // Channel amplitude ratio ("banana")
  name = fName + "_AmpRatio";
  fHampRatio = new TH2F(name.c_str(), name.c_str(), 100, -0.01, 1.99, 100, -0.01, 1.99);
  AddHistogram(fHampRatio);
  fHampRatio->GetXaxis()->SetTitle("a1/a2");
  fHampRatio->GetYaxis()->SetTitle("a0/a2");
  fHampRatio->SetOption("COLZ");

  // Channel vs strip amplitude
  name = fName + "_ch_vs_a2";
  fHchvsa2 = new TH2F(name.c_str(),name.c_str(), fNchan, -0.5, fNchan-0.5, 200, -5, 395);
  AddHistogram(fHchvsa2);
  fHchvsa2->SetOption("COLZ");

  // Channel vs common mode correction
  name = fName + "_chip_vs_cmc2";
  int chips = fNchan / 128; // one APV chip has 128 channels
  fHchipvscmc2 = new TProfile_Ref(name.c_str(), name.c_str(), chips, -0.5, chips-0.5);
  ((TProfile_Ref*)fHchipvscmc2)->SetReference(fReferenceDirectory);
  AddHistogram(fHchipvscmc2);
  fHchipvscmc2->SetStats(false);
  fHchipvscmc2->GetXaxis()->SetTitle("APV chip");
  fHchipvscmc2->GetYaxis()->SetTitle("Common mode correction");
  
  if (fExpertHistos) {
    // common mode correction value for each chip
    for (int i=0; i<chips; i++) {
      char tmp[5];
      sprintf(tmp, "%02i", i);
      name = fName + "_cmc2_chip" + tmp;
      TH1D* hChipCmc2 = new TH1D(name.c_str(), name.c_str(), 201, -100.5, 100.5);
      AddHistogram(hChipCmc2);
      vHchipcmc2s.push_back(hChipCmc2);
    }
  }

  // Single channel occupancies
  name = fName + "_occupancies";
  fHoccup = new TH1F_Ref(name.c_str(), name.c_str(), fVch->GetNbins(), fVch->GetMin(), fVch->GetMax(),fRateCounter,true);
  ((TH1F_Ref*)fHoccup)->SetReference(fReferenceDirectory);
  AddHistogram(fHoccup);
  fHoccup->SetStats(false);

  if (fExpertHistos) {
    // two histograms for the mean amplitude per pad
    // first the sum
    name = fName + "_sumAmp";
    fHsumAmp = new TH1F(name.c_str(), name.c_str(), fVch->GetNbins(), fVch->GetMin(), fVch->GetMax());
    // second the mean (sum divided by channels)
    name = fName + "_meanAmp";
    fHmeanAmp = new TH1F(name.c_str(), name.c_str(), fVch->GetNbins(), fVch->GetMin(), fVch->GetMax());
    AddHistogram(fHmeanAmp);
    fHmeanAmp->SetStats(false);
  }

  // cluster variable
  if ( fName.substr(0, 2)=="GP" ) // PixelGEM has two hemispheres, they are combined for clusters
    fVcPos  = AddVariable("Clu_Profile", fNchan/2, -0.5, (fNchan/2)-0.5, fNchan*fMAX_MULT);
  else
    fVcPos  = AddVariable("Clu_Profile",   fNchan, -0.5,     fNchan-0.5, fNchan*fMAX_MULT);
  fVcA2   = AddVariable("Clu_Amp2", 1000,    0., 4000., fNchan*fMAX_MULT);
  fVcSize = AddVariable("Clu_Size",   30,  -0.5,  29.5, fNchan*fMAX_MULT);
  fVcTime = AddVariable("Clu_Time",  100, -100.,  100., fNchan*fMAX_MULT);

  // Book histograms for GEM cluster data
  // Number of clusters
  std::string chitsname = fName + "Clu_Multiplicity";
  fHcHits=new TH1F(chitsname.c_str(),chitsname.c_str(),70,-0.5,69.5);
  AddHistogram(fHcHits);

  // Cluster profile
  std::string cchname = fName + fVcPos->GetName();
  fHcPos=new TH1F(cchname.c_str(), cchname.c_str(), fVcPos->GetNbins(), fVcPos->GetMin(), fVcPos->GetMax());
  AddHistogram(fHcPos);

  // Cluster amplitude 1
  name = fName + "Clu_Amp0";
  fHcA0=new TH1F(name.c_str(), name.c_str(), 1000, 0., 4000.);
  AddHistogram(fHcA0);

  // Cluster amplitude 2
  name = fName + "Clu_Amp1";
  fHcA1=new TH1F(name.c_str(), name.c_str(), 1000, 0., 4000.);
  AddHistogram(fHcA1);

  // Cluster amplitude 3
  std::string ca2name = fName + fVcA2->GetName();
  fHcA2=new TH1F(ca2name.c_str(), ca2name.c_str(), fVcA2->GetNbins(), fVcA2->GetMin(), fVcA2->GetMax());
  AddHistogram(fHcA2);

  //Cluster size
  std::string csname = fName + fVcSize->GetName();
  fHcSize=new TH1F(csname.c_str(), csname.c_str(), fVcSize->GetNbins(), fVcSize->GetMin(), fVcSize->GetMax());
  AddHistogram(fHcSize);

  // Cluster amplitude ratio ("banana")
  name = fName + "Clu_AmpRatio";
  fHcAmpRatio = new TH2F(name.c_str(), name.c_str(), 100, -0.01, 1.99, 100, -0.01, 1.99);
  AddHistogram(fHcAmpRatio);
  fHcAmpRatio->SetStats(false);
  fHcAmpRatio->GetXaxis()->SetTitle("A1/A2");
  fHcAmpRatio->GetYaxis()->SetTitle("A0/A2");
  fHcAmpRatio->SetOption("COLZ");

  // Cluster time
  std::string ctname = fName + fVcTime->GetName();
  fHcTime = new TH1F(ctname.c_str(), ctname.c_str(), fVcTime->GetNbins(), fVcTime->GetMin(), fVcTime->GetMax());
  AddHistogram(fHcTime);
  fHcTime->GetXaxis()->SetTitle("time / ns");

  if (fExpertHistos) {
    // two histograms for the mean amplitude per cluster
    // first the sum
    name = fName + "Clu_SumAmp";
    fHcSumAmp = new TH1F(name.c_str(), name.c_str(), fVcPos->GetNbins(), fVcPos->GetMin(), fVcPos->GetMax());
    // second the mean (sum divided by channels)
    name = fName + "Clu_MeanAmp";
    fHcMeanAmp = new TH1F(name.c_str(), name.c_str(), fVcPos->GetNbins(), fVcPos->GetMin(), fVcPos->GetMax());
    AddHistogram(fHcMeanAmp);
    fHcMeanAmp->SetStats(false);
  }

  // histograms after time cut
  // Number of clusters after time cut
  name = fName + "TClu_Multiplicity";
  fHcHitsT=new TH1F(name.c_str(),name.c_str(),70,-0.5,69.5);
  AddHistogram(fHcHitsT);

  // Cluster profile after time cut
  name = fName + "T" + fVcPos->GetName();
  fHcPosT=new TH1F(name.c_str(),name.c_str(), fVcPos->GetNbins(), fVcPos->GetMin(), fVcPos->GetMax());
  AddHistogram(fHcPosT);

  // Cluster amplitude 3 after time cut
  name = fName + "T" + fVcA2->GetName();
  fHcA2T=new TH1F(name.c_str(), name.c_str(), fVcA2->GetNbins(), fVcA2->GetMin(), fVcA2->GetMax());
  AddHistogram(fHcA2T);

  //Cluster size after time cut
  name = fName+"TClu_Size";
  fHcSizeT=new TH1F(name.c_str(),  name.c_str(), fVcSize->GetNbins(), fVcSize->GetMin(), fVcSize->GetMax());
  AddHistogram(fHcSizeT);

  // Cluster amplitude ratio ("banana") after time cut
  name = fName + "TClu_AmpRatio";
  fHcAmpRatioT = new TH2F(name.c_str(), name.c_str(), 100, -0.05, 1.95, 100, -0.05, 1.95);
  AddHistogram(fHcAmpRatioT);
  fHcAmpRatioT->SetStats(false);
  fHcAmpRatioT->GetXaxis()->SetTitle("A1/A2");
  fHcAmpRatioT->GetYaxis()->SetTitle("A0/A2");
  fHcAmpRatioT->SetOption("COLZ");

  // Cluster time after time cut
  name = fName + "T" + fVcTime->GetName();
  fHcTimeT = new TH1F(name.c_str(), name.c_str(), fVcTime->GetNbins(), fVcTime->GetMin(), fVcTime->GetMax());
  AddHistogram(fHcTimeT);
  fHcTimeT->GetXaxis()->SetTitle("time / ns");

  if (fExpertHistos) {
    // two histograms for the mean amplitude per cluster
    // first the sum
    name = fName + "TClu_SumAmp";
    fHcSumAmpT = new TH1F(name.c_str(), name.c_str(), fVcPos->GetNbins(), fVcPos->GetMin(), fVcPos->GetMax());
    // second the mean (sum divided by channels)
    name = fName + "TClu_MeanAmp";
    fHcMeanAmpT = new TH1F(name.c_str(), name.c_str(), fVcPos->GetNbins(), fVcPos->GetMin(), fVcPos->GetMax());
    AddHistogram(fHcMeanAmpT);
    fHcMeanAmpT->SetStats(false);
  }

  // Create tree leaves
  if(tree) {
    std::string chitsleavlist = chitsname + "/I";
    std::string cchleavlist   = cchname + "[" + chitsname +"]/F";
    std::string ca2leavlist   = ca2name + "[" + chitsname +"]/F";
    std::string csleavlist    = csname  + "[" + chitsname +"]/F";
    std::string ctleavlist    = ctname  + "[" + chitsname +"]/F";

    tree->Branch(chitsname.c_str(), &fNclustKept,          chitsleavlist.c_str(), 32000);
    tree->Branch(ca2name.c_str()  ,  fVcA2->GetValues(),   ca2leavlist.c_str(),   32000);
    tree->Branch(cchname.c_str()  ,  fVcPos->GetValues(),  cchleavlist.c_str(),   32000);
    tree->Branch(csname.c_str()   ,  fVcSize->GetValues(), csleavlist.c_str(),    32000);
    tree->Branch(ctname.c_str()   ,  fVcTime->GetValues(), ctleavlist.c_str(),    32000);
  }
}

void PlaneGEM::Reset() {
  PlaneAPV::Reset();
  fDigits.clear();

  if (fGeom) {
    GeomPlaneGEM *GEMGeom = dynamic_cast<GeomPlaneGEM*>(fGeom);
    if (GEMGeom) GEMGeom->ResetPlane();
  }

  fNclustKept = 0;
}

void PlaneGEM::EndEvent(const CS::DaqEvent &event) {
  // need geometry object for calculation of hits and cuts and ...
  GeomPlaneGEM *GEMGeom(0);
  if (fGeom) GEMGeom = dynamic_cast<GeomPlaneGEM*>(fGeom);

  if (thr_flag) TThread::Lock();

  // Get TCS phase
  try {
    TCSphase = event.GetTT().GetPhaseTCS();
  } catch (std::logic_error& e) {
    std::cerr << "PlaneGEM::EndEvent: TCS phase not available: " << e.what() << std::endl;
    TCSphase = 0;
  }

  // save number of hits
  fNhitsKept = 0;

  // only draw the common mode noise correction per chip once per event
  int chips = fNchan / 128;
  bool* filledCMC = new bool[chips];
  for (int i=0; i<chips; i++) filledCMC[i] = false;

  // loop over digits
  for (std::map<int,CDigit1>::const_iterator id = fDigits.begin(); id!=fDigits.end(); id++) {
    const CDigit1& digit = id->second;

    int   channel = digit.ch;
    int   hem     = (int) digit.dt[9];
    int   chip    = (int) digit.dt[10];

    double a0      = digit.dt[0];
    double a1      = digit.dt[1];
    double a2      = digit.dt[2];
    double a02     = (a2>0.) ? a0/a2 : 0.;
    double a12     = (a2>0.) ? a1/a2 : 0.;
    double cmc2    = digit.dt[7];

    if (hem) // for GM hem is always 0, for GP either +/-1, so this is only done for GP
      channel = channel + (1-hem)*(fNchan/4);
      
    // store data into variables
    fVch->Store(channel);
    fVa0->Store(a0);
    fVa1->Store(a1);
    fVa2->Store(a2);

    // save into histogram
    fHch->Fill(channel);
    if (fExpertHistos)
      fHsumAmp->Fill(channel, a2);
    fHa0->Fill(a0);
    fHa1->Fill(a1);
    fHa2->Fill(a2);
    fHchvsa2->Fill(channel, a2);

    // fill CMC per chip, only once per event
    if (channel/128 < chips) {
      if (filledCMC[channel/128]== false) {
        filledCMC[channel/128] = true;
#if USE_DATABASE==1
        if (maxBaselines.count(chip)>0) {
          fHchipvscmc2->Fill(channel/128, cmc2-(1024.-maxBaselines[chip]));
          if (fExpertHistos)
            vHchipcmc2s[channel/128]->Fill(cmc2-(1024.-maxBaselines[chip]));
        } else {
#endif
          fHchipvscmc2->Fill(channel/128, cmc2);
          if (fExpertHistos)
            vHchipcmc2s[channel/128]->Fill(cmc2);
#if USE_DATABASE==1
          if (!fPrintedMissCalib) {
            std::cerr << "PlaneGEM::EndEvent: No mean baseline for APV chips of plane " << fName << "..." << std::endl;
            fPrintedMissCalib=true;
          }
        }
#endif
      }
    } else
      std::cerr << "Strange channel number..." << std::endl;

    fHampRatio->Fill(a12,a02);

    double time;
    if ( GEMGeom && GEMGeom->CalcTime(digit, time) ) { // is there a valid time, variable has still to be filled
      fHt->Fill(TCSphase - time);
      fVt->Store(TCSphase - time);
    } else
      fVt->Store(101.);

    fNhitsKept++;
  }
  fHhits->Fill(fNhitsKept);

  if (!(fRateCounter%25)) {
    fHoccup->Reset("ICE");
    fHoccup->Add(fHch, 1/double(fRateCounter));

    if (fExpertHistos)
      fHmeanAmp->Divide(fHsumAmp, fHch);
  }

  // cleanup check of CMC
  delete [] filledCMC;

  if (thr_flag) TThread::UnLock();
}

void PlaneGEM::Clusterize() {
  // No geometry
  if (!fGeom) return;

  if (thr_flag) TThread::Lock();

  // Do clustering
  std::set<CCluster1>& clusters = fGeom->Clusterize(fDigits);

  // save number of clusters (before and after time cut)
  unsigned int NrClu(0), NrTClu(0);

  // this is were the real important stuff happens - cluster handling
  for (std::set<CCluster1>::iterator cluit=clusters.begin(); cluit!=clusters.end(); cluit++) {
    assert( cluit->size ); // otherwise he have forgotten to remove something above

    // get data from clusters
    double strip = cluit->dt[3];
    double amp0  = cluit->dt[0];
    double amp1  = cluit->dt[1];
    double amp2  = cluit->dt[2];
    int size     = cluit->size;

    // store in variables
    fVcPos ->Store(cluit->pos);
    fVcA2  ->Store(amp2);
    fVcSize->Store(size);

    // fill histograms
    fHcPos     ->Fill(strip);
    if (fExpertHistos)
      fHcSumAmp->Fill(strip, amp2);
    fHcA0      ->Fill(amp0);
    fHcA1      ->Fill(amp1);
    fHcA2      ->Fill(amp2);
    fHcSize    ->Fill(size);
    fHcAmpRatio->Fill(amp1/amp2, amp0/amp2);

    // clusters within tight time window (cut on amplitude ratio)
    const float a02min=0.1;
    const float a02max=0.6;
    const float a12min=0.5;
    const float a12max=0.95;

    float a02=amp0/amp2;
    float a12=amp1/amp2;

    if ( a02>=a02min && a02<=a02max && a12>=a12min && a12<=a12max ) {
      fHcPosT     ->Fill(strip);
      if (fExpertHistos)
        fHcSumAmpT->Fill(strip, amp2);
      fHcA2T      ->Fill(amp2);
      fHcSizeT    ->Fill(size);
      fHcAmpRatioT->Fill(amp1/amp2, amp0/amp2);

      // increase number of clusters after time cut
      NrTClu++;
    }

    // time information
    if ( cluit->dt.size() > 4 ) { // there is a valid time
      fHcTime->Fill(TCSphase - cluit->dt[4]);
      if ( a02>=a02min && a02<=a02max && a12>=a12min && a12<=a12max )
        fHcTimeT->Fill(TCSphase - cluit->dt[4]);

      fVcTime->Store(TCSphase - cluit->dt[4]);
    } else // we still have to fill the variable which is saved into the tree, histogram is not filled
      fVcTime->Store(101.);

    // increase number of clusters
    NrClu++;
  }
  fHcHits ->Fill(NrClu);
  fHcHitsT->Fill(NrTClu);
  fNclustKept = NrClu;

  if (!(fRateCounter%25)) {
    if (fExpertHistos) {
      fHcMeanAmp ->Divide(fHcSumAmp,  fHcPos);
      fHcMeanAmpT->Divide(fHcSumAmpT, fHcPosT);
    }
  }

  if (thr_flag) TThread::UnLock();
}

#if USE_DATABASE == 1
void PlaneGEM::ReadCalib(const tm &t)
{
  //  Read calibration data for each strip
  try{
    ReadFromDataBase(calib_data, t);
    // std::cout<<"PlaneGEM::ReadCalib() ==> "
    // <<this->GetName()<<" calibrations are found !"<<std::endl;

    maxBaselines.clear();

    unsigned int calib_size=0;

    std::map<unsigned int, PlaneGEM::APVcalib>::iterator it;
    for (it=calib_data.begin(); it!=calib_data.end(); it++) {
      calib_size += it->second.channel.size();

      float max=   0.;
      std::vector<APVchannel>::iterator itchan;
      for (itchan=it->second.channel.begin(); itchan!=it->second.channel.end(); itchan++) {
        if (max<itchan->ped) max = itchan->ped;
      }
      if (it->second.channel.size()>0) {
        maxBaselines[it->second.channel.begin()->chipId] = max;
      }
    }

    if(calib_size != (unsigned) fNchan) {
      std::cerr<<GetName()<<": Size of Calibration File is not correct ! Should be : "
               <<fNchan<<" Is "<<calib_size<<" "
               <<t.tm_mday<<"."<<t.tm_mon+1<<"."<<t.tm_year+1900<<" "
               <<t.tm_hour<<":"<<t.tm_min<<":"<<t.tm_sec<<std::endl;
      calib_data.clear();
    }
    else
      fUseCalib = true;
  }

  catch(CS::Exception& e) {
    std::cerr<<e.what()<<std::endl;
  }
  catch(const std::exception &e) {
    std::cerr<<e.what()<<std::endl;
  }
  catch(...) {
    std::cout<<"PlaneGEM::ReadCalib() ==> "<<GetName()
	<<" calibrations, valid for ";
    std::cout<<t.tm_mday<<"."<<t.tm_mon+1<<"."<<t.tm_year+1900<<" "
	<<t.tm_hour<<":"<<t.tm_min<<":"<<t.tm_sec
	<<", not found in DB"<<std::endl;
  }

  //clear current time calibration
  calib_time.clear();

  // Read calibration data for conversion of amplitude->time for each plane
  if (fDataBase==NULL)
    throw CS::Exception("PlaneGEM::ReadCalib():  data base is not opened.");
  try{
    //    std::cout << "Reading timing calibrations for " << GetName() << std::endl;
    struct tm tt(t);
    CDB::Time tp(mktime(&tt),0);
    std::string strdata("");
    fDataBase->read(fName,strdata,tp,"timing");
    if (strdata == "") {
      std::cerr << "PlaneGEM::ReadCalib() "<<GetName()<<" : no timing calibration file found"<<std::endl;
      return;
    }
    std::istringstream istrdata(strdata.c_str());
    //    std::cout << strdata << std::endl;
    istrdata >> calib_time;
  }
  catch(CS::Exception& e) {
    std::cout<<"rethrowing"<<std::endl;
    throw;
  }

  if (calib_time.size()==10) { // old calibration format, used including 2006
    // Just debug printouts
    if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_GM")!=0) {
      const float& r1_min = this->calib_time[0];
      const float& r1_max = this->calib_time[1];
      const float& t0_1 = this->calib_time[2];
      const float& sl_1 = this->calib_time[3];
      const float& tc_1 = this->calib_time[4];
      const float& r2_min = this->calib_time[5];
      const float& r2_max = this->calib_time[6];
      const float& t0_2 = this->calib_time[7];
      const float& sl_2 = this->calib_time[8];
      const float& tc_2 = this->calib_time[9];
      std::cout<<std::endl<<"Timing calibrations read in for "<<GetName()<<std::endl;
      std::cout<<"r1_min = "<<r1_min<<std::endl;
      std::cout<<"r1_max = "<<r1_max<<std::endl;
      std::cout<<"r2_min = "<<r2_min<<std::endl;
      std::cout<<"r2_max = "<<r2_max<<std::endl;
      std::cout<<"t0_1   = "<<t0_1<<std::endl;
      std::cout<<"t0_2   = "<<t0_2<<std::endl;
      std::cout<<"sl_1   = "<<sl_1<<std::endl;
      std::cout<<"sl_2   = "<<sl_2<<std::endl;
      std::cout<<"tc_1   = "<<tc_1<<std::endl;
      std::cout<<"tc_2   = "<<tc_2<<std::endl;
      std::cout<<std::endl;
    }
  } else if (calib_time.size() == 24) { // new calibration format, used from 2007 on
    // Just debug printouts
    if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_GM")!=0) {
      std::cout<<std::endl<<"Timing calibrations read in for "<<GetName()<<std::endl;
      // omitting the output of the covariance matrix which is used for error calculation
      std::cout<<"Ratio 02: a0="<<calib_time[ 0]<<", t0="<<calib_time[ 1]<<", r0="<<calib_time[ 2]<<std::endl;
      std::cout<<"Ratio 12: a0="<<calib_time[12]<<", t0="<<calib_time[13]<<", r0="<<calib_time[14]<<std::endl;
      std::cout<<std::endl;
    }
  } else {
    calib_time.clear(); // remove whatever has been loaded

    std::cout<<"PlaneGEM::ReadCalib() ==> "<<GetName()<<" timing calibrations, valid for ";
    std::cout<<t.tm_mday<<"."<<t.tm_mon+1<<"."<<t.tm_year+1900<<" "
             <<t.tm_hour<<":"<<t.tm_min<<":"<<t.tm_sec
             <<", not found in DB"<<std::endl;
  }
}

istream& operator>>(istream& in, std::map<unsigned int, PlaneGEM::APVcalib> &c) {
  char firstChar;
  char chipname[10];
  int ldcId, srcId, adcId, chipId;
  bool ignore=true;
  c.clear();

  while (in) {
    in >> firstChar;

    if (in.eof()) // firstchar could not be read, as the file end was reached
      continue;

    if (firstChar == '#') {
      in >> chipname >> ldcId >> srcId >> adcId  >> chipId;
      if ( c.count(chipId)==0 ) { // first entry found for this chip ID
        c.insert( std::pair<unsigned int, PlaneGEM::APVcalib>(chipId, PlaneGEM::APVcalib()) );
        ignore=false;
      } else { // this seems to be a broken pedestal file, two times the same APV ID... Ignore input till next chip
        std::cerr << "Broken pedestal file found for GEM. The same APV ID was found several times." << std::endl;
        ignore=true;
      }
    } else {
      int f = atoi(&firstChar);
      float p, s, cp, cs;
      in >> p >> s >> cp >> cs;
      if (!ignore)
        c[chipId].channel.push_back(PlaneGEM::APVchannel(f, p, s, cp, cs, ldcId, srcId, adcId, chipId));
    }
  }

  return in;
}

#endif //USE_DATABASE

void PlaneGEM::TextOutput(ostream& out) {

  // Write maximum of Landau of cluster amplitude 3 to file
  if(float gain = Gain())
    out<<"\tgain "<<gain<<std::endl;
  else
    out<<"\tgain ***"<<std::endl;
}

float PlaneGEM::Gain() {
  static const int statcut = 1000;
  if(fHcA2T->GetEntries() < statcut) return 0;

  TF1 *lan = new TF1("lan","landau",40.,4000.);
  fHcA2T->Fit(lan,"RQ0");
  if (lan->GetChisquare()==0.) // fit gone bad
    return 0;

  return lan->GetParameter(1);
}

