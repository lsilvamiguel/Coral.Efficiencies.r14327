#include "PlaneScifiJ.h"
#include "TriggerTime.h"
#include "GeomPlaneScifiJ.h"

ClassImp(PlaneScifiJ);

void PlaneScifiJ::Init(TTree* tree) {
  Plane1V::Init(tree);

  fFillfDigits = false;

  // Cluster histograms
  float posmax = 50;

  if(fGeom){
    //    GeomPlaneScifiJ* geompl = (GeomPlaneScifiJ*) fGeom;
    posmax = fGeom->GetHalfSize();
  }

  fVcch = AddVariable("_clusterPosition",100,-posmax,posmax,fNchan*fMAX_MULT);
  fVcs = AddVariable("_clusterSize",30, -0.5, 29.5, fNchan*fMAX_MULT);
  fVtc = AddVariable("_tcorr",100,-1000,1000,fNchan*fMAX_MULT);

  std::string chitsname = fName + "_clusterMultiplicity";

  std::string tcname = fName+"_tcorr";
  fHtc=new TH1F(tcname.c_str(),tcname.c_str(), 400, -80, 80);
  AddHistogram(fHtc);
  fHtc->GetXaxis()->SetTitle("corrected time (ns)");
  std::string tcleavlist = tcname + "[" + chitsname +"]/F";

  std::string tcvcname = fName+"_tcorr_vs_ch";
  fHtcVSch=new TH2F(tcvcname.c_str(),tcvcname.c_str(), 200, -20, 20,
		    fVch->GetNbins(), fVch->GetMin(), fVch->GetMax());
  fHtcVSch->SetOption("col");
  AddHistogram(fHtcVSch);

  // clusters ----------------------------------------------------------

  //multiplicity within the time cut on the whole detector
  //std::string chitsname = fName + "_clusterMultiplicity";
  fHchits=new TH1F(chitsname.c_str(),chitsname.c_str(),70,-0.5,69.5);
  AddHistogram(fHchits);
  std::string chitsleavlist = chitsname + "/I";

  //branch : cluster map
  std::string cchname = fName + fVcch->GetName();
  fHcch=new TH1F(cchname.c_str(),cchname.c_str(),
                 fVcch->GetNbins(), fVcch->GetMin(), fVcch->GetMax());
  AddHistogram(fHcch);
  std::string cchleavlist = cchname + "[" + chitsname +"]/F";

  //Cluster size
  std::string csname = fName+"_clusterSize";
  fHcs=new TH1F(csname.c_str(),csname.c_str(),
                fVcs->GetNbins(), fVcs->GetMin(), fVcs->GetMax());
  AddHistogram(fHcs);
  std::string csleavlist = csname + "[" + chitsname +"]/F";

  fHchVsTis=new TH2F((fName+"_ch_vs_tis").c_str(), (fName+" channel versus time in spill").c_str(),
                     fVch->GetNbins(), fVch->GetMin(), fVch->GetMax(),
                     100, 0., 15.);
  fHchVsTis->SetOption("COLZ");
  AddHistogram(fHchVsTis);

}




void PlaneScifiJ::EndEvent(const CS::DaqEvent &event) {

  if (thr_flag) TThread::Lock();

  fDigits.clear();

//   const float reftime = fIsHR ? event.GetTT().GetTimeHigh():
//     event.GetTT().GetTimeNorm();
//   const float factor = fIsHR ? CS::ChipF1::GetUnitHigh() :
//     CS::ChipF1::GetUnitNorm();

//   for (int i=0; i<fNhits; i++) {
//
//     double time = CS::ChipF1::TimeDifference(fTime[i],reftime);
//     int channel=fChannel[i];
  typedef std::list<CS::Chip::Digit*>::iterator lDIter;
  for (lDIter ii = lDigits.begin(); ii != lDigits.end(); ii++) {
    CS::ChipF1::Digit* iii = dynamic_cast<CS::ChipF1::Digit*> (*ii);
    if (!iii) {
      std::cerr<<"PlaneScifiJ::EndEvent: a digit is not a F1 one, strange...\n";
      continue;
    }
    register float factor = iii->GetTimeUnit();
    register double time = (double) (iii->GetTimeDecoded() / factor);
    register int channel = iii->GetChannel();


    float timeT0 = -8000;
    if (fUseCalib && channel<(int)calib_data.size())
      timeT0 = calib_data[channel].t0;
    double tcorr = (time - timeT0)*factor;

    std::vector<double> data;
    data.push_back(tcorr);
    fDigits.insert(make_pair(channel,CDigit1(channel,data)));	

    if( fVch->Test(channel) ) {

      fHtc->Fill(tcorr); //corrected time
      fHtcVSch->Fill(tcorr, channel); //corrected time

      if( fVt->Test(time) ) {
	fVch->Store(channel);
	fVt->Store(time);
	fHch->Fill(channel);
	fHt->Fill(time);
	fHtvsch->Fill(channel,time);
	fHchVsTis->Fill(channel, event.GetTT().GetTimeInSpill());
     	fNhitsKept++;
      }

      if( fOnTrigTVarID->Test(time) )  {
        fOnTrigTVarID->Store(time);
	fHonTT->Fill(time);        // on trigger time distribution
        fHonTP->Fill(channel);  // on trigger profile
        fNhitsKeptOnTrig++;
      }

      if( fOffTrigTVarID->Test(time) )  {
        fOffTrigTVarID->Store(time);
        fHoffTT->Fill(time);        // off trigger time distribution
        fHoffTP->Fill(channel);  // off trigger profile
        fNhitsKeptOffTrig++;
      }
    }
  }
  fHhits->Fill(fNhitsKept);

  if((fRateCounter%Plane1V::fRATE_UPDATE)==0) {

    // corrected on-trigger profile
    if( fOffTrigTVarID->GetMax() != fOffTrigTVarID->GetMin() ) {
      float scale = (float) ( fOnTrigTVarID->GetMax() - fOnTrigTVarID->GetMin() )/( fOffTrigTVarID->GetMax() - fOffTrigTVarID->GetMin() );
      for(register int i=1; i<=fHonTP->GetNbinsX(); i++)
	fHoncorrTP->SetBinContent(i, (float) fHonTP->GetBinContent(i) - scale*fHoffTP->GetBinContent(i) );
    }

    // rates histogram
    float timewin = ( fOffTrigTVarID->GetMax() - fOffTrigTVarID->GetMin() )*fF1_TICK*fRateCounter;
    if( timewin > 0 )
      for(register int i=1; i<=fHonTP->GetNbinsX(); i++) {
	fHrates->SetBinContent(i,(float) fHoffTP->GetBinContent(i)/timewin);
	fHrates->SetBinError(i,(float) fHoffTP->GetBinError(i)/timewin);
      }
  }
  if (thr_flag) TThread::UnLock();
}



void PlaneScifiJ::Clusterize() {
  if (!fGeom) return;
  std::set<CCluster1>& clusters = fGeom->Clusterize(fDigits);

  if (thr_flag) TThread::Lock();
  int nclusters=0;
  typedef std::set<CCluster1>::iterator IC;
  for (IC i = clusters.begin(); i!=clusters.end(); i++) {

    float cpos = i->pos;
    float size = i->size;

    if (fVcch->Test(cpos) &&
        fVcs->Test(size)) {

      fHcs->Fill(size);
      fHcch->Fill(cpos);
      fVcch->Store(cpos);
      nclusters++;
    }
    else {
      //update set of clusters
//      clusters.erase(i);
    }
  }
  fNclustKept = nclusters;
  fHchits->Fill(nclusters);

  if (thr_flag) TThread::UnLock();
}


#if USE_DATABASE == 1
void PlaneScifiJ::ReadCalib(const tm &t)
{
  // read-in corresponding calibration constants
  try{
    ReadFromDataBase(calib_data,t);
    //   std::cout<<"PlaneScifiJ::ReadCalib() ==> "
    //	<<this->GetName()<<" calibrations are found !"<<std::endl;

    if(calib_data.size() <  (unsigned) fNchan) {
      //    std::cerr<<"Size of Calibration File is not correct ! Should be : "
      //	  <<fNchan<<" Is "<<calib_data.size()<<" "
      //	  <<t.tm_mday<<"."<<t.tm_mon+1<<"."<<t.tm_year+1900<<" "
      //  <<t.tm_hour<<":"<<t.tm_min<<":"<<t.tm_sec<<std::endl;
    }
    else
      fUseCalib = true;

    //is not ok!!!!!!!
    fUseCalib=true;
  }

  catch(CS::Exception& e) {
    std::cerr<<e.what()<<std::endl;
  }
  catch(const std::exception &e) {
    std::cerr<<e.what()<<std::endl;
  }
  catch(...) {
    std::cout<<"PlaneScifiJ::ReadCalib() ==> "<<GetName()
	<<" calibrations, valid for ";
    std::cout<<t.tm_mday<<"."<<t.tm_mon+1<<"."<<t.tm_year+1900<<" "
	<<t.tm_hour<<":"<<t.tm_min<<":"<<t.tm_sec
	<<", not found in DB"<<std::endl;
  }

}
#endif //USE_DATABASE






