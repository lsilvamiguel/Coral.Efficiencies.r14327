#include "PlanePGEM.h"
#include "GeomPlanePGEM.h"
#include "ChipAPV.h"

#include <cmath>
#include <stdexcept>

ClassImp(PlanePGEM);

PlanePGEM::PlanePGEM(const char *detname,int nchan, int center, int width, bool pixel)
  : PlaneGEM(detname,nchan,center,width,pixel) {
  INeedGeom();
}

void PlanePGEM::Init(TTree* tree) {
  if (fPixel) { // pixel plane
    // override the default settings for hittimes from PlaneAPV before
    // calling the init method, as there the histograms are created
    fVt->SetRange(-100., 100.);
    fVt->SetNbins(200);

    PlaneAPV::Init(tree);

    // Geometry and calibrations
    if ( fGeom ) {
      GeomPlanePGEM* geompl = (GeomPlanePGEM*) fGeom;
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

    // hit pads
    name = fName + "_HitMap";
    fHhitMapPix = new TH2F(name.c_str(), name.c_str(), 32, -0.5, 31.5, 32, -0.5, 31.5);
    AddHistogram(fHhitMapPix);
    fHhitMapPix->SetStats(false);
    fHhitMapPix->SetOption("COLZ");

    // pad occupancies
    name = fName + "_Occupancy";
    fHoccupPix = new TH2F(name.c_str(), name.c_str(), 32, -0.5, 31.5, 32, -0.5, 31.5);
    AddHistogram(fHoccupPix);
    fHoccupPix->SetStats(false);
    fHoccupPix->SetOption("COLZ");

    if (fExpertHistos) {
      // two histograms for the mean amplitude per pad
      // first the sum
      name = fName + "_HitMapSumAmp";
      fHhitMapSumAmp = new TH2F(name.c_str(), name.c_str(), 32, -0.5, 31.5, 32, -0.5, 31.5);
      // second the mean (sum divided by HitMap)
      name = fName + "_HitMapMeanAmp";
      fHhitMapMeanAmp = new TH2F(name.c_str(), name.c_str(), 32, -0.5, 31.5, 32, -0.5, 31.5);
      AddHistogram(fHhitMapMeanAmp);
      fHhitMapMeanAmp->SetStats(false);
      fHhitMapMeanAmp->SetOption("COLZ");
    }

    // cluster variables
    fVcPosPixX = AddVariable("Clu_PosPixX", 32,  -0.5,  31.5, fNchan*fMAX_MULT);
    fVcPosPixY = AddVariable("Clu_PosPixY", 32,  -0.5,  31.5, fNchan*fMAX_MULT);
    fVcA2      = AddVariable("Clu_Amp2",  1000,    0., 4000., fNchan*fMAX_MULT);
    fVcSize    = AddVariable("Clu_Size",    30,  -0.5,  29.5, fNchan*fMAX_MULT);
    fVcTime    = AddVariable("Clu_Time",   100, -100.,  100., fNchan*fMAX_MULT);

    // cluster histograms
    // Number of clusters
    std::string chitsname = fName + "Clu_Multiplicity";
    fHcHits=new TH1F(chitsname.c_str(),chitsname.c_str(),70,-0.5,69.5);
    AddHistogram(fHcHits);

    // Cluster profile
    std::string cchname = fName + "Clu_PosPix";
    fHcPosPix=new TH2F(cchname.c_str(), cchname.c_str(), fVcPosPixX->GetNbins(), fVcPosPixX->GetMin(), fVcPosPixX->GetMax(),
                                                         fVcPosPixY->GetNbins(), fVcPosPixY->GetMin(), fVcPosPixY->GetMax());
    fHcPosPix->SetStats(false);
    fHcPosPix->SetOption("COLZ");
    AddHistogram(fHcPosPix);

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

    // Cluster size
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
      fHcSumAmpPix = new TH2F(name.c_str(), name.c_str(), 32, -0.5, 31.5, 32, -0.5, 31.5);
      // second the mean (sum divided by PosPix)
      name = fName + "Clu_MeanAmp";
      fHcMeanAmpPix = new TH2F(name.c_str(), name.c_str(), 32, -0.5, 31.5, 32, -0.5, 31.5);
      AddHistogram(fHcMeanAmpPix);
      fHcMeanAmpPix->SetStats(false);
      fHcMeanAmpPix->SetOption("COLZ");
    }

    // histograms after time cut
    // Number of clusters after time cut
    name = fName + "TClu_Multiplicity";
    fHcHitsT=new TH1F(name.c_str(),name.c_str(),70,-0.5,69.5);
    AddHistogram(fHcHitsT);

    // Cluster profile after time cut
    name = fName + "TClu_PosPix";
    fHcPosPixT=new TH2F(name.c_str(), name.c_str(), fVcPosPixX->GetNbins(), fVcPosPixX->GetMin(), fVcPosPixX->GetMax(),
                                                    fVcPosPixY->GetNbins(), fVcPosPixY->GetMin(), fVcPosPixY->GetMax());
    fHcPosPixT->SetStats(false);
    fHcPosPixT->SetOption("COLZ");
    AddHistogram(fHcPosPixT);

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
      fHcSumAmpPixT = new TH2F(name.c_str(), name.c_str(), 32, -0.5, 31.5, 32, -0.5, 31.5);
      // second the mean (sum divided by PosPix)
      name = fName + "TClu_MeanAmp";
      fHcMeanAmpPixT = new TH2F(name.c_str(), name.c_str(), 32, -0.5, 31.5, 32, -0.5, 31.5);
      AddHistogram(fHcMeanAmpPixT);
      fHcMeanAmpPixT->SetStats(false);
      fHcMeanAmpPixT->SetOption("COLZ");
    }

    // create tree leaves
    if (tree) {
      std::string chitsleavlist = chitsname + "/I";
      std::string cpxleavlist   = cchname + "X[" + chitsname +"]/F";
      std::string cpyleavlist   = cchname + "Y[" + chitsname +"]/F";
      std::string ca2leavlist   = ca2name +  "[" + chitsname +"]/F";
      std::string csleavlist    = csname  +  "[" + chitsname +"]/F";
      std::string ctleavlist    = ctname  +  "[" + chitsname +"]/F";

      tree->Branch(chitsname    .c_str(), &fNclustKept,             chitsleavlist.c_str(), 32000);
      tree->Branch(ca2name      .c_str(),  fVcA2->GetValues(),      ca2leavlist.c_str(),   32000);
      tree->Branch((cchname+"X").c_str(),  fVcPosPixX->GetValues(), cpxleavlist.c_str(),   32000);
      tree->Branch((cchname+"Y").c_str(),  fVcPosPixY->GetValues(), cpyleavlist.c_str(),   32000);
      tree->Branch(csname       .c_str(),  fVcSize->GetValues(),    csleavlist.c_str(),    32000);
      tree->Branch(ctname       .c_str(),  fVcTime->GetValues(),    ctleavlist.c_str(),    32000);
    }
  } else { // strip plane of PixelGEM detector
    PlaneGEM::Init(tree);
  }
}

void PlanePGEM::EndEvent(const CS::DaqEvent &event) {
  if (fPixel) { // pixel plane
    // need geometry object for calculation of hits and cuts and ...
    GeomPlanePGEM *pixelGEMGeom(0);
    if (fGeom) pixelGEMGeom = dynamic_cast<GeomPlanePGEM*>(fGeom);

    if (thr_flag) TThread::Lock();

    // Get TCS phase
    try {
      TCSphase = event.GetTT().GetPhaseTCS();
    } catch (std::logic_error& e) {
      std::cerr << "PlanePGEM::EndEvent: TCS phase not available: " << e.what() << std::endl;
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
      int   pixX    = (int) digit.dt[0];
      int   pixY    = (int) digit.dt[1];
      int   chip    = (int) digit.dt[10];

      double a0      = digit.dt[2];
      double a1      = digit.dt[3];
      double a2      = digit.dt[4];
      double a02     = (a2>0.) ? a0/a2 : 0.;
      double a12     = (a2>0.) ? a1/a2 : 0.;
      double cmc2    = digit.dt[9];

      // store data into variables
      fVch->Store(channel);
      fVa0->Store(a0);
      fVa1->Store(a1);
      fVa2->Store(a2);

      // fill histograms
      fHch->Fill(channel);
      fHhitMapPix->Fill(pixX, pixY);
      if (fExpertHistos)
        fHhitMapSumAmp->Fill(pixX, pixY, a2);
      fHa0->Fill(a0);
      fHa1->Fill(a1);
      fHa2->Fill(a2);
      fHchvsa2->Fill(channel, a2);

      // fill CMC per chip, only once per event
      if (channel/128 < chips) {
        if (filledCMC[channel/128]== false) {
          filledCMC[channel/128] = true;
#if USE_DATABASE==1
          if (maxBaselines.count(chip)>0)  {
            fHchipvscmc2->Fill(channel/128, cmc2-(1024-maxBaselines[chip]));
            if (fExpertHistos)
              vHchipcmc2s[channel/128]->Fill(cmc2-(1024-maxBaselines[chip]));
          } else {
#endif
            fHchipvscmc2->Fill(channel/128, cmc2);
            if (fExpertHistos)
              vHchipcmc2s[channel/128]->Fill(cmc2);
#if USE_DATABASE==1
            if (!fPrintedMissCalib) {
              std::cerr << "PlanePGEM::EndEvent: No mean baseline for APV chips of plane " << fName << "..." << std::endl;
              fPrintedMissCalib=true;
            }
          }
#endif
        }
      } else
        std::cerr << "Strange channel number..." << std::endl;

      fHampRatio->Fill(a12,a02);
 
      double time;
      if ( pixelGEMGeom && pixelGEMGeom->CalcTime(digit, time) ) { // is there a valid time, variable still has to be filled
        fHt->Fill(TCSphase - time);
        fVt->Store(TCSphase - time);
      } else
        fVt->Store(101.);

      fNhitsKept++;
    }
    fHhits->Fill(fNhitsKept);

    if (!(fRateCounter%25)) { // update occupancy histogram
      fHoccupPix->Reset("ICE");
      fHoccupPix->Add(fHhitMapPix, 1/double(fRateCounter));

      if (fExpertHistos)
        fHhitMapMeanAmp->Divide(fHhitMapSumAmp, fHhitMapPix);
    }

    // cleanup check of CMC
    delete [] filledCMC;

    if (thr_flag) TThread::UnLock();
  } else { // strip plane of PixelGEM detector
    PlaneGEM::EndEvent(event);
  }
}

void PlanePGEM::Clusterize() {
  if (fPixel) { // pixel plane
    // No geometry
    if (!fGeom) return;

    if (thr_flag) TThread::Lock();

    // Do clustering
    std::set<CCluster1>& clusters = fGeom->Clusterize(fDigits);

    // save number of clusters (before and after time cut)
    unsigned int NrClu(0), NrTClu(0);

    // this is were the real important stuff happens - cluster handling
    for (std::set<CCluster1>::iterator cluit=clusters.begin(); cluit!=clusters.end(); cluit++) {
      // get data from clusters
      double pixX  = cluit->dt[3];
      double posX  = cluit->er[3];
      double pixY  = cluit->dt[4];
      double posY  = cluit->er[4];
      double amp0  = cluit->dt[0];
      double amp1  = cluit->dt[1];
      double amp2  = cluit->dt[2];
      int size     = cluit->size;

      // store in variables
      fVcPosPixX->Store(posX);
      fVcPosPixY->Store(posY);
      fVcA2     ->Store(amp2);
      fVcSize   ->Store(size);

      // fill histograms
      fHcPosPix     ->Fill(pixX, pixY);
      if (fExpertHistos)
        fHcSumAmpPix->Fill(pixX, pixY, amp2);
      fHcA0         ->Fill(amp0);
      fHcA1         ->Fill(amp1);
      fHcA2         ->Fill(amp2);
      fHcSize       ->Fill(size);
      fHcAmpRatio   ->Fill(amp1/amp2, amp0/amp2);

      // clusters within tight time window (cut on amplitude ratio)
      const float a02min=0.1;
      const float a02max=0.6;
      const float a12min=0.5;
      const float a12max=0.95;

      float a02=amp0/amp2;
      float a12=amp1/amp2;

      if ( a02>=a02min && a02<=a02max && a12>=a12min && a12<=a12max ) {
        fHcPosPixT     ->Fill(pixX, pixY);
        if (fExpertHistos)
          fHcSumAmpPixT->Fill(pixX, pixY, amp2);
        fHcA2T         ->Fill(amp2);
        fHcSizeT       ->Fill(size);
        fHcAmpRatioT   ->Fill(amp1/amp2, amp0/amp2);
 
        // increase number of clusters after time cut
        NrTClu++;
      }

      // time information
      if ( cluit->dt.size() > 5 ) { // there is a valid time
        fHcTime->Fill(TCSphase - cluit->dt[5]);
        if ( a02>=a02min && a02<=a02max && a12>=a12min && a12<=a12max )
          fHcTimeT->Fill(TCSphase - cluit->dt[5]);

        fVcTime->Store(TCSphase - cluit->dt[5]);
      } else // we still have to fill the variable which is saved into the tree, histogram is not filled
        fVcTime->Store(101.);

      // increase number of clusters
      NrClu++;
    }
    fHcHits ->Fill(NrClu);
    fHcHitsT->Fill(NrTClu);
    fNclustKept = NrClu;

    if (!(fRateCounter%25)) { // update mean amplitudes histogram
      if (fExpertHistos) {
        fHcMeanAmpPix ->Divide(fHcSumAmpPix,  fHcPosPix);
        fHcMeanAmpPixT->Divide(fHcSumAmpPixT, fHcPosPixT);
      }
    }

    if (thr_flag) TThread::UnLock();
  } else { // strip plane of PixelGEM detector
    PlaneGEM::Clusterize();
  }
}










