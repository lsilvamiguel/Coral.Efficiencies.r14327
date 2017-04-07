#include "PlanePMumega.h"
#include "GeomPlanePMumega.h"
#include "ChipAPV.h"
#include "ChipAPVRICH.h"
#include "PlaneAPV.h"

#include <TF1.h>
#include <cmath>
#include <stdexcept>
//#include <algo.h>

ClassImp (PlanePMumega);


PlanePMumega::PlanePMumega(const char *detname, int nchan, int center, int width, int pixeltype)
  : PlaneAPV(detname,nchan, center, width, pixeltype==1) {

  fPixelMM=(pixeltype==2);
  fPixelMMsimplified = false;
  if (nchan == 768) fPixelMMsimplified = true;
  INeedGeom();
  fVa0->SetRange(-100, 900);
  fVa1->SetRange(-100, 900);
  fVa2->SetRange(-100, 900);

#if USE_DATA == 1
  calib_arr = 0;
#endif

}

void PlanePMumega::Init(TTree* tree) {
  // override the default settings for hittimes from planeAPV before
  // calling the init method, ad there the histograms are created
  fVt->SetRange(-100., 100.);
  fVt->SetNbins(200);

  PlaneAPV::Init(tree);

    // Geometry and calibrations
  if (fGeom) {
    GeomPlanePMumega* geompl = (GeomPlanePMumega*) fGeom; // GeomPlanePMumega needs coding
    geompl->SetParameters();
#if USE_DATABASE==1
    geompl->SetCalibrations(calib_data, fRunMaps, calib_time);
#else
    geompl->SetCalibrations(fRunMaps);
#endif
  }

  for (register int iapv = 0; iapv < MAX_NBAPV_MumegaAPV; iapv++) {
    fCMapv[iapv][0] = 0.;
    fCMapv[iapv][1] = 0.;
    fCMapv[iapv][2] = 0.;
  }

  for (register int ii = 0; ii < MAX_STAT_MumegaAPV; ii++) {
    fStat[ii].sum = 0.; fStat[ii].sum2 = 0.; fStat[ii].nb =0;
    fStatCM[ii].sum = 0.; fStatCM[ii].sum2 = 0.; fStatCM[ii].nb = 0;
  }

  std::string name;

  // Time from a1/a2 vs tcsphase
  name = fName + "_a1overa2_vs_TCSphase";
  fHtvsTCSph=new TH2F(name.c_str(),name.c_str(), 40, -10, 30,
                      100, 0. ,2.);
  fHtvsTCSph->SetOption("colz");
  AddHistogram(fHtvsTCSph);

  // a1/a2 amplitude ratio
  name = fName + "_a1overa2";
  fHa1a2Ratio = new TH1F_Ref(name.c_str(),name.c_str(),
                             100, 0., 2.,fRateCounter);
  ((TH1F_Ref*)fHa1a2Ratio)->SetReference(fReferenceDirectory);
  AddHistogram(fHa1a2Ratio);
  fHa1a2Ratio->GetXaxis()->SetTitle("a1/a2");

  // Channel amplitude ratio
  name = fName + "_AmpRatio";
  fHampRatio = new TH2F(name.c_str(), name.c_str(), 100, -0.01, 1.99, 100, -0.01, 1.99);
  fHampRatio->SetOption("colz");
  AddHistogram(fHampRatio);
  fHampRatio->GetXaxis()->SetTitle("a1/a2");
  fHampRatio->GetYaxis()->SetTitle("a0/a2");

  // Channel vs strip amplitude
  name = fName + "_ch_vs_a2";
  fHchvsa2 = new TH2F(name.c_str(),name.c_str(), fNchan, -0.5, fNchan-0.5, 500, -100, 900);
  fHchvsa2->SetOption("colz");
  AddHistogram(fHchvsa2);

  // Amplitude 3rd sample
  std::string a2name = fName + fVa2->GetName() + "_vs_spill";
  fHa2s = new TH2F(a2name.c_str(), a2name.c_str(), 200, 0., 200.,
                   fVa2->GetNbins(),
                   fVa2->GetMin(),
                   fVa2->GetMax());
  AddHistogram(fHa2s);

  if (fExpertHistos) {

    // Amplitude 1st sample
    std::string a0name = fName + fVa0->GetName() + "_vs_spill";
    fHa0s = new TH2F(a0name.c_str(), a0name.c_str(), 200, 0., 200.,
                     fVa0->GetNbins(),
                     fVa0->GetMin(),
                     fVa0->GetMax());
    AddHistogram(fHa0s);

    // Amplitude 2nd sample
    std::string a1name = fName + fVa1->GetName() + "_vs_spill";
    fHa1s = new TH2F(a1name.c_str(), a1name.c_str(), 200, 0., 200.,
                     fVa1->GetNbins(),
                     fVa1->GetMin(),
                     fVa1->GetMax());
    AddHistogram(fHa1s);

    // a1/a2 versus spill
    std::string a1a2name = fName + "_a1overa2_vs_spill";
    fHa1a2s = new TH2F(a1a2name.c_str(), a1a2name.c_str(), 200, 0., 200., 100, 0., 2.);
    AddHistogram(fHa1a2s);
    fHa1a2s->GetXaxis()->SetTitle("number of spill");
    fHa1a2s->GetYaxis()->SetTitle("a1/a2");

    // a1/a2 CM versus spill
    std::string a1a2CMname = fName + "_a1overa2CM_vs_spill";
    fHa1a2CMs = new TH2F(a1a2CMname.c_str(), a1a2CMname.c_str(), 200, 0., 200., 100, 0., 2.);
    AddHistogram(fHa1a2CMs);
    fHa1a2s->GetXaxis()->SetTitle("number of spill");
    fHa1a2s->GetYaxis()->SetTitle("a1/a2_CM");

    // a1/a2 CM amplitude ratio
    name = fName + "_a1overa2CM";
    fHa1a2RatioCM = new TH1F_Ref(name.c_str(),name.c_str(),
                               100, 0., 2.,fRateCounter);
    ((TH1F_Ref*)fHa1a2RatioCM)->SetReference(fReferenceDirectory);
    AddHistogram(fHa1a2RatioCM);
    fHa1a2RatioCM->GetXaxis()->SetTitle("a1/a2_CM");

    // Times from channel amplitudes ratio
    name = fName + "_a0/a2_vs_a1/a2_a2_gt_20";
    fHtimeRatio = new TH2F(name.c_str(),name.c_str(),
                           150, 0., 3., 150, 0., 3.);
    AddHistogram(fHtimeRatio);
    fHtimeRatio->GetXaxis()->SetTitle("a1/a2");
    fHtimeRatio->GetYaxis()->SetTitle("a0/a2");
    fHtimeRatio->SetOption("colz");

    // Channel vs strip amplitude
    name = fName + "_ch_vs_a2ped";
    fHchvsa2ped = new TH2F(name.c_str(),name.c_str(),
                           fNchan, -0.5, fNchan-0.5,
                           500, -100, 900);
    AddHistogram(fHchvsa2ped);

    // Channel vs strip amplitude
    name = fName + "_ch_vs_a2CM";
    fHchvsa2CM = new TH2F(name.c_str(),name.c_str(),
                          fNchan, -0.5, fNchan-0.5,
                          400, -400, 400);
    AddHistogram(fHchvsa2CM);

    // a2 CM average profile
    name = fName + "_a2_avg_CM";
    fHavgampCM = new TH1D(name.c_str(),name.c_str(),
                      fVch->GetNbins(),
                      fVch->GetMin(), fVch->GetMax());
    AddHistogram(fHavgampCM);

    // a2 CM sigma profile
    name = fName + "_a2_sigma_CM";
    fHsigampCM = new TH1D(name.c_str(),name.c_str(),
                        fVch->GetNbins(),
                        fVch->GetMin(), fVch->GetMax());
    AddHistogram(fHsigampCM);

  }

  // Single channel occupancies
  name = fName + "_occupancies";
  fHoccup = new TH1F_Ref(name.c_str(), name.c_str(), fVch->GetNbins(), fVch->GetMin(), fVch->GetMax(),fRateCounter,true);
  ((TH1F_Ref*)fHoccup)->SetReference(fReferenceDirectory);
  AddHistogram(fHoccup);
  fHoccup->SetStats(false);

  // a2 average profile
  name = fName + "_a2_avg";
  fHavgamp=new TH1D(name.c_str(),name.c_str(),
                    fVch->GetNbins(),
                    fVch->GetMin(), fVch->GetMax());
  AddHistogram(fHavgamp);

  // a2 sigma profile
  name = fName + "_a2_sigma";
  fHsigamp = new TH1D(name.c_str(),name.c_str(),
                    fVch->GetNbins(),
                    fVch->GetMin(), fVch->GetMax());
  AddHistogram(fHsigamp);

  if (fPixel) { //pixel plane

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
  }

  else if (fPixelMM) {
    // 2D hits counting
    name=fName+"_HitMapMM";
    register int nvert = 40;
    if  (fPixelMMsimplified) nvert = 12;
    fHhitMapPix = new TH2F(name.c_str(), name.c_str(), 128, -25.6, 25.6, nvert, -25.0, 25.0);
    AddHistogram(fHhitMapPix);
    fHhitMapPix->SetStats(false);
    fHhitMapPix->SetOption("COLZ");

    // 2D pad occupancies
    name=fName +"_OccupancyMM";
    fHoccupPix = new TH2F(name.c_str(), name.c_str(), 128, -25.6, 25.6, nvert, -25.0, 25.0);
    AddHistogram(fHoccupPix);
    fHoccupPix->SetStats(false);
    fHoccupPix->SetOption("COLZ");

    // 2D normalized hits counting
    name=fName+"_HitMapMMNorm";
    fHhitMapPixNorm = new TH2F(name.c_str(), name.c_str(), 128, -25.6, 25.6, nvert, -25.0, 25.0);
    AddHistogram(fHhitMapPixNorm);
    fHhitMapPixNorm->SetStats(false);
    fHhitMapPixNorm->SetOption("COLZ");

    // 2D normalized pad occupancies
    name=fName +"_OccupancyMMNorm";
    fHoccupPixNorm = new TH2F(name.c_str(), name.c_str(), 128, -25.6, 25.6, nvert, -25.0, 25.0);
    AddHistogram(fHoccupPixNorm);
    fHoccupPixNorm->SetStats(false);
    fHoccupPixNorm->SetOption("COLZ");

    // hits counting vs pixel number
    name=fName+"_hitVsPix";
    fHhitVsPix = new TH1D(name.c_str(), name.c_str(), 1280, 0., 1280.);
    AddHistogram(fHhitVsPix);

    // pad occupancies vs pixel number
    name=fName +"_occupancyVsPix";
    fHoccupVsPix = new TH1D(name.c_str(), name.c_str(), 1280, 0., 1280);
    AddHistogram(fHoccupVsPix);
  }


  fVcPos=0;
  fVcPosPixX = 0;
  fVcPosPixY = 0;
  // cluster variables
  if (fPixel) { //pixel plane
    fVcPosPixX = AddVariable("Clu_PosPixX", 32, -0.5, 31.5, fNchan*fMAX_MULT);
    fVcPosPixY = AddVariable("Clu_PosPixY", 32, -0.5, 31.5, fNchan*fMAX_MULT);
  }
  else if(fPixelMM)
  {
    fVcPosPixX = AddVariable("Clu_PosPixX",128, -25.6, 25.6,fNchan*fMAX_MULT);
    fVcPosPixY = AddVariable("Clu_PosPixY", 40, -25.0, 25.0,fNchan*fMAX_MULT);
  }
  else { // strip plane of PixelMicromegas detector
    if ( fName.substr(0, 2)=="GP" ) // PixelGEM has two hemispheres, they are combined for clusters
      fVcPos  = AddVariable("Clu_Profile", fNchan/2, -0.5, (fNchan/2)-0.5, fNchan*fMAX_MULT);
    else
      fVcPos  = AddVariable("Clu_Profile",   fNchan, -0.5,     fNchan-0.5, fNchan*fMAX_MULT);
  }
  fVcA2 = AddVariable("Clu_Amp2", 1000, 0., 4000., fNchan*fMAX_MULT);
  fVcSize = AddVariable("Clu_Size", 30, -0.5, 29.5, fNchan*fMAX_MULT);
  fVcTime = AddVariable("Clu_Time", 100, -100., 100., fNchan*fMAX_MULT);

  // cluster histograms
  // Number of clusters
  std::string chitsname = fName + "Clu_Multiplicity";
  fHcHits = new TH1F(chitsname.c_str(), chitsname.c_str(), 70, -0.5, 69.5);
  AddHistogram(fHcHits);

  // Cluster profile
  std::string cchname;
  if (fPixel) { //pixel plane
    cchname = fName + "Clu_PosPix";
    fHcPosPix = new TH2F(cchname.c_str(), cchname.c_str(), fVcPosPixX->GetNbins(),\
                         fVcPosPixX->GetMin(), fVcPosPixX->GetMax(), fVcPosPixY->GetNbins(),\
                         fVcPosPixY->GetMin(), fVcPosPixY->GetMax());
    fHcPosPix->SetStats(false);
    fHcPosPix->SetOption("COLZ");
    AddHistogram(fHcPosPix);
  } else if(fPixelMM){
    cchname = fName + "Clu_PosPixMM";
    fHcPosPix = new TH2F(cchname.c_str(), cchname.c_str(), fVcPosPixX->GetNbins(),\
                         fVcPosPixX->GetMin(), fVcPosPixX->GetMax(), fVcPosPixY->GetNbins(),\
                         fVcPosPixY->GetMin(), fVcPosPixY->GetMax());
    fHcPosPix->SetStats(false);
    fHcPosPix->SetOption("COLZ");
    AddHistogram(fHcPosPix);
  } else { // strip plane of PixelMicromegas detector
    cchname = fName + fVcPos->GetName();
    fHcPos=new TH1F(cchname.c_str(), cchname.c_str(), fVcPos->GetNbins(), fVcPos->GetMin(), fVcPos->GetMax());
    AddHistogram(fHcPos);
  }

  // Cluster amplitude 1
  name = fName + "Clu_Amp0";
  fHcA0 = new TH1F(name.c_str(), name.c_str(), 1000, 0., 4000.);
  AddHistogram(fHcA0);

  // Cluster amplitude 2
  name = fName + "Clu_Amp1";
  fHcA1 = new TH1F(name.c_str(), name.c_str(), 1000, 0., 4000.);
  AddHistogram(fHcA1);

  // Cluster amplitude 3
  std::string ca2name = fName + fVcA2->GetName();
  fHcA2 = new TH1F(ca2name.c_str(), ca2name.c_str(), fVcA2->GetNbins(), fVcA2->GetMin(),\
                   fVcA2->GetMax());
  AddHistogram(fHcA2);

  // Cluster size
  std::string csname = fName + fVcSize->GetName();
  fHcSize = new TH1F(csname.c_str(), csname.c_str(), fVcSize->GetNbins(), fVcSize->GetMin(),\
                     fVcSize->GetMax());
  AddHistogram(fHcSize);

  // Cluster amplitude ratio ("banana")
  name = fName + "Clu_AmpRatio";
  fHcAmpRatio = new TH2F(name.c_str(), name.c_str(), 100, -0.01, 1.99, 100, -0.01, 1.99);
  AddHistogram(fHcAmpRatio);
  fHcAmpRatio->SetStats(false);
  fHcAmpRatio->GetXaxis()->SetTitle("A1/A2");
  fHcAmpRatio->GetYaxis()->SetTitle("A0/A2");

  // Cluster time
  std::string ctname = fName + fVcTime->GetName();
  fHcTime = new TH1F(ctname.c_str(), ctname.c_str(), fVcTime->GetNbins(), fVcTime->GetMin(),\
                     fVcTime->GetMax());
  AddHistogram(fHcTime);
  fHcTime->GetXaxis()->SetTitle("time / ns");

  // histograms after time cut
  // Number of clusters after time cut
  name = fName + "TClu_Multiplicity";
  fHcHitsT = new TH1F(name.c_str(), name.c_str(), 70, -0.5, 69.5);
  AddHistogram(fHcHitsT);

  // Cluster profile after time cut
  if (fPixel) { //pixel plane
    name = fName + "TClu_PosPix";
    fHcPosPixT=new TH2F(name.c_str(), name.c_str(), fVcPosPixX->GetNbins(), fVcPosPixX->GetMin(), fVcPosPixX->GetMax(),	fVcPosPixY->GetNbins(), fVcPosPixY->GetMin(), fVcPosPixY->GetMax());
    fHcPosPixT->SetStats(false);
    fHcPosPixT->SetOption("COLZ");
    AddHistogram(fHcPosPixT);
  } else if(fPixelMM) {
    name = fName + "TClu_PosPixMM";
    fHcPosPixT=new TH2F(name.c_str(), name.c_str(), fVcPosPixX->GetNbins(), fVcPosPixX->GetMin(), fVcPosPixX->GetMax(),	fVcPosPixY->GetNbins(), fVcPosPixY->GetMin(), fVcPosPixY->GetMax());
    fHcPosPixT->SetStats(false);
    fHcPosPixT->SetOption("COLZ");
    AddHistogram(fHcPosPixT);
  } else { // strip plane of PixelMicromegas detector
    name = fName + "T" + fVcPos->GetName();
    fHcPosT=new TH1F(name.c_str(),name.c_str(), fVcPos->GetNbins(), fVcPos->GetMin(),\
                     fVcPos->GetMax());
    AddHistogram(fHcPosT);
  }

  // Cluster amplitude 3 after time cut
  name = fName + "T" + fVcA2->GetName();
  fHcA2T=new TH1F(name.c_str(), name.c_str(), fVcA2->GetNbins(), fVcA2->GetMin(),\
                  fVcA2->GetMax());
  AddHistogram(fHcA2T);

  // Cluster size after time cut
  name = fName + "TClu_Size";
  fHcSizeT = new TH1F(name.c_str(), name.c_str(), fVcSize->GetNbins(),\
                      fVcSize->GetMin(), fVcSize->GetMax());
  AddHistogram(fHcSizeT);

  // Cluster amplitude ratio ("banana") after time cut
  name = fName + "TClu_AmpRatio";
  fHcAmpRatioT = new TH2F(name.c_str(), name.c_str(), 100, -0.05, 1.95, 100, -0.05, 1.95);
  AddHistogram(fHcAmpRatioT);
  fHcAmpRatioT->SetStats(false);
  fHcAmpRatioT->GetXaxis()->SetTitle("A1/A2");
  fHcAmpRatioT->GetYaxis()->SetTitle("A0/A2");

  // Cluster time after time cut
  name = fName + "T" + fVcTime->GetName();
  fHcTimeT = new TH1F(name.c_str(), name.c_str(), fVcTime->GetNbins(),\
                      fVcTime->GetMin(), fVcTime->GetMax());
  AddHistogram(fHcTimeT);
  fHcTimeT->GetXaxis()->SetTitle("time / ns");

  // create tree leaves
  if (tree) {
    std::string chitsleavlist = chitsname + "/I";
    std::string cpxleavlist   = cchname + "X[" + chitsname +"]/F";
    std::string cpyleavlist   = cchname + "Y[" + chitsname +"]/F";
    std::string cchleavlist   = cchname + "[" + chitsname +"]/F";
    std::string ca2leavlist   = ca2name +  "[" + chitsname +"]/F";
    std::string csleavlist    = csname  +  "[" + chitsname +"]/F";
    std::string ctleavlist    = ctname  +  "[" + chitsname +"]/F";

    tree->Branch(chitsname.c_str(), &fNclustKept,          chitsleavlist.c_str(), 32000);
    tree->Branch(ca2name.c_str()  ,  fVcA2->GetValues(),   ca2leavlist.c_str(),   32000);
    if (fPixel) { //pixel plane
      tree->Branch((cchname+"X").c_str(), fVcPosPixX->GetValues(), cpxleavlist.c_str(), 32000);
      tree->Branch((cchname+"Y").c_str(), fVcPosPixY->GetValues(), cpyleavlist.c_str(), 32000);
    } else if(fPixelMM) {
      tree->Branch((cchname+"XMM").c_str(), fVcPosPixX->GetValues(), cpxleavlist.c_str(), 32000);
      tree->Branch((cchname+"YMM").c_str(), fVcPosPixY->GetValues(), cpyleavlist.c_str(), 32000);
    } else { // strip plane of PixelMicromegas detector
      tree->Branch(cchname.c_str()  ,  fVcPos->GetValues(),  cchleavlist.c_str(),   32000);
    }
    tree->Branch(csname.c_str()   ,  fVcSize->GetValues(), csleavlist.c_str(),    32000);
    tree->Branch(ctname.c_str()   ,  fVcTime->GetValues(), ctleavlist.c_str(),    32000);
  }
}

void PlanePMumega::Reset() {
  PlaneAPV::Reset();

  if (fGeom) {
    GeomPlanePMumega *PMumegaGeom = dynamic_cast<GeomPlanePMumega*>(fGeom);
    if (PMumegaGeom) PMumegaGeom->ResetPlane();
  }

  fNclustKept = 0;
  lDigits.clear();
  fDigits.clear();
}


void PlanePMumega::ResetHistograms() {
  PlaneAPV::ResetHistograms();

  for (register int ii = 0; ii < MAX_STAT_MumegaAPV; ii++) {
    fStat[ii].sum = 0.; fStat[ii].sum2 = 0.; fStat[ii].nb = 0;
    fStatCM[ii].sum = 0.; fStatCM[ii].sum2 = 0.; fStatCM[ii].nb = 0;
  }
}



void PlanePMumega::StoreDigit(float x, float y, int channel, int amp0, int amp1, int amp2, int cmc0, int cmc1, int
cmc2, int chip, int pix, int connb) {
  if (fNhits==0) fDigits.clear();
  float a12ratio = (amp2==0) ? 0 : float(amp1)/float(amp2);
  float a02ratio = (amp2==0) ? 0 : float(amp0)/float(amp2);
  std::vector<double> data;
  data.push_back(x);
  data.push_back(y);
  data.push_back(amp0);
  data.push_back(amp1);
  data.push_back(amp2);
  data.push_back(a12ratio);
  data.push_back(a02ratio);
  data.push_back(cmc0);
  data.push_back(cmc1);
  data.push_back(cmc2);
  data.push_back(chip);
  data.push_back(pix);
  data.push_back(connb);
  //data.push_back(a02ratio);
  int chanchip = channel + connb*MAX_APV_NBCH;
  fDigits.insert(make_pair(chanchip,CDigit1(chanchip,data)));
  fNhits++;

// if (fPixelMMsimplified) cerr<<"chanchip "<<chanchip<<" channel "<<channel<<" connb "<<connb<<" x "<<x<<" y "<<y<<" amp2 "<<amp2<<" pix "<<pix<<" chip "<<chip<<std::endl;
}



void PlanePMumega::StoreDigit(CS::Chip::Digit* digit) {

  if (fPixelMM) {
   std:: vector<float> data=digit->GetNtupleData();
   if(data.size()>=8) {
    this->StoreDigit( (float)data[0], (float)data[1], (int)data[2],
		      (int)data[5],   (int)data[6],   (int)data[7],
		      (int)data[12],  (int)data[13],  (int)data[14],
		      (int) data[3],  (int)data[15],  (int)data[16]);
   }

  } else PlaneAPV::StoreDigit(digit);

//   lDigits.push_back(digit);
}



void PlanePMumega::EndEvent(const CS::DaqEvent &event) {
  bool sparsefg = true, ampTotCleaned = false;

  // need geometry object for calculation of hits and cuts and ...
  GeomPlanePMumega *pixelMumegaGeom;
  if (fGeom) pixelMumegaGeom = dynamic_cast<GeomPlanePMumega*>(fGeom);
//   GeomPlaneMumega *MumegaGeom(0);
//   if (fGeom) MumegaGeom = dynamic_cast<GeomPlaneMumega*>(fGeom);

  if (thr_flag) TThread::Lock();

  register bool rateflag = false;
  if ((fRateCounter&0x3f) == 0) {
    rateflag = true;
  }

  // Get TCS phase
  try {
    TCSphase = event.GetTT().GetPhaseTCS() - 40.; // minus 40. to stay in agreement
                                                  // with previous definition of TCS
                                                  // phase in COOOL
  } catch (std::logic_error& e) {
    std::cerr << "PlanePMumega::EndEvent: TCS phase not available: " << e.what() << std::endl;
    TCSphase = 0;
  }

//     if ( ! sparsefg) {
//       memset(&ampTot, 0, sizeof(ampTot));
//     }

  // save number of hits
  fNhitsKept = 0;


  // loop over digits
  for (std::map<int,CDigit1>::const_iterator id = fDigits.begin(); id!=fDigits .end(); id++) {
    const CDigit1& digit = id->second;

    int channel = (int )digit.ch;
    float pixX = (float) digit.dt[0];
    float pixY = (float) digit.dt[1];
    int   hem     = (int) digit.dt[9];
    int pix = -1;

    double a0, a1, a2;
    if (fPixel) { // pixel plane
      a0 = digit.dt[2];
      a1 = digit.dt[3];
      a2 = digit.dt[4];
    }
    else if (fPixelMM)
    {
        a0 = digit.dt[2];
        a1 = digit.dt[3];
        a2 = digit.dt[4];
        pix =(int) digit.dt[11];

    }
    else { // strip plane of PixelMicromegas detector
      a0      = digit.dt[0];
      a1      = digit.dt[1];
      a2      = digit.dt[2];
    }
    double a02 = (a2>0.) ? a0/a2 : 0.;
    double a12 = (a2>0.) ? a1/a2 : 0.;
    if ((!fPixel) && (!fPixelMM) && hem) // for GM hem is always 0, for GP either +/-1, so this is only done for GP
      channel += (1-hem)*(fNchan/4);

// if (fPixelMMsimplified) cerr<<"channel "<<channel<<" pixX "<<pixX<<" pixY "<<pixY<<" hem "<<hem<<" pix "<<pix<<" conn "<<digit.dt[12]<<" a2 "<<a2<<endl;

    // store data into variables
    fVch->Store(channel);
    fVa0->Store(a0);
    fVa1->Store(a1);
    fVa2->Store(a2);

    // applys cuts on variables
    if (! ( fVch->Test(channel) && fVa0->Test(a0) && fVa1->Test(a1) && fVa2->Test(a2) ) ) continue;

    // fill histograms
    fHch->Fill(channel);
    if (fPixel) { // pixel plane
      fHhitMapPix->Fill(pixX, pixY);
    } else if (fPixelMM) {
      if (fPixelMMsimplified) {
        if(pixX>-12.8 && pixX <12.8 && pixY >-12.5 &&pixY <12.5) {
	  fHhitMapPix->Fill(pixX, pixY);
	  fHhitMapPixNorm->Fill(pixX, pixY);
        } else {
	  fHhitMapPix->Fill(pixX, pixY);
	  fHhitMapPix->Fill(pixX, pixY+4.16666667);
	  fHhitMapPix->Fill(pixX, pixY-4.16666667);
	  fHhitMapPixNorm->Fill(pixX, pixY, 1/3.);
	  fHhitMapPixNorm->Fill(pixX, pixY+4.16666667, 1/3.);
	  fHhitMapPixNorm->Fill(pixX, pixY-4.16666667, 1/3.);
        }
      } else {
        if(pixX>-12.8 && pixX <12.8 && pixY >-12.5 &&pixY <12.5) {
	  fHhitMapPix->Fill(pixX, pixY+0.625);
	  fHhitMapPix->Fill(pixX, pixY-0.625);
	  fHhitMapPixNorm->Fill(pixX, pixY+0.625);
	  fHhitMapPixNorm->Fill(pixX, pixY-0.625);
        } else  {
	  fHhitMapPix->Fill(pixX, pixY);
	  fHhitMapPix->Fill(pixX, pixY+1.25);
	  fHhitMapPix->Fill(pixX, pixY+2.5);
	  fHhitMapPix->Fill(pixX, pixY-1.25);
	  fHhitMapPix->Fill(pixX, pixY-2.5);
	  fHhitMapPixNorm->Fill(pixX, pixY, 2/5.);
	  fHhitMapPixNorm->Fill(pixX, pixY+1.25, 2/5.);
	  fHhitMapPixNorm->Fill(pixX, pixY+2.5, 2/5.);
	  fHhitMapPixNorm->Fill(pixX, pixY-1.25, 2/5.);
	  fHhitMapPixNorm->Fill(pixX, pixY-2.5, 2/5.);
        }
      }

      fHhitVsPix->Fill(pix);

//     std::cout<<"channel="<<channel<<" pix="<<pix<<" pixX ="<<pixX<<" pixY ="<<pixY<<std::endl;

    }
//     else { // strip plane of PixelMicromegas detector
//     }
    fHa0->Fill(a0);
    fHa1->Fill(a1);
    fHa2->Fill(a2);

    fHchvsa2->Fill(channel, a2);
    fHa1a2Ratio->Fill(a12);
    fHtvsTCSph->Fill(TCSphase, a12);
    fHampRatio->Fill(a12, a02);
    if (fExpertHistos) {
      fHa0s->Fill(event.GetBurstNumber(), a0);
      fHa1s->Fill(event.GetBurstNumber(), a1);
      fHa2s->Fill(event.GetBurstNumber(), a2);
      fHa1a2s->Fill(event.GetBurstNumber(), a12);
      if (a2>20) fHtimeRatio->Fill(float(a12), float(a02));
    }
    if (!(fRateCounter%100)) { // update occcupancy histogram
      fHoccup->Reset("ICE");
      fHoccup->Add(fHch, 1/double(fRateCounter));
      if (fPixel) { // pixel plane
        fHoccupPix->Reset("ICE");
        fHoccupPix->Add(fHhitMapPix, 1/double(fRateCounter));
      } else if (fPixelMM) {
        fHoccupPix->Reset("ICE");
        fHoccupPix->Add(fHhitMapPix, 1/double(fRateCounter));
        fHoccupPixNorm->Reset("ICE");
        fHoccupPixNorm->Add(fHhitMapPixNorm, 1/double(fRateCounter));
        fHoccupVsPix->Reset("ICE");
        fHoccupVsPix->Add(fHhitVsPix, 1/double(fRateCounter));
      }
    }


    fStat[channel].sum += a2;
    fStat[channel].sum2 += a2*a2;
    fStat[channel].nb++;
    if (fPixel) { // pixel plane
      fStat[channel].xpx = short(pixX);
      fStat[channel].ypx = short(pixY);
    }


    //      double time;
    //     if ( pixelMumegaGeom && pixelMumegaGeom->CalcTime(digit, time) ) { // is there a valid time, variable still has tio be filled
    //        fHt->Fill(TCSphase - time);
    //        fVt->Store(TCSphase - time);
    //      } else
    //    if ( MumegaGeom && MumegaGeom->CalcTime(digit, time) ) { // is there a valid time, variable has still to be filled
    //      fHt->Fill(TCSphase - time);
    //      fVt->Store(TCSphase - time);
    //    } else
   // fVt->Store(101.);
    fVt->Store(a12); // temporary solution before to take TCSPhase into account
    fHt->Fill(a12);

    fNhitsKept++;



    if ( 1) {
      if (!ampTotCleaned) {
        memset(&ampTot, 0, sizeof(ampTot));
        ampTotCleaned = true;
      }
      register int apvnb = (channel & ~0x7f) >> 7;
      register int chnb = channel & 0x7f;
      ampTot[apvnb].ampSort[0].amp[chnb] = a0;
      ampTot[apvnb].ampSort[1].amp[chnb] = a1;
      ampTot[apvnb].ampSort[2].amp[chnb] = a2;
      ampTot[apvnb].ampFlag[chnb] = 1;
    }

  }


  if (rateflag) {
    for (register int ii=0; ii<MAX_STAT_MumegaAPV; ii++) {
      register int nb = fStat[ii].nb;
      if (nb == 0) continue;
      register double avg = fStat[ii].sum / nb;
      register double sigma = sqrt(fStat[ii].sum2/nb - fStat[ii].sum*fStat[ii].sum/nb/nb);
      fHavgamp->SetBinContent(fHavgamp->GetXaxis()->FindBin(ii),avg);
      fHsigamp->SetBinContent(fHsigamp->GetXaxis()->FindBin(ii),sigma);
    }
  }


  if (fExpertHistos) {
    CalcCommonModeCorrection(MAX_AMP_DEVIATION);
    // loop over digits, only used for expert histos
    for (register int iapv = 0; iapv < fNchan/128; iapv++) {
      // Loop over channels
      for ( register int ich = 0; ich < MAX_APV_NBCH; ich++ ) {

        register int ch = ich + iapv*MAX_APV_NBCH;

        if ( 1) {
          register double ped= 0.;

          register double a2CM = -(ampTot[iapv].ampSort[2].amp[ich] - fCMapv[iapv][2] - ped);
          register double a1CM = -(ampTot[iapv].ampSort[1].amp[ich] - fCMapv[iapv][1] - ped);
  //         if (fExpertHistos) {
            fHchvsa2ped->Fill(ch, ampTot[iapv].ampSort[2].amp[ich]-ped);
            fHchvsa2CM->Fill(ch, a2CM);
  //         }

          fStatCM[ch].sum += a2CM;
          fStatCM[ch].sum2 += a2CM*a2CM;
          fStatCM[ch].nb++;
        }

      }

    }
  }




  if (fExpertHistos && rateflag) {
    for (register int ii=0; ii<MAX_STAT_MumegaAPV; ii++) {
      register int nb = fStatCM[ii].nb;
      if (nb == 0) continue;
      register double avgCM = fStatCM[ii].sum / nb;
      register double sigmaCM = sqrt(fStatCM[ii].sum2/nb - fStatCM[ii].sum*fStatCM[ii].sum/nb/nb);
      fHavgampCM->SetBinContent(fHavgampCM->GetXaxis()->FindBin(ii),avgCM);
      fHsigampCM->SetBinContent(fHsigampCM->GetXaxis()->FindBin(ii),sigmaCM);
    }
  }



  typedef std::list<CS::Chip::Digit*>::iterator lDIter;
  //fprintf(stderr, "fDigits size %d lDigits size %d \n", fDigits.size(), lDigits.size());

  // only used for expert histos
  if (fExpertHistos) for (lDIter ii = lDigits.begin(); ii != lDigits.end(); ii++) {

    CS::ChipAPV::Digit* idgt = dynamic_cast<CS::ChipAPV::Digit*> (*ii);

    if (!idgt) {
      std::cerr<<"PlanePMumega::EndEvent: a digit is not a APVPMumega one, strange...\n";
      (*ii)->Print();
      continue;
    }


    int channel = idgt->GetChannel();
    float a0 = float (idgt->GetAmplitude()[0]);
    float a1 = float(idgt->GetAmplitude()[1]);
    float a2 = float(idgt->GetAmplitude()[2]);
    //      sparsefg = idgt->IsSparsifed();
    sparsefg = false;

    if ( ! sparsefg) {
      if (!ampTotCleaned) {
        memset(&ampTot, 0, sizeof(ampTot));
        ampTotCleaned = true;
      }
      register int apvnb = (channel & ~0x7f) >> 7;
      register int chnb = channel & 0x7f;
      ampTot[apvnb].ampSort[0].amp[chnb] = a0;
      ampTot[apvnb].ampSort[1].amp[chnb] = a1;
      ampTot[apvnb].ampSort[2].amp[chnb] = a2;
      ampTot[apvnb].ampFlag[chnb] = 1;
    }

    if ( ! sparsefg) {
      CalcCommonModeCorrection(MAX_AMP_DEVIATION);

      for (register int ch = 0; ch < fNchan; ch++) {
        register int apvnb = (ch & ~0x7f) >> 7;
        register int chnb = ch & 0x7f;
        register double ped;

        if (ampTot[apvnb].ampFlag[chnb] != 1 ) continue;

	if (fUseCalib) {
#if USE_DATABASE == 1
	  ped = calib_arr[ch].pedestal;
#else
	  ped = 0;
#endif
	} else {
          ped = 0;
	}

	register double a2CM = -(ampTot[apvnb].ampSort[2].amp[chnb] - fCMapv[apvnb][2] - ped);
        register double a1CM = -(ampTot[apvnb].ampSort[1].amp[chnb] - fCMapv[apvnb][1] - ped);

        double a12CM = (a2CM>0.) ? a1CM/a2CM : 0.;
        if (fHchvsa2ped) fHchvsa2ped->Fill(ch, ampTot[apvnb].ampSort[2].amp[chnb]-ped);
        if (fHchvsa2CM) fHchvsa2CM->Fill(ch, a2CM);
        if (fHa1a2CMs) fHa1a2CMs->Fill(event.GetBurstNumber(), a12CM);
        if (fHa1a2RatioCM) fHa1a2RatioCM->Fill(a12CM);


        fStatCM[ch].sum += a2CM;
        fStatCM[ch].sum2 += a2CM*a2CM;
        fStatCM[ch].nb++;

        if (rateflag) {
	  for (register int ii=0; ii<MAX_STAT_MumegaAPV; ii++) {
	    register int nb = fStatCM[ii].nb;
	    if (nb == 0) continue;
	    register double avgCM = fStatCM[ii].sum / nb;
	    register double sigmaCM = sqrt(fStatCM[ii].sum2/nb - fStatCM[ii].sum*fStatCM[ii].sum/nb/nb);
	    if (fHavgampCM) fHavgampCM->SetBinContent(fHavgampCM->GetXaxis()->FindBin(ii),avgCM);
	    if (fHsigampCM) fHsigampCM->SetBinContent(fHsigampCM->GetXaxis()->FindBin(ii),sigmaCM);
	  }
        }
      }
    }

  }

  fHhits->Fill(fNhitsKept);

  if (thr_flag) TThread::UnLock();
}


void PlanePMumega::Clusterize() {

  // No geometry
  if (!fGeom) return;

  if (thr_flag) TThread::Lock();

  // Do clustering
  std::set<CCluster1>& clusters = fGeom->Clusterize(fDigits);

  // save number of clusters (before and after time cut)
  unsigned int NrClu(0), NrTClu(0);

  // this is where the real important stuff happens - cluster handling
  for (std::set<CCluster1>::iterator cluit = clusters.begin(); cluit != clusters.end(); cluit++) {
    assert( cluit->size ); // otherwise we have forgotten to remove something above

    // get data from clusters
    double pixX = cluit->dt[3];
    double posX = cluit->er[3];
    double pixY = cluit->dt[4];
    double posY = cluit->er[4];
    double amp0 = cluit->dt[0];
    double amp1 = cluit->dt[1];
    double amp2 = cluit->dt[2];
    double strip = cluit->dt[3];
    int size = cluit->size;
    double pos = cluit->pos;

    // store in variables
    if (fVcPosPixX) fVcPosPixX->Store(posX);
    if (fVcPosPixY) fVcPosPixY->Store(posY);
    if (fVcPos) fVcPos->Store(pos);
    fVcA2->Store(amp2);
    fVcSize->Store(size);

      // fill histograms
    if (fPixel) { // pixel plane
      fHcPosPix->Fill(pixX, pixY);
    } else if(fPixelMM) {
      fHcPosPix->Fill(pixX, pixY);
    } else  { // strip plane of PixelMicromegas detector
      fHcPos->Fill(strip);
    }
    fHcA0->Fill(amp0);
    fHcA1->Fill(amp1);
    fHcA2->Fill(amp2);
    fHcSize->Fill(size);
    fHcAmpRatio->Fill(amp1/amp2, amp0/amp2);

    // clusters within tight time window (cut on amplitude ratio)
    const float a02min = 0.0;
    const float a02max = 0.5;
    const float a12min = 0.2;
    const float a12max = 0.95;

    float a02 = amp0/amp2;
    float a12 = amp1/amp2;

    if ( a02>=a02min && a02<=a02max && a12 >= a12min && a12<=a12max ) {
      if (fPixel) { // pixel plane
        fHcPosPixT->Fill(pixX, pixY);
      } else if(fPixelMM) {
        fHcPosPixT->Fill(pixX, pixY);
      } else  { // strip plane of PixelMicromegas detector
        fHcPosT->Fill(strip);
      }
      fHcA2T->Fill(amp2);
      fHcSizeT->Fill(size);
      fHcAmpRatioT->Fill(amp1/amp2, amp0/amp2);

      // increase number of clusters after time cut
      NrTClu++;
    }

    // time information
    if (fPixel) { // pixel plane
      if (cluit->dt.size() > 5) { // there is a valid time
	fHcTime->Fill(TCSphase - cluit->dt[5]);
	if ( a02 >= a02min && a02 <= a02max &&a12 >= a12min && a12 <= a12max )
	  fHcTimeT->Fill(TCSphase - cluit->dt[5]);

	fVcTime->Store(TCSphase - cluit->dt[5]);
      } else // we still have to fill the variable which is saved into the tree, histogram is not filled
	fVcTime->Store(101.);
    } else if(fPixelMM) {
      if (cluit->dt.size() > 5) { // there is a valid time
// 	fHcTime->Fill(TCSphase - cluit->dt[5]);
	fHcTime->Fill(cluit->dt[5]);
	if ( a02 >= a02min && a02 <= a02max &&a12 >= a12min && a12 <= a12max )
// 	  fHcTimeT->Fill(TCSphase - cluit->dt[5]);
	  fHcTimeT->Fill(cluit->dt[5]);

	fVcTime->Store(TCSphase - cluit->dt[5]);
      } else // we still have to fill the variable which is saved into the tree, histogram is not filled
	fVcTime->Store(101.);
    } else { // strip plane of PixelMicromegas detector
      if ( cluit->dt.size() > 4 ) { // there is a valid time
// 	fHcTime->Fill(TCSphase - cluit->dt[4]);
	fHcTime->Fill(cluit->dt[4]);
	if ( a02 >= a02min && a02 <= a02max && a12 >= a12min && a12 <= a12max )
// 	  fHcTimeT->Fill(TCSphase - cluit->dt[4]);
	  fHcTimeT->Fill(cluit->dt[4]);

	fVcTime->Store(TCSphase - cluit->dt[4]);
      } else // we still have to fill the variable which is saved into the tree, histogram is not filled
	fVcTime->Store(101.);
    }

    // increase number of clusters
    NrClu++;
  }

  fHcHits->Fill(NrClu);
  fHcHitsT->Fill(NrTClu);
  fNclustKept=NrClu;

  if (thr_flag) TThread::UnLock();
}


#if USE_DATABASE == 1
void PlanePMumega::ReadCalib(const tm &t)
{
  //  Read calibration data for each strip
  try{
    ReadFromDataBase(calib_data, t);
    // std::cout<<"PlaneGEM::ReadCalib() ==> "
    // <<this->GetName()<<" calibrations are found !"<<std::endl;

    unsigned int calib_size=0;

    std::map<unsigned int, PlanePMumega::APVcalib>::iterator it;
    for (it=calib_data.begin(); it!=calib_data.end(); it++)
      calib_size += it->second.channel.size();

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
    std::cout<<"PlanePMumega::ReadCalib() ==> "<<GetName()
	     <<" calibrations, valid for ";
    std::cout<<t.tm_mday<<"."<<t.tm_mon+1<<"."<<t.tm_year+1900<<" "
	     <<t.tm_hour<<":"<<t.tm_min<<":"<<t.tm_sec
	     <<", not found in DB"<<std::endl;
  }

  //clear current time calibration
  calib_time.clear();

  // Read calibration data for conversion of amplitude->time for each plane
  if (fDataBase==NULL)
    throw CS::Exception("PlanePMumega::ReadCalib():  data base is not opened.");
  try{
    //    std::cout << "Reading timing calibrations for " << GetName() << std::endl;
    struct tm tt(t);
    CDB::Time tp(mktime(&tt),0);
    std::string strdata("");
    fDataBase->read(fName,strdata,tp,"timing");
    if (strdata == "") {
      std::cerr << "PlanePMumega::ReadCalib() "<<GetName()<<" : no timing calibration file found"<<std::endl;
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
    if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_MP")!=0) {
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
    if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_MP")!=0) {
      std::cout<<std::endl<<"Timing calibrations read in for "<<GetName()<<std::endl;
      // omitting the output of the covariance matrix which is used for error calculation
      std::cout<<"Ratio 02: a0="<<calib_time[ 0]<<", t0="<<calib_time[ 1]<<", r0="<<calib_time[ 2]<<std::endl;
      std::cout<<"Ratio 12: a0="<<calib_time[12]<<", t0="<<calib_time[13]<<", r0="<<calib_time[14]<<std::endl;
      std::cout<<std::endl;
    }
  } else {
    calib_time.clear(); // remove whatever has been loaded

    std::cout<<"PlanePMumega::ReadCalib() ==> "<<GetName()<<" timing calibrations, valid for ";
    std::cout<<t.tm_mday<<"."<<t.tm_mon+1<<"."<<t.tm_year+1900<<" "
	     <<t.tm_hour<<":"<<t.tm_min<<":"<<t.tm_sec
	     <<", not found in DB"<<std::endl;
  }
}


istream& operator>>(istream& in, std::map<unsigned int, PlanePMumega::APVcalib> &c) {
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
      if ( c.count(chipId+16*adcId)==0 ) { // first entry found for this chip ID
	c.insert( std::pair<unsigned int, PlanePMumega::APVcalib>(chipId+16*adcId, PlanePMumega::APVcalib()) );
	ignore=false;
      } else { // this seems to be a broken pedestal file, two times the same APV ID... Ignore input till next chip
	std::cerr << "Broken pedestal file found for pixelMM. The same APV ID "<<chipId<<" ADC "<<adcId<<", sID "<<srcId<<" was found several times." << std::endl;
	ignore=true;
      }
    } else {
      int f = atoi(&firstChar);
      float p, s, cp, cs;
      in >> p >> s >> cp >> cs;
      if (!ignore)
	c[chipId+16*adcId].channel.push_back(PlanePMumega::APVchannel(f, p, s, cp, cs, ldcId, srcId, adcId, chipId));
    }
  }

  return in;
}

#endif //USE_DATABASE

void PlanePMumega::TextOutput(ostream& out) {

  // Write maximum of Landau of cluster amplitude 3 to file
  if(float gain = Gain())
    out<<"\tgain "<<gain<<std::endl;
  else
    out<<"\tgain ***"<<std::endl;
}


float PlanePMumega::Gain() {
  static const int statcut = 1000;
  if(fHcA2T->GetEntries() < statcut) return 0;

  TF1 *lan = new TF1("lan","landau",40.,4000.);
  fHcA2T->Fit(lan,"RQ0");
  if (lan->GetChisquare()==0.) // fit gone bad
    return 0;

  return lan->GetParameter(1);
}


void PlanePMumega::CalcCommonModeCorrection(double cut) {

  // Local array for sorted amplitudes of active channels
  ampS3 ampAPV;

  for (register int iapv = 0; iapv < fNchan/128; iapv++) {
    // Loop over channels
    int nb_activech = 0;
    register int iactch = 0;
    memset(&ampAPV, 0, sizeof(ampAPV));
    for ( register int ich = 0; ich < MAX_APV_NBCH; ich++ ) {

      register int ch = ich + iapv*MAX_APV_NBCH;
      // Consider only active channels
#if USE_DATABASE == 1
      if ( fUseCalib ) {
	if ( calib_arr[ch].flag == 1 ) {

          // Copy active channels into new array for sorting
//          ampSort[0].amp[iactch] = fVa0->GetValues()[ch] - calib_arr[ch].pedestal;
//          ampSort[1].amp[iactch] = fVa1->GetValues()[ch] - calib_arr[ch].pedestal;
//          ampSort[2].amp[iactch] = fVa2->GetValues()[ch] - calib_arr[ch].pedestal;
          ampAPV.ampSort[0].amp[iactch] = ampTot[iapv].ampSort[0].amp[ich] - calib_arr[ch].pedestal;
          ampAPV.ampSort[1].amp[iactch] = ampTot[iapv].ampSort[1].amp[ich] - calib_arr[ch].pedestal;
          ampAPV.ampSort[2].amp[iactch] = ampTot[iapv].ampSort[2].amp[ich] - calib_arr[ch].pedestal;
          iactch++;
	}
      } else {

          // Copy active channels into new array for sorting
          ampAPV.ampSort[0].amp[iactch] = ampTot[iapv].ampSort[0].amp[ich];
          ampAPV.ampSort[1].amp[iactch] = ampTot[iapv].ampSort[1].amp[ich];
          ampAPV.ampSort[2].amp[iactch] = ampTot[iapv].ampSort[2].amp[ich];
          iactch++;
      }
#endif
    }
    nb_activech = iactch;

    for (register int ismp=0; ismp<3; ismp++) {

      // Sort in ascending order
      //      double *ampEnd = ampAPV.ampSort[ismp].amp + nb_activech;
      //      sort(ampAPV.ampSort[ismp].amp, ampEnd);
      //	copy(ampSort, ampEnd, ostream_iterator<double>(cout, " "));

      // Calculate median
      register double median;
      register int ctr;
      if ( nb_activech%2 == 1 ) {
        ctr = (nb_activech-1)/2;
        median = ampAPV.ampSort[ismp].amp[ctr];
      }
      else {
        ctr = nb_activech/2;
        median = (ampAPV.ampSort[ismp].amp[ctr-1]+ampAPV.ampSort[ismp].amp[ctr])/2.;
      }

      // Loop over active channels
      register double ampMean = 0;
      register int n = 0;
      for ( register int i=0; i<nb_activech; i++ ) {

        // Sum up channels with amplitudes around median
        if ( ((median-ampAPV.ampSort[ismp].amp[i]) < cut) &&
	     ((ampAPV.ampSort[ismp].amp[i]-ampAPV.ampSort[ismp].amp[ctr+i-n]) < cut) ) {
          ampMean += ampAPV.ampSort[ismp].amp[i];
          n++;
        }
      }

      // Calculate mean of channels without signal for CM correction
      fCMapv[iapv][ismp] = ampMean / n;
    }
  }

}
