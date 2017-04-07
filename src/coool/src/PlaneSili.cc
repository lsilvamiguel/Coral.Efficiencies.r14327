#include <cmath>
#include "PlaneSili.h"
#include <algorithm>
#include <stdexcept>
#include "GeomPlaneSili.h"
#include "ChipAPV.h"

#define SILAPV_A2_CUT 10

ClassImp(PlaneSili);

bool PlaneSili::fXYinit = false;
PlaneSili::xy_s PlaneSili::fXY[50];

double eSi(double s, double *p) {
  double apc=p[2]+p[3];
  double amc=p[2]-p[3];
  double f=s-p[5]+0.5*p[4]*amc;
  double ac=p[2]*p[3];
  double x=0.5/ac*(apc*f-amc*sqrt(f*f+p[4]*p[4]*ac));
  return x+p[1];
}
double APV_t(double x, double *p) {
  return eSi(-log(-log(x/p[0])), p);
}

void PlaneSili::Init(TTree* tree) {
 
  PlaneAPV::Init(tree);

  fHchVsTis=0;
  fHctime1=fHctime2=0;
  fHa0s=fHa1s=fHa2s=fHa0CMs=fHa1CMs=fHa2CMs=0;
  fHavgamp=fHsigamp=fHavgampCM=fHsigampCM=0;
  fHavgamp2D=fHsigamp2D=fHavgampCM2D=fHsigampCM2D=fHevtampCM2D=fHevtampCMcut2D=0;
  fHampSubRatio=0;
  fHampSubRap=0;
  fHa0CM=fHa1CM=fHa2CM=0;


  for (register int iapv = 0; iapv < MAX_NBAPV; iapv++) {
    cCMC_fCMapv[iapv][0] = 0.;
    cCMC_fCMapv[iapv][1] = 0.;
    cCMC_fCMapv[iapv][2] = 0.;
  }

  for (register int ii =0; ii < MAX_STAT; ii++) {
    fStat[ii].sum = 0.; fStat[ii].sum2 = 0.; fStat[ii].nb = 0;
    fStatCM[ii].sum = 0.; fStatCM[ii].sum2 = 0.; fStatCM[ii].nb = 0;
  }


  // create histograms specific to Silicon detectors here
  
  std::string name5 = fName + "_timing"; 
  fHtiming = new TH2F(name5.c_str(),"timing plot",
		      100, -0.2, 2, 100, -0.2, 2.0);
  AddHistogram(fHtiming);
  fHtiming->GetXaxis()->SetTitle("a1/a2");
  fHtiming->GetYaxis()->SetTitle("a0/a2");
  
  
  std::string name6 = fName + "_Ch_vs_Ampl"; 
  fChAmp =  new TH2F(name6.c_str(),"chan-ampl",
		     fNchan, 0, fNchan,
		     1025, -0.5, 1024.5);
  AddHistogram(fChAmp);
  
  // Channel vs common mode correction
  std::string namecmc = fName + "_ch_vs_cmc2";
  fHchvscmc2 = new TH2F(namecmc.c_str(),namecmc.c_str(),
		    fNchan, -0.5, fNchan-0.5,
		    200, -0.5 , 399.5);
  AddHistogram(fHchvscmc2);

  std::string name7 = fName + "_occupancies"; 
  fHoccup = new TH1F_Ref(name7.c_str(), name7.c_str(), fVch->GetNbins(),
			 fVch->GetMin(), fVch->GetMax(),fRateCounter,true);  
  ((TH1F_Ref*)fHoccup)->SetReference(fReferenceDirectory);
  AddHistogram(fHoccup);
  fHoccup->SetStats(false);
  
  // Cluster histograms
  float posmax = 50;
  
  if(fGeom)
    posmax = fGeom->GetHalfSize();
  
  fVcch = AddVariable("_clusterPosition",300,-posmax,posmax,fNchan*fMAX_MULT);
  fVca2 = AddVariable("_clusterSumAmp", 1000, 0, 1000, fNchan*fMAX_MULT);
  fVcs = AddVariable("_clusterSize",30, -0.5, 29.5, fNchan*fMAX_MULT);  
  
  // clusters ----------------------------------------------------------
  
  //multiplicity within the (time) amplitude cut on the whole detector
  std::string chitsname = fName + "_clusterMultiplicity";
  fHchits=new TH1F(chitsname.c_str(),chitsname.c_str(),70,-0.5,69.5);
  AddHistogram(fHchits);
  std::string chitsleavlist = chitsname + "/I";
  
  //branch : cluster map
  std::string cchname = fName + fVcch->GetName();
  fHcch=new TH1F(cchname.c_str(),cchname.c_str(),
		 fVcch->GetNbins(), fVcch->GetMin(), fVcch->GetMax());
  AddHistogram(fHcch);
  std::string cchleavlist = cchname + "[" + chitsname +"]/F";
  
  
  //branch : cluster amplitude
  std::string ca2name = fName + fVca2->GetName();
  fHca2=new TH1F(ca2name.c_str(),ca2name.c_str(),
		 fVca2->GetNbins(), fVca2->GetMin(), fVca2->GetMax());
  AddHistogram(fHca2);
  
  std::string ca2leavlist = ca2name + "[" + chitsname +"]/F";
  
  std::string catcsname = fName + "Amp vs. TCS";
  fHAmpTCStime = new TH2F(catcsname.c_str(),catcsname.c_str(),
			  fVca2->GetNbins(), fVca2->GetMin(), fVca2->GetMax(),
			  50, -5, 30);
  AddHistogram(fHAmpTCStime);
  
  std::string catcsleavlist = catcsname + "[" + chitsname +"]/F";
  
  //Cluster size
  std::string csname = fName+"_clusterSize"; 
  fHcs=new TH1F(csname.c_str(),csname.c_str(),
		fVcs->GetNbins(), fVcs->GetMin(), fVcs->GetMax());
  AddHistogram(fHcs);
  std::string csleavlist = csname + "[" + chitsname +"]/F";
  
  
  
  std::string name8 = fName + "_clusterTiming"; 
  fHctiming = new TH2F(name8.c_str(), "cluster timing plot",
		       100, -0.2, 2, 100, -0.2, 2.0);
  AddHistogram(fHctiming);
  fHctiming->SetStats(false);
  fHctiming->GetXaxis()->SetTitle("a1/a2");
  fHctiming->GetYaxis()->SetTitle("a0/a2");
  
  
  if (fExpertHistos) {
    std::string name9 = fName + "_clusterT1"; 
    fHctime1 = new TH1F(name9.c_str(), "cluster time 1", 50, -100, 100);
    AddHistogram(fHctime1);
    fHctime1->GetXaxis()->SetTitle("time from a0");


    std::string name10 = fName + "_clusterT2"; 
    fHctime2 = new TH1F(name10.c_str(), "cluster time 2", 50, -100, 100);
    AddHistogram(fHctime2);
    fHctime2->GetXaxis()->SetTitle("time from a1");
  }
  
  std::string name11 = fName + "_clusterTCorrelation"; 
  fHctimeCorr = new TH2F(name11.c_str(), "cluster time Correlation",
                         150, -45, 30, 150, -25, 50);
  AddHistogram(fHctimeCorr);
  fHctimeCorr->SetStats(false);
  fHctimeCorr->GetXaxis()->SetTitle("time from a0");
  fHctimeCorr->GetYaxis()->SetTitle("time from a1");
  
  
  std::string name13 = fName + "_TCluster-TCS"; 
  fHctimeCorrClusterTCS = new TH2F(name13.c_str(), "cluster time tcsphase",
				   40, -5, 30, 160, -40, 40);
  AddHistogram(fHctimeCorrClusterTCS);
  fHctimeCorrClusterTCS->SetStats(false);
  fHctimeCorrClusterTCS->GetXaxis()->SetTitle("from tcsphase");
  fHctimeCorrClusterTCS->GetYaxis()->SetTitle("time from clustering");
  fHctimeCorrClusterTCS->SetOption("col");
  
  std::string name13a = fName + "_ratio0-TCS"; 
  fHcRatio0TCS = new TH2F(name13a.c_str(), "a0/a2 vs. tcsphase",
			  40, -5, 30, 130, -0.1, 2.5);
  AddHistogram(fHcRatio0TCS);
  fHcRatio0TCS->SetStats(false);
  fHcRatio0TCS->GetXaxis()->SetTitle("tcsphase");
  fHcRatio0TCS->GetYaxis()->SetTitle("a0/a2");
  fHcRatio0TCS->SetOption("col");
  
  std::string name13b = fName + "_ratio1-TCS"; 
  fHcRatio1TCS = new TH2F(name13b.c_str(), "a1/a2 vs. tcsphase",
			  40, -5, 30, 130, -0.1, 2.5);
  AddHistogram(fHcRatio1TCS);
  fHcRatio1TCS->SetStats(false);
  fHcRatio1TCS->GetXaxis()->SetTitle("from tcsphase");
  fHcRatio1TCS->GetYaxis()->SetTitle("time from clustering");
  fHcRatio1TCS->SetOption("col");
  
  std::string name13c = fName + "_tcalib-tcorrected"; 
  fHcTiTi = new TH2F(name13c.c_str(), "t - t",
		     70, -60, 80, 70, -60, 80);
  AddHistogram(fHcTiTi);
  fHcTiTi->SetStats(true);
  fHcTiTi->GetXaxis()->SetTitle("sergei");
  fHcTiTi->GetYaxis()->SetTitle("correct");
  
  std::string name12 = fName + "_clusterTime"; 
  fHctime = new TH1F(name12.c_str(), "cluster time",
		     100, -50, 50);
  AddHistogram(fHctime); 
  fHctime->GetXaxis()->SetTitle("time (ns) from clusters");
  
  std::string name14 = fName + "_correctedTime"; 
  fHcRightTime = new TH1F(name14.c_str(), "time corrected",
			  160,-80, 80);
  AddHistogram(fHcRightTime); 
  fHcRightTime->GetXaxis()->SetTitle("corrected time (ns)");
  std::string t_cor_leavlist =  name14 + "[" + chitsname + "]/F"; 
  ft_cor = AddVariable("_correctedTime",200,-100,1000,fNchan*fMAX_MULT); 
  
  std::string name15 = fName + "_calibTime"; 
  fHcCalibTime = new TH1F(name15.c_str(), "time corrected",
			  160,-80, 80);
  AddHistogram(fHcCalibTime); 
  fHcCalibTime->GetXaxis()->SetTitle("calibrated time (ns)");  

  fHchVsTis = new TH2F((fName + "_chVsTis").c_str(), (fName + " channel versus time in spill").c_str(),
                       fVch->GetNbins(), fVch->GetMin(), fVch->GetMax(),
                       100, 0., 15.);
  fHchVsTis->SetOption("COLZ");
  AddHistogram(fHchVsTis);

  if (fExpertHistos) {
      // Amplitude 1st sample
    std::string a0name = fName + fVa0->GetName() + "_vs_spill";
    fHa0s=new TH2F(a0name.c_str(),a0name.c_str(), 200, 0., 200.,
		     fVa0->GetNbins(),
		     fVa0->GetMin(),
		     fVa0->GetMax());
    AddHistogram(fHa0s);

    // Amplitude 2nd sample
    std::string a1name = fName + fVa1->GetName() + "_vs_spill";
    fHa1s=new TH2F(a1name.c_str(),a1name.c_str(), 200, 0., 200.,
		     fVa1->GetNbins(),
		     fVa1->GetMin(),
		     fVa1->GetMax());
    AddHistogram(fHa1s);

    // Amplitude 3rd sample
    std::string a2name = fName + fVa2->GetName() + "_vs_spill";
    fHa2s=new TH2F(a2name.c_str(),a2name.c_str(), 200, 0., 200.,
		     fVa2->GetNbins(),
		     fVa2->GetMin(),
		     fVa2->GetMax());
    AddHistogram(fHa2s);

    // Amplitude 1st sample CM
    std::string name16 = fName + fVa0->GetName() + "CM_vs_spill";
    fHa0CMs=new TH2F(name16.c_str(),name16.c_str(), 200, 0., 200.,
 		     fVa0->GetNbins(),
 		     fVa0->GetMin(),
 		     fVa0->GetMax());
    AddHistogram(fHa0CMs);

    // Amplitude 2nd sample CM
    std::string name17 = fName + fVa1->GetName() + "CM_vs_spill";
    fHa1CMs=new TH2F(name17.c_str(),name17.c_str(), 200, 0., 200.,
 		     fVa1->GetNbins(),
 		     fVa1->GetMin(),
 		     fVa1->GetMax());
    AddHistogram(fHa1CMs);

    // Amplitude 3rd sample CM
    std::string name18 = fName + fVa2->GetName() + "CM_vs_spill";
    fHa2CMs=new TH2F(name17.c_str(),name17.c_str(), 200, 0., 200.,
 		     fVa2->GetNbins(),
 		     fVa2->GetMin(),
 		     fVa2->GetMax());
    AddHistogram(fHa2CMs);
  
    // a2 average profile  
    std::string name21 = fName + "_ch_avg";
    fHavgamp=new TH1D(name21.c_str(),name21.c_str(),
		      fVch->GetNbins(),
		      fVch->GetMin(), fVch->GetMax());
    AddHistogram(fHavgamp);

    // a2 sigma profile
    std::string name22 = fName + "_ch_sigma";
    fHsigamp=new TH1D(name22.c_str(),name22.c_str(),
		      fVch->GetNbins(),
		      fVch->GetMin(), fVch->GetMax());
    AddHistogram(fHsigamp);

    // 2D a2 average profile
    std::string name23 = fName + "_2D_ch_avg";
    fHavgamp2D=new TH2D(name23.c_str(),name23.c_str(),18,0.,18.,
			72, 0., 72.);
    AddHistogram(fHavgamp2D);

    // 2D a2 sigma profile
    std::string name24 = fName + "_2D_ch_sigma";
    fHsigamp2D=new TH2D(name24.c_str(),name24.c_str(),18,0.,18.,
			72, 0., 72.);
    AddHistogram(fHsigamp2D);

    // a2CM average profile
    std::string name25 = fName + "_ch_CM_avg";
    fHavgampCM=new TH1D(name25.c_str(),name25.c_str(),
			 fVch->GetNbins(),
			fVch->GetMin(), fVch->GetMax());
    AddHistogram(fHavgampCM);

    // a2CM sigma profile
    std::string name26 = fName + "_ch_CM_sigma";
    fHsigampCM=new TH1D(name26.c_str(),name26.c_str(),
			 fVch->GetNbins(),
			fVch->GetMin(), fVch->GetMax());
    AddHistogram(fHsigampCM);

    // 2D a2CM average profile
    std::string name27= fName + "_2D_ch_CM_avg";
    fHavgampCM2D=new TH2D(name27.c_str(),name27.c_str(),18,0.,18.,
			  72, 0., 72.);
    AddHistogram(fHavgampCM2D);

    // 2D a2CM sigma profile
    std::string name28 = fName + "_2D_ch_CM_sigma";
    fHsigampCM2D=new TH2D(name28.c_str(),name28.c_str(),18,0.,18.,
			  72, 0., 72.);
    AddHistogram(fHsigampCM2D);

    // 2D event with a2 CM amplitude
    std::string name29 = fName + "_2D_ch_CM_event";
    fHevtampCM2D=new TH2D(name29.c_str(),name29.c_str(),18,0.,18.,
			  72, 0., 72.);
    AddHistogram(fHevtampCM2D);

    // 2D event with a2 CM amplitude with cut
    std::string name30 = fName + "_2D_ch_CM_cut_event";
    fHevtampCMcut2D=new TH2D(name30.c_str(),name30.c_str(),18,0.,18.,
			     72, 0., 72.);
    AddHistogram(fHevtampCMcut2D);

    // Channel amplitudes substraction correlation
    std::string name31 = fName + "_a2-a1_vs_a1-a0";
    fHampSubRatio = new TH2F(name31.c_str(),name31.c_str(),
			  200, -100., 100., 200, -100., 100.);
    AddHistogram(fHampSubRatio);
    fHampSubRatio->GetXaxis()->SetTitle("a1-a0");
    fHampSubRatio->GetYaxis()->SetTitle("a2-a1");

    // Channel amplitudes substraction ratio
    std::string name32 = fName + "_a2-a1_dividedby_a1-a0";
    fHampSubRap = new TH1F(name32.c_str(),name32.c_str(), 400, -3., 5.);
    AddHistogram(fHampSubRap);
    fHampSubRap->GetXaxis()->SetTitle("(a2-a1)/(a1-a0)");
    //-
    // Amplitude 1st sample
    std::string a0CMname = fName + fVa0->GetName() + "_CMcorr";
    fHa0CM=new TH1F(a0CMname.c_str(),a0CMname.c_str(),
		     fVa0->GetNbins(),
		     fVa0->GetMin(),
		     fVa0->GetMax()); 
    AddHistogram(fHa0CM);

    // Amplitude 2nd sample
    std::string a1CMname = fName + fVa1->GetName() + "_CMcorr";
    fHa1CM=new TH1F(a1CMname.c_str(),a1CMname.c_str(),
		    fVa1->GetNbins(),
		    fVa1->GetMin(),
		    fVa1->GetMax()); 
    AddHistogram(fHa1CM);

    // Amplitude 3rd sample
    std::string a2CMname = fName + fVa2->GetName() + "_CMcorr";
    fHa2CM=new TH1F(a2CMname.c_str(),a2CMname.c_str(),
		     fVa2->GetNbins(),
		     fVa2->GetMin(),
		     fVa2->GetMax()); 
    AddHistogram(fHa2CM);
  }

  //-
  if(tree) {
    tree->Branch(chitsname.c_str(),&fNclustKept,
		 chitsleavlist.c_str(),32000);
    tree->Branch(cchname.c_str(),fVcch->GetValues(),
		 cchleavlist.c_str(),32000);
    tree->Branch(ca2name.c_str(),fVca2->GetValues(),
		 ca2leavlist.c_str(),32000);
    tree->Branch(csname.c_str(),fVcs->GetValues(),
		 csleavlist.c_str(),32000);
    tree->Branch(name14.c_str(), ft_cor->GetValues(), t_cor_leavlist.c_str(), 32000);
  }  
  timing_initialised=false;

}

void PlaneSili::EndEvent(const CS::DaqEvent &event) {
 
  if (thr_flag) TThread::Lock();
  // fRateCounter++; //removed, done directly by Monitor.cc
  register bool rateflag = false;
  
  //Random numbers to remove artefacts
  register double ran0;
  register double ran1;
  register double ran2;

  if((fRateCounter&0x3f)==0) {
    rateflag = true;
    if (fHevtampCM2D) fHevtampCM2D->Reset("ICE");
    if (fHevtampCMcut2D) fHevtampCMcut2D->Reset("ICE");
  }
  

  // Get TCS phase
  try {
    tcsphase = event.GetTT().GetPhaseTCS() - 40.; // minus 40. to stay in agreement
                                                  // with previous definition of TCS
                                                  // phase in COOOL
  } catch (std::logic_error& e) {
    std::cerr << "PlaneSili::EndEvent: TCS phase not available: " << e.what() << std::endl;
    tcsphase = 0;
  }

  for (std::map<int,CDigit1>::const_iterator id = fDigits.begin(); 
       id!=fDigits.end(); id++) {

    const CDigit1& digit = id->second;

    int channel=digit.ch;
    int t = 0;  // no time information yet at this level
    int a0=(int)digit.dt[0];
    int a1=(int)digit.dt[1];
    int a2=(int)digit.dt[2];
    float a12 = (a2==0) ? 0 : float(a1)/float(a2);
    float a02 = (a2==0) ? 0 : float(a0)/float(a2);
    int cmc0 = (int)digit.dt[5];
    int cmc1 = (int)digit.dt[6];
    int cmc2 = (int)digit.dt[7];
    //    std::cout<<"Common mode correction "<<cmc0<<" "<<cmc1<<" "<<cmc2<<std::endl;

//     int channel=digit.ch;
//     int a0=(int)digit.dt[0];
//     int a1=(int)digit.dt[1];
//     int a2=(int)digit.dt[2];
//     int t =(int)digit.dt[3];

    if(fVch->Test(channel) &&
       fVa0->Test(a0) &&
       fVa1->Test(a1) &&
       fVa2->Test(a2)) {

      
      fVch->Store(channel);
      fVa0->Store(a0);
      fVa1->Store(a1);
      fVa2->Store(a2);
      fVt->Store(t);

      fHch->Fill(channel);
      fHa0->Fill(a0);
      fHa1->Fill(a1);
      fHa2->Fill(a2);
      if (fHa0s) fHa0s->Fill(event.GetBurstNumber(),a0);
      if (fHa1s) fHa1s->Fill(event.GetBurstNumber(),a1);
      if (fHa2s) fHa2s->Fill(event.GetBurstNumber(),a2);
      fHt->Fill(t);
      fChAmp->Fill(channel, a2);
      fHchvscmc2->Fill(channel, cmc2);
      fHchVsTis->Fill(channel, event.GetTT().GetTimeInSpill());
     
      //Fill random numbers to remove artefacts
      ran0=gRandom->Rndm(0)-1;
      ran1=gRandom->Rndm(1)-1;
      ran2=gRandom->Rndm(2)-1;

      fHtiming->Fill(float(a1+ran1)/float(a2+ran2), float(a0+ran0)/float(a2+ran2));

      if (fHampSubRatio) fHampSubRatio->Fill(a1-a0,a2-a1);
      if (fHampSubRap) if (a1-a0) fHampSubRap->Fill((a2-a1)*1./(a1-a0));
     
      fNhitsKept++;
      if (!(fRateCounter%25)) {
	fHoccup->Reset("ICE");
	fHoccup->Add(fHch, 1/double(fRateCounter));
      }

      fStat[channel].sum += a2;
      fStat[channel].sum2 += a2*a2;
      fStat[channel].nb++;
    }
  }
  
  fHhits->Fill(fNhitsKept);

  if (fUseCalib) {
    calcCommonModeCorrection(MAX_AMP_DEVIATION);
    for (int ch = 0; ch < fNchan; ch++) {
   
      // int apvnb = (ch & ~0x7f) >> 7;
      int apvnb =ch/128;
      //std::cerr << "apvnb" << apvnb << std::endl;
      float ped;
      if (calib_arr) 
        ped = calib_arr[ch].pedestal;
      else
        ped = 0;
      double a0CM = -(fVa0->GetValues()[ch] - cCMC_fCMapv[apvnb][0] - ped);
      double a1CM = -(fVa1->GetValues()[ch] - cCMC_fCMapv[apvnb][1] - ped);
      double a2CM = -(fVa2->GetValues()[ch] - cCMC_fCMapv[apvnb][2] - ped);

      if (fHa0CM) fHa0CM->Fill(a0CM);
      if (fHa1CM) fHa1CM->Fill(a1CM);
      if (fHa2CM) fHa2CM->Fill(a2CM);
      if (fHa0CMs) fHa0CMs->Fill(event.GetBurstNumber(),a0CM);
      if (fHa1CMs) fHa1CMs->Fill(event.GetBurstNumber(),a1CM);
      if (fHa2CMs) fHa2CMs->Fill(event.GetBurstNumber(),a2CM);

      fStatCM[ch].sum += a2CM;
      fStatCM[ch].sum2 += a2CM*a2CM;
      fStatCM[ch].nb++;
       if (rateflag) {
         register xy_s *xys = ch2xy(ch);
         if (fHevtampCM2D) fHevtampCM2D->Fill(xys->x, xys->y, a2CM);
         if (fHevtampCMcut2D) if (a2CM > SILAPV_A2_CUT) fHevtampCMcut2D->Fill(xys->x, xys->y, a2CM);
       }
   }
  }


   if (rateflag) {
     for (register int ii=0; ii<MAX_STAT; ii++) {
       register int nb = fStat[ii].nb;
       register int nbCM = fStat[ii].nb;
       if (nb == 0) continue;
       register xy_s *xys = ch2xy(ii);
       register double avg = fStat[ii].sum / nb;
       register double sigma = sqrt(fStat[ii].sum2/nb - fStat[ii].sum*fStat[ii].sum/nb/nb);

       if (fHavgamp) fHavgamp->SetBinContent(fHavgamp->GetXaxis()->FindBin(ii),avg);
       if (fHsigamp) fHsigamp->SetBinContent(fHsigamp->GetXaxis()->FindBin(ii),sigma);
       if (fHavgamp2D) fHavgamp2D->SetBinContent(fHavgamp2D->GetXaxis()->FindBin(xys->x),
                                 fHavgamp2D->GetYaxis()->FindBin(xys->y),avg);
       if (fHsigamp2D) fHsigamp2D->SetBinContent(fHsigamp2D->GetXaxis()->FindBin(xys->x),
                                 fHsigamp2D->GetYaxis()->FindBin(xys->y),sigma);
       if (nbCM) {
         register double avgCM = fStatCM[ii].sum / nbCM;
         register double sigmaCM = sqrt(fStatCM[ii].sum2/nb - fStatCM[ii].sum*fStatCM[ii].sum/nbCM/nbCM);
         if (fHavgampCM) fHavgampCM->SetBinContent(fHavgampCM->GetXaxis()->FindBin(ii),avgCM);
         if (fHsigampCM) fHsigampCM->SetBinContent(fHsigampCM->GetXaxis()->FindBin(ii),sigmaCM);
         if (fHavgampCM2D) fHavgampCM2D->SetBinContent(fHavgampCM2D->GetXaxis()->FindBin(xys->x),
                                   fHavgampCM2D->GetYaxis()->FindBin(xys->y),avgCM);
         if (fHsigampCM2D) fHsigampCM2D->SetBinContent(fHsigampCM2D->GetXaxis()->FindBin(xys->x),
                                   fHsigampCM2D->GetYaxis()->FindBin(xys->y),sigmaCM);
       }
     }
   }




  if (thr_flag) TThread::UnLock();
}


//-----------------------------------------------------------------------------


PlaneSili::xy_s *PlaneSili::ch2xy (int ch) {

  static xy_s xys;

  xys.x=-1; xys.y=-1;

  register short apvnb = (ch & ~0x7f) >> 7;
  register short chinapv = (ch & 0x7f);

  register short chincon = -1, connb = -1;

  register short apvnbc = apvnb & 0x3;
  register short cardnb = (apvnb & ~0x3) >> 2;
  switch (apvnbc) {
   case 3:
    if (chinapv<12) { connb = 6; chincon = chinapv + 36; break; }
    if (chinapv<60) { connb = 7; chincon = chinapv - 12; break; }
    if (chinapv<108) { connb = 8; chincon = chinapv - 60; break; }
    break;
   case 2:
    if (chinapv<24) { connb = 4; chincon = chinapv + 24; break; }
    if (chinapv<72) { connb = 5; chincon = chinapv - 24; break; }
    if (chinapv<108) { connb = 6; chincon = chinapv - 72; break; }
    break;
   case 1:
    if (chinapv<36) { connb = 2; chincon = chinapv + 12; break; }
    if (chinapv<84) { connb = 3; chincon = chinapv - 36; break; }
    if (chinapv<108) { connb = 4; chincon = chinapv - 84; break; }
    break;
   case 0:
    if (chinapv<48) { connb = 0; chincon = chinapv + 0; break; }
    if (chinapv<96) { connb = 1; chincon = chinapv - 48; break; }
    if (chinapv<108) { connb = 2; chincon = chinapv - 96; break; }
    break;
  }

  if (chincon >=0 && connb >=0 && fXY[chincon].x >=0 && fXY[chincon].x >=0) {
    xys.x = fXY[chincon].x + cardnb*6;
    xys.y = fXY[chincon].y + connb*8;
  }


  return &xys;
}


//-----------------------------------------------------------------------------


void PlaneSili::initconndata (void) {

  if (fXYinit) return;

  register int ii = 0;
  ii = 0; fXY[ii].x = 3; fXY[ii].y = 0;
  ii = 2; fXY[ii].x = 3; fXY[ii].y = 1;
  ii = 4; fXY[ii].x = 4; fXY[ii].y = 0;
  ii = 6; fXY[ii].x = 5; fXY[ii].y = 0;
  ii = 8; fXY[ii].x = 4; fXY[ii].y = 1;
  ii = 10; fXY[ii].x = 3; fXY[ii].y = 2;
  ii = 12; fXY[ii].x = 5; fXY[ii].y = 1;
  ii = 14; fXY[ii].x = 4; fXY[ii].y = 2;
  ii = 16; fXY[ii].x = 5; fXY[ii].y = 2;
  ii = 18; fXY[ii].x = 5; fXY[ii].y = 3;
  ii = 20; fXY[ii].x = 4; fXY[ii].y = 3;
  ii = 22; fXY[ii].x = 3; fXY[ii].y = 3;
  ii = 24; fXY[ii].x = 3; fXY[ii].y = 4;
  ii = 26; fXY[ii].x = 4; fXY[ii].y = 4;
  ii = 28; fXY[ii].x = 5; fXY[ii].y = 4;
  ii = 30; fXY[ii].x = 5; fXY[ii].y = 5;
  ii = 32; fXY[ii].x = 4; fXY[ii].y = 5;
  ii = 34; fXY[ii].x = 5; fXY[ii].y = 6;
  ii = 36; fXY[ii].x = 3; fXY[ii].y = 5;
  ii = 38; fXY[ii].x = 4; fXY[ii].y = 6;
  ii = 40; fXY[ii].x = 5; fXY[ii].y = 7;
  ii = 42; fXY[ii].x = 4; fXY[ii].y = 7;
  ii = 44; fXY[ii].x = 3; fXY[ii].y = 6;
  ii = 46; fXY[ii].x = 3; fXY[ii].y = 7;

  for ( register int jj = 1; jj < 48; jj+=2 ) {
    fXY[jj].x = 5-fXY[jj-1].x;
    fXY[jj].y = fXY[jj-1].y;
  }

  fXYinit = true;
}

//-----------------------------------------------------------------------------



void PlaneSili::Clusterize() {
  
  if (!fGeom) return;
  std::set<CCluster1>& clusters = fGeom->Clusterize(fDigits);
  
  if (thr_flag) TThread::Lock();
  int nclusters=0;
  typedef std::set<CCluster1>::iterator IC;
 
  for (IC i = clusters.begin(); i!=clusters.end(); i++) {
    float cpos = i->pos;
    float size = i->size;
    float amp2= /*i->dt[0] + i->dt[1] +*/ i->dt[2];
 
    float r0 = i->dt[3]; //a0/a2
    float r1 = i->dt[4]; //a1/a2
    fHcRatio0TCS->Fill(tcsphase, r0);
    fHcRatio1TCS->Fill(tcsphase, r1);
    
    if (fVca2->Test(amp2) &&
      	fVcch->Test(cpos) &&
      	fVcs->Test(size)) {
      
      fHcs->Fill(size);
      fHcch->Fill(cpos);
      fHctiming->Fill(i->dt[4], i->dt[3]);

      // Sergei's timing stuff (using calibration data)
      
      // Unpack time calibrations
      // ========================


      if (fUseCalib) {  // added by Damien because this part crashes when no calib are read  

      // range of validity of linear region of r0=a0/a2
      const float& r0_min = this->calib_time[0]; //[0]; 
      const float& r0_max = this->calib_time[1]; //[1];
      // range of validity of linear region of r1=a1/a2
      const float& r1_min = this->calib_time[2]; //[2];
      const float& r1_max = this->calib_time[3]; //[3];
      // T0s, taken from the fit of "time= f(r)" at r = 0 
      // (r is a0/a2 or a1/a2, time is track time from SciFis) 
      // Currently is common for all SI planes. 
      //To be tuned to compensate some global timing changes.
      const float& t0_0 = this->calib_time[4];
      const float& t0_1 = this->calib_time[5];
      // fitted slopes of linear regions of "time = f(r)" dependence
      const float& sl_0 = this->calib_time[6];
      const float& sl_1 = this->calib_time[7];
      // fine time corrections, obtained by 
      //"Track time - Calculated Time" distribution fit
      // (as fitted T0s are not precise enought).
      const float& tc_0 = this->calib_time[8];
      const float& tc_1 = this->calib_time[9];
      
      // const float TCS_T0 = 40.0;
      // Calibration constant (Maximum of TCS distr hist.)
      const float time_resol = 3.0;         
      // for weighted mean calculation. 
      // (has to be also subject for calibrations)
      float hit_time_0 =  1000; float hit_time_err_0 = 1.E+10;
      float hit_time_1 =  1000; float hit_time_err_1 = 1.E+10;
      float hit_time   =  1000; float hit_time_err   = 1.E+10;

      double tcs_cor = tcsphase;
      // -TCS_T0; compansated by tcsphaseshift in coool
    
      if(r0_min < r0 && r0 < r0_max)  {
	hit_time_0 = tcs_cor+t0_0 + sl_0*r0 - tc_0;
	hit_time_err_0 = time_resol;
      }
      if(r1_min < r1 && r1 < r1_max)  {
	hit_time_1 = tcs_cor+t0_1 + sl_1*r1 - tc_1;
	hit_time_err_1 = time_resol; 
      }
      
      float w0= 1./(hit_time_err_0* hit_time_err_0);
      float w1= 1./(hit_time_err_1* hit_time_err_1);
      hit_time = (hit_time_0*w0 + hit_time_1*w1)/(w0+w1);
      hit_time_err = sqrt(1./(w0+w1));

      fHcCalibTime->Fill(hit_time);
      // const_cast<CCluster1&>(*i).dt.push_back(hit_time);

      //
      // Final timing stuff (using calibration data)
      //
      // range of validity for r0=a0/a2

      if (this->calib_time.size()>25) {
	const float& fR0_min = this->calib_time[10];
	const float& fR0_max = this->calib_time[11];
	// range of validity for r1=a1/a2
	const float& fR1_min = this->calib_time[18];
	const float& fR1_max = this->calib_time[19];
	if (!timing_initialised) {
	  fPar0[0]=(double)this->calib_time[12];
	  fPar0[1]=(double)this->calib_time[13];
	  fPar0[2]=(double)this->calib_time[14];
	  fPar0[3]=(double)this->calib_time[15];
	  fPar0[4]=(double)this->calib_time[16];
	  fPar0[5]=(double)this->calib_time[17];
	  fPar1[0]=(double)this->calib_time[20];
	  fPar1[1]=(double)this->calib_time[21];
	  fPar1[2]=(double)this->calib_time[22];
	  fPar1[3]=(double)this->calib_time[23];
	  fPar1[4]=(double)this->calib_time[24];
	  fPar1[5]=(double)this->calib_time[25];
	  timing_initialised=true;
	}
	const float f_time_resol = 3.0;         
	float f_hit_time_0 =  1000; float f_hit_time_err_0 = 1.E+10;
	float f_hit_time_1 =  1000; float f_hit_time_err_1 = 1.E+10;
	float f_hit_time   =  1000; float f_hit_time_err   = 1.E+10;
	
	float f_tcs_cor = tcsphase;
	
	if(fR0_min < r0 && r0 < fR0_max)  {
	  f_hit_time_0 = -APV_t(r0, fPar0) + f_tcs_cor;
	  f_hit_time_err_0 = f_time_resol;
	}
	if(fR1_min < r1 && r1 < fR1_max)  {
	  f_hit_time_1 = -APV_t(r1, fPar1) + f_tcs_cor;
	  if (r1>1.0) f_hit_time_err_1 = 25.0;
	  else        f_hit_time_err_1 =  5.0;
	}
	
	float f_w0= 1./(f_hit_time_err_0* f_hit_time_err_0);
	float f_w1= 1./(f_hit_time_err_1* f_hit_time_err_1);
	f_hit_time = (f_hit_time_0*f_w0 + f_hit_time_1*f_w1)/(f_w0+f_w1);
	f_hit_time_err = sqrt(1./(f_w0+f_w1));
	
	fHctime ->Fill(f_hit_time   - f_tcs_cor);
	if (fHctime1) fHctime1->Fill(f_hit_time_0 - f_tcs_cor);
	if (fHctime2) fHctime2->Fill(f_hit_time_1 - f_tcs_cor);
	fHctimeCorr->Fill(f_hit_time_0 - f_tcs_cor, f_hit_time_1 - f_tcs_cor);
	fHctimeCorrClusterTCS->Fill(tcsphase, f_hit_time-f_tcs_cor);

	fHcRightTime->Fill(f_hit_time);
	ft_cor->Store(f_hit_time);
	if (fabs(f_hit_time)<10) {
	  fHca2->Fill(amp2);
	  fHAmpTCStime->Fill(amp2, tcsphase);;
	}
	fHcTiTi->Fill(hit_time, f_hit_time);
	
	const_cast<CCluster1&>(*i).dt.push_back(f_hit_time_0);
	const_cast<CCluster1&>(*i).dt.push_back(f_hit_time_1);
	const_cast<CCluster1&>(*i).dt.push_back(f_hit_time);
	
      }
      }  // end of if (usecalib)
      fVcch->Store(cpos);
      fVca2->Store(amp2);
      
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

//-----------------------------------------------------------------------------
// Calculate common mode noise correction
// (from gemMonitor written by B. Ketzer)
//-----------------------------------------------------------------------------
void PlaneSili::calcCommonModeCorrection(double cut) {
  
  // Local array for sorted amplitudes of active channels
  typedef struct { double amp[MAX_APV_NBCH]; } ampS_s;
  ampS_s ampSort[3];
  
  for (register int iapv = 0; iapv < fNchan/128; iapv++) {

    // Loop over channels
    int nb_activech = 0;
    register int iactch = 0;
    for ( register int ich = 0; ich < MAX_APV_NBCH-1; ich++ ) {
      register int ch = ich + iapv*MAX_APV_NBCH;

      // calib_arr     = new calib_s[fNchan];

      // Consider only active channels

      if (calib_arr && (calib_arr[ch].flag == 1) ) {
	
        // Copy active channels into new array for sorting
        ampSort[0].amp[iactch] = fVa0->GetValues()[ch] - calib_arr[ch].pedestal;
        ampSort[1].amp[iactch] = fVa1->GetValues()[ch] - calib_arr[ch].pedestal;
        ampSort[2].amp[iactch] = fVa2->GetValues()[ch] - calib_arr[ch].pedestal;
//	std::cerr << "0: " << ampSort[0].amp[iactch] << " 1: " << ampSort[1].amp[iactch] << "2: " << ampSort[2].amp[iactch] << endl;
        iactch++;
      }
      nb_activech = iactch;
 
      for (register int ismp=0; ismp<3; ismp++) {
	
	// Sort in ascending order
	double *ampEnd = ampSort[ismp].amp + nb_activech;
	std::sort(ampSort[ismp].amp, ampEnd);
	//	copy(ampSort, ampEnd, ostream_iterator<double>(cout, " "));
 
	// Calculate median
	register double median;
	register int ctr;
	if ( nb_activech%2 == 1 ) {
	  ctr = (nb_activech-1)/2;
	  median = ampSort[ismp].amp[ctr];

	}
	else {

	  ctr = nb_activech/2;
	  median = (ampSort[ismp].amp[ctr-1]+ampSort[ismp].amp[ctr])/2.;
	}
	
	// Loop over active channels 
	register double ampMean = 0;
	register int n = 0;
	//	cerr << "nb_activech: " << nb_activech  << endl;
	for ( register int i=0; i<nb_activech; i++ ) {

	  //cerr << "median:"  << median-ampSort[ismp].amp[i]  << " cut:"  << cut
	  //   << "ampSort:"  <<  ampSort[ismp].amp[i]-ampSort[ismp].amp[ctr+i-n]  << endl;
	  // Sum up channels with amplitudes around median
	  if ( ((median-ampSort[ismp].amp[i]) < cut) && 
	       ((ampSort[ismp].amp[i]-ampSort[ismp].amp[ctr+i-n]) < cut) ) {
	    ampMean += ampSort[ismp].amp[i];
	    n++;
	  }
	}	
	// Calculate mean of channels without signal for CM correction
	cCMC_fCMapv[iapv][ismp] = ampMean / n;
      }
    }
  }
}
 

//-----------------------------------------------------------------------------


  
#if USE_DATABASE == 1
void PlaneSili::ReadCalib(const tm &t)
  {
    // read-in corresponding calibration constants
    try{
//      std::cout<<"Trying to read from data base for " <<fName<<std::endl;
      ReadFromDataBase(calib_data,t);
//      std::cout<<"PlaneSili::ReadCalib() ==> "
//	       <<this->GetName()<<" calibrations are found !"<<std::endl;
      
      if(calib_data.size() <  (unsigned) fNchan) {
	//    std::cerr<<"Size of Calibration File is not correct ! Should be : "
	//          <<fNchan<<" Is "<<calib_data.size()<<" "
	//          <<t.tm_mday<<"."<<t.tm_mon+1<<"."<<t.tm_year+1900<<" "
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
    catch(...) {      std::cout<<"PlaneSili::ReadCalib() ==> "<<GetName()
	       <<" calibrations, valid for ";
      std::cout<<t.tm_mday<<"."<<t.tm_mon+1<<"."<<t.tm_year+1900<<" "
	       <<t.tm_hour<<":"<<t.tm_min<<":"<<t.tm_sec
	       <<", not found in DB"<<std::endl;
    }
    
    // Read calibration data for conversion of amplitude->time for each plane
    if (fDataBase==NULL)
      throw CS::Exception("PlaneSili::ReadCalib():  data base is not opened.");
    try{
      std::string tcaln("do_timing_calib/");
      tcaln += GetName();
      ifstream tcalf(tcaln.c_str());
      std::string strdata("");
      if (!tcalf) {
	struct tm tt(t);
	CDB::Time tp(mktime(&tt),0);
	fDataBase->read(fName,strdata,tp,"timing");
	if (strdata == "") {
	  std::cerr << "PlaneSili::ReadCalib() "<<GetName()
		    <<" : no timing calibration file found"<<std::endl;
	  return;
	}
      }
      else {
	std::string tcalf_;
	while (tcalf >> tcalf_) strdata += tcalf_+" ";
	std::cerr << "READING FROM FILE!!!!" << std::endl;
	std::cerr << "^^^^^^^^^^^^^^^^^^^^^" << std::endl;
	std::cerr << "Read:" << strdata.c_str() << std::endl;
      }
      //    std::cout << "Reading timing calibrations for " << GetName() << std::endl;
      std::istringstream istrdata(strdata.c_str());
      //    std::cout << strdata << std::endl;
      istrdata >> calib_time;
    }
    catch(CS::Exception& e) {
      std::cout<<"rethrowing"<<std::endl;
      throw;
    }
  
    if (calib_time.size()!=0) {
      
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
    } else {
      std::cout<<"PlaneSili::ReadCalib() ==> "<<GetName()
	       <<" timing calibrations, valid for ";
      std::cout<<t.tm_mday<<"."<<t.tm_mon+1<<"."<<t.tm_year+1900<<" "
	       <<t.tm_hour<<":"<<t.tm_min<<":"<<t.tm_sec
	       <<", not found in DB"<<std::endl;
    }
    //-
    if (fUseCalib) {
      
      calib_arr     = new calib_s[fNchan];
      
      // need apvRich.xml map to sort calibration data
      typedef CS::Chip::Maps::const_iterator m_it;
      typedef std::vector<ChannelCalib>::const_iterator vec_APVcal;
      
      // assume same source id for all apvs/channels from one RiAPV
      unsigned int srcId = 0;
      for (vec_APVcal ical = calib_data.begin(); ical != (calib_data.end()); ical++) {
	const ChannelCalib& cal0 = *ical;
	if (&(cal0.calib_data[0])) {
        const ChannelCalib& ch0 = cal0.calib_data[0];
        srcId = ch0.srcId;
        break;
	}
      }
      
      if (srcId) {
	m_it m_low = fRunMaps->lower_bound(uint64(srcId)<<48);
	// fRunMaps->Print(std::cout, "");
	for( m_it cc=m_low; ((*cc).first>>48)==srcId; cc++ ) {
	  const CS::ChipAPV::Digit *m =
	    dynamic_cast<const CS::ChipAPV::Digit*>((*cc).second);
	  if( m==NULL ) {
	    std::cerr << "ChipAPV wrong map.\n";
	    continue;
	  }
	  if (strcmp(m->GetDetID().GetName().c_str(), GetName()))
	    continue;
	  
	  register int strip = m->GetChannel();
	  CS::ChipAPV::DataID dataid(m->GetDataID());
	  // there are 2 different adcId for a given APVRich TBname 
	  const ChannelCalib& a_cal = calib_data[m->GetChip()+6*((dataid.u.s.adc_id - 1)&1)-1];
	  if (calib_arr) {
            if (&(a_cal.calib_data[0])) {
	      const ChannelCalib a_ch = a_cal.calib_data[m->GetChipChannel()];
	      calib_arr[strip].flag     = a_ch.flag;
	      calib_arr[strip].pedestal = a_ch.ped;
	      calib_arr[strip].sigma    = a_ch.sigma;
	    } else {
	      calib_arr[strip].flag     = 0;
	      calib_arr[strip].pedestal = 0;
	      calib_arr[strip].sigma    = 0;
	    }
          }
	}
	if (getenv("COND_DB_DEBUG")!=0 || getenv("COND_DB_DEBUG_RA")!=0) {
	  std::cout << "Calibration values for APV "<<fName<<":\n";
	  std::cout << "  ch  flag   ped   sigma\n";
	  if (calib_arr) for (register int ii=0; ii<fNchan; ii++) {
	    std::cout<<"  "<<ii<<"  "<<calib_arr[ii].flag<<"   "<<calib_arr[ii].pedestal<<"   "<<calib_arr[ii].sigma<<std::endl;
	  }
	}
      }
    }
  }


















#endif //USE_DATABASE



