#include <cmath>

#include "PlaneECAL2.h"
#include <exception>
#include "TProfile.h"
#include "ChipSADC.h"
#include "PlaneECAL2Panel.h"
#include <vector>
#include <stdexcept>
using namespace std;

#define nSADC_Sampl 32       // number of SADC samples
#define Dpedstart  0         // Phys.ev. pedestal window
#define Dpedend    5
#define Dpulsstart 6         // Phys.ev. puls window
#define Dpulsend   31
#define Dcalpedstart 0       // LED pedestal window
#define Dcalpedend  7
#define Dledstart   8        // LED puls  window
#define Dledend     31
int SADCoverflow=4092;   // overflow for mSADC channels

ClassImp(PlaneECAL2);


PlaneECAL2::PlaneECAL2(const char *detname,int ncol, int nrow, int center, int width)
  : Plane(detname),fNrows(nrow),fNcols(ncol),fNchan(nrow*ncol), hLed(0), hLedChnl(0),
  hoverflow(0), hxynoise(0), hampnoise(0), hrnd(0), hledlen(0), hLedTm(0), hPhysTm(0),
  hmnled(0), hmnledrms(0), hmnledamp(0), hmnledr(0), hmnped(0), hmnpedr(0),hTimech(0),
    hTimechprof(0),hLEDTimech(0), hLedProf(0), hPedProf(0), href(0), ledcalibCFD(0) {

  nSADCchan=nb_ECAL2_SADCcnl;
  fNsamples=nSADC_Sampl;
  fNSADChits=0;
//
  SADCampcut = 20; SADCLEDcut = 10;  // SADC  amplitude/LED  cuts
//
  href=NULL, hprofx=NULL, hprofy=NULL, hrefXY=NULL; hxynoisecut=NULL; hampnoisecut=NULL;

  wwwled_cent=0.;
  wwwled_perif=0.;
 
  

  std::string name;
  pedstart=Dpedstart;
  pedend=Dpedend;
  name=fName+"_Phys_ped_window";
  fVped=AddVariable(name.c_str(),100,pedstart,pedend,5);
//
  pulsstart=Dpulsstart;
  pulsend=Dpulsend;
  name=fName+"_Phys_puls_window";
  fVpuls=AddVariable(name.c_str(),100,pulsstart,pulsend,5);
//
  calpedstart=Dcalpedstart;
  calpedend=Dcalpedend;
  name=fName+"_LED_ped_window";
  fVledped=AddVariable(name.c_str(),100,calpedstart,calpedend,5);
//
  ledstart=Dledstart;
  ledend=Dledend;
  name=fName+"_LED_puls_window";
  fVled=AddVariable(name.c_str(),100,ledstart,ledend,5);
//
  name=fName+"_SADC_amp/led_cuts";
  fVscuts=AddVariable(name.c_str(),100,SADCampcut,SADCLEDcut,5);

  std::string rowname=fName+"_row";
  std::string colname=fName+"_col";
  std::string ampname=fName+"_amp";

  fVrow=AddVariable(rowname.c_str(),fNrows,0,fNrows ,nSADCchan);
  fVcol=AddVariable(colname.c_str(),fNcols,0,fNcols ,nSADCchan);
  fVamp=AddVariable(ampname.c_str(),100,center-width,center+width,nSADCchan);
  
}


void PlaneECAL2::Init(TTree* tree) {
  std::string myname = "EC02_";
  std::string name, title;
  char scutamp[10],scutled[10];
  sprintf(scutamp,"%3.1i",SADCampcut);
  sprintf(scutled,"%3.1i",SADCLEDcut);
   fTcalib.resize(3072);

  name = myname + "SADC_hits#";
  hSADChit=new TH1F_Ref(name.c_str(),name.c_str(),300,0.,300., fRateCounter);
  ((TH1F_Ref*)hSADChit)->SetReference(fReferenceDirectory);
  AddHistogram(hSADChit);

  name = myname + "SADC_hits#(amp>thr)";
  title = myname + "SADC_hits#,amp>"+scutamp;
  hSADCcuthits=new TH1F_Ref(name.c_str(),title.c_str(),200,0.,200.,fRateCounter);
  ((TH1F_Ref*)hSADCcuthits)->SetReference(fReferenceDirectory);
  AddHistogram(hSADCcuthits);

  name = myname+"Nhits_vs_XY";
  title = fName+" Number of hits vs XY";
  fHrc1=new TH2F(name.c_str(),title.c_str(),fNcols,0.,fNcols,fNrows,0.,fNrows);
  fHrc1->SetOption("col");
  AddHistogram(fHrc1);

  name = myname+"Nhits_vs_XY(amp>"+scutamp+")";
  title = fName+" Number of hits vs XY, amp>"+scutamp;
  hnhitcut=new TH2F(name.c_str(),title.c_str(),fNcols,0.,fNcols,fNrows,0.,fNrows);
  hnhitcut->SetOption("col");
  AddHistogram(hnhitcut);

  name = myname + "Ampl_vs_XY";
  title = fName + " Cell energy depositions";
  fHrca=new TH2F(name.c_str(),title.c_str(),fNcols,0.,fNcols,fNrows,0.,fNrows);
  fHrca->SetOption("col");
  AddHistogram(fHrca);

  name = myname + "Ampl_vs_XY(amp>"+scutamp+")";
  title = myname + "Amplitude vs XY, amp>"+scutamp;
  fHxyac = new TH2F(name.c_str(), title.c_str(),fNcols, 0, fNcols, fNrows, 0, fNrows);
  fHxyac->SetOption("col");
  AddHistogram(fHxyac);

  name = myname + "Ampl_vs_chnl#";
  fHavsadr = new TH2F(name.c_str(), name.c_str(),
		fNrows*fNcols, 0, fNrows*fNcols,
		200, 0.,4100.);
  fHavsadr->SetOption("col");
  AddHistogram(fHavsadr);

  name = myname + "SADC_amp";
  title = fName + "SADC amplitude";
  fHa=new TH1F_Ref(name.c_str(),title.c_str(),
	       410,0.,4100., fRateCounter);
  ((TH1F_Ref*)fHa)->SetReference(fReferenceDirectory);
  AddHistogram(fHa);

  name = myname + "SADC_SumAmpl";
  title = fName + " SADC Sum Amplitude";
  hSADCsum=new TH1F_Ref(name.c_str(),title.c_str(),200,0.,20000., fRateCounter);
  ((TH1F_Ref*)hSADCsum)->SetReference(fReferenceDirectory);
  AddHistogram(hSADCsum);

  name = myname + "SADC_SumAmp(amp>thr)";
  title = fName + " SADC Sum amplitude, cell amp>"+scutamp;
  hSADCcutsum=new TH1F_Ref(name.c_str(),title.c_str(),200,0.,20000., fRateCounter);
  ((TH1F_Ref*)hSADCcutsum)->SetReference(fReferenceDirectory);
  AddHistogram(hSADCcutsum);

//  name = myname + "Ampl_vs_chnl#(amp>thr)";
//  title = fName + "Cell energy deposition vs chnl# (cut ampl>thr)";
//  fHchac = new TH2F(name.c_str(), title.c_str(),
//		    fNrows*fNcols, 0, fNrows*fNcols,
//		    fVamp->GetNbins(), fVamp->GetMin(),
//		    fVamp->GetMax());
//  AddHistogram(fHchac);

  name = myname + "SADC_peds";
  title = fName + " SADC pedestals";
  hPed=new TH1F_Ref(name.c_str(),title.c_str(), 110, 0., 440., fRateCounter);
//  ((TH1F_Ref*)hPed)->SetReference(fReferenceDirectory);
  AddHistogram(hPed);

//  name = myname + "SADC_pedRMS";
//  title = fName + " SADC pedRMS (counts)";
//  hPedRMS=new TH1F_Ref(name.c_str(),title.c_str(), 100, 0., 25., fRateCounter);
//  ((TH1F_Ref*)hPedRMS)->SetReference(fReferenceDirectory);
//  AddHistogram(hPedRMS);

  name = myname + "Peds_vs_XY";
  title = fName + " Peds vs XY";
  hPedXY=new TH2F(name.c_str(),title.c_str(),fNcols,0.,fNcols,fNrows,0.,fNrows);
  hPedXY->SetOption("col");
  AddHistogram(hPedXY);

  name = myname + "Peds/RMS_vs_chnl#";
  title = fName + " Pedestal/RMS vs chnl#";
  hPedProf=new TProfile(name.c_str(),title.c_str(),fNrows*fNcols,0.,fNrows*fNcols,"s");
  AddHistogram(hPedProf);

  name = myname + "Led_ev_length";
  title = fName + "Led event length";
  hledlen=new TH1F_Ref(name.c_str(),title.c_str(),310,0,3100, fRateCounter);
  ((TH1F_Ref*)hledlen)->SetReference(fReferenceDirectory);
  AddHistogram(hledlen);

  name = myname + "Leds_vs_XY";
  title = fName +  "Leds, Led>"+scutled;
  hLed=new TH2F(name.c_str(),title.c_str(),fNcols,0.,fNcols,fNrows,0.,fNrows);
  hLed->SetOption("col");
  AddHistogram(hLed);

  
  name = myname + "LedReference";
  title = fName + "Led deviation from reference";
  href= new TProfile(name.c_str(), title.c_str(), fNchan,0.,fNchan,"s");
  AddHistogram(href);

  name = myname+"BadLedRef_vs_XY";
  title = fName+"Bad Ref.Led vs XY";
  hrefXY=new TH2F(name.c_str(),title.c_str(),fNcols,0.,fNcols,fNrows,0.,fNrows);
  hrefXY->SetOption("col");
  AddHistogram(hrefXY);

  name = myname + "Led/RMS_vs_chnl#";
  title = fName + "LEDS/RMS vs channel#";
  hLedProf=new TProfile(name.c_str(),title.c_str(), fNrows*fNcols, 0., fNrows*fNcols,"s");
  AddHistogram(hLedProf);

  name = myname + "X_hitprofile";
  title = fName + " hit X-profile";
  hprofx=new TH1F_Ref(name.c_str(),title.c_str(), fNcols, 0.,fNcols, fRateCounter);
  ((TH1F_Ref*)hprofx)->SetReference(fReferenceDirectory);
  AddHistogram(hprofx);

  name = myname + "Y_hitprofile";
  title = fName + " hit Y-profile";
  hprofy=new TH1F_Ref(name.c_str(),title.c_str(), fNrows, 0.,fNrows, fRateCounter);
  ((TH1F_Ref*)hprofy)->SetReference(fReferenceDirectory);
  AddHistogram(hprofy);

  name = myname + "MnLed_vs_chnl#";
  title = fName + "Mean LED vs channel#";
  hmnled=new TH1F_Ref(name.c_str(),title.c_str(),
                              fNrows*fNcols, 0., fNrows*fNcols, fRateCounter);
  ((TH1F_Ref*)hmnled)->SetReference(fReferenceDirectory);
  AddHistogram(hmnled);
//
    name = myname + "MeanPedRMS";
    title = fName + " SADC mean  pedestal RMS";
    hmnpedr=new TH1F_Ref(name.c_str(),title.c_str(), 100, 0., 20., fRateCounter);
//    ((TH1F_Ref*)hmnpedr)->SetReference(fReferenceDirectory);
    AddHistogram(hmnpedr);
    
    name    = myname + "Time_vs_chnl";
    title   = fName  + "Time vs channel#";
    hTimech = new TH2F(name.c_str(),title.c_str(),fNchan,0,fNchan,700,-350,350);
    hTimech->SetOption("col");
    AddHistogram(hTimech);
    
    name    = myname + "Time_vs_chnlprof";
    title   = fName  + "Time vs channelprof#";
    hTimechprof = new TProfile2D(name.c_str(),title.c_str(),64,0,64,48,0,48);
    hTimechprof->SetOption("colz");
    AddHistogram(hTimechprof);

    name    = myname + "LED_Time_chnl";
    title   = fName  + "LED Timechannel#";
    hLEDTimech = new TH2F(name.c_str(),title.c_str(),fNchan,0,fNchan,400,-50,350);
    hLEDTimech->SetOption("col");
    AddHistogram(hLEDTimech);
    



  if (fExpertHistos) {

  name = myname + "Nhits_vs_chnl#";
  fHch=new TH1F_Ref(name.c_str(),name.c_str(),fNrows*fNcols,0,fNrows*fNcols,fRateCounter);
  ((TH1F_Ref*)fHch)->SetReference(fReferenceDirectory);
  AddHistogram(fHch);

    name = myname + "Leds";
    title = fName + "Mean Led amplitudes";
    hmnledamp=new TH1F_Ref(name.c_str(),title.c_str(),400,0.,4200., fRateCounter);
    ((TH1F_Ref*)hmnledamp)->SetReference(fReferenceDirectory);
    AddHistogram(hmnledamp);

    name = myname + "LedsRMS";
    title = fName + " Led RMSs(%),Led>"+scutled;
    hmnledr=new TH1F_Ref(name.c_str(),title.c_str(),400,0.,100., fRateCounter);
    ((TH1F_Ref*)hmnledr)->SetReference(fReferenceDirectory);
    AddHistogram(hmnledr);

    name = myname + "Led_vs_chnl#";
    title = fName + "LED vs channel#";
    hLedChnl=new TH2F(name.c_str(),title.c_str(),
                              fNrows*fNcols, 0., fNrows*fNcols,200,0.,4100.);
    AddHistogram(hLedChnl);

    name = myname + "overflow";
    title = fName + "ADC overfloved channels";
    hoverflow=new TH2F(name.c_str(),title.c_str(),fNcols,0.,fNcols,fNrows,0.,fNrows);
    hoverflow->SetOption("col");
    AddHistogram(hoverflow);

    name = myname + "Random_ev_length";
    title = fName + "Random trigger event length";
    hrnd=new TH1F_Ref(name.c_str(),title.c_str(),100,0,300, fRateCounter);
    ((TH1F_Ref*)hrnd)->SetReference(fReferenceDirectory);
    AddHistogram(hrnd);

    name = myname + "Nhits_vs_XY_rnd_trig";
    title = fName + " hits vs XY, random trigger";
    hxynoise = new TH2F(name.c_str(), title.c_str(),fNcols, 0, fNcols, fNrows, 0, fNrows);
    hxynoise->SetOption("col");
    AddHistogram(hxynoise);

    name = myname + "Nhits_vs_XY_rnd_trig_cut";
    title = fName + " hits vs XY, random trigger, amp cut";
    hxynoisecut = new TH2F(name.c_str(), title.c_str(),fNcols, 0, fNcols, fNrows, 0, fNrows);
    hxynoisecut->SetOption("col");
    AddHistogram(hxynoisecut);

    name = myname + "Amp_vs_XY_rnd_trig";
    hampnoise = new TH2F(name.c_str(),name.c_str(),fNcols, 0, fNcols, fNrows, 0, fNrows);
    hampnoise->SetOption("col");
    AddHistogram(hampnoise);

    name = myname + "Amp_vs_XY_rnd_trig_cut";
    hampnoisecut = new TH2F(name.c_str(),name.c_str(),fNcols, 0, fNcols, fNrows, 0, fNrows);
    hampnoisecut->SetOption("col");
    AddHistogram(hampnoisecut);

    name = myname + "SADC_Led_timing";
    title = fName + "SADC Led max.ampl vs sample#,amp>"+scutled;
    hLedTm=new TH1F_Ref(name.c_str(),title.c_str(),32,0.,32., fRateCounter);
//    ((TH1F_Ref*)hLedTm)->SetReference(fReferenceDirectory);
    AddHistogram(hLedTm);

    name = myname + "SADC_Phys_ev_timing";
    title = fName + "SADC,Phys ev. max.ampl vs sample#,amp>"+scutamp;
    hPhysTm=new TH1F_Ref(name.c_str(),title.c_str(),32,0.,32., fRateCounter);
//    ((TH1F_Ref*)hPhysTm)->SetReference(fReferenceDirectory);
    AddHistogram(hPhysTm);

    name = myname + "MeanPed_vs_chnl#";
    title = fName + "Mean pedestal vs channel#";
    hmnped=new TH1F_Ref(name.c_str(),title.c_str(),
                              fNrows*fNcols, 0., fNrows*fNcols, fRateCounter);
    ((TH1F_Ref*)hmnped)->SetReference(fReferenceDirectory);
    AddHistogram(hmnped);

    name = myname + "MnLedRMS_vs_chnl#";
    title = fName + "Mean LED RMS(%) vs channel#, Led>"+scutled;
    hmnledrms=new TH1F_Ref(name.c_str(),title.c_str(),
                              fNrows*fNcols, 0., fNrows*fNcols, fRateCounter);
    ((TH1F_Ref*)hmnledrms)->SetReference(fReferenceDirectory);
    AddHistogram(hmnledrms);

    name = myname + "LedSADC_samples";
    hSADCsmpl=new TH1F(name.c_str(),name.c_str(),32,0,32);
    AddHistogram(hSADCsmpl);
  }

  std::string hitsname = fName + "_hits";
  std::string hitsleavlist = hitsname + "/I";
  std::string leavlist0 =fVrow->GetName()  + "[" + hitsname + "]/F";
  std::string leavlist1 =fVcol->GetName()  + "[" + hitsname + "]/F";
  std::string leavlist2 =fVamp->GetName()  + "[" + hitsname + "]/F";

  if(tree) {
    fIsInTree = true;
    tree->Branch(hitsname.c_str(),&fNhitsKept,
		 hitsleavlist.c_str(),32000);
    tree->Branch(fVrow->GetName().c_str(),fVrow->GetValues(),
		 leavlist0.c_str(),32000);
    tree->Branch(fVcol->GetName().c_str(),fVcol->GetValues(),
		 leavlist1.c_str(),32000);
    tree->Branch(fVamp->GetName().c_str(),fVamp->GetValues(),
		 leavlist2.c_str(),32000);

  }

  if(ledcalibCFD==NULL)
    ledcalibCFD=new vector<Float_t>(3072,0);
  
}


//PlaneECAL2::~PlaneECAL2() {}

void PlaneECAL2::StoreDigit(int col, int row, std::vector<float>& data) {
//
  int sadclen=data.size()-2;
  if(sadclen<1) return;
//
  int nhit=fNSADChits;
  if (nhit< nSADCchan) {
    fSADCRow[nhit]=row;
    fSADCCol[nhit]=col;
    int len=fNsamples;
    if(sadclen<len){ len=sadclen;
                    for(int i=len;i<fNsamples;i++) fSADCev[nhit][i]=0; }
    for(int i=0;i<len;i++) fSADCev[nhit][i]=(int)data[i+2];
    fNhits++;
    fNSADChits++;  

  }
}


void PlaneECAL2::StoreDigit(int col, int row, int amp) {
}


void PlaneECAL2::StoreDigit(CS::Chip::Digit* digit) {
  if(fNhits==0){ fNSADChits=0;}        // clear
//
  std::vector<float> data=digit->GetNtupleData();
  if(data.size()<3) return;
  const CS::ChipSADC::Digit* sadcdig = dynamic_cast<const CS::ChipSADC::Digit*>(digit);
  if(sadcdig != NULL ) {
         int ix=sadcdig->GetX(); int iy=sadcdig->GetY();
         if(ix < 0||iy<0||ix>63||iy>47) return;
         this->StoreDigit(ix,iy,data);
  } else {
          int ix=(int) data[0]; int iy=(int) data[1];
//          if(ix < 0||iy<8||ix>63||iy>39) return;
          if(ix < 0||iy<0||ix>63||iy>47) return;
          this->StoreDigit(ix,iy,(int) data[2]);
  }
}


void PlaneECAL2::EndEvent(const CS::DaqEvent &event) {
  int evtype=event.GetType();
  if(!(evtype==7||evtype==8))  return;
  int trmsk=event.GetTrigger();
  const CS::uint32 Atr = event.GetHeader().GetTypeAttributes()[0]&0xf8000000;
  if(evtype==8&&Atr!=0x68000000) return;      // select only calorimeters calib.ev.
  if (thr_flag) TThread::Lock();

  int rndtr=trmsk&0x800;                          // random trigger
  if(hrnd && evtype==7&&rndtr)  hrnd->Fill(fNhits);
  if(hledlen && evtype==8) hledlen->Fill(fNhits);
//
//==================================== SADC ===================================
//
  double sum=0., time=0.;
  double SADCsum=0., SADCcutsum=0.;
  int SADCcuthits=0;
  if(evtype==8) {  hPedXY->Reset(); hLed->Reset();
                   wwwled_perif=0.; wwwled_cent=0.;
  }      
  // Get TCS phase
  double TCSphase;
  try {
    TCSphase = event.GetTT().GetPhaseTCS();
  } catch (std::logic_error& e) {
    std::cerr << "PlaneECAL2::EndEvent: TCS phase not available: " << e.what() << std::endl;
    TCSphase = 0;
  }                            
  
  for ( int i=0; i<fNSADChits; i++) {            // SADC
    int col=fSADCCol[i];
    int row=fSADCRow[i];
    unsigned int adr=row + col*fNrows; 
//
    double RMS=0.,ped=0.,puls=0.;
    if(evtype==8) {                                                 // Leds
         ped=SADCped(calpedstart,calpedend,&fSADCev[i][0],RMS);
         puls=SADCpuls(ledstart,ledend,&fSADCev[i][0], ped,sum,time);
         if(puls<0.) puls=0.;
         hPed->Fill(ped);
//         hPedRMS->Fill(RMS);
         hPedProf->Fill(adr,ped);
         hPedXY->Fill(col,row,ped);
         hLedProf->Fill(adr,puls);
         //Fill timing histogramm with an additional offset to the official timing due to some difference in algorithm
         double cfdtime=-999;
         cfdtime=CFDtime(pulsstart,pulsend,&fSADCev[i][0]);
         if(cfdtime!=-999 && (int)adr<fNchan && (ledcalibCFD->size()>0)&& (adr<ledcalibCFD->size())) {
           hLEDTimech->Fill(adr,(12.86*cfdtime)-ledcalibCFD->at(adr));
         }
         if(col>15&&col<48) wwwled_cent+=puls; else wwwled_perif+=puls;
         if(puls>SADCLEDcut)  hLed->Fill(col,row,puls);
         if (fExpertHistos) {
           hLedChnl->Fill(adr,puls);
           if(puls>SADCLEDcut) hLedTm->Fill(time);
         }
         continue;
    }
//
    if (evtype==7&&hoverflow) for(int j=0;j<fNsamples;j++){
       if(fSADCev[i][j]>SADCoverflow) hoverflow->Fill(col,row);
    }
//
    ped=SADCped(pedstart, pedend, &fSADCev[i][0],RMS);
    puls=SADCpuls(pulsstart, pulsend,&fSADCev[i][0], ped,sum,time);
    if(puls<0.) puls=0.;
    if(evtype==7&&rndtr)  {   //random trigger
      if(fExpertHistos )  {   
         hxynoise->Fill(col,row);
         hampnoise->Fill(col,row,puls);
         if(puls>SADCampcut){
            hxynoisecut->Fill(col,row);
            hampnoisecut->Fill(col,row,puls);
	 } 
      }
         continue;
    }
    fHa->Fill(puls); 
    fHavsadr->Fill(adr,puls);
    fHrca->Fill(col,row,puls);
    fHrc1->Fill(col,row);
    
        //Fill timing histogramm with an additional offset to the official timing due to some difference in algorithm
    double cfdtime=-999;
    //cfdtime=CFDtime(pulsstart,pulsend,&fSADCev[i][0])-1;
    cfdtime=Hmaxtime(pulsstart,pulsend,&fSADCev[i][0]);
    if(cfdtime!=-999 && adr<(unsigned)fNchan) {
        hTimech->Fill(adr,(12.86*cfdtime+TCSphase)-fTcalib.at(adr));
        if(abs((12.86*cfdtime+TCSphase)-fTcalib.at(adr))<20)
          hTimechprof->Fill(col,row,(12.86*cfdtime+TCSphase)-fTcalib.at(adr));
    }
    
    
    
    if(fExpertHistos )  {   
            fHch->Fill(adr);
    }
    SADCsum+=puls;
    if(puls>SADCampcut){
       if (hPhysTm) hPhysTm->Fill(time);
       if(row==23||row==24)  hprofx->Fill(col,puls);
       if(col==34||col==35)  hprofy->Fill(row,puls);
       SADCcuthits++; SADCcutsum+=puls;
//       fHchac->Fill(adr,puls);
       fHxyac->Fill(col,row,puls);
       hnhitcut->Fill(col,row);
   }
    if(fVrow->Test(row) &&
       fVcol->Test(col) &&
       fVamp->Test(puls)) {

      fVrow->Store(row);
      fVcol->Store(col);
      fVamp->Store(puls);
      fNhitsKept++;
    }
 }
  if(evtype==7&&!rndtr) {
     hSADCcuthits->Fill(SADCcuthits);
     hSADChit->Fill(fNSADChits);
     hSADCsum->Fill(SADCsum);
     hSADCcutsum->Fill(SADCcutsum);

  }
  if(evtype==8) {
    if (hLedProf) { makeProfiles(hLedProf, hmnled, hmnledrms, hmnledamp,hmnledr,NULL);
            makeReference(hLedProf,hmnled,href);
    }
    if (hPedProf) makeProfiles(hPedProf, hmnped, NULL, NULL, NULL, hmnpedr);
  }

  if (thr_flag) TThread::UnLock();
}


void PlaneECAL2::ControlPanel(const TGWindow *p, const TGWindow *main) {

  if (!fControlPanel) fControlPanel = new PlaneECAL2Panel(p, main, 100, 100, this);
}


void PlaneECAL2::AmpSpectrum() {

  std::string hdname = "EC02_1ch";
  if(fCurChan<fNchan) {
    //fHistList.pop_back();
    fHavsadr->ProjectionY(hdname.c_str(),fCurChan+1,fCurChan+1,"")->Draw();
  }
}

//==============================================================================

double PlaneECAL2::SADCped(int ipbg, int  ipend,int *data, double &RMS) {
   double mnped=0.;
   RMS=0.;
   for (int i=ipbg;i<=ipend;i++)  mnped+=data[i];
   mnped/=(ipend-ipbg+1);
   for (int i=ipbg;i<=ipend;i++)   { RMS+=(data[i]-mnped)*(data[i]-mnped);
}
   RMS=sqrt(RMS/(ipend-ipbg+1));
 return mnped;
}

//==============================================================================

double PlaneECAL2::SADCpuls(int ipbg, int  ipend,int *data,
                         double ped,double &sum,double &time) {
   double puls=0., pulspos=0.;
   sum=0.; time=0.;
   for (int i=ipbg;i<=ipend;i++) {
              double amp=data[i]-ped;
              if(amp>puls) {puls=amp; pulspos=i;}
              if(amp<0.) amp=0.;
              sum+=amp;  time+=amp*i;
   }
   if(sum>0.) time=time/sum;
   time=pulspos;
 return puls;
}

//==============================================================================

void PlaneECAL2::makeProfiles(TProfile *prof, TH1F *mean, TH1F *rms,
                         TH1F *integral,TH1F *rmsint,TH1F *rmsint1)
{
   int nbins=prof->GetNbinsX();
   if(mean!=NULL) mean->Reset();
   if(rms!=NULL) rms->Reset();
   if(integral!=NULL) integral->Reset();
   if(rmsint!=NULL) rmsint->Reset();
   if(rmsint1!=NULL) rmsint1->Reset();
   for (int i=1;i<nbins+1;i++) {
       int cnl=i-1;
       int ix=cnl/fNrows; int iy=cnl-ix*fNrows;
//       if(ix>31&&ix<42&&iy>19&&iy<29) continue;
       if(ix>33&&ix<36&&iy>22&&iy<25) continue;
       double Mean=prof->GetBinContent(i);
       double x=prof->GetBinCenter(i);
       if(mean!=NULL) mean->Fill(x,Mean);
       if(integral!=NULL) integral->Fill(Mean);
       double RMS1=prof->GetBinError(i);
       double RMS=0.;
       if(Mean>10.) RMS=RMS1/Mean*100.;
       if(rms!=NULL) rms->Fill(x,RMS);
       if(rmsint!=NULL) rmsint->Fill(RMS);
       if(rmsint1!=NULL) rmsint1->Fill(RMS1);
   }
}

//==============================================================================

void PlaneECAL2::makeReference(TProfile *prof, TH1F *mean, TProfile *ref)
{
//   float RefLedMin=50., RefLedMax=4000.;
   float RefLedMin=10., RefLedMax=3800.;
   if(ref!=NULL) ref->Reset(); else return;
   if(hrefXY!=NULL) hrefXY->Reset(); else return;
   if(mean==NULL) return;
   if(prof==NULL) return;
   int nbins=prof->GetNbinsX();
   const float* myref=((TH1F_Ref*)mean)->GetRef();
   if(myref==NULL) return;
   for (int i=1;i<nbins+1;i++) {
       double Mean=prof->GetBinContent(i);
       double x=prof->GetBinCenter(i);
       double Ref=myref[i-1];
       double Diff=200.;
       if(Ref>RefLedMin&&Ref<RefLedMax) Diff=(Mean/Ref)*100.;
       if(Diff<5.) Diff=5.;  
       ref->Fill(x,Diff);
       int cnl=i-1;
       int ix=cnl/fNrows; int iy=cnl-ix*fNrows;
       if((Diff<80||Diff>120)&& Diff!=200.)  hrefXY->Fill(ix,iy,Diff);
  }
}

//==============================================================================
//calculates the signal time using a CFD, output has to be corrected by
//calibration constants and TCSPhase
//==============================================================================
double PlaneECAL2::CFDtime(int ipbg, int ipend, int *data) {
  double data_delayed=0;
  double data_orig=0;
  double diff=0; 
  double diff_delayed=0;

  //currently some hardcoded parts, because the dynamic one is not tested yet
  //  float  basl1=data[0];
  //  float  basl2=data[1];
  for(  int i=0; i<ipend;i++) {
    if(i>2) {

      data_delayed=(data[i-2]-50)*2;
      data_orig=(data[i]-50);

      diff_delayed=diff;
      diff= data_orig-data_delayed;
      if(i>3) {
        if(diff_delayed -diff >10 &&diff<0 && diff_delayed>=0) {
          return i-5+(float)((diff_delayed)/(diff_delayed-diff)); 
        }

      }
    }
  }
  return -999;
};

double PlaneECAL2::Hmaxtime(int ipbg, int ipend, int *data) {
  int maxsample=0;
  double max=0.;
  double pedeven=data[0];
  double pedodd=data[1];
  for(  int i=0; i<ipend;i++) {
    if(i%2==0) {
      if(data[i]-pedeven>max){
        maxsample =i;
        max=data[i]-pedeven;
      }
    }
    else {
      if(data[i]-pedodd>max){
        maxsample =i;
        max=data[i]-pedodd;
      }
    }
  }
  if(max<20)
    return -99999;
  for(  int i=1; i<maxsample;i++) {
    if(i%2==0) {
      if(((data[i-1]-pedodd) <= (max/2.)) &&  ((data[i]-pedeven) >= (max/2.)) ) {
        return ( i-4 + (max/2. - (data[i-1]-pedodd))/(data[i]-pedeven-data[i-1]+pedodd) );
      }
    }
    else {
      if(((data[i-1]-pedeven) <= (max/2.)) &&  ((data[i]-pedodd) >= (max/2.)) ) {
        return ( i-4 + (max/2. - (data[i-1]-pedeven))/(data[i]-pedodd-data[i-1]+pedeven) );
      }

    }
  }
  return -99999;
};



//==============================================================================
void PlaneECAL2::TextOutput(ostream& out) {
    wwwled_perif/=(32*48);
    if(wwwled_perif<0.) wwwled_perif=0.; if(wwwled_perif>4000.) wwwled_perif=4000.;
    wwwled_cent/=(32*48);
    if(wwwled_cent<0.) wwwled_cent=0.; if(wwwled_cent>4000.) wwwled_cent=4000.;
    out<<"\tled_c "<<wwwled_cent<<std::endl;
    out<<"\tled_p "<<wwwled_perif<<std::endl;
    wwwled_cent=0.;
    wwwled_perif=0.;
}

#if USE_DATABASE == 1
void PlaneECAL2::ReadCalib(const tm &t)
{
  readledcalibCFD();
  // read-in corresponding calibration constants
  std::vector<TimeCalib> pre_tcalib_data;
  try{
    ReadFromDataBase(pre_tcalib_data,t,"_TimeCALIB");

//     printf("\nT0:  %f\n",fT0);
//     for(unsigned int i =0;i<pre_tcalib_data.size();i++) 
//       pre_tcalib_data.at(i).Print();

    if(pre_tcalib_data.size() < (unsigned) fNchan ){
      fUseCalib = false;
    }
    else    
      fUseCalib = true;

    fTcalib.resize(pre_tcalib_data.size());
//     if(fUseCalib =true)  // WTF ??? DN 14/10/2014
    if (fUseCalib == true) {
      for(unsigned int i=0;i<pre_tcalib_data.size();i++) {
        fTcalib.at( pre_tcalib_data.at(i).x*fNrows + pre_tcalib_data.at(i).y) = pre_tcalib_data.at(i).time + fT0;
//         cout<< fTcalib.at( pre_tcalib_data.at(i).x*fNrows + pre_tcalib_data.at(i).y)<<endl;
      }
    }
    else {
      for(int i=0;i<fNchan;i++) {
        fTcalib.at( pre_tcalib_data.at(i).x*fNrows + pre_tcalib_data.at(i).y) = 0;
      }
    }

  }
  

  catch(CS::Exception& e) {
    std::cerr<<__FILE__<<"  "<<__LINE__<<"  "<<e.what()<<std::endl;
  }
  catch(const std::exception &e) {
    std::cerr<<e.what()<<std::endl;
  }
  catch(...) {
    std::cout<<"PlaneECAL2::ReadCalib() ==> "<<GetName()
      <<" calibrations, valid for ";
    std::cout<<t.tm_mday<<"."<<t.tm_mon+1<<"."<<t.tm_year+1900<<" "
      <<t.tm_hour<<":"<<t.tm_min<<":"<<t.tm_sec
      <<", not found in DB"<<std::endl;
  }
}

template <class T>
void PlaneECAL2::ReadFromDataBase(T& v, const tm& t,const char* typecalib)
{
  if( fDataBase==NULL )
    throw CS::Exception("Plane::ReadFromDataBase():  data base is not opened.");
  try{
    struct tm tt(t);
    CDB::Time tp(mktime(&tt),0);
    std::string strdata("");
    fDataBase->read(fName,strdata,tp,typecalib);
    if (strdata == "") {
      std::cerr << "Plane::ReadFromDataBase() "<<GetName()<<" : no calibration file found"<<std::endl;      
      return; 
    }
    char tmp[256];
    char tmp2[256];
    std::istringstream istrdata(strdata.c_str());
    
    
    istrdata.getline(tmp,256);
    istrdata.getline(tmp,256);
    istrdata.getline(tmp,256);
    sscanf(tmp,"%s%f",tmp2,&fT0) ;
    
    istrdata >> v;
  }
  catch(CS::Exception& e) {
    std::cout<<"rethrowing"<<std::endl;
    throw;
  }
}

#endif //USE_DATABASE

void PlaneECAL2::readledcalibCFD() {
  try {
    if(ledcalibCFD==NULL)
      ledcalibCFD=new vector<Float_t>(3072,0);
    ledcalibCFD->clear();
    string filename="/online/detector/tum/ledcalibCFD";
    ifstream istrdata;
    istrdata.open(filename.c_str());
    if(!istrdata.is_open()) {
      cerr<<"Error opening file"<<filename<<endl;
      for(int i =0;i<3072;i++) 
        ledcalibCFD->push_back(0);
      return;
    }   
    char tmp[256];
    while(!istrdata.eof()) {

      istrdata.getline(tmp,256);
      ledcalibCFD->push_back(atoi(tmp));

    }
    istrdata.close();
  }
  catch (exception& e) {
    cout<<e.what()<<endl;
  }

}

