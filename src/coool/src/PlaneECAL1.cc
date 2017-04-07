#include <cmath>

#include "PlaneECAL1.h"
#include "PlaneFEM.h"
#include "TProfile.h"
#include "ChipSADC.h"
#include "PlaneECAL1Panel.h"

#define nSADC_Sampl 32
#define Dpedstart  0    
#define Dpedend    5
#define Dpulsstart 6
#define Dpulsend   31
#define Dcalpedstart 0
//#define Dcalpedend  11
//#define Dledstart   12
#define Dcalpedend  6
#define Dledstart   7
#define Dledend     31

const int SADCoverflow=1020;   // overflow for ADC channels
const int PlaneECAL1::fMAX_MULT = 1;

ClassImp(PlaneECAL1); \

PlaneECAL1::PlaneECAL1(const char *detname,int ncol, int nrow, int center, int width)
  : Plane(detname),fNrows(nrow),fNcols(ncol),fNchan(nrow*ncol), hLed(0), hLedChnl(0),
  hoverflow(0), hxynoise(0), hampnoise(0), hrnd(0), hledlen(0), hLedTm(0), hPhysTm(0),
  hmnled(0), hmnledrms(0), hmnledamp(0), hmnledr(0), hmnped(0), hmnpedr(0), hLedProf(0),
  hLasPuls(nrow*ncol), hLasAmpl(nrow*ncol) {

  nSADCchan=nb_ECAL1_SADCcnl;
  fNsamples=nSADC_Sampl;
  fNhits=0;

  SADCampcut = 10;
  SADCLEDcut = 20;
  RefLedMin=50;
  RefLedMax=1000;
  href=NULL, hprofx=NULL, hprofy=NULL, hrefXY=NULL;
  wwwled=0.;
//
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
  std::string ampcutname=fName+"_amp/led_cuts";
  fVampcut=AddVariable(ampcutname.c_str(),100,SADCampcut,SADCLEDcut,5);
//
  nCounterX=0;
  nCounterY=16;
//  name=fName+"_Led_counter";
//  fNcounter=AddVariable(name.c_str(),100,nCounterX,nCounterY,5);
//
  fRow = new int[fNchan*fMAX_MULT];
  fCol = new int[fNchan*fMAX_MULT];
  fAmp = new int[fNchan*fMAX_MULT];
  fAmpsum = new int[fNchan*fMAX_MULT];

  fPulseAmp = new double[fNchan];

  std::string rowname=fName+"_row";
  std::string colname=fName+"_col";
  std::string ampname=fName+"_amp";
  std::string sumname=fName+"_sum";

  fVrow=AddVariable(rowname.c_str(),fNrows,0,fNrows ,fNchan*fMAX_MULT);
  fVcol=AddVariable(colname.c_str(),fNcols,0,fNcols ,fNchan*fMAX_MULT);
  fVamp=AddVariable(ampname.c_str(),100,center-width,center+width,
		    fNchan*fMAX_MULT);
//  fVamp=AddVariable(sumname.c_str(),100,center-width,center+width,
//		    fNchan*fMAX_MULT);
}

void PlaneECAL1::Init(TTree* tree) {
  std::string myname;
  if(fName=="EC01P00") myname+= "EC01G_";
  else if(fName=="EC01P01") myname+= "EC01M_";
  else if(fName=="EC01P02") myname+= "EC01O_";
  
  bool noupdate=true;

  std::string name,title;
  char cutamp[10],cutled[10];
  sprintf(cutamp,"%3.1i",SADCampcut);
  sprintf(cutled,"%3.1i",SADCLEDcut);

  name = myname + "SADC_hits#";
  hSADChit=new TH1F_Ref(name.c_str(),name.c_str(),100,0,100, fRateCounter);
  ((TH1F_Ref*)hSADChit)->SetReference(fReferenceDirectory);
  AddHistogram(hSADChit);

  name = myname + "SADC_hits#(amp>thr)";
  title = myname + "SADC_hits#,amp>"+cutamp;
  hSADCcuthits=new TH1F_Ref(name.c_str(),title.c_str(), 30, 0, 30, fRateCounter);   
  ((TH1F_Ref*)hSADCcuthits)->SetReference(fReferenceDirectory);
  AddHistogram(hSADCcuthits);

  name = myname+"Nhits_vs_XY";
  title = fName+" Number of hits vs XY";
  fHrc1=new TH2F(name.c_str(),title.c_str(),fNcols,0.,fNcols,fNrows,0.,fNrows);
  fHrc1->SetOption("box");
  AddHistogram(fHrc1);

  name = myname+"Nhits_vs_XY(amp>thr)";
  title = fName+" Number of hits vs XY, amp>"+cutamp;
  hnhitcut=new TH2F(name.c_str(),title.c_str(),fNcols,0.,fNcols,fNrows,0.,fNrows);
  hnhitcut->SetOption("box");
  AddHistogram(hnhitcut);

  name = myname + "Nhits_vs_chnl#";
  fHch=new TH1F_Ref(name.c_str(),name.c_str(), fNrows*fNcols, 0,
                     fNrows*fNcols, fRateCounter);
  ((TH1F_Ref*)fHch)->SetReference(fReferenceDirectory);
  AddHistogram(fHch);

  name = myname + "Ampl_vs_XY";
  title = fName + "Amplitude vs XY";
  fHrca=new TH2F(name.c_str(),title.c_str(),fNcols,0.,fNcols,fNrows,0.,fNrows);
  fHrca->SetOption("col");
  AddHistogram(fHrca);

  name = myname + "Ampl_vs_XY(amp>thr)";
  title = myname + "Ampl. vs XY,amp>"+cutamp;
  fHxyac = new TH2F(name.c_str(), title.c_str(),
                      fNcols, 0, fNcols, fNrows, 0, fNrows);
  fHxyac->SetOption("col");
  AddHistogram(fHxyac);

  name = myname + "Ampl_vs_chnl#";
  fHavsadr = new TH2F(name.c_str(), name.c_str(),
		fNrows*fNcols, 0, fNrows*fNcols, 110, 0.,1100.);
  AddHistogram(fHavsadr);

  name = myname + "SADC_amp";
  title = fName + "SADC max. amplitude for each chnl.";
  fHa=new TH1F_Ref(name.c_str(),title.c_str(), 225,0., 1100, fRateCounter);
  ((TH1F_Ref*)fHa)->SetReference(fReferenceDirectory);
  AddHistogram(fHa);

  name = myname + "SumAmpl";
  title = fName + " SADC Sum Amplitude";
  hSADCsum=new TH1F_Ref(name.c_str(),title.c_str(),200,0.,2000., fRateCounter);
  ((TH1F_Ref*)hSADCsum)->SetReference(fReferenceDirectory);
  AddHistogram(hSADCsum);

  name = myname + "SumAmp(amp>thr)";
  title = fName + "SADC Sum Amplitude with amp>"+cutamp;
  hSADCcutsum=new TH1F_Ref(name.c_str(),title.c_str(), 200,0.,2000., fRateCounter);
  ((TH1F_Ref*)hSADCcutsum)->SetReference(fReferenceDirectory);
  AddHistogram(hSADCcutsum);
  
  name = myname + "Peds";
  title = fName + "SADC pedestals";
  hPed=new TH1F_Ref(name.c_str(),title.c_str(), 110, 0., 1100., fRateCounter);
  ((TH1F_Ref*)hPed)->SetReference(fReferenceDirectory);
  AddHistogram(hPed);

  name = myname + "PedRMS";
  title = fName + " SADC pedestal RMS (counts)";
  hPedRMS=new TH1F_Ref(name.c_str(),title.c_str(), 100, 0., 20., fRateCounter);
  ((TH1F_Ref*)hPedRMS)->SetReference(fReferenceDirectory);
  AddHistogram(hPedRMS);

  name = myname + "Peds_vs_XY";
  title = fName + "Peds vs XY";
  hPedXY=new TH2F(name.c_str(),title.c_str(),fNcols,0.,fNcols,fNrows,0.,fNrows);
  hPedXY->SetOption("col");
  AddHistogram(hPedXY);

  name = myname + "Peds/RMS_vs_chnl#";
  title = fName + " SADC pedestal/RMS vs chnl#";
  hPedProf=new TProfile(name.c_str(),title.c_str(),fNchan, 0., fNchan,"s");
  AddHistogram(hPedProf);
  
  name = myname + "Leds_vs_XY";
  title = fName + " Leds vs XY, led>"+cutled;
  hLed=new TH2F(name.c_str(),title.c_str(),fNcols,0.,fNcols,fNrows,0.,fNrows);
  hLed->SetOption("col");
  AddHistogram(hLed);

  name = myname + "Led_ev_length";
  title = fName + "Led event length";
  hledlen=new TH1F_Ref(name.c_str(),title.c_str(),305,0,610, fRateCounter);
  ((TH1F_Ref*)hledlen)->SetReference(fReferenceDirectory);
  AddHistogram(hledlen);

  name = myname + "LedReference";
  title = fName + " Led deviation from reference(%)";
  href= new TProfile(name.c_str(), title.c_str(), fNchan,0.,fNchan,"s");
  AddHistogram(href);

  name = myname+"BadLedRef_vs_XY";
  title = fName+"Bad Led Ref. vs XY";
  hrefXY=new TH2F(name.c_str(),title.c_str(),fNcols,0.,fNcols,fNrows,0.,fNrows);
  hrefXY->SetOption("col");
  AddHistogram(hrefXY);

  name = myname + "Led/RMS_vs_chnl#";
  title = fName + "LEDS/RMS vs channel#";
  hLedProf=new TProfile(name.c_str(),title.c_str(),fNrows*fNcols, 0., fNrows*fNcols,"s");
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
                              fNrows*fNcols, 0., fNrows*fNcols, fRateCounter, noupdate);
  ((TH1F_Ref*)hmnled)->SetReference(fReferenceDirectory);
  AddHistogram(hmnled);
//
//
  if (fExpertHistos) {
  
    for ( int icol=0; icol<fNcols; icol++) {
    	for( int irow=0; irow<fNrows; irow++ ) {
	int ch=irow + icol*fNrows;
 	char b1[222],b2[222];
	
	sprintf(b1,"%s_%d_%d_Las_Pulse_vs_sample#",myname.c_str(),icol,irow);
	sprintf(b2,"%s %d %d Las Pulse vs sample#",fName.c_str(),icol,irow);
	hLasPuls[ch]=new TProfile(b1,b2,fNsamples, 0., fNsamples,"s");
	AddHistogram(hLasPuls[ch]);
	
	}
    }
  
    for ( int icol=0; icol<fNcols; icol++) {
    	for( int irow=0; irow<fNrows; irow++ ) {
	int ch=irow + icol*fNrows;
 	char b1[222],b2[222];
	
	sprintf(b1,"%s_%d_%d_Las_Ampl",myname.c_str(),icol,irow);
	sprintf(b2,"%s %d %d Las Amplitude corr ped",fName.c_str(),icol,irow);
	hLasAmpl[ch]=new TH1F_Ref(b1,b2,256, 0., 1024., fRateCounter);
	((TH1F_Ref*)hLasAmpl[ch])->SetReference(fReferenceDirectory);
	AddHistogram(hLasAmpl[ch]);
	
	}
    }
  
    name = myname + "MeanPedRMS";
    title = fName + " SADC mean  pedestal RMS";
    hmnpedr=new TH1F(name.c_str(),title.c_str(), 100, 0., 20.);
    AddHistogram(hmnpedr);

    name = myname + "Leds";
    title = fName + "Mean Led amplitudes";
    hmnledamp=new TH1F(name.c_str(),title.c_str(),110,0.,1100.);
    AddHistogram(hmnledamp);

    name = myname + "LedsRMS";
    title = fName + " Led RMSs(%),led>"+cutled;
    hmnledr=new TH1F(name.c_str(),title.c_str(),400,0.,100.);
    AddHistogram(hmnledr);

    name = myname + "Led_vs_chnl#";
    title = fName + "LED vs channel#";
    hLedChnl=new TH2F(name.c_str(),title.c_str(), 
                              fNrows*fNcols, 0., fNrows*fNcols,204,0.,1020.);
    AddHistogram(hLedChnl);

    name = myname + "overflow";
    title = fName + "ADC overfloved channels";
    hoverflow=new TH2F(name.c_str(),title.c_str(),fNcols,0.,fNcols,fNrows,0.,fNrows);
    hoverflow->SetOption("col");
    AddHistogram(hoverflow);

    name = myname + "Random_ev_length";
    title = fName + "Random trigger event length";
    hrnd=new TH1F(name.c_str(),title.c_str(),50,0,50);
    AddHistogram(hrnd);

    name = myname + "Nhits_vs_XY_rnd_trig";
    title = fName + " hits vs XY, Random trigger";
    hxynoise = new TH2F(name.c_str(), title.c_str(),fNcols, 0, fNcols, fNrows, 0, fNrows);
    hxynoise->SetOption("box");
    AddHistogram(hxynoise);

    name = myname + "SADC_Led_timing";
    title = fName + "SADC Led max.ampl vs sample#,amp>"+cutled;
    hLedTm=new TH1F(name.c_str(),title.c_str(),32,0.,32.);
    AddHistogram(hLedTm);

    name = myname + "Phys_ev_timing";
    title = fName + "SADC Phys event max.ampl vs sample#,amp>"+cutamp;
    hPhysTm=new TH1F(name.c_str(),title.c_str(),32,0.,32.);
    AddHistogram(hPhysTm);

    name = myname + "MeanPed_vs_chnl#";
    title = fName + "Mean SADC pedestal vs channel#";
    hmnped=new TH1F(name.c_str(),title.c_str(), 
                              fNrows*fNcols, 0., fNrows*fNcols);
    AddHistogram(hmnped);

    name = myname + "MnLedRMS_vs_chnl#";
    title = fName + "Mean LED RMS(%) vs channel#,Led>"+cutled;
    hmnledrms=new TH1F(name.c_str(),title.c_str(), 
                              fNrows*fNcols, 0., fNrows*fNcols);
    AddHistogram(hmnledrms);

    char cntrnb[10];
    sprintf(cntrnb,"%2.2i-%2.2i",nCounterX,nCounterY);
    name = myname + "LED samples";
    title = fName + "LED event,cntr="+cntrnb;
    hampnoise=new TH2F(name.c_str(),title.c_str(),nSADC_Sampl,0.,nSADC_Sampl,1024,0.,1024);
    hampnoise->SetOption("col");
    AddHistogram(hampnoise);
  }

  std::string hitsname = fName + "_hits";
  std::string hitsleavlist = hitsname + "/I";
  if(tree) {
    fIsInTree = true;
    std::string leavlist0 =fVrow->GetName()  + "[" + hitsname + "]/F";
    std::string leavlist1 =fVcol->GetName()  + "[" + hitsname + "]/F";
    std::string leavlist2 =fVamp->GetName()  + "[" + hitsname + "]/F";
    tree->Branch(hitsname.c_str(),&fNhitsKept,
		 hitsleavlist.c_str(),32000);
    tree->Branch(fVrow->GetName().c_str(),fVrow->GetValues(),
		 leavlist0.c_str(),32000);
    tree->Branch(fVcol->GetName().c_str(),fVcol->GetValues(),
		 leavlist1.c_str(),32000);
    tree->Branch(fVamp->GetName().c_str(),fVamp->GetValues(),
		 leavlist2.c_str(),32000);

  }
}

void PlaneECAL1::StoreDigit(int col, int row, std::vector<float>& data) {
//
  int sadclen=data.size()-2;
  if(sadclen<1) return;
//
  if (fNhits < nSADCchan) {
    fSADCRow[fNhits]=row;
    fSADCCol[fNhits]=col;
    int len=fNsamples;
    if(sadclen<len){ len=sadclen;
                    for(int i=len;i<fNsamples;i++) fSADCev[fNhits][i]=0; }
    for(int i=0;i<len;i++) fSADCev[fNhits][i]=(int)data[i+2];
    fNhits++;
  }
}

void PlaneECAL1::StoreDigit(CS::Chip::Digit* digit) {
  std::vector<float> data=digit->GetNtupleData();
  if(data.size()<3) return;
  const CS::ChipSADC::Digit* sadcdig = dynamic_cast<const CS::ChipSADC::Digit*>(digit);
  if(sadcdig != NULL ) {
         int ix=sadcdig->GetX(); int iy=sadcdig->GetY();
         if(ix < 0||iy<0||ix>fNcols||iy>fNrows) return; 
         this->StoreDigit(ix,iy,data);
  }
}

void PlaneECAL1::EndEvent(const CS::DaqEvent &event) {
  int evtype=event.GetType();
  if(!(evtype==7||evtype==8))  return;
  int trmsk=event.GetTrigger();   
  const CS::uint32 Atr = event.GetHeader().GetTypeAttributes()[0]&0xf8000000;
  if(evtype==8&&Atr!=0x68000000) return;        // select calorimeter calib.ev.

  if (thr_flag) TThread::Lock();

  int rndtr=trmsk&0x800;     // random trigger 
  if(hrnd && evtype==7&&rndtr)  hrnd->Fill(fNhits);
  if(hledlen && evtype==8) hledlen->Fill(fNhits);
  if(evtype==8) wwwled=0.;

  double sum=0., time=0.;
  double SADCsum=0., SADCcutsum=0.;
  int SADCcuthits=0; int amp;
  for ( int i=0; i<fNhits; i++) {
    int col=fSADCCol[i];
    int row=fSADCRow[i];
    int adr=row + col*fNrows;
//
    if(evtype==7) {                // Overflow
      for(int j=0;j<fNsamples;j++){
       if(hoverflow && fSADCev[i][j]>SADCoverflow) hoverflow->Fill(col,row);
      }
     }
//     
    double RMS=0.,ped=0.,puls=0.;
    if(evtype==8) {                // Leds
           if (fExpertHistos) for( register int j=0;j<fNsamples;j++){hLasPuls[adr]->Fill(j,fSADCev[i][j]);}
         ped=SADCped(calpedstart,calpedend,&fSADCev[i][0],RMS);
         puls=SADCpuls(ledstart,ledend,&fSADCev[i][0], ped,sum,time);
	 if (adr < fNchan) {
	  fPulseAmp[adr] = puls;
	  if (fExpertHistos) hLasAmpl[adr]->Fill(puls);
	 }
         if(puls<0.) puls=0.;
         wwwled+=puls;
         hPed->Fill(ped);
         hPedRMS->Fill(RMS);
         hPedProf->Fill(adr,ped);
         hPedXY->Fill(col,row,ped);
         hLedProf->Fill(adr,puls);
//	 std::cout<<"Ampl FEM 3: "<<PlaneFEM::GetFEMAmp()<<std::endl;
         if(puls>SADCLEDcut) hLed->Fill(col,row,puls);
         if (fExpertHistos) {
           hLedChnl->Fill(adr,puls);
           if(col==nCounterX&&row==nCounterY) {
             for(int j=0;j<fNsamples;j++){
                 hampnoise->Fill(j,fSADCev[i][j]);
	   
             }
           }
           if(puls>SADCLEDcut) {
                 hLedTm->Fill(time);
           }
         }
         amp=(int)puls;
         if(fVrow->Test(row) &&
            fVcol->Test(col) &&
            fVamp->Test(amp)) {
              fVrow->Store(row);
              fVcol->Store(col);
              fVamp->Store(amp);
//            fVsum->Store(sum);
              fNhitsKept++;
         }
         continue;
    }


    ped=SADCped(pedstart, pedend, &fSADCev[i][0],RMS);
    puls=SADCpuls(pulsstart, pulsend,&fSADCev[i][0], ped,sum,time);
    if(puls<0.) puls=0.;
//
    amp=(int)puls;
    if(fVrow->Test(row) &&
       fVcol->Test(col) &&
       fVamp->Test(amp)) {
           fVrow->Store(row);
           fVcol->Store(col);
           fVamp->Store(amp);
//           fVsum->Store(sum);
       fNhitsKept++;
    }
//
    if(evtype==7&&rndtr)  {
      if(fExpertHistos)  {
         hxynoise->Fill(col,row);
//         hampnoise->Fill(col,row,puls);
       }
         continue;
    }                                       //random trigger
    fHavsadr->Fill(adr,puls);
    fHrca->Fill(col,row,puls);
    fHrc1->Fill(col,row);
    fHch->Fill(adr);
    fHa->Fill(puls);
    SADCsum+=puls;
    if(puls>SADCampcut){
       hprofx->Fill(col,puls);
       hprofy->Fill(row,puls);
       if (hPhysTm) hPhysTm->Fill(time);
       SADCcuthits++;
       SADCcutsum+=puls;
//       fHchac->Fill(adr,puls);
       fHxyac->Fill(col,row,puls);
       hnhitcut->Fill(col,row);
   }
 }
  if(evtype==7&&!rndtr){
     hSADCcuthits->Fill(SADCcuthits);
     hSADChit->Fill(fNhits);
     hSADCsum->Fill(SADCsum);
     hSADCcutsum->Fill(SADCcutsum);
  }
  if(evtype==8) { 
     if (hLedProf) {
        makeProfiles(hLedProf, hmnled, hmnledrms, hmnledamp,hmnledr,NULL);
        makeReference(hLedProf,hmnled,href);
     }
     if (hPedProf) makeProfiles(hPedProf, hmnped, NULL, NULL, NULL, hmnpedr);
  }

  if (thr_flag) TThread::UnLock();
}


void PlaneECAL1::ControlPanel(const TGWindow *p, const TGWindow *main) {

  if (!fControlPanel) fControlPanel = new PlaneECAL1Panel(p, main, 100, 100, this);
}


void PlaneECAL1::AmpSpectrum() {

  std::string hdname = fName+"_1ch";
  if(fCurChan<fNchan) {
    //fHistList.pop_back();
    fHavsadr->ProjectionY(hdname.c_str(),fCurChan+1,fCurChan+1,"")->Draw();
  }
}

//==============================================================================

double PlaneECAL1::SADCped(int ipbg, int  ipend,int *data, double &RMS) {
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

double PlaneECAL1::SADCpuls(int ipbg, int  ipend,int *data,
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

void PlaneECAL1::makeProfiles(TProfile *prof, TH1F *mean, TH1F *rms,
                         TH1F *integral,TH1F *rmsint,TH1F *rmsint1)
{

   int nbins=prof->GetNbinsX();
   if(mean!=NULL) mean->Reset();
   if(rms!=NULL) rms->Reset();
   if(integral!=NULL) integral->Reset();
   if(rmsint!=NULL) rmsint->Reset();
   if(rmsint1!=NULL) rmsint1->Reset();
   for (int i=1;i<nbins+1;i++) {
    if(fName=="EC01P00") {
       int cnl=i-1;
       int ix=cnl/fNrows; int iy=cnl-ix*fNrows;
//       if(ix>7&&ix<36&&iy>3&&iy<20) continue;
       if(ix>11&&ix<32&&iy>5&&iy<18) continue;                    // new 2012 GAMS central hole
    }
       double Mean=prof->GetBinContent(i);
       double x=prof->GetBinCenter(i);
       if(mean!=NULL) mean->Fill(x,Mean);
       if(integral!=NULL) integral->Fill(Mean);
       double RMS1=prof->GetBinError(i);
       double RMS=0.;
       if(Mean>SADCLEDcut) RMS=RMS1/Mean*100.; else RMS=0.;
       if(rms!=NULL) rms->Fill(x,RMS);
       if(rmsint!=NULL) rmsint->Fill(RMS);
       if(rmsint1!=NULL) rmsint1->Fill(RMS1);
  }
}

//==============================================================================

void PlaneECAL1::makeReference(TProfile *prof, TH1F *mean, TProfile *ref)
{
   if(!(fName=="EC01P00"||fName=="EC01P01"||fName=="EC01P02")) return;
   if(ref!=NULL) ref->Reset(); else return;
   if(hrefXY!=NULL) hrefXY->Reset(); else return;
   if(mean==NULL) return;
   if(prof==NULL) return;
   int nbins=prof->GetNbinsX();
   const float* myref=((TH1F_Ref*)mean)->GetRef();
   if(myref==NULL) return;
   for (int i=1;i<nbins+1;i++) {
     if(fName=="EC01P00") {
       int cnl=i-1;
       int ix=cnl/fNrows; int iy=cnl-ix*fNrows;
//       if(ix>7&&ix<36&&iy>3&&iy<20) continue;
       if(ix>11&&ix<32&&iy>5&&iy<18) continue;                    //  new 2012 central GAMS hole
     }
       double Mean=prof->GetBinContent(i);
       double x=prof->GetBinCenter(i);
       double Ref=myref[i-1];
       double Diff=200.;
       if(Ref>RefLedMin&&Ref<RefLedMax) Diff=(Mean/Ref)*100.;
       if(Diff<2.) Diff=2.;  
       ref->Fill(x,Diff);
       int cnl=i-1;
       int ix=cnl/fNrows; int iy=cnl-ix*fNrows;
       if((Diff<80||Diff>120)&& Diff!=200.) hrefXY->Fill(ix,iy,Diff);
    }
}
PlaneECAL1::~PlaneECAL1() {
  delete [] fRow; delete [] fCol; delete [] fAmp; delete [] fPulseAmp;
}
//==============================================================================
void PlaneECAL1::TextOutput(ostream& out) {
    std::string name;
    double ncn;
//    if(fName=="EC01P00") { ncn=608.; name="ecal1_G ";}
    if(fName=="EC01P00") { ncn=816.; name="ecal1_G ";}                               //  new 2012 central GAMS hole
    else if(fName=="EC01P01") {ncn=572.; name="ecal1_M ";}
    else if(fName=="EC01P02") {ncn=320.; name="ecal1_O ";}
    else return;
    wwwled/=ncn;
    if(wwwled<0.) wwwled=0.; if(wwwled>4000.) wwwled=4000.;
    out<<"\tled_c "<<wwwled<<std::endl;
    wwwled=0.;
}




