#include <cmath>

#include "PlaneFEM.h"
#include "PlaneFEMPanel.h"
#include "TProfile.h"
#include "ChipSADC.h"

#define nSADC_Sampl 32
#define Dpedstart  0    
#define Dpedend    2
#define Dpulsstart 3
#define Dpulsend   31
#define Dcalpedstart 0
#define Dcalpedend  2
#define Dledstart   3
#define Dledend     31

const int PlaneFEM::fMAX_MULT = 1;

ClassImp(PlaneFEM); \

PlaneFEM::PlaneFEM(const char *detname,int ncol, int nrow, int center, int width)
  : Plane(detname),fNrows(nrow),fNcols(ncol),fNchan(nrow*ncol), hLed(0), hLedChnl(0),
  hxynoise(0), hampnoise(0), hLedTm(0), hPhysTm(0),
  hmnled(0), hmnledrms(0), hmnledamp(0), hmnledr(0), hmnped(0), hmnpedr(0), hLedProf(0),
  hFemPuls(nrow), hFemAmpl(nrow), AmpLedF(nrow) {

  nSADCchan=nb_FEM_SADCcnl;
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

void PlaneFEM::Init(TTree* tree) {
  std::string myname;
  if(fName=="EC01FEM") myname+= "EC01F_";

  bool noupdate=true;

  std::string name,title;
  char cutamp[10],cutled[10];
  sprintf(cutamp,"%3.1i",SADCampcut);
  sprintf(cutled,"%3.1i",SADCLEDcut);

  name = myname + "Ampl_vs_chnl#";
  fHavsadr = new TH2F(name.c_str(), name.c_str(),
		fNrows*fNcols, 0, fNrows*fNcols, 110, 0.,1100.);
  AddHistogram(fHavsadr);

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
  

    for( int ch=1; ch<=nSADCchan; ch++ ) {
 	char b1[222],b2[222];
	
	sprintf(b1,"%sFEM%d_Pulse_vs_sample#",myname.c_str(),ch);
	sprintf(b2,"%s FEM %d Pulse vs sample#",fName.c_str(),ch);
	hFemPuls[ch-1]=new TProfile(b1,b2,fNsamples, 0., fNsamples,"s");
	AddHistogram(hFemPuls[ch-1]);
	}
	
    for( int ch=1; ch<=nSADCchan; ch++ ) {
 	char b1[222],b2[222];
	
	sprintf(b1,"%sFEM%d_Ampl",myname.c_str(),ch);
	sprintf(b2,"%s FEM %d Amplitude corr ped",fName.c_str(),ch);
	hFemAmpl[ch-1]=new TH1F_Ref(b1,b2,256, 0., 1024.,fRateCounter);
	((TH1F_Ref*)hFemAmpl[ch-1])->SetReference(fReferenceDirectory);
	AddHistogram(hFemAmpl[ch-1]);
	
	}

//
//
  if (fExpertHistos) {
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

void PlaneFEM::StoreDigit(int col, int row, std::vector<float>& data) {
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
// Modif for polarity change in SADC (630,7)
//    for(int i=0;i<len;i++) fSADCev[fNhits][i]=(int)data[i+2];
    for(int i=0;i<len;i++) fSADCev[fNhits][i]=1023-(int)data[i+2];
    fNhits++;
  }
}

void PlaneFEM::StoreDigit(CS::Chip::Digit* digit) {
  std::vector<float> data=digit->GetNtupleData();
  if(data.size()<3) return;
  const CS::ChipSADC::Digit* sadcdig = dynamic_cast<const CS::ChipSADC::Digit*>(digit);
  if(sadcdig != NULL ) {
         int ix=sadcdig->GetX(); int iy=sadcdig->GetY();
         if(ix < 0||iy<0||ix>fNcols||iy>fNrows) return; 
         this->StoreDigit(ix,iy,data);
  }
}

void PlaneFEM::EndEvent(const CS::DaqEvent &event) {
  int evtype=event.GetType();
  if(!(evtype==7||evtype==8))  return;
  int trmsk=event.GetTrigger();   
  const CS::uint32 Atr = event.GetHeader().GetTypeAttributes()[0]&0xf8000000;
  if(evtype==8&&Atr!=0x68000000) return;        // select calorimeter calib.ev.

  for( register int j=0;j<8;j++){
    AmpLedF[j]=0;
  }


  if (thr_flag) TThread::Lock();

  int rndtr=trmsk&0x800;     // random trigger 
  if(evtype==8) wwwled=0.;

  double sum=0., time=0.;
  double SADCsum=0., SADCcutsum=0.;
  int SADCcuthits=0; int amp;
  for ( register int i=0; i<fNhits; i++) {
    int col=fSADCCol[i];
    int row=fSADCRow[i];
    int adr=row + col*fNrows;
//
// Histos FEM ampl/sample
//

	if (col==0&&!rndtr&&evtype==8 && (row >= 0) && (row <= 7)) {
	  for( register int j=0;j<fNsamples;j++){hFemPuls[row]->Fill(j,fSADCev[i][j]);}
	}

//     
    double RMS=0.,ped=0.,puls=0.;
    if(evtype==8) {                // Leds
         ped=SADCped(calpedstart,calpedend,&fSADCev[i][0],RMS);
         puls=SADCpuls(ledstart,ledend,&fSADCev[i][0], ped,sum,time);
	 if (col == 0 && row >= 0 && row < 8) {
	   hFemAmpl[row]->Fill(puls);
	   AmpLedF[row]=puls;
	 }
         if(puls<0.) puls=0.;
         wwwled+=puls;
         hPed->Fill(ped);
         hPedRMS->Fill(RMS);
         hPedProf->Fill(adr,ped);
         hPedXY->Fill(col,row,ped);
         hLedProf->Fill(adr,puls);
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
    SADCsum+=puls;
    if(puls>SADCampcut){
       hprofx->Fill(col,puls);
       hprofy->Fill(row,puls);
       if (hPhysTm) hPhysTm->Fill(time);
       SADCcuthits++;
       SADCcutsum+=puls;
   }
 //end of loop on Hits
   
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

void PlaneFEM::ControlPanel(const TGWindow *p, const TGWindow *main) {

  if (!fControlPanel) fControlPanel = new PlaneFEMPanel(p, main, 100, 100, this);
}

void PlaneFEM::AmpSpectrum() {

  std::string hdname = fName+"_1ch";
  if(fCurChan<fNchan) {
    //fHistList.pop_back();
    fHavsadr->ProjectionY(hdname.c_str(),fCurChan+1,fCurChan+1,"")->Draw();
  }
}


//==============================================================================

double PlaneFEM::SADCped(int ipbg, int  ipend,int *data, double &RMS) {
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

double PlaneFEM::SADCpuls(int ipbg, int  ipend,int *data,
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

void PlaneFEM::makeProfiles(TProfile *prof, TH1F *mean, TH1F *rms,
                         TH1F *integral,TH1F *rmsint,TH1F *rmsint1)
{

   int nbins=prof->GetNbinsX();
   if(mean!=NULL) mean->Reset();
   if(rms!=NULL) rms->Reset();
   if(integral!=NULL) integral->Reset();
   if(rmsint!=NULL) rmsint->Reset();
   if(rmsint1!=NULL) rmsint1->Reset();
   for (int i=1;i<nbins+1;i++) {
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

void PlaneFEM::makeReference(TProfile *prof, TH1F *mean, TProfile *ref)
{
   if(!(fName=="EC01FEM")) return;
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
       if(Diff<2.) Diff=2.;  
       ref->Fill(x,Diff);
       int cnl=i-1;
       int ix=cnl/fNrows; int iy=cnl-ix*fNrows;
       if((Diff<80||Diff>120)&& Diff!=200.) hrefXY->Fill(ix,iy,Diff);
    }
}
PlaneFEM::~PlaneFEM() {
  delete [] fRow; delete [] fCol; delete [] fAmp;
}
//==============================================================================
void PlaneFEM::TextOutput(ostream& out) {
    std::string name;
    double ncn;
    if(fName=="EC01FEM") { ncn=8.; name="FEM ";}
    else return;
    wwwled/=ncn;
    if(wwwled<0.) wwwled=0.; if(wwwled>4000.) wwwled=4000.;
    out<<"\tled_c "<<wwwled<<std::endl;
    wwwled=0.;
}




