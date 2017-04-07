//official version 26 August 2009
//official version 20 September 2009
// 3 points around i_max 26.oct.09 
// parabola around i_max 28.oct.09
// 
// ECAL0      29 Sept 2012
// ECAL0SUM   01 Oct  2012
// ECAL0FEM   01 Oct  2012
// created by Gavrishchuk O. 10 October 2012
//
#include <cmath>
  
#include "PlaneECAL0.h"
#include "ChipSADC.h"
#include "GroupDAQ.h"

#include "DaqEvent.h"

#include <string.h>
#include <cstring>
#include <cstdlib>
#include <sstream>

//#include "TColor.h"
#include <cstdio>
#include "Plane2V.h"
#include "PlaneECAL0.h"
#include <strstream>
#include <iomanip>
#include <iostream>
 
#include <TProfile.h>
#define  fNsamples 32
#define  pedlen 5
 
using namespace std;

ClassImp(PlaneECAL0);
const int PlaneECAL0::fAmpCut = 1;
//--------------------------------------------------------------------------------------------  
void PlaneECAL0::Init(TTree* tree) {
//--------------------------------------------------------------------------------------------  
  Plane2V::Init(tree);
                
      if(strncmp(GetName(),"EC00P1__",8)== 0) {fNrows = 27; fNcols=33;}
      if(strncmp(GetName(),"EC00FEM_",8)== 0) {fNrows = 2; fNcols=2;}
      if(strncmp(GetName(),"EC00SUM_",8)== 0) {fNrows = 2; fNcols=2;}

  char cell_name[132];

  std::string amprowcol;
  std::string name, title;

//  cout << "All fNameNew : " << fName.c_str()<< endl;
  
  std::string ampname = fName + "_Sumampcut";
  fHsac=new TH1F(ampname.c_str(),ampname.c_str(),
		 fVamp->GetNbins(), fVamp->GetMin(),
		 fVamp->GetMax());
  fHsac->SetFillColor(8);
  AddHistogram(fHsac);

  std::string adrname = fName + "_adr";
  fHch=new TH1F_Ref( adrname.c_str(),adrname.c_str(), fNrows*fNcols, 0,
                     fNrows*fNcols, fRateCounter);
  ((TH1F_Ref*)fHch)->SetReference(fReferenceDirectory);
  
  fHch->SetFillColor(3);
  AddHistogram(fHch);
   
  std::string aadrcut = fName + "_adrVSampcut";
  fHchac = new TH2F(aadrcut.c_str(), aadrcut.c_str(), 
		    fNrows*fNcols, 0, fNrows*fNcols,
		    fVamp->GetNbins(), fVamp->GetMin(), 
		    fVamp->GetMax());
               fHchac->SetOption("boxcol");
  AddHistogram(fHchac);
  
  // 7    
  std::string ampmaxname = fName + "_maxamp";
  fHma=new TH1F_Ref(ampmaxname.c_str(),ampmaxname.c_str(),
                       fVamp->GetNbins(), fVamp->GetMin(),
                       fVamp->GetMax(), fRateCounter);
  ((TH1F_Ref*)fHma)->SetReference(fReferenceDirectory);
  fHma->SetFillColor(5);
  AddHistogram(fHma);

  // 8    
  std::string hitsthname = fName + "_hitsth";
  fHhth=new TH1F_Ref(hitsthname.c_str(),hitsthname.c_str(), 50, 0, 50, fRateCounter); 
  //   fNrows*fNcols, 0, fNrows*fNcols);
  ((TH1F_Ref*)fHhth)->SetReference(fReferenceDirectory);
  fHhth->SetFillColor(7);
  AddHistogram(fHhth);
  
  // 9
  std::string xyacut = fName + "_xVSyAcut";
  fHxyac = new TH2F(xyacut.c_str(), xyacut.c_str(),
                      fNcols, 0, fNcols, fNrows, 0, fNrows);
               fHxyac->SetOption("boxcolz");
               fHxyac->SetFillColor(4);
  AddHistogram(fHxyac);
  
// 10 
  std::string xynoise = fName + "_xVSynoise";
  fHxyn = new TH2F(xynoise.c_str(), xynoise.c_str(),fNcols, 0, fNcols, fNrows, 0, fNrows);
  fHxyn->SetOption("boxcolz");
  AddHistogram(fHxyn);

  
   
name = fName + "_LED_sum";
  fSum_led = new TH1F(name.c_str(),name.c_str(), 512, 0, 4096);
  fSum_led->SetFillColor(30);
  AddHistogram(fSum_led);
//  
name = fName + "_LED_xy";
    fAmpLed = new TH2F(name.c_str(),name.c_str(),fNcols,0,fNcols,fNrows,0,fNrows);
    fAmpLed->SetOption("colz");
//  fAmpLed->SetFillColor(2);
//  fAmpLed->SetMinimum(0.0);
  AddHistogram(fAmpLed);
//  
//name = fName + "_LED_profile";
//   fProf = new TProfile(name.c_str(),name.c_str(),fNrows*fNcols,0,fNrows*fNcols);
//          fProf->SetMaximum(200.0);
//            fProf->SetMinimum(0.0);
//          fProf->SetDrawOption("E1");
//          fProf->SetOption("E1BARHIST");
//            fProf->SetDrawOption("E1BAR");
//          fProf->SetFillColor(30);
//            fProf->SetMarkerColor(2);
//            fProf->SetMarkerSize(0.4);
//            fProf->SetMarkerStyle(20);
// AddHistogram(fProf);
//
name = fName + "_LED_adr";
   fHchampled = new TH2F(name.c_str(),name.c_str(),fNrows*fNcols,0,fNrows*fNcols,512,0,4096);
//  Int_t ci;   // for color index setting
//  ci = TColor::GetColor("#c82dd2");
//   fHchampled->SetMarkerColor(ci);
//   fHchampled->SetMarkerColor(1);
//   fHchampled->SetMarkerStyle(20);
//   fHchampled->SetMarkerSize(0.4);
   fHchampled->SetOption("boxcol");
   fHchampled->SetMinimum(0.0);
   AddHistogram(fHchampled);
//
name = fName + "adr_RMS_Ped";
   fRMS_Ped_ch = new TH2F(name.c_str(),name.c_str(),fNrows*fNcols,0,fNrows*fNcols,20,0,20);
   fRMS_Ped_ch->SetOption("boxcol");
  AddHistogram(fRMS_Ped_ch);

name = fName + "adr_Ped";
fPed_ch = new TH2F(name.c_str(),name.c_str(),fNrows*fNcols,0,fNrows*fNcols,100,0,100);
   fPed_ch->SetOption("boxcol");
   fPed_ch->SetMinimum(0.0);
   AddHistogram(fPed_ch);

name = fName + "_xy_Ped";
fPed_xy = new TH2F(name.c_str(),name.c_str(),fNcols,0,fNcols,fNrows,0,fNrows);
fPed_xy->SetOption("boxcolz");
     AddHistogram(fPed_xy);

name = fName + "_xy_RMS_Ped";
fRMS_Ped_xy = new TH2F(name.c_str(),name.c_str(),fNcols,0,fNcols,fNrows,0,fNrows);
     fRMS_Ped_xy->SetOption("boxcolz");
     AddHistogram(fRMS_Ped_xy);

name = fName + "xy_Time";
fTime_xy = new TH2F(name.c_str(),name.c_str(),fNcols,0,fNcols,fNrows,0,fNrows);
     fTime_xy->SetOption("boxcolz");
     AddHistogram(fTime_xy);

name = fName + "sum_time";
   fSum_time = new TH1F(name.c_str(),name.c_str(),320,0,32);
  AddHistogram(fSum_time);

name = fName + "_sum_Ped";
fPed_sum = new TH1F(name.c_str(),name.c_str(),100,0,100);
   AddHistogram(fPed_sum);

name = fName + "_sum_RMS_Ped";
fRMS_Ped_sum = new TH1F(name.c_str(),name.c_str(),200,0,20);
   AddHistogram(fRMS_Ped_sum);

//  
//------------------------------------------------------------------------------
//

if (fExpertHistos) {
name = fName + "_Bad";
   fBad = new TH1F(name.c_str(),name.c_str(),1024,0,4096);
  AddHistogram(fBad);
   
name = fName + "_Good";
   fGood = new TH1F(name.c_str(),name.c_str(),1024,0,4096);
  AddHistogram(fGood);
//  
//------------------------------------------------------------------------------
//
  for( int x = 0; x < fNcols; x++) {
  for( int y = 0; y < fNrows; y++) {
if(strncmp(GetName(),"EC00P1__",8)== 0){

        if (y >  5  && y < 21 && x >   5 && x < 27) continue; // ECAL0 window
        
        if (y >= 0  && y <  6 && x >=  0 && x <  3) continue; // ECAL0 low left
        if (y >= 0  && y <  6 && x >  29 && x < 33) continue; // ECAL0 low right
        if (y > 20  && y < 27 && x >=  0 && x <  3) continue; // ECAL0  up left 
        if (y > 20  && y < 27 && x >  29 && x < 33) continue; // ECAL0  up right
        
}
//----------------------------------------------------------- 
sprintf(cell_name,"Time_X_%d_Y_%d",x,y);
name  = fName + cell_name;
title = fName + cell_name;
  fCH_time[x][y] = new TH1F(name.c_str(), title.c_str(),32,0,32);
  AddHistogram(fCH_time[x][y]);
//-----------------------------------------------------------
sprintf(cell_name,"Ampl_X_%d_Y_%d",x,y);
name  = fName + cell_name;
title = fName + cell_name;
  fCH_amp[x][y] = new TH1F(name.c_str(), title.c_str(),1024,0,4096);
  AddHistogram(fCH_amp[x][y]);
//-----------------------------------------------------------
sprintf(cell_name,"Samp_X_%d_Y_%d",x,y);
name  = fName + cell_name;
title = fName + cell_name;
    fSProfile[x][y] = new TH2F(name.c_str(),title.c_str(),32,0,32,256,0,4096);
    fSProfile[x][y]->SetMinimum(0.0);  
    AddHistogram(fSProfile[x][y]);
//-------------------------------------------------
   }}
   
  } // expert
 
} // end void P1aneECAL0::Init(TTree* tree)

//fSProfile
//fSProfile
                                
//-------------------------------------------------------
//-------------------------------------------------------
//-------------------------------------------------------
//-------------------------------------------------------

void PlaneECAL0::EndEvent(const CS::DaqEvent &event) {
//Stat_t bb1;                         
//Stat_t ee1;                         
//Stat_t limit;                            


if(fNhits==0) return;
//if(fNhits>1) return;
        
//cout<<"    fNhits   =   "<<fNhits<<endl;

if (thr_flag) TThread::Lock();
  
      if(strncmp(GetName(),"EC00P1__",8)== 0) {fNrows = 27; fNcols=33;}
      if(strncmp(GetName(),"EC00FEM_",8)== 0) {fNrows = 2; fNcols=2;}
      if(strncmp(GetName(),"EC00SUM_",8)== 0) {fNrows = 2; fNcols=2;}

         const CS::DaqEvent::Header head = event.GetHeader();
         const CS::DaqEvent::EventType eventType = head.GetEventType();
const CS::uint32 Atr = head.GetTypeAttributes()[0]&0xf8000000;
//cout<<"eventType  = "<<eventType<<endl;
        if(eventType == 8 && Atr==0x68000000) {   // LEDs = calibration event ------

                
                
                //int tt = head.GetTriggerNumber();
//cout<<"tt= "<<tt<<endl;
//     if(strncmp(GetName(),"EC00P1__",8)== 0) 
       {
for (int i=0; i<fNhits; i++) {
    int row=fRow[i];
    int col=fCol[i];
    int amp=fAmp[i];
    int adr=fRow[i]+fCol[i]*fNrows; 
//if(fVrow->Test(row) && fVcol->Test(col) && fVamp->Test(amp))
{
//cout<<"ECAL0 LED --> x = "<<col<<"  y = "<<row<<"  led = "<<ECAL0_Led[col][row]<<endl; 
//--------------------------------------------------------------------------------------
//
//   ----------------- LED filling -----------------------------------------------------
//
                           fHchampled->Fill(adr,amp); // original LED values
	                     fSum_led->Fill(amp); 
//                                fProf->Fill(adr,amp);
if (amp>5)                   fAmpLed ->Fill(col,row,amp);                       
//-------------------------------------------------
// bb1 = fProf->GetBinContent(adr);                           
// ee1 = fProf->GetBinError(adr);                           
// ee1 = fProf->GetEntries()/476.;                           
 }//if(fVrow->Test(row) && fVco
 }//for (int i=0; i<fNhits; i++
 }//if(strncmp(GetName(),"HC01P1__",8)== 0)
               
} // if(eventType == 8 && Atr==0x68000000) // LEDs = calibration event -----------------
else 
{

//  int adrmax=0;
  int hcut = 0;
  int sumcut = 0; 
  int ampmax = 0; 
  int xampmax = -1;
  int yampmax = -1;
//--------------------------------------------------------------------------------------

for (int i=0; i<fNhits; i++) {
    int row=fRow[i];    int col=fCol[i];    
    int amp=fAmp[i];

//cout<<"amp = "<<amp<<endl;

///////////////if (amp==0) continue;

if(fVrow->Test(row) && fVcol->Test(col) && fVamp->Test(amp)) {
      int adr=fRow[i]+fCol[i]*fNrows; 

//if(strncmp(GetName(),"HC01P1S1",8)== 0) cerr<<"col= "<<col<<" row= "<<row<<" adr= "<<adr<<" fNrows= "<<fNrows<<endl;

      if (fAmp[i] > ampmax)  {
        ampmax = fAmp[i];
        xampmax = col;
        yampmax = row;
      }
      if (fAmp[i] > fAmpCut) {
        hcut++;     
        sumcut+=fAmp[i];
         fHchac->Fill(adr,amp);
         fHxyac->Fill(col,row,amp);
         fHch->Fill(adr);
      }    
       fHrca->Fill(col,row,amp);
       fHrc1->Fill(col,row);
       fHa->Fill(amp);
      
    fHavsadr->Fill(adr,amp);
      
      fVrow->Store(row);
      fVcol->Store(col);
      fVamp->Store(amp);

      fNhitsKept++;
    }
  //cout<<"amp = "<<amp<<endl;
}

  if(sumcut>0) fHsac->Fill(sumcut);
  if(ampmax>0) fHma->Fill(ampmax);
  fHhth->Fill(hcut);
  fHhit->Fill(fNhitsKept);
  if (ampmax < fAmpCut) fHxyn->Fill(xampmax,yampmax);

}
  if (thr_flag) TThread::UnLock();

}
void PlaneECAL0::StoreDigit(int col, int row, std::vector<float>& data) {

              if(strncmp(GetName(),"EC00P1__",8)== 0) {fNrows = 27; fNcols=33;}
              if(strncmp(GetName(),"EC00FEM_",8)== 0) {fNrows = 2; fNcols=2;}
              if(strncmp(GetName(),"EC00SUM_",8)== 0) {fNrows = 2; fNcols=2;}

  int sadclen=data.size()-2;

// cout<<"sadclen  = "<< sadclen<<endl;
  
    if(sadclen<1) return;
//
 double ped,rms_ped;
 int aaa;
 double amp;
int min_amp;
 int max_amp;
 int i_max;
int i_min;
 int sadc_amp;
 double time_mean,time_mean0,time_sigm;

double a_par,b_par,c_par,y_max;
double x1,x2,x3;
double y1,y2,y3;
int i,i1,i2,i3;

char cell_m[100];

  if (fNhits < fNchan*fMAX_MULT) {
    fRow[fNhits]=row;
    fCol[fNhits]=col;
    int fSADCev[fNsamples];
    int len=fNsamples;
    if(sadclen<len){ len=sadclen;
                    for(i=len;i<fNsamples;i++) fSADCev[i]=0; }
    for(i=0;i<len;i++) {fSADCev[i]=(int)data[i+2];}
//
//------------------------ pedestalls  --------------------------------------------    
//
    ped=0;
    for(i=0;i<pedlen;i++) ped+=fSADCev[i];  
    ped=ped/pedlen;
//cout<<"ped = "<<ped<<endl;    

aaa=row+fNrows*col;
fPed_sum->Fill(ped);    //  total pedestall filling
fPed_ch->Fill(aaa,ped); // pedestall filling
fPed_xy->Fill(col,row,ped); // pedestall x-y filling
sprintf(cell_m,"cell_%d",aaa);
rms_ped=fPed_ch->ProjectionY(cell_m,aaa,aaa)->GetRMS(1);

fRMS_Ped_ch  -> Fill(aaa,rms_ped); // RMS ped filling
fRMS_Ped_sum -> Fill(rms_ped);
fRMS_Ped_xy  -> Fill(col,row,rms_ped);
//
//----------------------------- maximum & minimum samle search -------------
//
   y_max=0;
 min_amp=10000;
 max_amp=0;
   i_max=40;
   i_min=-1;
for(i=pedlen; i<len-pedlen; i++) {     
 sadc_amp=fSADCev[i];
 if(sadc_amp < min_amp )  {min_amp = sadc_amp;  i_min = i;}
 if(sadc_amp > max_amp )  {max_amp = sadc_amp;  i_max = i;}
} 
//
//--------------------------- parabola to obtain time and amp -------------------------
//
// y=a*x**2 + b*x +c ;
// x_max=-b/(2*a);
// y_max=(4*a*c-b*b)/(4*a) ;

i1=i_max-1; i2=i_max;  i3=i_max   + 2;

x1=(double)i1; x2=(double)i2; x3=(double)i3;
y1=(double)fSADCev[i1]; y2=(double)fSADCev[i2]; y3=(double)fSADCev[i3];

a_par=(y3-(x3*(y2-y1)+x2*y1-x1*y2)/(x2-x1))/(x3*(x3-x1-x2)+x1*x2);
b_par=a_par*(x1+x2)-(y2-y1)/(x2-x1);
c_par=(x1*y2-x2*y1)/(x1-x2) + a_par*x2*x1 ;

time_mean0=0;

if(a_par != 0) {
time_mean0=b_par/(2.*a_par);
y_max=(4.0*a_par*c_par - b_par*b_par)/(4.0*a_par) ;
}
//------------------------------------------------------

                fSum_time->Fill(time_mean0);
                fTime_xy->Fill(col,row,time_mean0); 

amp=y_max-ped;

//cout<<"amp = "<<amp<<"   y_max = "<<y_max-ped<<"   c_par = "<<c_par<<endl;
   
//-----------------------------------------------------------------------
 
if (fExpertHistos)
         { 
                
fCH_time[col][row]->Fill(time_mean0);    // fill time_mean0
time_mean=fCH_time[col][row]->GetMean(); // new time arounn mean
time_sigm=fCH_time[col][row]->GetRMS(); // new time arounn mean
 
for(i=0; i<fNsamples; i++)  fSProfile[col][row]->Fill(i,fSADCev[i]-ped);  

  
if(time_mean0<5 ) 
{
fAmp[fNhits]=4000; 
fBad->Fill(amp); 
 return; 
} 
else 
{ 
fGood->Fill(amp); 
fCH_amp[col][row] ->Fill(amp);  
} // time conditions ------------------- t
 
} // end Expert histos  

fAmp[fNhits]=(int)amp;

fNhits++;

  }//if (fNhits <
}// end void

void PlaneECAL0::StoreDigit(CS::Chip::Digit* digit) {
  std::vector<float> data=digit->GetNtupleData();
//  cout<<" data.size()   = "<<data.size()<<endl;
  if(data.size()<3) return;
  const CS::ChipSADC::Digit* sadcdig = dynamic_cast<const CS::ChipSADC::Digit*>(digit);
//  cout<<"  sadcdig  = "<<sadcdig<<endl;
  if(sadcdig != NULL ) {
         int ix=sadcdig->GetX(); int iy=sadcdig->GetY();
         if(ix < 0||iy<0||ix>fNcols||iy>fNrows) return; 
         this->StoreDigit(ix,iy,data);
  }
}
