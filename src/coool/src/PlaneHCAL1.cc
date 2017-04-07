#include <cmath>
  
#include "PlaneHCAL1.h"
#include "ChipSADC.h"
#include "GroupDAQ.h"

#include "DaqEvent.h"

#include <string.h>
#include <cstring>
#include <cstdlib>
#include <sstream>


#include <cstdio>
#include "Plane2V.h"
#include "PlaneHCAL1.h"
#include <iomanip>
#include <iostream>

#include <TProfile.h>
#define  fNsamples 32
#define  pedlen 4

using namespace std;

ClassImp(PlaneHCAL1);
const int PlaneHCAL1::fAmpCut = 30;
//--------------------------------------------------------------------------------------------  
void PlaneHCAL1::Init(TTree* tree) {
//--------------------------------------------------------------------------------------------  
  Plane2V::Init(tree);

  char cell_name[132];
  char cell_name2[132];
  std::string amprowcol;

  fHsac = fHch = fHma = fHhth = fSum_led = NULL;

  fHchac = fHxyac = fHxyn = fHrefled = NULL;
  
  fProf = NULL;

  memset(fCHamp,0,sizeof(fCHamp));

//  cout << "***********************************************************************************************All fNameNew : " << fName.c_str()<< endl;

  if(strncmp(GetName(),"HC01P1__",8)== 0) {
   PlaneHCAL1::Calib_coef();
//  int x,y;
//
  for( int x = 0;x < 28; x++) {
  for( int y = 0;y < 20; y++) {

if (y > 7 && y < 12 && x > 9 && x < 18) continue;
if ( (x==0&&y==0)||(x==0&&y==1)||(x==0&&y==2)||(x==0&&y==3) ) continue;
if ( (x==1&&y==0)||(x==1&&y==1)||(x==1&&y==2)||(x==1&&y==3) ) continue;
if ( (x==2&&y==0)||(x==2&&y==1) ) continue;
if ( (x==3&&y==0)||(x==3&&y==1) ) continue;

if ( (x==24&&y==0)||(x==24&&y==1) ) continue;
if ( (x==25&&y==0)||(x==25&&y==1) ) continue;
if ( (x==26&&y==0)||(x==26&&y==1)||(x==26&&y==2)||(x==26&&y==3) ) continue;
if ( (x==27&&y==0)||(x==27&&y==1)||(x==27&&y==2)||(x==27&&y==3) ) continue;

if ( (x==0&&y==19)||(x==0&&y==18)||(x==0&&y==17)||(x==0&&y==16) ) continue;
if ( (x==1&&y==19)||(x==1&&y==18)||(x==1&&y==17)||(x==1&&y==16) ) continue;
if ( (x==2&&y==19)||(x==2&&y==18) ) continue;
if ( (x==3&&y==19)||(x==3&&y==18) ) continue;

if ( (x==24&&y==19)||(x==24&&y==18) ) continue;
if ( (x==25&&y==19)||(x==25&&y==18) ) continue;
if ( (x==26&&y==19)||(x==26&&y==18)||(x==26&&y==17)||(x==26&&y==16) ) continue;
if ( (x==27&&y==19)||(x==27&&y==18)||(x==27&&y==17)||(x==27&&y==16) ) continue;

 int adr=y + x*20;
sprintf(cell_name2,"%d_Ch",adr);
sprintf(cell_name,"X_%d_Y_%d",x,y);
  fCHamp[x][y] = new TH1F(cell_name2, cell_name, 200, -100, 100);
// fprintf(stderr, "PlaneHCAL1::Init: GetName %s, x %d, y %d, fCHamp %p\n", GetName(), x, y, fCHamp[x][y]);
    AddHistogram(fCHamp[x][y]);
   }}
//
  fSum_led = new TH1F("Sum_led", "Sum_led", 200, -100, 100);
  fSum_led->SetFillColor(30);
  AddHistogram(fSum_led);
//  
  fHrefled = new TH2F("RefLed", "RefLed", 28, 0,28, 20, 0, 20);
   fHrefled->SetOption("BOX");
   fHrefled->SetFillColor(2);
  AddHistogram(fHrefled);
//  
  fProf = new TProfile("prof", "prof", 560, 0, 560);
                       fProf->SetMaximum(100.0);
                        fProf->SetMinimum(-100.0);
            fProf->SetMarkerColor(2);
            fProf->SetMarkerSize(0.4);
            fProf->SetMarkerStyle(20);
  AddHistogram(fProf);
//---------------------------------------------------------------------------------
} // end HC01P1__
  
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
  fHhth=new TH1F_Ref(hitsthname.c_str(),hitsthname.c_str(), 50, 0, 50, fRateCounter);     // like in first histo
  //   fNrows*fNcols, 0, fNrows*fNcols);
  ((TH1F_Ref*)fHhth)->SetReference(fReferenceDirectory);
  fHhth->SetFillColor(7);
  AddHistogram(fHhth);
  
  // 9
  std::string xyacut = fName + "_xVSyAcut";
  fHxyac = new TH2F(xyacut.c_str(), xyacut.c_str(),
                      fNcols, 0, fNcols, fNrows, 0, fNrows);
               fHxyac->SetDrawOption("box");
               fHxyac->SetFillColor(4);
  AddHistogram(fHxyac);

  
  // 10 
  std::string xynoise = fName + "_xVSynoise";
  fHxyn = new TH2F(xynoise.c_str(), xynoise.c_str(),
                      fNcols, 0, fNcols, fNrows, 0, fNrows);
  fHxyn->SetDrawOption("box");
  AddHistogram(fHxyn);

}  // end void PlaneHCAL1::Init(TTree* tree)
//--------------------------------------------------------------------------------------------------------
  void  PlaneHCAL1::Calib_coef() {
//--------------------------------------------------------------------------------------------------------
float hcal_led[28][20]={ {   -1  ,   -1  ,   -1  ,   -1  ,  390.2,  708.2,  431.3,  364.2,  597.9,  378.3,  433.4,  431.7,  725.8,  349.7,  777.7,  763.9,   -1  ,   -1  ,   -1  ,   -1   },
                         {   -1  ,   -1  ,   -1  ,   -1  ,  659.7,  796.6,  159.1,  938.3,  429.5,  727.9,  501.9,  326.3,  522.3,  650.1, 1403.0,  663.9,   -1  ,   -1  ,   -1  ,   -1   },
                         {   -1  ,   -1  ,  515.9,  980.5,  718.3,  394.4,  885.8,  482.4,  602.4,  713.2, 1234.0,  823.1,  356.2,  669.8,  854.6,  303.7,  731.5, 1164.0,   -1  ,   -1   },
                         {   -1  ,   -1  ,  119.6,  794.5,  636.0,  375.7,  467.1,  248.4,  838.2,  989.8,  797.2, 1328.0,  435.6,  197.2,  847.2,  773.0,  935.9,  163.9,   -1  ,   -1   },
                         {  846.3,  415.1, 1095.0,  995.6,  430.5,  397.0,  350.2,  173.5, 1176.0, 1410.0, 1132.0,  888.0, 1051.0, 1271.0,  648.4,  845.0,  740.9, 1566.0,  869.9,  835.8 },
                         {  662.5,  951.2,  493.4,  795.5,  267.8,  223.7,  309.8,  363.8,  831.8,  896.9,  621.9,  399.6,  328.2,  678.1,  685.2,  769.4, 1029.0, 1945.0,  644.3,  733.7 },
                         {  731.3, 1349.0,  763.5,  791.2,  342.6,  692.4,  658.5,  592.0,  489.8,  484.0,  573.2,  572.2,  116.9,  518.2,  285.4,  574.3,  387.6, 1264.0, 1100.0,  609.6 },
                         {  527.7, 1490.0,  827.4,  992.3,  483.6,  425.8,  522.0,  281.9, 1900.0, 1040.0,  675.8, 1266.0,  511.1,  574.1,  303.3,  328.3,  773.9,  678.5,  646.9,  741.3 },
                         {  746.9,  548.2, 2181.0, 1827.0,  452.2,  502.4,  694.6,  561.6,  943.9, 1650.0,  948.4,  266.1,  157.5,  394.1,  687.0,  663.9,  499.8,  635.7,  868.7,  327.1 },
                         {  831.3, 1387.0,  568.0,  686.0,  414.4,  458.7,  846.5,  690.8,  485.7,  460.4,  426.1,  676.1,  482.2,  351.9,  653.8,  500.4,  517.5,  730.5,  347.4,  788.4 },
                         { 1117  ,  233.9,  744.4,  656.1,  617.8,  647.3,  967.5,  967.5,   -1  ,   -1  ,   -1  ,   -1  ,   91.28, 485.9,  171.0,  342.2,  278.0,  327.6,  246.0,  460.1 },
                         {  774.1,  975.3,  451.7,  953.5,  749.5,  393.9,  668.0,  358.6,   -1  ,   -1  ,   -1  ,   -1  ,  419.8,  296.4,  302.1,  550.3,  691.3,  319.9,  420.5,  453.8 },
                         {  610.3,  579.7,  857.0, 1401.0,  457.9,  838.1,  575.2,  896.4,   -1  ,   -1  ,   -1  ,   -1  ,  306.4,  567.7,  263.9,  637.2,  577.8,  338.0,  287.3,  605.3 },
                         {  809.7,  801.6,  769.3,  761.3,  675.3,  960.6,  697.4,  612.9,   -1  ,   -1  ,   -1  ,   -1  ,  481.4,  617.8,  270.1,  292.4,  242.3,  245.6,  519.9,  191.4 },
                         {  214.1,  153.1,  435.7,  184.9,  312.2,  699.9,  494.0,  310.2,   -1  ,   -1  ,   -1  ,   -1  ,  298.4,  262.7,  374.2,  433.9,  313.8,  333.8,  959.3,  749.4 },
                         {  289.7,  385.5,  346.0,  861.8,  372.8,  511.0,  609.6,  470.0,   -1  ,   -1  ,   -1  ,   -1  ,  195.2,  187.3,  590.1,  175.1,  249.8,  483.6,  290.7,  678.4 },
                         {  306.2,  661.9,  765.8,  402.2,  499.7,  331.8,  485.6,  348.1,   -1  ,   -1  ,   -1  ,   -1  ,  360.0,  315.7,  187.1,  203.9,  390.0,  525.4,  580.6,  445.8 },
                         {  645.2,  464.1,  581.6,  592.0,  865.8,  725.2,  456.6,  133.7,   -1  ,   -1  ,   -1  ,   -1  ,  217.7,  277.6,  325.9,  267.8,  593.8,  519.4,  393.8,    0   },
                         {  330.7, 1232.0, 1181.0, 1503.0, 1482.0,  942.6, 1186.0, 1203.0,  683.0,  881.3,  812.8,  953.4,  147.1,  286.0,  264.3,  575.8,  474.3, 1157.0,  873.2,  983.8 },
                         { 1874.0, 1409.0,  853.5, 1554.0, 1084.0, 1019.0,  869.5,  484.9, 1221.0,  721.1,  538.4, 1134.0, 1510.0,  458.0,  756.3,  781.9,  403.9,  986.6,  439.1,  103.3 },
                         {  730.5, 1788.0, 1156.0, 1309.0, 1862.0, 1062.0, 2829.0, 1134.0,  227.2, 1599.0,  791.4,  623.6,  429.4,  953.0,  962.3,  695.4,  690.6,  631.9,  539.4,  622.4 },
                         {    0.0, 1520.0, 1174.0, 1086.0, 1577.0,  701.8, 1064.0,  732.2,  900.9, 1559.0,  924.0,  686.1,  484.1, 1047.0,  159.9,  573.6,  991.7,  571.2,  890.9,  678.1 },
                         {  888.2,  485.4,  434.1,  449.7,  509.9,  431.9,  198.0,  151.5,  300.2,  561.5,  414.0,  369.9,  736.5,  877.3,  941.7,  726.2,  814.5,  550.5,  552.1,  542.5 },
                         {  203.5,  686.1,  781.2, 1177.0,  332.5,  515.3,   79.28, 386.5,  263.9,  623.0,  213.7,  382.4, 1206.0, 1511.0,  640.5,  470.8,  872.7,  463.7,  529.5,  197.6 },
                         {   -1  ,   -1  ,  759.7, 1064.0,  189.2,  279.3,  185.2,  248.0,  367.5,  474.1,  143.4,  421.8,  954.6,  800.4,  932.6,  501.5,  294.1,  549.9,   -1  ,   -1   },
                         {   -1  ,   -1  ,  965.7,  695.6,  382.1,  436.0,  375.8,  289.5,  701.9,  531.7,  328.4,  877.0, 1751.0, 1346.0, 1906.0, 1348.0, 1061.0,  415.3,   -1  ,   -1   },
                         {   -1  ,   -1  ,   -1  ,   -1  , 1632.0,  929.6, 1156.0,  515.8,  959.9,  632.2,  790.3,  449.4,  479.2,  372.1, 1147.0,  348.2,   -1  ,   -1  ,   -1  ,   -1   },
                         {   -1  ,   -1  ,   -1  ,   -1  ,  640.5,  677.7,  346.2,  281.9,  717.7,  378.2, 1097.0,  462.5,  493.2,  422.3,  961.7,  755.7,   -1  ,   -1  ,   -1  ,   -1   } };
  double pin_norm,led_norm;
  pin_norm=1.;
  led_norm=1.;
  int m,n;
   char a[256];
   float a_amp,s1,s2,Ecoef,s4;

//   FILE *fi = fopen("/online/mnt/detector/calo/hcal1/ref_led/led_50049.dat", "r");
// THIS HAS TO GO TO A CALIBRATION FILE !!!
  FILE *fi = fopen("/online/mnt/detector/calo/hcal1/ref_led/ref_led.dat", "r");

if(!fi ) {
cerr <<"================ File HCAL_LED.dat not found =====" << "\n";
cout <<"================ LED spectra - From Table   =====" << "\n" << endl;
 for(int i=0; i<28; i++){
 for(int j=0; j<20; j++){
    HCAL1_Led[i][j]  = hcal_led[i][j]; 
//cout<<"------------- HCALl LEDs  ------>  "<<HCAL1_Led[i][j]<<endl;    
                                      }}
} else {

cout<<"================ File led_50049.dat from "<<"/online/mnt/detector/calo/hcal1/ref_led/"<<endl;


 fscanf(fi,"%s %s %s %s %s %s %s", a ,a ,a ,a ,a ,a ,a  );
// cout<<a1<<endl;
 fscanf(fi,"%s %s", a ,a  );
// cout<<a1<<endl;
 fscanf(fi,"%s %s %s %s %s %s %s", a ,a ,a ,a ,a ,a ,a  );
// cout<<a1<<endl;
 fscanf(fi,"%s %s", a ,a  );
// cout<<a1<<endl;
 fscanf(fi,"%s %s %s", a ,a ,a  );
// cout<<a1<<endl;
 
 for(int i=0; i<28; i++){
 for(int j=0; j<20; j++){
 fscanf(fi,"%d %d %f %f %f %f %f", &m ,&n ,&a_amp ,&s1 ,&s2 ,&Ecoef ,&s4 );

//cout<<"&m ,&n ,&a_amp ,&s1 ,&s2 ,&Ecoef ,&s4 ="<<m<<" "<<n<<" "<<a_amp<<" "<<s1<<" "<<s2<<" "<<Ecoef<<" "<<s4<<endl;

    HCAL1_Led[i][j]= a_amp; 
//cout<<"------------- HCAL1 LEDs ------>  "<<a_amp<<endl;    
    
    }}
}   
 }




void PlaneHCAL1::EndEvent(const CS::DaqEvent &event) {
float ampl = -1;

if(fNhits==0) return;
  if (thr_flag) TThread::Lock();

      if(strncmp(GetName(),"HC01P1S1",8)== 0) fNrows = 5;
      if(strncmp(GetName(),"HC01P1S2",8)== 0) fNrows = 4;
      if(strncmp(GetName(),"HC01P1S3",8)== 0) fNrows = 5;
      if(strncmp(GetName(),"HC01P1S4",8)== 0) fNrows = 4;



         const CS::DaqEvent::Header head = event.GetHeader();
         const CS::DaqEvent::EventType eventType = head.GetEventType();
const CS::uint32 Atr = head.GetTypeAttributes()[0]&0xf8000000;

        if(eventType == 8 && Atr==0x68000000) { 




//int tt = head.GetTriggerNumber();
//cout<<"tt= "<<tt<<endl;

      if(strncmp(GetName(),"HC01P1__",8)== 0)  {
for (int i=0; i<fNhits; i++) {
    int row=fRow[i];    int col=fCol[i];    float amp=fAmp[i];

if(fVrow->Test(row) && fVcol->Test(col) && fVamp->Test(amp)) {
      int adr=fRow[i] + fCol[i]*fNrows;
      if (! fCHamp[col][row]) continue; // added by DN 29/6/2011 as some strange col/row values can appear...
if (row == 0  && col == 26) continue;
    if ( HCAL1_Led[col][row] > 0.0 ) ampl = 10.0*(amp - HCAL1_Led[col][row])/HCAL1_Led[col][row] - 21.2;
// else ampl=-70.0; // adr=359 - is sweetched  of!
//cerr<<"col= "<<col<< "row= "<<row<<" adr= "<<adr<<endl;
// fprintf(stderr, "PlaneHCAL1::Init: col %d, row %d, ampl %f, fCHamp %p\n", col, row, ampl, fCHamp[col][row]);
                          fCHamp[col][row] ->Fill(ampl);        
	                        fSum_led->Fill(ampl); 
        float Mean = fCHamp[col][row]->GetMean();
if (Mean>=100.0) Mean = 99.0; // Too big ampl
if (Mean<=-100.0) Mean = -99.0; // HV is switched off!
// Fill the profile hist
//    if(k==1) { k=-1;
//                       fProf->SetMaximum(100.0);
//                        fProf->SetMinimum(-100.0);
             
//            fProf->SetMarkerColor(2);
//            fProf->SetMarkerSize(0.4);
//            fProf->SetMarkerStyle(23);
//             }            
fProf->Fill(adr,Mean);






 if (fabs(Mean)>20) {
                    fHrefled->Fill(col,row);
            //         fHrefled->SetMarkerColor(2);
            //         fHrefled->SetMarkerSize(0.2);
            //         fHrefled->SetMarkerStyle(23);
              }                     
                    
                                                             }
                            }
                                                }
               
                            } else {

//  int adrmax=0;
  int hcut = 0;
 float sumcut = 0; 
 float ampmax = 0; 
  int xampmax = -1;
  int yampmax = -1;
//------------------------------------------------------------------------------------------------

for (int i=0; i<fNhits; i++) {
    int row=fRow[i];    int col=fCol[i];    float amp=fAmp[i];

if(fVrow->Test(row) && fVcol->Test(col) && fVamp->Test(amp)) {
      int adr=fRow[i] + fCol[i]*fNrows; 



//      if(strncmp(GetName(),"HC01P1S1",8)== 0) cerr<<"col= "<<col<<" row= "<<row<<" adr= "<<adr<<" fNrows= "<<fNrows<<endl;






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
  }
  fHsac->Fill(sumcut);
  fHma->Fill(ampmax);
  fHhth->Fill(hcut);
  fHhit->Fill(fNhitsKept);
  if (ampmax < fAmpCut) fHxyn->Fill(xampmax,yampmax);

}
   if (thr_flag) TThread::UnLock();

}
void PlaneHCAL1::StoreDigit(int col, int row, std::vector<float>& data) {
//
  int sadclen=data.size()-2;
  if(sadclen<1) return;
//
  if (fNhits < fNchan*fMAX_MULT) {
    fRow[fNhits]=row;
    fCol[fNhits]=col;
    int fSADCev[fNsamples];
    int len=fNsamples;
    if(sadclen<len){ len=sadclen;
                    for(int i=len;i<fNsamples;i++) fSADCev[i]=0; }
    for(int i=0;i<len;i++) fSADCev[i]=(int)data[i+2];
    double ped=0;
    for(int i=0;i<pedlen;i++) ped+=fSADCev[i];   ped=ped/pedlen;
    double amp=0;
    for(int i=pedlen;i<len;i++) amp+=fSADCev[i]-ped;
    fAmp[fNhits]=(int)amp;
    fNhits++;
  }
}

void PlaneHCAL1::StoreDigit(CS::Chip::Digit* digit) {
  std::vector<float> data=digit->GetNtupleData();
  if(data.size()<3) return;
  const CS::ChipSADC::Digit* sadcdig = dynamic_cast<const CS::ChipSADC::Digit*>(digit);
  if(sadcdig != NULL ) {
         int ix=sadcdig->GetX(); int iy=sadcdig->GetY();
         if(ix < 0||iy<0||ix>fNcols||iy>fNrows) return; 
         this->StoreDigit(ix,iy,data);
  }
}
