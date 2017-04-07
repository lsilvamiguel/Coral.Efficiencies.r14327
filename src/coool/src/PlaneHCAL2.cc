#include "PlaneHCAL2.h"
#include "ChipSADC.h"
#define  fNsamples 32
#define  pedlen 4

ClassImp(PlaneHCAL2);

const int PlaneHCAL2::fAmpCut = 20;

void PlaneHCAL2::Init(TTree* tree) {
  
  Plane2V::Init(tree);
  
  std::string ampname = fName + "_Sumampcut";
  fHsac=new TH1F(ampname.c_str(),ampname.c_str(),
		 fVamp->GetNbins(), fVamp->GetMin(),
		 fVamp->GetMax());
  AddHistogram(fHsac);

  std::string adrname = fName + "_adr";
  fHch=new TH1F( adrname.c_str(),adrname.c_str(), fNrows*fNcols, 0, fNrows*fNcols);
  AddHistogram(fHch);
   
  std::string aadrcut = fName + "_adrVSampcut";
  fHchac = new TH2F(aadrcut.c_str(), aadrcut.c_str(), 
		      fNrows*fNcols, 0, fNrows*fNcols,
		      fVamp->GetNbins(), fVamp->GetMin(), 
		      fVamp->GetMax());
  AddHistogram(fHchac);
  
  // 7    
  std::string ampmaxname = fName + "_maxamp";
  fHma=new TH1F(ampmaxname.c_str(),ampmaxname.c_str(),
                       fVamp->GetNbins(), fVamp->GetMin(),
                       fVamp->GetMax());
  AddHistogram(fHma);

  // 8    
  std::string hitsthname = fName + "_hitsth";
  fHhth=new TH1F(hitsthname.c_str(),hitsthname.c_str(), 50, 0, 50);     // like in first histo
  //   fNrows*fNcols, 0, fNrows*fNcols);
  AddHistogram(fHhth);
  
  // 9
  std::string xyacut = fName + "_xVSyAcut";
  fHxyac = new TH2F(xyacut.c_str(), xyacut.c_str(),
                      fNcols, 0, fNcols, fNrows, 0, fNrows);
  AddHistogram(fHxyac);

  
  // 10 
  std::string xynoise = fName + "_xVSynoise";
  fHxyn = new TH2F(xynoise.c_str(), xynoise.c_str(),
                      fNcols, 0, fNcols, fNrows, 0, fNrows);
  AddHistogram(fHxyn);

  std::string name, title;
  name = fName + "Leds_vs_XY";
  title = fName + " Leds vs XY, amp>thr";
  hLed=new TH2F(name.c_str(),title.c_str(),fNcols,0.,fNcols,fNrows,0.,fNrows);
  hLed->SetOption("col");
  AddHistogram(hLed);

  name = fName + "Led_vs_chnl#";
  title = fName + "LED vs channel#";
  hLedChnl=new TH2F(name.c_str(),title.c_str(), 
                            fNrows*fNcols, 0., fNrows*fNcols,200,0.,4100.);
  AddHistogram(hLedChnl);

}

void PlaneHCAL2::EndEvent( const CS::DaqEvent &event) {

  int evtype=event.GetType();
  if(!(evtype==7||evtype==8))  return;
//  if(evtype!=7) std::cout<<"===== "<<evtype<<std::endl;

  if (thr_flag) TThread::Lock();

  int hcut = 0;
  int sumcut = 0; 
  int ampmax = 0; 
  int xampmax = -1;
  int yampmax = -1;
  
  
  for (int i=0; i<fNhits; i++) {
    
    int row=fRow[i];
    int col=fCol[i];
    int amp=fAmp[i];
    int adr=fRow[i] + fCol[i]*fNrows; 

    if(evtype==8) {                                  // Leds
         if(amp>fAmpCut) hLed->Fill(col,row,amp);
         hLedChnl->Fill(adr,amp);
         continue;
    }

    if(fVrow->Test(row) &&
       fVcol->Test(col) &&
       fVamp->Test(amp)) {    

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
      }
    
      fHrca->Fill(col,row,amp);
      fHrc1->Fill(col,row);
      fHa->Fill(amp);

      fHch->Fill(adr);
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
  
  if (thr_flag) TThread::UnLock();
}
void PlaneHCAL2::StoreDigit(int col, int row, std::vector<float>& data) {
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
    for(int i=0;i<len;i++) {fSADCev[i]=(int)data[i+2];}
    double ped=0;
    for(int i=0;i<pedlen;i++) ped+=fSADCev[i];   ped=ped/pedlen;
    double amp=0;
    for(int i=pedlen;i<len;i++) amp+=fSADCev[i]-ped;
    fAmp[fNhits]=(int)amp;
    fNhits++;
  }
}

void PlaneHCAL2::StoreDigit(CS::Chip::Digit* digit) {
  std::vector<float> data=digit->GetNtupleData();
  if(data.size()<3) return;
  const CS::ChipSADC::Digit* sadcdig = dynamic_cast<const CS::ChipSADC::Digit*>(digit);
  if(sadcdig != NULL ) {
         int ix=sadcdig->GetX(); int iy=sadcdig->GetY();
         if(ix < 0||iy<0||ix>fNcols||iy>fNrows) return; 
         this->StoreDigit(ix,iy,data);
  }
}

