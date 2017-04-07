 //  version 29 Sept 2012 for ECAL0
 //
#ifndef __PlaneECAL0__
#define __PlaneECAL0__

#include "Plane2V.h"
#include "PlaneECAL0.h"

#include "TTree.h"

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
  
#include <TProfile.h>
   
   
/// Plane class for ECAL0

class PlaneECAL0 : public  Plane2V {
   
 private:
  
  static const int fAmpCut;
  double ECAL0_Led[33][27];

  TH1F *fHsac, *fHch, *fHma, *fHhth, *fBad, *fGood;
  TH2F *fHchac, *fHxyac, *fHxyn;

  TH1F *fCH_amp[33][27], *fSum_led, *fSum_time;
  TH1F *fCH_time[33][27];
     
  //TH2F *fHrefled_low,*fHrefled_hig;
  TH2F *fAmpLed,*fHchampled;
  TH2F *fPed_xy;  
  TH2F *fTime_xy; 
  TH2F *fPed_ch; 
  TH2F *fRMS_Ped_xy;  
  TH2F *fRMS_Ped_ch;

  TH1F *fPed_sum;  
  TH1F *fRMS_Ped_sum;  
   
  TProfile *fProf;
//  TProfile *fSProfile[33][27];
  TH2F *fSProfile[33][27];
   
 public:
  
    PlaneECAL0(const char *detname,int nrow, int ncol, int center, int width):
    Plane2V(detname,nrow, ncol, center, width) {}

  ~PlaneECAL0() {}

  void Init(TTree* tree =0);
//  void  Calib_coef();
  void StoreDigit(int col, int row,  std::vector<float>& data);

  #ifndef  __CINT__
  void EndEvent(const CS::DaqEvent &event);
  void StoreDigit(CS::Chip::Digit* digit);
#endif


  ClassDef(PlaneECAL0,1)
};

#endif
