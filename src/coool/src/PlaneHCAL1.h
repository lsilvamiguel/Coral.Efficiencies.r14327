//
#ifndef __PlaneHCAL1__
#define __PlaneHCAL1__

#include "Plane2V.h"
#include "PlaneHCAL1.h"

#include "TTree.h"

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>

#include <TProfile.h>


/// Plane class for HCAL1

class PlaneHCAL1 : public  Plane2V {
  
 private:
//  double pin_norm,led_norm;

  static const int fAmpCut;
  double HCAL1_Led[28][20];
//  double pin_norm,led_norm;

  TH1F *fHsac, *fHch, *fHma, *fHhth;
  TH2F *fHchac, *fHxyac, *fHxyn;

  TH1F *fCHamp[28][20], *fSum_led;
  TH2F *fHrefled;
   TProfile *fProf;
 public:
  
    PlaneHCAL1(const char *detname,int nrow, int ncol, int center, int width):
    Plane2V(detname,nrow, ncol, center, width) {}

  ~PlaneHCAL1() {}

  void Init(TTree* tree =0);
  void  Calib_coef();

#ifndef  __CINT__
  void EndEvent(const CS::DaqEvent &event);
  void StoreDigit(CS::Chip::Digit* digit);
#endif
  void StoreDigit(int col, int row,  std::vector<float>& data);

  ClassDef(PlaneHCAL1,1)
};

#endif




