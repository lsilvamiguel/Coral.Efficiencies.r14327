#ifndef __PlaneHCAL2__
#define __PlaneHCAL2__

#include "Plane2V.h"
#include "PlaneHCAL2.h"

#include "TTree.h"

/// Plane class for HCAL2

class PlaneHCAL2 : public  Plane2V {
  
 private:
  static const int fAmpCut;

  TH1F *fHsac, *fHch, *fHma, *fHhth;
  TH2F *fHchac, *fHxyac, *fHxyn;
  TH2F *hLed, *hLedChnl;

 public:
  
  PlaneHCAL2(const char *detname,int nrow, int ncol, int center, int width):
    Plane2V(detname,nrow, ncol, center, width) {}


  ~PlaneHCAL2() {}

  void Init(TTree* tree =0);

#ifndef  __CINT__
//  void EndEvent(int reftime, const CS::DaqEvent &event);
  void EndEvent(const CS::DaqEvent &event);
  void StoreDigit(CS::Chip::Digit* digit);
#endif
  void StoreDigit(int col, int row,  std::vector<float>& data);

  ClassDef(PlaneHCAL2,1)
};

#endif




