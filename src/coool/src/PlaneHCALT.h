#ifndef __PlaneHCALT__
#define __PlaneHCALT__

#include "Plane2V.h"

#include "TTree.h"

/// Plane class for HCAL1

class PlaneHCALT : public  Plane2V {

 public:
  
  PlaneHCALT(const char *detname,int nrow, int ncol, int center, int width):
    Plane2V(detname,nrow, ncol, center, width) {}

  ~PlaneHCALT() {}

  void Init(TTree* tree =0);
#ifndef  __CINT__
  void StoreDigit(CS::Chip::Digit* digit); 
  void EndEvent(const CS::DaqEvent &event);
#endif

  ClassDef(PlaneHCALT,1)
};

#endif




