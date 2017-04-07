#ifndef __PlaneHCAL2Sum__
#define __PlaneHCAL2Sum__

#include "Plane1V.h"

#include "TTree.h"

/// Plane for beam momentum stations

class PlaneHCAL2Sum : public  Plane1V {
  
 public:
  
  PlaneHCAL2Sum(const char *detname,int nchan, int center, int width):
    Plane1V(detname,nchan,center,width) {}
  
  ~PlaneHCAL2Sum() {}

  ClassDef(PlaneHCAL2Sum,1)
};

#endif




