#ifndef __PlaneVBHodo__
#define __PlaneBVHodo__

#include "Plane1V.h"

#include "TTree.h"

/// Veto Box inner Hodoscope plane

class PlaneVBHodo : public  Plane1V {
  
 public:
  
  PlaneVBHodo(const char *detname,int nchan, int center, int width):
    Plane1V(detname,nchan,center,width) {}
  
  ~PlaneVBHodo() {}
  
  ClassDef(PlaneVBHodo,1)
};

#endif




