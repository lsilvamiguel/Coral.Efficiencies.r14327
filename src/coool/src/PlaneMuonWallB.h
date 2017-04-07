#ifndef __PlaneMuonWallB__
#define __PlaneMuonWallB__

#include "Plane1V.h"

#include "TTree.h"

/// Plane for beam momentum stations

class PlaneMuonWallB : public  Plane1V {
  
 public:
  
  PlaneMuonWallB(const char *detname,int nchan, int center, int width):
    Plane1V(detname,nchan,center,width) {}
  
  ~PlaneMuonWallB() {}

  ClassDef(PlaneMuonWallB,1)
};

#endif




