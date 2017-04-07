#ifndef __PlaneMuonWall__
#define __PlaneMuonWall__

#include "Plane1V.h"

#include "TTree.h"

/// Plane for beam momentum stations

class PlaneMuonWall : public  Plane1V {
  
 public:
  
  PlaneMuonWall(const char *detname,int nchan, int center, int width):
    Plane1V(detname,nchan,center,width) {}
  
  ~PlaneMuonWall() {}

  ClassDef(PlaneMuonWall,1)
};

#endif




