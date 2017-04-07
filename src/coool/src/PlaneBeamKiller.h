#ifndef __PlaneBeamKiller__
#define __PlaneBeamKiller__

#include "Plane1V.h"

#include "TTree.h"

/// Beam Killer plane

class PlaneBeamKiller : public  Plane1V {
  
 public:
  
  PlaneBeamKiller(const char *detname,int nchan, int center, int width):
    Plane1V(detname,nchan,center,width) {}
  
  ~PlaneBeamKiller() {}
  
  ClassDef(PlaneBeamKiller,1)
};

#endif




