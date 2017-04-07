#ifndef __PlaneRCath__
#define __PlaneRCath__

#include "Plane2V.h"

#include "TTree.h"

/// Plane for RICH cathodes

class PlaneRCath : public  Plane2V {
  
 public:
  
  PlaneRCath(const char *detname,int nrow, int ncol, int center, int width):
    Plane2V(detname,nrow, ncol, center, width) {}
  
  ~PlaneRCath() {}

  void Init(TTree* tree =0);

#ifndef __CINT__
  void EndEvent(const CS::DaqEvent &event);
#endif
  ClassDef(PlaneRCath,1)
};

#endif




