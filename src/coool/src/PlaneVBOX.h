#ifndef __PlaneVBOX__
#define __PlaneVBOX__

#include "PlaneHCAL1.h"
#include "PlaneVBOX.h"

#include "TTree.h"

/// Plane class for VetoBox

class PlaneVBOX : public PlaneHCAL1 {
  
 private:
  int lAmpHistOffset;

 public:
  
  PlaneVBOX(const char *detname,int nrow, int ncol, int center, int width):
    PlaneHCAL1(detname,nrow, ncol, center, width) {}

  ~PlaneVBOX() {}

  void Init(TTree* tree =0);

#ifndef  __CINT__
  void EndEvent(const CS::DaqEvent &event);
#endif

  ClassDef(PlaneVBOX,1)
};

#endif




