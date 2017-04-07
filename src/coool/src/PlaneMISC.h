#ifndef __PlaneMISC__
#define __PlaneMISC__

#include "Plane1V.h"

#include <vector>
#include "TF1.h"
#include "TTree.h"

/// Trigger Hodoscope plane

class PlaneMISC : public  Plane1V {

 public:

  PlaneMISC(const char *detname,int nchan, int center, int width);
  void Init(TTree *tree);
  void Reset();

#ifndef __CINT__
  void EndEvent(const CS::DaqEvent &event);
#endif

  ~PlaneMISC() {};

  ClassDef(PlaneMISC,1)
};

#endif


















