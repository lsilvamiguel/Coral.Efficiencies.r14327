#ifndef __PlaneTrigHodo__
#define __PlaneTrigHodo__

#include "Plane1V.h"

#include "TTree.h"

/// Trigger Hodoscope plane

class PlaneTrigHodo : public  Plane1V {
  
 public:
  
  PlaneTrigHodo(const char *detname,int nchan, int center, int width):
    Plane1V(detname,nchan,center,width) {}
  
  ~PlaneTrigHodo() {}
  
  ClassDef(PlaneTrigHodo,1)
};

#endif




