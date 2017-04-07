#ifndef __PlanePrimakoffHodo__
#define __PlanePrimakoffHodo__

#include "Plane1V.h"

#include "TTree.h"

/// Primakoff Hodoscope plane

class PlanePrimakoffHodo : public  Plane1V {
  
 public:
  
  PlanePrimakoffHodo(const char *detname,int nchan, int center, int width);
  
  ~PlanePrimakoffHodo() {}

  ClassDef(PlanePrimakoffHodo,1)
};

#endif




