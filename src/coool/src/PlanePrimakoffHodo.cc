#include "PlanePrimakoffHodo.h"

ClassImp(PlanePrimakoffHodo);

PlanePrimakoffHodo::PlanePrimakoffHodo(const char *detname,int nchan, int center, int width)
  :
  Plane1V(detname,nchan,center,width) {

  // change of ranges in time histogram, requested by Andrea Ferrero
  delete(fVt);
  fVt = AddVariable("_t",1000,-4500,-500,fNchan*fMAX_MULT);
}
