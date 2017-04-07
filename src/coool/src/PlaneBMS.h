#ifndef __PlaneBMS__
#define __PlaneBMS__

#include "Plane1V.h"

#include <sstream>
#include <string>

#include "TTree.h"

/// Plane for beam momentum stations

// maximum number of triggers for the trigger specific histograms

class PlaneBMS : public  Plane1V {
 private:
 static const int fMaxTriggerNumber=12;
 TH1F *fHtrigger_spec_time[fMaxTriggerNumber];
 TH1F *fBMS3_repaired_timing; 


 public:
  
  PlaneBMS(const char *detname,int nchan, int center, int width):
    Plane1V(detname,nchan,center,width) {}
  
  ~PlaneBMS() {}
  void Init(TTree* tree =0);
#ifndef  __CINT__
  void EndEvent(const CS::DaqEvent &event);
#endif

   ClassDef(PlaneBMS,1)
};

#endif




