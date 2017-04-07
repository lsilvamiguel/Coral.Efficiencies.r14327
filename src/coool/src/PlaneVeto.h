#ifndef __PlaneVeto__
#define __PlaneVeto__

#include "Plane1V.h"

#include "TTree.h"

/// Plane for beam momentum stations

// maximum number of triggers for the trigger specific histograms

class PlaneVeto : public  Plane1V {
 private:
 static const int fMaxTriggerNumber=11;
  
 public:
  
  PlaneVeto(const char *detname,int nchan, int center, int width):
    Plane1V(detname,nchan,center,width) {}
  
  ~PlaneVeto() {}
  void Init(TTree* tree =0);
#ifndef  __CINT__
  void EndEvent(const CS::DaqEvent &event);
#endif

  TH1F *fHtrigger_spec_time[fMaxTriggerNumber];
  TH2F *fHtrigger_spec_timeVSch[fMaxTriggerNumber];
  

  ClassDef(PlaneVeto,1)
};

#endif




