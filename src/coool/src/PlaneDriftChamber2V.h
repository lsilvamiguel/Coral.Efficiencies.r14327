#ifndef __PlaneDriftChamber2V__
#define __PlaneDriftChamber2V__

#include "Plane1V.h"

#include "TTree.h"
#include "TThread.h"

/// Plane for drift chambers
class PlaneDriftChamber2V : public Plane1V {
 
 private:
   Variable *fVtemperature;
   Variable *fVport;
   TH2F *fHtemperature;

 public:
  
  PlaneDriftChamber2V( const char *detname,int nchan, int center, int width );

  void NewRun(int runnumber);
  virtual void ControlPanel(const TGWindow* p, const TGWindow* main); 

  void Init( TTree* tree =0 );
  
#ifndef __CINT__
  void StoreDigit(CS::Chip::Digit* digit);

  void EndEvent( const CS::DaqEvent &event );
#endif

  /// Sinica TDC tick in ms
  static const float fSINICA_TICK;
  
  ClassDef(PlaneDriftChamber2V,1)
};

#endif




