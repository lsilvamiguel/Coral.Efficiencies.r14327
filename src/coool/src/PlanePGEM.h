#ifndef __PlanePGEM__
#define __PlanePGEM__

#include "PlaneAPV.h"
#include "TThread.h"

#include "PlaneGEM.h"

/// Plane for PixelGEM detectors (readout with APV chips).

class PlanePGEM : public PlaneGEM {
  
 private:
  // hit histograms
  TH2F* fHhitMapPix;                        // hit map in case of pixel plane
  TH2F* fHoccupPix;                         // occupancy in case of pixel plane
  TH2F* fHhitMapSumAmp;
  TH2F* fHhitMapMeanAmp;                    // mean hit amplitude per pixel

  // cluster variables
  Variable* fVcPosPixX;                     // cluster position in case of pixel plane
  Variable* fVcPosPixY;

  // cluster histograms
  TH2F *fHcPosPix,     *fHcPosPixT;         // cluster position in case of pixel plane
  TH2F *fHcSumAmpPix,  *fHcSumAmpPixT;
  TH2F *fHcMeanAmpPix, *fHcMeanAmpPixT;     // mean cluster amplitude

 public:
  
  PlanePGEM(const char *detname,int nchan, int center, int width, bool pixel);
  
#ifndef  __CINT__
  virtual void EndEvent(const CS::DaqEvent &event);
#endif 

  virtual void Init(TTree* tree =0);

  virtual void Clusterize();
  
  ClassDef(PlanePGEM,0)
};

#endif










