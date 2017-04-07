#ifndef __GroupGEM__
#define __GroupGEM__

#include "Group.h"
#include "PlaneTrigHodo.h"
#include "PlaneScifiG.h"

#if USE_TRACK == 1
#include "TrackMaker.h"
#endif

class GroupGEM : public Group {
  
 public:
  GroupGEM(const char* name): Group(name) {}
  
  void Init();

#ifndef __CINT__
  void EndEvent(const CS::DaqEvent &event);
#endif
  void ControlPanel(const TGWindow *p, const TGWindow *main);

  ClassDef(GroupGEM,0)
};

const int GEMi=9;
const float GEMcti[9] = {-200, -75,  -50,  -25,  0.0 ,   12,   25,  50 , 200};
const float GEMcr0[9] = {2.00, 1.68, 1.40, 1.15, 0.35, 0.20, 0.10, 0.05, 0.0};
const float GEMcr1[9] = {1.70, 1.25, 1.15, 1.10, 0.80, 0.65, 0.50, 0.25, 0.0};

class GroupGEMTrack : public Group {

 private:
#if USE_TRACK == 1
  float tcsJitter, triggerTime, triggerTimeHR;
  float scifical[2][200];
  int nSciFi[2];
  float pos0[2], slop[2];
  char sciName[2][10];
  PlaneScifiG* scifiTimePlane[2];
  TH1* hSciTime[2];
  TH2* hSciTimeVsChan[2];
  PlaneTrigHodo *jitterClock;
  TH1 *hTCSjitter,
    *hGEMTime02, *hGEMTime12, *hChi22, *hdTSci;
  TH2 *hGEMRatio02VsSci, *hGEMRatio12VsSci,
    *hGEMTime02VsSci, *hGEMTime12VsSci, *hGEMTimeCorr;
  float dtdr0[GEMi], dtdr1[GEMi];
  TrackMaker *fTrackMaker;
  int maxNumberofClusters;
#endif
  
 public:
  GroupGEMTrack(const char* name): Group(name) {}
  
  void Init();

  void ResetHistograms();
  void Tracking();

#ifndef __CINT__
  void EndEvent(const CS::DaqEvent &event);
#endif
  void ControlPanel(const TGWindow *p, const TGWindow *main);
  
  ClassDef(GroupGEMTrack,0)
};

#endif










