#ifndef __GroupSili__
#define __GroupSili__

#include "Group.h"
#include "PlaneMISC.h"
#include "PlaneScifiJ.h"

#if USE_TRACK == 1
#include "TrackMaker.h"
#endif

class GroupSili : public Group {

 public:

  GroupSili(const char* name): Group(name) {}
  
  void Init();

#ifndef __CINT__
  void EndEvent(const CS::DaqEvent &event);
#endif
  void ControlPanel(const TGWindow *p, const TGWindow *main);

  ClassDef(GroupSili,0)
};

class GroupSiliTrack : public Group {

 private:
  bool firstread;
#if USE_TRACK == 1
  double tcsphase;
  //TH1 *hChi22;
  TH1 *hdTSci;
  TH2 **hSiliR0VsSciTCS, **hSiliR1VsSciTCS, **hSiliR0VsR1; 
  TH2 *hHitMap, *hHitDMap;
  TrackMaker *fTrackMaker;
  int maxNumberofClusters;
  double scifiTimeDiffCut;
  double scifiTriggerCut;
  double chi2Cut;
  int siliconTimeCutMode;
  int siliconActiveTimeCut;
  int deactivatePlane;
  char *deactivatedPlaneName;
  Tracker **deactivatedPlane;
  std::vector<const char*> silName;
  std::vector<const char*> sciName;
 

  TTree* tsit;
  int   *silMult;    // cluster multiplicity[12] changed to[20]
  float *silResi;  // cluster residual in current track
  float *silPosi;  // cluster position
  float *silAmpl;  // cluster amplitude
  float *silRat0;  // cluster r0
  float *silRat1;  // cluster r1
  float *silT0;    // time from r0
  float *silT1;    // time from r1
  float *silTa;    //
  float TimeFiber;    //
  float TimePhase;    //
  float chi2;         //
  float chi2p;        

#endif
 public:
  GroupSiliTrack(const char* name): Group(name) {firstread=true;}
  
  void Init();

  void ResetHistograms();
  void Tracking();

#ifndef __CINT__
  void EndEvent(const CS::DaqEvent &event);
#endif
  void ControlPanel(const TGWindow *p, const TGWindow *main);
  
  int ntrack;
  double crossing_point[3];
  double slopes[2];
  double crossing_point2[3];
  double slopes2[2];
  
  ClassDef(GroupSiliTrack,0)
};

#endif










