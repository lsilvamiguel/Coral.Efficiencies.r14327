#ifndef __GroupMwpc__
#define __GroupMwpc__

#include "Group.h"

#if USE_TRACK == 1
#include "TrackMaker.h"
#endif

class TH1F;
class TH2F;
class TH1;

class GroupMwpc : public Group {

 private:
  static const int fGATE;
  int fMM;

  int fNevent;
  int fNstation;
  TH2F* fHxvsy;
  TH1F* fHrecErr, *fHnRecEv, *fHx, *fHy;
  TH1F* fHxEff, *fHuEff, *fHvEff;
  TH1*  fEvDisplay;
#if USE_TRACK == 1
  TrackMaker *fTrackMaker;
#endif

 public:
  GroupMwpc(const char* name): Group(name) {}

  void Init();
  void ResetHistograms();
  void Tracking();

  ClassDef(GroupMwpc,0)
};

#endif
