#ifndef __GroupMuonWallA__
#define __GroupMuonWallA__

#include "Group.h"

#if USE_TRACK == 1
#include "TrackMaker.h"
#endif

class GroupMuonWallA : public Group {
 private:
  static const int fGATE;
  int fMM;
  
#if USE_TRACK == 1
  //    TrackMaker *fTrackMaker;
#endif

 public:
  GroupMuonWallA(const char* name);

#if USE_TRACK == 1
  //    void Init();

  //void ResetHistograms();

  //void Tracking();
#endif

  void ControlPanel(const TGWindow *p, const TGWindow *main);

  ClassDef(GroupMuonWallA,0)
};

#endif

