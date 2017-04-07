#ifndef __GroupMumega__
#define __GroupMumega__

#include "Group.h"

#if USE_TRACK == 1
#include "TrackMaker.h"
#endif

class GroupMumega : public Group {

 private:
  static const int fGATE;
  int fMM;
  
#if USE_TRACK == 1
  TrackMaker *fTrackMaker;
#endif

 public:
  GroupMumega(const char* name);

  void Init();

  void ResetHistograms();

  void Tracking();

  void ControlPanel(const TGWindow *p, const TGWindow *main);

  ClassDef(GroupMumega,0)
};

#endif
