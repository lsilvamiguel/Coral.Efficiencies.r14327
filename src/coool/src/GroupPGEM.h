#ifndef __GroupPGEM__
#define __GroupPGEM__

#include "Group.h"
#include "PlaneTrigHodo.h"
#include "PlaneScifiG.h"

#if USE_TRACK == 1
#include "TrackMaker.h"
#endif

class GroupPGEM : public Group {
  
 public:
  GroupPGEM(const char* name): Group(name) {}
  
  void Init();

#ifndef __CINT__
  void EndEvent(const CS::DaqEvent &event);
#endif
  void ControlPanel(const TGWindow *p, const TGWindow *main);

  ClassDef(GroupPGEM,0)
};

#endif










