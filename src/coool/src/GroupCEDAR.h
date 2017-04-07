#ifndef __GroupCEDAR__
#define __GroupCEDAR__

#include "Group.h"
#include "PlaneCEDAR.h"

class GroupCEDAR : public Group {

 public:
  GroupCEDAR(const char* name): Group(name) {}
  void Init();
#ifndef __CINT__
  void EndEvent(const CS::DaqEvent &event);
#endif
  void ControlPanel(const TGWindow *p, const TGWindow *main);
  ClassDef(GroupCEDAR,0)
};
#endif










