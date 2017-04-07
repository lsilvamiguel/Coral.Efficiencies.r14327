#ifndef __GroupBeamKiller__
#define __GroupBeamKiller__

#include "Group.h"
#include "PlaneBeamKiller.h"

class GroupBeamKiller : public Group {

 public:
  GroupBeamKiller(const char* name): Group(name) {}
  void Init();
#ifndef __CINT__
  void EndEvent(const CS::DaqEvent &event);
#endif
  void ControlPanel(const TGWindow *p, const TGWindow *main);
  ClassDef(GroupBeamKiller,0)
};
#endif










