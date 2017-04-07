#ifndef __GroupScifiJ__
#define __GroupScifiJ__

#include "Group.h"

class GroupScifiJ : public Group {

 public:
  GroupScifiJ(const char* name): Group(name) {}
  void Init();

#ifndef __CINT__
  void EndEvent(const CS::DaqEvent &event);
#endif
  void ControlPanel(const TGWindow *p, const TGWindow *main);

  ClassDef(GroupScifiJ,0)
};

#endif
