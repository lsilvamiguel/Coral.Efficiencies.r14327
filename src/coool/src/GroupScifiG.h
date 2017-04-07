#ifndef __GroupScifiG__
#define __GroupScifiG__

#include "Group.h"

class GroupScifiG : public Group {

 public:
  GroupScifiG(const char* name): Group(name) {}
  void Init();

#ifndef __CINT__
  void EndEvent(const CS::DaqEvent &event);
#endif
  void ControlPanel(const TGWindow *p, const TGWindow *main);

  ClassDef(GroupScifiG,0)
};

#endif
