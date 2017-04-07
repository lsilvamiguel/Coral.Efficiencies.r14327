#ifndef __GroupDAQ__
#define __GroupDAQ__

#include "Group.h"

class GroupDAQ : public Group {

 public:
  GroupDAQ(const char* name) : Group(name) {}
  void Init();

#ifndef __CINT__
  void EndEvent(const CS::DaqEvent &event);
#endif
  ClassDef(GroupDAQ,0)
};

#endif
