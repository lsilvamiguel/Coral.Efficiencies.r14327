#ifndef __GroupRW__
#define __GroupRW__

#include "Group.h"

class TH1F;

class GroupRichWall : public Group {

 private:
  TH1F* fHHitMult;
 public:
  GroupRichWall(const char* name);

  void Init();

#ifndef __CINT__
  /// Fills the histograms
  virtual void EndEvent(const CS::DaqEvent &event);
#endif

  ClassDef(GroupRichWall,1)
};

#endif
