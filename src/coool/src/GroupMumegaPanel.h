#ifndef __GroupMumegaPanel__
#define __GroupMumegaPanel__

#include "GroupMumega.h"
#include "GroupPanel.h"

#include <RQ_OBJECT.h>

#include "TThread.h"

/// GroupMumega control panel 

class GroupMumegaPanel : public GroupPanel {

 private:

  /// associated Group
  GroupMumega  *fGroupMumega;

 public:
  GroupMumegaPanel(const TGWindow *p,const TGWindow *main,
		 UInt_t w, UInt_t h, GroupMumega *Group);
  virtual ~GroupMumegaPanel();
  
  /// Creates one tab for various control buttons 
  void CreateControlTab();

#if ROOT_VERSION_CODE >= ROOT_30200
  RQ_OBJECT("GroupMumegaPanel")  
#else
  RQ_OBJECT()
#endif
  ClassDef(GroupMumegaPanel,0)
};


#endif





