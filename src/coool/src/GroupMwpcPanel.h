#ifndef __GroupMwpcPanel__
#define __GroupMwpcPanel__

#include "GroupMwpc.h"
#include "GroupPanel.h"

#include <RQ_OBJECT.h>

#include "TThread.h"

/// GroupMwpc control panel 

class GroupMwpcPanel : public GroupPanel {

 private:

  /// associated Group
  GroupMwpc  *fGroupMwpc;

 public:
  GroupMwpcPanel(const TGWindow *p,const TGWindow *main,
		 UInt_t w, UInt_t h, GroupMwpc *Group);
  virtual ~GroupMwpcPanel();
  
  /// Creates one tab for various control buttons 
  void CreateControlTab();

#if ROOT_VERSION_CODE >= ROOT_30200
  RQ_OBJECT("GroupMwpcPanel")  
#else
  RQ_OBJECT()
#endif
  ClassDef(GroupMwpcPanel,0)
};


#endif





