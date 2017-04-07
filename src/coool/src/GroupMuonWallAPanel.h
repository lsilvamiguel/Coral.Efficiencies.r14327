#ifndef __GroupMuonWallAPanel__
#define __GroupMuonWallAPanel__

#include "GroupMuonWallA.h"
#include "GroupPanel.h"

#include <RQ_OBJECT.h>

#include "TThread.h"

/// GroupMuonWallA control panel 

class GroupMuonWallAPanel : public GroupPanel {

 private:

  /// associated Group
  GroupMuonWallA  *fGroupMuonWallA;

 public:
  GroupMuonWallAPanel(const TGWindow *p,const TGWindow *main,
		 UInt_t w, UInt_t h, GroupMuonWallA *Group);
  virtual ~GroupMuonWallAPanel();
  
  /// Creates one tab for various control buttons 
  void CreateControlTab();

#if ROOT_VERSION_CODE >= ROOT_30200
  RQ_OBJECT("GroupMuonWallAPanel")
#else
  RQ_OBJECT()
#endif
  ClassDef(GroupMuonWallAPanel,0)
};


#endif





