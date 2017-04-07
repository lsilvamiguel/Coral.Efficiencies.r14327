#ifndef __GroupScifiJPanel__
#define __GroupScifiJPanel__

#include "GroupScifiJ.h"
#include "GroupPanel.h"

#include <RQ_OBJECT.h>

#include "TThread.h"

/// GroupScifiJ control panel 

class GroupScifiJPanel : public GroupPanel {

 private:

  /// associated Group
  GroupScifiJ  *fGroupScifiJ;

 public:
  GroupScifiJPanel(const TGWindow *p,const TGWindow *main,
		 UInt_t w, UInt_t h, GroupScifiJ *Group);
  virtual ~GroupScifiJPanel();
  
  /// Creates one tab for various control buttons 
  void CreateControlTab();

#if ROOT_VERSION_CODE >= ROOT_30200
  RQ_OBJECT("GroupScifiJPanel")  
#else
  RQ_OBJECT()
#endif
  ClassDef(GroupScifiJPanel,0)
};


#endif





