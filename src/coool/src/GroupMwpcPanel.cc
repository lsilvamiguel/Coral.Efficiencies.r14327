#include "GroupMwpcPanel.h"

ClassImp(GroupMwpcPanel);

GroupMwpcPanel:: GroupMwpcPanel(const TGWindow *p,const TGWindow *main,
			      UInt_t w, UInt_t h, GroupMwpc *Group)
  : GroupPanel(p, main, w, h, Group), fGroupMwpc(Group) {

  CreateControlTab();
}

GroupMwpcPanel::~GroupMwpcPanel() {
}

void GroupMwpcPanel::CreateControlTab() {


  //control tab
  TGCompositeFrame *tf=fTab->AddTab("Control");
  
  fMain->MapSubwindows();  
  fMain->Resize(fMain->GetDefaultSize());
  fMain->Layout();  
}
