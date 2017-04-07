#include "GroupMumegaPanel.h"

ClassImp(GroupMumegaPanel);

GroupMumegaPanel:: GroupMumegaPanel(const TGWindow *p,const TGWindow *main,
			      UInt_t w, UInt_t h, GroupMumega *Group)
  : GroupPanel(p, main, w, h, Group), fGroupMumega(Group) {

  CreateControlTab();
}

GroupMumegaPanel::~GroupMumegaPanel() {
}

void GroupMumegaPanel::CreateControlTab() {


  //control tab
//   TGCompositeFrame *tf=fTab->AddTab("Control");
  fTab->AddTab("Control");

  fMain->MapSubwindows();  
  fMain->Resize(fMain->GetDefaultSize());
  fMain->Layout();  
}
