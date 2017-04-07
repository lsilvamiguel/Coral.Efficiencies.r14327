#include "GroupMuonWallAPanel.h"

ClassImp(GroupMuonWallAPanel);

GroupMuonWallAPanel:: GroupMuonWallAPanel(const TGWindow *p,const TGWindow *main,
			      UInt_t w, UInt_t h, GroupMuonWallA *Group)
  : GroupPanel(p, main, w, h, Group), fGroupMuonWallA(Group) {

  CreateControlTab();
}

GroupMuonWallAPanel::~GroupMuonWallAPanel() {
}

void GroupMuonWallAPanel::CreateControlTab() {

  //control tab
//   TGCompositeFrame *tf=fTab->AddTab("Control");
  fTab->AddTab("Control");

  fMain->MapSubwindows();  
  fMain->Resize(fMain->GetDefaultSize());
  fMain->Layout();  
}




