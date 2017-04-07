#include "GroupScifiJPanel.h"

ClassImp(GroupScifiJPanel);

GroupScifiJPanel:: GroupScifiJPanel(const TGWindow *p,const TGWindow *main,
			      UInt_t w, UInt_t h, GroupScifiJ *Group)
  : GroupPanel(p, main, w, h, Group), fGroupScifiJ(Group) {

  CreateControlTab();
}

GroupScifiJPanel::~GroupScifiJPanel() {
}

void GroupScifiJPanel::CreateControlTab() {


  //control tab
//   TGCompositeFrame *tf=fTab->AddTab("Control");
  fTab->AddTab("Control");

  fMain->MapSubwindows();  
  fMain->Resize(fMain->GetDefaultSize());
  fMain->Layout();  
}
