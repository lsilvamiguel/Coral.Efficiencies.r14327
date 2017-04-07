#ifndef __GroupPanel__
#define __GroupPanel__

#include "Group.h"

#include "TThread.h"

#include <RQ_OBJECT.h>
#include <TGFrame.h>
#include <TGButton.h>
#include <TGTextEntry.h>
#include <TGTab.h>
#include <TGLabel.h>
#include <TGIcon.h>
#include <TGMsgBox.h>

/// base class for Group panels
class GroupPanel : public TObject {

 protected:

  /// Associated Group
  Group *fGroup; 

  TGTransientFrame    *fMain;
  TGTab               *fTab;
  TGHorizontalFrame *fFrame1;
  TGButton *fCutOkButton, *fHistOkButton, *fCancelButton;
  TGLayoutHints *fL1, *fL2, *fL3, *fL4, *fL5, *fL6;
  
  TGLabel *fCutVar;
  TGLabel *fCutMin;
  TGLabel *fCutMax;

  TGLabel *fHistVar;
  TGLabel *fHistBin;
  TGLabel *fHistMin;
  TGLabel *fHistMax;

  std::vector<TGLabel*> fCutLabel, fHistLabel;
  std::vector<TGTextEntry*> fCutEntryMin,fHistEntryMin;
  std::vector<TGTextEntry*> fCutEntryMax,fHistEntryMax;
  std::vector<TGTextEntry*> fCutEntryBin,fHistEntryBin;
  TGCompositeFrame *fHistFrame,*fCutFrame;

 public:

  GroupPanel(const TGWindow *p,const TGWindow *main,
		 UInt_t w, UInt_t h, Group *Group);
  virtual ~GroupPanel();
  
  // Slots  
  /// fMain->CloseWindow() connected here                             
  void CloseWindow() {fGroup->PanelClosed(); delete this;}
  
  /// all basic TGButton::Clicked() connected here
  void HandleButtons(int id=-1);             
  
#if ROOT_VERSION_CODE >= ROOT_30200
  RQ_OBJECT("GroupPanel")
#else
  RQ_OBJECT()
#endif
  ClassDef(GroupPanel,0)
};

#endif





