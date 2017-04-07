#ifndef __PlanePropFrame__
#define __PlanePropFrame__

#include <TGFrame.h>
#include <TGButton.h>
#include <TGTextEntry.h>
#include <TGTab.h>
#include <TGLabel.h>
#include <TGIcon.h>
#include <TGMsgBox.h>

#include "Plane.h"

enum  EPropCommandIdentifiers {
  M_PROPERTIES_CANCEL,
  M_PROPERTIES_OK
};

class PlanePropFrame : public TGTransientFrame  {

 private:
  const TGWindow *fMain;
  Plane *fPlane;
  TGTab               *fTab;
  TGHorizontalFrame *fFrame1;
  TGButton *fOkButton, *fCancelButton;
  TGLayoutHints *fL1, *fL2, *fL3, *fL4, *fL5;

  TGLabel *fVar;
  TGLabel *fBin;
  TGLabel *fMin;
  TGLabel *fMax;
  
  vector<TGLabel*> fCutLabel;

  vector<TGTextEntry*> fCutEntryMin;
  vector<TGTextEntry*> fCutEntryMax;
  vector<TGTextEntry*> fCutEntryBin;

 public:
  PlanePropFrame(const TGWindow *p,const TGWindow *main,
		 UInt_t w, UInt_t h, Plane *plane);
  virtual ~PlanePropFrame();
  
  virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t);
  void CloseWindow();
};

#endif


