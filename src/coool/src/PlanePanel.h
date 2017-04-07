#ifndef __PlanePanel__
#define __PlanePanel__

#include "Plane.h"

#include "TThread.h"

#include <RQ_OBJECT.h>
#include <TGFrame.h>
#include <TGButton.h>
#include <TGTextEntry.h>
#include <TGTab.h>
#include <TGLabel.h>
#include <TGIcon.h>
#include <TGMsgBox.h>
#include <TGSlider.h>
#include <TRootEmbeddedCanvas.h>
#include <TGComboBox.h>

#include "BitButton.h"

/// base class for Plane panels
class PlanePanel : public TObject {

 private:
  unsigned int fTMPattern;

 protected:

  /// Associated Plane
  Plane *fPlane;

  TGTransientFrame    *fMain;
  TGTab               *fTab;
  TGHorizontalFrame *fFrame1, *fFrame2, *fFrame3, *fFrame4;
  TGButton *fCutOkButton, *fHistOkButton, *fTrigOkButton, *fCancelButton;
  TGLayoutHints *fL1, *fL2, *fL3, *fL4, *fL5, *fL6, *fL7, *fL8, *fL9;

  TGCanvas *fHistCanvas, *fCutCanvas;

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

  TGComboBox           *fTTBox;
  TGTextEntry          *fHexaMask;
  TGCompositeFrame     *fTrigMask1;
  TGGroupFrame         *fTrigMask;
  TGGroupFrame         *fTrigType;
  std::vector<TGGroupFrame*> fByte;
  std::vector<BitButton*>    fBitButton;

 public:

  PlanePanel(const TGWindow *p,const TGWindow *main,
		 UInt_t w, UInt_t h, Plane *plane);
  virtual ~PlanePanel();

  // Slots
  /// fMain->CloseWindow() connected here
  void CloseWindow() {fPlane->PanelClosed(); delete this;}

  /// all basic TGButton::Clicked() connected here
  void HandleButtons(int id=-1);

  /// fHexaMask ReturnPressed() connected here
  bool ConvertHexa2Bin();

  /// fBitButton Clicked() connected here
  bool ConvertBin2Hexa();

  /// fTTBox Selected(int) connected here
  void TTSelected(int id);

#if ROOT_VERSION_CODE >= ROOT_30200
  RQ_OBJECT("PlanePanel")
#else
  RQ_OBJECT()
#endif
  ClassDef(PlanePanel,1)
};

#endif





