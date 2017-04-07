#ifndef __PlaneECAL1Panel__
#define __PlaneECAL1Panel__

#include "PlaneECAL1.h"
#include "PlanePanel.h"

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

#include "TThread.h"

/// PlaneECAL1 control panel

class PlaneECAL1Panel : public PlanePanel {

 private:

  /// Associated Plane
  //PlaneECAL1 *fPlane;

  //TGTransientFrame    *fMain;
  //TGTab               *fTab;
  //TGHorizontalFrame *fFrame1;
  //TGButton *fCutOkButton, *fHistOkButton, *fCancelButton;
  //TGLayoutHints *fL1, *fL2, *fL3, *fL4, *fL5, *fL6;

  //TGLabel *fCutVar;
  //TGLabel *fCutMin;
  //TGLabel *fCutMax;

  //TGLabel *fHistVar;
  //TGLabel *fHistBin;
  //TGLabel *fHistMin;
  //TGLabel *fHistMax;

  //vector<TGLabel*> fCutLabel, fHistLabel;
  //vector<TGTextEntry*> fCutEntryMin,fHistEntryMin;
  //vector<TGTextEntry*> fCutEntryMax,fHistEntryMax;
  //vector<TGTextEntry*> fCutEntryBin,fHistEntryBin;
  //TGCompositeFrame *fHistFrame,*fCutFrame;

  PlaneECAL1  *fPlaneECAL1;

  TGSlider *fChanSlider;
  TGCompositeFrame *fSliderFrame;
  TGTextEntry *fSliderText;
  TRootEmbeddedCanvas *fChanDisplay;

 public:
  PlaneECAL1Panel(const TGWindow *p,const TGWindow *main,
		 UInt_t w, UInt_t h, PlaneECAL1 *plane);
  virtual ~PlaneECAL1Panel();

  /// Creates tab for scan of the channels time spectra
  void CreateControlTab();

  /// \return slider position
  int GetChannel() {return fChanSlider->GetPosition();}

  // Slots
  /// fMain->CloseWindow() connected here
  //virtual void CloseWindow() {delete this;}

  /// all basic TGButton::Clicked() connected here
  //void HandleButtons(int id=-1);

  ///fSlider->PositionChanged(int) connected here
  void SliderMoved(int pos);

  /// Draw spectrum for a given channel
  void Draw();

  /// fSliderText->TextChanged(const char*) connected here
  void TextChanged(const char* text);

#if ROOT_VERSION_CODE >= ROOT_30200
  RQ_OBJECT("PlaneECAL1Panel")
#else
  RQ_OBJECT()
#endif
  ClassDef(PlaneECAL1Panel,1)
};


#endif
