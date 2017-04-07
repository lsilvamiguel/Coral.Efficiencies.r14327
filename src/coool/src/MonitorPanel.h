#ifndef __MonitorPanel__
#define __MonitorPanel__

#include "Monitor.h"

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
//#include <TGComboBox.h>

#include "BitButton.h"

///  class for Monitor panel
class MonitorPanel : public TObject {

 private:
  unsigned int fTMPattern;
  bool         fTMLowerOnly;
  bool         fTMStrict;

  std::set<int> fEvtTypesChecked;

  /// Associated Monitor
  Monitor *fMonitor; 

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

  //TGComboBox           *fTTBox;
  std::vector<TGCheckButton*>  fEvtTypeChecks;
  TGTextEntry          *fHexaMask;
  TGCompositeFrame     *fTrigMask1;
  TGGroupFrame         *fTrigMask;  
  TGGroupFrame         *fTrigType;
  std::vector<TGGroupFrame*> fByte;
  std::vector<BitButton*>    fBitButton;
  TGCheckButton*        fCBLower;
  TGCheckButton*        fCBStrict;

 public:

  MonitorPanel(const TGWindow *p,const TGWindow *main,
	       UInt_t w, UInt_t h, Monitor *monitor, unsigned int trigmask,
           bool lowerOnly, bool strict);
  ~MonitorPanel();
  
  // Slots  
  /// fMain->CloseWindow() connected here                             
  void CloseWindow() {fMonitor->PanelClosed(); delete this;}
  
  /// all basic TGButton::Clicked() connected here
  void HandleButtons(int id=-1);             
  
  /// fHexaMask ReturnPressed() connected here
  bool ConvertHexa2Bin();

  /// fBitButton Clicked() connected here
  bool ConvertBin2Hexa();
 
  /// all evt type check buttons connected here
  void EvtTypeChecked(int id=-1);
  
  /// same as clicking on the button
  void CheckEvtType(int id);

  /// all trigger mask check buttons connected here
  void TrigMaskChecked(int id=-1);

#if ROOT_VERSION_CODE >= ROOT_30200
  RQ_OBJECT("MonitorPanel")
#else
  RQ_OBJECT()
#endif
  ClassDef(MonitorPanel,0)
};

#endif





