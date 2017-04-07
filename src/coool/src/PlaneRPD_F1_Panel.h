#ifndef __PlaneRPD_F1_Panel__
#define __PlaneRPD_F1_Panel__

#include "PlaneRPD_F1.h"
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

/// PlaneRPD_F1 control panel 

class PlaneRPD_F1_Panel : public PlanePanel {

 private:

  PlaneRPD_F1  *fPlaneRPD_F1;

  TGSlider *fChanSlider;
  TGCompositeFrame *fSliderFrame;
  TGTextEntry *fSliderText;
  TRootEmbeddedCanvas *fChanDisplay;

/*   TGSlider *fChanSlider_2; */
/*   TGCompositeFrame *fSliderFrame_2; */
/*   TGTextEntry *fSliderText_2; */
/*   TRootEmbeddedCanvas *fChanDisplay_2; */

 public:
  PlaneRPD_F1_Panel(const TGWindow *p,const TGWindow *main,
		 UInt_t w, UInt_t h, PlaneRPD_F1 *plane);
  virtual ~PlaneRPD_F1_Panel();
  
  /// Creates tab for scan of the channels time spectra
  virtual void CreateControlTab();

  /// \return slider position
  int GetChannel() {return fChanSlider->GetPosition();}

  /// \return slider position
/*   int GetChannel_2() {return fChanSlider_2->GetPosition();} */

  ///fSlider->PositionChanged(int) connected here
  void SliderMoved(int pos);                 

  ///fSlider_2->PositionChanged(int) connected here
/*   void SliderMoved_2(int pos); */

  /// Draw tdc spectrum for a given channel
  void DrawTdc();                         

  /// Draw multiplicity spectrum for a given channel
/*   void DrawMultiplicity(); */
  
  /// fSliderText->TextChanged(const char*) connected here    
  void TextChanged(const char* text);        

  /// fSliderText_2->TextChanged(const char*) connected here 
/*   void TextChanged_2(const char* text); */

#if ROOT_VERSION_CODE >= ROOT_30200
  RQ_OBJECT("PlaneRPD_F1_Panel")  
#else
  RQ_OBJECT()
#endif
  ClassDef(PlaneRPD_F1_Panel,1)
};


#endif
