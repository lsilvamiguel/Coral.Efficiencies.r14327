#ifndef __PlaneRPD_SADC_Panel__
#define __PlaneRPD_SADC_Panel__

#include "PlaneRPD_SADC.h"
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

/// PlaneRPD_SADC control panel 

class PlaneRPD_SADC_Panel : public PlanePanel
{
  
 private:
    
  PlaneRPD_SADC  *fPlaneRPD_SADC;
  
  TGSlider *fChanSlider;
  TGCompositeFrame *fSliderFrame;
  TGTextEntry *fSliderText;
  TRootEmbeddedCanvas *fChanDisplay;
  
  TGSlider *fChanSlider_2;
  TGCompositeFrame *fSliderFrame_2;
  TGTextEntry *fSliderText_2;
  TRootEmbeddedCanvas *fChanDisplay_2;
  
  TGSlider *fChanSlider_3;
  TGCompositeFrame *fSliderFrame_3;
  TGTextEntry *fSliderText_3;
  TRootEmbeddedCanvas *fChanDisplay_3;
  
  TGSlider *fChanSlider_4;
  TGCompositeFrame *fSliderFrame_4;
  TGTextEntry *fSliderText_4;
  TRootEmbeddedCanvas *fChanDisplay_4;
  
 public:
  PlaneRPD_SADC_Panel(const TGWindow *p,const TGWindow *main,
		      UInt_t w, UInt_t h, PlaneRPD_SADC *plane);
  virtual ~PlaneRPD_SADC_Panel();
  
  /// Creates tab for scan of the channels time spectra
  virtual void CreateControlTab();
    
    /// \return slider position
  int GetChannel() {return fChanSlider->GetPosition();}

  /// \return slider position
  int GetChannel_2() {return fChanSlider_2->GetPosition();}

  /// \return slider position
  int GetChannel_3() {return fChanSlider_3->GetPosition();}

  /// \return slider position
  int GetChannel_4() {return fChanSlider_4->GetPosition();}

  ///fSlider->PositionChanged(int) connected here
  void SliderMoved(int pos);                 

  ///fSlider->PositionChanged(int) connected here
  void SliderMoved_2(int pos);                 

  ///fSlider->PositionChanged(int) connected here
  void SliderMoved_3(int pos);

  ///fSlider->PositionChanged(int) connected here
  void SliderMoved_4(int pos);

  /// Draw the amplitude spectrum for a given channel
  void DrawAmplitude();                         

  /// Draw adc spectrum for a given channel
  void DrawAdc();                         

  /// Draw pedestals spectrum for a given channel
  void DrawPedestals();

  /// Draw beginning of signal spectrum for a given channel
  void DrawBeginningOfSignal();

  /// fSliderText->TextChanged(const char*) connected here    
  void TextChanged(const char* text);        

  /// fSliderText->TextChanged(const char*) connected here    
  void TextChanged_2(const char* text);

  /// fSliderText->TextChanged(const char*) connected here    
  void TextChanged_3(const char* text);

  /// fSliderText->TextChanged(const char*) connected here    
  void TextChanged_4(const char* text);

#if ROOT_VERSION_CODE >= ROOT_30200
  RQ_OBJECT("PlaneRPD_SADC_Panel")  
#else
  RQ_OBJECT()
#endif
  ClassDef(PlaneRPD_SADC_Panel,1)
};


#endif
