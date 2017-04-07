#include "PlaneRPD_SADC_Panel.h"
#include "TPad.h"
#include "TCanvas.h"

ClassImp(PlaneRPD_SADC_Panel);

PlaneRPD_SADC_Panel:: PlaneRPD_SADC_Panel(const TGWindow *p,const TGWindow *main,
			      UInt_t w, UInt_t h, PlaneRPD_SADC *plane)
  : PlanePanel(p, main, w, h, plane), fPlaneRPD_SADC(plane)
{
  CreateControlTab();
}

PlaneRPD_SADC_Panel::~PlaneRPD_SADC_Panel()
{
  delete fSliderText;
  delete fSliderFrame;
  delete fChanSlider;
  delete fChanDisplay;
  delete fSliderText_2;
  delete fSliderFrame_2;
  delete fChanSlider_2;
  delete fChanDisplay_2;
  delete fSliderText_3;
  delete fSliderFrame_3;
  delete fChanSlider_3;
  delete fChanDisplay_3;
  delete fSliderText_4;
  delete fSliderFrame_4;
  delete fChanSlider_4;
  delete fChanDisplay_4;
}

void PlaneRPD_SADC_Panel::CreateControlTab()
{

  // Here create the tab for amplitude
  TGCompositeFrame *tf=fTab->AddTab("ControlAmplitude");
  tf->SetLayoutManager(new TGVerticalLayout(tf));

  fSliderFrame=new TGCompositeFrame(tf, 60, 20, kHorizontalFrame);
  
  // slider
  fChanSlider = new TGHSlider(fSliderFrame, 300, kSlider1 | kScaleBoth);
  fChanSlider->Connect("PositionChanged(int)","PlaneRPD_SADC_Panel",this,
		       "SliderMoved(int)");
  fChanSlider->Connect("Released()","PlaneRPD_SADC_Panel",this,
		       "DrawAmplitude()");

  fChanSlider->SetRange(0,fPlaneRPD_SADC->GetNchannels()-1);
  
  fSliderFrame->AddFrame(fChanSlider);
  
  // channel number display   
  fSliderText = new TGTextEntry(fSliderFrame, new TGTextBuffer(10));
  fSliderText->Connect("TextChanged(const char*)","PlaneRPD_SADC_Panel",this, 
		       "TextChanged(const char*)");
  fSliderText->Connect("ReturnPressed()","PlaneRPD_SADC_Panel",this, 
		       "DrawAmplitude()");
  fSliderFrame->AddFrame(fSliderText,fL4);
  fSliderText->Resize(50,fSliderText->GetDefaultHeight());

  tf->AddFrame(fSliderFrame,fL1);
  
  //canvas
  fChanDisplay = new TRootEmbeddedCanvas("fChanDisplay",tf, 200, 200);
  tf->AddFrame(fChanDisplay,fL3);

  fMain->MapSubwindows();  
  fMain->Resize(fMain->GetDefaultSize());
  fMain->Layout();  

  // And here create the tab for Adc
  TGCompositeFrame *tf_2=fTab->AddTab("ControlAdc");
  tf_2->SetLayoutManager(new TGVerticalLayout(tf_2));

  fSliderFrame_2=new TGCompositeFrame(tf_2, 60, 20, kHorizontalFrame);

  fChanSlider_2 = new TGHSlider(fSliderFrame_2, 300, kSlider1 | kScaleBoth);
  fChanSlider_2->Connect("PositionChanged(int)","PlaneRPD_SADC_Panel",this,
			 "SliderMoved_2(int)");
  fChanSlider_2->Connect("Released()","PlaneRPD_SADC_Panel",this,
			 "DrawAdc()");

  fChanSlider_2->SetRange(0,fPlaneRPD_SADC->GetNchannels()-1);
  
  fSliderFrame_2->AddFrame(fChanSlider_2);

  fSliderText_2 = new TGTextEntry(fSliderFrame_2, new TGTextBuffer(10));
  fSliderText_2->Connect("TextChanged(const char*)","PlaneRPD_SADC_Panel",this, 
			 "TextChanged_2(const char*)");
  fSliderText_2->Connect("ReturnPressed()","PlaneRPD_SADC_Panel",this, 
			 "DrawAdc()");
  fSliderFrame_2->AddFrame(fSliderText_2,fL4);
  fSliderText_2->Resize(50,fSliderText_2->GetDefaultHeight());

  tf_2->AddFrame(fSliderFrame_2,fL1);
  
  fChanDisplay_2 = new TRootEmbeddedCanvas("fChanDisplay_2",tf_2, 200, 200);
  tf_2->AddFrame(fChanDisplay_2,fL3);
  
  fMain->MapSubwindows();  
  fMain->Resize(fMain->GetDefaultSize());
  fMain->Layout();  


  // And here create the tab for Pedestals
  TGCompositeFrame *tf_3=fTab->AddTab("ControlPedestals");
  tf_3->SetLayoutManager(new TGVerticalLayout(tf_3));

  fSliderFrame_3=new TGCompositeFrame(tf_3, 60, 20, kHorizontalFrame);

  fChanSlider_3 = new TGHSlider(fSliderFrame_3, 300, kSlider1 | kScaleBoth);
  fChanSlider_3->Connect("PositionChanged(int)","PlaneRPD_SADC_Panel",this,
			 "SliderMoved_3(int)");
  fChanSlider_3->Connect("Released()","PlaneRPD_SADC_Panel",this,
			 "DrawPedestals()");

  fChanSlider_3->SetRange(0,fPlaneRPD_SADC->GetNchannels()-1);

  fSliderFrame_3->AddFrame(fChanSlider_3);

  fSliderText_3 = new TGTextEntry(fSliderFrame_3, new TGTextBuffer(10));
  fSliderText_3->Connect("TextChanged(const char*)","PlaneRPD_SADC_Panel",this, 
			 "TextChanged_3(const char*)");
  fSliderText_3->Connect("ReturnPressed()","PlaneRPD_SADC_Panel",this, 
			 "DrawPedestals()");
  fSliderFrame_3->AddFrame(fSliderText_3,fL4);
  fSliderText_3->Resize(50,fSliderText_3->GetDefaultHeight());

  tf_3->AddFrame(fSliderFrame_3,fL1);

  fChanDisplay_3 = new TRootEmbeddedCanvas("fChanDisplay_3",tf_3, 200, 200);
  tf_3->AddFrame(fChanDisplay_3,fL3);

  fMain->MapSubwindows();  
  fMain->Resize(fMain->GetDefaultSize());
  fMain->Layout();  

  // And here create the tab for Pedestals
  TGCompositeFrame *tf_4=fTab->AddTab("ControlBeginningOfSignal");
  tf_4->SetLayoutManager(new TGVerticalLayout(tf_4));

  fSliderFrame_4=new TGCompositeFrame(tf_4, 60, 20, kHorizontalFrame);

  fChanSlider_4 = new TGHSlider(fSliderFrame_4, 300, kSlider1 | kScaleBoth);
  fChanSlider_4->Connect("PositionChanged(int)","PlaneRPD_SADC_Panel",this,
			 "SliderMoved_4(int)");
  fChanSlider_4->Connect("Released()","PlaneRPD_SADC_Panel",this,
			 "DrawBeginningOfSignal()");

  fChanSlider_4->SetRange(0,fPlaneRPD_SADC->GetNchannels()-1);

  fSliderFrame_4->AddFrame(fChanSlider_4);

  fSliderText_4 = new TGTextEntry(fSliderFrame_4, new TGTextBuffer(10));
  fSliderText_4->Connect("TextChanged(const char*)","PlaneRPD_SADC_Panel",this, 
			 "TextChanged_4(const char*)");
  fSliderText_4->Connect("ReturnPressed()","PlaneRPD_SADC_Panel",this, 
			 "DrawBeginningOfSignal()");
  fSliderFrame_4->AddFrame(fSliderText_4,fL4);
  fSliderText_4->Resize(50,fSliderText_4->GetDefaultHeight());

  tf_4->AddFrame(fSliderFrame_4,fL1);

  fChanDisplay_4 = new TRootEmbeddedCanvas("fChanDisplay_4",tf_4, 200, 200);
  tf_4->AddFrame(fChanDisplay_4,fL3);

  fMain->MapSubwindows();  
  fMain->Resize(fMain->GetDefaultSize());
  fMain->Layout();  

}

void PlaneRPD_SADC_Panel::SliderMoved(int pos) {
  char buf[10];
  sprintf(buf, "%d", pos);
  fSliderText->GetBuffer()->Clear();
  fSliderText->GetBuffer()->AddText(0, buf);
  gClient->NeedRedraw(fSliderText);
  fPlaneRPD_SADC->SetChannel(pos);
}

void PlaneRPD_SADC_Panel::SliderMoved_2(int pos) {
  char buf[10];
  sprintf(buf, "%d", pos);
  fSliderText_2->GetBuffer()->Clear();
  fSliderText_2->GetBuffer()->AddText(0, buf);
  gClient->NeedRedraw(fSliderText_2);
  fPlaneRPD_SADC->SetChannel_2(pos);
}

void PlaneRPD_SADC_Panel::SliderMoved_3(int pos)
{
  char buf[10];
  sprintf(buf, "%d", pos);
  fSliderText_3->GetBuffer()->Clear();
  fSliderText_3->GetBuffer()->AddText(0, buf);
  gClient->NeedRedraw(fSliderText_3);
  fPlaneRPD_SADC->SetChannel_3(pos);
}

void PlaneRPD_SADC_Panel::SliderMoved_4(int pos)
{
  char buf[10];
  sprintf(buf, "%d", pos);
  fSliderText_4->GetBuffer()->Clear();
  fSliderText_4->GetBuffer()->AddText(0, buf);
  gClient->NeedRedraw(fSliderText_4);
  fPlaneRPD_SADC->SetChannel_4(pos);
}

void PlaneRPD_SADC_Panel::DrawAmplitude() {
  fChanDisplay->GetCanvas()->cd(0);
  fPlaneRPD_SADC->MaxAmplitudeSpectrum();
  TPad::Pad()->Modified();
  TPad::Pad()->Update();
}

void PlaneRPD_SADC_Panel::DrawAdc() {
  fChanDisplay_2->GetCanvas()->cd(0);
  fPlaneRPD_SADC->AdcSpectrum();
  TPad::Pad()->Modified();
  TPad::Pad()->Update();
}

void PlaneRPD_SADC_Panel::DrawPedestals()
{
  fChanDisplay_3->GetCanvas()->cd(0);
  fPlaneRPD_SADC->PedestalsSpectrum();
  TPad::Pad()->Modified();
  TPad::Pad()->Update();
}

void PlaneRPD_SADC_Panel::DrawBeginningOfSignal()
{
  fChanDisplay_4->GetCanvas()->cd(0);
  fPlaneRPD_SADC->BeginningOfSignalSpectrum();
  TPad::Pad()->Modified();
  TPad::Pad()->Update();
}

void PlaneRPD_SADC_Panel::TextChanged(const char* text) {
  fChanSlider->SetPosition(atoi(fSliderText->GetBuffer()->GetString()));
  fPlaneRPD_SADC->SetChannel(atoi(fSliderText->GetBuffer()->GetString()));
}
 
void PlaneRPD_SADC_Panel::TextChanged_2(const char* text) {
  fChanSlider_2->SetPosition(atoi(fSliderText_2->GetBuffer()->GetString()));
  fPlaneRPD_SADC->SetChannel_2(atoi(fSliderText_2->GetBuffer()->GetString()));
}

void PlaneRPD_SADC_Panel::TextChanged_3(const char* text)
{
  fChanSlider_3->SetPosition(atoi(fSliderText_3->GetBuffer()->GetString()));
  fPlaneRPD_SADC->SetChannel_3(atoi(fSliderText_3->GetBuffer()->GetString()));
}

void PlaneRPD_SADC_Panel::TextChanged_4(const char* text)
{
  fChanSlider_4->SetPosition(atoi(fSliderText_4->GetBuffer()->GetString()));
  fPlaneRPD_SADC->SetChannel_4(atoi(fSliderText_4->GetBuffer()->GetString()));
}
