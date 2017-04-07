#include "PlaneRPD_F1_Panel.h"
#include "TPad.h"
#include "TCanvas.h"

ClassImp(PlaneRPD_F1_Panel);

PlaneRPD_F1_Panel:: PlaneRPD_F1_Panel(const TGWindow *p,const TGWindow *main,
			      UInt_t w, UInt_t h, PlaneRPD_F1 *plane)
  : PlanePanel(p, main, w, h, plane), fPlaneRPD_F1(plane)
{
  CreateControlTab();
}

PlaneRPD_F1_Panel::~PlaneRPD_F1_Panel() {

  delete fSliderText;
  delete fSliderFrame;
  delete fChanSlider;
  delete fChanDisplay;
//   delete fSliderText_2;
//   delete fSliderFrame_2;
//   delete fChanSlider_2;
//   delete fChanDisplay_2;
}

void PlaneRPD_F1_Panel::CreateControlTab()
{
  //control tab for tdc
  TGCompositeFrame *tf=fTab->AddTab("ControlTdc");
  tf->SetLayoutManager(new TGVerticalLayout(tf));

  fSliderFrame=new TGCompositeFrame(tf, 60, 20, kHorizontalFrame);
  
  // slider
  fChanSlider = new TGHSlider(fSliderFrame, 300, kSlider1 | kScaleBoth);
  fChanSlider->Connect("PositionChanged(int)","PlaneRPD_F1_Panel",this,
		       "SliderMoved(int)");
  fChanSlider->Connect("Released()","PlaneRPD_F1_Panel",this,
		       "DrawTdc()");

  fChanSlider->SetRange(0,fPlaneRPD_F1->GetNchannels()-1);
  
  fSliderFrame->AddFrame(fChanSlider);
  
  // channel number display   
  fSliderText = new TGTextEntry(fSliderFrame, new TGTextBuffer(10));
  fSliderText->Connect("TextChanged(const char*)","PlaneRPD_F1_Panel",this, 
		       "TextChanged(const char*)");
  fSliderText->Connect("ReturnPressed()","PlaneRPD_F1_Panel",this, 
		       "DrawTdc()");
  fSliderFrame->AddFrame(fSliderText,fL4);
  fSliderText->Resize(50,fSliderText->GetDefaultHeight());

  tf->AddFrame(fSliderFrame,fL1);
  
  //canvas
  fChanDisplay = new TRootEmbeddedCanvas("fChanDisplay",tf, 200, 200);
  tf->AddFrame(fChanDisplay,fL3);
  
  fMain->MapSubwindows();  
  fMain->Resize(fMain->GetDefaultSize());
  fMain->Layout();  

  //control tab for multiplicity
//   TGCompositeFrame *tf_2=fTab->AddTab("ControlMultiplicity");
//   tf_2->SetLayoutManager(new TGVerticalLayout(tf_2));

//   fSliderFrame_2=new TGCompositeFrame(tf_2, 60, 20, kHorizontalFrame);

//   // slider
//   fChanSlider_2 = new TGHSlider(fSliderFrame_2, 300, kSlider1 | kScaleBoth);
//   fChanSlider_2->Connect("PositionChanged(int)","PlaneRPD_F1_Panel",this,
// 			 "SliderMoved_2(int)");
//   fChanSlider_2->Connect("Released()","PlaneRPD_F1_Panel",this,
// 			 "DrawMultiplicity()");

//   fChanSlider_2->SetRange(0,fPlaneRPD_F1->GetNchannels()-1);

//   fSliderFrame_2->AddFrame(fChanSlider_2);

//   // channel number display   
//   fSliderText_2 = new TGTextEntry(fSliderFrame_2, new TGTextBuffer(10));
//   fSliderText_2->Connect("TextChanged(const char*)","PlaneRPD_F1_Panel",this, 
// 			 "TextChanged_2(const char*)");
//   fSliderText_2->Connect("ReturnPressed()","PlaneRPD_F1_Panel",this, 
// 			 "DrawMultiplicity()");
//   fSliderFrame_2->AddFrame(fSliderText_2,fL4);
//   fSliderText_2->Resize(50,fSliderText_2->GetDefaultHeight());

//   tf_2->AddFrame(fSliderFrame_2,fL1);

//   //canvas
//   fChanDisplay_2 = new TRootEmbeddedCanvas("fChanDisplay_2",tf_2, 200, 200);
//   tf_2->AddFrame(fChanDisplay_2,fL3);

//   fMain->MapSubwindows();  
//   fMain->Resize(fMain->GetDefaultSize());
//   fMain->Layout();  
}

void PlaneRPD_F1_Panel::SliderMoved(int pos)
{
  char buf[10];
  sprintf(buf, "%d", pos);
  fSliderText->GetBuffer()->Clear();
  fSliderText->GetBuffer()->AddText(0, buf);
  gClient->NeedRedraw(fSliderText);
//   std::cout<<pos<<"\n";
  fPlaneRPD_F1->SetChannel(pos);
}

// void PlaneRPD_F1_Panel::SliderMoved_2(int pos)
// {
//   char buf[10];
//   sprintf(buf, "%d", pos);
//   fSliderText_2->GetBuffer()->Clear();
//   fSliderText_2->GetBuffer()->AddText(0, buf);
//   gClient->NeedRedraw(fSliderText_2);
//   fPlaneRPD_F1->SetChannel_2(pos);
// }

void PlaneRPD_F1_Panel::DrawTdc()
{
  fChanDisplay->GetCanvas()->cd(0);
  fPlaneRPD_F1->TdcSpectrum();
  TPad::Pad()->Modified();
  TPad::Pad()->Update();
}

// void PlaneRPD_F1_Panel::DrawMultiplicity()
// {
//   fChanDisplay_2->GetCanvas()->cd(0);
//   fPlaneRPD_F1->MultiplicitySpectrum();
//   TPad::Pad()->Modified();
//   TPad::Pad()->Update();
// }

void PlaneRPD_F1_Panel::TextChanged(const char* text) {
  fChanSlider->SetPosition(atoi(fSliderText->GetBuffer()->GetString()));
  fPlaneRPD_F1->SetChannel(atoi(fSliderText->GetBuffer()->GetString()));
}

// void PlaneRPD_F1_Panel::TextChanged_2(const char* text) {
//   fChanSlider_2->SetPosition(atoi(fSliderText_2->GetBuffer()->GetString()));
//   fPlaneRPD_F1->SetChannel_2(atoi(fSliderText_2->GetBuffer()->GetString()));
// }
