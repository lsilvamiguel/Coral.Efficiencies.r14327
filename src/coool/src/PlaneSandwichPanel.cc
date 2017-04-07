#include "PlaneSandwichPanel.h"
#include "TPad.h"

ClassImp(PlaneSandwichPanel);

PlaneSandwichPanel:: PlaneSandwichPanel(const TGWindow *p,const TGWindow *main,
			      UInt_t w, UInt_t h, PlaneSandwich *plane)
  : PlanePanel(p, main, w, h, plane), fPlaneSandwich(plane) {

  //cout<<"this is a PlaneSandwichPanel"<<endl;
  CreateControlTab();
}

PlaneSandwichPanel::~PlaneSandwichPanel() {

  delete fSliderText;
  delete fSliderFrame;
  delete fChanSlider;
  delete fChanDisplay;
}

void PlaneSandwichPanel::CreateControlTab() {


  //cout<<"creating control tab"<<endl;

  //control tab
  TGCompositeFrame *tf=fTab->AddTab("Control");
  tf->SetLayoutManager(new TGVerticalLayout(tf));

  fSliderFrame=new TGCompositeFrame(tf, 60, 20, kHorizontalFrame);
  
  // slider
  fChanSlider = new TGHSlider(fSliderFrame, 300, kSlider1 | kScaleBoth);
  fChanSlider->Connect("PositionChanged(int)","PlaneSandwichPanel",this,
		       "SliderMoved(int)");
  fChanSlider->Connect("Released()","PlaneSandwichPanel",this,
		       "Draw()");

  fChanSlider->SetRange(0,fPlaneSandwich->GetNchannels()-1);
  
  fSliderFrame->AddFrame(fChanSlider);
  
  // channel number display   
  fSliderText = new TGTextEntry(fSliderFrame, new TGTextBuffer(10));
  fSliderText->Connect("TextChanged(const char*)","PlaneSandwichPanel",this, 
		       "TextChanged(const char*)");
  fSliderText->Connect("ReturnPressed()","PlaneSandwichPanel",this, 
		       "Draw()");
  fSliderFrame->AddFrame(fSliderText,fL4);
  fSliderText->Resize(50,fSliderText->GetDefaultHeight());

  tf->AddFrame(fSliderFrame,fL1);
  
  //canvas
  fChanDisplay = new TRootEmbeddedCanvas("fChanDisplay",tf, 200, 200);
  tf->AddFrame(fChanDisplay,fL3);
  
  fMain->MapSubwindows();  
  fMain->Resize(fMain->GetDefaultSize());
  fMain->Layout();  
}

void PlaneSandwichPanel::SliderMoved(int pos) {
  char buf[10];
  sprintf(buf, "%d", pos);
  fSliderText->GetBuffer()->Clear();
  fSliderText->GetBuffer()->AddText(0, buf);
  gClient->NeedRedraw(fSliderText);
  //fPlaneSandwich->SetChannel(pos);
}

void PlaneSandwichPanel::Draw() {
  //fPlaneSandwich->TimeSpectrum();
  TPad::Pad()->Modified();
  TPad::Pad()->Update();
}

void PlaneSandwichPanel::TextChanged(const char* text) {
  fChanSlider->SetPosition(atoi(fSliderText->GetBuffer()->GetString()));
  //fPlaneSandwich->SetChannel(atoi(fSliderText->GetBuffer()->GetString()));
}
