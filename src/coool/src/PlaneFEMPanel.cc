#include "PlaneFEMPanel.h"
#include "TPad.h"

ClassImp(PlaneFEMPanel);

PlaneFEMPanel:: PlaneFEMPanel(const TGWindow *p,const TGWindow *main,
			      UInt_t w, UInt_t h, PlaneFEM *plane)
  : PlanePanel(p, main, w, h, plane), fPlaneFEM(plane) {

  //cout<<"this is a PlaneFEMPanel"<<endl;
  CreateControlTab();
}

PlaneFEMPanel::~PlaneFEMPanel() {

  delete fSliderText;
  delete fSliderFrame;
  delete fChanSlider;
  delete fChanDisplay;
}

void PlaneFEMPanel::CreateControlTab() {


  //cout<<"creating control tab"<<endl;

  //control tab
  TGCompositeFrame *tf=fTab->AddTab("Control");
  tf->SetLayoutManager(new TGVerticalLayout(tf));

  fSliderFrame=new TGCompositeFrame(tf, 60, 20, kHorizontalFrame);

  // slider
  fChanSlider = new TGHSlider(fSliderFrame, 300, kSlider1 | kScaleBoth);
  fChanSlider->Connect("PositionChanged(int)","PlaneFEMPanel",this,
		       "SliderMoved(int)");
  fChanSlider->Connect("Released()","PlaneFEMPanel",this,
		       "Draw()");

  fChanSlider->SetRange(0,fPlaneFEM->GetNchannels()-1);

  fSliderFrame->AddFrame(fChanSlider);

  // channel number display
  fSliderText = new TGTextEntry(fSliderFrame, new TGTextBuffer(10));
  fSliderText->Connect("TextChanged(const char*)","PlaneFEMPanel",this,
		       "TextChanged(const char*)");
  fSliderText->Connect("ReturnPressed()","PlaneFEMPanel",this,
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

void PlaneFEMPanel::SliderMoved(int pos) {
  char buf[10];
  sprintf(buf, "%d", pos);
  fSliderText->GetBuffer()->Clear();
  fSliderText->GetBuffer()->AddText(0, buf);
  gClient->NeedRedraw(fSliderText);
  fPlaneFEM->SetChannel(pos);
}

void PlaneFEMPanel::Draw() {
  fPlaneFEM->AmpSpectrum();
  TPad::Pad()->Modified();
  TPad::Pad()->Update();
}

void PlaneFEMPanel::TextChanged(const char* text) {
  fChanSlider->SetPosition(atoi(fSliderText->GetBuffer()->GetString()));
  fPlaneFEM->SetChannel(atoi(fSliderText->GetBuffer()->GetString()));
}



