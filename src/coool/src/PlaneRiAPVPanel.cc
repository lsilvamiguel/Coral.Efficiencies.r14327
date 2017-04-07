#include "PlaneRiAPVPanel.h"
#include "TPad.h"

ClassImp(PlaneRiAPVPanel);


PlaneRiAPVPanel:: PlaneRiAPVPanel(const TGWindow *p,const TGWindow *main,
			      UInt_t w, UInt_t h, PlaneRiAPV *plane)
  : PlanePanel(p, main, w, h, plane), fPlaneRiAPV(plane) {

  //cout<<"this is a PlaneRiAPVPanel"<<endl;
  CreateControlTab();
}


PlaneRiAPVPanel::~PlaneRiAPVPanel() {

  delete fSliderText;
  delete fSliderFrame;
  delete fChanSlider;
  delete fChanDisplay;
}


void PlaneRiAPVPanel::CreateControlTab() {

  //cout<<"creating control tab"<<endl;

  //control tab
  TGCompositeFrame *tf=fTab->AddTab("Control");
  tf->SetLayoutManager(new TGVerticalLayout(tf));

  fSliderFrame=new TGCompositeFrame(tf, 60, 20, kHorizontalFrame);

  // slider
  fChanSlider = new TGHSlider(fSliderFrame, 300, kSlider1 | kScaleBoth);
  fChanSlider->Connect("PositionChanged(int)","PlaneRiAPVPanel",this,
		       "SliderMoved(int)");
  fChanSlider->Connect("Released()","PlaneRiAPVPanel",this,
		       "Draw()");

  fChanSlider->SetRange(0,fPlaneRiAPV->GetNchannels()-1);

  fSliderFrame->AddFrame(fChanSlider);

  // channel number display
  fSliderText = new TGTextEntry(fSliderFrame, new TGTextBuffer(10));
  fSliderText->Connect("TextChanged(const char*)","PlaneRiAPVPanel",this,
		       "TextChanged(const char*)");
  fSliderText->Connect("ReturnPressed()","PlaneRiAPVPanel",this,
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


void PlaneRiAPVPanel::SliderMoved(int pos) {
  char buf[10];
  sprintf(buf, "%d", pos);
  fSliderText->GetBuffer()->Clear();
  fSliderText->GetBuffer()->AddText(0, buf);
  gClient->NeedRedraw(fSliderText);
  fPlaneRiAPV->SetChannel(pos);
}


void PlaneRiAPVPanel::Draw() {
  fPlaneRiAPV->ChannelA2Spectrum();
  TPad::Pad()->Modified();
  TPad::Pad()->Update();
}


void PlaneRiAPVPanel::TextChanged(const char* text) {
  fChanSlider->SetPosition(atoi(fSliderText->GetBuffer()->GetString()));
  fPlaneRiAPV->SetChannel(atoi(fSliderText->GetBuffer()->GetString()));
}



