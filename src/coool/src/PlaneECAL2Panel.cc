#include "PlaneECAL2Panel.h"
#include "TPad.h"

ClassImp(PlaneECAL2Panel);

PlaneECAL2Panel:: PlaneECAL2Panel(const TGWindow *p,const TGWindow *main,
			      UInt_t w, UInt_t h, PlaneECAL2 *plane)
  : PlanePanel(p, main, w, h, plane), fPlaneECAL2(plane) {

  //cout<<"this is a PlaneECAL2Panel"<<endl;
  CreateControlTab();
}


PlaneECAL2Panel::~PlaneECAL2Panel() {

  delete fSliderText;
  delete fSliderFrame;
  delete fChanSlider;
  delete fChanDisplay;
}


void PlaneECAL2Panel::CreateControlTab() {

  //cout<<"creating control tab"<<endl;

  //control tab
  TGCompositeFrame *tf=fTab->AddTab("Control");
  tf->SetLayoutManager(new TGVerticalLayout(tf));

  fSliderFrame=new TGCompositeFrame(tf, 60, 20, kHorizontalFrame);

  // slider
  fChanSlider = new TGHSlider(fSliderFrame, 300, kSlider1 | kScaleBoth);
  fChanSlider->Connect("PositionChanged(int)","PlaneECAL2Panel",this,
		       "SliderMoved(int)");
  fChanSlider->Connect("Released()","PlaneECAL2Panel",this,
		       "Draw()");

  fChanSlider->SetRange(0,fPlaneECAL2->GetNchannels()-1);

  fSliderFrame->AddFrame(fChanSlider);

  // channel number display
  fSliderText = new TGTextEntry(fSliderFrame, new TGTextBuffer(10));
  fSliderText->Connect("TextChanged(const char*)","PlaneECAL2Panel",this,
		       "TextChanged(const char*)");
  fSliderText->Connect("ReturnPressed()","PlaneECAL2Panel",this,
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


void PlaneECAL2Panel::SliderMoved(int pos) {
  char buf[10];
  sprintf(buf, "%d", pos);
  fSliderText->GetBuffer()->Clear();
  fSliderText->GetBuffer()->AddText(0, buf);
  gClient->NeedRedraw(fSliderText);
  fPlaneECAL2->SetChannel(pos);
}


void PlaneECAL2Panel::Draw() {
  fPlaneECAL2->AmpSpectrum();
  TPad::Pad()->Modified();
  TPad::Pad()->Update();
}


void PlaneECAL2Panel::TextChanged(const char* text) {
  fChanSlider->SetPosition(atoi(fSliderText->GetBuffer()->GetString()));
  fPlaneECAL2->SetChannel(atoi(fSliderText->GetBuffer()->GetString()));
}



