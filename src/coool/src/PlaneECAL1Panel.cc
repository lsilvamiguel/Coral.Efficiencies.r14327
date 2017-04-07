#include "PlaneECAL1Panel.h"
#include "TPad.h"

ClassImp(PlaneECAL1Panel);

PlaneECAL1Panel:: PlaneECAL1Panel(const TGWindow *p,const TGWindow *main,
			      UInt_t w, UInt_t h, PlaneECAL1 *plane)
  : PlanePanel(p, main, w, h, plane), fPlaneECAL1(plane) {

  //cout<<"this is a PlaneECAL1Panel"<<endl;
  CreateControlTab();
}

PlaneECAL1Panel::~PlaneECAL1Panel() {

  delete fSliderText;
  delete fSliderFrame;
  delete fChanSlider;
  delete fChanDisplay;
}

void PlaneECAL1Panel::CreateControlTab() {


  //cout<<"creating control tab"<<endl;

  //control tab
  TGCompositeFrame *tf=fTab->AddTab("Control");
  tf->SetLayoutManager(new TGVerticalLayout(tf));

  fSliderFrame=new TGCompositeFrame(tf, 60, 20, kHorizontalFrame);

  // slider
  fChanSlider = new TGHSlider(fSliderFrame, 300, kSlider1 | kScaleBoth);
  fChanSlider->Connect("PositionChanged(int)","PlaneECAL1Panel",this,
		       "SliderMoved(int)");
  fChanSlider->Connect("Released()","PlaneECAL1Panel",this,
		       "Draw()");

  fChanSlider->SetRange(0,fPlaneECAL1->GetNchannels()-1);

  fSliderFrame->AddFrame(fChanSlider);

  // channel number display
  fSliderText = new TGTextEntry(fSliderFrame, new TGTextBuffer(10));
  fSliderText->Connect("TextChanged(const char*)","PlaneECAL1Panel",this,
		       "TextChanged(const char*)");
  fSliderText->Connect("ReturnPressed()","PlaneECAL1Panel",this,
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

void PlaneECAL1Panel::SliderMoved(int pos) {
  char buf[10];
  sprintf(buf, "%d", pos);
  fSliderText->GetBuffer()->Clear();
  fSliderText->GetBuffer()->AddText(0, buf);
  gClient->NeedRedraw(fSliderText);
  fPlaneECAL1->SetChannel(pos);
}

void PlaneECAL1Panel::Draw() {
  fPlaneECAL1->AmpSpectrum();
  TPad::Pad()->Modified();
  TPad::Pad()->Update();
}

void PlaneECAL1Panel::TextChanged(const char* text) {
  fChanSlider->SetPosition(atoi(fSliderText->GetBuffer()->GetString()));
  fPlaneECAL1->SetChannel(atoi(fSliderText->GetBuffer()->GetString()));
}



