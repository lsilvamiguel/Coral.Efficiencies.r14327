#include "Plane1VPanel.h"
#include "TPad.h"

ClassImp(Plane1VPanel);

Plane1VPanel:: Plane1VPanel(const TGWindow *p,const TGWindow *main,
			      UInt_t w, UInt_t h, Plane1V *plane)
  : PlanePanel(p, main, w, h, plane), fPlane1V(plane) {

  //cout<<"this is a Plane1VPanel"<<endl;
  CreateControlTab();
}

Plane1VPanel::~Plane1VPanel() {

  delete fSliderText;
  delete fSliderFrame;
  delete fChanSlider;
  delete fChanDisplay;
}

void Plane1VPanel::CreateControlTab() {


  //cout<<"creating control tab"<<endl;

  //control tab
  TGCompositeFrame *tf=fTab->AddTab("Control");
  tf->SetLayoutManager(new TGVerticalLayout(tf));

  fSliderFrame=new TGCompositeFrame(tf, 60, 20, kHorizontalFrame);
  
  // slider
  fChanSlider = new TGHSlider(fSliderFrame, 300, kSlider1 | kScaleBoth);
  fChanSlider->Connect("PositionChanged(int)","Plane1VPanel",this,
		       "SliderMoved(int)");
  fChanSlider->Connect("Released()","Plane1VPanel",this,
		       "Draw()");

  fChanSlider->SetRange(0,fPlane1V->GetNchannels()-1);
  
  fSliderFrame->AddFrame(fChanSlider);
  
  // channel number display   
  fSliderText = new TGTextEntry(fSliderFrame, new TGTextBuffer(10));
  fSliderText->Connect("TextChanged(const char*)","Plane1VPanel",this, 
		       "TextChanged(const char*)");
  fSliderText->Connect("ReturnPressed()","Plane1VPanel",this, 
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

void Plane1VPanel::SliderMoved(int pos) {
  char buf[10];
  sprintf(buf, "%d", pos);
  fSliderText->GetBuffer()->Clear();
  fSliderText->GetBuffer()->AddText(0, buf);
  gClient->NeedRedraw(fSliderText);
  fPlane1V->SetChannel(pos);
}

void Plane1VPanel::Draw() {
  fPlane1V->TimeSpectrum();
  TPad::Pad()->Modified();
  TPad::Pad()->Update();
}

void Plane1VPanel::TextChanged(const char* text) {
  fChanSlider->SetPosition(atoi(fSliderText->GetBuffer()->GetString()));
  fPlane1V->SetChannel(atoi(fSliderText->GetBuffer()->GetString()));
}



