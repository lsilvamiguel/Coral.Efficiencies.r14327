#include "Plane2VPanel.h"
#include "TPad.h"

ClassImp(Plane2VPanel);

Plane2VPanel:: Plane2VPanel(const TGWindow *p,const TGWindow *main,
			      UInt_t w, UInt_t h, Plane2V *plane)
  : PlanePanel(p, main, w, h, plane), fPlane2V(plane) {

  //cout<<"this is a Plane2VPanel"<<endl;
  CreateControlTab();
}

Plane2VPanel::~Plane2VPanel() {

  delete fSliderText;
  delete fSliderFrame;
  delete fChanSlider;
  delete fChanDisplay;
}

void Plane2VPanel::CreateControlTab() {


  //cout<<"creating control tab"<<endl;

  //control tab
  TGCompositeFrame *tf=fTab->AddTab("Control");
  tf->SetLayoutManager(new TGVerticalLayout(tf));

  fSliderFrame=new TGCompositeFrame(tf, 60, 20, kHorizontalFrame);
  
  // slider
  fChanSlider = new TGHSlider(fSliderFrame, 300, kSlider1 | kScaleBoth);
  fChanSlider->Connect("PositionChanged(int)","Plane2VPanel",this,
		       "SliderMoved(int)");
  fChanSlider->Connect("Released()","Plane2VPanel",this,
		       "Draw()");

  fChanSlider->SetRange(0,fPlane2V->GetNchannels()-1);
  
  fSliderFrame->AddFrame(fChanSlider);
  
  // channel number display   
  fSliderText = new TGTextEntry(fSliderFrame, new TGTextBuffer(10));
  fSliderText->Connect("TextChanged(const char*)","Plane2VPanel",this, 
		       "TextChanged(const char*)");
  fSliderText->Connect("ReturnPressed()","Plane2VPanel",this, 
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

void Plane2VPanel::SliderMoved(int pos) {
  char buf[10];
  sprintf(buf, "%d", pos);
  fSliderText->GetBuffer()->Clear();
  fSliderText->GetBuffer()->AddText(0, buf);
  gClient->NeedRedraw(fSliderText);
  fPlane2V->SetChannel(pos);
}

void Plane2VPanel::Draw() {
  fPlane2V->AmpSpectrum();
  TPad::Pad()->Modified();
  TPad::Pad()->Update();
}

void Plane2VPanel::TextChanged(const char* text) {
  fChanSlider->SetPosition(atoi(fSliderText->GetBuffer()->GetString()));
  fPlane2V->SetChannel(atoi(fSliderText->GetBuffer()->GetString()));
}



