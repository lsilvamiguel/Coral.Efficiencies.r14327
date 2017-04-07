#include "PlanePropFrame.h"


///////////////////////       PlanePropFrame /////////////////////////////

PlanePropFrame::PlanePropFrame(const TGWindow *p, const TGWindow *main, 
			       UInt_t w, UInt_t h, Plane *plane) 
  : TGTransientFrame(p,main,w,h),
    fMain(main),
    fPlane(plane) { 

  

  fFrame1 = new TGHorizontalFrame(this, 60, 20, kFixedWidth);
  
  fOkButton = new TGTextButton(fFrame1, "&Ok", M_PROPERTIES_OK);
  fOkButton->Associate(this);
  fCancelButton = new TGTextButton(fFrame1, "&Cancel", M_PROPERTIES_CANCEL);
  fCancelButton->Associate(this);
  
  fL1 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,
			  2, 2, 2, 2);
  fL2 = new TGLayoutHints(kLHintsBottom | kLHintsRight, 2, 2, 5, 1);
  fL3 = new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
					 kLHintsExpandY, 2, 2, 5, 1);  

  fL4 = new TGLayoutHints(kLHintsRight, 2, 2, 5, 1);
  fL5 = new TGLayoutHints(kLHintsLeft, 2, 2, 5, 1);

  fFrame1->AddFrame(fOkButton, fL1);
  fFrame1->AddFrame(fCancelButton, fL1);
  
  fFrame1->Resize(150, fOkButton->GetDefaultHeight());
  AddFrame(fFrame1, fL2);
  
  //--------- create Tab widget and some composite frames for Tab testing
  fTab = new TGTab(this, 300, 300);  

  //Tab for Cuts 
  TGCompositeFrame *tf=fTab->AddTab("Variables");
  tf->SetLayoutManager(new TGMatrixLayout(tf,0,4));
  vector<TH1*>& vars=fPlane->GetVariables(); 

  fVar=new TGLabel(tf,"");
  fBin=new TGLabel(tf,"Number of bins");
  fMin=new TGLabel(tf,"Minimum");
  fMax=new TGLabel(tf,"Maximum");
  tf->AddFrame(fVar);
  tf->AddFrame(fBin);
  tf->AddFrame(fMin);
  tf->AddFrame(fMax);
 
  for(unsigned int i=0;i<vars.size();i++) {

    char bin[20];
    char min[20];
    char max[20];

    sprintf(bin,"%d",vars[i]->GetNbinsX());
    sprintf(min,"%d",(int) vars[i]->GetXaxis()->GetXmin());
    sprintf(max,"%d",(int) vars[i]->GetXaxis()->GetXmax());


    TGLabel *chb=new TGLabel(tf,vars[i]->GetName());
    fCutLabel.push_back(chb);
    tf->AddFrame(chb);
    //chb->Resize(100,chb->GetDefaultHeight());
    
    TGTextEntry *tebin=new TGTextEntry(tf,new TGTextBuffer(1034));
    fCutEntryBin.push_back(tebin);
    tf->AddFrame(tebin);
    tebin->Resize(100,tebin->GetDefaultHeight());
    tebin->SetText(bin);
   
    TGTextEntry *temin=new TGTextEntry(tf,new TGTextBuffer(1034));
    fCutEntryMin.push_back(temin);
    tf->AddFrame(temin);
    temin->Resize(100,temin->GetDefaultHeight());
    temin->SetText(min);

    TGTextEntry *temax=new TGTextEntry(tf,new TGTextBuffer(1034));
    fCutEntryMax.push_back(temax);
    tf->AddFrame(temax);
    temax->Resize(100,temax->GetDefaultHeight());    
    temax->SetText(max);
  }

  //Tab for general options
  tf=fTab->AddTab("General");
  


  //the end
  TGLayoutHints *fL3 = new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
					 kLHintsExpandY, 2, 2, 5, 1);
  AddFrame(fTab, fL3);

  // position relative to the parent's window
  Window_t wdum;
  int ax, ay;
  gVirtualX->TranslateCoordinates(main->GetId(), GetParent()->GetId(),
			   (((TGFrame *) main)->GetWidth() - fWidth) >> 1,
			   (((TGFrame *) main)->GetHeight() - fHeight) >> 1,
                           ax, ay, wdum);
  Move(ax, ay);

  this->SetWindowName((fPlane->GetName()).c_str());
  this->MapSubwindows();
  this->Resize(this->GetDefaultSize());
  this->MapWindow();    
}

PlanePropFrame::~PlanePropFrame() {


  delete fVar;
  delete fBin;
  delete fMin;
  delete fMax;

  fCutLabel.clear();
  fCutEntryMin.clear();
  fCutEntryMax.clear();
  fCutEntryBin.clear();
  delete fTab;
  delete fL1; delete fL2; delete fL3; delete fL4; delete fL5;
  delete fCancelButton; delete fOkButton;
  delete fFrame1;
  
}

Bool_t PlanePropFrame::ProcessMessage(Long_t msg, Long_t parm1, Long_t) {

  switch (GET_MSG(msg)) {
  
  case kC_COMMAND:

    switch (GET_SUBMSG(msg)) {
    case kCM_BUTTON:

      switch(parm1) {
      case M_PROPERTIES_OK:
	{
	  if(fPlane->IsInTree()) {
	    new TGMsgBox(gClient->GetRoot(),this,
			 "message",
			 "This detector being in the tree, the whole tree will be reset. Continue ?",
			 kMBIconExclamation,
			 kMBCancel |  kMBOk);
	  }
	  vector<TH1*>& vars=fPlane->GetVariables(); 
	  for(unsigned int i=0;i<vars.size();i++) {
	    int bin=atoi(fCutEntryBin[i]->GetText());
	    int min=atoi(fCutEntryMin[i]->GetText());
	    int max=atoi(fCutEntryMax[i]->GetText());

	    vars[i]->SetBins(bin,min,max);
	    vars[i]->Reset();
	  }
	  break;
	}
      case M_PROPERTIES_CANCEL:
	CloseWindow();
	break;
      default:
	break;
      }

      break;
    default:
      break;
    }
  
  default:
    break;
  }
  return true;
}

void PlanePropFrame::CloseWindow() {

  delete this;
}




