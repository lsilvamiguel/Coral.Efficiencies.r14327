
#include <iostream>
#include "PlanePanel.h"

ClassImp(PlanePanel);


enum  EPropCommandIdentifiers {
  M_PROPERTIES_CANCEL,
  M_HIST_OK,
  M_CUT_OK,
  M_TRIG_OK,
  M_TT_COMBO
};


PlanePanel:: PlanePanel(const TGWindow *p,const TGWindow *main,
			      UInt_t w, UInt_t h, Plane *plane) :
  TObject(),fPlane(plane) { 

  fMain = new TGTransientFrame(p,main,w,h);
  fMain->Connect("CloseWindow()", "PlanePanel", this, 
  		 "CloseWindow()");
   
  fL1 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,
			  2, 2, 2, 2);
  fL2 = new TGLayoutHints(kLHintsBottom | kLHintsRight, 2, 4, 5, 2);
  fL3 = new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
					 kLHintsExpandY, 2, 2, 5, 1);  

  fL4 = new TGLayoutHints(kLHintsRight, 2, 2, 5, 1);
  fL5 = new TGLayoutHints(kLHintsLeft, 2, 2, 5, 1);
  fL6 = new TGLayoutHints(kLHintsCenterY | kLHintsCenterX, 2, 2, 5, 1);
  fL7 = new TGLayoutHints(kLHintsBottom | kLHintsRight, 0, 0, 0, 0);
  fL8 = new TGLayoutHints(kLHintsTop | kLHintsExpandX |
					 kLHintsExpandY, 2, 2, 5, 1);  
  fL9 = new TGLayoutHints(kLHintsBottom | kLHintsRight, 2, 2, 5, 1);
 
  fFrame1 = new TGHorizontalFrame(fMain, 60, 20, kFixedWidth);
  fCancelButton = new TGTextButton(fFrame1, "&Close", M_PROPERTIES_CANCEL);
  fCancelButton->Connect("Clicked()","PlanePanel",this,"HandleButtons()");
  fFrame1->AddFrame(fCancelButton, fL1);  
  fFrame1->Resize(150, fCancelButton->GetDefaultHeight());
  fMain->AddFrame(fFrame1, fL2);
    

  //--------- create Tab widget and some composite frames for Tab testing
  fTab = new TGTab(fMain, 300, 300);  

  //Tab for histograms 
  TGCompositeFrame *tf=fTab->AddTab("Histograms");
  tf->SetLayoutManager(new TGVerticalLayout(tf));
  fHistCanvas = new TGCanvas(tf, 10, 10, kSunkenFrame | kDoubleBorder);

  fHistFrame=new TGCompositeFrame(fHistCanvas->GetViewPort(), 60, 20);
  fHistCanvas->SetContainer(fHistFrame);

  fHistFrame->SetLayoutManager(new TGMatrixLayout(fHistFrame,0,4));
  std::vector<TH1*>& hists=fPlane->GetHistoList(); 

  fHistVar=new TGLabel(fHistFrame,"");
  fHistBin=new TGLabel(fHistFrame,"Number of bins");
  fHistMin=new TGLabel(fHistFrame,"Minimum");
  fHistMax=new TGLabel(fHistFrame,"Maximum");
  fHistFrame->AddFrame(fHistVar);
  fHistFrame->AddFrame(fHistBin);
  fHistFrame->AddFrame(fHistMin);
  fHistFrame->AddFrame(fHistMax);
  
  for(unsigned int i=0;i<hists.size();i++) {

    for(int j=0; j<hists[i]->GetDimension();j++) {
      
      TAxis *axis=0;
      std::string suffix;
      if(j==0) {axis=hists[i]->GetXaxis(); suffix="_x"; }
      if(j==1) {axis=hists[i]->GetYaxis(); suffix="_y"; }
      if(j==2) {axis=hists[i]->GetZaxis(); suffix="_z"; }

      char bin[20];
      char min[20];
      char max[20];

      sprintf(bin,"%d",axis->GetNbins());
      sprintf(min,"%d",(int)axis->GetXmin());
      sprintf(max,"%d",(int)axis->GetXmax());
      
      std::string name=hists[i]->GetName() + suffix;
      
      TGLabel *chb=new TGLabel(fHistFrame,name.c_str());
      fHistLabel.push_back(chb);
      fHistFrame->AddFrame(chb);
      //chb->Resize(100,chb->GetDefaultHeight());
      
      TGTextEntry *tebin=new TGTextEntry(fHistFrame,new TGTextBuffer(1034));
      fHistEntryBin.push_back(tebin);
      fHistFrame->AddFrame(tebin);
      tebin->Resize(100,tebin->GetDefaultHeight());
      tebin->SetText(bin);
   
      TGTextEntry *temin=new TGTextEntry(fHistFrame,new TGTextBuffer(1034));
      fHistEntryMin.push_back(temin);
      fHistFrame->AddFrame(temin);
      temin->Resize(100,temin->GetDefaultHeight());
      temin->SetText(min);

      TGTextEntry *temax=new TGTextEntry(fHistFrame,new TGTextBuffer(1034));
      fHistEntryMax.push_back(temax);
      fHistFrame->AddFrame(temax);
      temax->Resize(100,temax->GetDefaultHeight());    
      temax->SetText(max);
    }
  }
  //  fHistCanvas->AddFrame(fHistFrame);
  tf->AddFrame(fHistCanvas,fL8);
  fHistCanvas->Resize(300,300);

  fFrame2 = new TGHorizontalFrame(tf, 60, 20, kFixedWidth);
  fHistOkButton = new TGTextButton(fFrame2, "   &Ok  ", M_HIST_OK);
  fHistOkButton->Connect("Clicked()","PlanePanel",this,"HandleButtons()");
  fFrame2->AddFrame(fHistOkButton, fL1);  
  tf->AddFrame(fFrame2, fL9);
  fFrame2->Resize(150, fHistOkButton->GetDefaultHeight());

  // tab for cuts
  tf=fTab->AddTab("Cuts");
  tf->SetLayoutManager(new TGVerticalLayout(tf));
  fCutCanvas = new TGCanvas(tf, 10, 10, kSunkenFrame | kDoubleBorder);

  fCutFrame=new TGCompositeFrame(fCutCanvas->GetViewPort(), 60, 20);
  fCutFrame->SetLayoutManager(new TGMatrixLayout(fCutFrame,0,3));
  fCutCanvas->SetContainer(fCutFrame);

  fCutVar=new TGLabel(fCutFrame,"");
  fCutMin=new TGLabel(fCutFrame,"Minimum");
  fCutMax=new TGLabel(fCutFrame,"Maximum");

  fCutFrame->AddFrame(fCutVar);
  fCutFrame->AddFrame(fCutMin);
  fCutFrame->AddFrame(fCutMax);
 
  std::vector<Variable*>& vars=fPlane->GetVariables();
  for(unsigned int i=0;i<vars.size();i++) {

    char min[20];
    char max[20];

    sprintf(min,"%d",(int) vars[i]->GetMin());
    sprintf(max,"%d",(int) vars[i]->GetMax());


    TGLabel *chb=new TGLabel(fCutFrame,(vars[i]->GetName()).c_str());
    fCutLabel.push_back(chb);
    fCutFrame->AddFrame(chb);
    //chb->Resize(100,chb->GetDefaultHeight());
   
    TGTextEntry *temin=new TGTextEntry(fCutFrame,new TGTextBuffer(1034));
    fCutEntryMin.push_back(temin);
    fCutFrame->AddFrame(temin);
    temin->Resize(100,temin->GetDefaultHeight());
    temin->SetText(min);

    TGTextEntry *temax=new TGTextEntry(fCutFrame,new TGTextBuffer(1034));
    fCutEntryMax.push_back(temax);
    fCutFrame->AddFrame(temax);
    temax->Resize(100,temax->GetDefaultHeight());    
    temax->SetText(max);
  }
  tf->AddFrame(fCutCanvas,fL8);
  fCutCanvas->Resize(300,300);

  fFrame3 = new TGHorizontalFrame(tf, 60, 20, kFixedWidth);
  fCutOkButton = new TGTextButton(fFrame3, "   &Ok  ", M_CUT_OK);
  fCutOkButton->Connect("Clicked()","PlanePanel",this,"HandleButtons()");
  fFrame3->AddFrame(fCutOkButton, fL1);  
  tf->AddFrame(fFrame3, fL9);
  fFrame3->Resize(150, fCutOkButton->GetDefaultHeight());

  // tab for trigger
  tf=fTab->AddTab("Trigger");
  tf->SetLayoutManager(new TGVerticalLayout(tf));

  // trigger type
  fTrigType = new TGGroupFrame(tf, "Trigger Type", kVerticalFrame);
  fTTBox = new TGComboBox(fTrigType, M_TT_COMBO);
  fTTBox->Connect("Selected(int)","PlanePanel",this,"TTSelected(int)");
    
  std::string prefix = "Trigger type ";
  char suffix[2];
  for (register int i=0; i<4; i++) {
    sprintf(suffix,"%d",i+1);
    std::string tt = prefix + suffix; 
    fTTBox->AddEntry(tt.c_str(),i);
  }
  fTTBox->Select(0);
  fTrigType->AddFrame(fTTBox,fL6);
  fTTBox->Resize(150 ,30);
  
  // trigger mask
  fTrigMask = new TGGroupFrame(tf, "Trigger Mask", kVerticalFrame);
  fTrigMask1 = new TGCompositeFrame(fTrigMask,60,20);
  fTrigMask1->SetLayoutManager(new TGHorizontalLayout(fTrigMask1));

  fByte.push_back(new TGGroupFrame(fTrigMask1, "8", kHorizontalFrame));
  fByte.push_back(new TGGroupFrame(fTrigMask1, "16", kHorizontalFrame));
  fByte.push_back(new TGGroupFrame(fTrigMask1, "24", kHorizontalFrame));
  fByte.push_back(new TGGroupFrame(fTrigMask1, "32", kHorizontalFrame));

  for(int i=0; i<32; i++) {
    int bytei=i/8;
    //cout<<bytei<<endl;
    BitButton *button=new BitButton(fByte[bytei],i);
    button->Connect("Clicked()","PlanePanel",this,"ConvertBin2Hexa()");
    fBitButton.push_back(button);
    fByte[bytei]->AddFrame(button,fL7);
  }
  for(int i=3; i>=0; i--) {
    fTrigMask1->AddFrame(fByte[i]);
  }

  fHexaMask = new TGTextEntry(fTrigMask,"0x00000000");
  fHexaMask->Connect("ReturnPressed()","PlanePanel",this, 
		       "ConvertHexa2Bin()");
  fTrigMask->AddFrame(fTrigMask1);
  fTrigMask->AddFrame(fHexaMask,fL6);

  tf->AddFrame(fTrigType,fL1);
  tf->AddFrame(fTrigMask,fL1);

  fFrame4 = new TGHorizontalFrame(tf, 60, 20, kFixedWidth);
  fTrigOkButton = new TGTextButton(fFrame4, "   &Ok  ", M_TRIG_OK);
  fTrigOkButton->Connect("Clicked()","PlanePanel",this,"HandleButtons()");
  fFrame4->AddFrame(fTrigOkButton, fL1);  
  tf->AddFrame(fFrame4, fL9);
  fFrame4->Resize(150, fCutOkButton->GetDefaultHeight()); 

  //the end
  TGLayoutHints *fL3 = new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
					 kLHintsExpandY, 2, 2, 5, 1);
  fMain->AddFrame(fTab, fL3);

  // position relative to the parent's window
  Window_t wdum;
  int ax, ay;
  gVirtualX->TranslateCoordinates(fMain->GetId(), fMain->GetParent()->GetId(),
				  (((TGFrame *) fMain)->GetWidth() - 
				   fMain->GetWidth()) >> 1,
				  (((TGFrame *) fMain)->GetHeight() - 
				   fMain->GetHeight()) >> 1,
				  ax, ay, wdum);
  fMain->Move(ax, ay);

  fMain->SetWindowName(fPlane->GetName());
  fMain->MapSubwindows();
  fMain->Resize(fMain->GetDefaultSize());
  fMain->MapWindow();     

}

PlanePanel::~PlanePanel() {

  delete fHexaMask;
  fByte.clear();
  fBitButton.clear();
  delete fTrigMask1;
  delete fTrigMask;
  delete fTrigType;

  delete fCutVar;
  delete fCutMin;
  delete fCutMax;
  delete fHistVar;
  delete fHistBin;
  delete fHistMin;
  delete fHistMax;

  fCutLabel.clear();
  fCutEntryMin.clear();
  fCutEntryMax.clear();
  fCutEntryBin.clear();
  fHistEntryMin.clear();
  fHistEntryMax.clear();
  fHistEntryBin.clear();

  delete fCutCanvas;

  delete fHistOkButton; delete fCutOkButton; delete fTrigOkButton;
  delete fFrame4;
  delete fFrame3;
  delete fFrame2;
  delete fHistFrame;
  delete fHistCanvas;
  delete fCutFrame;

  delete fTab;
  delete fCancelButton; 
  delete fFrame1;
  delete fL1; delete fL2; delete fL3; delete fL4; delete fL5; delete fL6;
  delete fL7; delete fL8; delete fL9;
}

void PlanePanel::HandleButtons(int id) {

  if (id == -1) {
    TGButton *btn = (TGButton *) gTQSender;
    id = btn->WidgetId();
  }

  switch(id) {
  case M_HIST_OK:
    {
      std::vector<TH1*>& hists=fPlane->GetHistoList();

      if (thr_flag) TThread::Lock();
      int row=0;
      for(unsigned int i=0;i<hists.size();i++) {
	int bin[3];
	int min[3];
	int max[3];
	for(int j=0; j<hists[i]->GetDimension();j++) {
	  bin[j]=atoi(fHistEntryBin[row]->GetText());
	  min[j]=atoi(fHistEntryMin[row]->GetText());
	  max[j]=atoi(fHistEntryMax[row]->GetText());
	  row++;
	}

	switch(hists[i]->GetDimension()) {
	case 1:
	  hists[i]->SetBins(bin[0],min[0],max[0]);
	  break;
	case 2:
	  hists[i]->SetBins(bin[0],min[0],max[0],
			    bin[1],min[1],max[1]);	  
	  break;
	case 3:
	  hists[i]->SetBins(bin[0],min[0],max[0],
			    bin[1],min[1],max[1],
			    bin[2],min[2],max[2]);	
	  break;
	default:
	  break;
	}
	hists[i]->Reset();
      }
      if (thr_flag) TThread::UnLock();
      
      break;
    }
    
  case M_CUT_OK:
    {
      if(fPlane->IsInTree()) {
	new TGMsgBox(gClient->GetRoot(),fMain,
		     "message",
		     "This detector being in the tree, the whole tree will be reset. Continue ?",
		     kMBIconExclamation,
		     kMBCancel |  kMBOk);
      }
      std::vector<Variable*>& vars=fPlane->GetVariables(); 
      //vector<TH1*>& hists=fPlane->GetHistoList();
      for(unsigned int i=0;i<vars.size();i++) {
	int min=atoi(fCutEntryMin[i]->GetText());
	int max=atoi(fCutEntryMax[i]->GetText());

	vars[i]->Reset();
	vars[i]->SetRange(min,max);
      }
                        
      fPlane->ResetHistograms();
      fPlane->Modified(true);      
       
      break;
    }

  case M_TRIG_OK:
    std::cout<<"trigger tab OK button clicked"<<std::endl;
    fPlane->SetTriggerMask(fTMPattern);
    break;
  case M_PROPERTIES_CANCEL:
    fMain->SendCloseMessage();
    break;
  default:
    break;
  }
}

bool PlanePanel::ConvertHexa2Bin() {

  char **endptr=0;
  unsigned int tm = strtoul(fHexaMask->GetText(),endptr,16);

  if(strcmp(*endptr,"\0")==0 && strcmp(fHexaMask->GetText(),"\0") !=0 ) {
    int mask=1;
    for(size_t i=0;i<sizeof(int)*8;i++) {
      if(i == fBitButton.size()) break;
      fBitButton[i]->SetValue((tm & mask) >> i);
      mask=mask<<1;
    }
    return true;
  }
  else 
    return false;
}

bool PlanePanel::ConvertBin2Hexa() {
  
  fTMPattern=0;
  
  for(size_t i=0;i<fBitButton.size();i++) {
    if(i == sizeof(int)*8) break;
    fTMPattern += fBitButton[i]->GetValue() << i;
  } 
  
  char txt[10];
  sprintf(txt,"%x",fTMPattern);
  std::string txthm = "0x";
  txthm += txt;
  
  fHexaMask->GetBuffer()->Clear();
  fHexaMask->GetBuffer()->AddText(0, txthm.c_str());
  gClient->NeedRedraw(fHexaMask);
  
  return true;
}

void PlanePanel::TTSelected(int id) {
  
  std::cout<<"trigger type "<<id<<" selected for plane "<<fPlane->GetName()<<std::endl;
}



