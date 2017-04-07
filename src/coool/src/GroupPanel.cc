#include "GroupPanel.h"

ClassImp(GroupPanel);


enum  EPropCommandIdentifiers {
  M_PROPERTIES_CANCEL,
  M_HIST_OK,
  M_CUT_OK,
};


GroupPanel:: GroupPanel(const TGWindow *p,const TGWindow *main,
			      UInt_t w, UInt_t h, Group *Group) :
  TObject(),fGroup(Group) {

  fMain = new TGTransientFrame(p,main,w,h);
  fMain->Connect("CloseWindow()", "GroupPanel", this,
  		 "CloseWindow()");

  fFrame1 = new TGHorizontalFrame(fMain, 60, 20, kFixedWidth);


  fCancelButton = new TGTextButton(fFrame1, "&Cancel", M_PROPERTIES_CANCEL);
  fCancelButton->Connect("Clicked()","GroupPanel",this,"HandleButtons()");
  fL1 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,
			  2, 2, 2, 2);
  fL2 = new TGLayoutHints(kLHintsBottom | kLHintsRight, 2, 2, 5, 1);
  fL3 = new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
					 kLHintsExpandY, 2, 2, 5, 1);

  fL4 = new TGLayoutHints(kLHintsRight, 2, 2, 5, 1);
  fL5 = new TGLayoutHints(kLHintsLeft, 2, 2, 5, 1);
  fL6 = new TGLayoutHints(kLHintsCenterY | kLHintsCenterX, 2, 2, 5, 1);



  fFrame1->AddFrame(fCancelButton, fL1);

  fFrame1->Resize(150, fCancelButton->GetDefaultHeight());
  fMain->AddFrame(fFrame1, fL2);

  //--------- create Tab widget and some composite frames for Tab testing
  fTab = new TGTab(fMain, 300, 300);

  //Tab for histograms
  TGCompositeFrame *tf=fTab->AddTab("Histograms");
  tf->SetLayoutManager(new TGVerticalLayout(tf));

  fHistFrame=new TGCompositeFrame(tf, 60, 20);
  fHistFrame->SetLayoutManager(new TGMatrixLayout(fHistFrame,0,4));
  std::vector<TH1*>& hists=fGroup->GetHistoList();

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

      TAxis *axis = 0;
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
  tf->AddFrame(fHistFrame);

  fHistOkButton = new TGTextButton(tf, "   &Ok  ", M_HIST_OK);
  fHistOkButton->Connect("Clicked()","GroupPanel",this,"HandleButtons()");

  tf->AddFrame(fHistOkButton, fL6);

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

  fMain->SetWindowName(fGroup->GetName());
  fMain->MapSubwindows();
  fMain->Resize(fMain->GetDefaultSize());
  fMain->MapWindow();

}


GroupPanel::~GroupPanel() {

  delete fHistVar;
  delete fHistBin;
  delete fHistMin;
  delete fHistMax;

  fHistEntryMin.clear();
  fHistEntryMax.clear();
  fHistEntryBin.clear();

  delete fHistOkButton;
  delete fHistFrame;

  delete fTab;
  delete fCancelButton;
  delete fFrame1;
  delete fL1; delete fL2; delete fL3; delete fL4; delete fL5; delete fL6;
}


void GroupPanel::HandleButtons(int id) {

  if (id == -1) {
    TGButton *btn = (TGButton *) gTQSender;
    id = btn->WidgetId();
  }

  switch(id) {
  case M_HIST_OK:
    {
      std::vector<TH1*>& hists=fGroup->GetHistoList();

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

  case M_PROPERTIES_CANCEL:
    fMain->SendCloseMessage();
    break;
  default:
    break;
  }
}

