#include "MonitorPanel.h"

ClassImp(MonitorPanel);


enum  EPropCommandIdentifiers {
  M_PROPERTIES_CANCEL,
  M_TRIG_OK,
  M_TT_COMBO
};


MonitorPanel:: MonitorPanel(const TGWindow *p,const TGWindow *main,
			    UInt_t w, UInt_t h, Monitor *monitor,
			    unsigned int trigmask, bool lowerOnly, bool strict) :
  TObject(), fTMPattern(trigmask), fTMLowerOnly(lowerOnly),
  fTMStrict(strict), fMonitor(monitor) {

  fMain = new TGTransientFrame(p,main,w,h);
  fMain->Connect("CloseWindow()", "MonitorPanel", this,
  		 "CloseWindow()");

  fL1 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,
			  2, 2, 2, 2);
  fL2 = new TGLayoutHints(kLHintsBottom | kLHintsRight, 2, 4, 5, 2);
  fL3 = new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
					 kLHintsExpandY, 2, 2, 5, 1);

  fL4 = new TGLayoutHints(kLHintsCenterX | kLHintsExpandX, 2, 2, 5, 1);
  fL5 = new TGLayoutHints(kLHintsLeft, 2, 2, 5, 1);
  fL6 = new TGLayoutHints(kLHintsCenterY | kLHintsCenterX, 2, 2, 5, 1);
  fL7 = new TGLayoutHints(kLHintsBottom | kLHintsRight, 0, 0, 0, 0);
  fL8 = new TGLayoutHints(kLHintsTop | kLHintsExpandX |
					 kLHintsExpandY, 2, 2, 5, 1);
  fL9 = new TGLayoutHints(kLHintsBottom | kLHintsRight, 2, 2, 5, 1);

  fFrame1 = new TGHorizontalFrame(fMain, 60, 20, kFixedWidth);
  fCancelButton = new TGTextButton(fFrame1, "&Close", M_PROPERTIES_CANCEL);
  fCancelButton->Connect("Clicked()","MonitorPanel",this,"HandleButtons()");
  fFrame1->AddFrame(fCancelButton, fL1);
  fFrame1->Resize(150, fCancelButton->GetDefaultHeight());
  fMain->AddFrame(fFrame1, fL2);


  //--------- create Tab widget and some composite frames for Tab testing
  fTab = new TGTab(fMain, 300, 300);


  // tab for trigger
  TGCompositeFrame *tf=fTab->AddTab("Trigger");
  tf->SetLayoutManager(new TGVerticalLayout(tf));

  // trigger type
  fTrigType = new TGGroupFrame(tf, "Trigger Type", kVerticalFrame);
  const std::map<int,std::string>& evtypes = fMonitor->GetEvtTypes();
  typedef std::map<int,std::string>::const_iterator IT;
  for(IT it=evtypes.begin(); it!=evtypes.end(); it++) {
    fEvtTypeChecks.push_back(new TGCheckButton(fTrigType, it->second.c_str()
					       ,it->first));
    fEvtTypeChecks.back()->Connect("Clicked()","MonitorPanel",this,"EvtTypeChecked(int)");
    fTrigType->AddFrame(fEvtTypeChecks.back());
  }
  
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
    //std::cout<<bytei<<std::endl;
    BitButton *button=new BitButton(fByte[bytei],i);
    button->SetValue((fTMPattern>>i) & 1);
    button->Connect("Clicked()","MonitorPanel",this,"ConvertBin2Hexa()");
    fBitButton.push_back(button);
    fByte[bytei]->AddFrame(button,fL7);
  }
  for(int i=3; i>=0; i--) {
    fTrigMask1->AddFrame(fByte[i]);
  }

  fHexaMask = new TGTextEntry(fTrigMask,"0x00000000");
  fHexaMask->Connect("ReturnPressed()","MonitorPanel",this,
		       "ConvertHexa2Bin()");

  fCBLower = new TGCheckButton(fTrigMask, "only compare lower 16 bits", 1);
  fCBLower->Connect("Clicked()","MonitorPanel",this,"TrigMaskChecked(int)");
  if (fTMLowerOnly)
    fCBLower->SetState(kButtonDown);
  else
    fCBLower->SetState(kButtonUp);

  fCBStrict = new TGCheckButton(fTrigMask, "strict comparison", 2);
  fCBStrict->Connect("Clicked()","MonitorPanel",this,"TrigMaskChecked(int)");
  if (fTMStrict)
    fCBStrict->SetState(kButtonDown);
  else
    fCBStrict->SetState(kButtonUp);

  fTrigMask->AddFrame(fTrigMask1);
  fTrigMask->AddFrame(fHexaMask,fL6);
  fTrigMask->AddFrame(fCBLower);
  fTrigMask->AddFrame(fCBStrict);

  tf->AddFrame(fTrigType,fL1);
  tf->AddFrame(fTrigMask,fL1);

  fFrame4 = new TGHorizontalFrame(tf, 60, 20, kFixedWidth);
  fTrigOkButton = new TGTextButton(fFrame4, "   &Ok  ", M_TRIG_OK);
  fTrigOkButton->Connect("Clicked()","MonitorPanel",this,"HandleButtons()");
  fFrame4->AddFrame(fTrigOkButton, fL1);
  tf->AddFrame(fFrame4, fL9);
  fFrame4->Resize(150, fTrigOkButton->GetDefaultHeight());

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
  fMain->SetWindowName(fMonitor->GetName());
  fMain->MapSubwindows();
  fMain->Resize(fMain->GetDefaultSize());
  fMain->MapWindow();
  ConvertBin2Hexa();
}

MonitorPanel::~MonitorPanel() {

  delete fHexaMask;
  fByte.clear();
  fBitButton.clear();
  delete fTrigMask1;
  delete fTrigMask;
  delete fTrigType;

  delete fTrigOkButton;
  delete fFrame4;



  delete fTab;
  delete fCancelButton;
  delete fFrame1;
  delete fL1; delete fL2; delete fL3; delete fL4; delete fL5; delete fL6;
  delete fL7; delete fL8; delete fL9;
}

void MonitorPanel::HandleButtons(int id) {

  if (id == -1) {
    TGButton *btn = (TGButton *) gTQSender;
    id = btn->WidgetId();
  }

  switch(id) {

  case M_TRIG_OK:
    if (thr_flag) TThread::Lock();
    fMonitor->SetTriggerMask(fTMPattern);
    fMonitor->SetTriggerMaskLowerOnly(fTMLowerOnly);
    fMonitor->SetTriggerMaskStrict(fTMStrict);
    fMonitor->CheckEvtTypes(fEvtTypesChecked);
    if (thr_flag) TThread::UnLock();
    break;
  case M_PROPERTIES_CANCEL:
    fMain->SendCloseMessage();
    break;
  default:
    break;
  }
}

bool MonitorPanel::ConvertHexa2Bin() {

  char **endptr=0;
  fTMPattern = strtoul(fHexaMask->GetText(),endptr,16);
  
  int mask=1;
  for(size_t i=0;i<sizeof(int)*8;i++) {
    if(i == fBitButton.size()) break;
    fBitButton[i]->SetValue((fTMPattern & mask) >> i);
    mask=mask<<1;
  }
  return true;
}

bool MonitorPanel::ConvertBin2Hexa() {

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

void MonitorPanel::EvtTypeChecked(int id) {
  TGButton *btn = (TGButton *) gTQSender;
  id = btn->WidgetId();
  
  if(btn->GetState()==kButtonDown) {      
    fEvtTypesChecked.insert(id);
  }
  else {
    std::set<int>::iterator it = fEvtTypesChecked.find(id);
    if(it!=fEvtTypesChecked.end()) {
      fEvtTypesChecked.erase(it);
    }
    else 
      std::cerr<<"MonitorPanel::EvtTypeChecked ERROR !"<<std::endl;
  }  
}

void MonitorPanel::CheckEvtType(int id) {
  for(unsigned i=0; i<fEvtTypeChecks.size(); i++) {
    if(fEvtTypeChecks[i]->WidgetId() == id) {
      fEvtTypeChecks[i]->SetState(kButtonDown);
      fEvtTypesChecked.insert(id);
      break;
    } 
  }
}

void MonitorPanel::TrigMaskChecked(int id) {
  TGButton *btn = (TGButton *) gTQSender;
  id = btn->WidgetId();
  
  bool state(false);
  if (btn->GetState()==kButtonDown)
    state=true;

  switch (id) {
    case 1:
      fTMLowerOnly = state;
      break;
    case 2:
      fTMStrict = state;
      break;
    default:
      std::cerr << "MonitorPanel::TriggerMaskChecked: called with wrong ID: " << id << std::endl;
  }
}





