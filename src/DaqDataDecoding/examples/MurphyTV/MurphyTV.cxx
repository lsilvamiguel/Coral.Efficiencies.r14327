
#include "MurphyTV.h"


/*****************************************************************************************************
*****************************************************************************************************/
#include "errsdescr.h"
bool MurphyTV::alive=false;
  
MurphyTV::MurphyTV(CollectErrs *ce,TGWindow *main, bool calibrationTrigger, bool small)
 :TGTransientFrame(gClient->GetRoot(),main,(1280-small*620),(960-small*480)){
  alive=true;
  trash=ce;
   
  if (small)  SetWMSizeHints(600,400,1800,2600,5,5);
  else        SetWMSizeHints(800,600,1800,2800,5,5);
  

  fTab = new TGTab(this, 800, 600);

  fL1 = new TGLayoutHints(kLHintsBottom | kLHintsCenterX , 0, 0, 8, 8);
  fL2 = new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0);
  fL3 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2);

  TGCompositeFrame *tf = fTab->AddTab("Catches");  
  
  catchlist=new CatchList(tf,300,300,trash,CATCHLIST,calibrationTrigger, small);
  tf->AddFrame(catchlist,fL3);
  catchlist->Associate(this);

  tf = fTab->AddTab("Errors");  
  errslist =new ErrsList(tf,300,300,trash,ERRSLIST,small);
  tf->AddFrame(errslist,fL3);
  errslist->Associate(this);  

  tf = fTab->AddTab("Monitor");
  monitorhist=new HistMonitor(tf,300,300,trash);
  tf->AddFrame(monitorhist, fL3);  
   
  tf = fTab->AddTab("History");
  history = new HistHist(tf,300,300,trash);
  tf->AddFrame(history,fL3);  
  
  tf = fTab->AddTab("Event Sizes");  
  sizeshist = new HistSizes(tf,300,300,trash);
  tf->AddFrame(sizeshist,fL3);
  
  tf = fTab->AddTab("Data Flow");  
  dsizeshist = new HistHistD(tf,300,300,trash);
  tf->AddFrame(dsizeshist,fL3);  
  
  AddFrame(fTab, fL3);
  
  fTBClose=new TGTextButton(this,"Close",TB_CLOSE);
  AddFrame(fTBClose,fL1);

  

  MapSubwindows();

  Layout();
  SetWindowName("MurphyTV - Online COMPASS data monitoring");
  SetIconName("MurphyTV");
  MapWindow();

  fTimer=new TTimer(500);
  fTimer->SetObject(this);
  fTimer->TurnOn(); 

}

MurphyTV::~MurphyTV(){

  fTimer->TurnOff();

  delete fTimer;
  
  delete fTBClose;
  delete catchlist;
  delete errslist;
  delete monitorhist;
  delete history;
  delete sizeshist;
  delete dsizeshist;

  delete fL1;
  delete fL2;
  delete fL3;
  delete fTab;

  alive=false;

}

Bool_t MurphyTV::HandleTimer(TTimer *timer)
{
  tdiv=(tdiv+1) % 60;
  
  switch(fTab->GetCurrent()){
    case 0:
      catchlist->Update(tdiv == 0, tdiv == 0);
      break;
    case 1:
      errslist->Update(tdiv == 0);
      break;
    case 2:
      if ((tdiv % 6)==0)
      monitorhist->Update();
      break;
    case 3:
      if ((tdiv % 6)==0)
      history->Update();
      break;
    case 4:
      if ((tdiv % 6)==0)
      sizeshist->Update();
      break;
    case 5:
      if ((tdiv % 6)==0)
      dsizeshist->Update();
      break;
  
  }
  return kTRUE;
}


  
Bool_t MurphyTV::ProcessMessage(Long_t msg,Long_t parm1,Long_t parm2){
  switch(GET_MSG(msg)) {
    case kC_COMMAND:
      switch(GET_SUBMSG(msg)){

        case kCM_BUTTON:
          switch(parm1) {
            case CATCHLIST:
              errslist->ChoseSource(parm2);
              fTab->SetTab(1);
              break;
            case ERRSLIST:
              catchlist->SelectSpecial((DaqErrorType)parm2);
              fTab->SetTab(0);
              break;
            case TB_CLOSE:
              CloseWindow();
	            break;
            default:
              break;
          }
          break;
        default:
	        break;
      }
      break;
    default:
      break;
  }
  return kTRUE;
}
void MurphyTV::CloseWindow()
{
  delete this;
}


/*****************************************************************************************************
*****************************************************************************************************/

SorTable::SorTable(TGWindow *p,UInt_t wi, UInt_t h,UInt_t opt,Int_t id, bool small)
:TGCompositeFrame(p, wi, h, opt | kVerticalFrame,fgWhitePixel)
{

  MsgWindow=p;
  sortby=0;
  
  
  if (small) fnt=gClient->GetFontByName("fixed");
  else                 fnt=gClient->GetFontByName("-*-courier-bold-r-*-*-20-*");
  cwidth=gVirtualX->TextWidth(fnt, "0", 1);
  lastclick=gSystem->Now();
  lastclicked=0xffffffff; 
  
  gClient->GetColorByName("red",red);
  gClient->GetColorByName("yellow",yellow);
  gClient->GetColorByName("orange",orange);
  gClient->GetColorByName("white",white);
 
  
  fLH = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0);
  fLH2 = new TGLayoutHints(kLHintsExpandY|kLHintsLeft ,0, 0, 0, 0); 
  fLH1 = new TGLayoutHints(kLHintsExpandX|kLHintsTop  ,0, 0, 0, 0);     
  

  fCF= new TGCompositeFrame(this, 300, 20,kHorizontalFrame|kFixedHeight );
  fTBFill=0;
  AddFrame(fCF,fLH1);
  
  fLB=new TGListBox(this,id,0);
  fLB->IntegralHeight(kFALSE);
  fLB->Resize(100,100);

  tggc=new TGGC(TGTextLBEntry::GetDefaultGC());
  
  tggc->SetFont(gVirtualX->GetFontHandle(fnt));
  fLB->AddEntry(new TGTextLBEntry(fLB->GetContainer(),new TGString("empty"),1,tggc->GetGC(),fnt),
                                  new TGLayoutHints(kLHintsTop|kLHintsLeft,0,0,0,0)); 

  AddFrame(fLB,fLH);    
}

SorTable::~SorTable()
{
  delete fLB;  
  for (vector<HCollumn *>::iterator it=HCols.begin();it!=HCols.end();it++) delete (*it);
  delete fTBFill; 
  delete fCF;
  delete fLH;
  delete fLH1;
  delete fLH2;
  delete tggc;
}

void SorTable::AddCollumn(char *title,uint32 width,bool sortdir,char *tooltip) {
  uint32 col=HCols.size();
  HCollumn *hc=new HCollumn(fCF,col,title,width,cwidth,sortdir,tooltip);
  HCols.push_back(hc);
  fCF->AddFrame(hc,fLH2);
  hc->Associate(this);
}

void SorTable::LastCollumn(){
  fTBFill=new TGTextButton(fCF,"",0);
  fCF->AddFrame(fTBFill,fLH);
  Layout();
}

void SorTable::AddRow(Int_t Id) {
  Row r(HCols.size(),Id);
  Data.push_back(r);
}

SorTable::Entry &SorTable::GetEntry(uint32 col,uint32 row) {
  if (row<Data.size()) if (col<Data[row].size()) return Data[row][col]; 
  return defdum;
}
void SorTable::SortBy(uint32 col,bool dir) {
  if (col>=HCols.size()) return;
  HCols[col]->sortdir=dir;
  sortby=col;
}



int SorTable::c_sort_h(const void *a1,const void *a2){
  Entry &e1=*((scs *)a1)->val;
  Entry &e2=*((scs *)a2)->val;
  if (e1.Disp > e2.Disp ) return 1;
  else if (e1.Disp < e2.Disp ) return -1;
  else {
    if (e1.Disp==Entry::INT) {
      if (e1.IntVal > e2.IntVal) return 1; else if (e1.IntVal < e2.IntVal) return -1;
    }
    else if (e1.Disp==Entry::DOUBLE) {
      if (e1.DoubleVal > e2.DoubleVal) return 1; else if (e1.DoubleVal < e2.DoubleVal) return -1;
    }
    else if (e1.Disp==Entry::TEXT) return e1.TextVal.CompareTo(e2.TextVal);
  }
  return 0;
}
 
void SorTable::Update()
{
  uint32 a,b,c,d;
  uint32 nr=Data.size();
  Int_t lsel=fLB->GetSelected();
  fLB->RemoveEntries(-32767,32768);  
     
  
  if (nr==0 || HCols.size()==0) {
     fLB->AddEntry(new TGTextLBEntry(fLB->GetContainer(),new TGString("empty   "),0,tggc->GetGC(),fnt),
                                      new TGLayoutHints(kLHintsTop|kLHintsLeft,0,0,0,0)); 
  }   
  else {
    scs s[nr];
    char buf[5000];
    
    if (sortby>=HCols.size()) sortby=0;
    for (a=0;a<nr;a++) {
      s[a].nr=a;
      s[a].val=& Data[a][sortby];
    }
    qsort(&s[0],nr,sizeof(scs),c_sort_h);

    for (a=0;a<nr;a++){    
      if (HCols[sortby]->sortdir) d=s[a].nr; else d=s[nr-a-1].nr;    
      for (b=c=0;b<HCols.size();b++) {
        Entry &e=Data[d][b];
        int w=HCols[b]->Width;
        color=white;
	if (Data[d][2].DoubleVal>=0.005) color=yellow;
	if (Data[d][2].DoubleVal>=0.05) color=orange;
	if (Data[d][2].DoubleVal>=0.50) color=red;
        if (e.Disp==Entry::TEXT) sprintf(buf+c,"%-*s    ",w-3,e.TextVal.Data());
        else if (e.Disp==Entry::DOUBLE&&HCols.size()>9&&b==2) sprintf(buf+c,"%*.1f%%    ",w-3,e.DoubleVal*100.0);
        else if (e.Disp==Entry::DOUBLE) sprintf(buf+c,"%*.2f    ",w-3,e.DoubleVal);
        else if (e.Disp==Entry::INT) sprintf(buf+c,"%*d    ",w-3,e.IntVal);
        else for (int fb=0;fb<w+3;fb++) buf[c+fb]=' ';
        c+=w;
      } 
      fLB->AddEntry(new TGTextLBEntry(fLB->GetContainer(),new TGString(buf),d,tggc->GetGC(),fnt,0,color),
                                        new TGLayoutHints(kLHintsTop|kLHintsLeft,0,0,0,0));    									                            
    } 
    ((TGLBContainer *)(fLB->GetContainer()))->Select(lsel);
  }
  fLB->MapSubwindows();
  fLB->Layout();
}



Bool_t SorTable::ProcessMessage(Long_t msg,Long_t parm1,Long_t parm2){

  if (GET_MSG(msg)==kC_COMMAND){
    if ( GET_SUBMSG(msg)==kCM_BUTTON){
      if (0<=parm1 && parm1<(Long_t)HCols.size()) {
        if ((Long_t)sortby==parm1) HCols[sortby]->sortdir=!HCols[sortby]->sortdir;
        sortby=parm1;
        Update();
      } 
    }
    if ( GET_SUBMSG(msg)==kCM_LISTBOX && 0<=parm2 && parm2<(Long_t)Data.size()) {
      if ((long)(gSystem->Now()-lastclick)>250 || lastclicked!=parm2) SendMessage(MsgWindow, msg,parm1, Data[parm2].Id);
      else SendMessage(MsgWindow,MK_MSG(kC_COMMAND,kCM_BUTTON), parm1,Data[parm2].Id); 
      lastclick=gSystem->Now();
      lastclicked=parm2;
    }
  }        
  return kTRUE;
}

SorTable::HCollumn::HCollumn(TGWindow *p,uint32 id,char *txt,uint32 w,uint32 cw,bool sdir,char *tooltip) 
:TGCompositeFrame(p, w*cw, 20, kFixedWidth )
{
  Width=w;
  sortdir=sdir;
  Button=new TGTextButton(this,new TGHotString(txt),id);
  fLH = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY ,0, 0, 0, 0);
  if (tooltip) Button->SetToolTipText(tooltip,1000);
  AddFrame(Button,fLH);
} 

SorTable::Entry::Entry(Int_t val,Int_t norm,bool normit) {
  if (norm==0) {TextVal="no data"; Disp=TEXT;}
  else if (normit) {DoubleVal=(Double_t)val/(Double_t)norm; Disp=DOUBLE;}
  else {IntVal=val; Disp=INT;}
}

/*****************************************************************************************************
*****************************************************************************************************/

MurphyTV::CatchList::CatchList(TGWindow *p,UInt_t w, UInt_t
h,CollectErrs*ce,Int_t id, bool calibrationTrigger, bool small)
:TGCompositeFrame(p, w, h, kVerticalFrame)
{
  uint32 a;
  MsgWindow=p;
  Id=id;
  trash=ce;
  special=DaqError::EVENT_MISSING_SRCID;
  asource=0;
  dispsrc=0xffffff;
  
  fLH = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0);  
  fLH1 = new TGLayoutHints(kLHintsCenterY | kLHintsLeft ,5, 5, 5, 5);  
  fLH2 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX ,0,0,0,0); //2, 2, 2, 2); 
  fLH3 = new TGLayoutHints(kLHintsCenterY | kLHintsExpandX ,5, 5, 5, 5);
  fLH4 = new TGLayoutHints(kLHintsExpandY | kLHintsLeft ,0, 0, 0, 0);
  fLH5 = new TGLayoutHints(kLHintsLeft | kLHintsExpandY ,80-small*70, 0, 0, 0);
  
  
  fCF = new TGCompositeFrame(this, 100, 25, kHorizontalFrame);
 
  fLabel = new TGLabel(fCF,new TGHotString("Special="));
  fCF->AddFrame(fLabel,fLH1);
  fCombo = new TGComboBox(fCF,CB_SPECIAL);
  fCF->AddFrame(fCombo,fLH1);
  fCheck = new TGCheckButton(fCF,"Normalize");
  fCheck->SetState(kButtonDown);
  fCheck->SetToolTipText("Norm everything on the number of recorded events",1000);
  fCheck->Associate(this);
  fCF->AddFrame(fCheck,fLH1);
 
  AddFrame(fCF,fLH2);

  for (a=0 ;a<ErrZoo::GetNrErrs();a++ ) 
    fCombo->AddEntry(ErrZoo::GetTitle(a), ErrZoo::GetType(a));
  fCombo->Select(special);
  fCombo->Resize(180, 20);
  fCombo->Associate(this);
  
  if (calibrationTrigger) {
  fCalibrationTrigger=new TGLabel(this, "Calibration Trigger ONLY");
  gClient->GetColorByName("red",red);
  fCalibrationTrigger->SetTextJustify(kTextLeft);
  fCalibrationTrigger->ChangeBackground(red);
  fCF->AddFrame(fCalibrationTrigger,fLH5);
  }
  
 
    
  fRunNumber =  new TGLabel(fCF,"Run number: n.a.");
  fRunNumber->SetTextJustify(kTextLeft);
  fCF->AddFrame(fRunNumber, fLH5);
  
  
  
  fSpillNumber =  new TGLabel(fCF,"Spill number: n.a.");
  fSpillNumber->SetTextJustify(kTextLeft);
  fCF->AddFrame(fSpillNumber, fLH5);
  
     
  fLastUpdate =  new TGLabel(fCF,"Last event: n.a.");
  fLastUpdate->SetTextJustify(kTextLeft);
  fCF->AddFrame(fLastUpdate, fLH5);
 
                             
  catchtab=new SorTable(this,w,h,kSunkenFrame|kDoubleBorder,LB_CATCHES, small);
 
  catchtab->AddCollumn("SourceID",10,true);
  catchtab->AddCollumn("Type",10,true,"Detector Type");
  catchtab->AddCollumn("BadEvents",10,false,"Number of events which had errors on this SrcID");
  catchtab->AddCollumn("#Header",10,true,"Number of header words (or 0 for GEM)");
  catchtab->AddCollumn("#Data",10,true,"Number of data words");
  catchtab->AddCollumn("#Errors",10,false,"Number of errors on SourceID");
  catchtab->AddCollumn("Special",10,false,"How often the error selected in \"special\" combo box occured");
  catchtab->AddCollumn("Special at",10,true,"Trigger# when error \"Special\" occured the first time");
  catchtab->AddCollumn("#Spills",10,true,"Nb of spills where this SrcID is affected by an error");
  catchtab->AddCollumn("Last spill",10,true,"Last spill where this SrcID has got an error");
  catchtab->LastCollumn();
  catchtab->SortBy(2,false);

  AddFrame(catchtab,fLH); 


  fCF2= new TGCompositeFrame(this,300,120,kHorizontalFrame );
  
  fGF=new TGGroupFrame(fCF2,"Frontends on selected SourceID");  

  porttab=new SorTable(fGF,300,120,kSunkenFrame|kDoubleBorder,LB_PORTS, small);
  porttab->AddCollumn("Port#",10,true,"Frontend input port number on CATCH");
  porttab->AddCollumn("GeoID/Value",10,true,"Frontend geographicID/Value of error word");
  porttab->AddCollumn("#Errors",10,false,"Number of errors on frontend");
  porttab->AddCollumn("#Header",10,true,"Number of header words from frontend");
  porttab->AddCollumn("#Data",10,true,"Number of data words from frontend");
  porttab->AddCollumn("Special",10,false,"How often the error selected in \"special\" combo box occured");
  porttab->LastCollumn();
  porttab->SortBy(2,false);
  
  fGF->AddFrame(porttab,fLH);
  fCF2->AddFrame(fGF,fLH);
  fGF2=new TGGroupFrame(fCF2,"Attached Detectors");
  fTVDets=new TGTextView(fGF2,150,120,-1,0,fGF->GetDefaultFrameBackground());
  fTVDets->ChangeOptions(0);
  fGF2->AddFrame(fTVDets,fLH4);
  
  fCF2->AddFrame(fGF2,fLH4);
  
  AddFrame(fCF2,fLH);
     
  Update(true);


}


MurphyTV::CatchList::~CatchList()
{ 
  delete fLastUpdate;
  delete catchtab;
  delete porttab;


  delete fTVDets;
  delete fGF;
  delete fGF2;
  
  delete fCheck;
  delete fCombo;
  delete fLabel;
  delete fCF;
  delete fCF2;
  delete fLH5;
  delete fLH4;
  delete fLH3;
  delete fLH2;
  delete fLH1;
  delete fLH;
}

const char *MurphyTV::CatchList::SrcType(uint32 src) {
  struct Alias {uint32 min ;uint32 max ;char txt[16]; };
  static const Alias atb[]={{0,0,"DAQ"},
     {1,1,"Filter"},{2,9,"Master-T"},{10,20,"Scaler"},{21,22,"ScalerFi15"},{23,60,"Scaler"},{61,74,"Trigger"},
     {75,79,"Veto"},{80,99,"Trigger"},{100,101,"RPD"},{102,102,"Munich"},
     {103,103,"RPD"},{110,110,"Munich"},
     {128,129,"SciFi15"},{144,146,"SciFi03"},{147,149,"SciFi04"},
     {150,219,"SciFi-D"},{220,239,"SciFi-W"},
     {240,249,"BMS"},{250,250,"Cedar"},
     {251,270,"DC"},{271,319,"W45"},
     {320,332,"Straw"},{380,383,"PMM"},
     {384,400,"MicroM"},{401,430,"MW2"},{431,440,"RW"},{441,479,"MWPC"},
     {480,499,"MW1"},{500,539,"RICH-MAPMT"},{540,600,"RICH-APV"},
     {601,610,"HCAL"},{611,639,"ECAL"},
     {640,699,"Silicon"},{700,700,"CedarADC"},{701,749,"GEM"},{750,799,"PGEM"},
     {800,800,"MT-GADC"},{801,801,"MT-GTDC"},
     {810,811,"CedarADC"},
     {820,825,"LaCameraA"},{830,835,"LaCameraB"},
     {840,841,"SciFi-ADC"},{850,857,"SciFi35"},
     {860,861,"SciFi01"},{870,871,"ScalerFI01"},
     {880,880,"TIGER-CA"},{881,881,"TIGER-FI"},
     {882,889,"DC05"},
     {896,1023,"SMUX"}
    };
  for (uint32 a=0 ;a<sizeof(atb)/sizeof(Alias) ;a++) if (atb[a].min<=src && atb[a].max>=src) return atb[a].txt;
  return "- ??? -";

}

void MurphyTV::CatchList::SrcInfos(uint32 src) {
  if (asource!=dispsrc) {
    set <string> names;
    trash->GetDetsAtCatch(asource,names);
    fTVDets->Clear();
    for (set <string>::iterator it = names.begin();it!=names.end();it++) 
      fTVDets->AddLine(it->c_str());
    names.clear();
    dispsrc=asource;
  }
}

void MurphyTV::CatchList::Update(bool redraw,bool redraw2){
  int a,b;
  int NrSrc=0;
  static uint32 RunNr=0;
  
  catchtab->ClearRows();
  porttab->ClearRows();
  
  trash->Lock();
  time_t EvntTime=trash->GetEventTime();
  
  bool dnorm=( fCheck->GetState()==kButtonDown ); 
  uint32 norm=trash->GetNrRecEvents();
  for (a=b=0;a<CollectErrs::MAX_SOURCE;a++){
    if (trash->CatchExists(a)){ 
      NrSrc++;
      catchtab->AddRow(a); 
      catchtab->GetEntry(0,b)=SorTable::Entry(a);
      catchtab->GetEntry(1,b)=SorTable::Entry(SrcType(a));
      catchtab->GetEntry(2,b)=SorTable::Entry(trash->GetCatchNrBadEvnts(a),norm,dnorm);
      catchtab->GetEntry(3,b)=SorTable::Entry(trash->GetNrHeaderOnCatch(a),norm,dnorm);
      catchtab->GetEntry(4,b)=SorTable::Entry(trash->GetNrDataOnCatch(a),norm,dnorm);
      catchtab->GetEntry(5,b)=SorTable::Entry(trash->GetCatchErrSum(a),norm,dnorm);
      catchtab->GetEntry(6,b)=SorTable::Entry(trash->GetErrOnCatch(a,special),norm,dnorm);
      catchtab->GetEntry(7,b)=SorTable::Entry((Int_t)trash->GetFirstErrPosOnCatch(a,special));
      catchtab->GetEntry(8,b)=SorTable::Entry(trash->GetNbSpillWithErrGen(a),norm,kFALSE);
      catchtab->GetEntry(9,b)=SorTable::Entry(trash->GetLastSpillWithErrGen(a),norm,kFALSE);
      b++;
    }
  }
  
  
  map <uint16,uint32> prtocc;
  map <uint16,uint32> geoocc;  

  const map <CollectErrs::VPortID,CollectErrs::VPortStat> &pst=trash->GetVPortStats(asource);
  for (CollectErrs::VPortStatMap::const_iterator it=pst.begin();it!=pst.end();it++) {
    uint16 ap=it->first.Port();
    uint16 ag=it->first.GeoID();
    if (ap!=CollectErrs::INVALID && ag!=CollectErrs::INVALID) {
      prtocc[ap]++;
      geoocc[ag]++;
    }
  }
 
  b=0;
  for (CollectErrs::VPortStatMap::const_iterator it=pst.begin();it!=pst.end();it++) {
    uint16 ap=it->first.Port();
    uint16 ag=it->first.GeoID();
    uint32 apo=prtocc[ap];
    uint32 ago=geoocc[ag];

    if (ap==CollectErrs::INVALID && ago==1) continue;
    if (ag==CollectErrs::INVALID && apo==1) continue;

    CollectErrs::VPortStat vps=it->second;

    if (ap!=CollectErrs::INVALID && ag!=CollectErrs::INVALID){
      if (apo==1) {
        CollectErrs::VPortStatMap::const_iterator it2=pst.find(CollectErrs::VPortID(ap,CollectErrs::INVALID));
        if (it2!=pst.end()) vps+=it2->second;
      }
      if (ago==1) {
        CollectErrs::VPortStatMap::const_iterator it2=pst.find(CollectErrs::VPortID(CollectErrs::INVALID,ag));
        if (it2!=pst.end()) vps+=it2->second;
      }
    }

    porttab->AddRow(b);
    if (ap==CollectErrs::INVALID) porttab->GetEntry(0,b)=SorTable::Entry("   ");
    else porttab->GetEntry(0,b)=SorTable::Entry(ap);
    if (ag==CollectErrs::INVALID) porttab->GetEntry(1,b)=SorTable::Entry("   ");
    else porttab->GetEntry(1,b)=SorTable::Entry(ag);
    porttab->GetEntry(2,b)=SorTable::Entry(vps.NumAnyErr,norm,dnorm);
    porttab->GetEntry(3,b)=SorTable::Entry(vps.NumHeader,norm,dnorm);
    porttab->GetEntry(4,b)=SorTable::Entry(vps.NumData,norm,dnorm);
    porttab->GetEntry(5,b)=SorTable::Entry(vps.GetNumErr(special),norm,dnorm);
    b++;
  }

  trash->Unlock();
  
  if (redraw) catchtab->Update();
  if (redraw2) {
    porttab->Update();  
    if (asource!=dispsrc) {
      set <string> names;
      trash->GetDetsAtCatch(asource,names);
      fTVDets->Clear();
      if (names.empty()) fTVDets->AddLine("not available");
      for (set <string>::iterator it = names.begin();it!=names.end();it++) 
        fTVDets->AddLine(it->c_str());
        
      names.clear();
      dispsrc=asource;
      
    }
  }
  
  char buf[100];
  char buf2[100];
  ULong_t red, green, orange, grey;
  gClient->GetColorByName("red",red);
  gClient->GetColorByName("green",green);
  gClient->GetColorByName("orange",orange);
  grey=fCF->GetDefaultFrameBackground();
  
  uint32 nre=trash->GetNrRecEvents();
  uint32 SpillNr=trash->GetEventSpillNb();
   
  if (nre!=0) {
    sprintf(buf,"                                    ");
    strncpy(buf,ctime(&EvntTime),strlen(ctime(&EvntTime))-1);
    sprintf(buf2,"Last event: %s",buf);
    fLastUpdate->SetText(buf2);
    if (lastSpillNr!=SpillNr) {lastSpillNr=SpillNr;lastSpillChange=time(0);}
    if (abs(time(0)-EvntTime)>60||abs(time(0)-lastSpillChange)>60) {
        count_of_no_data =(count_of_no_data+1) % 600;
        if (count_of_no_data==0) 
                system("ssh -f pccorc21 play /online/soft/DaqDataDecoding/examples/MurphyTV/somethingwrong.wav&");
        if (mycolor!=red) {mycolor=red;mycolor2=grey;} else
                    	  {mycolor=grey;mycolor2=red;}
     	fLastUpdate->ChangeBackground(mycolor);
        fRunNumber->ChangeBackground(mycolor);
        fSpillNumber->ChangeBackground(mycolor2);
    } else {
        count_of_no_data=0;    
        fLastUpdate->ChangeBackground(green);
        fRunNumber->ChangeBackground(green);
        fSpillNumber->ChangeBackground(green);
    }
    sprintf(buf,"Run number: %u",RunNr);
    fRunNumber->SetText(buf);
    
    sprintf(buf,"Spill number: %3.3d",SpillNr);
    fSpillNumber->SetText(buf);
    
	
  } else {
    sprintf(buf,"Last event: n.a.                                     ");
    fLastUpdate->SetText(buf);
    sprintf(buf,"Run number: n.a.        ");
    fRunNumber->SetText(buf);
    sprintf(buf,"Spill number: n.a. ");
    fSpillNumber->SetText(buf);
  } 
  
  if (RunNr!=trash->GetRunNr()) {
  	RunNr=trash->GetRunNr();
	catchtab->SortBy(2, false);
	porttab->SortBy(2,false);
        Update(true);
  }
}

void MurphyTV::CatchList::SelectSpecial(DaqErrorType type){
  fCombo->Select(type);
  special=type;
  catchtab->SortBy(6,false);
  porttab->SortBy(5,false);
  Update(true);
}

Bool_t MurphyTV::CatchList::ProcessMessage(Long_t msg,Long_t parm1,Long_t parm2){
  uint32 ns;
  switch(GET_MSG(msg)) {
    case kC_COMMAND:
      switch(GET_SUBMSG(msg)){
        case kCM_COMBOBOX:
          if (parm1==CB_SPECIAL) SelectSpecial((DaqErrorType)parm2);          
          break;
        case kCM_BUTTON:
          if (parm1==LB_CATCHES)
            SendMessage(MsgWindow,MK_MSG(kC_COMMAND,kCM_BUTTON), Id,parm2);
          break;
        case kCM_CHECKBUTTON:
          Update(true);
          break;
        case kCM_LISTBOX:
          switch(parm1){
            case LB_CATCHES:
              ns=parm2;
              if (ns>CollectErrs::MAX_SOURCE || asource==ns) break;              
              asource=ns;
              Update(false,true);
              SrcInfos(asource);
              break;
            default:
              break;
          }
          break;
        default:
	  break;
      }
      break;
    default:
      break;
  }
  return kTRUE;
}


/*****************************************************************************************************
*****************************************************************************************************/

MurphyTV::ErrsList::ErrsList(TGWindow *p,UInt_t w, UInt_t h,CollectErrs*ce,Int_t id, bool small)
:TGCompositeFrame(p, w, h, kVerticalFrame)
{
  MsgWindow=p;
  Id=id;
  trash=ce;

  asource=2;
  disperr=0xffffff;
  
  fCF = new TGCompositeFrame(this, 100, 25, kHorizontalFrame);
  fLH3 = new TGLayoutHints(kLHintsCenterX | kLHintsCenterY ,5, 5, 5, 5);
  fLH1 = new TGLayoutHints(kLHintsCenterY | kLHintsLeft ,5, 5, 5, 5);
  fLabel = new TGLabel(fCF,new TGHotString("SrcID="));
  fCF->AddFrame(fLabel,fLH1);
  TGTextBuffer *tbuf = new TGTextBuffer(10);
  tbuf->AddText(0, "2");
  fTESrc = new TGTextEntry(fCF, tbuf,TE_SRC);
  fTESrc->SetToolTipText("Show information about this SourceId in the \"-On SrcID\" collumns",1000);
  fTESrc->Resize(50, fTESrc->GetDefaultHeight());
  fTESrc->Associate(this);
  fCF->AddFrame(fTESrc,fLH1);
  fCBHide = new TGCheckButton(fCF,"Show only errors on SrcID");
  fCBHide->SetState(kButtonUp);
  fCBHide->SetToolTipText("Show only errors which occured on the selected SourceID",1000);
  fCBHide->Associate(this);
  fCF->AddFrame(fCBHide,fLH3);
  fCBNorm = new TGCheckButton(fCF,"Normalize");
  fCBNorm->SetState(kButtonDown);
  fCBNorm->SetToolTipText("Norm everything on the number of recorded events",1000);
  fCBNorm->Associate(this);
  fCF->AddFrame(fCBNorm,fLH3);

  fLH2 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX , 2, 2, 2, 2);  
  AddFrame(fCF,fLH2);
                             
  errtab=new SorTable(this,w,h,kSunkenFrame|kDoubleBorder,LB_ERRS, small);
  errtab->AddCollumn("Error Type",30,true);
  errtab->AddCollumn("Count",12,false,"How often the error occured at all");
  errtab->AddCollumn("-On SrcID",12,false,"How often the error occured on SourceID SrcID");
  errtab->AddCollumn("#BadEvents",12,false,"In how many events the error occured");
  errtab->AddCollumn("-On SrcID",12,false,"In how many events the error occured on SourceID SrcID");
  errtab->AddCollumn("#BadSrcIDs",12,false,"On how many sourceIDs the error occured");
  errtab->AddCollumn("#Spills",12,false,"Nb of spills where this SrcID is affected by this error");
  errtab->AddCollumn("Last spill",12,false,"Last spill where this SrcID is affected by this error");
  errtab->LastCollumn();
  errtab->SortBy(1,false);
  fLH = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2);
  AddFrame(errtab,fLH); 
  fGF= new TGGroupFrame(this,"Error Description :");
  
  fTVDescr=new TGTextView(fGF,300,80,-1,0,fGF->GetDefaultFrameBackground());
  fTVDescr->ChangeOptions(0);
  fGF->AddFrame(fTVDescr,fLH2);
  
  AddFrame(fGF,fLH2);
  
  Update(true);
}


MurphyTV::ErrsList::~ErrsList()
{
  delete errtab;
  delete fTESrc;
  delete fLabel;
  delete fCF;
  delete fTVDescr;
  delete fGF;
  delete fLH3;
  delete fLH2;
  delete fLH1;
  delete fLH;
}


void MurphyTV::ErrsList::Update(bool redraw){
  uint32 a,b,c,d;

  bool dnorm=( fCBNorm->GetState()==kButtonDown ); 

  errtab->ClearRows(); 
  trash->Lock();
  
  uint32 norm=trash->GetNrRecEvents();  
  for (a=b=0;a<ErrZoo::GetNrErrs();a++) {
    DaqErrorType et=ErrZoo::GetType(a);
    uint32 ne=trash->GetErr(et);
    if (ne && ( trash->GetErrOnCatch(asource,et) || fCBHide->GetState()==kButtonUp )) {
      errtab->AddRow(a);
      errtab->GetEntry(0,b)=SorTable::Entry(ErrZoo::GetTitle(a));
      errtab->GetEntry(1,b)=SorTable::Entry(ne,norm,dnorm);
      errtab->GetEntry(2,b)=SorTable::Entry(trash->GetErrOnCatch(asource,et),norm,dnorm);
      errtab->GetEntry(3,b)=SorTable::Entry(trash->GetNrEvntsWithErr(et),norm,dnorm);
      errtab->GetEntry(4,b)=SorTable::Entry(trash->GetEvntsWithErrOnCatch(asource,et),norm,dnorm);
      for (c=d=0;c<CollectErrs::MAX_SOURCE;c++) if (trash->GetErrOnCatch(c,et)) d++; 
      errtab->GetEntry(5,b)=SorTable::Entry((Int_t)d);     
      errtab->GetEntry(6,b)=SorTable::Entry(trash->GetNbSpillWithErr(asource,et),norm,kFALSE);
      errtab->GetEntry(7,b)=SorTable::Entry(trash->GetLastSpillWithErr(asource,et),norm,kFALSE);
      b++;
    }
  }

  trash->Unlock();  
  if (redraw) errtab->Update();
  Int_t sel=errtab->GetSelected();  
  if (disperr!=sel) {
    fTVDescr->LoadBuffer(ErrZoo::GetDescription(sel));
    disperr=sel;
  }
}

void MurphyTV::ErrsList::ChoseSource(uint32 src) {
  char buf[100];
  if (trash->CatchExists(src)) asource=src;
  sprintf(buf,"%d",asource);
  fTESrc->SetText(buf);
  errtab->SortBy(2,false);
  fCBHide->SetState(kButtonDown);
  Update();
}

Bool_t MurphyTV::ErrsList::ProcessMessage(Long_t msg,Long_t parm1,Long_t parm2){
  int a;
  switch(GET_MSG(msg)) {
    case kC_COMMAND:
      switch(GET_SUBMSG(msg)){
        case kCM_BUTTON:
          if (parm1==LB_ERRS)
            SendMessage(MsgWindow,MK_MSG(kC_COMMAND,kCM_BUTTON), Id,ErrZoo::GetType(parm2));
          break;
        case kCM_CHECKBUTTON:
          Update();
          break;
        case kCM_LISTBOX:
          switch(parm1){
            case LB_ERRS:
              fTVDescr->LoadBuffer(ErrZoo::GetDescription(parm2));
              disperr=parm2;
              break;
            default:
              break;
          }
          break;
        default:
	        break;
      }
      break;
    case kC_TEXTENTRY:
      if (GET_SUBMSG(msg)==kTE_ENTER){
        switch(parm1){
          case TE_SRC:            
            if (1==sscanf(fTESrc->GetText(),"%d",&a)) ChoseSource(a);
            else ChoseSource(asource);                       
            break;
          default:
            break;
        }  
      }
      break;
    default:
      break;
  }
  return kTRUE;
}


/*****************************************************************************************************
*****************************************************************************************************/

MurphyTV::HistMonitor::HistMonitor(TGWindow *p,UInt_t w,UInt_t h,CollectErrs *ce)
:TGCompositeFrame(p, w, h, kVerticalFrame)
{
  trash=ce;
  DRun=trash->GetRunCount();
  Ackno=true;
  
  fL1 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 0, 0, 2, 2);
  fL2 = new TGLayoutHints(kLHintsExpandY | kLHintsRight , 2, 0, 2, 2);
  fL4 = new TGLayoutHints(kLHintsCenterY | kLHintsExpandX , 0, 0, 2, 2);
  fL3 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 0, 0, 2, 2);

  fGF = new TGGroupFrame(this,"General Information",kHorizontalFrame);

  for (uint32 a=0;a<3;a++) {
    fCFc[a] = new TGCompositeFrame(fGF, 180, 200, kVerticalFrame);
    for (uint32 b=0;b<4;b++) {
      fLInfos[a][b]=new TGLabel(fCFc[a],"                               ");
      fLInfos[a][b]->SetTextJustify(kTextLeft);
      fCFc[a]->AddFrame(fLInfos[a][b],fL1);
    }    
    fGF->AddFrame(fCFc[a],fL3);
  }
  AddFrame(fGF,fL1);
  
  fLMsg= new TGLabel(this,"Messages:");
  fLMsg->SetTextJustify(kTextLeft);
  AddFrame(fLMsg,fL1);
    
  fTVMsg = new TGTextView(this,300,300);
  AddFrame(fTVMsg,fL3);
  fCF2= new TGCompositeFrame(this, 180, 200, kHorizontalFrame);
  fnt=gClient->GetFontByName("-adobe-courier-bold-r-normal--34-240-100-100-m-200-iso8859-1");
  tggc= new TGGC(TGLabel::GetDefaultGC());
  tggc->SetFont(gVirtualX->GetFontHandle(fnt));
  fLShit = new TGLabel(fCF2,"no news",tggc->GetGC(),fnt,kSunkenFrame);
  fCF2->AddFrame(fLShit,fL4);
  fTBAha = new TGTextButton(fCF2,"Alert off till next change",TB_AHA);
  fTBAha->SetToolTipText("Turn the alert on the left off till something changes again",1000);
  fTBAha->Associate(this);
  fCF2->AddFrame(fTBAha,fL2);
  AddFrame(fCF2,fL1);
  
  fLInfos[0][1]->SetText("Current Event time stamp :");
}


MurphyTV::HistMonitor::~HistMonitor()
{
  delete fLShit;
  delete tggc;
  delete fTBAha;
  delete fCF2;
  for (uint32 a=0;a<3;a++) {
    for (uint32 b=0;b<4;b++) delete fLInfos[a][b];
    delete fCFc[a];
  }
  
  delete fGF;
  delete fL1;
  delete fL2; 
  delete fL3; 
  delete fL4;
}

void MurphyTV::HistMonitor::Update()
{
  char buf[200];
  TGText nt;

  trash->Lock();
  
  uint32 RunNr=trash->GetRunNr();
 
  if (DRun!=trash->GetRunCount() ) {
    fTVMsg->Clear();
    sprintf(buf,"== New Run %d started ==============================================",RunNr);
    fTVMsg->AddLine(buf);
    fTVMsg->AddLine("   ");
    repstat.clear();
    DRun=trash->GetRunCount();
  }

  uint32 nre=trash->GetNrRecEvents();
  uint32 nem=0;
  if (nre>100) for (uint32 a=0;a<CollectErrs::MAX_SOURCE;a++) {
    bool ce[ErrZoo::GetNrErrs()];
    if (!trash->CatchExists(a)) continue;
    for (uint32 b=0;b<ErrZoo::GetNrErrs();b++) ce[b]=false;
    RepStat::iterator it=repstat.find(a);
    if (it !=repstat.end()) {
      for (vector <msg>::iterator er=it->second.begin();er!=it->second.end();er++) {
        ce[er->errind]=true;
        float se=(float)trash->GetEvntsWithErrOnCatch(a,ErrZoo::GetType(er->errind));
        float pe=100.0*se/(float)nre;
        se=pe/sqrt(se);
        if (er->pbad<10.0) {
          if (se<3.0) se=3.0;
        }
        else if (se<10.0) se=10.0;
        if ((pe>se + er->pbad ) || (pe+se<er->pbad )) {
          sprintf(buf,"On SourceID %4d : Error \"%s\" now in %5.2f %% of the events",a,ErrZoo::GetTitle(er->errind),pe);
          nt.InsLine(nt.RowCount(),buf);
          nem++;
          er->pbad=pe;
        }
      }
    }
    for (uint32 b=0;b<ErrZoo::GetNrErrs();b++){
      if (ce[b]) continue;
      float se=(float)trash->GetEvntsWithErrOnCatch(a,ErrZoo::GetType(b));
      float pe=100.0*se/(float)nre;
      se=pe/sqrt(se);
      if (pe>1.0+se) {
        sprintf(buf,"On SourceID %4d : new Error \"%s\" in %5.2f %% of the events",a,ErrZoo::GetTitle(b),pe);
        nt.InsLine(nt.RowCount(),buf);
        nem++;
        repstat[a].push_back(msg(b,pe));
      }
    } 
    if (nem>200) {
      nt.InsLine(nt.RowCount(),"***************************************************");
      nt.InsLine(nt.RowCount(),"*** Attention !! Too many new error messages !! ***");
      nt.InsLine(nt.RowCount(),"*** Only first ~200 are displayed.              ***");
      nt.InsLine(nt.RowCount(),"***************************************************");
      break;
    }
  } 

  time_t EvntTime=trash->GetEventTime();
  uint32 ErrSum=trash->GetErrSum();
  uint32 PSize=trash->GetParsedSize();
  uint32 NrSrc=0;
  uint32 NrGood=0;
  uint32 NrGone=0;
  for (uint32 a=0;a<CollectErrs::MAX_SOURCE;a++) if (trash->CatchExists(a)) {
    NrSrc++;
    if (trash->GetCatchErrSum(a)==0) NrGood++;
    if (trash->GetErrOnCatch(a,DaqError::EVENT_MISSING_SRCID)) NrGone++; 
  }   
  
  trash->Unlock();
  
  if (nem) {
    if (fTVMsg->GetText()->RowCount()+nem>260) {
      fTVMsg->Clear();
      sprintf(buf,"== New Run %d started ==============================================",RunNr);
      fTVMsg->AddLine(buf);
      fTVMsg->AddLine("( Too many messages. Old ones are deleted. )");
      fTVMsg->AddLine("----------------------------------------------------------------------------");
    }
    fTVMsg->AddText(&nt); 
    fLShit->SetText("You got new errors.");
    Ackno=false;
    ULong_t clr=0;
    gClient->GetColorByName("Red",clr);
    fLShit->ChangeBackground(clr);
  } 

  if (nre) {
    sprintf(buf,"RunNr : %d",RunNr);
    fLInfos[0][0]->SetText(buf);
    fLInfos[0][2]->SetText(ctime(&EvntTime));
    sprintf(buf,"Rec. events : %d",nre);
    fLInfos[1][0]->SetText(buf);
    sprintf(buf,"#Errors : %d",ErrSum);
    fLInfos[1][1]->SetText(buf);
    sprintf(buf,"per event : %f",(float)(ErrSum)/(float)(nre));
    fLInfos[1][2]->SetText(buf);  
    sprintf(buf,"Event size : %d",PSize/nre);
    fLInfos[1][3]->SetText(buf);
    sprintf(buf,"#SourceIDs : %d",NrSrc);
    fLInfos[2][0]->SetText(buf);
    sprintf(buf,"missing : %d",NrGone);
    fLInfos[2][1]->SetText(buf);
    sprintf(buf,"error-free : %d",NrGood);
    fLInfos[2][2]->SetText(buf);
  }
  else {
    fLInfos[0][0]->SetText("RunNr : n.a");
    fLInfos[0][2]->SetText("n.a");
    fLInfos[0][3]->SetText("Trigger# : n.a");
    fLInfos[1][0]->SetText("Rec. events : n.a");
    fLInfos[1][1]->SetText("#Errors : n.a");
    fLInfos[1][2]->SetText("per event : n.a");
    fLInfos[1][3]->SetText("Event size : n.a");  
    fLInfos[2][0]->SetText("#SourceIDs : n.a");
    fLInfos[2][1]->SetText("missing : n.a");
    fLInfos[2][2]->SetText("error-free : n.a");
  }  
  
}

Bool_t MurphyTV::HistMonitor::ProcessMessage(Long_t msg,Long_t parm1,Long_t parm2){
  switch(GET_MSG(msg)) {
    case kC_COMMAND:
      switch(GET_SUBMSG(msg)){
        case kCM_BUTTON:
          if (parm1==TB_AHA) {
            if (Ackno) break;
            fLShit->SetText("no news");
            fLShit->ChangeBackground(TGFrame::GetDefaultFrameBackground());
            fTVMsg->AddLine("----------------------------------------------------------------------------");
            Ackno=true;
          }
          break;
        default:
          break;
      }
      break;
    default:
      break;
  }
  return kTRUE;
}



/*****************************************************************************************************
*****************************************************************************************************/

MurphyTV::HistHist::HistHist(TGWindow *p,UInt_t w,UInt_t h,CollectErrs *ce)
:TGCompositeFrame(p, w, h, kVerticalFrame)
{
  trash=ce;
  errkat=-1;
  
  fL3 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2);

  gStyle->SetCanvasColor(kWhite);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.15);    
  gStyle->SetPadBottomMargin(0.15);
  
  fEc2 = new TRootEmbeddedCanvas("ec2", this, 300, 300);
  AddFrame(fEc2, fL3); 
  fEc2->GetCanvas()->SetBorderMode(0);
  fF2 = new TGCompositeFrame(this, 100, 30, kHorizontalFrame);  
  fL5 = new TGLayoutHints(kLHintsCenterY | kLHintsLeft ,5, 5, 5, 5);  
  fLabel2 = new TGLabel(fF2,"Display error:");
  fF2->AddFrame(fLabel2,fL5);
  fCombo2 = new TGComboBox(fF2,CB_ERRKAT2);
  fF2->AddFrame(fCombo2, fL5);
  fL4 = new TGLayoutHints(kLHintsBottom | kLHintsLeft | kLHintsExpandX , 2, 2, 2, 2);  
  AddFrame(fF2, fL4);
  fCombo2->AddEntry("All", -1); 
  for (uint32 a = 0; a < ErrZoo::GetNrErrs(); a++) fCombo2->AddEntry(ErrZoo::GetTitle(a), a);  

  fCombo2->Select(errkat);
  fCombo2->Resize(200, 20);
  fCombo2->Associate(this);

  fEc2->GetCanvas()->cd();
//  fEc2->GetCanvas()->SetLogy();
  fEc2->GetCanvas()->SetGridy(1);

  gStyle->SetNdivisions(0,"X");  
  fEc2->GetCanvas()->cd();
  fHist2=new TH1F("fHist2","",trash->GetHistBufSize(),0,trash->GetHistBufSize());
  fHist2->SetMaximum(1E2);
  fHist2->SetMinimum(1E-4);
  fHist2->Scale(1);
  fHist2->SetStats(kFALSE);
  fHist2->GetYaxis()->SetLabelFont(132);
  fHist2->GetYaxis()->SetLabelSize(0.05);
  fHist2->GetYaxis()->SetTitleFont(42);
  fHist2->GetYaxis()->SetTitleSize(0.05);
  fHist2->GetYaxis()->SetTitleOffset(1);

  fHist2->GetYaxis()->SetTitle("Errors per rec. event");
  fEc2->GetCanvas()->cd();
  fHist2->Draw();

  fXaxis2 = new TGaxis(0,0, trash->GetHistBufSize(), 0,0,
                       trash->GetHistBufSize()*trash->GetHistBufRes()/1E3,510,"-=");
  fXaxis2->SetTitle("Events * 10^3");
  fXaxis2->SetTitleOffset(1.25);
  fXaxis2->SetLabelFont(132);
  fXaxis2->SetLabelOffset(-0.075);
  fXaxis2->SetLabelSize(0.045);
  fXaxis2->SetTitleSize(0.045);
 
  fXaxis2->Draw();

  ChoseTEHist(errkat);
}

MurphyTV::HistHist::~HistHist()
{
  delete fHist2;
  delete fXaxis2;
  delete fEc2;
  delete fCombo2;
  delete fLabel2;
  delete fL3;
}

float MurphyTV::HistHist::GetErrKatRate(int bin,int kat){
  float erg=0;
  uint32 nre=trash->GetRecEventsHist(bin);
  if (nre==0) return 0;  
  if (kat>=0) {
    erg=trash->GetErrHist(ErrZoo::GetType(kat),bin);
  }
  else {
    for (uint32 a=0;a<ErrZoo::GetNrErrs();a++) erg+=trash->GetErrHist(ErrZoo::GetType(a),bin);       
  }  
  return erg/(float)nre;
}

void MurphyTV::HistHist::Update()
{
  for (uint32 b=0;b<trash->GetHistBufSize();b++) fHist2->SetBinContent(b+1,GetErrKatRate(b,errkat));  
  fEc2->GetCanvas()->Modified();
  fEc2->GetCanvas()->Update();
}

void MurphyTV::HistHist::ChoseTEHist(int kat) {
  char H2Title[200];
  errkat=kat;   
  if (errkat>=0) sprintf(H2Title,"History of %s error",ErrZoo::GetTitle(errkat)); 
  else sprintf(H2Title,"History of total error number"); 
  fHist2->SetNameTitle("fHist2",H2Title);
  Update();
}


Bool_t MurphyTV::HistHist::ProcessMessage(Long_t msg,Long_t parm1,Long_t parm2){
  switch(GET_MSG(msg)) {
    case kC_COMMAND:
      switch(GET_SUBMSG(msg)){
        case kCM_COMBOBOX:
	        switch(parm1) {
	          case CB_ERRKAT2:
              if ((parm2>=(Long_t)ErrZoo::GetNrErrs())||(parm2<-1)) {
                parm2=0;
                fCombo2->Select(parm2);
              }
              ChoseTEHist(parm2);
	            break;
	          default:
	            break;
	        }
          break;
        default:
	        break;
      }
      break;
    default:
      break;
  }
  return kTRUE;
}


/*****************************************************************************************************
*****************************************************************************************************/

MurphyTV::HistSizes::HistSizes(TGWindow *p,UInt_t w,UInt_t h,CollectErrs *ce)
:TGCompositeFrame(p, w, h, kVerticalFrame)
{

  trash=ce;

  SourceID1=0x02;
  
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.1);    
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetOptStat(1001100);
  
  fEc1 = new TRootEmbeddedCanvas("ec4", this, 300, 300);
  fL3 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2);
  AddFrame(fEc1,fL3);
  //fEc1->GetCanvas()->SetEditable(kFALSE);
  fEc1->GetCanvas()->SetBorderMode(0);

  fF1 = new TGCompositeFrame(this, 100, 30, kHorizontalFrame);
  
  fL5 = new TGLayoutHints(kLHintsCenterY | kLHintsLeft ,5, 5, 5, 5);
  fLabel1 = new TGLabel(fF1,new TGHotString("SourceID :"));
  fF1->AddFrame(fLabel1,fL5);
  TGTextBuffer *tbuf = new TGTextBuffer(10);
  tbuf->AddText(0, "2");
  fEdit1 = new TGTextEntry(fF1, tbuf,TE_CATCHSEL1);
  fEdit1->Resize(50, fEdit1->GetDefaultHeight());
  fEdit1->Associate(this);
  fF1->AddFrame(fEdit1,fL5);
  
  fL4 = new TGLayoutHints(kLHintsBottom | kLHintsLeft | kLHintsExpandX , 2, 2, 2, 2);  
  AddFrame(fF1, fL4);


  fEc1->GetCanvas()->cd();
  fEc1->GetCanvas()->SetLogy();
  fEc1->GetCanvas()->SetGridy(1);

  gStyle->SetNdivisions(510,"X");
  fEc1->GetCanvas()->cd();
  fHist1=new TH1F("fHist4","Event size Distribution",trash->GetSizeHistSize(),
                  0,trash->GetSizeHistRes()*trash->GetSizeHistSize());

  fHist1->SetMaximum(1E5);
  fHist1->SetMinimum(1);
  fHist1->SetStats(kTRUE);

  fHist1->GetYaxis()->SetTitle("#Events");
  fHist1->GetXaxis()->SetTitle("Event Size/Byte");

  fEc1->GetCanvas()->cd();
  fHist1->Draw("");
}


MurphyTV::HistSizes::~HistSizes()
{
  delete fHist1;
  delete fEc1;
  delete fL3;
  delete fL4;
  delete fL5;
  delete fEdit1;
  delete fLabel1;  
}


void MurphyTV::HistSizes::Update()
{
  uint32 a;  
  for (a=0;a<trash->GetSizeHistSize();a++) fHist1->SetBinContent(a,trash->GetSizeHistOfCatch(SourceID1,a));
  fEc1->GetCanvas()->Modified();
  fEc1->GetCanvas()->Update();
}

Bool_t MurphyTV::HistSizes::ProcessMessage(Long_t msg,Long_t parm1,Long_t parm2){
  int a;
  char buf[20];

  switch(GET_MSG(msg)) {
    case kC_TEXTENTRY:
      if (GET_SUBMSG(msg)==kTE_ENTER){
        switch(parm1){
          case TE_CATCHSEL1:
            if (1==sscanf(fEdit1->GetText(),"%d",&a)) if (trash->CatchExists(a)) SourceID1=a;
            sprintf(buf,"%d",SourceID1);
            fEdit1->SetText(buf);                        
            break;
          default:
            break;
        }  
      }
      break;
    default:
      break;
  }
  return kTRUE;
}


/*****************************************************************************************************
*****************************************************************************************************/

MurphyTV::HistHistD::HistHistD(TGWindow *p,UInt_t w,UInt_t h,CollectErrs *ce)
:TGCompositeFrame(p, w, h, kVerticalFrame)
{
  trash=ce;
  
  fL3 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2);

  gStyle->SetCanvasColor(kWhite);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.15);    
  gStyle->SetPadBottomMargin(0.15);
  
  fEc2 = new TRootEmbeddedCanvas("ec2", this, 300, 300);
  AddFrame(fEc2, fL3);  
  
  fEc2->GetCanvas()->SetBorderMode(0);
  fEc2->GetCanvas()->cd();
 // fEc2->GetCanvas()->SetLogy();
  fEc2->GetCanvas()->SetGridy(1);

  gStyle->SetNdivisions(0,"X");  
  fEc2->GetCanvas()->cd();
  fHist2=new TH1F("fHist2D","History of event size",trash->GetHistBufSize(),0,trash->GetHistBufSize());
  fHist2->SetMaximum(1E5);
  fHist2->SetMinimum(0);
  fHist2->Scale(1);
  fHist2->SetStats(kFALSE);
  fHist2->GetYaxis()->SetLabelFont(132);
  fHist2->GetYaxis()->SetLabelSize(0.05);
  fHist2->GetYaxis()->SetTitleFont(42);
  fHist2->GetYaxis()->SetTitleSize(0.05);
  fHist2->GetYaxis()->SetTitleOffset(1.5);
  fHist2->GetYaxis()->SetTitle("Byte per rec. event");
  fEc2->GetCanvas()->cd();
  fHist2->Draw();

  fXaxis2 = new TGaxis(0,0, trash->GetHistBufSize(), 0,0,
                       trash->GetHistBufSize()*trash->GetHistBufRes()/1E3,510,"-=");
  fXaxis2->SetTitle("Events * 10^3");
  fXaxis2->SetTitleOffset(1.25);
  fXaxis2->SetLabelFont(132);
  fXaxis2->SetLabelOffset(-0.075);
  fXaxis2->SetLabelSize(0.045);
  fXaxis2->SetTitleSize(0.045);
 
  fXaxis2->Draw();
}

MurphyTV::HistHistD::~HistHistD()
{
  delete fHist2;
  delete fXaxis2;
  delete fEc2;
  delete fL3;
}


void MurphyTV::HistHistD::Update()
{
  for (uint32 b=0;b<trash->GetHistBufSize();b++)  {   
    if (trash->GetRecEventsHist(b))
      fHist2->SetBinContent(b+1,(float)trash->GetEventSizeHist(b)/(float)trash->GetRecEventsHist(b));
    else fHist2->SetBinContent(b+1,0);
  }
  fEc2->GetCanvas()->Modified();
  fEc2->GetCanvas()->Update();

}

/*****************************************************************************************************
*****************************************************************************************************/

HelpDlg::HelpDlg(TGWindow *main,HelpDlg **me,char *title,char *txt):TRootHelpDialog(main,title,800,400) {
  Me=me;
  SetText(txt);
  Popup();
}

  

