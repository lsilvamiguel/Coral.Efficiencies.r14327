
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <pthread.h>
#include <time.h>

#include "MurphyTV.h"
#include "monitor.h"

using CS::uint32;
using CS::DaqError;
using CS::DaqOption;
using CS::ObjectXML;
using CS::Chip;
using namespace std;

/*****************************************************************
 class to execute function DoSomething() in a background thread
 DoSomething is to be rewritten in child classes.
*****************************************************************/

class BasicThread {  
  public:
    enum tStatus {RUNNING,TERMINATING,STOPPED};
    
    BasicThread(bool);
    virtual ~BasicThread();
    void RunThread();
    void StopThread();
    virtual void DoSomething();
    tStatus GetStatus() {return status;}
    void Done() {status=TERMINATING;}
  private:
    pthread_t TheThread;
    tStatus status;
    static void *TheThreadFunc(BasicThread *that);
};

BasicThread::BasicThread(bool run){
	status=STOPPED;
  if (run) RunThread();
}

BasicThread::~BasicThread(){
  StopThread();
}

void BasicThread::RunThread(){
  if (status==STOPPED) {
    if (0==pthread_create( &TheThread, NULL, (void *(*)(void *))&TheThreadFunc,(void *) this)) {
      status=RUNNING;
    }
  } 
}

void BasicThread::StopThread(){
  if (status!=STOPPED) {
    status=TERMINATING;
    pthread_join ( TheThread, NULL );
    status=STOPPED;
  }
}

void BasicThread::DoSomething(){}
void *BasicThread::TheThreadFunc(BasicThread *that){
  while (TERMINATING!=that->GetStatus()) that->DoSomething();    
  return 0;
}


/*****************************************************************
class to get events from DATE and hand them to CollectErrs
in background thread
*****************************************************************/

class ErrThread:public BasicThread{
  public:
    ErrThread(char *mfile=0);
    ~ErrThread();
    void DoSomething();
    bool Connect(const char *source);
    void Disconnect();
    CollectErrs *GetTrash() {return trash;}
    bool IsConnected() {return connected;}
    void CalibrationTrigger(bool ct) { calibrationTrigger=ct; }
  private:
    CollectErrs *trash;
    bool connected;
    bool calibrationTrigger;
};

ErrThread::ErrThread(char *mfile):BasicThread(false){
  trash=new CollectErrs(mfile); 
}

ErrThread::~ErrThread(){
  StopThread();
   delete trash;
}

void ErrThread::DoSomething(){ 
  int st;
  char *ptr;
  DaqEvent *evnt;
	
  st=monitorGetEventDynamic((void **)&ptr); //get raw event from DATE 
  if (0==st) {
    connected=true;
    if (0==ptr) {
    	usleep(100000); // 100 msec. sleep, prevent too high CPU usage when idle
	return;}  
    evnt=new DaqEvent(ptr);  //prepare event for use with decoding library

    trash->DecodeEvent(evnt);    //decode event and scan for errors
    trash->Lock();               //prevent collision with other thread (GUI)
    trash->HandleNewEvent(evnt); //update error counters 
    trash->Unlock();
    delete evnt;
    free(ptr);   
  } else { 
    connected=false;
    printf("Date error:%s \n",monitorDecodeError(st));
    Done();
  }
}

bool ErrThread::Connect(const char *source=0){
  int err=0;
  static char *table[6]={"All","no","PHY","yes",0};
  static char *tableCalibrationTrigger[6]={"All","no","CAL","yes",0};
  static char datasource[1000]="";  
  
  if (source) {
      if (strcmp(source,"auto")==0 && strncmp(datasource,"@",1)!=0) return false; 
      if (strcmp(source,"auto")!=0) strncpy(datasource,source,1000); 
  }
  StopThread();
  connected=false;
    if (strlen(datasource)) {
      err=monitorSetDataSource(datasource);   //tell DATE where to get the events
      if (err) {printf("Date error:%s \n",monitorDecodeError(err)); return false;}
    } else return false;
  
  monitorSetNowait();
  err=monitorDeclareMp("MurphyTV");
  if (err) {printf("Date error:%s \n",monitorDecodeError(err)); return false;}
  
  if (calibrationTrigger) 
      //We only want calibration triggers	
      err=monitorDeclareTable(&tableCalibrationTrigger[0]);
  else  
      //We only want physics events
      err=monitorDeclareTable(&table[0]);
  
  if (err) {printf("Date error:%s \n",monitorDecodeError(err)); return false;}
      
  //do a reset after all connects if not in autoconnect mode  
  if (strcmp(source,"auto")!=0) trash->ResetAtNextEvent(); 
  RunThread();       
  return true;
}

void ErrThread::Disconnect() {
  StopThread();
  monitorLogout();  //say goodbye to DATE
  connected=false;
}


/*****************************************************************************
this class provides the little control window, using ROOT GUI classes
*****************************************************************************/

class MainFrame : public TGMainFrame {
  enum CommandIdentifiers {
     TE_SOURCE,TE_MAPS ,TE_CATCH,
     TB_EXIT,TB_CON,TB_DISCON,TB_FILE,TB_DISPLAY,TB_HELP,TB_LOAD,TB_LOADMAPS
  };
private:
  ErrThread *thread;
  MurphyTV *myfr;      

  HelpDlg *HelpGuiDlg;
  TGTextEntry *fEdit1,*fTEMaps;
  TGLabel *fLabel1,*fLMaps;
  TGCompositeFrame *fCF,*fCF2,*fCFMaps;
  TGTextButton *fTB1,*fTB2,*fTB3,*fTB4,*fTB5,*fTB6,*fTBMaps;
  TGStatusBar *fSB;
  TGLayoutHints *fLH1,*fLH2,*fLH3,*fLH4, *fLH5;
  TTimer *fTimer;

  uint32 lastrec;
  TTime lasttime;
  TTime lastDataArrivalTime;
  static char *HelpGuiTxt;

public:
  MainFrame(bool display=false,char *dfile=0,char *mfile=0,
            bool calibrationTriggerOnly=false, bool smallfont=false);
  virtual ~MainFrame();
  void Load();
  bool small;
  bool calibrationTrigger;
  void LoadMaps();
  Bool_t HandleTimer(TTimer *timer);
  virtual void CloseWindow();
  virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t);
};

#include "helpgui.h"

MainFrame::MainFrame(bool display,char *dfile,char *mfile,
                     bool calibrationTriggerOnly, bool smallfont) 
  :TGMainFrame(gClient->GetRoot(), 401, 151,kVerticalFrame)
{
  small=smallfont;
  calibrationTrigger=calibrationTriggerOnly;
  thread=new ErrThread(mfile);
  thread->CalibrationTrigger(calibrationTrigger); 
  myfr=0;
  HelpGuiDlg=0;

  fLH2 = new TGLayoutHints(kLHintsBottom | kLHintsExpandX  ,1, 1, 1, 1);
  fLH1 = new TGLayoutHints(kLHintsExpandX |kLHintsExpandY |kLHintsTop ,2, 2, 2, 2);  
  fLH3 = new TGLayoutHints(kLHintsCenterY | kLHintsLeft,4,2,2,2);
  fLH4 = new TGLayoutHints(kLHintsCenterY | kLHintsExpandX,2,2,2,4);

  fSB=new TGStatusBar(this,300,30);
  Int_t parts[2]={20,80};
  fSB->SetParts(parts,2);
  AddFrame(fSB,fLH2);
      
  fCF = new TGCompositeFrame(this, 300, 25,kHorizontalFrame );
   
  fTB1=new TGTextButton(fCF,new TGHotString("E&xit"),TB_EXIT);
  fCF->AddFrame(fTB1,fLH1);
  fTB1->Associate(this);
  fTB2=new TGTextButton(fCF,new TGHotString("&Connect"),TB_CON);
  fCF->AddFrame(fTB2,fLH1);
  fTB2->Associate(this);
  fTB3=new TGTextButton(fCF,new TGHotString("Di&sconnect"),TB_DISCON);
  fCF->AddFrame(fTB3,fLH1);  
  fTB3->Associate(this);
  fTB4=new TGTextButton(fCF,new TGHotString("&Display"),TB_DISPLAY);
  fCF->AddFrame(fTB4,fLH1);
  fTB4->Associate(this);
  fTB6=new TGTextButton(fCF,new TGHotString("&Help"),TB_HELP);
  fCF->AddFrame(fTB6,fLH1);
  fTB6->Associate(this);
  AddFrame(fCF,fLH1);
  
  fCF2 = new TGCompositeFrame(this, 300, 30, kHorizontalFrame);
  
  fLabel1 = new TGLabel(fCF2,new TGHotString("Data source: "));
  fCF2->AddFrame(fLabel1,fLH3);

  TGTextBuffer *tbuf = new TGTextBuffer(999);
  if (dfile) tbuf->AddText(0,dfile);
  fEdit1 = new TGTextEntry(fCF2, tbuf,TE_SOURCE);
  fEdit1->Resize(150, fEdit1->GetDefaultHeight());
  fEdit1->Associate(this);
  fCF2->AddFrame(fEdit1,fLH4);
  fTB5 = new TGTextButton(fCF2,"...",TB_LOAD);
  fTB5->Associate(this);
  fCF2->AddFrame(fTB5,fLH3);
  
  AddFrame(fCF2,fLH2);
  
  fCFMaps = new TGCompositeFrame(this, 300, 30,kHorizontalFrame );
  
  fLMaps = new TGLabel(fCFMaps,new TGHotString("Mapping file: "));
  fCFMaps->AddFrame(fLMaps,fLH3);

  TGTextBuffer *tbuf2 = new TGTextBuffer(999);
  if (mfile) tbuf2->AddText(0,mfile);
  fTEMaps = new TGTextEntry(fCFMaps, tbuf2, TE_MAPS);
  fTEMaps->Resize(150, fEdit1->GetDefaultHeight());
  fTEMaps->Associate(this);
  fCFMaps->AddFrame(fTEMaps,fLH4);
  fTBMaps = new TGTextButton(fCFMaps,"...",TB_LOADMAPS);
  fTBMaps->Associate(this);
  fCFMaps->AddFrame(fTBMaps,fLH3);
  
  AddFrame(fCFMaps,fLH2);  
  
  SetWindowName("MurphyTV");

  MapSubwindows();

  Resize(450, 140);
  MapWindow();
  
  fTimer=new TTimer(500);
  fTimer->SetObject(this);
  lastrec=thread->GetTrash()->GetNrRecEvents();
  lasttime=gSystem->Now();
  fTimer->TurnOn(); 
  if (dfile) thread->Connect(dfile);
  if (display) myfr=new MurphyTV(thread->GetTrash(),this, calibrationTrigger, small);
  lastDataArrivalTime=gSystem->Now();
}

 MainFrame::~MainFrame()
{
  fTimer->TurnOff();
  delete fTimer;
  
  delete fSB;
  delete fTB1;
  delete fTB2;
  delete fTB3;
  delete fTB4;
  delete fTB5;
  delete fTB6;
  delete fEdit1;
  delete fLabel1;
  delete fCF;
  delete fCF2;
  
  delete fTBMaps;
  delete fTEMaps;
  delete fLMaps;
  delete fCFMaps;
  
  delete fLH1;
  delete fLH2;
  delete fLH3;
  delete fLH4;
  delete HelpGuiDlg;
	delete thread;
  if (MurphyTV::IsOpen()) delete myfr;
}

Bool_t MainFrame::HandleTimer(TTimer *timer)
{
  const char *st;
  char buf[50];
  uint32 a;
  long dt;
  long timeWithoutData;
  TTime now;
  
  if (thread->IsConnected()) st="Connected";
  else st="Not connected";
  fSB->SetText(st,0); 
  a=thread->GetTrash()->GetNrRecEvents();
  now=gSystem->Now();
  dt=long(now-lasttime);
  if (dt>0) sprintf(buf,"Got %d events (%4ld per sec)",a,(1000*(a-lastrec))/dt);
  else sprintf(buf,"???");
  fSB->SetText(buf,1);
  if ((a-lastrec)>0) lastDataArrivalTime=now;
  timeWithoutData=now-lastDataArrivalTime;  
  lastrec=a;
  lasttime=now;
  if (timeWithoutData>120000) { // reconnect if online data and no data for 2 minute
  		lastDataArrivalTime=now;
  		thread->Connect("auto");  
		}
  return kTRUE;
}

void MainFrame::Load(){
  const char *filetypes[] = {"DATE files","*.raw","DATE files","*.dat","All files","*",0,0 };
  TGFileInfo fi;
  fi.fFileTypes = filetypes;
  new TGFileDialog(gClient->GetRoot(), gClient->GetRoot(), kFDOpen,&fi);
  if (fi.fFilename) {
    fEdit1->SetText(fi.fFilename);
    thread->Connect(fi.fFilename);
  }
}
void MainFrame::LoadMaps(){
  const char *filetypes[] = {"Map files","*.xml","All files","*",0,0 };
  TGFileInfo fi;
  fi.fFileTypes = filetypes;
  new TGFileDialog(gClient->GetRoot(), gClient->GetRoot(), kFDOpen,&fi);
  if (fi.fFilename) {
    fTEMaps->SetText(fi.fFilename);
    thread->Disconnect();
    thread->GetTrash()->ReadMaps(fi.fFilename);
  }
}
void MainFrame::CloseWindow()
{
   TGMainFrame::CloseWindow();
   delete this;
   gApplication->Terminate(0);
}

Bool_t MainFrame::ProcessMessage(Long_t msg, Long_t parm1, Long_t)
{
  switch (GET_MSG(msg)) {

    case kC_COMMAND:
      switch (GET_SUBMSG(msg)) {
        case kCM_BUTTON:
          switch (parm1) {
            case TB_EXIT:
	      thread->Disconnect();
	      exit(0); 
              break;
            case TB_DISPLAY:
              if (!MurphyTV::IsOpen()) myfr=new MurphyTV(thread->GetTrash(),this, calibrationTrigger, small);
              break;
            case TB_CON:
	      thread->Disconnect();
	      thread->GetTrash()->ReadMaps();
              thread->Connect(fEdit1->GetText());        
              break;
            case TB_DISCON:
              thread->Disconnect();
              break; 
            case TB_HELP:
              if (0==HelpGuiDlg) HelpGuiDlg= new HelpDlg(this,&HelpGuiDlg,"Help on MurphyTV",HelpGuiTxt);
              break;  
            case TB_LOAD:
              Load();
              break;   
            case TB_LOADMAPS:
              LoadMaps();
              break;  
            default:
              break;              
          }
          break;
        default:
          break;
      }
      
    case kC_TEXTENTRY:
      if (GET_SUBMSG(msg)==kTE_ENTER){
        switch(parm1){
          case TE_SOURCE:
            thread->Connect(fEdit1->GetText());        
            break;
          case TE_MAPS:
	    thread->Disconnect();
            thread->GetTrash()->ReadMaps(fTEMaps->GetText());
            break;
	    case TE_CATCH:
	    thread->Disconnect();
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





extern void InitGui();
VoidFuncPtr_t initfuncs[]={InitGui,0};
TROOT root("GUI","gte",initfuncs);
int main(int argc,char **argv){ 
  char *dfile=0;
  char *mfile=0;
  bool display=false;
  bool small=false;
  bool calibrationTrigger=false;
 
  bool se=false;
  for (int a=1;a<argc;a++) {
    size_t al=strlen(argv[a]);
    char *ra=argv[a];
    char b='s';
    
    if (argv[a][0]=='-') {
      if (al<2) {se=true;break;}
      if (argv[a][1]=='d') {
        if (al==2) display=true; else {se=true;break;}
        continue;
      }
      if (argv[a][1]=='C') {
        if (al==2) calibrationTrigger=true;     
	else {se=true;break;}
        continue;
      }
      if (argv[a][1]=='f') {
        if (al==2) small=true;     
	else {se=true;break;}
        continue;
      }
      b=argv[a][1];
      if (al<=2) {
        a++;
	if (a>=argc) {se=true;break;}
        ra=argv[a];
      } else ra=argv[a]+2;
    }      
   
    if (b=='s') dfile=ra; 
    else if (b=='m') mfile=ra;
    else {se=true;break;}
  }
  
  if (se) { 
    cout << "\nusage: " 
         << "MurphyTV \n" 
         << "  [-s <data source>] \n" 
	 << "  [-m <mapping directory>]\n"
	 << "  [-d open display]\n"
	 << "  [-C analyse only calibration events]\n"
	 << "  [-f use small fonts]\n"
	 << endl;
    return 0;
  }

 
  Chip::SetFillStatistics(false);
  char *d1[2]={"MurphyTV",0};
  int d2=1;
  TApplication app("MurphyTV",&d2,d1);

  new MainFrame(display,dfile,mfile, calibrationTrigger, small);
  app.Run();
  
}
