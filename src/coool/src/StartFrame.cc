#include "StartFrame.h"
#include "config.h"
#include <iostream>
#include "compass_xpm.h"

const char* StartFrame::maptype[] = { "map files",   "*.xml",
				      "All files",   "*.*",
				       0,               0 };
const char* StartFrame::datatype[] = { "data files",   "*.dat",
				       "data files",   "*.raw",
				       "All files",   "*.*",
		 		       0,               0 };
const char* StartFrame::geomtype[] = { "data files",   "*.dat",
				       "All files",   "*.*",
				       0,               0 };


ClassImp(StartFrame)

namespace { // local routine
void usage(char* pgm_name, const char* err_string) {
  std::cerr << "\n\nWarning !! the arguments of the monitoring program have changed !!\n\n";
  if (err_string) std::cerr << err_string <<std::endl;
  std::cerr << "usage: " << gApplication->Argv(0) 
       << " [-map map_file] [-root root_file]"
       << " [-geom detector_file]"
       << " [-group groups_file] [-cfg user_configuration_file]" 
       << " [-workenv working environment string]"
       << " [data file or DATE source]\n";
  exit(2);
}
}


StartFrame::StartFrame(const TGWindow *p,MainFrame *mainf, UInt_t w, 
		       UInt_t h) 
  : TObject(), fMainframe(mainf) {

  fRootFile = "";
  fGroupFile = "/afs/cern.ch/compass/detector/monitor/groups.xml";

//  fGeomFile = "/afs/cern.ch/compass/detector/geometry/2002/detectors.22018.minus.dat";
  fGeomFile = "/afs/cern.ch/compass/detector/geometry/2004/detectors.for.coool.dat";
//   fMapFile = "/afs/cern.ch/compass/detector/maps/2003.xml";
  fMapFile = "/afs/cern.ch/compass/detector/maps";
  //   fDataFile = "/afs/cern.ch/compass/scratch/d12/data/cdr11002-18201.dat";
  fDataFile = "@pccore12:";
  fParamFile = "/afs/cern.ch/compass/detector/monitor/default_params";
  fConfigFile = "";
  std::string wkenv = "AFS";

  int nbarg = mainf->Argc();
  register int pargc = 1;

  while (pargc < nbarg) {
    char *pargv, *pargvn;
    pargv =  mainf->Argv(pargc);
    pargvn = 0;
    if ((pargc + 1) < nbarg) pargvn =  mainf->Argv(pargc+1);
    
    if (pargv[0] == '-') {   // it is an option
      if (strcmp(pargv, "-map") == 0) {  // option to give Map file
        if (pargvn && (pargvn[0] != '-')) fMapFile = pargvn;
          else usage(mainf->Argv(0), "badly formed -m option");
        pargc+=2;
        continue;
      }
      if (strcmp(pargv, "-geom") == 0) {  // option to give Geom file
        if (pargvn && (pargvn[0] != '-')) fGeomFile = pargvn;
          else usage(mainf->Argv(0), "badly formed -geom option");
        pargc+=2;
        continue;
      }
      if (strcmp(pargv, "-root") == 0) {  // option to give root file
        if (pargvn)
	  if (pargvn[0] != '-') {
	    fRootFile = pargvn;
	    pargc+=2;
	  }
	  else usage(mainf->Argv(0), "badly formed -root option, root file beginning by -");
	else usage(mainf->Argv(0), "badly formed -root option");
        continue;
      }
      if (strcmp(pargv, "-group") == 0) {  // option to give groups file
        if (pargvn && (pargvn[0] != '-')) fGroupFile = pargvn;
          else usage(mainf->Argv(0), "badly formed -group option");
        pargc+=2;
        continue;
      }
      if (strcmp(pargv, "-cfg") == 0) {  // option to give user config file
        if (pargvn && (pargvn[0] != '-')) fConfigFile = pargvn;
          else usage(mainf->Argv(0), "badly formed -cfg option");
        pargc+=2;
        continue;
      }
      if (strcmp(pargv, "-workenv") == 0) {  // option to give working environment string
        if (pargvn && (pargvn[0] != '-')) wkenv = pargvn;
          else usage(mainf->Argv(0), "badly formed -workenv option");
        pargc+=2;
        continue;
      }
      if (strcmp(pargv, "--help") == 0) {  // option to request help
        usage(mainf->Argv(0), 0);
      }
      // unknown option found
      char stmp[200];
      sprintf(stmp, "unknown option %s", pargv);
      usage(mainf->Argv(0), stmp);
    }
    else {   // data file name
      if ((pargc + 1) == nbarg) fDataFile = pargv;
      else usage(mainf->Argv(0),
                 "data file name not at the end of command line or more than one data file name");
      pargc++;
      continue;
    }
  }

//   cerr << "fMapFile " << fMapFile << " fRootFile " << fRootFile << " fGroupFile " << fGroupFile;
//   cerr << " fDataFile " << fDataFile <<endl;

  fMonitor=fMainframe->GetMonitor();
  fMonitor->SetWkEnvStr(wkenv.c_str());

  fMain = new TGTransientFrame(p, mainf->GetFrame(), w, h);
  fMain->Connect("CloseWindow()", "StartFrame", this, 
  		 "CloseWindow()");

  fL1 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandY | kLHintsExpandX, 0,0,0,0);
  fL2 = new TGLayoutHints(kLHintsExpandX | kLHintsBottom,
			  2, 2, 2, 2);
  fL3 = new TGLayoutHints(kLHintsBottom | kLHintsLeft,
			  2, 2, 2, 2);
  fL4 = new TGLayoutHints( kLHintsBottom | kLHintsRight,
			  2, 1, 1, 1);
  fL5 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 0,0,0,0);
  
  // Upper frame (Monitored detectors & Tasks)
  fUpFrame = new TGCompositeFrame(fMain, 60, 20, kHorizontalFrame);

  // Check frame
  fCheckFrame = new TGGroupFrame(fUpFrame, "Monitored detectors", kHorizontalFrame);
  fCheckFrame->SetLayoutManager(new TGHorizontalLayout(fCheckFrame));
  
  fMonDetFrame = new TGVerticalFrame(fCheckFrame, 60, 40, kVerticalFrame);
  std::map<std::string,bool>& detintree=fMonitor->GetDetInTree();
  typedef std::map<std::string,bool>::iterator IM;
  int n=0;
  for(IM i=detintree.begin(); i!=detintree.end(); i++) {
    TGCheckButton *check=new TGCheckButton(fMonDetFrame,
					   (i->first).c_str(),n);
    if(i->second)
      check->SetState(kButtonDown);
    else
      check->SetState(kButtonUp);
    
    fCheckButton.push_back(check);
    check->Connect("Clicked()","StartFrame",this,"Checked()");
    fMonDetFrame->AddFrame(check);
    n++;
  }
  fCheckFrame->AddFrame(fMonDetFrame,fL5);
  fMonDetFrame->Resize(75, (unsigned int) n*fCheckButton[0]->GetDefaultHeight());

  fAllNoneFrame = new TGVerticalFrame(fCheckFrame, 60, 20, kFixedWidth);

  fDefButton = new TGTextButton(fAllNoneFrame," &Default ",1);
  fDefButton->Connect("Clicked()","StartFrame",this,"DefClicked()");
  fAllNoneFrame->AddFrame(fDefButton,fL2);  

  fAllButton = new TGTextButton(fAllNoneFrame," &All ",1);
  fAllButton->Connect("Clicked()","StartFrame",this,"AllClicked()");
  fAllNoneFrame->AddFrame(fAllButton,fL2);  

  fNoneButton = new TGTextButton(fAllNoneFrame,"&None",1);
  fNoneButton->Connect("Clicked()","StartFrame",this,"NoneClicked()");
  fAllNoneFrame->AddFrame(fNoneButton,fL2);
  
  fCheckFrame->AddFrame(fAllNoneFrame,fL4);
  fAllNoneFrame->Resize(75, (unsigned int) 1.5*fAllButton->GetDefaultHeight());


  // Tasks frame
  fTaskFrame = new TGGroupFrame(fUpFrame, "Tasks", kVerticalFrame); 
  fTaskFrame->SetLayoutManager(new TGVerticalLayout(fTaskFrame));

  fCheckTree = new TGCheckButton(fTaskFrame,"Fill Tree",n);
  fCheckTree->Connect("Clicked()","StartFrame",this,"CheckedTree()");
  fCheckTree->SetToolTipText("fill a tree with the Planes variables in the root file");
  fTaskFrame->AddFrame(fCheckTree);

  fCheckText = new TGCheckButton(fTaskFrame,"Text Output",n);
  fCheckText->Connect("Clicked()","StartFrame",this,"CheckedText()");
  fCheckText->SetToolTipText("create at the end a text file with some numerical values computed from histograms");
  fTaskFrame->AddFrame(fCheckText);

  fCheckTrack = new TGCheckButton(fTaskFrame,"Tracking",n);
  fCheckTrack->Connect("Clicked()","StartFrame",this,"CheckedTrack()");
  fTaskFrame->AddFrame(fCheckTrack);
  
  fCheckCluster = new TGCheckButton(fTaskFrame,"Clustering",n);
  fCheckCluster->Connect("Clicked()","StartFrame",this,"CheckedCluster()");
  fCheckCluster->SetState(kButtonDown);
  fTaskFrame->AddFrame(fCheckCluster);
    
  fCheckRef = new TGCheckButton(fTaskFrame,"Show ref. plots",n);
  fCheckRef->Connect("Clicked()","StartFrame",this,"CheckedRef()");
  fCheckRef->SetToolTipText("show reference histos in red in some histograms");
  fCheckRef->SetState(kButtonDown);
  fTaskFrame->AddFrame(fCheckRef);

  fCheckCalib = new TGCheckButton(fTaskFrame,"Use Calibrations",n);
  fCheckCalib->Connect("Clicked()","StartFrame",this,"CheckedCalib()");
  fCheckCalib->SetToolTipText("read and use the calibration data");
  fCheckCalib->SetState(kButtonDown);
  fTaskFrame->AddFrame(fCheckCalib);

  fCheckReadFailedEvt = new TGCheckButton(fTaskFrame,"Read evt with decoding partly failed",n);
  fCheckReadFailedEvt->Connect("Clicked()","StartFrame",this,"CheckedReadFailedEvt()");
  fCheckReadFailedEvt->SetToolTipText("this happens for example for some calibration events with no trigger time");
  fCheckReadFailedEvt->SetState(kButtonDown);
  fTaskFrame->AddFrame(fCheckReadFailedEvt);

  fCheckExpert = new TGCheckButton(fTaskFrame,"Activate expert histos",n);
  fCheckExpert->Connect("Clicked()","StartFrame",this,"CheckedExpert()");
  fCheckExpert->SetToolTipText("activate additional histograms normal peoples do not care about");
  fTaskFrame->AddFrame(fCheckExpert);


  // Input frame
  fInputFrame = new TGGroupFrame(fMain,"Input files", kVerticalFrame);
  fInputFrame->SetLayoutManager(new TGVerticalLayout(fInputFrame));
  
  fRawFrame = new TGCompositeFrame(fInputFrame, 60, 20, kHorizontalFrame);  
  fRawLabel = new TGLabel(fRawFrame,new TGHotString(fDataFile.c_str()));
  fRawButton = new TGTextButton(fRawFrame, "&Raw data:", 100);
  fRawButton->Connect("Clicked()","StartFrame",this,"HandleButtons()");
  fRawButton->SetToolTipText("Open raw data file");
  fRawFrame->AddFrame(fRawButton,fL3);
  fRawFrame->AddFrame(fRawLabel,fL3);

  fMapFrame = new TGCompositeFrame(fInputFrame, 60, 20, kHorizontalFrame);
  fMapLabel = new TGLabel(fMapFrame,new TGHotString(fMapFile.c_str()));  
  fMapButton = new TGTextButton(fMapFrame, "&Map file :", 101);
  fMapButton->SetToolTipText("Open mapping file");
  fMapButton->Connect("Clicked()","StartFrame",this,"HandleButtons()");
  fMapFrame->AddFrame(fMapButton,fL3);
  fMapFrame->AddFrame(fMapLabel,fL3);

  fGeomFrame = new TGCompositeFrame(fInputFrame, 60, 20, kHorizontalFrame);
  fGeomLabel = new TGLabel(fGeomFrame,new TGHotString(fGeomFile.c_str()));  
  fGeomButton = new TGTextButton(fGeomFrame, "&Geom file :", 102);
  fGeomButton->SetToolTipText("Open mapping file");
  fGeomButton->Connect("Clicked()","StartFrame",this,"HandleButtons()");
  fGeomFrame->AddFrame(fGeomButton,fL3);
  fGeomFrame->AddFrame(fGeomLabel,fL3);

  fInputFrame->AddFrame(fRawFrame);
  fInputFrame->AddFrame(fMapFrame);
  fInputFrame->AddFrame(fGeomFrame);

  // down frame
  fDownFrame = new TGCompositeFrame(fMain, 60, 20, kHorizontalFrame);  

  // decoding Frame
  fDecodingFrame = new TGGroupFrame(fDownFrame,"Decoding", kVerticalFrame);
  
  fSetupButton = new TGTextButton(fDecodingFrame, "&Setup :", 103);
  fSetupButton->SetToolTipText("Open control panel");
  fSetupButton->Connect("Clicked()","StartFrame",this,"HandleButtons()");
  
  fDecodingFrameDown = new TGCompositeFrame(fDecodingFrame,60,20, kHorizontalFrame);
  fEvtMax = new TGTextEntry(fDecodingFrameDown,new TGTextBuffer(10));
  fEvtLabel = new TGLabel(fDecodingFrameDown, "Number of events :");
  fDecodingFrameDown->AddFrame(fEvtLabel,fL3);
  fDecodingFrameDown->AddFrame(fEvtMax,fL2);

  fMinSpacingEntry = new TGTextEntry(fDecodingFrameDown,new TGTextBuffer(10));
  fMinSpacingLabel = new TGLabel(fDecodingFrameDown, "Min spacing of events :");
  fDecodingFrameDown->AddFrame(fMinSpacingLabel,fL3);
  fDecodingFrameDown->AddFrame(fMinSpacingEntry,fL2);

  // Entry to set reference run number
  fDecodingFrameDown2 = new TGCompositeFrame(fDecodingFrame,60,20, kHorizontalFrame);
  fRefEntry = new TGTextEntry(fDecodingFrameDown2,new TGTextBuffer(20));
  fRefLabel = new TGLabel(fDecodingFrameDown2, "Reference run or file :");
  fDecodingFrameDown2->AddFrame(fRefLabel,fL3);
  fDecodingFrameDown2->AddFrame(fRefEntry,fL2);

  fDecodingFrame->AddFrame(fSetupButton,fL3);
  fDecodingFrame->AddFrame(fDecodingFrameDown,fL5);
  fDecodingFrame->AddFrame(fDecodingFrameDown2,fL5);

  // start button 
#if ROOT_VERSION_CODE >= ROOT_VERSION(4,4,2)
  // cheating with const...
  register void* adr_compass_xpm = (void*) &compass_xpm[0];
  char** compass_xpm_noconst =  (char**) adr_compass_xpm;
  fStartButton = new TGPictureButton(fDownFrame, gClient->GetPicturePool()->GetPicture("compass_xpm", compass_xpm_noconst), 102);
#else
  fStartButton = new TGTextButton(fDownFrame,"&Start", 102);
#endif
  fStartButton->SetToolTipText("Start Decoding");
  fStartButton->Connect("Clicked()","MainFrame",fMainframe,"Init()");

  fUpFrame->AddFrame(fCheckFrame,fL5);
  fUpFrame->AddFrame(fTaskFrame);

  fDownFrame->AddFrame(fDecodingFrame,fL1);
  fDownFrame->AddFrame(fStartButton,fL4);

  fMain->AddFrame(fUpFrame,fL5);
  fMain->AddFrame(fInputFrame);
  fMain->AddFrame(fDownFrame,fL5);

  // position relative to the parent's window
  Window_t wdum;
  int ax, ay;
  TGMainFrame *main=fMainframe->GetFrame();
  gVirtualX->TranslateCoordinates(main->GetId(), fMain->GetParent()->GetId(),
				  (((TGFrame *) main)->GetWidth() 
				   + fMain->GetWidth()) >> 1,
				  (((TGFrame *) main)->GetHeight() 
				   + fMain->GetHeight()) >> 1,
                           ax, ay, wdum);
  fMain->Move(ax, ay);


  fMain->SetWindowName("Input");
  fMain->MapSubwindows();
  fMain->Resize(fMain->GetDefaultSize());
  fMain->MapWindow();    


}


StartFrame::~StartFrame() {

  delete fRefEntry;
  delete fRefLabel;
  delete fEvtMax;
  delete fEvtLabel;
  delete fMinSpacingEntry;
  delete fMinSpacingLabel;
  delete fSetupButton;
  delete fDecodingFrameDown;
  delete fDecodingFrame;
  delete fStartButton;
  delete fDownFrame;
  delete fMapButton;
  delete fMapLabel;
  delete fMapFrame;
  delete fRawButton;
  delete fRawLabel;
  delete fRawFrame;
  delete fInputFrame;
  fCheckButton.clear();
  delete fDefButton;
  delete fAllButton;
  delete fNoneButton;
  delete fMonDetFrame;
  delete fAllNoneFrame;
  delete fCheckFrame;
  delete fCheckTree;
  delete fCheckText;
  delete fCheckTrack;
  delete fCheckCluster;
  delete fCheckExpert;
  delete fTaskFrame;
  delete fUpFrame;
#if ROOT_VERSION_CODE >= ROOT_VERSION(4,4,2)
  gClient->FreePicture(gClient->GetPicture("compass_xpm"));
#else
  #if AFS_ACCESS == 1
    gClient->FreePicture(gClient->GetPicture("/afs/cern.ch/compass/detector/monitor/compass.xpm"));
  #endif
#endif
  delete fL5;
  delete fL4;
  delete fL3;
  delete fL2;
  delete fL1;

  //cout<<"deleting fMain"<<endl;
  delete fMain;
  //cout<<"StartFrame deleted"<<endl;
}


void StartFrame::Checked(int id) {

  TGButton *btn = (TGButton *) gTQSender;
  if (id == -1) {
    id = btn->WidgetId();
  }

  std::map<std::string,bool>& detintree=fMonitor->GetDetInTree();
  typedef std::map<std::string,bool>::iterator IM;
  int i=0;
  for(IM mi=detintree.begin();mi!=detintree.end();mi++) {
    if(i==id) {
      if(btn->GetState()==kButtonDown) { 
	mi->second=true;
	//cout<<"true"<<endl;
      }
      else if(btn->GetState()==kButtonUp) { 
	mi->second=false;
	//cout<<"false"<<endl;
      }
    }
    i++;
  }
}

void StartFrame::HandleButtons(int id) {

  if (id == -1) {
    TGButton *btn = (TGButton *) gTQSender;
    id = btn->WidgetId();
  }
  TGFileInfo fi;
  switch(id) {
  case 100 :
    fi.fFileTypes = CHAR_TGFILEINFO(datatype);
    new TGFileDialog(gClient->GetRoot(), fMain, kFDOpen,&fi);	    
    if(fi.fFilename!=0) {
      std::string name=fi.fFilename;
      unsigned int pos=name.find("@",0,1);
      if(pos<name.size()) {
	name.erase(0,pos);
	fDataFile=name.c_str();
      }
      else 
	fDataFile=fi.fFilename;
      TGString *l=new TGString(fDataFile.c_str());
      fRawLabel->SetText(l);
      fRawLabel->Resize(fRawLabel->GetDefaultSize());
      fMain->MapSubwindows();
      fMain->Layout();
      fMain->Resize(fMain->GetDefaultSize());
      fMain->MapWindow(); 
    }
    break;
  case 101 :
    fi.fFileTypes = CHAR_TGFILEINFO(maptype);
    new TGFileDialog(gClient->GetRoot(), fMain, kFDOpen,&fi);	
    if(fi.fFilename!=0) {
      fMapFile=fi.fFilename; 
      TGString *l=new TGString(fMapFile.c_str());
      fMapLabel->SetText(l);
      fMapLabel->Resize(fMapLabel->GetDefaultSize());
      fMain->MapSubwindows();
      fMain->Layout();
      fMain->Resize(fMain->GetDefaultSize());
      fMain->MapWindow();      
    }
    break;    
  case 102 :
    fi.fFileTypes = CHAR_TGFILEINFO(geomtype);
    new TGFileDialog(gClient->GetRoot(), fMain, kFDOpen,&fi);	
    if(fi.fFilename!=0) {
      fGeomFile=fi.fFilename; 
      TGString *l=new TGString(fGeomFile.c_str());
      fGeomLabel->SetText(l);
      fGeomLabel->Resize(fGeomLabel->GetDefaultSize());
      fMain->MapSubwindows();
      fMain->Layout();
      fMain->Resize(fMain->GetDefaultSize());
      fMain->MapWindow();      
    }
    break;    
  case 103:
    fMonitor->ControlPanel(gClient->GetRoot(),fMain);
    break;
  default:
    break;
  }
}


void StartFrame::CheckedTree() {
  TGButton *btn = (TGButton *) gTQSender;
  if(btn->GetState()==kButtonDown) { 
    fMonitor->DoTree(true);
  }
  else if(btn->GetState()==kButtonUp) { 
    fMonitor->DoTree(false);
  } 
}


void StartFrame::CheckedText() {
  TGButton *btn = (TGButton *) gTQSender;
  if(btn->GetState()==kButtonDown) { 
    fMonitor->DoTextOutput(true);
  }
  else if(btn->GetState()==kButtonUp) { 
    fMonitor->DoTextOutput(false);
  } 
}


void StartFrame::CheckedTrack() {
  TGButton *btn = (TGButton *) gTQSender;
  if(btn->GetState()==kButtonDown) { 
    fMonitor->DoTracking(true);
  }
  else if(btn->GetState()==kButtonUp) { 
    fMonitor->DoTracking(false);
  } 
}


void StartFrame::CheckedCluster() {
  TGButton *btn = (TGButton *) gTQSender;
  if(btn->GetState()==kButtonDown) { 
    fMonitor->DoClustering(true);
  }
  else if(btn->GetState()==kButtonUp) { 
    fMonitor->DoClustering(false);
  } 
}


void StartFrame::CheckedRef() {
  TGButton *btn = (TGButton *) gTQSender;
  if(btn->GetState()==kButtonDown) { 
    fMonitor->DoShowRef(true);
  }
  else if(btn->GetState()==kButtonUp) { 
    fMonitor->DoShowRef(false);
  } 
}


void StartFrame::CheckedCalib() {
  TGButton *btn = (TGButton *) gTQSender;
  if(btn->GetState()==kButtonDown) { 
    fMonitor->DoCalib(true);
  }
  else if(btn->GetState()==kButtonUp) { 
    fMonitor->DoCalib(false);
  } 
}


void StartFrame::CheckedReadFailedEvt() {
  TGButton *btn = (TGButton *) gTQSender;
  if(btn->GetState()==kButtonDown) { 
    fMonitor->DoReadFailedEvt(true);
  }
  else if(btn->GetState()==kButtonUp) { 
    fMonitor->DoReadFailedEvt(false);
  } 
}


void StartFrame::CheckedExpert() {
  TGButton *btn = (TGButton *) gTQSender;
  if(btn->GetState()==kButtonDown) {
    fMonitor->DoExpertHisto(true);
  }
  else if(btn->GetState()==kButtonUp) {
    fMonitor->DoExpertHisto(false);
  } 
}


void StartFrame::DefClicked() {
  for(unsigned i=0; i<fCheckButton.size(); i++ ) {
    fCheckButton[i]->SetState(kButtonDown);
  }
  std::map<std::string,bool>& detintree=fMonitor->GetDetInTree();
  typedef std::map<std::string,bool>::iterator IM;
  for(IM mi=detintree.begin();mi!=detintree.end();mi++) {
    mi->second=true;
  }
  IM im;
  register int i=0;
  for(IM mi=detintree.begin(); mi!=detintree.end(); mi++) {
// exclude always   
     if (   (mi->first == "RICH_BORA")
         || (mi->first == "VetoBox")
// for non-Primakoff hadron runs
//         || (mi->first == "BMS")
// for hadron runs
//         || (mi->first == "SciFis(Warsaw)")
// for muon runs
//         || (mi->first == "RPD")
//         || (mi->first == "Sandwich")
//         || (mi->first == "CEDARs")
// for Drell-Yan runs:
         || (mi->first == "Camera")
         || (mi->first == "ECAL0")
         || (mi->first == "RPD")
         || (mi->first == "Sandwich")
         || (mi->first == "SciFis(Warsaw)")
         || (mi->first == "FEM")
         || (mi->first == "Micromegas")
         || (mi->first == "RICH_APV")
         || (mi->first == "Silicons")
        ) {
      mi->second=false;
      fCheckButton[i]->SetState(kButtonUp);
    }
    i++;
  }
}


void StartFrame::AllClicked() {
  for(unsigned i=0; i<fCheckButton.size(); i++ ) {
    fCheckButton[i]->SetState(kButtonDown);
  }
  std::map<std::string,bool>& detintree=fMonitor->GetDetInTree();
  typedef std::map<std::string,bool>::iterator IM;
  for(IM mi=detintree.begin();mi!=detintree.end();mi++) {
    mi->second=true;
  }
}


void StartFrame::NoneClicked() {
  for(unsigned i=0; i<fCheckButton.size(); i++ ) {
    fCheckButton[i]->SetState(kButtonUp);
  }
  std::map<std::string,bool>& detintree=fMonitor->GetDetInTree();
  typedef std::map<std::string,bool>::iterator IM;
  for(IM mi=detintree.begin();mi!=detintree.end();mi++) {
    mi->second=false;
  }
}

