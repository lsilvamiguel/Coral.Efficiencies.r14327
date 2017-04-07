#ifndef __Monitor__
#define __Monitor__

// stl
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include <set>
#include <string>
#include <map>

// compass decoding library

#ifndef __CINT__
#include "DaqEventsManager.h"
#include "config.h"
#include "utils.h"
// #include "Chip.h"
#include "DaqEvent.h"
// #include "ChipF1.h"
// #include "ChipADC.h"
// #include "ChipAPV.h"
#include "DaqOption.h"
#endif

// monitoring headers

#include "threadstate.h"

#include "Group.h"

#include "Plane.h"

// root

#include "TTree.h"
#include "TROOT.h"
#include "TFile.h"
#include "TThread.h"
#include "TH1.h"
#include "TStopwatch.h"

// calib database

#if USE_DATABASE == 1
#include "MySQLDB.h"
#endif

class MonitorPanel;

/*! \brief Central monitoring class.

  This class creates a set of Plane objects from what it reads
  in the mapping file.

  Decoding then starts as a separate thread, where digits are passed
  to the Plane's.
*/
class Monitor : public TObject {

 private:

  /// decoding thread
  TThread *fTh;

  /// true if the thread is running
  bool fThreadRun;

  /// true if coool is closing (to avoid mutex lock conflicts)
  bool fClosingWindow;

  /// Flag to reset histos at each begin of run
  bool  fRunClear;

  /// current run number
  unsigned int fRunNumber;

  /// time of the first event of the run
  tm  fRunTime;

#ifndef __CINT__

  /// mapping file
  //AZ//string fMapfile;

  /// decoding manager
  CS::DaqEventsManager *fManager;

  /// decodes a given chip. digits are passed to Planes here.
//   static void         Decode(CS::Chip& chip);
#endif

  /// map of detector planes
  static std::map<std::string,Plane*> fDetMap;
  static std::map<std::string,Plane*> fDetMapRegistered;
  static std::map<std::string,Plane*> fDetMapNotRegistered;
  static std::vector<Plane*>       fPlanes;

  /// map of Plane Groups
  static std::map<std::string, Group*> fGroupMap;

  static const char* fAllPlanesName;

  /// filename for personnal options
  static const char* fPersoFile;

  /// detector types included in the tree
  std::map<std::string,bool> fDetInTree;

  /// global tree
  TTree              *fTree;

  /// if true, a tree is created and filled
  bool                fDoTree;

  /// if true, tracking is performed
  bool                fDoTracking;

  /// if true, clustering is performed
  bool                fDoClustering;

  /// if true, text output
  bool                fDoTextOutput;

  /// if true, calibrations will be used
  bool                fDoCalib;

  /// if true, decoding will stop if it meets an end of run event
  bool                fStopAtEndRun;

  /// if true, read event even if decoding failed (for calib events)
  bool                fReadFailedEvt;

  /// if true, activate expert histograms
  bool                fExpertHistos;

  /// if true, event reading is suspended, restart when come back to false
  bool                fStopGo;

  /// output root file
  TFile              *fRootFile;

  /// flag set if rootfile has been flushed to disk (but not really closed)
  bool       fRootFileFlushed;

  /// output root dir
  std::string              fRootDir;

  /// output text file name
  std::string              fTextFile;

  /// output text dir
  std::string              fTextDir;

  /// number of root files opened so far
  int                 fNfiles;

  /// stream from DATE file
  ifstream           *fRawData;

  /// stream from mapping file
  //ifstream           *fMapXml;

  /// event number
  int                 fNevent;

  /// number of events to be decoded
  int                 fNeventMax;

  /// evt-in-spill number in spill of last decoded Event
  int                 fLastDecodedEventnum;

  /// minimum difference between two accepted event numbers (to avoid sampling only the start of spill)
  int                 fMinEvtspacing;


  /// if false, ends decoding
  bool       flag_end;

  /// if true, the decoding thread is finished
  bool       fThreadFinished;

  /// true if Init has been called
  bool                fIsInitialized;

  /// spill number
  int                 fSpill;

  /// event number in run
  int                 fEvtinrun;

  /// event number in spill
  int                 fEvtinspl;

  /// time in second
  int                 fTimesec;

  /// time in microsecond
  int                 fTimemic;

  /// current event type
  int                 fEvtType;

  /// all event types
  std::map<int,std::string>     fEvtTypes;

  std::set<int>            fEvtTypesOK;

  /// trigger pattern
  unsigned int        fPattern;
  /// trigger attributes,
  int                 fAttributes;

  /// Online Filter bits
  int                 fFilterBits;

  /// Associated control panel
  MonitorPanel        *fControlPanel;

  /// trigger mask
  unsigned int        fTMPattern;

  /// only consider lower 16 bits of the trigger mask
  bool                fTMLowerOnly;

  /// strictly compare the trigger mask
  bool                fTMStrict;

#if USE_DATABASE == 1
  /// pointer to database (if USE_DATABASE defined)
  MySQLDB             *fDataBase;
#endif

  /// benchmark for decoding
  TStopwatch *fBenchDecoding;

  /// benchmark for decoding
  TStopwatch *fBenchClustering;

  /// benchmark for decoding
  TStopwatch *fBenchTracking;

  /// fills detectors list from information read in mapping file
  void                FillDetList();

  /// read mapping
  void                ReadMaps(int run,const char *mapfile);

  /// creates a NewDetector or returns a pointer to it if it exists
  Plane*              NewDetector(const char* detname, unsigned int detnum);
  
  /// working environment string which tags the environment where Coool should run
  std::string fWkEnvStr;



 public:

  void Run();


  Monitor();
/*    Monitor(string mapfile, string rootfile,string datafile, string groupfile,string geomfile); */
  virtual ~Monitor();

  /// Opens the I/O files and fills the detectors list
  void  Init(std::string mapfile, std::string rootfile, std::list<std::string>& datafilelist,
             std::string groupfile, std::string geomfile);

  /// Opens the I/O files and fills the detectors list for only one datafile
  void  Init(std::string mapfile, std::string rootfile, std::string datafile,
             std::string groupfile, std::string geomfile);

  void ResetTree();

  /// \return Plane map
  std::map<std::string,Plane*>& GetDetMap() const {return fDetMap;}

  /// \return Group map
  std::map<std::string,Group*>& GetGroupMap() const {return fGroupMap;}

  /// \return Plane with name detname
  Plane* GetPlane(const char* detname) const {
    std::string sdetname=detname;
    std::map<std::string, Plane*>::iterator gi=fDetMap.find(sdetname);
    if(gi!=fDetMap.end()) return gi->second;
    else return 0;
  }

  /// \return Group with name grpname
  Group* GetGroup(const char* grpname) const {
    std::string sgrpname=grpname;
    std::map<std::string, Group*>::iterator gi=fGroupMap.find(sgrpname);
    if(gi!=fGroupMap.end()) return gi->second;
    else return 0;
  }

  /// \return map of "IsInTree flags" for all Planes
  std::map<std::string,bool>& GetDetInTree() {return fDetInTree;}

  /// \return event number
  int GetEventNumber() const {return fNevent;}
  unsigned int GetRunNumber() const {return fRunNumber;}
  bool GetThreadState() const {return fThreadRun;}

  /// sets the number of events to be decoded
  void SetEventNumber(int nevent) {
    if(fIsInitialized && nevent)  fManager->SetEventsMax(nevent);
    fNeventMax = nevent;
  }

  /// sets the minimum event spacing
  void SetEventSpacing(int nevent) {
    fMinEvtspacing = nevent;
  }

  /// sets root output directory
  void SetRootDir(const char* rootdir) {fRootDir = rootdir;}

  /// sets text output file
  void SetTextFile(const char* textfile) {fTextFile = textfile;}

  /// sets text output directory
  void SetTextDir(const char* textdir) {fTextDir = textdir;}

  /// reset all the plane and group histos
  void ResetHistos();

  /// called at new run
  void NewRun(int runnumber);

  /// Start decoding thread
  int ThreadStart();

  /// Stop decoding thread
  int ThreadStop();

  /// Request to stop decoding thread
  void ThreadStopRequest();

  /// close the root file
  void                CloseRootFile();

  /// find an histo from all the planes
  TH1* GetHisto(const std::string hName);

  /// decoding Thread
  static void Thread0(void* arg);

  /// Add a new group
  Group*  NewGroup(const char* grpname, const char* grptype = 0);

  /// Load the groups description from a xml file
  Bool_t LoadGroups(const char* filename);

  /// Sets the flag for stopping at EOR event
  void StopAtEndRun(bool state) {fStopAtEndRun=state;}

  /// Sets the flag for activating expert histograms
  void DoExpertHisto(bool state) {fExpertHistos=state;}

  /// Sets the flag for stopping at EOR event
  void DoReadFailedEvt(bool state) {fReadFailedEvt=state;}

  /// Sets the flag for tree filling
  void DoTree(bool state) {fDoTree=state;}

  /// Sets the flag for tracking
  void DoTracking(bool state) {fDoTracking=state;}

  /// Sets the flag for clustering
  void DoClustering(bool state) {fDoClustering=state;}

  /// Sets the flag in TH1F_Ref for showing reference plots
  void DoShowRef(bool state) { TH1F_Ref::SetShowRef(state); }

  /// Sets the global flag for using or not calibrations
  void DoCalib(bool state) {fDoCalib = state;}

  /// Sets the flag for text ouput
  void DoTextOutput(bool state) {fDoTextOutput=state;}

  /// Creates the control panel
  void ControlPanel(const TGWindow *p, const TGWindow *main);

  /// Called when the control panel is closed
  void PanelClosed() {fControlPanel=0;}

  /// Sets the trigger mask
  void SetTriggerMask(unsigned int tmask) {
    fTMPattern = tmask;
  }

  /// Set the trigger mask consider only lowest 16 bits flag
  void SetTriggerMaskLowerOnly(bool lowerOnly) {
    fTMLowerOnly = lowerOnly;
  }

  /// Set the trigger mask strict comparison flag
  void SetTriggerMaskStrict(bool strict) {
    fTMStrict = strict;
  }

  /// sets accepted event types
  void CheckEvtTypes(const std::set<int>& evttypesok) {
    fEvtTypesOK = evttypesok;
  }

  /// sets accepted event types from a string "type1,type2,..."
  void CheckEvtTypes(const char* evttypesstr);


#ifndef __CINT__
  /// gets list of event types;
  const std::map<int, std::string>& GetEvtTypes() const {
    return fEvtTypes;
  }
#endif

  /// Set and Get the run clear flag
  void SetRunClear(bool f)  { fRunClear = f; }
  bool GetRunClear()  { return fRunClear; }

  /// Set and Get the stopgo flag (suspend and restart events reading)
  void SetStopGo(bool f)  { fStopGo = f; }
  bool GetStopGo()  { return fStopGo; }

  bool GetThreadFinished()  { return fThreadFinished; }

  void SetClosingWindow(bool f)  { fClosingWindow = f; }

  /// opens and write in text output file
  void TextOutput(int runnumber);

  /// set the working environment string
  void SetWkEnvStr(const char* wkEnvStrChr)
    { fWkEnvStr = wkEnvStrChr; }

  /// return rootfile TFile* pointer
  TFile* GetRootFile()  { return fRootFile; }


  ClassDef(Monitor,0)

};

#endif






