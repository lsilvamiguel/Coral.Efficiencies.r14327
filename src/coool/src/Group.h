#ifndef __Group__
#define __Group__

#include <string>
#include <vector>
#include "Plane.h"
#include "Reference.h"

#ifndef __CINT__
#include "DaqEvent.h"
#endif

class GroupPanel;

class Group {

 private:

  /// default reference directory is a fairly good run
  static std::string defaultRefFileName;

  /// forced reference file is given by user
  static std::string forcedRefFileName;
  
  /// working environment string which tags the environment where Coool should run
  static std::string wkEnvStr;

 protected:

  /// where reference histos can be found
  TDirectory *fReferenceDirectory;

  /// Group name
  std::string fName;

  /// Group type 
  std::string fType;

  /// vector of Plane's belonging to the group
  std::vector<const Plane*> fPlanes;

  /// vector of histograms declared in this group
  std::vector<TH1*> fHistList;

  /// Associated GroupPanel
  GroupPanel *fControlPanel;

  /// Number of events seen by Monitor
  int fRateCounter;

  /// if true, activate expert histograms
  static bool fExpertHistos;

  /// do we do tracking ? needed for initializations
  static bool fDoTracking;

#if USE_DATABASE == 1
  /// pointer to database
  MySQLDB* fDataBase;

  /// returns name of the reference file from the database
  const char* GetRefDirNameDB();
#endif

  /// give a directory pointer to the default reference file
  TDirectory *GetDefaultRefDir() {return fReferenceDirectory;}

  /// opens root reference file
  virtual void OpenReference();

  /// returns name of the reference file
  const char* GetRefDirName();

  /// add one histogram to the histogram list
  void AddHistogram(TH1* hist);


 public:
  Group(const char* name):
    fReferenceDirectory(0), fName(name), fType("NOTYPE"), fControlPanel(0),
      fRateCounter(0) {}
  virtual ~Group();

  /// must be called after group is complete
  virtual void Init() {  fRateCounter = 0; }

  /// sets group type (yes, this should be in the constructors. yes I'm lazy)
  void SetType(const char* type) {fType = type;}

  /// resets all histograms
  virtual void ResetHistograms();

  /// resets all histograms if one plane was modified
  void ResetIfNeeded();

  /// called at new run
  virtual void NewRun(int runnumber) {}

  /// Add a plane to the list of contained Planes
  void Add(const Plane* plane);

  /// Fills the histograms
  virtual void EndEvent() {};    

#ifndef __CINT__
  /// Fills the histograms
  virtual void EndEvent(const CS::DaqEvent &event) {};
#endif

  /// Do the tracking
  virtual void Tracking() {}

  /// \return vector of Plane's
  std::vector<const Plane*>& GetPlane() {return fPlanes;}

  /// \return vector oh histograms
  std::vector<TH1*>& GetHistoList() {return fHistList;}

  /// \return name of the group
  const char* GetName() const {return fName.c_str();}

  /// \return type of the group
  const char* GetType() const;

  /// Creates the GroupPanel
  virtual void ControlPanel(const TGWindow *p, const TGWindow *main);

  /// called by the associated GroupPanel when it is closed
  void PanelClosed() {fControlPanel=0;}

  /// update fRateCounter (increase it by 1 basically...)
  void UpdateRateCounter() { fRateCounter++; }

  /// find an histo from the histo list
  TH1* GetHisto( const std::string hName );

#if USE_DATABASE == 1
  void setDBpt(MySQLDB* dbpt)  { fDataBase = dbpt; }
#endif

  /// gives a file name which will be used as reference instead of DB ones
  static void ForceRefFileName(const char* filename)
    { forcedRefFileName = filename; }

  /// set the working environment string
  static void setWkEnvStr(const char* wkEnvStrChr)
    { wkEnvStr = wkEnvStrChr; }

  /// set expert histo flag
  static void setExpertHistoFlag(bool experthisto)
    { fExpertHistos = experthisto; }

  /// set dotracking flag
  static void setDoTrackingFlag(bool state)
    { fDoTracking = state; }

  ClassDef(Group,0)

};

#endif




