#ifndef __Plane__
#define __Plane__

#include <string>
#include <vector>
#include <ctime>

#include "TFile.h"
#include "TWebFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"
#include "TGWindow.h"

#include "threadstate.h"

#ifndef __CINT__
#include "config.h"
#include "Chip.h"
#include "DaqEvent.h"
#include "DaqEventsManager.h"
#endif

#include "Variable.h"
#include "Containers.h"
#include "Reference.h"


#if USE_DATABASE == 1
#include "MySQLDB.h"
#endif

#include "GeomPlane.h"

#define ROOT_30200 197120
#if ROOT_VERSION_CODE >= 197120  /* 3.02.00 */
    #define CHAR_TGFILEINFO(a) a
#else
    #define CHAR_TGFILEINFO(a) const_cast<char**>(a)
#endif

using namespace std;

class PlanePanel;

/// abstract Plane base class
class Plane : public TObject {

 protected:
  TDirectory *fReferenceDirectory;

  /// Number of events seen by Monitor
  int fRateCounter;


 private:

  /// map<TBName, pointer to Plane instance>
  static std::map<std::string,const Plane*> sInstances;

  /// default reference file is a fairly good run
  static std::string defaultRefFileName;

  /// forced reference file is given by user
  static std::string forcedRefFileName;

  /// working environment string which tags the environment where Coool should run
  static std::string wkEnvStr;

  /// vector of histos
  std::vector<TH1*> fHistList;

  /// vector of variables (subset of fHistList)
  std::vector<Variable*> fVariables;

  /// set to true if one plane has been modified
  static bool sModified;

  /// true if this plane is specified to be modified
  bool fModified;

  /// do we need geometry ?
  bool fNeedGeom;


 protected:

  /// do we do tracking ? needed for initializations
  static bool fDoTracking;

  /// if true, activate expert histograms
  static bool fExpertHistos;

  /// geometry
  GeomPlane *fGeom;

  /// name of the detector
  std::string fName;

  /// True if the detector is in the global tree. 2bremoved
  bool fIsInTree;

  /// number of hits
  int fNhits;

  /// number of hits passing the cuts
  int fNhitsKept;

#ifndef __CINT__
  /// list of Digit kept for this plane
  std::list<CS::Chip::Digit*> lDigits;

//   const CS::DaqEvent *fEvent;  // not used
#endif

  /// control panel
  PlanePanel *fControlPanel;

  /// trigger mask
  unsigned int fTMPattern;

  /// calibration availability flag
  bool  fUseCalib;

  /// tree
  TTree *fTree;

#if USE_DATABASE == 1
  /// pointer to database
  MySQLDB* fDataBase;

  /// reads calibration data from the database
  template <class T>
  void ReadFromDataBase(T& v, const tm& t) {
    if( fDataBase==NULL )
      throw CS::Exception("Plane::ReadFromDataBase():  data base is not opened.");
    try{
      struct tm tt(t);
      CDB::Time tp(mktime(&tt),0);
      std::string strdata("");
      fDataBase->read(fName,strdata,tp);
      if (strdata == "") {
        std::cerr << "Plane::ReadFromDataBase() "<<GetName()<<" : no calibration file found"<<std::endl;
        return;
      }
      std::istringstream istrdata(strdata.c_str());
      istrdata >> v;
    }
    catch(CS::Exception& e) {
      std::cout<<"rethrowing"<<std::endl;
      throw;
    }
  }


  /// returns name of the reference file from the database
  const char* GetRefDirNameDB();
#endif

  /// returns name of the reference file
  const char* GetRefDirName();

  /// create one variable and add it to the variable list
  Variable* AddVariable(const char* name, int nbins, float min, float max,
			int maxsize);

  /// add one variable to the variable list
  void AddVariable(Variable* var);

  /// add one histogram to the histogram list
  void AddHistogram(TH1* hist);

 public:

  Plane(const char* name, TTree *tree=0);

  virtual ~Plane();

  /// Resets the plane. Must be called at the beginning of each event
  virtual void Reset();

  /// Resets all histograms
  virtual void ResetHistograms();

  /// tells Plane if this Plane has been modified
  void Modified(bool state) {
    fModified = state;
    if(state) sModified = true;
  }

  #ifndef __CINT__
  const std::list<CS::Chip::Digit*>& GetDigits(void) const {return lDigits;}
  #endif

  bool IsModified() const {return fModified;}

  /// semaphore : tells if one Plane has been modified
  static bool ModifiedSemaphore() {return Plane::sModified;}

  /// aknowledge ModifiedSem()
  static void ModifiedAcknowledge() {Plane::sModified = false;}

  /// Called at new run
  virtual void NewRun(int runnumber) {}

#ifndef __CINT__
  virtual void StoreDigit(CS::Chip::Digit* digit)=0;

  /*! \brief Called at the end of the event
    Applies the cuts, fills histograms and tree.
    This function must be called at the end of the event
  */
  virtual void EndEvent(const CS::DaqEvent &event) {}
#endif

  /// do the clustering
  virtual void Clusterize() {}

  /// Books histograms and branchs tree
  virtual void Init(TTree* tree)=0;

  /// Books histograms and branchs tree
  virtual void Init() {
    Init(fTree);
  }

  /// Sets tree pointer
  void SetTree(TTree *tree) {fTree=tree;}

  /// Creates the PlanePanel
  virtual void ControlPanel(const TGWindow *p, const TGWindow *main) {}

  /// Sets the trigger mask
  void SetTriggerMask(unsigned int tmask) {fTMPattern = tmask;}

  /// \return detector name
  const char* GetName() const {return fName.c_str();}

  /// \return detector type
  const char* GetType() const;

  /// give a directory pointer to the default reference file
  TDirectory *GetDefaultRefDir() {return fReferenceDirectory;}

  /// opens root reference file
  virtual void OpenReference();

  /// \return vector of histograms
  std::vector<TH1*>& GetHistoList() {return fHistList;}

  /// \return variable named vName if any or 0
  TH1* GetHisto( const std::string hName );

  /// \return vector of variables
  std::vector<Variable*>& GetVariables() {return fVariables;}
  const std::vector<Variable*>& GetVariables() const {return fVariables;}

  /// \return variable named vName if any or 0
  Variable* GetVariable( const std::string vName );

  /// true if the detector is in the tree
  bool          IsInTree() {return fIsInTree;}

  /// \returns number of hits passing the cuts
  int           GetNhits() {return fNhitsKept;}

  /// \returns associated geomplane
  const GeomPlane* GetGeometry() const {return fGeom;}

  /// called by PlanePanel::CloseWindow()
  void PanelClosed() {fControlPanel=0;}

  /// read the calibration data
  virtual void ReadCalib(const tm& t) {};

  /// copy constructor
  bool operator== (const Plane&) const;

  /// write Plane output to a stream
  virtual void TextOutput(ostream& out) {}

  /// called at the end, before writing histograms and text output
  virtual void End() {}

#if USE_DATABASE == 1
  void setDBpt(MySQLDB* dbpt)  { fDataBase = dbpt; }
#endif

  void SetGeometry(GeomPlane* geom);

  /// formatted writing of a set of data.
  friend ostream& operator<<(ostream& out, const std::set<int>& dataset);

  /// tells if geometry is needed by this plane
  bool NeedGeom() const {return fNeedGeom;}

  /// geometry is needed by this plane
  void INeedGeom() {fNeedGeom = true;}

  /// update fRateCounter (increase it by 1 basically...)
  void UpdateRateCounter() { fRateCounter++; }

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

  /// returns const pointer to another plane instance
  const Plane* GetPlane(const char* name) const;

  ClassDef(Plane,0)

};

////////////////////////////////////////////////////////////////////////////////

template <class T>
istream& operator>>(istream &in,std::vector<T> &v)
{
  v.clear();

  do
    v.push_back(T());
  while( in>>v.back() );

  v.pop_back();  // Remove last element becase it was not read.

  if( in.eof() )
    in.clear();

  return in;
}

// added for gcc 3.2.2
// template <class _T1, class _T2>
// inline std::pair<_T1, _T2> make_pair(const _T1& __x, const _T2& __y)
// {
//   return std::pair<_T1, _T2>(__x, __y);
// }

#endif









