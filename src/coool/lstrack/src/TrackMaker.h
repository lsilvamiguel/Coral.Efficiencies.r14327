#ifndef __TrackMaker
#define __TrackMaker

// uncomment for using the FastSkip class (see below)
#define __FASTSKIP__

#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <unistd.h>

#include "TObject.h"
#include "TFile.h"
#include "TMatrix.h"
#include "TVector.h"
#include "TGeometry.h"
#include "TNode.h"
#include "TBRIK.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TStopwatch.h"

//#include "Geometry.h"

#include "Tracker.h"
#include "Doublet.h"
#include "Station.h"

// class FastSkip accelerates finding tracks with
// good chi2 by skipping the calculation of chi2 for tracks
// where a good chi2 can be excluded from a previous (bad)
// chi2 calculation. 
// developed and coded 11/2002 Jan Friedrich, E18, TU Munich

#ifdef __FASTSKIP__
class FastSkip {
  int ntrackers;
  int lat;
  int combMult;
  int *multiplik;
  int *max_index;
  int *fastCount;
  bool* momCom;
  int *mComb;
  bool skipPast;
  short *jmin, *jmax;

  int shMult;
  int *shMultiplik;
  unsigned short* skipShort;
  
 public:
  FastSkip(int ntr);
  void Reset(int *cluMult);
  bool CheckCluster(int *cl);
  
  bool ShortLT(int ti, int i, int *cl);
  
  void LoopTracker(int ti, int i, int *cl);
  void AddSkip(int *cl);
  void PrintStatistics();
#ifdef __FASTSKIP_TIMER__
  TStopwatch *fBenchLoopTracker, 
    *fBenchCheckCluster;
#endif
};
#endif    

class TrackMaker : public TObject {
  
 private:
  enum Modes {SP=1, LSF=2};
  
  static const unsigned int fCnofDoublets;
  static const unsigned int fMaxTracks = 100;
  
  //name
  std::string fName;
  
  //const Geometry* fGeometry; 
  TFile *fFile; 
  
  // cuts
  int fMinSize;
  float fMinpchi2;
  float fMaxchi2;
  float fR2min;
  float fR2max;
  float fYpmin;
  float fYpmax;
  float fZpmin;
  float fZpmax;

  // beginning and end of zone
  float fXmin;
  float fXmax;

  // abscissa of reconstruction plane
  float fXRecon;

  // arrays for track block in output tree
  float *fX;
  float *fY;
  float *fZ;
  float *fYp;
  float *fZp;
  float *fChi2;

  float fChi2Sum;
  int fChi2N;

  int fEventMask;   // selects entries for tree filling

  bool fBadEvent; // Set at clustering stage if the number of clusters is not 
                  // correct
  
  Int_t fMode;   // 1: space points,  2: least squares fit
  unsigned fNcomb;  // number of combinations
  int fNcombOk; // number of correct combinations

  //  map<int, TMatrix*> fCovMats;  Not yet implemented
  TMatrix *fCovMat;       // covariance matrix
  TVector *fConstVec;     // vector of constants

  std::vector<Tracker*> fTrackers;       // all trackers
  std::vector<Tracker*> fTrackersSel;    // selected trackers
  std::vector<Tracker*> fTrackersID;     // all trackers, indexed by ID
                                    // -> full of holes ! 

  std::vector<Tracker*> fCheckTrackers;  // used only for space points
  std::vector<Doublet*> fDoublets;       // used only for space points
  std::vector<Station*> fStations;       // used only for space points
  
  std::vector<Cluster*> fCurComb;     // if no cluster, fCurComb[i]=0
  int              fCurCombSize;

#ifdef __FASTSKIP__
  int*             fCurClusterInds;
  int*             fCurClusterMult;
#endif

  std::vector<Track*> fTracks;   // owns the tracks
  std::vector<Track*> fTracksOk; 

  // event display 
  TGeometry *gGeometry;
  TBRIK *fWorldBrik;
  TNode *fWorldNode;


  void AddTracker(Tracker* tracker);

 public:

  TH1F *fHyp;
  TH2F *fHypvsy;
  TH1F *fHzp;
  TH2F *fHzpvsz;
  TH1F *fHchi2;
  TH1F *fHlogchi2;
  TH1F *fHNtracks;
  TH1F *fHprobchi2;
  TH2F *fHzvsy;

  std::vector<TH1*> fHists;


  TrackMaker(const char* name, Int_t mode,float maxchi2);
  TrackMaker(const char* name, Int_t mode,float maxchi2, std::ifstream &in, 
	     TFile *file);
  // TrackMaker(Int_t mode,float maxchi2, const Geometry *geom, TFile *file);
  ~TrackMaker();

  void Default(); // default system

  Tracker* AddTracker(Int_t id, const char *name, Int_t nwires, Float_t x, 
		      Float_t y, Float_t z, float dx, float dy, float dz, 
		      Float_t angle, Float_t inpitch, float outpitch);

  Tracker* AddTracker(Int_t id, const char *name, Int_t nwires, Float_t x, 
		      Float_t y, Float_t z, float dx, float dy, float dz,
		      Float_t angle, Float_t inpitch, float outpitch, 
		      float ctmin, float ctmax, float ctotmin);

  void AddTracker(const Tracker* tracker); 

  // Necessary only for SP
  Doublet* CreateDoublet(const char* name, unsigned int id1, unsigned int id2);
  Station* CreateStation(const char* name, Doublet* d1, Doublet *d2);

  void Init();      // called when the system is declared
  void Reset();     // called before each event
  void ResetHistos(); // resets all histograms

  // add digit to plane id
  void AddDigit(unsigned int id, int channel, int t, int tot);

  // add cluster to plane id (cm, origin at the center)
  void AddClusterWRS(unsigned int id, float pos, int size, float res);

  // add cluster to plane id (cm, origin at the center)
  void AddClusterWRS(unsigned int id, float pos, float time, int size, float res);

  // add cluster to plane id (unit is channel number)
  void AddCluster(unsigned int id, float pos, float t, float tot, float s);

  // clusterize all planes
  Int_t Clusterize(int nclustperplane);      // clustering for all planes

  // finds all cluster combinations
  Int_t Candidates();

  // recursive cluster association                      
  int   NextTracker(); 

  // 2b removed
  void  LsfTrack();
  
  // create a track from current cluster combination (LSF)
  void  MakeTrack();

  // calculate efficiency for all tracks
  void  Efficiency();

  // checks candidate tracks using check planes (Space Points)
  Int_t CheckCandidates(Float_t checkwidth, unsigned int okthr);

  // calculates residuals on all planes
  void  Residuals();
  
  // draw full system
  void Draw(bool plane,bool station, bool doublet, bool tracks);

  // prints information
  void Dump();
  
  // prints geometry
  void OutConfig(const char *filename);
  
  

  // plots chi2 probability
  void PlotProbchi2();

  // plots yp
  void  PlotYp();

  // plots zp
  void  PlotZp();  

  // plots all residuals vs pos
  void PlotResVsPos();

  // plots all residuals
  void PlotRes();  

  // plots all cluster times
  void PlotCTimes();

 
  // plots profile @ x=0
  void PlotProfile();

  // adds tracks block to the tree
  void OutputTree(TTree *tree);

  // returns TrackMaker name
  const char* GetName() const {return fName.c_str();}

  // get number of good combinations
  Int_t GetNCombOk() {return fNcombOk;}

  // returns plane with index i
  Tracker* GetTracker(Int_t i) {return fTrackersID[i];}

  // returns vector of all planes
  std::vector<Tracker*>& GetTrackers() {return fTrackers;}

  // returns vector of selected planes. (? wildcard)
  std::vector<Tracker*>& GetTrackers(const char* selection);
  
  // returns doublet with index i
  Doublet* GetDoublet(int i) {return fDoublets[i];}

  // returns station with index i
  Station* GetStation(int i) {return fStations[i];}

  // returns event mask
  int GetMask() const {return fEventMask;}

  // returns vector of good tracks (SP) -> 2b modified
  const std::vector<Track*>& GetTracks() {return fTracksOk;}

  // returns vector of good tracks (SP) -> 2b modified
  const std::vector<Track*>& GetAllTracks() {return fTracks;}
  
  // returns vector of histograms
  std::vector<TH1*>& GetHistograms() {return fHists;}

  // returns mean chi2
  float GetMeanChi2() {return fChi2Sum/fChi2N;}
  
  // sets min number of planes for tracking (LSF)
  void SetNTrackers(int nplanes) {fMinSize=nplanes;}

  // sets chi2 cut
  void SetPChi2(float pchi2) {fMinpchi2=pchi2;}

  // sets chi2 cut
  void SetChi2(float chi2) {fMaxchi2=chi2;}

  void TimeCuts(float min, float max);

  void TimeCuts(int nsigma);

  // removes all cuts
  void NoCuts() {
    fMinpchi2=0;
    fMaxchi2=0;
  }

  // sets R range @ x=0. Should be done on the considered plane 
  void SetRRange(float rmin, float rmax) {
    fR2min=rmin*rmin; 
    fR2max=rmax*rmax;
  }

  // sets range in tany
  void SetYpRange(float ypmin, float ypmax) {
    fYpmin=ypmin;
    fYpmax=ypmax;
  }

  // sets range in tanz
  void SetZpRange(float zpmin, float zpmax) {
    fZpmin=zpmin;
    fZpmax=zpmax;
  }

#ifdef __FASTSKIP__
  FastSkip* fastSkip;
#endif
  
  bool string_match(const char* str, const char* pattern);

  ClassDef(TrackMaker,0)
};

#endif





