// $Id: TEv.h 14069 2015-09-17 20:44:46Z lsilva $

#ifndef TEv_h
#define TEv_h

/*!
  \class TEv
  \brief Event-for-tracking

  This class provides
  an interface to pattern recognition
  and track fit algorithms

  \warning Only one instance the class is allowed
  \author Sergei.Gerassimov@cern.ch
*/

/*
  (Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/TEv.h":
  i) Data members "eventTime,eventTRef" and related methods "SetEventTime",
  "SetEventTime" and "UpdateDrifts".
  ii) Methods specific to lattice: "Quadruples","CleanTrackList","DumpEvent",
  "ForeTrack2*","BackTrack*", "getAssociatedHit", etc... And class "TSpacePt"
  used by some of them.
  iii) Data members "reTracking" and "BMSSmearing" (cf. "TEv::TracksFit2") and
  related accessors and mutators.
*/

#include <cassert>
#include <CsSTD.h>
#include "CsVertex.h"
#include "THit.h"
#include "THitMC.h"
#include "TKine.h"
#include "TVtxMC.h"
#include "TTrack.h"

class CsEvent;
class CsCluster;
class CsZone;
class CsTrack;

class TSpacePt {
public:
  TSpacePt() {
    y=z = 0; nHits = 0; hPat = 0;
    for (int i = 0; i<8; i++) hs[i] = 0;
    chi2 = 0;
    sc2 = 0; ssc = 0; ss2 = 0;
    scu = 0; ssu = 0; su2 = 0;
  }
  TSpacePt(const TSpacePt &spt) {
    y = spt.y; z = spt.z; nHits = spt.nHits; hPat = spt.hPat;
    for (int i = 0; i<8; i++) hs[i] = spt.hs[i];
    chi2 = spt.chi2;
    sc2 = spt.sc2; ssc = spt.ssc; ss2 = spt.ss2;
    scu = spt.scu; ssu = spt.ssu; su2 = spt.su2;
  }
  void reset() {
    nHits = 0; hPat = 0; chi2 = 0;
    sc2 = 0; ssc = 0; ss2 = 0;
    scu = 0; ssu = 0; su2 = 0;
  }
  float y, z;
  int nHits; unsigned int hPat;
  THit *hs[8]; // *hZ, *hZp, *hY, *hYp, *hU, *hUp, *hV, *hVp
  double chi2, sc2, ssc, ss2, scu, ssu, su2;
};

class TTrackPair2 {
public:
  double Chi2;            //!< Reduced chi2
  short int IFit;         //!< Fit type bit pat (KFback:2,QN:8)
  unsigned int NDics;     //!< Number of hits taken into account in the QN fit
  float  Chi2tot;         //!< Total chi2
  THlx   Hfirst;          //!< Track parameters in starting point
  //  double Pinv;                //!< track candidate q/P
  //  double ErrP;                //!< track candidate Sigma(q/P)^2
  std::list<TTrack>::iterator iL;  //!< "upstream"   track segment iterator
  std::list<TTrack>::iterator iR;  //!< "downstream" track segment iterator
  //! "Less" operator
  bool operator < (const TTrackPair2& tp) const { return (Chi2 < tp.Chi2);}
  int WorstPl;            //!< Plane w/ worst chi2 increment.
};

class TEv {
public:
  static unsigned int NobjCreated;
  static unsigned int NobjDestructed;

  TEv();
  ~TEv();

  // Accessors
  static TEv* Ptr();   //!< Returns pointer to this class
  static TEv& Ref();   //!< Returns reference to this class

  CsEvent *ptrEvt();   //!< Get pointer to CsEvent object

  bool LRSelect;       //!< Enable LR selection

  //! Returns "true" in case of MC event
  bool IsMC() const {return isMC;}

  //! Returns Run #
  unsigned int Run()   const {return run;}

  //! Returns Event #
  unsigned int Event() const {return event;}

  //! Returns Event #
  unsigned int EventInBurst() const {return ev_in_burst;}

  //! Returns event trigger mask
  unsigned int TrigMask()       const;
  std::string       TrigMaskString() const;

  //! vector of the TKine objects
  const std::vector<TKine>& vKine()   { return vecKine; }

  //! Accessor to i-th TKine object of the vKine
  const TKine& vKine(int i) const { return vecKine[i]; }

  //! Vector of the TVtxMC objects
  const std::vector<TVtxMC>& vVtxMC() { return vecVtxMC; }

  //! Accessor to i-th TVtxMC object of the vVtxMC
  const TVtxMC& vVtxMC(int i) const { return vecVtxMC[i]; }

  //! Vector of the THitMC objects
  const std::vector<THitMC>& vHitMC() { return vecHitMC; }

  //! Accessor to i-th THitMC object of the vHitMC
  const THitMC& vHitMC(int i) const { return vecHitMC[i]; }

  //! Vector of the THit objects
  const std::vector<THit>&   vHit()   { return vecHit;   }

  //! Accessor to i-th THit object of the vHit
  const THit& vHit(int i) const { return vecHit[i]; }

  //! List of the TTrack objects
  const std::list<TTrack>&   lTrack()   { return listTrack; }

  //! Import CORAL clusters. Depending upon "beamTelescope": all of them (=0), or beam telescope alone (=1), or all but beam telescope (=-1).
  void ImportClusters(const std::list<CsCluster*> &lCsCl,
		      int beamTelescope = 0);

  //! Returns clusters to  CORAL. opt = "all" or "unused"
  std::list<CsCluster*> ExportClusters(std::string opt);

  //! Export found || briged || fitted tracks to CORAL
  void ExportTracks(std::list<CsTrack*>& lCsTrk);

  //! Import CORAL tracks for brigding or fitting. Importing mode = "cp" (copy) or "mv" (move)
  void ImportTracks(std::list<CsTrack*>& lCsTrk, std::string mode);

  //! Simple fit of track segments (to be done before bridging)
  void FitSegments();

  //! Pre-pattern recognition in the CORAL zones
  void PrePattern (CsZone* zone);
  void PrePattern1(CsZone* zone); //!< alternative Pre-pattern
  void PrePattern2(CsZone* zone); //!< alternative Pre-pattern
  void PrePattern3(CsZone* zone); //!< alternative Pre-pattern
  void ParaxialPR();

  //! Bridging, i.e. connecting track segments from different detector groups
  void BridgeSegments ();
  void BridgeSegments1(); //!< Alternative #1
  void BridgeSegments2(); //!< Alternative #2
  void BridgeSegments3(); //!< Alternative #3

  //! Subroutines of alternative bridging #2: Extending track by backtracking
  void BackTrackZ2();               //!< Backtrack 0x4-TTrack's into zone 0x2
  void BackTrackSAS();              //!< Backtrack SAS-TTrack's into zone 0x1
  void BackTrackZ1(int selection);  //!< Backtrack (selected) 0x2-TTrack's into zone 0x1
  void BackTrackVD();               //!< Backtrack &0x2-TTrack's to the DY vertex detector
  bool BackTrack2SI(TTrack &t);     //!< Re-assess SIs in argument track by backtracking its non-SI piece
  //! Subroutines of alternative bridging #2 (cont'd)
  void BridgeOverTarget(std::list<TTrackPair2> &lTrackPair0,
			TTrack &t, int nDFs0x1, double chi20x1);
  bool UpdateKF(int imag, int targetField,
		TTrack &tx, TTrackPair2 &tp,
		int &prvRId, TTrack &prvRT, double &prvCop, double *Hnoise,
		TTrackPair2 *minTP, double &minChi2,
		double &chi2);
  void DoubleBridge();              //!< Bridge over both SM1 and SM2

  //! Extend track by foretracking
  void ForeTrack2RICHWall();  //!< Foretrack 0x3-TTrack's to RICHWall
  int getAssociatedHit(float yp, float zp,
		       float c, float s, const TPlane *p, THit *h,
		       const TPlane *pp, THit **hp);  //!< Called by "ForeTrack2RICHWall"
  void ForeTrack2MA();        //!< Foretrack to muWallA
  int findMWPAHits(TTrack &t, const TStation *sPA, THit **hsPA, int nPAs);  //!< Called by "ForeTrack2MA"
  void ForeTrack2Hs();        //!< Foretrack 0x6-TTrack's to Hodos
  void UpdateHitStatus(int ipli, int iplf, unsigned long long plPat,
		       std::list<unsigned int> &fishyTrks, int &nHits);  //!< Called by "ForeTrack2RICHWall/MA"


  //! Tracks' Fit with quality cuts
  void TracksFit ();
  void TracksFit1(); //!< Alternative track fit
  void TracksFit2(); //!< Alternative track fit
  void BeamsFit2();  //!< Fit beam tracks (beamTelescope piece) w/ P taken from option. In MC, then overwrite P w/ smeared true value.
  void TracksFit3(); //!< Alternative track fit

  //! Connect track segments through muon wall (not used in TraFDic)
  void BridgeMuons();

  //! Dealing w/ yoke tracks: split them and add tail back, after head is modified 
  void SplitYokeTr(TTrack &t, TTrack *&yokeTr, const THlx *&yokeTlHlast);
  void AddYokeTrTail(TTrack &t, TTrack *yokeTr);

  // Reconstruction monitoring
  void Monitor();

  int Quadruples(int ipl, float *yhit, float *zhit, int *yref, int *zref,
		 int ifl[]);

  //! Information consistency checks
  void Checks();

  //! Simple V0 finding
  void SearchV0();

  //! Print Monte-Carlo information (if exists)
  void PrintMC(int mode=0);

  //! Print reconstructed tracks
  void PrintRecTracks(int mode=0);

  //! Clean track list
  void CleanTrackList(int mode = 0);
  void CleanTrackListSI();

  //! Dump for debug purpose
  void DumpEvent();
  void DumpTrackList();
  void DumpHits();

  //! Get array of incident tracks, w/, partial, rejection of those w/ continuation into the spectrometer
  void GetIncidentTracks(int nBeamsMx, int &nBeams, TTrack **beams);
  //! Build ``space point'' consistent w/ TTrack t's Haux based on Z-proj THit hZ
  void BuildSpacePt(TTrack &t, THit *hZ, unsigned int *usedPlPatZ,
		    float dy2, float dyz, float dz2,
		    std::vector<TSpacePt> &spts,
		    const THlx *Haux = 0, bool free = true);

  //! Set "eventTime" from beam-tracks' time, if the latter can be uniquely defined.
  void SetEventTime();
  //! Accessor for "eventTime"
  double GetEventTime() const { return eventTime+eventTRef; }
  //! Update hits from drift-like detectors w/ "TEv::eventTime", prior to track finding in spectrometer, and hence while propagation time is not yet known.
  void UpdateDrifts();
  //! Refit tracks using info from vertexing, viz., so far, event's time. Returns true if any actual refit took place (whether it modified anything or not...)
  bool TracksRefit(CsVertex *pVertex, std::list<CsTrack*> &tracksToBeDeleted);

  //! Flag retracking, after 1st pass of tracking and vertexing.
  void FlagReTracking() { reTracking = true; }

  //! BMS smearing. Used in MC for debugging purposes. Accessor and mutator.
  void   SetBMSSmearing(double value) { BMSSmearing = value; }
  double GetBMSSmearing() const { return BMSSmearing; }

  friend class TDisplay;

private:
  static TEv* address;  //!< pointer to itself
  bool         isMC;        //!< MC data.      Retrieved from parent "CsEvent".
  bool         hadronJob;   //!< Hadron job.   Retrieved from parent "CsEvent".
  unsigned int run;         //!< Run #.        Retrieved from parent "CsEvent".
  unsigned int event;       //!< Event #.      Retrieved from parent "CsEvent".
  unsigned int ev_in_burst; //!< # in burst.   Retrieved from parent "CsEvent".
  unsigned int trig_mask;   //!< Trigger mask. Retrieved from parent "CsEvent".
  double       eventTime;   //!< Event time. Default = 0. Corresponds to the time diff between current T0 used in drifts and the actual (or best estimate of) T0.
  double       eventTRef;   //!< Reference time = back-up of event time after the latter been subtracted from the T0 of drift detectors. Default = 0. Updated via "UpdateDrifts".
  bool    evtTConsidered;   //!< Flag set when event time has been considered, whether it was actually taken into account or not.
  bool         reTracking;  //!< Flags retracking, after 1st pass of tracking and vertexing.
  double       BMSSmearing; //!< Used for debugging MC, cf. "TEv::TracksFit2".

  std::vector<TKine>   vecKine;
  std::vector<TVtxMC>  vecVtxMC;
  std::vector<THitMC>  vecHitMC;
  std::vector<THit>    vecHit;
  std::list  <TTrack>  listTrack;

  std::map<CsCluster*,int,std::less<CsCluster*> > mapCsCl2Hit;  // CsCluster* -> THit index map
  std::map<CsMCHit*,int,std::less<CsMCHit*> > mapCsMCHit2HitMC; // CsMCHit*   -> THitMC index map
  std::map<int,int,std::less<int> > mapCsTrId2TrId;             // CsTrack ID -> TTrack ID map

  void GetMCInfo();        //!< Get current event MC information

  //! Correct cluster positions on MM for magnetic field effects and geometry of the detector
  void CorrectMMClusters(TTrack& track);

  //! Reconstruction quality monitoring histograms (for MC data)
  void MCMonitor();
  //! Monitoring histograms for real data
  void RDMonitor();
};
#endif //TEv_h
