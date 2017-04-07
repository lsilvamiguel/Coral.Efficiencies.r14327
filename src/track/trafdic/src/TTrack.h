// $Id: TTrack.h 14069 2015-09-17 20:44:46Z lsilva $

#ifndef TTrack_h
#define TTrack_h

#include "CsSTD.h"
#include <set>
#include <map>
#include <cassert>
#include "TSetup.h"
#include "THlx.h"
#include "THit.h"
#include "TKine.h"

class CsTrack;

/*!
  \brief Track for tracking :-)

  Main working class of the track finding.
  In the succesive steps of the tracking, its various data members are filled.
  All data members are private. Class TEv is declared "friend" of all tracking
  algorithms in member functions of the class TEv.

  \warning Only for internal use in pattern recognition and track fit.
  \author Sergei.Gerassimov@cern.ch
*/

/*
  Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/TTrack":
   i) Methods: "SubHit", "SubHits", "GetId", "UpdateProjs".
  ii) Member objects "IFit", "NDFs".
 iii) Member objects related to dico fit ("Haux" (auxilliary helix), "vGuests"
  (guestimates) and "GuestsGroups", "Chi2aux") and method "QNewtonFit"
  iv) Associated reconstructed track.
   v) Alternative KF.
  */

class TTrack {

public:
  static unsigned int NobjCreated;
  static unsigned int NobjDestructed;

  TTrack();                 //!< Constructor
  TTrack(const TTrack& t);  //!< Copy constructor
  ~TTrack();                //!< Destructor

  TTrack& operator = (const TTrack& t);   //!< = operators
  
  //        ********** METHODS **********

  void Print(int level, std::string com="") const; //!< Print track contents
  void UpdateProjs();                  //!< Update set of proj. "sProj"

  // Beware! Following methods do not update track's type
  void AddHit(const THit &h);          //!< Add hit "h"
  bool Merge(const TTrack &ti, const TTrack &tj);  //!< Requires no common plane, return false if not
  void AnnexHit(const THit &h);        //!< Add hit "h", while plane "h.IPlane" may already have been inserted
  const THit *ReplaceHit(const THit *h);       //!< Replace pre-existing THit w/ argument one on corresponding plane. If there is indeed a pre-existing THit, a pointer to it is returned. If not, 0 is returned. 
  void RemoveHit(const THit &h);       //!< Remove hit "h"
  void SubHit(const THit &h);          //!< Erase hit "h"
  int  SubHit(int ipl);                //!< Erase hit from plane "ipl", return index of erased hit.
  void SubHits(int ipli, int iplf);    //!< Erase hits from planes in "[ipli,iplf]"
  int  CancelHit(int ipl,              //!< Erase hit from plane "ipl" leaving plane inserted, if already, and then return hit#. The hit->track reference may be left unchanged.
		 bool updateTID = true);
  int  CancelHit_Clip(int ipl,         //!< Same as CancelHit, but clip the list of plane ref. if cancelled hit turns out to be first or last.
		 bool updateTID = true);
  void RestoreHit(int ipl, int ihit,   //!< Restore hit #"ihit" on plane #"ipl". The hit->track reference may be left unchanged. 
		 bool updateTID = true);
  void SubstituteHit(const THit &h1, const THit &h2);  //!< Substitute hit h1 by hit h2 of the same plane
  void Shorten(int ipl);               //!< Erase all hits down to and including plane "ipl"
  void Behead(int ipl);                //!< Erase all hits up   to and including plane "ipl"
  void Clip(bool downstream = false);  //!< Erase planes inserted up(resp. down)stream of first (resp. last) hit
  bool BendingInfo();                  //!< True if >1 hits in bending projections upstream of SM2
  void Evaluate(int ipli, int iplf,    //!< Evaluate # of space points and projections in track, between TPlanes "ipli" and "iplf"
		unsigned int mode,
		int &nSpacePts, int &nAllProjs);
  bool CrossRPipe();                   //!< Return whether track crosses RICH pipe, as determined by straight extrapolation from its first helix.
  void WorstHit(std::map<int,double> &mChi2,//!< Return plane#, zone# and value of worst chi2 increment 
		int &ipl, int &igr, double &incr);
  int  SingleDiff(const TTrack *t,     //!< Check there's a single diff w.r.t. argument TTrack, return its index (or -1), and type of diff.
		  int &typeDiff) const;
  void Append(TTrack &t, std::string mode = "");    //!< Append track "t" to "this" track.
  int  CommonHits(const TTrack &t);            //!< Count number of common hits for track t and "this" track
  void FindKine();                             //!< Find corresponding MC track
  void InsertMissedPlanes();                   //!< Insert missed planes in the list of plane referencies
  bool UseHitTime(int clean = 0);              //!< Track's time. Cleaning upon non null arg, Returns false if track to be erased following cleaning.
  bool QuickKF(int dir, int mode = 0);         //!< Simplified Kalman fit
  bool FullKF (int dir, int ipl = -1,
	       bool scaleUp = true);           //!< Kalman fit w/ Multi-Scattering (detectors' radiation length + material maps or ROOTGeometry), starting w/ plane "ipl", scaling up initial cov. matrix.
  bool KFwExtraMS(double radLenFr);            //!< Alternative Kalman fit, w/ extra MS
  bool Refine (int iter);                      //!< Misc. track improvments
  double GetSmoothed(THlx& h, int ipl, bool f=true) const; //!< get "smoothed" helix on the plane# ipl. Returns Chi2/ndf
  double GetSmoothed(THlx& h, double X) const; //!< Get "smoothed" helix at given X. Returns Chi2/ndf

  bool QNewtonFit(int dir, int mode = 0);      //!< Quasi-Newton fit (w/ propagator or Dico)
  int Fit2IncidentTrack(int nBeams, TTrack **beams, bool mode); //!< QN Fit w/ Y vertex fixed @ intersection w/ incident track in Z dimension

  bool Pickup(const std::list<int> &ipls);     //!< Extrapolate to and pickup from argument planes

  // Accessors
  float RadLenFraction() const { return radLenFr; }  //!< Returns total fraction of radiation length cumulated during last extrapolation.
  float ELoss() const { return eLoss; } //!< Returns total energy loss cumulated during last extrapolation.

  //! returns 1 if this track has a segment in the detectors' group igr, otherwise, retuns 0
  short unsigned int InGroup(int igr) const;

  //! Returns TTrack' type, i.e. bit pattern of zones 
  short unsigned int GetType() const { return(Type); }

  //! returns number of detectors' groups where this track had been found
  short unsigned int NGroups() const;

  //! returns track type bit map (Type) in printable form 
  short unsigned int BitMap() const;

  //! returns pointer to corresponding CsTrack object or NULL
  CsTrack* PtrTrk() const;

  //! returns q/P, taken from the first THlx
  const double& Pinv() const;

  //! Retuns helix (at upstream 'u' or downstream 'd' end)
  const THlx& H(char c) const {
    if(c == 'u') return Hfirst;
    if(c == 'd') return Hlast;
    assert(false);
    return(Hfirst);  // just to get rid of compiler warnings
  }

  //! Returns list of hits
  const std::list<int>& lHPat() const { return lHitPat; }
  //! Returns list of planes
  const std::list<int>& lPRef() const { return lPlnRef; }

  //! Returns Id
  unsigned int GetId() const { return Id; }

  //! Returns associated TKine (i.e. MC) track# "IKine"
  int GetIKine() const { return IKine; }
  //! Returns number of hits from the same MC track "NHsame"
  int GetNHsame() const { return NHsame; }

  //! Dump hit list to screen (for debug purposes).
  void DumpHits(int);

  // START: bg 2006/03/15
  //friend int  TDisplay::DrawTracks(int);
  //friend void TDisplay::TrackInfo (int);
  friend class TDisplay;
  // END: bg 2006/03/15
  
private:

  static unsigned int TrackCounter;    //!< to store unique track id 

  unsigned int Id;    // Unique track ID
  short    int Type;  // Pattern of zones where tracks is found.
  int          Scifi; // Pattern of zones for various attributes: 0xff: scifi-(almost)only, 0xff00: VSAT-(almost)only, 0xff0000: ILMO w/ 1st hodo TStation alone
  short    int IMark; // Flag. If =-1, flags a track (expected to be a single-zone segment) which has been merged w/ (appended to) another track, signalling that it has to be erased (which is done, in principle, once the merging is validated). If >0, flags the segments of a yoke track (i.e. a track through SM2 yoke; both upstream and downstream segments being flagged). The flag is also used in "TEv::BackTrackZ2" to, temporarily, tag tracks that are deemed unreliable, and which hits can be freed.
  short    int IFit;  // Fit type bit pat (Straight:0x1,KFback:0x2,KFfore:0x4,QN:0x8,Very good track:0x10,Full KF:0x60,QN w/ Y-vertex fixed:0x80,beam cont'd into spectro:100) 
  bool fFitDone;      // "true" if forward   FullKF done
  bool bFitDone;      // "true" if bnackward FullKF done
  bool insPlaneDone;  // "true" if "missing" planes have been inserted into hit maps
  int          IKine;   // Associated TKine (i.e. MC) track# (MC only) 
  unsigned int NHits;   // Number of associated hits
  unsigned int NDFs;    // NDFs = #hits + #pixel_hits
  unsigned int NDics;   // Number of associated hits w/in dico
  unsigned int NHsame;  // Number of hits from the same MC track 
  int          IKsame;  // Reference to above mentioned MC track 
  double       Chi2tot; // Total track's chi2
  float  MeanTime;      // Mean time (ns)
  float  SigmaTime;     // Mean time error (<0 if no time measuring hit)
  float  DispTime;      // Time dispersion
  float  radLenFr;      // X/X0 (total fraction of radiation length passed), cumulated during last extrapolation.
  float  eLoss;         // Energy loss, cumulated during last extrapolation.
  CsTrack *ptrTrk;      // Pointer to CsTrack
  int    NShowerHits;   // #hits from preshower (RICHwall) detector
  int    HasShower;     // =1: #hits > cut, =2: &= d#hits/dZ > cut
  int    Associate;     // Associated reconstructed track (defaults to -1). Association (so far) concerns tracks (presumed) non-interacting and continued beyond SM1, i.e. 0x1X (X=3,7,f), and their 0xX and 0x10 pieces: 0x10.Associate = 0x1X.

  THlx Hfirst;          // Helix (i.e. track's parameters) at first measured point
  THlx Hlast;           // Helix (i.e. track's parameters) at last  measured point
  std::list<THlx> lHsm; // Helices - results of Kalman "smoothing"
  THlx Haux;            // Auxilliary storage for 1st helix (used by dico)

  std::set<int> sProj;  // Set of projections the track was built from (Nota bene: this member is not updated by TraFDic) 
  std::set<int> sIHit;  // Set of hit indices

  std::list<int> lPlnRef;    // Plane numbers where track was followed
  std::list<int> lHitPat;    // Hit pattern ( -3 - plane is OFF, -2 - out of active area, -1 - hit not found, >=0 - ref. to vHit. 

  // following datamembers are NOT in copy constructor as
  // it's temporary storage of information for Kalman smoothing.  
  std::map<int,THlx> mHef;   // Helices extrapolated in forward  fit on plane #ipl
  std::map<int,THlx> mHeb;   // Helices extrapolated in backward fit on plane #ipl
  std::map<int,THlx> mHuf;   // Helices updated      in forward  fit on plane #ipl
  std::map<int,THlx> mHub;   // Helices updated      in backward fit on plane #ipl

  std::vector<float> vGuests;// Guestimates from fit
  int GuestsGroups;     // Pattern of groups where guestimates
  double Chi2aux;       // Auxilliary Chi2 (for TLattice)

  std::map<int,double> mChi2;// Chi2 increments from hit on plane #ipl (NOT in copy constructor)

  
  friend class TEv; 
  friend class EqId;

};


inline short unsigned int TTrack::InGroup(int igr) const
{
  const TSetup& setup = TSetup::Ref();
  return( Type>>igr & 1 );
}

inline short unsigned int TTrack::NGroups() const
{
  const TSetup& setup = TSetup::Ref();
  short unsigned int ngr=0;
  for(int igr=0; igr < int(setup.vIplFirst().size()); igr++) ngr += Type>>igr & 1;
  return(ngr);
}

inline short unsigned int TTrack::BitMap() const
{
  short unsigned int ityp(0);
  const TSetup& setup = TSetup::Ref();
  for(int ig = 0; ig < int(setup.vIplFirst().size()); ig++){
    ityp += int(pow(10.,ig)) * (Type>>ig & 1);
  }
  return(ityp);
}

inline CsTrack* TTrack::PtrTrk() const
{
  return( ptrTrk );
}

inline const double& TTrack::Pinv()const
{
  return(Hfirst(5));
}
#endif
