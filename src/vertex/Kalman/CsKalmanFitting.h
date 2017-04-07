/*!
   \file    CsKalmanFitting.h
   \brief   Compass Vertex Parameters Fitting Class.
   \author  Alexandre Korzenev
   \version $Revision: 1.27 $
   \date    $Date: 2010/10/26 13:01:53 $ 

*/

#ifndef CsKalmanFitting_h
#define CsKalmanFitting_h

#include "CsVrtFitting.h"
#include "CsHistograms.h"
#include "THlx.h"
#include "CsVTrack.h"

/*! \class CsKalmanFitting
    \brief Coral interface to Kalman vertex fit.
*/

// ***** MASS CUTS *****
const int V0vvCut = 5, V0ZsCut = 400;
const double V0pTCut = .025, V0cthCut = .9999;
const double M_K   = 0.493677;
const double M_K0  = 0.497672;
const double M_Lam = 1.115684;
const double M_D0  = 1.8645;
const double D0LowCut = -.6,    D0UpCut = .4;
const double DSLowCut = -.0025, DSUpCut = .0038;
const double dm_phi = 0.080, dm_K0 = 0.060, dm_Lam = 0.040;

class CsKalmanFitting : public CsVrtFitting, public CsEndOfJob {

 private:
  
  int    hist_;             //!< 0 - Off, 1 - ON.
 
  int    PrimReduce_;
  int    Schema_;       //!< 0 - IKF, 1 - DKF.
  int    Covar_;
  bool   Agm_;          //!< Track's error matrix augmenting
  bool   Specials_;     //!< 0 - no specials, 1 - specials are kept in fit
  double CosPrimMM_;
  double InvChi2Trk_;   //!< Chi2 cut on track in IKF. 
  double DirChi2Trk_;   //!< Chi2 cut on track in DKF.
  double Chi2VCut_;     //!< Vertex with chi2 greater than "Chi2VCut" is removed.
  double ScalMS_;       //!< Factor used to scale down a MS component of an error matrix.
  double RefPlane_;     //!< Reference plane for spectro tracks (>=1e6 means no reference).
  double RefMargin_;    //!< Reference is varied until vertex<RefPlane-Margin.
  double RefBeam_;      //!< Beam Reference plane (>=1e6 means no reference).
  double CurrentRefPlane_, CurrentRefBeam_;  //!< References for current vertex.
  bool BeamHelixFirst_; //!< Use helix @ 1st measured point for beam track.
  double DistInit_;     //!< if dist between initial and final vert position is bigger - fit redone.
  double ELcorr_;       //!< coefficient for energy losses.
  int    Print_[6];     //!< Print info.
  int    BeamP0_;       //!< Beam momentum reference (for histogramming purposes)
  CsVertex *BestPVertex_;//!< Pointer to Best Primary Vertex. Returns null pointer, if none. N.B.: Only meaningful after "doFitting" has been executed.

  CsHist1D *hPullsSm[5], *hPullsMC[5], *hPullsMCvrt[3];
  CsHist1D *hDeltaKalman[3];
  CsHist1D *hdE, *hm_KK, *hm_Kpi, *hm_pipi, *hm_ppi, *hm_pip, *hm_mumu;
  CsHist1D *hm_pi3m, *hm_pi3e, *hm_pi5m, *hm_pi5e;
  CsHist1D *hSecDeltaKalman[3];
  CsHist1D *hSecPullsMCvrt[3];

  CsHist1F *OM_pipi,*OM_Zvtx,*OM_NTrkM,*OM_BeamP,*OM_BeamPz,*OM_BeamPx,*OM_BeamPy,*OM_BeamTx,*OM_BeamTy;
  CsHist2F *OM_ZvsX,*OM_ZvsY;

  CsHist1F *OM_ybj,*OM_zhad,*OM_mult;

  int statistics_[20];

 private:

  //!< Build the CsVtrack's from args vertex and lists of CsTrack's and specials
  void setTrks(TMtx &XVtx, std::list<CsTrack*> &trks, bool isPrimary, std::map<CsTrack*,bool> &specials, std::list<CsVTrack> &vtrks);
  //!< Reset the CsVTrack's given argument vertex "XVtx"
  void resetTrks(TMtx &XVtx, std::list<CsVTrack> &vtrks);
  void Kalman( TMtx &Xn, TMtx &Cn, std::list<CsVTrack>& vtrks, double& chi2 );
  void GetDerivs(THlx* HelRef, TMtx &Xe,TMtx &Ake, TMtx &Bke, TMtx &Qke, TMtx &Hke, TMtx &COV);
  bool Globalfit(TMtx &CO, TMtx &X0, std::list<CsVTrack> &vVertTrk);
  bool Smoother(TMtx &C0,TMtx &Cn, TMtx &X0,TMtx &Xn, double &Chi2tot, std::list<CsVTrack> &vVertTrk );
  bool InvFilter(TMtx &Cn, TMtx &Xn, CsVTrack &it, double &, TMtx *Xnew,TMtx *Cnew);
  bool CovMatrix( TMtx &Cn, std::list<CsVTrack> &vVertTrk );
  bool DirFilter(TMtx &Cn, TMtx &Xn, CsVTrack &it, double & );
  void ELossCorrection( TMtx& Xn, std::list<CsVTrack>& vVertTrk );
  bool Pulls( std::list<CsVTrack>& vVertTrk, TMtx& Cn, TMtx& Xn );
  bool PullsMC( std::list<CsVTrack>& vVertTrk );
  // Best primary Vertex
  bool   BpVDY_;     //!< Use special DY algorithm to determine Best pVertex
  double BpVTimeW_;  //!< Weight of time chi2 (summed over mu-tracks) in BpV criterion 
  void getBestVertex(CsVertex *vrt, CsVertex *&bestV,
		     int &nonInteractingBV,
		     bool hasMup,      bool &hasMupBV,
		     bool hasBMS,      bool &hasBMSBV,
		     int piNCandidate, int &piNBV);
  void getBestVertexDY(CsVertex *vrt, CsVertex *&bestV,
		       int &nGoodTrksBV, int &nMuswTBV,
		       double &tChi2BV, // Defined upon return only if nMuswTBV> 0
		       bool hasMup,      bool &hasMupBV,
		       bool hasBMS,      bool &hasBMSBV,
		       int piNCandidate, int &piNBV);

 public:

  CsKalmanFitting();              //!< constructor
  ~CsKalmanFitting() {};          //!< destructor

  /*! \fn virtual bool doFitting(list<CsVertex*> &vrts, map<CsTrack*,bool> &specials, bool reTrackingON, double &T0)
      \brief Performs the vertex/track parameters fit for the given lists of vertex/tracks. 
      \param \e vrts = Reference to list of CsVertex objects to be fitted.
      \param \e specials = Reference to map of scattered-mu ID.
      \param \e reTrackingON: true if expecting to be in the context of a 2nd pass of vertexing, after a re-tracking has taken place.
      \param \e T0: If not a null pointer, on input it holds the T0 used in tracking: if this turns out to deviate too much from eventual best vertex' time, and \e reTrackingON is non true, then doFitting aborts and returns false so that a re-tracking be requested. On output, \e T0 holds the Best Vertex time, which is the recommended reference time to use as T0 in the re-tracking.  
   */ 
  bool doFitting(std::list<CsVertex*> &vrts, std::map<CsTrack*,bool> &specials,
		 bool reTrackingON, double *T0);
  double rebuildTracksCut_; //!< Cut on bestVertex beam track time conditioning the request for re-tracking.

  /*! \fn virtual const CsVertex *getBestPVertex() const
    \brief Returns a point to the Best Primary Vertex, if defined.
   */
  const CsVertex *getBestPVertex() const { return BestPVertex_; }

  //! "End of job" method  
  bool end();
};

#endif //CsKalmanFitting_h
