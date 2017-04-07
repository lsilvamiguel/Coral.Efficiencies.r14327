/*!
   \file    CsAverPattern.h
   \brief   Compass Vertex Pattern Class.
   \author  Alexandre Korzenev
   \version $Revision: 1.22 $
   \date    $Date: 2010/11/11 01:09:30 $ 

*/

#ifndef CsAverPattern_h
#define CsAverPattern_h



/*! \class CsVtrPattern 
    \brief Coral interface to Kalman Vertex Reconstruction Pattern.
*/

#include "CsVrtPattern.h"
#include "CsHistograms.h"
#include "THlx.h"

class CsAverPattern : public CsVrtPattern, public CsEndOfJob {
 
 private:
  
  bool  Print_[6];          //!< Print info.
  int    hist_;             //!< 0 - Off, 1 - ON.
  /*! \fn bool mode_;
    \brief For debugging.
    
    \arg \a 0 : standard
    \arg \a 1 : only correct beam track is involved (for MC only)
    \arg \a 2 : all involved tracks are correct (for MC only)

    The \a "1" and \a "2" modes serves for debugging: one uses a set of tracks 
    created in primary vertex (it is known from MC association) and
    just does a fit. If correct MC association for beam track is not 
    found an event is skiped.
  */
  int    mode_;
  bool   NSpec_;          //!< true - skip ev if mu' is not found.
    
  double CosPrimMM_;        //!< Cut on cos(mu,mu'): below cut, the initial guess for Z of vertex is derived from mu' only.
  double CosPrim_;          //!< Cut on cos(beam,track): above cut, Z of CDA is not checked.
  double AcceptTill_;       //!< The track is accepted for Kalman Fit if its momentum is smaller than \a AcceptTill_.
 
  double MomCut_;           //!< Cut for accepted momentum.
  double BeamCand_;         //!< Do loop over all beam candidates inside of interval \a BeamCand_.
  
  double LinDist_;          //!< Cut on \e z in case of liniar extrapolation in prefilter.
  double HelDist_;          //!< Cut on distance to the beam track in case of helix extrapolation.
 
  double CnSigmas_[3];      //!< Diagonal elements of covariance matrix

  double TimePrimCut_;  //!< Cut on time for tracks in primary vertex
  double TimeSecCut_;   //!< Cut on time difference between 2 tracks

  bool findPrim_;           //!< true - finding of primary vertex is ON, false - OFF.
  bool DDistUseMMap_;       //!< Uncertainty on helix Distance (in FindPrimary) takes material (be it from map or ROOTG) into account
  bool findSec_;            //!< true - finding of secondaries is ON, false - OFF.
  bool refitTracks_;        //!< true - Refitting tracks is ON, false - OFF.
  bool rebuildTracks_;      //!< true - Requesting a re-tracking is ON, false - OFF.
  double rebuildTracksCut_; //!< Cut on the beam track time conditioning the requesting a re-tracking.

  double SecDist_;          //!< accepted distance between two tracks of different signed.
  double SecXX0_;           //!< The number of radiation lengths allowed by tracks in V0.

  CsVertex *T0SettingVertex_; //!< Pointer to preliminary Best Primary Vertex, used to set the T0 of tracking, if any. N.B.: Only meaningful after "doPattern" has been executed.

  
  CsHist1D *hDeltaPref[3];
  CsHist1D *hSecDeltaPref[3];

  std::list<CsTrack*>  tracks_;
  std::list<CsVertex*> vrts_;
  std::map<CsTrack*,bool> specials_; //!< Is needed to select special tracks for primary vertex.

  int statistics_[20];
  
 private:

  /*! \fn bool FindPrimary( THlx** HTrMom, int* vIact, CsTrack** TrkRef, int tn ); 
    \brief Fills vrts_ with primary vertex. 
  */
  bool FindPrimary( THlx** HTrMom, int* vIact, CsTrack** TrkRef, int tn );

  /*! \fn bool FindSecondaries( THlx** HTrMom, int* vIact, CsTrack** TrkRef, int tn );
    \brief Fills vrts_ with secondary vertices. 
  */
  bool FindSecondaries( THlx** HTrMom, int* vIact, CsTrack** TrkRef, int tn );

 
 public:

  CsAverPattern();          //!< constructor
  ~CsAverPattern();         //!< destructor

  /*! \fn bool doPattern( vector<CsParticle*> &particles ) = 0;
    \brief This method performs the vertex pattern recognition on the argument set of particles. Returns \c true if the operation ended correctly; returns \c false otherwise.
    \param particles = Vector of pointers to the particles.
    \param reTrackT0 = if non null on input, allows "doPattern" to request a re-tracking, by returning on output the event time to be used in the re-tracking.
  */
  bool doPattern(std::vector<CsParticle*> &particles, double *reTrackT0 = 0);

  /*! \fn bool getPatterns( list<CsVertex*>& vrts, map<CsTrack*,bool>& specials )
      \brief This method returns the list of found vertex candidates after the pattern recognition.
      \param vrts = List of vertex candidates.
      \param specials = Map of special tracks (incident, scattered mu...).
   */
  bool getPatterns( std::list<CsVertex*>& vrts, std::map<CsTrack*,bool>& specials );

  /*! \fn virtual const CsVertex *getT0SettingVertex() const
    \brief Returns a pointer to the vertex selected by the PR as the potential Best Primary Vertex and hence assigning the T0 used in a possible refit of all tracks. 
   */
  const CsVertex *getT0SettingVertex() const { return T0SettingVertex_; }

  //! "End of job" method 
  bool end();
};

#endif //CsAverPattern_h
