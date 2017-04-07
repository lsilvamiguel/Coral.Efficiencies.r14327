// $Id: CsMCUtils.h,v 1.9 2004/01/13 01:09:35 ybedfer Exp $

/*!
   \file    CsMCUtils.h
   \brief   Compass Montecarlo Utilities Class.
   \author  Benigno Gobbo
   \version $Revision: 1.9 $
   \date    $Date: 2004/01/13 01:09:35 $
   \par History:
   20030401 gcc 3 compliant <br>
   20030401 speed up procedures <br>
*/

#ifndef CsMCUtils_h
#define CsMCUtils_h

#include "CsSTD.h"
#include "CsTrack.h"
#include "CsMCTrack.h"

/*! \class CsMCUtils 
    \brief Monte Carlo Utilities Class.

    Contains some method to easy compare MC True with MC Reco infos.
*/

class CsMCUtils{

 public:

  CsMCUtils();  //!< Constructor

  /*! \fn static std::list<CsMCHits*> getAssociatedHits( const CsCluster* cluster );
     \brief Returns the list of MC hits associated to a cluster
     \param cluster The input cluster
  */
  static std::list<CsMCHit*> getAssociatedHits( const CsCluster* cluster );

  /*! \fn static std::list<CsClusters*> getAssociatedClusters( const CsMCHit* hit );
     \brief Returns the list of clusters associated to a MC hit
     \param hit The input hit
  */
  static std::list<CsCluster*> getAssociatedClusters( const CsMCHit* hit );

  /*! \fn static std::list<CsMCTrack*> getAssociatedMCTracks( const CsTrack* track, const float minHits );
     \brief Returns the list of Monte Carlo True tracks associated to 
            the given track. 
     \param track the Reconstructed track
     \param minHits the % of MC Hits that must be in common between true and
            reconstructed track. Association in dove via: "good" clusters
            -> digits -> "origin" track hits
  */
  static std::list<CsMCTrack*> getAssociatedMCTracks( const CsTrack* track, 
					 const float minHits );  

  /*! \fn static std::list<CsMCTrack*> getAssociatedMCTracksRECON( const CsTrack* track, const float minHits );
     \brief Returns the list of Monte Carlo True tracks associated to 
            the given track. Number of common hits is calculated
            according RECON algorithm.
     \param track the Reconstructed track
     \param minHits the % of MC Hits that must be in common between true and
            reconstructed track. Association in dove via: "good" clusters
            -> digits -> "origin" track hits
  */
  static std::list<CsMCTrack*> getAssociatedMCTracksRECON( const CsTrack* track, 
					 const float minHits );  

  /*! \fn static CsMCTrack* getAssociatedMCTrack( const CsTrack* track, int& nhits )
     \brief Returns the Monte Carlo True track with the maximum number
            of hits in common with track.
     \param track the Reconstructed track
     \param nhits Return value: the number of hits in common
     \warning Return pointer could be null!
  */
  static CsMCTrack* getAssociatedMCTrack( const CsTrack* track, int& nhits );

  /*! \fn static CsMCTrack* getAssociatedMCTrackRECON( const CsTrack* track, int& nhits )
     \brief Returns the Monte Carlo True track with the maximum number
            of hits in common with track. Number of common hits is calculated
            according RECON algorithm.
     \param track the Reconstructed track
     \param nhits Return value: the number of hits in common
     \warning Return pointer could be null!
  */
  static CsMCTrack* getAssociatedMCTrackRECON( const CsTrack* track, int& nhits );

  /*! \fn static std::list<CsTrack*> getAssociatedTracks( const CsMCTrack* mctrack, const float minHits );
     \brief Returns the list of Reconstructed tracks associated to 
            the given Monte Carlo True track. 
     \param mcTrack The Monte Carlo true Track 
     \param minHits the % of MC Hits that must be in common between true and
            reconstructed tracks. Association in dove via: "good" clusters
            -> digits -> "origin" track hits
  */
  static std::list<CsTrack*> getAssociatedTracks( const CsMCTrack* mcTrack, 
				      const float minHits );  

  /*! \fn static CsTrack* getAssociatedTrack( const CsMCTrack* mcTrack, int& nhits )
     \brief Returns the Reconstructed track with the maximum number
            of hits in common with mcTrack.
     \param mcTrack The Monte Carlo true Track 
     \param nhits Return value: the number of hits in common
     \warning Return pointer could be null!
  */
  static CsTrack* getAssociatedTrack( const CsMCTrack* mcTrack, int& nhits );

  /*! \fn static std::list<CsTrack*> getAffiliatedTracks( const CsMCTrack* mctrack, const float minClusters );
     \brief Returns a list of Reconstructed tracks associated to
            the given Monte Carlo True track. Nota Bene: Association criterion
	    differs from that of "getAssociatedTracks". (This particular
	    function has been introduced to meet the requirements of
	    phast (At least for phast's version <= 6.008))
     \param mcTrack The Monte Carlo true Track 
     \param minHits the % of reco'd track's Clusters that must be in common
            between MC and reco'd tracks. Association in dove via:
	    "good" clusters -> digits -> "origin" track hits
  */
  static std::list<CsTrack*> getAffiliatedTracks( const CsMCTrack* mcTrack, 
						  const float minClusters);

  /*! \fn static bool isAGoodCluster( const CsCluster* cluster )
     \brief Returns \c true if the cluster is a real signal cluster;
            returns \c false if the cluster is from noise digits or 
            if it is a mirror cluster in drift-like detectors.   
     \param cluster The Cluster
  */
  static bool isAGoodCluster( const CsCluster* cluster );
  
  /*! \fn static bool isACorrelatedTrack( const CsMCTrack* mcTrack )
     \brief returns \c true if the track is correlated in time with the 
            trigger, returns \c false if the track originates from pileup
     \param mcTrack The track
  */
  static bool isACorrelatedTrack( const CsMCTrack* mcTrack );

  /*! \fn static bool isACorrelatedHit( const CsMCHit* mcHit )
     \brief returns \c true if the hit belong to a track correlated 
            in time with the trigger, 
     	    returns \c false if the hit is due to a pileup track
     \param mcHit The hit
  */
  static bool isACorrelatedHit( const CsMCHit* mcHit );

 private:

  static unsigned int _currentEvent;
  static std::multimap<const CsMCHit*, const CsCluster*> _hc; 
  static void _buildHCMultimap( void );

};

#endif // CsMCUtils_h
