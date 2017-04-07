// $Id: CsMCUtils.cc,v 1.23 2010/06/02 10:56:12 tnagel Exp $

/*!
   \file    CsMCUtils.cc
   \brief   Compass Montecarlo Utilities Class.
   \author  Benigno Gobbo
   \version $Revision: 1.23 $
   \date    $Date: 2010/06/02 10:56:12 $
*/

#include "CsMCUtils.h"

#include "CsEvent.h"
#include "CsGeant3.h"
#include "CsErrLog.h"
#include "CsCluster.h"
#include "CsDigit.h"
#include "CsMCDigit.h"
#include "CsMCHit.h"
#include "CsMCTrkHit.h"
#include "CsDetector.h"
#include "CsOpt.h"

#include <CLHEP/Matrix/Matrix.h>

# include <algorithm>

using namespace std;
using namespace CLHEP;

unsigned int CsMCUtils::_currentEvent = 0;
std::multimap<const CsMCHit*, const CsCluster*> CsMCUtils::_hc;


CsMCUtils::CsMCUtils() {
}

std::list<CsMCTrack*> CsMCUtils::getAssociatedMCTracks( const CsTrack* track, 
							const float minHits ) {

  std::list<CsMCTrack*> mcTracks; mcTracks.clear();

  // check if MC run...
  if( ! (CsEvent::Instance()->isAMonteCarloEvent()) ) {
    // MC data only...
    CsErrLog::mes( elError, "This method works only on MC data." );
    return mcTracks;
  }

  // find all hits associated to the track
  CsTrack theTrack = *track;
  std::list<CsMCHit*>   trackHits, clusterHits;
  std::list<CsCluster*> clusters = theTrack.getClusters(); // track clusters
  std::list<CsCluster*>::iterator Ic;
  for( Ic=clusters.begin(); Ic!=clusters.end(); Ic++ ) {
    if( isAGoodCluster( *Ic ) ) {
      clusterHits.clear();
      std::list<CsDigit*> digits = (*Ic)->getDigitsList(); // cluster digits
      std::list<CsDigit*>::iterator Id;
      for( Id=digits.begin(); Id!=digits.end(); Id++ ) {
	CsMCDigit* dig = dynamic_cast<CsMCDigit*>( *Id );
	if( dig != 0 ) {
	  std::list<CsMCHit*> hits = dig->getHits(); // digit hits 
	  if( !hits.empty() ) {
	    std::list<CsMCHit*>::iterator Ih; 
	    for( Ih=hits.begin(); Ih!=hits.end(); Ih++ ) {
	      if( (*Ih)->getOrigin() == 0 ) {
		// if two or more digits are originated from one MC hit.
		if( std::find(clusterHits.begin(),clusterHits.end(),*(Ih)) == clusterHits.end())
		  clusterHits.push_back( *(Ih) );
	      }
	    }
	  }
	}
      }
      if(clusterHits.size()==1) trackHits.push_back(clusterHits.front());
    }
  }

  // find all MC tracks with more than minHits % hits in common with track
  int NtrackHits = theTrack.getClusters().size();

  std::list<CsMCTrack*> allMCTracks = CsGeant3::Instance()->getMCTracks();
  std::list<CsMCTrack*>::iterator It; 

  for( It=allMCTracks.begin(); It!=allMCTracks.end(); It++ ) {
    int NcommonHits = 0;
    std::list<CsMCHit*> hits = (*It)->getMCHits();
    std::list<CsMCHit*>::iterator Ih; 
    for( Ih=hits.begin(); Ih!=hits.end(); Ih++ ) {
      if( std::find( trackHits.begin(),trackHits.end(),*Ih ) != trackHits.end() ) {
	NcommonHits++;
      }      
    }
    if( (float(NcommonHits)/float(NtrackHits)*100) > minHits ) {
      mcTracks.push_back( *(It) );
    }
  }

  return( mcTracks );

}

std::list<CsMCTrack*> CsMCUtils::getAssociatedMCTracksRECON( const CsTrack* track, 
						   const float minHits ) {

  std::list<CsMCTrack*> mcTracks;
  mcTracks.clear();

  // check if MC run...
  if( ! (CsEvent::Instance()->isAMonteCarloEvent()) ) {
    // MC data only...
    CsErrLog::mes( elError, "This method works only on MC data." );
    return mcTracks;
  }

  // find all hits associated to the track
  CsTrack theTrack = *track;
  std::list<CsMCHit*>   trackHits, clusterHits;
  std::list<CsCluster*> clusters = theTrack.getClusters(); // track clusters
  std::list<CsCluster*>::iterator Ic;
  for( Ic=clusters.begin(); Ic!=clusters.end(); Ic++ ) {
    
    clusterHits.clear();
    std::list<CsDigit*> digits = (*Ic)->getDigitsList(); // cluster digits
    std::list<CsDigit*>::iterator Id;
    for( Id=digits.begin(); Id!=digits.end(); Id++ ) {
      CsMCDigit* dig = dynamic_cast<CsMCDigit*>( *Id );
      if( dig != 0 ) {
	std::list<CsMCHit*> hits = dig->getHits(); // digit hits 
	if( !hits.empty() ) {
	  std::list<CsMCHit*>::iterator Ih; 
	  for( Ih=hits.begin(); Ih!=hits.end(); Ih++ ) {
	    if( (*Ih)->getOrigin() == 0 ) {
	      // if two or more digits are originated from one MC hit.
	      if( std::find(clusterHits.begin(),clusterHits.end(),*(Ih)) == clusterHits.end())
		clusterHits.push_back( *(Ih) );
	    }
	  }
	}
      }
    }
    std::list<CsMCHit*>::iterator Ih; 
    for( Ih=clusterHits.begin(); Ih!=clusterHits.end(); Ih++ ) 
      trackHits.push_back( *(Ih) );
    
  }
  
  // find all MC tracks with more than minHits % hits in common with track
  int NtrackHits = theTrack.getClusters().size();

  std::list<CsMCTrack*> allMCTracks = CsGeant3::Instance()->getMCTracks();
  std::list<CsMCTrack*>::iterator It; 

  for( It=allMCTracks.begin(); It!=allMCTracks.end(); It++ ) {
    int NcommonHits = 0;
    std::list<CsMCHit*> hits = (*It)->getMCHits();
    std::list<CsMCHit*>::iterator Ih; 
    for( Ih=hits.begin(); Ih!=hits.end(); Ih++ ) {
      if( std::find( trackHits.begin(),trackHits.end(),*Ih ) != trackHits.end() ) {
	NcommonHits++;
      }      
    }
    if( (float(NcommonHits)/float(NtrackHits)*100) > minHits ) {
      mcTracks.push_back( *(It) );
    }
  }

  return( mcTracks );

}

CsMCTrack* CsMCUtils::getAssociatedMCTrack( const CsTrack* track, int& nhits ) {

  CsMCTrack* TheMCTrack = 0;

  // check if MC run...
  if( ! (CsEvent::Instance()->isAMonteCarloEvent()) ) {
    // MC data only...
    CsErrLog::mes( elError, "This method works only on MC data." );
    return TheMCTrack;
  }

  // find all hits associated to the track
  CsTrack theTrack = *track;
  std::list<CsMCHit*>   trackHits, clusterHits;
  std::list<CsCluster*> clusters = theTrack.getClusters(); // track clusters
  std::list<CsCluster*>::iterator Ic;
  for( Ic=clusters.begin(); Ic!=clusters.end(); Ic++ ) {
    if( isAGoodCluster( *Ic ) ) {
      clusterHits.clear();
      std::list<CsDigit*> digits = (*Ic)->getDigitsList(); // cluster digits
      std::list<CsDigit*>::iterator Id;
      for( Id=digits.begin(); Id!=digits.end(); Id++ ) {
	CsMCDigit* dig = dynamic_cast<CsMCDigit*>( *Id );
	if( dig != 0 ) {
	  std::list<CsMCHit*> hits = dig->getHits(); // digit hits 
	  if( !hits.empty() ) {
	    std::list<CsMCHit*>::iterator Ih; 
	    for( Ih=hits.begin(); Ih!=hits.end(); Ih++ ) {
	      if( (*Ih)->getOrigin() == 0 ) {
		// if two or more digits are originated from one MC hit.
		if( std::find(clusterHits.begin(),clusterHits.end(),*(Ih)) == clusterHits.end())
		  clusterHits.push_back( *(Ih) );
	      }
	    }
	  }
	}
      }
      if(clusterHits.size()==1) trackHits.push_back(clusterHits.front());
    }
  }

  // find the MC track with the maximum number of hits in common with track
  int MaxHits = 0; 

  std::list<CsMCTrack*> allMCTracks = CsGeant3::Instance()->getMCTracks();
  std::list<CsMCTrack*>::iterator It; 

  for( It=allMCTracks.begin(); It!=allMCTracks.end(); It++ ) {
    int NcommonHits = 0;
    std::list<CsMCHit*> hits = (*It)->getMCHits();
    std::list<CsMCHit*>::iterator Ih; 
    for( Ih=hits.begin(); Ih!=hits.end(); Ih++ ) {
      if( std::find( trackHits.begin(),trackHits.end(),*Ih ) != trackHits.end() ) {
	NcommonHits++;
      }      
    }
    if( NcommonHits > MaxHits ) {
      MaxHits = NcommonHits;
      TheMCTrack =  *It;
    }
  }

  nhits = MaxHits;
  return( TheMCTrack );

}


CsMCTrack* CsMCUtils::getAssociatedMCTrackRECON( const CsTrack* track, int& nhits ) {

  CsMCTrack* TheMCTrack = 0;

  // check if MC run...
  if( ! (CsEvent::Instance()->isAMonteCarloEvent()) ) {
    // MC data only...
    CsErrLog::mes( elError, "This method works only on MC data." );
    return TheMCTrack;
  }

  // find all hits associated to the track
  CsTrack theTrack = *track;
  std::list<CsMCHit*>   trackHits, clusterHits;
  std::list<CsCluster*> clusters = theTrack.getClusters(); // track clusters
  std::list<CsCluster*>::iterator Ic;
  for( Ic=clusters.begin(); Ic!=clusters.end(); Ic++ ) {

    clusterHits.clear();
    std::list<CsDigit*> digits = (*Ic)->getDigitsList(); // cluster digits
    std::list<CsDigit*>::iterator Id;
    for( Id=digits.begin(); Id!=digits.end(); Id++ ) {
      CsMCDigit* dig = dynamic_cast<CsMCDigit*>( *Id );
      if( dig != 0 ) {
	std::list<CsMCHit*> hits = dig->getHits(); // digit hits 
	if( !hits.empty() ) {
	  std::list<CsMCHit*>::iterator Ih; 
	  for( Ih=hits.begin(); Ih!=hits.end(); Ih++ ) {
	    if( (*Ih)->getOrigin() == 0 ) {
	      // if two or more digits are originated from one MC hit.
	      if( std::find(clusterHits.begin(),clusterHits.end(),*(Ih)) == clusterHits.end())
		clusterHits.push_back( *(Ih) );
	    }
	  }
	}
      }
    }
    std::list<CsMCHit*>::iterator Ih; 
    for( Ih=clusterHits.begin(); Ih!=clusterHits.end(); Ih++ ) 
      trackHits.push_back( *(Ih) );
  }

  // find the MC track with the maximum number of hits in common with track
  int MaxHits = 0; 

  std::list<CsMCTrack*> allMCTracks = CsGeant3::Instance()->getMCTracks();
  std::list<CsMCTrack*>::iterator It; 

  for( It=allMCTracks.begin(); It!=allMCTracks.end(); It++ ) {
    int NcommonHits = 0;
    std::list<CsMCHit*> hits = (*It)->getMCHits();
    std::list<CsMCHit*>::iterator Ih; 
    for( Ih=hits.begin(); Ih!=hits.end(); Ih++ ) {
      if( std::find( trackHits.begin(),trackHits.end(),*Ih ) != trackHits.end() ) {
	NcommonHits++;
      }      
    }
    if( NcommonHits > MaxHits ) {
      MaxHits = NcommonHits;
      TheMCTrack =  *It;
    }
  }

  nhits = MaxHits;
  return( TheMCTrack );

}


std::list<CsTrack*> CsMCUtils::getAssociatedTracks( const CsMCTrack* mcTrack, 
					       const float minHits ) {
  std::list<CsTrack*> tracks;
  tracks.clear();

  // check if MC run...
  if( ! (CsEvent::Instance()->isAMonteCarloEvent()) ) {
    // MC data only...
    CsErrLog::mes( elError, "This method works only on MC data." );
    return tracks;
  }

  // get all hits associated to the MC track
  CsMCTrack theMCTrack = *mcTrack;
  std::list<CsMCHit*> mctrackHits = theMCTrack.getMCHits();

  // find all tracks with more than minHits % hits in common with MC track
  int NmctrackHits = mctrackHits.size();

  std::list<CsTrack*> allTracks = CsEvent::Instance()->getTracks();
  std::list<CsTrack*>::iterator It; 

  for( It=allTracks.begin(); It!=allTracks.end(); It++ ) {
    int NcommonHits = 0;
    std::list<CsCluster*> clusters = (*It)->getClusters(); // track clusters
    std::list<CsCluster*>::iterator Ic;
    for( Ic=clusters.begin(); Ic!=clusters.end(); Ic++ ) {
      if( isAGoodCluster( *Ic ) ) {
	std::list<CsDigit*> digits = (*Ic)->getDigitsList(); // cluster digits
	std::list<CsDigit*>::iterator Id;
	for( Id=digits.begin(); Id!=digits.end(); Id++ ) {
	  CsMCDigit* dig = dynamic_cast<CsMCDigit*>( *Id );
	  if( dig != 0 ) {
	    std::list<CsMCHit*> hits = dig->getHits(); // digit hits 
	    if( !hits.empty() ) {
	      std::list<CsMCHit*>::iterator Ih; 
	      for( Ih=hits.begin(); Ih!=hits.end(); Ih++ ) {
		if( (*Ih)->getOrigin() == 0 && 
		   std::find( mctrackHits.begin(), mctrackHits.end(), *Ih ) 
		    != mctrackHits.end() ) {
		  NcommonHits++;
		}      
	      }
	    }
	  }
	}
      }
    }
    if( (float(NcommonHits)/float(NmctrackHits)*100) > minHits ) {
      tracks.push_back( *(It) );
    }
  }

  return( tracks );

}


CsTrack* CsMCUtils::getAssociatedTrack( const CsMCTrack* mcTrack, int& nhits ) {

  CsTrack* TheTrack = 0;
  // check if MC run...

  if( ! (CsEvent::Instance()->isAMonteCarloEvent()) ) {
    // MC data only...
    CsErrLog::mes( elError, "This method works only on MC data." );
    return TheTrack;
  }

  // get all hits associated to the MC track
  CsMCTrack theMCTrack = *mcTrack;
  std::list<CsMCHit*> mctrackHits = theMCTrack.getMCHits();

  // find the track withthe maximum number of hits in common with the MC track
  int MaxHits = 0;

  std::list<CsTrack*> allTracks = CsEvent::Instance()->getTracks();
  std::list<CsTrack*>::iterator It; 

  for( It=allTracks.begin(); It!=allTracks.end(); It++ ) {
    int NcommonHits = 0;
    std::list<CsCluster*> clusters = (*It)->getClusters(); // track clusters
    std::list<CsCluster*>::iterator Ic;
    for( Ic=clusters.begin(); Ic!=clusters.end(); Ic++ ) {
      if( isAGoodCluster( *Ic ) ) {
	std::list<CsDigit*> digits = (*Ic)->getDigitsList(); // cluster digits
	std::list<CsDigit*>::iterator Id;
	for( Id=digits.begin(); Id!=digits.end(); Id++ ) {
	  CsMCDigit* dig = dynamic_cast<CsMCDigit*>( *Id );
	  if( dig != 0 ) {
	    std::list<CsMCHit*> hits = dig->getHits(); // digit hits 
	    if( !hits.empty() ) {
	      std::list<CsMCHit*>::iterator Ih; 
	      for( Ih=hits.begin(); Ih!=hits.end(); Ih++ ) {
		if( (*Ih)->getOrigin() == 0 && 
		    std::find( mctrackHits.begin(), mctrackHits.end(), *Ih ) 
		    != mctrackHits.end() ) {
		  NcommonHits++;
		}      
	      }
	    }
	  }
	}
      }
    }
    if( NcommonHits > MaxHits ) {
      MaxHits = NcommonHits;
      TheTrack = *It;
    }
  }

  nhits = MaxHits;
  return( TheTrack );

}

std::list<CsTrack*> CsMCUtils::getAffiliatedTracks( const CsMCTrack* mcTrack, 
						    const float minClusters) {
  std::list<CsTrack*> tracks;
  tracks.clear();

  // check if MC run...
  if( ! (CsEvent::Instance()->isAMonteCarloEvent()) ) {
    // MC data only...
    CsErrLog::mes( elError, "This method works only on MC data." );
    return tracks;
  }

  // get all hits associated to the MC track
  CsMCTrack theMCTrack = *mcTrack;
  std::list<CsMCHit*> mctrackHits = theMCTrack.getMCHits();

  // find all tracks with more than minClusters % of their clusters in common with MC

  std::list<CsTrack*> allTracks = CsEvent::Instance()->getTracks();
  std::list<CsTrack*>::iterator It; 

  for( It=allTracks.begin(); It!=allTracks.end(); It++ ) {
    int NcommonClusters = 0;
    std::list<CsCluster*> clusters = (*It)->getClusters(); // track clusters
    std::list<CsCluster*>::iterator Ic;
    for( Ic=clusters.begin(); Ic!=clusters.end(); Ic++ ) {
      if( isAGoodCluster( *Ic ) ) {
	std::list<CsDigit*> digits = (*Ic)->getDigitsList(); // cluster digits
	std::list<CsDigit*>::iterator Id;
	bool hasHitinCommon = false;
	for( Id=digits.begin(); Id!=digits.end(); Id++ ) {
	  CsMCDigit* dig = dynamic_cast<CsMCDigit*>( *Id );
	  if( dig != 0 ) {
	    std::list<CsMCHit*> hits = dig->getHits(); // digit hits 
	    if( !hits.empty() ) {
	      std::list<CsMCHit*>::iterator Ih; 
	      for( Ih=hits.begin(); Ih!=hits.end(); Ih++ ) {
		if( (*Ih)->getOrigin() == 0 && 
		   std::find( mctrackHits.begin(), mctrackHits.end(), *Ih )
		    != mctrackHits.end() ) {
		  hasHitinCommon = true; break;
		}
	      }
	    }
	  }
	}
	if (hasHitinCommon) NcommonClusters++;
      }
    }
    if( (float(NcommonClusters)/clusters.size()*100) > minClusters ) {
      tracks.push_back( *(It) );
    }
  }

  return( tracks );

}

std::list<CsCluster*> CsMCUtils::getAssociatedClusters( const CsMCHit* hit ) {

  std::list<CsCluster*> myClusters;
  myClusters.clear();

  // check if MC run...
  if( ! (CsEvent::Instance()->isAMonteCarloEvent()) ) {
    // MC data only...
    CsErrLog::mes( elError, "This method works only on MC data." );
    return myClusters;
  }

  CsMCUtils::_buildHCMultimap();

  for( std::pair<std::multimap<const CsMCHit*,const CsCluster*>::iterator,
	 std::multimap<const CsMCHit*,const CsCluster*>::iterator> Im = 
	 _hc.equal_range( hit ); Im.first != Im.second; ++Im.first ) {
    myClusters.push_back( const_cast<CsCluster*>((*Im.first).second) );
  }

  return( myClusters );

}

std::list<CsMCHit*> CsMCUtils::getAssociatedHits( const CsCluster* cluster ) {

  std::list<CsMCHit*> myHits;
  myHits.clear();

  // check if MC run...
  if( ! (CsEvent::Instance()->isAMonteCarloEvent()) ) {
    // MC data only...
    CsErrLog::mes( elError, "This method works only on MC data." );
    return myHits;
  }

  const std::list<CsDigit*> digits = cluster->getDigitsList(); // cluster digits
  std::list<CsDigit*>::const_iterator Id;
  for( Id=digits.begin(); Id!=digits.end(); Id++ ) {
    CsMCDigit* dig = dynamic_cast<CsMCDigit*>( *Id );
    if( dig != 0 ) {
      std::list<CsMCHit*> hits = dig->getHits(); // digit hits 
      if( !hits.empty() ) {
	std::list<CsMCHit*>::iterator Ih; 
	for( Ih=hits.begin(); Ih!=hits.end(); Ih++ ) {
	  if( std::find( myHits.begin(), myHits.end(), *Ih ) == myHits.end() ) {
	    myHits.push_back( *Ih );
	  }      
	}
      }
    }
  }
  
  return( myHits );

}


bool CsMCUtils::isAGoodCluster(const CsCluster* cluster) {

  //  "Good" cluster means: not the ghost mirror in a pair of left/right
  // CsClusters from a drift detector.
  //  In all other circumstances (not a drift, drift not processed as such, not
  // clusterised as such), the method returns true.
  //  The algotithm compares the positions of the "CsCluster" and of the MC hit
  // it originates from, w.r.t. the central wire of the drift cell.
  //  !! WARNING !!: This algorithm is flawed: it does not provide for the MC
  // hits that generate several CsCluster's in the several drift cells they
  // traverse via the "multiHitDecoding" option of the various "makeMCDecoding"
  // methods. A more accurate version of it would use the additional info
  // in CsMCHit's CsMCTrkHit sub-class. But it would be quite complex.
  //  The best way to get the genuine vs. ghost info for drift clusters is to
  // access the "isGenuine" attribute of CsCluster. N.B.: this attribute is
  // only defined in the particular case of MC drfit CsCluster's. In other
  // cases it's undetermined.

  //            *************** METHODE? ***************
  static bool first = true; static int methode = 0; // Default = 0
  if (first) {
    CsOpt::Instance()->getOpt("MCUtils","isAGood",methode); first = false; 
  }

  //         *************** CHECK THAT MC DATA ***************
  if (!(CsEvent::Instance()->isAMonteCarloEvent())) {
    CsErrLog::mes( elError, "This method works only on MC data." );
    return false;
  }
  // ********** IF "Exact" OR "Smeared" CLUSTERISATION MODE... **********
  string clusterMakerType = CsEvent::Instance()->getClusterizationType();
  if (clusterMakerType=="MCExact" || clusterMakerType=="MCSmeared") {
    return true;                                        // ...RETURN ALWAYS TRUE
  }

  //      ********** GET CORRESPONDING CELL# (whenceforth wire#) **********
  // (Retrieving the cell# from the CsCluster's U and associated CsDetector's
  // pitch and offset, is not straightforward: may depend upon the details of
  // the CsDetector irregular spacing, cf. CsRichWallDetector.
  //  => get cell# from associated CsDigit.)
  const std::list<CsDigit*> digits = cluster->getDigitsList();
  if (digits.empty()) // This should never happen...
    return false;     // ...if it does though...
  // ...else
  int wire = digits.front()->getAddress();
  std::list<CsMCHit*> hits = getAssociatedHits( cluster );  // ***** GET MC HITS
  static double umirror; if (methode) {              // ***** GET MIRROR CLUSTER
    if (cluster->hasMirrorCluster()) {
      umirror = (cluster->getMirrorCluster())->getU();
    }
    else return true;
  }
  double uclus = cluster->getU();                     // ***** GET CsCluster's U

  std::list<CsMCHit*>::iterator Ih;
  for (Ih = hits.begin(); Ih!=hits.end(); Ih++) {
  
    const CsDetector* det = dynamic_cast<CsDetector*>((*Ih)->getDet());
    if (det==NULL) continue;
    if (det->hasDrift()) {

      HepMatrix rotM(3,3); rotM = det->getRotWRS();
      double wirD               = det->getWirD();
      double wirP               = det->getWirP();

      double x = (*Ih)->getX(), y = (*Ih)->getY(), z = (*Ih)->getZ();

      int err; HepMatrix iRotM(3,3); iRotM = rotM.inverse( err );
      double uhit  = iRotM(1,1)*x + iRotM(1,2)*y + iRotM(1,3)*z;
      double uwire = wirD + det->Wire2Pos(wire);
  
      if (methode) {// ***** ALTERNATIVE METHOD allowing for COALESCENCE *****
	// (- In close enough mirrors pairs, both mirrors can be considered
	//    genuine.
	//  - TraFDic is itself coalescing such pairs into a single cluster
	//   (a "THit" in TraFFic terms), cf. lattice's "TEv::ImportClusters",
	//   to which the genuine mirror is associated. Therefore, when TraFDic
	//   is in use, the present method may not be so usefull.)
        static float coalesce = 4.;
	HepMatrix ccov = cluster->getCov();
	float sigu = sqrt(ccov(1,1));	
	// do not check if mirror cluster exists, done above
	if( fabs( uclus - uhit ) <= fabs( umirror - uhit ) || 
            fabs(uclus-umirror) < coalesce*sigu ) {
	  return true; 
	} 
      }
      else {
	if( fabs(uhit-uwire)  <= wirP/2 &&
	    fabs(uclus-uwire) <= wirP/2 &&  
	    ( ( uhit<uwire && uclus<uwire ) || 
	      ( uhit>uwire && uclus>uwire ) ) ) { 
	  return true;
	} 
      }
    }
    else return true;   // If NOT a DRIFT DETECTOR: RETURN TRUE
  }
  // Reaching in here means "hits.empty()" or all hits w/o associated CsDetector
  return false;                               // => RETURN FALSE
}

bool CsMCUtils::isACorrelatedTrack( const CsMCTrack* mcTrack ) {
  
  // check if MC run...
  if( ! (CsEvent::Instance()->isAMonteCarloEvent()) ) {
    // MC data only...
    CsErrLog::mes( elError, "This method works only on MC data." );
    return false;
  }

  const CsMCTrack* motherTrack;
  const CsMCVertex* inVertex = mcTrack->getInVertex();
  int i;
  for (i=0;i<100;i++) {
    motherTrack = inVertex->getInTrack();
    if (motherTrack == NULL)	// found original vertex
      return (inVertex->getGnum() == 1);
    inVertex= motherTrack->getInVertex();
  }
  string str = "does not find the original vertex";	
  CsErrLog::mes( elFatal, str );
  return false;
}
  
bool CsMCUtils::isACorrelatedHit( const CsMCHit* mcHit ) {

  const CsMCTrack* mcTrack = mcHit->getMCTrack();
  return isACorrelatedTrack(mcTrack);
} 

void CsMCUtils::_buildHCMultimap() {
  
  // should never happen...
  if( ! (CsEvent::Instance()->isAMonteCarloEvent()) ) return;

  if( _currentEvent == CsEvent::Instance()->getEventNumberInRun() ) return;

  _currentEvent = CsEvent::Instance()->getEventNumberInRun();

  _hc.clear();

  std::list<CsMCHit*> myHits;
  std::list<CsCluster*> clusters = CsEvent::Instance()->getClusters();
  for( std::list<CsCluster*>::iterator Ic=clusters.begin();
       Ic!=clusters.end(); Ic++ ) {
    myHits.clear();
    std::list<CsDigit*> digits = (*Ic)->getDigitsList(); // cluster digits
    for( std::list<CsDigit*>::iterator Id=digits.begin(); 
	 Id!=digits.end(); Id++ ) {
      CsMCDigit* dig = dynamic_cast<CsMCDigit*>( *Id );
      if( dig != 0 ) {
	std::list<CsMCHit*> hits = dig->getHits(); // digit hits 
	if( !hits.empty() ) {
	  for( std::list<CsMCHit*>::iterator Ih=hits.begin(); 
	       Ih!=hits.end(); Ih++ ) {
	    if( std::find( myHits.begin(), myHits.end(), 
			   *Ih ) == myHits.end() ) {
	      myHits.push_back( *Ih );
	    }      
	  }
	}
      }
    }
    for( std::list<CsMCHit*>::iterator Ih=myHits.begin(); 
	 Ih!=myHits.end(); Ih++ ) {
      _hc.insert( std::multimap<CsMCHit*,CsCluster*>::value_type(*Ih,*Ic) );
    }
  }
}
