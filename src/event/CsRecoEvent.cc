// $Id: CsRecoEvent.cc,v 1.40 2010/08/26 17:25:19 tnagel Exp $

/*!
   \file    CsRecoEvent.h
   \brief   Compass Class CsRecoEvent
   \author  Benigno Gobbo
   \version $Revision: 1.40 $
   \date    $Date: 2010/08/26 17:25:19 $

   This will contain all reconstructed transient objects. 
*/

#include "CsRecoEvent.h"
#ifdef COMPASS_USE_OSPACE_STD
# include <ospace/std/algorithm>
# include <ospace/std/functional>
#else
# include <algorithm>
# include <functional>
#endif
#include "CsMCDigit.h"
#include "CsGeom.h"
#include "CsBeamRecons.h"
#include "Reco/CalorimeterParticle.h"

// Static Members definition
int CsRecoEvent::_nextTrkId = 0;
// end of static members definition

CsRecoEvent::CsRecoEvent() {
  clear();
}  

CsRecoEvent::~CsRecoEvent() {
  clear();
}  

void CsRecoEvent::clear() {
  clearTracks();
  clearDigits();
  clearClusters();
  clearCalObjs();
  clearBeamTracksList();
  clearVertices();
  clearParticles();
  resetTriggerTime();
  for( int i=0; i<8; i++ ) _scalers[i]=0;
  _hodoData.clear();
  for( int i=0; i<10; i++ ) _eventHeader[i] = 0;
  _timeinspill = 0;
  _incidentFlux[0] = 0;
  _incidentFlux[1] = 0;
}

void CsRecoEvent::reset() {
  clearTracks();
  clearBeamTracksList();
  clearCalObjs();
  clearParticles();
  clearVertices();
}

void CsRecoEvent::setParticles( const std::vector<CsParticle*> &particles ) {
  _particles = particles;
}

void CsRecoEvent::setTracks( const std::list<CsTrack*> &tracks ) {
  if( !_tracks.empty() ) {
    for( std::list<CsTrack*>:: iterator I=_tracks.begin(); I!=_tracks.end(); I++ ) {
      delete *I;
    }
    _tracks.clear();
  }

  for( std::list<CsTrack*>::const_iterator I=tracks.begin(); I!=tracks.end(); I++ )
    _tracks.push_back( new CsTrack(**I) );
}

void CsRecoEvent::clearTracks() {
  if( !_tracks.empty() ) {
    std::list<CsTrack*>:: iterator i;
    for( i=_tracks.begin(); i!=_tracks.end(); i++ ) {
      delete *i;
    }
    _tracks.clear();
  }
  _nextTrkId = 0;
}

void CsRecoEvent::setVertices( const std::list<CsVertex*> &vertices ) {
  _vertices = vertices;
}

void CsRecoEvent::clearVertices() {
  if( !_vertices.empty() ) {
    std::list<CsVertex*>:: iterator i;
    for( i=_vertices.begin(); i!=_vertices.end(); i++ ) {
      delete *i;
    }
    _vertices.clear();
  }
}

void CsRecoEvent::clearDigits() {
  if( !_digits.empty() ) {
    std::list<CsDigit*>:: iterator i;
    for( i=_digits.begin(); i!=_digits.end(); i++ ) {
      CsMCDigit* mdig = dynamic_cast<CsMCDigit*>(*i);
      if( mdig == NULL ) {
	delete *i;
      }
      else {
	delete mdig;
      }
    }
    _digits.clear();
  }
}

void CsRecoEvent::clearClusters() {
  if( !_clusters.empty() ) {
    std::list<CsCluster*>:: iterator i;
    for( i=_clusters.begin(); i!=_clusters.end(); i++ ) {
      delete *i;
    }
    _clusters.clear();
		
    // empty all detectors list of clusters
    std::list<CsDetector*> d = CsGeom::Instance()->getDetectors();
    std::list<CsDetector*>::iterator Id; 
    for( Id=d.begin(); Id != d.end(); Id++ ) (*Id)->clearClusterList(); 
  }
}

void CsRecoEvent::clearBeamTracks() {
  if( !_beamTracks.empty() ) {
    std::list<CsBeam*>:: iterator i;
    for( i=_beamTracks.begin(); i!=_beamTracks.end(); i++ ) {
      delete *i;
    }
  }
}

void CsRecoEvent::subtractCluster(void) { 
  delete *(_clusters.rbegin());
  _clusters.pop_back(); 
}

struct CsRecoEvent::_sortDigits : 
  public std::binary_function<CsDigit*, CsDigit*, bool> {
  bool operator() (CsDigit* d1, CsDigit* d2) { 

    double z1 = 0;
    double z2 = 0;
    CsDetector* det1 = dynamic_cast<CsDetector*>(d1->getDet());
    CsDetector* det2 = dynamic_cast<CsDetector*>(d2->getDet());
    if( det1 != 0 )  z1 = det1->getZcm();
    if( det2 != 0 )  z2 = det2->getZcm();

    int addr1 = d1->getAddress();
    int addr2 = d2->getAddress();

    if( z1 < z2 ) {
      return( true );
    }
    else if( z1 == z2 ) {
      if( addr1 < addr2 ) {
	return( true );
      }
      else {
	return( false );
      }
    }
    else {
      return( false );
    }
  }

};

struct CsRecoEvent::_sortClusters : 
  public std::binary_function<CsCluster*, CsCluster*, bool> {
  bool operator() (CsCluster* c1, CsCluster* c2) { 
    double u1 = c1->getU();
    double u2 = c2->getU();
    double v1 = c1->getV();
    double v2 = c2->getV();
    double w1 = c1->getW();
    double w2 = c2->getW();
  
    if( w1 < w2 ) {
      return( true );
    }
    else if( w1 == w2 ) {
      if( u1 < u2 ) {
	return( true );
      }
    else if( u1 == u2 ) {
      if( v1 < v2 ) {
	return( true );
      }
      else {
	return( false );
      }
    }
      else {
	return( false );
      }
    }
    else {
      return( false );
    }
  }
};

void CsRecoEvent::sortDigits() { 
  _digits.sort( _sortDigits() ); 
}

void CsRecoEvent::sortClusters() { 
  _clusters.sort( _sortClusters() ); 
}

void CsRecoEvent::clearParticles(void) {
  if( !_particles.empty() ) {
    for( unsigned int i=0; i<_particles.size(); i++ ) {
      delete _particles[i];
    }
    _particles.clear();
  }
}

void CsRecoEvent::setCalObjsVector( 
        const std::vector<Reco::CalorimeterParticle> calobjs ) {
  _calObjs.resize( calobjs.size() );
  for( unsigned int i=0; i<calobjs.size(); i++ ) {
    Reco::CalorimeterParticle* calobj = 
      new Reco::CalorimeterParticle( calobjs[i] );
    _calObjs[i] = calobj;
  }
}    

void CsRecoEvent::clearCalObjs(void) {
  if( !_calObjs.empty() ) {
    for( unsigned int i=0; i<_calObjs.size(); i++ ) {
      delete _calObjs[i];
    }
    _calObjs.clear();
  }
}


void CsRecoEvent::setScalers( const int* scalers ) {
  for( int i=0; i<8; i++ ) {
    _scalers[i] = scalers[i];
  }
}

void CsRecoEvent::setEventHeader( const unsigned int datum, const unsigned int position ) {
  if( position<10 ) {
    _eventHeader[position] = datum;
  }
}

void CsRecoEvent::setEventHeader( const unsigned int* data ) {
  for( unsigned i=0; i<10; i++ ) {
    _eventHeader[i] = data[i];
  }
}
