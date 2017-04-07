// $Id: CsTrack.cc,v 1.41 2010/06/23 09:49:58 tnagel Exp $

/*!
   \file    CsTrack.cc
   \brief   Compass Track Class.
   \author  Benigno Gobbo
   \version $Revision: 1.41 $
   \date    $Date: 2010/06/23 09:49:58 $
*/

#include "coral_config.h"

#include "CsTrack.h"
#include "CsEvent.h"
#ifdef COMPASS_USE_OSPACE_STD
# include <ospace/std/algorithm>
#else
# include <algorithm>
#endif
#include <string.h>
#include "CsErrLog.h"
#include "CsGeom.h"

using namespace std;

CsTrack::~CsTrack()
{
}

CsTrack::CsTrack() { 
#ifndef USE_Migration
  _id        = CsEvent::Instance()->getTrackId();
#else
  _id = 0;
#endif
  _chi2      = 0;
  _ndfs      = 0;
  _time      = 0;
  _timeError = -1;
  _hasTime   = false;
  for( int i=0; i<CSTRACK_RICHDATASIZE; i++ ) _rich1Probs[i] = 0;
  _hasRich1  = false;
  _hasCEDAR  = false;
  _XX0       = 0;
  for( int i=0; i<CSTRACK_MAPSIZE; i++ ) { _expectedDets[i]=0; _firedDets[i]=0; }
  _associatedMCTrack = 0;
  _associatedTrack = -1;
  _nShowerHits = 0;
  _hasShower = 0;
}

CsTrack::CsTrack( const list<CsCluster*> &clusters ) :
  _clusters(clusters)
{
  _id        = CsEvent::Instance()->getTrackId();
  _chi2      = 0;
  _ndfs      = (int)clusters.size();
  _time      = 0;
  _timeError = -1;
  _hasTime   = false;
  for( int i=0; i<CSTRACK_RICHDATASIZE; i++ ) _rich1Probs[i] = 0;
  _hasRich1  = false;
  _hasCEDAR  = false;
  _XX0       = 0;
  for( int i=0; i<CSTRACK_MAPSIZE; i++ ) { _expectedDets[i]=0; _firedDets[i]=0; }
  _associatedMCTrack = 0;
  _associatedTrack = -1;
  _nShowerHits = 0;
  _hasShower = 0;
}

CsTrack::CsTrack( const list<CsCluster*> &clusters, const list<CsZone*> &zones ) :
  _clusters(clusters), _zones(zones)
{
  _id        = CsEvent::Instance()->getTrackId();
  _chi2      = 0;
  _ndfs      = (int)clusters.size();
  _time      = 0;
  _timeError = -1;
  _hasTime   = false;
  for( int i=0; i<CSTRACK_RICHDATASIZE; i++ ) _rich1Probs[i] = 0;
  _hasRich1  = false;
  _hasCEDAR  = false;
  _XX0       = 0;
  for( int i=0; i<CSTRACK_MAPSIZE; i++ ) { _expectedDets[i]=0; _firedDets[i]=0; }
  _associatedMCTrack = 0;
  _associatedTrack = -1;
  _nShowerHits = 0;
  _hasShower = 0;
}

CsTrack::CsTrack( const vector<CsHelix> &helices, 
		  const list<CsCluster*> &clusters,
		  const list<CsZone*> &zones ) :
  _helices(helices), _clusters(clusters), _zones(zones)
{
  _id        = CsEvent::Instance()->getTrackId();
  _chi2      = 0;
  _ndfs      = (int)clusters.size();
  _time      = 0;
  _timeError = -1;
  _hasTime   = false;
  for( int i=0; i<CSTRACK_RICHDATASIZE; i++ ) _rich1Probs[i] = 0;
  _hasRich1  = false;
  _hasCEDAR  = false;
  _XX0       = 0;
  for( int i=0; i<CSTRACK_MAPSIZE; i++ ) { _expectedDets[i]=0; _firedDets[i]=0; }
  _associatedMCTrack = 0;
  _associatedTrack = -1;
  _nShowerHits = 0;
  _hasShower = 0;
}

CsTrack::CsTrack( const vector<CsHelix> &helices, 
		  const list<CsCluster*> &clusters,
		  const list<CsZone*> &zones, double chi2 ) :
  _chi2(chi2), _helices(helices), _clusters(clusters), _zones(zones)
{
  _id        = CsEvent::Instance()->getTrackId();
  _ndfs      = (int)clusters.size();
  _time      = 0;
  _timeError = -1;
  _hasTime   = false;
  for( int i=0; i<CSTRACK_RICHDATASIZE; i++ ) _rich1Probs[i] = 0;
  _hasRich1  = false;
  _hasCEDAR  = false;
  _XX0       = 0;
  for( int i=0; i<CSTRACK_MAPSIZE; i++ ) { _expectedDets[i]=0; _firedDets[i]=0; }
  _associatedMCTrack = 0;
  _associatedTrack = -1;
  _nShowerHits = 0;
  _hasShower = 0;
}

CsTrack::CsTrack( const CsTrack& track ) :
  _id(track._id), _chi2(track._chi2), _ndfs(track._ndfs),
  _helices(track._helices), 
  _clusters(track._clusters), _zones(track._zones), 
  _vertices(track._vertices), _time(track._time), 
  _timeError(track._timeError), _hasTime(track._hasTime),
  _hasRich1(track._hasRich1), _cedarInfo(track._cedarInfo), _hasCEDAR(track._hasCEDAR),
  _XX0(track._XX0), _associatedMCTrack(track._associatedMCTrack),
  _associatedTrack(track._associatedTrack),
  _nShowerHits(track._nShowerHits), _hasShower(track._hasShower)
{
  memcpy( _rich1Probs, track._rich1Probs, sizeof( _rich1Probs ) ); 
  memcpy( _expectedDets, track._expectedDets, sizeof( _expectedDets ) ); 
  memcpy( _firedDets, track._firedDets, sizeof( _firedDets ) ); 
}

CsTrack& CsTrack::operator=( const CsTrack& track ) {
  if( this    != &track ) {
    _id        = track._id;
    _chi2      = track._chi2;
    _ndfs      = track._ndfs;
    _helices   = track._helices;
    _clusters  = track._clusters;
    _zones     = track._zones;
    _vertices  = track._vertices;
    _time      = track._time;
    _timeError = track._timeError;
    _hasTime   = track._hasTime;
    memcpy( _rich1Probs, track._rich1Probs, sizeof( _rich1Probs ) ); 
    _hasRich1  = track._hasRich1;
    _hasCEDAR  = track._hasCEDAR;
    _cedarInfo = track._cedarInfo;
    _XX0       = track._XX0;
    memcpy( _expectedDets, track._expectedDets, sizeof( _expectedDets ) ); 
    memcpy( _firedDets, track._firedDets, sizeof( _firedDets ) ); 
    _associatedMCTrack = track._associatedMCTrack;
    _associatedTrack = track._associatedTrack;
    _nShowerHits = track._nShowerHits;
    _hasShower   = track._hasShower;
  }
  return( *this );
}

bool CsTrack::operator==( const CsTrack& track ) {
  if (_id         == track._id         &&
      _chi2       == track._chi2       &&
      _ndfs       == track._ndfs       &&
      _helices    == track._helices    &&
      _clusters   == track._clusters   &&
      _zones      == track._zones      &&
      _vertices   == track._vertices   &&
      _time       == track._time       &&
      _timeError  == track._timeError  &&
      _hasTime    == track._hasTime    &&
      memcmp(_rich1Probs,  track._rich1Probs,  sizeof(_rich1Probs)) &&
      _hasRich1   == track._hasRich1   &&
      _hasCEDAR   == track._hasCEDAR   &&
      _cedarInfo  == track._cedarInfo  &&
      _XX0        == track._XX0        &&
      memcmp(_expectedDets,track._expectedDets,sizeof(_expectedDets)) &&
      memcmp(_firedDets,   track._firedDets,   sizeof(_firedDets)) &&
      _nShowerHits == track._nShowerHits &&
      _hasShower == track._hasShower) {
    return true ;
  }
  else
    return false;

}

void CsTrack::addHelix( const CsHelix &helix ) {
  if( find( _helices.begin(), _helices.end(), helix ) == _helices.end() ) {
    _helices.push_back( helix );
  }
}

void CsTrack::addCluster( CsCluster& cluster ) {
  if( find(_clusters.begin(), _clusters.end(), &cluster)==_clusters.end() ) {
    _clusters.push_back( &cluster );
  }
}

void CsTrack::addZone( CsZone& zone ) {
  if( find( _zones.begin(), _zones.end(), &zone ) == _zones.end() ) {
    _zones.push_back( &zone );
  }
}

CsVertex* CsTrack::getFirstVertex() const { 
  if( !_vertices.empty() ) {
    return( _vertices.front() ); 
  }
  else {
    return( NULL );
  }
}

double CsTrack::getMeanTime() const {
  if( _hasTime ) {
    return( _time );
  }
  else {
    CsErrLog::mes( elError, "Mean time not set for this track." );
    // should end processing...
    throw( "Mean time not set for this track." );
    return( 0 );
  }
}

double CsTrack::getMeanTimeError() const {
  if( _hasTime ) {
    return( _timeError );
  }
  else {
    CsErrLog::mes( elError, "Mean time error not set for this track." );
    // should end processing...
    throw( "Mean time not set for this track." );
    return( -1 );
  }
}

void CsTrack::addCEDARInfo() {

  double zTarget = CsGeom::Instance()->getTargetCenter();  // in mm
  CsHelix* foundHelix = NULL;
  //get helix downstream of the target
  for(std::vector<CsHelix>::iterator It = _helices.begin(); It != _helices.end(); ++It) {
    if(It->getZ() < zTarget)
      if(!foundHelix || It->getZ() > foundHelix->getZ())
        foundHelix = &(*It);
  }

  if(foundHelix) {
    bool _hasCEDAR = true;
    double dxCEDAR, dyCEDAR;
    //get beam divergence between the two cedars, according to Laus transport matrix
    //as the parameters are read from a config file a cedar instance has to be accessed
    //the parameters are identical for both cedars (only one call necessary)
    if(CsGeom::Instance()->getCEDARs().size() != 0)
      CsGeom::Instance()->getCEDARs()[0]->getDXDY(foundHelix, dxCEDAR, dyCEDAR);

    for(unsigned int ce = 0; ce < CsGeom::Instance()->getCEDARs().size(); ++ce) {
      CsGeom::Instance()->getCEDARs()[ce]->getLikelihoods(_cedarInfo, dxCEDAR, dyCEDAR);
    }
  }
}

void CsTrack::setRich1Probs( const double* probs ) {
  memcpy( _rich1Probs, probs, sizeof(_rich1Probs) );
  _hasRich1  = true;
}


const double* CsTrack::getRich1Probs() const {
  if( _hasRich1 ) {
    return( _rich1Probs );
  }
  else {
    CsErrLog::mes( elError, "Rich probabilities not set for this track." );
    // should end processing...
    throw( "Rich probabilities not set for this track." );
    return( NULL );
  } 
}

double CsTrack::BkgLikelihood() {
  if( _hasRich1 ) {
    return( _rich1Probs[0] );
  }
  else {
    CsErrLog::mes( elError, "Rich probabilities not set for this track." );
    // should end processing...
    throw( "Rich probabilities not set for this track." );
    return( 0 );
  } 
}

double CsTrack::PionLikelihood() {
  if( _hasRich1 ) {
    return( _rich1Probs[1] );
  }
  else {
    CsErrLog::mes( elError, "Rich probabilities not set for this track." );
    // should end processing...
    throw( "Rich probabilities not set for this track." );
    return( 0 );
  } 
}

double CsTrack::KaonLikelihood() {
  if( _hasRich1 ) {
    return( _rich1Probs[2] );
  }
  else {
    CsErrLog::mes( elError, "Rich probabilities not set for this track." );
    // should end processing...
    throw( "Rich probabilities not set for this track." );
    return( 0 );
  } 
}

double CsTrack::ProtonLikelihood() {
  if( _hasRich1 ) {
    return( _rich1Probs[3] );
  }
  else {
    CsErrLog::mes( elError, "Rich probabilities not set for this track." );
    // should end processing...
    throw( "Rich probabilities not set for this track." );
    return( 0 );
  } 
}

double CsTrack::PionDLikeDIndex() {
  if( _hasRich1 ) {
    return( _rich1Probs[4] );
  }
  else {
    CsErrLog::mes( elError, "Rich probabilities not set for this track." );
    // should end processing...
    throw( "Rich probabilities not set for this track." );
    return( 0 );
  } 
}

double CsTrack::KaonDLikeDIndex() {
  if( _hasRich1 ) {
    return( _rich1Probs[5] );
  }
  else {
    CsErrLog::mes( elError, "Rich probabilities not set for this track." );
    // should end processing...
    throw( "Rich probabilities not set for this track." );
    return( 0 );
  } 
}

double CsTrack::ProtonDLikeDIndex() {
  if( _hasRich1 ) {
    return( _rich1Probs[6] );
  }
  else {
    CsErrLog::mes( elError, "Rich probabilities not set for this track." );
    // should end processing...
    throw( "Rich probabilities not set for this track." );
    return( 0 );
  } 
}

double CsTrack::Rich1ThetaMaxLike() {
  if( _hasRich1 ) {
    return( _rich1Probs[7] );
  }
  else {
    CsErrLog::mes( elError, "Rich probabilities not set for this track." );
    // should end processing...
    throw( "Rich probabilities not set for this track." );
    return( 0 );
  } 
}

double CsTrack::Rich1Theta() {
  if( _hasRich1 ) {
    return( _rich1Probs[8] );
  }
  else {
    CsErrLog::mes( elError, "Rich response not available for this track." );
    // should end processing...
    throw( "Rich response not available for this track." );
    return( 0 );
  } 
}

int CsTrack::Rich1NbPhot() {
  if( _hasRich1 ) {
    return( int( _rich1Probs[9] ) );
  }
  else {
    CsErrLog::mes( elError, "Rich response not available for this track." );
    // should end processing...
    throw( "Rich response not available for this track." );
    return( 0 );
  } 
}

double CsTrack::Rich1ThetaFit() {
  if( _hasRich1 ) {
    return( _rich1Probs[10] );
  }
  else {
    CsErrLog::mes( elError, "Rich probabilities not set for this track." );
    // should end processing...
    throw( "Rich probabilities not set for this track." );
    return( 0 );
  } 
}

double CsTrack::Rich1RingChi2() {
  if( _hasRich1 ) {
    return( _rich1Probs[11] );
  }
  else {
    CsErrLog::mes( elError, "Rich probabilities not set for this track." );
    // should end processing...
    throw( "Rich probabilities not set for this track." );
    return( 0 );
  } 
}

double CsTrack::PionChi2() {
  if( _hasRich1 ) {
    return( _rich1Probs[12] );
  }
  else {
    CsErrLog::mes( elError, "Rich probabilities not set for this track." );
    // should end processing...
    throw( "Rich probabilities not set for this track." );
    return( 0 );
  } 
}

double CsTrack::KaonChi2() {
  if( _hasRich1 ) {
    return( _rich1Probs[13] );
  }
  else {
    CsErrLog::mes( elError, "Rich probabilities not set for this track." );
    // should end processing...
    throw( "Rich probabilities not set for this track." );
    return( 0 );
  } 
}

double CsTrack::ProtonChi2() {
  if( _hasRich1 ) {
    return( _rich1Probs[14] );
  }
  else {
    CsErrLog::mes( elError, "Rich probabilities not set for this track." );
    // should end processing...
    throw( "Rich probabilities not set for this track." );
    return( 0 );
  } 
}

double CsTrack::ElectronLikelihood( void ) {
  if( _hasRich1 && CSTRACK_RICHDATASIZE > 15 ) {
    return( _rich1Probs[15] );
  }
  else {
    CsErrLog::mes( elError, "Rich probabilities not set for this track." );
    // should end processing...
    throw( "Rich probabilities not set for this track." );
    return( 0 );
  } 
}

double CsTrack::MuonLikelihood( void ){
  if( _hasRich1 && CSTRACK_RICHDATASIZE > 16 ) {
    return( _rich1Probs[16] );
  }
  else {
    CsErrLog::mes( elError, "Rich probabilities not set for this track." );
    // should end processing...
    throw( "Rich probabilities not set for this track." );
    return( 0 );
  } 
}

double CsTrack::ElectronDLikeDIndex( void ){
  if( _hasRich1 && CSTRACK_RICHDATASIZE > 17 ) {
    return( _rich1Probs[17] );
  }
  else {
    CsErrLog::mes( elError, "Rich probabilities not set for this track." );
    // should end processing...
    throw( "Rich probabilities not set for this track." );
    return( 0 );
  } 
}

double CsTrack::MuonDLikeDIndex( void ){
  if( _hasRich1 && CSTRACK_RICHDATASIZE > 18 ) {
    return( _rich1Probs[18] );
  }
  else {
    CsErrLog::mes( elError, "Rich probabilities not set for this track." );
    // should end processing...
    throw( "Rich probabilities not set for this track." );
    return( 0 );
  } 
}

double CsTrack::ElectronChi2( void ){
  if( _hasRich1 && CSTRACK_RICHDATASIZE > 19 ) {
    return( _rich1Probs[19] );
  }
  else {
    CsErrLog::mes( elError, "Rich probabilities not set for this track." );
    // should end processing...
    throw( "Rich probabilities not set for this track." );
    return( 0 );
  } 
}

double CsTrack::MuonChi2( void ){
  if( _hasRich1 && CSTRACK_RICHDATASIZE > 20 ) {
    return( _rich1Probs[20] );
  }
  else {
    CsErrLog::mes( elError, "Rich probabilities not set for this track." );
    // should end processing...
    throw( "Rich probabilities not set for this track." );
    return( 0 );
  } 
}

bool CsTrack::setExpectedDetsBitmap( const unsigned int* bitmap ) {
  memcpy( _expectedDets, bitmap, sizeof(_expectedDets) );
  return true;
}

bool CsTrack::setFiredDetsBitmap( const unsigned int* bitmap ) {
  memcpy( _firedDets, bitmap, sizeof(_firedDets) );
  return true;
}

const unsigned int* CsTrack::getExpectedDetsBitmap() const {
  return( _expectedDets );
}

const unsigned int* CsTrack::getFiredDetsBitmap() const {
  return( _firedDets );
}

int CsTrack::getNumberOfAssociatedClusters() const {
  int hits = _clusters.size();
  if( hits == 0 ) {
    for( int i=0; i<CSTRACK_MAPSIZE; i++ ){
      for( int j=0; j<32; j++ ) {
	hits += ((_firedDets[i]>>j)&1);
      }
    }
  }
  
  return hits;
}

const CsHelix* CsTrack::getHelixUpstreamOf( double z ) const {
  if ( _helices.size() == 0 )
    return NULL;

  // we cannot make any assumption on the ordering of helices
  vector<CsHelix>::const_iterator near = _helices.begin();
  for( vector<CsHelix>::const_iterator hlx = _helices.begin(); hlx != _helices.end(); hlx++ ) {
    // if initialised downstream of z, take any helix that is more upstream
    if ( near->getZ() > z && hlx->getZ() < near->getZ() )
      near = hlx;
    // else take only helices that are more downstream but still upstream of z
    else if ( hlx->getZ() > near->getZ() && hlx->getZ() <= z )
      near = hlx;
  }
  return &*near;
}
