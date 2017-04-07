// $Id: CsMCTrack.cc,v 1.10 2009/08/24 21:06:37 ybedfer Exp $

/*!
   \file    CsMCTrack.cc
   \brief   Compass Montecarlo Track Class.
   \author  Benigno
   \version $Revision: 1.10 $
   \date    $Date: 2009/08/24 21:06:37 $
*/

#include "CsMCTrack.h"

using namespace std;
using namespace CLHEP;

CsMCTrack::CsMCTrack() {
  Gnum_ = 0;
  p_ = HepLorentzVector( 0.0, 0.0, 0.0, 0.0 );
  particle_ = 0;
  inVertex_ = 0;
  outTracks_.clear();
  outVertices_.clear();
  recoTrackIDs_.clear();
}

CsMCTrack::CsMCTrack( int Gnum, double px, double py, double pz,
		      CsMCParticle particle, CsMCVertex& inVertex ) :
  Gnum_(Gnum), particle_(particle), inVertex_(&inVertex)
{
  double mass = particle_.getMass();
  double energy = px*px+py*py+pz*pz+mass*mass;
  energy = sqrt( energy );
  p_ = HepLorentzVector( px, py, pz, energy );
  outTracks_.clear();
  outVertices_.clear();
  recoTrackIDs_.clear();
}

CsMCTrack::CsMCTrack( int Gnum, Hep3Vector p,
		      CsMCParticle particle, CsMCVertex& inVertex ) :
  Gnum_(Gnum), particle_(particle), inVertex_(&inVertex)
{
  double mass = particle_.getMass();
  double mag  = p.mag(); 
  double energy = mag*mag+mass*mass;
  energy = sqrt( energy );
  p_ = HepLorentzVector( p, energy );
  outTracks_.clear();
  outVertices_.clear();
  recoTrackIDs_.clear();
}

CsMCTrack::CsMCTrack( const CsMCTrack& track ) {
  Gnum_        = track.Gnum_;
  p_           = track.p_;
  particle_    = track.particle_;
  MCHits_      = track.MCHits_;
  inVertex_    = track.inVertex_;
  outTracks_   = track.outTracks_;
  outVertices_ = track.outVertices_;
  recoTrackIDs_  = track.recoTrackIDs_;
}

CsMCTrack& CsMCTrack::operator=( const CsMCTrack& track ) {
  if( this != &track ) {
    Gnum_        = track.Gnum_;
    p_           = track.p_;
    particle_    = track.particle_;
    MCHits_      = track.MCHits_;
    inVertex_    = track.inVertex_;
    outTracks_   = track.outTracks_;
    outVertices_ = track.outVertices_;
    recoTrackIDs_  = track.recoTrackIDs_;
  }
  return( *this );
}

bool CsMCTrack::operator==( const CsMCTrack& track ) const {
  // Nota bene: the list of associated reco'd tracks is not checked
  // (Reason why:
  // i) Basically: At the time of introducing that part of the CsMCTrack's 
  //   code dealing with associated CsTrack's, I (Y.B.) do not know whether
  //   a change in the meaning of "operator==" might not affect any other
  //   part of CORAL.
  // ii) This sound reasonable: Cases where one wants to compare CsMCTracks
  //    are rather connected to event generation than to event reconstruction.)
  if( Gnum_        == track.Gnum_      &&
      p_           == track.p_         &&
      particle_    == track.particle_  &&
      inVertex_    == track.inVertex_  &&
      outTracks_   == track.outTracks_ && 
      outVertices_ == outVertices_     )
    return( true );
  else
    return( false );
}


bool CsMCTrack::operator<( const CsMCTrack& track ) const {
  if( Gnum_ < track.Gnum_ )
    return( true );
  else
    return( false );
}

list<CsMCHit*> CsMCTrack::getMCHits() const {
  return( MCHits_ );
}

void CsMCTrack::addMCHit( CsMCHit& MCHit ) {
if( MCHit.getMCTrack() == this )
  MCHits_.push_back( &MCHit );
}

list<CsMCTrack*> CsMCTrack::getOutTracks() const {
  return( outTracks_ );
}

list<CsMCVertex*> CsMCTrack::getOutVertices() const {
  return( outVertices_ );
}

void CsMCTrack::addOutTrack( CsMCTrack& MCTrack ) {
  outTracks_.push_back( & MCTrack );
}

void CsMCTrack::addOutVertex( CsMCVertex& MCVertex ) {
  outVertices_.push_back( & MCVertex );
}

int CsMCTrack::getGnum() const {
  return( Gnum_ );
}

double CsMCTrack::getPX() const {
  return( p_.px() );
}

double CsMCTrack::getPY() const {
  return( p_.py() );
}

double CsMCTrack::getPZ() const {
  return( p_.pz() );
}

double CsMCTrack::getE() const {
  return( p_.e() );
}

double CsMCTrack::getM() const {
  return( p_.m() );
}

HepLorentzVector CsMCTrack::getP() const {
  return( p_ );
}

CsMCParticle* CsMCTrack::getParticle() const { 
  return( const_cast<CsMCParticle*>( &particle_ ));
}

const CsMCVertex* CsMCTrack::getInVertex() const { 
  return( inVertex_ );
}
