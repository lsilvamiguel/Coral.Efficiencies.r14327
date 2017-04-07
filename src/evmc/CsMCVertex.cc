// $Id: CsMCVertex.cc,v 1.6 2003/04/09 15:48:30 benigno Exp $

/*!
   \file    CsMCVertex.cc
   \brief   Compass Monte Carlo Vertex Class.
   \author  Benigno Gobbo
   \version $Revision: 1.6 $
   \date    $Date: 2003/04/09 15:48:30 $
*/

#include "CsMCVertex.h"

using namespace std;

CsMCVertex::CsMCVertex() {
  x_ = 0;
  y_ = 0;
  z_ = 0;
  inTrack_ = 0;
  outTracks_.clear();
}

CsMCVertex::CsMCVertex( int Gnum, double x, double y, double z, double t ) :
  Gnum_(Gnum), x_(x), y_(y), z_(z), t_(t)
{
  inTrack_ = 0;
  outTracks_.clear();
}


CsMCVertex::CsMCVertex( int Gnum, double x, double y, double z, double t, 
			CsMCTrack& track ) :
   Gnum_(Gnum), x_(x), y_(y), z_(z), t_(t), inTrack_(&track)
{
  outTracks_.clear();
}

CsMCVertex::CsMCVertex( const CsMCVertex& vertex ) {
  Gnum_      = vertex.Gnum_;
  x_         = vertex.x_;
  y_         = vertex.y_;
  z_         = vertex.z_;
  t_         = vertex.t_;
  inTrack_   = vertex.inTrack_;
  outTracks_ = vertex.outTracks_;
}

CsMCVertex& CsMCVertex::operator=( const CsMCVertex& vertex ) {
  if( this != &vertex ) {
    Gnum_      = vertex.Gnum_;
    x_         = vertex.x_;
    y_         = vertex.y_;
    z_         = vertex.z_;
    t_         = vertex.t_;
    inTrack_   = vertex.inTrack_;
    outTracks_ = vertex.outTracks_;
  }
  return( *this );
}

bool CsMCVertex::operator==( const CsMCVertex& vertex ) const {
  if( x_ == vertex.x_ &&
      y_ == vertex.y_ &&
      z_ == vertex.z_ &&
      t_ == vertex.t_ &&
      Gnum_ == vertex.Gnum_ &&
      inTrack_ == vertex.inTrack_ &&
      outTracks_ == vertex.outTracks_ ) 
    return( true );
  else
    return( false );
}

bool CsMCVertex::operator<( const CsMCVertex& vertex ) const {
  if( Gnum_ < vertex.Gnum_ )
    return( true );
  else
    return( false );
}


void CsMCVertex::setInTrack( CsMCTrack& inTrack ) {
  inTrack_ = &inTrack;
}

void CsMCVertex::addOutTrack( CsMCTrack& outTrack ) {
  outTracks_.push_back( &outTrack );
}


int CsMCVertex::getGnum() const {
  return( Gnum_ );
}

double CsMCVertex::getX() const {
  return( x_ );
}

double CsMCVertex::getY() const {
  return( y_ );
}

double CsMCVertex::getZ() const {
  return( z_ );
}

double CsMCVertex::getT() const {
  return( t_ );
}

const CsMCTrack* CsMCVertex::getInTrack() const {
  return( inTrack_ );
}

list<CsMCTrack*> CsMCVertex::getOutTracks() {
  return( outTracks_ );
}
