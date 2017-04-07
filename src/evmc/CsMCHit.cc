// $Id: CsMCHit.cc,v 1.6 2001/11/29 18:27:54 bernet Exp $

/*!
   \file    CsMCHit.cc
   \brief   Compass Montecarlo Hits Class.
   \author  Guennadi Khaoustov
   \version $Revision: 1.6 $
   \date    $Date: 2001/11/29 18:27:54 $
*/

#include "CsMCHit.h"
#include "CsMCTrack.h"

//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CsMCHit::CsMCHit() {
  dtime_   = 0;
  MCTrack_ = 0;
  det_     = 0;
}
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CsMCHit::CsMCHit(double dtime, CsMCTrack& MCTrack, 
	         CsDet& det, int detid) 
  : dtime_(dtime), MCTrack_(&MCTrack), det_(&det), id_(detid) {}  
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CsMCHit::CsMCHit(double dtime, CsMCTrack& MCTrack, CsDet& det ) 
  : dtime_(dtime), MCTrack_(&MCTrack), det_(&det), id_(0) {}  
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CsMCHit::CsMCHit( const CsMCHit& hit ) :
  dtime_(hit.dtime_), MCTrack_(hit.MCTrack_),  det_(hit.det_), id_(hit.id_) 
{}
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CsMCHit& CsMCHit::operator=( const CsMCHit& hit ) {
  if( this != &hit ) {
    dtime_   = hit.dtime_;
    MCTrack_ = hit.MCTrack_;
    det_     = hit.det_;
    id_      = hit.id_;
  }
  return( *this );
}
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool CsMCHit::operator==( const CsMCHit& hit ) const {
  if( 
      dtime_   == hit.dtime_   &&
      MCTrack_ == hit.MCTrack_ &&
      det_     == hit.det_  &&
      id_      == hit.id_) 
    return( true );
  else
    return( false );
}
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool CsMCHit::operator<( const CsMCHit& hit ) const {

  if( MCTrack_->getGnum() < hit.MCTrack_->getGnum() ) {
    return( true );
  }
  else if( MCTrack_->getGnum() == hit.MCTrack_->getGnum() ) {
    if( dtime_ < hit.dtime_ ) {
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
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




