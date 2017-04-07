// $Id: CsMCTrkHit.cc,v 1.2 2001/11/29 18:27:55 bernet Exp $

/*!
   \file    CsMCTrkHit.cc
   \brief   Compass Montecarlo Hits Class.
   \author  Benigno Gobbo
   \version $Revision: 1.2 $
   \date    $Date: 2001/11/29 18:27:55 $
*/

#include "CsMCTrkHit.h"

using namespace CLHEP;

CsMCTrkHit::CsMCTrkHit() : CsMCHit() {
  x_       = 0;
  y_       = 0;
  z_       = 0;
  uin_     = 0;
  vin_     = 0;
  win_     = 0;
  uout_    = 0;
  vout_    = 0;
  wout_    = 0;
  elos_    = 0;
  eion_    = 0;
  orig_    = 0;
  p_.set(0., 0., 0.);
}
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CsMCTrkHit::CsMCTrkHit( double x, double y, double z, 
		  double uin, double vin, double win, 
		  double uout, double vout, double wout, 
		  double elos, double eion, double dtime,
		  Hep3Vector p, CsMCTrack& MCTrack, 
		  int orig, CsDetector& det, int detid ) :
  CsMCHit( dtime, MCTrack, det, detid),
  x_(x), y_(y), z_(z), 
  uin_(uin), vin_(vin), win_(win),  
  uout_(uout), vout_(vout), wout_(wout),  
  elos_(elos), eion_(eion), orig_(orig), p_(p) {
}  
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CsMCTrkHit::CsMCTrkHit( double x, double y, double z, 
		  double uin, double vin, double win, 
		  double uout, double vout, double wout, 
		  double elos, double eion, double dtime,
		  Hep3Vector p, CsMCTrack& MCTrack, 
		  int orig, CsDetector& det) :
  CsMCHit( dtime, MCTrack, det, det.GetID()),
  x_(x), y_(y), z_(z), 
  uin_(uin), vin_(vin), win_(win),  
  uout_(uout), vout_(vout), wout_(wout),  
  elos_(elos), eion_(eion), orig_(orig), p_(p) {
}  
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CsMCTrkHit::CsMCTrkHit( const CsMCTrkHit& hit ) :
  CsMCHit(hit),
  x_(hit.x_), y_(hit.y_), z_(hit.z_), 
  uin_(hit.uin_), vin_(hit.vin_), win_(hit.win_), 
  uout_(hit.uout_), vout_(hit.vout_), wout_(hit.wout_), 
  elos_(hit.elos_), eion_(hit.eion_),
  orig_(hit.orig_), p_(hit.p_)
{}
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CsMCTrkHit& CsMCTrkHit::operator=( const CsMCTrkHit& hit ) {
  if( this != &hit ) {
    CsMCHit(*this) = CsMCHit::operator=(hit);
    x_       = hit.x_;
    y_       = hit.y_;
    z_       = hit.z_;
    uin_     = hit.uin_;
    vin_     = hit.vin_;
    win_     = hit.win_;
    uout_    = hit.uout_; 
    vout_    = hit.vout_;
    wout_    = hit.wout_;
    elos_    = hit.elos_;
    eion_    = hit.eion_;
    orig_    = hit.orig_;
    p_       = hit.p_;
  }
  return( *this );
}
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool CsMCTrkHit::operator==( const CsMCTrkHit& hit ) const {
  if( CsMCHit::operator==(hit) &&
      x_       == hit.x_       && 
      y_       == hit.y_       &&
      z_       == hit.z_       &&
      uin_     == hit.uin_     && 
      vin_     == hit.vin_     &&
      win_     == hit.win_     &&
      uout_    == hit.uout_    && 
      vout_    == hit.vout_    &&
      wout_    == hit.wout_    &&
      elos_    == hit.elos_    &&
      eion_    == hit.eion_    &&
      orig_    == hit.orig_    &&
      p_       == hit.p_       ) {
    return( true );
  }
  else {
    return( false );
  }
}
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool CsMCTrkHit::operator<( const CsMCTrkHit& hit ) const {

  if( CsMCHit::operator<(hit) ) {
    return( true );
  }
  else if( CsMCHit::operator==(hit) ) {
    if( z_ < hit.z_ ) {
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
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++






