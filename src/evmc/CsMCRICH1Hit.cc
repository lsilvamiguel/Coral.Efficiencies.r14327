// $Id: CsMCRICH1Hit.cc,v 1.1 2000/08/10 19:04:27 khaustov Exp $

/*!
   \file    CsMCRICH1Hit.cc
   \brief   Compass Monte Carlo RICH1 Hit Class.
   \author  Guennadi Khaoustov
   \version $Revision: 1.1 $
   \date    $Date: 2000/08/10 19:04:27 $
*/

#include "CsMCRICH1Hit.h"

CsMCRICH1Hit::CsMCRICH1Hit() : CsMCHit() {
  xdet_       = 0;
  ydet_       = 0;
  zdet_       = 0;
  Vdet_       = 0;
  Wdet_       = 0;
  xprod_      = 0;
  yprod_      = 0;
  zprod_      = 0;
  xref_       = 0;
  yref_       = 0;
  zref_       = 0;
  photenergy_ = 0;
  cherangle_  = 0;
  cathode_    = 0;
  orig_       = 0;
  p_.set(0., 0., 0.);
}


CsMCRICH1Hit::CsMCRICH1Hit( const CsMCRICH1Hit& hit ) : CsMCHit(hit),
  xdet_(hit.xdet_), ydet_(hit.ydet_), zdet_(hit.zdet_),
  Vdet_(hit.Vdet_), Wdet_(hit.Wdet_),
  xprod_(hit.xprod_), yprod_(hit.yprod_), zprod_(hit.zprod_), 
  xref_(hit.xref_), yref_(hit.yref_), zref_(hit.zref_),
  photenergy_(hit.photenergy_),
  cherangle_(hit.cherangle_), cathode_(hit.cathode_),
  orig_(hit.orig_), p_(hit.p_)
{} 
 
CsMCRICH1Hit& CsMCRICH1Hit::operator=( const CsMCRICH1Hit& hit ) {
  if( this != &hit ) {
    CsMCHit(*this) = CsMCHit::operator=(hit);
    xdet_       = hit.xdet_;
    ydet_       = hit.ydet_;
    zdet_       = hit.zdet_;
    Vdet_       = hit.Vdet_;
    Wdet_       = hit.Wdet_;
    xprod_      = hit.xprod_;
    yprod_      = hit.yprod_;
    zprod_      = hit.zprod_;
    xref_       = hit.xref_;
    yref_       = hit.yref_;
    zref_       = hit.zref_;
    photenergy_ = hit.photenergy_;
    cherangle_  = hit.cherangle_;
    cathode_    = hit.cathode_;
    orig_    = hit.orig_;
    p_       = hit.p_;
  }
  return( *this );
}

bool CsMCRICH1Hit::operator==( const CsMCRICH1Hit& hit ) const {
  if( CsMCHit::operator==(hit) &&
      xdet_       == hit.xdet_       &&
      ydet_       == hit.ydet_       &&
      zdet_       == hit.zdet_       &&
      Vdet_       == hit.Vdet_       &&
      Wdet_       == hit.Wdet_       &&
      xprod_      == hit.xprod_      &&
      yprod_      == hit.yprod_      &&
      zprod_      == hit.zprod_      &&
      photenergy_ == hit.photenergy_ &&
      xref_       == hit.xref_       &&
      yref_       == hit.yref_       &&
      zref_       == hit.zref_       &&
      cherangle_  == hit.cherangle_  &&
      cathode_    == hit.cathode_    &&
      orig_    == hit.orig_    &&
      p_       == hit.p_ )
      return( true );
  else
    return( false );
}

bool CsMCRICH1Hit::operator<( const CsMCRICH1Hit& hit ) const {

    if( CsMCHit::operator<(hit)) return (true);
    if( CsMCHit::operator==(hit)) {
       if(cathode_ < hit.cathode_) return (true);
    }
    return (false);
}




