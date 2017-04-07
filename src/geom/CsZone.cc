// $Id: CsZone.cc,v 1.4 2003/04/17 09:20:21 benigno Exp $

/*!
   \file    CsZone.cc
   \brief   Compass Apparatus geometrical zone
   \author  Benigno Gobbo
   \version $Revision: 1.4 $
   \date    $Date: 2003/04/17 09:20:21 $
*/

#include "CsZone.h"
#include "CsDetector.h"

CsZone::CsZone( float Zmin, float Zmax, std::string name, 
		std::list<CsDetector*> alldets ) {

  zmin_ = Zmin;
  zmax_ = Zmax;
  name_ = name;
  
  //=== Set Zone unique identifier
  static unsigned int uid=0;
  id_ = uid; 
  uid++;
  
  firstDet_ = (*(alldets.begin()));
  lastDet_ = (*(alldets.begin()));
  std::list<CsDetector*>::iterator Id;
  for( Id=alldets.begin(); Id!=alldets.end(); Id++ ) {
    float z = (*Id)->getZcm();
    if( zmin_ < z && z < zmax_ ) {
      dets_.push_back( *Id );
      if( z < firstDet_->getZcm() ) firstDet_ = (*Id);
      if( z > lastDet_->getZcm() )  lastDet_  = (*Id);
    }
  }
}

CsZone::CsZone( const CsZone& zone ) :
  name_(zone.name_),
  zmin_(zone.zmin_), 
  zmax_(zone.zmax_), 
  dets_(zone.dets_), 
  firstDet_(zone.firstDet_), 
  lastDet_(zone.lastDet_) 
{}

CsZone& CsZone::operator=( const CsZone& zone ) {
  if( this   != &zone ) {
    name_     = zone.name_;
    zmin_     = zone.zmin_; 
    zmax_     = zone.zmax_; 
    dets_     = zone.dets_; 
    firstDet_ = zone.firstDet_; 
    lastDet_  = zone.lastDet_;
  }
  return( *this );
}

bool CsZone::operator==( const CsZone& zone ) {
  if( name_     == zone.name_      &&
      zmin_     == zone.zmin_      &&
      zmax_     == zone.zmax_      &&
      dets_     == zone.dets_      &&
      firstDet_ == zone.firstDet_  &&
      lastDet_  == zone.lastDet_   ) {
    return( true );
  }
  else {
    return( false );
  }
}

bool CsZone::operator<( const CsZone& zone ) {
  if( firstDet_ < zone.firstDet_  ) {
    return( true );
  }
  else {
    return( false );
  }
}

