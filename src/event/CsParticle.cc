// $Id: CsParticle.cc,v 1.10 2010/02/10 17:25:32 tnagel Exp $

/*!
   \file    CsParticle.cc
   \brief   Compass Particle Class.
   \author  Benigno Gobbo
   \version $Revision: 1.10 $
   \date    $Date: 2010/02/10 17:25:32 $
*/

#include "CsParticle.h"

CsParticle::CsParticle() {
  _type   = ORDINARY;
  _charge = -999;
  _track  = NULL;
  _calObjs.clear();
}

CsParticle::CsParticle( CsTrack* track ) {
  int nplus=0;
  int nminus=0;
//    vector<CsHelix> hlx = track->getHelices();
//    for( unsigned int i=0; i<hlx.size(); i++ ) {
//      if( hlx[i].getCop() > 0 ) {
//        nplus++;
//      }
//      else if( hlx[i].getCop() < 0 ) {
//        nminus++;
//      }
//    }
//    if( nplus > nminus ) {
//      _charge = 1;
//    }
//    else if( nplus < nminus ) {
//      _charge = -1;
//    }
//    else {
//      _charge = -999;
//    }
  if( !(track->getHelices()).empty() ) {
    double cop = (track->getHelices()).front().getCop();
    if( cop > 0 ) {
      _charge = 1;
    }
    else if( cop < 0 ) {
      _charge = -1;
    }
    else {
      _charge = -999;
    }
  }
  else {
    _charge = -999;
  }
  _track  = track;
  _type   = ORDINARY;
  _calObjs.clear();
}

CsParticle::CsParticle( Reco::CalorimeterParticle* calobj ) {
  _type   = ORDINARY;
  _charge = 0;
  _track  = NULL;
  _calObjs.push_back( calobj );
}


CsParticle::CsParticle( const CsParticle& part ) :
  _charge(part._charge), _track(part._track), _type(part._type), _calObjs(part._calObjs)
{
}

CsParticle& CsParticle::operator=( const CsParticle& part ) {
  if( this  != & part ) {
    _type    = part._type;
    _charge  = part._charge;
    _track   = part._track;
    _calObjs = part._calObjs;
  }
  return( *this );
}
 
