// $Id: CsCluster.cc,v 1.8 2006/11/29 06:12:52 ybedfer Exp $

/*!
   \file    CsCluster.cc
   \brief   Compass Single Coodinate Cluster Class.
   \author  Benigno Gobbo
   \version $Revision: 1.8 $
   \date    $Date: 2006/11/29 06:12:52 $
*/

#include "CsCluster.h"
#include "CsErrLog.h"

using namespace std;
using namespace CLHEP;

CsCluster::CsCluster() : 
  //x(3,0),
  hasAnalog_(false),  hasTime_(false),    
  time_(0), timeErr_(-1),
  hasAssociates_(false),
  hasMirror_(false),  myMirror_(NULL),
  LRProb_(1) //! By default, left/right probability is 1. The value is later reassessed for drift-like detectors. As to "isGenuine_" it's left undetermined and will only be set in the MC case and for the sole drift-like detectors.
{
  analog_.clear();
  analogErr_.clear();
  cov_ = HepMatrix(3,3,0);
  dets_.clear();
  digits_.clear();
  associates_.clear();
}

CsCluster::CsCluster( double u, double v, double w, const HepMatrix &cov  ) : 
  cov_(cov),
  hasAnalog_(false), hasTime_(false),
  time_(0), timeErr_(-1),
  hasAssociates_(false),
  hasMirror_(false),  myMirror_(NULL),
  LRProb_(1)
  //! by default, left/right probability is 1, for all non drift-like detectors.
{
  analog_.clear();
  analogErr_.clear();
  x.push_back(u);
  x.push_back(v);
  x.push_back(w);
  dets_.clear();
  digits_.clear();
  associates_.clear();
}

CsCluster::CsCluster( const vector<double> &xx, const HepMatrix &cov  ) : 
  x(xx),
  cov_(cov),
  hasAnalog_(false), hasTime_(false),
  time_(0),   timeErr_(-1),
  hasAssociates_(false),
  hasMirror_(false),  myMirror_(NULL),
  LRProb_(1)
  //! by default, left/right probability is 1, for all non drift-like detectors.
{
  analog_.clear();
  analogErr_.clear();
  dets_.clear();
  digits_.clear();
  associates_.clear();
}

CsCluster::CsCluster( const CsCluster& clus ) :
  x(clus.x), cov_(clus.cov_), 
  dets_(clus.dets_),             digits_(clus.digits_),
  hasAnalog_(clus.hasAnalog_),   analog_(clus.analog_), 
  analogErr_(clus.analogErr_),
  hasTime_(clus.hasTime_),       time_(clus.time_),
  timeErr_(clus.timeErr_),
  hasAssociates_(clus.hasAssociates_),
  associates_(clus.associates_),
  hasMirror_(clus.hasMirror_),  myMirror_(clus.myMirror_),
  LRProb_(clus.LRProb_)
{
}

CsCluster& CsCluster::operator=( const CsCluster& clus ) {
  if( this     != & clus ) {
    x           = clus.x; 
    cov_        = clus.cov_;
    dets_       = clus.dets_;
    digits_     = clus.digits_;
    hasAnalog_  = clus.hasAnalog_; 
    analog_     = clus.analog_;     
    analogErr_  = clus.analogErr_;
    hasTime_    = clus.hasTime_;
    time_       = clus.time_;
    timeErr_    = clus.timeErr_;
    hasAssociates_ = clus.hasAssociates_;
    associates_ = clus.associates_;
    hasMirror_  = clus.hasMirror_;
    myMirror_   = clus.myMirror_;
		LRProb_			= clus.LRProb_;
  }
  return( *this );
}
  
bool CsCluster::operator==( const CsCluster& clus ) const {
  if( x           == clus.x           && 
      dets_       == clus.dets_       &&
      digits_     == clus.digits_     &&
      hasAnalog_  == clus.hasAnalog_  && 
      analog_     == clus.analog_     &&
      analogErr_  == clus.analogErr_  &&
      hasTime_    == clus.hasTime_    &&
      time_       == clus.time_       &&
      timeErr_    == clus.timeErr_    &&
      hasAssociates_ == clus.hasAssociates_ &&
      associates_ == clus.associates_ &&
      hasMirror_  == clus.hasMirror_  && 
      myMirror_   == clus.myMirror_   &&
      LRProb_	  == clus.LRProb_ ) {
    return true;
  }
  else { 
    return false;
  }
}

bool CsCluster::operator<( const CsCluster& clus ) const {
  if( x.size()!=3 )
    return x<clus.x;

  if( x[2] < clus.x[2] ) {
    return( true );
  }
  else if( x[2] == clus.x[2] ) {
    if( x[0] < clus.x[0] ) {
      return( true );
    }
    else if( x[0] == clus.x[0] ) {
      if( x[1] < clus.x[1] ) {
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

void CsCluster::addDet( CsDetector& det ) {
  dets_.push_back( &det );
  dets_.sort();
}

void CsCluster::addDigit( CsDigit& digit ) {
  digits_.push_back( &digit );
}

list<CsZone*> CsCluster::getZonesList() {
  list <CsZone*> zones;
  list <CsDetector*>::iterator Id;
  for( Id=dets_.begin(); Id!=dets_.end(); Id++ ) {
    list <CsZone*> detZones = (*Id)->getMyZones();
    list <CsZone*>::iterator Iz;
    for( Iz=detZones.begin(); Iz!=detZones.end(); Iz++ ) {
      zones.push_back( *Iz );
    }
  }
  zones.sort();
  zones.unique();
  return( zones );
}

void CsCluster::setTime( const double time, const double timeError ) {
  hasTime_ = true;
  time_    = time;
  timeErr_    = timeError;
}

void CsCluster::setTimeError( const double timeError ) {
  timeErr_    = timeError;
}

bool CsCluster::getTime( double& time ) const {
  if( hasTime_ ) {
    time = time_;
    return( true );
  }
  else {
    time = 0;
    return( false );
  }
}

bool CsCluster::getTimeError( double& timeError ) const {
  if( hasTime_ ) {
    timeError = timeErr_;
    return( true );
  }
  else {
    timeError = -1;
    return( false );
  }
}

void CsCluster::setAnalog( const double analog, const double analogError ) {
  if(  analog_.size() == 0 ) {
    analog_.resize(1);
  }
  if(  analogErr_.size() == 0 ) {
    analogErr_.resize(1);
  }
  hasAnalog_    = true;
  analog_[0]    = analog;
  analogErr_[0] = analogError;
}

void CsCluster::setAnalogError( const double analogError ) {
  if(  analogErr_.size() == 0 ) {
    analogErr_.resize(1);
  }
  analogErr_[0] = analogError;
}

bool CsCluster::getAnalog( double& analog ) const {
  if( hasAnalog_ ) {
    analog = analog_[0];
    return( true );
  }
  else {
    analog = 0;
    return( false );
  }
}

bool CsCluster::getAnalogError( double& analogError ) const {
  if( hasAnalog_ ) {
    analogError = analogErr_[0];
    return( true );
  }
  else {
    analogError = -1;
    return( false );
  }
}

void CsCluster::addAnalogData( const double analogData, 
			       const double analogError ) {
  analog_.push_back( analogData );
  analogErr_.push_back( analogError );
}

bool CsCluster::updateAnalogDataVectElement( const double data, 
					     const unsigned int n ) {
  if( analog_.size() <= n ) {
    return false;
  } 
  else {
    analog_[n] = data;
    return true;
  }
}

bool CsCluster::updateAnalogDataErrorsVectElement( const double dataerr, 
						   const unsigned int n ){
  if( analogErr_.size() <= n ) {
    return false;
  } 
  else {
    analogErr_[n] = dataerr;
    return true;
  }
}

