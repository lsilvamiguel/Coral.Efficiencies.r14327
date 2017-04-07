// $Id: CsVertex.cc,v 1.9 2010/05/24 22:04:56 ybedfer Exp $

/*!
   \file    CsVertex.cc
   \brief   Compass Vertex Class.
   \author  Benigno Gobbo
   \version $Revision: 1.9 $
   \date    $Date: 2010/05/24 22:04:56 $
*/

#include "CsVertex.h"
#include "CsErrLog.h"

using namespace std;
using namespace CLHEP;

CsVertex::CsVertex() : type_(false), x_(0), y_(0), z_(0),
          Chi2tot_(-1), Cov_(0), Tracks_(0), isBestVertex_(false)
{
  trkPar_.clear();
}


CsVertex::CsVertex(double x, double y, double z) : 
  type_(false), x_(x), y_(y), z_(z), Chi2tot_(-1),
  Cov_(0), Tracks_(0), isBestVertex_(false)
{
  trkPar_.clear();
}


CsVertex::CsVertex(double x, double y, double z, CsTrack& track) :
  type_(false), x_(x), y_(y), z_(z), Chi2tot_(-1),
  Cov_(0), isBestVertex_(false)
{
  Tracks_.push_back( &track );
}


CsVertex::CsVertex(double x, double y, double z, const list<CsTrack*> &Tracks) :
  type_(false), x_(x), y_(y), z_(z), Chi2tot_(-1),
  Cov_(0), Tracks_(Tracks), isBestVertex_(false)
{
}


CsVertex::CsVertex(const CsVertex& vertex) : 
  type_(vertex.type_), x_(vertex.x_), y_(vertex.y_), z_(vertex.z_), 
  Chi2tot_(vertex.Chi2tot_),
  Tracks_(vertex.Tracks_), trkPar_(vertex.trkPar_),
  isBestVertex_(false) // Best vertex is unique: attribute not passed on
{
  Cov_.clear();
  for(unsigned int i=0; i<vertex.Cov_.size(); i++ ) {
    Cov_.push_back( new HepMatrix( *vertex.Cov_[i] ) );
  }
}


CsVertex::~CsVertex()
{
  vector<HepMatrix*>::iterator ic;
  for( ic = Cov_.begin(); ic != Cov_.end(); ic++ ) {
    if( *ic != 0 ) delete *ic;
    else CsErrLog::Instance()->mes(elError, "It is not possible to delete cov matrix." );
  }  
  Cov_.clear();
  Tracks_.clear();
  trkPar_.clear();
}


CsVertex& CsVertex::operator=( const CsVertex& vertex ) {
  if( this != &vertex ) {
    type_      = vertex.type_;
    x_         = vertex.x_;
    y_         = vertex.y_;
    z_         = vertex.z_;
    Tracks_    = vertex.Tracks_;
    Chi2tot_   = vertex.Chi2tot_;
    isBestVertex_ = false; // Best vertex is unique: attribute not passed on
    trkPar_    = vertex.trkPar_;
    Cov_.clear();
    for(unsigned int i=0; i<vertex.Cov_.size(); i++ )
      Cov_.push_back( new HepMatrix( *vertex.Cov_[i] ) );
  }
  return( *this );
}


bool CsVertex::operator==( const CsVertex& vertex ) const {
  if( type_   == vertex.type_    &&
      x_      == vertex.x_       &&
      y_      == vertex.y_       &&
      z_      == vertex.z_       &&
      Chi2tot_== vertex.Chi2tot_ &&
      Tracks_ == vertex.Tracks_  &&
      trkPar_== vertex.trkPar_ )//&& 
    //Cov_    == vertex.Cov_ )
    return( true );
  else
    return( false );
}


HepMatrix* CsVertex::getCov( int i ) const
{
  if( (unsigned int)i >= Cov_.size() ) {
    CsErrLog::Instance()->mes(elWarning, "Index of correlation matrix is out of range." );
    return 0;
  } 
  return( Cov_[i] ); 
}


void CsVertex::setCov( const vector<HepMatrix*> &cov )
{ 
  vector<HepMatrix*>::iterator im;
  for( im=Cov_.begin(); im!=Cov_.end(); im++ )
    if( (*im) != 0 ) delete *im; 
  copy( cov.begin(), cov.end(), Cov_.begin() );
}


void CsVertex::clearCov()
{
  vector<HepMatrix*>::iterator im;
  for( im=Cov_.begin(); im!=Cov_.end(); im++ ) 
    if( (*im) != 0 ) delete *im;
  
  Cov_.clear();
}


bool CsVertex::getPar( CsTrack* trk, Cs3Vector& par ) 
{
  map<CsTrack*,Cs3Vector>::iterator im = trkPar_.find( trk );
  if( im == trkPar_.end() ) return false;
  par = im->second;
  return true; 
}
