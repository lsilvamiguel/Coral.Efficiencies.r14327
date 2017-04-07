// $Id: CsSpacePoint.cc,v 1.9 2003/04/24 07:23:26 benigno Exp $

/*!
   \file    CsSpacePoint.cc 
   \brief   Space points bluid from clusters for a certain number of detectors.
   \author  Hugo Pereira
   \version $Revision: 1.9 $
   \date    $Date: 2003/04/24 07:23:26 $
*/

#include "CsSpacePoint.h"

#include "CLHEP/Matrix/Matrix.h"
#include "CsErrLog.h"
#include "CsCluster.h"

using namespace std;

//_____________________________________________________________________________
CsSpacePoint::CsSpacePoint( const CsDetFamily &dF, 
  list<CsCluster*> c, 
	double z, int mode) :
  dF_(&dF), 
  mode_(mode),
	minimised_Fast_(false),
  x_Fast_(0), y_Fast_(0), chi2_Fast_(-1),
	
  minimised_(false),
  x_(0), y_(0), z_(z), tx_(0), ty_(0), chi2_(-1)

{
	// fill clusters 
  list<CsCluster*>::iterator Ic;
	for( Ic = c.begin(); Ic != c.end(); Ic++ ) c_.push_back( *Ic );
	
}

//_____________________________________________________________________________
CsSpacePoint& CsSpacePoint::operator=( const CsSpacePoint &sp)
{
  if( this != &sp ) { *this = sp; }
  return *this;
}

//_______________________________________________________________________________
void CsSpacePoint::dump()
{
  cout << "CsSpacePoint::dump " << endl;
	
  list< CsCluster* >::iterator Ic;
  for( Ic = c_.begin(); Ic != c_.end(); Ic++ ) 
  cout << " CsCluster: u=" << (*Ic)->getU() << " w=" << (*Ic)->getW() << endl;
	
	if( minimised_Fast_ ) {
	  cout << "  Fast"
	    << " x="<< x_Fast_  << " y=" << y_Fast_ 
	    << " chi2=" << chi2_Fast_  << endl;
	}
	 
	if( minimised_ ) {
    cout << "  Full" 
      << " x=" << x_  << " y=" << y_ 
      << " tx=" << tx_ << " ty=" << ty_ 
      << " chi2=" << chi2_ << endl;
	}
	
	return;
}	
