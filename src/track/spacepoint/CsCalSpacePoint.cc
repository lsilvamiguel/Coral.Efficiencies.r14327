// $Id: CsCalSpacePoint.cc,v 1.3 2003/04/24 07:23:24 benigno Exp $

/*!
   \file    CsCalSpacePoint.cc 
   \brief   Calibrarion spacepoint, derived from CsSpacePoint Class.
   \author  Hugo Pereira
   \version $Revision: 1.3 $
   \date    $Date: 2003/04/24 07:23:24 $
*/

#include "CsCalSpacePoint.h"
#include "CsSpacePoint.h"

//_____________________________________________________________________________
CsCalSpacePoint::CsCalSpacePoint( const CsDetFamily &df, std::list<CsCluster*> c, double z, int mode):
  CsSpacePoint( df, c, z, mode ),
  detOff_( 0 ),
  found_( false ),
  clFound_( NULL )
{}

//_____________________________________________________________________________
CsCalSpacePoint& CsCalSpacePoint::operator=( const CsCalSpacePoint &sp)
{
  if( this != &sp ) { *this = sp; }
  return *this;
}

