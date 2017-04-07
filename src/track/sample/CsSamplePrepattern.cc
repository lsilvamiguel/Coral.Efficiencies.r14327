// $Id: CsSamplePrepattern.cc,v 1.4 2000/10/30 23:27:07 zvyagin Exp $

/*!
   \file    CsSamplePrepattern.cc
   \brief   Sample of track prepattern class.
   \author  Benigno Gobbo
   \version $Revision: 1.4 $
   \date    $Date: 2000/10/30 23:27:07 $

*/

#include "CsSamplePrepattern.h"
#include "CsErrLog.h"

CsSamplePrepattern::CsSamplePrepattern() {
}

CsSamplePrepattern::~CsSamplePrepattern() {
}

bool CsSamplePrepattern::doPrepattern( const list<CsCluster*> clusters, 
				       const list<CsZone*> zones) {
  CsErrLog::Instance()->mes( elError, "doPrepattern method not yet implemented" );
  return( false );
}

bool CsSamplePrepattern::getPatterns( list<CsTrack*>& tracks ) {
  CsErrLog::Instance()->mes( elError, "getPatterns method not yet implemented" );
  return( false );
}

const list<CsCluster*> &CsSamplePrepattern::getUnusedClusters(void) {
  list<CsCluster*> clusters;
  CsErrLog::Instance()->mes( elError, 
			     "getUnusedClusters method not yet implemented" );
  return( clusters );
}


