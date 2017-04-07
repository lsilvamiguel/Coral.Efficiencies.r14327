// $Id: CsSampleBridging.cc,v 1.3 2000/02/03 17:06:06 benigno Exp $

/*!
   \file    CsSampleBridging.cc
   \brief   Sample of track bridging class
   \author  Benigno Gobbo
   \version $Revision: 1.3 $
   \date    $Date: 2000/02/03 17:06:06 $

*/

#include "CsSampleBridging.h"
#include "CsErrLog.h"

CsSampleBridging::CsSampleBridging() {
}

CsSampleBridging::~CsSampleBridging() {
}

bool CsSampleBridging::doBridging( list<CsTrack*>& tracks, 
				   list<CsCluster*>& clusters ) {
  CsErrLog::Instance()->mes( elError, "doBridging method not yet implemented" );
  return( false );
}

