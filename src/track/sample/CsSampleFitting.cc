// $Id: CsSampleFitting.cc,v 1.3 2000/02/03 17:06:06 benigno Exp $

/*!
   \file    CsSampleFitting.cc
   \brief   Sample of track fitting class
   \author  Benigno Gobbo
   \version $Revision: 1.3 $
   \date    $Date: 2000/02/03 17:06:06 $

*/

#include "CsSampleFitting.h"
#include "CsErrLog.h"

CsSampleFitting::CsSampleFitting() {
}

CsSampleFitting::~CsSampleFitting() {
}

bool CsSampleFitting::doFitting( list<CsTrack*>& tracks,
				 const list<CsCluster*>& clusters ) {
  CsErrLog::Instance()->mes( elError, "doFitting method not yet implemented" );
  return( false );
}



