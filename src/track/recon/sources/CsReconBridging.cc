// $Id: CsReconBridging.cc,v 1.3 2003/04/22 12:02:38 benigno Exp $

/*!
   \file    CsReconBridging.cc
   \brief   Recon track bridging class
   \author  ...
   \version $Revision: 1.3 $
   \date    $Date: 

*/

#include "CsReconBridging.h"
#include "CsErrLog.h"
#include "RecOpt.h"

CsReconBridging::CsReconBridging() {
}

CsReconBridging::~CsReconBridging() {
}

bool CsReconBridging::doBridging( std::list<CsTrack*>& tracks, 
				   std::list<CsCluster*>& clusters ) {

if (RecOpt::Switch > 1 ) return( true); //do nothing	
  CsErrLog::Instance()->mes( elFatal, "Bridging is done together  with Prepattern in Recon package " );
  return( false );
}





