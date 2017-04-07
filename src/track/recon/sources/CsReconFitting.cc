// $Id: CsReconFitting.cc,v 1.3 2003/04/22 12:02:38 benigno Exp $

/*!
   \file    CsReconFitting.cc
   \brief   Recon track fitting class
   \author  ...
   \version $Revision: 1.3 $
   \date    $Date: 2003/04/22 12:02:38 $

*/

#include "CsReconFitting.h"
#include "CsErrLog.h"
#include "RecOpt.h"


CsReconFitting::CsReconFitting() {
}

CsReconFitting::~CsReconFitting() {
}

bool CsReconFitting::doFitting( std::list<CsTrack*>& tracks,
				 const std::list<CsCluster*>& clusters ) {

if (RecOpt::Switch > 1 ) return( true); //do nothing	
  CsErrLog::Instance()->mes( elFatal, "Fitting not supported on Recon package" );
  return( false );
}



