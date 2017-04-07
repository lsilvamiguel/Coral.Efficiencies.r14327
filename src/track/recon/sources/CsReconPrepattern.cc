
// $Id: CsReconPrepattern.cc,v 1.9 2003/04/22 12:02:39 benigno Exp $

/*!
   \file    CsReconPrepattern.h
   \brief   Recon track prepattern class.
   \author  ...
   \version $Revision: 1.9 $
   \date    $Date: 2003/04/22 12:02:39 $
*/
#include "Recon.h"
#include "Coral.h"
#include "CsReconPrepattern.h"
#include "CsErrLog.h"
#include "Recon.h"
#include "RecCall.h"
#include "RecOpt.h"
#include "CsStopwatch.h"

//Constructor
CsReconPrepattern::CsReconPrepattern() {
  new Recon;   //create  Recon package instance
  
}

CsReconPrepattern::~CsReconPrepattern() {}

bool CsReconPrepattern::doPrepattern( const std::list<CsCluster*> &clusters,
                                      const std::list<CsZone*> &zones) {
  if (RecOpt::Switch > 1 ) return( true); //do nothing	
 
  //CsErrLog::Instance()->mes( elFatal, "doPrepattern method not yet implemented for Recon" );
  
  //do nothing and return true,
  //doPrepattern is done inside getPatterns() method
  //
  return( true);
}

bool CsReconPrepattern::getPatterns( std::list<CsTrack*>& tracks ) {
  if (RecOpt::Switch > 1 ) return( true); //do nothing

    // CsStopwatch for elapsed time measurements...
    CsStopwatch stopwatches;
    int chrono = stopwatches.start(); // start a chronometer

  RecCall call;
  call.RecClus(tracks);
  //  printf(" CsReconPrepattern , jestem za RecCall \n"); 
   float time7 = stopwatches.stop(chrono); // stop a chronometer
   Recon::ref().countTime7(time7);

  return( true );
}

std::list<CsCluster*> CsReconPrepattern::getUnusedClusters(void) const {
    return Recon::ref().listUnusedClus;
// 
//   list<CsCluster*> clusters;
//  
//   
//   
//   clusters=Recon::ref().listUnusedClus;
//  
//  
//   if (RecOpt::Switch > 1 ) return(clusters); //do nothing
// 
// //  CsErrLog::Instance()->mes( elError, 
// //			     "getUnusedClusters method not yet implemented" );
//   return( clusters);
}









