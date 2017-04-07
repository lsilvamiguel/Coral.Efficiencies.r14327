// $Id: CleanDetFile.cc,v 1.15 2009/08/31 01:51:30 ybedfer Exp $
/*!
   \file    CleanDetFile.cc
   \brief   DetFileManager method to clean up detector file
   \author  Hugo Pereira
   \version $Revision: 1.15 $
   \date    $Date: 2009/08/31 01:51:30 $
*/
#include "DetFileManager.h"

//_______________________________________________________________________________
void DetFileManager::CleanDetFile( void )
{  
  //=== Sort Detectors
  Sort( "FI" ); 
  Sort( "SI" ); 
  Sort( "GM" );
  Sort( "GP" );
  Sort( "DC" );
  Sort( "MM" );
  Sort( "MP" );
  Sort( "PS" );
  Sort( "PB" );
  Sort( "PA" );
  Sort( "MA" );
  Sort( "H" ); 
  
  //=== fibers
  MatchCenters( "FI**X","FI**Y");
  MatchCenters( "FI**X","FI**V");

  //=== Silicons
  MatchCenters( "SI**X","SI**Y");
  MatchCenters( "SI**U","SI**V");
 
  //=== PS
//  MatchCenters( "PS**U","PB**V");
  MatchCenters( "PS**X","PS**Y");
  
  //=== PB
//  MatchCenters( "PB01U","PB02V");
//  MatchCenters( "PB01U","PB02V");
//  MatchCenters( "PB03U","PB04V");
//  MatchCenters( "PB05U","PB06V");
//  MatchCenters( "PB**X","PB**U");
  
  //=== PA
//  MatchCenters( "PA**U","PA**V");
//  MatchCenters( "PA**U","PA**X");

  //=== MA
  MatchCenters( "MA0*X*__", "MA0*Y*__");

  //=== DW
//  MatchCenters( "DW01X__", "DW01V__");
//  MatchCenters( "DW02Y__", "DW02U__");
    
  //=== GEM
  MatchCenters( "GM**X","GM**Y");
  MatchCenters( "GM**U","GM**V");

  //=== PixelGEM
  MatchCenters( "GP**X","GP**Y");
  MatchCenters( "GP**U","GP**V");
    
  //=== DC
  MatchCenters( "DC**X1","DC**Y1" );
  MatchCenters( "DC**X2","DC**Y2" );
  MatchCenters( "DC**U1","DC**V1" );
  MatchCenters( "DC**U2","DC**V2" );
  
  
  //=== MM
  MatchCenters( "MM**X1","MM**Y1" );
  MatchCenters( "MM**U1","MM**V1" );

  //=== MP
  MatchCenters( "MP**X1","MP**Y1" );
  MatchCenters( "MP**U1","MP**V1" );

}
