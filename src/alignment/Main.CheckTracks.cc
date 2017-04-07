// $Id: Main.CheckTracks.cc,v 1.7 2008/02/20 13:56:13 rgazda Exp $
/*!
   \file    Main.CheckTracks.cc
   \brief   Instanciate a CheckTracks object, draw all plots required in option file
   \author  Hugo Pereira
   \version $Revision: 1.7 $
   \date    $Date: 2008/02/20 13:56:13 $
*/

#include "CheckTracks.h"
#include <TROOT.h>
#include <TApplication.h>
using namespace std;
//! TROOT singleton
TROOT root("checkTracks", "automatic histogram selection for alignment output tree");

/*! 
  \fn int main( int argc, char *argv[] )
  \brief arguments are [-b] (batch mode) <option file>
*/
int main( int argc, char *argv[] )
{
  if(argc<2){
    //=== Parse command line arguments - Initialise CsOpt
    printf("Usage: checkTracks [-b] <OptionFile>\n");
    printf(" -b batch mode\n");
    return 0;
  }

  //=== Create default checkTrack object
  CheckTracks check;
  
  //==== Check if job is batch job
  for( int i=1; i<argc; i++ ) {
    if( string(argv[i])=="-b" ) {
      printf("checkTracks - INFO: this is a batch job.\n");
      check.isBatch_=true;
    }
  }
  char* file = argv[argc-1];
  
  //=== Initialise batch mode application
  TApplication theApp("App", &argc, argv);
  check.DrawFromOpt( file );
  return 0;
} 
