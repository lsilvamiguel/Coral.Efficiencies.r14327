// $Id: Main.display3D.cc,v 1.6 2010/09/07 18:26:59 tnagel Exp $
/*!
   \file    Main.display3D.cc
   \brief   Runs traffic/trafdic on modified option file to generate alignment tree
   \author  Hugo Pereira
   \version $Revision: 1.6 $
   \date    $Date: 2010/09/07 18:26:59 $
*/

#include "DaqDataDecoding/Exception.h"
#include "Coral.h"
#include "CsInit.h"
#include "CsOpt.h"
#include "DrawTracks3D.h"

#include <TROOT.h>
#include <TStyle.h>
#include <TApplication.h>

using namespace std;

/*! \fn int main( int argc, char *argv[] )
  \brief main file. 
  \param argc number of arguments
  \param argv arguments must at least contain option file
*/

int main( int argc, char *argv[] ) {


  try
  {                  
    
    // check options
    CsOpt *opt = CsOpt::Instance( argc, argv );
    
    // check that traffic event display is off
    vector< int > trafOpt;
    if( opt->getOpt( "TraF", "Graph", trafOpt ) &&
      trafOpt.size() && trafOpt[0] ) {
      cout << "Main - FATAL: cannot run together with traffic event display.\n"
        << "Switch of TraF Graph [0] in option file.\n";
      return 0;
    }
    
    // Check detector file
    string detFile;
    if( ! opt->getOpt("detector", "table", detFile) ) {
      cout << "Main - FATAL: no detector table found.\n";
      return 0;
    }
        
    // initialise coral    
    Coral::init( argc, argv ); 
       
    // initialise application 
    gROOT->Macro("rootlogon.C");
    char *root_argv[]={argv[0],NULL};
    int root_argc=1;
    TApplication application("App",&root_argc, root_argv);
    
    // check batch
    if( gROOT->IsBatch() ) {
      cout << "Main - FATAL: cannot run in batch mode. Abort.\n";
      return 0;
    }
        
    // Get first event
    if( !CsEvent::Instance()->getNextEvent() ) { 
      cout << "Main - FATAL: troubles getting first event.\n";
      return 0;
    }
    
    // Initialise DrawTracks3D
    DrawTracks3D dt( "", detFile.c_str(), true );
    application.Run();
  } 
     
  catch(std::exception &e ){ cerr << "Exception:\n" << e.what() << endl; }
  catch( ... ){ cerr << "Unknown exception!" << endl; }

  return 0;
}
    
