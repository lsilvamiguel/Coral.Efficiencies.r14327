// $Id: main.cc,v 1.7 2010/09/07 18:27:02 tnagel Exp $
#include "Coral.h"
#include "CsInit.h"
#include "CsEvent.h"
#include "CsGeom.h"
#include "CsStopwatch.h"
#include "CsRegistrySing.h"
#include "CsErrLog.h"

// Analysis Coral main program
int main( int argc, char *argv[] ) {

 try
 {

    // Package Initialization 
    CsInit::Instance( argc, argv );
    Coral::init( argc, argv );
    CsEvent* event = CsEvent::Instance();

    // Loop on events
    while( event->getNextEvent() ) {
      if( !(CsRegistrySing::Instance()->callEoeMethods()) ) {
        break;
      }
    }

    // End session
    CsRegistrySing::Instance()->callEndMethods();
    CsErrLog::Instance()->dump( elDebugging );
  }
  catch(std::exception &e )
 {
  cerr << "Exception:\n" << e.what() << endl;
 }
 catch( ... )
  {
   cerr << "Unknown exception!\n";
  }

  return 0;
}


