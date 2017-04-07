// $Id: main.cc,v 1.10 2010/09/07 18:27:02 tnagel Exp $

#include "CsInit.h"
#include "CsEvent.h"
#include "CsGeom.h"
#include "CsStopwatch.h"
#include "CsRegistrySing.h"
#include "CsErrLog.h"
#include "CoralUser.h"
#include "Coral.h"

// Analysis Coral main program
int main( int argc, char *argv[] ) {

  try
  {

    // Package Initialization 
#warning "Coral init is not needed. But some packages instantiate it. They should not."
    Coral::init( argc, argv ); 
    // this should be preferred...
    //CsInit::Instance( argc, argv );
    CsEvent* event = CsEvent::Instance();
    // User Initialization 
    CoralUserInit();

    int nevt=0;

    // Loop on events
    while( event->getNextEvent() ) {

      if( (++nevt)%1 == 0 ) cout << "Event: " << nevt << endl;
      CoralUserEvent();
      if( !(CsRegistrySing::Instance()->callEoeMethods()) ) {
        break;
      }
    }

    // End session
    CoralUserEnd();
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


