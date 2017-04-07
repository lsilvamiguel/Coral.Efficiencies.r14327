// $Id: main.cc,v 1.10 2010/09/07 18:27:02 tnagel Exp $

#include <csignal>
#include "CsInit.h"
#include "CsEvent.h"
#include "CsGeom.h"
#include "CsStopwatch.h"
#include "CsRegistrySing.h"
#include "CsErrLog.h"
#include "CoralUser.h"
#include "Coral.h"
#include <iostream>

static bool flag_end=false;

using namespace std;

// -----------------------------------------------------------------------------

void operation_system_signal(int n)
{
  cerr << "\n"
            "============================================\n"
            "=== The program has received signal: %3d ===\n"
    "============================================\n\n" << n << endl;
  if( flag_end )
  {
    cerr << "Forcing exit.\n\n";
    exit(1);
  }
  else
    flag_end = true;
}

// -----------------------------------------------------------------------------

// Analysis Coral main program
int main( int argc, char *argv[] ) {
  //--------------------------------------------------------------------------
  // Set signal handler to allow user to abort the program.
  (void) signal(0,      (void(*)(int))operation_system_signal); /* 0 */
  (void) signal(SIGHUP, (void(*)(int))operation_system_signal); /* 1 */
  (void) signal(SIGINT, (void(*)(int))operation_system_signal); /* 2 */
  (void) signal(SIGQUIT,(void(*)(int))operation_system_signal); /* 3 */
  (void) signal(SIGALRM,(void(*)(int))operation_system_signal); /* 4 */
  (void) signal(SIGTERM,(void(*)(int))operation_system_signal); /* 5 */
  //--------------------------------------------------------------------------

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
    while( !flag_end && event->getNextEvent() ) {
      nevt++; if( nevt%100 == 0 ) cout << "Event: " << nevt << endl;
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


