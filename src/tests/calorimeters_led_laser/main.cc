// $Id: main.cc,v 1.53 2010/09/07 18:27:02 tnagel Exp $

#include "DaqDataDecoding/Exception.h"
#include "CsInit.h"
#include "CsEvent.h"
#include "CsGeom.h"
#include "CsStopwatch.h"
#include "CsRegistrySing.h"
#include "CsErrLog.h"
#include "CoralUser.h"
#include "Coral.h"
#include "CsOpt.h"
#include "ObjectsCounter.h"

// Analysis Coral main program
int main( int argc, char *argv[] ) {

  try
  {  
    // Package Initialization 
    ///\todo "Coral init is not needed. But some packages instantiate it. They should not."
    Coral::init( argc, argv ); 
    // this should be preferred...
    //CsInit::Instance( argc, argv );
    CsInit* init = CsInit::Instance();
    CsOpt *opt = CsOpt::Instance(argc,argv); // Instantiate the Option Interpreter
    CsEvent* event = CsEvent::Instance();
    unsigned int print_each = 100; 
    std::string tag = "";  std::string key = "events print each"; int n;
    if( opt->getOpt( tag, key, n ) ) {
      print_each = (unsigned int) n;
    }
    // User Initialization 
    CoralUserInit();
    
    int nevt=0;
    
    // Loop on events    
    while( event->getNextEvent() ) { 
      nevt++;
      if((nevt)%print_each == 0 )
 	std::cout << "Event: " << nevt << std::endl;
      CoralUserEvent();

      if( nevt >= (int)init->getMaxEvents() ) break;

      if( !(CsRegistrySing::Instance()->callEoeMethods()) ) {
	break;
      }
    }
    // End session
    CoralUserEnd();
    CsRegistrySing::Instance()->callEndMethods();
    CsErrLog::Instance()->dump( elDebugging );
    
	               
  }  
  catch(std::exception &e ) { std::cerr << "Exception:\n" << e.what() << std::endl; }
  catch(std::string &e) { std::cerr << "Exception:\n" << e << "\n"; }
  catch(const char *e) { std::cerr << "Exception:\n" << e << "\n"; }
  catch( ... ) { std::cerr << "Unknown exception!!\n"; }

  //  CS::Exception::PrintStatistics();

  ObjectsCounterMaster::SetStream(&std::clog);
  ObjectsCounterMaster::Print();

  return 0;
}
