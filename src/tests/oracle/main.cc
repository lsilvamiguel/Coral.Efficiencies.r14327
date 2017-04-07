#include "CsInit.h"
#include "CsEvent.h"
#include "CsStopwatch.h"

#include "CsRegistrySing.h"
#include "CsErrLog.h"
#include "Coral.h"

#include "CoralUser.h"

using namespace std;

int main( int argc, char *argv[] ) 
{
  try
    {

      if(!CoralUserSetup(argc, argv))
	{
	  cerr << "main(): CoralUserSetup() failed." << endl;
	  return -1;
	}
      
      // Package Initialization 
      Coral::init( argc, argv ); 
      
      CsEvent* event = CsEvent::Instance();
      // User Initialization  
      if(!CoralUserInit())
	{
	  cerr << "main(): CoralUserInit() failed." << endl;
	  return -2;
	}

      int nevt=0;
      // Loop on events
      while( event->getNextEvent() ) 
	{
	  if( (++nevt)%10 == 0 ) 
	    {
	      cout << "Event: " << nevt << endl;
	    }

	  if(!CoralUserEvent())
	    {
	      cerr << "main(): CoralUserEvent() failed." << endl;
	      throw;
	    }

	  if( !(CsRegistrySing::Instance()->callEoeMethods()) ) 
	    {
	      break;
	    }
	  
	}
      // End session
      CsRegistrySing::Instance()->callEndMethods();
      if(!CoralUserEnd())
	{
	  cerr << "CoralUserEnd() failed." << endl;
	}
      CsErrLog::Instance()->dump( elDebugging );
    }
  catch(std::exception &e )
    {
      cerr << "Exception:\n" << e.what() << endl;
      return -3;
    }
  catch( const string &e )
    {
      cerr << e << "\n";
    }
  catch( const char *e )
    {
      cerr << e << "\n";
    }
  catch( ... )
    {
      cerr << "Unknown exception!\n";
      return -4;
    }
  
  CS::Exception::PrintStatistics();
  return 0;
}

