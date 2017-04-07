#include "Coral.h"
#include <iostream>
#include "Recon.h"
#include "RecCall.h"
#include "CsRegistry.h"
#include "CsEndOfJob.h"
#include "CsStopwatch.h"
#include "RecOpt.h"

using namespace std;

extern "C"
void reconend_() ;

Recon* Recon::address = 0; // init static pointer


 // CsStopwatch for elapsed time measurements... 
    CsStopwatch stopwatches;
    int chrono5 = stopwatches.start(); // start a chronometer

    
//! Constructor (Recon package initialisation)
Recon::Recon()
{
  if(address == 0){ //!< Protection against multiple instancies
    address = this;
    time1 =0;  
    time2 =0;
    time7 =0;

   // CsStopwatch for elapsed time measurements...
    CsStopwatch stopwatches;
    int chrono6 = stopwatches.start(); // start a chronometer
 
    //initialisation of Recon
    RecCall  ini;
    ini.ReconIni();
float time6 = stopwatches.stop(chrono6); // stop a chronometer
Recon::ref().countTime6(time6);
 
    // Register for "end of job" call
    CsRegistry reg;                
    reg.EOJRegistration( this );


    // Register for "end of event" call
    reg.EOERegistration( this );


  }  
}


//! Destructor
Recon::~Recon()
{
  address = 0;
}

// "End Of Event" method
bool Recon::eoe() {

  bool print = false;

  if(print){ // EOF statistics


}


  

  return( true );
}

 


// "End of Job" method 
bool Recon::end() {
   reconend_() ;
  cout<<"-----------> Recon End-Of-Job"<<endl;

float time5 = stopwatches.stop(chrono5); // stop a chronometer

 cout<<"Total time in Coral  = "<<time5 <<"  sec "<<endl;
 cout<<"Initialization time in Recon in Coral  = "<<time6 <<"  sec "<<endl;
 cout<<"Total time in Recon with intefaces = "<<time7 <<"  sec "<<endl;
 cout<<"Total time in Recon (MC)  = "<<time1 <<"  sec "<<endl;
 cout<<"Total time in Recon reconstruction module  = "<<time2 <<"  sec "<<endl;

 
 // if(address != NULL) delete this;   
    return( true );
}






