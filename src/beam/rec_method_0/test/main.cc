// $Id: main.cc,v 1.5 2000/12/05 08:02:30 khaustov Exp $
#include "coral_config.h"
#include "Coral.h"
#include "CsEvent.h"
#include "CsGeom.h"
#include "CsStopwatch.h"
#include "CsRegistrySing.h"
#include "CsBeamRecons.h"

/*!
  \fn int main( int argc, char *argv[] ) {
  \brief The HardCoral main 
  \author Benigno Gobbo
  \version $Revision: 1.5 $
  \date $Date: 2000/12/05 08:02:30 $
*/

int main( int argc, char *argv[] ) {

  // CsStopwatch for elapsed time measurements...
  CsStopwatch stopwatches;
  int chrono = stopwatches.start(); // start a chronometer
  double time;

  // Package Initialization 
  Coral* coral        = Coral::init( argc, argv );
  int nevent=0;
  list<CsBeam*>   bmtracks;      // beam track pointer list - main results of this staff
  int thisRun_= 0;
  int previousRun_=-1;
  CsBeamRecons::Instance()->bmrecons();
//
  while( coral->getNextEvent() ) {

//    thisRun_ = store_->getRun();
//    if( previousRun_ != thisRun_ ) {
//      previousRun_ = thisRun_;
//     CsRegistrySing::Instance()->callSorMethods();
//    }


  CsBeamRecons::Instance()->bmrecons();
  bmtracks = CsBeamRecons::Instance()->getBeam();
//  list<CsBeam*>:: iterator ib;
//   for(ib=bmtracks.begin();ib!=bmtracks.end(); ib++)
//      cout << "Time="<<(*ib)->getTrackTime()<< endl;
//  nevent++;
//  if(nevent>100) break;

  }

  cout << " That's all folks..." << endl;
  coral->end();

  return 0;
}


