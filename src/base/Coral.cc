// $Id: Coral.cc,v 1.17 2010/02/08 21:54:10 tnagel Exp $

/*!
   \file    Coral.cc
   \brief   Compass Reconstruction and AnaLysis Package.
   \author  Benigno Gobbo
   \version $Revision: 1.17 $
   \date    $Date: 2010/02/08 21:54:10 $
*/

#include "Coral.h"

using namespace std;

Coral* Coral::instance_ = 0;

Coral* Coral::init( int argc, char **argv ) {
 if( instance_ == 0 ) {
   instance_ = new Coral( argc, argv );
 }
 return instance_; 
}

Coral* Coral::Instance() {
 if( instance_ != 0 ) 
   return instance_;
 else {
   cerr << "Coral FATAL: wrong singleton usage." << endl;
   exit(1);
 }
}

Coral::Coral( int argc, char ** argv ) {

  // Instance the Compass Initializer and the Event
  CsInit::Instance( argc, argv );
  CsEvent::Instance();

}


