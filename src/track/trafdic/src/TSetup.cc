#include <iostream>
#include "TSetup.h"

/*!
  \file
  Costructor(s), destructor and trivial methods
  of the class TSetup
*/

TSetup* TSetup::address = 0; // init static pointer

/*! 
  TSetup constructor
  calls private member function TSetup::Init().
  This function put all necessary information into TSetup data members.

*/

TSetup::TSetup():
  NMags       (nMags),
  MagType     (magType),
  MagScale    (magScale),
  MagFlag1    (magFlag1),
  MagFlag2    (magFlag2),
  MagFieldInt (magFieldInt),
  BeamTAngle  (beamTAngle),
  TargetCenter(targetCenter),
  FirstMuWallA(firstMuWallA)
{
  NobjCreated++;
  if(address == 0){ //!< Protection against multiple instancies
    address = this;

    // Setup initialization  
    if( !Init() ) {
      std::cout<<"TSetup constructor ==> Initialisation failed "<<std::endl;
      assert(false);
    }

  } else {
    std::cout<<"Only one instance of the class TSetup is allowed"<<std::endl;
    assert(false);
  }
}

//! Destructor
TSetup::~TSetup()
{
  NobjDestructed++;
  address = 0;
}










