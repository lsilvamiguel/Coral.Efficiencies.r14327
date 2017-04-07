// $Id: CsRegistry.cc,v 1.3 2000/03/06 15:12:03 benigno Exp $

/*!
   \file    CsRegistry.cc
   \brief   Compass Registry Class
   \author  Benigno Gobbo 
   \version $Revision: 1.3 $
   \date    $Date: 2000/03/06 15:12:03 $
*/

#include "CsRegistry.h"

#include "CsRegistrySing.h"

CsRegistry::CsRegistry() {
}

bool CsRegistry::EOJRegistration( CsEndOfJob* ptr ) {
  CsRegistrySing* reg = CsRegistrySing::Instance();
  bool status = reg->EOJRegistration( ptr );
  return status;
}

bool CsRegistry::EOERegistration( CsEndOfEvent* ptr ) {
  CsRegistrySing* reg = CsRegistrySing::Instance();
  bool status = reg->EOERegistration( ptr );
  return status;
}

bool CsRegistry::SORRegistration( CsStartOfRun* ptr ) {
  CsRegistrySing* reg = CsRegistrySing::Instance();
  bool status = reg->SORRegistration( ptr );
  return status;
}





