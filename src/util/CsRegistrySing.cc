// $Id: CsRegistrySing.cc,v 1.8 2010/01/28 12:51:26 tnagel Exp $

/*!
   \file    CsRegistrySing.cc
   \brief   Compass Registry Singleton Class
   \author  Benigno Gobbo 
   \version $Revision: 1.8 $
   \date    $Date: 2010/01/28 12:51:26 $
*/

#include "CsRegistrySing.h"

using namespace std;

CsRegistrySing* CsRegistrySing::_instance = NULL;

CsRegistrySing* CsRegistrySing::Instance() {
 if( _instance == NULL ) {
   _instance = new CsRegistrySing();
 }
 return( _instance ); 
}

CsRegistrySing::CsRegistrySing() {
}

bool CsRegistrySing::EOJRegistration( CsEndOfJob* ptr ) {
  // Check if already registered...
  list<CsEndOfJob*>::iterator i;
  for( i=_EOJRegister.begin(); i!=_EOJRegister.end(); i++ ) {
    if( (*i) == ptr ) return false;
  }
  // OK: 
  _EOJRegister.push_back( ptr );
  return true;
}

bool CsRegistrySing::EOERegistration( CsEndOfEvent* ptr ) {
  // Check if already registered...
  list<CsEndOfEvent*>::iterator i;
  for( i=_EOERegister.begin(); i!=_EOERegister.end(); i++ ) {
    if( (*i) == ptr ) return false;
  }

  // OK: 
  _EOERegister.push_back( ptr );
  return true;
}


bool CsRegistrySing::SORRegistration( CsStartOfRun* ptr ) {
  // Check if already registered...
  list<CsStartOfRun*>::iterator i;
  for( i=_SORRegister.begin(); i!=_SORRegister.end(); i++ ) {
    if( (*i) == ptr ) return false;
  }

  // OK: 
  _SORRegister.push_back( ptr );
  return true;
}

bool CsRegistrySing::callEndMethods() {

  bool OK = true;

  list<CsEndOfJob*>::reverse_iterator i;
  if( ! _EOJRegister.empty() ) {
    for( i=_EOJRegister.rbegin(); i!=_EOJRegister.rend(); i++ ) {
      if( ! (*i)->end() ) OK = false;
    }
  }
  return OK;
}


bool CsRegistrySing::callEoeMethods() {

  bool OK = true;

  list<CsEndOfEvent*>::reverse_iterator i;
  if( ! _EOERegister.empty() ) {
    for( i=_EOERegister.rbegin(); i!=_EOERegister.rend(); i++ ) {
      if( ! (*i)->eoe() ) OK = false;
    }
  }
  return OK;
}


bool CsRegistrySing::callSorMethods() {

  bool OK = true;

  list<CsStartOfRun*>::reverse_iterator i;
  if( ! _SORRegister.empty() ) {
    for( i=_SORRegister.rbegin(); i!=_SORRegister.rend(); i++ ) {
      if( ! (*i)->sor() ) OK = false;
    }
  }
  return OK;
}



