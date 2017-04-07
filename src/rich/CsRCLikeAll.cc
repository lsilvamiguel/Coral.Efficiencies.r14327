
/*!
   \file    CsRCLikeAll.h
   \---------------------
   \brief   CsRCLikeAll class implementation
   \author  Paolo Schiavon
   \version 1.0
   \date    12 January 2003
*/

  #include <iostream>
  #include <ostream>
  #include <cstdio>


// --------------------------------
  #include "CsRCLikeAll.h"
// --------------------------------

  using namespace std;

//===========================================================================
  CsRCLikeAll::CsRCLikeAll() {}
//-----------------------------

//===========================================================================
  CsRCLikeAll::~CsRCLikeAll() {}
//------------------------------

//===========================================================================
  double CsRCLikeAll::normSignal( const double theIpo ) { return 0.; }
//--------------------------------------------------------------------

//===========================================================================
  double CsRCLikeAll::normBackgr( const double theIpo ) { return 0.; }
//--------------------------------------------------------------------

//===========================================================================
  double CsRCLikeAll::likeSignal( const CsRCPhoton* pPhot,
//--------------------------------------------------------
    const double theIpo ) { return 0.; }

//===========================================================================
  double CsRCLikeAll::likeBackgr( const CsRCPhoton* pPhot,
//--------------------------------------------------------
    const double theIpo ) { return 0.; }

//===========================================================================
  double CsRCLikeAll::getRingBackground( const double theReco ) {
//---------------------------------------------------------------
    return 0.; }
