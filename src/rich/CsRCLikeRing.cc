/*!
   \file    CsRCLikeRing.h
   \---------------------
   \brief   CsRCLikeRing class implementation
   \author  Paolo Schiavon
   \version 1.0
   \date    12 January 2003
*/

  #include <iostream>
  #include <ostream>
  #include <cstdio>


// --------------------------------
  #include "CsRCLikeRing.h"
// --------------------------------

  using namespace std;

//===========================================================================
  CsRCLikeRing::CsRCLikeRing() {}
//-------------------------------

//===========================================================================
  CsRCLikeRing::~CsRCLikeRing() {}
//--------------------------------

//===========================================================================
  double CsRCLikeRing::normSignal( const double theReco ) { return 0.; }
//----------------------------------------------------------------------

//===========================================================================
  double CsRCLikeRing::normBackgr( const double theReco,
//------------------------------------------------------
    const double theIpo ) { return 0.; }

//===========================================================================
  double CsRCLikeRing::likeSignal( const double thePhot,
//------------------------------------------------------
    const double theIpo, const double sigPhot ) { return 0. ; }

//===========================================================================
  double CsRCLikeRing::likeBackgr( CsRCPhoton* pPhot,
//---------------------------------------------------
    const double theIpo ) { return 0. ; }

//===========================================================================
  double CsRCLikeRing::getRingBackground( const double theReco ) {
//----------------------------------------------------------------
    return 0.; }
