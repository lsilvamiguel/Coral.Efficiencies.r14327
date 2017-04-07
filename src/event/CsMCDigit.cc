// $Id: CsMCDigit.cc,v 1.2 2003/04/10 09:38:58 benigno Exp $

/*!
   \file    CsMCDigit.cc
   \brief   Compass Monte Carlo Digit Base Class.
   \author  Benigno Gobbo
   \version $Revision: 1.2 $
   \date    $Date: 2003/04/10 09:38:58 $
*/

#include "CsMCDigit.h"

CsMCDigit::CsMCDigit() : CsDigit() {
  hits_.clear();
}

CsMCDigit::CsMCDigit(CsDet& det, int address ) : 
  CsDigit( det, address ) {
  hits_.clear();
}

CsMCDigit::CsMCDigit(CsDet& det, int address, double* data, int datasize) : 
  CsDigit( det, address, data, datasize ) {
  hits_.clear();
}

CsMCDigit::CsMCDigit( const CsMCDigit& digi ) :
  CsDigit( digi ) {
  hits_ = digi.hits_;
}


CsMCDigit& CsMCDigit::operator=( const CsMCDigit& digi ) {
  if( this != & digi ) {
    CsDigit(*this) = CsDigit::operator=(digi);
    hits_ = digi.hits_;
  }
  return( *this );
}

bool CsMCDigit::operator==( const CsMCDigit& digi ) const {
  if( CsDigit::operator==(digi) &&
      hits_ == digi.hits_ ) {
    return true;
  }
  else {
    return false;
  }
}

bool CsMCDigit::operator<( const CsMCDigit& digi ) const {
  if( CsDigit::operator<(digi) ) {
    return true;
  }
  else {
    return false;
  }
}

