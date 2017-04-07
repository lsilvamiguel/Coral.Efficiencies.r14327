// $Id: CsDigit.cc,v 1.2 2003/04/10 09:38:58 benigno Exp $

/*!
   \file    CsDigit.cc
   \brief   Compass Digit class.
   \author  Benigno Gobbo
   \version $Revision: 1.2 $
   \date    $Date: 2003/04/10 09:38:58 $
*/

#include "CsDigit.h"

#include <cstddef>

CsDigit::CsDigit() : 
  address_(0), 
  datasize_(0), data_(NULL),
  det_(0), 
  clusterized_(false) {
}

CsDigit::CsDigit( CsDet& det, int address ) 
  : address_(address), datasize_(0), data_(NULL), det_(&det), clusterized_(false)
   {
}
 
CsDigit::CsDigit( CsDet& det, int address, double* data, int datasize ) 
  : address_(address), datasize_(datasize), data_(NULL), det_(&det), clusterized_(false) 
{
  if( datasize_ > 0 ) {
    data_ = new double[datasize];
    for( int i=0; i<datasize_; i++ ) data_[i] = data[i];
  }
}

CsDigit::~CsDigit() {
  if( datasize_ > 0 ) delete []data_;
} 

CsDigit::CsDigit( const CsDigit& digi ) :
  address_(digi.address_), datasize_(digi.datasize_), data_(NULL),
  det_(digi.det_), clusterized_(digi.clusterized_)
{
  if( datasize_>0 )
  {
    data_ = new double[datasize_];
    for( int i=0; i<datasize_; i++ )
      data_[i] = digi.data_[i];
  }
}

CsDigit& CsDigit::operator=( const CsDigit& digi )
{
  if( this      != & digi )
  {
    det_         = digi.det_;
    address_     = digi.address_;
    datasize_    = digi.datasize_;
    clusterized_ = digi.clusterized_;
    delete [] data_;

    if( datasize_>0 )
    {
      data_ = new double[datasize_];
      for(int i=0; i<datasize_; i++ )
        data_[i] = digi.data_[i];
    }
    else
      data_ = NULL;
  }
  return( *this );
}
  
bool CsDigit::operator==( const CsDigit& digi ) const {
  if( det_         == digi.det_      &&
      address_     == digi.address_  &&
      datasize_    == digi.datasize_ &&
      clusterized_ == digi.clusterized_ ) {
    bool equal = true;
    for( int i=0; i<datasize_; i++ ) {
      equal = equal && ( data_[i] == digi.data_[i] );
    }
    if( equal ) {
      return true;
    }
    else {
      return false;
    }
  }
  else { 
    return false;
  }
}

bool CsDigit::operator<( const CsDigit& digi ) const {
  if( det_ < digi.det_ ) {
    return( true );
  }
  else if( det_ == digi.det_ ) {
    if( address_ < digi.address_ ) {
      return( true );
    }
    else {
      return( false );
    }
  }

  else {
    return( false );
  }
}

double CsDigit::getDatum() const { 
  if( datasize_ > 0 ) {
    return( data_[0] ); 
  }
  else {
    return 0;
  }
}

void CsDigit::replaceDatum( double datum ) {
  if( datasize_ == 0 ) {
    data_ = new double[1];
    datasize_ = 1;
  }
  data_[0] = datum;    
}
