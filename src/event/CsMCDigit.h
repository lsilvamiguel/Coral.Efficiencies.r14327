// $Id: CsMCDigit.h,v 1.4 2006/11/29 06:19:02 ybedfer Exp $

/*!
   \file    CsMCDigit.h 
   \brief   Compass Monte Carlo Digit Class.
   \author  Benigno Gobbo
   \version $Revision: 1.4 $
   \date    $Date: 2006/11/29 06:19:02 $
*/

#ifndef CsMCDigit_h
#define CsMCDigit_h

#include <list>
#include "CsDigit.h"
#include "CsMCHit.h"


/*! \class CsMCDigit 
    \brief Compass Monte Carlo Digit class.
*/

class CsMCDigit : public CsDigit {

 public:

  //! Default Constructor
  CsMCDigit(); 

  //! Constructor
  CsMCDigit( CsDet& det, int address ); 

  //! Constructor
  CsMCDigit( CsDet& det, int address, double* data, int datasize = 1 ); 

  //! Copy Constructor
  CsMCDigit( const CsMCDigit& );            

  //! Assignment Operator
  CsMCDigit& operator=( const CsMCDigit& ); 

  //! "equal to" Operator
  bool operator==( const CsMCDigit& ) const;  

  //! "less than" Operator
  bool operator<( const CsMCDigit& ) const;  

  //! Returns the list of pointers to the MC true Hits associated to the Digit.
  std::list<CsMCHit*> getHits() const { return( hits_ ); }

  //! add an Hit to this Digit
  void addHit( CsMCHit& hit ) { hits_.push_back( &hit ); }

  //! Sets left/right info.
  void setLR(int lr) { lr_ = lr; }
  //! Returns left/right info assigned to the hit.
  int getLR() const { return lr_; }

 private:

  std::list<CsMCHit*> hits_;  //!< list of pointers to the associated MC true Hits
  int lr_;  //!< Left/right info derived from MC truth bank.

};

#endif // CsMCDigit_h
