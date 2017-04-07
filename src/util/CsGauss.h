// $Id: CsGauss.h,v 1.2 2007/02/05 10:18:46 gobbo Exp $

/*!
   \file    CsGauss.h
   \brief   Gaussian-Distrubuted random number generator
   \author  Benigno Gobbo 
   \version $Revision: 1.2 $
   \date    $Date: 2007/02/05 10:18:46 $
*/

/*! \class CsGauss CsGauss.h
  \brief Gaussian-Distrubuted random number generator.

  This is a C++ porting of the MATHLIB RG32 function.
 */

class CsGauss {

public:

  //! creator
  inline CsGauss() { seed_ = 875949887; }

  //! set feed
  inline void setSeed( int seed ) { seed_ = seed; }

  //! get feed
  inline int  getSeed() { return( seed_ ); }

  //! get random number
  inline float random() {
    const float c = 1.1920928955078e-07;
    const int    m = 0x7fffffff; int j = 0;
    for( int i=0; i<12; i++ ){
      seed_ = (seed_ * 69069) & m;
      j = j + seed_/256;
    }
    return( c * ( ( j + 134 ) / 256 * 256 ) - 6.0 );
  }

private:

  int seed_;

};

