// $Id: CsRandom.h,v 1.7 2009/09/14 00:43:04 ybedfer Exp $

/*!
   \file    CsRandom.h
   \brief   Compass Random Number Generators Interface
   \author  Benigno Gobbo 
   \version $Revision: 1.7 $
   \date    $Date: 2009/09/14 00:43:04 $
*/

#ifndef CsRandom_h
#define CsRandom_h

#include <CLHEP/Random/Randomize.h>

/*! \class CsRandom 
    \brief   Compass Random Number Generators Interface

    Interface to CLHEP random numbers generators
*/

class CsRandom {

 public:

  //! singleton instantiation (but first).
  static CsRandom* Instance();

  //! single flat random number 
  static double flat();

  //! single gauss random number 
  static double gauss();

  //! Gaussian random numbers are generated two at the time. "falgGauss" will force generation of a new couple of values
  void flagGauss() { _gauss->setF(false); }

  //! single exponential random number 
  static double exp();

  //! single Breit Wigner random number 
  static double bw();

  //! single Poisson random number 
  static double poisson();

  //! single binomial random number 
  static double binomial();

  //! single Chi^2 random number 
  static double chi2();

  //! single gamma function random number 
  static double gamma();

  //! single StudentT random number 
  static double studT();

  //! set new seed
  void setSeed( long seed );

  //! get seed
  long getSeed() const { return( _seed ); }

 protected:

  //! constructor
  CsRandom();

 private:

  static CsRandom* _instance;

  CLHEP::HepRandomEngine* _theEngine;      //!< The Random Number Generator Engine
  long   _seed;                            //!< Engine seed (if needed)
  int    _lux;                             //!< Engine luxury (if needed)
  int    _index;                           //!< Engine index (if needed)
  
  static CLHEP::RandFlat* _flat;           //!< flat random number generator
  static CLHEP::RandGauss* _gauss;         //!< Gauss random number generator
  static CLHEP::RandExponential* _exp;     //!< exponential random number generator
  static CLHEP::RandBreitWigner* _bw;      //!< Breit Wigner random number generator
  static CLHEP::RandPoisson* _poisson;     //!< Poisson random number generator
  static CLHEP::RandBinomial* _binomial;   //!< binomial random number generator
  static CLHEP::RandChiSquare* _chi2;      //!< Chi^2 random number generator
  static CLHEP::RandGamma* _gamma;         //!< gamma random number generator
  static CLHEP::RandStudentT* _studT;      //!< StudentT random number generator
};

#endif // CsRandom_h
