// $Id: CsRandom.cc,v 1.9 2008/11/11 23:03:27 ybedfer Exp $

/*!
   \file    CsRandom.cc
   \brief   Compass Random Number Generators Interface
   \author  Benigno Gobbo 
   \version $Revision: 1.9 $
   \date    $Date: 2008/11/11 23:03:27 $
*/

#include <fstream>
#include <sstream>

#include "CsRandom.h"
#include "CsOpt.h"
#include "CsErrLog.h"

using namespace CLHEP;

CsRandom* CsRandom::_instance = NULL;

RandFlat*        CsRandom::_flat     = NULL;
RandGauss*       CsRandom::_gauss    = NULL;
RandExponential* CsRandom::_exp      = NULL;
RandBreitWigner* CsRandom::_bw       = NULL; 
RandPoisson*     CsRandom::_poisson  = NULL;
RandBinomial*    CsRandom::_binomial = NULL;
RandChiSquare*   CsRandom::_chi2     = NULL;
RandGamma*       CsRandom::_gamma    = NULL;
RandStudentT*    CsRandom::_studT    = NULL;

CsRandom* CsRandom::Instance() {
  if( _instance == NULL ) {
    _instance = new CsRandom();
  }
  return( _instance );
}

CsRandom::CsRandom( ) {


  CsOpt* opt = CsOpt::Instance();

  std::string tag = "random number";
  std::string key;
  std::string str;
  int n;

  key = "luxury";
  if( opt->getOpt( tag, key, n ) ) {
    _lux = n;
  }
  else {
    _lux = 3;
  }

  key = "index";
  if( opt->getOpt( tag, key, n ) ) {
    _index = n;
  }
  else {
    _index = 1;
  }

  key = "seed";
  if( opt->getOpt( tag, key, n ) ) {
    _seed = n;

    // special seed 0: following general convention, truely random seed is
    // chosen, read from Linux kernel entropy source
    if (n == 0) {
      std::ifstream f("/dev/random");
      f.read( (char *)(&_seed), sizeof(_seed) );
      if ( f.bad() )
        CsErrLog::mes( elFatal, "Error reading from /dev/random" ); 
      f.close();
      std::stringstream seedstr; seedstr << _seed;
      CsErrLog::mes( elError, std::string("seeding random generator ")
                     + "with entropy from /dev/random: " + seedstr.str()
                     + " (this is not an error)"); 
    }
  }
  else {
    _seed = 19990101;
  }

  key = "engine";
  if( opt->getOpt( tag, key, str ) ) {
    if( str == "JamesEngine" ) {
      // make sure that seed is in range [0,900000000] which is recommended
      // by CLHEP (cf. JamesRandom.cc in CLHEP versions 1.8.1.0, 2.0.3.2)
      if (_seed < 0) _seed = -_seed;
      while (_seed > 900000000) _seed >>= 1;
      _theEngine = new HepJamesRandom( _seed );
    }
    else if( str == "RandEngine" ) {
      _theEngine = new RandEngine( _seed );
    }
    else if( str == "DRand48Engine" ) {
      _theEngine = new DRand48Engine( _seed );
    }
    else if( str == "RanluxEngine" ) {
      _theEngine = new RanluxEngine( _seed, _lux );
    }
    else if( str == "RanecuEngine" ) {
      _theEngine = new RanecuEngine( _index );
      _theEngine->setSeed( _seed, _lux );
    }
    else if( str == "MTwistEngine" ) {
      _theEngine = new MTwistEngine( _seed );
    }
    else {
      CsErrLog::mes( elFatal, "Unknown random engine '" + str);
    }
  } else CsErrLog::mes( elFatal, "No random engine specified");

  _flat     = new RandFlat( _theEngine );
  _gauss    = new RandGauss( _theEngine );
  _exp      = new RandExponential( _theEngine );
  _bw       = new RandBreitWigner( _theEngine );
  _poisson  = new RandPoisson( _theEngine );
  _binomial = new RandBinomial( _theEngine );
  _chi2     = new RandChiSquare( _theEngine );
  _gamma    = new RandGamma( _theEngine );
  _studT    = new RandStudentT( _theEngine );

}

void CsRandom::setSeed( long seed ) {
  _seed = seed;
  _theEngine->setSeed( _seed, _lux );
}


double CsRandom::flat() { 
  if( _flat != NULL ) {
    return( _flat->fire() ); 
  }
  else {
    CsErrLog::mes( elFatal, "call before 1st Singleton instantiation" ); 
    return 0;
  }
}

double CsRandom::gauss() { 
  if( _gauss != NULL ) {
    return( _gauss->fire() ); 
  }
  else {
    CsErrLog::mes( elFatal, "call before 1st Singleton instantiation" ); 
    return 0;
  }
}

double CsRandom::exp() { 
  if( _exp != NULL ) {
    return( _exp->fire() ); 
  }  
  else {
    CsErrLog::mes( elFatal, "call before 1st Singleton instantiation" ); 
    return 0;
  }
}

double CsRandom::poisson() { 
  if( _poisson != NULL ) {
    return( _poisson->fire() ); 
  }
  else {
    CsErrLog::mes( elFatal, "call before 1st Singleton instantiation" ); 
    return 0;
  }
}

double CsRandom::binomial() { 
  if( _binomial != NULL ) {
    return( _binomial->fire() ); 
  }
  else {
    CsErrLog::mes( elFatal, "call before 1st Singleton instantiation" ); 
    return 0;
  }
}

double CsRandom::chi2() { 
  if( _chi2 != NULL ) {
    return( _chi2->fire() ); 
  }
  else {
    CsErrLog::mes( elFatal, "call before 1st Singleton instantiation" ); 
    return 0;
  }
}

double CsRandom::gamma() { 
  if( _gamma != NULL ) {
    return( _gamma->fire() ); 
  }
  else {
    CsErrLog::mes( elFatal, "call before 1st Singleton instantiation" ); 
    return 0;
  }
}

double CsRandom::studT() { 
  if( _studT != NULL ) {
    return( _studT->fire() ); 
  }
  else {
    CsErrLog::mes( elFatal, "call before 1st Singleton instantiation" ); 
    return 0;
  }
}
