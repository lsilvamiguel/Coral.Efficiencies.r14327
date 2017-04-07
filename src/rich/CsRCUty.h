  #ifndef  RCUTY_H
  #define  RCUTY_H

/*!
   \file    CsRCUty.h
   \brief   Execution Keys for CsRichOne.
   \author  Paolo Schiavon
   \version 0.02
   \date    June  2000
*/  


  #include <CLHEP/Vector/ThreeVector.h>
  #include "CLHEP/Matrix/Matrix.h"

  class CsRCUty {


    public:

  static CsRCUty* Instance();

  void print() const;

  double zRotMWR( const CLHEP::Hep3Vector&, const CLHEP::Hep3Vector& );
  double yRotMWR( const CLHEP::Hep3Vector&, const CLHEP::Hep3Vector& );

  void printVL( char*, int, int* );
  void printVL( char*, int, float* );
  void printVL( char*, int, double* );

  CLHEP::Hep3Vector rotfbcc( float, CLHEP::Hep3Vector&, CLHEP::Hep3Vector&, double&, double& );
  double patg( const double, const double );
  double psatg( const double, const double );


    protected:

  CsRCUty();

  ~CsRCUty();


    private:


  static CsRCUty* instance_;


  };

#endif
