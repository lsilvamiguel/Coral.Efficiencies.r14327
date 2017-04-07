/*!
   \file    CsRCUty.cc
   \brief   Utylities for CsRichOne.
   \author  Paolo Schiavon
   \version 0.02
   \date    June  2000
*/


#include <iostream>
#include <ostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "CsRCUty.h"

using namespace std;
using namespace CLHEP;

  CsRCUty* CsRCUty::instance_ = 0;

//===========================================================================
  CsRCUty* CsRCUty::Instance() {
//------------------------------
    if( instance_ == 0 ) instance_ = new CsRCUty();
    return instance_; 
  }

//===========================================================================
  CsRCUty::CsRCUty() { }
//----------------------

//===========================================================================
  void CsRCUty::print() const {
//-----------------------------
    cout << "   " << endl;
    cout << "Utylities :" << endl;
    cout << "-----------" << endl;
    cout << "   " << endl;
  }

//===========================================================================
  CsRCUty::~CsRCUty() { }
//-----------------------

//===========================================================================
  double CsRCUty::zRotMWR( const Hep3Vector &r,
//------------------------------------------------
                              const Hep3Vector &det ) {

//--- Paolo  -  May 1999.

      return  r.y() * det.y() + r.z() * det.z();  }



//===========================================================================
  double CsRCUty::yRotMWR( const Hep3Vector &r,
//------------------------------------------------
                              const Hep3Vector &det ) {

//--- Paolo  -  May 1999.

      return  - r.z() * det.y() + r.y() * det.z() ;  }


//===========================================================================
  void CsRCUty::printVL( char *name, int l, int *ve ) {
//-----------------------------------------------------
      printf("%s",name);
      for( int k=0; k < l; k++) { printf("   %i",ve[k]); }
      printf("\n");  }

//===========================================================================
  void CsRCUty::printVL( char *name, int l, float *ve ) {
//-------------------------------------------------------
      printf("%s",name);
      for( int k=0; k < l; k++) { printf("   %f",ve[k]); }
      printf("\n");  }

//===========================================================================
  void CsRCUty::printVL( char *name, int l, double *ve ) {
//--------------------------------------------------------
      printf("%s",name);
      for( int k=0; k < l; k++) { printf("   %f",ve[k]); }
      printf("\n");  }


//==========================================================================
  Hep3Vector CsRCUty::rotfbcc( float vers, Hep3Vector &vro, Hep3Vector &vin,
//-----------------------------------------------------------------------
                               double &theta, double &phi )             {

//- Paolo - november 2000

//- from drotzfb
//- vectors vro and vin need not to be unit vectors.

//- vers = +1 : direct rotation. :
//  rotates from reference S to S';
//  vin is the vector in S, vout in S';
//  vro is a vector (direction) defining (in S) the z'-axis of S';
//  the z-axis of S is rotated in two steps: a rotation around the
//  x-axis of S by an angle beta (S->S-temp) and then a rotation around
//  the y-axis of S-temp by an angle chi (S-temp->S');
//  this leaves the y-axis of S' in the yz-plane of S.

//- vers = -1 : inverse rotation. :
//  rotates from reference S' to S;
//  vin is the vector in S', vout in S;
//  vro is a vector (direction) defining (in S) the z'-axis of S';
//  the z'-axis of S' is rotated in two steps: a rotation around
//  the y'-axis of S' by an angle -chi (S'->S-temp) and then 
//  a rotation around the x-axis of S-temp by an angle -beta
//  (S-temp->S).


    double llo = 1.e+09;
    double mmo = 1.e+09;
    double nno = 1.e+09;
    theta = 1.e+09;
    phi   = 1.e+09;
    if( vro.mag2() != 0. ) {

    Hep3Vector vrou = vro.unit();
    double LL = vrou.x();
    double MM = vrou.y();
    double NN = vrou.z();
    double sqMN = sqrt( MM*MM + NN*NN );

    double ll = vin.x();
    double mm = vin.y();
    double nn = vin.z();
    double proj = ll*LL + mm*MM + nn*NN;

    if( vers == +1. ) {
      if( sqMN == 0 ) {
        llo = - nn;
        mmo =   mm;
        nno =   ll;
      } else {
        llo = (ll - LL*proj) / sqMN;
        mmo = (mm*NN - nn*MM) / sqMN;
        nno = proj;
      }
    }

    if( vers == -1. ) {
      if( sqMN == 0 ) {
        llo =   nn;
        mmo =   mm;
        nno = - ll;
      } else {
        llo = ll*sqMN + nn*LL;
        mmo = ( mm*NN - ll*LL*MM + nn*MM*sqMN) / sqMN;
        nno = (-mm*MM - ll*LL*NN + nn*NN*sqMN) / sqMN;
      }
    }
    double mag = sqrt( llo*llo + mmo*mmo + nno*nno );
    theta = acos( nno/mag );
    //phi = patg( llo, mmo );           //   ???  
    phi = patg( mmo, llo );             //   040420

    }
    Hep3Vector vout( llo, mmo, nno );

    return vout;

  }


//==========================================================================
  double CsRCUty::patg( const double sin, const double cos ) {
//------------------------------------------------------------

//- Paolo - november 2000

    static double PI = 3.1415926535;

    double angle = 0.;
    double sq = 0.;
    if( sin < 0. ) sq = 1.;
    if( cos >  0. ) { angle = atan( sin/cos ) + sq*2.*PI; }
    if( cos == 0. ) { angle = PI/2. + sq*PI; }
    if( cos <  0. ) { angle = atan( sin/cos ) + PI; }

//- angle in rads
    return angle;

  }


//==========================================================================
  double CsRCUty::psatg( const double sini, const double cosi ) {
//---------------------------------------------------------------

//- Paolo - April  2010

    static double PI = 3.1415926535;

    double angle = 0.;
    double norm = sqrt( sini*sini + cosi*cosi );
    double sin = sini/norm;
    double cos = cosi/norm;
    if( fabs( sin ) <= fabs( cos ) ) angle = atan( fabs( sin/cos ) );
    else  angle = PI/2. - atan( fabs( cos/sin ) );
    if( cos >  0.  &&  sin < 0. ) angle = 2.*PI - angle;
    if( cos <  0.  &&  sin > 0. ) angle = PI - angle;
    if( cos <  0.  &&  sin < 0. ) angle = PI + angle;

//- angle in rads
    return angle;

  }
