/*!
   \file    CsRCMIrrorNom.cc
   \------------------------
   \brief   CsRCMIrrorNom class implementation.
   \author  Paolo Schiavon
   \version 0.01
   \date    25 august 2000
*/

  #include <ostream>
  #include <cstdio>

//-----------------------------
  #include "CsRCMirrorNom.h"
  #include "CsRCMirrors.h"
//-----------------------------

  using namespace std;

//===========================================================================
  CsRCMirrorNom::CsRCMirrorNom() {
//--------------------------------
    name_ = " ";
    vC0_.set(0., 0., 0.);
    RRnom_ = 0.;
    RR_ = 0.;
    mirpo_ = ' ';
  }


//===========================================================================
  CsRCMirrorNom::CsRCMirrorNom( const int kMir, const string name,
//----------------------------------------------------------------
                                const double xC0, const double yC0,
				const double zC0, const double RR ) {

//- positions in mm, MRS
//  --------------------
    kMir_ = kMir;
    name_ = name;
    vC0nom_.setX( xC0 );
    vC0nom_.setY( yC0 );
    vC0nom_.setZ( zC0 );
    vC0_.setX( xC0 );
    vC0_.setY( yC0 );
    vC0_.setZ( zC0 );
    RRnom_ = RR;
    RR_ = RR;
    mirpo_ = name[0];

//- read average mirror centre value and compute new centre position
//  MAIN RS
//----------------------------------------------------------  030120
    double thetaAve = 0.074;
//@@-----------------------                      provisional
    double RRc = 0;
    double xC0c = 0.;
    double yC0c = 0.;
    double zC0c = 0.;
    CsOpt* opt = CsOpt::Instance();
    vector<float> vflo;
    bool boo = false;
    if( kMir == 0 ) {
      boo = opt->CsOpt::getOpt( "RICHONE", "MirrorUP", vflo );
      if( boo ) {
	RRc  = vflo[0];
	double xpos = xC0 + 0.;
	double ypos = yC0 + RR * sin(   thetaAve );
	double zpos = zC0 + RR * cos( - thetaAve );
	xC0c = xpos - 0;
	yC0c = ypos - RRc * sin(   thetaAve );
	zC0c = zpos - RRc * cos( - thetaAve );
	vC0_.setX( xC0c );
	vC0_.setY( yC0c );
	vC0_.setZ( zC0c );
	RR_ = RRc;
	cout << "RICHONE, CsRCMirrorNom::CsRCMirrorNom() :";
	cout << " average mirror UP   radius  "
	     << RR_ << endl;
      }
    }
    if( kMir == 1 ) {
      boo = opt->CsOpt::getOpt( "RICHONE", "MirrorDOWN", vflo );
      if( boo ) {
	RRc  = vflo[0];
	double xpos = xC0 + 0.;
	double ypos = yC0 + RR * sin(   thetaAve );
	double zpos = zC0 + RR * cos(   thetaAve );
	xC0c = xpos - 0;
	yC0c = ypos - RRc * sin(   thetaAve );
	zC0c = zpos - RRc * cos(   thetaAve );
	vC0_.setX( xC0c );
	vC0_.setY( yC0c );
	vC0_.setZ( zC0c );
	RR_ = RRc;
	cout << "RICHONE, CsRCMirrorNom::CsRCMirrorNom() :";
	cout << " average mirror DOWN radius  "
	     << RR_ << endl;
      }
    }

    //cout << name_ << "  " << RRnom_ << "  " << RR_ << "  " << vC0_ << endl;

  }

//===========================================================================
  CsRCMirrorNom::CsRCMirrorNom( const CsRCMirrorNom &mirr ) {
//-----------------------------------------------------------
    cout << "RICHONE : CsRCMirrorNom CopyConstructor" << endl;
    kMir_ = mirr.kMir_;
    name_ = mirr.name_;
    vC0nom_ = mirr.vC0nom_;
    vC0_ = mirr.vC0_;
    RRnom_ = mirr.RRnom_;
    RR_ = mirr.RR_;
    mirpo_ = mirr.mirpo_;
  }

//===========================================================================
  CsRCMirrorNom& CsRCMirrorNom::operator=( const CsRCMirrorNom &mirr ) {
//----------------------------------------------------------------------
    if( this != &mirr ) {
      kMir_ = mirr.kMir_;
      name_ = mirr.name_;
      vC0nom_ = mirr.vC0nom_;
      vC0_ = mirr.vC0_;
      RRnom_ = mirr.RRnom_;
      RR_ = mirr.RR_;
      mirpo_ = mirr.mirpo_;
    }
    return ( *this );
  }

//===========================================================================
  void CsRCMirrorNom::print() const {
//-----------------------------------
    static double rg = 180./3.1415926;
    cout << endl;
    cout << "Nominal mirror geometry :" << endl
         << "-------------------------" << endl;
    cout << "Mirror  " << name_ << " nr  " << kMir_
         << " :  centre = " << vC0nom_
         << " :  centre = " << vC0_
         << ",  nominal radius = " << RRnom_
         << ",  radius = " << RR_;
  }

//===========================================================================
  CsRCMirrorNom::~CsRCMirrorNom() { }
//-----------------------------------
