/*!
   \file    CsRCMirrorElem.cc
   \------------------------
   \brief   CsRCMirrorElem class implementation.
   \author  Paolo Schiavon
   \version 0.01
   \date    10 august 2000
*/

    #include <ostream>
    #include <cstdio>

//-----------------------------
    #include "CsRCMirrorElem.h"

    #include "CsRCMirrors.h"
    #include "CsRCRecConst.h"
//-----------------------------

    using namespace std;
    using namespace CLHEP;

//===========================================================================
    CsRCMirrorElem::CsRCMirrorElem() {
//------------------------------------
      kMir_ = -1;
      name_ = ' ';
      theta_ = 0.;
      phi_ = 0.;
      RR_ = 0.;
      deTheta_ = 0.;
      dePhi_ = 0.;
      delta_ = 0.;
      qfact_ = 0.;
      mirpo_ = ' ';
      mirNo_ = NULL;
      align_ = false;
      vpos_.set(0., 0., 0.);
      vC0_.set(0., 0., 0.);
      vtan_.set(0., 0., 0.);
    }


//===========================================================================
    CsRCMirrorElem::CsRCMirrorElem( const int kMir, const string name,
//-------------------------------------------------------------------
				    const double theta, const double phi,
                                    const double RR, const double deTheta,
				    const double dePhi,
                                    const double delta, const double qfact,
				    const int align ) {                      

//--- angles read in rads, positions in mm,
//    COMPASS MIRROR RS, with the origin in the centre of curvature
//    of the nominal mirrors ( up or down ), as from surveyors data :
//    ---------------------------------------------------------------

      CsRCRecConst *cons = CsRCRecConst::Instance();
      double rg = cons->RadDeg();

      kMir_ = kMir;
      name_ = name;
      theta_ = theta / rg;
      phi_ = phi / rg;
      RR_ = RR;
      deTheta_ = deTheta / rg;
      dePhi_ = dePhi / rg;
      delta_ = delta;
      qfact_ = qfact;
      align_ = false;
      //if( align == 0 ) align_ = true;
      //moved to CsRCMirrors::doAliMirrors 20/9/00

      char mirrUp = 'T';
      char mirrDw = 'B';
      CsRCMirrors* mirr = CsRCMirrors::Instance();
      list<CsRCMirrorNom*> lMirrNom = mirr->lMirrNom();
      list<CsRCMirrorNom*>::iterator in;
      if( name[1] == mirrUp ) {
        mirpo_ = 'U';
        in = lMirrNom.begin();
      }
      if( name[1] == mirrDw ) {
        mirpo_ = 'D';
        in = lMirrNom.begin(); in++;
      }
      mirNo_ = (*in);

      double RRave = (*in)->RR();
      double zpos = RRave * sin( theta_ ) * cos( phi_ );
      double xpos = RRave * sin( theta_ ) * sin( phi_ );
      double ypos = RRave * cos( theta_ );
//--- positions of the centres of the mirror element surfaces :
      vpos_.setX( xpos );
      vpos_.setY( ypos );
      vpos_.setZ( zpos );
      vpos_ += (*in)->vC0();
      //see vCePos in CsRCMirrors::doAliMirrors and 
      //              CsRCEventAnalysis::doMirrDetFit
      //cout << vpos_ << "  " << (*in)->vC0() << endl;

//--- position corrected by delta :
      double ddz = delta_ * sin( theta_ ) * cos( phi_ );
      double ddx = delta_ * sin( theta_ ) * sin( phi_ );
      double ddy = delta_ * cos( theta_ );
      Hep3Vector dpos( ddx, ddy, ddz );
      vpos_ += dpos;

//--- as from 020325
//--- warning : dePhi_(H) and deTheta_(V) are angular displacements of the
//    mirror axis
//    the mirror surveyor's RS z-axis always points up
//    ------------------------------------------------
      CsOpt* opt = CsOpt::Instance();
      bool boo = false;
      vector<float> vflo;
      boo = opt->CsOpt::getOpt( "RICHONE", name, vflo );
      if( boo ) {
	RR_ = (*in)->RRnom() + vflo[0];
	dePhi_ =  vflo[1] / 1000.;
	deTheta_ =  vflo[2] / 1000.;
#ifdef CsRCMirrorElem_PRINT_ALL
	cout << "RICHONE, CsRCMirrorElem() :";
	cout << " mirror element  " << name << " , radius and corr.s to ";
	cout << "phi and theta   " << RR_ << ", " << dePhi_ << ", "
	     << deTheta_ << endl;
        cout << "--------------------------------------------------------"
	     << "-----------------------------" << endl;
#endif
      }

//--- corrected positions of the centres of curvature of the mirror elements
//    ( MAIN RS ) :
      double theta_p = theta_ + deTheta_;                     //   020325
      double phi_p   = phi_ - dePhi_;                         //   020325
      double zc0 = vpos_.z() - RR_ * sin( theta_p ) * cos( phi_p );
      double xc0 = vpos_.x() - RR_ * sin( theta_p ) * sin( phi_p );
      double yc0 = vpos_.y() - RR_ * cos( theta_p );
      vC0_.setX( xc0 );
      vC0_.setY( yc0 );
      vC0_.setZ( zc0 );
      //cout << name_ << "  " << RRave << "  " << RR_ << "  " << vC0_ << endl;

//--- actual directions to the mirror elements ( MIRROR RS )
      double tana = vpos_.x() / vpos_.z();
      double tanb = vpos_.y() / vpos_.z();
      vtan_.setX( tana );
      vtan_.setY( tanb );
      vtan_.setZ(   1. );

    }

//===========================================================================
    CsRCMirrorElem::CsRCMirrorElem( const CsRCMirrorElem &mirr ) {
//----------------------------------------------------------------
      cout << "RICHONE : CsRCMirrorElem CopyConstructor" << endl;
      kMir_ = mirr.kMir_;
      name_ = mirr.name_;
      theta_ = mirr.theta_;
      phi_ = mirr.phi_;
      RR_ = mirr.RR_;
      deTheta_ = mirr.deTheta_;
      dePhi_ = mirr.dePhi_;
      delta_ = mirr.delta_;
      qfact_ = mirr.qfact_;
      mirpo_ = mirr.mirpo_;
      mirNo_ = mirr.mirNo_;
      align_ = mirr.align_;
      vpos_ = mirr.vpos_;
      vC0_ = mirr.vC0_;
      vtan_ = mirr.vtan_;
    }

//===========================================================================
    CsRCMirrorElem& CsRCMirrorElem::operator=( const CsRCMirrorElem &mirr ) {
//---------------------------------------------------------------------------
      if( this != &mirr ) {
	kMir_ = mirr.kMir_;
        name_ = mirr.name_;
        theta_ = mirr.theta_;
        phi_ = mirr.phi_;
        RR_ = mirr.RR_;
        deTheta_ = mirr.deTheta_;
        dePhi_ = mirr.dePhi_;
        delta_ = mirr.delta_;
        qfact_ = mirr.qfact_;
        mirpo_ = mirr.mirpo_;
        mirNo_ = mirr.mirNo_;
        align_ = mirr.align_;
        vpos_ = mirr.vpos_;
        vC0_ = mirr.vC0_;
        vtan_ = mirr.vtan_;
      }
      return ( *this );
    }

//===========================================================================
    void CsRCMirrorElem::print() const {
//--------------------------------------
      static double rg = 180./3.1415926;
      cout << endl;
      cout << "Mirror element geometry :" << endl
           << "-------------------------" << endl;
      cout << "Mirror  " << name_ << " nr  " << kMir_
           << " :  theta= " << theta_*rg
           << ",  phi= " << phi_*rg;
      cout << ",  RR= " << RR_ << ",  dTheta= " << deTheta_*rg
           << ",  dPhi= " << dePhi_*rg;
      cout << ",  delta= " << delta_ << ",  quality =" << qfact_
           << ",  align= " << align_ << endl;
    }

//===========================================================================
    CsRCMirrorElem::~CsRCMirrorElem() { }
//---------------------------------------
