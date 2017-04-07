/*!
   \file    CsRCPhotonDet.cc
   \------------------------
   \brief   CsRCPhotonDet class implementation.
   \author  Paolo Schiavon
   \version 0.01
   \date    28 august 2000, rev. 21 February 2001
*/


#    include <ostream>
#    include <cstdio>

#    include "CsErrLog.h"
#    include "CsOpt.h"

//-----------------------------
#    include "CsRCPhotonDet.h"

#    include "CsRCDetectors.h"
//-----------------------------

    using namespace std;
    using namespace CLHEP;

//===========================================================================
    CsRCPhotonDet::CsRCPhotonDet( int kDet, string name,
//------------------------------------------------------
                          double xDet0, double yDet0, double zDet0, 
	                  HepMatrix rotMatrix ) {

//--- positions in mm, MRS
//    --------------------
      kDet_ = kDet;
      name_ = name;
//--- detector centre
      vDet0In_.setX( xDet0 );
      vDet0In_.setY( yDet0 );
      vDet0In_.setZ( zDet0 );
//@@---------------------------------    //   020823
      CsOpt* opt = CsOpt::Instance();
      vector<float> vflo;
      bool boo = false;
      double corrx = 0.;
      double corry = 0.;
      double corrz = 0.;
      if( kDet == 0 ) {
	boo = opt->CsOpt::getOpt( "RICHONE", "PhotonDetectorUP", vflo );
	if( boo ) {
	  corrx = vflo[0];
	  corry = vflo[1];
	  corrz = vflo[2];
	  xDet0 += corrx;
	  yDet0 += corry;
	  zDet0 += corrz;
	  Hep3Vector vCorr( corrx, corry, corrz );
	  cout << "RICHONE, CsRCPhotonDet::CsRCPhotonDet() :";
	  cout << " corrections to photon detector UP   position  "
	       << vCorr << endl;
	}
      }
      if( kDet == 1 ) {
	boo = opt->CsOpt::getOpt( "RICHONE", "PhotonDetectorDOWN", vflo );
	if( boo ) {
	  corrx = vflo[0];
	  corry = vflo[1];
	  corrz = vflo[2];
	  xDet0 += corrx;
	  yDet0 += corry;
	  zDet0 += corrz;
	  Hep3Vector vCorr( corrx, corry, corrz );
	  cout << "RICHONE, CsRCPhotonDet::CsRCPhotonDet() :";
	  cout << " corrections to photon detector DOWN position  "
	       << vCorr << endl;
	}
      }
//@@---------------------------------    //   020823   !!!
      vDet0_.setX( xDet0 );
      vDet0_.setY( yDet0 );
      vDet0_.setZ( zDet0 );
//--- detector rot matrix, angle, direction, position
      rotMatrix_ = rotMatrix;
      detAng_ = asin( rotMatrix_[1][2] );
      vDcDet_.setX(  0.);
      vDcDet_.setY( -sin( detAng_ ));
      vDcDet_.setZ(  cos( detAng_ ));
      detpo_ = name[0];
    }


//===========================================================================
    CsRCPhotonDet::CsRCPhotonDet( const CsRCPhotonDet &det ) {
//------------------------------------------------------------
      cout << "RICHONE : CsRCPhotonDet CopyConstructor" << endl;
      kDet_ = det.kDet_;
      name_ = det.name_;
      vDet0_ = det.vDet0_;
      vDet0In_ = det.vDet0In_;
      rotMatrix_ = det.rotMatrix_;
      detAng_ = det.detAng_;
      vDcDet_ = det.vDcDet_;
      detpo_ = det.detpo_;
      vDetW_ = det.vDetW_;
      dyDRSMWR_ = det.dyDRSMWR_;
    }

//===========================================================================
    CsRCPhotonDet& CsRCPhotonDet::operator=( const CsRCPhotonDet &det ) {
//-----------------------------------------------------------------------
      cout << "RICHONE : CsRCPhotonDet Operator=" << endl;
      if( this != &det ) {
        kDet_ = det.kDet_;
        name_ = det.name_;
        vDet0_ = det.vDet0_;
	vDet0In_ = det.vDet0In_;
        rotMatrix_ = det.rotMatrix_;
        detAng_ = det.detAng_;
        vDcDet_ = det.vDcDet_;
        detpo_ = det.detpo_;
	vDetW_ = det.vDetW_;
	dyDRSMWR_ = det.dyDRSMWR_;
      }
      return ( *this );
    }

//===========================================================================
    void CsRCPhotonDet::print() const {
//-------------------------------------
      //      cout << endl;
      //      cout << "Photon detector geometry :" << endl
      //           << "--------------------------" << endl;
      cout << "Photon detector  " << name_ << "  nr " << kDet_ << endl;
      cout << "detector centres (MRS) = " << vDet0In_ << endl;
      cout << "detector centres corr. (MRS) = " << vDet0_ << endl;
      cout << "detector centres (MWR) = " << vDetW_ << endl;
      cout << "dyDRSMWR    " << dyDRSMWR_ << endl;
      cout << "angle (rad) = " << detAng_ << endl;
      cout << "rotation matrix (MRS) = ";
      for( int i=0; i<3; i++ ) for( int j=0; j<3; j++ ) 
        cout << rotMatrix_(i+1,j+1) << "  ";
      cout << endl;
      cout << endl;
      CsRCDetectors *dets = CsRCDetectors::Instance();
      cout << "zEntrWind  = " << dets->zEntrWind()
	   << "   zExitWind   = " << dets->zExitWind() << endl;
      cout << endl;
    }

//===========================================================================
    CsRCPhotonDet::~CsRCPhotonDet() { }
//-------------------------------------

