/*!
   \file    CsRCParticle.cc
   \-------------------------
   \brief   CsRCParticle class implementation.
   \author  Paolo Schiavon
   \version 0.01
   \date    1 October 2000
*/


  #include <ostream>
  #include <cstdio>


  #include "CsTrack.h"
  #include "CsMCTrack.h"

//-----------------------------
  #include "CsRCParticle.h"

  #include "CsRCExeKeys.h"
  #include "CsRCRecConst.h"
  #include "CsRCDetectors.h"
//-----------------------------

  #include <CLHEP/Vector/ThreeVector.h>
//-------------------------------------

  using namespace std;
  using namespace CLHEP;

//===========================================================================
  CsRCParticle::CsRCParticle() {
//------------------------------
    kPart_ = -1;
    pTrack_ = NULL;
    mom_ = 0.;
    flag_ = false;

    vPade_.clear();
    ddPaDet_.clear();
    thPade_.clear();
    kDetPart_ = -1;
    pMirPart_ = NULL;
    pRing_ = NULL;
    vPoPaMir0_.clear();
    pMCTrack_ = NULL;
    flagS_ = false;
  }

//- MC: myfile tracks (type 0 = OLD)
//===========================================================================
  CsRCParticle::CsRCParticle( int kPart, Hep3Vector vPosIn, Hep3Vector vDirIn,
//---------------------------------------------------------------------------
			      Hep3Vector vPosOut, double mom,
			      int iTrack, int iParTy, int ioExWd ) {
    kPart_ = kPart;
    pTrack_ = NULL;
    vPosIn_ = vPosIn;
    vDirIn_ = vDirIn;
    vPosOut_ = vPosOut;
    vPosEmP_.set(1.e+06, 0., 0.);
    vDirEmP_.set(1.e+06, 0., 0.);
    vPosExW_.set(1.e+06, 0., 0.);
    vDirExW_.set(1.e+06, 0., 0.);
    mom_ = mom;
    flag_ = true;

    thPamir_ = 1.e+06;
    vPade_.clear();
    ddPaDet_.clear();
    thPade_.clear();
    kDetPart_ = -1;
    pMirPart_ = NULL;
    pRing_ = NULL;

    vPoPaMir0_.clear();
    vCorPoPa0_.set(0., 0., 0.);            //   101220
    corPaRR_ = 0.;              //   101220

    pMCTrack_ = NULL;
    iTrack_ = iTrack;
    iPartTy_ = iParTy;
    ioExWd_ = ioExWd;
    charge_ = 0;
    for ( int kPaTy=2; kPaTy <= 14; kPaTy+=3)  {
      if ( iParTy == kPaTy   )  { charge_ = + 1;  break; }
      if ( iParTy == kPaTy+1 )  { charge_ = - 1;  break; }
    }
    flagS_ = true;

    nPhoCer_ = 0;
    thetaCer_ = 0.;
    MCmTime_ = 0.;

    zHelix0_ = 0.;
    zHelix1_ = 0.;
    nClus_ = 0;
    chiSq_ = 0.;
  }

//- MC: myfile tracks (type 1 = NEW or >1)
//===========================================================================
  CsRCParticle::CsRCParticle( int kPart, Hep3Vector vPosIn, Hep3Vector vDirIn,
//---------------------------------------------------------------------------
			      Hep3Vector vPosOut,
			      Hep3Vector vPosEmP, Hep3Vector vDirEmP,
			      Hep3Vector vPosExW, Hep3Vector vDirExW,
			      double mom,
			      int npho, double theta,
			      int iTrack, int iParTy, int ioExWd ) {
    kPart_ = kPart;
    pTrack_ = NULL;
    vPosIn_ = vPosIn;
    vDirIn_ = vDirIn;
    vPosOut_ = vPosOut;
    vPosEmP_ = vPosEmP;
    vDirEmP_ = vDirEmP;
    vPosExW_ = vPosExW;
    vDirExW_ = vDirExW;
    mom_ = mom;
    flag_ = true;

    thPamir_ = 1.e+06;
    vPade_.clear();
    ddPaDet_.clear();
    thPade_.clear();
    kDetPart_ = -1;
    pMirPart_ = NULL;
    pRing_ = NULL;

    vPoPaMir0_.clear();
    vCorPoPa0_.set(0., 0., 0.);            //   101220
    corPaRR_ = 0.;              //   101220

    nPhoCer_ = npho;
    thetaCer_ = theta;
    MCmTime_ = 0.;

    pMCTrack_ = NULL;
    iTrack_ = iTrack;
    iPartTy_ = iParTy;
    ioExWd_ = ioExWd;
    charge_ = 0;
    for ( int kPaTy=2; kPaTy <= 14; kPaTy+=3)  {
      if ( iParTy == kPaTy   )  { charge_ = + 1;  break; }
      if ( iParTy == kPaTy+1 )  { charge_ = - 1;  break; }
    }
    flagS_ = true;

    //nPhoCer_ = 0;             //   030812   !!!
    //thetaCer_ = 0.;           //   030812   !!!

    zHelix0_ = 0.;
    zHelix1_ = 0.;
    nClus_ = 0;
    chiSq_ = 0.;
  }

//- MC: MC tracks
//===========================================================================
  CsRCParticle::CsRCParticle( int kPart, CsMCTrack* pMCtrk,
//---------------------------------------------------------
			      Hep3Vector vPosIn, Hep3Vector vDirIn,
			      Hep3Vector vPosOut,
			      Hep3Vector vPosEmP, Hep3Vector vDirEmP,
			      Hep3Vector vPosExW, Hep3Vector vDirExW,
			      double mom,
                              int npho, double theta ) {
    kPart_ = kPart;
    pTrack_ = NULL;
    vPosIn_ = vPosIn;
    vDirIn_ = vDirIn;
    vPosOut_ = vPosOut;
    vPosEmP_ = vPosEmP;
    vDirEmP_ = vDirEmP;
    vPosExW_ = vPosExW;
    vDirExW_ = vDirExW;
    mom_ = mom;
    flag_ = true;

    thPamir_ = 1.e+06;
    vPade_.clear();
    ddPaDet_.clear();
    thPade_.clear();
    kDetPart_ = -1;
    pMirPart_ = NULL;
    pRing_ = NULL;

    vPoPaMir0_.clear();
    vCorPoPa0_.set(0., 0., 0.);            //   101220
    corPaRR_ = 0.;              //   101220

    pMCTrack_ = pMCtrk;
    charge_ = pMCTrack_->getParticle()->getCharge();
    iTrack_ = pMCTrack_->getGnum();
    iPartTy_ = pMCTrack_->getParticle()->getGeantNumber();
    ioExWd_ = 2;
    flagS_ = true;

    nPhoCer_ = npho;
    thetaCer_ = theta;
    MCmTime_ = 0.;

    zHelix0_ = 0.;
    zHelix1_ = 0.;
    nClus_ = 0;
    chiSq_ = 0.;
  }

//- MC: reconstructed tracks
//===========================================================================
  CsRCParticle::CsRCParticle( int kPart, CsTrack* ptrk, CsMCTrack* pMCtrk,
//------------------------------------------------------------------------
			      Hep3Vector vPosIn, Hep3Vector vDirIn,
			      Hep3Vector vPosOut,
			      Hep3Vector vPosEmP, Hep3Vector vDirEmP,
			      Hep3Vector vPosExW, Hep3Vector vDirExW,
			      double mom,
                              int npho, double theta ) {
    kPart_ = kPart;
    pTrack_ = ptrk;
    vPosIn_ = vPosIn;
    vDirIn_ = vDirIn;
    vPosOut_ = vPosOut;
    vPosEmP_ = vPosEmP;
    vDirEmP_ = vDirEmP;
    vPosExW_ = vPosExW;
    vDirExW_ = vDirExW;
    mom_ = mom;
    flag_ = true;

    thPamir_ = 1.e+06;
    vPade_.clear();
    ddPaDet_.clear();
    thPade_.clear();
    kDetPart_ = -1;
    pMirPart_ = NULL;
    pRing_ = NULL;

    vPoPaMir0_.clear();
    vCorPoPa0_.set(0., 0., 0.);            //   101220
    corPaRR_ = 0.;              //   101220

    pMCTrack_ = pMCtrk;
    charge_ = pMCTrack_->getParticle()->getCharge();
    iTrack_ = pMCTrack_->getGnum();
    iPartTy_ = pMCTrack_->getParticle()->getGeantNumber();
    ioExWd_ = 2;
    flagS_ = true;

    nPhoCer_ = npho;
    thetaCer_ = theta;
    MCmTime_ = 0.;

    zHelix0_ = 0.;
    zHelix1_ = 0.;
    nClus_ = 0;
    chiSq_ = 0.;
  }

//- Data: reconstructed tracks
//===========================================================================
  CsRCParticle::CsRCParticle( int kPart, CsTrack* ptrk,
//-----------------------------------------------------
			      Hep3Vector vPosIn, Hep3Vector vDirIn,
			      Hep3Vector vPosOut,
			      Hep3Vector vPosEmP, Hep3Vector vDirEmP,
			      Hep3Vector vPosExW, Hep3Vector vDirExW,
			      double mom, int charge ) {
    kPart_ = kPart;
    pTrack_ = ptrk;
    vPosIn_ = vPosIn;
    vDirIn_ = vDirIn;
    vPosOut_ = vPosOut;
    vPosEmP_ = vPosEmP;
    vDirEmP_ = vDirEmP;
    vPosExW_ = vPosExW;
    vDirExW_ = vDirExW;
    mom_ = mom;
    charge_ = charge;
    flag_ = true;

    thPamir_ = 1.e+06;
    vPade_.clear();
    ddPaDet_.clear();
    thPade_.clear();
    kDetPart_ = -1;
    pMirPart_ = NULL;
    pRing_ = NULL;

    vPoPaMir0_.clear();         //   101220
    vCorPoPa0_.set(0., 0., 0.); //   101220
    corPaRR_ = 0.;              //   101220

    pMCTrack_ = NULL;
    iTrack_ = 0;
    iPartTy_ = 0;
    ioExWd_ = 0;
    flagS_ = false;

    nPhoCer_ = 0;
    thetaCer_ = 0.;
    MCmTime_ = 0.;

    zHelix0_ = 0.;
    zHelix1_ = 0.;
    nClus_ = 0;
    chiSq_ = 0.;
  }

//- Copy constructor
//===========================================================================
  CsRCParticle::CsRCParticle( const CsRCParticle &part ) {
//--------------------------------------------------------
    //cout << "RICHONE : CsRCParticle CopyConstructor" << endl;
    kPart_ = part.kPart_;
    pTrack_ = part.pTrack_;
    vPosIn_ = part.vPosIn_;
    vDirIn_ = part.vDirIn_;
    vPosOut_ = part.vPosOut_;
    vPosEmP_ = part.vPosEmP_;
    vDirEmP_ = part.vDirEmP_;
    vPosExW_ = part.vPosExW_;
    vDirExW_ = part.vDirExW_;
    mom_ = part.mom_;
    charge_ = part.charge_;
    flag_ = part.flag_;

    pathLen_ = part.pathLen_;
    thPamir_ = part.thPamir_;
    vPade_ = part.vPade_;
    ddPaDet_ = part.ddPaDet_;
    thPade_ = part.thPade_;
    kDetPart_ = part.kDetPart_;
    pMirPart_ = part.pMirPart_;
    pRing_ = part.pRing_;

    vPoPaMir0_ = part.vPoPaMir0_;
    vCorPoPa0_ = part.vCorPoPa0_;       //   101220
    corPaRR_ = part.corPaRR_;           //   101220

    pMCTrack_ = part.pMCTrack_;
    iTrack_ = part.iTrack_;
    iPartTy_ = part.iPartTy_;
    ioExWd_ = part.ioExWd_;
    flagS_ = part.flagS_;

    nPhoCer_ = part.nPhoCer_;
    thetaCer_ = part.thetaCer_;
    MCmTime_ = part.MCmTime_;

    zHelix0_ = part.zHelix0_;
    zHelix1_ = part.zHelix1_;
    nClus_ = part.nClus_;
    chiSq_ = part.chiSq_;
  }

//===========================================================================
  CsRCParticle::~CsRCParticle() {
//-------------------------------
    vPade_.clear();
    ddPaDet_.clear();
    thPade_.clear();
  }

//===========================================================================
  double CsRCParticle::mass( int iPartTy ) {
//------------------------------------------
    double *massPartv = CsRCRecConst::Instance()->massPartv();
    return massPartv[iPartTy];
  }

//===========================================================================
  double CsRCParticle::paMass() {
//-------------------------------
    int kk = iPartTy_;
    if ( kk >= 200 ) { kk = kk - 200; }
    CsRCRecConst *cons = CsRCRecConst::Instance();
    double *massPartv = cons->massPartv();
    return massPartv[kk];
  }

//===========================================================================
  void CsRCParticle::print() const {
//----------------------------------
    cout << endl;
    cout << " Particle (input) " << kPart_ << " : " << endl;
    cout << "   pos. at window  in : " <<  vPosIn_ << endl;
    cout << setprecision( 5 );
    cout << "   dirc. at window in : " <<  vDirIn_ << endl;
    cout << setprecision( 2 );
    cout << "   pos. at window out : " <<  vPosOut_ << endl;
    cout << "   mom., charge and flag : " << mom_ << ",  " << charge_
	 << ",  " << flag_;
    //cout << "   pos. on det. : " << vPade_[0] << " " << vPade_[1]
    // << ",  angle with det. : " << thPade_[0] << " " << thPade_[1]
    // << ",  dist. from 'beam' on det. : " << ddPaDet_[0]
    // << " " << ddPaDet_[1] << endl;
    cout << "   z-start, z-end, clus, chi : " << zHelix0_ << ",  "
	 << zHelix1_ << ",  " << nClus_ << ",  " << chiSq_ << endl;
    cout << endl;
  }

//===========================================================================
  void CsRCParticle::printMC() const {
//------------------------------------
    cout << "   Geant track  " <<  pMCTrack_->getGnum()
	 << " , Geant particle   "
	 << pMCTrack_->getParticle()->getGeantNumber() << endl;
    cout << "   Cher.photons  " << nPhoCer_
         << "  angle  " << thetaCer_
	 << "  tkmTime  " << MCmTime_
	 << endl;
  }


//===========================================================================
  int CsRCParticle::getClosestCat() const {
//-----------------------------------------

//- Paolo - January  2007

    CsRCDetectors *dets = CsRCDetectors::Instance();

    int kCat = -1;
    double partX = vPade_[kDetPart_].x();
    double partY = vPade_[kDetPart_].y();

    list<CsRCCathode*> lCathodes = dets->lCathodes();
    list<CsRCCathode*>::iterator ic;
//- assign catNo closest to 'reflected' particle 'hit'
    double disttMn = 1000000.;
    for( ic=lCathodes.begin(); ic!=lCathodes.end(); ic++ ) {
      int kc = (*ic)->kCat();
      double distt = pow( partX - (*ic)->vOffCatW().x(), 2 ) +
        pow( partY - (*ic)->vOffCatW().y(), 2 );
      if( distt < disttMn ) {
        disttMn = distt;
        kCat = kc;
      }
    }

    return  kCat;
  }
