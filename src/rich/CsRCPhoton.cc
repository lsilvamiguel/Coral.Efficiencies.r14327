
/*!
   \file    CsRCPhoton.cc
   \----------------------
   \brief   CsRCPhoton class implementation.
   \author  Paolo Schiavon
   \version 0.01
   \date    1 October 2000
*/


  #include <iostream>
  #include <ostream>
  #include <cstdio>
  #include <cmath>

  #include <CLHEP/Vector/ThreeVector.h>
//-------------------------------------

  #include "CsErrLog.h"

// ---------------------------
  #include "CsRCPhoton.h"

  #include "CsRCDetectors.h"
  #include "CsRCPhotonDet.h"
  #include "CsRCMirrors.h"

  #include "CsRCParticle.h"
  #include "CsRCPad.h"
  #include "CsRCCluster.h"

  #include "CsRCPartPhotons.h"
  #include "CsRCEventPartPhotons.h"

  #include "CsRCExeKeys.h"
  #include "CsRCRecConst.h"
  #include "CsRCHistos.h"

  #include "CsRCUty.h"
// ---------------------------

  using namespace std;
  using namespace CLHEP;

//========================================================================
  CsRCPhoton::CsRCPhoton() {
//--------------------------
    kPhot_ = -1;
    the0_ = 0.;
    the_ = 0.;
    theNorm_ = 0.;
    phi_ = 0.;
    phiA_ = 0.;
    theB_ = 0.;
    theM_ = 0.;
    PH_ = 0.;
    ptToClu_ = NULL;
    ptToDet_ = NULL;
    pPartPhot_ = NULL;
    kDetClu_ = 0;
    likeFirst_ = false;
    flag_ = true;
  }

//========================================================================
  CsRCPhoton::CsRCPhoton( int kPhot, double the, double phi,
//----------------------------------------------------------
			  double phiA, double theB, double theM, double PH,
			  CsRCCluster* ptclu, CsRCPhotonDet* ptdet, 
			  CsRCPartPhotons* pPartPhot, int kdeclu ) {
    kPhot_ = kPhot;
    the0_ = the;
    the_ = the;
    theNorm_ = the;
    phi_ = phi;
    phiA_ = phiA;
    theB_ = theB;
    theM_ = theM;
    PH_ = PH;
    ptToClu_ = ptclu;
    ptToDet_ = ptdet;
    kDetClu_ = kdeclu;
    pPartPhot_ = pPartPhot;
    likeFirst_ = false;
    flag_ = true;

    CFRefInd_ = CsRCRecConst::Instance()->CFRefInd();
    for( int k=0; k<5; k++ ) {
      thetaIpo_[k] = pPartPhot->thetaIpo( k*3+2 );
    }

    if( isPMT() ) {
      theNorm_ = thetaVStoUV( the );
      //std::cout << the_ << "  " << theNorm_  << "  " << the_*facVStoUV()
      //	  << std::endl;
      CFRefInd_ = CsRCRecConst::Instance()->CFRefIndVS();
      for( int k=0; k<5; k++ ) {
	thetaIpo_[k] = pPartPhot->thetaIpoVS( k*3+2 );
      }
    }
  }

//========================================================================
  CsRCPhoton::CsRCPhoton( const CsRCPhoton &phot ) {
//--------------------------------------------------
    //cout << "RICHONE : CsRCPhoton CopyConstructor" << endl;
    kPhot_ = phot.kPhot_;
    the0_ = phot.the0_;
    the_ = phot.the_;
    theNorm_ = phot.theNorm_;
    phi_ = phot.phi_;
    phiA_ = phot.phiA_;
    theB_ = phot.theB_;
    theM_ = phot.theM_;
    PH_ = phot.PH_;
    ptToClu_ = phot.ptToClu_;
    ptToDet_ = phot.ptToDet_;
    pPartPhot_ = phot.pPartPhot_;
    kDetClu_ = phot.kDetClu_;
    CFRefInd_ = phot.CFRefInd_;
    facVStoUV_ = phot.facVStoUV_;
    for( int k=0; k<5; k++) thetaIpo_[k] = phot.thetaIpo_[k];
    likeFirst_ = phot.likeFirst_;
    flag_ = phot.flag_;
  }

//========================================================================
  void CsRCPhoton::print() const {
//--------------------------------
    cout << endl;
    cout << "   Photon  " << kPhot_
         << " :  theta  " << the_
         << ",  thetaNorm  " << theNorm_
         << ",  phi  " << phi_ 
         << ",  phiA " << phiA_ 
         << ",  theB " << theB_ 
         << ",  theM " << theM_ 
         << ",  PH " << PH_ << endl;
    cout << "                cluster  " << ptToClu_->kClu()
	 << ",  cathode   " << ptToClu_->ic()
	 << ",  detector  " << ptToDet_->kDet()
	 << ",  flag  " << flag_ << endl;
    list<CsRCPad*> lPads = ptToClu_->lPads();
    list<CsRCPad*>::iterator ia;
    cout << "                padPH  ";
    for( ia=lPads.begin(); ia!=lPads.end(); ia++ ) cout << (*ia)->PH() << " ";
    cout << endl;
  }

//========================================================================
  CsRCPhoton::~CsRCPhoton() {}
//----------------------------


//========================================================================
  double CsRCPhoton::sigmaPhoPid( const CsRCPartPhotons* papho ) const {
//----------------------------------------------------------------------

//--- interface function :
//    --------------------
//    Paolo  -  July  2001
//    Rev.      December  2002
//    Rev.      June  2006


      CsRCHistos& hist = CsRCHistos::Ref();
      double xh, yh, wh;

      double sigmaPhot = 1000000.;
      float sigma = 0.;

      if( isAPV() ) sigma = sigmaPhot6APV();
//@@                -----------------------
      else if( isPMT() ) sigma = sigmaPhot6PMT();
//@@                     -----------------------
      else sigma = sigmaPhot6();
//@@       --------------------

      //float corrBeta = papho->getCorrBeta( papho->sigmaPhot5( beta ) );
      double mom = papho->pPart()->mom();
      float corrMom = papho->getCorrMom( papho->sigmaPhot7( mom ) );
//@@---------------------------------------------------------------
      sigma *= corrMom;
//@@----------------------
      sigma *= papho->getCorrFactor();
//@@---------------------------------

      if( sigma > 0. ) sigmaPhot = sigma;

      //std::cout << "sigmaPhoPid  " << sigmaPhot << std::endl;

//--- monitor error values
      xh = phi_;
      yh = sigmaPhot;
      if( hist.hRC1702 ) hist.hRC1702->Fill( xh, yh );
//hh                     ----------------------------
      xh = mom;
      yh = sigmaPhot;
      if( hist.hRC1703 ) hist.hRC1703->Fill( xh, yh );
//hh                     ----------------------------

      return  sigmaPhot;

  }


//========================================================================
  double CsRCPhoton::sigmaPhot6() const {
//---------------------------------------


//  February  2000,  rev. April  2001


//- sigma single-photon from data (or MC) :
//  ---------------------------------------

    static float par[10];
    static double gr;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      CsOpt* opt = CsOpt::Instance();
      CsRCRecConst *cons = CsRCRecConst::Instance();

      gr = 0.0174533;

//--- defaults
//    ( from Dima's files ),
//    ( Paolo file /scra/paolo/rc5-78k-spho.hist )  22/11/99
//    --------------------------------------------
      //par[0] =   0.805;
      //par[1] =   0.0014;
      //par[2] = - 0.0689;
      //par[3] =   0.0580;
      //par[4] =   0.0533;

//    from Vadim's 6*3k ref files/lepto_full ( co-18k-15t.hist )  4/01
//    ----------------------------------------------------------
      //par[0] =   0.89;
      //par[1] =   0.0082;
      //par[2] = - 0.0516;
      //par[3] =   0.0548;
      //par[4] =   0.0537;

//    from data 2002 ( run 22642 ) 02/12/10
//    -------------------------------------
      par[0] =   1.16;
      par[1] = - 0.026;
      par[2] = - 0.0137;
      par[3] =   0.0143;
      par[4] =   0.0405;

//--- from rich1.options :
//    --------------------
      vector<float> vPar;
      bool boo = opt->CsOpt::getOpt( "RICHONE", "sigmaPhoPid", vPar );
      if( boo ) {
	for( unsigned int k=0; k<vPar.size(); k++ ) par[k] = vPar[k];

	if( cons->printConsts() ) {
	  cout << " RICHONE, CsRCPhoton::sigmaPhot6 :   ";
	  for( unsigned int k=0; k<vPar.size(); k++ ) cout << par[k] <<"  ";
	    cout << endl;
	}
      } else {
	//cout << " RICHONE, CsRCPhoton::sigmaPhot6 :   ";
	//cout << "NO parameters read, default used!" << endl;
        string mess = "RICHONE, CsRCPhoton::sigmaPhot6 : ";
	string err = "NO parameters read, default used!";
	mess.append( err );
	CsErrLog::Instance()->mes( elError, mess );
      }

    }
 

//- MC : from sigphopid.kumac/sigmapho6.f - histo.s : 3521-4 (+off)
//  ---------------------------------------------------------------
//- Data : from sigphoton-data.kumac  &  histo 3521-4 (+off)
//  --------------------------------------------------------

//- phi around particle, from normal to 'particle plane' (deg)

    float phiw = phi_ * gr;                    //   to rad

    float sigmaPhot = par[0] + 
                      par[1] * cos( phiw ) +
                      par[2] * cos( 2.*phiw ) +
                      par[3] * cos( 3.*phiw ) +
	              par[4] * cos( 4.*phiw );
    //cout << "sigmaPhot6  " << sigmaPhot << endl;

    return sigmaPhot;

  }


//========================================================================
  double CsRCPhoton::sigmaPhot6APV() const {
//------------------------------------------

//- Paolo  -  June  2006

//- sigma single-photon from data (or MC) for APV cathodes :
//  --------------------------------------------------------

    static float par[10];
    static double gr;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      CsOpt* opt = CsOpt::Instance();
      CsRCRecConst *cons = CsRCRecConst::Instance();

      gr = 0.0174533;

//    from data 2006 ( run 50704 ) 060908
//    -----------------------------------
      par[0] =   1.46;
      par[1] = - 0.0846;
      par[2] =   0.0568;
      par[3] =  -0.0076;
      par[4] = - 0.0240;

//--- from rich1.options :
//    --------------------
      vector<float> vPar;
      bool boo = opt->CsOpt::getOpt( "RICHONE", "sigmaPhoPidAPV", vPar );
      if( boo ) {
	for( unsigned int k=0; k<vPar.size(); k++ ) par[k] = vPar[k];

	if( cons->printConsts() ) {
	  std::cout << " RICHONE, CsRCPhoton::sigmaPhot6APV :   ";
	  std::cout << setprecision ( 3 );
	  for( unsigned int k=0; k<vPar.size(); k++ ) 
	    std::cout << par[k] <<"  ";
	    std::cout << endl;
	    std::cout << setprecision ( 2 );
	}
      } else {
	//std::cout << " RICHONE, CsRCPhoton::sigmaPhot6APV :   ";
	//std::cout << "NO parameters read, default used!" << endl;
        string mess = "RICHONE, CsRCPhoton::sigmaPhot6APV : ";
	string err = "NO parameters read, default used!";
	mess.append( err );
	CsErrLog::Instance()->mes( elError, mess );
      }
    }
 
//- Data : from sigphophi.kumac  &  histo 3521-4 (+2000)
//  ----------------------------------------------------

//- phi around particle, from normal to 'particle plane' (deg)

    float phiw = phi_ * gr;                    //   to rad

    float sigmaPhot = par[0] + 
                      par[1] * cos( phiw ) +
                      par[2] * cos( 2.*phiw ) +
                      par[3] * cos( 3.*phiw ) +
	              par[4] * cos( 4.*phiw );

//- Pro-memoria
    //^sigmaPhot = sigmaPhot6();

    return  sigmaPhot;

  }


//========================================================================
  double CsRCPhoton::sigmaPhot6PMT() const {
//------------------------------------------

//- Paolo  -  June  2006

//- sigma single-photon from data (or MC) for MAPMT cathodes :
//  ----------------------------------------------------------

    static float par[10];
    static double gr;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      CsOpt* opt = CsOpt::Instance();
      CsRCRecConst *cons = CsRCRecConst::Instance();

      gr = 0.0174533;

//    from data 2006 ( run 50704 ) 060908
//    -----------------------------------
      par[0] =   1.62;
      par[1] =   0.151;
      par[2] = - 0.0713;
      par[3] = - 0.0365;
      par[4] = - 0.0714;

//--- from rich1.options :
//    --------------------
      vector<float> vPar;
      bool boo = opt->CsOpt::getOpt( "RICHONE", "sigmaPhoPidPMT", vPar );
      if( boo ) {
	for( unsigned int k=0; k<vPar.size(); k++ ) par[k] = vPar[k];

	if( cons->printConsts() ) {
	  std::cout << " RICHONE, CsRCPhoton::sigmaPhot6PMT :   ";
	  std::cout << setprecision ( 3 );
	  for( unsigned int k=0; k<vPar.size(); k++ ) 
	    std::cout << par[k] <<"  ";
	    std::cout << endl;
	  std::cout << setprecision ( 2 );
	}
      } else {
	//std::cout << " RICHONE, CsRCPhoton::sigmaPhot6PMT :   ";
	//std::cout << "NO parameters read, default used!" << endl;
        string mess = "RICHONE, CsRCPhoton::sigmaPhot6PMT : ";
	string err = "NO parameters read, default used!";
	mess.append( err );
	CsErrLog::Instance()->mes( elError, mess );
      }
    }
 
//- Data : from sigphophi.kumac  &  histo 3521-4 (+2000)
//  ----------------------------------------------------

//- phi around particle, from normal to 'particle plane' (deg)

    float phiw = phi_ * gr;                    //   to rad

    float sigmaPhot = par[0] + 
                      par[1] * cos( phiw ) +
                      par[2] * cos( 2.*phiw ) +
                      par[3] * cos( 3.*phiw ) +
	              par[4] * cos( 4.*phiw );

//- Provisional
    //^sigmaPhot = sigmaPhot6();

    return  sigmaPhot;
  }


//========================================================================
  double CsRCPhoton::facVStoUV() const {
//--------------------------------------

//- Paolo
//- August  2005

    CsRCRecConst *cons = CsRCRecConst::Instance();

    double fact = 0.;
    if( the_ > 0.) {
      double fact0 = ( cons->CFRefIndUV() - cons->CFRefIndVS() ) /
                       cons->CFRefIndVS();
      fact0 /= (the_ * the_);
      fact0 *= 1000000.;

      fact = 1. + fact0 * ( 1. - 0.5* fact0 );

      //std::cout << setprecision( 4 );
      //std::cout << " facVStoUV   " << the_ << "  " << fact0 
      //  	  << "  " << 0.5*fact0*fact0 << std::endl;
    }

    return  fact;

  }


//========================================================================
  double CsRCPhoton::thetaVStoUV( const double the ) const {
//----------------------------------------------------------

//- Paolo
//- August  2005

//- angles in mrad

/*
    CsRCRecConst *cons = CsRCRecConst::Instance();

    double beta1 = cons->CFRefIndVS() * cos( the/1000. );
    double cos = 1.;
    if( cons->CFRefIndUV() > 0.) cos = beta1 / cons->CFRefIndUV();
    if( cos > 1.) cos = 1.;
    double theta = 1000. * acos( cos );
    //std::cout << "thetaVStoUV " << the << "  " << theta << std::endl;

    return  theta;
*/

//- January  2007

    return  pPartPhot_->thetaVStoUV( the );

  }


//========================================================================
  double CsRCPhoton::thetaUVtoVS( const double the ) const {
//----------------------------------------------------------

//- Paolo
//- August  2005

//- angles in mrad

/*
    CsRCRecConst *cons = CsRCRecConst::Instance();

    double beta1 = cons->CFRefIndUV() * cos( the/1000. );
    double cos = 1.;
    if( cons->CFRefIndVS() > 0.) cos = beta1 / cons->CFRefIndVS();
    if( cos > 1.) cos = 1.;
    double theta = 1000. * acos( cos );
    //std::cout << "thetaUVtoVS " << the << "  " << theta << std::endl;

    return  theta;
*/

//- January  2007

    return  pPartPhot_->thetaUVtoVS( the );

  }


//========================================================================
  double CsRCPhoton::getThetaIpo( const double theIpo ) const {
//-------------------------------------------------------------

//- Paolo
//- August  2005

    double theta = 0.;
    if( theIpo < 0. ) {
      int kHypo = ( int( abs( theIpo ) )-2 )/3 + 1;
      theta = thetaIpo( kHypo );
    } else {
      theta = theIpo;
      if( isPMT() ) theta = thetaUVtoVS( theIpo );
    }

    return  theta;
  }


//========================================================================
  void CsRCPhoton::flagAllPhotons( const CsRCPhoton* phot ) const {
//-----------------------------------------------------------------

//- Paolo
//- June  2007

    CsRCEventPartPhotons* paphos = CsRCEventPartPhotons::Instance();
    list<CsRCPartPhotons*> lPaPhotons = paphos->lPartPhotons();
    list<CsRCPartPhotons*>::iterator ipf;

    CsRCPhoton* phow = const_cast<CsRCPhoton*>(phot);
    phow->setFlag( false );
    CsRCCluster* clu = phot->ptToClu();

    for( ipf=lPaPhotons.begin(); ipf!=lPaPhotons.end(); ipf++ ) {
      if( !(*ipf)->flag() ) continue;
      if( (*ipf) == phot->pPartPhot() ) continue;
      list<CsRCPhoton*> lPhotons = (*ipf)->lPhotons();
      list<CsRCPhoton*>::iterator ih;
      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	if( (*ih)->ptToClu() == clu ) (*ih)->setFlag( false );
      }
    }

    return;
 }



//=========================================================================== 
  double CsRCPhoton::getPMTOptCorr() {
//------------------------------------

//- Paolo  -  December 2009

    CsRCDetectors *dets = CsRCDetectors::Instance();
    CsRCMirrors *mirr = CsRCMirrors::Instance();

//- As from 14/12/2009

    static double parCorrTheAve[49][6];

    parCorrTheAve[1][0] = 51.5818;
    parCorrTheAve[1][1] = 0.228304;
    parCorrTheAve[1][2] = 206.016;
    parCorrTheAve[1][3] = 51.3995;
    parCorrTheAve[1][4] = 1.0937;
    parCorrTheAve[1][5] = 230.674;

    parCorrTheAve[2][0] = 51.3608;
    parCorrTheAve[2][1] = 0.240205;
    parCorrTheAve[2][2] = 347.319;
    parCorrTheAve[2][3] = 51.4769;
    parCorrTheAve[2][4] = 0.762537;
    parCorrTheAve[2][5] = 254.133;

    parCorrTheAve[3][0] = 51.4884;
    parCorrTheAve[3][1] = 0.101601;
    parCorrTheAve[3][2] = 36.2885;
    parCorrTheAve[3][3] = 51.5525;
    parCorrTheAve[3][4] = 0.761132;
    parCorrTheAve[3][5] = 279.313;

    parCorrTheAve[4][0] = 51.3512;
    parCorrTheAve[4][1] = 0.208743;
    parCorrTheAve[4][2] = 93.0916;
    parCorrTheAve[4][3] = 51.5465;
    parCorrTheAve[4][4] = 0.902237;
    parCorrTheAve[4][5] = 306.337;

    parCorrTheAve[5][0] = 51.1471;
    parCorrTheAve[5][1] = 0.424416;
    parCorrTheAve[5][2] = 119.851;
    parCorrTheAve[5][3] = 51.5244;
    parCorrTheAve[5][4] = 1.64099;
    parCorrTheAve[5][5] = 318.808;

    parCorrTheAve[6][0] = 51.3689;
    parCorrTheAve[6][1] = 0.183587;
    parCorrTheAve[6][2] = 321.004;
    parCorrTheAve[6][3] = 51.4317;
    parCorrTheAve[6][4] = 0.710331;
    parCorrTheAve[6][5] = 204.493;

    parCorrTheAve[7][0] = 51.4679;
    parCorrTheAve[7][1] = 0.13017;
    parCorrTheAve[7][2] = 297.427;
    parCorrTheAve[7][3] = 51.5592;
    parCorrTheAve[7][4] = 0.676584;
    parCorrTheAve[7][5] = 232.357;

    parCorrTheAve[8][0] = 51.5022;
    parCorrTheAve[8][1] = 0.0980646;
    parCorrTheAve[8][2] = 321.705;
    parCorrTheAve[8][3] = 51.6034;
    parCorrTheAve[8][4] = 0.548216;
    parCorrTheAve[8][5] = 285.069;

    parCorrTheAve[9][0] = 51.4794;
    parCorrTheAve[9][1] = 0.118681;
    parCorrTheAve[9][2] = 356.123;
    parCorrTheAve[9][3] = 51.6032;
    parCorrTheAve[9][4] = 0.833881;
    parCorrTheAve[9][5] = 330.423;

    parCorrTheAve[10][0] = 51.448;
    parCorrTheAve[10][1] = 0.110625;
    parCorrTheAve[10][2] = 352.83;
    parCorrTheAve[10][3] = 51.5965;
    parCorrTheAve[10][4] = 1.05066;
    parCorrTheAve[10][5] = 344.283;

    parCorrTheAve[11][0] = 51.2593;
    parCorrTheAve[11][1] = 0.293669;
    parCorrTheAve[11][2] = 355.849;
    parCorrTheAve[11][3] = 51.5082;
    parCorrTheAve[11][4] = 0.790778;
    parCorrTheAve[11][5] = 173.695;

    parCorrTheAve[12][0] = 51.4839;
    parCorrTheAve[12][1] = 0.0760606;
    parCorrTheAve[12][2] = 286.274;
    parCorrTheAve[12][3] = 51.6331;
    parCorrTheAve[12][4] = 0.637741;
    parCorrTheAve[12][5] = 171.267;

    parCorrTheAve[13][0] = 51.4945;
    parCorrTheAve[13][1] = 0.0942979;
    parCorrTheAve[13][2] = 326.652;
    parCorrTheAve[13][3] = 51.6286;
    parCorrTheAve[13][4] = 0.251859;
    parCorrTheAve[13][5] = 43.9548;

    parCorrTheAve[14][0] = 51.5077;
    parCorrTheAve[14][1] = 0.132556;
    parCorrTheAve[14][2] = 335.564;
    parCorrTheAve[14][3] = 51.6276;
    parCorrTheAve[14][4] = 0.774712;
    parCorrTheAve[14][5] = 7.27865;

    parCorrTheAve[15][0] = 51.4392;
    parCorrTheAve[15][1] = 0.125418;
    parCorrTheAve[15][2] = 328.082;
    parCorrTheAve[15][3] = 51.576;
    parCorrTheAve[15][4] = 1.02161;
    parCorrTheAve[15][5] = 7.79996;

    parCorrTheAve[16][0] = 51.4039;
    parCorrTheAve[16][1] = 0.221805;
    parCorrTheAve[16][2] = 14.8365;
    parCorrTheAve[16][3] = 51.5689;
    parCorrTheAve[16][4] = 1.07227;
    parCorrTheAve[16][5] = 145.787;

    parCorrTheAve[17][0] = 51.6694;
    parCorrTheAve[17][1] = 0.158506;
    parCorrTheAve[17][2] = 101.496;
    parCorrTheAve[17][3] = 51.6886;
    parCorrTheAve[17][4] = 0.948467;
    parCorrTheAve[17][5] = 129.212;

    parCorrTheAve[18][0] = 51.6261;
    parCorrTheAve[18][1] = 0.144237;
    parCorrTheAve[18][2] = 59.4265;
    parCorrTheAve[18][3] = 51.6708;
    parCorrTheAve[18][4] = 0.858494;
    parCorrTheAve[18][5] = 87.5799;

    parCorrTheAve[19][0] = 51.6249;
    parCorrTheAve[19][1] = 0.191744;
    parCorrTheAve[19][2] = 40.8716;
    parCorrTheAve[19][3] = 51.6724;
    parCorrTheAve[19][4] = 1.01104;
    parCorrTheAve[19][5] = 50.3981;

    parCorrTheAve[20][0] = 51.4592;
    parCorrTheAve[20][1] = 0.139799;
    parCorrTheAve[20][2] = 121.43;
    parCorrTheAve[20][3] = 51.5656;
    parCorrTheAve[20][4] = 1.16214;
    parCorrTheAve[20][5] = 34.6579;

    parCorrTheAve[21][0] = 51.406;
    parCorrTheAve[21][1] = 0.538933;
    parCorrTheAve[21][2] = 103.828;
    parCorrTheAve[21][3] = 51.731;
    parCorrTheAve[21][4] = 2.63741;
    parCorrTheAve[21][5] = 130.436;

    parCorrTheAve[22][0] = 51.6713;
    parCorrTheAve[22][1] = 0.251763;
    parCorrTheAve[22][2] = 78.2018;
    parCorrTheAve[22][3] = 51.6924;
    parCorrTheAve[22][4] = 1.38444;
    parCorrTheAve[22][5] = 111.836;

    parCorrTheAve[23][0] = 51.7079;
    parCorrTheAve[23][1] = 0.265453;
    parCorrTheAve[23][2] = 67.5935;
    parCorrTheAve[23][3] = 51.7619;
    parCorrTheAve[23][4] = 1.32805;
    parCorrTheAve[23][5] = 90.2442;

    parCorrTheAve[24][0] = 51.6714;
    parCorrTheAve[24][1] = 0.3179;
    parCorrTheAve[24][2] = 48.7364;
    parCorrTheAve[24][3] = 51.6567;
    parCorrTheAve[24][4] = 1.46413;
    parCorrTheAve[24][5] = 65.7089;

    parCorrTheAve[25][0] = 51.5067;
    parCorrTheAve[25][1] = 0.688652;
    parCorrTheAve[25][2] = 244.82;
    parCorrTheAve[25][3] = 50.9209;
    parCorrTheAve[25][4] = 2.53101;
    parCorrTheAve[25][5] = 44.4146;


    const double DegRad = 0.0174533;
    const int ncTxy = 7;
    const double dDirx = 0.05;
    const double xDirLmn = -0.06 - dDirx/2.;
    const double xDirLmx =  0.24 + dDirx/2.;
    const double dDiry = 0.05;
    const double yDirLmn = -0.40 - dDiry/2.;
    const double yDirLmx = -0.10 + dDiry/2.;

    double dTheta = 0.;

    CsRCCluster* clu = ptToClu();
    double thePhoIn = theNorm();
    double phiPhoIn = phiA();

//- Photon Cherenkov angles
    int iCat = clu->ic();
    int kDetClu = dets->cathodePos( iCat );
    Hep3Vector vDcPartWw = pPartPhot()->vDcPartW()[kDetClu];
    double ll = sin( thePhoIn/1000. ) * cos( phiPhoIn*DegRad );
    double mm = sin( thePhoIn/1000. ) * sin( phiPhoIn*DegRad );
    double nn = cos( thePhoIn/1000. );
    Hep3Vector vDcPho( ll, mm, nn );
    vDcPho = vDcPho.unit();
    double thePho = 0.;
    double phiPho = 0.;
    Hep3Vector vDcPhoEmW = 
      CsRCUty::Instance()->rotfbcc( -1., vDcPartWw, vDcPho, thePho, phiPho );
//  ------------------------------------------------------------------------
//- Photon impact on mirror
    double RR = mirr->RRv( kDetClu );
    Hep3Vector vPoPhotWw = pPartPhot()->vPoPhotW()[kDetClu];
    Hep3Vector vPoC( 0., 0., 0. );
    Hep3Vector vPoPhoMir = mirr->vImpMir( vPoPhotWw, vDcPhoEmW, vPoC, RR );
//                         -----------------------------------------------
//- Normal to mirror at Photon impact
    Hep3Vector vDcNoPhoMir = (1./RR) * vPoPhoMir;
//- Photon reflected direction
    double cosPhoMir = vDcNoPhoMir * vDcPhoEmW;
    Hep3Vector vDcPhoRefl = 2.*cosPhoMir * vDcNoPhoMir - vDcPhoEmW;
//- Photon direction at detector (no QZW correction)
    Hep3Vector vDcPhoDetW = vDcPhoRefl.unit();

    int kQua = -1;
    Hep3Vector vDcPhoDetWw = vDcPhoDetW;
//- Mirroring for cathodes different from #03 (TopJura)
//  Optical Correction Table is for TopJura ONLY
//- TopJura
    if( iCat == dets->nCatPMT()[0] ) {
      kQua = 0;
    }
//- TopSaleve
    else if( iCat == dets->nCatPMT()[1] ) {
      vDcPhoDetWw.setX( -vDcPhoDetW.x() );
      kQua = 1;
    }
//- DownJura
    else if( iCat == dets->nCatPMT()[2] ) {
      vDcPhoDetWw.setY( -vDcPhoDetW.y() );
      kQua = 2;
    }
//- DownSaleve
    else if( iCat == dets->nCatPMT()[3] ) {
      vDcPhoDetWw.setX( -vDcPhoDetW.x() );
      vDcPhoDetWw.setY( -vDcPhoDetW.y() );
      kQua = 3;
    }
    double tgx = vDcPhoDetWw.x()/vDcPhoDetWw.z();
    double tgy = vDcPhoDetWw.y()/vDcPhoDetWw.z();

    int ktx, kty;
    ktx = -1;
    kty = -1;
    for( int ky=0; ky<ncTxy; ky++ ) {
      double yDirn = yDirLmn + ky * dDiry;
      double yDirx = yDirn + dDiry;
      if( tgy >= yDirn  &&  tgy < yDirx ) {
	kty = ky;
	ktx = -1;
	for( int kx=0; kx<ncTxy; kx++ ) {
	  double xDirn = xDirLmn + kx * dDirx;
	  double xDirx = xDirn + dDirx;
	  if( tgx >= xDirn  &&  tgx < xDirx ) {
	    ktx = kx;
	    break;
	  }
	}
	break;
      }
    }
    if( tgx < xDirLmn ) ktx = 0;
    if( tgx > xDirLmx ) ktx = ncTxy + 1;
    if( tgy < yDirLmn ) kty = 0;
    if( tgy > yDirLmx ) kty = ncTxy + 1;
    if( ktx <= 1 ) ktx = 1;
    if( ktx >= 5 ) ktx = 5;
    if( kty <= 1 ) kty = 1;
    if( kty >= 5 ) kty = 5;

    int kxy = (kty-1) * 5 + ktx;
    double r0t = parCorrTheAve[kxy][0];
    double rrt = parCorrTheAve[kxy][1];
    double f0t = parCorrTheAve[kxy][2];
    double r0p = parCorrTheAve[kxy][3];
    double rrp = parCorrTheAve[kxy][4];
    double f0p = parCorrTheAve[kxy][5];

    double thePhoCo = thePhoIn - 
      ( r0p + rrp * cos( (phiPhoIn - f0p)*DegRad ) - r0t );
    double phiPhoCo = phiPhoIn;
    dTheta = thePhoCo - thePhoIn;

    return  dTheta;
  }
