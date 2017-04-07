/*!
   \File    CsRCLikeRing02.cc
   \----------------------------
   \brief   CsRCLikeRing02 class implementation.
   \author  Paolo Schiavon
   \version 1.0
   \date    9 December 2002
*/

  #include <iostream>
  #include <ostream>
  #include <cstdio>

  #include "CsOpt.h"
  #include "CsHist.h"
  #include "CsHistograms.h"

//------------------------------
  #include "CsRCLikeRing.h"
  #include "CsRCLikeRing02.h"

  #include "CsRCParticle.h"
  #include "CsRCCluster.h"
  #include "CsRCEventPartPhotons.h"
  #include "CsRCPartPhotons.h"
  #include "CsRCPhoton.h"

  #include "CsRCRing.h"

  #include "CsRCHistos.h"

  #include "CsRCRecConst.h"
  #include "CsRCExeKeys.h"
//------------------------------

  using namespace std;

//===========================================================================
  CsRCLikeRing02::CsRCLikeRing02( const CsRCRing* ring ): CsRCLikeRing() {
//------------------------------------------------------------------------

    pPartPhot_ = ring->pPartPhot();
    theUpL_ = ring->theUpL();
    theLoL_ = ring->theLoL();
    //cout << pPartPhot_ << "  " << theUpL_ << "  " << theLoL_ << endl;
  }

//===========================================================================
  CsRCLikeRing02::~CsRCLikeRing02() {}
//------------------------------------


  extern "C" { float erf_( const float& ); }

//===========================================================================
  double CsRCLikeRing02::normSignal( const double theIpo ) {
//----------------------------------------------------------


//- compute Signal normalization for Likelihood
//  -------------------------------------------
//- Paolo  -  December 2002


//- average sigma-photon :
    //double momPart = pPartPhot_->pPart()->mom();
    //double sigPhoRec = pPartPhot_->sigmaPhoRec( momPart );
    CsRCPartPhotons* pt = const_cast<CsRCPartPhotons*>(pPartPhot_);
    double momPart = pt->pPart()->mom();
    double sigPhoRec = pt->sigmaPhoRec( momPart );
//@@                   --------------------------

    double normS = 0.5* ( erf_( (theUpL_-theIpo)/(sqrt( 2.)*sigPhoRec) )
			- erf_( (theLoL_-theIpo)/(sqrt( 2.)*sigPhoRec) ) );
    if( theIpo <= 0. ) normS = 0.;
    //cout << theUpL_ << "  " << theLoL_ << "  " << normS << "  "
    //     << theIpo << "  " << sigPhoRec << "  " << sqrt( 2.) << endl;
    return  normS;

  }


//===========================================================================
  double CsRCLikeRing02::normBackgr( const double theReco,
//--------------------------------------------------------
				     const double theIpo ) {


//- compute Background normalization for Likelihood
//  -----------------------------------------------
//- Paolo  -  December 2002


    CsRCRecConst * cons = CsRCRecConst::Instance();

    static float dBindTheta = 4.;
    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      CsRCHistos& hist = CsRCHistos::Ref();
//--- reference to histo 3561
      float XMin =  0.;
      float XMax = 70.;
      size_t nBin = 280;
      if( hist.hRC3561 ) nBin = hist.hRC3561->GetDim(1).GetNBins();
      if( hist.hRC3561 ) XMin = hist.hRC3561->GetDim(1).GetMin();
      if( hist.hRC3561 ) XMax = hist.hRC3561->GetDim(1).GetMax();
      dBindTheta = float( nBin ) / (XMax - XMin);
      //cout << "dBindTheta  " << dBindTheta << endl;
    }

//- average sigma-photon :
    //double momPart = pPartPhot_->pPart()->mom();
    //double sigPhoRec = pPartPhot_->sigmaPhoRec( momPart );
    CsRCPartPhotons* pt = const_cast<CsRCPartPhotons*>(pPartPhot_);
    double momPart = pt->pPart()->mom();
    double sigPhoRec = pt->sigmaPhoRec( momPart );
//                     --------------------------

//- from back4-fit, 2002, run 20330
//  -------------------------------
//  bkPar0 : (p0 + p1*theIpo),  bkPar1 : (p2 + p3*theIpo)
//  bkgr = bkPar0*the + bkPar1*the**2
//  -----------------------------------------------------

    double delta = theUpL_ - theLoL_;   //   useful wind.
//@@--------------------------------

    vector<float> par = CsRCEventPartPhotons::Instance()->vPara();
    //for( int k=0; k<4; k++ ) cout << par[k] << "  "; cout << endl;

    double bkPar0 = par[0] + par[1] * theIpo;
    double bkPar1 = par[2] + par[3] * theIpo;
    double normB = dBindTheta * delta *
      ( bkPar0* theReco + bkPar1*theReco*theReco + bkPar1* delta*delta/12.);

    return  normB;

  }


//===========================================================================
  double CsRCLikeRing02::likeSignal( const double thePhot,
//--------------------------------------------------------
				     const double theIpo,
				     const double sigPhot ) {


//- compute Signal value for Likelihood
//  -----------------------------------
//- Paolo  -  December 2002


    CsRCRecConst* cons = CsRCRecConst::Instance();
    static double TwoPI = cons->TwoPI();
    static double sq2pi = sqrt( TwoPI );

    double qSw = (thePhot - theIpo) / sigPhot;
    double signal = exp( - 0.5* qSw*qSw ) / (sigPhot*sq2pi);

    if( theIpo <= 0. ) signal = 0.;

    return  signal;

  }


//===========================================================================
  //double CsRCLikeRing02::likeBackgr( const double thePhot,
  double CsRCLikeRing02::likeBackgr( CsRCPhoton* pPhot,
//----------------------------------------------------------
				     const double theIpo ) {


//- compute Background value for Likelihood
//  ---------------------------------------
//- Paolo  -  December 2002

    double thePhot = pPhot->the();

    vector<float> par = CsRCEventPartPhotons::Instance()->vPara();
    //for( int k=0; k<4; k++ ) cout << par[k] << "  "; cout << endl;

    double bkPar0 = par[0] + par[1] * theIpo;
    double bkPar1 = par[2] + par[3] * theIpo;
    double backgr = bkPar0 * thePhot + bkPar1 * thePhot*thePhot;

    CsRCHistos& hist = CsRCHistos::Ref();
    double xh, yh, wh;
    if( backgr > 0. ) {
      xh = pPhot->ptToClu()->xc();
      yh = pPhot->ptToClu()->yc();
      wh = backgr;
      if( hist.hRC1530 ) hist.hRC1530->Fill( xh, yh, wh );
//    ---------------------------------------------------
      if( hist.hRC1540 ) hist.hRC1540->Fill( xh, yh );
//    -----------------------------------------------
    }

    return  backgr;

  }


//===========================================================================
  double CsRCLikeRing02::getRingBackground( const double theReco ) {
//------------------------------------------------------------------


//-  Paolo - September 2002
//    rev.  December 2002


    CsRCRecConst * cons = CsRCRecConst::Instance();

    static float dBindTheta = 4.;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

//--- reference to histo 3561
      CsRCHistos& hist = CsRCHistos::Ref();
      float XMin =  0.;
      float XMax = 70.;
      size_t nBin = 280;
      if( hist.hRC3561 ) nBin = hist.hRC3561->GetDim(1).GetNBins();
      if( hist.hRC3561 ) XMin = hist.hRC3561->GetDim(1).GetMin();
      if( hist.hRC3561 ) XMax = hist.hRC3561->GetDim(1).GetMax();
      dBindTheta = float( nBin ) / (XMax - XMin);

    }

    //double momPart = pPartPhot_->pPart()->mom();
    //double sigPhoRec = pPartPhot_->sigmaPhoRec( momPart );
    CsRCPartPhotons* pt = const_cast<CsRCPartPhotons*>(pPartPhot_);
    double momPart = pt->pPart()->mom();
    double sigPhoRec = pt->sigmaPhoRec( momPart );
//                     --------------------------
    double delta = theUpL_ - theLoL_;   //   useful window
//@@--------------------------------
    double normSK = 0.5* ( erf_( (theUpL_-theReco)/(sqrt( 2.)*sigPhoRec) )
			 - erf_( (theLoL_-theReco)/(sqrt( 2.)*sigPhoRec) ) );

    vector<float> par = CsRCEventPartPhotons::Instance()->vPara();
    //for( int k=0; k<4; k++ ) cout << par[k] << "  "; cout << endl;

    double bkPar0K = par[0] + par[1] * theReco;
    double bkPar1K = par[2] + par[3] * theReco;
    double normBK = dBindTheta * delta *
      (bkPar0K* theReco + bkPar1K*theReco*theReco + bkPar1K* delta*delta/12.);

    double corrBack = 0.;
    if( fabs( normSK + normBK ) > 0. ) corrBack = normSK / (normSK + normBK);

    //if( corrBack == 0. ) {
      //cout << theUpL_ << "  " << theLoL_ << "  " << theReco << "  "
      //     << sigPhoRec << "  " << corrBack << endl;
    //}

    return  corrBack;

  }
