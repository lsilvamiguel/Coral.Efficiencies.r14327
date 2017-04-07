/*!
   \File    CsRCLikeRing03.cc
   \----------------------------
   \brief   CsRCLikeRing03 class implementation.
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
  #include "CsRCLikeRing03.h"

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
  CsRCLikeRing03::CsRCLikeRing03( const CsRCRing* ring ): CsRCLikeRing() {
//------------------------------------------------------------------------

    pPartPhot_ = ring->pPartPhot();
    theUpL_ = ring->theUpL();
    theLoL_ = ring->theLoL();

  }

//===========================================================================
  CsRCLikeRing03::~CsRCLikeRing03() {}
//------------------------------------


  extern "C" { float erf_( const float& ); }

//===========================================================================
  double CsRCLikeRing03::normSignal( const double theIpo ) {
//----------------------------------------------------------


//- compute Signal normalization for Likelihood
//  -------------------------------------------
//- Paolo  -  December 2002
//  modified for TEST  January 2004


//- average sigma-photon :
    //double momPart = pPartPhot_->pPart()->mom();
    //double sigPhoRec = pPartPhot_->sigmaPhoRec( momPart );
    CsRCPartPhotons* pt = const_cast<CsRCPartPhotons*>(pPartPhot_);
    double momPart = pt->pPart()->mom();
    double sigPhoRec = pt->sigmaPhoRec( momPart );
//@@                   --------------------------

    double delta = theUpL_ - theLoL_;
    double theIpoL = theIpo - delta/2.;
    double theIpoU = theIpo + delta/2.;
    double theWL = theLoL_;
    if( theLoL_ < theIpoL ) theWL = theIpoL;
    double theWU = theUpL_;
    if( theUpL_ > theIpoU ) theWU = theIpoU;
    //std::cout << theLoL_ << "  " << theUpL_ << "  "
    //  << theIpoL << "  " << theIpoU << "  "
    //  << theWL << "  " << theWU << std::endl;
    double normS = 0.;
    if( theWU > theWL ) {
      normS = 0.5* ( erf_( (theWU-theIpo)/(sqrt( 2.)*sigPhoRec) )
  	      - erf_( (theWL-theIpo)/(sqrt( 2.)*sigPhoRec) ) );
    }
    //double normS = 0.5* ( erf_( (theUpL_-theIpo)/(sqrt( 2.)*sigPhoRec) )
    //		     - erf_( (theLoL_-theIpo)/(sqrt( 2.)*sigPhoRec) ) );

    if( theIpo <= 0. ) normS = 0.;
    //double normS = 0.;
    //std::cout << "normS  " << normS << std::endl;
    return  normS;

  }


//===========================================================================
  double CsRCLikeRing03::normBackgr( const double theReco,
//--------------------------------------------------------
				     const double theIpo ) {


//- compute Background normalization for Likelihood
//  -----------------------------------------------
//- Paolo  -  December 2002
//  modified for TEST  January 2004


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
//  ------------------------------
//  <bkgr> = par[0] + par[1]*the + par[2]*the**2 + par[3]*the**3
//  ------------------------------------------------------------

    double delta = theUpL_ - theLoL_;   //   useful wind.
//@@--------------------------------

    vector<float> par = CsRCEventPartPhotons::Instance()->vPara();
    double bkPar0 = par[0] + par[1] * theIpo;
    double bkPar1 = par[2] + par[3] * theIpo;

    double theIpoL = theIpo - delta/2.;
    double theIpoU = theIpo + delta/2.;
    double theWL = theLoL_;
    if( theLoL_ < theIpoL ) theWL = theIpoL;
    double theWU = theUpL_;
    if( theUpL_ > theIpoU ) theWU = theIpoU;
    delta = theWU - theWL;
    if( delta < 0. ) delta = 0.;

    double normB = dBindTheta * delta *
      ( bkPar0* theReco + bkPar1*theReco*theReco + bkPar1* delta*delta/12.);

    if( normB <= 0. ) normB = 1.;
    //double normB = 1.;
    //std::cout << "normB  " << normB << std::endl;
    return  normB;

  }


//===========================================================================
  double CsRCLikeRing03::likeSignal( const double thePhot,
//--------------------------------------------------------
				     const double theIpo,
				     const double sigPhot ) {


//- compute Signal value for Likelihood
//  -----------------------------------
//- Paolo  -  December 2002
//  modified for TEST  January 2004


    CsRCRecConst* cons = CsRCRecConst::Instance();
    static double TwoPI = cons->TwoPI();
    static double sq2pi = sqrt( TwoPI );

    double delta = theUpL_ - theLoL_;
    double theIpoL = theIpo - delta/2.;
    double theIpoU = theIpo + delta/2.;
    double theWL = theLoL_;
    if( theLoL_ < theIpoL ) theWL = theIpoL;
    double theWU = theUpL_;
    if( theUpL_ > theIpoU ) theWU = theIpoU;

    double signal = 0.;
    if( thePhot >= theWL  &&  thePhot <= theWU ) {
      double qSw = (thePhot - theIpo) / sigPhot;
      signal = exp( - 0.5* qSw*qSw ) / (sigPhot*sq2pi);
    }

    if( theIpo <= 0. ) signal = 0.;
    //double signal = 0.;
    //std::cout << "signal  " << signal << std::endl;
    return  signal;

  }


//===========================================================================
  //double CsRCLikeRing03::likeBackgr( const double thePhot,
  double CsRCLikeRing03::likeBackgr( CsRCPhoton* pPhot,
//----------------------------------------------------------
				     const double theIpo ) {


//- compute Background value for Likelihood
//  ---------------------------------------
//- Paolo  -  December 2002
//  modified for TEST  January 2004

    double thePhot = pPhot->the();

    vector<float> par = CsRCEventPartPhotons::Instance()->vPara();
    double bkPar0 = par[0] + par[1] * theIpo;
    double bkPar1 = par[2] + par[3] * theIpo;

    double delta = theUpL_ - theLoL_;
    double theIpoL = theIpo - delta/2.;
    double theIpoU = theIpo + delta/2.;
    double theWL = theLoL_;
    if( theLoL_ < theIpoL ) theWL = theIpoL;
    double theWU = theUpL_;
    if( theUpL_ > theIpoU ) theWU = theIpoU;

    double backgr = 0.;
    if( thePhot >= theWL  &&  thePhot <= theWU ) {
      backgr = bkPar0 * thePhot + bkPar1 * thePhot*thePhot;
    }
    //double backgr = 0.;
    //std::cout << "backgr  " << backgr << std::endl;
    return  backgr;

  }


//===========================================================================
  double CsRCLikeRing03::getRingBackground( const double theReco ) {
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

    double bkPar0K = par[0] + par[1] * theReco;
    double bkPar1K = par[2] + par[3] * theReco;
    double normBK = dBindTheta * delta *
      (bkPar0K* theReco + bkPar1K*theReco*theReco + bkPar1K* delta*delta/12.);

    double corrBack = 0.;
    if( fabs( normSK + normBK ) > 0. ) corrBack = normSK / (normSK + normBK);

    return  corrBack;

  }
