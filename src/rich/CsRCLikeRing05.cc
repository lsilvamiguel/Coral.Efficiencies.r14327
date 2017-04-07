/*!
   \File    CsRCLikeRing05.cc
   \----------------------------
   \brief   CsRCLikeRing05 class implementation.
   \author  Paolo Schiavon
   \version 1.0
   \date    May 2004
*/

  #include <iostream>
  #include <ostream>
  #include <cstdio>

  #include "CsOpt.h"
  #include "CsHist.h"
  #include "CsHistograms.h"

//------------------------------
  #include "CsRCLikeRing.h"
  #include "CsRCLikeRing05.h"

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
  CsRCLikeRing05::CsRCLikeRing05( const CsRCRing* ring ): CsRCLikeRing() {
//------------------------------------------------------------------------

    pPartPhot_ = ring->pPartPhot();
    theUpL_ = ring->theUpL();
    theLoL_ = ring->theLoL();

  }

//===========================================================================
  CsRCLikeRing05::~CsRCLikeRing05() {}
//------------------------------------


  extern "C" { float erf_( const float& ); }

//===========================================================================
  double CsRCLikeRing05::normSignal( const double theIpo ) {
//----------------------------------------------------------


//- compute Signal normalization for Likelihood
//  -------------------------------------------
//  dummy : assume likeRatio true
//- Paolo  -  May 2004


//- average sigma-photon :
    //CsRCPartPhotons* pt = const_cast<CsRCPartPhotons*>(pPartPhot_);
    //double momPart = pt->pPart()->mom();
    //double sigPhoRec = pt->sigmaPhoRec( momPart );
//@@                   --------------------------

    double normS = 0.;
    if( theIpo <= 0. ) normS = 0.;

    return  normS;

  }


//===========================================================================
  double CsRCLikeRing05::normBackgr( const double theReco,
//--------------------------------------------------------
				     const double theIpo ) {


//- compute Background normalization for Likelihood
//  -----------------------------------------------
//  dummy : assume likeRatio true
//- Paolo  -  May 2004


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
    //CsRCPartPhotons* pt = const_cast<CsRCPartPhotons*>(pPartPhot_);
    //double momPart = pt->pPart()->mom();
    //double sigPhoRec = pt->sigmaPhoRec( momPart );
//                     --------------------------

    //double delta = theUpL_ - theLoL_;   //   useful window
//@@--------------------------------

    double normB = 1.;

    return  normB;

  }


//===========================================================================
  double CsRCLikeRing05::likeSignal( const double thePhot,
//--------------------------------------------------------
				     const double theIpo,
				     const double sigPhot ) {


//- compute Signal value for Likelihood
//  -----------------------------------
//- Paolo  -  May 2004


    CsRCRecConst* cons = CsRCRecConst::Instance();
    static double TwoPI = cons->TwoPI();
    static double sq2pi = sqrt( TwoPI );

    double qSw = (thePhot - theIpo) / sigPhot;
    double signal = exp( - 0.5* qSw*qSw ) / (sigPhot*sq2pi);

    if( theIpo <= 0. ) signal = 0.;

    return  signal;

  }


//===========================================================================
  double CsRCLikeRing05::likeBackgr( CsRCPhoton* pPhot,
//----------------------------------------------------------
				     const double theIpo ) {


//- compute Background value for Likelihood
//  ---------------------------------------
//- Paolo  -  May 2004

//- TENTATIVE !


    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;
    }

//- from back-para.hist background map
//  ----------------------------------
    double thePhot = pPhot->the();

    double backgr = 0.;
    double backWgt = 0.;
    if( pPhot->ptToClu()->getBackWgt( backWgt ) ) {
//@@---------------------------------------------
      backgr = 0.095 * backWgt * thePhot;
    }

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
  double CsRCLikeRing05::getRingBackground( const double theReco ) {
//------------------------------------------------------------------


//-  Paolo - May 2004

//-  NOT YET IMPLEMENTED !


    CsRCRecConst * cons = CsRCRecConst::Instance();

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;
    }

    double corrBack = 0.;

    return  corrBack;

  }
