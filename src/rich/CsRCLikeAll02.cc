/*!
   \File    CsRCLikeAll02.cc
   \----------------------------
   \brief   CsRCLikeAll02 class implementation.
   \author  Paolo Schiavon
   \version 1.0
   \date    23 December 2002
*/


  #include <iostream>
  #include <ostream>
  #include <cstdio>

  #include "CsOpt.h"
  #include "CsHist.h"
  #include "CsHistograms.h"

//------------------------------
  #include "CsRCLikeAll.h"
  #include "CsRCLikeAll02.h"

  #include "CsRCParticle.h"
  #include "CsRCCluster.h"
  #include "CsRCEventPartPhotons.h"
  #include "CsRCPartPhotons.h"
  #include "CsRCPhoton.h"

  #include "CsRCHistos.h"

  #include "CsRCRecConst.h"
  #include "CsRCExeKeys.h"
//------------------------------

  using namespace std;

//===========================================================================
  CsRCLikeAll02::CsRCLikeAll02( const CsRCPartPhotons* papho ):
//-------------------------------------------------------------
    CsRCLikeAll() {

    pPartPhot_ = papho;

  }

//===========================================================================
  CsRCLikeAll02::~CsRCLikeAll02() {}
//----------------------------------


  extern "C" {
    float prob_( const float&, const int& );
    float erf_( const float& );
  }


//===========================================================================
  double CsRCLikeAll02::normSignal( const double theIpo ) {
//---------------------------------------------------------


//- compute Signal normalization for Likelihood
//  -------------------------------------------
//- Paolo  -  December 2002


    CsRCRecConst* cons = CsRCRecConst::Instance();

    double normS = 1.;

    return  normS;

  }


//===========================================================================
  double CsRCLikeAll02::normBackgr( const double theIpo ) {
//---------------------------------------------------------


//- compute Background normalization for Likelihood
//  -----------------------------------------------
//- Paolo  -  December 2002


    CsRCRecConst* cons = CsRCRecConst::Instance();
    static double theMaxRgQ = cons->theMaxRgQ();

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

    float theMaxRg = cons->xuScan();

//- from back4-fit, 2002, run 20330
//  -------------------------------
//  bkPar0 : (p0 + p1*theIpo),  bkPar1 : (p2 + p3*theIpo)
//  bkgr = bkPar0*the + bkPar1*the**2
//  -----------------------------------------------------

    vector<float> par = CsRCEventPartPhotons::Instance()->vPara();

    double bkPar0 = par[0] + par[1] * theIpo;
    double bkPar1 = par[2] + par[3] * theIpo;
    double normB = dBindTheta *
      ( bkPar0* theMaxRgQ /2. + bkPar1* theMaxRgQ*sqrt( theMaxRgQ )/3. );
//----------------------------------------------------------------------
    if( normB < 0. ) normB = 0.;

    return  normB;

  }


//===========================================================================
  //double CsRCLikeAll02::likeSignal( const double thePhot,
  double CsRCLikeAll02::likeSignal( const CsRCPhoton* pPhot,
//----------------------------------------------------------
				    const double theIpo ) {
  //const double sigPhot ) {


//- compute Signal value for Likelihood
//  -----------------------------------
//- Paolo  -  December 2002


    CsRCRecConst* cons = CsRCRecConst::Instance();
    static double TwoPI = cons->TwoPI();
    static double sq2pi = sqrt( TwoPI );

    double thePhot = pPhot->the();
    double sigPhot = pPhot->sigmaPhoPid( pPartPhot_ );

    double qSw = (thePhot - theIpo) / sigPhot;
    double signal = exp( - 0.5* qSw*qSw ) / (sigPhot*sq2pi);

    if( theIpo <= 0. ) signal = 0.;

    return  signal;

  }


//===========================================================================
  //double CsRCLikeAll02::likeBackgr( const double thePhot,
  double CsRCLikeAll02::likeBackgr( const CsRCPhoton* pPhot,
//---------------------------------------------------------
				    const double theIpo ) {

//- compute Background value for Likelihood
//  ---------------------------------------
//- Paolo  -  December 2002

    double thePhot = pPhot->the();

    vector<float> par = CsRCEventPartPhotons::Instance()->vPara();

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
