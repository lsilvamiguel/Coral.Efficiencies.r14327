/*!
   \File    CsRCLikeAll03.cc
   \----------------------------
   \brief   CsRCLikeAll03 class implementation.
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
  #include "CsRCLikeAll03.h"

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
  CsRCLikeAll03::CsRCLikeAll03( const CsRCPartPhotons* papho ):
//-------------------------------------------------------------
    CsRCLikeAll() {

    pPartPhot_ = papho;

  }

//===========================================================================
  CsRCLikeAll03::~CsRCLikeAll03() {}
//---------------------------------------------------------

  extern "C" {
    float prob_( const float&, const int& );
    float erf_( const float& );
  }


//===========================================================================
  double CsRCLikeAll03::normSignal( const double theIpo ) {
//---------------------------------------------------------


//- compute Signal normalization for Likelihood
//  -------------------------------------------
//- Paolo  -  December 2002
//  Rev. for TEST - March 2006         


    CsRCRecConst* cons = CsRCRecConst::Instance();
    float coo = cons->nZero() * 283. / 1000000.;

    //double normS = coo * theIpo*theIpo;
    double normS = 1;
//@@----------------

    return  normS;

  }


//===========================================================================
  double CsRCLikeAll03::normBackgr( const double theIpo ) {
//---------------------------------------------------------


//- compute Background normalization for Likelihood
//  -----------------------------------------------
//- Paolo  -  December 2002
//  Rev. for TEST - March 2006        


    CsRCRecConst* cons = CsRCRecConst::Instance();

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

//- from back5-fit, 2002, run 20330
//  -------------------------------
//  <bkgr> = par[0] + par[1]*the + par[2]*the**2 + par[3]*the**3
//  ------------------------------------------------------------

    vector<float> par = CsRCEventPartPhotons::Instance()->vPara();

    double normB = dBindTheta *
      ( par[0]* theMaxRg + par[1]* pow( theMaxRg, 2 ) +
        par[2]* pow( theMaxRg, 3 ) + par[3]* pow( theMaxRg, 4 ) );

    if( normB < 0. ) normB = 0.;
    normB = 0.;
//@@----------

    return  normB;

  }


//===========================================================================
  //double CsRCLikeAll03::likeSignal( const double thePhot,
  double CsRCLikeAll03::likeSignal( const CsRCPhoton* pPhot,
//----------------------------------------------------------
				    const double theIpo ) {
  //const double sigPhot ) {


//- compute Signal value for Likelihood
//  -----------------------------------
//- Paolo  -  December 2002
//  Rev. for TEST - March 2006        


    CsRCRecConst* cons = CsRCRecConst::Instance();
    static double TwoPI = cons->TwoPI();
    static double sq2pi = sqrt( TwoPI );
    float coo = cons->nZero() * 283. / 1000000.;

    double thePhot = pPhot->the();
    double sigPhot = pPhot->sigmaPhoPid( pPartPhot_ );

    double qSw = (thePhot - theIpo) / sigPhot;
    //double signal = coo * theIpo*theIpo * exp( - 0.5* qSw*qSw )
    //  / (sigPhot*sq2pi);
    double signal = exp( - 0.5* qSw*qSw ) / (sigPhot*sq2pi);
//@@-------------------------------------------------------

    if( theIpo <= 0. ) signal = 0.;

    return  signal;

  }


//===========================================================================
  //double CsRCLikeAll03::likeBackgr( const double thePhot,
  double CsRCLikeAll03::likeBackgr( const CsRCPhoton* pPhot,
//----------------------------------------------------------
				    const double theIpo ) {


//- compute Background value for Likelihood
//  ---------------------------------------
//- Paolo  -  December 2002
//  Rev. for TEST - March 2006        


    double thePhot = pPhot->the();

    vector<float> par = CsRCEventPartPhotons::Instance()->vPara();

    //double backgr = par[0] + par[1]* thePhot + par[2]* thePhot*thePhot +
    //  par[3]* thePhot*thePhot*thePhot;
    double backgr = 0.;
//@@------------------

    return  backgr;

  }


//===========================================================================
  double CsRCLikeAll03::getRingBackground( const double theReco ) {
//-----------------------------------------------------------------


//-  Paolo - January 2003


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

//- warning: theReco = thtaLikeMx!
//  ------------------------------
    //double momPart = pPart_->mom();
    //double sigPhoRec = sigmaPhoRec( momPart );
    CsRCPartPhotons* pt = const_cast<CsRCPartPhotons*>(pPartPhot_);
    double momPart = pt->pPart()->mom();
    double sigPhoRec = pt->sigmaPhoRec( momPart );
//                     --------------------------
    float dSigC = cons->sigCut() * sigPhoRec;
    double theLoL = theReco - dSigC;
    double theUpL = theReco + dSigC;
    double delta = theUpL - theLoL;   //   useful window
//@@------------------------------
    float coo = cons->nZero() * 283. / 1000000.;
    double normSK = 0.5* ( erf_( (theUpL-theReco)/(sqrt( 2.)*sigPhoRec) )
			 - erf_( (theLoL-theReco)/(sqrt( 2.)*sigPhoRec) ) );
    normSK *= coo * theReco*theReco;

    vector<float> par = CsRCEventPartPhotons::Instance()->vPara();
    //for( int k=0; k<4; k++ ) cout << par[k] << "  "; cout << endl;

    double normBK = dBindTheta * delta *
      ( par[0] + par[1]* theReco + par[2]* pow( theReco, 2 ) +
        par[3]* pow( theReco, 2 ) ) + 
      dBindTheta * pow( delta, 3 )/6. * ( par[2] + 3.* par[3]* theReco );

    double corrBack = 0.;
    if( fabs( normSK + normBK ) > 0. ) corrBack = normSK / (normSK + normBK);

    //cout << theUpL << "  " << theLoL << "  " << theReco << "  "
    //     << sigPhoRec << "  " << corrBack << endl;

    return  corrBack;

  }
