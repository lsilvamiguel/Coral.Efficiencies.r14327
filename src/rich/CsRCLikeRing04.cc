/*!
   \File    CsRCLikeRing04.cc
   \----------------------------
   \brief   CsRCLikeRing04 class implementation.
   \author  Paolo Schiavon
   \version 1.0
   \date    16 October 2003
*/

  #include <iostream>
  #include <ostream>
  #include <cstdio>

  #include "CsOpt.h"
  #include "CsHist.h"
  #include "CsHistograms.h"

//------------------------------
  #include "CsRCLikeRing.h"
  #include "CsRCLikeRing04.h"

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
  CsRCLikeRing04::CsRCLikeRing04( const CsRCRing* ring ): CsRCLikeRing() {
//------------------------------------------------------------------------

    pPartPhot_ = ring->pPartPhot();
    theUpL_ = ring->theUpL();
    theLoL_ = ring->theLoL();

  }

//===========================================================================
  CsRCLikeRing04::~CsRCLikeRing04() {}
//------------------------------------


  extern "C" { float erf_( const float& ); }

//===========================================================================
  double CsRCLikeRing04::normSignal( const double theIpo ) {
//----------------------------------------------------------


//- compute Signal normalization for Likelihood
//  -------------------------------------------
//  dummy : assume likeRatio true
//- Paolo  -  October 2003


//- average sigma-photon :
    CsRCPartPhotons* pt = const_cast<CsRCPartPhotons*>(pPartPhot_);
    //double momPart = pt->pPart()->mom();
    //double sigPhoRec = pt->sigmaPhoRec( momPart );
//@@                   --------------------------

    double normS = 0.;
    if( theIpo <= 0. ) normS = 0.;

    return  normS;

  }


//===========================================================================
  double CsRCLikeRing04::normBackgr( const double theReco,
//--------------------------------------------------------
				     const double theIpo ) {


//- compute Background normalization for Likelihood
//  -----------------------------------------------
//  dummy : assume likeRatio true
//- Paolo  -  October 2003


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
    CsRCPartPhotons* pt = const_cast<CsRCPartPhotons*>(pPartPhot_);
    //double momPart = pt->pPart()->mom();
    //double sigPhoRec = pt->sigmaPhoRec( momPart );
//                     --------------------------

    //double delta = theUpL_ - theLoL_;   //   useful wind.
//@@--------------------------------

    //vector<float> par = CsRCEventPartPhotons::Instance()->vPara();

    double normB = 1.;

    return  normB;

  }


//===========================================================================
  double CsRCLikeRing04::likeSignal( const double thePhot,
//--------------------------------------------------------
				     const double theIpo,
				     const double sigPhot ) {


//- compute Signal value for Likelihood
//  -----------------------------------
//- Paolo  -  October 2003


    CsRCRecConst* cons = CsRCRecConst::Instance();
    static double TwoPI = cons->TwoPI();
    static double sq2pi = sqrt( TwoPI );

    double qSw = (thePhot - theIpo) / sigPhot;
    double signal = exp( - 0.5* qSw*qSw ) / (sigPhot*sq2pi);

    if( theIpo <= 0. ) signal = 0.;

    return  signal;

  }


//===========================================================================
  //double CsRCLikeRing04::likeBackgr( const double thePhot,
  double CsRCLikeRing04::likeBackgr( CsRCPhoton* pPhot,
//----------------------------------------------------------
				     const double theIpo ) {


//- compute Background value for Likelihood
//  ---------------------------------------
//- Paolo  -  October 2003

    double thePhot = pPhot->the();

    static CsRCPartPhotons* pPartPhotPrec = NULL;
    static std::vector<float> par;
    static float paDeLm[7];
    static int nCo = 0;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      par = CsRCEventPartPhotons::Instance()->vPara();
      int i = 0;
      for( int k=0; k<12; k++ ) {
	std::cout << "      ";
	for( int j=0; j<6; j++ ) {
	  std::cout << par[i] << ",  ";
	  i++;
	}
	std::cout << std::endl;
      }

      paDeLm[0] =    0.;
      paDeLm[1] =  100.;
      paDeLm[2] =  200.;
      paDeLm[3] =  300.;
      paDeLm[4] =  400.;
      paDeLm[5] =  500.;
      paDeLm[6] = 5000.;

    }

//- from back5-fit, 2003, run 30580
//  ------------------------------
//  <bkgr> = par0 + par1*the + par2*the**2 + par3*the**3
//  ----------------------------------------------------
//  par0_3 are interpolated as a function of ddPaDet (6 points x 12)
//  and then                as a function of theIpo  (12 param.s )

    bool print = false;
    if( pPartPhot_ != pPartPhotPrec ) {
      //print = true;
      nCo = 0;
      pPartPhotPrec = const_cast<CsRCPartPhotons*>(pPartPhot_);
    }
    nCo++;
    //if( nCo == 300 ) print = true;

    int kDetPart = pPartPhot_->kDetPart();
    double ddPaD = pPartPhot_->pPart()->ddPaDet()[kDetPart];
    double theI = theIpo;
    double theP = thePhot;
    //ddPaD = 300.;
    //theI = 40.;
    //theP = 45.;

    double bkParD[12];
    for( int kp=0; kp<12; kp++ ) bkParD[kp] = 0.;
    for( int kp=0; kp<6; kp++ ) {
      if( ddPaD > paDeLm[kp]  &&  ddPaD <= paDeLm[kp+1] ) {
	for( int kq=0; kq<12; kq++ ) {
	  int kk = kp     * 12 + kq;
	  //int kj = (kp+1) * 12 + kq;
	  //bkParD[kq] = par[kk] + (par[kj]-par[kk]) * 
     	  //            (ddPaD-paDeLm[kp])/(paDeLm[kp+1]-paDeLm[kp]);
	  bkParD[kq] = par[kk];
	}
      }
    }
    if( print ) {
      //for( int kp=0; kp<12; kp++ ) std::cout << bkParD[kp] << "  ";
      //std::cout << std::endl;
    }
    double bkParT[4];
    for( int kt=0; kt<4; kt++ ) bkParT[kt] = 0.;
    for( int kt=0; kt<4; kt++ ) {
      bkParT[kt] = bkParD[3*kt] + bkParD[3*kt+1]*theI + 
	           bkParD[3*kt+2]*theI*theI;
    }
    if( print ) {
      //for( int kt=0; kt<4; kt++ ) std::cout << bkParT[kt] << "  ";
      //std::cout << std::endl;
    }

    double backgr = bkParT[0] + bkParT[1] * theP + bkParT[2] * theP*theP +
                    bkParT[3] * theP*theP*theP;
    if( print ) {
      //std:: cout << ddPaD << "  " << theI << "  " << backgr << std::endl;
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
  double CsRCLikeRing04::getRingBackground( const double theReco ) {
//------------------------------------------------------------------


//-  Paolo - October 2003


    CsRCRecConst * cons = CsRCRecConst::Instance();

    static float dBindTheta = 4.;
    static std::vector<float> par;
    static float paDeLm[7];

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

      par = CsRCEventPartPhotons::Instance()->vPara();
      //int i = 0;
      //for( int k=0; k<12; k++ ) {
      //  std::cout << "      ";
      //  for( int j=0; j<6; j++ ) {
      //    std::cout << par[i] << ",  ";
      //    i++;
      //  }
      //  std::cout << std::endl;
      //}

      paDeLm[0] =    0.;
      paDeLm[1] =  100.;
      paDeLm[2] =  200.;
      paDeLm[3] =  300.;
      paDeLm[4] =  400.;
      paDeLm[5] =  500.;
      paDeLm[6] = 5000.;

    }

    CsRCPartPhotons* pt = const_cast<CsRCPartPhotons*>(pPartPhot_);
    double momPart = pt->pPart()->mom();
    double sigPhoRec = pt->sigmaPhoRec( momPart );
//                     --------------------------
    double delta = theUpL_ - theLoL_;   //   useful window
//@@--------------------------------
    double normSK = 0.5* ( erf_( (theUpL_-theReco)/(sqrt( 2.)*sigPhoRec) )
			 - erf_( (theLoL_-theReco)/(sqrt( 2.)*sigPhoRec) ) );

    int kDetPart = pPartPhot_->kDetPart();
    double ddPaD = pPartPhot_->pPart()->ddPaDet()[kDetPart];
    double theI = theReco;
    double bkParD[12];
    for( int kp=0; kp<12; kp++ ) bkParD[kp] = 0.;
    for( int kp=0; kp<6; kp++ ) {
      if( ddPaD > paDeLm[kp]  &&  ddPaD <= paDeLm[kp+1] ) {
	for( int kq=0; kq<12; kq++ ) {
	  int kk = kp     * 12 + kq;
	  int kj = (kp+1) * 12 + kq;
	  bkParD[kq] = par[kk] + (par[kj]-par[kk]) * 
     	              (ddPaD-paDeLm[kp])/(paDeLm[kp+1]-paDeLm[kp]);
	}
      }
    }
    double bkParT[4];
    for( int kt=0; kt<4; kt++ ) bkParT[kt] = 0.;
    for( int kt=0; kt<4; kt++ ) {
      bkParT[kt] = bkParD[3*kt] + bkParD[3*kt+1]*theI + 
	           bkParD[3*kt+2]*theI*theI;
    }
//- approximate computation
    double normBK = dBindTheta * delta *
      ( bkParT[0] + bkParT[1] * theReco + bkParT[2] * theReco*theReco +
        bkParT[3] * theReco*theReco*theReco );

    double corrBack = 0.;
    if( fabs( normSK + normBK ) > 0. ) corrBack = normSK / (normSK + normBK);

    return  corrBack;

  }
