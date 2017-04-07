/*!
   \File    CsRCLikeAll04.cc
   \----------------------------
   \brief   CsRCLikeAll04 class implementation.
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
  #include "CsRCLikeAll.h"
  #include "CsRCLikeAll04.h"

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
  CsRCLikeAll04::CsRCLikeAll04( const CsRCPartPhotons* papho ):
//-------------------------------------------------------------
    CsRCLikeAll() {

    pPartPhot_ = papho;

  }

//===========================================================================
  CsRCLikeAll04::~CsRCLikeAll04() {}
//----------------------------------

  extern "C" {
    float prob_( const float&, const int& );
    float erf_( const float& );
  }


//===========================================================================
  double CsRCLikeAll04::normSignal( const double theIpo ) {
//---------------------------------------------------------


//- compute Signal normalization for Likelihood
//  -------------------------------------------
//  dummy : assume likeRatio true
//- Paolo  -  May 2004


    CsRCRecConst* cons = CsRCRecConst::Instance();

    double normS = 0.;
    if( theIpo <= 0. ) normS = 0.;

    return  normS;

  }


//===========================================================================
  double CsRCLikeAll04::normBackgr( const double theIpo ) {
//---------------------------------------------------------


//- compute Background normalization for Likelihood
//  -----------------------------------------------
//  dummy : assume likeRatio true
//- Paolo  -  May 2004


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

    double normB = 1.;

    return  normB;

  }


//===========================================================================
  //double CsRCLikeAll04::likeSignal( const double thePhot,
  double CsRCLikeAll04::likeSignal( const CsRCPhoton* pPhot,
//---------------------------------------------------------
				    const double theIpo ) {
 //const double sigPhot ) {


//- compute Signal value for Likelihood
//  -----------------------------------
//- Paolo  -  May 2004


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
  double CsRCLikeAll04::likeBackgr( const CsRCPhoton* pPhot,
//----------------------------------------------------------
				    const double theIpo ) {

//- compute Background value for Likelihood
//  ---------------------------------------
//- Paolo  -  May 2004

    double thePhot = pPhot->the();

    static std::vector<float> par;
    static float paDeLm[7];
    static int nCo = 0;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      par = CsRCEventPartPhotons::Instance()->vPara();
      int i = 0;
      for( int k=0; k<12; k++ ) {
	//std::cout << "      ";
	for( int j=0; j<6; j++ ) {
	  //std::cout << par[i] << ",  ";
	  i++;
	}
	//std::cout << std::endl;
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

    int kDetPart = pPartPhot_->kDetPart();
    double ddPaD = pPartPhot_->pPart()->ddPaDet()[kDetPart];
    double theI = theIpo;
    double theP = thePhot;

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
