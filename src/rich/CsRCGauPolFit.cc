/*!
   \file    CsRCGauPolFit.cc
   \-----------------------
   \brief   CsRCGauPolFit class implementation.
   \author  Paolo Schiavon
   \version 1.0
   \date    December 2003
*/


//----------------------------
  #include "CsRCChiSqFit.h"
  #include "CsRCGauPolFit.h"

  #include "CsRCHistos.h"

  #include "CsRCRecConst.h"
//--------------------------

  # include <ostream>

  #include "CLHEP/Matrix/Matrix.h"
  #include "CLHEP/Matrix/SymMatrix.h"
  #include "CLHEP/Matrix/DiagMatrix.h"
  #include "CLHEP/Matrix/Vector.h"

  using namespace std;
  using namespace CLHEP;

//===========================================================================
  CsRCGauPolFit::CsRCGauPolFit(): CsRCChiSqFit() {}
//-------------------------------------------------

//===========================================================================
  CsRCGauPolFit::CsRCGauPolFit( const int nPoint,
//-----------------------------------------------------------------------
				const double *xMeas, const double *yMeas,
                                const int nParam, const double *param,
				const int *iPaFit ):
    CsRCChiSqFit( nPoint, xMeas, yMeas, nParam, param, iPaFit ) {

    printKey_ = false;

    double error[nPoint];
    for( int kp=0; kp < nPoint; kp++ ) error[kp] = fabs( yMeas[kp] );

    getCovMax( error );
//  ------------------

    convTolSet_ = false;
  }


//===========================================================================
  CsRCGauPolFit::CsRCGauPolFit( const int nPoint,
//-----------------------------------------------------------------------
				const double *xMeas, const double *yMeas,
                                const double *error,
                                const int nParam, const double *param,
				const int *iPaFit ):
    CsRCChiSqFit( nPoint, xMeas, yMeas, nParam, param, iPaFit ) {

    printKey_ = false;
    //printKey_ = true;
    //std::cout << setprecision( 7 );

    getCovMax( error );
//  ------------------

    convTolSet_ = false;
  }


//===========================================================================
  CsRCGauPolFit::CsRCGauPolFit( const int nPoint,
//-----------------------------------------------------------------------
				const double *xMeas, const double *yMeas,
                                const double *error,
                                const int nParam, const double *param,
				const double *toll, const int *iPaFit ):
    CsRCChiSqFit( nPoint, xMeas, yMeas, nParam, param, iPaFit ) {

    printKey_ = false;
    //printKey_ = true;
    //std::cout << setprecision( 7 );

    getCovMax( error );
//  ------------------

    convTolSet_ = false;
    getConvTol( toll );
//  ------------------

  }


//===========================================================================
  void CsRCGauPolFit::print() const { CsRCChiSqFit::print(); }
//------------------------------------------------------------

//===========================================================================
  CsRCGauPolFit::~CsRCGauPolFit() {}
//----------------------------------


//===========================================================================
  void CsRCGauPolFit::getConvTol( const double *convTol ) {
//---------------------------------------------------------

//- Version : 1.0
//- Date : April 2006


    if( convTolSet_ ) return;
    HepVector VP( nPaF_, 0 );
    convTol_ = VP;
    int kPaFit = 0;
    for( int kParam=0; kParam < nParam_; kParam++ ) {
      if( iPaF_[kParam] == 0 ) {
	convTol_[kPaFit] = convTol[kParam];
	kPaFit++;
      }
    }
    //std::cout << "getConvTol" << setprecision(6) << std::endl;
    //for( int k=0; k<nPaF_; k++ ) std::cout << convTol_[k] << "  ";
    //std::cout << setprecision(2) << std::endl;

    convTolSet_ = true;

    return;
  }


//===========================================================================
  bool CsRCGauPolFit::doChiFit() {
//--------------------------------

//- Version : 1.0
//- Date : December 2003


//- Fit procedure :
//  ---------------
    if( !flagFit_ ) return  false;
    flagFit_ = true;

//- Start fitting procedure :
//    -------------------------
    if( !fitStart() ) return  false;
//      -----------
 
    int nIterMax = 10;

    bool doIter = true;
    while( doIter ) {

//--- Perform next iteration :
//    ------------------------
      nIter_++;

      if( !getDeriv() ) break;
//        -----------

      if( !fitIteration() ) break;
//        ---------------

      doIter = fitCheckConv();
//             --------------
      if( nIter_ >= nIterMax ) {
	doIter = false;
	flagFit_ = false;
      }
    }   /* end iteration loop */
    if( !flagFit_ ) return  flagFit_;

    if( !fitEnd() ) return  false;
//      ---------

    if( !fitPulls() ) return  false;
//      -----------

    return  flagFit_;

  }


//===========================================================================
  bool CsRCGauPolFit::fitStart() {
//--------------------------------

//--- Version : 1.0
//--- Date : December 2003


//--- Degrees of freedom of the fit :
//    -------------------------------
      nDegFree_ = nPoint_ - nPaF_;

      if( printKey_ ) {
	cout << endl;
        cout << "nPoint, nPara = " << nPoint_;
        cout << " , " << nPaF_ << endl;
      }

//--- Set convergence tolerancies on parameters to be fitted :
//    --------------------------------------------------------
      //static HepVector convTol( nParam_, 0 );
      double convTol[nParam_];
      convTol[0] = 0.1;
      convTol[1] = 0.00001;
      convTol[2] = 0.00001;
      convTol[3] = 0.1;
      convTol[4] = 0.00001;
      convTol[5] = 0.00001;
      convTol[6] = 0.1;
      convTol[7] = 0.1;
      convTol[8] = 0.1;
      convTol[9] = 0.1;

      getConvTol( convTol );
//-------------------------
      //HepVector VP( nPaF_, 0 );
      //convTol_ = VP;
      //int kPaFit = 0;
      //for( int kParam=0; kParam < nParam_; kParam++ ) {
      //  if( iPaF_[kParam] == 0 ) { 
      //    convTol_[kPaFit] = convTol[kParam];
      //    kPaFit++;
      //  }
      //}

//--- Compute initial residuals :
//    ---------------------------
      HepVector VY( nPoint_, 0 );
      Ya_ = VY;
      dYm_ = VY;
      double G1 = Pa_[0];
      double m1 = Pa_[1];
      double s1 = Pa_[2];
      double G2 = Pa_[3];
      double m2 = Pa_[4];
      double s2 = Pa_[5];
      double P0 = Pa_[6];
      double P1 = Pa_[7];
      double P2 = Pa_[8];
      double P3 = Pa_[9];
      for( int kp=0; kp < nPoint_; kp++ ) {
	double xms1 = 0.;
	if( s1 > 0. ) xms1 = (Xm_[kp] - m1) / s1;
        double expv1 = exp( - xms1*xms1 / 2. );
	double xms2 = 0.;
	if( s2 > 0. ) xms2 = (Xm_[kp] - m2) / s2;
        double expv2 = exp( - xms2*xms2 / 2. );
	Ya_[kp] = G1* expv1 + G2* expv2 + 
	  P0 + P1* Xm_[kp] + P2* pow( Xm_[kp], 2 ) + P3* pow( Xm_[kp], 3 );
	dYm_[kp] = Ya_[kp] - Ym_[kp];
      }

      if( printKey_ ) {
	cout << endl;
        cout << "Xm Start = " << Xm_.num_row();
	cout << Xm_.T() << endl;
	cout << "Ym Start = " << Ym_.num_row();
	cout << Ym_.T() << endl;
	cout << "dYm Start = " << dYm_.num_row();
        cout << dYm_.T() << endl;
        cout << "conv Toll = " << convTol_.num_row();
        cout << convTol_.T() << endl;
      }

      return  true;

  }


//===========================================================================
  bool CsRCGauPolFit::getDeriv() {
//--------------------------------

//--- Version : 1.0
//--- Date : December 2003


//--- Compute the derivative Matrix dYdX :
//    ------------------------------------
      int kPaFit = 0;
      for( int kParam=0; kParam < nParam_; kParam++ ) {
        if( iPaF_[kParam] == 0 ) { 
	  Pa_[kParam] = Xp_[kPaFit];
          kPaFit++;
        } else {
	  Pa_[kParam] = PaIn_[kParam];
	}
      }
      double G1 = Pa_[0];
      double m1 = Pa_[1];
      double s1 = Pa_[2];
      double G2 = Pa_[3];
      double m2 = Pa_[4];
      double s2 = Pa_[5];
      double P0 = Pa_[6];
      double P1 = Pa_[7];
      double P2 = Pa_[8];
      double P3 = Pa_[9];

      HepMatrix AA( nPoint_, nPaF_, 0 );
      mA_ = AA;
      for( int kp=0; kp < nPoint_; kp++ ) {
	double xms1 = 0.;
	if( s1 > 0. ) xms1 = (Xm_[kp] - m1) / s1;
        double expv1 = exp( - xms1*xms1 / 2. );
	double xms2 = 0.;
	if( s2 > 0. ) xms2 = (Xm_[kp] - m2) / s2;
        double expv2 = exp( - xms2*xms2 / 2. );
        int kq = 0;
	if( iPaF_[0] == 0 ) { mA_[kp][kq] = expv1;  kq++; }
	if( iPaF_[1] == 0 ) { mA_[kp][kq] = G1* expv1* xms1/ s1;  kq++; }
	if( iPaF_[2] == 0 ) { mA_[kp][kq] = G1* expv1* xms1*xms1/ s1;  kq++; }
	if( iPaF_[3] == 0 ) { mA_[kp][kq] = expv2;  kq++; }
	if( iPaF_[4] == 0 ) { mA_[kp][kq] = G2* expv2* xms2/ s2;  kq++; }
	if( iPaF_[5] == 0 ) { mA_[kp][kq] = G2* expv2* xms2*xms2/ s2;  kq++; }
	if( iPaF_[6] == 0 ) { mA_[kp][kq] = 1.;  kq++; }
	if( iPaF_[7] == 0 ) { mA_[kp][kq] = Xm_[kp];  kq++; }
	if( iPaF_[8] == 0 ) { mA_[kp][kq] = pow( Xm_[kp], 2 );  kq++; }
	if( iPaF_[9] == 0 ) { mA_[kp][kq] = pow( Xm_[kp], 3 );  kq++; }
      }

      if( printKey_ ) {
	cout << endl;
        cout << "derA = " << mA_.num_row() << " x " << mA_.num_col();
	cout << mA_ << endl;
      }

      return  true;

  }


//===========================================================================
  bool CsRCGauPolFit::fitIteration() {
//------------------------------------

//--- Version : 1.0
//--- Date : December 2003


      if( printKey_ ) {
        cout << endl << "------ Iter = " << nIter_ << endl;
      }

//--- Compute the covariance matrix on parameters Epf_
//    as Epf = inv(mAT * Wgm * mA) :
//    ------------------------------
      int invFlag;
      HepMatrix Wa = mA_.T() * Wgm_;
      Epf_ = (Wa * mA_).inverse( invFlag );
      if( invFlag != 0 ) flagFit_ = false;

      if( printKey_ ) {
	cout << endl;
        cout << "Epf = " << Epf_.num_row() << " x " << Epf_.num_col();
        cout << Epf_ << endl;
      }

      if( flagFit_ ) {

//----- Compute the vector of changes in parameters dXp_
//      as dXp = - Epf * mAT * Wgm * dYm :
//      ----------------------------------
        dXp_ = - (Epf_ * Wa * dYm_);
//----- and the vector of residuals dYm_, as dYm = dYm + mA * dXp :
//      -----------------------------------------------------------
        //dYm_ += mA_ * dXp_;
//----- and the new parameters Xp_, as Xp = Xp + dXp :
//      ----------------------------------------------
        Xp_ += dXp_;

//----- Compute the new full vector of parameters :
//      -------------------------------------------
	int kPaFit = 0;
        for( int kParam=0; kParam < nParam_; kParam++ ) {
          if( iPaF_[kParam] == 0 ) { 
	    Pa_[kParam] = Xp_[kPaFit];
            kPaFit++;
          } else {
	    Pa_[kParam] = PaIn_[kParam];
	  }
        }
        double G1 = Pa_[0];
        double m1 = Pa_[1];
        double s1 = Pa_[2];
        double G2 = Pa_[3];
        double m2 = Pa_[4];
        double s2 = Pa_[5];
        double P0 = Pa_[6];
        double P1 = Pa_[7];
        double P2 = Pa_[8];
        double P3 = Pa_[9];
//----- Compute the vector of new approximate var.s as Ya = F( Xp ) :
//      -------------------------------------------------------------
        for( int kp=0; kp < nPoint_; kp++ ) {
	  double xms1 = 0.;
	  if( s1 > 0. ) xms1 = (Xm_[kp] - m1) / s1;
          double expv1 = exp( - xms1*xms1 / 2. );
	  double xms2 = 0.;
	  if( s2 > 0. ) xms2 = (Xm_[kp] - m2) / s2;
          double expv2 = exp( - xms2*xms2 / 2. );
	  Ya_[kp] = G1* expv1 + G2* expv2 + 
	  P0 + P1* Xm_[kp] + P2* pow( Xm_[kp], 2 ) + P3* pow( Xm_[kp], 3 );
//----- and the vector of residuals dYm_, as dYm = Ya - Ym :
//      ----------------------------------------------------
	  dYm_[kp] = Ya_[kp] - Ym_[kp];
	}

        if( printKey_ )  { 
          cout << "dXp = " << dXp_.T() << endl;
	  cout << "Xp = " << Xp_.T() << endl;
          cout << "dYm = " << dYm_.T() << endl;
	}

      }

      return  flagFit_;

  }


//===========================================================================
  bool CsRCGauPolFit::fitCheckConv() {
//------------------------------------

//- Version : 1.0
//- Date : December 2003


    bool doIter = false;
    for( int kp=0; kp < nPaF_; kp++ )
      if( fabs( dXp_[kp] ) > convTol_[kp] ) doIter = true;

    return doIter;

  }


//===========================================================================
  bool CsRCGauPolFit::fitEnd() {
//------------------------------

//--- Version : 1.0
//--- Date : December 2003


//--- Compute Chi Square of the fit :
//    -------------------------------
      HepMatrix Ch = dYm_.T() * Wgm_ * dYm_;
      Ch /= nDegFree_;
      chiSquare_ = Ch[0][0];

//--- Output fitted parameters :
//    --------------------------
      int kPaFit = 0;
      for( int kParam=0; kParam < nParam_; kParam++ ) {
        if( iPaF_[kParam] == 0 ) { 
          Pa_[kParam] = Xp_[kPaFit];
          kPaFit++;
        } else {
	  Pa_[kParam] = PaIn_[kParam];
	}
      }
      Pf_ = Pa_;

      return  true;

  }


//===========================================================================
  bool CsRCGauPolFit::fitPulls() {
//--------------------------------

//--- Version : 1.0
//--- Date : December 2003


//--- Compute the fitted qq.s :
//    -------------------------
      Yf_ = Ym_ + dYm_;
      //cout << Yf_ << endl;

//--- Compute Covariance Matrix on fitted qq.s :
//    ------------------------------------------
      Erf_ = mA_ * Epf_ * mA_.T();
      //cout << Erf_ << endl;

//--- Compute Pulls :
//    ---------------
      HepVector VY( nPoint_, 0 );
      pully_ = VY;
      for( int kp=0; kp<nPoint_; kp++ ) {
	double den = Erm_[kp][kp] - Erf_[kp][kp];
      	pully_[kp] = 1000000.;
      	if( den > 0. ) pully_[kp] = (Ym_[kp] - Yf_[kp]) / sqrt( den );
	//cout << pully_[k] << "  ";
      }
      //cout << endl;

      return  true;

  }


//===========================================================================
  bool CsRCGauPolFit::doHist( const vector<CsHist2D*> vHisto ) {
//--------------------------------------------------------------

//--- Version : 1.0
//--- Date : May 2000, rev. 11/9/00, rev 19/4/01, rev 12/03 


    CsRCHistos& hist = CsRCHistos::Ref();
    float  xh, yh, wh;

    xh = Pa_[1];
    yh = Pa_[2];
    if( vHisto[0] ) vHisto[0]->Fill( xh, yh );
//                  -------------------------
    xh = PaIn_[0] - Pa_[0];
    yh = nPoint_;
    if( vHisto[1] ) vHisto[1]->Fill( xh, yh );
//                  -------------------------
    xh = chiSquare_;
    yh = nPoint_;
    if( vHisto[2] ) vHisto[2]->Fill( xh, yh );
//                  -------------------------
    xh = nIter_;
    yh = nPoint_;
    if( vHisto[3] ) vHisto[3]->Fill( xh, yh );
//                  -------------------------
    for( int kp=0; kp<nPoint_; kp++ ) {
      xh = pully_[kp];
      yh = nPoint_;
      if( vHisto[4] ) vHisto[4]->Fill( xh, yh );
//                    -------------------------
    }

    return  true;

  }
