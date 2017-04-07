/*!
   \file    CsRCCircleFit.cc
   \-----------------------
   \brief   CsRCCircleFit class implementation.
   \author  Paolo Schiavon
   \version 0.02,  rev. 20/6/00
   \date    May 2000
*/


//----------------------------
  #include "CsRCChiSqFit.h"
  #include "CsRCCircleFit.h"

  #include "CsRCHistos.h"

  #include "CsRCRecConst.h"
//--------------------------

  #include <ostream>
  #include <iomanip>

  #include "CLHEP/Matrix/Matrix.h"
  #include "CLHEP/Matrix/SymMatrix.h"
  #include "CLHEP/Matrix/DiagMatrix.h"
  #include "CLHEP/Matrix/Vector.h"

  using namespace std;
  using namespace CLHEP;

//===========================================================================
 CsRCCircleFit::CsRCCircleFit(): CsRCChiSqFit() {}
//------------------------------------------------

//===========================================================================
  CsRCCircleFit::CsRCCircleFit( const int nPoint,
//-----------------------------------------------------------------------
				const double *xMeas, const double *yMeas,
                                const double error,
                                const int nParam, const double *param,
				const int *iPaFit ):
    CsRCChiSqFit( nPoint, xMeas, yMeas, nParam, param, iPaFit ) {

    printKey_ = false;

    getCovMax( error );
//  ------------------

  }


//===========================================================================
  CsRCCircleFit::CsRCCircleFit( const int nPoint,
//-----------------------------------------------------------------------
				const double *xMeas, const double *yMeas,
                                const double *error,
                                const int nParam, const double *param,
				const int *iPaFit ):
    CsRCChiSqFit( nPoint, xMeas, yMeas, nParam, param, iPaFit ) {

    printKey_ = false;
    //printKey_ = true;

    getCovMax( error );
//  ------------------

  }


//===========================================================================
  void CsRCCircleFit::print() const { CsRCChiSqFit::print(); }
//------------------------------------------------------------

  CsRCCircleFit::~CsRCCircleFit() {}



//===========================================================================
  bool CsRCCircleFit::doChiFit() {
//--------------------------------

//- Version : 1.0
//- Date : May 2000


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
  bool CsRCCircleFit::fitStart() {
//--------------------------------

//--- Version : 1.0
//--- Date : May 2000


//--- Degrees of freedom of the fit :
//    -------------------------------
      nDegFree_ = nPoint_ - nPaF_;

      if( printKey_ ) {
	cout << endl;
        cout << "nPoint, nPara = " << nPoint_;
        cout << " , " << nPaF_ << endl;
      }

//--- Set convergenze tolerancies on parameters to be fitted :
//    --------------------------------------------------------
      static HepVector convTol( nParam_, 0 );
      convTol[0] = 0.001;
      convTol[1] = 0.001;
      convTol[2] = 0.001;
      HepVector VP( nPaF_, 0 );
      convTol_ = VP;
      int kPaFit = 0;
      for( int kParam=0; kParam < nParam_; kParam++ ) {
        if( iPaF_[kParam] == 0 ) { 
          convTol_[kPaFit] = convTol[kParam];
          kPaFit++;
        }
      }

//--- Compute initial residuals :
//    ---------------------------
      HepVector VY( nPoint_, 0 );
      dYm_ = VY;
      double RR = Pa_[0];
      double x0 = Pa_[1];
      double y0 = Pa_[2];
      for( int kp=0; kp < nPoint_; kp++ ) {
        dYm_[kp] = sqrt( (Xm_[kp]-x0)*(Xm_[kp]-x0) + 
                         (Ym_[kp]-y0)*(Ym_[kp]-y0) ) - RR;
      }
      if( printKey_ ) {
	cout << endl;
        cout << "dYm Start = " << dYm_.num_row();
        cout << dYm_.T() << endl;
	cout << setprecision( 4 );
        cout << "conv Toll = " << convTol_.num_row();
        cout << convTol_.T() << endl;
	cout << setprecision( 2 );
      }

      return  true;

  }


//===========================================================================
  bool CsRCCircleFit::getDeriv() {
//--------------------------------

//--- Version : 1.0
//--- Date : May 2000


//--- Compute the derivative Matrix dYdX :
//    ------------------------------------
      int kPaFit = 0;
      for( int kParam=0; kParam < nParam_; kParam++ ) {
        if( iPaF_[kParam] == 0 ) { 
          Pa_[kPaFit] = Xp_[kParam];
          kPaFit++;
        }
      }
      double RR = Pa_[0];
      double x0 = Pa_[1];
      double y0 = Pa_[2];

      HepMatrix AA( nPoint_, nPaF_, 0 );
      mA_ = AA;
      for( int kp=0; kp < nPoint_; kp++ ) {
        double drr = sqrt( (Xm_[kp]-x0)*(Xm_[kp]-x0) + 
                           (Ym_[kp]-y0)*(Ym_[kp]-y0) );
        int kq = 0;
	if( iPaF_[0] == 0 ) { mA_[kp][kq] = - 1.;  kq++; }
	if( iPaF_[1] == 0 ) { mA_[kp][kq] = - (Xm_[kp]-x0)/drr;  kq++; }
	if( iPaF_[2] == 0 ) { mA_[kp][kq] = - (Ym_[kp]-y0)/drr;  kq++; }
      }

      if( printKey_ ) {
	cout << endl;
        cout << "derA = " << mA_.num_row() << " x " << mA_.num_col();
	cout << mA_ << endl;
      }

      return  true;

  }


//===========================================================================
  bool CsRCCircleFit::fitIteration() {
//------------------------------------

//--- Version : 1.0
//--- Date : May 2000


      if( printKey_ ) {
        cout << endl << "------ Iter = " << nIter_ << endl;
      }

//--- Compute the covariance matrix on parameters Epf_
//    as Epf = inv(mAT * Wgm * mA) :
//    ------------------------------
      int invFlag;
      HepMatrix Wa = mA_.T() * Wgm_;
      //HepMatrix Ww = Wa * mA_;
      //double det = Ww.determinant();
      //if( det <= 0. ) return  false;
      //Epf_ = Ww.inverse( invFlag );
      Epf_ = (Wa * mA_).inverse( invFlag );
      if( invFlag != 0 ) flagFit_ = false;

      if( printKey_ ) {
	cout << endl;
        cout << "Epf = " << Epf_.num_row() << " x " << Epf_.num_col();
        cout << Epf_ << endl;
      }

      if( flagFit_ ) {

//--- Compute the vector of changes in parameters dXp_
//    as dXp = - Epf * mAT * Wgm * dYm :
//    ----------------------------------
        dXp_ = - (Epf_ * Wa * dYm_);
//--- and the vector of residuals dYm_, as dYm = dYm + mA * dXp :
//    -----------------------------------------------------------
	dYm_ += mA_ * dXp_;
//--- and the new parameters Xp_, as Xp = Xp + dXp :
//    ----------------------------------------------
        Xp_ += dXp_;

        if( printKey_ )  { 
	  cout << setprecision( 4 );
          cout << "dXp = " << dXp_.T() << endl;
	  cout << setprecision( 2 );
	  cout << "Xp = " << Xp_.T() << endl;
          cout << "dYm = " << dYm_.T() << endl;
	}

      }

      return  flagFit_;

  }


//===========================================================================
  bool CsRCCircleFit::fitCheckConv() {
//------------------------------------

//- Version : 1.0
//- Date : May 2000


    bool doIter = false;
    for( int kp=0; kp < nPaF_; kp++ )
      if( fabs( dXp_[kp] ) > convTol_[kp] ) doIter = true;

    return doIter;

  }


//===========================================================================
  bool CsRCCircleFit::fitEnd() {
//------------------------------

//--- Version : 1.0
//--- Date : May 2000


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
        }
      }
      Pf_ = Pa_;

      return  true;

  }


//===========================================================================
  bool CsRCCircleFit::fitPulls() {
//--------------------------------

//--- Version : 1.0
//--- Date : April 2001


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
      for( int k=0; k<nPoint_; k++ ) {
	double den = Erm_[k][k] - Erf_[k][k];
      	pully_[k] = 1000000.;
      	if( den > 0. ) pully_[k] = (Ym_[k] - Yf_[k]) / sqrt( den );
	//cout << pully_[k] << "  ";
      }
      //cout << endl;

      return  true;

  }


//===========================================================================
  bool CsRCCircleFit::doHist( const vector<CsHist2D*> vHisto ) {
//--------------------------------------------------------------

//--- Version : 1.0
//--- Date : May 2000, rev. 11/9/00, rev 19/401


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
    for( int k=0; k<nPoint_; k++ ) {
      xh = pully_[k];
      yh = nPoint_;
      if( vHisto[4] ) vHisto[4]->Fill( xh, yh );
//                    -------------------------
    }

    return  true;

  }
