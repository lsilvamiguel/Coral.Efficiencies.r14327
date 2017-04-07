/*!
   \file    CsRCEllipseFitTest.cc
   \-----------------------------
   \brief   CsRCEllipseFitTest class implementation.
   \author  Paolo Schiavon
   \version 0.02,  rev. 20/6/00
   \date    May 2000
*/


//----------------------------
  #include "CsRCChiSqFit.h"
  #include "CsRCEllipseFitTest.h"

  #include "CsRichOne.h"
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
  CsRCEllipseFitTest::CsRCEllipseFitTest(): CsRCChiSqFit() {}
//---------------------------------------------------


//===========================================================================
  CsRCEllipseFitTest::CsRCEllipseFitTest( const int nPoint,
//-------------------------------------------------------------------------
				  const double *xMeas, const double *yMeas,
                                  const double xError, const double yError,
                                  const int nParam, const double *param,
				  const int *iPaFit ):
    CsRCChiSqFit( nPoint, xMeas, yMeas, nPoint, nParam, param, iPaFit ) {

    printKey_ = false;

    getCovMax( xError, yError );
//  ---------------------------

  }


//===========================================================================
  CsRCEllipseFitTest::CsRCEllipseFitTest( const int nPoint,
//---------------------------------------------------------------------------
				  const double *xMeas, const double *yMeas,
				  const double *xError, const double *yError,
				  const int nParam, const double *param,
				  const int *iPaFit ):
    CsRCChiSqFit( nPoint, xMeas, yMeas, nPoint, nParam, param, iPaFit ) {

    printKey_ = false;
    printKey_ = true;//*
    //if( CsRichOne::Instance()->kEvent() == 241 )  printKey_ = true;

    getCovMax( xError, yError );
//  ---------------------------

  }


//===========================================================================
  void CsRCEllipseFitTest::print() const { CsRCChiSqFit::print(); }
//-------------------------------------------------------------

  CsRCEllipseFitTest::~CsRCEllipseFitTest() {}



//===========================================================================
  bool CsRCEllipseFitTest::doChiFit() {
//---------------------------------

//- Version : 1.0
//- Date : January  2003


//- This is a Type <3> Fit
//  ---------------------

//- Fit procedure :
//  ---------------
    if( !flagFit_ ) return  false;
    kFail_ = 0;
    flagFit_ = true;

//- Start fitting procedure :
//  -------------------------
    if( !fitStart() ) return  false;
//      -----------
    cout << "111 " << endl;//*
    int nIterMax = 20;
//@@-----------------

    bool doIter = true;
    while( doIter ) {

//--- Perform next iteration :
//    ------------------------
      nIter_++;

      if( !getDeriv() ) break;
//        -----------
      cout << "222 " << endl;//*
      if( !fitIteration() ) break;
//        ---------------

      doIter = fitCheckConv();
//             --------------
      if( nIter_ >= nIterMax ) {
	doIter = false;
	kFail_ = 9;
	flagFit_ = false;
      }

    }   /* end iteration */
    if( !flagFit_ ) return  flagFit_;

    if( !fitEnd() ) return  false;
//      ---------

    if( !fitPulls() ) return  false;
//      -----------

    return flagFit_;
  }


//===========================================================================
  bool CsRCEllipseFitTest::fitStart() {
//-------------------------------------

//--- Version : 1.0
//--- Date : January  2003


//--- this is a Type <3> Fit
//    ----------------------
      if( printKey_ ) cout << "--- fitStart" << endl;//*

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
      HepVector convTolP( nParam_, 0 );
//@@  -------------------
      convTolP[0] = 0.05;
      convTolP[1] = 0.05;
      //convTolP[2] = 0.5;
      convTolP[2] = 1.0;//*
      convTolP[3] = 0.025;
      convTolP[4] = 0.2;
//@@  -------------------
      HepVector VP( nPaF_, 0 );
      convTolP_ = VP;
      int kPaFit = 0;
      for( int kParam=0; kParam<nParam_; kParam++ ) {
        if( iPaF_[kParam] == 0 ) { 
          convTolP_[kPaFit] = convTolP[kParam];
          kPaFit++;
        }
      }
      //for( int k=0; k<10; k++ ) std::cout << convTolP_[k] << std::endl;//*

      //*HepVector convTolY( nPoint_, 0 );
      //*for( int kPoint=0; kPoint<nPoint_; kPoint++ ) convTolY[kPoint] = 0.05;
//@@                                                -----------------------
      //*convTolY_ = convTolY;

//--- Compute initial residuals :
//    ---------------------------
      Ya_ = Ym_;
      dYm_ = Ya_ - Ym_;

      for( int kParam=0; kParam<nParam_; kParam++ ) {
        if( iPaF_[kParam] == 0 ) {
          Xp_[kPaFit] = Pa_[kParam];
          kPaFit++;
        }
      }

      if( printKey_ ) {
	cout << endl;
        cout << "dYm Start = " << dYm_.num_row();
        cout << dYm_.T() << endl;
        //*cout << "Var conv Toll = " << convTolY_.num_row();
        //*cout << convTolY_.T() << endl;
	cout << "Par conv Toll = " << convTolP_.num_row();
        cout << convTolP_.T() << endl;
      }

      return  true;

  }


//===========================================================================
  bool CsRCEllipseFitTest::getDeriv() {
//-------------------------------------

//--- Version : 1.0
//--- Date : January  2003


//--- this is a Type <3> Fit
//    ----------------------
      if( printKey_ ) cout << "--- getDeriv" << endl;//*
      int kPaFit = 0;
      for( int kParam=0; kParam < nParam_; kParam++ ) {
        if( iPaF_[kParam] == 0 ) { 
          Pa_[kParam] = Xp_[kPaFit];
          kPaFit++;
        } else {
          Pa_[kParam] = PaIn_[kParam];
	}
      }
      double x0 = Pa_[0];
      double y0 = Pa_[1];
      double AA = Pa_[2];
      double eps = Pa_[3];
      double phi = Pa_[4];
      double AA2 = AA*AA;
      double eps2 = eps*eps;

      double cosPhi = cos( phi );
      double sinPhi = sin( phi );
      double cosPhi2 = cosPhi*cosPhi;
      double sinPhi2 = sinPhi*sinPhi;

//--- Compute current Constraint values :
//    -----------------------------------
      HepVector VC( nConst_, 0 );
      Co_ = VC;
      for( int kConst=0; kConst<nConst_; kConst++ ) {

	double deltaX = ( Ya_[2*kConst  ] - x0 );
	double deltaY = ( Ya_[2*kConst+1] - y0 );
	Co_[kConst] = 
	    deltaX*deltaX * ( eps2*cosPhi2 + sinPhi2 )
	  - deltaX*deltaY * ( eps2-1.)* 2.* sinPhi * cosPhi
	  + deltaY*deltaY * ( cosPhi2 + eps2*sinPhi2 )
	  - eps2 * AA2;

      }
      cout << "113 " << endl;//*
//--- Compute the derivative Matrix dC/dX :
//    -------------------------------------
      cout << nConst_ << "  " << nPaF_ << endl;//*
      //HepMatrix mAA( nConst_, nPaF_, 0 );//*
      HepMatrix mAA( nConst_, nPaF_ );//*
      cout << mAA << endl;//*
      mA_ = mAA;
      for( int kConst=0; kConst<nConst_; kConst++ ) {
      cout << "11300 " << endl;//*
	double deltaX = ( Ya_[2*kConst  ] - x0 );
	double deltaY = ( Ya_[2*kConst+1] - y0 );
      cout << "1130 " << endl;//*
        int kPara = 0;
	if( iPaF_[0] == 0 ) {
//------- dC/d(x0):
	  mA_[kConst][kPara] = 
            - 2.* deltaX * ( eps2*cosPhi2 + sinPhi2 )
	    +     deltaY * ( eps2-1.)* 2.* sinPhi * cosPhi;
	  kPara++;
	}
      cout << "1131 " << endl;//*
	if( iPaF_[1] == 0 ) {
//------- dC/d(y0):
	  mA_[kConst][kPara] =
	           deltaX * ( eps2-1.)* 2.* sinPhi * cosPhi
	    - 2.*  deltaY * ( cosPhi2 + eps2*sinPhi2 );
	  kPara++;
	}
      cout << "1132 " << endl;//*
	if( iPaF_[2] == 0 ) {
//------- dC/d(AA):
	  double mA = - eps2;
	  mA_[kConst][kPara] = mA * 2.* AA;
	  kPara++;
	}
      cout << "1133 " << endl;//*
	if( iPaF_[3] == 0 ) {
//------- dC/d(eps):
	  double mA =
              deltaX*deltaX * cosPhi2
	    - deltaX*deltaY * 2.* sinPhi * cosPhi;
	    + deltaY*deltaY * sinPhi2
	    - AA2;
          mA_[kConst][kPara] = mA * 2.* eps;
	  kPara++;
	}
      cout << "1134 " << endl;//*
	if( iPaF_[4] == 0 ) {
//------- dC/d(phi):
           mA_[kConst][kPara] = 
	       deltaY*deltaY - deltaX*deltaX * ( eps2-1.)* 2.* sinPhi * cosPhi
	     - deltaX*deltaY * ( eps2-1.)* 2.* ( cosPhi2 - sinPhi2 );
	  kPara++;
	}
      cout << "1135 " << endl;//*
      }
      cout << "114 " << endl;//*
//--- Compute the derivative Matrix dC/dY :
//    -------------------------------------
      HepMatrix mBB( nConst_, nPoint_, 0 );
      mB_ = mBB;
      for( int kConst=0; kConst<nConst_; kConst++ ) {

	double deltaX = ( Ya_[2*kConst  ] - x0 );
	double deltaY = ( Ya_[2*kConst+1] - y0 );
	int kCooX = 2* kConst;
        mB_[kConst][kCooX] = 
	  2.* deltaX * ( eps2*cosPhi2 + sinPhi2 )
	  -   deltaY * ( eps2-1.)* 2.* sinPhi * cosPhi;
	int kCooY = 2* kConst + 1;
        mB_[kConst][kCooY] = 
	 -     deltaX * ( eps2-1.)* 2.* sinPhi * cosPhi
	 + 2.* deltaY * ( cosPhi2 + eps2*sinPhi2 );

      }
      cout << "115 " << endl;//*
      if( printKey_ ) {
	cout << endl;
        cout << "constr = " << Co_ << endl;
	cout << endl;
        cout << "derA = " << mA_.num_row() << " x " << mA_.num_col();
	cout << mA_ << endl;
        cout << "derB = " << mB_.num_row() << " x " << mB_.num_col();
	cout << mB_ << endl;
      }

      return  true;

  }


//===========================================================================
  bool CsRCEllipseFitTest::fitIteration() {
//-----------------------------------------

//--- Version : 1.0
//--- Date : January  2003

//--- this is a Type <3> Fit
//    ----------------------

      if( printKey_ ) {
        cout << endl << "------ Iter = " << nIter_ << endl;
      }
      if( printKey_ ) cout << "--- fitIteration" << endl;//*

//--- Compute constraint residuals Ro_
//    as Ro = Co - mB * dYm :
//    -----------------------
      Ro_ = Co_ - mB_ * dYm_;

//--- Compute weight matrix on constraints Eco_
//    as Eco = inv(mB * Erm * mB.T()) :
//    ---------------------------------
      int invFlag = 0;
      Eco_ = (mB_ * Erm_ * mB_.T()).inverse( invFlag );
      if( invFlag != 0 ) {
	kFail_ = 2;
	flagFit_ = false;
	return  flagFit_;
      }

//--- Compute the covariance matrix on parameters Epf_
//    as Epf = inv(mAT * Eco * mA) :
//    ------------------------------
      HepMatrix Wa = mA_.T() * Eco_;
      Epf_ = (Wa * mA_).inverse( invFlag );
      if( invFlag != 0 ) {
	kFail_ = 3;
	flagFit_ = false;
	return  flagFit_;
      }

      if( printKey_ ) {
	cout << endl;
        cout << "Epf = " << Epf_.num_row() << " x " << Epf_.num_col();
        cout << Epf_ << endl;
      }

//--- Compute the vector of changes in parameters dXp_ :
//    --------------------------------------------------
      dXp_ = - (Epf_ * Wa * Ro_);

//--- Compute the vector of Lagrange mult. Lk_ :
//    ------------------------------------------
      Lk_ = Eco_ * (Ro_ + mA_ * dXp_);

//--- Compute the vector of changes in var.s dYw_ :
//    ---------------------------------------------
      dYw_ = - Erm_ * mB_.T() * Lk_;

//--- Compute the vector of changes in var.s dYw_ :
//    the changes in var.s and the new residuals       
//    ---------------------------------------------
      dYd_ = dYw_ - dYm_;
      dYm_ = dYw_;

//--- Compute the new parameters Xp_ and var.s Ya_ :
//    ----------------------------------------------
      Xp_ += dXp_;
      Ya_ = Ym_ + dYm_;

      if( printKey_ ) {
	cout << "dXp = " << dXp_.T() << endl;
	cout << "Xp = " << Xp_.T() << endl;
	cout << "dYm = " << dYm_.T() << endl;
	cout << "Ya = " << Ya_.T() << endl;
	cout << "dYd = " << dYd_.T() << endl;
      }

      return  flagFit_;

  }


//===========================================================================
  bool CsRCEllipseFitTest::fitCheckConv() {
//-----------------------------------------

//- Version : 1.0
//- Date : January  2003


    bool doIter = false;
    //*for( int kp=0; kp<nPoint_; kp++ )
      //*if( fabs( dYd_[kp] ) > convTolY_[kp] ) doIter = true;

    for( int kp=0; kp<nPaF_; kp++ )
      if( fabs( dXp_[kp] ) > convTolP_[kp] ) doIter = true;

    return  doIter;

  }


//===========================================================================
  bool CsRCEllipseFitTest::fitEnd() {
//-----------------------------------

//--- Version : 1.0
//--- Date : January  2003


//--- Compute Chi Square of the fit :
//    -------------------------------
      HepMatrix Ch = (mA_ * dXp_ + Ro_).T() * Lk_;
      if( nDegFree_ > 0 ) {
	Ch /= nDegFree_;
      } else  return  false;
      chiSquare_ = Ch[0][0];

      if( printKey_ ) {
	cout << "chisq = " << chiSquare_ << endl;
      }

//--- Output fitted parameters :
//    --------------------------
      int kPaFit = 0;
      for( int kParam=0; kParam<nParam_; kParam++ ) {
        if( iPaF_[kParam] == 0 ) { 
          Pa_[kParam] = Xp_[kPaFit];
          kPaFit++;
        }
      }
      Pf_ = Pa_;

      if( printKey_ ) {
	cout << "fit par.s = " << Pf_ << endl;
      }

      return  true;

  }


//===========================================================================
  bool CsRCEllipseFitTest::fitPulls() {
//---------------------------------

//--- Version : 1.0
//--- Date : January  2003


//--- Compute the fitted var.s :
//    --------------------------
      Yf_ = Ym_ + dYm_;
      //cout << Yf_ << endl;

//--- Compute Covariance Matrix on fitted var.s :
//    -------------------------------------------
      HepMatrix W1 = mB_ * Erm_;
      HepMatrix W2 = mA_.T() * Eco_ * W1;
      Erf_ = Erm_ - W1.T() * Eco_ * W1 + W2.T() * Epf_ * W2;
      //cout << Erf_ << endl;

      if( printKey_ ) {
	cout << "fit cov mx = " << Erf_ << endl;
      }

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

      if( printKey_ ) {
	cout << "pulls = " << pully_ << endl;
      }

      return  true;

  }


//===========================================================================
  bool CsRCEllipseFitTest::doHist( const vector<CsHist2D*> vHisto ) {
//---------------------------------------------------------------

//--- Version : 1.0
//--- Date : February  2003


    CsRCRecConst *cons = CsRCRecConst::Instance();
    CsRCHistos& hist = CsRCHistos::Ref();
    float  xh, yh, wh;

    xh = kFail_;
    yh = nIter_;
    if( vHisto[5] ) vHisto[5]->Fill( xh, yh );
//                  -------------------------

    if( !flagFit_ ) return  false;

    double xPade = PaIn_[5];
    double yPade = PaIn_[6];
    double ddPaDet = 0.;
    if( yPade >= 0. )
      ddPaDet = sqrt( xPade*xPade + (yPade-470.)*(yPade-470.) );
    if( yPade < 0. )
      ddPaDet = sqrt( xPade*xPade + (yPade+470.)*(yPade+470.) );

    xh = Pa_[0];
    yh = Pa_[1];
    if( vHisto[0] ) vHisto[0]->Fill( xh, yh );
//                  -------------------------
    xh = Pa_[2] - PaIn_[2];
    yh = ddPaDet;
    if( vHisto[1] ) vHisto[1]->Fill( xh, yh );
//                  -------------------------
    xh = Pa_[3];
    yh = (Pa_[4] - PaIn_[4]) * cons->RadDeg();
    if( vHisto[2] ) vHisto[2]->Fill( xh, yh );
//                  -------------------------
    xh = chiSquare_;
    yh = nPoint_;
    if( vHisto[3] ) vHisto[3]->Fill( xh, yh );
//                  -------------------------
    xh = nIter_;
    yh = nPoint_;
    if( vHisto[4] ) vHisto[4]->Fill( xh, yh );
//                  -------------------------
    for( int k=0; k<nPoint_; k++ ) {
      xh = pully_[k];
      yh = nPoint_;
      if( vHisto[6] ) vHisto[6]->Fill( xh, yh );
//                    -------------------------
    }
    double phiPa = 0.;
    if( xPade != 0. ) phiPa = atan( yPade/xPade );
      else  phiPa = 1.5708;
    xh = (Pa_[4] - phiPa) * cons->RadDeg();
    yh = phiPa * cons->RadDeg();
    if( vHisto[7] ) vHisto[7]->Fill( xh, yh );
//                  -------------------------

    if( printKey_ ) {
      cout << "histo.s = " << endl;
    }

    return  true;

  }
