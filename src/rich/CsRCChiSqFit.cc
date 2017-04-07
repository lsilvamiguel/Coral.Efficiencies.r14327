/*!
   \file    CsRCChiSqFit.cc
   \-----------------------
   \brief   CsRCChiSqFit class implementation.
   \author  Paolo Schiavon
   \version 0.02
   \date    May 2000
*/


//----------------------------
  #include "CsRCChiSqFit.h"
//----------------------------

  # include <ostream>

  #include "CLHEP/Matrix/Matrix.h"
  #include "CLHEP/Matrix/SymMatrix.h"
  #include "CLHEP/Matrix/DiagMatrix.h"
  #include "CLHEP/Matrix/Vector.h"
//------------------------------------

  using namespace std;
  using namespace CLHEP;

//===========================================================================
  CsRCChiSqFit::CsRCChiSqFit() {}
//-------------------------------


//===========================================================================
  CsRCChiSqFit::CsRCChiSqFit( const int nPoint,
//---------------------------------------------------------------------
			      const double *xMeas, const double *yMeas,
                              const int nParam, const double *param,
			      const int *iPaFit ) {

//--- Fit type <1>
//    ------------

    if( nPoint <= 0 ) {
      flagFit_ = false;
      return;
    }
    printKey_ = false;

    nIter_ = 0;
    chiSquare_ = 10000000.;
    flagFit_ = true;

    nPoint_ = nPoint;
    nParam_ = nParam;

    HepVector VY( nPoint, 0 );
    Ym_ = VY;
    Xm_ = VY;
    for( int kPoint=0; kPoint < nPoint; kPoint++ ) {
      Xm_[kPoint] = xMeas[kPoint];
      Ym_[kPoint] = yMeas[kPoint];
    }

    HepVector VP( nParam, 0 );
    PaIn_ = VP;
    Pa_ = VP;
    iPaF_ = VP;
    for( int kParam=0; kParam < nParam; kParam++ ) {
      PaIn_[kParam] = param[kParam];
      Pa_[kParam] = param[kParam];
      iPaF_[kParam] = iPaFit[kParam];
    }
    int kPaFit = 0;
    for( int kParam=0; kParam < nParam; kParam++ ) {
      if( iPaFit[kParam] == 0 ) kPaFit++;
    }
    nPaF_ = kPaFit;
    HepVector VX( nPaF_, 0 );
    Xp_ = VX;
    kPaFit = 0;
    for( int kParam=0; kParam < nParam; kParam++ ) {
      if( iPaFit[kParam] == 0 ) { 
        Xp_[kPaFit] = param[kParam];  kPaFit++;
      }
    }
  }


//===========================================================================
  CsRCChiSqFit::CsRCChiSqFit( const int nPoint, const double *yMeas,
//------------------------------------------------------------------
                              const int nConst ) {

//--- Fit type <2>
//    ------------

    if( nPoint <= 0 ) {
      flagFit_ = false;
      return;
    }
    printKey_ = false;

    nIter_ = 0;
    chiSquare_ = 10000000.;
    flagFit_ = true;

    nPoint_ = nPoint;
    nConst_ = nConst;

    HepVector VY( nPoint, 0 );
    Ym_ = VY;
    for( int kPoint=0; kPoint < nPoint; kPoint++ ) {
      Ym_[kPoint] = yMeas[kPoint];
    }
  }


//===========================================================================
  CsRCChiSqFit::CsRCChiSqFit( const int nPoint,
//---------------------------------------------------------------------
			      const double *xMeas, const double *yMeas,
                              const int nConst, const int nParam,
			      const double *param, const int *iPaFit ) {

//--- Fit type <3>
//    ------------

    if( nPoint <= 0 ) {
      flagFit_ = false;
      return;
    }
    printKey_ = false;

    nIter_ = 0;
    chiSquare_ = 10000000.;
    flagFit_ = true;

    nPoint_ = 2* nPoint;
    nConst_ = nConst;
    nParam_ = nParam;

    HepVector VY( nPoint_, 0 );
    Ym_ = VY;
    for( int kPoint=0; kPoint<nPoint; kPoint++ ) {
      Ym_[2*kPoint  ] = xMeas[kPoint];
      Ym_[2*kPoint+1] = yMeas[kPoint];
    }

    HepVector VP( nParam, 0 );
    PaIn_ = VP;
    Pa_ = VP;
    iPaF_ = VP;
    for( int kParam=0; kParam < nParam; kParam++ ) {
      PaIn_[kParam] = param[kParam];
      Pa_[kParam] = param[kParam];
      iPaF_[kParam] = iPaFit[kParam];
    }
    int kPaFit = 0;
    for( int kParam=0; kParam < nParam; kParam++ ) {
      if( iPaFit[kParam] == 0 ) kPaFit++;
    }
    nPaF_ = kPaFit;
    HepVector VX( nPaF_, 0 );
    Xp_ = VX;
    kPaFit = 0;
    for( int kParam=0; kParam < nParam; kParam++ ) {
      if( iPaFit[kParam] == 0 ) { 
        Xp_[kPaFit] = param[kParam];  kPaFit++;
      }
    }
  }


//===========================================================================
  void CsRCChiSqFit::getCovMax( const double error ) {
//----------------------------------------------------

    if( !flagFit_ ) return;
    if ( error <= 0. ) {
      flagFit_ = false;
      return;
    }
    // Diagonal ???
    HepDiagMatrix Er( nPoint_, 1 );
    Er *= error;
    Erm_ = Er;
    //int invFlag;                               //   040422
    //Wgm_ = Er.inverse( invFlag );
    //if( invFlag != 0 ) flagFit_ = false;
    HepDiagMatrix Wg( nPoint_, 1 );
    Wg *= 1./error;
    Wgm_ = Wg;

    if( printKey_ ) {
      cout << endl;
      cout << "Wgm = " << Wgm_.num_row() << " x " << Wgm_.num_col();
      //      cout << Wgm_ << endl;
      cout << endl;
      for( int k=0; k < nPoint_; k++) cout << Wgm_[k][k] << " ";
      cout << endl;
    }
  }


//===========================================================================
  void CsRCChiSqFit::getCovMax( const double *error ) {
//-----------------------------------------------------

    if( !flagFit_ ) return;
    HepMatrix Er( nPoint_, nPoint_, 0 );
    //for( int kd=0; kd < nPoint_; kd++) Er[kd][kd] = error[kd];
    for( int kd=0; kd < nPoint_; kd++) {
      if( error[kd] <= 0. ) {
	flagFit_ = false;
	return;
      }
      Er[kd][kd] = error[kd];
    }
    Erm_ = Er;
    //int invFlag;                               //   040422
    //Wgm_ = Er.inverse( invFlag );
    //if( invFlag != 0 ) flagFit_ = false;
    HepMatrix Wg( nPoint_, nPoint_, 0 );
    for( int kd=0; kd < nPoint_; kd++) Wg[kd][kd] = 1./error[kd];
    Wgm_ = Wg;

    if( printKey_ ) {
      cout << endl;
      cout << "Wgm = " << Wgm_.num_row() << " x " << Wgm_.num_col();
      //      cout << Wgm_ << endl;
      cout << endl;
      for( int k=0; k < nPoint_; k++ ) cout << Wgm_[k][k] << " ";
      cout << endl;
    }
  }


//===========================================================================
  void CsRCChiSqFit::getCovMax( const double xError, const double yError ) {
//--------------------------------------------------------------------------

    if( !flagFit_ ) return;
    if( xError <= 0.  ||  yError <= 0. ) {
      flagFit_ = false;
      return;
    }
    // Diagonal ???
    HepDiagMatrix Er( nPoint_, 1 );
    int nPoint = nPoint_/2;
    for( int kd=0; kd<nPoint; kd++) {
      int ke = 2*kd;
      Er[ke][ke] = xError;
      ke++;
      Er[ke][ke] = yError;
    }
    Erm_ = Er;
    //int invFlag;                               //   040422
    //Wgm_ = Erm_.inverse( invFlag );
    //if( invFlag != 0 ) flagFit_ = false;
    HepDiagMatrix Wg( nPoint_, 1 );
    for( int kd=0; kd<nPoint; kd++) {
      int ke = 2*kd;
      Wg[ke][ke] = 1./xError;
      ke++;
      Wg[ke][ke] = 1./yError;
    }
    Wgm_ = Wg;

    if( printKey_ ) {
      cout << endl;
      cout << "Erm = " << Erm_.num_row() << " x " << Erm_.num_col();
      //cout << Erm_ << endl;
      cout << endl;
      for( int k=0; k<nPoint_; k++ ) cout << Erm_[k][k] << " ";
      cout << endl;
      cout << "Wgm = " << Wgm_.num_row() << " x " << Wgm_.num_col();
      //cout << Wgm_ << endl;
      cout << endl;
      for( int k=0; k<nPoint_; k++ ) cout << Wgm_[k][k] << " ";
      cout << endl;
    }
  }




//===========================================================================
  void CsRCChiSqFit::getCovMax( const double *xError, const double *yError )
//--------------------------------------------------------------------------
  {

    if( !flagFit_ ) return;
    HepMatrix Er( nPoint_, nPoint_, 0 );
    int nPoint = nPoint_/2;
    for( int kd=0; kd<nPoint; kd++) {
      if( xError[kd] <= 0.  ||  yError[kd] <= 0. ) {
	flagFit_ = false;
	return;
      }
      int ke = 2*kd;
      Er[ke][ke] = xError[kd];
      ke++;
      Er[ke][ke] = yError[kd];
    }
    Erm_ = Er;
    //int invFlag;                               //   040422
    //Wgm_ = Erm_.inverse( invFlag );
    //if( invFlag != 0 ) flagFit_ = false;
    HepMatrix Wg( nPoint_, nPoint_, 0 );
    for( int kd=0; kd<nPoint; kd++) {
      int ke = 2*kd;
      Wg[ke][ke] = 1./xError[kd];
      ke++;
      Wg[ke][ke] = 1./yError[kd];
    }
    Wgm_ = Wg;


    if( printKey_ ) {
      cout << endl;
      cout << "Erm = " << Erm_.num_row() << " x " << Erm_.num_col();
      //cout << Erm_ << endl;
      cout << endl;
      for( int k=0; k<nPoint_; k++ ) cout << Erm_[k][k] << " ";
      cout << endl;
      cout << "Wgm = " << Wgm_.num_row() << " x " << Wgm_.num_col();
      //cout << Wgm_ << endl;
      cout << endl;
      for( int k=0; k<nPoint_; k++ ) cout << Wgm_[k][k] << " ";
      cout << endl;
    }
  }


//===========================================================================
  void CsRCChiSqFit::getCovMax( const HepMatrix& error ) {
//--------------------------------------------------------

    if( !flagFit_ ) return;
    double det = error.determinant();
    if( det == 0. ) {
      flagFit_ = false;
      return;
    }
    Erm_ = error;
    int invFlag;
    Wgm_ = Erm_.inverse( invFlag );
    if( invFlag != 0 ) flagFit_ = false;

    if( printKey_ ) {
      cout << endl;
      cout << "Wgm = " << Wgm_.num_row() << " x " << Wgm_.num_col();
      //      cout << Wgm_ << endl;
      cout << endl;
      for( int k=0; k < nPoint_; k++ ) cout << Wgm_[k][k] << " ";
      cout << endl;
    }
  }

  bool CsRCChiSqFit::doChiFit() { return false; }
  bool CsRCChiSqFit::fitStart() { return false; }
  bool CsRCChiSqFit::getDeriv() { return false; }
  bool CsRCChiSqFit::fitIteration() { return false; }
  bool CsRCChiSqFit::fitCheckConv() { return false; }
  bool CsRCChiSqFit::fitEnd() { return false; }
  bool CsRCChiSqFit::fitPulls() { return false; }

  bool CsRCChiSqFit::doHist( vector<CsHist2D*> vHisto ) { return false; }

  void CsRCChiSqFit::print() const {
    cout << "chi = " << chiSquare_ << ",  ";
    cout << "nit = " << nIter_ << ",  ";
    cout << "nu = " << nDegFree_ << ",  ";
    cout << "ParFit = ";
    for( int k=0; k<nParam_; k++ ) cout << Pf_[k] << "  ";
    cout << endl;
  }

  CsRCChiSqFit::~CsRCChiSqFit() {}
