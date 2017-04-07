/*!
   \file    CsRCCathode.cc
   \----------------------
   \brief   CsRCCathode class implementation.
   \author  Paolo Schiavon
   \version 0.01
   \date    21 February 2001, rev. August 2005
*/

#   include <ostream>
#   include <cstdio>
#   include <string>    

#   include "CsErrLog.h"

//-----------------------------
#   include "CsRCCathode.h"
#   include "CsRCDetectors.h"
#   include "CsRCPhotonDet.h"

#   include "CsRCRecConst.h"
#   include "CsRCExeKeys.h"
#   include "CsRCHistos.h"
//-----------------------------

    using namespace std;
	using namespace CLHEP;

//===========================================================================
    CsRCCathode::CsRCCathode( int kCat, int ID, string TBname, string name, 
              double xOffCat0, double yOffCat0, double zOffCat0,
              int nPadx, int nPady, double padx, double pady,
	      double ddQzW, double ddGap, HepMatrix corrRotMx ) {
//---------------------------------------------------------------

//--- positions in mm, MRS (MAIN REFERENCE SYSTEM)
//    --------------------
      kCat_ = kCat;
      ID_ = ID;
      TBname_ = TBname;
      name_ = name;
      char sU = 'U';
      char sD = 'D';
      string ss;
      ss = name;
      ss = ss.assign( ss, 1, 2 );
      number_ = atoi( ss.c_str() );
      if( name[0] == sU ) number_--;
      if( name[0] == sD ) number_ += 8 - 1;
      vOffCat0_.setX( xOffCat0 );
      vOffCat0_.setY( yOffCat0 );
      vOffCat0_.setZ( zOffCat0 );
      list<CsRCPhotonDet*> lPhoDets = CsRCDetectors::Instance()->lPhoDet();
      CsRCPhotonDet* det = NULL;
      if( name[0] == sU ) det = lPhoDets.front();
      if( name[0] == sD ) det = lPhoDets.back();
      vCat0_.setX( det->vDet0().x() + xOffCat0 );
      vCat0_.setY( det->vDet0().y() + yOffCat0 );
      vCat0_.setZ( det->vDet0().z() + zOffCat0 );
//--- number of pads per cathode, pad dimensions
      nPadx_ = nPadx;
      nPady_ = nPady;
      padx_ = padx;
      pady_ = pady;
      nPPadx_ = 0;
      nPPady_ = 0;
      ppadx_ = 0.;
      ppady_ = 0.;
      nPMTx_ = 0;
      nPMTy_ = 0;
      PMTx_ = 0.;
      PMTy_ = 0.;
//--- compute half cathode distance
      hCatx_ = padx * nPadx /2.;
      hCaty_ = pady * nPady /2.;
      xLimMnW_ = 0.;
      xLimMxW_ = 0.;
      yLimMnW_ = 0.;
      yLimMxW_ = 0.;
//--- Quartz window, gap
      ddQzW_ = ddQzW;
      ddGap_ = ddGap;
//--- cathode rot-matrix, angle, direction
      rotMatrix_ = det->rotMatrix() + corrRotMx;
//--- cathode type
      char sM = 'M';
      char sP = 'P';
      char sA = 'A';
      isPMT_ = false;
      //if( TBname[1] == sP ) isPMT_ = true;
      if( TBname[1] == sM  ||  TBname[1] == sP ) isPMT_ = true;
      isAPV_ = false;
      if( TBname[1] == sA ) isAPV_ = true;

//--- PMT-cathode local geometry
//    --------------------------
//    on PD plane, CRS (CATHODE REFERENCE SYSTEM), origin in x/yOffCatW
//    x-axis points always to Jura, y-axis always upwards.
//
      if( isPMT_ ) {
	nPPadx_ = nPadx_;
	nPPady_ = nPady_;
	nPMTx_ = 12;
	nPMTy_ = 12;
	PMTx_ = 48.2;
	PMTy_ = 46.7;
	int PPPMx = nPPadx_ / nPMTx_;
	int PPPMy = nPPady_ / nPMTy_;
	ppadx_ = PMTx_ / float( PPPMx );
	ppady_ = PMTy_ / float( PPPMy );

//----- T(top), B(bottom); J(Jura), S(Saleve)
//      TJ->3, TS->5, BJ->10. BS->12
        static const int TJ =  3;
        static const int TS =  5;
        static const int BJ = 10;
        static const int BS = 12;
        for( int k=0; k<nPMTx_; k++ ) xcePMTWC_[k] = 0.;
        for( int k=0; k<nPMTy_; k++ ) ycePMTWC_[k] = 0.;

//----- on PD plane, PRS (PMT REFERENCE SYSTEM), origin in x/ycePMTWC
//      x-axis points always to Jura, y-axis always upwards.
//----  2 3
//      0 1   everywhere
        static double xCornerPMT[4] = { 0., 0., 0., 0. };
        static double yCornerPMT[4] = { 0., 0., 0., 0. };
//----  20 21 22 23 24
//      15 16 17 18 19
//      10 11 12 13 14
//       5  6  7  8  9
//       0  1  2  3  4  everywhere
	static int const nCoPPad = (PPPMx+1) * (PPPMy+1);
	static const double x[25] = { -23.94, -11.97, 0., 11.97, 23.94,
				      -23.94, -11.97, 0., 11.97, 23.94,
				      -23.94, -11.97, 0., 11.97, 23.94,
				      -23.94, -11.97, 0., 11.97, 23.94,
				      -23.94, -11.97, 0., 11.97, 23.94, };
	static const double y1[5] = { -22.67, -22.96, -23.24, -23.53, -23.82};
	static const double y2[5] = { -11.04, -11.33, -11.62, -11.91, -12.20};
	static const double y3[5] = {   0.58,   0.29,   0.00, - 0.29, - 0.58};
	static const double y4[5] = {  12.20,  11.91,  11.62,  11.33,  11.04};
	static const double y5[5] = {  23.82,  23.63,  23.24,  22.96,  22.67};
        for( int k=0; k<nCoPPad; k++ ) xCoPPadPMT_[k] = 0.;
        for( int k=0; k<nCoPPad; k++ ) yCoPPadPMT_[k] = 0.;
	static int const ncePPad = PPPMx * PPPMy;
        for( int k=0; k<ncePPad; k++ ) xcePPadPMT_[k] = 0.;
        for( int k=0; k<ncePPad; k++ ) ycePPadPMT_[k] = 0.;

//----- TS   ( 5 )
	if( number_ == TS ) {
	  // x =  25.09 + (k-6) * padx; y = -25.56 - (k-6) * pady;
	  for( int k=0; k<nPMTx_; k++ ) xcePMTWC_[k] = -264.1 + k * PMTx_;
	  for( int k=0; k<nPMTy_; k++ ) ycePMTWC_[k] = -259.1 + k * PMTy_;
	  xCornerPMT[0] = -23.94;   yCornerPMT[0] = -22.67;
	  xCornerPMT[1] =  23.94;   yCornerPMT[1] = -23.82;
	  xCornerPMT[2] = -23.94;   yCornerPMT[2] =  23.82;
	  xCornerPMT[3] =  23.94;   yCornerPMT[3] =  22.67;
	  for( int k=0; k<25; k++ ) xCoPPadPMT_[k] = x[k];
	  for( int k= 0; k< 5; k++ ) yCoPPadPMT_[k] = y1[k   ];
	  for( int k= 5; k<10; k++ ) yCoPPadPMT_[k] = y2[k- 5];
	  for( int k=10; k<15; k++ ) yCoPPadPMT_[k] = y3[k-10];
	  for( int k=15; k<20; k++ ) yCoPPadPMT_[k] = y4[k-15];
	  for( int k=20; k<25; k++ ) yCoPPadPMT_[k] = y5[k-20];
	  for( int k=0; k<ncePPad; k++ ) {
	    int j = k + k/4;
            xcePPadPMT_[k] = ( xCoPPadPMT_[j] + xCoPPadPMT_[j+1]
			   + xCoPPadPMT_[j+5] + xCoPPadPMT_[j+6] ) /4.;
            ycePPadPMT_[k] = ( yCoPPadPMT_[j] + yCoPPadPMT_[j+1]
			   + yCoPPadPMT_[j+5] + yCoPPadPMT_[j+6] ) /4.;
	  }

	}
//----- TJ   ( 3 )
	if( number_ == TJ ) {
	  // x = -25.09 - (k-6) * padx; y = -25.56 - (k-6) * pady;
	  for( int k=0; k<nPMTx_; k++ ) xcePMTWC_[k] = -266.1 + k * PMTx_;
	  for( int k=0; k<nPMTy_; k++ ) ycePMTWC_[k] = -259.1 + k * PMTy_;
	  xCornerPMT[0] = -23.94;   yCornerPMT[0] = -23.82;
	  xCornerPMT[1] =  23.94;   yCornerPMT[1] = -22.67;
	  xCornerPMT[2] = -23.94;   yCornerPMT[2] =  22.67;
	  xCornerPMT[3] =  23.94;   yCornerPMT[3] =  23.82;
	  for( int k=0; k<25; k++ ) xCoPPadPMT_[k] = x[k];
	  for( int k= 0; k< 5; k++ ) yCoPPadPMT_[k] = - y5[k   ];
	  for( int k= 5; k<10; k++ ) yCoPPadPMT_[k] = - y4[k- 5];
	  for( int k=10; k<15; k++ ) yCoPPadPMT_[k] = - y3[k-10];
	  for( int k=15; k<20; k++ ) yCoPPadPMT_[k] = - y2[k-15];
	  for( int k=20; k<25; k++ ) yCoPPadPMT_[k] = - y1[k-20];
	  for( int k=0; k<ncePPad; k++ ) {
	    int j = k + k/4;
            xcePPadPMT_[k] = ( xCoPPadPMT_[j] + xCoPPadPMT_[j+1]
			   + xCoPPadPMT_[j+5] + xCoPPadPMT_[j+6] ) /4.;
            ycePPadPMT_[k] = ( yCoPPadPMT_[j] + yCoPPadPMT_[j+1]
			   + yCoPPadPMT_[j+5] + yCoPPadPMT_[j+6] ) /4.;
	  }
	}
//----- BS   ( 12 )
	if( number_ == BS ) {
	  // x =  25.09 + (k-6) * padx; y =  25.56 + (k-6) * pady;
	  for( int k=0; k<nPMTx_; k++ ) xcePMTWC_[k] = -264.1 + k * PMTx_;
	  for( int k=0; k<nPMTy_; k++ ) ycePMTWC_[k] = -254.6 + k * PMTy_;
	  xCornerPMT[0] = -23.94;   yCornerPMT[0] = -23.82;
	  xCornerPMT[1] =  23.94;   yCornerPMT[1] = -22.67;
	  xCornerPMT[2] = -23.94;   yCornerPMT[2] =  22.67;
	  xCornerPMT[3] =  23.94;   yCornerPMT[3] =  23.82;
	  for( int k=0; k<25; k++ ) xCoPPadPMT_[k] = x[k];
	  for( int k= 0; k< 5; k++ ) yCoPPadPMT_[k] = - y5[k   ];
	  for( int k= 5; k<10; k++ ) yCoPPadPMT_[k] = - y4[k- 5];
	  for( int k=10; k<15; k++ ) yCoPPadPMT_[k] = - y3[k-10];
	  for( int k=15; k<20; k++ ) yCoPPadPMT_[k] = - y2[k-15];
	  for( int k=20; k<25; k++ ) yCoPPadPMT_[k] = - y1[k-20];
	  for( int k=0; k<ncePPad; k++ ) {
	    int j = k + k/4;
            xcePPadPMT_[k] = ( xCoPPadPMT_[j] + xCoPPadPMT_[j+1]
			   + xCoPPadPMT_[j+5] + xCoPPadPMT_[j+6] ) /4.;
            ycePPadPMT_[k] = ( yCoPPadPMT_[j] + yCoPPadPMT_[j+1]
			   + yCoPPadPMT_[j+5] + yCoPPadPMT_[j+6] ) /4.;
	  }
	}
//----- BJ   ( 10 )
	if( number_ == BJ ) {
	  // x = -25.09 - (k-6) * padx; y =  25.56 + (k-6) * pady;
	  for( int k=0; k<nPMTx_; k++ ) xcePMTWC_[k] = -266.1 + k * PMTx_;
	  for( int k=0; k<nPMTy_; k++ ) ycePMTWC_[k] = -254.6 + k * PMTy_;
	  xCornerPMT[0] = -23.94;   yCornerPMT[0] = -22.67;
	  xCornerPMT[1] =  23.94;   yCornerPMT[1] = -23.82;
	  xCornerPMT[2] = -23.94;   yCornerPMT[2] =  23.82;
	  xCornerPMT[3] =  23.94;   yCornerPMT[3] =  22.67;
	  for( int k=0; k<25; k++ ) xCoPPadPMT_[k] = x[k];
	  for( int k= 0; k< 5; k++ ) yCoPPadPMT_[k] = y1[k   ];
	  for( int k= 5; k<10; k++ ) yCoPPadPMT_[k] = y2[k- 5];
	  for( int k=10; k<15; k++ ) yCoPPadPMT_[k] = y3[k-10];
	  for( int k=15; k<20; k++ ) yCoPPadPMT_[k] = y4[k-15];
	  for( int k=20; k<25; k++ ) yCoPPadPMT_[k] = y5[k-20];
	  for( int k=0; k<ncePPad; k++ ) {
	    int j = k + k/4;
            xcePPadPMT_[k] = ( xCoPPadPMT_[j] + xCoPPadPMT_[j+1]
			   + xCoPPadPMT_[j+5] + xCoPPadPMT_[j+6] ) /4.;
            ycePPadPMT_[k] = ( yCoPPadPMT_[j] + yCoPPadPMT_[j+1]
			   + yCoPPadPMT_[j+5] + yCoPPadPMT_[j+6] ) /4.;
	  }

	}

	plotGeo();
//---------------

      }      // end if isPMT_

    }

//===========================================================================
    CsRCCathode::CsRCCathode( const CsRCCathode &cat ) {
//------------------------------------------------------
      cout << "RICHONE : CsRCCathode CopyConstructor" << endl;
      kCat_ = cat.kCat_;
      ID_ = cat.ID_;
      TBname_ = cat.TBname_;
      name_ = cat.name_;
      vOffCat0_ = cat.vOffCat0_;
      vCat0_ = cat.vCat0_;
      vOffCatW_ = cat.vOffCatW_;
//--- number of pads per cathode, pad dimensions
      nPadx_ = cat.nPadx_;
      nPady_ = cat.nPady_;
      padx_ = cat.padx_;
      pady_ = cat.pady_;
      nPPadx_ = cat.nPPadx_;
      nPPady_ = cat.nPPady_;
      ppadx_ = cat.ppadx_;
      ppady_ = cat.ppady_;
      nPMTx_ = cat.nPMTx_;
      nPMTy_ = cat.nPMTy_;
      PMTx_ = cat.PMTx_;
      PMTy_ = cat.PMTy_;
      hCatx_ = cat.hCatx_;
      hCaty_ = cat.hCaty_;
      xLimMnW_ = cat.xLimMnW_;
      xLimMxW_ = cat.xLimMxW_;
      yLimMnW_ = cat.yLimMnW_;
      yLimMxW_ = cat.yLimMxW_;
      ddQzW_ = cat.ddQzW_;
      ddGap_ = cat.ddGap_;
//--- detector angles (to detector directions), detector centres
      rotMatrix_ = cat.rotMatrix_;

      isPMT_ = cat.isPMT_;
      isAPV_ = cat.isAPV_;
      for( int k=0; k<nPMTx_; k++ ) xcePMTWC_[k] = cat.xcePMTWC_[k];
      for( int k=0; k<nPMTy_; k++ ) ycePMTWC_[k] = cat.ycePMTWC_[k];
      for( int k=0; k<25; k++ ) xCoPPadPMT_[k] = cat.xCoPPadPMT_[k];
      for( int k=0; k<25; k++ ) yCoPPadPMT_[k] = cat.yCoPPadPMT_[k];
      for( int k=0; k<16; k++ ) xcePPadPMT_[k] = cat.xcePPadPMT_[k];
      for( int k=0; k<16; k++ ) ycePPadPMT_[k] = cat.ycePPadPMT_[k];
    }

//===========================================================================
    void CsRCCathode::print() const {
//-----------------------------------
      //      cout << endl;
      //      cout << "Cathode geometry :" << endl
      //           << "------------------" << endl;
      cout << "Cathode  " << name_ << "  nr " << kCat_
	   << "   ID " << ID_ << "  TBname  " << TBname_
           << ",  pads/cathode x/y = " << nPadx_ << " x " << nPady_
	   << ",  pad dimensions = " << padx_ << " x " << pady_ << endl;
      cout << "halfCathode x/y  = " << hCatx_ << "  " << hCaty_ << endl;
      cout << "cathode centre (MRS) = "<< vCat0_ << endl;
      cout << "cathode centre (MWR) = "<< vOffCatW_ << endl;
      cout << "quartz window width = " << ddQzW_ << ",  chamber gap = "
           << ddGap_ << endl;
      cout << "rotation matrix (MRS) = ";
      for( int i=0; i<3; i++ ) for( int j=0; j<3; j++ ) 
        cout << rotMatrix_(i+1,j+1) << "  ";
      cout << endl;
      cout << "isPMT = " << isPMT_ << "   isAPV = " << isAPV_ << endl;
      if( isPMT_ ) {
	cout << "PMT centre-x (CRS) = ";
	for( int k=0; k<nPadx_; k++ ) cout << xcePMTWC_[k] << ", ";
	cout << endl;
	cout << "PMT centre-y (CRS) = ";
	for( int k=0; k<nPady_; k++ ) cout << ycePMTWC_[k] << ", ";
	cout << endl;
	cout << "PMT grid-x (PRS) = " << endl;
	int i = 0;
	for( int k=0; k<5; k++ ) {
	  for( int j=0; j<5; j++ ) { cout << xCoPPadPMT_[i] << ", "; i++; }
	  cout << endl;
	}
	cout << "PMT grid-y (PRS) = " << endl;
	i = 0;
	for( int k=0; k<5; k++ ) {
	  for( int j=0; j<5; j++ ) { cout << yCoPPadPMT_[i] << ", "; i++; }
	  cout << endl;
	}
	cout << "PPad centre-x (PRS) = " << endl;
	i = 0;
	for( int k=0; k<4; k++ ) {
	  for( int j=0; j<4; j++ ) { cout << xcePPadPMT_[i] << ", "; i++; }
	  cout << endl;
	}
	cout << "PPad centre-y (PRS) = " << endl;
	i = 0;
	for( int k=0; k<4; k++ ) {
	  for( int j=0; j<4; j++ ) { cout << ycePPadPMT_[i] << ", "; i++; }
	  cout << endl;
	}
      }
      cout << endl;
    }

//===========================================================================
    CsRCCathode::~CsRCCathode() {}
//--------------------------------


//===========================================================================
    bool CsRCCathode::ccePMTWC( const int kPMT, double& xce, double& yce ) {
//--------------------------------------------------------------------------

//--- Paolo - August  2005


//--- Beware : kPMT from 0 to 143 (from Saleve to Jura  and  bottom to top)
//             ix, iy from 0 to 11 (same);   all cathodes
      xce = 0.;
      yce = 0.;
      if( !isPMT_ ) return  false;

      int iy = kPMT / 12;
      int ix = kPMT - iy * 12;
      xce = xcePMTWC_[ix];
      yce = ycePMTWC_[iy];

      return  true;

    }


//===========================================================================
    bool CsRCCathode::ccePPadWC( const int kPPadx, const int kPPady,
//------------------------------------------------------------------
				 double& xce, double& yce ) {

//--- Paolo - August  2005


//--- Beware : kPPadx/y from 0 to 47 (from Saleve to Jura  and  bottom to top)
//             all cathodes
      xce = 0.;
      yce = 0.;
      if( !isPMT_ ) return  false;

      int iPMTx = kPPadx / 4;
      int iPPdax = kPPadx - iPMTx * 4;
      int iPMTy = kPPady / 4;
      int iPPday = kPPady - iPMTy * 4;
      int iPPad = iPPday * 4 + iPPdax;
      xce = xcePMTWC_[iPMTx] + xcePPadPMT_[iPPad];
      yce = ycePMTWC_[iPMTy] + ycePPadPMT_[iPPad];

      return  true;

    }

//===========================================================================
    void CsRCCathode::plotGeo() const {
//-------------------------------------

//--- Paolo - March  2006

      //if( !CsRCHistos::Ref().bookHis()  ||
      //  CsRCHistos::Ref().levelBk() < 2 ) return;
//    --------------------------------------------

      CsRCRecConst *cons = CsRCRecConst::Instance();
      int kOffHi = cons->kOffHi();

      static std::vector<CsHist2D*> vRC1040;
      static int kCall = 0;

      static bool firstCall = true;
      if( firstCall ) {
        firstCall = false;

	CsHistograms::SetCurrentPath("/RICH");
	std::string hTitle = "cathode geo";
	for( int kH=0; kH<4; kH++ ) {
	  stringstream hN1040;
	  int kHist = kOffHi + 1040 + kH;
	  hN1040 << kHist;
          vRC1040.push_back( new CsHist2D( hN1040.str(), hTitle,
					   200, -50., 50., 200, -50., 50.) );
        }
	int kH;
	int kHist;
	kH = 4;
	int nch = 600;
	float xch = 300.;
	float ych = 300.;
	stringstream hN1044;
	kHist = kOffHi + 1040 + kH;
	hN1044 << kHist;
	vRC1040.push_back( new CsHist2D( hN1044.str(), hTitle,
					 nch, 0., xch, nch, 0., ych) );
        kH = 5;
	stringstream hN1045;
	kHist = kOffHi + 1040 + kH;
	hN1045 << kHist;
	vRC1040.push_back( new CsHist2D( hN1045.str(), hTitle,
					 nch, -xch, 0., nch, 0., ych) );
	kH = 6;
	stringstream hN1046;
	kHist = kOffHi + 1040 + kH;
	hN1046 << kHist;
	vRC1040.push_back( new CsHist2D( hN1046.str(), hTitle,
					 nch, 0., xch, nch, -ych, 0.) );
	kH = 7;
	stringstream hN1047;
	kHist = kOffHi + 1040 + kH;
	hN1047 << kHist;
	vRC1040.push_back( new CsHist2D( hN1047.str(), hTitle,
					 nch, -xch, 0., nch, -ych, 0.) );
        CsHistograms::SetCurrentPath("/");

      }

      int kHist = -1;
      if( number_ ==  3 ) kHist = 0;
      if( number_ ==  5 ) kHist = 1;
      if( number_ == 10 ) kHist = 2;
      if( number_ == 12 ) kHist = 3;
      if( kHist < 0 ) return;
      if( kHist >= int( vRC1040.size() ) ) return;      //   100726

      double xh, yh;
      for( int kc=0; kc<16; kc++ ) {
	xh = xcePPadPMT_[kc];
	yh = ycePPadPMT_[kc];
	if( vRC1040[kHist] ) vRC1040[kHist]->Fill( xh, yh );
//hh                         ------------------------------
      }
      kCall++;

      kHist += 4;
      if( kHist >= int( vRC1040.size() ) ) return;      //   100726
      float xOffCatW[16];
      float yOffCatW[16];
      float xOff = 288.;
      float yOff = 288.;
      xOffCatW[ 3] =  xOff;
      xOffCatW[ 5] = -xOff;
      xOffCatW[10] =  xOff;
      xOffCatW[12] = -xOff;
      yOffCatW[ 3] =  yOff;
      yOffCatW[ 5] =  yOff;
      yOffCatW[10] = -yOff;
      yOffCatW[12] = -yOff;
      for( int ixPMT=0; ixPMT<12; ixPMT++ ) {
	for( int iyPMT=0; iyPMT<12; iyPMT++ ) {
	  for( int iPPad=0; iPPad<16; iPPad++ ) {
	    xh = xcePPadPMT_[iPPad] + xcePMTWC_[ixPMT] + xOffCatW[number_];
	    yh = ycePPadPMT_[iPPad] + ycePMTWC_[iyPMT] + yOffCatW[number_];
	    if( vRC1040[kHist] ) vRC1040[kHist]->Fill( xh, yh );
//hh                             ------------------------------
	  }
	}
      }

      //std::cout << "CsRCCathode::plot " << kCall << std::endl;
      //if( kCall == 4 ) exit(0);
      
      return;
    }
