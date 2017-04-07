/*!
   \file    CsRCPad.cc
   \------------------
   \brief   CsRCPad class implementation.
   \author  Paolo Schiavon
   \version 0.01
   \date    October 2000, rev. September 2005
*/

    #include <iostream>
    #include <fstream>
    #include <iomanip>
    #include <sstream>
    #include <string>

    #include <ostream>
    #include <cstdio>

//----------------------------
    #include "CsRCPad.h"

    #include "CsRCDetectors.h"
//----------------------------

    #include "CsDigit.h"

    using namespace std;


//===========================================================================
    CsRCPad::CsRCPad() {}
//-----------------------

//===========================================================================
    CsRCPad::CsRCPad( int kPad, CsDigit* dig, int ic, int ix, int iy,
//-------------------------------------------------------------------
		      double PH ) {
      kPad_ = kPad;
      pDigit_ = dig;
      ic_ = ic;
      ix_ = ix;
      iy_ = iy;
      PH_ = PH;
      PHPack_ = 0.;
      PMTT0_ = 0.;
      PMTnum_ = 0;
      PMTcha_ = 0;
      flag_ = true;

      if( CsRCDetectors::Instance()->ptrToCat( ic )->isPMT() ) {
        PMTT0_ = PH_;
      }
      MCTime_ = 0.;
    }

//===========================================================================
    CsRCPad::CsRCPad( int kPad, int ic, int ix, int iy, double PH ) {
//-------------------------------------------------------------------
      kPad_ = kPad;
      pDigit_ = 0;
      ic_ = ic;
      ix_ = ix;
      iy_ = iy;
      PH_ = PH;
      PHPack_ = 0.;
      PMTT0_ = 0.;
      PMTnum_ = 0;
      PMTcha_ = 0;
      flag_ = true;

      CsRCDetectors* dets = CsRCDetectors::Instance();
      if( dets->ptrToCat( ic )->isPMT() ) {
        PMTT0_ = PH_;
      }
      MCTime_ = 0.;
    }

//===========================================================================
    CsRCPad::CsRCPad( const CsRCPad &pad ) {
//------------------------------------------
      //cout << "RICHONE : CsRCPad CopyConstructor" << endl;
      kPad_ = pad.kPad_;
      pDigit_ = pad.pDigit_;
      ic_ = pad.ic_;
      ix_ = pad.ix_;
      iy_ = pad.iy_;
      PH_ = pad.PH_;
      PHPack_ = pad.PHPack_;
      PMTT0_ = pad.PMTT0_;
      PMTnum_ = pad.PMTnum_;
      PMTcha_ = pad.PMTcha_;
      MCTime_ = pad.MCTime_;
      flag_ = pad.flag_;
      flagS_ = pad.flagS_;
    }

//===========================================================================
    void CsRCPad::print() const {
//-------------------------------
      cout << "   Pad  " << kPad_ << " : " << "  x-pos " << ix_
           << ",  y-pos " << iy_ << ",  cathode " << ic_
           << ",  PH " << PH_ << "  " << PHPack_ 
	   << ",  PMTT0 " << PMTT0_
	   << ",  PMTnum, cha " << PMTnum_ << "  " << PMTcha_
	   << ",  MCTime " << MCTime_
	   << ",  flag " << flag_ << endl;
    }

//===========================================================================
    CsRCPad::~CsRCPad() {}
//------------------------


//===========================================================================
  bool CsRCPad::setPMTpars() {
//----------------------------

//- Paolo  -  December 2006

    CsRCDetectors* dets = CsRCDetectors::Instance();
    CsRCCathode* cath = dets->ptrToCat( ic_ );
    if( !cath->isPMT() ) return  false;

//- PMT channel number is the same as the corresponding PPad number
//  (on the CsI plane)
    int iPPadx = ix_%4;
    int iPPady = iy_%4;
    int PMTcha = iPPady * 4 + iPPadx;
    if( PMTcha < 0  ||  PMTcha > 15 ) return  false;
    PMTcha_ = PMTcha;
    int iPMTx = int( ix_ / 4 );
    int iPMTy = int( iy_ / 4 );
    int PMTnum = iPMTy * cath->nPMTx() + iPMTx;
    //if( PMTnum < 0  ||  PMTnum > cath->nPMTx()*cath->nPMTy() )   // 080818
    //  return  false;
    if( PMTnum < 0  ||  PMTnum > (cath->nPMTx()*cath->nPMTy())-1 )
      return  false;
    PMTnum_ = PMTnum;

    //std::cout << "CsRCPad::setPMTpars " << ic_ << "  " << ix_
    //          << "  " << iy_ << "  " << iPPadx << "  " << iPPady
    //  	<< "  " << PMTnum << "  " << PMTcha << std::endl;

    return  true;

  }
