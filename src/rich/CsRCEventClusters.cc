
/*!
   \file    CsRCEventClusters.cc
   \----------------------------
   \brief   CsRCEventClusters class author.
   \implementation  Paolo Schiavon
   \version 0.1.1
   \date    October 2000, rev. August 2005
*/


  #include <iostream>
  #include <ostream>
  #include <cstdio>
  #include <cmath>
  #include <cstdlib>

  #include "CsErrLog.h"

  #include "CsMCHit.h"
  #include "CsMCDigit.h"
//------------------------------
  #include "CsRCEventClusters.h"
  #include "CsRCCluster.h"

  #include "CsRCEventPads.h"
  #include "CsRCPad.h"

  #include "CsRCHistos.h"

  #include "CsRCDetectors.h"
  #include "CsRCRecConst.h"
  #include "CsRCExeKeys.h"

  #include "CsRCChiSqFit.h"
  #include "CsRCCircleFit.h"
//------------------------------

  using namespace std;

  CsRCEventClusters* CsRCEventClusters::instance_ = 0;

//===========================================================================
  CsRCEventClusters* CsRCEventClusters::Instance() {
//--------------------------------------------------
    if( instance_ == 0 ) instance_ = new ( CsRCEventClusters );
    return instance_;
  }

//===========================================================================
  CsRCEventClusters::CsRCEventClusters() {
//----------------------------------------
  }

//===========================================================================
  CsRCEventClusters::CsRCEventClusters( const CsRCEventClusters &evclus ) {
//-------------------------------------------------------------------------
    cout << "RICHONE : CsRCEventClusters CopyConstructor" << endl;
    instance_ = evclus.instance_;
    lClusters_ = evclus.lClusters_;
    flag_ = evclus.flag_;
  }

//===========================================================================
  void CsRCEventClusters::clearEventClusters() {
//----------------------------------------------
    list<CsRCCluster*>::iterator ic;
    for( ic=lClusters_.begin(); ic!=lClusters_.end(); ic++ ) delete (*ic);
    lClusters_.clear();
  }

//===========================================================================
  void CsRCEventClusters::print() const {
//---------------------------------------
    cout << endl;
    cout << " Clusters : " << endl;
    list<CsRCCluster*>::const_iterator ic;
    for( ic=lClusters_.begin(); ic!=lClusters_.end(); ic++ ) (*ic)->print();
  }

//===========================================================================
  void CsRCEventClusters::printFull() const {
//-------------------------------------------
    cout << endl;
    cout << " Clusters : " << endl;
    list<CsRCCluster*>::const_iterator ic;
    for( ic=lClusters_.begin(); ic!=lClusters_.end(); ic++ ) {
      (*ic)->print();
      list<CsRCPad*> lPads = (*ic)->lPads();
      list<CsRCPad*>::iterator ia;
      for( ia=lPads.begin(); ia!=lPads.end(); ia++ ) (*ia)->print();
    }
  }

//===========================================================================
  void CsRCEventClusters::printFull( const int &nc ) const {
//----------------------------------------------------------
    cout << endl;
    cout << " Clusters : " << endl;
    list<CsRCCluster*>::const_iterator ic;
    int kc = 0;
    for( ic=lClusters_.begin(); ic!=lClusters_.end(); ic++ ) {
      (*ic)->print();
      list<CsRCPad*> lPads = (*ic)->lPads();
      list<CsRCPad*>::iterator ia;
      for( ia=lPads.begin(); ia!=lPads.end(); ia++ ) {
	cout << (*ia) << "  ";
	(*ia)->print();
      }
      if( kc++ >= nc ) break;
    }
  }

//===========================================================================
  CsRCEventClusters::~CsRCEventClusters() {
//-----------------------------------------
    clearEventClusters();
  }


//==========================================================================  
  void CsRCEventClusters::getEventClusters() {
//--------------------------------------------


//--- Paolo  -  May 2001


      CsRCExeKeys *key = CsRCExeKeys::Instance();
      CsRCHistos& hist = CsRCHistos::Ref();
      int kh;
      double xh, yh, wh;

      static bool killHalo = false;

      static bool firstCall = true;
      if( firstCall ) {
	firstCall = false;
      }

      xh = 30.5;
      if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                       ------------------------

      if( key->DoClu() ) doClustering();
//                       --------------
      else  getClusFromPads();
//          -----------------

      if( key->KillHaloClus() ) killHaloClus();
//                              --------------

      int nCluEv = lClusters_.size();

//--- monitoring histograms :
//    -----------------------
      xh = nCluEv;
      if( hist.hRC1526 ) hist.hRC1526->Fill( xh );
//hh                     ------------------------
      list<CsRCCluster*>::iterator cl;
      for( cl=lClusters_.begin(); cl!=lClusters_.end(); cl++ ) {
        xh = (*cl)->xc();
        yh = (*cl)->yc();
        if( hist.hRC1525 ) hist.hRC1525->Fill( xh, yh );
//hh                       ----------------------------
        xh = (*cl)->PH();
        yh = (*cl)->ic();
        if( hist.hRC1527 ) hist.hRC1527->Fill( xh, yh );
//hh                       ----------------------------
      }
      list<CsRCPad*> lPadAll = CsRCEventPads::Instance()->lPads();
      int nPadEv = lPadAll.size();
      xh = float( nCluEv ) / float( nPadEv );
      if( hist.hRC1528 ) hist.hRC1528->Fill( xh );
//hh                     ------------------------

      checkClusters();
//-------------------
      checkPMTClus();
//------------------

//--- conditional prints :
//    --------------------
      int kPrintEventClus = key->kPrintEventClus();
      if( kPrintEventClus == 1 ) {
        cout << endl;
        cout << " Clusters : " << nCluEv << endl;
      }
      if( kPrintEventClus == 2 ) print();

      xh = 31.5;
      if( hist.hRC1000 ) hist.hRC1000->Fill( xh, float( nPadEv ) );
//                       -----------------------------------------
      xh = 32.5;
      if( hist.hRC1000 ) hist.hRC1000->Fill( xh, float( nCluEv ) );
//                       -----------------------------------------
      if( lClusters_.empty() ) {
	xh = 33.5;
	if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                         ------------------------
      }

      xh = 39.5;
      if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                       ------------------------

      return;

  }


//==========================================================================
  void CsRCEventClusters::getClusFromPads() {
//-------------------------------------------

//--- Paolo
//--- Version 0.1.1 - October 2000


      CsRCExeKeys *key = CsRCExeKeys::Instance();

      CsRCHistos& hist = CsRCHistos::Ref();
      double xh, yh, wh;

      CsRCDetectors *dets = CsRCDetectors::Instance();
      list<CsRCCathode*> lCathodes = dets->lCathodes();

      static bool firstCall = true;
      if( firstCall ) {
        firstCall = false;
        key->acknoMethod( "CsRCEventClusters::getClusFromPads" );
      }

      flag_ = false;

      list<CsRCPad*> lPadAll = CsRCEventPads::Instance()->lPads();
      list<CsRCPad*>::iterator ip;
      for( ip=lPadAll.begin(); ip!=lPadAll.end(); ip++ ) {
	if( !(*ip)->flag() ) continue;

        int kClu = lClusters_.size();
        list<CsRCPad*> lPads;
        lPads.clear();
        lPads.push_back( (*ip) );
        int ic = (*ip)->ic();

	CsRCCathode* cat = dets->ptrToCat( ic );
	double padx = cat->padx();
	double pady = cat->pady();
	double hCatx = cat->hCatx();
	double hCaty = cat->hCaty();
	double xc = 0.;
	double yc = 0.;
	double PH = 0.;
	if( cat->isPMT() ) {
          xc = (*ip)->ix() + 0.5;
          yc = (*ip)->iy() + 0.5;
	  xc = (xc*padx - hCatx) + dets->vOffCatW( ic ).x();
          yc = (yc*pady - hCaty) + dets->vOffCatW( ic ).y();
          //---- alternate code obsolete
          //int ix = (*ip)->ix();
          //int iy = (*ip)->iy();
	  //int PPPMx = cat->nPPadx() / cat->nPMTx();
	  //int PPPMy = cat->nPPady() / cat->nPMTy();
	  //int ixPMT = ix / PPPMx;
	  //int ixPad = ix - ixPMT * PPPMx;
	  //int iyPMT = iy / PPPMy;
	  //int iyPad = iy - iyPMT * PPPMy;
	  //int iPPad = iyPad * PPPMx + ixPad;
	  //xc = cat->xcePMTWC()[ixPMT] + cat->xcePPadPMT()[iPPad] +
  	  //     dets->vOffCatW( ic ).x();
	  //yc = cat->ycePMTWC()[iyPMT] + cat->ycePPadPMT()[iPPad] +
  	  //     dets->vOffCatW( ic ).y();
          PH = (*ip)->PH();
	}
	else {
          xc = (*ip)->ix() + 0.5;
          yc = (*ip)->iy() + 0.5;
	  xc = (xc*padx - hCatx) + dets->vOffCatW( ic ).x();
          yc = (yc*pady - hCaty) + dets->vOffCatW( ic ).y();
          PH = (*ip)->PH();
	}

        lClusters_.push_back( new CsRCCluster( kClu, lPads,
					       ic, xc, yc, PH ) );

	xh = double( lPads.size() );
	if( hist.hRC1001 ) hist.hRC1001->Fill( xh );
//hh                       ------------------------
      }

      if( !lClusters_.empty() ) flag_ = true;
  }


//==========================================================================
  void CsRCEventClusters::getClusFromPads( int jCat ) {
//-----------------------------------------------------

//--- Paolo
//--- Version 1.0 - August 2005


      CsRCExeKeys *key = CsRCExeKeys::Instance();

      CsRCHistos& hist = CsRCHistos::Ref();
      double xh, yh, wh;

      CsRCDetectors *dets = CsRCDetectors::Instance();
      list<CsRCCathode*> lCathodes = dets->lCathodes();

      static bool firstCall = true;
      if( firstCall ) {
        firstCall = false;
        key->acknoMethod( "CsRCEventClusters::getClusFromPads(cat)" );
      }
      list<CsRCPad*> lPadAll = CsRCEventPads::Instance()->lPads();
      list<CsRCPad*>::iterator ip;
      for( ip=lPadAll.begin(); ip!=lPadAll.end(); ip++ ) {
	if( !(*ip)->flag() ) continue;

	if( (*ip)->ic() != jCat ) continue;

        int kClu = lClusters_.size();
        list<CsRCPad*> lPads;
        lPads.clear();
        lPads.push_back( (*ip) );
        int ic = (*ip)->ic();

	CsRCCathode* cat = dets->ptrToCat( ic );
	double padx = cat->padx();
	double pady = cat->pady();
	double hCatx = cat->hCatx();
	double hCaty = cat->hCaty();

	double xc = 0.;
	double yc = 0.;
	double PH = 0.;
	if( cat->isPMT() ) {
          double xc0 = (*ip)->ix() + 0.5;
          double yc0 = (*ip)->iy() + 0.5;
	  xc0 = (xc0*padx - hCatx) + dets->vOffCatW( ic ).x();
          yc0 = (yc0*pady - hCaty) + dets->vOffCatW( ic ).y();
	  //std::cout << "getClusFromPads : " << ic << "  " << xc0 << " "
	  //    << yc0 << " --- ";
	  if( !cat->ccePPadWC( (*ip)->ix(), (*ip)->iy(), xc, yc ) ) return;
	  xc += dets->vOffCatW( ic ).x();
	  yc += dets->vOffCatW( ic ).y();
	  //std::cout << xc << " " << yc << std::endl;
	  xh = xc - xc0;
	  yh = yc - yc0;
	  if( hist.hRC1022 ) hist.hRC1022->Fill( xh, yh );
//hh                         ----------------------------
          PH = (*ip)->PH();
          //---- alternate code obsolete
          //int ix = (*ip)->ix();
          //int iy = (*ip)->iy();
	  //int PPPMx = cat->nPPadx() / cat->nPMTx();
	  //int PPPMy = cat->nPPady() / cat->nPMTy();
	  //int ixPMT = ix / PPPMx;
	  //if( ixPMT < 0  || ixPMT > 11 ) {
	  //  std::cout << " CsRCEventClusters::getClusFromPads :"
	  //	      << " PMT out of range, ixPMT = " << ixPMT
	  //	      << std::endl;
	  //  continue;
	  //}
	  //int ixPad = ix - ixPMT * PPPMx;
	  //int iyPMT = iy / PPPMy;
	  //if( iyPMT < 0  || iyPMT > 11 ) {
          //  std::cout << " CsRCEventClusters::getClusFromPads :"
  	  //	      << " PMT out of range, iyPMT = " << iyPMT
	  //	      << std::endl;
	  //  continue;
	  //}
	  //int iyPad = iy - iyPMT * PPPMy;
	  //int iPPad = iyPad * PPPMx + ixPad;
	  //if( iPPad < 0  || iPPad > 15 ) {
          //  std::cout << " CsRCEventClusters::getClusFromPads :"
	  //	      << " PMT channel out of range, iPPad = " << iPPad
	  //	      << std::endl;
	  //  continue;
	  //}
	  //xc = cat->xcePMTWC()[ixPMT] + cat->xcePPadPMT()[iPPad] +
  	  //     dets->vOffCatW( ic ).x();
	  //yc = cat->ycePMTWC()[iyPMT] + cat->ycePPadPMT()[iPPad] +
  	  //     dets->vOffCatW( ic ).y();

//------- test 20061110
/*
	  float xH = 0.;
	  float yH = 0.;
	  float tH = 0.;
	  CsDigit* dg = (*ip)->pDigit();
	  list<CsMCHit*> padMCHits = 
	    dynamic_cast<CsMCDigit*>(dg)->getHits();
	  list<CsMCHit*>::iterator ih;
	  for( ih=padMCHits.begin(); ih!=padMCHits.end(); ih++ ) {
	    tH = (*ih)->getDTime();
	    xH = (*ih)->getX();
	    yH = (*ih)->getY();
	  }
	  std::cout << "getClusFromPads  "<< padMCHits.size();
	  std::cout << " clus " << xc << "  " << yc;
	  std::cout << " hits " << xH << "  " << yH << std::endl;
*/
	}
	else {
          xc = (*ip)->ix() + 0.5;
          yc = (*ip)->iy() + 0.5;
	  xc = (xc*padx - hCatx) + dets->vOffCatW( ic ).x();
          yc = (yc*pady - hCaty) + dets->vOffCatW( ic ).y();
          PH = (*ip)->PH();
	}
	//std::cout << "getClusFromPads : " << ic << "  " << cat->isPMT()
	//	    << std::endl;
        lClusters_.push_back( new CsRCCluster( kClu, lPads,
					       ic, xc, yc, PH ) );

	xh = double( lPads.size() );
	if( hist.hRC1001 ) hist.hRC1001->Fill( xh );
//hh                       ------------------------
      }

  }


//==========================================================================
  void CsRCEventClusters::doClustering() {
//----------------------------------------


//--- Paolo Schiavon.
//--- Version 1.0
//--- Date  August 2005


//--- SCAN all the cathodes separately :
//    ----------------------------------
      CsRCDetectors *dets = CsRCDetectors::Instance();
      int nCathode = dets->nCathode();
      flag_ = false;

      for( int jCat=0; jCat < nCathode; jCat++ ) {

	CsRCCathode* cat = dets->ptrToCat( jCat );

	if( cat->isPMT() ) getClusFromPads( jCat );
//                         -----------------------
	//else if( cat->isAPV() ) doClusteringAPV( jCat );
	else if( cat->isAPV() ) doClustering( jCat );
//                              --------------------
	else doClustering( jCat );
//           --------------------
      }
      if( !lClusters_.empty() ) flag_ = true;      //   ???

  }


//==========================================================================
  void CsRCEventClusters::doClustering( int jCat ) {
//--------------------------------------------------


//--- Paolo Schiavon.
//--- Version 0.1.3 - October 2000
//--- date  6/98
//--- for corrected x-y pad wires.
//--- from cluproxy version 3/a
//--- to C++ 22/6/99.
//--- in progress.

//--------------------------------------------------
//    Finds pad clusters from pad hits and ADC's,
//    for each cathode of the RICH photon chambers.
//--------------------------------------------------

      CsRCExeKeys *key = CsRCExeKeys::Instance();

      CsRCDetectors *dets = CsRCDetectors::Instance();
      list<CsRCCathode*> lCathodes = dets->lCathodes();

      static bool firstCall = true;
      if( firstCall ) {
        firstCall = false;
        key->acknoMethod( "CsRCEventClusters::doClustering" );
      }

      CsRCRecConst *cons = CsRCRecConst::Instance();
      int kOffHi = cons->kOffHi();
      float maxCluSize = cons->maxCluSize();

      CsRCHistos& hist = CsRCHistos::Ref();
      int kh;
      double xh, yh, wh;

      CsRCCathode* cat = dets->ptrToCat( jCat );
      double padx = cat->padx();
      double pady = cat->pady();
      int nPadx = cat->nPadx();
      int nPady = cat->nPady();
      double hCatx = cat->hCatx();
      double hCaty = cat->hCaty();

      float PHratio11 = 0.05;
//@@-------------------------
      float PHratio00x = 0.3;
      float PHratio00y = 0.6;
//@@-------------------------

      list<CsRCPad*> lPadAll = CsRCEventPads::Instance()->lPads();
      list<CsRCPad*>::iterator ip;
      int nPadEv = lPadAll.size();

//--- SCAN all the cathodes separately :
//    ----------------------------------

//----- note : jCat from 0 to 15 /  ixPad,iyPad from 0 to 71.
//      ----------------------------------------------------
	int kPadCat = 0;
	for( ip=lPadAll.begin(); ip!=lPadAll.end(); ip++ ) {
	  if( !(*ip)->flag() ) continue;
	  if( (*ip)->ic() == jCat ) kPadCat++;
	}
	int nPadCat = kPadCat;
        if( nPadCat == 0 ) return;

	int ixPad[nPadCat];
        int iyPad[nPadCat];
        float PHPad[nPadCat];
        CsRCPad *ptToPad[nPadCat];

        kPadCat = 0;
	for( ip=lPadAll.begin(); ip!=lPadAll.end(); ip++ ) {
	  if( !(*ip)->flag() ) continue;
	  if( (*ip)->ic() == jCat )  {
	    //ixPad[kPadCat] = (*ip)->ix();
	    //iyPad[kPadCat] = (*ip)->iy();
	    //PHPad[kPadCat] = (*ip)->PH();
	    ptToPad[kPadCat] = (*ip);
            kPadCat++;
	  }
	}

//----- ordering of pads by ascending ix (first) and iy (second)
//      --------------------------------------------------------
//----- added 020528 -----
	//cout << endl << jCat << "  " << nPadCat << endl;
        for( int kPad=0; kPad<nPadCat; kPad++ ) {
	  CsRCPad* pad = ptToPad[kPad];
	  int kord = pad->ix()*nPadx + pad->iy();
	  for( int jPad=0; jPad<kPad; jPad++ ) {
	    CsRCPad* pad = ptToPad[jPad];
	    int jord = pad->ix()*nPadx + pad->iy();
	    if( kord < jord ) {
	      CsRCPad* wptToPad = ptToPad[kPad];
	      for( int iPad=kPad; iPad>jPad; iPad-- ) {
		ptToPad[iPad] = ptToPad[iPad-1];
	      }
	      ptToPad[jPad] = wptToPad;
	      break;
	    }
	  }
	}
        for( int kPad=0; kPad<nPadCat; kPad++ ) {
	  CsRCPad* pad = ptToPad[kPad];
	  ixPad[kPad] = pad->ix();
	  iyPad[kPad] = pad->iy();
	  //PHPad[kPad] = pad->PH();
	  PHPad[kPad] = pad->PH() + 1.;                   //   020823

	  //cout << jCat << "  " << ixPad[kPad] << "  "
	  //     << iyPad[kPad] << "  " << PHPad[kPad] << endl;
	}

	bool flagPad[nPadCat];
        bool flagPadS[nPadCat];
        int ptToPadS[nPadCat];
        int kPadClu[nPadCat];

        int kCluCat = 0;
        float PHratio = 0.;

        for( int kPad=0; kPad<nPadCat; kPad++ ) flagPad[kPad] = true;

//----- loop :
//      ------
        int kPmax = 0;
        while ( kPmax >= 0 ) {

//------- find pad of max PH :
//        --------------------
          float PHmax = 0.;
          kPmax = -1;
          for( int kPad=0; kPad<nPadCat; kPad++ ) {
            if( flagPad[kPad] ) {
              if( PHPad[kPad] > PHmax ) {
		PHmax = PHPad[kPad];
		kPmax = kPad;
	      }
	    }
          }
          if( kPmax >= 0 ) {

	    int ixPmax = ixPad[kPmax];
            int iyPmax = iyPad[kPmax];
            flagPad[kPmax] = false;

//--------- look for 'connected' pads around maximum :
//          ------------------------------------------
            int kPadS = 0;
            ptToPadS[kPadS] = kPmax;
            flagPadS[kPadS] = true;
//--------- loop :
//          ------
            int iPad1 =  0;
            int iPad2 = -1;
            while ( kPadS > iPad2 ) {
              iPad1 = iPad2 + 1;
              iPad2 = kPadS;

              for( int kPa=iPad1; kPa<=iPad2; kPa++ ) {
		int kPmaxS = ptToPadS[kPa];
                int ixPmaxS = ixPad[kPmaxS];
                int iyPmaxS = iyPad[kPmaxS];

//------------- select pads along x or y only :
//              -------------------------------
//              take into account that pad index ix runs first, 
//              iy runs for fixed ix.
 
                for( int isgs=-1; isgs<=1; isgs+=2 ) {
		  int ix = ixPmaxS;
                  int iy = iyPmaxS + isgs;
                  int kc1 = kPmaxS + isgs;
                  int kc2 = 0;
                  if( isgs == 1 ) kc2 = nPadCat-1;
                  int kc = kc1;
                  if( kc >= 0  &&  kc < nPadCat ) {

                    if( ixPad[kc] == ix  &&  iyPad[kc] == iy )  {
		      if( flagPad[kc] ) {
		        flagPad[kc] = false;
                        kPadS++;
                        ptToPadS[kPadS] = kc;
                        flagPadS[kPadS] = true;
		      }
		    }
                    ix = ixPmaxS + isgs;
                    iy = iyPmaxS;
                    kc = kc1;
                    int iddc = isgs*(kc2-kc);
                    float ddxx = isgs*(ix-ixPad[kc]);
                    while ( iddc >= 0  &&  ddxx >= 0 )  {
		      if( ixPad[kc] == ix  &&  iyPad[kc] == iy ) {
			if( flagPad[kc] ) {
			  flagPad[kc] = false;
                          kPadS++;
                          ptToPadS[kPadS] = kc;
                          flagPadS[kPadS] = true;
			}
		      }
                      kc += isgs;
                      iddc = isgs*(kc2-kc);
                      if( iddc >= 0 ) { ddxx = isgs*(ix-ixPad[kc]); }
		    }
                  }   /* end if kc */

                }   /* end loop: isgs */

              }   /* end loop: kPa */

            }   /* end while: kPadS > iPad2 */

            int nPadS = ++kPadS;

            int offClu = 50;
            int ptToPadC[nPadS*offClu];

//--------- construct cluster around max :
//          ------------------------------
//          store first pad in first cluster :
//          ----------------------------------
            int kClu = 0;
            kPadS = 0;
            int kPoint  = ptToPadS[kPadS];
            int kSto = kClu * offClu + kPadS;
            ptToPadC[kSto] = kPoint;
            kPadClu[kClu] = 1;
            flagPadS[kPadS] = false;
            if( kPadS >= offClu ) {
              string str = 
	   "RICHONE, CsRCEventClusters::doClustering() : buffer overflow (1)";
              CsErrLog::Instance()->mes( elError, str );
	    }

//--------- loop :
//          ------
            int kPmaxB = 0;
            while ( kPmaxB >= 0 )  {

//----------- find secondary maximum :
//            ------------------------
              float PHmaxB = 0.;
              kPmaxB = -1;
              for( int kPadB=0; kPadB<nPadS; kPadB++ )  {
                if( flagPadS[kPadB] ) {
		  kPoint = ptToPadS[kPadB];
		  if( PHPad[kPoint] > PHmaxB ) {
		    kPmaxB = kPadB;
		    PHmaxB = PHPad[kPoint];
		  }
		}
              }
              if( kPmaxB >= 0 )  {
		flagPadS[kPmaxB] = false;
                kPoint  = ptToPadS[kPmaxB];
                int ixPmaxB = ixPad[kPoint];
                int iyPmaxB = iyPad[kPoint];

		bool bNewClu = false;
                int kCu = 0;
                while( kCu <= kClu  &&  !bNewClu ) {

//--------------- check if secondary maximum is a new cluster :
//                ---------------------------------------------
                  kSto = kCu * offClu;
                  kPoint = ptToPadC[kSto];
                  PHratio = PHmaxB / PHPad[kPoint];
		  //PHratio = 0.;                              //   !!!
                  int iddx = abs( ixPmaxB - ixPad[kPoint] );
                  int iddy = abs( iyPmaxB - iyPad[kPoint] );
//@@-------------------------------------------------------
                  if( ( iddx == 0 && iddy == 1 && PHratio > PHratio00y )
                      || ( iddy == 0 && iddx == 1 && PHratio > PHratio00x ) ) {
//@@--------------------------------------------------------
//----------------- store it :
//                  ----------
		    kClu++;
                    kSto = kClu * offClu;
                    ptToPadC[kSto] = ptToPadS[kPmaxB];
                    kPadClu[kClu] = 1;
		    bNewClu = true;
		  }
//--------------- check if secondary maximum adds to cluster :
//                --------------------------------------------
		  if( !bNewClu ) {
		    int iPa = kPadClu[kCu];
                    for( int kPa=0; kPa<iPa; kPa++ ) {
		      kSto = kCu * offClu + kPa;
                      kPoint = ptToPadC[kSto];
                      int iddx = abs( ixPmaxB - ixPad[kPoint] );
                      int iddy = abs( iyPmaxB - iyPad[kPoint] );
//@@-----------------------------------------------------------
		      if( iddx + iddy <= 1 ) {
//@@-----------------------------------------
			bNewClu = true;
			if( iPa >= offClu ) {
                          string str = 
	   "RICHONE, CsRCEventClusters::doClustering() : buffer overflow (2)";
                          CsErrLog::Instance()->mes( elError, str );
		          break;
		        }
//--------------------- add to cluster :
//                      ----------------
		        //kSto = kCu * offClu + iPa;
                        //ptToPadC[kSto] = ptToPadS[kPmaxB];
                        //kPadClu[kCu] += 1;
                        //break;
			kPadClu[kCu] += 1;               //   030922
//--------------------- check clu size :
//                      ----------------
			if( kPadClu[kCu] <= maxCluSize ) {
			  kSto = kCu * offClu + iPa;
			  ptToPadC[kSto] = ptToPadS[kPmaxB];
			  break;
//--------------------- otherwise return to pad pool :
//                      ------------------------------
			} else {
			  flagPad[ptToPadS[kPmaxB]] = true;
			  kPadClu[kCu] -= 1;
			  break;
			}
		      }
		    }
		  }   /* end if on bNewClu */

		  kCu++;
		}   /* end while on kCu... */

//------------- if none, store a new cluster :
//              ------------------------------
		//if( !bNewClu ) {
		//  kClu++;
                //  kSto = kClu * offClu;
                //  ptToPadC[kSto] = ptToPadS[kPmaxB];
                //  kPadClu[kClu] = 1;
		//}
//------------- if none, return to pad pool :               //   030922
//              -----------------------------
		if( !bNewClu ) flagPad[ptToPadS[kPmaxB]] = true;

              }   /* end if on kPmaxB >= 0 */

	    }   /* end while on kPmaxB >= 0 */


//--------- compute cluster cms :
//          ---------------------
	    for( int kCu=0; kCu<=kClu; kCu++ ) {

              int nPadClu = kPadClu[kCu];
              int kCluS = kCluCat + kCu;
              float xCluw = 0.;
              float yCluw = 0.;
              float PHCluw = 0.;
	      list<CsRCPad*> lPads;
              for( int kPa=0; kPa<nPadClu; kPa++ ) {
                kSto = kCu * offClu + kPa;
                kPoint = ptToPadC[kSto];
                xCluw += ixPad[kPoint] * PHPad[kPoint];
                yCluw += iyPad[kPoint] * PHPad[kPoint];
                PHCluw += PHPad[kPoint];
        	lPads.push_back( ptToPad[kPoint] );
              }
	      if( PHCluw <= 0. ) continue;              //   011128
              xCluw /= PHCluw;
              yCluw /= PHCluw;
              xCluw += 0.5;
              yCluw += 0.5;
//----------- cluster to MWR :
//            ----------------
	      xCluw = (xCluw*padx - hCatx) + dets->vOffCatW( jCat ).x();
	      yCluw = (yCluw*pady - hCaty) + dets->vOffCatW( jCat ).y();

	      //^PHCluw -= float(nPadClu);                //   020823
              PHCluw /= float(nPadClu);                //   061130

//----------- construct cluster object :
//            --------------------------
	      int jClu = lClusters_.size();
	      lClusters_.push_back( new CsRCCluster( jClu, lPads, jCat,
						     xCluw, yCluw, PHCluw ) );

              int kSto0 = kCu * offClu;
              int kPoint = ptToPadC[kSto0];
              for( int kCh=0; kCh<nPadClu; kCh++ ) {
                kSto = kSto0 + kCh;
                int kPoins = ptToPadC[kSto];
                float ddxx = ixPad[kPoins] - ixPad[kPoint];
                float ddyy = iyPad[kPoins] - iyPad[kPoint];
                xh = ddxx;
		yh = ddyy;
		if( hist.hRC1011 ) hist.hRC1011->Fill( xh, yh );
//hh                               ----------------------------
                //wh = PHPad[kPoins];
                wh = PHPad[kPoins] - 1.;                //   020823
		if( hist.hRC1021 ) hist.hRC1021->Fill( xh, yh, wh );
//hh                               --------------------------------
              }
              xh = nPadClu;
	      if( hist.hRC1001 ) hist.hRC1001->Fill( xh );
//hh                             ------------------------

	    }   /* end loop on kCu */

	    kCluCat += ++kClu;

	  }   /* end if on kPmax >= 0 */

	}   /* end of loop while on kPmax >= 0 */

	int nCluCat = kCluCat;

  }


//==========================================================================
  void CsRCEventClusters::doClusteringAPV( int jCat ) {
//-----------------------------------------------------
//- Pro memoria
  }


//==========================================================================  
  void CsRCEventClusters::killHaloClus() {
//----------------------------------------


//- Paolo  -  February  2003


    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;
      CsRCExeKeys* key = CsRCExeKeys::Instance();
      key->acknoMethod( "CsRCEventClusters::killHaloClus" );
    }

    CsRCHistos& hist = CsRCHistos::Ref();
    double xh, yh, wh;

    CsRCRecConst* cons = CsRCRecConst::Instance();
    double radDeg = cons->RadDeg();
    static float CFRefInd = cons->CFRefInd();

    double xPade =   0.;
    double yPade0 = 470.;
    double yPade  =   2.;
    double RR = 167.;
    double DeRR = 20.;

    int nClusters = lClusters_.size();
    double xClu[nClusters];
    double yClu[nClusters];
    double errClu[nClusters];
    //double sigma = 5.;
    double sigma = DeRR/3.46;

    list<CsRCCluster*>::iterator ic;
    int kClu = 0;
    for( ic=lClusters_.begin(); ic!=lClusters_.end(); ic++ ) {
      if( !(*ic)->flag() ) continue;

      double xCluw = (*ic)->xc();
      double yCluw = (*ic)->yc();
      if( yCluw >= 0. ) yPade =  yPade0;
      if( yCluw <= 0. ) yPade = -yPade0;
      xCluw -= xPade;
      yCluw -= yPade;
      double ddXXYY = sqrt( xCluw*xCluw + yCluw*yCluw );

      if( fabs( ddXXYY - RR ) < DeRR ) {
	xClu[kClu] = xCluw;
	yClu[kClu] = yCluw;
	errClu[kClu] = sigma*sigma;
	kClu++;
      }
    }
    int nClu = kClu;

//- fit :
//  -----
    bool exeFit = true;
    if( nClu < 5 ) exeFit = false;
//@@------------------------------------------
    if( exeFit ) {

      static vector<CsHist2D*> vRC3770;
      static int kOffHi = cons->kOffHi();
      static bool fitHaloDef = true;
      if( fitHaloDef ) {
	fitHaloDef = false;

	for( int j=0; j<10; j++) vRC3770.push_back( NULL );
	int level = 1;
        bool bk = CsRCHistos::Ref().bookHis()  &&
	  level <= CsRCHistos::Ref().levelBk();
        if( bk ) {
	  string hTitle;
	  vRC3770.clear();

          stringstream hN3771;
          hN3771 << kOffHi + 3771;
          hTitle = "c0fit - c0in";
          vRC3770.push_back( new CsHist2D( hN3771.str(), hTitle,
					   100, -10., 10., 100, -10., 10. ) );
          stringstream hN3772; 
          hN3772 << kOffHi + 3772;
          hTitle = "RRfit - RRin  vs  nDF";
          vRC3770.push_back( new CsHist2D( hN3772.str(), hTitle,
					   100, -10., 10., 50, 0., 50. ) );
	  stringstream hN3773;
	  hN3773 << kOffHi + 3773;
	  hTitle = "chisq / nu  vs  nDF";
	  vRC3770.push_back( new CsHist2D( hN3773.str(), hTitle,
					   100, 0., 10., 50, 0., 50. ) );
	  stringstream hN3774;
	  hN3774 << kOffHi + 3774;
	  hTitle = "n iter  vs  nDF";
	  vRC3770.push_back( new CsHist2D( hN3774.str(), hTitle,
					   20, 0., 20., 50, 0., 50. ) );
	  stringstream hN3775;
	  hN3775 << kOffHi + 3775;
	  hTitle = "pulls  vs  nDF";
          vRC3770.push_back( new CsHist2D( hN3775.str(), hTitle,
					   100, -5., 5., 50, 0., 50. ) );
	  stringstream hN3776;
	  hN3776 << kOffHi + 3776;
	  hTitle = "nClu  vs  nCluHa";
          vRC3770.push_back( new CsHist2D( hN3776.str(), hTitle,
					   100, 0., 100., 100, 0., 100. ) );
	  stringstream hN3777;
	  hN3777 << kOffHi + 3777;
	  hTitle = "yCluHa  vs  xCluHa";
          vRC3770.push_back( new CsHist2D( hN3777.str(), hTitle,
			     100, -200., 200., 100, -200., 200. ) );
	}
      }

      int nParam = 3;
      double param[nParam];
      param[0] = RR;
      param[1] = 0.;
      param[2] = 0.;
      int iPaFit[nParam];
      iPaFit[0] = 0;
      iPaFit[1] = 0;
      iPaFit[2] = 0;

      CsRCCircleFit oCirclev( nClu, xClu, yClu, errClu,
//    ---------------------------------------------------
			      nParam, param, iPaFit );
      //cout << "====== " << endl;
      if( oCirclev.doChiFit() ) {
//        -------------------

	oCirclev.doHist( vRC3770 );
//      --------------------------
	//cout << "------ " << oCirclev.para()[0] << "  " << RR << endl;
	//cout << oCirclev.para()[1] << "  " << oCirclev.para()[2] << endl;
	double xCircle = oCirclev.para()[1] + xPade;
	double yCircle = 0.;
	//double rrCircle = 167.;
	double rrCircle = oCirclev.para()[0];
	double DeCir = 4.;
//@@---------------------
	int kCluHa = 0;
	for( ic=lClusters_.begin(); ic!=lClusters_.end(); ic++ ) {
	  if( !(*ic)->flag() ) continue;

	  double xCluw = (*ic)->xc();
	  double yCluw = (*ic)->yc();
	  if( yCluw >= 0. ) yCircle = oCirclev.para()[2] + yPade0;
	  if( yCluw <= 0. ) yCircle = oCirclev.para()[2] - yPade0;
	  xCluw -= xCircle;
	  yCluw -= yCircle;
	  double ddXXYY = sqrt( xCluw*xCluw + yCluw*yCluw );
	  //cout << ddXXYY << "  ";
	  if( fabs( ddXXYY - rrCircle ) < DeCir ) {
	    (*ic)->setFlag( false );
	    kCluHa++;
	    xh = xCluw;
	    yh = yCluw;
	    if( vRC3770[6] ) vRC3770[6]->Fill( xh, yh );
//hh                         --------------------------
	  }
	}
	int nCluHa = kCluHa;
	xh = nClu;
	yh = nCluHa;
	if( vRC3770[5] ) vRC3770[5]->Fill( xh, yh );
//hh                     --------------------------
	//cout << endl;
	//cout << nClu << "  " << nCluHa << endl;
      }

    }

  }


//=========================================================================== 
  void CsRCEventClusters::checkClusters() {
//-----------------------------------------


//- Paolo - August  2003

//- WARNING : to be adapted to PMT !!!


    CsRCRecConst *cons = CsRCRecConst::Instance();
    int kOffHi = cons->kOffHi();
    float maxCluSize = cons->maxCluSize();

    CsRCDetectors *dets = CsRCDetectors::Instance();

    CsRCExeKeys *key = CsRCExeKeys::Instance();
    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      key->acknoMethod( "CsRCEventClusters::checkClusters" );
    }

    std::list<CsRCCluster*>::iterator ic;
    for( ic=lClusters_.begin(); ic!=lClusters_.end(); ic++ ) {
      if( !flag() ) continue;

      //if( (*ic)->isAPV() ) continue;
      //if( (*ic)->isPMT() ) continue;

      CsRCCathode* cath = dets->ptrToCat( (*ic)->ic() );
      double padx = cath->padx();
      double pady = cath->pady();
      int nPadx = cath->nPadx();
      int nPady = cath->nPady();
      double hCatx = padx * nPadx / 2.;
      double hCaty = pady * nPady / 2.;

      int kCat = (*ic)->ic();
      double xClu = (*ic)->xc();
      double yClu = (*ic)->yc();
      xClu = int( ( xClu - dets->vOffCatW( kCat ).x() + hCatx ) / padx );
      yClu = int( ( yClu - dets->vOffCatW( kCat ).y() + hCaty ) / pady );
      //std::cout << xClu << "  " << yClu << std::endl;
      std::list<CsRCPad*> lPads = (*ic)->lPads();
      std::list<CsRCPad*>::iterator ip;

      //if( lPads.size() < 5 ) continue;
      //if( !(lPads.size() == maxCluSize) ) continue;
      //if( (lPads.size() <= maxCluSize) ) continue;
//@@---------------------------------------------

      CsRCHistos& hist = CsRCHistos::Ref();
      int kh;
      double xh, yh, wh;
      static int kHHMx = 20;
      if( !CsRCHistos::Ref().bookHis() ||
	  CsRCHistos::Ref().levelBk() < 2 ) kHHMx = 0;
      static vector<CsHist2D*> vRC6500;
      static int kPlotP = 0;
      static bool bPlotP = true;
//@@---------------------------
      if( bPlotP  &&  kPlotP < kHHMx ) {
        CsHistograms::SetCurrentPath("/RICH");
        int khPlots = kOffHi + 6500 + kPlotP;
        int nCanPl = 25;
        float xCanPl = float( nCanPl );
        stringstream name0, title0;
        name0 << khPlots;
        title0 << "   ";
        vRC6500.push_back( new CsHist2D( name0.str(), title0.str(),
//hh    -----------------------------------------------------------
	  			         2*nCanPl, -xCanPl, xCanPl,
				         2*nCanPl, -xCanPl, xCanPl ) );

	CsHistograms::SetCurrentPath("");

	for( ip=lPads.begin(); ip!=lPads.end(); ip++ ) {
	  if( !flag() ) continue;

          xh = (*ip)->ix() - xClu;
	  yh = (*ip)->iy() - yClu;
          //wh = (*ip)->PH();
	  wh = 1.;
	  if( vRC6500[kPlotP] ) vRC6500[kPlotP]->Fill( xh, yh, wh );
//hh                            -----------------------------------
	}

      }
      kPlotP++;

    }

  }


//=========================================================================== 
  void CsRCEventClusters::checkPMTClus() {
//----------------------------------------


//- Paolo
//  August  2006


    if( lClusters_.size() == 0 ) return;

    CsRCExeKeys *key = CsRCExeKeys::Instance();
    CsRCDetectors *dets = CsRCDetectors::Instance();
    CsRCHistos& hist = CsRCHistos::Ref();
    double xh, yh, wh;
    CsRCRecConst *cons = CsRCRecConst::Instance();
    int kOffHi = cons->kOffHi();

    static std::vector<CsHist2D*> vRC8200;
    static int kHH = 0;
    static int kHHMx = 20;
    if( !CsRCHistos::Ref().bookHis() ||
	CsRCHistos::Ref().levelBk() < 2 ) kHHMx = 0;
    static int nCathode = 0;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;
      key->acknoMethod( "CsRCEventClusters::checkPMTClus()" );

      nCathode = dets->nCathode();
    }


    if( kHH >= kHHMx ) return;
    if( kHH < kHHMx ) {
      int kHist = kOffHi + 8200 + kHH;
      stringstream hN8200;
      hN8200 << kHist;
      string hTitle = "Clusters";
      CsHistograms::SetCurrentPath("/RICH");
      vRC8200.push_back( new CsHist2D( hN8200.str(), hTitle,
//hh  ------------------------------------------------------
			 268, -1608., 1608., 268, -1608., 1608. ) );
      CsHistograms::SetCurrentPath("/");

      list<CsRCCluster*>::iterator ic;
      for( ic=lClusters_.begin(); ic!=lClusters_.end(); ic++ ) {

	CsRCCathode* cath = dets->ptrToCat( (*ic)->ic() );
	//if( !cath->isPMT() ) continue;
	int nPPadx = cath->nPPadx();
	int nPPady = cath->nPPady();
	xh = (*ic)->xc();
	yh = (*ic)->yc();
	vRC8200[kHH]->Fill( xh, yh );
//      ----------------------------
      }
      kHH++;
    }

    return;
  }
