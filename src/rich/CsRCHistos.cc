
/*!
   \file    CsRCHistos.cc
   \brief   Histogram booking for CsRichOne.
   \author  Paolo Schiavon
   \version 0.01
   \date    $Date: 2011/02/04 11:10:12 $
*/


#include <iostream>
#include <ostream>
#include <cstdio>
#include <cstdlib>

#include <cmath>

// ---------------------------
#include "CsOpt.h"

// ---------------------------
#include "CsRCHistos.h"

#include "CsRCExeKeys.h"
#include "CsRCRecConst.h"
#include "CsRCDetectors.h"
#include "CsRCMirrors.h"
 // ---------------------------

  using namespace std;


  CsRCHistos* CsRCHistos::Ptr_ = 0;

//===========================================================================
  CsRCHistos& CsRCHistos::Ref() {
//-------------------------------
    if( Ptr_ == 0 ) Ptr_ = new CsRCHistos();
    return *Ptr_; 
  }

//===========================================================================
  void CsRCHistos::print() const {
//--------------------------------
    cout << endl;
    cout << "CsRCHistos pt = " << this << "  " << hRC3765 << endl;
    cout << endl;
  }

//===========================================================================
  CsRCHistos::~CsRCHistos() {
//---------------------------
  }


//===========================================================================
  CsRCHistos::CsRCHistos() {
//--------------------------


//--- July  2000, rev. January 2001


//- HBOOK Histogram booking :
//  -------------------------
//  Histo.s booked in ascending order
//  ---------------------------------

    static bool MCEvent = false;
    if( CsRCExeKeys::Instance()->MCarloEvent() ) MCEvent = true;
    static bool DataEvent = false;
    if( CsRCExeKeys::Instance()->DataEvent() ) DataEvent = true;
//- for MC run book all Data run histo.s (provisional)
//  --------------------------------------------------
    if( MCEvent ) DataEvent = true;
//@@------------------------------

    CsOpt* opt = CsOpt::Instance();
    bool boo;
    int kOpw;
    bool lprint = false;
    boo = opt->CsOpt::getOpt( "RICHONE", "PrintKeys", kOpw );
    if( boo ) if( kOpw > 0 ) lprint = true;

//- Histogram booking swhich
//  ------------------------
    string sOpw;
    static bool bookHis = false;
    boo = opt->CsOpt::getOpt( "RICHONE", "BookHistograms", sOpw );
    if( boo ) if( sOpw == "YES" ) bookHis = true;
    //bookHis_ = bookHis;
    int iOpw;
//- Histogram booking level: 0, 1, 2(=all)
//  --------------------------------------
    static int levelBk = 0;
    boo = opt->CsOpt::getOpt( "RICHONE", "BookHistoLevel", iOpw );
    if( boo ) levelBk = iOpw;
    //levelBk_ = levelBk;


//- Centralized histogram booking swhich and level
//  ----------------------------------------------
    string hLevStr = "none";
    if( CsOpt::Instance()->getOpt("RI01P", "hist level", hLevStr ) ||
	CsOpt::Instance()->getOpt("RI", "hist level", hLevStr ) ) {
      if( hLevStr == "none" ) bookHis = false;
      else if( hLevStr == "normal" ) { bookHis = true; levelBk = 0; }
      else if( hLevStr == "high" ) { bookHis = true; levelBk = 2; }
      else  bookHis = false;
      
    //} else  bookHis = false;
    }

    if( bookHis )
      cout << "RICHONE, CsRCHistos()  INFO: Histogram level for \"" << "RI01P"
	//<< "\" is " << hLevStr << " (" << levelBk+1 << ").\n";
	   << "\" is " << hLevStr << endl;

    bookHis_ = bookHis;
    levelBk_ = levelBk;

    //if( bookHis ) CsHistograms::SetCurrentPath("/RICH");
    CsHistograms::SetCurrentPath("/RICH");
//@@-------------------------------------

    vector<float> vflo;
    float momMin = 0.;
    float momMax = 0.;
    boo = opt->CsOpt::getOpt( "RICHONE", "MomentumRange", vflo );
    if( boo ) {
      momMin = vflo[0];
      momMax = vflo[1];
    }
    float momAcc = 220;
    if( momMax < 100. ) momAcc = 120.;
    int mBi = int( momMax - momMin );
    int mBA = int( momAcc );

    CsRCRecConst *cons = CsRCRecConst::Instance();
    static int kOffHi = cons->kOffHi();
    CsRCDetectors* dets = CsRCDetectors::Instance();
    static int nCathode = dets->nCathode();

    int k;
    string hTitle;
    bool bk = false;
    int level;
    int kHisBk = 0;


    hRC1000 = NULL;
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1000;
      hN1000 << kOffHi + 1000;
      hTitle = "rejections";
      hRC1000 = new CsHist1D( hN1000.str(), hTitle, 100, 0., 100. );
      kHisBk++;
    }

//- CsRCEventClusters::doClustering():
    hRC1001 = NULL;
    hRC1011 = NULL;
    hRC1021 = NULL;
    hRC1022 = NULL;
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1001;
      hN1001 << kOffHi + 1001;
      hTitle = "cluster size";
      hRC1001 = new CsHist1D( hN1001.str(), hTitle, 50, 0., 50. );
      kHisBk++;
    }
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1011;
      hN1011 << kOffHi + 1011;
      hTitle = "cluster shape around pad of max PH";
      hRC1011 = new CsHist2D( hN1011.str(), hTitle,
			      9, -4.5, 4.5, 9, -4.5, 4.5 );
      kHisBk++;
      stringstream hN1021;
      hN1021 << kOffHi + 1021;
      hTitle = "cluster PH shape around pad of max PH";
      hRC1021 = new CsHist2D( hN1021.str(), hTitle,
			      9, -4.5, 4.5, 9, -4.5, 4.5 );
      kHisBk++;
      stringstream hN1022;
      hN1022 << kOffHi + 1022;
      hTitle = "cluster PH shape around pad of max PH";
      hRC1022 = new CsHist2D( hN1022.str(), hTitle,
			      100, -50., 50., 100, -50., 50. );
      kHisBk++;
    }

//- from CsRichOne::printMemory():
    hRC1031 = NULL;
    level = 2;
    bk = bookHis && level<=levelBk;
    if( bk ) {
      stringstream hN1031;
      hN1031 << kOffHi + 1031;
      hTitle = "HEAP vs time";
      hRC1031 = new CsHist2D( hN1031.str(), hTitle,
			      50, 0., 50., 10, 0., 10. );
      kHisBk++;
    }

//- from CsRCEventParticles::getEventParticles():
    hRC1061 = NULL;
    hRC1062 = NULL;
    hRC1063 = NULL;
    hRC1066 = NULL;
    hRC1068 = NULL;
    level = 1;
    bk = bookHis && MCEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1061;
      hN1061 << kOffHi + 1061;
      hTitle = "thetaCerAv vs thetaCer - thetaCerAv";
      hRC1061 = new CsHist2D( hN1061.str(), hTitle,
			      100, -5., 5., 60, 0., 60. );
      kHisBk++;
      stringstream hN1062;
      hN1062 << kOffHi + 1062;
      hTitle = "Ephot vs thetaCerAv";
      hRC1062 = new CsHist2D( hN1062.str(), hTitle,
			      60, 0., 60., 50, 5.5, 8.0 );
      kHisBk++;
      stringstream hN1063;
      hN1063 << kOffHi + 1063;
      hTitle = "Ephot vs thetaCer - thetaCerAv";
      hRC1063 = new CsHist2D( hN1063.str(), hTitle,
			      100, -5., 5., 50, 5.5, 8.0 );
      kHisBk++;
      stringstream hN1066;
      hN1066 << kOffHi + 1066;
      hTitle = "Ephot vs nPhotons";
      hRC1066 = new CsHist2D( hN1066.str(), hTitle,
			      60, 0., 60., 50, 5.5, 8.0 );
      kHisBk++;
      stringstream hN1068;
      hN1068 << kOffHi + 1068;
      hTitle = "beta vs nPhotons";
      hRC1068 = new CsHist2D( hN1068.str(), hTitle,
			      60, 0., 60., 64, 0.9984, 1. );
      kHisBk++;
    }

//- booking in CsRCPartPhotons::getPartPhotons() test histo
//      CsHist2D* vRC1100(1-4);
//      CsHist2D* hRC1111;

//- CsRCEventPartPhotons::getEventPartPhotons():
    hRC1121 = NULL;
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1121;
      hN1121 << kOffHi + 1121;
      hTitle = "ParPhotons / event";
      hRC1121 = new CsHist1D( hN1121.str(), hTitle, 50, 0., 50. );
      kHisBk++;
    }

//- CsRCPartPhotons::getPartPhotons():
    hRC1131 = NULL;
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1131;
      hN1131 << kOffHi + 1131;
      hTitle = "Photons / PaPho";
      hRC1131 = new CsHist1D( hN1131.str(), hTitle, 500, 0., 500. );
      kHisBk++;
    }

//- CsRCEventParticles::setMCParticle()/setMCRecParticle()/setDataParticle():
    hRC1200 = NULL;
    hRC1201 = NULL;
    hRC1202 = NULL;
    hRC1203 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1200;
      hN1200 << kOffHi + 1200;
      hTitle = "y- vs x-track at PT centre";
      hRC1200 = new CsHist2D( hN1200.str(), hTitle,
			      100, -1000., 1000., 100, -1000., 1000. );
      kHisBk++;
      stringstream hN1201;
      hN1201 << kOffHi + 1201;
      hTitle = "y- vs x-track, entrance";
      hRC1201 = new CsHist2D( hN1201.str(), hTitle,
			      100, -1000., 1000., 100, -1000., 1000. );
      kHisBk++;
      stringstream hN1202;
      hN1202 << kOffHi + 1202;
      hTitle = "y- vs x-track, exit";
      hRC1202 = new CsHist2D( hN1202.str(), hTitle,
			      100, -1000., 1000., 100, -1000., 1000. );
      kHisBk++;
      stringstream hN1203;
      hN1203 << kOffHi + 1203;
      hTitle = "chsq-track vs SM";
      hRC1203 = new CsHist2D( hN1203.str(), hTitle,
			      100, 0., 500., 10, 0., 10. );
      kHisBk++;
    }
    hRC1204 = NULL;
    hRC1205 = NULL;
    hRC1206 = NULL;
    hRC1207 = NULL;
    hRC1208 = NULL;
    hRC1209 = NULL;
    level = 1;
    bk = bookHis && MCEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1204;
      hN1204 << kOffHi + 1204;
      hTitle = "Y- vs X-track pull, helix.0";
      hRC1204 = new CsHist2D( hN1204.str(), hTitle,
			      100, -5., 5., 100, -5., 5. );
      kHisBk++;
      stringstream hN1205;
      hN1205 << kOffHi + 1205;
      hTitle = "tgb- vs tga-track pull, helix.0";
      hRC1205 = new CsHist2D( hN1205.str(), hTitle,
			      100, -5., 5., 100, -5., 5. );
      kHisBk++;
      stringstream hN1206;
      hN1206 << kOffHi + 1206;
      hTitle = "X- vs Y-track pull, entrance w.";
      hRC1206 = new CsHist2D( hN1206.str(), hTitle,
			      100, -5., 5., 100, -5., 5. );
      kHisBk++;
      stringstream hN1207;
      hN1207 << kOffHi + 1207;
      hTitle = "tgb- vs tga-track pull, entrance w";
      hRC1207 = new CsHist2D( hN1207.str(), hTitle,
			      100, -5., 5., 100, -5., 5. );
      kHisBk++;
      stringstream hN1208;
      hN1208 << kOffHi + 1208;
      hTitle = "tgb- vs tga-track, helix.0";
      hRC1208 = new CsHist2D( hN1208.str(), hTitle,
			      100, -1., 1., 100, -1., 1. );
      kHisBk++;
      stringstream hN1209;
      hN1209 << kOffHi + 1209;
      hTitle = "tgb- vs tga-track, entrance w";
      hRC1209 = new CsHist2D( hN1209.str(), hTitle,
			      100, -1., 1., 100, -1., 1. );
      kHisBk++;
    }
//- CsRCEventParticles::setDataParticle():
    hRC1220 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1220;
      hN1220 << kOffHi + 1220;
      hTitle = "track mom. vs time";
      hRC1220 = new CsHist2D( hN1220.str(), hTitle,
      			      100, -100., 100., mBi, momMin, momMax );
      kHisBk++;
    }

//- booking in CsRCPartPhotons::checkLikeRatios() :
//     hRC1301-3, 5, 9-10

//- CsRCPartPhotons::nPhotExpctCat() :
    for( int kh=0; kh<nCathode; kh++ ) vRC1350.push_back( NULL );
    for( int kh=0; kh<nCathode; kh++ ) vRC1370.push_back( NULL );
    level = 2;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      vRC1350.clear();
      for( int kh=0; kh<nCathode; kh++ ){
	stringstream hN1350;
	hN1350 << kOffHi + 1350 + kh + 1;
	hTitle = "n-Photon Expected";
	vRC1350.push_back( new CsHist2D( hN1350.str(), hTitle, 100, 0., 100.,
					 6, 0., 6. ) );
	kHisBk++;
      }
      vRC1370.clear();
      for( int kh=0; kh<nCathode; kh++ ){
	stringstream hN1370;
	hN1370 << kOffHi + 1370 + kh + 1;
	hTitle = "n-Photon Expected, in vs out";
	vRC1370.push_back( new CsHist2D( hN1370.str(), hTitle, 50, 0., 100.,
					 50, 0., 100. ) );
	kHisBk++;
      }
    }

//- CsRCEventPads::checkPPads():
    for( int kh=0; kh<4; kh++ ) vRC1500.push_back( NULL );
    level = 2;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      vRC1500.clear();
      stringstream hN1501;
      hN1501 << kOffHi + 1501;
      hTitle = "norm of signal";
      vRC1500.push_back( new CsHist2D( hN1501.str(), hTitle, 48, 0., 48.,
				       48, 0., 48. ) );
      kHisBk++;
      stringstream hN1502;
      hN1502 << kOffHi + 1502;
      hTitle = "norm of signal";
      vRC1500.push_back( new CsHist2D( hN1502.str(), hTitle, 48, 0., 48.,
				       48, 0., 48. ) );
      kHisBk++;
      stringstream hN1503;
      hN1503 << kOffHi + 1503;
      hTitle = "norm of signal";
      vRC1500.push_back( new CsHist2D( hN1503.str(), hTitle, 48, 0., 48.,
				       48, 0., 48. ) );
      kHisBk++;
      stringstream hN1504;
      hN1504 << kOffHi + 1504;
      hTitle = "norm of signal";
      vRC1500.push_back( new CsHist2D( hN1504.str(), hTitle, 48, 0., 48.,
				       48, 0., 48. ) );
      kHisBk++;
    }

//- CsRCPartPhotons::nPhotExpct
    hRC1505 = NULL;
    hRC1506 = NULL;
    hRC1507 = NULL;
    hRC1508 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1505;
      hN1505 << kOffHi + 1505;
      hTitle = "n-phot expect";
      hRC1505 = new CsHist2D( hN1505.str(), hTitle, 60, -1500., 1500.,
			      60, -1500., 1500. );
      kHisBk++;
      stringstream hN1506;
      hN1506 << kOffHi + 1506;
      hTitle = "n-phot expect";
      hRC1506 = new CsHist2D( hN1506.str(), hTitle, 60, -1500., 1500.,
			      60, -1500., 1500. );
      kHisBk++;
      stringstream hN1507;
      hN1507 << kOffHi + 1507;
      hTitle = "fecUV";
      hRC1507 = new CsHist1D( hN1507.str(), hTitle, 50, 0., 1. );
      kHisBk++;
      stringstream hN1508;
      hN1508 << kOffHi + 1508;
      hTitle = "frcVS";
      hRC1508 = new CsHist1D( hN1508.str(), hTitle, 50, 0., 1. );
      kHisBk++;
    }

//- CsRCLikeAll05::normSignal, normBackgr(), likeSignal, likeBackgr():
    hRC1510 = NULL;
    hRC1511 = NULL;
    hRC1512 = NULL;
    hRC1513 = NULL;
    hRC1514 = NULL;
    hRC1515 = NULL;
    hRC1516 = NULL;
    hRC1517 = NULL;
    hRC1518 = NULL;
    hRC1519 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1510;
      hN1510 << kOffHi + 1510;
      hTitle = "norm of signal";
      hRC1510 = new CsHist2D( hN1510.str(), hTitle, 100, 0., 100.,
			      22, 0., 1.1 );
      kHisBk++;
      stringstream hN1511;
      hN1511 << kOffHi + 1511;
      hTitle = "norm of backgr";
      hRC1511 = new CsHist2D( hN1511.str(), hTitle, 100, 0., 100.,
			      100, 0., 100. );
      kHisBk++;
      stringstream hN1512;
      hN1512 << kOffHi + 1512;
      hTitle = "prob of signal";
      hRC1512 = new CsHist2D( hN1512.str(), hTitle, 100, 0., 100.,
			      11, 0., 1.1 );
      kHisBk++;
      stringstream hN1513;
      hN1513 << kOffHi + 1513;
      hTitle = "norm of split rings";
      hRC1513 = new CsHist2D( hN1513.str(), hTitle, 110, 0., 1.1,
			      10, 0., 10. );
      kHisBk++;
      stringstream hN1514;
      hN1514 << kOffHi + 1514;
      hTitle = "fracUse";
      hRC1514 = new CsHist2D( hN1514.str(), hTitle, 100, 0., 200.,
			      100, 0., 400. );
      kHisBk++;
      stringstream hN1515;
      hN1515 << kOffHi + 1515;
      hTitle = "fracUse";
      hRC1515 = new CsHist2D( hN1515.str(), hTitle, 100, 0., 200.,
			      100, 0., 400. );
      kHisBk++;
      stringstream hN1516;
      hN1516 << kOffHi + 1516;
      hTitle = "probSig";
      hRC1516 = new CsHist2D( hN1516.str(), hTitle, 100, 0., 200.,
			      100, 0., 400. );
      kHisBk++;
      stringstream hN1517;
      hN1517 << kOffHi + 1517;
      hTitle = "probSig";
      hRC1517 = new CsHist2D( hN1517.str(), hTitle, 100, 0., 200.,
			      100, 0., 400. );
      kHisBk++;
      stringstream hN1518;
      hN1518 << kOffHi + 1518;
      hTitle = "intBack";
      hRC1518 = new CsHist2D( hN1518.str(), hTitle, 100, 0., 500.,
			      100, 0., 100. );
      kHisBk++;
      stringstream hN1519;
      hN1519 << kOffHi + 1519;
      hTitle = "intBack";
      hRC1519 = new CsHist2D( hN1519.str(), hTitle, 100, 0., 500.,
			      100, 0., 100. );
      kHisBk++;
    }

//- CsRCEventPads::getEventPads():
    hRC1520 = NULL;
    hRC1521 = NULL;
    hRC1522 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1520;
      hN1520 << kOffHi + 1520;
      hTitle = "pads / event";
      hRC1520 = new CsHist2D( hN1520.str(), hTitle, 320, -160., 160.,
			      320, -160., 160. );
      //hRC1520 = new CsHist2D( hN1520.str(), hTitle, 300, -150., 150.,
      //			300, -150., 150. );
      kHisBk++;
    }
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1521;
      hN1521 << kOffHi + 1521;
      hTitle = "pads / event";
      hRC1521 = new CsHist1D( hN1521.str(), hTitle, 200, 0., 4000. );
      kHisBk++;
    }
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1522;
      hN1522 << kOffHi + 1522;
      hTitle = "pads / event vs cathode";
      hRC1522 = new CsHist2D( hN1522.str(), hTitle, 100, 0., 500.,
			      nCathode, 0., double( nCathode ) );
      kHisBk++;
    }

//- CsRCEventRings::checkNoise():
    hRC1523 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1523;
      hN1523 << kOffHi + 1523;
      hTitle = "noise pads / event vs cathode";
      hRC1523 = new CsHist2D( hN1523.str(), hTitle, 100, 0., 500.,
			      nCathode, 0., double( nCathode ) );
      kHisBk++;
    }

//- CsRCEventPads::getEventPads():
    hRC1524 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1524;
      hN1524 << kOffHi + 1524;
      hTitle = "pad PH vs cathode";
      hRC1524 = new CsHist2D( hN1524.str(), hTitle, 50, 0., 50.,
			      nCathode, 0., double( nCathode ) );
      kHisBk++;
    }

//- CsRCEventClusters::doClustering():
    hRC1525 = NULL;
    hRC1526 = NULL;
    hRC1527 = NULL;
    hRC1528 = NULL;
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1525;
      hN1525 << kOffHi + 1525;
      hTitle = "y-clu vs x-clu on detectors";
      hRC1525 = new CsHist2D( hN1525.str(), hTitle,
			      200, -1600., 1600., 200, -1600., 1600. );
      kHisBk++;
    }
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1526;
      hN1526 << kOffHi + 1526;
      hTitle = "clusters / event";
      hRC1526 = new CsHist1D( hN1526.str(), hTitle, 200, 0., 4000. );
      kHisBk++;
    }
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1527;
      hN1527 << kOffHi + 1527;
      hTitle = "pad PH vs cathode";
      hRC1527 = new CsHist2D( hN1527.str(), hTitle, 50, 0., 1000.,
			      nCathode, 0., double( nCathode ) );
      kHisBk++;
      stringstream hN1528;
      hN1528 << kOffHi + 1528;
      hTitle = "clusters / pad";
      hRC1528 = new CsHist1D( hN1528.str(), hTitle, 60, 0.5, 1.1 );
      kHisBk++;
    }

//- CsRCEventRings::checkRingPH():
    hRC1529 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1529;
      hN1529 << kOffHi + 1529;
      hTitle = "ring pad PH vs angle";
      hRC1529 = new CsHist2D( hN1529.str(), hTitle, 50, 0., 50.,
			      60, 0., 60. );
      kHisBk++;
    }
//- CsRCLikeRing0X::likeBackgr() and CsRCCluster::getBackWgt():
    hRC1530 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1530;
      hN1530 << kOffHi + 1530;
      hTitle = "weight vs x-clu, y-clu on detectors";
      hRC1530 = new CsHist2D( hN1530.str(), hTitle,
			      200, -1600., 1600., 200, -1600., 1600. );
      kHisBk++;
    }

//- CsRCEventParticles::getEventParticles():
    hRC1531 = NULL;
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1531;
      hN1531 << kOffHi + 1531;
      hTitle = "particles / event - input";
      hRC1531 = new CsHist1D( hN1531.str(), hTitle, 50, 0., 50. );
      kHisBk++;
    }

//- CsRCEventParticles::partAnalysis():
    hRC1533 = NULL;
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1533;
      hN1533 << kOffHi + 1533;
      hTitle = "particles / event - filter";
      hRC1533 = new CsHist1D( hN1533.str(), hTitle, 50, 0., 50. );
      kHisBk++;
    }

//- CsRCEventRings::getEventRings():
    hRC1535 = NULL;
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1535;
      hN1535 << kOffHi + 1535;
      hTitle = "rings / event - reconstr.";
      hRC1535 = new CsHist1D( hN1535.str(), hTitle, 50, 0., 50. );
      kHisBk++;
    }
    hRC1536 = NULL;
    level = 2;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1536;
      hN1536 << kOffHi + 1536;
      hTitle = "ring time diff.          ";
      hRC1536 = new CsHist1D( hN1536.str(), hTitle, 100, -25., 25. );
      kHisBk++;
    }

//- CsRCEventAnalysis::dataMonitor() 
    hRC1537 = NULL;
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1537;
      hN1537 << kOffHi + 1537;
      hTitle = "rings / event - accepted";
      hRC1537 = new CsHist1D( hN1537.str(), hTitle, 50, 0., 50. );
      kHisBk++;
    }

//- CsRCEventAnalysis::partIdent() 
    hRC1539 = NULL;
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1539;
      hN1539 << kOffHi + 1539;
      hTitle = "rings / event - identified";
      hRC1539 = new CsHist1D( hN1539.str(), hTitle, 50, 0., 50. );
      kHisBk++;
    }

//- CsRCLikeRing0X::likeBackgr() and CsRCCluster::getBackWgt():
    hRC1540 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1540;
      hN1540 << kOffHi + 1540;
      hTitle = "counts vs x-clu, y-clu on detectors";
      hRC1540 = new CsHist2D( hN1540.str(), hTitle,
			      200, -1600., 1600., 200, -1600., 1600. );
      kHisBk++;
    }

//- CsRCRing::getRingPk():
    hRC1545 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1545;
      hN1545 << kOffHi + 1545;
      hTitle = "photons in peak, Peak";
      hRC1545 = new CsHist1D( hN1545.str(), hTitle, 100, 0., 100. );
      kHisBk++;
    }
//- CsRCRing::getRingBf()
    hRC1547 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1547;
      hN1547 << kOffHi + 1547;
      hTitle = "photons in peak, BackFilter";
      hRC1547 = new CsHist1D( hN1547.str(), hTitle, 100, 0., 100. );
      kHisBk++;
    }
//- CsRCRing::getRingWa()
    hRC1548 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1548;
      hN1548 << kOffHi + 1548;
      hTitle = "photons in peak, 3Sigma";
      hRC1548 = new CsHist1D( hN1548.str(), hTitle, 100, 0., 100. );
      kHisBk++;
    }
//- CsRCRing::getRingBf():
    hRC1551 = NULL;
     level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1551;
      hN1551 << kOffHi + 1551;
      hTitle = "peak window vs betarec";
      hRC1551 = new CsHist2D( hN1551.str(), hTitle,
			      20, 0., 20., 64, 0.9984, 1. );
      kHisBk++;
    }

//- CsRCEventPartPhotons::moniCFRefInd():
    hRC1556 = NULL;
    hRC1557 = NULL;
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1556;
      hN1556 << kOffHi + 1556;
      hTitle = "n-1 (VS)";
      hRC1556 = new CsHist2D( hN1556.str(), hTitle, 200, 0.0002, 0.0022,
			      100, 0., 100. );
                              //100, 70., 170. );
      kHisBk++;
      stringstream hN1557;
      hN1557 << kOffHi + 1557;
      hTitle = "n-1 (UV)";
      hRC1557 = new CsHist2D( hN1557.str(), hTitle, 200, 0.0002, 0.0022,
			      100, 0., 100. );
                              //100, 70., 170. );
      kHisBk++;
    }
    hRC1558 = NULL;
    hRC1559 = NULL;
    level = -1;      //      book always
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1558;
      hN1558 << kOffHi + 1558;
      hTitle = "n-1 (VS)";
      hRC1558 = new CsHist1D( hN1558.str(), hTitle, 200, 0.0002, 0.0022 );
      kHisBk++;
      stringstream hN1559;
      hN1559 << kOffHi + 1559;
      hTitle = "n-1 (UV)";
      hRC1559 = new CsHist1D( hN1559.str(), hTitle, 200, 0.0002, 0.0022 );
      kHisBk++;
    }

//- CsRCEventAnalysis::partIdent():
    hRC1568 = NULL;
    hRC1569 = NULL;
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1568;
      hN1568 << kOffHi + 1568;
      hTitle = "chi-square";
      hRC1568 = new CsHist1D( hN1568.str(), hTitle, 100, 0., 10. );
      kHisBk++;
    }
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1569;
      hN1569 << kOffHi + 1569;
      hTitle = "chi-square vs n photons";
      hRC1569 = new CsHist2D( hN1569.str(), hTitle,
			      100, 0., 100., 24, 6., 30. );
      kHisBk++;
    }

//- CsRCPartPhotons::getLikeAll(), nPhotExpct(): 
    hRC1571 = NULL;
    hRC1572 = NULL;
    hRC1573 = NULL;
//- CsRCLikeAll05::normSignal(), normBackgr(), likeSignal(), likeBackgr()
    hRC1574 = NULL;
    hRC1575 = NULL;
    hRC1576 = NULL;
    hRC1577 = NULL;
    hRC1578 = NULL;
    hRC1579 = NULL;
    level = 2;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1571;
      hN1571 << kOffHi + 1571;
      hTitle = "nPhotons-nExpect";
      hRC1571 = new CsHist1D( hN1571.str(), hTitle, 200, -100., 100. );
      kHisBk++;
      stringstream hN1572;
      hN1572 << kOffHi + 1572;
      hTitle = "nPhoExpect PMT";
      hRC1572 = new CsHist1D( hN1572.str(), hTitle, 100, 0., 100. );
      kHisBk++;
      stringstream hN1573;
      hN1573 << kOffHi + 1573;
      hTitle = "nPhoExpect APV";
      hRC1573 = new CsHist1D( hN1573.str(), hTitle, 100, 0., 100. );
      kHisBk++;
      stringstream hN1574;
      hN1574 << kOffHi + 1574;
      hTitle = "signal on PMT.s";
      hRC1574 = new CsHist1D( hN1574.str(), hTitle, 200, 0., 1. );
      kHisBk++;
      stringstream hN1575;
      hN1575 << kOffHi + 1575;
      hTitle = "signal on APV.s";
      hRC1575 = new CsHist1D( hN1575.str(), hTitle, 200, 0., 1. );
      kHisBk++;
      stringstream hN1576;
      hN1576 << kOffHi + 1576;
      hTitle = "backgr on PMT.s";
      hRC1576 = new CsHist1D( hN1576.str(), hTitle, 200, 0., 1. );
      kHisBk++;
      stringstream hN1577;
      hN1577 << kOffHi + 1577;
      hTitle = "normSignal";
      hRC1577 = new CsHist1D( hN1577.str(), hTitle, 200, 0., 200. );
      kHisBk++;
      stringstream hN1578;
      hN1578 << kOffHi + 1578;
      hTitle = "normBackgr";
      hRC1578 = new CsHist1D( hN1578.str(), hTitle, 200, 0., 200. );
      kHisBk++;
      stringstream hN1579;
      hN1579 << kOffHi + 1579;
      hTitle = "backgr on APV";
      hRC1579 = new CsHist1D( hN1579.str(), hTitle, 200, 0., 1. );
      kHisBk++;
    }

//- CsRCRing::histCntRatio():
    hRC1581 = NULL;
    hRC1582 = NULL;
    hRC1583 = NULL;
    hRC1584 = NULL;
    hRC1585 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1581;
      hN1581 << kOffHi + 1581;
      hTitle = "PH peak max";
      hRC1581 = new CsHist1D( hN1581.str(), hTitle, 100, 0., 100. );  //phtw=1
      kHisBk++;
      stringstream hN1582;
      hN1582 << kOffHi + 1582;
      hTitle = "PH peak average";
      hRC1582 = new CsHist1D( hN1582.str(), hTitle, 100, 0., 100. );  //phtw=1
      kHisBk++;
      stringstream hN1583;
      hN1583 << kOffHi + 1583;
      hTitle = "PH peak average/max";
      hRC1583 = new CsHist1D( hN1583.str(), hTitle, 100, 0., 2. );
      kHisBk++;
      stringstream hN1584;
      hN1584 << kOffHi + 1584;
      hTitle = "peak shape";
      hRC1584 = new CsHist1D( hN1584.str(), hTitle, 100, 0., 100. );
      kHisBk++;
      stringstream hN1585;
      hN1585 << kOffHi + 1585;
      hTitle = "PH peak total";
      hRC1585 = new CsHist1D( hN1585.str(), hTitle, 100, 0., 100. );
      kHisBk++;
    }

//- CsRCEventAnalysis::partIdent() (cont.):
    for( int kh=0; kh<5; kh++ ) vRC1590.push_back( NULL );
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      vRC1590.clear();
      stringstream hN1591;
      hN1591 << kOffHi + 1591;
      hTitle = "reconstructed mass vs likelihood ratio - elctron";
      vRC1590.push_back( new CsHist2D( hN1591.str(), hTitle, 
				       500, 0.0, 1.0, 100, 1., 6. ) );
      stringstream hN1592;
      hN1592 << kOffHi + 1592;
      hTitle = "reconstructed mass vs likelihood ratio - muon";
      vRC1590.push_back( new CsHist2D( hN1592.str(), hTitle, 
				       500, 0.0, 1.0, 100, 1., 6. ) );
      stringstream hN1593;
      hN1593 << kOffHi + 1593;
      hTitle = "reconstructed mass vs likelihood ratio - pion";
      vRC1590.push_back( new CsHist2D( hN1593.str(), hTitle, 
				       500, 0.0, 1.0, 100, 1., 6. ) );
      stringstream hN1594;
      hN1594 << kOffHi + 1594;
      hTitle = "reconstructed mass vs likelihood ratio - kaon";
      vRC1590.push_back( new CsHist2D( hN1594.str(), hTitle, 
				       500, 0.0, 1.0, 100, 1., 6. ) );
      stringstream hN1595;
      hN1595 << kOffHi + 1595;
      hTitle = "reconstructed mass vs likelihood ratio - proton";
      vRC1590.push_back( new CsHist2D( hN1595.str(), hTitle, 
				       500, 0.0, 1.0, 100, 1., 6. ) );
    }

//- CsRCEventPartPhotons::moniCFRefInd() :
    for( int kh=0; kh<nCathode; kh++ ) vRC1600.push_back( NULL );
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      vRC1600.clear();
      for( int kh=0; kh<nCathode; kh++ ) {
	stringstream hN1600;
	hN1600 << kOffHi + 1600 + kh;
	hTitle = "n-1 (VS)";
	vRC1600.push_back( new CsHist1D( hN1600.str(), hTitle,
					 200, 0.0002, 0.0022 ) );
	kHisBk++;
      }
    }

//- CsRCLikeAll05::normSignal, likeSignal:
    hRC1623 = NULL;
    hRC1624 = NULL;
    hRC1625 = NULL;
    hRC1626 = NULL;
    hRC1627 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1623;
      hN1623 << kOffHi + 1623;
      hTitle = "ring ang. vs cath.";
      hRC1623 = new CsHist2D( hN1623.str(), hTitle, 90, 0., 360.,
			      10, 0., 10. );
      kHisBk++;
      stringstream hN1624;
      hN1624 << kOffHi + 1624;
      hTitle = "the0 vs rr0, fracUse = 0";
      hRC1624 = new CsHist2D( hN1624.str(), hTitle, 100, 0., 200.,
			      100, 0., 400. );
      kHisBk++;
      stringstream hN1625;
      hN1625 << kOffHi + 1625;
      hTitle = "the0 vs rr0, fracUse = 1";
      hRC1625 = new CsHist2D( hN1625.str(), hTitle, 100, 0., 200.,
			      100, 0., 400. );
      kHisBk++;
      stringstream hN1626;
      hN1626 << kOffHi + 1626;
      hTitle = "the0 vs rr0, fracUse = 0";
      hRC1626 = new CsHist2D( hN1626.str(), hTitle, 100, 0., 200.,
			      100, 0., 400. );
      kHisBk++;
      stringstream hN1627;
      hN1627 << kOffHi + 1627;
      hTitle = "the0 vs rr0, fracUse = 1";
      hRC1627 = new CsHist2D( hN1627.str(), hTitle, 100, 0., 200.,
			      100, 0., 400. );
      kHisBk++;
    }

//- CsRCEventAnalysis::partIdent() (cont.):
//- CsRCEventRings::partChiSqIdent() : 1630, 1631
//- CsRCRing:: getThetaWgav() : 1632
//- CsRCEventRings::partRingAllIdent() : 1633
    hRC1621 = NULL;
    hRC1630 = NULL;
    hRC1631 = NULL;
    hRC1632 = NULL;
    hRC1633 = NULL;
    hRC1634 = NULL;
    hRC1635 = NULL;
    hRC1636 = NULL;
    hRC1637 = NULL;
    hRC1638 = NULL;
    hRC1639 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1621;
      hN1621 << kOffHi + 1621;
      hTitle = "n-1 Ipo";
      hRC1621 = new CsHist1D( hN1621.str(), hTitle, 200, 1.3, 1.7 );
      kHisBk++;
      stringstream hN1630;
      hN1630 << kOffHi + 1630;
      hTitle = "theta rec - theta Fit";
      hRC1630 = new CsHist1D( hN1630.str(), hTitle, 200, -2., 2. );
      kHisBk++;
      stringstream hN1631;
      hN1631 << kOffHi + 1631;
      hTitle = "theta rec - theta Fit";
      hRC1631 = new CsHist1D( hN1631.str(), hTitle, 200, -2., 2. );
      kHisBk++;
      stringstream hN1632;
      hN1632 << kOffHi + 1632;
      hTitle = "pulls theta WgAv";
      hRC1632 = new CsHist1D( hN1632.str(), hTitle, 100, -5., 5. );
      kHisBk++;
      stringstream hN1633;
      hN1633 << kOffHi + 1633;
      hTitle = "theta rec - theta Like";
      hRC1633 = new CsHist1D( hN1633.str(), hTitle, 200, -2., 2. );
      kHisBk++;
      //stringstream hN1634;
      //hN1634 << kOffHi + 1634;
      //hTitle = "n-Photons, corr back";
      //hRC1634 = new CsHist1D( hN1634.str(), hTitle, 100, 0., 100. );
      //kHisBk++;
      stringstream hN1635;
      hN1635 << kOffHi + 1635;
      hTitle = "theta rec - theta Ide";
      hRC1635 = new CsHist2D( hN1635.str(), hTitle, 100, -2.5, 2.5,
			      5, 0., 5. );
      kHisBk++;
    }
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1634;
      hN1634 << kOffHi + 1634;
      hTitle = "theta MC UV - theta Ide  vs theta Ide";
      hRC1634 = new CsHist2D( hN1634.str(), hTitle, 100, -0.25, 0.25,
			      54, 4., 58.);
      kHisBk++;
      stringstream hN1636;
      hN1636 << kOffHi + 1636;
      hTitle = "theta rec - theta Ide ([b] 0.99995-1)";
      hRC1636 = new CsHist1D( hN1636.str(), hTitle, 100, -2.5, 2.5 );
      kHisBk++;
      stringstream hN1637;
      hN1637 << kOffHi + 1637;
      hTitle = "theta MC VS - theta Ide  vs theta Ide";
      hRC1637 = new CsHist2D( hN1637.str(), hTitle, 100, -0.25, 0.25,
			      54, 4., 58.);
      kHisBk++;
    }
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1638;
      hN1638 << kOffHi + 1638;
      hTitle = "nPhotons RingCorr - Expected";
      hRC1638 = new CsHist2D( hN1638.str(), hTitle, 50, -25., 25.,
			      5, 0., 5. );
      kHisBk++;
      stringstream hN1639;
      hN1639 << kOffHi + 1639;
      hTitle = "nZero";
      hRC1639 = new CsHist1D( hN1639.str(), hTitle, 100, 0., 50. );
      kHisBk++;
    }

    hRC1640 = NULL;
    hRC1641 = NULL;
    hRC1642 = NULL;
    hRC1643 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1640;
      hN1640 << kOffHi + 1640;
      hTitle = "TheReco vs momMeas";
      hRC1640 = new CsHist2D( hN1640.str(), hTitle, 
			      100, 0., 50., 70, 0., 70. );
      kHisBk++;
      stringstream hN1641;
      hN1641 << kOffHi + 1641;
      hTitle = "Pion Like vs momMeas";
      hRC1641 = new CsHist2D( hN1641.str(), hTitle, 
			      100, 0., 50., 100, 0., 5. );
      kHisBk++;
      stringstream hN1642;
      hN1642 << kOffHi + 1642;
      hTitle = "Kaon Like vs momMeas";
      hRC1642 = new CsHist2D( hN1642.str(), hTitle, 
			      100, 0., 50., 100, 0., 5. );
      kHisBk++;
      stringstream hN1643;
      hN1643 << kOffHi + 1643;
      hTitle = "Proton Like vs momMeas";
      hRC1643 = new CsHist2D( hN1643.str(), hTitle, 
			      100, 0., 50., 100, 0., 5. );
      kHisBk++;
    }

    for( int kh=0; kh<4; kh++ ) vRC1650.push_back( NULL );
    hRC1655 = NULL;
    //for( int kh=0; kh<4; kh++ ) vRC1655.push_back( NULL );
    level = 1;
    bk = bookHis && MCEvent && level<=levelBk;
    if( bk ) {
      vRC1650.clear();
      stringstream hN1651;
      hN1651 << kOffHi + 1651;
      hTitle = "pion mass";
      vRC1650.push_back( new CsHist1D( hN1651.str(), hTitle, 500, 0.0, 1.0 ) );
      kHisBk++;
      stringstream hN1652;
      hN1652 << kOffHi + 1652;
      hTitle = "Kaon mass";
      vRC1650.push_back( new CsHist1D( hN1652.str(), hTitle, 500, 0.0, 1.0 ) );
      kHisBk++;
      stringstream hN1653;
      hN1653 << kOffHi + 1653;
      hTitle = "proton mass";
      vRC1650.push_back( new CsHist1D( hN1653.str(), hTitle, 500, 0.0, 1.0 ) );
      kHisBk++;
      stringstream hN1655;
      hN1655 << kOffHi + 1655;
      hTitle = "mass reconst.d vs mass MC";
      hRC1655 = new CsHist2D( hN1655.str(), hTitle,
			      10, 0., 10., 10, 0., 10. );
    }

//- 020916
    for( int kh=0; kh<4; kh++ ) vRC1655.push_back( NULL );
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      kHisBk++;
      vRC1655.clear();
      stringstream hN1656;
      hN1656 << kOffHi + 1656;
      hTitle = "mom vs dmom/mom, pion";
      vRC1655.push_back( new CsHist2D( hN1656.str(), hTitle,
				     100, -0.25, 0.25, mBi, momMin, momMax ) );
      kHisBk++;
      stringstream hN1657;
      hN1657 << kOffHi + 1657;
      hTitle = "mom vs dmom/mom, Kaon";
      vRC1655.push_back( new CsHist2D( hN1657.str(), hTitle,
				     100, -0.25, 0.25, mBi, momMin, momMax ) );
      kHisBk++;
      stringstream hN1658;
      hN1658 << kOffHi + 1658;
      hTitle = "proton mass";
      vRC1655.push_back( new CsHist2D( hN1658.str(), hTitle,
				       100, -0.25, 0.25, 30, 0., 30. ) );
      kHisBk++;
    }

    hRC1661 = NULL;
    hRC1662 = NULL;
    hRC1663 = NULL;
    hRC1664 = NULL;
    hRC1665 = NULL;
    hRC1666 = NULL;
    hRC1667 = NULL;
    hRC1668 = NULL;
    hRC1669 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1661;
      hN1661 << kOffHi + 1661;
      hTitle = "chi-square minimum, reconstructed";
      hRC1661 = new CsHist2D( hN1661.str(), hTitle, 100, 0., 10., 5, 0., 5. );
      kHisBk++;
      stringstream hN1662;
      hN1662 << kOffHi + 1662;
      hTitle = "chi-square > minimum";
      hRC1662 = new CsHist1D( hN1662.str(), hTitle, 100, 0., 10. );
      kHisBk++;
      stringstream hN1663;
      hN1663 << kOffHi + 1663;
      hTitle = "chi-square - chi-square-minimum vs chi-square";
      hRC1663 = new CsHist2D( hN1663.str(), hTitle,
			    100, 0., 2.5, 100, 0., 10. );
      kHisBk++;
      stringstream hN1664;
      hN1664 << kOffHi + 1664;
      hTitle = "second to max likelihood ratio vs ident mass";
      hRC1664 = new CsHist2D( hN1664.str(), hTitle,
			      100, 0., 1., 5, 0., 5. );
      kHisBk++;
      stringstream hN1665;
      hN1665 << kOffHi + 1665;
      hTitle = "likelihood probability vs 2nd probability";
      hRC1665 = new CsHist2D( hN1665.str(), hTitle,
			      100, 0., 0.2, 100, 0., 0.2 );
      kHisBk++;
      stringstream hN1666;
      hN1666 << kOffHi + 1666;
      hTitle = "chi-square probability, max";
      hRC1666 = new CsHist1D( hN1666.str(), hTitle, 100, 0., 1. );
      kHisBk++;
      stringstream hN1667;
      hN1667 << kOffHi + 1667;
      hTitle = "chi-square probability, < max";
      hRC1667 = new CsHist1D( hN1667.str(), hTitle, 100, 0., 1. );
      kHisBk++;
      stringstream hN1668;
      hN1668 << kOffHi + 1668;
      hTitle = "chi-square probability, max-less";
      hRC1668 = new CsHist1D( hN1668.str(), hTitle, 100, 0., 0.1 );
      kHisBk++;
    }
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1669;
      hN1669 << kOffHi + 1669;
      hTitle = "likelihood probability vs part type";
    //hRC1669 = new CsHist2D( hN1669.str(), hTitle, 100, 0., 0.2, 5, 0., 5. );
      hRC1669 = new CsHist2D( hN1669.str(), hTitle, 200, 0., 1., 5, 0., 5. );
      kHisBk++;
    }

    //for( int kh=0; kh<12; kh++ ) vRC1670.push_back( NULL );
    for( int kh=0; kh<10; kh++ ) vRC1670.push_back( NULL );
    level = 1;
    bk = bookHis && MCEvent && level<=levelBk;
    if( bk ) {
      vRC1670.clear();
      //for( int kh=0; kh<12; kh++ ) {
      for( int kh=0; kh<10; kh++ ) {
        stringstream hN1670;
        hN1670 << kOffHi + 1670 + kh;
        hTitle = "mass reconst.d vs mass MC in mom. bins";
        vRC1670.push_back( new CsHist2D( hN1670.str(), hTitle,
					 5, 0., 5., 5, 0., 5. ) );
        kHisBk++;
      }
    }

//- CsRCEventAnalysis::partIdent() (cont.):
    for( int kh=0; kh<5; kh++ ) vRC1680.push_back( NULL );
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      vRC1680.clear();
      stringstream hN1681;
      hN1681 << kOffHi + 1681;
      hTitle = "theta rec - theta Hypo  vs theta rec";
      vRC1680.push_back( new CsHist2D( hN1681.str(), hTitle,
				       100, -2.5, 2.5, 54, 4., 58.) );
      kHisBk++;
      stringstream hN1682;
      hN1682 << kOffHi + 1682;
      hTitle = "theta rec - theta Hypo  vs theta rec";
      vRC1680.push_back( new CsHist2D( hN1682.str(), hTitle,
				       100, -2.5, 2.5, 54, 4., 58.) );
      kHisBk++;
      stringstream hN1683;
      hN1683 << kOffHi + 1683;
      hTitle = "theta rec - theta Hypo  vs theta rec";
      vRC1680.push_back( new CsHist2D( hN1683.str(), hTitle,
				       100, -2.5, 2.5, 54, 4., 58.) );
      kHisBk++;
      stringstream hN1684;
      hN1684 << kOffHi + 1684;
      hTitle = "theta rec - theta Hypo  vs theta rec";
      vRC1680.push_back( new CsHist2D( hN1684.str(), hTitle,
				       100, -2.5, 2.5, 54, 4., 58.) );
      kHisBk++;
      stringstream hN1685;
      hN1685 << kOffHi + 1685;
      hTitle = "theta rec - theta Hypo  vs theta rec";
      vRC1680.push_back( new CsHist2D( hN1685.str(), hTitle,
				       100, -2.5, 2.5, 54, 4., 58.) );
      kHisBk++;
    }

//- CsRCEventAnalysis::partIdent() (cont.):
    hRC1691 = NULL;
    hRC1692 = NULL;
    hRC1693 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1691;
      hN1691 << kOffHi + 1691;
      hTitle = "likelihood less than backgr vs mom";
      hRC1691 = new CsHist2D( hN1691.str(), hTitle,
			      100, 0., 0.2, 100, 0., 60. );
      kHisBk++;
    }
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1692;
      hN1692 << kOffHi + 1692;
      hTitle = "reconstructed mass vs mom";
      hRC1692 = new CsHist2D( hN1692.str(), hTitle, 500, 0.0, 1.0,
			      60, 0., 60. );
      kHisBk++;
    }
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1693;
      hN1693 << kOffHi + 1693;
      hTitle = "reconstructed mass vs part. probability";
      hRC1693 = new CsHist2D( hN1693.str(), hTitle, 500, 0., 1., 5, 0., 5. );
      kHisBk++;
    }
    hRC1694 = NULL;
    hRC1695 = NULL;
    hRC1696 = NULL;
    level = 1;
    bk = bookHis && MCEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1694;
      hN1694 << kOffHi + 1694;
      hTitle = "reconstructed mass vs theReco, no mu.s";
      hRC1694 = new CsHist2D( hN1694.str(), hTitle, 500, 0.0, 1.0,
			      50, 0., 60. );
      kHisBk++;
      stringstream hN1695;
      hN1695 << kOffHi + 1695;
      hTitle = "reconstructed mass vs part. probability";
      hRC1695 = new CsHist2D( hN1695.str(), hTitle, 500, 0., 1., 5, 0., 5. );
      kHisBk++;
      stringstream hN1696;
      hN1696 << kOffHi + 1696;
      hTitle = "reconstructed mass vs theReco, no mu.s";
      hRC1696 = new CsHist2D( hN1696.str(), hTitle, 500, 0.0, 1.0,
			      60, 50., 56. );
      kHisBk++;
    }

//- CsRCPartPhotons::sigmaPhoRec(), CsRCPhoton::sigmaPhoPid(),
//  CsRCRing::sigmaRing(), CsRCRing::getRingPk() :
    hRC1701 = NULL;
    hRC1702 = NULL;
    hRC1703 = NULL;
    hRC1705 = NULL;
    hRC1708 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN1701;
      hN1701 << kOffHi + 1701;
      hTitle = "sigma Photon for Rec vs momentum";
      hRC1701 = new CsHist2D( hN1701.str(), hTitle, 
			      60, 0., 60., 100, 0., 3. );
      kHisBk++;
      stringstream hN1702;
      hN1702 << kOffHi + 1702;
      hTitle = "sigma Photon for PID vs momentum";
      hRC1702 = new CsHist2D( hN1702.str(), hTitle, 
			      90, 0., 360., 100, 0., 3. );
      kHisBk++;
      stringstream hN1703;
      hN1703 << kOffHi + 1703;
      hTitle = "sigma Photon for PID vs phi Photon";
      hRC1703 = new CsHist2D( hN1703.str(), hTitle, 
			      60, 0., 60., 100, 0., 3. );
      kHisBk++;
      stringstream hN1705;
      hN1705 << kOffHi + 1705;
      hTitle = "sigma Ring vs beta-1";
      hRC1705 = new CsHist2D( hN1705.str(), hTitle, 
			      64, -0.0016, 0., 100, 0, 3. );
      kHisBk++;
      stringstream hN1708;
      hN1708 << kOffHi + 1708;
      hTitle = "nMore  vs  nCan";
      hRC1708 = new CsHist2D( hN1708.str(), hTitle, 
			      10, 0., 10., 10, 0., 10. );
      kHisBk++;
    }

//- booking moved to CsRCMirrors::doAliMirrAll():
//        vRC1800.push_back( ...
//        hRC1980 ...

//- CsRCPartPhotons::doPMTOptCorr():
    for( int kh=0; kh<4; kh++ ) vRC2100.push_back( NULL );
    for( int kh=0; kh<4; kh++ ) vRC2105.push_back( NULL );
    for( int kh=0; kh<4; kh++ ) vRC2110.push_back( NULL );
    for( int kh=0; kh<12; kh++ ) vRC2120.push_back( NULL );
    level = 2;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      vRC2100.clear();
      vRC2105.clear();
      vRC2110.clear();
      vRC2120.clear();
      for( int kh=0; kh<4; kh++ ) {
	stringstream hN2100;
	hN2100 << kOffHi + 2100 + kh + 1;
	hTitle = "xPMTcorr-xPMT vs tgX";
	vRC2100.push_back( new CsHist2D( hN2100.str(), hTitle,
					 //100, -0.2, 0.2, 100, -10., 10.) );
					 100, -0.5, 0.5, 100, -10., 10.) );
	kHisBk++;
      }
      for( int kh=0; kh<4; kh++ ) {
	stringstream hN2105;
	hN2105 << kOffHi + 2105 + kh + 1;
	//hTitle = "yPMTcorr-yPMT vs xPMTcorr-xPMT";
	hTitle = "yPMTcorr vs xPMTcorr";
	vRC2105.push_back( new CsHist2D( hN2105.str(), hTitle,
					 //100, -10., 10., 100, -10., 10.) );
					 200, -50., 50., 200, -50., 50.) );
	kHisBk++;
      }
      for( int kh=0; kh<4; kh++ ) {
	stringstream hN2110;
	hN2110 << kOffHi + 2110 + kh + 1;
	hTitle = "yPMTcorr-yPMT vs tgY";
	vRC2110.push_back( new CsHist2D( hN2110.str(), hTitle,
					 100, -0.5, 0.5, 100, -10., 10.) );
	kHisBk++;
      }
      for( int kh=0; kh<12; kh++ ) {
	stringstream hN2120;
	hN2120 << kOffHi + 2120 + kh + 1;
	hTitle = "PMT cha";
	vRC2120.push_back( new CsHist1D( hN2120.str(), hTitle, 16, 0., 16.) );
	kHisBk++;
      }
    }

//- CsRCEventAnalysis::MCMonitor():
    hRC3001 = NULL;
    hRC3002 = NULL;
    level = 1;
    bk = bookHis && MCEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3001;
      hN3001 << kOffHi + 3001;
      hTitle = "thphrec-thphcal, photon";
      hRC3001 = new CsHist1D( hN3001.str(), hTitle, 100, -10., 10. );
      kHisBk++;
      stringstream hN3002;
      hN3002 << kOffHi + 3002;
      hTitle = "thphrec-thphcal, ring";
      hRC3002 = new CsHist1D( hN3002.str(), hTitle, 100, -2.5, 2.5 );
      kHisBk++;
    }

//- CsRCEventAnalysis::dataMonitor():
    hRC3006 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3006;
      hN3006 << kOffHi + 3006;
      hTitle = "photon spectrum";
      hRC3006 = new CsHist1D( hN3006.str(), hTitle, 60, 0., 60. );
      kHisBk++;
    }

//- CsRCEventAnalysis::MCMonitor() (cont.):
    hRC3008 = NULL;
    hRC3009 = NULL;
    hRC3010 = NULL;
    level = 1;
    bk = bookHis && MCEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3008;
      hN3008 << kOffHi + 3008;
      hTitle = "thphcal vs thphrec-thphcal, accepted";
      hRC3008 = new CsHist2D( hN3008.str(), hTitle,
			      100, -5., 5., 60, 0., 60. );
      kHisBk++;
      stringstream hN3009;
      hN3009 << kOffHi + 3009;
      hTitle = "mom. vs thphrec-thphcal, accepted";
      hRC3009 = new CsHist2D( hN3009.str(), hTitle,
			      100, -5., 5., mBi, momMin, momMax );
      kHisBk++;
      stringstream hN3010;
      hN3010 << kOffHi + 3010;
      hTitle = "betacal vs thphrec-thphcal, accepted";
      hRC3010 = new CsHist2D( hN3010.str(), hTitle,
			      100, -5., 5., 80, 0.9984, 1. );
      kHisBk++;
    }
    hRC3011 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3011;
      hN3011 << kOffHi + 3011;
      hTitle = "n photons in ring vs betarec";
      hRC3011 = new CsHist2D( hN3011.str(), hTitle,
			      50, 0., 50., 80, 0.9984, 1. );
      kHisBk++;
    }
    hRC3012 = NULL;
    hRC3013 = NULL;
    hRC3014 = NULL;
    hRC3015 = NULL;
    hRC3016 = NULL;
    hRC3017 = NULL;
    hRC3019 = NULL;
    hRC3020 = NULL;
    hRC3022 = NULL;
    hRC3023 = NULL;
    hRC3025 = NULL;
    hRC3026 = NULL;
    hRC3028 = NULL;
    hRC3029 = NULL;
    hRC3031 = NULL;
    hRC3034 = NULL;
    level = 1;
    bk = bookHis && MCEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3012;
      hN3012 << kOffHi + 3012;
      hTitle = "n pho vs thphrec-thphcal, reconst.d";
      hRC3012 = new CsHist2D( hN3012.str(), hTitle,
			      100, -5., 5., 40, 0., 40. );
      kHisBk++;
      stringstream hN3013;
      hN3013 << kOffHi + 3013;
      hTitle = "betacal vs thphrec-thphcal, reconst.d";
      hRC3013 = new CsHist2D( hN3013.str(), hTitle,
			      100, -5., 5., 80, 0.9984, 1. );
      kHisBk++;
      stringstream hN3014;
      hN3014 << kOffHi + 3014;
      hTitle = "betarec vs betacal, accepted";
      hRC3014 = new CsHist2D( hN3014.str(), hTitle,
			      80, 0.9984, 1., 80, 0.9984, 1. );
      kHisBk++;
      stringstream hN3015;
      hN3015 << kOffHi + 3015;
      hTitle = "n photon vs betacal, reconst.d";
      hRC3015 = new CsHist1D( hN3015.str(), hTitle, 64, 0.9984, 1. );
      kHisBk++;
      stringstream hN3016;
      hN3016 << kOffHi + 3016;
      hTitle = "part. accepted vs betacal";
      hRC3016 = new CsHist1D( hN3016.str(), hTitle, 64, 0.9984, 1. );
      kHisBk++;
      stringstream hN3017;
      hN3017 << kOffHi + 3017;
      hTitle = "part. recont.d vs betacal";
      hRC3017 = new CsHist1D( hN3017.str(), hTitle, 64, 0.9984, 1. );
      kHisBk++;
      stringstream hN3019;
      hN3019 << kOffHi + 3019;
      hTitle = "cos * beta * n vs betacal";
      hRC3019 = new CsHist1D( hN3019.str(), hTitle, 64, 0.9984, 1. );
      kHisBk++;
      stringstream hN3020;
      hN3020 << kOffHi + 3020;
      hTitle = "thphrec vs betacal";
      hRC3020 = new CsHist1D( hN3020.str(), hTitle, 64, 0.9984, 1. );
      kHisBk++;
      stringstream hN3022;
      hN3022 << kOffHi + 3022;
      hTitle = "part. accepted vs mom.";
      hRC3022 = new CsHist1D( hN3022.str(), hTitle, mBi, momMin, momMax );
      kHisBk++;
      stringstream hN3023;
      hN3023 << kOffHi + 3023;
      hTitle = "part. reconst.d vs mom.";
      hRC3023 = new CsHist1D( hN3023.str(), hTitle, mBi, momMin, momMax );
      kHisBk++;
      stringstream hN3025;
      hN3025 << kOffHi + 3025;
      hTitle = "part. accepted vs n photons";
      hRC3025 = new CsHist1D( hN3025.str(), hTitle, 30, 0., 60. );
      kHisBk++;
      stringstream hN3026;
      hN3026 << kOffHi + 3026;
      hTitle = "part. reconst.d vs n photons";
      hRC3026 = new CsHist1D( hN3026.str(), hTitle, 30, 0., 60. );
      kHisBk++;
      stringstream hN3028;
      hN3028 << kOffHi + 3028;
      hTitle = "betacal vs reconst. cut, accepted";
      hRC3028 = new CsHist2D( hN3028.str(), hTitle,
			      60, 0., 6., 64, 0.9984, 1. );
      kHisBk++;
      stringstream hN3029;
      hN3029 << kOffHi + 3029;
      hTitle = "particle-type vs mom, rejected";
      hRC3029 = new CsHist2D( hN3029.str(), hTitle,
			      mBi, momMin, momMax, 15, 1., 15. );
      kHisBk++;
      stringstream hN3031;
      hN3031 << kOffHi + 3031;
      hTitle = "thphrec-thphcal, beta = 1";
      hRC3031 = new CsHist1D( hN3031.str(), hTitle, 100, -2.5, 2.5  );
      kHisBk++;

      stringstream hN3034;
      hN3034 << kOffHi + 3034;
      hTitle = "thetarec vs theta calc";
      hRC3034 = new CsHist2D( hN3034.str(), hTitle,
			      108, 4., 58., 108, 4., 58. );
      kHisBk++;
    }

//- CsRCEventAnalysis::dataMonitor() (cont.):
    hRC3036 = NULL;
    hRC3037 = NULL;
    hRC3042 = NULL;
    hRC3043 = NULL;
    hRC3045 = NULL;
    hRC3046 = NULL;
    hRC3047 = NULL;
    hRC3050 = NULL;
    hRC3051 = NULL;
    hRC3052 = NULL;
    hRC3053 = NULL;
    hRC3054 = NULL;
    hRC3055 = NULL;
    hRC3056 = NULL;
    hRC3057 = NULL;
    hRC3058 = NULL;
    hRC3059 = NULL;
    hRC3060 = NULL;
    hRC3061 = NULL;
    hRC3062 = NULL;
    hRC3063 = NULL;
    hRC3064 = NULL;
    hRC3065 = NULL;
    hRC3067 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3036;
      hN3036 << kOffHi + 3036;
      hTitle = "particles useful vs betarec";
      hRC3036 = new CsHist1D( hN3036.str(), hTitle, 64, 0.9984, 1. );
      kHisBk++;
      stringstream hN3037;
      hN3037 << kOffHi + 3037;
      hTitle = "particles reconst.d vs betarec";
      hRC3037 = new CsHist1D( hN3037.str(), hTitle, 64, 0.9984, 1. );
      kHisBk++;
    }
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3042;
      hN3042 << kOffHi + 3042;
      hTitle = "particles useful vs momentum";
      hRC3042 = new CsHist1D( hN3042.str(), hTitle, mBi, momMin, momMax );
      kHisBk++;
      stringstream hN3043;
      hN3043 << kOffHi + 3043;
      hTitle = "particles reconst.d vs momentum";
      hRC3043 = new CsHist1D( hN3043.str(), hTitle, mBi, momMin, momMax );
      kHisBk++;
    }
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3045;
      hN3045 << kOffHi + 3045;
      hTitle = "particles useful vs n photon";
      //hRC3045 = new CsHist1D( hN3045.str(), hTitle, 30, 0., 60. );
      hRC3045 = new CsHist1D( hN3045.str(), hTitle, 100, 0., 100. );
      kHisBk++;
    }
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3046;
      hN3046 << kOffHi + 3046;
      hTitle = "particles reconst.d vs n photon";
      //hRC3046 = new CsHist1D( hN3046.str(), hTitle, 30, 0., 60. );
      hRC3046 = new CsHist1D( hN3046.str(), hTitle, 100, 0., 100. );
      kHisBk++;
    }
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3047;
      hN3047 << kOffHi + 3047;
      hTitle = "reconst.d ring ringness";
      hRC3047 = new CsHist1D( hN3047.str(), hTitle, 100, 0., 20. );
      kHisBk++;
      stringstream hN3050;
      hN3050 << kOffHi + 3050;
      hTitle = "reconst.d cluster multipl.";
      hRC3050 = new CsHist1D( hN3050.str(), hTitle, 20, 0., 20. );
      kHisBk++;
      stringstream hN3051;
      hN3051 << kOffHi + 3051;
      hTitle = "thphrec-thphipo, beta = 1";
      hRC3051 = new CsHist1D( hN3051.str(), hTitle, 100, -2.5, 2.5 );
      kHisBk++;
      stringstream hN3052;
      hN3052 << kOffHi + 3052;
      hTitle = "thphrec, beta = 1";
      hRC3052 = new CsHist1D( hN3052.str(), hTitle, 140, 0., 70. );
      kHisBk++;
      stringstream hN3053;
      hN3053 << kOffHi + 3053;
      hTitle = "betacal vs thphrec-thphcalw";
      hRC3053 = new CsHist2D( hN3053.str(), hTitle,
			      100, -5., 5., 80, 0.9984, 1. );
      kHisBk++;
    }
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3054;
      hN3054 << kOffHi + 3054;
      hTitle = "y-part vs x-part on mirror up, reconst.d";
      hRC3054 = new CsHist2D( hN3054.str(), hTitle,
			      100, -2000., 2000., 100, -2000., 2000. );
      kHisBk++;
      stringstream hN3055;
      hN3055 << kOffHi + 3055;
      hTitle = "y-part vs x-part on mirror down, reconst.d";
      hRC3055 = new CsHist2D( hN3055.str(), hTitle,
			      100, -2000., 2000., 100, -2000., 2000. );
      kHisBk++;
      stringstream hN3056;
      hN3056 << kOffHi + 3056;
      hTitle = "y-part vs x-part on detectors, reconst.d";
      hRC3056 = new CsHist2D( hN3056.str(), hTitle,
			      100, -1600., 1600., 100, -1600., 1600. );
      kHisBk++;
    }
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3057;
      hN3057 << kOffHi + 3057;
      hTitle = "thphoton, beta = 1";
      hRC3057 = new CsHist1D( hN3057.str(), hTitle, 140, 0., 70. );
      kHisBk++;
      stringstream hN3058;
      hN3058 << kOffHi + 3058;
      hTitle = "thphoton-thphrec, beta = 1";
      hRC3058 = new CsHist1D( hN3058.str(), hTitle, 100, -10., 10. );
      kHisBk++;
    }
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3059;
      hN3059 << kOffHi + 3059;
      hTitle = "n-1";
      //hRC3059 = new CsHist1D( hN3059.str(), hTitle, 100, 0.0002, 0.00245 );
      //hRC3059 = new CsHist1D( hN3059.str(), hTitle, 200, 0.00050, 0.00250 );
      hRC3059 = new CsHist1D( hN3059.str(), hTitle, 400, 0.0002, 0.0022 );
      kHisBk++;
      stringstream hN3060;
      hN3060 << kOffHi + 3060;
      hTitle = "n-1 PMT only";
      hRC3060 = new CsHist1D( hN3060.str(), hTitle, 400, 0.0002, 0.0022 );
      kHisBk++;
    }
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3061;
      hN3061 << kOffHi + 3061;
      hTitle = "photons/ring, det up";
      hRC3061 = new CsHist1D( hN3061.str(), hTitle, 30, 0., 60. );
      kHisBk++;

      stringstream hN3062;
      hN3062 << kOffHi + 3062;
      hTitle = "photons/ring, split";
      hRC3062 = new CsHist1D( hN3062.str(), hTitle, 30, 0., 60. );
      kHisBk++;

      stringstream hN3063;
      hN3063 << kOffHi + 3063;
      hTitle = "photons/ring, det down";
      hRC3063 = new CsHist1D( hN3063.str(), hTitle, 30, 0., 60. );
      kHisBk++;

      stringstream hN3064;
      hN3064 << kOffHi + 3064;
      hTitle = "photons/ring, det down vs det up";
      hRC3064 = new CsHist2D( hN3064.str(), hTitle, 60, 0., 60.,
			      60, 0., 60. );
      kHisBk++;

      stringstream hN3065;
      hN3065 << kOffHi + 3065;
      hTitle = "photons/ring vs theta rec";
      hRC3065 = new CsHist1D( hN3065.str(), hTitle, 120, 0., 60. );
      kHisBk++;
    }
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3067;
      hN3067 << kOffHi + 3067;
      hTitle = "theta rec";
      hRC3067 = new CsHist1D( hN3067.str(), hTitle, 120, 0., 60. );
      kHisBk++;
    }

//- CsRCEventAnalysis::partIdent() (cont.):
    hRC3072 = NULL;
    hRC3073 = NULL;
    hRC3074 = NULL;
    hRC3075 = NULL;
    hRC3076 = NULL;
    hRC3077 = NULL;
    hRC3078 = NULL;
    hRC3079 = NULL;
    level = 1;
    bk = bookHis && MCEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3072;
      hN3072 << kOffHi + 3072;
      hTitle = "pions accepted vs mom.";
      hRC3072 = new CsHist1D( hN3072.str(), hTitle, mBi, momMin, momMax );
      kHisBk++;
      stringstream hN3073;
      hN3073 << kOffHi + 3073;
      hTitle = "kaons accepted vs mom.";
      hRC3073 = new CsHist1D( hN3073.str(), hTitle, mBi, momMin, momMax );
      kHisBk++;
      stringstream hN3074;
      hN3074 << kOffHi + 3074;
      hTitle = "protons accepted vs mom.";
      hRC3074 = new CsHist1D( hN3074.str(), hTitle, mBi, momMin, momMax );
      kHisBk++;
      stringstream hN3075;
      hN3075 << kOffHi + 3075;
      hTitle = "part. accepted vs n photons";
      hRC3075 = new CsHist1D( hN3075.str(), hTitle, 30, 0., 60. );
      kHisBk++;
      stringstream hN3076;
      hN3076 << kOffHi + 3076;
      hTitle = "part. identif.d vs mom.";
      hRC3076 = new CsHist1D( hN3076.str(), hTitle, mBi, momMin, momMax );
      kHisBk++;
      stringstream hN3077;
      hN3077 << kOffHi + 3077;
      hTitle = "part. identif.d vs mom.";
      hRC3077 = new CsHist1D( hN3077.str(), hTitle, mBi, momMin, momMax );
      kHisBk++;
      stringstream hN3078;
      hN3078 << kOffHi + 3078;
      hTitle = "part. identif.d vs mom.";
      hRC3078 = new CsHist1D( hN3078.str(), hTitle, mBi, momMin, momMax );
      kHisBk++;
      stringstream hN3079;
      hN3079 << kOffHi + 3079;
      hTitle = "part. identif.d vs n photons";
      hRC3079 = new CsHist1D( hN3079.str(), hTitle, 30, 0., 60. );
      kHisBk++;
    }

//- booking in CsRCRing::checkLikelihood()
//      CsHist1D* vRC3070;
//      vector<CsHist1D*> vRC3080;

//- CsRCEventAnalysis::MCMonitor() (cont.):
    hRC3101 = NULL;
    hRC3111 = NULL;
    level = 1;
    bk = bookHis && MCEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3101;
      hN3101 << kOffHi + 3101;
      hTitle = "hits MC vs Hits Ring - Accepted";
      hRC3101 = new CsHist2D( hN3101.str(), hTitle,
			    60, 0., 60., 60, 0., 60. );
      kHisBk++;
      stringstream hN3111;
      hN3111 << kOffHi + 3111;
      hTitle = "hits MC vs Hits Ring - Reconstructed";
      hRC3111 = new CsHist2D( hN3111.str(), hTitle,
			    60, 0., 60., 60, 0., 60. );
      kHisBk++;
    }

//- booking in CsRCRing::checkChiSquare()
//      vector<CsHist1D*> vRC3120;
//- booking in CsRCRing::checkLikeVsIndex()
//      vector<CsHist1D*> vRC3160;

//- CsRCEventAnalysis::cluStructure():
    hRC3250 = NULL;
    hRC3252 = NULL;
    hRC3255 = NULL;
    hRC3256 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3250;
      hN3250 << kOffHi + 3250;
      hTitle = "reconst.d cluster multipl.";
      hRC3250 = new CsHist1D( hN3250.str(), hTitle, 20, 0., 20. );
      kHisBk++;

      stringstream hN3252;
      hN3252 << kOffHi + 3252;
      hTitle = "ring clu PH vs phi";
      hRC3252 = new CsHist1D( hN3252.str(), hTitle, 90, 0., 360. );
      kHisBk++;

      stringstream hN3255;
      hN3255 << kOffHi + 3255;
      hTitle = "cluster shape around pad of max PH, accpt rings";
      hRC3255 = new CsHist2D( hN3255.str(), hTitle,
			      9, -4.5, 4.5, 9, -4.5, 4.5 );
      kHisBk++;
      stringstream hN3256;
      hN3256 << kOffHi + 3256;
      hTitle = "cluster PH shape around pad of max PH, accpt rings";
      hRC3256 = new CsHist2D( hN3256.str(), hTitle,
			      9, -4.5, 4.5, 9, -4.5, 4.5 );
      kHisBk++;
    }

//- booking in CsRCPartPhotons::checkLikeReso() {
//      std::vector<CsHist2D*> vRC3270(0-6);

//- booking in CsRCRing::checkChiReso() {
//      std::vector<CsHist2D*> vRC3280(0-6);

//- booking in CsRCPartPhotons::getLikeProb() {
//      std::vector<CsHist2D*> vRC3290(0-2);

//- booking in CsRCCluster::getBackWgt
//      vector<CsHist2D*> vRC3300[0-3]

//- booking in CsRCPartPhotons::GetPMTPhotonPosition
//      vector<CsHist2D*> vRC3310[0-3] (3311-3314)
//- booking in CsRCPartPhotons::doPMTOptCorr
//      vector<CsHist2D*> vRC3320[0-63] (3321-3384)
//- booking in CsRCEventPartPhotons::checkAMPSCorr()
//      vector<CsHist2D*> vRC3390[0-3] (3391-3394)
//- booking in CsRCEventRings::checkAMPSCorr()
//      vector<CsHist2D*> vRC3395[0-3] (3396-3399)

//- booking in CsRCEventRings::checkPhotonAngle()
//      vector<CsHist2D*> vRC3400[0-3], vRC3405[0-3], vRC3410[0-5],
//                        vRC3420[0-5], vRC3430[0-3]

//- booking in CsRCEventAnalysis::checkThetaLikeMax()
//      vector<CsHist1D*> vRC3440[0-4],

//- booking in CsRCEventAnalysis::checkLikeDistr()
//      vector<CsHist2D*> vRC3450[0-3], vRC3460[0-3], vRC3470[0-3],
//                        vRC3480[0-3], vRC3490[0-5]



//- CsRCPartPhotons::getPartPhotons():
    hRC3501 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3501;
      hN3501 << kOffHi + 3501;
      hTitle = "part. path length vs mom";
      hRC3501 = new CsHist2D( hN3501.str(), hTitle,
			      100, 2600., 3100., 60, 0., 60. );
      kHisBk++;
    }
//- CsRCPartPhotons::getPartPhotons():
    hRC3502 = NULL;
    hRC3503 = NULL;
    hRC3504 = NULL;
    level = 2;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3502;
      hN3502 << kOffHi + 3502;
      hTitle = "phiA vs phi";
      hRC3502 = new CsHist2D( hN3502.str(), hTitle,
			      90, 0., 360., 90, 0., 360. );
      kHisBk++;
      stringstream hN3503;
      hN3503 << kOffHi + 3503;
      hTitle = "det - mirror geometry";
      hRC3503 = new CsHist2D( hN3503.str(), hTitle, 
			      200, 5500., 9500., 150, 0., 3000. );
      kHisBk++;
      stringstream hN3504;
      hN3504 << kOffHi + 3504;
      hTitle = "theClu(Qz) - theClu";
      hRC3504 = new CsHist2D( hN3504.str(), hTitle, 
			      100, -2., 2., 100, 0., 360. );
      kHisBk++;
    }

//- CsRCMirrors::doSelMirrors():
    hRC3505 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3505;
      hN3505 << kOffHi + 3505;
      hTitle = "y- vs x-particle on mirrors";
      hRC3505 = new CsHist2D( hN3505.str(), hTitle,
			      200, -2000., 2000., 200, -2000., 2000. );
      kHisBk++;
    }

//- CsRCPartPhotons::getPartPhotons():
    hRC3506 = NULL;
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3506;
      hN3506 << kOffHi + 3506;
      hTitle = "y-clu vs x-clu (processed) on detectors";
      hRC3506 = new CsHist2D( hN3506.str(), hTitle,
			      200, -1600., 1600., 200, -1600., 1600. );
      kHisBk++;
    }

//- CsRCMirrors::doAliMirrors():
    hRC3507 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3507;
      hN3507 << kOffHi + 3507;
      hTitle = "y- vs x-particle (processed) on mirrors";
      hRC3507 = new CsHist2D( hN3507.str(), hTitle,
			      200, -2000., 2000., 200, -2000., 2000. );
      kHisBk++;
    }

//- CsRCPartPhotons::getSplitLimit():
    hRC3508 = NULL;
    hRC3509 = NULL;
    level = 2;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3508;
      hN3508 << kOffHi + 3508;
      hTitle = "split ring cut down vs up";
      hRC3508 = new CsHist2D( hN3508.str(), hTitle,
			      60, 200., 500., 60, -500., -200. );
      kHisBk++;
      stringstream hN3509;
      hN3509 << kOffHi + 3509;
      hTitle = "paDet down vs up";
      hRC3509 = new CsHist2D( hN3509.str(), hTitle,
			      100, 200., 1000., 100, -1000., -200. );
      kHisBk++;
    }

//- CsRCPartPhotons::getPartPhotons():
    hRC3510 = NULL;
    level = 2;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3510;
      hN3510 << kOffHi + 3510;
      hTitle = "theta corr mirr - theta Nocorr vs phi";
      hRC3510 = new CsHist2D( hN3510.str(), hTitle,
			      120, -3., 3., 90, 0., 360. );
      kHisBk++;
    }

//- CsRCEventAnalysis::dataMonitor() (cont.):
    hRC3519 = NULL;
    hRC3520 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3519;
      hN3519 << kOffHi + 3519;
      hTitle = "theta part-det vs theta part";
      hRC3519 = new CsHist2D( hN3519.str(), hTitle,
			      80, 0., 40., 90, 0., 30. );
      kHisBk++;
      stringstream hN3520;
      hN3520 << kOffHi + 3520;
      hTitle = "theta part-det vs phi part";
      hRC3520 = new CsHist2D( hN3520.str(), hTitle,
			      80, 0., 40., 90, 0., 360. );
      kHisBk++;
    }

//- CsRCEventRings::singlePhoton():
    for( int kh=0; kh<4; kh++ ) vRC3520.push_back( NULL );
    hRC3525 = NULL;
    hRC3526 = NULL;
    hRC3527 = NULL;
    hRC3528 = NULL;
    hRC3529 = NULL;
    hRC3530 = NULL;
    hRC3531 = NULL;
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      vRC3520.clear();
      stringstream hN3521;
      hN3521 << kOffHi + 3521;
      hTitle = "phi vs thphoton-thphrec, theta part-det = 7.5";
      vRC3520.push_back( new CsHist2D( hN3521.str(), hTitle,
				       100, -10., 10., 60, 0., 360. ) );
      kHisBk++;
      stringstream hN3522;
      hN3522 << kOffHi + 3522;
      hTitle = "phi vs thphoton-thphrec, theta part-det = 12.5";
      vRC3520.push_back( new CsHist2D( hN3522.str(), hTitle,
				       100, -10., 10., 60, 0., 360. ) );
      kHisBk++;
      stringstream hN3523;
      hN3523 << kOffHi + 3523;
      hTitle = "phi vs thphoton-thphrec, theta part-det = 17.5";
      vRC3520.push_back( new CsHist2D( hN3523.str(), hTitle,
				       100, -10., 10., 60, 0., 360. ) );
      kHisBk++;
      stringstream hN3524;
      hN3524 << kOffHi + 3524;
      hTitle = "phi vs thphoton-thphrec, theta part-det = 22.5";
      vRC3520.push_back( new CsHist2D( hN3524.str(), hTitle,
				       100, -10., 10., 60, 0., 360. ) );
      kHisBk++;
      stringstream hN3525;
      hN3525 << kOffHi + 3525;
      hTitle = "theta photon - theta ring, beta = 1";
      hRC3525 = new CsHist1D( hN3525.str(), hTitle, 100, -10., 10. );
      kHisBk++;
      stringstream hN3526;
      hN3526 << kOffHi + 3526;
      hTitle = "phi photon vs theta photon - theta ring, beta = 1 (Nzero)";
      hRC3526 = new CsHist2D( hN3526.str(), hTitle,
			      100, -10., 10., 60, 0., 360. );
      kHisBk++;
    }
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3527;
      hN3527 << kOffHi + 3527;
      hTitle = "theta photon vs theta photon - theta ring";
      hRC3527 = new CsHist2D( hN3527.str(), hTitle,
			      100, -10., 10., 60, 0., 60. );
      kHisBk++;
    }
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3528;
      hN3528 << kOffHi + 3528;
      hTitle = "theta photon vs theta photon - theta ring";
      hRC3528 = new CsHist2D( hN3528.str(), hTitle,
			      100, -10., 10., 72, 12., 16. );
      kHisBk++;
      stringstream hN3529;
      hN3529 << kOffHi + 3529;
      hTitle = "theta particle-detector";
      hRC3529 = new CsHist1D( hN3529.str(), hTitle, 80, 0., 40. );
      kHisBk++;
      stringstream hN3530;
      hN3530 << kOffHi + 3530;
      hTitle = "phi photon vs thetaBM";
      hRC3530 = new CsHist2D( hN3530.str(), hTitle,
			      100, -10., 10., 60, 0., 360. );
      kHisBk++;
      stringstream hN3531;
      hN3531 << kOffHi + 3531;
      hTitle = "phiA photon vs theta photon - theta ring, beta = 1";
      hRC3531 = new CsHist2D( hN3531.str(), hTitle,
			      100, -10., 10., 60, 0., 360. );
      kHisBk++;
    }

//- CsRCEventPartPhotons::checkSignal():
    hRC3533 = NULL;
    for( int kh=0; kh<4; kh++ ) vRC3535.push_back( NULL );
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3533;
      hN3533 << kOffHi + 3533;
      hTitle = "mom. vs thphoton - theReco";
      hRC3533 = new CsHist2D( hN3533.str(), hTitle,
			      100, -10., 10., 60, 0., 60. );
      kHisBk++;

      vRC3535.clear();
      stringstream hN3536;
      hN3536 << kOffHi + 3536;
      hTitle = "mom. vs thphoton - theIpo, pion hypo.";
      vRC3535.push_back( new CsHist2D( hN3536.str(), hTitle,
				       100, -10., 10., 60, 0., 60. ) );
      kHisBk++;
      stringstream hN3537;
      hN3537 << kOffHi + 3537;
      hTitle = "mom. vs thphoton - theIpo, kaon hypo.";
      vRC3535.push_back( new CsHist2D( hN3537.str(), hTitle,
				       100, -10., 10., 60, 0., 60. ) );
      kHisBk++;
      stringstream hN3538;
      hN3538 << kOffHi + 3538;
      hTitle = "mom. vs thphoton - theIpo, proton hypo.";
      vRC3535.push_back( new CsHist2D( hN3538.str(), hTitle,
				       100, -10., 10., 60, 0., 60. ) );
      kHisBk++;
      stringstream hN3539;
      hN3539 << kOffHi + 3539;
      hTitle = "mom. counts for thphoton - theIpo";
      vRC3535.push_back( new CsHist2D( hN3539.str(), hTitle,
				       3, 0., 3., 60, 0., 60. ) );
      kHisBk++;
    }

//- CsRCEventRings::checkSignal():
    for( int kh=0; kh<4; kh++ ) vRC3540.push_back( NULL );
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      vRC3540.clear();
      stringstream hN3541;
      hN3541 << kOffHi + 3541;
      hTitle = "mom. vs thphoton - theIpo, pion hypo.";
      vRC3540.push_back( new CsHist2D( hN3541.str(), hTitle,
				       100, -10., 10., 60, 0., 60. ) );
      kHisBk++;
      stringstream hN3542;
      hN3542 << kOffHi + 3542;
      hTitle = "mom. vs thphoton - theIpo, kaon hypo.";
      vRC3540.push_back( new CsHist2D( hN3542.str(), hTitle,
				       100, -10., 10., 60, 0., 60. ) );
      kHisBk++;
      stringstream hN3543;
      hN3543 << kOffHi + 3543;
      hTitle = "mom. vs thphoton - theIpo, proton hypo.";
      vRC3540.push_back( new CsHist2D( hN3543.str(), hTitle,
				       100, -10., 10., 60, 0., 60. ) );
      kHisBk++;
      stringstream hN3544;
      hN3544 << kOffHi + 3544;
      hTitle = "mom. counts for thphoton - theIpo";
      vRC3540.push_back( new CsHist2D( hN3544.str(), hTitle,
				       3, 0., 3., 60, 0., 60. ) );
      kHisBk++;
    }

//- CsRCEventPartPhotons::checkSignal():
    for( int kh=0; kh<4; kh++ ) vRC3545.push_back( NULL );
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      vRC3545.clear();
      stringstream hN3546;
      hN3546 << kOffHi + 3546;
      hTitle = "ddPaDet vs thphoton - theIpo, pion hypo.";
      vRC3545.push_back( new CsHist2D( hN3546.str(), hTitle,
				       100, -10., 10., 100, 0., 1000. ) );
      kHisBk++;
      stringstream hN3547;
      hN3547 << kOffHi + 3547;
      hTitle = "ddPaDet vs thphoton - theIpo, kaon hypo.";
      vRC3545.push_back( new CsHist2D( hN3547.str(), hTitle,
				       100, -10., 10., 100, 0., 1000. ) );
      kHisBk++;
      stringstream hN3548;
      hN3548 << kOffHi + 3548;
      hTitle = "ddPaDet vs thphoton - theIpo, proton hypo.";
      vRC3545.push_back( new CsHist2D( hN3548.str(), hTitle,
				       100, -10., 10., 100, 0., 1000. ) );
      kHisBk++;
      stringstream hN3549;
      hN3549 << kOffHi + 3549;
      hTitle = "ddPaDet counts for thphoton - theIpo";
      vRC3545.push_back( new CsHist2D( hN3549.str(), hTitle,
				       3, 0., 3., 100, 0., 1000. ) );
      kHisBk++;
    }

//- CsRCEventPartPhotons::bkgrPhotons():
    hRC3550 = NULL;
    hRC3551 = NULL;
    level = 1;
    bk = bookHis && MCEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3550;
      hN3550 << kOffHi + 3550;
      hTitle = "counts theta ring vs theta photon";
      hRC3550 = new CsHist1D( hN3550.str(), hTitle, 54, 4., 58. );
      kHisBk++;
      stringstream hN3551;
      hN3551 << kOffHi + 3551;
      hTitle = "theta ring vs theta photon";
      hRC3551 = new CsHist2D( hN3551.str(), hTitle,
      //		      216, 4., 58., 54, 4., 58. );
			      280, 0., 70., 54, 4., 58. );
      kHisBk++;
    }

//- CsRCEventAnalysis::ringSelection():
    hRC3553 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3553;
      hN3553 << kOffHi + 3553;
      hTitle = "x-p on dets theta ring";
      hRC3553 = new CsHist2D( hN3553.str(), hTitle, 
			      140, 0., 70., 100, -500., 500. );
      kHisBk++;
    }

//- CsRCEventPartPhotons::rawSignal():
    hRC3554 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3554;
      hN3554 << kOffHi + 3554;
      hTitle = "raw photon signal, theta";
      hRC3554 = new CsHist2D( hN3554.str(), hTitle, 100, 0., 100.,
			      2, 0., 2. );
      kHisBk++;
    }
    hRC3555 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3555;
      hN3555 << kOffHi + 3555;
      hTitle = "raw signal on detectors, R";
      hRC3555 = new CsHist1D( hN3555.str(), hTitle, 100, 0., 250. );
      kHisBk++;
    }
    hRC3556 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3556;
      hN3556 << kOffHi + 3556;
      hTitle = "raw signal on detectors, ring";
      hRC3556 = new CsHist2D( hN3556.str(), hTitle, 
			      100, -250., 250., 100, -250., 250. );
      kHisBk++;
    }

//- CsRCEventAnalysis::ringSelection():
    hRC3557 = NULL;
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3557;
      hN3557 << kOffHi + 3557;
      hTitle = "raw signal on detectors, R";
      hRC3557 = new CsHist1D( hN3557.str(), hTitle, 100, 0., 250. );
      kHisBk++;
    }
    hRC3558 = NULL;
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3558;
      hN3558 << kOffHi + 3558;
      hTitle = "raw signal on detectors, ring";
      hRC3558 = new CsHist2D( hN3558.str(), hTitle, 
			      100, -250., 250., 100, -250., 250. );
      kHisBk++;
    }

//- CsRCEventRings::singlePhoton():
    hRC3559 = NULL;
    level = -1;      //      book always
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3559;
      hN3559 << kOffHi + 3559;
      hTitle = "n-1";
      hRC3559 = new CsHist1D( hN3559.str(), hTitle, 200, 0.0002, 0.0022 );
      kHisBk++;
    }

//- CsRCEventRings::ringSignal():
    for( int kh=0; kh<4; kh++ ) vRC3570.push_back( NULL );
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      vRC3570.clear();
      stringstream hN3571;
      hN3571 << kOffHi + 3571;
      hTitle = "integrated ring signal";
      vRC3570.push_back( new CsHist2D( hN3571.str(), hTitle,
				       60, -60., 60., 60, -60., 60. ) );
      stringstream hN3572;
      hN3572 << kOffHi + 3572;
      hTitle = "integrated ring signal";
      vRC3570.push_back( new CsHist2D( hN3572.str(), hTitle,
				       60, -60., 60., 60, -60., 60. ) );
      stringstream hN3573;
      hN3573 << kOffHi + 3573;
      hTitle = "integrated ring signal";
      vRC3570.push_back( new CsHist2D( hN3573.str(), hTitle,
				       60, -60., 60., 60, -60., 60. ) );
      stringstream hN3574;
      hN3574 << kOffHi + 3574;
      hTitle = "integrated ring signal";
      vRC3570.push_back( new CsHist2D( hN3574.str(), hTitle,
				       60, -60., 60., 60, -60., 60. ) );
      kHisBk = kHisBk + 4;
    }

//- CsRCEventAnalysis::ringSelection():
    for( int kh=0; kh<4; kh++ ) vRC3580.push_back( NULL );
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      vRC3580.clear();
      stringstream hN3581;
      hN3581 << kOffHi + 3581;
      hTitle = "integrated ring signal";
      vRC3580.push_back( new CsHist2D( hN3581.str(), hTitle,
				       60, -60., 60., 60, -60., 60. ) );
      stringstream hN3582;
      hN3582 << kOffHi + 3582;
      hTitle = "integrated ring signal";
      vRC3580.push_back( new CsHist2D( hN3582.str(), hTitle,
				       60, -60., 60., 60, -60., 60. ) );
      stringstream hN3583;
      hN3583 << kOffHi + 3583;
      hTitle = "integrated ring signal";
      vRC3580.push_back( new CsHist2D( hN3583.str(), hTitle,
				       60, -60., 60., 60, -60., 60. ) );
      stringstream hN3584;
      hN3584 << kOffHi + 3584;
      hTitle = "integrated ring signal";
      vRC3580.push_back( new CsHist2D( hN3584.str(), hTitle,
				       60, -60., 60., 60, -60., 60. ) );
      kHisBk = kHisBk + 4;
    }

//- CsRCEventRings::singlePhoton():
    hRC3560 = NULL;
    hRC3561 = NULL;
    hRC3562 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3560;
      hN3560 << kOffHi + 3560;
      hTitle = "counts theta ring vs theta photon";
      hRC3560 = new CsHist1D( hN3560.str(), hTitle, 54, 4., 58. );
      kHisBk++;
      stringstream hN3561;
      hN3561 << kOffHi + 3561;
      hTitle = "theta ring vs theta photon";
      hRC3561 = new CsHist2D( hN3561.str(), hTitle,
      //		      216, 4., 58., 54, 4., 58. );
			      280, 0., 70., 54, 4., 58. );
      kHisBk++;
      stringstream hN3562;
      hN3562 << kOffHi + 3562;
      hTitle = "counts, beta = 1";
      hRC3562 = new CsHist1D( hN3562.str(), hTitle, 10, 0., 10. );
      kHisBk++;
    }

//- CsRCEventPartPhotons::checkSignal()
    hRC3563 = NULL;
    level = 2;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3563;
      hN3563 << kOffHi + 3563;
      hTitle = "APV background tests";
      hRC3563 = new CsHist2D( hN3563.str(), hTitle,
			      200, -100., 100., 40, 0., 40. );
      kHisBk++;
    }

//- CsRCEventRings::singlePhoton(): APV tests
    for( int kh=0; kh<4; kh++ ) vRC3565.push_back( NULL );
    hRC3569 = NULL;
    level = 2;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      vRC3565.clear();
      stringstream hN3565;
      hN3565 << kOffHi + 3565;
      hTitle = "theta ring vs theta photon";
      vRC3565.push_back( new CsHist2D( hN3565.str(), hTitle,
				       280, 0., 70., 54, 4., 58. ) );
      stringstream hN3566;
      hN3566 << kOffHi + 3566;
      hTitle = "theta ring vs theta photon";
      vRC3565.push_back( new CsHist2D( hN3566.str(), hTitle,
				       280, 0., 70., 54, 4., 58. ) );
      stringstream hN3567;
      hN3567 << kOffHi + 3567;
      hTitle = "theta ring vs theta photon";
      vRC3565.push_back( new CsHist2D( hN3567.str(), hTitle,
				       280, 0., 70., 54, 4., 58. ) );
      stringstream hN3568;
      hN3568 << kOffHi + 3568;
      hTitle = "theta ring vs theta photon";
      vRC3565.push_back( new CsHist2D( hN3568.str(), hTitle,
				       280, 0., 70., 54, 4., 58. ) );
      stringstream hN3569;
      hN3569 << kOffHi + 3569;
      hTitle = "counts theta ring vs theta photon";
      hRC3569 = new CsHist2D( hN3569.str(), hTitle,
			      54, 4., 58., 4, 0., 4. );
    }

//- CsRCEventRings::singlePhoton():
    for( int kh=0; kh<6; kh++ ) vRC3590.push_back( NULL );
    hRC3597 = NULL;
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      vRC3590.clear();
      stringstream hN3591;
      hN3591 << kOffHi + 3591;
      hTitle = "theta ring vs theta photon";
      vRC3590.push_back( new CsHist2D( hN3591.str(), hTitle,
				       280, 0., 70., 54, 4., 58. ) );
      stringstream hN3592;
      hN3592 << kOffHi + 3592;
      hTitle = "theta ring vs theta photon";
      vRC3590.push_back( new CsHist2D( hN3592.str(), hTitle,
				       280, 0., 70., 54, 4., 58. ) );
      stringstream hN3593;
      hN3593 << kOffHi + 3593;
      hTitle = "theta ring vs theta photon";
      vRC3590.push_back( new CsHist2D( hN3593.str(), hTitle,
				       280, 0., 70., 54, 4., 58. ) );
      stringstream hN3594;
      hN3594 << kOffHi + 3594;
      hTitle = "theta ring vs theta photon";
      vRC3590.push_back( new CsHist2D( hN3594.str(), hTitle,
				       280, 0., 70., 54, 4., 58. ) );
      stringstream hN3595;
      hN3595 << kOffHi + 3595;
      hTitle = "theta ring vs theta photon";
      vRC3590.push_back( new CsHist2D( hN3595.str(), hTitle,
				       280, 0., 70., 54, 4., 58. ) );
      stringstream hN3596;
      hN3596 << kOffHi + 3596;
      hTitle = "theta ring vs theta photon";
      vRC3590.push_back( new CsHist2D( hN3596.str(), hTitle,
				       280, 0., 70., 54, 4., 58. ) );
      kHisBk = kHisBk + 6;
      stringstream hN3597;
      hN3597 << kOffHi + 3597;
      hTitle = "counts theta ring vs theta photon";
      hRC3597 = new CsHist2D( hN3597.str(), hTitle, 
			      54, 4., 58., 6, 0., 600. );
      kHisBk++;
    }

//- CsRCEventAnalysis::dataMonitor() (cont.):
    hRC3601 = NULL;
    hRC3602 = NULL;
    hRC3603 = NULL;
    hRC3604 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3601;
      hN3601 << kOffHi + 3601;
      hTitle = "y-part vs x-part, entr. wind.";
      hRC3601 = new CsHist2D( hN3601.str(), hTitle,
			      100, -1000., 1000., 100, -1000., 1000. );
      kHisBk++;
      stringstream hN3602;
      hN3602 << kOffHi + 3602;
      hTitle = "m-part vs l-part, entr. wind.";
      hRC3602 = new CsHist2D( hN3602.str(), hTitle,
			      100, -0.1, 0.1, 100, -0.1, 0.1 );
      kHisBk++;
      stringstream hN3603;
      hN3603 << kOffHi + 3603;
      hTitle = "l-part vs x-part, entr. wind.";
      hRC3603 = new CsHist2D( hN3603.str(), hTitle,
			      100, -1000., 1000., 100, -0.1, 0.1 );
      kHisBk++;
      stringstream hN3604;
      hN3604 << kOffHi + 3604;
      hTitle = "m-part vs y-part, entr. wind.";
      hRC3604 = new CsHist2D( hN3604.str(), hTitle,
			      100, -1000., 1000., 100, -0.1, 0.1 );
      kHisBk++;
    }

//- CsRCEventParticles::exitWindow():
    hRC3605 = NULL;
    hRC3606 = NULL;
    hRC3607 = NULL;
    hRC3608 = NULL;
    hRC3609 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3605;
      hN3605 << kOffHi + 3605;
      hTitle = "dRR-part vs mom., exit wind.";
      hRC3605 = new CsHist2D( hN3605.str(), hTitle,
			      mBi*2, momMin, momMax, 200, 0., 100. );
      kHisBk++;
      stringstream hN3606;
      hN3606 << kOffHi + 3606;
      hTitle = "dx-part vs dy-part, exit wind.";
      hRC3606 = new CsHist2D( hN3606.str(), hTitle,
			      100, -10., 10., 100, -10., 10. );
      kHisBk++;
      stringstream hN3607;
      hN3607 << kOffHi + 3607;
      hTitle = "dRR-normalized vs mom., exit wind.";
      hRC3607 = new CsHist2D( hN3607.str(), hTitle,
			      mBi*2, momMin, momMax, 100, 0., 10. );
      kHisBk++;
      stringstream hN3608;
      hN3608 << kOffHi + 3608;
      hTitle = "dx-part vs mom., exit wind.";
      hRC3608 = new CsHist2D( hN3608.str(), hTitle,
			      mBi*2, momMin, momMax, 200, -10., 10. );
      kHisBk++;
      stringstream hN3609;
      hN3609 << kOffHi + 3609;
      hTitle = "dy-part vs mom., exit wind.";
      hRC3609 = new CsHist2D( hN3609.str(), hTitle,
			      mBi*2, momMin, momMax, 200, -10., 10. );
      kHisBk++;
    }

//- CsRCEventAnalysis::dataMonitor() (cont.):
    hRC3610 = NULL;
    hRC3611 = NULL;
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3610;
      hN3610 << kOffHi + 3610;
      hTitle = "particle momentum, accepted";
      hRC3610 = new CsHist1D( hN3610.str(), hTitle, mBA, 0., momAcc );
      kHisBk++;
      stringstream hN3611;
      hN3611 << kOffHi + 3611;
      hTitle = "particle momentum, all";
      hRC3611 = new CsHist1D( hN3611.str(), hTitle, mBA, 0., momAcc );
      kHisBk++;
    }

//- CsRCPartPhotons::mcsPartCorr():
    hRC3613 = NULL;
    hRC3614 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3613;
      hN3613 << kOffHi + 3613;
      hTitle = "photon emiss pos corr Y vs X";
      hRC3613 = new CsHist2D( hN3613.str(), hTitle,
			      100, -1., 1., 100, -1., 1. );
      kHisBk++;
      stringstream hN3614;
      hN3614 << kOffHi + 3614;
      hTitle = "photon emiss direc corr Y vs X";
      hRC3614 = new CsHist2D( hN3614.str(), hTitle,
			      100, -1., 1., 100, -1., 1. );
      kHisBk++;
    }

//- CsRCEventAnalysis::MCMonitor() (cont.):
    hRC3615 = NULL;
    level = 1;
    bk = bookHis && MCEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3615;
      hN3615 << kOffHi + 3615;
      hTitle = "particle type, all";
      hRC3615 = new CsHist2D( hN3615.str(), hTitle, 30, 0., 30.,
			      10, 0., 10. );
      kHisBk++;
    }

//- CsRCEventParticles::partAnalysis():
    hRC3617 = NULL;
    level = 1;
    bk = bookHis && MCEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3617;
      hN3617 << kOffHi + 3617;
      hTitle = "particle type, rejected";
      hRC3617 = new CsHist1D( hN3617.str(), hTitle, 30, 0., 30. );
      kHisBk++;
    }
    hRC3618 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3618;
      hN3618 << kOffHi + 3618;
      hTitle = "particle ring vs dist.";
      hRC3618 = new CsHist2D( hN3618.str(), hTitle, 60, 0., 60.,
			      100, 0., 2000. );
      kHisBk++;
    }

//- CsRCEventAnalysis::MCMonitor() (cont.):
    for( int j=0; j<5; j++) vRC3620.push_back( NULL );
    hRC3631 = NULL;
    hRC3632 = NULL;
    level = 1;
    bk = bookHis && MCEvent && level<=levelBk;
    if( bk ) {
      vRC3620.clear();
      stringstream hN3623;
      hN3623 << kOffHi + 3623;
      hTitle = "e+- momentum, accepted";
      vRC3620.push_back( new CsHist1D( hN3623.str(), 
				       hTitle, mBA, 0., momAcc ) );
      kHisBk++;
      stringstream hN3624;
      hN3624 << kOffHi + 3624;
      hTitle = "mu+- momentum, accepted";
      vRC3620.push_back( new CsHist1D( hN3624.str(), 
				       hTitle, mBA, 0., momAcc ) );
      kHisBk++;
      stringstream hN3625;
      hN3625 << kOffHi + 3625;
      hTitle = "pi+- momentum, accepted";
      vRC3620.push_back( new CsHist1D( hN3625.str(), 
				       hTitle, mBA, 0., momAcc ) );
      kHisBk++;
      stringstream hN3626;
      hN3626 << kOffHi + 3626;
      hTitle = "K+- momentum, accepted";
      vRC3620.push_back( new CsHist1D( hN3626.str(), 
				       hTitle, mBA, 0., momAcc ) );
      kHisBk++;
      stringstream hN3627;
      hN3627 << kOffHi + 3627;
      hTitle = "ppbar+- momentum, accepted";
      vRC3620.push_back( new CsHist1D( hN3627.str(), 
				       hTitle, mBA, 0., momAcc ) );
      kHisBk++;

      stringstream hN3631;
      hN3631 << kOffHi + 3631;
      hTitle = "theta MC - theta Wgav";
      hRC3631 = new CsHist1D( hN3631.str(), hTitle, 200, -5., 5. );
      kHisBk++;
      stringstream hN3632;
      hN3632 << kOffHi + 3632;
      hTitle = "theta MC - theta Fit";
      hRC3632 = new CsHist2D( hN3632.str(), hTitle, 200, -5., 5.,
			      60, 0., 60. );
      kHisBk++;
    }

//- CsRCEventParticles::exitWindow():
    hRC3635 = NULL;
    hRC3636 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3635;
      hN3635 << kOffHi + 3635;
      hTitle = "Y- vs X-part at exit window";
      hRC3635 = new CsHist2D( hN3635.str(), hTitle,
			      200, -2000., 2000., 200, -2000., 2000.);
      kHisBk++;
      stringstream hN3636;
      hN3636 << kOffHi + 3636;
      hTitle = "Y- vs X-part at exit window - cut";
      hRC3636 = new CsHist2D( hN3636.str(), hTitle,
			      200, -2000., 2000., 200, -2000., 2000.);
      kHisBk++;
    }

//- CsRCEventParticles::partCorr():
    hRC3637 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3637;
      hN3637 << kOffHi + 3637;
      hTitle = "dtana-part vs dtanb-part";
      hRC3637 = new CsHist2D( hN3637.str(), hTitle,
			      100, -0.025, 0.025, 100, -0.025, 0.025 );
      kHisBk++;
    }

//- CsRCEventRings::singlePhoton():
    hRC3640 = NULL;
    for( int j=0; j<5; j++) vRC3640.push_back( NULL );
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3640;
      hN3640 << kOffHi + 3640;
      hTitle = "theta part-lens";
      hRC3640 = new CsHist1D( hN3640.str(), hTitle, 80, 0., 40.);
      kHisBk++;
      vRC3640.clear();
      stringstream hN3641;
      hN3641 << kOffHi + 3641;
      hTitle = "phi vs thphoton-thphrec, theta part-lens";
      vRC3640.push_back( new CsHist2D( hN3641.str(), hTitle,
				       100, -10., 10., 60, 0., 360. ) );
      kHisBk++;
      stringstream hN3642;
      hN3642 << kOffHi + 3642;
      hTitle = "phi vs thphoton-thphrec, theta part-lens";
      vRC3640.push_back( new CsHist2D( hN3642.str(), hTitle,
				       100, -10., 10., 60, 0., 360. ) );
      kHisBk++;
      stringstream hN3643;
      hN3643 << kOffHi + 3643;
      hTitle = "phi vs thphoton-thphrec, theta part-lens";
      vRC3640.push_back( new CsHist2D( hN3643.str(), hTitle,
				       100, -10., 10., 60, 0., 360. ) );
      kHisBk++;
      stringstream hN3644;
      hN3644 << kOffHi + 3644;
      hTitle = "phi vs thphoton-thphrec, theta part-lens";
      vRC3640.push_back( new CsHist2D( hN3644.str(), hTitle,
				       100, -10., 10., 60, 0., 360. ) );
      kHisBk++;
    }

//- CsRCEventRings::singlePhoton():
    for( int kh=0; kh<16; kh++ ) vRC3660.push_back( NULL );
    hRC3677 = NULL;
    hRC3678 = NULL;
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      vRC3660.clear();
      for( int kh=0; kh<16; kh++ ) {
	stringstream hN3660;
	hN3660 << kOffHi + 3661 + kh;
	hTitle = "theta ring vs theta photon";
	vRC3660.push_back( new CsHist2D( hN3660.str(), hTitle,
					 280, 0., 70., 54, 4., 58. ) );
	kHisBk++;
      }
      stringstream hN3677;
      hN3677 << kOffHi + 3677;
      hTitle = "counts theta ring vs theta photon";
      hRC3677 = new CsHist2D( hN3677.str(), hTitle, 
			      54, 4., 58., 16, 0., 16. );
      kHisBk++;
      stringstream hN3678;
      hN3678 << kOffHi + 3678;
      hTitle = "counts theta ring vs theta photon";
      hRC3678 = new CsHist2D( hN3678.str(), hTitle, 
			      54, 4., 58., 16, 0., 16. );
      kHisBk++;
    }

//- CsRCEventRings::singlePhoton():
    for( int kh=0; kh<6; kh++ ) vRC3680.push_back( NULL );
    hRC3687 = NULL;
    hRC3688 = NULL;
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      vRC3680.clear();
      stringstream hN3681;
      hN3681 << kOffHi + 3681;
      hTitle = "theta ring vs theta photon";
      vRC3680.push_back( new CsHist2D( hN3681.str(), hTitle,
				       280, 0., 70., 54, 4., 58. ) );
      stringstream hN3682;
      hN3682 << kOffHi + 3682;
      hTitle = "theta ring vs theta photon";
      vRC3680.push_back( new CsHist2D( hN3682.str(), hTitle,
				       280, 0., 70., 54, 4., 58. ) );
      stringstream hN3683;
      hN3683 << kOffHi + 3683;
      hTitle = "theta ring vs theta photon";
      vRC3680.push_back( new CsHist2D( hN3683.str(), hTitle,
				       280, 0., 70., 54, 4., 58. ) );
      stringstream hN3684;
      hN3684 << kOffHi + 3684;
      hTitle = "theta ring vs theta photon";
      vRC3680.push_back( new CsHist2D( hN3684.str(), hTitle,
				       280, 0., 70., 54, 4., 58. ) );
      stringstream hN3685;
      hN3685 << kOffHi + 3685;
      hTitle = "theta ring vs theta photon";
      vRC3680.push_back( new CsHist2D( hN3685.str(), hTitle,
				       280, 0., 70., 54, 4., 58. ) );
      stringstream hN3686;
      hN3686 << kOffHi + 3686;
      hTitle = "theta ring vs theta photon";
      vRC3680.push_back( new CsHist2D( hN3686.str(), hTitle,
				       280, 0., 70., 54, 4., 58. ) );
      kHisBk = kHisBk + 6;
      stringstream hN3687;
      hN3687 << kOffHi + 3687;
      hTitle = "counts theta ring vs theta photon";
      hRC3687 = new CsHist2D( hN3687.str(), hTitle, 
			      54, 4., 58., 6, 0., 600. );
      kHisBk++;
      stringstream hN3688;
      hN3688 << kOffHi + 3688;
      hTitle = "counts theta ring vs theta photon";
      hRC3688 = new CsHist2D( hN3688.str(), hTitle, 
			      100, -1500., 1500., 100, -1500., 1500. );
      kHisBk++;
    }
//- CsRCEventRings::singlePhoton():
    for( int kh=0; kh<6; kh++ ) vRC3690.push_back( NULL );
    hRC3697 = NULL;
    hRC3698 = NULL;
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      vRC3690.clear();
      stringstream hN3691;
      hN3691 << kOffHi + 3691;
      hTitle = "theta ring vs theta photon";
      vRC3690.push_back( new CsHist2D( hN3691.str(), hTitle,
				       280, 0., 70., 54, 4., 58. ) );
      stringstream hN3692;
      hN3692 << kOffHi + 3692;
      hTitle = "theta ring vs theta photon";
      vRC3690.push_back( new CsHist2D( hN3692.str(), hTitle,
				       280, 0., 70., 54, 4., 58. ) );
      stringstream hN3693;
      hN3693 << kOffHi + 3693;
      hTitle = "theta ring vs theta photon";
      vRC3690.push_back( new CsHist2D( hN3693.str(), hTitle,
				       280, 0., 70., 54, 4., 58. ) );
      stringstream hN3694;
      hN3694 << kOffHi + 3694;
      hTitle = "theta ring vs theta photon";
      vRC3690.push_back( new CsHist2D( hN3694.str(), hTitle,
				       280, 0., 70., 54, 4., 58. ) );
      stringstream hN3695;
      hN3695 << kOffHi + 3695;
      hTitle = "theta ring vs theta photon";
      vRC3690.push_back( new CsHist2D( hN3695.str(), hTitle,
				       280, 0., 70., 54, 4., 58. ) );
      stringstream hN3696;
      hN3696 << kOffHi + 3696;
      hTitle = "theta ring vs theta photon";
      vRC3690.push_back( new CsHist2D( hN3696.str(), hTitle,
				       280, 0., 70., 54, 4., 58. ) );
      kHisBk = kHisBk + 6;
      stringstream hN3697;
      hN3697 << kOffHi + 3697;
      hTitle = "counts theta ring vs theta photon";
      hRC3697 = new CsHist2D( hN3697.str(), hTitle, 
			      54, 4., 58., 6, 0., 600. );
      kHisBk++;
      stringstream hN3698;
      hN3698 << kOffHi + 3698;
      hTitle = "counts theta ring vs theta photon";
      hRC3698 = new CsHist2D( hN3698.str(), hTitle, 
			      100, -1500., 1500., 100, -1500., 1500. );
      kHisBk++;
    }

//- CsRCCircleFit::doHist():
    for( int j=0; j<5; j++) vRC3700.push_back( NULL );
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      vRC3700.clear();
      int npo = 100;
      stringstream hN3701;
      hN3701 << kOffHi + 3701;
      hTitle = "c0in - c0fit";
      vRC3700.push_back( new CsHist2D( hN3701.str(), hTitle, 100, -2.5, 2.5,
			               100, -2.5, 2.5 ) );
      kHisBk++;
      stringstream hN3702; 
      hN3702 << kOffHi + 3702;
      hTitle = "thphrec - thphfit";
      vRC3700.push_back( new CsHist2D( hN3702.str(), hTitle, 100, -1.0, 1.0,
			               npo, 0., float( npo ) ) );
      kHisBk++;
      stringstream hN3703;
      hN3703 << kOffHi + 3703;
      hTitle = "chisq / nu";
      vRC3700.push_back( new CsHist2D( hN3703.str(), hTitle, 100, 0., 10.,
			               npo, 0., float( npo ) ) );
      kHisBk++;
      stringstream hN3704;
      hN3704 << kOffHi + 3704;
      hTitle = "n iter";
      vRC3700.push_back( new CsHist2D( hN3704.str(), hTitle, 20, 0., 20.,
                                       npo, 0., float( npo ) ) );
      kHisBk++;
      stringstream hN3705;
      hN3705 << kOffHi + 3705;
      hTitle = "pulls";
      vRC3700.push_back( new CsHist2D( hN3705.str(), hTitle, 100, -5., 5.,
			               npo, 0., float( npo ) ) );
      kHisBk++;
    }

    for( int j=0; j<7; j++) vRC3710.push_back( NULL );
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      vRC3710.clear();
      int npo = 100;
      stringstream hN3711;
      hN3711 << kOffHi + 3711;
      hTitle = "c0in - c0fit";
      vRC3710.push_back( new CsHist2D( hN3711.str(), hTitle, 100, -10., 10 ,
			               100, -10., 10. ) );
      kHisBk++;
      stringstream hN3712; 
      hN3712 << kOffHi + 3712;
      hTitle = "radrec - radfit";
      vRC3710.push_back( new CsHist2D( hN3712.str(), hTitle, 100, -10., 10.,
			               npo, 0., float( npo ) ) );
      kHisBk++;
      stringstream hN3713;
      hN3713 << kOffHi + 3713;
      hTitle = "chisq / nu";
      vRC3710.push_back( new CsHist2D( hN3713.str(), hTitle, 100, 0., 10.,
			               npo, 0., float( npo ) ) );
      kHisBk++;
      stringstream hN3714;
      hN3714 << kOffHi + 3714;
      hTitle = "n iter";
      vRC3710.push_back( new CsHist2D( hN3714.str(), hTitle, 20, 0., 20.,
                                       npo, 0., float( npo ) ) );
      kHisBk++;
      stringstream hN3715;
      hN3715 << kOffHi + 3715;
      hTitle = "pulls";
      vRC3710.push_back( new CsHist2D( hN3715.str(), hTitle, 100, -5., 5.,
			               npo, 0., float( npo ) ) );
      kHisBk++;
      stringstream hN3716; 
      hN3716 << kOffHi + 3716;
      hTitle = "c0in - c0fit --- det UP";
      vRC3710.push_back( new CsHist2D( hN3716.str(), hTitle, 100, -10., 10.,
			               100, -10., 10. ) );
      kHisBk++;
      stringstream hN3717; 
      hN3717 << kOffHi + 3717;
      hTitle = "c0in - c0fit --- det DOWN";
      vRC3710.push_back( new CsHist2D( hN3717.str(), hTitle, 100, -10., 10.,
			               100, -10., 10. ) );
      kHisBk++;
    }

    for( int j=0; j<8; j++ ) vRC3720.push_back( NULL );
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      vRC3720.clear();
      for( int kh=0; kh<8; kh++ ) {
	stringstream hN3720;
	hN3720 << kOffHi + 3720 + kh;
	hTitle = "c0in - c0fit";
	vRC3720.push_back( new CsHist2D( hN3720.str(), hTitle, 100, -10., 10 ,
					 100, -10., 10. ) );
	kHisBk++;
      }
    }
    for( int j=0; j<8; j++ ) vRC3730.push_back( NULL );
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      vRC3730.clear();
      for( int kh=0; kh<8; kh++ ) {
	stringstream hN3730;
	hN3730 << kOffHi + 3730 + kh;
	hTitle = "c0in - c0fit";
	vRC3730.push_back( new CsHist2D( hN3730.str(), hTitle, 100, -10., 10 ,
					 100, -10., 10. ) );
	kHisBk++;
      }
    }

//- booking in CsRCRing::getDetEFit():
//        static vector<CsHist2D*> RC3740

//- CsRCMirrors::doAliMirrors():
    hRC3750 = NULL;
    hRC3751 = NULL;
    hRC3760 = NULL;
    hRC3765 = NULL;
    level = 1;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      stringstream hN3750;
      hN3750 << kOffHi + 3750;
      hTitle = "alignment data";
      hRC3750 = new CsHist1D( hN3750.str(), hTitle, 1000, 0., 1000. );
      kHisBk++;
      stringstream hN3751;
      hN3751 << kOffHi + 3751;
      hTitle = "yc-mirr.elem. vs xc-mirr.elem";
      hRC3751 = new CsHist2D( hN3751.str(), hTitle,
			      300, -3000., 3000., 300, -3000., 3000. );
      kHisBk++;
      stringstream hN3760;
      hN3760 << kOffHi + 3760;
      hTitle = "alignment fit data";
      hRC3760 = new CsHist1D( hN3760.str(), hTitle, 600, 0., 600. );
      kHisBk++;
      stringstream hN3765;
      hN3765 << kOffHi + 3765;
      hTitle = "particles on mirror element";
      hRC3765 = new CsHist2D( hN3765.str(), hTitle,
			      60, -3000., 3000., 60, -3000., 3000. );
      kHisBk++;
    }

//- booking in CsRCEventClusters::killHaloClus():
//        static vector<CsHist2D*> RC3770

//- booking moved to CsRCMirrors::doAliMirrors(): (20/9/00)
//        vRC3800.push_back( ... ->116
//        hRC3980 ...

//- booking in CsRCEventDisplay::doEveDisplay() 
//      vector<CsHist2D*> vRC4000;
//      vector<CsHist2D*> vRC4300;
//      vector<CsHist1D*> vRC4600;
//      vector<CsHist1D*> vRC4900;

//- booking in CsRCEventAnalysis::hitDisplay()
//      static vector<CsHist2D*> vRC5000; (0->19)
//      static vector<CsHist2D*> vRC5100; (0->19)
//      static vector<CsHist2D*> vRC5200; (0->19)
//- booking in CsRCEventAnalysis::BMDisplay()
//      static vector<CsHist2D*> vRC5300; (0->19)
//- obsolete and conflicting with the following...

//- booking in CsRCEventDisplay::sideDisplay() 
//      vector<CsHist1D*> vRC5300;
//- booking in CsRCEventDisplay::sideDisplay() 
//      vector<CsHist1D*> vRC5900;

//- booking in CsRCRing::peakScan() 
//      vector<CsHist1D*> vRC6000;
//      vector<CsHist1D*> vRC6100;
//- booking in CsRCRing::peakChiSearch() 
//      vector<CsHist1D*> vRC6200;
//- booking in CsRCRing::peakMassSearch() 
//      vector<CsHist1D*> vRC6200;
//- booking in CsRCRing::histPeakSearch() 
//      vector<CsHist1D*> vRC6300;

//- booking in CsRCEventPads::setThreshold()
//      vector<CsHist2D*> vRC6400;
//      vector<CsHist2D*> vRC6420;
//      vector<CsHist2D*> vRC6440;

//- booking in CsRCEventClusters::checkClusters()
//      vector<CsHist2D*> vRC6500; (0->19)

//- CsRCEventRings::singlePhoton():
    for( int kh=0; kh<15; kh++ ) vRC6550.push_back( NULL );
    level = 0;
    bk = bookHis && DataEvent && level<=levelBk;
    if( bk ) {
      vRC6550.clear();
      for( int kh=0; kh<15; kh++ ) {
	stringstream hN6550;
	hN6550 << kOffHi + 6550 + kh + 1;
	hTitle = "phi vs thphoton-thphorec";
        vRC6550.push_back( new CsHist2D( hN6550.str(), hTitle,
					 100, -10., 10., 60, 0., 360. ) );
	kHisBk++;
      }
    }

//- booking in CsRCEventPartPhotons::checkAMPSCorr()
//      vector<CsHist2D*> vRC6630[0-4] (6631-6634)
//      vector<CsHist2D*> vRC6635[0-4] (6636-6639)
//      vector<CsHist2D*> vRC6640[0-4] (6641-6644)

//- booking in CsRCEventAnalysis::moniCFRefInf()
//      vector<CsHist2D*> vRC6660[0-7] (6661-6668)
//      vector<CsHist2D*> vRC6670[0-3] (6671-6674)
//- booking in CsRCEventAnalysis::moniCFRefInfR()
//      vector<CsHist2D*> vRC6680[0-19] (6681-6700)

//- booking in CsRCEventRings::checkAMPSCorr()
//      vector<CsHist2D*> vRC6700[0-...] (6701-6705)
//      vector<CsHist2D*> vRC6710[0-...] (6711-6715)
//- booking in CsRCPartPhotons::checkLikeCorr()
//      vector<CsHist2D*> vRC6750[0-9] (6751-6760)
//      vector<CsHist2D*> vRC6760[0-8] (6761-6769)
//      vector<CsHist2D*> vRC6770[0-8] (6771-6779)
//      vector<CsHist2D*> vRC6780[0-8] (6781-6789)
//      vector<CsHist2D*> vRC6790[0-8] (6791-6799)

//- booking in CsRCEventRings::checkRichMom()
//      vector<CsHist2D*> vRC6800[0-9] (6801-6809)
//- booking in CsRCEventRings::checkPhiTheta()
//      vector<CsHist2D*> vRC6850[0-9] (6851-6865)

//- booking in CsRCEventAnalysis::moniCFRefInfR()
//      vector<CsHist2D*> vRC6900[0-37] (6901-6938)

//- booking in CsRCEventAnalysis::checkRingEff
//        vRC7000.push_back( ... ->50
//        vRC7050.push_back( ... ->50

//- booking moved to CsRCMirrors::doAliMirrPhi(): (6/4/04)
//        vRC7100.push_back( ... ->116
//        vRC7300.push_back( ... ->116
//        hRC7450, hRC7451
//        vRC7500.push_back( ... ->116
//        hRC7650, hRC7651
//        vRC7700.push_back( ... ->116
//        hRC7850

//- booking in CsRCPartPhotons::checkLikeReso() {
//      std::vector<CsHist2D*> vRC8000(0-19);
//      std::vector<CsHist2D*> vRC8020(0-19);
//      std::vector<CsHist2D*> vRC8050(0-19);

//- booking in CsRCRing::checkChiReso() {
//      std::vector<CsHist2D*> vRC8070(0-19);

//- booking in CsRCEventPads::checkPPads() {
//      std::vector<CsHist2D*> vRC8100(0-19);

//- booking in CsRCEventPads::setDataPads() {
//      std::vector<CsHist1D*> vRC8120(0-1)->(0-2)
//      std::vector<CsHist1D*> vRC8130(0-3)

//- booking in CsRCRing::setTime() {
//      std::vector<CsHist1D*> vRC8140(0-1)

//- booking in CsRCEventClusters::checkPMTClus() {
//      static std::vector<CsHist2D*> vRC8200(1-20)

//- booking in CsRCEventRings::checkRingFit() {
//      static std::vector<CsHist2D*> vRC8300(1-54)
//- booking in CsRCEventRings::checkRingFit() {
//      static std::vector<CsHist2D*> vRC8400(1-54)
//- booking in CsRCEventAnalysis::checkOptCorr() {
//      static std::vector<CsHist2D*> vRC8500(1-100)
//      static std::vector<CsHist2D*> vRC8600(1-100)
//      static std::vector<CsHist2D*> vRC8700(1-100)
//- booking in CsRCEventAnalysis::singlePhotonCAT() {
//      static std::vector<CsHist2D*> vRC8800(1-100)
//- booking in CsRCEventAnalysis::checkLikeDisTheta() {
//      static std::vector<CsHist2D*> vRC8900(1-20)
//- booking in CsRCEventAnalysis::singlePhotonCAT() {
//      static std::vector<CsHist2D*> vRC8950(1-10)

//    print();

    CsHistograms::SetCurrentPath("/");
//@@---------------------------------

    //if( lprint && bookHis ) {
    if( bookHis ) {
      cout << "RICHONE, CsRCHistos() : Histograms required"
	   << ", level " << levelBk << " : " << kHisBk
	   << " histograms booked" << endl;
      cout << "-------------------------------------------"
	   << "---------------" 
	   << "------------------" << endl;
    } else {
      cout << "RICHONE, CsRCHistos() : NO histograms booked" << endl;
      cout << "--------------------------------------------" << endl;
      cout << "nHist = " << kHisBk << endl;
    }

  }
