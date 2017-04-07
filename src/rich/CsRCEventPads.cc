/*!
   \file    CsRCEventPads.cc
   \------------------------
   \brief   CsRCEventPads class implementation.
   \author  Paolo Schiavon
   \version 0.03
   \date    October 2000, rev. August 2005
*/


    #include <ostream>
    #include <cstdio>

    #include <unistd.h>
    #include <fcntl.h>
    #include <sys/types.h>
    #include <sys/stat.h>

    #include "Coral.h"

    #include "CsDigit.h"
    #include "CsMCDigit.h"
    #include "CsRICH1Detector.h"
    #include "CsErrLog.h"

//----------------------------
    #include "CsRCEventPads.h"
    #include "CsRCPad.h"

    #include "CsRCDetectors.h"
    #include "CsRCEventParticles.h"

    #include "CsRCExeKeys.h"
    #include "CsRCRecConst.h"
    #include "CsRCHistos.h"

    #include "CsRCMirrors.h"
    #include "CsRCParticle.h"
    #include "CsRCUty.h"

    #include "../evmc/CsMCRICH1Hit.h"
    #include "../geom/CsRICH1UpGrade.h"
//----------------------------

  using namespace std;
  using namespace CLHEP;

  CsRCEventPads* CsRCEventPads::instance_ = 0;

//==========================================================================
  CsRCEventPads* CsRCEventPads::Instance() {
//------------------------------------------
    if( instance_ == 0 ) instance_ = new CsRCEventPads();
    return instance_;
  }

//==========================================================================
  CsRCEventPads::CsRCEventPads() { }
//----------------------------------

//==========================================================================
  CsRCEventPads::CsRCEventPads( const CsRCEventPads &evpads ) {
//-------------------------------------------------------------
    cout << "RICHONE : CsRCEventPads CopyConstructor" << endl;
    instance_ = evpads.instance_;
    lPads_ = evpads.lPads_;
    flag_ = evpads.flag_;
  }

//==========================================================================
  CsRCEventPads& CsRCEventPads::operator=( const CsRCEventPads &evpads ) {
//------------------------------------------------------------------------
    if( this != &evpads ) {
      cout << "RICHONE : CsRCEventPads Operator=" << endl;
      instance_ = evpads.instance_;
      lPads_ = evpads.lPads_;
      flag_ = evpads.flag_;
    }
    return ( *this );
  }

//==========================================================================
  void CsRCEventPads::clearEventPads() {
//--------------------------------------
    list<CsRCPad*>::iterator ia;
    for( ia=lPads_.begin(); ia!=lPads_.end(); ia++ ) delete (*ia);
    lPads_.clear();
  }

//==========================================================================
  void CsRCEventPads::clearEventPads( list<CsRCPad*> &lPads ) {
//-------------------------------------------------------------
    list<CsRCPad*>::const_iterator ia;
    for( ia=lPads.begin(); ia!=lPads.end(); ia++ ) delete (*ia);
    lPads.clear();
  }

//==========================================================================
  void CsRCEventPads::print() const {
//-----------------------------------
    cout << endl;
    cout << " Pads : " << endl;
    list<CsRCPad*>::const_iterator ia;
    for( ia=lPads_.begin(); ia!=lPads_.end(); ia++ ) {
      (*ia)->print();
    }
  }

//==========================================================================
  void CsRCEventPads::print( const int &np ) const {
//--------------------------------------------------
    cout << endl;
    cout << " Pads : " << endl;
    list<CsRCPad*>::const_iterator ia;
    int ka = 0;
    for( ia=lPads_.begin(); ia!=lPads_.end(); ia++ ) {
      cout << (*ia) << "  ";
      (*ia)->print();
      if( ka++ >= np ) break;
    }
  }

//==========================================================================
  CsRCEventPads::~CsRCEventPads() {
//---------------------------------
    clearEventPads();
  }



#include "coral_config.h"
#if USE_RFIO
#  include <shift.h>
#endif



//==========================================================================
  bool CsRCEventPads::getEventPads() {
//------------------------------------


//--- Paolo
//    May  2006


      CsRCHistos& hist = CsRCHistos::Ref();
      double xh, yh, wh;

      CsRICH1Detector* rich = CsGeom::Instance()->getRich1Detector();

      CsRCExeKeys *key = CsRCExeKeys::Instance();
      static bool readMyFile = key->readMyFile();
//@@--------------------------------------------

      CsRCDetectors *dets = CsRCDetectors::Instance();
      list<CsRCCathode*> lCathodes = dets->lCathodes();
      static int nCathode = 0;

      static bool firstCall = true;
      if( firstCall ) {
        firstCall = false;
        key->acknoMethod( "CsRCEventPads::getEventPads" );

	if( !readMyFile  &&  !(rich->borasMultiarrayOK()) ) {
	  key->setUsePadThresh( false );
//@@-----------------------------------
	  if( key->UsePadThresh() ) {
	    CsErrLog::Instance()->mes( elWarning,
            "RICH pad thresholds NOT available; cathode thresholds used" );
	    std::cout << " RICHONE, CsRCEventPads::getEventPads() : "
		      << "RICH pad thresholds file NOT available; "
		      << "cathode thresholds used" << std::endl;
	  }
	}

	nCathode = dets->nCathode();
      }

      int nPadEv = 0;

      if( key->MCarloEvent() ) {

        if( readMyFile ) {

//------- Myfile MC event :
//        -----------------
          if( !setMyPads() ) return  false;
//             -----------
	} else {

//------- COMGeant event :
//        ----------------
	  if( !setMCRecPads() ) return  false;
//             --------------
	}
      }

      if( key->DataEvent() ) {

        if( readMyFile ) {

//------- Myfile Data event :
//        -------------------
          if( !setMyPads() ) return  false;
//             -----------
        } else {

//------- Raw DATA event :
//        ----------------
	  if( !setDataPads() ) return  false;
//             -------------
	}
      }


//--- monitoring histograms :
//    -----------------------
      nPadEv = lPads_.size();
      xh = nPadEv;
      if( hist.hRC1521 ) hist.hRC1521->Fill( xh );
//hh                     ------------------------

      if( lPads_.size() == 0 ) return  false;

      checkPPads();
//    ------------

      int nPadCat[nCathode];
      for( int kc=0; kc<nCathode; kc++ ) nPadCat[kc] = 0;
      list<CsRCPad*>::iterator pd;
      double padx0 = dets->ptrToCat( lPads_.front()->ic() )->padx();
      double pady0 = dets->ptrToCat( lPads_.front()->ic() )->pady();
      int ip = 0;
      for( pd=lPads_.begin(); pd!=lPads_.end(); pd++ ) {
	int ic = (*pd)->ic();
	CsRCCathode* cat = dets->ptrToCat( ic );
	double padx = cat->padx();
	int nHCatx = cat->nPadx()/2;
	float facx = padx / padx0;
	float offx = dets->vOffCatW( ic ).x() / padx - nHCatx;
	offx *= facx;
	if( offx >= 0 ) offx += 0.5;
	if( offx <  0 ) offx -= 0.5;
	int ioffx = int( offx );
	double pady = cat->pady();
	int nHCaty = cat->nPady()/2;
	float facy = pady / pady0;
	float offy = dets->vOffCatW( ic ).y() / pady - nHCaty;
	offy *= facy;
	if( offy >= 0 ) offy -= 25. + 0.5;
	if( offy <  0 ) offy += 25. - 0.5;
	int ioffy = int( offy );
        xh = (*pd)->ix() *facx + ioffx + 0.5;
        yh = (*pd)->iy() *facy + ioffy + 0.5;
	if( hist.hRC1520 ) hist.hRC1520->Fill( xh, yh );
//hh                       ----------------------------
	//wh = (*pd)->PH();;
        //if( hist.hRC1520 ) hist.hRC1520->Fill( xh, yh, wh );

        xh = (*pd)->PH();
        yh = ic;
        if( hist.hRC1524 ) hist.hRC1524->Fill( xh, yh );
//hh                       ----------------------------
	nPadCat[ic]++;
	ip++;
      }
      for( int kc=0; kc<nCathode; kc++ ) {
        xh = nPadCat[kc];
        yh = kc;
        if( hist.hRC1522 ) hist.hRC1522->Fill( xh, yh );
//hh                       ----------------------------
      }
//--- conditional prints :
//    --------------------
      int kPrintEventPads = key->kPrintEventPads();
      if( kPrintEventPads == 1 ) {
	cout << endl;
	cout << " Pads : " << nPadEv << endl;
      }

      if( kPrintEventPads == 2 ) print();

      return true;

  }


//===========================================================================
  bool CsRCEventPads::killHalo( int& cath, int& ix, int& iy, double& PH ) {
//-------------------------------------------------------------------------


//- Paolo  -  June 2002.


//- rejection of 'beam halo' pads - myfile :
//  ----------------------------------------
//  provisional!

//- WARNING : NOT adapted to PMT RICH !!!

//- defaults
//  --------
    static float RR = 18.;
    static float dRRout = 4.;
    static float dRRin  = 4.;

    static float RR1, RR2;

//- from rich1.options :
//  --------------------
    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;
      CsRCExeKeys::Instance()->acknoMethod( "CsRCEventPads::killHalo" );

      vector<float> vPar;
      bool boo = CsOpt::Instance()->getOpt( "RICHONE", "killHaloLim", vPar );
      if( boo ) {
        dRRout = vPar[0];
        dRRin  = vPar[1];
      }
      float CFRefInd = CsRCRecConst::Instance()->CFRefInd();
      if( CFRefInd > 1. ) {
	RR = sqrt( 2.* (CFRefInd-1.) );
	RR *= 3300. / 8.;
      }
      RR1 = RR - dRRin;
      RR2 = RR + dRRout;
      //cout << RR << "  " << dRRout << "  " << dRRin << endl;
    }

    bool kill = false;

    int ixc = 0;
    int iyc = 0;
    if( cath ==  3 ) { ixc = -3; iyc = 27; }
    if( cath ==  5 ) { ixc = 74; iyc = 27; }
    if( cath == 10 ) { ixc = -3; iyc = 45; }
    if( cath == 12 ) { ixc = 74; iyc = 45; }
    float rr = sqrt( double( (ix - ixc)*(ix - ixc) + (iy - iyc)*(iy - iyc) ) );

    if( rr >= RR1 && rr <= RR2 ) kill = true;

    return  kill;

  }


//===========================================================================
  bool CsRCEventPads::setThreshold( int& cath, int& ix, int& iy, double& PH )
//---------------------------------------------------------------------------
  {


//- Paolo  -  August 2003.


//- sets software threshold for each pad :
//  --------------------------------------
//  in progress!

//--- WARNING : to be adapted to PMT RICH !!!


    CsRCExeKeys *key = CsRCExeKeys::Instance();
    static bool readMyFile = key->readMyFile();
//@@------------------------------------------

    CsRCRecConst *cons = CsRCRecConst::Instance();
    int kOffHi = cons->kOffHi();

    CsRCHistos& hist = CsRCHistos::Ref();
    double xh, yh, wh;

    CsRICH1Detector* rich = CsGeom::Instance()->getRich1Detector();

    static float *cathPadThre;

    static int nCathode = 0;
    static int nPadx = 0;
    static int nPady = 0;

    static std::vector<CsHist2D*> vRC6400;
    static std::vector<CsHist2D*> vRC6420;
    static std::vector<CsHist2D*> vRC6440;
    static std::vector<CsHist1D*> vRC6460;

    iboramultiarray& padThresh = rich->getThresholds();
    fboramultiarray& padPedest = rich->getPedestals();
    fboramultiarray& padSigma  = rich->getSigmas();
    //for( int i=0; i<16; i++ ) for( int j=0; j<72; j++ )
    //  for( int k=0; k<72; k++ ) std::cout << padPedest[i][j][k] << "  ";
    //std::cout << std::endl;
    //exit(1);

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      key->acknoMethod( "CsRCEventPads::setThreshold" );

      CsRCDetectors *dets = CsRCDetectors::Instance();
      std::list<CsRCCathode*> lCathodes = dets->lCathodes();
      nCathode = dets->nCathode();
      nPadx = lCathodes.front()->nPadx();
      nPady = lCathodes.front()->nPady();

      CsOpt* opt = CsOpt::Instance();
      if( readMyFile ) {
	int myFileRun = CsRCEventParticles::Instance()->myFileRun();
	if( myFileRun == 0  ||  !rich->readBorasTables( myFileRun ) ) {
//                               ----------------------------------
	  key->setUsePadThresh( false );
//@@-----------------------------------
	  CsErrLog::Instance()->mes( elWarning,
          "RICH pad thresholds NOT available; cathode thresholds used" );
	  std::cout << "RICHONE, CsRCEventPads::getEventPads() : "
		    << "RICH pad thresholds file NOT available; "
		    << "cathode thresholds used" << std::endl;
	  return  false;
	}
      } else {
	if( !(rich->borasMultiarrayOK()) ) {
	  key->setUsePadThresh( false );
//@@-----------------------------------
	  CsErrLog::Instance()->mes( elWarning,
          "RICH pad thresholds NOT available; cathode thresholds used" );
	  std::cout << "RICHONE, CsRCEventPads::getEventPads() : "
		    << "RICH pad thresholds file NOT available; "
		    << "cathode thresholds used" << std::endl;
	  return  false;
	}
      }

      if( hist.bookHis()  &&  hist.levelBk() == 2 ) {
	CsHistograms::SetCurrentPath("/RICH");
	int nChan = nPadx * nPady;
	for( int kh=0; kh<nCathode; kh++ ) {
	  int khPlots = kOffHi + 6400 + kh;
	  stringstream name0, title0;
	  name0 << khPlots;
	  title0 << "   ";
	  vRC6400.push_back( new CsHist2D( name0.str(), title0.str(),
//hh      -----------------------------------------------------------
					   nPadx, 0., double( nPadx ),
					   nPady, 0., double( nPady ) ) );
	  //nChan, 0., double( nChan ) ) );
	  khPlots = kOffHi + 6420 + kh;
	  stringstream name1, title1;
	  name1 << khPlots;
	  title1 << "   ";
	  vRC6420.push_back( new CsHist2D( name1.str(), title1.str(),
//hh      -----------------------------------------------------------
					   nPadx, 0., double( nPadx ),
					   nPady, 0., double( nPady ) ) );
	  //nChan, 0., double( nChan ) ) );
	  khPlots = kOffHi + 6440 + kh;
	  stringstream name2, title2;
	  name2 << khPlots;
	  title2 << "   ";
	  vRC6440.push_back( new CsHist2D( name2.str(), title2.str(),
//hh      -----------------------------------------------------------
					   nPadx, 0., double( nPadx ),
					   nPady, 0., double( nPady ) ) );
	  //nChan, 0., double( nChan ) ) );
	  khPlots = kOffHi + 6460 + kh;
	  stringstream name3, title3;
	  name3 << khPlots;
	  title3 << "   ";
	  vRC6460.push_back( new CsHist1D( name3.str(), title3.str(),
//hh      -----------------------------------------------------------
					   100, -10., 10. ) );
	}
	CsHistograms::SetCurrentPath("");
      }

      //cathPadThre = CsRCRecConst::Instance()->cathodeThre();
      cathPadThre = CsRCRecConst::Instance()->cathPadThre();

      for( int kca=0; kca<nCathode; kca++ ) {
	for( int kx=0; kx<nPadx; kx++ ) {
	  for( int ky=0; ky<nPady; ky++ ) {
	    float thresh = padThresh[kca][kx][ky];
	    float pedest = padPedest[kca][kx][ky];
	    float sigma  = padSigma[kca][kx][ky];
	    //std::cout << thresh << "  " << pedest << "  " 
            //          << sigma << std::endl;
	    float thrUse = pedest + 2.5 * sigma;
	    xh = double( kx ) + 0.5;
	    yh = double( ky ) + 0.5;
	    wh = thrUse - thresh;
	    if( vRC6440[kca] ) vRC6440[kca]->Fill( xh, yh, wh );
//hh                           --------------------------------
	    xh = thrUse - thresh;
	    if( vRC6460[kca] ) vRC6460[kca]->Fill( xh );
//hh                           ------------------------
	  }
	}
      }

    }

    bool bAccp = false;
    //std::cout << cath << "  " << ix << "  " << iy 
    //          << "  " << PH << std::endl;

    //iboramultiarray& padThresh = rich->getThresholds();
    //fboramultiarray& padPedest = rich->getPedestals();
    //fboramultiarray& padSigma  = rich->getSigmas();
    //for( int i=0; i<16; i++ ) for( int j=0; j<72; j++ )
    //  for( int k=0; k<72; k++ ) std::cout << padPedest[i][j][k] << "  ";
    //std::cout << std::endl;
    //exit(1);

    float thresh = padThresh[cath][ix][iy];
    float pedest = padPedest[cath][ix][iy];
    float sigma  = padSigma[cath][ix][iy];
    float PHin = PH;
    float ThrUse = 0.;
    if( !(thresh < 0. || pedest < 0. || sigma < 0.) ) {
      PHin = PH + thresh;
      ThrUse = pedest + cathPadThre[cath] * sigma;
    }
    //std::cout << PH << "  " << thresh << "  " << pedest << "  "
    //    << cathPadThre[cath] << "  " << sigma << std::endl;
    float ddPH = PHin - ThrUse;
    if( ddPH  >= 0. ) bAccp = true;
//    ---------------
    float ddTh = ThrUse - thresh;

    static bool bCount = true;
    if( bCount ) {
      //bCount = false;
      //xh = float( (iy * nPadx) + ix ) + 0.5;
      xh = double( ix ) + 0.5;
      yh = double( iy ) + 0.5;
      wh = ddPH;
      //wh = ddPH / sigma;
      //wh = ddTh;
      if( vRC6400[cath] ) vRC6400[cath]->Fill( xh, yh, wh );
//hh                      ---------------------------------
      wh = 1.;
      if( vRC6420[cath] ) vRC6420[cath]->Fill( xh, yh, wh );
//hh                      ---------------------------------
    }

    return  bAccp;

  }


  extern "C" { double rndm_( double& ); }

//==========================================================================
  void CsRCEventPads::genMCRings() {
//----------------------------------

//- Paolo
//  August   2005
//  Very Provisional Tool!

    CsRCRecConst *cons = CsRCRecConst::Instance();
    CsRCDetectors *dets = CsRCDetectors::Instance();
    CsRCEventParticles* evparts = CsRCEventParticles::Instance();
    std::list<CsRCParticle*> lParticles = evparts->lParticles();

    lPads_.clear();

    int kPart = 0;
    std::list<CsRCParticle*>::iterator ia;
    for( ia=lParticles.begin(); ia!=lParticles.end(); ia++ ) {
      CsRCParticle* part = (*ia);
      //if( kPart > 0 ) break;
//@@-----------------------
    Hep3Vector vPoPart0 = part->vPosIn();
    Hep3Vector vDcPart0 = part->vDirIn();
    Hep3Vector vPoPaDetW = dets->vImpDetW( vPoPart0, vDcPart0 );
    static const float ddYY = 204;
    if( fabs( vPoPaDetW.y() ) < 288.+ddYY ) continue;     // NO split rings!
    if( fabs( vPoPaDetW.y() ) < 600.+ddYY  &&
        fabs( vPoPaDetW.x() ) < 600.) continue;           // PMT&APV cat
//@@------------------------------------------------
    double mom = part->mom();
//- assume pion mass !
//  ------------------
    double betaIpo = mom / sqrt( mom*mom + 0.01948 );
    double cosThe = 0.;
    cosThe = 1./ (betaIpo * cons->CFRefIndUV() );
    if( cosThe > 1.) cosThe = 1.;
    double thetaUV = acos( cosThe );
    cosThe = 1./ (betaIpo * cons->CFRefIndVS() );
    if( cosThe > 1.) cosThe = 1.;
    double thetaVS = acos( cosThe );
    //std::cout << "genMCRings : " << kDetPart << "  " << mom
    //          << "  " << vPoPaDetW( 1 ) << "  " << vPoPaDetW( 2 )
    //          << std::endl;
    CsRCCathode* cat = NULL;
    //double xc0 = vPoPaDetW( 1 );
    //double yc0 = vPoPaDetW( 2 );
    int nPhotons = 50;
//@@-----------------
    for( int kf=0; kf<nPhotons; kf++ ) {
      bool inCat = false;
      double ww;
      double thew, phiw;
      double phi = rndm_( ww ) * 6.28318;
      double zPho0 = dets->zEntrWind() + rndm_( ww ) * 
	( dets->zExitWind() - dets->zEntrWind() );
      double xPho0 = vPoPart0.x() + vDcPart0.x()*( zPho0-vPoPart0.z() );
      double yPho0 = vPoPart0.y() + vDcPart0.y()*( zPho0-vPoPart0.z() );
      Hep3Vector vPoPho0( xPho0, yPho0, zPho0 );
      double the = thetaUV;
      Hep3Vector vDcPhoPUV( sin( the )*cos( phi ), sin( the )*sin( phi ),
		  	    cos( the ) );
      vDcPhoPUV = vDcPhoPUV.unit();
      Hep3Vector vDcPho0UV = CsRCUty::Instance()->
	rotfbcc( -1., vDcPart0, vDcPhoPUV, thew, phiw );
      Hep3Vector vPoPhoDetWUV = dets->vImpDetW( vPoPho0, vDcPho0UV );
      double xcUV = vPoPhoDetWUV.x();
      double ycUV = vPoPhoDetWUV.y();
      the = thetaVS;
      Hep3Vector vDcPhoPVS( sin( the )*cos( phi ), sin( the )*sin( phi ),
		 	    cos( the ) );
      vDcPhoPVS = vDcPhoPVS.unit();
      Hep3Vector vDcPho0VS = CsRCUty::Instance()->
	rotfbcc( -1., vDcPart0, vDcPhoPVS, thew, phiw );
      Hep3Vector vPoPhoDetWVS = dets->vImpDetW( vPoPho0, vDcPho0VS );
      double xcVS = vPoPhoDetWVS.x();
      double ycVS = vPoPhoDetWVS.y();
      //std::cout << "genMCRings : " << xcVS << "  " << ycVS
      //       	  << "  " << xcUV << "  " << ycUV << "  " << phi << std::endl;
      for( int kc=0; kc<16; kc++ ) {
	cat = dets->ptrToCat( kc );
	double limXMn = cat->xLimMnW();
	double limYMx = cat->yLimMxW();
	double limYMn = cat->yLimMnW();
	double limXMx = cat->xLimMxW();
	//std::cout << "genMCRings : " << kc << "  " << limXMn 
	//          << "  " << limXMx
	//          << "  " << limYMn << "  " << limYMx << std::endl;
	bool bSave = false;
	bool bPush = true;
	std::list<CsRCPad*>::iterator ip;
	int kPad = 0;
	int kPadx = 0;
	int kPady = 0;
	double PHPad = 1.;
        if( xcVS > limXMn  &&  xcVS < limXMx  &&
	    ycVS > limYMn  &&  ycVS < limYMx ) {
	  if( cat->isPMT() ) {
	    inCat = true;
	    //std::cout << "genMCRingsVS: ------ " << std::endl;
	    double xx = xcVS - dets->vOffCatW( kc ).x() + cat->hCatx();
	    double dxx = cat->ppadx();
	    kPadx = int( xx / dxx );
	    double yy = ycVS - dets->vOffCatW( kc ).y() + cat->hCaty();
	    double dyy = cat->ppady();
	    //std::cout << "genMCRingsVS: " << dxx << " " << dyy << std::endl;
	    kPady = int( yy / dyy );
	    //std::cout << "genMCRingsVS: " << kc << " " << xx << " " << kPadx
	    //	  << " " << yy << " " << kPady << std::endl;
  	    kPad = lPads_.size();
	    bPush = true;
	    for( ip=lPads_.begin(); ip!=lPads_.end(); ip++ ) {
	      if( kc == (*ip)->ic()  &&  kPadx == (*ip)->ix()  &&
		  kPady == (*ip)->iy() ) {
		bPush = false;
		//std::cout << "genMCRings VS : double count" << std::endl;
		break;
	      }
	    }
	    if( bPush ) {
	      lPads_.push_back( new CsRCPad( kPad, kc, kPadx, kPady,
					     PHPad ) );
	      bSave = true;
	    }
	  }
	}
	if( bSave ) continue;
        if( xcUV > limXMn  &&  xcUV < limXMx  &&
	    ycUV > limYMn  &&  ycUV < limYMx ) {
	  //std::cout << "genMCRings : " << kc << "  " << limXMn 
	  //        << "  " << limXMx
	  //        << "  " << limYMn << "  " << limYMx << std::endl;
	  if( !cat->isPMT() ) {
	    inCat = true;
	    //std::cout << "genMCRingsUV: " << std::endl;
	    double xx = xcUV - dets->vOffCatW( kc ).x() + cat->hCatx();
	    double dxx = cat->padx();
	    kPadx = int( xx / dxx );
	    double yy = ycUV - dets->vOffCatW( kc ).y() + cat->hCaty();
	    double dyy = cat->pady();
	    //std::cout << "genMCRingsUV: " << dxx << " " << dyy << std::endl;
	    kPady = int( yy / dyy );
	    //std::cout << "genMCRingsUV: " << kc << " " << xx << " " << kPadx
	    //	  << " " << yy << " " << kPady << std::endl;
	    kPad = lPads_.size();
	    bPush = true;
	    for( ip=lPads_.begin(); ip!=lPads_.end(); ip++ ) {
	      if( kc == (*ip)->ic()  &&  kPadx == (*ip)->ix()  &&
		  kPady == (*ip)->iy() ) {
		bPush = false;
		//std::cout << "genMCRings UV : double count" << std::endl;
		break;
	      }
	    }
	    if( bPush ) {
	      lPads_.push_back( new CsRCPad( kPad, kc, kPadx, kPady,
					     PHPad ) );
	      bSave = true;
	    }
	  }
	}
      }
      //if( !inCat ) std::cout << "genMCRings : out Cat" << std::endl;
    }

    kPart++;
    }

  }


//===========================================================================
  bool CsRCEventPads::setMyPads() {
//---------------------------------


//- Paolo
//  May  2006


    CsRCExeKeys *key = CsRCExeKeys::Instance();
    CsRCDetectors *dets = CsRCDetectors::Instance();

    int bytePerInt = sizeof(int);
    int bytePerFloat = sizeof(float);
    static int fh;
    int iret;
    int pSw;

    int ibuff[10];
    float fbuff[2000];

    int nPadEv = 0;
    static int nCathode = 0;
    static float *cathodeThre;

    //bool dumpBf = true;
    bool dumpBf = false;
    if( dumpBf ) std::cout << std::endl;
    if( dumpBf ) std::cout << "CsRCEventPads::setMyPads" << std::endl;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;
      key->acknoMethod( "CsRCEventPads::setMyPads" );

      nCathode = dets->nCathode();
      cathodeThre = CsRCRecConst::Instance()->cathodeThre();
    }

    fh = CsRCEventParticles::Instance()->myFile();
//- read nPadEv
    pSw = bytePerInt;
    iret = read( fh, ibuff, pSw );
//  -----------------------------
    nPadEv = ibuff[0];
    if( dumpBf ) std::cout << "nPadEv " << nPadEv << std::endl;
    if( nPadEv > 5000 ) {
      std::cout << std::endl;
      std::cout << " RICHONE, CsRCEventPads::setMyPads() :"
      	        << " stop reading myFile, nPadEv > 5000  "
		<< nPadEv << std::endl;
      CsRichOne::Instance()->setEndMyFile( true );
      return  false;
    }

    for( int kPad=0; kPad<nPadEv; kPad++ ) {

      pSw = 3 * bytePerInt;
      iret = read( fh, ibuff, pSw );
//    -----------------------------
//--- read  iCatw from 1 to 16,  iXPadw,iYPadw from 1 to 72
//    -----------------------------------------------------
      int iCatw = ibuff[0];
      int iXPadw = ibuff[1];
      int iYPadw = ibuff[2];
//--- go to iCatw from 0 to 15,  iXPadw,iYPadw from 0 to 71
//    -----------------------------------------------------
//    (myfile-data is compatible with myfile-MC)
      //cout << nPadEv << "  " << iCatw << "  " << iXPadw << "  "
      //     << iYPadw << "  ";
      iCatw--;
      iXPadw--;   iYPadw--;

      CsRCCathode* cath = dets->ptrToCat( iCatw );
      if( cath == NULL ) {
	std::cout << " RICHONE, CsRCEventPads::setMyPads() :"
		  << " stop reading myFile" << std::endl;
	CsRichOne::Instance()->setEndMyFile( true );
	return  false;
      }
//--- APV :
      if( cath->isAPV() ) {

//----- OUT for RUN from 51908 on !!!
/*	int cathode = iCatw;
//----- swap cathode number according to Sergei
	if( cathode == 0) cathode = 6;
	else if( cathode ==  6 ) cathode =  0;
	else if( cathode ==  1 ) cathode =  7;
	else if( cathode ==  7 ) cathode =  1;
	else if( cathode ==  2 ) cathode =  4;
	else if( cathode ==  4 ) cathode =  2;
	else if( cathode ==  3 ) cathode =  5;
	else if( cathode ==  5 ) cathode =  3;
	else if( cathode ==  8 ) cathode = 14;
	else if( cathode == 14 ) cathode =  8;
	else if( cathode ==  9 ) cathode = 15;
	else if( cathode == 15 ) cathode =  9;
	else if( cathode == 10 ) cathode = 12;
	else if( cathode == 12 ) cathode = 10;
	else if( cathode == 11 ) cathode = 13;
	else if( cathode == 13 ) cathode = 11;
	iCatw = cathode;
*/
      }

      pSw = bytePerFloat;
      iret = read( fh, fbuff, pSw );
//---------------------------------
      double PHPack = fbuff[0];
      double PHPadw = 0.;
//--- mod for packing A2, A0, A1 into a float (see->myfile)
      double mod = 300.;
      if( cath->isAPV() ) {
	double A2 = int(PHPack/(mod*mod));
	PHPadw = A2;
      }
      else if( cath->isPMT() ) {
	PHPadw = PHPack;
      }
      else {
	PHPadw = PHPack;
      }

//--- check skip event
//    ----------------
      if( CsRichOne::Instance()->skipMyEvent() ) continue;
//    ---------------------------------------------------


//--- accept pad :
//    ------------
      bool accPad = true;

      if( iCatw >= nCathode ) return false;
      int nPadx = dets->ptrToCat( iCatw )->nPadx();
      if( iXPadw >= nPadx ) return false;
      int nPady = dets->ptrToCat( iCatw )->nPady();
      if( iYPadw >= nPady ) return false;

      bool bSelCath = dets->ptrToCat( iCatw )->isPMT();
      if( key->selPMTonly() ) if( !bSelCath ) accPad = false;
      if( key->selAPVonly() ) if( bSelCath ) accPad = false;
//@@--------------------------------------------------------

      if( !cath->isAPV()  &&  !cath->isPMT() ) {
	if( key->UsePadThresh() ) {
//------- set individual thresholds in Digit PH for each pad : (030825)
//        ----------------------------------------------------
//------- obsolete :
//        ----------
          if( !setThreshold( iCatw, iXPadw, iYPadw, PHPadw ) ) accPad = false;
//             ---------------------------------------------
	}
	else {
//----- set thresholds in Digit PH for each cathode : (011010)
//      ---------------------------------------------
          if( PHPadw < cathodeThre[iCatw] ) accPad = false;
	}
      }

//--- kill 'halo' pads on detector : (0206011)
//    ------------------------------
      if( key->KillHaloPads() ) {
	if( killHalo( iCatw, iXPadw, iYPadw, PHPadw ) ) accPad = false;
//@@        -----------------------------------------
      }

      if( accPad ) {

//----- construct the Pad object :
//      --------------------------
	if( key->MCarloEvent() ) {
  	  int jPad = lPads_.size();
	  lPads_.push_back( new CsRCPad( jPad, iCatw, iXPadw, iYPadw,
		  	  	         PHPadw ) );
	}
	if( key->DataEvent() ) {
	  int jPad = lPads_.size();
	  lPads_.push_back( new CsRCPad( jPad, NULL, iCatw, iXPadw, iYPadw,
					 PHPadw ) );
	}

	if( cath->isPMT() ) lPads_.back()->setPMTpars();
//-----------------------------------------------------

/*
//----- TEST 081111 ------------------------------
	if( cath->isPMT() ) {
	  //bool goPad = true;
	  bool goPad = false;
	  CsRCPad* pad = lPads_.back();
	  for( int kp=0; kp<16; kp++ ) {
	    CsRCPad* pad = lPads_.back();
	    //if( pad->PMTcha() <  4 ) goPad = false;
	    //if( pad->PMTcha() > 11 ) goPad = false;
	    if( pad->PMTcha() ==  5 ) goPad = true;
	    if( pad->PMTcha() ==  6 ) goPad = true;
	    if( pad->PMTcha() ==  9 ) goPad = true;
	    if( pad->PMTcha() == 10 ) goPad = true;
	  }
	  //if( !goPad ) pad->setFlag( false );
	  if( goPad ) pad->setFlag( false );
//@@----------------------------------------
	}
//----- TEST 081111 ------------------------------
*/

	if( cath->isAPV() ) {
	  if( !checkAPVPad( PHPack, mod ) ) lPads_.back()->setFlag( false );
//@@-----------------------------------------------------------------------
	}
	if( cath->isPMT() ) {
	  double Time = PHPadw;
	  if( !checkPMTPPad( Time ) ) lPads_.back()->setFlag( false );
//@@-----------------------------------------------------------------
	}
      }

    }   /* end loop on pads: kPad */

//- skip myFile record trailer
    int nSkip = 1;
    if( key->readGf5() ) nSkip = 2;
    pSw = nSkip * bytePerInt;
    iret = read( fh, ibuff, pSw );
//-------------------------------
    if( dumpBf ) std::cout << "myFile record length " << ibuff[0] 
			   << std::endl << std::endl;


    if( CsRichOne::Instance()->skipMyEvent() ) return  false;
//  --------------------------------------------------------
    return  true;

  }


//===========================================================================
  bool CsRCEventPads::setMCRecPads() {
//------------------------------------


//- Paolo
//  May  2006


    CsRCExeKeys *key = CsRCExeKeys::Instance();
    CsRCDetectors *dets = CsRCDetectors::Instance();

    int nPadEv = 0;
    static int nCathode = 0;
    static float *cathodeThre;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;
      key->acknoMethod( "CsRCEventPads::setMCRecPads" );

      nCathode = dets->nCathode();
      cathodeThre = CsRCRecConst::Instance()->cathodeThre();
    }

//- get digit pointers for this events
//  ----------------------------------
    CsRICH1Detector* rich = CsGeom::Instance()->getRich1Detector();
    list<CsDigit*> mydigits;
    mydigits = rich->getMyDigits();
    //mydigits = dynamic_cast<CsRICH1UpGrade*>(rich)->getMyDigits();
    //cout << endl << " Digits : " << mydigits.size() << endl;

    bool bOnePart = false;
//@@---------------------
    CsMCTrack* mctrackSel = NULL;

    bool bNoNoise = false;
//@@---------------------

//- loop over digits
//  ----------------
    //std::cout << "setMCRecPads ------  " << std::endl;
    list<CsDigit*>::iterator id;
    int kDig = 0;
    for( id=mydigits.begin(); id!=mydigits.end(); id++ ) {

//------------------------------------------------------------------------
//--- suppress detector noise, for background studies (9/01) :
//    (but it is better to use the corresponding OPTION)
//--- WARNING : can not be used from 'myfile'
      if( bNoNoise ) {
	if( dynamic_cast<CsMCDigit*>(*id)->getHits().empty() ) continue;
      }
//------------------------------------------------------------------------

//--- flag pad :
//    ----------
      bool flagPad = true;

      //int cathode, ix, iy;
      //int addr = (*id)->getAddress();
//--- //provisional; what about APV's and PMT's ???
      //double PH = (*id)->getDatum();
      //double PHPack = PH;
      //cathode = rich->getCathode(addr);
      //ix = rich->getPadX(addr);
      //iy = rich->getPadY(addr);
//--- //go to cathode from 0 to 15,  ix, iy from 0 to 71
//    ------------------------------------------------
      //cathode--;
      //ix--;   iy--;

//--- Revised 090216
      int addr = (*id)->getAddress();
      int cathode = rich->getCathode(addr)-1;
      int ix = 0;
      int iy = 0;
      double PH = 0.;
      double PHPack = 0.;
      double PMTT0 = 0.;
      double Time = 0.;
      CsRCCathode* cath = dets->ptrToCat( cathode );
      ix = rich->getPadX(addr)-1;
      iy = rich->getPadY(addr)-1;
      int nData = (*id)->getDataSize();
      double* data = (*id)->getData();
      //std::cout << "setMCRecPads " << cathode << " " << ix << " " << iy
      //	<< "  " << nData << std::endl;
//--- APV :
      if( cath->isAPV() ) {
	double A2 = data[0];
	double A0 = data[2];
	double A1 = data[3];
//----- mod and lim's for packing A2, A0, A1 into a float (PHPack->myfile)
	double mod = 300.;
	double lim2  = 999.;
	double lim01 = mod - 1.;
	if( A2 > lim2 ) A2 = lim2;
	if( A0 > lim01 ) A0 = lim01;
	if( A1 > lim01 ) A1 = lim01;
	PHPack = A2*mod*mod + A0*mod + A1;
	//std::cout << "APV :  ";
	//for( int k=0; k<nData; k++ ) std::cout << data[k] << "  ";
	//std::cout << std::endl;

	if( !checkAPVPad( PHPack, mod ) ) flagPad = false;
//           --------------------------
      }
//--- MAPMT :
      else if( cath->isPMT() ) {
	Time = data[1];
	PH   = Time;
	PHPack = PH;
//----- //Time offset for each PMT channel from CondDB
	bool t0Cal = false;
	t0Cal= rich->T0_flag();
	if( t0Cal  &&  nData >= 3 ) {
	  double T0 = data[2];
	  PMTT0 = T0;
	  //Time += PMTT0;
	  //std::cout << "t0  " << PMTT0 << "  " << Time << std::endl;
	}
	//std::cout << "PMT :  ";
	//for( int k=0; k<nData; k++ ) std::cout << data[k] << "  ";
	//std::cout << std::endl;
	//std::cout << "t0  " << PMTT0 << " " << Time << " "
	//<< PH << std::endl;

	if( !checkPMTPPad( Time ) ) flagPad = false;
//           --------------------
      }
//--- Gassiplex :
      else {
	PH = (*id)->getDatum();
	PHPack = PH;
      }

//--- padMCTime : 040127
//    ------------------
      list<CsMCHit*> padMCHs = dynamic_cast<CsMCDigit*>(*id)->getHits();
      //std::cout << "setMCRecPads  " << cathode;
      //std::cout << " setMCRecPads  " << padMCHs.size() << std::endl;
      list<CsMCHit*>::iterator ih;
      float padMCTime = 0.;
      for( ih=padMCHs.begin(); ih!=padMCHs.end(); ih++ ) {
	float time = (*ih)->getDTime();
	//std::cout << "getEventPads " << time << "  ";
        if( fabs( time ) >= 1. ) {      //   ???
	  padMCTime = time;
	  break;
	}
      }
      //std::cout << std::endl;
      //std::cout << "----- " << padMCTime << std::endl;

//--- Test : 061025
//    -------------
      if( dets->ptrToCat( cathode )->isPMT() ) {
	//std::cout << "setMCRecPads --- " << cathode << std::endl;
        //list<CsMCHit*> padMCHits = dynamic_cast<CsMCDigit*>((*id))->getHits();
	//std::cout << "setMCRecPads  " << padMCHits.size() << std::endl;
        //for( ih=padMCHits.begin(); ih!=padMCHits.end(); ih++ ) {
	  //CsMCRICH1Hit* hit = dynamic_cast<CsMCRICH1Hit*>((*ih));
          //float tH = (*ih)->getDTime();
	  //float xH = (*ih)->getX();
	  //float yH = (*ih)->getY();
	//}
        //std::vector<CsMCRICH1Hit*> padMCHits =
	//dynamic_cast<CsRICH1UpGrade*>(rich)->GetMCHits();
	//std::cout << "setMCRecPads  " << padMCHits.size() << std::endl;
        //for( unsigned int kh=0; kh<padMCHits.size(); kh++ ) {
	//std::cout << "digs " << cathode << " " << ix << " " << iy
	//  << "  " << data[1] << std::endl;
	//int kh = kDig;
	//CsMCRICH1Hit* hit = padMCHits[kh];
	//float tH = hit->getDTime();
	//float xH = hit->getX();
	//float yH = hit->getY();
	//std::cout << "hits " << tH << " " << xH << " " << yH << std::endl;
	//}
      }

//--- accept pad :
//    ------------
      bool accPad = true;

      bool bSelCath = dets->ptrToCat( cathode )->isPMT();
      if( key->selPMTonly() ) if( !bSelCath ) accPad = false;
      if( key->selAPVonly() ) if( bSelCath ) accPad = false;
//@@--------------------------------------------------------

      if( !cath->isAPV()  &&  !cath->isPMT() ) {
	if( key->UsePadThresh() ) {
//------- set individual thresholds in Digit PH for each pad : (030825)
//        ----------------------------------------------------
//------- obsolete :
//        ---------- 
	  if( !setThreshold( cathode, ix, iy, PH ) ) accPad = false;
//             -----------------------------------
	}
	else {
//------- set thresholds in Digit PH for each cathode : (011010)
//        ---------------------------------------------
          if( PH < cathodeThre[cathode] ) accPad = false;
	}
      }

//--- kill 'halo' pads on detector : (0206011)
//    ------------------------------
      if( key->KillHaloPads() ) {
	if( killHalo( cathode, ix, iy, PH ) ) accPad = false;
//@@        -------------------------------
      }

      if( accPad ) {

//----- construct the Pad object :
//      --------------------------

	int jPad = lPads_.size();
	lPads_.push_back( new CsRCPad( jPad, (*id),
				       cathode, ix, iy, PH ) );

	lPads_.back()->setPHPack( PHPack );
//----------------------------------------
	lPads_.back()->setMCTime( padMCTime );
//-------------------------------------------

	if( cath->isPMT() ) lPads_.back()->setPMTpars();
//-----------------------------------------------------

	if( !flagPad ) lPads_.back()->setFlag( false );
//@@--------------------------------------------------
      }

      kDig++;

//------------------------------------------------------------------------
//--- select only one particle, for background studies (9/01) :
//    use together with the corresponding CsMCRICH1Detector.cc code
//    -------------------------------------------------------------
      //if( mctrackSel != NULL )
      //cout << mctrackSel->getParticle()->getGeantNumber() << "  ";
      if( bOnePart ) {
	if( mctrackSel == NULL ) {
	  CsDigit* digit = (*id);
	  list<CsMCHit*> mcHits = dynamic_cast<CsMCDigit*>(digit)->getHits();
	  list<CsMCHit*>::iterator ih;
	  for( ih=mcHits.begin(); ih!=mcHits.end(); ih++ ) {
	    mctrackSel = (*ih)->getMCTrack();
	    //if( mctrackSel != NULL ) cout << "--- " 
	    //<< mctrackSel->getParticle()->getGeantNumber() << "  ";
	    break;
	  }
	}
      }
      //cout << endl;
//------------------------------------------------------------------------
    }

//------------------------------------------------------------------------
    if( bOnePart ) {
      //cout << mctrackSel << endl;
      list<CsRCParticle*> lParticles =
        CsRCEventParticles::Instance()->lParticles();
      list<CsRCParticle*>::iterator ip;
      for( ip=lParticles.begin(); ip!=lParticles.end(); ip++ ) {
	(*ip)->setFlag( false );
//      -----------------------
        CsMCTrack* mctrack = (*ip)->pMCTrack();
        if( mctrack == mctrackSel ) (*ip)->setFlag( true );
//                                  ----------------------
        //if( mctrackSel != NULL ) {
        //  cout << mctrackSel->getParticle()->getGeantNumber()
        //       << "  " << mctrackSel << "  ";
        //  cout << mctrack->getParticle()->getGeantNumber()
        //       << "  " << mctrack << "  " << (*ip)->flag() << endl;
        //}
      }
      //cout << endl;
    }
//------------------------------------------------------------------------

    return  true;

  }


//===========================================================================
  bool CsRCEventPads::setDataPads() {
//-----------------------------------


//- Paolo
//  May  2006


    CsRCExeKeys *key = CsRCExeKeys::Instance();
    CsRCDetectors *dets = CsRCDetectors::Instance();
    CsRCRecConst* cons = CsRCRecConst::Instance();
    CsRCHistos& hist = CsRCHistos::Ref();
    double xh, yh, wh;

    int nPadEv = 0;
    static int nCathode = 0;
    static float *cathodeThre;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;
      key->acknoMethod( "CsRCEventPads::setDataPads" );

      nCathode = dets->nCathode();
      cathodeThre = cons->cathodeThre();
    }

    xh = 20.5;
    if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                     ------------------------

//- get digit pointers for this events
//  ----------------------------------
    CsRICH1Detector* rich = CsGeom::Instance()->getRich1Detector();
    list<CsDigit*> mydigits;
    mydigits = rich->getMyDigits();

//- loop over digits
//  ----------------
    list<CsDigit*>::iterator id;
    int kDig = 0;
    for( id=mydigits.begin(); id!=mydigits.end(); id++ ) {

//--- flag pad :
//    ----------
      bool flagPad = true;

      xh = 21.5;
      if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                       ------------------------

// OLD
//- as from Sasha : (010905)
//  ---------------
//    double *data = (*id)->getData();
//    int cathode = int ( data[0] );
//    int ix = int( data[1] );
//    int iy = int( data[2] );
//    double PH = data[3];

//- NEW, as from Kolosov (051202)
//  ----------------------------- 
      //int addr = (*id)->getAddress();
      //int cathode = rich->getCathode(addr)-1;
      //int ix = rich->getPadX(addr)-1;
      //int iy = rich->getPadY(addr)-1;
      //double PH = (*id)->getDatum();

//- NEW, as from Kolosov (051202) and Paolo (060808)
//  ------------------------------------------------
      int addr = (*id)->getAddress();
      int cathode = rich->getCathode(addr)-1;
      int ix = 0;
      int iy = 0;
      double PH = 0.;
      double PHPack = 0.;
      double PMTT0 = 0.;
      double Time = 0.;
      CsRCCathode* cath = dets->ptrToCat( cathode );
//--- APV :
      if( cath->isAPV() ) {
	int nData = (*id)->getDataSize();
	double* data = (*id)->getData();
	ix = rich->getPadX(addr)-1;
	iy = rich->getPadY(addr)-1;
//----- as from Kolosov (061121)
	double A2 = data[0];
	double A0 = data[2];
	double A1 = data[3];
	PH = A2;
	//A2 = 666.;
	//A0 = 199.;
	//A1 = 199.;
//----- mod and lim's for packing A2, A0, A1 into a float (PHPack->myfile)
	double mod = 300.;
	double lim2  = 999.;
	double lim01 = mod - 1.;
	if( A2 > lim2 ) A2 = lim2;
	if( A0 > lim01 ) A0 = lim01;
	if( A1 > lim01 ) A1 = lim01;
	PHPack = A2*mod*mod + A0*mod + A1;
	//double mmod = mod*mod;
	//double a2 = int(PHPack/mmod);
	//double a0 = int((PHPack - a2*mmod)/mod);
	//double a1 = int(PHPack - a2*mmod - a0*mod);
	//std::cout << PHPack << std::endl;
	//std::cout << A2 << "  " << A0 << "  " << A1 << std::endl;
	//std::cout << a2 << "  " << a0 << "  " << a1 << std::endl;

	if( !checkAPVPad( PHPack, mod ) ) flagPad = false;
//           --------------------------
      }
//--- MAPMT :
      else if( cath->isPMT() ) {
	int nData = (*id)->getDataSize();
	double* data = (*id)->getData();
	ix = rich->getPadX(addr)-1;
	iy = rich->getPadY(addr)-1;
	//---- swap in X - obsolete
	//ix = cath->nPPadx()-1 - ix;
	//---- swap inside PMTs - obsolete
	//swapPMTChs( ix, iy );
        //---- alternate code - obsolete
	//int iPMT = rich->getPadX(addr)-1;
	//int iCha = rich->getPadY(addr)-1;
	//int PPPMx = cath->nPPadx()/cath->nPMTx();
	//ix = (iPMT-1 % cath->nPMTx())* PPPMx + (iCha-1 % PPPMx) + 1;
	//int PPPMy = cath->nPPady()/cath->nPMTy();
	//iy = ((iPMT-1)/cath->nPMTx())* PPPMy + ((iCha-1)/PPPMx) + 1;
	Time = data[1];
	PH   = Time;
	PHPack = PH;
	//---- Time offset for each PMT channel from CondDB
	bool t0Cal = false;
	t0Cal= rich->T0_flag();
	if( t0Cal  &&  nData >= 3 ) {
	  double T0 = data[2];
	  PMTT0 = T0;
	  //Time += PMTT0;
	  //std::cout << "t0  " << PMTT0 << "  " << Time << std::endl;
	}
	//std::cout << "t0  " << PMTT0 << " " << Time << " "
	//<< PH << std::endl;
	if( !checkPMTPPad( Time ) ) flagPad = false;
//           --------------------
      }
//--- Gassiplex :
      else {
	ix = rich->getPadX(addr)-1;
	iy = rich->getPadY(addr)-1;
	PH = (*id)->getDatum();
	PHPack = PH;
	//std::cout << cathode << "  " << ix << "  " << iy << "  "
        //        << PH << std::endl;
      }

//--- VERY OLD chamber (=cathode) swapping : (010905)
//    --------------------------------------
//--- cathode = 15 - cathode
//    ----------------------
      //int kChamber = cathode / 2;
      //int plusOne = cathode - kChamber * 2;
      //int kChamberSw = 7 -  kChamber;
      //cathode = kChamberSw * 2 + 1 - plusOne;
//--- corrected by Sasha in ChipGassiplex.cc : (020517)
//--------------------------------------------

//--- cathode is already from 0 to 15,  ix, iy from 0 to 71
//                                          or from 0 to 47
//    -----------------------------------------------------

//--- accept pad :
//    ------------
      bool accPad = true;

      bool bSelCath = dets->ptrToCat( cathode )->isPMT();
      if( key->selPMTonly() ) if( !bSelCath ) accPad = false;
      if( key->selAPVonly() ) if( bSelCath ) accPad = false;
//@@--------------------------------------------------------

      if( !cath->isAPV()  &&  !cath->isPMT() ) {
	if( key->UsePadThresh() ) {
//------- set individual thresholds in Digit PH for each pad : (030825)
//        ----------------------------------------------------
//------- obsolete :
//        ----------
	  if( !setThreshold( cathode, ix, iy, PH ) ) accPad = false;
//             -----------------------------------
	}
	else {
//------- set thresholds in Digit PH for each cathode : (011010)
//        ---------------------------------------------
          if( PH < cathodeThre[cathode] ) accPad = false;
	}
      }

//--- kill 'halo' pads on detector : (0206011) (never used)
//    ------------------------------
      if( key->KillHaloPads() ) {
	if( killHalo( cathode, ix, iy, PH ) ) accPad = false;
//@@        -------------------------------
      }

      if( accPad ) {

//----- construct the Pad object :
//      --------------------------

	int jPad = lPads_.size();
	lPads_.push_back( new CsRCPad( jPad, (*id),
				       cathode, ix, iy, PH ) );

	lPads_.back()->setPHPack( PHPack );
//----------------------------------------
	if( cath->isPMT() ) lPads_.back()->setPMTT0( PMTT0 );
//----------------------------------------------------------

	if( cath->isPMT() ) lPads_.back()->setPMTpars();
//-----------------------------------------------------

	xh = 22.5;
	if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                         ------------------------

	if( !flagPad ) lPads_.back()->setFlag( false );
//@@--------------------------------------------------
      }

      kDig++;

    }

    if( lPads_.empty() ) {
      xh = 23.5;
      if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                       ------------------------
    }

    xh = 29.5;
    if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                     ------------------------

    return  true;

  }


//===========================================================================
  void CsRCEventPads::checkPPads() {
//----------------------------------


//- Paolo
//  August  2006


    if( !CsRichOne::Instance()->UpRICHJob() ) return;
//  ------------------------------------------------

    if( lPads_.size() == 0 ) return;

    if( !CsRCHistos::Ref().bookHis()  ||
        CsRCHistos::Ref().levelBk() < 2 ) return;
//  --------------------------------------------

    CsRCExeKeys *key = CsRCExeKeys::Instance();
    CsRCDetectors *dets = CsRCDetectors::Instance();
    CsRCHistos& hist = CsRCHistos::Ref();
    double xh, yh, wh;
    CsRCRecConst *cons = CsRCRecConst::Instance();
    int kOffHi = cons->kOffHi();

    int nPadEv = lPads_.size();

    static std::vector<CsHist2D*> vRC8100;
    static int kHH = 0;
    static int kHHMx = 20;
    if( !CsRCHistos::Ref().bookHis() ||
	CsRCHistos::Ref().levelBk() < 2 ) kHHMx = 0;
    static int nCathode = 0;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;
      key->acknoMethod( "CsRCEventPads::checkPPads" );

      nCathode = dets->nCathode();
    }

    list<CsRCPad*>::iterator ip;
    for( ip=lPads_.begin(); ip!=lPads_.end(); ip++ ) {
      //^if( !(*ip)->flag() ) continue;

      CsRCCathode* cath = dets->ptrToCat( (*ip)->ic() );
      if( !cath->isPMT() ) continue;

      int khh = -1;
      int ic = (*ip)->ic();
      for( int kc=0; kc<int(dets->nCatPMT().size()); kc++ ) {
	if( ic == dets->nCatPMT()[kc] ) khh = kc;
      }
      if( khh >= 0  &&  khh < int(dets->nCatPMT().size()) ) {
	xh = (*ip)->ix();
	yh = (*ip)->iy();
	if( hist.vRC1500[khh] ) hist.vRC1500[khh]->Fill( xh, yh );
//      ---------------------------------------------------------
      }

    }

    if( kHH >= kHHMx ) return;
    if( kHH < kHHMx ) {
      int kHist = kOffHi + 8100 + kHH;
      stringstream hN8100;
      hN8100 << kHist;
      string hTitle = "PPads";
      CsHistograms::SetCurrentPath("/RICH");
      vRC8100.push_back( new CsHist2D( hN8100.str(), hTitle,
//hh  ------------------------------------------------------
				       96, 0., 96., 96, 0., 96. ) );
      CsHistograms::SetCurrentPath("/");

      for( ip=lPads_.begin(); ip!=lPads_.end(); ip++ ) {
	//^if( !(*ip)->flag() ) continue;

	CsRCCathode* cath = dets->ptrToCat( (*ip)->ic() );
	if( !cath->isPMT() ) continue;
	int nPPadx = cath->nPPadx();
	int nPPady = cath->nPPady();
	int ic = (*ip)->ic();
	int ix = (*ip)->ix();
	int iy = (*ip)->iy();
	if( ic == dets->nCatPMT()[0] ) { ix += nPPadx; iy += nPPady; }
	if( ic == dets->nCatPMT()[1] ) { iy += nPPady; }
	if( ic == dets->nCatPMT()[2] ) { ix += nPPadx; }
	if( ic == dets->nCatPMT()[3] ) { }
	vRC8100[kHH]->Fill( float( ix ), float( iy ) );
//      ----------------------------------------------
      }
      kHH++;
    }

    return;
  }


//===========================================================================
  void CsRCEventPads::swapPMTChs( int& ix, int& iy ) {
//----------------------------------------------------


//- Paolo
//  August  2006

    int ixPMT = ix / 4;
    int ixPT = ix - ixPMT * 4;
    int iyPMT = iy / 4;
    int iyPT = iy - iyPMT * 4;
    int icPT = ixPT + iyPT * 4;
    int icwPT = 15 - icPT;
    int iywPT = icwPT / 4;
    int ixwPT = icwPT - iywPT * 4;
    ix = ixwPT + ixPMT * 4;
    iy = iywPT + iyPMT * 4;

    return;

  }


//===========================================================================
  bool CsRCEventPads::checkPMTPPad( const double Time ) {
//-------------------------------------------------------

//- Paolo
//  November  2006

    CsRCExeKeys *key = CsRCExeKeys::Instance();
    CsRCRecConst* cons = CsRCRecConst::Instance();
    CsRICH1Detector* rich = CsGeom::Instance()->getRich1Detector();

    static float offset, cutLow, cutHigh;
    static std::vector<CsHist1D*> vRC8120;
    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;
      key->acknoMethod( "CsRCEventPads::checkPMTPPad" );

      //---- Default values
      offset  =    -1681.;
      cutLow  = -1000000.;
      cutHigh =  1000000.;
      CsOpt* opt = CsOpt::Instance();
      bool boo;
      vector<float> vflo;
      boo = opt->CsOpt::getOpt( "RICHONE", "PMTTimeCuts", vflo );
      if( boo ) {
        offset  = vflo[0];
	cutLow  = vflo[1];
	cutHigh = vflo[2];
      }
      if( cons->printConsts() ) {
        std::cout << " RICHONE, CsRCEventPads::checkPMTPPad : " 
             << "  offset = " << offset
             << "  cutLow = " << cutLow
             << "  cutHigh = " << cutHigh << std::endl;
      }

      if( CsRCHistos::Ref().bookHis() && CsRCHistos::Ref().levelBk() >= 2 ) {
	int kOffHi = cons->kOffHi();
        CsHistograms::SetCurrentPath("/RICH");
        int kHist = 0;
        string hTitle;
        kHist = kOffHi + 8121;
        stringstream hN8121;
	hN8121 << kHist;
        hTitle = "PPad Time";
        vRC8120.push_back( new CsHist1D( hN8121.str(), hTitle,
//hh    ------------------------------------------------------
	  	  	  	         200, -100., 100. ) );
        kHist = kOffHi + 8122;
        stringstream hN8122;
	hN8122 << kHist;
        hTitle = "PPad Time Cut";
        vRC8120.push_back( new CsHist1D( hN8122.str(), hTitle,
//hh    ------------------------------------------------------
	  	  	  	         200, -100., 100. ) );
//----- Added 1/12/2010
        kHist = kOffHi + 8123;
        stringstream hN8123;
	hN8123 << kHist;
        hTitle = "PPad abs Time";
        vRC8120.push_back( new CsHist1D( hN8123.str(), hTitle,
//hh    ------------------------------------------------------
	  	  	  	         200, -2000., 2000. ) );
	CsHistograms::SetCurrentPath("/");
      }
    }

    double xh, yh, wh;
    int kh;

    kh = 2;
    if( int( vRC8120.size() ) > kh ) {
      xh = Time;
      if( vRC8120[kh] ) vRC8120[kh]->Fill( xh );
//hh  -----------------------------------------
    }

    bool t0Cal = false;
    if( !key->readMyFile() ) t0Cal = rich->T0_flag();
    if( t0Cal ) offset = 0.;
    //std::cout << "t0Cal  " << t0Cal << std::endl;

    float xTime = Time - offset;
    //std::cout << "t0  " << offset << "  " << Time << "  "
    //          << xTime << std::endl;
    kh = 0;
    if( int( vRC8120.size() ) > kh ) {
      xh = xTime;
      if( vRC8120[kh] ) vRC8120[kh]->Fill( xh );
//hh  -----------------------------------------
    }

//- APPLY CUTS
//  ----------
    if( xTime < cutLow  ||  xTime > cutHigh ) return  false;
//  -------------------------------------------------------

    kh = 1;
    if( int( vRC8120.size() ) > kh ) {
      xh = xTime;
      if( vRC8120[kh] ) vRC8120[kh]->Fill( xh );
//hh  -----------------------------------------
    }

    return true;
  }


//===========================================================================
  bool CsRCEventPads::checkAPVPad( const double PH , const double mod ) {
//-----------------------------------------------------------------------

//- Paolo
//  November  2006

    CsRCExeKeys *key = CsRCExeKeys::Instance();
    CsRCRecConst* cons = CsRCRecConst::Instance();

    static float cutLow, cutHigh;
    static std::vector<CsHist1D*> vRC8130;
    static bool doCuts = false;
    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;
      key->acknoMethod( "CsRCEventPads::checkAPVPad" );

      cutLow  = -1000000.;
      cutHigh =  1000000.;
      doCuts = false;
      CsOpt* opt = CsOpt::Instance();
      bool boo;
      vector<float> vflo;
      boo = opt->CsOpt::getOpt( "RICHONE", "APVAmpliCuts", vflo );
      if( boo ) {
	cutLow = vflo[0];
	cutHigh = vflo[1];
	doCuts = true;
      }
      if( cons->printConsts() ) {
        std::cout << " RICHONE, CsRCEventPads::checkAPVPad :   " 
             << "cutLow = " << cutLow
             << "  cutHigh = " << cutHigh << std::endl;
      }

      if( CsRCHistos::Ref().bookHis() && CsRCHistos::Ref().levelBk() >= 2 ) {
	int kOffHi = cons->kOffHi();
        CsHistograms::SetCurrentPath("/RICH");
        int kHist = 0;
        string hTitle;
        kHist = kOffHi + 8131;
        stringstream hN8131;
	hN8131 << kHist;
        hTitle = "APV A2";
        vRC8130.push_back( new CsHist1D( hN8131.str(), hTitle,
//hh    ------------------------------------------------------
	  	  	  	         300, 0., 300. ) );
        kHist = kOffHi + 8132;
        stringstream hN8132;
	hN8132 << kHist;
        hTitle = "APV A0";
        vRC8130.push_back( new CsHist1D( hN8132.str(), hTitle,
//hh    ------------------------------------------------------
	  	  	  	         150, 0., 300. ) );
        kHist = kOffHi + 8133;
        stringstream hN8133;
	hN8133 << kHist;
        hTitle = "APV A1";
        vRC8130.push_back( new CsHist1D( hN8133.str(), hTitle,
//hh    ------------------------------------------------------
	  	  	  	         150, 0., 300. ) );
        kHist = kOffHi + 8134;
        stringstream hN8134;
	hN8134 << kHist;
        hTitle = "APV A2 Cut";
        vRC8130.push_back( new CsHist1D( hN8134.str(), hTitle,
//hh    ------------------------------------------------------
	  	  	  	         300, 0., 300. ) );
//----- Added 11/12/2010
        kHist = kOffHi + 8135;
        stringstream hN8135;
	hN8135 << kHist;
        hTitle = "APV A1-A0";
        vRC8130.push_back( new CsHist1D( hN8135.str(), hTitle,
//hh    ------------------------------------------------------
	  	  	  	         300, -300., 300. ) );
        kHist = kOffHi + 8136;
        stringstream hN8136;
	hN8136 << kHist;
        hTitle = "APV A2-A1";
        vRC8130.push_back( new CsHist1D( hN8136.str(), hTitle,
//hh    ------------------------------------------------------
	  	  	  	         300, -300., 300. ) );
	CsHistograms::SetCurrentPath("/");
      }
    }

    double mmod = mod*mod;
    double A2 = int(PH/mmod);
    double A0 = int((PH - A2*mmod)/mod);
    double A1 = int(PH - A2*mmod - A0*mod);
    //std::cout << "-- " << PH << std::endl;
    //std::cout << A2 << "  " << A0 << "  " << A1 << std::endl;

    double xh, yh, wh;
    int kh;

    kh = 0;
    if( int( vRC8130.size() ) > kh ) {
      xh = A2;
      if( vRC8130[kh] ) vRC8130[kh]->Fill( xh );
//hh  ---------------------------------------
    }
    kh = 1;
    if( int( vRC8130.size() ) > kh ) {
      xh = A0;
      if( vRC8130[kh] ) vRC8130[kh]->Fill( xh );
//hh  ---------------------------------------
    }
    kh = 2;
    if( int( vRC8130.size() ) > kh ) {
      xh = A1;
      if( vRC8130[kh] ) vRC8130[kh]->Fill( xh );
//hh  ---------------------------------------
    }
    kh = 4;
    if( int( vRC8130.size() ) > kh ) {
      xh = A1 - A0;
      if( vRC8130[kh] ) vRC8130[kh]->Fill( xh );
//hh  ---------------------------------------
    }
    kh = 5;
    if( int( vRC8130.size() ) > kh ) {
      xh = A2 - A1;
      if( vRC8130[kh] ) vRC8130[kh]->Fill( xh );
//hh  ---------------------------------------
    }

//- APPLY CUTS
//  ----------
    if( doCuts ) {
      if( A2 < cutLow  ||  A2 >= cutHigh ) return  false;
//    --------------------------------------------------
      if( !( A0 <= A1  &&  A1 < A2 ) ) return  false;
//    ----------------------------------------------
    }

    kh = 3;
    if( int( vRC8130.size() ) > kh ) {
      xh = A2;
      if( vRC8130[kh] ) vRC8130[kh]->Fill( xh );
//hh  ---------------------------------------
    }

    return true;
  }
