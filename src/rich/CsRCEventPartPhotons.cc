
/*!
   \file    CsRCEventPartPhotons.cc
   \-------------------------------
   \brief   CsRCEventPartPhotons class implementation.
   \author  Paolo Schiavon
   \version 0.01
   \date    1 October 2000
*/


  #include <iostream>
  #include <ostream>
  #include <cstdio>
  #include <cmath>

  #include <CLHEP/Vector/ThreeVector.h>
//-------------------------------------

  #include "CsErrLog.h"

// --------------------------------
  #include "CsRCEventPartPhotons.h"

  #include "CsRichOne.h"

  #include "CsRCPartPhotons.h"
  #include "CsRCPhoton.h"

  #include "CsRCEventClusters.h"
  #include "CsRCCluster.h"

  #include "CsRCDetectors.h"
  #include "CsRCCathode.h"

  #include "CsRCMirrors.h"
  #include "CsRCEventParticles.h"
  #include "CsRCParticle.h"

  #include "CsRCExeKeys.h"
  #include "CsRCRecConst.h"
  #include "CsRCHistos.h"
// ---------------------------

  using namespace std;
  using namespace CLHEP;

  CsRCEventPartPhotons* CsRCEventPartPhotons::instance_ = 0;

//==========================================================================
  CsRCEventPartPhotons* CsRCEventPartPhotons::Instance() {
//--------------------------------------------------------
    if( instance_ == 0 ) instance_ = new CsRCEventPartPhotons();
    return instance_;
  }

//==========================================================================
  CsRCEventPartPhotons::CsRCEventPartPhotons() {
//----------------------------------------------
    lPartPhotons_.clear();
    flagPartPhoSel_ = false;
  }

//==========================================================================
  CsRCEventPartPhotons::CsRCEventPartPhotons( 
//-------------------------------------------
                        const CsRCEventPartPhotons &partphot ) {
    cout << "RICHONE : CsRCEventPartPhotons CopyConstructor" << endl;
    instance_ = partphot.instance_;
    lPartPhotons_ = partphot.lPartPhotons_;
    flag_ = partphot.flag_;
    flagPartPhoSel_ = partphot.flagPartPhoSel_;
  }

//==========================================================================
  void CsRCEventPartPhotons::clearEventPartPhotons() {
//----------------------------------------------------
    list<CsRCPartPhotons*>::iterator ip;
    for( ip=lPartPhotons_.begin(); ip!=lPartPhotons_.end(); ip++ ) delete *ip;
    lPartPhotons_.clear();
  }

//==========================================================================
  void CsRCEventPartPhotons::print() const {
//------------------------------------------
    list<CsRCPartPhotons*>::const_iterator ipf;
    for( ipf=lPartPhotons_.begin(); ipf!=lPartPhotons_.end(); ipf++ ) {
      (*ipf)->print();
    }
  }

//==========================================================================
  CsRCEventPartPhotons::~CsRCEventPartPhotons() {
//-----------------------------------------------
    clearEventPartPhotons();
  }


//==========================================================================  
  void CsRCEventPartPhotons::getEventPartPhotons() {
//--------------------------------------------------


//--- Paolo  -  October 2000


      CsRCExeKeys *key = CsRCExeKeys::Instance();
      CsRCRecConst *cons = CsRCRecConst::Instance();
      CsRCHistos& hist = CsRCHistos::Ref();
      double xh, yh, wh;

      xh = 40.5;
      if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                       ------------------------

      flag_ = false;

//--- LOOP over all PARTICLES:
//    ------------------------
      CsRCEventParticles* parts = CsRCEventParticles::Instance();
      list<CsRCParticle*> lParticles = parts->lParticles();
      list<CsRCParticle*>::iterator ip;
      for( ip=lParticles.begin(); ip!=lParticles.end(); ip++ ) {
	if( !(*ip)->flag() ) continue;

	xh = 41.5;
	if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                         ------------------------

//----- create part-photons object :
//      ----------------------------
        if( CsRCMirrors::Instance()->doSelMirrors( (*ip) ) ) {
//      ----------------------------------------------------

	  int kPaPhot = lPartPhotons_.size();
	  lPartPhotons_.push_back( new CsRCPartPhotons( kPaPhot, (*ip) ) );
//                                 -------------------------------------

	  xh = 44.5;
	  if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                           ------------------------
	  xh = 45.5;
	  if( hist.hRC1000 ) 
	    hist.hRC1000->Fill( xh, lPartPhotons_.back()->lPhotons().size() );
//          -----------------------------------------------------------------
	  if( lPartPhotons_.back()->lPhotons().empty() ) {
	    xh = 42.5;
	    if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                             ------------------------
	  }
	  if( !lPartPhotons_.back()->flag() ) {
	    xh = 43.5;
	    if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                             ------------------------
	  }
	}
      }

      if( !lPartPhotons_.empty() )  flag_ = true;
      if( !flag_ ) {
	if( key->kPrintRejections() == 1 ) {
          std::cout << "RICHONE - getEventPartPhotons : Ev "
                    << CsRichOne::Instance()->kEvent() << "  n-PaPhos "
		    << lPartPhotons_.size() << "  mom "
		    << (*ip)->mom() << std::endl;
	}
      }

      if( lPartPhotons_.empty() ) {
	xh = 46.5;
	if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                         ------------------------
      }

      //^if( key->MCarloEvent() ) bkgrPhotons();
//                             -------------

//--- monitor raw signal on detectors :
//    ---------------------------------
      //^rawSignal();
//    -----------

//--- select 'useful' part-photons :
//    -----------------------------
      bool PartPhoSel = key->PartPhoSelect();
      flagPartPhoSel_ = false;
      if( PartPhoSel )  { flagPartPhoSel_ = true;  partPhoSelection(); }
//                                                 ------------------
      list<CsRCPartPhotons*>::iterator ipf;
      for( ipf=lPartPhotons_.begin(); ipf!=lPartPhotons_.end(); ipf++ ) {
	if( !(*ipf)->flag() ) {
	  xh = 47.5;
	  if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                           ------------------------
	}
      }

//    monitor the C4F10 refractive index :
//    ------------------------------------
      moniCFRefInd();
//    --------------

//--- monitor raw signal :
//    --------------------
      checkSignal();
//    -------------

//--- All photons particle identification :
//    -------------------------------------
      setBackgrType();
//    ---------------
      if( cons->likeType() == "ALL" ) partAllIdent();
//                                    --------------
      else if( cons->ringDefMode() == "MAXLIKE" ) partAllIdent();
//                                                --------------

//--- checks :
//    --------
      bool checkAMPS = true;
      if( checkAMPS ) checkAMPSCorr();
//                    ---------------

      bool checkLIKE = true;
      if( checkLIKE ) checkLikeCorr();
//                     ---------------

//--- conditional prints :
//    --------------------
      int nPaPhot = lPartPhotons_.size();

      int kPrintPartPhotons = key->kPrintPartPhotons();
      if( kPrintPartPhotons == 1 ) {
        cout << endl;
        cout << " PartPhotons : " << nPaPhot << endl;
      }
      if( kPrintPartPhotons >= 2 ) print();


//--- monitoring histograms :
//    -----------------------
      xh = nPaPhot;
      if( hist.hRC1121 ) hist.hRC1121->Fill( xh );
//hh                     ------------------------

      xh = 49.5;
      if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                       ------------------------

      return;

  }


  extern "C" { double rndm_( double& ); }

//===========================================================================
  void CsRCEventPartPhotons::bkgrPhotons() {
//------------------------------------------


//--- Paolo  -  September 2001


      CsRCHistos& hist = CsRCHistos::Ref();
      int kh;
      double xh, yh, wh;

//--- background vs theta MC, evaluation with fit :
//    ---------------------------------------------
      list<CsRCPartPhotons*>::iterator ipf;
      for( ipf=lPartPhotons_.begin(); ipf!=lPartPhotons_.end(); ipf++ ) {
	if( !(*ipf)->flag() ) continue;
	double thetaCer = (*ipf)->pPart()->thetaCer();
	//double thetaCer = (*ipf)->pPart()->mom();           // !!!
	//double ww;
	//double thetaCer = 4. + 54.* rndm_( ww );            // !!!
        list<CsRCPhoton*> lPhotons = (*ipf)->lPhotons();
        list<CsRCPhoton*>::iterator ih;
        if( lPhotons.size() > 1 )  {   
          for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
//--------- theta photon vs theta MC ( Like background ):
//          ---------------------------------------------
            xh = (*ih)->the();
            yh = thetaCer;
            if( hist.hRC3551 ) hist.hRC3551->Fill( xh, yh );
//hh                           ----------------------------
          }   /* end of loop on photons */

//------- counts for theta photon vs theta MC :
//        -------------------------------------
	  xh = thetaCer;
	  if( hist.hRC3550 ) hist.hRC3550->Fill( xh );
//hh                         ------------------------
	}
      }

  }


//===========================================================================
  void CsRCEventPartPhotons::rawSignal() {
//----------------------------------------


//--- Paolo  -  October 2001
//      rev.    December 2002


      CsRCHistos& hist = CsRCHistos::Ref();
      int kh;
      double xh, yh, wh;


//--- raw signal on detector planes :
//    -------------------------------
      list<CsRCPartPhotons*>::iterator ipf;
      for( ipf=lPartPhotons_.begin(); ipf!=lPartPhotons_.end(); ipf++ ) {
	if( !(*ipf)->flag() ) continue;

        list<CsRCPhoton*> lPhotons = (*ipf)->lPhotons();
        list<CsRCPhoton*>::iterator ih;
        if( lPhotons.size() > 1 )  {   
          for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {

	    xh = (*ih)->the();
	    yh = -1.;
	    if( (*ih)->isPMT() ) yh = 0.5;
	    if( (*ih)->isAPV() ) yh = 1.5;
            if( hist.hRC3554 ) hist.hRC3554->Fill( xh, yh );
//hh                           ----------------------------

	    double xPade = (*ipf)->pPart()->vPade()[(*ih)->kDetClu()].x();
	    double yPade = (*ipf)->pPart()->vPade()[(*ih)->kDetClu()].y();
	    double xClu = (*ih)->ptToClu()->xc();
	    double yClu = (*ih)->ptToClu()->yc();
	    double dd = sqrt( (xClu-xPade)*(xClu-xPade) + 
			      (yClu-yPade)*(yClu-yPade) );
            xh = dd;
            if( hist.hRC3555 ) hist.hRC3555->Fill( xh );
//hh                           ------------------------
	    xh = xClu - xPade;
	    yh = yClu - yPade;
            if( hist.hRC3556 ) hist.hRC3556->Fill( xh, yh );
//hh                           ----------------------------
          }   /* end of loop on photons */
	}
      }

  }


//===========================================================================
  void CsRCEventPartPhotons::partPhoSelection() {
//-----------------------------------------------


//--- part-photon selection :
//    -----------------------
//--- Paolo  -  November 2002
//      rev.    December 2002


      CsRCExeKeys *key = CsRCExeKeys::Instance();
      CsRCRecConst *cons = CsRCRecConst::Instance();

      CsRCHistos& hist = CsRCHistos::Ref();
      double xh, yh, wh;

      static bool firstCall = true;
      if( firstCall ) { 
        firstCall = false;

        key->acknoMethod( "CsRCEventPartPhotons::PartPhoSelection" );
      }

      float momMinProc = cons->momMinProc();
      float momMaxProc = cons->momMaxProc();
      int nPhotMinRing = cons->nPhotMinRing();

//--- loop over part-photons :
//    ------------------------
      CsRCEventPartPhotons* paphos = CsRCEventPartPhotons::Instance();
      list<CsRCPartPhotons*> lPaPhotons = paphos->lPartPhotons();
      list<CsRCPartPhotons*>::iterator ipf;

      for( ipf=lPaPhotons.begin(); ipf!=lPaPhotons.end(); ipf++ ) {
	if( !(*ipf)->flag() ) continue;

	CsRCPartPhotons* papho = (*ipf);
	if( !papho ) continue;

        bool bSele = true;

//----- one part-photon only :
	//if( !( lPaPhotons.size() == 1 ) ) bSele = false;
//@@---------------------------------------------------

//----- particle momentum :
	CsRCParticle* part = papho->pPart();
        float momMeas = part->mom();
        if( !( momMeas > momMinProc ) ) bSele = false;
	if( !( momMeas < momMaxProc ) ) bSele = false;
//@@-------------------------------------------------

//----- 'photons' in use :
	std::list<CsRCPhoton*> lPhotons = papho->lPhotons();
        if( int(lPhotons.size()) < nPhotMinRing ) bSele = false;
//@@-----------------------------------------------------------

//----- particle on detectors :
        int kDetPart = papho->kDetPart();
	double xPade = papho->vPoPaDetW()[kDetPart].x();
	//if( !( xPade < 0. ) ) bSele = false;
	//if( !( xPade > 0. ) ) bSele = false;
//@@---------------------------------------
	double yPade = papho->vPoPaDetW()[kDetPart].y();
	//if( !( yPade < 0. ) ) bSele = false;      //   down
	//if( !( yPade > 0. ) ) bSele = false;      //   up
//@@---------------------------------------
	//if( fabs( yPade ) < 600.+204.) bSele = false;     //   cat PMT&APV
//@@------------------------------------------------

	double ddPaDet = part->ddPaDet()[kDetPart];
	//if( !(ddPaDet > 400. ) )  bSele = false;
	//if( !(ddPaDet >   0.  &&  ddPaDet < 100.) )  bSele = false;
	//if( !(ddPaDet > 100.  &&  ddPaDet < 200.) )  bSele = false;
	//if( !(ddPaDet > 200.  &&  ddPaDet < 300.) )  bSele = false;
	//if( !(ddPaDet > 300.  &&  ddPaDet < 400.) )  bSele = false;
	//if( !(ddPaDet > 400.  &&  ddPaDet < 500.) )  bSele = false;
	//if( !(ddPaDet > 500.) )  bSele = false;
//@@-------------------------------------------

//----- particle at entrance window :
	double ppxy;
        ppxy = sqrt( pow( part->vPosIn().x(), 2 ) +
      		     pow( part->vPosIn().y(), 2 ) );
	//if( !(ppxy > 50.) ) bSele = false;
//@@-------------------------------------
	double ttxy;
        ttxy = sqrt( pow( part->vDirIn().x()/part->vDirIn().z(), 2 ) +
      		     pow( part->vDirIn().y()/part->vDirIn().z(), 2 ) );
	double thepa = atan( ttxy );
	ttxy *= 1000.;
	thepa *= 1000.;
	//if( !(ttxy > 100.) ) bSele = false;
	//if( !(thepa < 5.) ) bSele = false;
//@@-------------------------------------

//----- photons on det up or down :
	//list<CsRCPhoton*> lPhotons = papho->lPhotons();
	list<CsRCPhoton*>::iterator ih;
	int nPhoUp = 0;
        int nPhoDw = 0;
        for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
          int kDet = (*ih)->kDetClu();
          if( kDet == 0 ) nPhoUp++;
          if( kDet == 1 ) nPhoDw++;
        }
	//if( !( nPhoUp < nPhoDw ) ) bSele = false;      //   down
	//if( !( nPhoDw < nPhoUp ) ) bSele = false;      //   up
//@@--------------------------------------------

	bool exeTest = false;
	if( exeTest ) {
//------- process selected particles to check MAPMT-APV boundary :
//        --------------------------------------------------------
	  bool partUse = false;
	  CsRCDetectors *dets = CsRCDetectors::Instance();
	  int kDetPart = papho->kDetPart();
          Hep3Vector vPade = part->vPade()[kDetPart];
          //float DD = 240.;
          float DD = 290.;
//@@---------------------
	  int kCa = -1;
	  if( kDetPart == 0 ) {
	    kCa = 3;
	    float yLim = dets->vOffCatW( kCa ).y() + DD;
	    if( vPade.y() > yLim ) partUse = true;
	    float xLim = dets->vOffCatW( kCa ).x() + DD;
	    if( vPade.x() > xLim ) partUse = true;
	    kCa = 5;
	    xLim = dets->vOffCatW( kCa ).x() - DD;
	    if( vPade.x() < xLim ) partUse = true;
	  }
	  if( kDetPart == 1 ) {
	    kCa = 10;
	    float yLim = dets->vOffCatW( kCa ).y() - DD;
	    if( vPade.y() < yLim ) partUse = true;
	    float xLim = dets->vOffCatW( kCa ).x() + DD;
	    if( vPade.x() > xLim ) partUse = true;
	    kCa = 12;
	    xLim = dets->vOffCatW( kCa ).x() - DD;
	    if( vPade.x() < xLim ) partUse = true;
	  }
	  if( !partUse ) bSele = false;
//@@----------------------------------
	}


	if( !bSele ) papho->setFlag( false );
//                   -----------------------
      }

  }


//===========================================================================
  void CsRCEventPartPhotons::checkSignal() {
//------------------------------------------


//- Paolo  -  September 2002
//    rev.    December 2002
//    rev.    September 2005
//    rev.    July 2008


    CsRCRecConst* cons = CsRCRecConst::Instance();
    int kOffHi = cons->kOffHi();
    float CFRefInd = cons->CFRefInd();
    double* massPartv = cons->massPartv();

    static std::vector<CsHist1D*> vRC3390;
    CsRCHistos& hist = CsRCHistos::Ref();
    int kh;
    double xh, yh, wh;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;
      for( int kh=0; kh<4; kh++ ) vRC3390.push_back( NULL );
      if( CsRCHistos::Ref().bookHis() &&
	  CsRCHistos::Ref().levelBk() >= 2 ) {
	CsHistograms::SetCurrentPath("/RICH");
	vRC3390.clear();
	int kHist = 0;
	string hTitle;
        hTitle = "thePho-thePi qua-1";
	stringstream hN3391;
	kHist = kOffHi + 3391;
	hN3391 << kHist;
	vRC3390.push_back( new CsHist1D( hN3391.str(), hTitle,
					 100, -10., 10. ) );
        hTitle = "thePho-thePi qua-2";
	stringstream hN3392;
	kHist = kOffHi + 3392;
	hN3392 << kHist;
	vRC3390.push_back( new CsHist1D( hN3392.str(), hTitle,
					 100, -10., 10. ) );
        hTitle = "thePho-thePi qua-3";
	stringstream hN3393;
	kHist = kOffHi + 3393;
	hN3393 << kHist;
	vRC3390.push_back( new CsHist1D( hN3393.str(), hTitle,
					 100, -10., 10. ) );
        hTitle = "thePho-thePi qua-4";
	stringstream hN3394;
	kHist = kOffHi + 3394;
	hN3394 << kHist;
	vRC3390.push_back( new CsHist1D( hN3394.str(), hTitle,
					 100, -10., 10. ) );
        CsHistograms::SetCurrentPath("/");
      }
    }

//- LOOP over PART-PHOTONS :
//  ------------------------
    CsRCEventPartPhotons* paphos = CsRCEventPartPhotons::Instance();
    list<CsRCPartPhotons*> lPaPhotons = paphos->lPartPhotons();
    list<CsRCPartPhotons*>::iterator ipf;

    for( ipf=lPaPhotons.begin(); ipf!=lPaPhotons.end(); ipf++ ) {
      CsRCPartPhotons* papho = (*ipf);
      if( !papho ) continue;    
      if( !papho->flag() ) continue;
//    -----------------------------

      CsRCParticle* part = papho->pPart();
      float momPart = part->mom();

      int kDetPart = papho->kDetPart();
      double ddPaDet = part->ddPaDet()[kDetPart];

      list<CsRCPhoton*> lPhotons = papho->lPhotons();
      list<CsRCPhoton*>::iterator ih;
      if( lPhotons.size() <= 2 ) continue;
//@@-------------------------------------
      //bool split = false;                                            //
      //for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {         //
      //  if( (*ih)->kDetClu() != kDetPart ) { split = true; break; }  //
      //}                                                              //
      //if( !split ) return;                                           //
//    //-------------------                                            //

//--- check for pion, kaon or proton only
//    -----------------------------------
      int kIpo = 0;
      for( int kPaTy=8; kPaTy<=14; kPaTy+=3 ) {

	double theIpo = (*ipf)->thetaIpo( kPaTy );
	double theIpoVS = (*ipf)->thetaIpoVS( kPaTy );
	if( theIpo <= 0. ) continue;

	for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
          double thePhot = (*ih)->the();

          //if( (*ih)->kDetClu() != kDetPart ) continue;             //
          //if( (*ih)->kDetClu() == kDetPart ) continue;             //
//        -------------------------------------------
	  double theIpow = theIpo;
          if( (*ih)->isPMT() ) theIpow = theIpoVS;
          if( theIpow <= 0. ) continue;

          xh = thePhot - theIpow;
          yh = momPart;
          if( hist.vRC3535[kIpo] ) hist.vRC3535[kIpo]->Fill( xh, yh );
//hh                               ----------------------------------
          xh = thePhot - theIpow;
          yh = ddPaDet;
          if( hist.vRC3545[kIpo] ) hist.vRC3545[kIpo]->Fill( xh, yh );
//hh                               ----------------------------------
//------- 21/07/2008 ---
	  CsRCDetectors *dets = CsRCDetectors::Instance();
	  if( dets->nCatPMT().size() <= 0 ) continue;
	  if( kPaTy == 8 ) {
	    int cClu = (*ih)->ptToClu()->ic();
	    xh = thePhot - theIpow;
	    if( cClu == dets->nCatPMT()[0] ) {
	      int kHH = 0;
	      if( vRC3390[kHH] ) vRC3390[kHH]->Fill( xh );
//hh                             ------------------------
	    }
	    if( cClu == dets->nCatPMT()[1] ) {
	      int kHH = 1;
	      if( vRC3390[kHH] ) vRC3390[kHH]->Fill( xh );
//hh                             ------------------------
	    }
	    if( cClu == dets->nCatPMT()[2] ) {
	      int kHH = 2;
	      if( vRC3390[kHH] ) vRC3390[kHH]->Fill( xh );
//hh                             ------------------------
	    }
	    if( cClu == dets->nCatPMT()[3] ) {
	      int kHH = 3;
	      if( vRC3390[kHH] ) vRC3390[kHH]->Fill( xh );
//hh                             ------------------------
	    }
	  }
//------- 21/07/2008 ---
        }
        xh = float( kIpo ) + 0.5;
        yh = momPart;
        if( hist.vRC3535[3] ) hist.vRC3535[3]->Fill( xh, yh );
//hh                          -------------------------------
        xh = float( kIpo ) + 0.5;
        yh = ddPaDet;
        if( hist.vRC3545[3] ) hist.vRC3545[3]->Fill( xh, yh );
//hh                          -------------------------------
        kIpo++;
      }   /* end of loop on kPaTy */

    }


//-------------------------------------------------------------------
//- 'Correlated' background tests : (050304)
//  -------------------------------
//  WORNING : this tests use the SAME histo.s as the APV tests in
//- CsRCEventRings::singlePhoton() : ensure that those tests are OFF!
//
    bool exeTest = false;
    //bool exeTest = true;
//@@--------------------
//- WARNING : to be adapted to PMT RICH !!!
    if( exeTest ) {
      for( ipf=lPaPhotons.begin(); ipf!=lPaPhotons.end(); ipf++ ) {
        CsRCPartPhotons* papho = (*ipf);
        if( !papho ) continue;    
        if( !papho->flag() ) continue;
//      -----------------------------
	CsRCParticle* part = papho->pPart();
	float momPart = part->mom();
	//if( part->charge() < 0. ) continue;
//@@--------------------------------------
	int kDetPart = papho->kDetPart();
	double ddPaDet = part->ddPaDet()[kDetPart];
	int kPaTy = 8;                     //   assume all pions!
	float massIpo = massPartv[kPaTy];
	double betaIpo = momPart / sqrt( momPart*momPart + massIpo*massIpo );
	double cosTheW = 1./ (betaIpo * CFRefInd);
	if( cosTheW > 1.) continue;
	double focLen = papho->pPart()->pMirPart()->RR() / 2.;
	double rrIpo = acos( cosTheW ) * focLen;
	int caPart = papho->iCaPa();

	// --- check cathode pos.
	CsRCDetectors *dets = CsRCDetectors::Instance();
	std::list<CsRCCathode*> lCathodes = dets->lCathodes();
	//std::list<CsRCCathode*>::iterator ith;
	//for( ith=lCathodes.begin(); ith!=lCathodes.end(); ith++ ) {
	//  if( (*ith)->number() == caPart ) {
	//std::cout << caPart << "  " << (*ith)->vCat0() << std::endl;
        //  }
	//}

	// part counts :
	int khC = caPart;
	int khCs = lCathodes.size() - 1;
	//if( khC < 0  ||  khC > 15 ) continue;
	if( khC < 0  ||  khC > khCs ) continue;
	xh = 0.5;
	yh = 32 + 0.5;
	if( hist.hRC3563 ) hist.hRC3563->Fill( xh, yh );
//hh                       ----------------------------
	//std::list<CsRCPhoton*> lPhotons = papho->lPhotons();
	//std::list<CsRCPhoton*>::iterator ih;
	//for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	//  int caClu = (*ih)->ptToClu()->ic();
	//  if( caClu != caPart ) continue;
	//  int kDetClu = (*ih)->kDetClu();
	//  double xClu = (*ih)->ptToClu()->xc();
	//  double yClu = (*ih)->ptToClu()->yc();
	CsRCEventClusters* evclus = CsRCEventClusters::Instance();
        std::list<CsRCCluster*> lClusters = evclus->lClusters();
        std::list<CsRCCluster*>::iterator cl;
	for( cl=lClusters.begin(); cl!=lClusters.end(); cl++ ) {
	  int caClu = (*cl)->ic();
	  int kDetClu = 0;
	  if( caClu >= 8 ) kDetClu = 1;
	  double xClu = (*cl)->xc();
	  double yClu = (*cl)->yc();
	  double xPade = papho->pPart()->vPade()[kDetClu].x();
	  double yPade = papho->pPart()->vPade()[kDetClu].y();
	  double rr = sqrt( (xClu-xPade)*(xClu-xPade) +
			    (yClu-yPade)*(yClu-yPade) );
	  xh = (rr - rrIpo);
	  yh = float( khC+16 ) + 0.5;
	  if( hist.hRC3563 ) hist.hRC3563->Fill( xh, yh );
//hh                         ----------------------------
	  if( caClu == caPart ) {
	    yh = float( khC ) + 0.5;
	    if( hist.hRC3563 ) hist.hRC3563->Fill( xh, yh );
//hh                           ----------------------------
	  }
	}
      }
    }


  }


//===========================================================================
  void CsRCEventPartPhotons::checkAMPSCorr() {
//--------------------------------------------

//- Paolo  -  July 2008


    if( !CsRichOne::Instance()->UpRICHJob() ) return;
    if( !CsRCHistos::Ref().bookHis()  ||
        CsRCHistos::Ref().levelBk() < 2 ) return;
//  --------------------------------------------

    CsRCRecConst* cons = CsRCRecConst::Instance();
    int kOffHi = cons->kOffHi();
    float CFRefInd = cons->CFRefInd();
    double* massPartv = cons->massPartv();
    CsRCDetectors *dets = CsRCDetectors::Instance();

    static std::vector<CsHist1D*> vRC6630;
    static std::vector<CsHist2D*> vRC6635;
    static std::vector<CsHist1D*> vRC6640;
    static std::vector<CsHist2D*> vRC6650;

    CsRCHistos& hist = CsRCHistos::Ref();
    int kh;
    double xh, yh, wh;
    int kHH;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;
      for( int kh=0; kh<4; kh++ ) vRC6630.push_back( NULL );
      for( int kh=0; kh<4; kh++ ) vRC6635.push_back( NULL );
      for( int kh=0; kh<4; kh++ ) vRC6640.push_back( NULL );
      for( int kh=0; kh<9; kh++ ) vRC6650.push_back( NULL );
      if( CsRCHistos::Ref().bookHis() &&
	  CsRCHistos::Ref().levelBk() >= 2 ) {
	CsHistograms::SetCurrentPath("/RICH");
	vRC6630.clear();
	int kHist = 0;
	string hTitle;
        hTitle = "thePho-thePi JuUp";
	stringstream hN6631;
	kHist = kOffHi + 6631;
	hN6631 << kHist;
	vRC6630.push_back( new CsHist1D( hN6631.str(), hTitle,
					 100, -10., 10. ) );
        hTitle = "thePho-thePi SaUp";
	stringstream hN6632;
	kHist = kOffHi + 6632;
	hN6632 << kHist;
	vRC6630.push_back( new CsHist1D( hN6632.str(), hTitle,
					 100, -10., 10. ) );
        hTitle = "thePho-thePi JuDw";
	stringstream hN6633;
	kHist = kOffHi + 6633;
	hN6633 << kHist;
	vRC6630.push_back( new CsHist1D( hN6633.str(), hTitle,
					 100, -10., 10. ) );
        hTitle = "thePho-thePi SaDw";
	stringstream hN6634;
	kHist = kOffHi + 6634;
	hN6634 << kHist;
	vRC6630.push_back( new CsHist1D( hN6634.str(), hTitle,
					 100, -10., 10. ) );
	vRC6635.clear();
        hTitle = "ID vs mom";
	stringstream hN6636;
	kHist = kOffHi + 6636;
	hN6636 << kHist;
	vRC6635.push_back( new CsHist2D( hN6636.str(), hTitle,
					 7, 0., 7., 50, 0., 50. ) );
	stringstream hN6637;
	kHist = kOffHi + 6637;
	hN6637 << kHist;
	vRC6635.push_back( new CsHist2D( hN6637.str(), hTitle,
					 7, 0., 7., 50, 0., 50. ) );
	stringstream hN6638;
	kHist = kOffHi + 6638;
	hN6638 << kHist;
	vRC6635.push_back( new CsHist2D( hN6638.str(), hTitle,
					 100, -10., 10., 50, 0., 50. ) );
	stringstream hN6639;
	kHist = kOffHi + 6639;
	hN6639 << kHist;
	vRC6635.push_back( new CsHist2D( hN6639.str(), hTitle,
					 100, -10., 10., 50, 0., 50. ) );
	vRC6640.clear();
	stringstream hN6641;
	kHist = kOffHi + 6641;
	hN6641 << kHist;
	vRC6640.push_back( new CsHist1D( hN6641.str(), hTitle,
					 100, -10., 10. ) );
	stringstream hN6642;
	kHist = kOffHi + 6642;
	hN6642 << kHist;
	vRC6640.push_back( new CsHist1D( hN6642.str(), hTitle,
					 100, -10., 10. ) );
	stringstream hN6643;
	kHist = kOffHi + 6643;
	hN6643 << kHist;
	vRC6640.push_back( new CsHist1D( hN6643.str(), hTitle,
					 100, 0., 100. ) );
	stringstream hN6644;
	kHist = kOffHi + 6644;
	hN6644 << kHist;
	vRC6640.push_back( new CsHist1D( hN6644.str(), hTitle,
					 100, 0., 100. ) );
	vRC6650.clear();
	stringstream hN6651;
	kHist = kOffHi + 6651;
	hN6651 << kHist;
	vRC6650.push_back( new CsHist2D( hN6651.str(), hTitle,
					 100, -10., 10., 80, 0., 400. ) );
	stringstream hN6652;
	kHist = kOffHi + 6652;
	hN6652 << kHist;
	vRC6650.push_back( new CsHist2D( hN6652.str(), hTitle,
					 100, -10., 10., 80, 0., 400. ) );
	stringstream hN6653;
	kHist = kOffHi + 6653;
	hN6653 << kHist;
	vRC6650.push_back( new CsHist2D( hN6653.str(), hTitle,
					 100, -10., 10., 80, 0., 400. ) );
	stringstream hN6654;
	kHist = kOffHi + 6654;
	hN6654 << kHist;
	vRC6650.push_back( new CsHist2D( hN6654.str(), hTitle,
					 100, -400., 400., 80, 0., 400. ) );
	stringstream hN6655;
	kHist = kOffHi + 6655;
	hN6655 << kHist;
	vRC6650.push_back( new CsHist2D( hN6655.str(), hTitle,
					 100, -400., 400., 80, 0., 400. ) );
	stringstream hN6656;
	kHist = kOffHi + 6656;
	hN6656 << kHist;
	vRC6650.push_back( new CsHist2D( hN6656.str(), hTitle,
					 100, -400., 400., 80, 0., 400. ) );
	stringstream hN6657;
	kHist = kOffHi + 6657;
	hN6657 << kHist;
	vRC6650.push_back( new CsHist2D( hN6657.str(), hTitle,
					 100, -400., 400., 80, 0., 400. ) );
	stringstream hN6658;
	kHist = kOffHi + 6658;
	hN6658 << kHist;
	vRC6650.push_back( new CsHist2D( hN6658.str(), hTitle,
					 100, -400., 400., 80, 0., 400. ) );
	stringstream hN6659;
	kHist = kOffHi + 6659;
	hN6659 << kHist;
	vRC6650.push_back( new CsHist2D( hN6659.str(), hTitle,
					 100, -400., 400., 80, 0., 400. ) );
        CsHistograms::SetCurrentPath("/");
      }
    }

    static int nTot = 0;
    static int nPion = 0;
    static int nTotp = 0;
    double* probaLK;

//- LOOP over PART-PHOTONS
//  ----------------------
    list<CsRCPartPhotons*>::iterator ipf;
    for( ipf=lPartPhotons_.begin(); ipf!=lPartPhotons_.end(); ipf++ ) {
      if( !(*ipf)->flag() ) continue;
      CsRCPartPhotons* papho = (*ipf);

      CsRCParticle* part = papho->pPart();
      float momPart = part->mom();

      bool bAcc = true;

      double likeBkg = papho->probaLKBgAll();
      probaLK = papho->probaLKAll();
      double likeElec = probaLK[ 2];
      double likeMuon = probaLK[ 5];
      double likePion = probaLK[ 8];
      double likeKaon = probaLK[11];
      double likeProton = probaLK[14];
//--- Determine ID particle
      double likew[6];
      likew[0] = likeBkg;
      likew[1] = likeElec;
      likew[2] = likeMuon;
      likew[3] = likePion;
      likew[4] = likeKaon;
      likew[5] = likeProton;
      int kID = -1;
      bool bID = false;
      for( int kHy=0; kHy<6; kHy++ ) {
	if( likew[kHy] > 0. ) {
	  for( int kp=0; kp<6; kp++ ) {
	    if( kp == kHy ) continue; 
	    if( likew[kp] > likew[kHy] ) { bID = false; break; }
	    bID = true;
	  }
	  if( bID ) { kID = kHy; break; }
	}
      }
      kHH = 0;
      xh = float( kID ) + 0.5;
      yh = momPart;
      if( vRC6635[kHH] ) vRC6635[kHH]->Fill( xh, yh );
//hh                     ----------------------------
      xh = 6.5;
      if( vRC6635[kHH] ) vRC6635[kHH]->Fill( xh, yh );
//hh                     ----------------------------
      kID = -1;
      bID = false;
//--- Select ID pions (analysis way)
      if( likePion > 0.  &&  likePion > likeBkg  &&
	  likePion > likeKaon  &&  likePion > likeProton ) bID = true;
      if( momPart < 8.  &&  likePion < likeElec ) bID = false;
      if( bID ) kID = 3;
//--- Select ID kaons (analysis way)
      bID = false;
      if( likeKaon > 0.  &&  likeKaon > likeBkg  &&
	  likeKaon > likePion  &&  likeKaon > likeProton ) bID = true;
      if( momPart < 8.  &&  likeKaon < likeElec ) bID = false;
      if( bID ) kID = 4;
      kHH = 1;
      if( kID > 0 ) {
	xh = float( kID ) + 0.5;
	yh = momPart;
	if( vRC6635[kHH] ) vRC6635[kHH]->Fill( xh, yh );
//hh                       ----------------------------
      }
      xh = 6.5;
      if( vRC6635[kHH] ) vRC6635[kHH]->Fill( xh, yh );
//hh                     ----------------------------

      list<CsRCPhoton*> lPhotons = papho->lPhotons();
//@@------------------------------------------------
      if( kID == 3  ||  kID == 4 ) {
        if( papho->thetaLikeSet() ) {
	  double theLike = papho->thetaLike();
	  if( theLike > 0. ) {
	    if( lPhotons.size() > 5 ) {
//@@---------------------------------
  	      int kPaTy = -1;
	      if( kID == 3 ) {
	        kPaTy =  8;           // Pion
	        kHH = 2;
	      }
  	      if( kID == 4 ) {
	        kPaTy = 11;           // Kaon
	        kHH = 3;
	      }
	      double theIpo = 0.;
	      if( kPaTy > 0 ) theIpo = papho->thetaIpo( kPaTy );
	      if( theIpo > 0. ) {
	        xh = theLike - theIpo;
	        yh = momPart;
	        if( vRC6635[kHH] ) vRC6635[kHH]->Fill( xh, yh );
//hh                               ----------------------------
	      }
	    }
	  }
        }
      }


//--- Select ID pions
      bool IDpion = false;
      if( likePion > 0.  &&  likePion > likeBkg  &&
	  likePion > likeElec  &&  likePion > likeMuon  &&
	  likePion > likeKaon  &&  likePion > likeProton ) IDpion = true;
      if( !IDpion ) bAcc = false;
//@@  --------------------------
      nTot++;
      if( IDpion ) nPion++;
      if( nTot%10000000 == 0 ) std::cout << "checkAMPSCorr : " << nTot << "  "
				     << nPion << std::endl;

      //list<CsRCPhoton*> lPhotons = papho->lPhotons();
//@@------------------------------------------------
      list<CsRCPhoton*>::iterator ih;
      if( lPhotons.size() <= 2 ) bAcc = false;
//@@-----------------------------------------

      if( !bAcc ) continue;

      //double theReco = ring->the();

//--- check for one among electron, muon, pion, kaon or proton
//    --------------------------------------------------------
      int kIpo = 0;
//    select pion
      int kPaTy = 8;
//@@---------------
      double theIpo = papho->thetaIpo( kPaTy );
      double theIpoVS = papho->thetaIpoVS( kPaTy );
      if( theIpo <= 0. ) continue;

      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
        double theIpow = theIpo;
        if( (*ih)->isPMT() ) theIpow = theIpoVS;
        //double theIpow = theReco;
//@@----------------------------
        if( theIpow <= 0. ) continue;
	if( dets->nCatPMT().size() > 0 ) {
          double thePhot = (*ih)->the();
          int cClu = (*ih)->ptToClu()->ic();
          xh = thePhot - theIpow;
          if( cClu == dets->nCatPMT()[0] ) {
            kHH = 0;
            if( vRC6630[kHH] ) vRC6630[kHH]->Fill( xh );
//hh                           ------------------------
          }
          if( cClu == dets->nCatPMT()[1] ) {
            kHH = 1;
            if( vRC6630[kHH] ) vRC6630[kHH]->Fill( xh );
//hh                           ------------------------
          }
          if( cClu == dets->nCatPMT()[2] ) {
            kHH = 2;
            if( vRC6630[kHH] ) vRC6630[kHH]->Fill( xh );
//hh                           ------------------------
          }
          if( cClu == dets->nCatPMT()[3] ) {
            kHH = 3;
            if( vRC6630[kHH] ) vRC6630[kHH]->Fill( xh );
//hh                           ------------------------
          }
	}
      }
      kIpo++;

    }

    for( ipf=lPartPhotons_.begin(); ipf!=lPartPhotons_.end(); ipf++ ) {
      if( !(*ipf)->flag() ) continue;
      CsRCPartPhotons* papho = (*ipf);

      list<CsRCPhoton*> lPhotons = papho->lPhotons();
      list<CsRCPhoton*>::iterator ih;
      int nPMT = 0;
      int nAPV = 0;
      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	if( (*ih)->isPMT() ) nPMT++;
	if( (*ih)->isAPV() ) nAPV++;
      }

      CsRCParticle* part = papho->pPart();
      float momPart = part->mom();
//--- particle at entrance window :
      double ttxy;
      ttxy = sqrt( pow( part->vDirIn().x()/part->vDirIn().z(), 2 ) +
    		   pow( part->vDirIn().y()/part->vDirIn().z(), 2 ) );
      double thePart = atan( ttxy );
      thePart *= 1000.;

//--- check for pion only
      int kPaTy = 8;
//@@---------------
      double theIpo = papho->thetaIpo( kPaTy );
      //double theIpoVS = papho->thetaIpoVS( kPaTy );
      if( theIpo <= 0. ) continue;

      double likeBkg = papho->probaLKBgAll();
      probaLK = papho->probaLKAll();
      double likeElec = probaLK[ 2];
      double likeMuon = probaLK[ 5];
      double likePion = probaLK[ 8];
      double likeKaon = probaLK[11];
      double likeProton = probaLK[14];

      if( papho->thetaLikeSet() ) {
//----- Select ID pions
        bool IDpion = false;
        if( likePion > 0.  &&  likePion > likeBkg  &&
	    likePion > likeElec  &&  likePion > likeMuon  &&
	    likePion > likeKaon  &&  likePion > likeProton ) IDpion = true;
        if( IDpion ) {

	  double theLike = papho->thetaLike();
	  //std::cout << "checkAMPSCorr: theLike/Ipo " << theLike << "  "
	  //          << theIpo << "  " << thePart << std::endl;
	  if( theLike > 0. ) {

	    if( lPhotons.size() > 5 ) {
//@@-----------------------------------
	      xh = theLike - theIpo;
	      yh = thePart;
	      kHH = 0;
	      if( vRC6650[kHH] ) vRC6650[kHH]->Fill( xh, yh );
//hh                             ----------------------------

	      int nPMT = 0;
	      int nAPV = 0;
	      int nCsI = 0;
	      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
		if( (*ih)->isPMT() ) nPMT++;
		if( (*ih)->isAPV() ) nAPV++;
		nCsI++;
	      }
	      if( !CsRichOne::Instance()->UpRICHJob() ) nPMT = nCsI;

	      if( nAPV == 0 ) {
		kHH = 1;
		if( vRC6650[kHH] ) vRC6650[kHH]->Fill( xh, yh );
//hh                               ----------------------------
	      }
	      if( nPMT == 0 ) {
		kHH = 2;
		if( vRC6650[kHH] ) vRC6650[kHH]->Fill( xh, yh );
//hh                               ----------------------------
	      }
	      int kDetPart = papho->kDetPart();
	      //if( kDetPart == 0 ) {
	      if( kDetPart == 0  &&  nPMT != 0 ) {
	        //xh = 1000.*part->thPade()[kDetPart];
	        //yh = thePart;
	        //kHH = 3;
	        //if( vRC6650[kHH] ) vRC6650[kHH]->Fill( xh, yh );
//hh            //                   ----------------------------
	        xh = 1000.*atan( papho->vDcPaReflW()[kDetPart].x() /
				 papho->vDcPaReflW()[kDetPart].z() );
	        if( papho->vPoPaDetW()[kDetPart].x() > 0. ) {
		  xh -= 87.5;
		  kHH = 3;
		}
	        if( papho->vPoPaDetW()[kDetPart].x() < 0. ) {
		  xh += 87.5;
		  kHH = 4;
		}
	        yh = thePart;
	        if( vRC6650[kHH] ) vRC6650[kHH]->Fill( xh, yh );
//hh                               ----------------------------
	        xh = 1000.*atan( papho->vDcPaReflW()[kDetPart].y() /
				 papho->vDcPaReflW()[kDetPart].z() );
	        xh += 270.5;
	        yh = thePart;
	        kHH = 5;
	        if( vRC6650[kHH] ) vRC6650[kHH]->Fill( xh, yh );
//hh                               ----------------------------
	      }
	      //if( kDetPart == 1 ) {
	      if( kDetPart == 1  &&  nPMT != 0 ) {
	        xh = 1000.*atan( papho->vDcPaReflW()[kDetPart].x() /
				 papho->vDcPaReflW()[kDetPart].z() );
	        if( papho->vPoPaDetW()[kDetPart].x() > 0. ) {
		  xh -= 87.5;
		  kHH = 6;
		}
	        if( papho->vPoPaDetW()[kDetPart].x() < 0. ) {
		  xh += 87.5;
		  kHH = 7;
		}
	        yh = thePart;
	        if( vRC6650[kHH] ) vRC6650[kHH]->Fill( xh, yh );
//hh                               ----------------------------
	        xh = 1000.*atan( papho->vDcPaReflW()[kDetPart].y() /
				 papho->vDcPaReflW()[kDetPart].z() );
	        xh -= 270.5;
	        yh = thePart;
	        kHH = 8;
	        if( vRC6650[kHH] ) vRC6650[kHH]->Fill( xh, yh );
//hh                               ----------------------------
	      }
	    }
	  }
	}
      }

      double theAvePMT = 0.;
      double theAveAPV = 0.;
      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	if( (*ih)->isPMT() ) theAvePMT += (*ih)->the0();
	if( (*ih)->isAPV() ) theAveAPV += (*ih)->the0();
      }
      if( lPhotons.size() < 30 ) continue;
//@@-------------------------------------
      if( nPMT == 0 ) continue;
//@@--------------------------
      if( nAPV == 0 ) continue;
//@@--------------------------
      if( float( nAPV )/float( nPMT ) < 0.3 ) continue;
      if( float( nPMT )/float( nAPV ) < 0.3 ) continue;
//@@--------------------------------------------------
      if( nPMT > 0 ) theAvePMT /= float( nPMT );
      if( nAPV > 0 ) theAveAPV /= float( nAPV );
      double chisqPMT = 0.;
      double chisqAPV = 0.;
      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	if( (*ih)->isPMT() ) {
	  chisqPMT += ((*ih)->the0()-theAvePMT)*((*ih)->the0()-theAvePMT);
	  kHH = 0;
	  xh = (*ih)->the0() - theAvePMT;
	  if( vRC6640[kHH] ) vRC6640[kHH]->Fill( xh );
//hh                         ------------------------
	}
	if( (*ih)->isAPV() ) {
	  chisqAPV += ((*ih)->the0()-theAveAPV)*((*ih)->the0()-theAveAPV);
	  kHH = 1;
	  xh = (*ih)->the0() - theAveAPV;
	  if( vRC6640[kHH] ) vRC6640[kHH]->Fill( xh );
//hh                         ------------------------
	}
      }
      if( nPMT > 0 ) chisqPMT /= float( nPMT );
      if( nAPV > 0 ) chisqAPV /= float( nAPV );
      //std::cout << "checkAMPSCorr: chisqPMT/APV  " << chisqPMT << "  "
      //   	  << chisqAPV << "  " << nPMT << "  " << nAPV << std::endl;
      if( chisqPMT > 0. ) {
        kHH = 2;
        xh = chisqPMT;
        if( vRC6640[kHH] ) vRC6640[kHH]->Fill( xh );
//hh                       ------------------------
      }
      if( chisqAPV > 0. ) {
        kHH = 3;
        xh = chisqAPV;
        if( vRC6640[kHH] ) vRC6640[kHH]->Fill( xh );
//hh                       ------------------------
      }

    }

    return;
  }


//===========================================================================
  void CsRCEventPartPhotons::setBackgrType() {
//--------------------------------------------


//--- Paolo  -  Deceember 2002


    CsRCRecConst* cons = CsRCRecConst::Instance();
    CsOpt* opt = CsOpt::Instance();
    string sflo = "   ";

    if( cons->backgrType() == "02" ) {
      if( !getBackgrParam02() ) {
//         ------------------
	cout << " RICHONE, CsRCEventPartPhotons::setBackgrType() :   ";
	cout << "NO background parametrization (02) available!";
	cout << "  Default values used";
	cout << endl;
	string mess = "RICHONE, CsRCEventPartPhotons::setBackgrType() : ";
	string err = "NO background parametrization available!";
	mess.append( err );
	CsErrLog::Instance()->mes( elError, mess );
	//CsErrLog::Instance()->mes( elFatal, mess );
      }
    }

    else if( cons->backgrType() == "03" ) {
      if( !getBackgrParam03() ) {
//         ------------------
	cout << " RICHONE, CsRCEventPartPhotons::setBackgrType() :   ";
	cout << "NO background parametrization (03) available!";
	cout << "  Default values used";
	cout << endl;
	string mess = "RICHONE, CsRCEventPartPhotons::setBackgrType() : ";
	string err = "NO background parametrization available!";
	mess.append( err );
	CsErrLog::Instance()->mes( elError, mess );
	//CsErrLog::Instance()->mes( elFatal, mess );
      }
    }

    else if( cons->backgrType() == "04" ) {
      if( !getBackgrParam04() ) {
//         ------------------
	cout << " RICHONE, CsRCEventPartPhotons::setBackgrType() :   ";
	cout << "NO background parametrization (04) available!";
	cout << "  Default values used";
	cout << endl;
	string mess = "RICHONE, CsRCEventPartPhotons::setBackgrType() : ";
	string err = "NO background parametrization available!";
	mess.append( err );
	CsErrLog::Instance()->mes( elError, mess );
	//CsErrLog::Instance()->mes( elFatal, mess );
      }
    }

    else if( cons->backgrType() == "05" ) {
      CsRCEventClusters* evclus = CsRCEventClusters::Instance();
      if( evclus->lClusters().size() == 0 ) return;
//    --------------------------------------------
      CsRCCluster* pClu = evclus->lClusters().front();
      double backWgt = 0.;
      if( !pClu->getBackWgt( backWgt ) ) {
	cout << " RICHONE, CsRCEventPartPhotons::setBackgrType() :   ";
	cout << "NO background parametrization (05) available!";
	cout << endl;
	string mess = "RICHONE, CsRCEventPartPhotons::setBackgrType() : ";
	string err = "NO background parametrization available!";
	mess.append( err );
	//CsErrLog::Instance()->mes( elError, mess );
	CsErrLog::Instance()->mes( elFatal, mess );
      }
    }

//map 050204 (Stefano) OBSOLETE!
    else if( cons->backgrType() == "BKGMAP" ) {
      bool getBackgrParamMap = true;
      if( !getBackgrParamMap ) {
//         ------------------
	cout << " RICHONE, CsRCEventPartPhotons::setBackgrType() :   ";
	cout << "NO background map available!";
	cout << "  Default values used";
	cout << endl;
	string mess = "RICHONE, CsRCEventPartPhotons::setBackgrType() : ";
	string err = "NO background map available!";
	mess.append( err );
	CsErrLog::Instance()->mes( elError, mess );
	//CsErrLog::Instance()->mes( elFatal, mess );
      }
    }
//map

    else {
      cout << " RICHONE, CsRCEventPartPhotons::setBackgrType() :   ";
      cout << "NO background parametrization " << cons->backgrType()
	   << " available!";
      //cout << "  Default values used";
      cout << endl;
      string mess = "RICHONE, CsRCEventPartPhotons::setBackgrType() : ";
      string err = "NO background parametrization available!";
      mess.append( err );
      //CsErrLog::Instance()->mes( elError, mess );
      CsErrLog::Instance()->mes( elFatal, mess );
    }

  }


//===========================================================================
  bool CsRCEventPartPhotons::getBackgrParam02() {
//-----------------------------------------------

//--- Paolo  -  December 2002


      CsRCExeKeys *key = CsRCExeKeys::Instance();
      CsRCRecConst *cons = CsRCRecConst::Instance();

      static bool done = false;

      static bool firstCall = true;
      if( firstCall ) {
	firstCall = false;

        key->acknoMethod( "CsRCEventPartPhotons::getBackgrParam02" );

//----- defaults :
//      ----------
//----- Normal Intensity runs (2002)
//      ---------------------
	vPara_.push_back(  0.007 );
	vPara_.push_back( -0.00013 );
	vPara_.push_back( -0.00008 );
	vPara_.push_back(  0.0000017 );

//----- from rich1.options :
//      --------------------
        CsOpt* opt = CsOpt::Instance();
	bool boo = false;
        vector<float> vPar;
        boo = opt->CsOpt::getOpt( "RICHONE", "BackgrParam02", vPar );
	int kP = vPar.size();
        if( boo ) {
	  vPara_.clear();
	  for( int k=0; k<kP; k++ ) vPara_.push_back( vPar[k] );
	  done = true;
        }
	if( done ) {
	  if( cons->printConsts() ) {
	    //cout << endl;
	    cout << " RICHONE, CsRCEventPartPhotons::getBackgrParam02 :  ";
	    for( int k=0; k<kP; k++ ) cout << vPara_[k] << ",  ";
	    cout << endl;
	  }
	}
	else {
	  //cout << endl;
	  //cout << " RICHONE, CsRCEventPartPhotons::getBackgrParam02 :  ";
	  //cout << "NO parameters read, DEFAULT used !!!" << endl;
	  string mess = "RICHONE, CsRCEventPartPhotons::getBackgrParam02 : ";
	  string err = "NO parameters read, DEFAULT used !!!";
	  mess.append( err );
	  CsErrLog::Instance()->mes( elError, mess );
	}
      }

      return done;

  }


//===========================================================================
  bool CsRCEventPartPhotons::getBackgrParam03() {
//-----------------------------------------------

//--- Paolo  -  December 2002


      CsRCExeKeys *key = CsRCExeKeys::Instance();
      CsRCRecConst *cons = CsRCRecConst::Instance();

      static bool done = false;

      static bool firstCall = true;
      if( firstCall ) {
	firstCall = false;

        key->acknoMethod( "CsRCEventPartPhotons::getBackgrParam03" );

//----- defaults :
//      ----------
//----- Normal Intensity runs (2002)
//      ---------------------
	vPara_.push_back(  0.01825 );
	vPara_.push_back(  0.00305 );
	vPara_.push_back(  0.00034 );
	vPara_.push_back( -0.0000045 );

//----- Low Intensity runs
//      ------------------
	//vPara_.push_back( 0.0035 );
	//vPara_.push_back( 0.00075 );
	//vPara_.push_back( 0.00006 );
	//vPara_.push_back( -0.0000008 );

//----- from rich1.options :
//      --------------------
        CsOpt* opt = CsOpt::Instance();
	bool boo = false;
        vector<float> vPar;
        boo = opt->CsOpt::getOpt( "RICHONE", "BackgrParam03", vPar );
	int kP = vPar.size();
        if( boo ) {
	  vPara_.clear();
	  for( int k=0; k<kP; k++ ) vPara_.push_back( vPar[k] );
	  done = true;
        }
	if( done ) {
	  if( cons->printConsts() ) {
	    //cout << endl;
	    cout << " RICHONE, CsRCEventPartPhotons::getBackgrParam03 :  ";
	    for( int k=0; k<kP; k++ ) cout << vPara_[k] << ",  ";
	    cout << endl;
	  }
	}
	else {
	  //cout << endl;
	  //cout << " RICHONE, CsRCEventPartPhotons::getBackgrParam03 :  ";
	  //cout << "NO parameters read, DEFAULT used !!!" << endl;
	  string mess = "RICHONE, CsRCEventPartPhotons::getBackgrParam03 : ";
	  string err = "NO parameters read, DEFAULT used !!!";
	  mess.append( err );
	  CsErrLog::Instance()->mes( elError, mess );
	}
      }

      return done;

  }


//===========================================================================
  bool CsRCEventPartPhotons::getBackgrParam04() {
//-----------------------------------------------

//--- Paolo  -  October 2003


      CsRCExeKeys *key = CsRCExeKeys::Instance();
      CsRCRecConst *cons = CsRCRecConst::Instance();

      static bool done = false;

      static double paraDf[72];

      static bool firstCall = true;
      if( firstCall ) {
	firstCall = false;

        key->acknoMethod( "CsRCEventPartPhotons::getBackgrParam04" );

//----- defaults :
//      ----------
//----- Normal Intensity runs (30580 - 2003)
//      ---------------------
	static double paraDf[72] = {
	  -0.0309219, 0.00158499, -9.02395E-06,
 	   0.0227756, -0.000936618, 8.37554E-06,
	  -0.000370517, 3.42803E-05, -3.22077E-07,
	   1.57102E-06, -2.50466E-07, 2.27486E-09,
	  -0.0238598, 0.00125499, -1.57299E-05,
	   0.0457975, -0.00117096, 8.35176E-06,
	  -0.00156255, 5.44526E-05, -4.84769E-07,
	   1.45287E-05, -5.63393E-07, 5.45586E-09,
	   0.00408633, 5.29808E-05, -1.04503E-06,
	   0.018217, -0.000492701, 5.33781E-06,
	  -0.000521644, 2.8075E-05, -3.59568E-07, 
	   4.52405E-06, -3.1444E-07, 4.23362E-09,
	  -0.0453265, 0.00243912, -2.77743E-05,
	   0.0155947, -0.000744271, 8.72587E-06,
	  -0.000404629, 2.8496E-05, -3.57915E-07,
	   3.59428E-06, -2.78464E-07, 3.67996E-09,
	   0.00921289, -0.000886046, 1.53291E-05, 
	   0.00407539, -4.86444E-05, -2.99829E-07,
	  -7.58437E-05, 2.11586E-06, -1.83189E-09, 
	   8.98362E-07, -3.31804E-08, 2.9512E-10,
	  -0.00897714, 0.000164209, 1.93331E-06,
	   0.00595999, -0.00019179, 1.54872E-06,
	  -0.000168463, 7.07223E-06, -5.6664E-08, 
	   1.96138E-06, -9.3458E-08, 9.08867E-10
	};
	for( int kPa=0; kPa<72; kPa++ ) vPara_.push_back( paraDf[kPa] );

//----- Low Intensity runs
//      ------------------
//      none...


//----- from rich1.options :
//      --------------------
        CsOpt* opt = CsOpt::Instance();
	bool boo = false;
        vector<float> vPar;
        boo = opt->CsOpt::getOpt( "RICHONE", "BackgrParam04", vPar );
	int kP = vPar.size();
        if( boo ) {
	  vPara_.clear();
	  for( int k=0; k<kP; k++ ) vPara_.push_back( vPar[k] );
	  done = true;
        }
	if( done ) {
	  if( cons->printConsts() ) {
	    std::cout << setprecision( 7 );
	    //cout << endl;
	    cout << " RICHONE, CsRCEventPartPhotons::getBackgrParam04 :  ";
	    cout << endl;
	    //int i = 0;
	    //for( int k=0; k<4; k++ ) {
	    //  cout << "      ";
	    //  for( int j=0; j<6; j++ ) {
	    //    cout << vPara_[i] << ",  ";
	    //    i++;
	    //  }
	    //  cout << endl;
	    //}
	    int i = 0;
	    for( int k=0; k<12; k++ ) {
	      cout << "      ";
	      for( int j=0; j<6; j++ ) {
		cout << vPara_[i] << ",  ";
		i++;
	      }
	      cout << endl;
	    }
	    std::cout << setprecision( 2 );
	  }
	}
	else {
	  //cout << endl;
	  //cout << " RICHONE, CsRCEventPartPhotons::getBackgrParam04 :  ";
	  //cout << "NO parameters read, DEFAULT used !!!" << endl;
	  string mess = "RICHONE, CsRCEventPartPhotons::getBackgrParam04 : ";
	  string err = "NO parameters read, DEFAULT used !!!";
	  mess.append( err );
	  CsErrLog::Instance()->mes( elError, mess );
	}
      }

      return done;

  }


//===========================================================================
  void CsRCEventPartPhotons::partAllIdent() {
//-------------------------------------------


//--- Part=photon particle identification
//    -----------------------------------
//--- Paolo  -  December 2002
//    rev. August 2005          


      CsRCExeKeys* key = CsRCExeKeys::Instance();

      static bool firstCall = true;
      if( firstCall ) { 
        firstCall = false;

        key->acknoMethod( "CsRCEventPartPhotons::PartIdent" );

	if( key->likeONLY() ) {
	  std::cout << std::endl;
	  std::cout << "RICHONE, CsRCEventPartPhotons() : WARNING,"
	            << " tests on Likelihood ACTIVE!" << std::endl;
	  std::cout << "------------------------------------------"
	            << "----------------------------" << std::endl;
	}
      }

      CsRCRecConst* cons = CsRCRecConst::Instance();
      int kOffHi = cons->kOffHi();
      int nPhotMinRing = cons->nPhotMinRing();

      int nProb = cons->outBufferSize();

      CsRCHistos& hist = CsRCHistos::Ref();
      double xh, yh, wh;


      double likeDV = cons->likeDefVa();
//@@-----------------------------------

//--- PARTICLE IDENTIFICATION :
//    -------------------------
      int nPIDEv = 0;

//--- loop over part-photons :
//    ------------------------
      bool isFirst = true;
      list<CsRCPartPhotons*> lPartPhotons = lPartPhotons_;
      list<CsRCPartPhotons*>::iterator ipf;
      for( ipf=lPartPhotons.begin(); ipf!=lPartPhotons.end(); ipf++ ) {
        if( !(*ipf)->flag() ) continue;

//----- define SELECTED part-photons :
//      -----------------------------
        bool bSele = true;

//----- enough photons :
//      ----------------
	list<CsRCPhoton*> lPhotons = (*ipf)->lPhotons();
        list<CsRCPhoton*>::iterator ih;
        //^int nPhoRing = lPhotons.size();
	int nPhotons = lPhotons.size();
	if( nPhotons == 0 ) bSele = false;
//@@-------------------------------------

//----- particle momentum over muon threshold (1.91 GeV/c) :
//      ----------------------------------------------------
	CsRCParticle* part = (*ipf)->pPart();
        float momMeas = part->mom();
        //if( momMeas < 1.91 ) bSele = false;
//@@--------------------------------------

	for( int kPro=0; kPro<nProb; kPro++) (*ipf)->setPartProb( kPro,  0.);

//----- process selected part-photons :
//      -------------------------------
        if( bSele ) {

	  if( isFirst ) (*ipf)->setLikeFirst( true );
	  isFirst = false;

//------- compute maximum likelihood and its angle :
//        ------------------------------------------
	  double theRecoLkR = 0.;
	  long double likeMax = -100000.;
	  double theRecoLw = 0.;
	  long double likeMw = likeDV;
	  if( key->doThetaLikeMax()  &&  !key->likeONLY() ) {

    	    //if( (*ipf)->getThetaLikeMax( likeMw, theRecoLw ) ) {
	      //likeMax = likeMw;
	      //theRecoLkR = theRecoLw;
	      //(*ipf)->setThetaLike( theRecoLkR );

	    if( (*ipf)->getThetaLikeMax() ) {
//              -------------------------
	      likeMax = (*ipf)->pLikeMax();
	      theRecoLkR = (*ipf)->thetaLike();
	      //std::cout << theRecoLkR << "  " << momMeas << std::endl;
	    }
	  }

	  (*ipf)->setPartProb( 7, theRecoLkR );
//@@------------------------------------------
	  //cout << "     " << likeMax << "  " << theRecoLkR << endl;

//------- Cerenkov angle :
//        ----------------
	  double theReco = 0.;
	  if( (*ipf)->thetaLikeSet() ) theReco = theRecoLkR;
//@@                                   --------------------

	  static bool checkLk = true;
	  if( !key->doCheckLike()  ||  key->likeONLY() ) checkLk = false;
	  //if( checkLk ) (*ipf)->checkLikelihood( theRecoLkR );

	  int kTyLL = 2;
//@@-------------------
          double probaLK[31];
          double probaLKAPV[31];
          double probaLKPMT[31];
          double derivLK[31];
          double derivLKAPV[31];
          double derivLKPMT[31];
          double probaChi[31];
	  int nPhotAPV = 0;
	  int nPhotPMT = 0;
//------- look for all(5) particles :
//        ---------------------------
          for( int kPaTy=kTyLL; kPaTy<=14; kPaTy+=3 ) {
//        -------------------------------------------

	    probaLK[kPaTy] = likeDV;
	    probaLKAPV[kPaTy] = likeDV;
	    probaLKPMT[kPaTy] = likeDV;
	    derivLK[kPaTy] =  0.;
	    derivLKAPV[kPaTy] =  0.;
	    derivLKPMT[kPaTy] =  0.;
	    probaChi[kPaTy] = 10000000.;

//--------- check particle threshold :
//          --------------------------
	    double theIpo = (*ipf)->thetaIpo( kPaTy );
//- 081016
	    double theIpoUV = (*ipf)->thetaIpoUV( kPaTy );
	    double theIpoVS = (*ipf)->thetaIpoVS( kPaTy );
	    bool compLike = false;
	    if( CsRichOne::Instance()->UpRICHJob() ) {
//----------- ASSUME particle momentum above both (UV ans VS) thresholds
	      if( theIpoUV > 0.  &&  theIpoVS > 0. ) compLike = true;
	    } else {
//----------- Particle momentum above UV threshold
	      ///if( theIpo >= 0. ) {      //   081016
	      if( theIpo > 0. ) compLike = true;
	    }
	    if( compLike ) {

	      //if( theIpoUV > 0. && theIpoVS < 0. ) 
	      //std::cout << "theIpo  " << theIpoUV << "  " << theIpoVS
	      //  	  << std::endl;
//- 081016

//----------- PART-PHOTON LIKELIHOOD :
//            ------------------------
//----------- compute likelihood for the given mass :
//            ---------------------------------------
//----------- compute Likelihood using measured angle and
//            theIpo computed from 'right' refr. index

	      //theIpo = double( - kPaTy );

	      (*ipf)->setLikeONLY( true );
//@@-------------------------------------
	      (*ipf)->setPionONLY( false );
	      if( kPaTy ==  8 ) (*ipf)->setPionONLY( true );
	      (*ipf)->setKaonONLY( false );
	      if( kPaTy == 11 ) (*ipf)->setKaonONLY( true );
//@@-------------------------------------------------------
	      int nPhot = 0;
	      long double pLike = (*ipf)->getLikelihood( theIpo, nPhot );
//                                --------------------------------------
	      probaLK[kPaTy] = pLike;
//@@          ----------------------
	      (*ipf)->setLikeONLY( false );
	      (*ipf)->setPionONLY( false );
	      (*ipf)->setKaonONLY( false );


//----------- compute likelihood DERIVATIVES :
//            -------------------------------- (080506)
	      double pLikeAPV = 0.;
	      double pLikePMT = 0.;
	      double derAPV = 0.;
	      double derPMT = 0.;
	      if( pLike != likeDV ) {
		if( !key->likeONLY() ) {
	          (*ipf)->getLikeDeriv( kPaTy, theIpo, pLike,
//                -------------------------------------------
					derAPV, derPMT,
					pLikeAPV, pLikePMT,
					nPhotAPV, nPhotPMT );
		}
	      }
	      probaLKAPV[kPaTy] = pLikeAPV;
	      probaLKPMT[kPaTy] = pLikePMT;
	      derivLK[kPaTy] = derAPV;
	      derivLKAPV[kPaTy] = derAPV;
              derivLKPMT[kPaTy] = derPMT;

/*	      if( pLike != likeDV ) {
		double DCFRefInd = 0.000025;
		double cosTheWD = 1./ ( betaIpo * (CFRefInd+DCFRefInd) );
		if( cosTheWD > 1.) cosTheWD = 1.;
		double theIpoD = acos( cosTheWD ) * 1000.;
                long double pLikeD = 0.;
//------------- compute Likelihood using measured angle and
//              theIpo normalized to ref. index
		if( !key->likeONLY() ) {
	          pLikeD = (*ipf)->getLikelihood( theIpoD, nPhot );
//                         -----------------------------------------
		  double likeDer = (pLikeD - pLike) / DCFRefInd;
		  derivLK[kPaTy] = likeDer;
		}
		//cout << kPaTy << "  " << pLike << "  " << pLikeD << "  "
		//     << likeDer << endl;
	      }*/


//----------- PART-PHOTON QSQUARE :
//            ---------------------
//----------- compute qsquare for the given mass :
//            ------------------------------------
//----------- compute qsquare using measured angle and
//            theIpo computed from 'right' refr. index

	      //theIpo = double( - kPaTy );      // pro-memoria

	      double qSqua = (*ipf)->getQsquare( theIpo, nPhot );
//                           -----------------------------------
	      //probaChi[kPaTy] = exp( - qSqua);
	      probaChi[kPaTy] = qSqua;

            }   // end if on theIpo >= 0.

          }   // end loop on particle type: kPaTy


//------- set variables to part-photon :
//        ------------------------------
	  (*ipf)->setProbaLKAll( &probaLK[0] );
	  (*ipf)->setDerivLKAll( &derivLK[0] );
	  (*ipf)->setProbaChiAll( &probaChi[0] );
	  (*ipf)->setProbaLKAllAPV( &probaLKAPV[0] );
	  (*ipf)->setProbaLKAllPMT( &probaLKPMT[0] );
	  (*ipf)->setDerivLKAllAPV( &derivLKAPV[0] );
	  (*ipf)->setDerivLKAllPMT( &derivLKPMT[0] );
//@@------------------------------------------------

//------- define partProbs :
//        ------------------
	  (*ipf)->setPartProb( 1,  probaLK[ 8] );
	  (*ipf)->setPartProb( 2,  probaLK[11] );
	  (*ipf)->setPartProb( 3,  probaLK[14] );
	  (*ipf)->setPartProb( 4,  derivLK[ 8] );
	  (*ipf)->setPartProb( 5,  derivLK[11] );
	  (*ipf)->setPartProb( 6,  derivLK[14] );
//@@      --------------------------------------
	  if( nProb > 15 ) {
    	    (*ipf)->setPartProb( 15, probaLK[ 2] );
	    (*ipf)->setPartProb( 16, probaLK[ 5] );
	    (*ipf)->setPartProb( 17, derivLK[ 2] );
	    (*ipf)->setPartProb( 18, derivLK[ 5] );
//@@        --------------------------------------
	  }
	  //^(*ipf)->setPartProb(  9, double( nPhoRing ) );
	  (*ipf)->setPartProb(  9, double( nPhotons ) );
	  (*ipf)->setPartProb( 12,  probaChi[ 8] );
	  (*ipf)->setPartProb( 13,  probaChi[11] );
	  (*ipf)->setPartProb( 14,  probaChi[14] );
//@@      ----------------------------------------
	  if( nProb > 15 ) {
 	    (*ipf)->setPartProb( 19,  probaChi[ 2] );
	    (*ipf)->setPartProb( 20,  probaChi[ 5] );
//@@        ----------------------------------------
	  }
//------- (080506)
	  if( nProb > 21 ) {
 	    (*ipf)->setPartProb( 23,  derivLKPMT[ 8] );
	    (*ipf)->setPartProb( 24,  derivLKPMT[11] );
 	    (*ipf)->setPartProb( 25,  derivLKPMT[14] );
	    (*ipf)->setPartProb( 26,  derivLKPMT[ 2] );
	    (*ipf)->setPartProb( 27,  derivLKPMT[ 5] );
	    (*ipf)->setPartProb( 28, double( nPhotAPV ) );
	    (*ipf)->setPartProb( 29, double( nPhotPMT ) );
 	    //# (*ipf)->setPartProb( 30,  probaLKPMT[ 8] );
	    //# (*ipf)->setPartProb( 31,  probaLKPMT[11] );
 	    //# (*ipf)->setPartProb( 32,  probaLKPMT[14] );
	    //# (*ipf)->setPartProb( 33,  probaLKPMT[ 2] );
	    //# (*ipf)->setPartProb( 34,  probaLKPMT[ 5] );
//          ---------------------------------------------
	  }
          //std::cout << "lk    " << probaLK[ 8] << "  " << probaLK[11]
	  //          << "  " << probaLK[14] << std::endl;
	  //std::cout << "chi   " << probaChi[ 8] << "  " << probaChi[11]
	  //          << "  " << probaChi[14] << std::endl;

//------- compute likelihood of 'background' :
//        ------------------------------------
	  int nPhotBk = 0;
	  double theIpoBk = 0.;
	  long double  pLikeBk = (*ipf)->getLikelihood( theIpoBk, nPhotBk );
//                               ------------------------------------------
	  //cout << "     " << pLikeBk << endl;


	  (*ipf)->setProbaLKBgAll( pLikeBk );
//@@----------------------------------------
	  (*ipf)->setPartProb( 0, pLikeBk );
//@@---------------------------------------

	  //^(*ipf)->setNPhotAll( nPhoRing );
	  (*ipf)->setNPhotAll( nPhotons );
//------- (080506)
	  (*ipf)->setNPhotAllAPV( nPhotAPV );
	  (*ipf)->setNPhotAllPMT( nPhotPMT );
//@@----------------------------------------

	  (*ipf)->setPartProbsSet();
//@@-------------------------------

	  if( checkLk ) (*ipf)->checkLikeReso();
//                      -----------------------
	  if( checkLk ) (*ipf)->getLikeProb();
//                      ---------------------

	  (*ipf)->setLikeFirst( false );
	}

      }

  }


//===========================================================================
  void CsRCEventPartPhotons::moniCFRefInd() {
//-------------------------------------------


//- monitor for C4F10 refractive Index
//  ----------------------------------
//- Paolo  -  May 2005, rev. August 2005


    CsRCExeKeys *key = CsRCExeKeys::Instance();
    CsRCHistos& hist = CsRCHistos::Ref();
    double xh, yh, wh;

    CsOpt* opt = CsOpt::Instance();
    static float momLimit = 50.;
//@@---------------------------

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      key->acknoMethod( "CsRCEventPartPhotons::moniCFRefInd" );

      vector<float> vPar;
      bool boo = opt->CsOpt::getOpt( "RICHONE", "n-1Limits", vPar );
      if( boo ) momLimit  = vPar[0];

      CsRCRecConst* cons = CsRCRecConst::Instance();
      double Nminus1VS = (cons->CFRefIndVS() - 1.)* 100.;
      if( hist.hRC1556 ) 
	hist.hRC1556->Fill( 0.000205, 0.5, Nminus1VS );
//hh    ----------------------------------------------
      double Nminus1UV = (cons->CFRefIndUV() - 1.)* 100.;
      if( hist.hRC1557 ) 
	hist.hRC1557->Fill( 0.000205, 0.5, Nminus1UV );
//hh    ----------------------------------------------
    }

    list<CsRCPartPhotons*>::iterator ipf;
    for( ipf=lPartPhotons_.begin(); ipf!=lPartPhotons_.end(); ipf++ ) {
      if( !(*ipf)->flag() ) continue;

      CsRCParticle* part = (*ipf)->pPart();
      if( part->mom() > momLimit ) continue;
//@@---------------------------------------

      if( key->MCarloEvent() ) {
	if( part->iPartT() != 8  &&  part->iPartT() != 9 ) return;
//@@-------------------------------------------------------------
      }

      list<CsRCPhoton*> lPhotons = (*ipf)->lPhotons();
      list<CsRCPhoton*>::iterator ih;
      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
        if( !(*ih)->flag() ) continue;

        double thePho = (*ih)->the0();

//----- ( n-1 ) as 1/(beta*cos(theta-photon)) - 1
//      -----------------------------------------
//      assume PION mass
//      ----------------
        float nn = fabs( part->mom() * cos( thePho/1000. ) );
        float Nminus1 = 0.;
        if( nn > 0. ) {
//------- Pion mass
          Nminus1 = sqrt( 0.01948 + part->mom()*part->mom() ) / nn;
//------- Muon mass (test)
          //Nminus1 = sqrt( 0.011172 + part->mom()*part->mom() ) / nn;
          Nminus1 = Nminus1 - 1.;
        }
        xh = Nminus1;
	yh = part->mom();
	if( (*ih)->isPMT() ) {
	  if( hist.hRC1556 ) hist.hRC1556->Fill( xh, yh );
//hh                         ----------------------------
	  if( hist.hRC1558 ) hist.hRC1558->Fill( xh );
//hh                         ------------------------
	} else {
	  if( hist.hRC1557 ) hist.hRC1557->Fill( xh, yh );
//hh                         ----------------------------
	  if( hist.hRC1559 ) hist.hRC1559->Fill( xh );
//hh                         ------------------------
	}
        int ic = (*ih)->ptToClu()->ic();
	if( hist.vRC1600[ic] ) hist.vRC1600[ic]->Fill( xh );
      }

    }

    return;
  }


//===========================================================================
  void CsRCEventPartPhotons::checkLikeCorr() {
//--------------------------------------------

//- Paolo  -  October  2008


    if( !CsRCHistos::Ref().bookHis()  ||
        CsRCHistos::Ref().levelBk() < 2 ) return;
//  --------------------------------------------

    CsRCExeKeys* key = CsRCExeKeys::Instance();
    CsRCRecConst* cons = CsRCRecConst::Instance();

    CsRCHistos& hist = CsRCHistos::Ref();
    int kOffHi = cons->kOffHi();
    static std::vector<CsHist2D*> vRC6750;
    static std::vector<CsHist2D*> vRC6760;
    static std::vector<CsHist2D*> vRC6770;
    static std::vector<CsHist2D*> vRC6780;
    static std::vector<CsHist2D*> vRC6790;
    double xh, yh, wh;

    bool print = false;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      for( int kh=0; kh<10; kh++ ) vRC6750.push_back( NULL );
      for( int kh=0; kh<10; kh++ ) vRC6760.push_back( NULL );
      for( int kh=0; kh<10; kh++ ) vRC6770.push_back( NULL );
      for( int kh=0; kh<10; kh++ ) vRC6780.push_back( NULL );
      for( int kh=0; kh<10; kh++ ) vRC6790.push_back( NULL );
      if( CsRCHistos::Ref().bookHis() &&
	  CsRCHistos::Ref().levelBk() >= 2 ) {
	CsHistograms::SetCurrentPath("/RICH");
	vRC6750.clear();
	for( int kh=0; kh<10; kh++ ) {
	  string hTitle = " ";
	  stringstream hN6750;
	  int kHist = kOffHi + 6750 + kh + 1;
	  hN6750 << kHist;
	  vRC6750.push_back( new CsHist2D( hN6750.str(), hTitle,
					   100, -0.1,  0.1, 50, 0., 50. ) );
	                                   //200, -2.,  2., 50, 0., 50. ) );
	}
	vRC6760.clear();
	for( int kh=0; kh<5; kh++ ) {
	  string hTitle = " ";
	  stringstream hN6760;
	  int kHist = kOffHi + 6760 + kh + 1;
	  hN6760 << kHist;
	  vRC6760.push_back( new CsHist2D( hN6760.str(), hTitle,
					   //100, -0.1,  0.1, 50, 0., 50. ) );
	                                   200, -2.,  2., 50, 0., 2. ) );
	}
	for( int kh=5; kh<9; kh++ ) {
	  string hTitle = " ";
	  stringstream hN6760;
	  int kHist = kOffHi + 6760 + kh + 1;
	  hN6760 << kHist;
	  vRC6760.push_back( new CsHist2D( hN6760.str(), hTitle,
					   //100, -0.1,  0.1, 50, 0., 50. ) );
	                                   200, -1.,  1., 50, 0., 50. ) );
	}
	vRC6770.clear();
	for( int kh=0; kh<5; kh++ ) {
	  string hTitle = " ";
	  stringstream hN6770;
	  int kHist = kOffHi + 6770 + kh + 1;
	  hN6770 << kHist;
	  vRC6770.push_back( new CsHist2D( hN6770.str(), hTitle,
	                                   100, -1.,  3., 100, -1., 3. ) );
	}
	for( int kh=0; kh<4; kh++ ) {
	  string hTitle = " ";
	  stringstream hN6770;
	  int kHist = kOffHi + 6770 + kh + 6;
	  hN6770 << kHist;
	  vRC6770.push_back( new CsHist2D( hN6770.str(), hTitle,
	                                   100, -500., 500., 50, 0., 50. ) );
	}
	vRC6780.clear();
	for( int kh=0; kh<9; kh++ ) {
	  string hTitle = " ";
	  stringstream hN6780;
	  int kHist = kOffHi + 6780 + kh + 1;
	  hN6780 << kHist;
	  vRC6780.push_back( new CsHist2D( hN6780.str(), hTitle,
	                                   100, -1.,  3., 100, -1., 3. ) );
	}
	vRC6790.clear();
	string hTitle = " ";
        int kHist = 0;
	kHist = kOffHi + 6790 + 1;
        stringstream hN6791;
        hN6791 << kHist;
        vRC6790.push_back( new CsHist2D( hN6791.str(), hTitle,
                                         101, 1.0014, 1.0016, 20, 0., 20. ) );
	kHist = kOffHi + 6790 + 2;
        stringstream hN6792;
	hN6792 << kHist;
        vRC6790.push_back( new CsHist2D( hN6792.str(), hTitle,
                                         101, 10., 60., 20, 0., 20. ) );
	for( int kh=0; kh<3; kh++ ) {
	  kHist = kOffHi + 6790 + kh + 3;
          stringstream hN6793;
	  hN6793 << kHist;
          vRC6790.push_back( new CsHist2D( hN6793.str(), hTitle,
                                           100, -200., 200., 50, 0., 50. ) );
	}
	CsHistograms::SetCurrentPath("/");
      }
    }

    int nProb = cons->outBufferSize();
    double likeDV = cons->likeDefVa();
//- order of mag. of n :     1.001500 (UV)   1.001350 (VS)
//- UV refr. index
    float CFRefIndP = cons->CFRefInd();
    //static double dRefIndT = 0.000000;
    //static double dRefIndT = 0.000015;      // to get 'true' like 1.1% (VS)
    static double dRefIndT = 0.000035;      // to get 'true' like 2.6% (VS)
    //static double dRefIndT = -0.000035;      // to get 'true' like 2.6% (VS)
    //static double dRefIndT = 0.000050;      // to get 'true' like 3.7% (VS)
    //static double dRefIndT = 0.000065;      // to get 'true' like 4.8% (VS)
    //static double dRefIndT = 0.000100;      // to get 'true' like 7.4% (VS)
//@@---------------------------------
    double CFRefIndT = CFRefIndP + dRefIndT;
//- VS refr. index
    float CFRefIndPVS = cons->CFRefIndVS();
    //static double dRefIndTVS = 0.000015;
    static double dRefIndTVS = 0.000025;
    //static double dRefIndTVS = 0.000035;
    //static double dRefIndTVS = 0.000050;
    //static double dRefIndTVS = 0.000065;
//@@-----------------------------------
    double CFRefIndTVS = CFRefIndPVS + dRefIndTVS;

    static double dRefIndD = 0.000050;      // for deriv. comput. 3.7% (VS)
//@@---------------------------------

//- loop over part-photons :
//  ------------------------
    bool isFirst = true;
    list<CsRCPartPhotons*> lPartPhotons = lPartPhotons_;
    list<CsRCPartPhotons*>::iterator ipf;
    for( ipf=lPartPhotons.begin(); ipf!=lPartPhotons.end(); ipf++ ) {
      if( !(*ipf)->flag() ) continue;
      CsRCPartPhotons* papho = (*ipf);

      float momPart = papho->pPart()->mom();

      double* probaLK;
      probaLK = papho->probaLKAll();
//--- Production Likelihoods
      double likeBkgP = papho->probaLKBgAll();
      double likeElecP = probaLK[ 2];
      double likeMuonP = probaLK[ 5];
      double likePionP = probaLK[ 8];
      double likeKaonP = probaLK[11];
      double likeProtonP = probaLK[14];
      if( likeElecP == likeBkgP ) continue;      // NO photons from signal
//--- Select ID pions (analysis way) Production Likelihoods 
      bool bIDpionP = false;
      if( likePionP > likeDV  &&  likePionP > likeBkgP  &&
          likePionP > likeKaonP  &&  likePionP > likeProtonP )
	bIDpionP = true;
      if( momPart < 8.  &&  likePionP < likeElecP ) bIDpionP = false;
//--- Select ID kaons (analysis way) Production Likelihoods
      bool bIDkaonP = false;
      if( likeKaonP > likeDV  &&  likeKaonP > likeBkgP  &&
          likeKaonP > likePionP  &&  likeKaonP > likeProtonP )
	bIDkaonP = true;
      if( momPart < 8.  &&  likeKaonP < likeElecP ) bIDkaonP = false;

      double probaLKT[31];
      double probaLKC[31];
      //for( int kPaTy=1; kPaTy<=31; kPaTy++ ) {
//--- corrected   100201   cppcheck
      for( int kPaTy=0; kPaTy<31; kPaTy++ ) {
	probaLKT[kPaTy] = likeDV;
	probaLKC[kPaTy] = likeDV;
      }

//--- ASSUME production Index CFRefIndP 'WRONG' and
//                            CFRefIndT the GOOD one ('TRUE')

//--- For all(5) particles :
//    ----------------------
//    Compute the 'TRUE' Likelihoods
      for( int kPaTy=2; kPaTy<=14; kPaTy+=3 ) {
	double pLikeProd = probaLK[kPaTy];
	if( pLikeProd == likeDV ) continue;      //   under thresh.
	double massPart = papho->pPart()->mass( kPaTy );
	double theIpoT = papho->getThetaIpo( CFRefIndT, momPart, massPart );
//                              -------------------------------------------
//----- Set modified = 'TRUE' refr. indexes
	cons->setCFRefInd( CFRefIndT );
	cons->setCFRefIndVS( CFRefIndTVS );
	int nPhotT = 0;
	probaLKT[kPaTy] = papho->getLikelihood( theIpoT, nPhotT );
//                               --------------------------------
//----- Restore Production refr. indexes
	cons->setCFRefInd( CFRefIndP );
	cons->setCFRefIndVS( CFRefIndPVS );
      }
      int nPhotB = 0;
      double theIpoB = 0.;
      double probaBkT = papho->getLikelihood( theIpoB, nPhotB );
//                             --------------------------------
//--- 'TRUE' Likelihoods
      double likeBkgT = probaBkT;
      double likeElecT = probaLKT[ 2];
      double likeMuonT = probaLKT[ 5];
      double likePionT = probaLKT[ 8];
      double likeKaonT = probaLKT[11];
      double likeProtonT = probaLKT[14];

//--- Production Derivatives : UpRICHJob
      double* derLKAPV;
      derLKAPV = papho->derivLKAllAPV();
      double* derLKPMT;
      derLKPMT = papho->derivLKAllPMT();
      int nPhotAPV = int( papho->partProbs()[28] );
      int nPhotPMT = int( papho->partProbs()[29] );
//--- Production Derivatives : NOT UpRICHJob
      double* derLK;
      derLK = papho->derivLKAll();
      int nPhotTOT = int( papho->partProbs()[9] );

      for( int kPaTy=2; kPaTy<=14; kPaTy+=3 ) {

	double pLikeProd = probaLK[kPaTy];
	if( pLikeProd == likeDV ) continue;      //   under thresh.
	double pLikeTrue = likeDV;
	double pLikeCorr = likeDV;
	double lCorr = 0.;
	pLikeTrue = probaLKT[kPaTy];
        if( pLikeTrue == pLikeProd ) continue;   // thresh. probl. OBSOLETE

	if( CsRichOne::Instance()->UpRICHJob()  &&  nProb > 21 ) {

//------- Compute Like Corr from derivatives
          //lCorr = derLKAPV[kPaTy]*dRefIndT + derLKPMT[kPaTy]*dRefIndT;
          lCorr = derLKAPV[kPaTy]*dRefIndT + derLKPMT[kPaTy]*dRefIndTVS;
          pLikeCorr = pLikeProd + lCorr;
          probaLKC[kPaTy] = pLikeCorr;

//------- Check APPROXIMATE intermediate proc.
	  double pLikeCorrA = 0.;
          bool checkApp = true;
          if( checkApp ) {
	    bool twoPtDer = true;
	    double massPart = papho->pPart()->mass( kPaTy );
            double DCFRefInd = 0.000050;
	    if( twoPtDer ) DCFRefInd /= 2.;
            double CFRefIndDP = cons->CFRefInd() + DCFRefInd;
	    double CFRefIndDM = cons->CFRefInd() - DCFRefInd;
            double theIpoDP = 0.;
	    double theIpoDM = 0.;
	    theIpoDP = papho->getThetaIpo( CFRefIndDP, momPart, massPart );
//                            --------------------------------------------
	    theIpoDM = papho->getThetaIpo( CFRefIndDM, momPart, massPart );
//                            --------------------------------------------
	    double likeDer = 0.;
	    int nPhotDP = 0;
            long double pLikeDP = likeDV;
            long double pLikeDM = likeDV;
            if( theIpoDP > 0. ) {
              pLikeDP = papho->getLikelihood( theIpoDP, nPhotDP );
//                             ----------------------------------
	      likeDer = (pLikeDP - pLikeProd) / DCFRefInd;
	    }
	    if( twoPtDer  &&   theIpoDP > 0.  &&  theIpoDM > 0. ) {
	      pLikeDM = papho->getLikelihood( theIpoDM, nPhotDP );
//                             ----------------------------------
	      likeDer = (pLikeDP - pLikeDM) / (2.*DCFRefInd);
	    }
	    //likeDer = derLKAPV[kPaTy] + derLKPMT[kPaTy];
	    //lCorr = likeDer*dRefIndT;
	    //lCorr = likeDer*dRefIndTVS;
	    double fUV = float( nPhotAPV )/float( nPhotTOT );
	    double fVS = float( nPhotPMT )/float( nPhotTOT );
	    //lCorr = likeDer * (fUV*dRefIndT + fVS*dRefIndTVS);
	    pLikeCorrA = pLikeProd + lCorr;
	    if( nPhotDP > 0 ) {
	      //xh = pLikeCorr;
	      xh = pLikeTrue;
	      yh = pLikeCorrA;
	      int kHH = (kPaTy-2)/3;
              if( kHH >= 0  &&  kHH < int( vRC6770.size() ) ) {
       	        if( vRC6770[kHH] ) vRC6770[kHH]->Fill( xh, yh );
//hh                               ----------------------------
	      }
	      if( nPhotAPV == 0 ) {
	        int kHH = (kPaTy-8)/3;
                if( kHH >= 0  &&  kHH < int( vRC6780.size() ) ) {
       	          if( vRC6780[kHH] ) vRC6780[kHH]->Fill( xh, yh );
//hh                                 ----------------------------
	        }
	      }
	      else if( nPhotPMT == 0 ) {
	        int kHH = (kPaTy-8)/3;
		if( kHH >=0 ) {
		  kHH += 6;
                  if( kHH >= 0  &&  kHH < int( vRC6780.size() ) ) {
       	            if( vRC6780[kHH] ) vRC6780[kHH]->Fill( xh, yh );
//hh                                   ----------------------------
		  }
	        }
	      } else {
	        int kHH = (kPaTy-8)/3;
		if( kHH >=0 ) {
		  kHH += 3;
                  if( kHH >= 0  &&  kHH < int( vRC6780.size() ) ) {
       	            if( vRC6780[kHH] ) vRC6780[kHH]->Fill( xh, yh );
//hh                                   ----------------------------
		  }
	        }
	      }
	      xh = likeDer - (derLKAPV[kPaTy]+derLKPMT[kPaTy]);
	      yh = momPart;
	      kHH = (kPaTy-8)/3;
	      if( kHH >=0 ) {
		kHH += 2;
		if( kHH >= 0  &&  kHH < int( vRC6790.size() ) ) {
		  if( vRC6790[kHH] ) vRC6790[kHH]->Fill( xh, yh );
//hh                                 ----------------------------
		}
	      }
	    }
          }      // if   checkApp
//------- histo's :
          //xh = (pLikeCorr - pLikeProd) / pLikeProd;
	  //xh = (pLikeCorr - pLikeTrue) / (pLikeProd - pLikeTrue);  // cccc
	  //xh = (pLikeCorr - pLikeTrue) / pLikeProd;                // d
	  xh = (pLikeCorrA - pLikeTrue) / pLikeProd;
          yh = momPart;
          if( nPhotAPV > 0  &&  nPhotPMT > 0 ) {
            int kHH = (kPaTy-2)/3;
            if( kHH >= 0  &&  kHH < int( vRC6750.size() ) ) {
       	      if( vRC6750[kHH] ) vRC6750[kHH]->Fill( xh, yh );
//hh                             ----------------------------
	    }
	  } else {
	    int kHH = (kPaTy-2)/3 + 5;
            if( kHH >= 0  &&  kHH < int( vRC6750.size() ) ) {
      	      if( vRC6750[kHH] ) vRC6750[kHH]->Fill( xh, yh );
//hh                             ----------------------------
	    }
	  }
	  yh = float( nPhotAPV )/float( nPhotPMT );
          int kHH = (kPaTy-2)/3;
          if( kHH >= 0  &&  kHH < int( vRC6760.size() ) ) {
            if( vRC6760[kHH] ) vRC6760[kHH]->Fill( xh, yh );
//hh                           ----------------------------
          }
	  print = false;
	  if( print ) {
	    if( nPhotAPV > 0  &&  nPhotPMT > 0 ) {
      	      std::cout << "  " <<  setprecision(8);
              std::cout << kPaTy << " - " << nPhotAPV << "  " << nPhotPMT
                        << " - " << pLikeProd << "  " << pLikeTrue << " - "
                        << pLikeTrue << "  "
      	                //<< pLikeAPV << "  " << pLikePMT << " - "
                        << lCorr << "  " << pLikeCorr << "  " << xh
                        << std::endl;
	    }
	  }

	} else {

//------- Compute Like Corr from derivatives
          lCorr = derLK[kPaTy]*dRefIndT;
          pLikeCorr = pLikeProd + lCorr;
          probaLKC[kPaTy] = pLikeCorr;

//------- histo's :
          xh = (pLikeCorr - pLikeProd) / pLikeProd;
	  //xh = (pLikeCorr - pLikeTrue) / (pLikeProd - pLikeTrue);
          yh = momPart;
          if( nPhotTOT > 0 ) {
            int kHH = (kPaTy-2)/3;
            if( kHH >= 0  &&  kHH < int( vRC6750.size() ) ) {
       	      if( vRC6750[kHH] ) vRC6750[kHH]->Fill( xh, yh );
//hh                             ----------------------------
	    }
	  }

        }      //   if/else  UpRICHJob()

      }      //   for  kPaTy
      print = false;
      if( print ) {
	if( nPhotAPV > 0  &&  nPhotPMT > 0 ) {
          std::cout << std::endl;
          std::cout << setprecision(8);
          std::cout << dRefIndT << "  " << dRefIndTVS << std::endl;
          for( int k=2; k<=14; k+=3 ) std::cout << probaLK[k] << "  ";
          std::cout << std::endl;
          for( int k=2; k<=14; k+=3 ) std::cout << probaLKT[k] << "  ";
          std::cout << std::endl;
          for( int k=2; k<=14; k+=3 ) std::cout << probaLKC[k] << "  ";
          std::cout << std::endl;
	}
      }

      double derLikePionP = 0.;
      double derLikeKaonP = 0.;
      int nPhotH = 0;
      if( CsRichOne::Instance()->UpRICHJob()  &&  nProb > 21 ) {
	derLikePionP = papho->derivLKAllPMT()[8];
	derLikeKaonP = papho->derivLKAllPMT()[11];
	nPhotH = nPhotPMT;
      } else {
	derLikePionP = papho->derivLKAll()[8];
	derLikeKaonP = papho->derivLKAll()[11];
	nPhotH = nPhotTOT;
      }
//--- Histo Pion (ID or notID) Production (PMT) derivatives
      if( bIDpionP  &&  probaLK[8] > likeDV  &&  nPhotH > 0 ) {
        xh = derLikePionP;
        yh = momPart;
        int kHH = 5;
        if( kHH >= 0  &&  kHH < int( vRC6770.size() ) ) {
          if( vRC6770[kHH] ) vRC6770[kHH]->Fill( xh, yh );
//hh                         ----------------------------
        }
      }
      if( !bIDpionP  &&  probaLK[8] > likeDV  &&  nPhotH > 0 ) {
        xh = derLikePionP;
        yh = momPart;
        int kHH = 6;
        if( kHH >= 0  &&  kHH < int( vRC6770.size() ) ) {
          if( vRC6770[kHH] ) vRC6770[kHH]->Fill( xh, yh );
//hh                         ----------------------------
        }
      }
//--- Histo Kaon (ID or notID) Production derivatives
      if( bIDkaonP  &&  probaLK[11] > likeDV  &&  nPhotH > 0 ) {
        xh = derLikeKaonP;
        yh = momPart;
        int kHH = 7;
        if( kHH >= 0  &&  kHH < int( vRC6770.size() ) ) {
          if( vRC6770[kHH] ) vRC6770[kHH]->Fill( xh, yh );
//hh                         ----------------------------
        }
      }
      if( !bIDkaonP  &&  probaLK[11] > likeDV  &&  nPhotH > 0 ) {
        xh = derLikeKaonP;
        yh = momPart;
        int kHH = 8;
        if( kHH >= 0  &&  kHH < int( vRC6770.size() ) ) {
          if( vRC6770[kHH] ) vRC6770[kHH]->Fill( xh, yh );
//hh                         ----------------------------
        }
      }


/*
//--- Check Pion and Kaon ID inversion
      double pLikeCorr;
//--- Select ID pions (analysis way) 'TRUE' likelihoods
      bool bIDpionT = false;
      if( likePionT > likeDV  &&  likePionT > likeBkgT  &&
          likePionT > likeKaonT  &&  likePionT > likeProtonT )
	bIDpionT = true;
      if( momPart < 8.  &&  likePionT < likeElecT ) bIDpionT = false;
//--- Histo Pion ID inversion
      pLikeCorr = probaLKC[8];
      //pLikeCorr = probaLK[8];      //
      if( bIDpionT  &&  probaLKT[5] > likeDV  &&  pLikeCorr > likeDV ) {
        xh = pLikeCorr - probaLKT[5];      // Like ID-pion corr - muon true
        yh = momPart;
        int kHH = 5;
        if( kHH >= 0  &&  kHH < int( vRC6760.size() ) ) {
      	  if( vRC6760[kHH] ) vRC6760[kHH]->Fill( xh, yh );
//hh                         ----------------------------
        }
      }
      if( bIDpionT  &&  probaLKT[11] > likeDV  &&  pLikeCorr > likeDV ) {
        xh = pLikeCorr - probaLKT[11];      // Like ID-pion corr - kaon true
        yh = momPart;
        int kHH = 6;
        if( kHH >= 0  &&  kHH < int( vRC6760.size() ) ) {
          if( vRC6760[kHH] ) vRC6760[kHH]->Fill( xh, yh );
//hh                         ----------------------------
        }
      }
//--- Select ID kaons (analysis way) 'TRUE' likelihoods
      bool bIDkaonT = false;
      if( likeKaonT > likeDV  &&  likeKaonT > likeBkgT  &&
          likeKaonT > likePionT  &&  likeKaonT > likeProtonT )
	bIDkaonT = true;
      if( momPart < 8.  &&  likeKaonT < likeElecT ) bIDkaonT = false;
//--- Histo Kaon ID inversion
      pLikeCorr = probaLKC[11];
      //pLikeCorr = probaLK[11];      //
      if( bIDkaonT  &&  probaLKT[8] > likeDV  &&  pLikeCorr > likeDV ) {
        xh = pLikeCorr - probaLKT[8];      // Like ID-kaon corr - pion true
        yh = momPart;
        int kHH = 7;
        if( kHH >= 0  &&  kHH < int( vRC6760.size() ) ) {
          if( vRC6760[kHH] ) vRC6760[kHH]->Fill( xh, yh );
//hh                         ----------------------------
        }
      }
      if( bIDkaonT  &&  probaLKT[14] > likeDV  &&  pLikeCorr > likeDV ) {
        xh = pLikeCorr - probaLKT[14];      // Like ID-kaon corr - prot true
        yh = momPart;
        int kHH = 8;
        if( kHH >= 0  &&  kHH < int( vRC6760.size() ) ) {
      	  if( vRC6760[kHH] ) vRC6760[kHH]->Fill( xh, yh );
//hh                         ----------------------------
        }
      }
*/

//--- 081105
//--- ASSUME production Index CFRefIndP 'GOOD' and
//                            CFRefIndT the wrong one ('WRONG')
//--- USE ID Pions AND KAONS (analysis way) Production likelihoods
//--- Check Pion and Kaon ID inversion
      double derLKAPVW = 0.;
      double derLKPMTW = 0.;
      double pLikeAPVW = likeDV;
      double pLikePMTW = likeDV;
      double pLikeCorrW[31];
      for( int k=0; k<31; k++ ) pLikeCorrW[k] = likeDV;
      double probaLKW[31];
      for( int k=0; k<31; k++ ) probaLKW[k] = probaLKT[k];
      int nPhotAPVW = 0;
      int nPhotPMTW = 0;
      double lCorr = 0.;
      for( int kPaTy=5; kPaTy<=14; kPaTy+=3 ) {
        double massPart = papho->pPart()->mass( kPaTy );
        double theIpoW = papho->getThetaIpo( CFRefIndT, momPart, massPart );
//                              -------------------------------------------
        double pLikeProd = probaLK[kPaTy];
        double pLikeWrong = probaLKW[kPaTy];
        papho->getLikeDeriv( kPaTy, theIpoW, pLikeWrong, derLKAPVW, derLKPMTW,
//             ---------------------------------------------------------------
		  	     pLikeAPVW, pLikePMTW, nPhotAPVW, nPhotPMTW );
        lCorr = derLKAPVW*dRefIndT + derLKPMTW*dRefIndTVS;      //
        double likeDer = derLKAPVW + derLKPMTW;
        //lCorr = likeDer*dRefIndT;
        //lCorr = likeDer*dRefIndTVS;                             //
        int nPhotTOTW = nPhotAPVW + nPhotPMTW;
        double fUV = float( nPhotAPVW )/float( nPhotTOTW );
        double fVS = float( nPhotPMTW )/float( nPhotTOTW );
        //lCorr = likeDer * (fUV*dRefIndT + fVS*dRefIndTVS);      //
        pLikeCorrW[kPaTy] = pLikeWrong - lCorr;
      }
      double pLikeCorr = likeDV;
//--- Histo Pion ID inversion
      pLikeCorr = pLikeCorrW[8];
      //pLikeCorr = probaLKW[8];      //
      int kPaTy = 5;
      if( bIDpionP  &&  probaLK[kPaTy] > likeDV  &&  pLikeCorr > likeDV ) {
        xh = pLikeCorr - probaLK[kPaTy];      // Like ID-pion corr - muon true
        yh = momPart;
        int kHH = 5;
        if( kHH >= 0  &&  kHH < int( vRC6760.size() ) ) {
      	  if( vRC6760[kHH] ) vRC6760[kHH]->Fill( xh, yh );
//hh                         ----------------------------
        }
      }
      kPaTy = 11;
      if( bIDpionP  &&  probaLK[kPaTy] > likeDV  &&  pLikeCorr > likeDV ) {
        xh = pLikeCorr - probaLK[kPaTy];     // Like ID-pion corr - kaon true
        yh = momPart;
        int kHH = 6;
        if( kHH >= 0  &&  kHH < int( vRC6760.size() ) ) {
          if( vRC6760[kHH] ) vRC6760[kHH]->Fill( xh, yh );
//hh                         ----------------------------
        }
      }
//--- Histo Kaon ID inversion
      pLikeCorr = pLikeCorrW[11];
      //pLikeCorr = probaLKW[11];      //
      kPaTy = 8;
      if( bIDkaonP  &&  probaLK[kPaTy] > likeDV  &&  pLikeCorr > likeDV ) {
        xh = pLikeCorr - probaLK[kPaTy];     // Like ID-kaon corr - pion true
        yh = momPart;
        int kHH = 7;
        if( kHH >= 0  &&  kHH < int( vRC6760.size() ) ) {
          if( vRC6760[kHH] ) vRC6760[kHH]->Fill( xh, yh );
//hh                         ----------------------------
        }
      }
      kPaTy = 14;
      if( bIDkaonP  &&  probaLK[kPaTy] > likeDV  &&  pLikeCorr > likeDV ) {
        xh = pLikeCorr - probaLK[kPaTy];     // Like ID-kaon corr - prot true
        yh = momPart;
        int kHH = 8;
        if( kHH >= 0  &&  kHH < int( vRC6760.size() ) ) {
      	  if( vRC6760[kHH] ) vRC6760[kHH]->Fill( xh, yh );
//hh                         ----------------------------
        }
      }


      bool checkProf = true;
      if( checkProf ) {
	static int nHMax = 20;
	static int nPro = 0;
//----- Check Like profile vs refr. index
	double massPart = 0.;
        if( bIDpionP ) massPart = papho->pPart()->mass( 8 );
	if( bIDkaonP ) massPart = papho->pPart()->mass( 11 );
//----- ID pions (analysis way)
        if( bIDpionP ) {
        //if( !bIDpionP ) {
//----- ID kaons (analysis way)
	//if( bIDkaonP ) {
	  if( papho->likeProSet() ) {
  	    double theMin = papho->theLikeProMn();
	    double theMax = papho->theLikeProMx();
	    int nIpo = papho->nLikePro();
            double dThe = (theMax - theMin)/float( nIpo-1 );
            double thek[nIpo];
            double pLikek[nIpo];
            for( int kIpo=0; kIpo<nIpo; kIpo++ ) {
              thek[kIpo] = papho->theLikePro()[kIpo];
              pLikek[kIpo] = papho->pLikePro()[kIpo];
            }
	    if( nPro < nHMax ) {
	      double CFRefIndW1 = CFRefIndP - 3.*dRefIndT;
	      double CFRefIndW2 = CFRefIndP + 3.*dRefIndT;
	      //std::cout << setprecision(8);
	      //std::cout << CFRefIndP << "  " << CFRefIndW1 << "  "
	      //	  << CFRefIndW2 << std::endl;
	      int nPo = 50;
	      double CFRefIndW = 0.;
	      double dCFRefIndW = (CFRefIndW2 - CFRefIndW1)/float( nPo );
	      xh = CFRefIndW1 + dCFRefIndW/2.;
	      yh = double( nPro ) + 0.5;
	      wh = momPart;
	      int kHH = 0;
	      if( kHH >= 0  &&  kHH < int( vRC6790.size() ) ) {
		if( vRC6790[kHH] ) vRC6790[kHH]->Fill( xh, yh, wh );
//hh                               --------------------------------
	      }
	      for( int kPo=1; kPo<nPo; kPo++ ) {
		CFRefIndW = CFRefIndW1 + kPo*dCFRefIndW;
		double pLike = 0.;
		double theta =
		  papho->getThetaIpo( CFRefIndW, momPart, massPart );
//                       -------------------------------------------
		for( int kIpo=0; kIpo<nIpo-1; kIpo++ ) {
		  if( theta >= thek[kIpo]  &&  theta < thek[kIpo+1] ) {
		    pLike = pLikek[kIpo] + (theta - thek[kIpo])/dThe *
		      (pLikek[kIpo+1] - pLikek[kIpo]);
		    break;
		  }
		}
		//xh = theta;
		xh = CFRefIndW;
		yh = double( nPro ) + 0.5;
		wh = pLike;
		int kHH = 0;
		if( kHH >= 0  &&  kHH < int( vRC6790.size() ) ) {
		  if( vRC6790[kHH] ) vRC6790[kHH]->Fill( xh, yh, wh );
//hh                                 --------------------------------
		}
	      }

	      int kIpo = 0;
	      xh = thek[kIpo];
	      yh = double( nPro ) + 0.5;
	      wh = momPart;
	      kHH = 1;
	      if( kHH >= 0  &&  kHH < int( vRC6790.size() ) ) {
		if( vRC6790[kHH] ) vRC6790[kHH]->Fill( xh, yh, wh );
//hh                             ------------------------------
	      }
	      CFRefIndW1 = CFRefIndP + dRefIndD/2.;
	      CFRefIndW2 = CFRefIndP - dRefIndD/2.;
	      double CFRefIndW3 = CFRefIndP + dRefIndD;
	      double theta0 =
	        papho->getThetaIpo( CFRefIndP, momPart, massPart );
//                     -------------------------------------------
	      kIpo = 1;
	      xh = thek[kIpo];
	      yh = double( nPro ) + 0.5;
	      wh = theta0;
	      kHH = 1;
	      if( kHH >= 0  &&  kHH < int( vRC6790.size() ) ) {
		if( vRC6790[kHH] ) vRC6790[kHH]->Fill( xh, yh, wh );
//hh                             ------------------------------
	      }
	      double theta1 =
	        papho->getThetaIpo( CFRefIndW1, momPart, massPart );
//                     --------------------------------------------
	      kIpo = 2;
	      xh = thek[kIpo];
	      yh = double( nPro ) + 0.5;
	      wh = theta1;
	      kHH = 1;
	      if( kHH >= 0  &&  kHH < int( vRC6790.size() ) ) {
		if( vRC6790[kHH] ) vRC6790[kHH]->Fill( xh, yh, wh );
//hh                             ------------------------------
	      }
	      double theta2 =
	        papho->getThetaIpo( CFRefIndW2, momPart, massPart );
//                     --------------------------------------------
	      kIpo = 3;
	      xh = thek[kIpo];
	      yh = double( nPro ) + 0.5;
	      wh = theta2;
	      kHH = 1;
	      if( kHH >= 0  &&  kHH < int( vRC6790.size() ) ) {
		if( vRC6790[kHH] ) vRC6790[kHH]->Fill( xh, yh, wh );
//hh                             ------------------------------
	      }
	      double theta3 =
	        papho->getThetaIpo( CFRefIndW3, momPart, massPart );
//                     --------------------------------------------
	      kIpo = 4;
	      xh = thek[kIpo];
	      yh = double( nPro ) + 0.5;
	      wh = theta3;
	      kHH = 1;
	      if( kHH >= 0  &&  kHH < int( vRC6790.size() ) ) {
		if( vRC6790[kHH] ) vRC6790[kHH]->Fill( xh, yh, wh );
//hh                             ------------------------------
	      }
	      //std::cout << theta0 << "  " << theta1 << "  " << theta2
	      //	  << std::endl;
	      for( int kIpo=5; kIpo<nIpo; kIpo++ ) {
		xh = thek[kIpo];
		yh = double( nPro ) + 0.5;
		wh = pLikek[kIpo];
		int kHH = 1;
		if( kHH >= 0  &&  kHH < int( vRC6790.size() ) ) {
		  if( vRC6790[kHH] ) vRC6790[kHH]->Fill( xh, yh, wh );
//hh                               ------------------------------
		}
	      }
	      nPro++;
	    }
          }
        }
      }      // if  checkProf

    }      //   for  lPartPhotons

    return;
  }


