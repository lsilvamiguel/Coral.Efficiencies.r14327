/*!
   \file    CsRCEventRings.cc
   \-------------------------
   \brief   CsRCEventRings class implementation.
   \author  Paolo Schiavon
   \version 0.02,  rev. 20/6/00,  rev. October 2000
   \date    25 June 1999
*/

//- revised  August 2005


  #include <iostream>
  #include <ostream>
  #include <cstdio>
  #include <cmath>
  #include <cstdlib>

  #include <CLHEP/Vector/ThreeVector.h>
  #include "CLHEP/Matrix/Matrix.h"

  #include "CsErrLog.h"
  #include "CsTrack.h"
//-------------------------------------


// ----------------------------
  #include "CsRCEventRings.h"

  #include "CsRCMirrors.h"
  #include "CsRCDetectors.h"

  #include "CsRCEventParticles.h"
  #include "CsRCParticle.h"
  #include "CsRCEventPads.h"
  #include "CsRCPad.h"
  #include "CsRCCluster.h"

  #include "CsRCEventPartPhotons.h"
  #include "CsRCPartPhotons.h"
  #include "CsRCPhoton.h"

  #include "CsRCLikeAll.h"
  #include "CsRCLikeAll02.h"
  #include "CsRCLikeAll03.h"

  #include "CsRCRing.h"
  #include "CsRCLikeRing.h"
  #include "CsRCLikeRing02.h"
  #include "CsRCLikeRing03.h"
  #include "CsRCLikeRing04.h"

  #include "CsRCEventAnalysis.h"

  #include "CsRCChiSqFit.h"
  #include "CsRCCircleFit.h"

  #include "CsRCExeKeys.h"
  #include "CsRCRecConst.h"
  #include "CsRCHistos.h"
// ----------------------------

  using namespace std;
  using namespace CLHEP;

  CsRCEventRings* CsRCEventRings::instance_ = 0;

//===========================================================================
  CsRCEventRings* CsRCEventRings::Instance() {
//--------------------------------------------
    if( instance_ == 0 ) instance_ = new CsRCEventRings();
    return instance_; 
  }

//===========================================================================
  CsRCEventRings::CsRCEventRings() {
//----------------------------------
    lRings_.clear();
    flagRingSel_ = false;
  }

//===========================================================================
  CsRCEventRings::CsRCEventRings( const CsRCEventRings &evring ) {
//----------------------------------------------------------------
    cout << "RICHONE : CsRCEventRings CopyConstructor" << endl;
    instance_ = evring.instance_;
    lRings_ = evring.lRings_;
    flag_ = evring.flag_;
    flagRingSel_ = evring.flagRingSel_;
  }

//===========================================================================
  void CsRCEventRings::clearEventRings() {
//----------------------------------------
    list<CsRCRing*>::iterator ir;
    for( ir=lRings_.begin(); ir!=lRings_.end(); ir++ ) delete *ir;
    lRings_.clear();
  }

//===========================================================================
  void CsRCEventRings::print() const {
//------------------------------------
    list<CsRCRing*>::const_iterator ir;
    for( ir=lRings_.begin(); ir!=lRings_.end(); ir++ ) (*ir)->print();
  }

//===========================================================================
  CsRCEventRings::~CsRCEventRings() {
//-----------------------------------
    clearEventRings();
  }


//===========================================================================
  void CsRCEventRings::getEventRings() {
//--------------------------------------

//--- Procedure for ring recognition.
//    -------------------------------
//--- Paolo  -  May 1999
//--- from :  ringrec-3.f,  RingRec-2.cc,  RingRec3.cc
//            ----------------------------------------
//--- from :  ringrec-4.f
//            -----------
//--- rev.  October  2000
//--- rev.  December  2002
//--- rev.  August  2005


//--- from "CsRCExecKeys"
      CsRCExeKeys *key = CsRCExeKeys::Instance();
      bool DoClu = key->DoClu();
      bool DoBkgFil = key->DoBkgFil();
      bool DoWgtAv = key->DoWgtAv();
      bool SigmaCut = key->SigmaCut();
      bool AliMirror = key->AliMirror();
      int kPrint = key->kPrintEventRings();

//--- from "CsRCRecConst"
      CsRCRecConst *cons = CsRCRecConst::Instance();
      int kOffHi = cons->kOffHi();

      CsRCHistos& hist = CsRCHistos::Ref();
      int kh;
      double xh, wh;

      xh = 50.5;
      if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                       ------------------------


//--------- RING RECOGNITION ----------
//-------------------------------------

//--- LOOP over PART-PHOTONS :
//    ------------------------
      CsRCEventPartPhotons* paphos = CsRCEventPartPhotons::Instance();
      list<CsRCPartPhotons*> lPaPhotons = paphos->lPartPhotons();
      list<CsRCPartPhotons*>::iterator ipf;

      for( ipf=lPaPhotons.begin(); ipf!=lPaPhotons.end(); ipf++ ) {
	if( !(*ipf)->flag() ) continue;

	xh = 51.5;
	if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                         ------------------------

//----- change photon the() angle definition :
//      --------------------------------------
	(*ipf)->setPhoTheNorm();
//@@---------------------------

//----- first order ring definition :
//      -----------------------------
	CsRCRing* ring = new CsRCRing();
        int kRing = lRings_.size();

	if( cons->ringDefMode() == "MAXLIKE" ) {
          ring->getRingMaxLike( kRing, (*ipf) );
//        -------------------------------------
	} else {
          ring->getRingPk( kRing, (*ipf) );
//        --------------------------------
	}
	if( ring->flag() ) {

	  lRings_.push_back( ring );

	  xh = 52.5;
	  if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                           ------------------------

//------- monitoring photons in ring moved to getRingPk
	  if( kPrint >= 3 ) {
  	    cout << endl << " Ring Pk :";
	    ring->print();
	    ring->printPhotons();
	    ring->printPhoFlag();
	  }

	} else { delete ring; }

      }   /* end loop on Part-photons: ipf */

      if( lRings_.empty() ) {
	xh = 53.5;
	if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                         ------------------------
      }


//--- flag photons (clusters) of ring 'signal', for all part-photons :
//    ----------------------------------------------------------------
      //^if( DoBkgFil )  flagPhoSignal();
//                    ---------------

//--- LOOP over RINGS found (Pk):
//    ---------------------------
      list<CsRCRing*>::iterator ir;
      for( ir=lRings_.begin(); ir!=lRings_.end(); ir++ ) {
        if( !(*ir)->flag() ) continue;

//----- ring MEANTIME (PMT only) :
//      --------------------------
        (*ir)->setTime();
//      ---------------- 

//----- BACKGROUND filter :
//      -------------------
        if( !(*ir)->flag() ) continue;
        if( DoBkgFil )  {
	  if( kPrint >= 3 ) {
	    cout << endl << " Ring Pk,  ";
	    (*ir)->printPhoFlag();
	  }

          //^(*ir)->getRingBf();
//        ------------------

//------- monitoring histograms moved to getRingBf
          if( kPrint >= 3 ) {
            cout << endl << " Ring Bf :";
            (*ir)->print();
          }
        }   /* end if on DoBkgFil */


//----- Weighted PHI average :
//      ----------------------
        if( !(*ir)->flag() ) continue;
        if( DoWgtAv )  {

          //^(*ir)->getRingWa();
//        ------------------
          if( kPrint >= 3 ) {
            cout << endl << " Ring Wa :";
            (*ir)->print();
          }
        }   /* end if on DoWgtAv */


//----- three sigma CUT around peak :
//      -----------------------------
        if( !(*ir)->flag() ) continue;
        if( SigmaCut ) {

          //^(*ir)->getRing3S();
//        ------------------
//------- monitoring photons in ring moved to getRing3S
          if( kPrint >= 3 ) {
            cout << endl << " Ring 3S :";
            (*ir)->print();
          }
        }   /* end if on SigmaCut */

      }   /* end loop on Rings : ir */


//--- monitoring histogram :
//    ----------------------
      int nRingEv = 0;
      for( ir=lRings_.begin(); ir!=lRings_.end(); ir++ ) {
        if( (*ir)->flag() ) nRingEv++;
      }
      xh = nRingEv;
      if( hist.hRC1535 ) hist.hRC1535->Fill( xh );
//hh                     ------------------------

//--- select 'useful' rings :
//    -----------------------
      flagRingSel_ = false;
      bool RingSel = key->RingSelect();
      if( RingSel )  { flagRingSel_ = true;  ringSelection(); }
//                                           ---------------

//--- sigma single-photon analysis :
//    ------------------------------
      singlePhoton();
//    --------------
      singlePhotonPMT();
//    -----------------

      checkSignal();
//    -------------

      bool bCheckAMPS = true;
      if( bCheckAMPS ) checkAMPSCorr();
//                     ---------------

      checkPhotonAngle();
//    ------------------

      checkPMTTimeDif();
//    -----------------

//--- monitoring histogram :
//    ----------------------
      int nRingSel = 0;
      for( ir=lRings_.begin(); ir!=lRings_.end(); ir++ ) {
        if( (*ir)->flag() ) nRingSel++;
      }
      xh = nRingSel;
      if( hist.hRC1537 ) hist.hRC1537->Fill( xh );
//hh                     ------------------------

      //^checkRingPH();
//    -------------
      //^checkNoise();
//    ------------

      //^ringSignal();
//    ------------

//--- set Ring Theta angle type :
//    ---------------------------
      for( ir=lRings_.begin(); ir!=lRings_.end(); ir++ ) {
	if( (*ir)->flag() ) (*ir)->setThetaType();
//                          ---------------------
      }
//--- Ring photons particle identification :
//    --------------------------------------
      //setBackgrType();
//    ---------------
      partChiSqIdent();
//    ----------------
      if( cons->likeType() == "ALL" ) partRingAllIdent();
//                                    ------------------
      if( cons->likeType() == "RING" ) partRingIdent();
//                                     ---------------

      //checkEventSignal();
//    ------------------

      bool bCheckRichMom = true;
      if( bCheckRichMom ) checkRichMom();
//    ----------------------------------

      bool bCheckPhiTheta = true;
      if( bCheckPhiTheta ) checkPhiTheta();
//    ------------------------------------

      bool bCheckRingFit = true;
      if( bCheckRingFit ) checkRingFit();
//    ----------------------------------

//--- conditional prints :
//    --------------------
      if( kPrint == 1 ) {
        cout << endl;
        cout << " Rings : " << nRingEv << endl;
      }
      if( kPrint == 2 ) print();

      xh = 59.5;
      if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                       ------------------------

      return;

  }


//===========================================================================
  void CsRCEventRings::flagPhoSignal() {
//--------------------------------------


//--- Paolo  -  May 1999,  rev. October 2000


//--- from "CsRCExecKeys"
      CsRCExeKeys *key = CsRCExeKeys::Instance();
      bool DoClu = key->DoClu();
      int kDoClu = 0;
      if( DoClu ) kDoClu = 0;
      if( !DoClu ) kDoClu = 1;

//--- from "CsRCRecConst"
      CsRCRecConst *cons = CsRCRecConst::Instance();
      int nPhotMin = cons->nPhotMin();
      float *nSigmaBf = cons->nSigmaBf();
      float sigBfCut = cons->sigBfCut();


//--- flag clusters of ring 'signal', for all particles :
//    ---------------------------------------------------
      float Cut = nSigmaBf[kDoClu] * sigBfCut;

//--- loop over part-photons :
      CsRCEventPartPhotons* paphos = CsRCEventPartPhotons::Instance();
      list<CsRCPartPhotons*> lPaPhotons = paphos->lPartPhotons();
      list<CsRCPartPhotons*>::iterator ipf;

      for( ipf=lPaPhotons.begin(); ipf!=lPaPhotons.end(); ipf++ ) {
	if( !(*ipf)->flag() ) continue;

        CsRCRing* ring = (*ipf)->pRing();
	if( ring == NULL ) continue;
        float theMN = ring->the() - Cut;
        float theMX = ring->the() + Cut;

//----- loop over 'photons' of each part-photons :
	list<CsRCPhoton*> lPhotons = (*ipf)->lPhotons();
	list<CsRCPhoton*>::iterator ih;
//$$        if( lPhotons.size() > nPhotMin )  {  /*  pro-memoria  */
        for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
          double thew = (*ih)->the();
          if( thew > theMN  &&  thew < theMX )  {
	    (*ipf)->setCluSignal( (*ih)->ptToClu() );
	  }
	}   /* end loop on photons */
//$$        }

      }   /* end loop on part-photons */

//--- loop over part-photons :
      list<CsRCPartPhotons*>::iterator ipg;
      for( ipf=lPaPhotons.begin(); ipf!=lPaPhotons.end(); ipf++ ) {
	if( !(*ipf)->flag() ) continue;

//----- loop over 'photons' of each part-photons :
	list<CsRCPhoton*> lPhotons = (*ipf)->lPhotons();
	list<CsRCPhoton*>::iterator ih;
        for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {

          bool bBreak = false;
          for( ipg=lPaPhotons.begin(); ipg!=lPaPhotons.end(); ipg++ ) {
	    if( !(*ipg)->flag() ) continue;

	    if( (*ipg) != (*ipf) ) {
	      list<CsRCCluster*> lCluSignal = (*ipg)->lCluSignal();
	      list<CsRCCluster*>::iterator ic;
	      for( ic=lCluSignal.begin(); ic!=lCluSignal.end(); ic++ ) {
		if( (*ih)->ptToClu() == (*ic) ) { bBreak = true;   break; }
              }
	    }
	    if( bBreak ) break;

	  }
          if( bBreak ) (*ih)->setFlag( false );

        }   /* end loop on photons */

      }   /* end loop on Part-photons */

  }


//===========================================================================
  void CsRCEventRings::ringSelection() {
//--------------------------------------


//--- ring selection (for re-analysis)
//    --------------------------------
//--- Paolo  -  November 2001
//      rev.    December 2002


      CsRCExeKeys *key = CsRCExeKeys::Instance();
      CsRCRecConst *cons = CsRCRecConst::Instance();

      CsRCHistos& hist = CsRCHistos::Ref();
      double xh, yh, wh;

      static bool firstCall = true;
      if( firstCall ) { 
        firstCall = false;

        key->acknoMethod( "CsRCEventRings::RingSelection" );
      }

      int nPhotMinRing = cons->nPhotMinRing();
      float momMinProc = cons->momMinProc();
      float momMaxProc = cons->momMaxProc();
      float CFRefInd = cons->CFRefInd();

//--- loop over rings :
//    -----------------
      CsRCEventRings* evrings = CsRCEventRings::Instance();
      list<CsRCRing*> lRings = evrings->lRings();
      list<CsRCRing*>::iterator ir;
      for( ir=lRings.begin(); ir!=lRings.end(); ir++ ) {
        if( !(*ir)->flag() ) continue;

        bool bSele = true;

//----- one ring only :
	//if( !( lRings.size() == 1 ) ) bSele = false;
//@@-----------------------------------------------

//----- particle momentum :
	CsRCParticle* part = (*ir)->pPartPhot()->pPart();
        float momMeas = part->mom();
        if( !( momMeas > momMinProc ) ) bSele = false;
//@@-------------------------------------------------
        //if( !( momMeas > 60.) ) bSele = false;
        //if( !( momMeas < 20.) ) bSele = false;
	//if( !( momMeas < 100.) ) bSele = false;
	//if( !( momMeas < 50.) ) bSele = false;
	if( !( momMeas < momMaxProc ) ) bSele = false;
//@@-------------------------------------------------

//----- photons in the ring :
	list<CsRCPhoton*> lPhotons = (*ir)->lPhotons();
        list<CsRCPhoton*>::iterator ih;
        int nPhoRing = lPhotons.size();
        //if( !( nPhoRing >= 4 ) ) bSele = false;
        //if( !( nPhoRing >= 10 ) ) bSele = false;
        if( !( nPhoRing >= nPhotMinRing ) ) bSele = false;
//@@-----------------------------------------------------
        //if( !( nPhoRing < 30 ) ) bSele = false;
//@@------------------------------------------

//----- particle on detectors :
        int kDetPart = (*ir)->pPartPhot()->kDetPart();
	double xPade = (*ir)->pPartPhot()->vPoPaDetW()[kDetPart].x();
	//if( !( xPade < 0. ) ) bSele = false;
	//if( !( xPade > 0. ) ) bSele = false;
	//if( !( xPade > 200.  &&  xPade < 300. ) ) bSele = false;
//@@-----------------------------------------------------------
	double yPade = (*ir)->pPartPhot()->vPoPaDetW()[kDetPart].y();
	//if( !( yPade < 0. ) ) bSele = false;      //   down
	//if( !( yPade > 0. ) ) bSele = false;      //   up
//@@---------------------------------------

	double ddPaDet = part->ddPaDet()[kDetPart];
	//if( !(ddPaDet > 400. ) )  bSele = false;
//@@-------------------------------------------

//----- particle at PT centre :
//----- WARNING: needs TRAFFIC initialization! (Extrapolate...)
	static double cov[15];
	for( int k=0; k<15; k++ ) cov[k] = 0.;
	int j = 1;
	for( int k=0; k<15; k+=j ) { cov[k] = 1.; j++; }
        CsHelix hxInW = CsHelix( part->vPosIn().x(), part->vPosIn().y(),
				 part->vPosIn().z(),
				 part->vDirIn().x()/part->vDirIn().z(),
				 part->vDirIn().y()/part->vDirIn().z(),
			 	 part->charge()/momMeas, cov );
	CsHelix hxPTg;
	double xPTg = 0.;
	double yPTg = 0.;
	double xyPTg = 0.;
	//if( hxInW.Extrapolate( -350., hxPTg ) ) {
	//  xPTg = hxPTg.getX();
	//  yPTg = hxPTg.getY();
	//cout << xPTg << "  " << yPTg << endl; 
        //  xyPTg = pow( xPTg/ 100., 2 ) + pow( yPTg/ 100., 2 );
	//}
	//if( xyPTg > 1. ) continue;
//----------------------------------------------------

//----- particle at entrance window :
	double ppxy;
        //ppxy = ( pow( part->vPosIn().x()/ 200., 2 ) +
      	//	   pow( part->vPosIn().y()/ 200., 2 ) );
	//if( !(ppxy > 1.) ) bSele = false;
        ppxy = sqrt( pow( part->vPosIn().x(), 2 ) +
      		     pow( part->vPosIn().y(), 2 ) );
	//if( !(ppxy > 50.) ) bSele = false;
//@@-------------------------------------
	double ttxy;
        //ttxy = ( pow( part->vDirIn().x()/part->vDirIn().z()/ 0.050, 2 ) +
      	//	   pow( part->vDirIn().y()/part->vDirIn().z()/ 0.050, 2 ) );
	//if( !(ttxy > 1.) ) bSele = false;
        ttxy = sqrt( pow( part->vDirIn().x()/part->vDirIn().z(), 2 ) +
      		     pow( part->vDirIn().y()/part->vDirIn().z(), 2 ) );
	ttxy *= 1000.;
	//if( !(ttxy > 100.) ) bSele = false;
//@@-------------------------------------

//----- divergent tracks :
        //if( !( part->vPosIn().y() > 0.  &&
        //       part->vDirIn().y()/part->vDirIn().z() > 0.  ||
        //       part->vPosIn().y() < 0.  &&
        //       part->vDirIn().y()/part->vDirIn().z() < 0. ) ) bSele = false;
//@@-----------------------------------------------------------------------

//----- track mom from ? :
//      NOT AVAILABLE on myfiles!
        CsTrack* track = part->pTrack();
	double zHelix0 = 0.;
	double zHelix1 = 0.;
        if( track != NULL ) {
          vector<CsHelix> vHelix = track->getHelices();
          zHelix0 = vHelix[0].getZ();
          zHelix1 = 0.;
          if( vHelix.size() == 2 ) zHelix1 = vHelix[1].getZ();
	}
        int kSM12 = -1;
	if( zHelix0 < 3500.  &&  zHelix1 < 18000.) kSM12 = 1;
	if( zHelix0 > 3500.  &&  zHelix1 > 18000.) kSM12 = 2;
	if( zHelix0 < 3500.  &&  zHelix1 > 18000.) kSM12 = 3;
        //if( !(kSM12 == 2) ) bSele = false;
//@@-------------------------------------

//----- photons on det up or down :
	int nPhoUp = 0;
        int nPhoDw = 0;
        for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
          int kDet = (*ih)->kDetClu();
          if( kDet == 0 ) nPhoUp++;
          if( kDet == 1 ) nPhoDw++;
//------- check on photon phi :
	  double phiDet = atan( (*ih)->ptToClu()->yc() / 
				(*ih)->ptToClu()->xc() );
	  phiDet *= 180./ 3.1415926;
	  if( phiDet < 0.) phiDet += 360.;
	  //cout << (*ih)->kPhot() << "  " << phiDet << "  "
	  //     << (*ih)->phiA() << "  " << (*ih)->phi() << endl; 
        }
	//if( !( nPhoUp < nPhoDw ) ) bSele = false;      //   down
	//if( !( nPhoDw < nPhoUp ) ) bSele = false;      //   up
//@@--------------------------------------------

//----- theta ring :
	//if( !((*ir)->the() >= 48. && (*ir)->the() <= 53.) ) bSele = false;
	//if( !((*ir)->the() >= 30. && (*ir)->the() <= 53.) ) bSele = false;
//@@---------------------------------------------------------------------

//----- select mass region :
	double cosReco = cos( (*ir)->the()/1000. );
	double recoMass = CFRefInd*cosReco*CFRefInd*cosReco - 1.;
	if( recoMass < 0. ) recoMass = 0.;
	recoMass = sqrt( recoMass );
	recoMass *= momMeas;
	//if( !(recoMass > 0.25  &&  recoMass < 0.40) ) bSele = false;
//@@---------------------------------------------------------------

	bool bThre = true;
	for( int kPaTy=2; kPaTy<=14; kPaTy+=3 ) {
	  double mass = cons->massPartv()[kPaTy];
	  double momThresh = mass / sqrt( CFRefInd*CFRefInd - 1.);
	  if( momMeas < momThresh ) {
	    bThre = false;
	    break;
	  }
	}
	//if( !bThre ) if( lPhotons.size() <= 6.) bSele = false;
//@@---------------------------------------------------------


	if( !bSele ) (*ir)->setFlag( false );
//                   -----------------------

	if( bSele ) {
	  xh = (*ir)->the();
	  yh = xPade;
	  if( hist.hRC3553 ) hist.hRC3553->Fill( xh, yh );
//                           ----------------------------
	}

	//^ringSignal( (*ir) );
//      -------------------

      }

//----- test 070828 ----------------------------------
//^^	testPhiKzero();      NOT WORKING!
//      --------------
//----- test 070828 ----------------------------------

  }


//===========================================================================
  void CsRCEventRings::singlePhoton() {
//-------------------------------------


//--- Paolo  -  December 1999,   rev.  October 2000
//              moved from CsRCEventRings.cc and rev. August 2002


//--- sigma single-photon analysis :
//    ------------------------------
//--- background vs theta, evaluation with fit :
//    ------------------------------------------

//--- from "CsRCRecConst"
      CsRCRecConst *cons = CsRCRecConst::Instance();
      int kOffHi = cons->kOffHi();
      double RadDeg = cons->RadDeg();
      float CFRefInd = cons->CFRefInd();

      CsRCHistos& hist = CsRCHistos::Ref();
      int kh;
      double xh, yh, wh;


      CsOpt* opt = CsOpt::Instance();
      static float momLimit = 50.;

//--- from rich1.options :
//    --------------------
      static bool firstCall = true;
      if( firstCall ) {
	firstCall = false;

	vector<float> vPar;
	bool boo = opt->CsOpt::getOpt( "RICHONE", "n-1Limits", vPar );
	if( boo ) momLimit  = vPar[0];
      }


//--- loop over rings :
//    -----------------
      CsRCEventRings* evrings = CsRCEventRings::Instance();
      list<CsRCRing*> lRings = evrings->lRings();
      list<CsRCRing*>::iterator ir;
      for( ir=lRings.begin(); ir!=lRings.end(); ir++ ) {     //   020814

        CsRCRing* ring = (*ir);                              //   020814
        if( !ring ) continue;                                //   020814
        CsRCPartPhotons* papho = ring->pPartPhot();          //   020814
        if( !papho ) continue;                               //   020814
//----- use selected rings only
	if( !(*ir)->flag() ) continue;                       //   021203
//      -----------------------------

        int kDetPart = papho->kDetPart();
        CsRCParticle* part = papho->pPart();
        double thePaDet = part->thPade()[kDetPart];
        thePaDet *= RadDeg;   //   to deg

        list<CsRCPhoton*> lPhotons = papho->lPhotons();
        int nPhoRing = ring->lPhotons().size();

        double theReco = ring->the();
        double betaReco = 0.;
	double betaw = cos( theReco/1000. ) * CFRefInd;
	if( betaw > 0.) betaReco = 1./ betaw;
        if( betaReco > 1. ) betaReco = 1.;

        if( betaReco > 0.99995 ) if( hist.hRC3562 ) hist.hRC3562->Fill( 0.5 );
//hh                                                -------------------------
        if( nPhoRing > 1 )  {   //   to avoid peaks at zero
        //if( nPhoRing >= 8 )  {   //   to work with a clean sample

          xh = thePaDet;
          if( hist.hRC3529 ) hist.hRC3529->Fill( xh );
//hh                         ------------------------

          kh = -1;
          if( thePaDet >  5.  &&  thePaDet <= 10. ) kh = 0;
          if( thePaDet > 10.  &&  thePaDet <= 15. ) kh = 1;
          if( thePaDet > 15.  &&  thePaDet <= 20. ) kh = 2;
          if( thePaDet > 20.  &&  thePaDet <= 25. ) kh = 3;
	  double ddPaDet = part->ddPaDet()[kDetPart];
	  int khP = -1;
          if( ddPaDet >   0.  &&  ddPaDet <= 100.) khP = 0;
          if( ddPaDet > 100.  &&  ddPaDet <= 200.) khP = 1;
          if( ddPaDet > 200.  &&  ddPaDet <= 300.) khP = 2;
          if( ddPaDet > 300.  &&  ddPaDet <= 400.) khP = 3;
          if( ddPaDet > 400.  &&  ddPaDet <= 500.) khP = 4;
          if( ddPaDet > 500.  &&  ddPaDet <= 600.) khP = 5;

	  double xPade = part->vPade()[part->kDetPart()].x();
	  double yPade = part->vPade()[part->kDetPart()].y();
//------------------------------------------------------------
//------- tests APV - 050220  &&  correlated backg : (50225)
//        ------------------------------------------
	  int khC = -1;
	  bool exeTest = false;
	  //bool exeTest = true;
//@@----------------------------
	  if( exeTest ) {
            const float xC0 = 317.5;   //   MWR
	    const float yC0 = 557.4;   //   MWR 
            const float dCx = 100.;
            const float dCy = 100.;
//@@------------------------------
	    bool cathPP = false;
            if( ( xPade >= xC0-dCx  &&  xPade <  xC0+dCx )  &&
                ( yPade >= yC0-dCy  &&  yPade <  yC0+dCy ) ) {
	      khC = 0;
	      cathPP = true;      // PP
	    }
	    bool cathNP = false;
            if(( xPade <= -xC0+dCx  &&  xPade > -xC0-dCx )  && 
               ( yPade >=  yC0-dCy  &&  yPade <  yC0+dCy ) ) {
	      khC = 1;
	      cathNP = true;      // NP
	    }
	    bool cathNN = false;
            if(( xPade <= -xC0+dCx  &&  xPade > -xC0-dCx )  && 
               ( yPade <= -yC0+dCy  &&  yPade > -yC0-dCy ) ) {
	      khC = 2;
	      cathNN = true;      // NN
	    }
  	    bool cathPN = false;
            if(( xPade >=  xC0-dCx  &&  xPade <  xC0+dCx )  && 
               ( yPade <= -yC0+dCy  &&  yPade > -yC0-dCy ) ) {
	      khC = 3;
	      cathPN = true;      // PN
	    }
	  }
//------------------------------------------------------------
	  bool ringPMTonly = ring->ringPMTonly();
	  bool ringAPVonly = ring->ringAPVonly();

//------- test 070801 ----------
          //std::cout << std::endl;
	  static const double TwoPi = cons->TwoPI();
	  CsRCDetectors* dets = CsRCDetectors::Instance();
	  static int nCathode = dets->nCathode();
	  int kCatHit[nCathode];
	  int nCatHit[nCathode];
	  double aCatHit[nCathode];
	  bool catHit = false;
	  if( theReco > 0.) {
	    papho->setLikeONLY( false );
	    catHit = papho->getHitCathodes( theReco, &kCatHit[0],
//          -----------------------------------------------------
		     &nCatHit[0], &aCatHit[0] );
	    if( catHit ) {
	      list<CsRCCathode*> lCathodes = dets->lCathodes();
	      list<CsRCCathode*>::iterator ic;
	      double normTo = 0.;
	      int nNorm = 0;
	      for( int kCat=0; kCat<nCathode; kCat++ ) {
		double theNorm = aCatHit[kCat]/TwoPi;
		normTo += theNorm;
		//if( nCatHit[kCat] > 2 ) 
		//std::cout << kCat << "  " << theNorm << std::endl;
		if( kCatHit[kCat] >= 0 ) {
		  nNorm++;
		  //std::cout << kCat << "  " << theNorm << std::endl;
		  xh = theReco;
		  yh = float( kCat );
		  wh = theNorm;
		  if( hist.hRC3677 ) hist.hRC3677->Fill( xh, yh, wh );
//hh                                 --------------------------------
		  //break;
		}
		if( kCatHit[kCat] == 0 ) {
		  //if( nCatHit[kCat] == 0 ) 
		  //std::cout << kCat << "  " << theNorm << std::endl;
		  xh = theReco;
		  yh = float( kCat );
		  if( hist.hRC3678 ) hist.hRC3678->Fill( xh, yh );
//hh                                 ----------------------------
		}
	      }
	      //if( nNorm > 1 ) {
	      //std::cout << normTo << " - ";
	      //for( int k=0; k<16; k++) if( kCatHit[k] >= 0 ) 
	      //std::cout << aCatHit[k]/TwoPi << "  " << nCatHit[k] << " - ";
	      //std::cout << std::endl;
	      //}
	    }
	  }
//------- test 070801 ----------

          list<CsRCPhoton*>::iterator ih;
          for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {

  	    double thePho = (*ih)->the();
	    double phiPho = (*ih)->phi();
	    double phiPhoA = (*ih)->phiA();
	    double thetaBM = (*ih)->theB() - (*ih)->theM();

//--------- (theta photon-theta ring) vs phi & thPade :
//          -------------------------------------------
            xh = thePho - theReco;
            yh = phiPho;
            if( kh >= 0  &&  kh <= 3 ) {
              if( hist.vRC3520[kh] ) hist.vRC3520[kh]->Fill( xh, yh );
//hh                                 --------------------------------
	    }

//--------- Added   081216
//--------- (theta photon-theta ring) vs phi & thPade :
//          -------------------------------------------
            xh = thePho - theReco;
            yh = phiPho;
	    int kHH = kh;
            if( kHH >= 0  &&  kHH < int( hist.vRC6550.size() ) ) {
              if( hist.vRC6550[kHH] ) hist.vRC6550[kHH]->Fill( xh, yh );
//hh                                  ---------------------------------
	    }
	    if( (*ih)->isPMT() ) {
	      int kHH = kh + 4;
              if( kHH >= 0  &&  kHH < int( hist.vRC6550.size() ) ) {
                if( hist.vRC6550[kHH] ) hist.vRC6550[kHH]->Fill( xh, yh );
//hh                                    ---------------------------------
	      }
	    }
	    if( (*ih)->isAPV() ) {
	      int kHH = kh + 8;
              if( kHH >= 0  &&  kHH < int( hist.vRC6550.size() ) ) {
                if( hist.vRC6550[kHH] ) hist.vRC6550[kHH]->Fill( xh, yh );
//hh                                    ---------------------------------
	      }
	    }
//---------------   081216

//--------- (theta photon-theta ring),  beta close to 1 :
//          ---------------------------------------------
            if( betaReco > 0.99995 ) {

              xh = thePho - theReco;
              if( hist.hRC3525 ) hist.hRC3525->Fill( xh );
//hh                             ------------------------
//----------- (theta photon-theta ring) vs phi photon :
//            -----------------------------------------
              xh = thePho - theReco;
              yh = phiPho;
              if( hist.hRC3526 ) hist.hRC3526->Fill( xh, yh );
//hh                             ----------------------------
//----------- (theta photon-theta ring) vs phiA photon :
//            ------------------------------------------
              xh = thePho - theReco;
              yh = phiPhoA;
              if( hist.hRC3531 ) hist.hRC3531->Fill( xh, yh );
//hh                             ----------------------------
	    }

//--------- (theta photon-theta ring) vs theta ring :
//          -----------------------------------------
            xh = thePho - theReco;
            yh = theReco;
            if( hist.hRC3527 ) hist.hRC3527->Fill( xh, yh );
//hh                           ----------------------------
//--------- (theta photon-theta ring) vs theta part-mirror :
//          ------------------------------------------------
            xh = thePho - theReco;
            yh = part->thPamir() * RadDeg;
            if( hist.hRC3528 ) hist.hRC3528->Fill( xh, yh );
//hh                           ----------------------------

//--------- phi photon vs thetaBM :
//          -----------------------
	    xh = thetaBM;
            yh = phiPho;
            if( hist.hRC3530 ) hist.hRC3530->Fill( xh, yh );
//hh                           ----------------------------

//--------- theta ring vs theta photon ( Like BACKGROUND ):
//          -----------------------------------------------
            xh = thePho;
            yh = theReco;
            if( hist.hRC3561 ) hist.hRC3561->Fill( xh, yh );
//hh                           ----------------------------

	    if( CsRichOne::Instance()->UpRICHJob() ) {
  	      if( ringAPVonly  &&  (*ih)->isAPV() ) {
                xh = thePho;
                yh = theReco;
	        if( khP >= 0  &&  khP <= 5 ) {
	          if( hist.vRC3680[khP] ) hist.vRC3680[khP]->Fill( xh, yh );
//hh                                      ---------------------------------
		}
	      }
	      if( ringPMTonly  &&  (*ih)->isPMT() ) {
                xh = thePho;
                yh = theReco;
	        if( khP >= 0  &&  khP <= 5 ) {
	          if( hist.vRC3690[khP] ) hist.vRC3690[khP]->Fill( xh, yh );
//hh                                      ---------------------------------
		}
	      }
//----------- test 070801 ----------
	      if( catHit ) {
		int kCat = (*ih)->ptToClu()->ic();
		xh = thePho;
		yh = theReco;
		if( hist.vRC3660[kCat] ) hist.vRC3660[kCat]->Fill( xh, yh );
//hh                                     ----------------------------------
	      }
//----------- test 070801 ----------
	    } else {
              xh = thePho;
              yh = theReco;
	      if( khP >= 0  &&  khP <= 5 ) {
	        if( hist.vRC3590[khP] ) hist.vRC3590[khP]->Fill( xh, yh );
//hh                                    ---------------------------------
	      }
	    }

//--------- theta photon vs theta ring ( tests APV ): 050220
//          -----------------------------------------
            xh = thePho;
            yh = theReco;
	    if( khC >= 0  &&  khC <= 3 ) {
	      if( hist.vRC3565[khC] ) hist.vRC3565[khC]->Fill( xh, yh );
//hh                                  ---------------------------------
	    }

//--------- ( n-1 ) as 1/(beta*cos(theta-photon)) - 1
//          -----------------------------------------
//          assume pion mass (and particle momentum limit)
	    if( part->mom() < momLimit ) {
//@@--------------------------------------
              //^double thePhow = (*ih)->the();
  	      double thePhow = (*ih)->the0();
	      float nm1 = fabs( part->mom() * cos( thePhow/1000. ) );
	      if( nm1 > 0. )
	        nm1 = sqrt( 0.01948 + part->mom()*part->mom() ) / nm1;
	      xh = nm1 - 1.;
	      if( hist.hRC3559 ) hist.hRC3559->Fill( xh );
//hh                             ------------------------
	    }

          }   /* end of loop on photons */

//------- counts for theta photon vs theta ring :
//        ---------------------------------------
	  xh = theReco;
	  if( hist.hRC3560 ) hist.hRC3560->Fill( xh );
//hh                         ------------------------

	  // Infra gets coral to crash (Y.B.)
          // Bug corrected and fill restored (P.S.) 040227
	  if( CsRichOne::Instance()->UpRICHJob() ) {
	    if( ringAPVonly ) {
	      xh = theReco;
	      yh = ddPaDet;
              if( hist.hRC3687 ) hist.hRC3687->Fill( xh, yh );
//hh                             ----------------------------
	      xh = xPade;
	      yh = yPade;
              if( hist.hRC3688 ) hist.hRC3688->Fill( xh, yh );
//hh                             ----------------------------
	    }
	    if( ringPMTonly ) {
	      xh = theReco;
	      yh = ddPaDet;
              if( hist.hRC3697 ) hist.hRC3697->Fill( xh, yh );
//hh                             ----------------------------
	      xh = xPade;
	      yh = yPade;
              if( hist.hRC3698 ) hist.hRC3698->Fill( xh, yh );
//hh                             ----------------------------
	    }
	  } else {
	    xh = theReco;
	    yh = ddPaDet;
	    if( hist.hRC3597 ) hist.hRC3597->Fill( xh, yh );
//hh                           ----------------------------
	  }

//------- counts for theta photon vs theta ring ( tests APV ): 050220
//        ----------------------------------------------------
	  if( khC >= 0  &&  khC <= 3 ) {
	    xh = theReco;
	    yh = float( khC ) + 0.5;
	    if( hist.hRC3569 ) hist.hRC3569->Fill( xh, yh );
//hh                           ----------------------------
	  }

        }   /* End of if on nPhoRing */

      }


//-------------------------------------------------------------------
//--- 'Correlated' background tests : (050225)
//    -------------------------------
//    WORNING : this tests use the same histo.s as the APV tests!
//--- loop over part-photons :
//    ------------------------
      bool exeTest = false;
      //bool exeTest = true;
//@@----------------------
      if( exeTest ) {

        CsRCEventPartPhotons* evepapho = CsRCEventPartPhotons::Instance();
        std::list<CsRCPartPhotons*> lPartPhotons = evepapho->lPartPhotons();
        std::list<CsRCPartPhotons*>::iterator ipf;
        for( ipf=lPartPhotons.begin(); ipf!=lPartPhotons.end(); ipf++ ) {
	  if( !(*ipf)->flag() ) continue;
	  int khC = -1;
	  double theReco = 0.;
          int kDetPart = (*ipf)->kDetPart();
	  double ddPaDet = (*ipf)->pPart()->ddPaDet()[kDetPart];
          if( ddPaDet > 180. ) {
            int kPart = (*ipf)->pPart()->kPart();
            if( kPart >= 100 ) {
              list<CsRCPartPhotons*>::iterator ipg;
	      for( ipg=lPartPhotons.begin(); ipg!=lPartPhotons.end(); ipg++) {
		if( !(*ipg)->flag() ) continue;
	        if( (*ipg)->pPart()->kPart() == kPart-100 ) {
		  if( (*ipg)->pRing() ) theReco = (*ipg)->pRing()->the();
		  //if( (*ipg)->thetaLikeSet() ) theReco =(*ipg)->thetaLike();
                  khC = 1;
	          break;
	        }
	      }
	    }
	    else  {
	      if( (*ipf)->pRing() ) theReco = (*ipf)->pRing()->the();
  	      //if( (*ipf)->thetaLikeSet() ) theReco = (*ipf)->thetaLike();
	      khC = 0;
	    }
          }
	  if( khC >= 0  &&  khC <= 1 ) {
	    xh = theReco;
	    yh = float( khC ) + 0.5;
	    if( hist.hRC3569 ) hist.hRC3569->Fill( xh, yh );
//hh                           ----------------------------
	    list<CsRCPhoton*> lPhotons = (*ipf)->lPhotons();
	    list<CsRCPhoton*>::iterator ih;
	    for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	      double thePho = (*ih)->the();
              xh = thePho;
              yh = theReco;
              if( hist.vRC3565[khC] ) hist.vRC3565[khC]->Fill( xh, yh );
//hh                                  ---------------------------------
	    }
	  }
	}

      }
//-------------------------------------------------------------------

  }


//===========================================================================
  void CsRCEventRings::checkSignal() {
//------------------------------------

//- Paolo  -  September 2002
//  revised   August 2005


    CsRCRecConst* cons = CsRCRecConst::Instance();
    int kOffHi = cons->kOffHi();
    float CFRefInd = cons->CFRefInd();
    double* massPartv = cons->massPartv();

    CsRCHistos& hist = CsRCHistos::Ref();
    int kh;
    double xh, yh, wh;

//- LOOP over RINGS :
//  -----------------
    CsRCEventRings* evrings = CsRCEventRings::Instance();
    list<CsRCRing*> lRings = evrings->lRings();
    list<CsRCRing*>::iterator ir;
    for( ir=lRings.begin(); ir!=lRings.end(); ir++ ) {

      CsRCRing* ring = (*ir);
      if( !ring ) continue;
      CsRCPartPhotons* papho = ring->pPartPhot();
      if( !papho ) continue;

//    check on ring flag
      if( !ring->flag() ) continue;
//    ----------------------------

      CsRCParticle* part = ring->pPartPhot()->pPart();
      float momPart = part->mom();

      //^list<CsRCPhoton*> lPhotons = papho->lPhotons();
      list<CsRCPhoton*> lPhotons = ring->lPhotons();
      list<CsRCPhoton*>::iterator ih;
      if( lPhotons.size() <= 2 ) continue;
//@@-------------------------------------

      double theReco = ring->the();
      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	double thePhot = (*ih)->the();
	xh = thePhot - theReco;
	yh = momPart;
	if( hist.hRC3533 ) hist.hRC3533->Fill( xh, yh );
//hh                       ----------------------------
      }

//--- check for pion, kaon or proton only
//    -----------------------------------
      int kIpo = 0;
      for( int kPaTy=8; kPaTy<=14; kPaTy+=3 ) {

	double theIpo = papho->thetaIpo( kPaTy );
	if( theIpo > 0. ) {

	  for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	    double thePhot = (*ih)->the();
	    xh = thePhot - theIpo;
	    yh = momPart;
	    if( hist.vRC3540[kIpo] ) hist.vRC3540[kIpo]->Fill( xh, yh );
//hh                                 ----------------------------------
	  }
	  xh = float( kIpo ) + 0.5;
	  yh = momPart;
	  if( hist.vRC3540[3] ) hist.vRC3540[3]->Fill( xh, yh );
//hh                            -------------------------------
	}
	kIpo++;

      }   /* end of loop on kPaTy */

    }

  } 


//===========================================================================
  void CsRCEventRings::checkAMPSCorr() {
//--------------------------------------

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

    static std::vector<CsHist1D*> vRC3395;
    static std::vector<CsHist2D*> vRC6700;
    static std::vector<CsHist2D*> vRC6710;
    static std::vector<CsHist2D*> vRC6720;

    CsRCHistos& hist = CsRCHistos::Ref();
    int kh;
    double xh, yh, wh;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;
      for( int kh=0; kh<4; kh++ ) vRC3395.push_back( NULL );
      for( int kh=0; kh<9; kh++ ) vRC6700.push_back( NULL );
      for( int kh=0; kh<9; kh++ ) vRC6710.push_back( NULL );
      for( int kh=0; kh<9; kh++ ) vRC6720.push_back( NULL );
      if( CsRCHistos::Ref().bookHis() &&
	  CsRCHistos::Ref().levelBk() >= 2 ) {
	CsHistograms::SetCurrentPath("/RICH");
	vRC3395.clear();
	int kHist = 0;
	string hTitle;
        hTitle = "thePho-thePi JuUp";
	stringstream hN3396;
	kHist = kOffHi + 3396;
	hN3396 << kHist;
	vRC3395.push_back( new CsHist1D( hN3396.str(), hTitle,
					 100, -10., 10. ) );
        hTitle = "thePho-thePi SaUp";
	stringstream hN3397;
	kHist = kOffHi + 3397;
	hN3397 << kHist;
	vRC3395.push_back( new CsHist1D( hN3397.str(), hTitle,
					 100, -10., 10. ) );
        hTitle = "thePho-thePi JuDw";
	stringstream hN3398;
	kHist = kOffHi + 3398;
	hN3398 << kHist;
	vRC3395.push_back( new CsHist1D( hN3398.str(), hTitle,
					 100, -10., 10. ) );
        hTitle = "thePho-thePi SaDw";
	stringstream hN3399;
	kHist = kOffHi + 3399;
	hN3399 << kHist;
	vRC3395.push_back( new CsHist1D( hN3399.str(), hTitle,
					 100, -10., 10. ) );

	vRC6700.clear();
	hTitle = " ";
	for( int kh=0; kh<5; kh++ ) {
	  kHist = kOffHi + 6700 + kh + 1;
	  stringstream hN6700;
	  hN6700 << kHist;
	  vRC6700.push_back( new CsHist2D( hN6700.str(), hTitle,
					   100, -10., 10., 90, 0., 360. ) );
	}
	vRC6710.clear();
	hTitle = " ";
	for( int kh=0; kh<5; kh++ ) {
	  kHist = kOffHi + 6710 + kh + 1;
	  stringstream hN6710;
	  hN6710 << kHist;
	  vRC6710.push_back( new CsHist2D( hN6710.str(), hTitle,
					   100, -10., 10., 90, 0., 360. ) );
	}
	vRC6720.clear();
	hTitle = " ";
	for( int kh=0; kh<8; kh++ ) {
	  kHist = kOffHi + 6720 + kh + 1;
	  stringstream hN6720;
	  hN6720 << kHist;
	  vRC6720.push_back( new CsHist2D( hN6720.str(), hTitle,
					   100, -10., 10., 90, 0., 360. ) );
	}
        CsHistograms::SetCurrentPath("/");
      }
    }

    double* probaLK;

//- LOOP over RINGS :
//  -----------------
    CsRCEventRings* evrings = CsRCEventRings::Instance();
    list<CsRCRing*> lRings = evrings->lRings();
    list<CsRCRing*>::iterator ir;
    for( ir=lRings.begin(); ir!=lRings.end(); ir++ ) {

      CsRCRing* ring = (*ir);
      if( !ring ) continue;
      CsRCPartPhotons* papho = ring->pPartPhot();
      if( !papho ) continue;
//    check on ring flag
      if( !ring->flag() ) continue;
//    ----------------------------

      CsRCParticle* part = ring->pPartPhot()->pPart();
      float momPart = part->mom();

//--- Select ID pions
      double likeBkg = papho->probaLKBgAll();
      probaLK = papho->probaLKAll();
      double likeElec = probaLK[ 2];
      double likeMuon = probaLK[ 5];
      double likePion = probaLK[ 8];
      double likeKaon = probaLK[11];
      double likeProton = probaLK[14];
      bool IDpion = false;
      if( likePion > 0.  &&  likePion > likeBkg  &&
	  likePion > likeElec  &&  likePion > likeMuon  &&
	  likePion > likeKaon  &&  likePion > likeProton ) IDpion = true;
      if( !IDpion ) continue;
//@@  ----------------------

//--- Select full (or split) rings
      int kRingLoc = ring->kRingLoc();
      //if( kRingLoc == 0 ) continue;
      //if( kRingLoc != 0 ) continue;
//@@  ----------------------------

      //list<CsRCPhoton*> lPhotons = papho->lPhotons();
      list<CsRCPhoton*> lPhotons = ring->lPhotons();
//@@-----------------------------------------------
      list<CsRCPhoton*>::iterator ih;
      if( lPhotons.size() <= 2 ) continue;
//@@-------------------------------------

      double theReco = ring->the();

//--- check for one from elctron, muon, pion, kaon or proton
//    ------------------------------------------------------
      int kIpo = 0;
      for( int kPaTy=2; kPaTy<=14; kPaTy+=3 ) {
	if( kPaTy != 8 ) continue;               // select pion

	double theIpo = papho->thetaIpo( kPaTy );
	//double theIpoVS = papho->thetaIpoVS( kPaTy );   // NO!
	if( theIpo <= 0. ) continue;

	for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	  //double theIpow = theIpo;
	  double theIpow = theReco;
//@@------------------------------
          //if( (*ih)->isPMT() ) theIpow = theIpoVS;   // NO!
          if( theIpow <= 0. ) continue;
	  if( dets->nCatPMT().size() <= 0 ) continue;

          double thePhot = (*ih)->the();
	  int cClu = (*ih)->ptToClu()->ic();
          xh = thePhot - theIpow;
          if( cClu == dets->nCatPMT()[0] ) {
            int kHH = 0;
            if( vRC3395[kHH] ) vRC3395[kHH]->Fill( xh );
//hh                           ------------------------
          }
          if( cClu == dets->nCatPMT()[1] ) {
            int kHH = 1;
            if( vRC3395[kHH] ) vRC3395[kHH]->Fill( xh );
//hh                           ------------------------
          }
          if( cClu == dets->nCatPMT()[2] ) {
            int kHH = 2;
            if( vRC3395[kHH] ) vRC3395[kHH]->Fill( xh );
//hh                           ------------------------
          }
          if( cClu == dets->nCatPMT()[3] ) {
            int kHH = 3;
            if( vRC3395[kHH] ) vRC3395[kHH]->Fill( xh );
//hh                           ------------------------
	  }
        }
        kIpo++;
      }   /* end of loop on kPaTy */

    }

//- Added 080918
    for( ir=lRings.begin(); ir!=lRings.end(); ir++ ) {

      CsRCRing* ring = (*ir);
      if( !ring ) continue;
      CsRCPartPhotons* papho = ring->pPartPhot();
      if( !papho ) continue;
//    check on ring flag
      if( !ring->flag() ) continue;
//    ----------------------------

      CsRCParticle* part = ring->pPartPhot()->pPart();
      float momPart = part->mom();

      list<CsRCPhoton*> lPhotons = ring->lPhotons();
//@@-----------------------------------------------
      list<CsRCPhoton*>::iterator ih;
      if( lPhotons.size() <= 2 ) continue;
//@@-------------------------------------
      double mom0 = 0.;
      double dmom = 10.;
      int kHH = -1;
      for( int kp=0; kp<5; kp++ ) {
	if( momPart >= mom0  &&  momPart < (mom0+dmom) ) {
	  kHH = kp;
	  break;
	}
	mom0 += dmom;
      }
      double xPade = part->vPade()[part->kDetPart()].x();
      double yPade = part->vPade()[part->kDetPart()].y();
      int kHK = -1;
      double tgxPade = papho->vDcPaReflW()[part->kDetPart()].x()/
	papho->vDcPaReflW()[part->kDetPart()].z();
      double tgx0 = -0.40;
      double dtgx =  0.10;
      for( int kt=0; kt<8; kt++ ) {
	if( tgxPade >= tgx0  &&  tgxPade < (tgx0+dtgx) ) {
	  kHK = kt;
	  break;
	}
	tgx0 += dtgx;
      }
      //double theReco = papho->thetaUVtoVS( ring->the() );
      double theReco = 0.;
      int nPMT = 0;
      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	if( !(*ih)->isPMT() ) continue;
	theReco += (*ih)->theB();
	nPMT++;
      }
      if( nPMT > 0 ) theReco /= nPMT;
      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	if( !(*ih)->isPMT() ) continue;
	double thePhoP = (*ih)->theB();
	double thePhoPO = (*ih)->theM();
	double phiPho = (*ih)->phi();
	double phiPhoA = (*ih)->phiA();
	//if( thePhoP == thePhoPO )
	//std::cout << thePhoP << "  " << thePhoPO << std::endl;
	//std::cout << phiPho << "  " << phiPhoA << std::end;
	if( kHH >= 0  &&  kHH < int( vRC6700.size() ) ) {
	  xh = thePhoPO - thePhoP;
	  //yh = phiPho;
	  yh = phiPhoA;
	  if( vRC6700[kHH] ) vRC6700[kHH]->Fill( xh, yh );
//hh                         ----------------------------
	}
	if( kHH >= 0  &&  kHH < int( vRC6710.size() ) ) {
	  xh = thePhoP - theReco;
	  //xh = thePhoPO - theReco;
	  //yh = phiPho;
	  yh = phiPhoA;
	  if( vRC6710[kHH] ) vRC6710[kHH]->Fill( xh, yh );
//hh                         ----------------------------
	}
	if( kHK >= 0  &&  kHK < int( vRC6720.size() ) ) {
	  xh = thePhoP - theReco;
	  //xh = thePhoPO - theReco;
	  //yh = phiPho;
	  yh = phiPhoA;
	  if( vRC6720[kHK] ) vRC6720[kHK]->Fill( xh, yh );
//hh                         ----------------------------
	}
      }

    }

  }


//===========================================================================
  void CsRCEventRings::checkNoise() {
//-----------------------------------


//- Paolo  -  June 2001


    int nCathode = CsRCDetectors::Instance()->nCathode();

    //list<CsRCPad*> lNoisePads = CsRCEventPads::Instance()->lPads();
    list<CsRCPad*> lPadsAll = CsRCEventPads::Instance()->lPads();

//- LOOP over RINGS found :
//  -----------------------
    int nSgnPadCat[nCathode];
    for( int k=0; k<nCathode; k++ ) nSgnPadCat[k] = 0;
    list<CsRCRing*>::iterator ir;
    for( ir=lRings_.begin(); ir!=lRings_.end(); ir++ ) {
      if( !(*ir)->flag() ) continue;

      list<CsRCPhoton*> lPhotons = (*ir)->lPhotons();
      list<CsRCPhoton*>::iterator ih;
      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {

	list<CsRCPad*> lPads = (*ih)->ptToClu()->lPads();
        list<CsRCPad*>::iterator ip;
        for( ip=lPads.begin(); ip!=lPads.end(); ip++ ) {
	  nSgnPadCat[(*ip)->ic()]++;
	  //const CsRCPad pad = *(*ip);
	  //lNoisePads.remove( pad );
	}
      }
    }

    int nTotPadCat[nCathode];
    for( int k=0; k<nCathode; k++ ) nTotPadCat[k] = 0;
    //int nNoiPadCat[nCathode];
    //for( int k=0; k<nCathode; k++ ) nNoiPadCat[k] = 0;
    list<CsRCPad*>::iterator ip;
    //for( ip=lNoisePads.begin(); ip!=lNoisePads.end(); ip++ ) {
    //  nNoiPadCat[(*ip)->ic()]++;
    //}
    for( ip=lPadsAll.begin(); ip!=lPadsAll.end(); ip++ ) {
      if( !(*ip)->flag() ) continue;
      nTotPadCat[(*ip)->ic()]++;
    }

    CsRCHistos& hist = CsRCHistos::Ref();
    for( int kCat=0; kCat<nCathode; kCat++ ) {
      float xh = nTotPadCat[kCat] - nSgnPadCat[kCat];
      float yh = kCat;
      if( hist.hRC1523 ) hist.hRC1523->Fill( xh, yh );
    }

  }


//===========================================================================
  void CsRCEventRings::checkRingPH() {
//------------------------------------


//- Paolo  -  October 2001


    CsRCHistos& hist = CsRCHistos::Ref();

//- LOOP over RINGS found :
//  -----------------------
    list<CsRCRing*>::iterator ir;
    for( ir=lRings_.begin(); ir!=lRings_.end(); ir++ ) {
      if( !(*ir)->flag() ) continue;
      double theta = (*ir)->the();

      list<CsRCPhoton*> lPhotons = (*ir)->lPhotons();
      list<CsRCPhoton*>::iterator ih;
      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	//if( !(*ih)->flag() ) continue;

	list<CsRCPad*> lPads = (*ih)->ptToClu()->lPads();
        list<CsRCPad*>::iterator ip;
        for( ip=lPads.begin(); ip!=lPads.end(); ip++ ) {

	  float xh = (*ip)->PH();
	  float yh = theta;
	  if( hist.hRC1529 ) hist.hRC1529->Fill( xh, yh );
//        -----------------------------------------------

	}
      }
    }

  }


//===========================================================================
  void CsRCEventRings::ringSignal() {
//-----------------------------------


//--- Paolo  - April  2002


      CsRCHistos& hist = CsRCHistos::Ref();
      int kh;
      double xh, yh, wh;


//--- integrated ring signal :
//    ------------------------
      list<CsRCRing*>::iterator ir;
      for( ir=lRings_.begin(); ir!=lRings_.end(); ir++ ) {
	if( !(*ir)->flag() ) continue;

	CsRCParticle* part = (*ir)->pPartPhot()->pPart();
	double xPade = part->vPade()[part->kDetPart()].x();
	double yPade = part->vPade()[part->kDetPart()].y();
	double radDeg = CsRCRecConst::Instance()->RadDeg();
//----- provisional cuts for test (020414):
//      -----------------------------------
        //if( part->mom() > 20. ) continue;
        //if( part->mom() < 60. ) continue;
	//if( yPade > -450. ) continue;         //   lower RICH half
	//if( yPade <  450. ) continue;         //   upper RICH half
	//double rPade = sqrt( xPade*xPade + 
	//		     (fabs(yPade)-250.)*(fabs(yPade)-250.) );
	//if( rPade < 300. ) continue;
	//cout << part->mom() << "  " << yPade << endl;
//      -----------------------------------
	int kHist = -1;
	if( xPade > -400.  &&  xPade <= -200.) kHist = 0;
	if( xPade > -200.  &&  xPade <=    0.) kHist = 1;
	if( xPade >    0.  &&  xPade <=  200.) kHist = 2;
	if( xPade >  200.  &&  xPade <=  400.) kHist = 3;
	if( kHist < 0 ) continue;
        list<CsRCPhoton*> lPhotons = (*ir)->lPhotons();
        list<CsRCPhoton*>::iterator ih;
        if( lPhotons.size() > 1 )  {   
          for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
            xh = (*ih)->the() * cos( (*ih)->phiA()*radDeg );
            yh = (*ih)->the() * sin( (*ih)->phiA()*radDeg );
            if( hist.vRC3570[kHist] ) hist.vRC3570[kHist]->Fill( xh, yh );
//hh                                  -----------------------------------
          }   /* end of loop on photons */
	}
      }

  }


//===========================================================================
  void CsRCEventRings::ringSignal( CsRCRing* ring ) {
//---------------------------------------------------


//--- Paolo  - April  2002


      CsRCHistos& hist = CsRCHistos::Ref();
      int kh;
      double xh, yh, wh;

      if( !ring->flag() ) return;

      CsRCParticle* part = ring->pPartPhot()->pPart();
      double xPade = part->vPade()[part->kDetPart()].x();
      double yPade = part->vPade()[part->kDetPart()].y();
      double radDeg = CsRCRecConst::Instance()->RadDeg();

      list<CsRCPhoton*>::iterator ih;

//--- raw signal on detector planes :
//    -------------------------------
      CsRCPartPhotons* papho = ring->pPartPhot();
      if( papho ) {
        list<CsRCPhoton*> lPhotons = papho->lPhotons();
        for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
          double xPade = papho->pPart()->vPade()[(*ih)->kDetClu()].x();
          double yPade = papho->pPart()->vPade()[(*ih)->kDetClu()].y();
          double xClu = (*ih)->ptToClu()->xc();
          double yClu = (*ih)->ptToClu()->yc();
          double dd = sqrt( (xClu-xPade)*(xClu-xPade) + 
      		            (yClu-yPade)*(yClu-yPade) );
          xh = dd;
          if( hist.hRC3557 ) hist.hRC3557->Fill( xh );
//hh                         ------------------------
          xh = xClu - xPade;
          yh = yClu - yPade;
          if( hist.hRC3558 ) hist.hRC3558->Fill( xh, yh );
//hh                         ----------------------------
        }   /* end of loop on photons */
      }

//--- integrated ring signal :
//    ------------------------
      int kHist = -1;
      if( xPade > -400.  &&  xPade <= -200.) kHist = 0;
      if( xPade > -200.  &&  xPade <=    0.) kHist = 1;
      if( xPade >    0.  &&  xPade <=  200.) kHist = 2;
      if( xPade >  200.  &&  xPade <=  400.) kHist = 3;
      if( kHist >= 0 ) { 
        list<CsRCPhoton*> lPhotons = ring->lPhotons();
        if( lPhotons.size() >= 5 )  {   
          for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
            xh = (*ih)->the() * cos( (*ih)->phiA()*radDeg );
            yh = (*ih)->the() * sin( (*ih)->phiA()*radDeg );
            if( hist.vRC3580[kHist] ) hist.vRC3580[kHist]->Fill( xh, yh );
//hh                                  -----------------------------------
          }   /* end of loop on photons */
	}
      }

  }


//===========================================================================
  void CsRCEventRings::setBackgrType() {
//--------------------------------------


//--- Paolo  -  Deceember 2002

//--- not used


  }


  extern "C" {
    float prob_( const float&, const int& );
    float erf_( const float& );
  }

//===========================================================================
  void CsRCEventRings::partChiSqIdent() {
//---------------------------------------


//--- Ring particle identification
//    ----------------------------
//--- Paolo  -  December 2002
//    revised   August 2005


      CsRCExeKeys* key = CsRCExeKeys::Instance();

      static bool firstCall = true;
      if( firstCall ) { 
        firstCall = false;

        key->acknoMethod( "CsRCEventRings::PartChiSqIdent" );
      }

      CsRCRecConst* cons = CsRCRecConst::Instance();
      int kOffHi = cons->kOffHi();
      float CFRefInd = cons->CFRefInd();
      double* massPartv = cons->massPartv();
      int nPhotMinRing = cons->nPhotMinRing();

      int nProb = cons->outBufferSize();

      CsRCHistos& hist = CsRCHistos::Ref();
      double xh, yh, wh;


//--- CHISQ PARTICLE IDENTIFICATION :
//    -------------------------------
      int nPIDEv = 0;

//--- loop over rings :
//    -----------------
      list<CsRCRing*> lRings = lRings_;
      list<CsRCRing*>::iterator ir;
      for( ir=lRings.begin(); ir!=lRings.end(); ir++ ) {
        if( !(*ir)->flag() ) continue;

//----- define SELECTED rings :
//      -----------------------
        bool bSele = true;

//----- enough photons in the ring :
//      ---------------------------
	list<CsRCPhoton*> lPhotons = (*ir)->lPhotons();
        list<CsRCPhoton*>::iterator ih;
        int nPhoRing = lPhotons.size();
        if( nPhoRing < nPhotMinRing ) bSele = false;
//@@                                  -------------

//----- particle momentum over muon threshold (1.91 GeV/c) :
//      ----------------------------------------------------
	CsRCParticle* part = (*ir)->pPartPhot()->pPart();
        float momMeas = part->mom();
        //if( momMeas < 1.91 ) bSele = false;
//@@--------------------------------------

//----- reconstructed Cerenkov angle :
//      ------------------------------
	double theReco = (*ir)->theReco();

	double theRecoR = (*ir)->the();
	double theRecoWg = (*ir)->getThetaWgav();
	double theRecoFt = (*ir)->thetaRFit();

	xh = theRecoR - theRecoWg;
	if( hist.hRC1630 ) hist.hRC1630-> Fill( xh );
//hh                       -------------------------
	if( theRecoR != theRecoFt ) {
  	  xh = theRecoR - theRecoFt;
	  if( hist.hRC1631 ) hist.hRC1631-> Fill( xh );
//hh                         -------------------------
	}

//----- process selected rings :
//      ------------------------
        if( bSele ) {

	  static bool checkChi = true;
	  if( checkChi ) (*ir)->checkChiSquare( theReco );
//                       --------------------------------

	  int kTyLL = 2;
//@@-------------------
          int kFigMe = 0;
	  double qSquare[31];
          double probaQS[31];
          double fgMeritMx = 10000000.;
          double fgMeritMn = fgMeritMx;
          int kTyMin = 30;
	  int nPhot = 0;
//------- look for pion, kaon or proton only :
//        ------------------------------------
          for( int kPaTy=kTyLL; kPaTy<=14; kPaTy+=3 ) {

            qSquare[kPaTy] = fgMeritMx;
	    probaQS[kPaTy] = -1.;

//	    Remember at thresh -> 0 phot <=> backgr ...
//@@----------------------------------------------------
//--------- check particle threshold :
//          --------------------------
	    double theIpo = (*ir)->pPartPhot()->thetaIpo( kPaTy );
	    if( theIpo >= 0. ) {

//----------- CHISQUARE :
//            -----------
//----------- compute qsquare for the given mass :
//            ------------------------------------
	      double qSqua = (*ir)->getQsquare( theIpo, nPhot );
//                           ----------------------------------
	      qSquare[kPaTy] = qSqua;
	      probaQS[kPaTy] = 0.;
	      if( qSqua > 0. ) probaQS[kPaTy] = prob_( qSqua, nPhot );


//----------- choose figure-of-merit type for PID :
//            -------------------------------------
//            kFigMe = 1: chi-square;  kFigMe = 2: likelihood :
//            -------------------------------------------------
              double fgMerit = qSqua;

//----------- define selected minimum (chisq) = maximum (like)
              if( fgMerit < fgMeritMn ) {
                fgMeritMn = fgMerit;
                kTyMin = kPaTy;
	      }

            }   // end if on cosTheW

          }   // end loop on particle type: kPaTy
	  if( kTyMin < 30 ) nPIDEv++;

//------- set variables to ring :
//        -----------------------
	  (*ir)->setnPhotQsQ( nPhoRing );
	  (*ir)->setQSquareRing( &qSquare[0] );
//@@------------------------------------------
//------- define partProbs :
//        ------------------
	  (*ir)->setPartProb( 9, double( nPhoRing ) );
	  (*ir)->setPartProb( 12,  qSquare[ 8] );
	  (*ir)->setPartProb( 13,  qSquare[11] );
	  (*ir)->setPartProb( 14,  qSquare[14] );
//@@      --------------------------------------
	  if( nProb > 15 ) {
 	    (*ir)->setPartProb( 19,  qSquare[ 2] );
	    (*ir)->setPartProb( 20,  qSquare[ 5] );
//@@        --------------------------------------
	  }
	  //cout << "  " << partProbs[12] << "  " << partProbs[13]
	  //     << "  " << partProbs[14] << endl;


//------- reconstructed ring chi-square :
//        -------------------------------
//        warning : what is theReco?
	  int nPhoQsq = 0;
	  double qSquaR = (*ir)->getRingQsq( theReco, nPhoQsq );
//                        -------------------------------------
          if( nPhoQsq > 0 ) {
	    xh = qSquaR*nPhoQsq;
            yh = nPhoQsq;
            if( hist.hRC1569 ) hist.hRC1569->Fill( xh, yh );
//hh                           ----------------------------
	  }
          xh = qSquaR;
          if( hist.hRC1568 ) hist.hRC1568->Fill( xh );
//hh                         ------------------------

	  (*ir)->setRingQsQ( qSquaR );
//@@---------------------------------
	  (*ir)->setPartProb( 11, qSquaR );
//@@--------------------------------------

	}

      }

  }


//===========================================================================
  void CsRCEventRings::partRingAllIdent() {
//-----------------------------------------


//--- Ring particle identification
//    ----------------------------
//--- Paolo  -  January 2003


      CsRCExeKeys* key = CsRCExeKeys::Instance();

      static bool firstCall = true;
      if( firstCall ) { 
        firstCall = false;

        key->acknoMethod( "CsRCEventRings::PartRingAllIdent" );
      }

      CsRCRecConst* cons = CsRCRecConst::Instance();
      int kOffHi = cons->kOffHi();
      float CFRefInd = cons->CFRefInd();
      double* massPartv = cons->massPartv();
      int nPhotMinRing = cons->nPhotMinRing();

      CsRCHistos& hist = CsRCHistos::Ref();
      double xh, yh, wh;


//--- ALL+RING PARTICLE IDENTIFICATION :
//    ----------------------------------
      int nPIDEv = 0;

//--- loop over rings :
//    -----------------
      list<CsRCRing*> lRings = lRings_;
      list<CsRCRing*>::iterator ir;
      for( ir=lRings.begin(); ir!=lRings.end(); ir++ ) {
        if( !(*ir)->flag() ) continue;

//----- define SELECTED rings :
//      -----------------------
        bool bSele = true;

//----- enough photons in the ring :
//      ---------------------------
	list<CsRCPhoton*> lPhotons = (*ir)->lPhotons();
        list<CsRCPhoton*>::iterator ih;
        int nPhoRing = lPhotons.size();
        if( nPhoRing < nPhotMinRing ) bSele = false;
//@@                                  -------------

//----- particle momentum over muon threshold (1.91 GeV/c) :
//      ----------------------------------------------------
	CsRCParticle* part = (*ir)->pPartPhot()->pPart();
        float momMeas = part->mom();
        //if( momMeas < 1.91 ) bSele = false;
//@@--------------------------------------

//----- process selected rings :
//      ------------------------
        if( bSele ) {

//------- reconstructed Cerenkov angle :
//        ------------------------------
	  double theReco = (*ir)->theReco();
	  double theRecoR = (*ir)->the();

	  double theRecoLkR = 0.;
	  if( (*ir)->pPartPhot()->thetaLikeSet() ) {
	    theRecoLkR = (*ir)->pPartPhot()->thetaLike();
            xh = theRecoR - theRecoLkR;
	    if( hist.hRC1633 ) hist.hRC1633-> Fill( xh );
//hh                           -------------------------
	  }

	  CsRCLikeAll* likeAll = NULL;
	  if( cons->backgrType() == "02" )
	    likeAll = new CsRCLikeAll02( (*ir)->pPartPhot() );
	  if( cons->backgrType() == "03" )
	    likeAll = new CsRCLikeAll03( (*ir)->pPartPhot() );
//@@      ---------------------------------------------------
	  if( likeAll ) {
	    (*ir)->setRingBack( likeAll->getRingBackground( theReco ) );
//@@        -----------------------------------------------------------
	  }
	  delete likeAll;

	}

      }

  }


//===========================================================================
  void CsRCEventRings::partRingIdent() {
//--------------------------------------


//--- Ring particle identification
//    ----------------------------
//--- Paolo  -  December 2002


      CsRCExeKeys* key = CsRCExeKeys::Instance();

      static bool firstCall = true;
      if( firstCall ) { 
        firstCall = false;

        key->acknoMethod( "CsRCEventRings::PartRingIdent" );

	if( key->likeONLY() ) {
	  std::cout << std::endl;
	  std::cout << "RICHONE, CsRCEventRings() : WARNING,"
	            << " tests on Likelihood ACTIVE!" << std::endl;
	  std::cout << "------------------------------------"
	            << "----------------------------" << std::endl;
	}
      }

      CsRCRecConst* cons = CsRCRecConst::Instance();
      int kOffHi = cons->kOffHi();
      float CFRefInd = cons->CFRefInd();
      double* massPartv = cons->massPartv();
      int nPhotMinRing = cons->nPhotMinRing();

      int nProb = cons->outBufferSize();

      CsRCHistos& hist = CsRCHistos::Ref();
      double xh, yh, wh;


//--- PARTICLE IDENTIFICATION :
//    -------------------------
      int nPIDEv = 0;

//--- loop over rings :
//    -----------------
      list<CsRCRing*> lRings = lRings_;
      list<CsRCRing*>::iterator ir;
      for( ir=lRings.begin(); ir!=lRings.end(); ir++ ) {
        if( !(*ir)->flag() ) continue;

//----- define SELECTED rings :
//      -----------------------
        bool bSele = true;

//----- enough photons in the ring :
//      ---------------------------
	list<CsRCPhoton*> lPhotons = (*ir)->lPhotons();
        list<CsRCPhoton*>::iterator ih;
        int nPhoRing = lPhotons.size();
        if( nPhoRing < nPhotMinRing ) bSele = false;
//@@                                  -------------

//----- particle momentum over muon threshold (1.91 GeV/c) :
//      ----------------------------------------------------
	CsRCParticle* part = (*ir)->pPartPhot()->pPart();
        float momMeas = part->mom();
        //if( momMeas < 1.91 ) bSele = false;
//@@--------------------------------------

//----- reconstructed Cerenkov angle :
//      ------------------------------
	double theReco = (*ir)->theReco();

	double theRecoR = (*ir)->the();
	double theRecoFt = (*ir)->thetaRFit();

	//xh = theRecoR - theRecoFt;
	//if( hist.hRC1631 ) hist.hRC1631-> Fill( xh );
//hh                       -------------------------

//----- process selected rings :
//      ------------------------
        if( bSele ) {

//------- compute maximum likelihood and its angle :
//        ------------------------------------------
	  double theRecoLkR = 0.;
	  long double likeMax = -1.;
	  double theRecoLw = 0.;
	  long double likeMw = -1.;
	  if(  key->doThetaLikeMax()  &&  !key->likeONLY() ) {
	    if( (*ir)->getThetaLikeMax( likeMw, theRecoLw ) ) {
//              -------------------------------------------
	      likeMax = likeMw;
	      theRecoLkR = theRecoLw;
	      (*ir)->setThetaLike( theRecoLkR );
//            ---------------------------------
	    }    
	  }

	  //(*ir)->setThetaLike( theRecoLkR );
//        //---------------------------------
	  (*ir)->setPartProb( 7, theRecoLkR );
//@@-----------------------------------------
	  //cout << "     " << likeMax << "  " << theRecoLkR << endl;

	  if( cons->thetaType() == "LIKE" ) theReco = theRecoLkR;
//@@                                        --------------------

	  if( theRecoLkR > 0. ) {
	    xh = theRecoR - theRecoLkR;
	    if( hist.hRC1633 ) hist.hRC1633-> Fill( xh );
//hh                           -------------------------
	  }

	  static bool checkLk = true;
	  if( !key->doCheckLike()  ||  key->likeONLY() ) checkLk = false;
	  if( checkLk ) (*ir)->checkLikelihood( theReco );
//                      ---------------------------------
	  static bool checkLRf = true;
	  if( !key->doCheckLike()  ||  key->likeONLY() ) checkLRf = false;
	  if( checkLRf ) (*ir)->checkLikeVsIndex( theReco );
//                       ----------------------------------

	  int kTyLL = 2;
//@@-------------------
          int kFigMe = 0;
          double probaLK[31];
          double derivLK[31];
          double deriv2LK[31];
//------- look for all(5) particles :
//        ---------------------------
          for( int kPaTy=kTyLL; kPaTy<=14; kPaTy+=3 ) {
//        -------------------------------------------

            float massIpo = part->mass( kPaTy );
            double betaIpo = momMeas /
              sqrt( momMeas*momMeas + massIpo*massIpo );
//	    Remember at thresh -> 0 phot <=> backgr ...
//@@---------------------------------------------------

	    probaLK[kPaTy] = -1.;
	    derivLK[kPaTy] =  0.;
	    deriv2LK[kPaTy] =  0.;

//--------- check particle threshold :
//          --------------------------
	    double theIpo = (*ir)->pPartPhot()->thetaIpo( kPaTy );
	    if( theIpo >= 0. ) {

//----------- RING LIKELIHOOD :
//            -----------------
//----------- compute likelihood for the given mass :
//            ---------------------------------------
//----------- compute Likelihood using measured angle and
//            theIpo computed from 'right' ref index

	      //theIpo = double( - kPaTy );

	      int nPhot = 0;
	      long double pLike = -1.;
	      pLike = (*ir)->getLikelihood( theReco, theIpo, nPhot );
//            ------------------------------------------------------
	      //pLike /= likeMax;
//@@----------------------------
	      //if( pLike > 1. ) std::cout<< pLike << std::endl;
//                             ------------------------------
	      if( pLike >= 0. ) probaLK[kPaTy] = pLike;
//@@                            ----------------------

//----------- compute likelihood derivatives :
//            --------------------------------
//----------- compute Likelihood using measured angle and
//            theIpo normalized to ref. index
	      if( pLike >= 0. ) {
		//double DCFRefInd = 0.000025;
		double DCFRefInd = 0.000050;
//@@---------------------------------------
		//double cosTheWD = 1./ ( betaIpo * (CFRefInd+DCFRefInd) );
		//double theIpoD = acos( cosTheWD ) * 1000.;
		//long double pLikeD = -1.;
	        //pLikeD = (*ir)->getLikelihood( theReco, theIpoD, nPhot );
//              //--------------------------------------------------------
		//double likeDerX = (pLikeD - pLike) / DCFRefInd;

		double cosTheWP = 1./ ( betaIpo * (CFRefInd + DCFRefInd) );
		if( cosTheWP > 1.) cosTheWP = 1.;
		double theIpoP = acos( cosTheWP ) * 1000.;
		//long double pLikeP = -1.;
		long double pLikeP = 0.;
		if( !key->likeONLY() ) {
	          pLikeP = (*ir)->getLikelihood( theReco, theIpoP, nPhot );
//                --------------------------------------------------------
		}
		double cosTheWM = 1./ ( betaIpo * (CFRefInd - DCFRefInd) );
		if( cosTheWM > 1.) cosTheWM = 1.;
		double theIpoM = acos( cosTheWM ) * 1000.;
		//long double pLikeM = -1.;
		long double pLikeM = 0.;
		if( !key->likeONLY() ) {
	          pLikeM = (*ir)->getLikelihood( theReco, theIpoM, nPhot );
//                --------------------------------------------------------
		}
		double likeDer = 0.5 * (pLikeP - pLikeM) / DCFRefInd;
		derivLK[kPaTy] = likeDer;

		double likeDer2 = (pLikeP - 2.*pLike + pLikeM) / 
		  (DCFRefInd*DCFRefInd);
		deriv2LK[kPaTy] = likeDer2;
		//cout << kPaTy << "  " << pLikeP << "  " << pLike << "  "
		//     << pLikeM << endl;
		//cout << "  " << likeDer << "  " << likeDer << "  "
		//     << likeDer2 << endl;

	      }

            }   // end if on cosTheW

          }   // end loop on particle type: kPaTy

//------- set variables to ring :
//        -----------------------
	  (*ir)->setProbaLKRing( &probaLK[0] );
	  (*ir)->setDerivLKRing( &derivLK[0] );
//@@------------------------------------------

//------- define partProbs :
//        ------------------
	  (*ir)->setPartProb( 1,  probaLK[ 8] );
	  (*ir)->setPartProb( 2,  probaLK[11] );
	  (*ir)->setPartProb( 3,  probaLK[14] );
	  (*ir)->setPartProb( 4,  derivLK[ 8] );
	  (*ir)->setPartProb( 5,  derivLK[11] );
	  (*ir)->setPartProb( 6,  derivLK[14] );
//@@      -------------------------------------
	  if( nProb > 15 ) {
    	    (*ir)->setPartProb( 15, probaLK[ 2] );
	    (*ir)->setPartProb( 16, probaLK[ 5] );
	    (*ir)->setPartProb( 17, derivLK[ 2] );
	    (*ir)->setPartProb( 18, derivLK[ 5] );
//@@        -------------------------------------
	  }
//@@      prov. 030411 ------------------------
	  //(*ir)->setPartProb( 12,  deriv2LK[ 8] );   //   !!! WARNING
	  //(*ir)->setPartProb( 13,  deriv2LK[11] );   //   !!! WARNING
	  //(*ir)->setPartProb( 14,  deriv2LK[14] );   //   !!! WARNING
//@@      -------------------------------------

//------- compute likelihood of 'background' :
//        ------------------------------------
	  int nPhotBk = 0;
	  double theIpoBk = 0.;
	  long double pLikeBk = -1.;
       	  pLikeBk = (*ir)->getLikelihood( theReco, theIpoBk, nPhotBk );
//        ------------------------------------------------------------
          //cout << "     " << pLikeBk << endl;
	  //pLikeBk /= likeMax;
//@@--------------------------
	  //if( pLikeBk > 1. ) std::cout<< pLikeBk << std::endl;
//                           --------------------------------
	  if( pLikeBk > 1. ) pLikeBk = 1.;

	  (*ir)->setProbaLKBgRing( pLikeBk );
//@@----------------------------------------
	  (*ir)->setPartProb( 0, pLikeBk );
//@@--------------------------------------

//------- compute ring background :
//        -------------------------
	  string backType = cons->backgrType();
	  CsRCLikeRing* likeRing = NULL;
	  if( backType == "02" ) likeRing = new CsRCLikeRing02( (*ir) );
	  if( backType == "03" ) likeRing = new CsRCLikeRing03( (*ir) );
	  if( backType == "04" ) likeRing = new CsRCLikeRing04( (*ir) );
//@@                             --------------------------------------
	  if( likeRing ) {
	    (*ir)->setRingBack( likeRing->getRingBackground( theReco ) );
//@@        ------------------------------------------------------------
	  }
	  delete likeRing;


	  (*ir)->setPartProbsSet();
//@@------------------------------

	}

      }

  }


//===========================================================================
  void CsRCEventRings::checkPhotonAngle() {
//-----------------------------------------


//- Paolo  - Dcember  2005


    if( !CsRichOne::Instance()->UpRICHJob() ) return;
//  ------------------------------------------------

    if( !CsRCHistos::Ref().bookHis()  ||
        CsRCHistos::Ref().levelBk() < 2 ) return;
//  --------------------------------------------

    CsRCExeKeys* key = CsRCExeKeys::Instance();
    CsRCRecConst *cons = CsRCRecConst::Instance();
    static int kOffHi = cons->kOffHi();
    static double radDeg = cons->RadDeg();
    CsRCMirrors *mirr = CsRCMirrors::Instance();
    CsRCDetectors *dets = CsRCDetectors::Instance();
    CsRCHistos& hist = CsRCHistos::Ref();
    int kh;
    double xh, yh, wh;

    static std::vector<CsHist2D*> hRCPaIn;
    static std::vector<CsHist2D*> hRCPhoAd;
    static std::vector<CsHist2D*> hRCPhoPa;
    static std::vector<CsHist2D*> hRCPhoAc;
    static std::vector<CsHist2D*> hRCPaAc;
    
    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      key->acknoMethod( "CsRCEventRings::checkPhotonAngle" );

      std::string hTitle = "tgy vs tgx";
      int kName = 0;
      for( int kH=0; kH<4; kH++ ) {
	hTitle = "PaThx vs PaX";
	kName = kOffHi + 3400 + kH;
	stringstream hNPaIn;
	hNPaIn << kName;
	hRCPaIn.push_back( new CsHist2D( hNPaIn.str(), hTitle,
		  	    100, -1000., 1000., 100, -250., 250. ) );
	hTitle = "thy vs thx";
	kName = kOffHi + 3405 + kH;
	stringstream hNPhoAd;
	hNPhoAd << kName;
	hRCPhoAd.push_back( new CsHist2D( hNPhoAd.str(), hTitle,
		  	    100, -300., 300., 100, -400., 400. ) );
	hTitle = "thxPho vs thxPa";
	kName = kOffHi + 3430 + kH;
	stringstream hNPhoPa;
	hNPhoPa << kName;
	hRCPhoPa.push_back( new CsHist2D( hNPhoPa.str(), hTitle,
		  	    100, -400., 400., 100, -400., 400. ) );
      }
      for( int kH=0; kH<6; kH++ ) {
	hTitle = "thy vs thx";
	kName = kOffHi + 3410 + kH;
	stringstream hNPhoAc;
	hNPhoAc << kName;
	hRCPhoAc.push_back( new CsHist2D( hNPhoAc.str(), hTitle,
		  	    100, -300., 300., 100, -400., 400. ) );
	hTitle = "Pay vs Pax";
	kName = kOffHi + 3420 + kH;
	stringstream hNPaAc;
	hNPaAc << kName;
	hRCPaAc.push_back( new CsHist2D( hNPaAc.str(), hTitle,
		  	   100, 0., 600., 100, 200., 800. ) );
      }
    }


//- all rings :
//  -----------
    list<CsRCRing*>::iterator ir;
    for( ir=lRings_.begin(); ir!=lRings_.end(); ir++ ) {
      if( !(*ir)->flag() ) continue;

      CsRCPartPhotons* papho = (*ir)->pPartPhot();
      CsRCParticle* part = papho->pPart();
      double xPa = part->vPosIn().x();
      double tgxPa = part->vDirIn().x()/part->vDirIn().z();;
      double xPade = part->vPade()[part->kDetPart()].x();
      double yPade = part->vPade()[part->kDetPart()].y();
      double tgxPaRW = papho->vDcPaReflW()[part->kDetPart()].x()/
	papho->vDcPaReflW()[part->kDetPart()].z();

      int kCat = papho->iCaPa();
      int kHist = -1;
      for( int kc=0; kc<int(dets->nCatPMT().size()); kc++ ) {
	if( kCat == dets->nCatPMT()[kc] ) kHist = kc;
      }
      //^if( kCat ==  3 ) kHist = 0;
      //^if( kCat ==  5 ) kHist = 1;
      //^if( kCat == 10 ) kHist = 2;
      //^if( kCat == 12 ) kHist = 3;
      //^if( kHist >= 0 ) {
      if( kHist >= 0  &&  kHist < int(dets->nCatPMT().size()) ) {
	xh = xPa;
	yh = atan( tgxPa ) * 1000.;
	if( hRCPaIn[kHist] ) hRCPaIn[kHist]->Fill( xh, yh );
//hh                         ------------------------------
      }

      list<CsRCPhoton*> lPhotons = (*ir)->lPhotons();
      list<CsRCPhoton*>::iterator ih;
      if( lPhotons.size() < 5 ) continue;
//@@------------------------------------
      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {

	int kDetClu = (*ih)->kDetClu();
	float RR = mirr->RRv( kDetClu );
	Hep3Vector vPoPhotW = papho->vPoPhotW()[kDetClu];
        double al = sin((*ih)->the()/1000.) * cos((*ih)->phiA()*radDeg);
        double am = sin((*ih)->the()/1000.) * sin((*ih)->phiA()*radDeg);
        double an = cos((*ih)->the()/1000.);
        Hep3Vector vDcPhoEmW( al, am, an );
        vDcPhoEmW = vDcPhoEmW.unit();
//----- assume photon emitted from average point on part. traj.
//      -------------------------------------------------------
//----- 'photon' impact on mirror (MWR) :
//      ---------------------------------
        Hep3Vector vPoC( 0., 0., 0. );
        Hep3Vector vPoPhoMir =
          mirr->vImpMir( vPoPhotW, vDcPhoEmW, vPoC, RR );
//        ----------------------------------------------
//----- normal to mirror at 'photon' impact :
//      -------------------------------------
        Hep3Vector vDcNoPhoMir = (1./RR) * vPoPhoMir;
//----- 'photon' reflected direction :
//      ------------------------------
        double cosPhoMir = vDcNoPhoMir * vDcPhoEmW;
        Hep3Vector vDcPhoRefl = 2.*cosPhoMir * vDcNoPhoMir - vDcPhoEmW;
//----- 'photon' incidence angles on detector :
//      ---------------------------------------
        double tgxPho = vDcPhoRefl.x() / vDcPhoRefl.z();
        double tgyPho = vDcPhoRefl.y() / vDcPhoRefl.z();
        //std::cout << tgxPho << "  " << tgyPho << std::endl;

        int kCat = (*ih)->ptToClu()->ic();
        int kHist = -1;
	for( int kc=0; kc<int(dets->nCatPMT().size()); kc++ ) {
	  if( kCat == dets->nCatPMT()[kc] ) kHist = kc;
	}
        //^if( kCat ==  3 ) kHist = 0;
        //^if( kCat ==  5 ) kHist = 1;
        //^if( kCat == 10 ) kHist = 2;
        //^if( kCat == 12 ) kHist = 3;
	xh = 0.;
	yh = 0.;
        //^if( kHist >= 0 ) {
	if( kHist >= 0  &&  kHist < int(dets->nCatPMT().size()) ) {
	  xh = atan( tgxPho ) * 1000.;
	  yh = atan( tgyPho ) * 1000.;
          if( hRCPhoAd[kHist] ) hRCPhoAd[kHist]->Fill( xh, yh );
//hh                            -------------------------------
	  xh = atan( tgxPa ) * 1000.;
	  //xh = atan( tgxPaRW ) * 1000.;
	  yh = atan( tgxPho ) * 1000.;
          if( hRCPhoPa[kHist] ) hRCPhoPa[kHist]->Fill( xh, yh );
//hh                            -------------------------------
        }
        int kZone = -1;
        //^if(  kCat ==  3 ) {
	if(  kCat == dets->nCatPMT()[0] ) {
          float ddYY = 204;
          //double xc = (*ih)->ptToClu()->xc();
          //double yc = (*ih)->ptToClu()->yc();
          double xc = xPade;
          double yc = yPade;
          if( yc < 290.+ddYY ) {
            if( xc < 190. ) kZone = 0;
            if( xc > 190.  &&  xc <= 380.) kZone = 1;
            if( xc > 380.) kZone = 2;
          }
          if( yc > 290.+ddYY ) {
            if( xc < 190. ) kZone = 3;
            if( xc > 190.  &&  xc <= 380.) kZone = 4;
            if( xc > 380.) kZone = 5;
          }
          if( kZone >= 0 ) {
	    xh = atan( tgxPho ) * 1000.;
	    yh = atan( tgyPho ) * 1000.;
            if( hRCPhoAc[kZone] ) hRCPhoAc[kZone]->Fill( xh, yh );
//hh                              -------------------------------
	    if( ih == lPhotons.begin() ) {
	      xh = xPade;
	      yh = yPade;
	      if( hRCPaAc[kZone] ) hRCPaAc[kZone]->Fill( xh, yh );
//hh                               ------------------------------
	    }
	  }
	}
      }

    }

  }


//===========================================================================
  void CsRCEventRings::singlePhotonPMT() {
//----------------------------------------

//--- Paolo  -  April 2006


      if( !CsRichOne::Instance()->UpRICHJob() ) return;
//    ------------------------------------------------

      CsRCExeKeys* key = CsRCExeKeys::Instance();
      CsRCRecConst *cons = CsRCRecConst::Instance();
      int kOffHi = cons->kOffHi();
      double RadDeg = cons->RadDeg();
      float CFRefInd = cons->CFRefInd();
      CsRCDetectors *dets = CsRCDetectors::Instance();

      CsRCHistos& hist = CsRCHistos::Ref();
      int kh;
      double xh, yh, wh;

      static bool firstCall = true;
      if( firstCall ) {
	firstCall = false;

	key->acknoMethod( "CsRCEventRings::singlePhotonPMT" );
      }

//@@-----------------------------------------------------------
      static Hep3Vector
	vNormToLens3W( tan( 5.5/RadDeg), tan(-15./RadDeg), 1.);
      vNormToLens3W = vNormToLens3W.unit();
      static Hep3Vector
	vNormToLens5W( tan(-5.5/RadDeg), tan(-15./RadDeg), 1.);
      vNormToLens5W = vNormToLens5W.unit();
      static Hep3Vector
	vNormToLens10W( tan( 5.5/RadDeg), tan( 15./RadDeg), 1.);
      vNormToLens10W = vNormToLens10W.unit();
      static Hep3Vector
	vNormToLens12W( tan(-5.5/RadDeg), tan( 15./RadDeg), 1.);
      vNormToLens12W = vNormToLens12W.unit();
      //std::cout << vNormToLens3W << "  " << vNormToLens5W << "  "
      //        << vNormToLens10W << "  " << vNormToLens12W << std::endl;
//@@-----------------------------------------------------------

//--- loop over rings :
//    -----------------
      CsRCEventRings* evrings = CsRCEventRings::Instance();
      list<CsRCRing*> lRings = evrings->lRings();
      list<CsRCRing*>::iterator ir;
      for( ir=lRings.begin(); ir!=lRings.end(); ir++ ) {

        CsRCRing* ring = (*ir);
        if( !ring ) continue;
        CsRCPartPhotons* papho = ring->pPartPhot();
        if( !papho ) continue;
	if( !(*ir)->flag() ) continue;

        int kDetPart = papho->kDetPart();
	int iCaPa = papho->iCaPa();
	Hep3Vector vDcPaReflW = papho->vDcPaReflW()[kDetPart];
	Hep3Vector vNormToLens;
	//^if( iCaPa == 3 ) vNormToLens = vNormToLens3W;
	//^else if( iCaPa == 5 ) vNormToLens = vNormToLens5W;
	//^else if( iCaPa == 10 ) vNormToLens = vNormToLens10W;
	//^else if( iCaPa == 12 ) vNormToLens = vNormToLens12W;
	if( iCaPa == dets->nCatPMT()[0] ) vNormToLens = vNormToLens3W;
	else if( iCaPa == dets->nCatPMT()[1] ) vNormToLens = vNormToLens5W;
	else if( iCaPa == dets->nCatPMT()[2] ) vNormToLens = vNormToLens10W;
	else if( iCaPa == dets->nCatPMT()[3] ) vNormToLens = vNormToLens12W;
	else  continue;
	//std::cout << vNormToLens3W << "  " << vDcPaReflW << std::endl;
	double thePaLens = acos( vDcPaReflW.dot( vNormToLens ) );
	//std::cout << thePaLens << std::endl;
        thePaLens *= RadDeg;   //   to deg

	xh = thePaLens;
	if( hist.hRC3640 ) hist.hRC3640->Fill( xh );
//hh                       ------------------------
	kh = -1;
	if( thePaLens >  0.  &&  thePaLens <=  3. ) kh = 0;
	if( thePaLens >  3.  &&  thePaLens <=  6. ) kh = 1;
	if( thePaLens >  6.  &&  thePaLens <=  9. ) kh = 2;
	if( thePaLens >  9.  &&  thePaLens <= 12. ) kh = 3;
	if( kh >= 0  &&  kh <= 3 ) {
          double theReco = ring->the();
          list<CsRCPhoton*> lPhotons = papho->lPhotons();
	  list<CsRCPhoton*>::iterator ih;
	  for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
  	    double thePho = (*ih)->the();
	    double phiPho = (*ih)->phi();
	    if( !(*ih)->isPMT() ) continue;

  	    xh = thePho - theReco;
	    yh = phiPho;
	    if( hist.vRC3640[kh] ) hist.vRC3640[kh]->Fill( xh, yh );
//hh                               --------------------------------
	  }
	}

      }

      return;
  }


//===========================================================================
  void CsRCEventRings::checkPMTTimeDif() {
//----------------------------------------

//- Paolo  -  March 2007


    if( !CsRichOne::Instance()->UpRICHJob() ) return;
//  ------------------------------------------------

    CsRCExeKeys* key = CsRCExeKeys::Instance();
    CsRCRecConst *cons = CsRCRecConst::Instance();
    int kOffHi = cons->kOffHi();

    CsRCHistos& hist = CsRCHistos::Ref();
    int kh;
    double xh, yh, wh;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      key->acknoMethod( "CsRCEventRings::checkPMTTimeDif" );
    }

    if( lRings_.size() > 0 ) {
      bool set = false;
      double time0 = 0.;
      list<CsRCRing*>::iterator ir;
      for( ir=lRings_.begin(); ir!=lRings_.end(); ir++ ) {
	if( !(*ir)->flag() ) continue;
	if( !lRings_.front()->mTimeSet() ) continue;
	if( !set ) {
	  set = true;
	  time0 = (*ir)->mTime();
	  continue;
	}
	if( set ) {
          double dTime = (*ir)->mTime() - time0;
          //std::cout << dTime << std::endl;
          xh = dTime;
          if( hist.hRC1536 ) hist.hRC1536->Fill( xh );
//hh                         ------------------------
        }
      }
    }

    return;
  }


//===========================================================================
  void CsRCEventRings::checkEventSignal() {
//-----------------------------------------

//- Paolo  -  March 2007

//  WARNING : DO NOT FORGET the 70. mrad fiducial radius!


    if( !CsRichOne::Instance()->UpRICHJob() ) return;
//  ------------------------------------------------

    CsRCExeKeys* key = CsRCExeKeys::Instance();
    CsRCRecConst *cons = CsRCRecConst::Instance();
    int kOffHi = cons->kOffHi();

    CsRCHistos& hist = CsRCHistos::Ref();
    int kh;
    double xh, yh, wh;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      key->acknoMethod( "CsRCEventRings::checkEventSignal" );
    }
    //std::cout << std::endl;
    //std::cout << lRings_.size() << std::endl;
    if( lRings_.size() > 0 ) {
      list<CsRCRing*>::iterator ir;
      for( ir=lRings_.begin(); ir!=lRings_.end(); ir++ ) {
	if( !(*ir)->flag() ) continue;

	double nPhotRg = (*ir)->lPhotons().size();
	double nPhotExRg = (*ir)->pPartPhot()->normSg();
	double nPhotBk = 0;
	double nPhotExBk = 0;
	list<CsRCRing*>::iterator jr;
	for( jr=lRings_.begin(); jr!=lRings_.end(); jr++ ) {
	  if( !(*jr)->flag() ) continue;
	  if( (*jr) == (*ir) ) continue;
	  nPhotBk += (*jr)->lPhotons().size();
	  nPhotExBk += (*jr)->pPartPhot()->normSg();
	}
	//std::cout << (*ir)->kRing() << "  " << nPhotRg << "  " << nPhotBk
	//          << "  " << nPhotExRg << "  " << nPhotExBk << std::endl;

      }
    }

    return;
  }


//===========================================================================
  void CsRCEventRings::testPhiKzero() {
//-------------------------------------

//- Paolo  -  August 2007

//- WARNING : NOT working!(?)

    CsRichOne* richone = CsRichOne::Instance();
std::cout << "testPhiKzero  " << std::endl;
std::cout << richone->evGood() << " 1 " << richone->kMuonPart()
          << "  " << richone->zVertex()
          << "  " << richone->phiMass() << std::endl;
    int kPartCode = richone->kMuonPart();
    int kMuPart = 0;
    int kK1Part = 0;
    int kK2Part = 0;
    kK2Part = kPartCode / 10000;
    kK1Part = (kPartCode - kK2Part*10000) / 100;
    kMuPart = kPartCode - kK2Part*10000 - kK1Part*100;

    CsRCEventParticles* evparts = CsRCEventParticles::Instance();
    list<CsRCParticle*> lParticles = evparts->lParticles();
    list<CsRCParticle*>::iterator ipk;
    CsRCParticle* pMuPart = NULL;
    CsRCParticle* pK1Part = NULL;
    CsRCParticle* pK2Part = NULL;
    int kPart = 1;
    for( ipk=lParticles.begin(); ipk!=lParticles.end(); ipk++ ) {
      if( kPart == kMuPart ) pMuPart = (*ipk);
      if( kPart == kK1Part ) pK1Part = (*ipk);
      if( kPart == kK2Part ) pK2Part = (*ipk);
      kPart++;
    }
//std::cout << pMuPart << "  " << pK1Part << "  "
//    	<< pK2Part << std::endl;

//- particle at PT centre :
//- WARNING: needs TRAFFIC initialization! (Extrapolate...)
    static double cov[15];
    for( int k=0; k<15; k++ ) cov[k] = 0.;
    int j = 1;
    for( int k=0; k<15; k+=j ) { cov[k] = 1.; j++; }

    double zVertex = richone->zVertex();
    double dXdZ = 0.;
    double dYdZ = 0.;
    CsRCParticle* part = NULL;
    CsHelix hxInW1, hxVtx1;
    double mom1 = 0.;
    double momVx1 = 0.;
    Hep3Vector vDirVx1( 0. );
    part = pK1Part;
    if( part ) {
      hxInW1 = CsHelix( part->vPosIn().x(), part->vPosIn().y(),
    	  		part->vPosIn().z(),
		        part->vDirIn().x()/part->vDirIn().z(),
	       	        part->vDirIn().y()/part->vDirIn().z(),
	       	        part->charge()/part->mom(), cov );
      if( hxInW1.Extrapolate( zVertex, hxVtx1 ) ) {
        momVx1 = fabs(1./hxVtx1.getCop());
	dXdZ = hxVtx1.getDXDZ();
	dYdZ = hxVtx1.getDYDZ();
	vDirVx1.setX( dXdZ );
	vDirVx1.setY( dYdZ );
	vDirVx1.setZ( 1. );
	vDirVx1 = vDirVx1.unit();
//std::cout << part->charge() << " 1 ";
//std::cout << " 1 " << hxVtx1.getX() << "  " << hxVtx1.getY()
//  	    << "  " << hxVtx1.getZ() << std::endl;
      }
    mom1 = part->mom();
    }
    CsHelix hxInW2, hxVtx2;
    double mom2 = 0.;
    double momVx2 = 0.;
    Hep3Vector vDirVx2( 0. );
    part = pK2Part;
    if( part ) {
      hxInW2 = CsHelix( part->vPosIn().x(), part->vPosIn().y(),
			part->vPosIn().z(),
		        part->vDirIn().x()/part->vDirIn().z(),
	      	        part->vDirIn().y()/part->vDirIn().z(),
	       	        part->charge()/part->mom(), cov );
      if( hxInW2.Extrapolate( zVertex, hxVtx2 ) ) {
        momVx2 = fabs(1./hxVtx2.getCop());
	dXdZ = hxVtx2.getDXDZ();
	dYdZ = hxVtx2.getDYDZ();
	vDirVx2.setX( dXdZ );
	vDirVx2.setY( dYdZ );
	vDirVx2.setZ( 1. );
	vDirVx2 = vDirVx2.unit();
//std::cout << part->charge() << " 2 " << std::endl;
      }
      mom2 = part->mom();
    }
if( pK1Part  &&  pK2Part ) {
  //std::cout << std::endl;
std::cout << 1./hxInW1.getCop() << "  " << 1./hxVtx1.getCop() << "  " 
          << 1./hxInW2.getCop() << "  " << 1./hxVtx2.getCop()
          << " -    " << zVertex << std::endl;
std::cout << setprecision( 5 );
std::cout << pK1Part->vDirIn() << "  " << vDirVx1 << std::endl;
std::cout << setprecision( 2 );
std::cout << setprecision( 5 );
std::cout << pK2Part->vDirIn() << "  " << vDirVx2 << std::endl;
std::cout << setprecision( 2 );
}

    if( pK1Part  &&  pK2Part ) {
      Hep3Vector v1( momVx1*vDirVx1 );
      Hep3Vector v2( momVx2*vDirVx2 );
      Hep3Vector v0( v1 + v2 );
      double KMass = 0.49368;
      double E1 = sqrt( KMass*KMass + v1.mag2() );
      double E2 = sqrt( KMass*KMass + v2.mag2() );
      double E0 = E1 + E2;
      double PhiMass = sqrt( E0*E0 - v0.mag2() );
std::cout << setprecision( 5 );
std::cout << v1 << "  " << v2 << "  " << v0 << std::endl;
std::cout << E1 << "  " << E2 << "  " << E0 << std::endl;
std::cout << richone->phiMass() << "  " << PhiMass << std::endl;
std::cout << setprecision( 2 );
    }

    return;
  }



//===========================================================================
  void CsRCEventRings::checkRichMom() {
//-------------------------------------

//- Paolo  -  November  2008


    if( !CsRCHistos::Ref().bookHis()  ||
        CsRCHistos::Ref().levelBk() < 2 ) return;
//  --------------------------------------------

    CsRCExeKeys* key = CsRCExeKeys::Instance();
    CsRCRecConst* cons = CsRCRecConst::Instance();

    CsRCHistos& hist = CsRCHistos::Ref();
    int kOffHi = cons->kOffHi();
    static std::vector<CsHist2D*> vRC6800;
    double xh, yh, wh;

    bool print = false;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      for( int kh=0; kh<12; kh++ ) vRC6800.push_back( NULL );
      if( CsRCHistos::Ref().bookHis() &&
	  CsRCHistos::Ref().levelBk() >= 2 ) {
	CsHistograms::SetCurrentPath("/RICH");
	vRC6800.clear();
	for( int kh=0; kh<12; kh++ ) {
	  string hTitle = " ";
	  stringstream hN6800;
	  int kHist = kOffHi + 6800 + kh + 1;
	  hN6800 << kHist;
	  vRC6800.push_back( new CsHist2D( hN6800.str(), hTitle,
					   100, -0.25, 0.25, 100, 0., 50. ) );
	}
	CsHistograms::SetCurrentPath("/");
      }
    }

    double likeDV = cons->likeDefVa();

//- loop over part-photons :
//  ------------------------
    bool isFirst = true;
    CsRCEventPartPhotons* paphos = CsRCEventPartPhotons::Instance();
    list<CsRCPartPhotons*> lPartPhotons = paphos->lPartPhotons();
    list<CsRCPartPhotons*>::iterator ipf;
    for( ipf=lPartPhotons.begin(); ipf!=lPartPhotons.end(); ipf++ ) {
      if( !(*ipf)->flag() ) continue;
      CsRCPartPhotons* papho = (*ipf);
      CsRCRing* ring = papho->pRing();
      if( !ring ) continue;

      float momPart = papho->pPart()->mom();

      double* probaLK;
      probaLK = papho->probaLKAll();
//--- Production Likelihoods
      double likeBkg = papho->probaLKBgAll();
      double likeElec = probaLK[ 2];
      double likeMuon = probaLK[ 5];
      double likePion = probaLK[ 8];
      double likeKaon = probaLK[11];
      double likeProton = probaLK[14];
      bool bIDpart[4];
      if( likeElec == likeBkg ) continue;      // NO photons from signal

//--- Select ID pions (analysis way)
      bool bIDpion = false;
      if( likePion > likeDV  &&  likePion > likeBkg  &&
          likePion > likeKaon  &&  likePion > likeProton )
	bIDpion = true;
      if( momPart < 8.  &&  likePion < likeElec ) bIDpion = false;
      bIDpart[0] = bIDpion;
//--- Select ID kaons (analysis way)
      bool bIDkaon = false;
      if( likeKaon > likeDV  &&  likeKaon > likeBkg  &&
          likeKaon > likePion  &&  likeKaon > likeProton )
	bIDkaon = true;
      if( momPart < 8.  &&  likeKaon < likeElec ) bIDkaon = false;
      bIDpart[1] = bIDkaon;
//--- Select ID protons (analysis way) ?
      bool bIDproton = false;
      if( likeProton > likeDV  &&  likeProton > likeBkg  &&
          likeProton > likePion  &&  likeProton > likeKaon )
	bIDproton = true;
      if( momPart < 8.  &&  likeProton < likeElec ) bIDproton = false;
      bIDpart[2] = bIDproton;
//--- Select ID muons (analysis way) ???
      bool bIDmuon = false;
      if( likeMuon > likeDV  &&  likeMuon > likeBkg  &&
          likeMuon > likePion  &&
          likeMuon > likeKaon  &&  likeMuon > likeProton )
	bIDmuon = true;
      if( momPart < 8.  &&  likeMuon < likeElec ) bIDmuon = false;
      bIDpart[3] = bIDmuon;

      double theReco = ring->theReco();
      double theRecoFt = ring->thetaRFit();
      double theLike = 0.;
      if( papho->thetaLikeSet() ) theLike = papho->thetaLike();
      double theRing = 0.;

      double CFRefInd = cons->CFRefInd();
      double gammaPart = 0;
      double betaPart = 0.;
      double massPart = 0.;
      double momReco = 0.;

      int kPart;
//--- ID particles
      theRing = theReco;
      if( theRing > 0. ) {
        betaPart = 1./( cos( theRing/1000. ) * CFRefInd );
        if( betaPart < 1.) {
	  gammaPart = 1./ sqrt(1.- betaPart*betaPart);
        } else  { gammaPart = 0.; }
        kPart = 0;
        for( int kPaTy=5; kPaTy<=14; kPaTy+=3 ) {
	  massPart = papho->pPart()->mass( kPaTy );
	  momReco = gammaPart*betaPart * massPart;
	  xh = (momReco - momPart) / momPart;
	  yh = momPart;
	  //int kHH = kPart;
	  int jPart = kPart - 1;
	  if(  jPart < 0 ) jPart = 3;
	  int kHH = jPart;
	  if( bIDpart[jPart] ) {
            if( kHH >= 0  &&  kHH < int( vRC6800.size() ) ) {
	      if( vRC6800[kHH] ) vRC6800[kHH]->Fill( xh, yh );
//hh                             ----------------------------
	    }
	  }
	  kPart++;
	}
      }
      theRing = theRecoFt;
      if( theRing > 0. ) {
        betaPart = 1./( cos( theRing/1000. ) * CFRefInd );
        if( betaPart < 1.) {
	  gammaPart = 1./ sqrt(1.- betaPart*betaPart);
        } else  { gammaPart = 0.; }
        kPart = 0;
        for( int kPaTy=5; kPaTy<=14; kPaTy+=3 ) {
	  massPart = papho->pPart()->mass( kPaTy );
	  momReco = gammaPart*betaPart * massPart;
	  xh = (momReco - momPart) / momPart;
	  yh = momPart;
	  //int kHH = kPart + 4;
	  int jPart = kPart - 1;
	  if(  jPart < 0 ) jPart = 3;
	  int kHH = jPart + 4;
	  if( bIDpart[jPart] ) {
            if( kHH >= 0  &&  kHH < int( vRC6800.size() ) ) {
	      if( vRC6800[kHH] ) vRC6800[kHH]->Fill( xh, yh );
//hh                             ----------------------------
	    }
	  }
	  kPart++;
	}
      }
      theRing = theLike;
      if( theRing > 0. ) {
        betaPart = 1./( cos( theRing/1000. ) * CFRefInd );
        if( betaPart < 1.) {
	  gammaPart = 1./ sqrt(1.- betaPart*betaPart);
        } else  { gammaPart = 0.; }
        kPart = 0;
        for( int kPaTy=5; kPaTy<=14; kPaTy+=3 ) {
	  massPart = papho->pPart()->mass( kPaTy );
	  momReco = gammaPart*betaPart * massPart;
	  xh = (momReco - momPart) / momPart;
	  yh = momPart;
	  //int kHH = kPart + 8;
	  int jPart = kPart - 1;
	  if(  jPart < 0 ) jPart = 3;
	  int kHH = jPart + 8;
	  if( bIDpart[jPart] ) {
            if( kHH >= 0  &&  kHH < int( vRC6800.size() ) ) {
	      if( vRC6800[kHH] ) vRC6800[kHH]->Fill( xh, yh );
//hh                             ----------------------------
	    }
	  }
	  kPart++;
	}
      }

    }

    return;
  }


//===========================================================================
  void CsRCEventRings::checkPhiTheta() {
//--------------------------------------

//- Paolo  -  November  2008


    if( !CsRCHistos::Ref().bookHis()  ||
        CsRCHistos::Ref().levelBk() < 2 ) return;
//  --------------------------------------------

    CsRCExeKeys* key = CsRCExeKeys::Instance();
    CsRCRecConst* cons = CsRCRecConst::Instance();

    CsRCHistos& hist = CsRCHistos::Ref();
    int kOffHi = cons->kOffHi();
    static std::vector<CsHist2D*> vRC6850;
    double xh, yh, wh;

    bool print = false;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      for( int kh=0; kh<12; kh++ ) vRC6850.push_back( NULL );
      if( CsRCHistos::Ref().bookHis() &&
	  CsRCHistos::Ref().levelBk() >= 2 ) {
	CsHistograms::SetCurrentPath("/RICH");
	vRC6850.clear();
	for( int kh=0; kh<12; kh++ ) {
	  string hTitle = " ";
	  stringstream hN6850;
	  int kHist = kOffHi + 6850 + kh + 1;
	  hN6850 << kHist;
	  vRC6850.push_back( new CsHist2D( hN6850.str(), hTitle,
					   100, -10., 10., 100, 0., 360. ) );
	}
	for( int kh=12; kh<15; kh++ ) {
	  string hTitle = " ";
	  stringstream hN6850;
	  int kHist = kOffHi + 6850 + kh + 1;
	  hN6850 << kHist;
	  vRC6850.push_back( new CsHist2D( hN6850.str(), hTitle,
					   100, -1200., 1200.,
					   100, -1200., 1200. ) );
	}
	CsHistograms::SetCurrentPath("/");
      }
    }

    double likeDV = cons->likeDefVa();

//- loop over part-photons :
//  ------------------------
    bool isFirst = true;
    CsRCEventPartPhotons* paphos = CsRCEventPartPhotons::Instance();
    list<CsRCPartPhotons*> lPartPhotons = paphos->lPartPhotons();
    list<CsRCPartPhotons*>::iterator ipf;
    for( ipf=lPartPhotons.begin(); ipf!=lPartPhotons.end(); ipf++ ) {
      if( !(*ipf)->flag() ) continue;
      CsRCPartPhotons* papho = (*ipf);
      CsRCRing* ring = papho->pRing();
      if( !ring ) continue;

      CsRCParticle* part = papho->pPart();
      float momPart = part->mom();

      double* probaLK;
      probaLK = papho->probaLKAll();
//--- Production Likelihoods
      double likeBkg = papho->probaLKBgAll();
      double likeElec = probaLK[ 2];
      double likeMuon = probaLK[ 5];
      double likePion = probaLK[ 8];
      double likeKaon = probaLK[11];
      double likeProton = probaLK[14];
      bool bIDpart[4];
      if( likeElec == likeBkg ) continue;      // NO photons from signal

//--- Select ID pions (analysis way)
      bool bIDpion = false;
      if( likePion > likeDV  &&  likePion > likeBkg  &&
          likePion > likeKaon  &&  likePion > likeProton )
	bIDpion = true;
      if( momPart < 8.  &&  likePion < likeElec ) bIDpion = false;
      bIDpart[0] = bIDpion;
//--- Select ID kaons (analysis way)
      bool bIDkaon = false;
      if( likeKaon > likeDV  &&  likeKaon > likeBkg  &&
          likeKaon > likePion  &&  likeKaon > likeProton )
	bIDkaon = true;
      if( momPart < 8.  &&  likeKaon < likeElec ) bIDkaon = false;
      bIDpart[1] = bIDkaon;
//--- Select ID protons (analysis way) ?
      bool bIDproton = false;
      if( likeProton > likeDV  &&  likeProton > likeBkg  &&
          likeProton > likePion  &&  likeProton > likeKaon )
	bIDproton = true;
      if( momPart < 8.  &&  likeProton < likeElec ) bIDproton = false;
      bIDpart[2] = bIDproton;

      double theReco = ring->theReco();
      double theRecoFt = ring->thetaRFit();
      double theLike = 0.;
      if( papho->thetaLikeSet() ) theLike = papho->thetaLike();

      double theRing = theReco;
      //double theRing = theRecoFt;
      //double theRing = theLike;
//@-----------------------------

      bool bSele = true;

//    using FIXED limits
//--- mom lower limit
      double momLLimit =  5.;       // GeV/c
      if( momPart < momLLimit ) bSele = false;
//@@-----------------------------------------
//--- mom upper limit
      double momULimit = 50.;       // GeV/c
      if( momPart > momULimit ) bSele = false;
//@@-----------------------------------------
//--- particle angle at entrance window :
      double tgxy;
      tgxy = sqrt( pow( part->vDirIn().x()/part->vDirIn().z(), 2 ) +
  	  	   pow( part->vDirIn().y()/part->vDirIn().z(), 2 ) );
      tgxy *= 1000.;
//--- tgxy lower limit
      double tgxyLLimit = 20.;      // ~mrad
      if( tgxy < tgxyLLimit ) bSele = false;
//@@---------------------------------------
//--- nPhoton lower limit
      int nPhoLLimit = 10;
      //int nPhoLLimit =  5;
      if( int( ring->lPhotons().size() ) < nPhoLLimit ) bSele = false;
//@@-----------------------------------------------------------------

      if( bSele ) {

//----- loop over 'photons' of each RING :
        list<CsRCPhoton*> lPhotons = ring->lPhotons();
        list<CsRCPhoton*>::iterator ih;
	bool isPMT = false;
	bool isAPV = false;
        for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	  if( (*ih)->isPMT() ) isPMT = true;
	  if( (*ih)->isAPV() ) isAPV = true;
	}
	bool PMTonly = false;
	bool APVonly = false;
	if( !isAPV  &&  isPMT ) PMTonly = true;
	if( isAPV  &&  !isPMT ) APVonly = true;
	int kDetPart = papho->kDetPart();
	double xPade = part->vPade()[kDetPart].x();
	double yPade = part->vPade()[kDetPart].y();
	int kHQ = 0;
	if( xPade < 0.  &&  yPade > 0. ) kHQ = 0;
	if( xPade > 0.  &&  yPade > 0. ) kHQ = 1;
	if( xPade < 0.  &&  yPade < 0. ) kHQ = 2;
	if( xPade > 0.  &&  yPade < 0. ) kHQ = 3;
	int kHP = 0;
	if( PMTonly ) kHP = 12;
	else if( APVonly ) kHP = 14;
	else  kHP = 13;
	xh = xPade;
	yh = yPade;
	if( kHP >= 0  &&  kHP < int( vRC6850.size() ) ) {
	  if( vRC6850[kHP] ) vRC6850[kHP]->Fill( xh, yh );
//hh                         ----------------------------
	}
	int kHH = 0;
	if( PMTonly ) kHH = 3*kHQ;
	else if( APVonly ) kHH = 3*kHQ + 2;
	else  kHH = 3*kHQ + 1;
        for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
          double theta = (*ih)->the();
	  double phi = (*ih)->phi();
	  //double phi = (*ih)->phiA();
//@-----------------------------------
	  xh = theta - theRing;
	  yh = phi;
          if( kHH >= 0  &&  kHH < int( vRC6850.size() ) ) {
	    if( vRC6850[kHH] ) vRC6850[kHH]->Fill( xh, yh );
//hh                           ----------------------------
	  }
	}
      }

    }

    return;
  }


//===========================================================================
  void CsRCEventRings::checkRingFit() {
//-------------------------------------

//- Paolo  -  March  2009


    if( !CsRCHistos::Ref().bookHis()  ||
        CsRCHistos::Ref().levelBk() < 2 ) return;
//  --------------------------------------------

    CsRCExeKeys* key = CsRCExeKeys::Instance();
    CsRCRecConst* cons = CsRCRecConst::Instance();
    double radDeg = cons->RadDeg();
    CsRCDetectors *dets = CsRCDetectors::Instance();
    CsRCHistos& hist = CsRCHistos::Ref();
    int kOffHi = cons->kOffHi();
    double xh, yh, wh;
    int kHH;

    static std::vector<CsHist2D*> vRC8300;
    static std::vector<CsHist2D*> vRC8400;

    bool print = false;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      for( int kh=0; kh<60; kh++ ) vRC8300.push_back( NULL );
      for( int kh=0; kh<60; kh++ ) vRC8400.push_back( NULL );
      if( CsRCHistos::Ref().bookHis() &&
	  CsRCHistos::Ref().levelBk() >= 2 ) {
	CsHistograms::SetCurrentPath("/RICH");
	int kHist = 0;
	string hTitle = " ";
	int khh = 0;
	khh = 8300;
	vRC8300.clear();
	for( int kc=0; kc<16; kc++ ) {
	  kHist = kOffHi + khh + 5*kc + 1;
	  stringstream hN8301;
          hN8301 << kHist;
          vRC8300.push_back( new CsHist2D( hN8301.str(), hTitle,
					   100, -4., 4., 100, -4., 4. ) );
	  kHist = kOffHi + khh + 5*kc + 2;
	  stringstream hN8302;
          hN8302 << kHist;
          vRC8300.push_back( new CsHist2D( hN8302.str(), hTitle,
					   100, -1.0, 1.0, 100, 0., 100. ) );
	  kHist = kOffHi + khh + 5*kc + 3;
	  stringstream hN8303;
          hN8303 << kHist;
          vRC8300.push_back( new CsHist2D( hN8303.str(), hTitle,
					   100, -5., 5., 100, 0., 100. ) );
	  kHist = kOffHi + khh + 5*kc + 4;
	  stringstream hN8304;
          hN8304 << kHist;
          vRC8300.push_back( new CsHist2D( hN8304.str(), hTitle,
					   100, -4., 4., 100, 0., 100. ) );
	  kHist = kOffHi + khh + 5*kc + 5;
	  stringstream hN8305;
          hN8305 << kHist;
          vRC8300.push_back( new CsHist2D( hN8305.str(), hTitle,
					   100, -4., 4., 100, 0., 100. ) );
	}
	for( int kk=0; kk<2; kk++ ) {
	  kHist = kOffHi + khh + 80 + 2*kk + 1;
	  stringstream hN8351;
          hN8351 << kHist;
          vRC8300.push_back( new CsHist2D( hN8351.str(), hTitle,
					   100, 0., 10., 1, 0., 1. ) );
	  kHist = kOffHi + khh + 80 + 2*kk + 2;
	  stringstream hN8352;
          hN8352 << kHist;
          vRC8300.push_back( new CsHist2D( hN8352.str(), hTitle,
					   20, 0., 20., 1, 0., 1. ) );
	}
	khh = 8400;
	vRC8400.clear();
	for( int kc=0; kc<16; kc++ ) {
	  kHist = kOffHi + khh + 5*kc + 1;
	  stringstream hN8401;
          hN8401 << kHist;
          vRC8400.push_back( new CsHist2D( hN8401.str(), hTitle,
					   100, -10., 10., 100, -10., 10. ) );
	  kHist = kOffHi + khh + 5*kc + 2;
	  stringstream hN8402;
          hN8402 << kHist;
          vRC8400.push_back( new CsHist2D( hN8402.str(), hTitle,
					   100, -10., 10., 100, 0., 100. ) );
	  kHist = kOffHi + khh + 5*kc + 3;
	  stringstream hN8403;
          hN8403 << kHist;
          vRC8400.push_back( new CsHist2D( hN8403.str(), hTitle,
					   100, -5., 5., 100, 0., 100. ) );
	  kHist = kOffHi + khh + 5*kc + 4;
	  stringstream hN8404;
          hN8404 << kHist;
          vRC8400.push_back( new CsHist2D( hN8404.str(), hTitle,
					   100, -10., 10., 100, 0., 100. ) );
	  kHist = kOffHi + khh + 5*kc + 5;
	  stringstream hN8405;
          hN8405 << kHist;
          vRC8400.push_back( new CsHist2D( hN8405.str(), hTitle,
					   100, -10., 10., 100, 0., 100. ) );
	}
	for( int kk=0; kk<2; kk++ ) {
	  kHist = kOffHi + khh + 80 + 2*kk + 1;
	  stringstream hN8451;
          hN8451 << kHist;
          vRC8400.push_back( new CsHist2D( hN8451.str(), hTitle,
					   100, 0., 10., 1, 0., 1. ) );
	  kHist = kOffHi + khh + 80 + 2*kk + 2;
	  stringstream hN8452;
          hN8452 << kHist;
          vRC8400.push_back( new CsHist2D( hN8452.str(), hTitle,
					   20, 0., 20., 1, 0., 1. ) );
	}
	CsHistograms::SetCurrentPath("/");
      }
    }

    CsRCEventPartPhotons* evepapho = CsRCEventPartPhotons::Instance();
    list<CsRCPartPhotons*> lPaPhotons = evepapho->lPartPhotons();
    list<CsRCPartPhotons*>::iterator ipf;
    for( ipf=lPaPhotons.begin(); ipf!=lPaPhotons.end(); ipf++ ) {
      if( !(*ipf) ) continue;
      CsRCPartPhotons* papho = (*ipf);
      if( !papho->flag() ) continue;

      CsRCRing* ring = papho->pRing();
      if( !ring ) continue;
      if( !ring->flag() ) continue;

      list<CsRCPhoton*> lPhotons = ring->lPhotons();
      list<CsRCPhoton*>::iterator ih;
      bool PMTonly = true;
      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
        if( !(*ih)->flag() ) continue;
	if( (*ih)->isPMT() ) continue;
	PMTonly = false;
	break;
      }
      bool APVonly = true;
      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
        if( !(*ih)->flag() ) continue;
	if( (*ih)->isAPV() ) continue;
	APVonly = false;
	break;
      }

//    determine if RING photons are ONE Cathode only
      int cClu = -1;
      int kCat = -1;
      bool PMTOneCat = false;
      if( PMTonly ) {
        for( int kc=0; kc<16; kc++ ) {
	  kCat = kc;
	  for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	    if( !(*ih)->flag() ) continue;
	    cClu = (*ih)->ptToClu()->ic();
	    if( cClu == kCat ) continue;
	    kCat = -1;
	    break;
	  }
	  if( kCat >= 0 ) {
	    PMTOneCat = true;
	    break;
	  }
        }
      }
      bool APVOneCat = false;
      if( APVonly ) {
        for( int kc=0; kc<16; kc++ ) {
	  kCat = kc;
	  for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	    if( !(*ih)->flag() ) continue;
	    cClu = (*ih)->ptToClu()->ic();
	    if( cClu == kCat ) continue;
	    kCat = -1;
	    break;
	  }
	  if( kCat >= 0 ) {
	    APVOneCat = true;
	    break;
	  }
        }
      }
      if( PMTOneCat ) APVOneCat = false;

//--- n-Photon/Ring lower limit
      int nPhoLLimit = 0;
      if( PMTonly ) nPhoLLimit = 10;
      if( APVonly ) nPhoLLimit =  4;
      int nPhotons = lPhotons.size();
      if( nPhotons < nPhoLLimit ) continue;
//@@--------------------------------------

//--- Ring fit in S.Y. plane :
//    ------------------------
      bool exeFit = true;
      if( !PMTOneCat  &&  !APVOneCat ) exeFit = false;
      int kDE = 0;
      if( PMTOneCat ) kDE = 0;
      if( APVOneCat ) kDE = 1;
      if( exeFit ) {
	double theRing = ring->the();
	int nParam = 3;
	double param[nParam];
	param[0] = theRing;
	param[1] = 0.;
	param[2] = 0.;
	int iPaFit[nParam];
	iPaFit[0] = 0;
	iPaFit[1] = 0;
	iPaFit[2] = 0;

	double momPart = papho->pPart()->mom();
	double sigPhoRec = papho->sigmaPhoRec( momPart );
//                         -----------------------------
	double dSigC = cons->sigCut() * sigPhoRec;
//      -----------------------------------------
	double theLoL = theRing - dSigC;
	double theUpL = theRing + dSigC;
	int kPhot = 0;
	double thePhot[nPhotons];
	double phiPhot[nPhotons];
	double sigPhot[nPhotons];
	for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
  	  double theW = (*ih)->the();
          if( theW >= theLoL  &&  theW <= theUpL ) {
	    double phiW = (*ih)->phi();
	    double sigW = (*ih)->sigmaPhoPid( papho );
//                        ---------------------------
	    thePhot[kPhot] = theW * cos( phiW/radDeg );
	    phiPhot[kPhot] = theW * sin( phiW/radDeg );
	    sigPhot[kPhot] = sigW*sigW;
	    kPhot++;
	  }
	}
	int nPoint = kPhot;

        CsRCCircleFit oCirclev( nPoint, thePhot, phiPhot, sigPhot,
//      ----------------------------------------------------------
                                nParam, param, iPaFit );
        if( oCirclev.doChiFit() ) {
//          -------------------

	  double theFit = oCirclev.para()[0];
	  double xFit = oCirclev.para()[1];
	  double yFit = oCirclev.para()[2];
	  double chiSq = oCirclev.chiSquare();
	  int nIter = oCirclev.nIter();
	  HepVector pull = oCirclev.pull();
	  xh = xFit;
	  yh = yFit;
	  kHH = 5*kCat + 0;
	  if( kHH >= 0  &&  kHH < int( vRC8300.size() ) ) {
	    if( vRC8300[kHH] ) vRC8300[kHH]->Fill( xh, yh );
//hh                           ----------------------------
	  }
	  xh = theFit - theRing;
	  yh = nPoint;
	  kHH = 5*kCat + 1;
	  if( kHH >= 0  &&  kHH < int( vRC8300.size() ) ) {
	    if( vRC8300[kHH] ) vRC8300[kHH]->Fill( xh, yh );
//hh                           ----------------------------
	  }
	  xh = xFit;
	  yh = nPoint;
	  kHH = 5*kCat + 3;
	  if( kHH >= 0  &&  kHH < int( vRC8300.size() ) ) {
	    if( vRC8300[kHH] ) vRC8300[kHH]->Fill( xh, yh );
//hh                           ----------------------------
	  }
	  xh = yFit;
	  yh = nPoint;
	  kHH = 5*kCat + 4;
	  if( kHH >= 0  &&  kHH < int( vRC8300.size() ) ) {
	    if( vRC8300[kHH] ) vRC8300[kHH]->Fill( xh, yh );
//hh                           ----------------------------
	  }
	  xh = chiSq;
	  yh = 0.5;
	  kHH = 80 + 2*kDE;
	  if( kHH >= 0  &&  kHH < int( vRC8300.size() ) ) {
	    if( vRC8300[kHH] ) vRC8300[kHH]->Fill( xh, yh );
//hh                           ----------------------------
	  }
	  xh = nIter;
	  yh = 0.5;
	  kHH = 80 + 2*kDE + 1;
	  if( kHH >= 0  &&  kHH < int( vRC8300.size() ) ) {
	    if( vRC8300[kHH] ) vRC8300[kHH]->Fill( xh, yh );
//hh                           ----------------------------
	  }
	  kHH = 5*kCat + 2;
	  if( kHH >= 0  &&  kHH < int( vRC8300.size() ) ) {
	    for( int k=0; k<nPoint; k++ ) {
	      xh = pull[k];
	      yh = nPoint;
	      if( vRC8300[kHH] ) vRC8300[kHH]->Fill( xh, yh );
//hh                             ----------------------------
	    }
	  }
        }

      }

//--- Ring fit in Detector plane :
//    ----------------------------
      exeFit = true;
      if( !PMTOneCat  &&  !APVOneCat ) exeFit = false;
      kDE = 0;
      if( PMTOneCat ) kDE = 0;
      if( APVOneCat ) kDE = 1;
      if( exeFit ) {
	double theRing = ring->the();
	int nParam = 3;
	double param[nParam];
	param[0] = theRing * 3.20;
	param[1] = 0.;
	param[2] = 0.;
	int iPaFit[nParam];
	iPaFit[0] = 0;
	iPaFit[1] = 0;
	iPaFit[2] = 0;

	int kDetPart = papho->kDetPart();
	double xPade = papho->vPoPaDetW()[kDetPart].x();
	double yPade = papho->vPoPaDetW()[kDetPart].y();
	double momPart = papho->pPart()->mom();
	double sigPhoRec = papho->sigmaPhoRec( momPart );
//                         -----------------------------
        double sigC = cons->nSigAliMir();
//      --------------------------------
        double dSigC = sigC * sigPhoRec;
	double theLoL = theRing - dSigC;
	double theUpL = theRing + dSigC;
        double sigma = 5.;
//@@---------------------
	int kPhot = 0;
	double xClu[nPhotons];
	double yClu[nPhotons];
	double errClu[nPhotons];
	for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
  	  double theW = (*ih)->the();
          if( theW >= theLoL  &&  theW <= theUpL ) {
	    xClu[kPhot] = (*ih)->ptToClu()->xc() - xPade;
	    yClu[kPhot] = (*ih)->ptToClu()->yc() - yPade;
	    errClu[kPhot] = sigma*sigma;
	    kPhot++;
	  }
	}
	int nPoint = kPhot;

        CsRCCircleFit oCirclev( nPoint, xClu, yClu, errClu,
//      ---------------------------------------------------
                                nParam, param, iPaFit );
        if( oCirclev.doChiFit() ) {
//          -------------------

	  double RFit = oCirclev.para()[0];
	  double xFit = oCirclev.para()[1];
	  double yFit = oCirclev.para()[2];
	  double chiSq = oCirclev.chiSquare();
	  int nIter = oCirclev.nIter();
	  HepVector pull = oCirclev.pull();
	  xh = xFit;
	  yh = yFit;
	  kHH = 5*kCat + 0;
	  if( kHH >= 0  &&  kHH < int( vRC8400.size() ) ) {
	    if( vRC8400[kHH] ) vRC8400[kHH]->Fill( xh, yh );
//hh                           ----------------------------
	  }
	  xh = RFit - theRing * 3.20;
	  yh = nPoint;
	  kHH = 5*kCat + 1;
	  if( kHH >= 0  &&  kHH < int( vRC8400.size() ) ) {
	    if( vRC8400[kHH] ) vRC8400[kHH]->Fill( xh, yh );
//hh                           ----------------------------
	  }
	  xh = xFit;
	  yh = nPoint;
	  kHH = 5*kCat + 3;
	  if( kHH >= 0  &&  kHH < int( vRC8400.size() ) ) {
	    if( vRC8400[kHH] ) vRC8400[kHH]->Fill( xh, yh );
//hh                           ----------------------------
	  }
	  xh = yFit;
	  yh = nPoint;
	  kHH = 5*kCat + 4;
	  if( kHH >= 0  &&  kHH < int( vRC8400.size() ) ) {
	    if( vRC8400[kHH] ) vRC8400[kHH]->Fill( xh, yh );
//hh                           ----------------------------
	  }
	  xh = chiSq;
	  yh = 0.5;
	  kHH = 80 + 2*kDE;
	  if( kHH >= 0  &&  kHH < int( vRC8400.size() ) ) {
	    if( vRC8400[kHH] ) vRC8400[kHH]->Fill( xh, yh );
//hh                           ----------------------------
	  }
	  xh = nIter;
	  yh = 0.5;
	  kHH = 80 + 2*kDE + 1;
	  if( kHH >= 0  &&  kHH < int( vRC8400.size() ) ) {
	    if( vRC8400[kHH] ) vRC8400[kHH]->Fill( xh, yh );
//hh                           ----------------------------
	  }
	  kHH = 5*kCat + 2;
	  if( kHH >= 0  &&  kHH < int( vRC8400.size() ) ) {
	    for( int k=0; k<nPoint; k++ ) {
	      xh = pull[k];
	      yh = nPoint;
	      if( vRC8400[kHH] ) vRC8400[kHH]->Fill( xh, yh );
//hh                             ----------------------------
	    }
	  }
        }

      }

    }

    return;
  }

