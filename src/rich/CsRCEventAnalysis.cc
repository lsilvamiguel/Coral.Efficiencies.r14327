
/*!
   \File    CsRCEventAnalysis.cc
   \----------------------------
   \brief   CsRCEventAnalysis class implementation.
   \author  Paolo Schiavon
   \version 0.10,  rev. 28/6/00
   \date    6 December 1999,  rev.  October  2000
*/


  #include <iostream>
  #include <ostream>
  #include <cstdio>

  #include "CsHist.h"
  #include "CsHistograms.h"
  #include "CsMCHit.h"
  #include "CsMCRICH1Hit.h"
  #include "CsMCDigit.h"
  #include "CsGeom.h"
  #include "CsRICH1Detector.h"
  #include "CsMCTrack.h"
  #include "CsTrack.h"

//------------------------------
  #include "CsRCEventAnalysis.h"

  #include "CsRichOne.h"

  #include "CsRCDetectors.h"
  #include "CsRCMirrors.h"

  #include "CsRCEventPads.h"
  #include "CsRCPad.h"
  #include "CsRCEventClusters.h"
  #include "CsRCCluster.h"

  #include "CsRCEventParticles.h"

  #include "CsRCEventPartPhotons.h"
  #include "CsRCPartPhotons.h"
  #include "CsRCPhoton.h"

  #include "CsRCEventRings.h"
  #include "CsRCRing.h"

  #include "CsRCChiSqFit.h"
  #include "CsRCCircleFit.h"
  #include "CsRCEllipseFit.h"
  #include "CsRCEllipseFitTest.h"   //*
  #include "CsRCGauPolFit.h"

  #include "CsRCExeKeys.h"
  #include "CsRCRecConst.h"
  #include "CsRCHistos.h"

  #include "CsRCLikeAll.h"
  #include "CsRCLikeAll05.h"

  #include "CsRCUty.h"
//------------------------------

  using namespace std;
  using namespace CLHEP;

  CsRCEventAnalysis* CsRCEventAnalysis::instance_ = 0;

//===========================================================================
  CsRCEventAnalysis* CsRCEventAnalysis::Instance() {
//--------------------------------------------------
    if( instance_ == 0 ) instance_ = new CsRCEventAnalysis();
    return instance_;
  }

//===========================================================================
  CsRCEventAnalysis::CsRCEventAnalysis() {
//----------------------------------------
    nMCProcPart_ = 0;
    nDaProcPart_ = 0;
    nMCRecPart_ = 0;
    nDaRecPart_ = 0;
    flag_ = false;
    flagMCMon_ = false;
    flagDataMon_ = false;
    flagPID_ = false;
    flagAliM_ = false;
  }

//===========================================================================
  CsRCEventAnalysis::CsRCEventAnalysis( const CsRCEventAnalysis &evea ) {
//-----------------------------------------------------------------------
    cout << "RICHONE : CsRCEventAnalysis CopyConstructor" << endl;
    instance_ = evea.instance_;

    nMCProcPart_ = evea.nMCProcPart_;
    nDaProcPart_ = evea.nDaProcPart_;
    nMCRecPart_ = evea.nMCRecPart_;
    nDaRecPart_ = evea.nDaRecPart_;

    flag_ = evea.flag_;
    flagMCMon_ = evea.flagMCMon_;
    flagDataMon_ = evea.flagDataMon_;
    flagPID_ = evea.flagPID_;
    flagAliM_ = evea.flagAliM_;
  }

//===========================================================================
  void CsRCEventAnalysis::print() {
//---------------------------------
    CsRCRecConst *cons = CsRCRecConst::Instance();

    double likeDV = cons->likeDefVa();
//@@---------------------------------

    if( flagPID_ ) {

      if( cons->likeType() == "ALL" ) {
	cout << endl;
	cout << " Particle likelihoods on accepted particles :" << endl;
	CsRCEventPartPhotons* evepapho = CsRCEventPartPhotons::Instance();
	list<CsRCPartPhotons*> lPaPhotons = evepapho->lPartPhotons();
	list<CsRCPartPhotons*>::iterator ipf;
	for( ipf=lPaPhotons.begin(); ipf!=lPaPhotons.end(); ipf++ ) {
	  if( !(*ipf) ) continue;
	  if( !(*ipf)->flag() ) continue;

	  //for( int k=0; k<4; k++ ) cout << (*ipf)->partProbs()[k] << "  ";
	  //cout << endl;

	  std::cout << setprecision( 4 );
          cout << "   PartPhot  " << (*ipf)->kPaPhot() << ",  Particle  "
               << (*ipf)->pPart()->kPart();
	  double probp = (*ipf)->partProbs()[1];
	  if( probp != likeDV ) cout << ",  like Pion  " << probp;

	  double probk = (*ipf)->partProbs()[2];
	  if( probk != likeDV ) cout << ",  like Kaon  " << probk;

	  double probr = (*ipf)->partProbs()[3];
	  if( probr != likeDV ) cout << ",  like Proton  " << probr;

	  double probb = (*ipf)->partProbs()[0];
	  if( probb != likeDV ) cout << ",  like Backgr.  " << probb;

	  std::cout << "    fr. " << (*ipf)->fracUse();
	  cout << endl;
	  cout << "                               ";
	  double dprob = (*ipf)->partProbs()[4];
	  if( probp != likeDV ) cout << "dL/dn Pion  " << dprob;

	  dprob = (*ipf)->partProbs()[5];
	  if( probk != likeDV ) cout << ",  dL/dn Kaon  " << dprob;

	  dprob = (*ipf)->partProbs()[6];
	  if( probr != likeDV ) cout << ",  dL/dn Proton  " << dprob;
	  cout << endl;

          std::cout << setprecision( 2 );
	  CsRCRing* ring = (*ipf)->pRing();
	  if( ring  &&  ring->flag() ) {
	    cout << "                               ";
	    double prob = ring->partProbs()[12];
	    if( prob < 10000000. ) cout << "chisq Pion  " << prob;

	    prob = ring->partProbs()[13];
            if( prob < 10000000. ) cout << ",  chisq Kaon  " << prob;

	    prob = ring->partProbs()[14];
            if( prob < 10000000. ) cout << ",  chisq Proton  " << prob;
	    cout << endl;

	    cout << "                               ";
	    cout << "thetaRec = " <<  ring->partProbs()[8];
	    //cout << ",  thetaLike = " <<  ring->partProbs()[7];
	    cout << ",  thetaLike = " <<  (*ipf)->partProbs()[7];
	    cout << ",  thetaFit = " <<  ring->partProbs()[10];
            std::cout << setprecision( 0 );
	    cout << ",  photons = " <<  ring->partProbs()[9];
            std::cout << setprecision( 2 );
	    cout << ",  ring chisq = " <<  ring->partProbs()[11];
	    cout << endl;
	  } else {
	    cout << "                               ";
	    cout << "thetaLike = " <<  (*ipf)->partProbs()[7];
	    cout << endl;
	  }
	}
      }

      if( cons->likeType() == "RING" ) {
	cout << endl;
	cout << " Particle likelihoods on accepted rings :" << endl;
	CsRCEventRings* evrings = CsRCEventRings::Instance();
	list<CsRCRing*> lRings = evrings->lRings();
	list<CsRCRing*>::iterator ir;
	for( ir=lRings.begin(); ir!=lRings.end(); ir++ ) {
	  if( !(*ir) ) continue;
	  if( !(*ir)->flag() ) continue;

          cout << "   Ring  " << (*ir)->kRing() << ",  Particle  "
               << (*ir)->pPartPhot()->pPart()->kPart();
	  double probp = (*ir)->partProbs()[1];
	  if( probp >= 0. ) cout << ",  like Pion  " << probp;
	  double probk = (*ir)->partProbs()[2];
	  if( probk >= 0. ) cout << ",  like Kaon  " << probk;
	  double probr = (*ir)->partProbs()[3];
	  if( probr >= 0. ) cout << ",  like Proton  " << probr;
	  double probb = (*ir)->partProbs()[0];
	  if( probb >= 0. ) cout << ",  like Backgr.  " << probb;
	  cout << endl;
	  cout << "                           ";
	  double dprob = (*ir)->partProbs()[4];
	  if( probp >= 0. ) cout << "dL/dn Pion  " << dprob;
	  dprob = (*ir)->partProbs()[5];
	  if( probk >= 0. ) cout << ",  dL/dn Kaon  " << dprob;
	  dprob = (*ir)->partProbs()[6];
	  if( probr >= 0. ) cout << ",  dL/dn Proton  " << dprob;
	  cout << endl;

	  cout << "                           ";
	  double prob = (*ir)->partProbs()[12];
	  if( prob >= 0. ) cout << "chisq Pion  " << prob;
	  prob = (*ir)->partProbs()[13];
	  if( prob >= 0. ) cout << ",  chisq Kaon  " << prob;
	  prob = (*ir)->partProbs()[14];
	  if( prob >= 0. ) cout << ",  chisq Proton  " << prob;
	  cout << endl;
	  cout << "                           ";
	  cout << "thetaRec = " <<  (*ir)->partProbs()[8];
	  cout << ",  thetaLike = " <<  (*ir)->partProbs()[7];
	  cout << ",  thetaFit = " <<  (*ir)->partProbs()[10];
	  cout << ",  photons = " <<  (*ir)->partProbs()[9];
	  cout << ",  ring chisq = " <<  (*ir)->partProbs()[11];
	  cout << endl;
	}
      }

      if( cons->likeType() == "CHISQ" ) {
	cout << endl;
	cout << " Particle chi-squares on accepted rings :" << endl;
	CsRCEventRings* evrings = CsRCEventRings::Instance();
	list<CsRCRing*> lRings = evrings->lRings();
	list<CsRCRing*>::iterator ir;
	for( ir=lRings.begin(); ir!=lRings.end(); ir++ ) {
	  if( !(*ir) ) continue;
	  if( !(*ir)->flag() ) continue;

          cout << "   Ring  " << (*ir)->kRing() << ",  Particle  "
               << (*ir)->pPartPhot()->pPart()->kPart();
	  double prob = (*ir)->partProbs()[12];
	  if( prob >= 0. ) cout << "   chisq Pion  " << prob;
	  prob = (*ir)->partProbs()[13];
	  if( prob >= 0. ) cout << ",  chisq Kaon  " << prob;
	  prob = (*ir)->partProbs()[14];
	  if( prob >= 0. ) cout << ",  chisq Proton  " << prob;
	  cout << endl;
	  cout << "                           ";
	  cout << "thetaRec = " <<  (*ir)->partProbs()[8];
	  //cout << ",  thetaLike = " <<  (*ir)->partProbs()[7];
	  cout << ",  thetaFit = " <<  (*ir)->partProbs()[10];
	  cout << ",  photons = " <<  (*ir)->partProbs()[9];
	  cout << ",  ring chisq = " <<  (*ir)->partProbs()[11];
	  cout << endl;
	}
      }
    }
  }

//===========================================================================
  CsRCEventAnalysis::~CsRCEventAnalysis() {}
//------------------------------------------


  extern "C" {
    float prob_( const float&, const int& );
    float erf_( const float& );
  }

//===========================================================================
  void CsRCEventAnalysis::doEveAnalysis() {
//-----------------------------------------


//--- data monitoring and ring analysis.
//    ----------------------------------
//--- Paolo  -  December 1999,  rev. October 2000


//--- from "CsRCExecKeys"
      CsRCExeKeys *key = CsRCExeKeys::Instance();
      bool MCMoni = key->MCMoni();
      bool DataMoni = key->DataMoni();
      bool PartIdent = key->PartIdent();
      bool AliMirror = key->AliMirror();

//--- from "CsRCRecConst"
      CsRCRecConst *cons = CsRCRecConst::Instance();
      int kOffHi = cons->kOffHi();

      CsRCEventParticles* evparts = CsRCEventParticles::Instance();

      bool MCarloEvent = key->MCarloEvent();

      static bool firstCall = true;
      if( firstCall ) {
	firstCall = false;

        key->acknoMethod( "CsRCEventAnalysis" );
      }

//--- monitor of MC variables :
//    -------------------------
      flagMCMon_ = false;
      if( MCMoni  &&  MCarloEvent )  { flagMCMon_ = true;  MCMonitor(); }
//                                                         -----------   

//--- monitor of data variables :
//    ---------------------------
      flagDataMon_ = false;
      if( DataMoni )  { flagDataMon_ = true;  dataMonitor(); }
//                                            -------------

//--- monitor of particle identification :
//    ------------------------------------
      flagPID_ = false;
      if( PartIdent )  { flagPID_ = true;  PIDMonitor(); }
//                                         ------------

//--- mirror element alignment :
//    --------------------------
      flagAliM_ = false;
      if( AliMirror ) { 
	flagAliM_ = true;
	//CsRCMirrors::Instance()->doAliMirrors();
	CsRCMirrors::Instance()->doAliMirrPhi();
//                               --------------
	//CsRCMirrors::Instance()->doAliMirrAll();
//                               --------------
	doMirrDetFit();
//      --------------
      }

//--- Output from RICH :
//    ------------------
//--- moved to CsRichOne::doRichOne() : 070430
      //setProbTsoTrack();
//    -----------------

      flag_ = true;

//--- conditional prints :
//    --------------------
      int kPrintEventAnalysis = key->kPrintEventAnalysis();
      if( kPrintEventAnalysis == 1 ) print();

  }


//===========================================================================
  void CsRCEventAnalysis::MCMonitor() {
//-------------------------------------


//--- monitor with MC variables
//    -------------------------
//--- Paolo  -  December 1999,  rev. October  2000



      CsRCExeKeys *key = CsRCExeKeys::Instance();

      static bool readMyFile;

      static bool firstCall = true;
      if( firstCall ) {
	firstCall = false;

        key->acknoMethod( "CsRCEventAnalysis::MCMonitor" );

	readMyFile = key->readMyFile();
      }

//--- from "CsRCRecConst"
      CsRCRecConst *cons = CsRCRecConst::Instance();
      int kOffHi = cons->kOffHi();
      float CFRefInd = cons->CFRefInd();
      double *massPartv = cons->massPartv();
      int nPhotMinRing = cons->nPhotMinRing();

      CsRCHistos& hist = CsRCHistos::Ref();
      int kh;
      double xh, yh, wh;

      bool MCarloEvent = key->MCarloEvent();

//--- 090217 ------
      checkMCPMTs();
//-----------------

//--- loop over rings :
//    -----------------
      CsRCEventRings* evrings = CsRCEventRings::Instance();
      list<CsRCRing*> lRings = evrings->lRings();
      list<CsRCRing*>::iterator ir;
      for( ir=lRings.begin(); ir!=lRings.end(); ir++ ) {
        if( !(*ir)->flag() ) continue;

//----- select particle type (MC) :
//      ---------------------------
        CsRCEventParticles* evparts = CsRCEventParticles::Instance();
	CsRCParticle* part = (*ir)->pPartPhot()->pPart();
	if( !evparts->partSelect( part ) ) continue;
//@@                 ------------------

        float momPart = part->mom();

//----- photons in ring :
//      -----------------
        int nPhoRing = (*ir)->lPhotons().size();

//----- reconstructed Cerenkov angle :
//      ------------------------------
	double theReco = (*ir)->theReco();
//                       ----------------
	double theRecoWg = (*ir)->getThetaWgav();
//                         ---------------------
	double theRecoFt = (*ir)->getThetaRFit();
//                         ---------------------
	double theMC = part->thetaCer();
	xh = theMC - theRecoWg;
	if( hist.hRC3631 ) hist.hRC3631-> Fill( xh );
//hh                       -------------------------
	//xh = theMC - theRecoFt;
	//yh = theMC;
	//if( hist.hRC3632 ) hist.hRC3632-> Fill( xh, yh );
//hh                       -----------------------------
	if( (*ir)->pPartPhot()->howPMT() > 0.8 ) {
	  //list<CsRCPhoton*> lPhotons = (*ir)->pPartPhot()->lPhotons();
	  list<CsRCPhoton*> lPhotons = (*ir)->lPhotons();
	  list<CsRCPhoton*>::iterator ih;
	  double theRing = 0.;
	  int kPhot = 0;
	  for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	    if( !(*ih)->flag() ) continue;
	    theRing += (*ih)->the0();
	    kPhot++;
	  }
	  if( kPhot > 0 ) {
	    theRing /= float( kPhot );
	    xh = theMC - theRing;
	    yh = theMC;
	    if( hist.hRC3632 ) hist.hRC3632-> Fill( xh, yh );
//hh                           -----------------------------
	  }
	}

        double betaReco = 1./( cos( theReco/1000. ) * CFRefInd );
        if( betaReco > 1.) betaReco = 1.;

//----- correct reconstructed theta for E losses :
//      ------------------------------------------
	//double theRecoCo = theReco - offTheRing( betaReco );
//@@-------------------------------------------------------

        int iPartT = part->iPartT();
        if( iPartT > 200) iPartT -= 200;
	if( iPartT > 30 ) iPartT = 0;               //   wrong myFile?
        double massPart = massPartv[iPartT];

        //xh = iPartT;
        //if( hist.hRC3615 ) hist.hRC3615->Fill( xh );
//hh                       ------------------------


//----- define ACCEPTED rings :
//      -----------------------
        bool bAccp = true;
//@@---------------------

//----- select particle charge ( MC ) :
//      -------------------------------
        int chgPart = part->charge();
	//        if( chgPart != +1 ) bAccp = false;
	//        if( chgPart != -1 ) bAccp = false;

//----- particle up to the exit window : ( for compatibility with myFile )
//      --------------------------------
        if( part->ioExWd() != 2 ) bAccp = false;
//@@---                           -------------

//----- enough clusters in the ring :
//      -----------------------------
	// !!!        if( nPhoRing <= 0 ) bAccp = false;
        if( nPhoRing < nPhotMinRing ) bAccp = false;
//@@---                               -------------

        if( bAccp ) {

//------- plot (accepted) particle momenta :
//        ----------------------------------
          xh = momPart;
          int khw = -1;
          if( iPartT ==  2  || iPartT ==  3 ) khw = 0;
          if( iPartT ==  5  || iPartT ==  6 ) khw = 1;
          if( iPartT ==  8  || iPartT ==  9 ) khw = 2;
          if( iPartT == 11  || iPartT == 12 ) khw = 3;
          if( iPartT == 14  || iPartT == 15 ) khw = 4;
          if( khw >= 0 ) { 
	    if( hist.vRC3620[khw] ) hist.vRC3620[khw] ->Fill( xh );
//hh                                ------------------------------
	  }

//------- enough photons in the ring :
//        ----------------------------
          if( nPhoRing > 0 ) {
	  // !!!          if( nPhoRing >= nPhotMinRing ) {
//@@-------------------------------------
            double betaCalc = momPart /
              sqrt( momPart*momPart + massPart*massPart );
            double cosTheCa = 1./ (betaCalc * CFRefInd);

//--------- only particle over threshold :
//          ------------------------------
            if( cosTheCa <= 1. ) {

	      (*ir)->setFlagOverThrs( true );

              double theCalc = acos( cosTheCa ) * 1000.;
//@@---------------------------------------------------
	      // !!!	      double theCalc = part->thetaCer();

              float dTheta = theReco - theCalc;
              double cosTheRe = cos( theReco/1000. );

              nMCProcPart_++;

              xh = dTheta;
              yh = theCalc;
              if( hist.hRC3008 ) hist.hRC3008->Fill( xh, yh );
//hh                             ----------------------------
              xh = dTheta;
              yh = momPart;
              if( hist.hRC3009 ) hist.hRC3009->Fill( xh, yh );
//hh                             ----------------------------
              xh = dTheta;
              yh = betaCalc;
              if( hist.hRC3010 ) hist.hRC3010->Fill( xh, yh );
//hh                             ----------------------------
              xh = betaCalc;
              if( hist.hRC3016 ) hist.hRC3016->Fill( xh );
//hh                             ------------------------
              xh = betaCalc;
              yh = betaReco;
              if( hist.hRC3014 ) hist.hRC3014->Fill( xh, yh );
//hh                             ----------------------------
              xh = momPart;
              if( hist.hRC3022 ) hist.hRC3022->Fill( xh );
//hh                             ------------------------
              xh = nPhoRing;
              if( hist.hRC3025 ) hist.hRC3025->Fill( xh );
//hh                             ------------------------
              xh = theCalc;
	      yh = theReco;
              if( hist.hRC3034 ) hist.hRC3034->Fill( xh, yh );
//hh                             ----------------------------

//----------- check track hits vs ring hits :
//            -------------------------------
	      bool testHit = false;
	      //bool testHit = true;
//@@------------------------------
	      if( readMyFile ) testHit = false;
//-------------------------------------------------------------------
	      list<CsMCHit*> lRingHits;
              list<CsMCHit*> lTrackHits;
	      int nRHt;
	      int nEqu;
	      if( testHit ) {

		chkHitLists( (*ir), nRHt, nEqu,
			     lTrackHits, lRingHits );
		chkHitPrint( (*ir), nRHt, nEqu,
			     lTrackHits, lRingHits );
	
                xh = lRingHits.size();
                yh = lTrackHits.size(); 
                if( hist.hRC3101 ) hist.hRC3101->Fill( xh, yh );
//hh                               ----------------------------
	        hitDisplay( 5100, lTrackHits );
	        hitDisplay( 5200, lRingHits );
	      }
//-------------------------------------------------------------------


//----------- define 'reconstructed' rings :
//            ------------------------------
              bool bReco = false;
//@-----------------------------

              float recCut = (*ir)->recoCut( betaCalc );
//@@---------------------------------------------------
              if( fabs( theReco - theCalc ) < recCut ) bReco = true;
//@@---------                                          ------------

              xh = recCut;
              yh = betaCalc;
              if( hist.hRC3028 ) hist.hRC3028->Fill( xh, yh );
//hh                             ----------------------------

              if( bReco ) {

		(*ir)->setFlagReco( true );
//              --------------------------
		nMCRecPart_++;

                xh = dTheta;
                yh = nPhoRing;
                if( hist.hRC3012 ) hist.hRC3012->Fill( xh, yh );
//hh                               ----------------------------
                xh = dTheta;
                yh = betaCalc;
                if( hist.hRC3013 ) hist.hRC3013->Fill( xh, yh );
//hh                               ----------------------------
                xh = betaCalc;
                wh = nPhoRing;
	        if( hist.hRC3015 ) hist.hRC3015->Fill( xh, wh );
//hh                               ----------------------------
                xh = betaCalc;
       	        if( hist.hRC3017 ) hist.hRC3017->Fill( xh );
//hh                               ------------------------
                xh = betaCalc;
                wh = cosTheRe*betaCalc*CFRefInd;
	        if( hist.hRC3019 ) hist.hRC3019->Fill( xh, wh );
//hh                               ----------------------------
                xh = betaCalc;
                wh = theReco;
	        if( hist.hRC3020 ) hist.hRC3020->Fill( xh, wh );
//hh                               ----------------------------
                xh = momPart;
       	        if( hist.hRC3023 ) hist.hRC3023->Fill( xh );
//hh                               ------------------------
                xh = nPhoRing;
       	        if( hist.hRC3026 ) hist.hRC3026->Fill( xh );
//hh                               ------------------------

                if( betaCalc >= 0.99995 )  {
                  xh = dTheta;
                  if( hist.hRC3031 ) hist.hRC3031->Fill( xh );
//hh                                 ------------------------
                }
	        list<CsRCPhoton*> lPhotons = (*ir)->lPhotons();
	        list<CsRCPhoton*>::iterator ih;
		for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
                  double thePhot = (*ih)->the();
                  double phiPhot = (*ih)->phi();
                  xh = thePhot - theCalc;
                  if( hist.hRC3001 ) hist.hRC3001->Fill( xh );
//hh                                 ------------------------
                }
                xh = dTheta;
                if( hist.hRC3002 ) hist.hRC3002->Fill( xh );
//hh                               ------------------------
                xh = nPhoRing;
                yh = betaCalc;
                if( hist.hRC3011 ) hist.hRC3011->Fill( xh, yh );
//hh                               ----------------------------

//----------------------------------------------------------------------
		if( testHit ) {
		  chkHitPrint( (*ir), nRHt, nEqu,
			       lTrackHits, lRingHits );

                  xh = lRingHits.size();
                  yh = lTrackHits.size(); 
                  if( hist.hRC3111 ) hist.hRC3111->Fill( xh, yh );
//hh                                 ----------------------------
		}
//----------------------------------------------------------------------

              }   /* end if on bReco */

            }   /* end of if on cosTheCa */

          }   /* end if on nPhoRing */

        }  else  {     /* else on bAccp */

 
          xh = momPart;
          yh = iPartT;
          if( hist.hRC3029 ) hist.hRC3029->Fill( xh, yh );
//hh                         ----------------------------

        }   /* end if on bAccp */

      }   /* end of loop on Rings */

  }


//===========================================================================
  void CsRCEventAnalysis::dataMonitor() {
//---------------------------------------


//--- monitor of data
//    ---------------
//--- Paolo  -  December 1999,  rev.  October  2000



      CsRCExeKeys *key = CsRCExeKeys::Instance();

      static bool firstCall = true;
      if( firstCall ) {
	firstCall = false;

        key->acknoMethod( "CsRCEventAnalysis::DataMonitor" );
      }

      int kTyLL = 2;
//@@---------------

//--- from "ReconConst.h"
      CsRCRecConst* cons = CsRCRecConst::Instance();
      int kOffHi = cons->kOffHi();
      float CFRefInd = cons->CFRefInd();
      double TwoPI = cons->TwoPI();
      double RadDeg = cons->RadDeg();
      double *massPartv = cons->massPartv();
      int nPhotMinRing = cons->nPhotMinRing();
      float momMinProc = cons->momMinProc();

      CsRCHistos& hist = CsRCHistos::Ref();
      int kh;
      double xh, yh, wh;

//--- 090225
      ///checkDataPMTs();
//-------------------

//--- loop over rings :
//    -----------------
      int nRingEv = 0;
      CsRCEventRings* evrings = CsRCEventRings::Instance();
      list<CsRCRing*> lRings = evrings->lRings();
      list<CsRCRing*>::iterator ir;
      for( ir=lRings.begin(); ir!=lRings.end(); ir++ ) {
        if( !(*ir)->flag() ) continue;

//----- reconstructed Cerenkov angle :
//      ------------------------------
	double theReco = (*ir)->theReco();
//                       ----------------
	double theRecoWg = (*ir)->getThetaWgav();
//                         ---------------------
	double theRecoFt = (*ir)->getThetaRFit();
//                         ---------------------

        double betaReco = 1./( cos( theReco/1000. ) * CFRefInd );
        if( betaReco > 1.) betaReco = 1.;

//----- correct reconstructed theta for E losses (not used) :
//      -----------------------------------------------------
	//double theRecoCo = theReco - offTheRing( betaReco );
//@@-------------------------------------------------------

	CsRCParticle* part = (*ir)->pPartPhot()->pPart();
        float momPart = part->mom();

        xh = momPart;
        if( hist.hRC3611 ) hist.hRC3611->Fill( xh );
//hh                       ------------------------

//----- photons in ring :
//      -----------------
        int nPhoRing = (*ir)->lPhotons().size();


//----- define ACCEPTED rings :
//      -----------------------
        bool bAccp = true;
//@@---------------------

//----- particle mom. over pion threshold :
//      -----------------------------------
        if( momPart < momMinProc ) bAccp = false;
//@@---                            -------------

//----- enough photons in the ring :
//      ----------------------------
	// !!!	if( nPhoRing <= 0 ) bAccp = false;
        if( nPhoRing < nPhotMinRing ) bAccp = false;
//@@---                               -------------

        if( bAccp )  {

	  nRingEv++;

//------- plot (accepted) particle at entrance window :
//        ---------------------------------------------
          xh = part->xa();
          yh = part->ya();
          if( hist.hRC3601 ) hist.hRC3601->Fill( xh, yh );
//hh                         ----------------------------
          xh = part->ld();
          yh = part->md();
          if( hist.hRC3602 ) hist.hRC3602->Fill( xh, yh );
//hh                         ----------------------------
          xh = part->xa();
          yh = part->ld();
          if( hist.hRC3603 ) hist.hRC3603->Fill( xh, yh );
//hh                         ----------------------------
          xh = part->ya();
          yh = part->md();
          if( hist.hRC3604 ) hist.hRC3604->Fill( xh, yh );
//hh                         ----------------------------

//------- plot (accepted) particle momentum :
//        -----------------------------------
          xh = momPart;
          if( hist.hRC3610 ) hist.hRC3610->Fill( xh );
//hh                         ------------------------

//------- plot (accepted) angle particle-detector :
//        -----------------------------------------
	  int kDetPart = (*ir)->pPartPhot()->kDetPart();
          double thePade = part->thPade()[kDetPart];
          double thePart = acos( part->nd() );
          double phiPart = atan2( part->md(), part->ld() );
          if( phiPart < 0. ) phiPart += TwoPI;

          xh = thePade * RadDeg;
          yh = thePart * RadDeg;
          if( hist.hRC3519 ) hist.hRC3519->Fill( xh, yh );
//hh                         ----------------------------
          xh = thePade * RadDeg;
          yh = phiPart * RadDeg;
          if( hist.hRC3520 ) hist.hRC3520->Fill( xh, yh );
//hh                         ----------------------------

//------- photon-in-ring spectrum :
//        -------------------------
          xh = nPhoRing;
          if( hist.hRC3006 ) hist.hRC3006->Fill( xh );
//hh                         ------------------------

//------- cluster-in-ring multiplicity :
//        ------------------------------
	  list<CsRCPhoton*> lPhotons = (*ir)->lPhotons();
	  list<CsRCPhoton*>::iterator ih;
	  for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	    int nPadClu = (*ih)->ptToClu()->lPads().size();

            xh = nPadClu;
            if( hist.hRC3050 ) hist.hRC3050->Fill( xh );
//hh                           ------------------------
	  }
	  //^cluStructure( (*ir) );
//        ---------------------

//------- photon-in-ring spectrum for split rings :
//        -----------------------------------------
	  int kRingLoc = (*ir)->kRingLoc();
	  if( kRingLoc == +1 ) {
            xh = nPhoRing;
            if( hist.hRC3061 ) hist.hRC3061->Fill( xh );
//hh                           ------------------------
	  }
	  if( kRingLoc ==  0 ) {
            xh = nPhoRing;
            if( hist.hRC3062 ) hist.hRC3062->Fill( xh );
//hh                           ------------------------
          }
	  if( kRingLoc == -1 ) {
            xh = nPhoRing;
            if( hist.hRC3063 ) hist.hRC3063->Fill( xh );
//hh                           ------------------------
          }

	  int nPhoUp = 0;
	  int nPhoDw = 0;
	  for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	    int kDet = (*ih)->kDetClu();
	    if( kDet == 0 ) nPhoUp++;
	    if( kDet == 1 ) nPhoDw++;
	  }
          xh = nPhoUp;
          yh = nPhoDw;
          if( hist.hRC3064 ) hist.hRC3064->Fill( xh, yh );
//hh                         ----------------------------

          nDaProcPart_++;

          xh = betaReco;
          if( hist.hRC3036 ) hist.hRC3036->Fill( xh );
//hh                         ------------------------
          xh = momPart;
          if( hist.hRC3042 ) hist.hRC3042->Fill( xh );
//hh                         ------------------------
          xh = nPhoRing;
          if( hist.hRC3045 ) hist.hRC3045->Fill( xh );
//hh                         ------------------------
	  xh = theReco;
          wh = nPhoRing;
	  if( hist.hRC3065 ) hist.hRC3065->Fill( xh, wh );
//hh                         ----------------------------
          xh = theReco;
       	  if( hist.hRC3067 ) hist.hRC3067->Fill( xh );
//hh                         ------------------------

//------- data analysis histo.s (on rings) :
//        ----------------------------------
	  if( betaReco >= 0.99995 ) {
	    xh = theReco;
	    if( hist.hRC3052 ) hist.hRC3052->Fill( xh );
//hh                           ------------------------
	    for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	      xh = (*ih)->the();
	      if( hist.hRC3057 ) hist.hRC3057->Fill( xh );
//hh                             ------------------------
	      xh = (*ih)->the() - theReco;
	      if( hist.hRC3058 ) hist.hRC3058->Fill( xh );
//hh                             ------------------------
	    }
	  }

//------- n - 1 :
//        -------
	  moniCFRefInf( (*ir) );                   //   020830
//        ---------------------

//------- define 'reconstructed' rings :
//        ------------------------------
          bool bReco = false;
//@@------------------------

//------- check for pion, kaon or proton only (input = all) :
//        ---------------------------------------------------
          float recCut;
	  for( int kPaTy=kTyLL; kPaTy<=14; kPaTy+=3 ) {
//        -------------------------------------------

            float massIpo = massPartv[kPaTy];
            double betaIpo = momPart /
              sqrt( momPart*momPart + massIpo*massIpo );
	    double theIpo = (*ir)->pPartPhot()->thetaIpo( kPaTy );
	    if( theIpo > 0. ) {

              double dTheta = theReco - theIpo;

              xh = dTheta;
              yh = betaReco;
              if( hist.hRC3053 ) hist.hRC3053->Fill( xh, yh );
//hh                             ----------------------------
              if( betaIpo >= 0.99995 ) {
                xh = dTheta;
                if( hist.hRC3051 ) hist.hRC3051->Fill( xh );
//hh                               ------------------------
	      }

              recCut = (*ir)->recoCut( betaIpo );
//@@                   -------------------------
              if( fabs( dTheta ) < recCut ) bReco = true;
//                                          ------------

            }   /* end of if on theIpo */

          }   /* end of loop on kPaTy */

          if( bReco ) {

	    if( !key->MCarloEvent() ) (*ir)->setFlagReco( true );
//                                    --------------------------
	    nDaRecPart_++;

            xh = betaReco;
            if( hist.hRC3037 ) hist.hRC3037->Fill( xh );
//hh                           ------------------------
            xh = momPart;
            if( hist.hRC3043 ) hist.hRC3043->Fill( xh );
//hh                           ------------------------
            xh = nPhoRing;
            if( hist.hRC3046 ) hist.hRC3046->Fill( xh );
//hh                           ------------------------

	    bool exeFits = true;
	    //^bool exeFits = false;
//@@---------------------------
	    if( exeFits ) {

	    CsRCCircleFit oCircle;
	    if( (*ir)->getDetRFit( hist.vRC3710, oCircle ) );
//              ------------------------------------------
	    CsRCCircleFit* pCircle = &oCircle;

	    vector<CsHist2D*> vDummy;
	    //^CsRCEllipseFit oEllipse;
	    //^if( (*ir)->getDetEFit( vDummy, oEllipse ) );
//              -------------------------------------
	    //^CsRCEllipseFit* pEllipse = &oEllipse;

	    if( pCircle->flag()  &&  pCircle->nIter() < 10 ) {

	      int kDetPart = (*ir)->pPartPhot()->kDetPart();
	      double xPade = (*ir)->pPartPhot()->vPoPaDetW()[kDetPart].x();
	      double yPade = (*ir)->pPartPhot()->vPoPaDetW()[kDetPart].y();
	      double pass = 200.;
	      xPade += 4.* pass;;
	      int kPosX = int( xPade / pass );
	      if( kPosX > 7 ) kPosX = 7;
	      if( kPosX < 0 ) kPosX = 0;
	      if( yPade > 0. ) yPade -= 400.;
	      if( yPade < 0. ) yPade += 400.;
	      yPade += 4.* pass;
	      int kPosY = int( yPade / pass );
	      if( kPosY > 7 ) kPosY = 7;
	      if( kPosY < 0 ) kPosY = 0;

	      xh = pCircle->para()[1];
	      yh = pCircle->para()[2];
	      if( hist.vRC3720[kPosX] ) hist.vRC3720[kPosX]->Fill( xh, yh );
//                                      -----------------------------------
	      xh = pCircle->para()[1];
	      yh = pCircle->para()[2];
	      if( hist.vRC3730[kPosY] ) hist.vRC3730[kPosY]->Fill( xh, yh );
//                                      -----------------------------------
	    }
	    }

	    double ringNess = 0;
	    //if( (*ir)->getRingness( ringNess ) ) {
	    if( (*ir)->getRingness( (*ir)->lPhotons(), ringNess ) ) {
//              -------------------------------------------------
	      xh = ringNess;
	      if( hist.hRC3047 ) hist.hRC3047->Fill( xh );
//hh                             ------------------------
	    }

          }   /* end of if on bReco */


          xh = part->vPoPaMir0()[kDetPart].x();
          yh = part->vPoPaMir0()[kDetPart].y();
          if( kDetPart == 0 ) { 
            if( hist.hRC3054 ) hist.hRC3054->Fill( xh, yh );
//hh                           ----------------------------
	  }
          if( kDetPart == 1 ) { 
            if( hist.hRC3055 ) hist.hRC3055->Fill( xh, yh );
//hh                           ----------------------------
          }

          xh = (*ir)->pPartPhot()->vPoPaDetW()[kDetPart].x();
          yh = (*ir)->pPartPhot()->vPoPaDetW()[kDetPart].y();
          if( hist.hRC3056 ) hist.hRC3056->Fill( xh, yh );
//hh                         ----------------------------

        }   /* end of if on bAccp */

      }   /* end of loop on Rings */

//--- moved to CsRCEventRings
      //xh = nRingEv;
      //if( hist.hRC1537 ) hist.hRC1537->Fill( xh );
//hh                     ------------------------

      moniCFRefIndR();
//    ---------------

      checkRingEff();
//    --------------

      checkOptCorr();
//    --------------

      //singlePhotonCAT();
//    -----------------


      return;
  }


//===========================================================================
  void CsRCEventAnalysis::cluStructure( CsRCRing* ring ) {
//--------------------------------------------------------


//- Paolo  -  October 2001

    CsRCHistos& hist = CsRCHistos::Ref();
    double xh, yh, wh;

    list<CsRCPhoton*> lPhotons = ring->lPhotons();
    list<CsRCPhoton*>::iterator ih;
    for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
      int nPadClu = (*ih)->ptToClu()->lPads().size();
      list<CsRCPad*> lPads = (*ih)->ptToClu()->lPads();
      list<CsRCPad*>::iterator id;
      //float padPHmx = 0.;
      float padPHmx = -1.;                          //   020823
      CsRCPad* padmx = NULL;
      for( id=lPads.begin(); id!=lPads.end(); id++ ) {
	if( (*id)->PH() > padPHmx ) {
	  padPHmx = (*id)->PH();
	  padmx = (*id);
	}
	xh = (*ih)->phi();
	wh = (*id)->PH();
	if( hist.hRC3252 ) hist.hRC3252->Fill( xh, wh );
//                         ----------------------------
	float ddxx = (*id)->ix() - lPads.front()->ix();
	float ddyy = (*id)->iy() - lPads.front()->iy();
	xh = ddxx;
	yh = ddyy;
	if( hist.hRC3255 ) hist.hRC3255->Fill( xh, yh );
//                         ----------------------------
	wh = (*id)->PH();
	if( hist.hRC3256 ) hist.hRC3256->Fill( xh, yh, wh );
//                         --------------------------------
      }
      if( padmx != lPads.front() ) continue;
      if( nPadClu > 1 ) {
	//cout << nPadClu << "  " << (*ih)->PH() << "  " << padPHmx 
	//     << "  " << ((*ih)->PH()-padPHmx)/padPHmx << endl;
      }

      xh = nPadClu;
      if( hist.hRC3250 ) hist.hRC3250->Fill( xh );
//hh                     ------------------------
    }

  }


//===========================================================================
  void CsRCEventAnalysis::PIDMonitor() {
//--------------------------------------


//--- particle identification monitor
//    -------------------------------
//--- Paolo  -  December 1999,  rev. June 2000,  
//                              rev. October 2000
//                              rev. January 2001
//                              rev. June 2001
//                              rev. December 2002


      CsRCExeKeys* key = CsRCExeKeys::Instance();

      static bool firstCall = true;
      if( firstCall ) { 
        firstCall = false;

        key->acknoMethod( "CsRCEventAnalysis::PIDMonitor" );
      }

      CsRCRecConst* cons = CsRCRecConst::Instance();
      int kOffHi = cons->kOffHi();
      float CFRefInd = cons->CFRefInd();
      double TwoPI = cons->TwoPI();
      double* massPartv = cons->massPartv();
      unsigned int nPhotMinRing = cons->nPhotMinRing();

      static double sq2pi = sqrt( TwoPI );

      CsRCHistos& hist = CsRCHistos::Ref();
      double xh, yh, wh;

      double likeDV = cons->likeDefVa();
//@@-----------------------------------


//--- PARTICLE IDENTIFICATION monitor :
//    ---------------------------------

//--- loop over rings :
//    -----------------
      CsRCEventRings* evrings = CsRCEventRings::Instance();
      list<CsRCRing*> lRings = evrings->lRings();
      list<CsRCRing*>::iterator ir;
      for( ir=lRings.begin(); ir!=lRings.end(); ir++ ) {
        if( !(*ir)->flag() ) continue;

//----- monitor input rings (particles) (MC only) :
//      -------------------------------------------
        if( key->MCarloEvent() ) {
	  if( (*ir)->pPartPhot()->pPart()->flag() ) {
	    int iPartT = (*ir)->pPartPhot()->pPart()->iPartT();
	    if( iPartT > 200) iPartT -= 200;
	    if( iPartT > 30 ) iPartT = 0;               //   wrong myFile?
	    if( hist.hRC3615 ) hist.hRC3615->Fill( double( iPartT ), 3.5 );
//hh                           -------------------------------------------
	  }
	}
      }


      double recoMass;
      double probMass;
      int nPIDEv = 0;

//--- LOOP over PART-PHOTONS :
//    ------------------------
      CsRCEventPartPhotons* evepapho = CsRCEventPartPhotons::Instance();
      list<CsRCPartPhotons*> lPaPhotons = evepapho->lPartPhotons();
      list<CsRCPartPhotons*>::iterator ipf;

      for( ipf=lPaPhotons.begin(); ipf!=lPaPhotons.end(); ipf++ ) {
	if( !(*ipf) ) continue;
	CsRCPartPhotons* papho = (*ipf);
	if( !papho->flag() ) continue;

	CsRCRing* ring = papho->pRing();
	bool RING = false;
	if( ring != NULL  &&  
	    ring->flag()  &&
	    ring->lPhotons().size() >= nPhotMinRing ) RING = true;
//@@    ---------------------------------------------------------

	double* probaLK;
	double* derivLK;
	double* qSquare;
	double qSquaR = 0.;
	int nPhoRing = 0;
	double lZero[31];
	for( int k=0; k<31; k++ ) lZero[k] = 0.;
	probaLK = lZero;
	derivLK = lZero;
	qSquare = lZero;
	double likePion = likeDV;
	double likeKaon = likeDV;
	double likeProton = likeDV;
	double likeBkgr = likeDV;
	bool PIDset = false;

	double theRecoR = 0.;
	double theReco = 0.;
	double theRecoWg = 0.;
	double theRecoFt = 0.;

	if( RING ) {
//------- reconstructed Cerenkov angle :
//        ------------------------------
	  theRecoR = ring->the();
	  theReco = ring->theReco();
	  theRecoWg = ring->getThetaWgav();
	  theRecoFt = ring->thetaRFit();

//------- from QSQUARE :
//        --------------
	  nPhoRing = ring->nPhotQsQ();
	  qSquaR = ring->ringQsQ();
	  qSquare = ring->qSquareRing();
	}

	double theRecoLkR = 0.;
//----- PART-PHOTONS  particle identification :
//      ---------------------------------------
	if( cons->likeType() == "ALL" ) {

	  if( papho->thetaLikeSet() ) theRecoLkR = papho->thetaLike();
	  probaLK = papho->probaLKAll();
	  derivLK = papho->derivLKAll();
	  likeBkgr = papho->probaLKBgAll();
	  likePion = probaLK[ 8];
	  likeKaon = probaLK[11];
	  likeProton = probaLK[14];

	  PIDset = true;
	}

//----- RING  particle identification :
//      -------------------------------
	else if( cons->likeType() == "RING"  &&  RING ) {

	  if( ring->thetaLikeSet() ) theRecoLkR = ring->thetaLike();
	  probaLK = ring->probaLKRing();
	  derivLK = ring->derivLKRing();
	  likeBkgr = ring->probaLKBgRing();
	  likePion = probaLK[ 8];
	  likeKaon = probaLK[11];
	  likeProton = probaLK[14];

	  PIDset = true;
	}

//----- CHISQ  particle identification :
//      --------------------------------
	else if( cons->likeType() == "CHISQ"  &&  RING ) {

	  PIDset = true;
	}

	else  continue;
	if( !PIDset ) continue;


//----- define SELECTED part-photons or rings :
//      ---------------------------------------
        bool bSele = true;

//----- //particle momentum over muon threshold (1.91 GeV/c) :
//      ----------------------------------------------------
	CsRCParticle* part = papho->pPart();
        float momMeas = part->mom();
        //if( momMeas < 1.91 ) bSele = false;
//@@--------------------------------------

	if( RING ) BMDisplay( ring );
//@@--------------------------------


//----- process selected part-photons or rings :
//      ----------------------------------------
        if( bSele ) {

//------- monitor selected rings (particles) (MC) :
//        -----------------------------------------
          if( key->MCarloEvent() ) {
	    if( papho->pPart()->flag() ) {
	      int iPartT = papho->pPart()->iPartT();
	      if( iPartT > 200) iPartT -= 200;
	      if( iPartT > 30 ) iPartT = 0;   //   wrong myFile?
	      if( hist.hRC3615 ) hist.hRC3615->Fill( double( iPartT ), 4.5 );
//hh                             -------------------------------------------
	    }
	  }

          double betaReco = 1./( cos( theReco/1000. ) * CFRefInd );
          float gammaReco;
          if( betaReco > 1.) betaReco = 1.;
          if( betaReco < 1.) {
            gammaReco = 1./ sqrt(1.- betaReco*betaReco);
	  } else  { gammaReco = 0.; }

//------- correct reconstructed theta for E losses : (currently not used)
//        ------------------------------------------
	  //double theRecoCo = theReco - offTheRing( betaReco );
//@@---------------------------------------------------------

	  int kTyLL = 2;
//@@-------------------
          int kFigMe = 0;
          double probaQS[31];
          bool kOvThre[31];
          double fgMeritMx = 10000000.;
          double fgMeritMn = fgMeritMx;
          int kTyMin = 30;
          float recCut = 0.;
//------- look for electron, muon, pion, kaon or proton only :
//        ----------------------------------------------------
          for( int kPaTy=kTyLL; kPaTy<=14; kPaTy+=3 ) {
//        ---------------------------------------------

	    kOvThre[kPaTy] = false;
	    kOvThre[kPaTy+1] = false;
	    probaQS[kPaTy] = -1.;

            float massIpo = part->mass( kPaTy );
            double betaIpo = momMeas /
              sqrt( momMeas*momMeas + massIpo*massIpo );
	    //betaIpo *= 1.0000;   //   mom res at thresh(?)
	    //remember at thresh -> 0 phot <=> backgr ...
//@@--------------------------------------------------------

//--------- check particle threshold :
//          --------------------------
	    double theIpo = papho->thetaIpo( kPaTy );
	    if( theIpo > 0. ) {

              kOvThre[kPaTy] = true;
              kOvThre[kPaTy+1] = true;

	      double pLike = probaLK[kPaTy];

	      double qSqua = 0.;
	      int nPhot = 0;
	      if( RING ) {
	        qSqua = qSquare[kPaTy];
	        nPhot = nPhoRing;
		probaQS[kPaTy] = 0.;
		if( qSqua > 0. ) probaQS[kPaTy] = prob_( qSqua, nPhot );
	      }

//----------- choose figure-of-merit type for PID :
//            -------------------------------------
//            kFigMe = 1: chi-square;  kFigMe = 2: likelihood :
//            -------------------------------------------------
	      if( cons->likeType() == "CHISQ"  &&  RING ) kFigMe = 1;
//@@----------------------------------------------------------------
	      else if( cons->likeType() == "RING"  &&  RING ) kFigMe = 2;
//@@--------------------------------------------------------------------
	      else if( cons->likeType() == "ALL" ) kFigMe = 2;
//@@---------------------------------------------------------
	      else  kFigMe = 0;

//----------- define selected minimum (chisq) = maximum (like) :
//            --------------------------------------------------
              double fgMerit = 0;
	      if( kFigMe == 1 ) fgMerit = qSqua;
              if( kFigMe == 2 ) fgMerit = 1./ pLike;
              if( fgMerit < fgMeritMn ) {
                fgMeritMn = fgMerit;
                kTyMin = kPaTy;
//@@            --------------
                if( RING ) recCut = ring->recoCut( betaIpo );
//                                  ------------------------
	      }

            }   // end if on theIpo

          }   // end loop on particle type: kPaTy

	  if( kTyMin < 30 ) {
	    if( cons->likeType() == "ALL" ) nPIDEv++;
	    else if( cons->likeType() == "RING"  &&  RING ) nPIDEv++;
	    else if( cons->likeType() == "CHISQ"  &&  RING ) nPIDEv++;
	  } else {
	    continue;
	  }


//------- LOCAL ANALYSIS :
//        ----------------
//------- assign mass and probability :
//        -----------------------------
	  float probCut = 0.;                   //   no cut
//@@------------------------
          float chiCut = 1000.;                 //   no cut
//@@--------------------------

          recoMass = -1.;
          probMass = -1.;
          if( kTyMin >= kTyLL  &&  kTyMin <= 14 ) {
//        -----------------------------------------
            if( kOvThre[kTyMin] ) {
              recoMass = massPartv[kTyMin];
              if( kFigMe == 1 ) probMass = probaQS[kTyMin];
              if( kFigMe == 2 ) probMass = probaLK[kTyMin];
	    }
          }

	  if( RING ) {
//--------- theta resolution histos :
//          -------------------------
            float massIde = massPartv[kTyMin];
            double betaIde = momMeas /
              sqrt( momMeas*momMeas + massIde*massIde );
 
	    double theIde = papho->thetaIpo( kTyMin );
	    //double theIde = papho->thetaIpoVS( kTyMin );
	    if( theIde > 0. ) {
	    if( kTyMin >= kTyLL  &&  kTyMin <= 14 ) {
//          -----------------------------------------
	      double dTheta = theReco - theIde;
              //double dTheta = getThetaWgav( (*ir) ) - theIde;

              xh = dTheta;
              yh = (kTyMin-kTyLL)/3 + 0.5;
//                 ----------------
              if( hist.hRC1635 ) hist.hRC1635->Fill( xh, yh );
//hh                             ----------------------------
              if( betaIde >= 0.99995 ) {
                xh = dTheta;
                if( hist.hRC1636 ) hist.hRC1636->Fill( xh );
//hh                               ------------------------
	      }
	      if( key->MCarloEvent() ) {
		if( papho->howPMT() < 0.2 ) {
		  xh = part->thetaCer() - papho->thetaIpo( part->iPartT() );
		  yh = papho->thetaIpo( part->iPartT() );
		  if( hist.hRC1634 ) hist.hRC1634->Fill( xh, yh );
//hh                                 ----------------------------
		}
		if( papho->howPMT() > 0.8 ) {
		  xh = part->thetaCer() - papho->thetaIpoVS( part->iPartT() );
		  yh = papho->thetaIpoVS( part->iPartT() );
		  if( hist.hRC1637 ) hist.hRC1637->Fill( xh, yh );
//hh                                 ----------------------------
		}
	      }
              xh = dTheta;
	      yh = theIde;
              int kHHy = (kTyMin-kTyLL)/3;
//                       ----------------
              if( kHHy >= 0  &&  kHHy < 5 ) {
		if( hist.vRC1680[kHHy] ) hist.vRC1680[kHHy]->Fill( xh, yh );
//hh                                     ----------------------------------
	      }

	      if( ring->flagBack() ) {                 //   ???
		double corrRing = ring->ringBack();
//                                ----------------
		int nPhoRingCo = int( nPhoRing * corrRing );
		int nPhoExp = int( ring->getNPhotExpected( theIde ) );
//                                 --------------------------------
		xh = nPhoRingCo - nPhoExp;
		yh = (kTyMin-kTyLL)/3 + 0.5;
//                   ----------------
		if( hist.hRC1638 ) hist.hRC1638-> Fill( xh, yh );
//hh                               -----------------------------
//------------- OBSO!
		double nPhoRingN0 = nPhoRing * corrRing;
		//xh = nPhoRingN0;
		//if( hist.hRC1634 ) hist.hRC1634-> Fill( xh );
//hh                               -------------------------

		double pathLen = papho->pPart()->pathLen();
	        double nZero = nPhoRingN0 /
		  ( pathLen/10.* (theIde/1000.)*(theIde/1000.) );
		xh = nZero;
		if( hist.hRC1639 ) hist.hRC1639-> Fill( xh );
//hh                               -------------------------
	      }
	    }
	    }
	  }

	  if( RING ) {
//--------- chi-square histos :
//          -------------------
            if( kTyMin >= 5 && kTyMin <= 14 ) {
//          ---------------------------------
              if( qSquare[kTyMin] < fgMeritMx ) {

		int kChi = (kTyMin-5)/3;
                xh = qSquare[kTyMin];
		yh = float( kChi ) + 0.5;
          	if( xh > 0. ) 
		  if( hist.hRC1661 ) hist.hRC1661-> Fill( xh, yh );
//hh                                 -----------------------------

		xh = probaQS[kTyMin];
		if( hist.hRC1666 ) hist.hRC1666-> Fill( xh );
//hh                               -------------------------

                for( int kPaTy=8; kPaTy<=14; kPaTy+=3 ) {
//              ---------------------------------------
                  if( kPaTy != kTyMin )  {
                    xh = qSquare[kPaTy];
                    if( hist.hRC1662 ) hist.hRC1662->Fill( xh );
//hh                                   ------------------------
                    xh = fabs( qSquare[kPaTy] - qSquare[kTyMin] );
                    yh = qSquare[kTyMin];
                    if( hist.hRC1663 ) hist.hRC1663->Fill( xh, yh );
//hh                                   ----------------------------
                    xh = probaQS[kPaTy];
                    if( hist.hRC1667 ) hist.hRC1667->Fill( xh );
//hh                                   ------------------------
                    xh = fabs( probaQS[kPaTy] - probaQS[kTyMin] );
                    if( hist.hRC1668 ) hist.hRC1668->Fill( xh );
//hh                                   ------------------------
                  }   // end of if on kTyMin
                }   // end of loop on kPaTy

              }   // end of if on fgMeritMx
            }   // end of if on particle type: kTyMin
	  }


//------- likelihood histos :
//        -------------------
	  double partProbMx = 0.;
          for( int kPaTy=kTyLL; kPaTy<=14; kPaTy+=3 ) {
//        -------------------------------------------
	    if( probaLK[kPaTy] > partProbMx ) partProbMx = probaLK[kPaTy];
	  }
	  if( likeBkgr >= partProbMx ) {
	    xh = partProbMx;
	    yh = momMeas;
	    if( hist.hRC1691 ) hist.hRC1691->Fill( xh, yh );
//hh                           ----------------------------
	  }

          xh = momMeas;
	  yh = theReco;
          if( hist.hRC1640 ) hist.hRC1640->Fill( xh, yh );
//hh                         ----------------------------
	  double probRef = likeBkgr;
	  double prPion = likePion - probRef;
	  double prKaon = likeKaon - probRef;
	  double prProton = likeProton - probRef;
	  if( prPion >= 0. && prPion > prKaon && prPion > prProton ) {
	    xh = momMeas;
	    yh = likePion / likeBkgr;
	    if( hist.hRC1641 ) hist.hRC1641->Fill( xh, yh );
//hh                           ----------------------------
	  }
	  if( prKaon >= 0. && prKaon > prPion && prKaon > prProton ) {
	    xh = momMeas;
	    yh = likeKaon / likeBkgr;
	    if( hist.hRC1642 ) hist.hRC1642->Fill( xh, yh );
//hh                           ----------------------------
	  }
	  if( prProton >= 0. && prProton > prPion && prProton > prKaon ) {
	    xh = momMeas;
	    yh = likeProton / likeBkgr;
	    if( hist.hRC1643 ) hist.hRC1643->Fill( xh, yh );
//hh                           ----------------------------
	  }


//------- likelihood histos (pions, kaons and protons only):
//        --------------------------------------------------
          if( kTyMin >= 8  &&  kTyMin <= 14 ) {
//        -----------------------------------
	    double probaLKMx = probaLK[kTyMin];
	    double diffProMin = 1000000.;
	    int kPaTy2 = 0;
            for( int kPaTy=8; kPaTy<=14; kPaTy+=3 ) {
//          ---------------------------------------
              if( kPaTy == kTyMin ) continue;
	      if( probaLK[kPaTy] < 0. ) continue;
	      double diffPro = fabs( probaLKMx - probaLK[kPaTy] );
	      if( diffPro < diffProMin ) {
	        diffProMin = diffPro;
		kPaTy2 = kPaTy;
	      }
	    }
	    if( probaLKMx > 0.  &&  kPaTy2 > 0 ) {
              xh = probaLK[kPaTy2] / probaLKMx;
	      yh = float( (kTyMin-5)/3 ) + 0.5;
//                 ---------------------
	      if( hist.hRC1664 ) hist.hRC1664->Fill( xh, yh );
//hh                             ----------------------------
            }
            xh = probaLKMx;
	    yh = probaLK[kPaTy2];
	    if( hist.hRC1665 ) hist.hRC1665->Fill( xh, yh );
//hh                           ----------------------------
	    xh = probaLK[kTyMin];
            yh = float( (kTyMin-5)/3 ) + 0.5;
//               ---------------------
            if( hist.hRC1669 ) hist.hRC1669->Fill( xh, yh );
//hh                           ----------------------------
	    xh = likeBkgr;
            yh = 0.5;
            if( hist.hRC1669 ) hist.hRC1669->Fill( xh, yh );
//hh                           ----------------------------
	  }

  	  //checkLikeDistr();
//        ----------------


//------- reconstructed mass spectrum histos :
//        ------------------------------------

	  bool bReco = false;
	  for( int kPaTy=kTyLL; kPaTy<=14; kPaTy+=3 ) {
//        -------------------------------------------
	    float massIde = massPartv[kPaTy];
            double betaIde = momMeas /
              sqrt( momMeas*momMeas + massIde*massIde );
	    double theIde = papho->thetaIpo( kPaTy );
	    if( theIde < 0.) continue;
	    if( fabs( theReco - theIde ) < ring->recoCut( betaIde ) ) {
	      bReco = true;
	      break;
	    }
	  }
	  double cosReco = cos( theReco/1000. );
	  double measMass = CFRefInd*cosReco*CFRefInd*cosReco - 1.;
	  if( measMass < 0. ) measMass = 0.;
	  measMass = sqrt( measMass );
	  measMass *= momMeas;

	  //if( bReco ) {
	    xh = measMass;
	    yh = momMeas;
	    if( hist.hRC1692 ) hist.hRC1692->Fill( xh, yh );
//hh                           ----------------------------
	  //}
	  int kProbMx = -1;
	  int kChiMn = -1;
	  if( kFigMe == 1 ) {
	    double partChiMn = 1000000.;
	    int kChim = 0;
	    for( int kPaTy=8; kPaTy<=14; kPaTy+=3 ) {
//          ---------------------------------------
	      if( qSquare[kPaTy] < partChiMn ) {
		partChiMn = qSquare[kPaTy];
		kChim = kPaTy;
	      }
	    }
	    int kChiMn = (kChim-5)/3;
//                       -----------
	    if( partChiMn < chiCut ) {
	      xh = measMass;
	      yh = kChiMn + 0.5;
	      if( hist.hRC1693 ) hist.hRC1693->Fill( xh, yh );
//hh                             ----------------------------
	    }
	    kProbMx = kChiMn;
	  }

	  if( kFigMe == 2 ) {
	    double partProbMx = 0.;
	    kProbMx = -1;
	    for( int kPaTy=8; kPaTy<=14; kPaTy+=3 ) {
//          ---------------------------------------
	      if( probaLK[kPaTy] > partProbMx ) {
		partProbMx = probaLK[kPaTy];
		kProbMx = kPaTy;
	      }
	    }
	    if( kProbMx >= 0 ) {
	      //if( bReco ) {
	      //if( measMass > 0.440  &&  measMass < 0.560 ) {
	      xh = measMass;
	      yh = probaLK[kProbMx]/likeBkgr;
	      int kHistMx = (kProbMx-5)/3 + 1;
//                          -----------------
	      if( kHistMx > 0  &&  kHistMx < 5 ) {
	        if( hist.vRC1590[kHistMx] )
		  hist.vRC1590[kHistMx]->Fill( xh, yh );
//hh              -------------------------------------
	      }
	      //}
	      if( probaLK[kProbMx]/likeBkgr > probCut ) {
	        xh = measMass;
	        yh = (kProbMx-5)/3 + 0.5;
//                   -------------
	        if( hist.hRC1693 ) hist.hRC1693->Fill( xh, yh );
//hh                               ----------------------------
	      }
	      //}
	    }
	  }


//------- momentum resolution from RICH (data) :
//        --------------------------------------
	  if( RING ) {
	    int kHistMx = (kProbMx-5)/3;
//                        -------------
	    if( kHistMx > 0  &&  kHistMx < 4 ) {
	      if( nPhoRing >= 10 ) {
	        int iPartTw = kProbMx;
	        float momReco = gammaReco*betaReco * part->mass( iPartTw );
		xh = (momReco - momMeas) / momMeas;
		yh = momMeas;
		int khw = kHistMx - 1;
	        if( khw >= 0 ) {
		  if( hist.vRC1655[khw] ) hist.vRC1655[khw]->Fill( xh, yh );
//hh                                      ---------------------------------
	        }
	      }
	    }
	  }


//------- MONTECARLO analysis :
//        ---------------------
	  if( RING  &&  key->MCarloEvent() ) {

            int iPartTw = part->iPartT();
            if( iPartTw > 200 ) iPartTw -= 200;
	    if( iPartTw > 30 ) iPartTw = 0;               //   wrong myFile?
//@@--------------------------------------

//--------- mass spectrum (NO muons):
//          -------------------------
            if( !(iPartTw ==  5  || iPartTw ==  6 ) ) {
	      xh = measMass;
	      yh = theReco;
	      if( hist.hRC1694 ) hist.hRC1694->Fill( xh, yh );
//hh                             ----------------------------
	      xh = measMass;
      	      //yh = kProbMx + 0.5;
	      int kHistMx = (kProbMx-5)/3;
//                          -------------
      	      yh = kHistMx + 0.5;
              if( hist.hRC1695 ) hist.hRC1695->Fill( xh, yh );
//hh                             ----------------------------
	      xh = measMass;
	      yh = theReco;
	      if( hist.hRC1696 ) hist.hRC1696->Fill( xh, yh );
//hh                             ----------------------------
	    }

//--------- pions, kaons or protons :
//          -------------------------
            int khw = -1;
            if( iPartTw ==  8  || iPartTw ==  9 ) khw = 0;
            if( iPartTw == 11  || iPartTw == 12 ) khw = 1;
            if( iPartTw == 14  || iPartTw == 15 ) khw = 2;
	    // !!!           xh = recoMass;
	    xh = measMass;
            if( khw >= 0 ) { 
	      if( hist.vRC1650[khw] ) hist.vRC1650[khw]->Fill( xh );
//hh                                  -----------------------------
	    }

//--------- PID table :
//          -----------
            //if( probMass >= 0. ) {                //   030326   !!!

//----------- monitor ID rings (particles) (MC) :
//            -----------------------------------
	      if( hist.hRC3615 ) hist.hRC3615->Fill( double( iPartTw ), 5.5 );
//hh                             --------------------------------------------
              xh = 0.;
              yh = 0.;
              if( iPartTw ==  2  || iPartTw ==  3 ) xh = 0.5;
	      else if( iPartTw ==  5  || iPartTw ==  6 ) xh = 1.5;
              else if( iPartTw ==  8  || iPartTw ==  9 ) xh = 2.5;
              else if( iPartTw == 11  || iPartTw == 12 ) xh = 3.5;
              else if( iPartTw == 14  || iPartTw == 15 ) xh = 4.5;
	      else xh = 5.5;
	      if( kOvThre[iPartTw] ) {
                if( recoMass == part->mass(  2 ) ) yh = 0.5;
                else if( recoMass == part->mass(  5 ) ) yh = 1.5;
	        else if( recoMass == part->mass(  8 ) ) yh = 2.5;
                else if( recoMass == part->mass( 11 ) ) yh = 3.5;
                else if( recoMass == part->mass( 14 ) ) yh = 4.5;
	      }
	      else if( !kOvThre[iPartTw] ) {
		if( recoMass == part->mass(  2 ) ) yh = 5.5;
                else if( recoMass == part->mass(  5 ) ) yh = 6.5;
	        else if( recoMass == part->mass(  8 ) ) yh = 7.5;
                else if( recoMass == part->mass( 11 ) ) yh = 8.5;
                else if( recoMass == part->mass( 14 ) ) yh = 9.5;
	      }
	      //else yh = 6.5;
              if( xh > 0. &&  yh > 0. ) {
	        if( hist.hRC1655 ) hist.hRC1655->Fill( xh, yh );
//Hh                               ----------------------------
	      }
              int kTab = int( momMeas / 5. );
              //if( kTab > 11 ) kTab = 11;
//----------- limit to 10 slices to avoid hist owerwriting
	      if( kTab > 9 ) kTab = 9;
              if( xh > 0. &&  yh > 0. ) { 
                if( hist.vRC1670[kTab]) hist.vRC1670[kTab]->Fill( xh, yh ); 
//hh                                    ----------------------------------
              }
	    //}   // end if on probMass

//--------- momentum resolution from RICH :    (moved to data 020916)
//          -------------------------------
	    //float momReco = gammaReco*betaReco * part->mass( iPartTw );
	    //xh = (momReco - momMeas) / momMeas;
            //yh = momMeas;
            //khw = -1;
            //if( iPartTw ==  8  || iPartTw ==  9 ) khw = 0;
            //if( iPartTw == 11  || iPartTw == 12 ) khw = 1;
            //if( iPartTw == 14  || iPartTw == 15 ) khw = 2;
            //if( khw >= 0 ) {
	      //if( hist.vRC1655[khw] ) hist.vRC1655[khw]->Fill( xh, yh );
//hh                                  ---------------------------------
	    //}

//--------- likelihood particle identification efficiency :
//          -----------------------------------------------
            //if( kOvThre[iPartTw] ) {                //   030326   !!!
              khw = -1;
              if( iPartTw ==  8  || iPartTw ==  9 ) khw = 0;
              if( iPartTw == 11  || iPartTw == 12 ) khw = 1;
              if( iPartTw == 14  || iPartTw == 15 ) khw = 2;
              xh = momMeas;
              if( khw == 0 ) if( hist.hRC3072 ) hist.hRC3072->Fill( xh );
              if( khw == 1 ) if( hist.hRC3073 ) hist.hRC3073->Fill( xh );
              if( khw == 2 ) if( hist.hRC3074 ) hist.hRC3074->Fill( xh );
//hh                                            ------------------------
              xh = nPhoRing;
              if( khw >= 0 ) if( hist.hRC3075 ) hist.hRC3075->Fill( xh );
//hh                                            ------------------------
              double partProbMx = 0.;
              int kProbMx = -1;
	      int iProb = 0;
	      int iProbMx = -1;
	      //for( int kPaTy=8; kPaTy<=14; kPaTy+=3 ) {     //   !!!
	      for( int kPaTy=kTyLL; kPaTy<=14; kPaTy+=3 ) {
//            ---------------------------------------------
		if( probaLK[kPaTy] > partProbMx ) {
		  partProbMx = probaLK[kPaTy];
		  kProbMx = kPaTy;
		  iProbMx = iProb;
		}
		iProb++;
	      }
	      iProbMx -= 2;                           //   !!!
              if( kProbMx >= 0 && iProbMx == khw ) {
                xh = momMeas;
                if( khw == 0 ) if( hist.hRC3076 ) hist.hRC3076->Fill( xh );
                if( khw == 1 ) if( hist.hRC3077 ) hist.hRC3077->Fill( xh );
                if( khw == 2 ) if( hist.hRC3078 ) hist.hRC3078->Fill( xh );
//hh                                              ------------------------
                xh = nPhoRing;
                if( khw >= 0 ) if( hist.hRC3079 ) hist.hRC3079->Fill( xh );
//hh                                              ------------------------
	      }
	      //}

          }   // end if on MCarloEvent

        }   // end if on bSele


      }   // end loop on part-photons

      checkLikeDistr();
//    ----------------

      checkThetaLikeMax();
//    -------------------

      checkLikeDisTheta();
//    -------------------


//--- moved to CsRCEventAnalysis::doEveAnalysis 040513
      //setProbsToTrack();
//@@-------------------

      xh = nPIDEv;
      if( hist.hRC1539 ) hist.hRC1539->Fill( xh );
//hh                     ------------------------


  }


//===========================================================================
  void CsRCEventAnalysis::chkHitLists( CsRCRing* ring, int &nRHt, int &nEqu,
//--------------------------------------------------------------------------
				       list<CsMCHit*>& lTrackHits,
				       list<CsMCHit*>& lRingHits ) {


//--- Paolo  -  January 2001

      list<CsRCPhoton*> lPhotons = ring->lPhotons();
      list<CsRCPhoton*>::iterator ih;
      int nPhoRing = lPhotons.size();
 
      lRingHits.clear();
      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	list<CsRCPad*> lPads = (*ih)->ptToClu()->lPads();
	list<CsRCPad*>::iterator ia;
	for( ia=lPads.begin(); ia!=lPads.end(); ia++ ) {
	  CsMCDigit* dig =
	    dynamic_cast<CsMCDigit*>((*ia)->pDigit());
	  list<CsMCHit*> lHits = dig->getHits();
	  list<CsMCHit*>::iterator ii;
	  for( ii=lHits.begin(); ii!=lHits.end(); ii++ ) {
	    lRingHits.push_back( (*ii) );
	  }
	}
      }
      nRHt = lRingHits.size();
      lRingHits.sort();
      lRingHits.unique();

      lTrackHits.clear();
      CsRICH1Detector* rich =
	CsGeom::Instance()->getRich1Detector();
      list<CsMCHit*> lMChits =
	ring->pPartPhot()->pPart()->pMCTrack()->getMCHits();
      list<CsMCHit*>::iterator ii;
      for( ii=lMChits.begin(); ii!=lMChits.end(); ii++ ) {
	if( (*ii)->getOrigin() != 0 ) continue;
	CsDetector* det =
          dynamic_cast<CsDetector*>((*ii)->getDet());
        if( det == dynamic_cast<CsDetector*>( rich ) ) {
          lTrackHits.push_back( (*ii) );
        }
      }
      lTrackHits.sort();
      lTrackHits.unique();
      nEqu = 0;
      for( ii=lTrackHits.begin(); ii!=lTrackHits.end(); ii++ ) {
        list<CsMCHit*>::iterator ik;
        for( ik=lRingHits.begin(); ik!=lRingHits.end(); ik++ ) {
          if( (*ik) == (*ii) ) { nEqu++; break; }
        }
      }

  }


//===========================================================================
  void CsRCEventAnalysis::chkHitPrint( CsRCRing* ring, int &nRHt, int &nEqu,
//--------------------------------------------------------------------------
				       list<CsMCHit*>& lTrackHits,
				       list<CsMCHit*>& lRingHits ) {

//--- Paolo  -  January 2001

      list<CsRCPhoton*> lPhotons = ring->lPhotons();
      list<CsRCPhoton*>::iterator ih;
      int nPhoRing = lPhotons.size();
      list<CsMCHit*>::iterator ii;

      cout << endl << "=========   : " << endl;
      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	list<CsRCPad*> lPads = (*ih)->ptToClu()->lPads();
	cout << lPads.size() << " (";
	list<CsRCPad*>::iterator ia;
	for( ia=lPads.begin(); ia!=lPads.end(); ia++ ) {
	  CsMCDigit* dig = 
	    dynamic_cast<CsMCDigit*>((*ia)->pDigit());
	  list<CsMCHit*> lHits = dig->getHits();
	  cout << lHits.size() << " ";
	  list<CsMCHit*>::iterator ii;
	  for( ii=lHits.begin(); ii!=lHits.end(); ii++ ) {
	    cout << (*ii) << " ";
	  }
	}
	cout << ") ";
      }
      cout << endl;
      for( ii=lRingHits.begin(); ii!=lRingHits.end(); ii++ ) {
	cout << (*ii) << " ";
      }
      cout << endl;
      for( ii=lTrackHits.begin(); ii!=lTrackHits.end(); ii++ ) {
	cout << (*ii) << " ";
      }
      cout << endl;
      cout << ""
	   << nPhoRing << "   "
	   << lRingHits.size() << "(" << nEqu << "=)" 
           << "(" << nRHt << ")   "
	   << lTrackHits.size() << "   "
           << endl;
  }


//===========================================================================
  void CsRCEventAnalysis::hitDisplay( int histN, list<CsMCHit*> lHits ) {
//-----------------------------------------------------------------------


//--- data monitoring and ring analysis.
//    ----------------------------------
//--- Paolo  -  November 2000


      CsRCRecConst *cons = CsRCRecConst::Instance();
      int kOffHi = cons->kOffHi();

      CsRCDetectors *dets = CsRCDetectors::Instance();

      static vector<CsHist2D*> vRC5000;
      static int kHH = 0;
      static int kHHMx = 20;
      if( !CsRCHistos::Ref().bookHis() ||
	  CsRCHistos::Ref().levelBk() < 2 ) kHHMx = 0;
      static int kH5000 = 0;
      static int kH5100 = 0;
      static int kH5200 = 0;

      static bool firstCall = true;
      if( firstCall ) {
        firstCall = false;
      }

      int kh;
      double xh, yh, wh;

      int kH = 0;
      if( kHH < kHHMx ) {
	if( histN == 5000 ) kH = kH5000;
	if( histN == 5100 ) kH = kH5100;
	if( histN == 5200 ) kH = kH5200;
        int khDisp = kOffHi + histN + kH;
        stringstream name0, title0;
        name0 << khDisp;
        title0 << "   ";
	CsHistograms::SetCurrentPath("/RICH");
        vRC5000.push_back( new CsHist2D( name0.str(), title0.str(),
//hh    -----------------------------------------------------------
        310, -1240., 1240., 310, -1240., 1240. ) );
        CsHistograms::SetCurrentPath("/");

	list<CsMCHit*>::iterator ii;
	for( ii=lHits.begin(); ii!=lHits.end(); ii++ ) {
	  CsMCRICH1Hit* hit = dynamic_cast<CsMCRICH1Hit*>( (*ii) );
	  int ic = hit->getCathode() - 1;
	  xh = hit->getYdetDRS() + dets->vOffCatW( ic ).x();
	  yh = hit->getZdetDRS() + dets->vOffCatW( ic ).y();
	  wh  = 1.;
	  vRC5000[kHH]->Fill( xh, yh, wh );
//hh      --------------------------------
	}
      }
      kHH++;
      if( histN == 5000 ) kH5000++;
      if( histN == 5100 ) kH5100++;
      if( histN == 5200 ) kH5200++;
  }


//===========================================================================
  void CsRCEventAnalysis::BMDisplay( CsRCRing* ring ) {
//-----------------------------------------------------


//--- data monitoring and ring analysis
//    ---------------------------------
//--- Paolo  -  March 2001


      CsRCRecConst *cons = CsRCRecConst::Instance();
      int kOffHi = cons->kOffHi();

      CsRCDetectors *dets = CsRCDetectors::Instance();

      static vector<CsHist2D*> vRC5300;
      static int kHH = 0;
      static int kHHMx = 20;
      if( !CsRCHistos::Ref().bookHis() ||
	  CsRCHistos::Ref().levelBk() < 2 ) kHHMx = 0;
      static int histN = 5300;
      static int kH5300 = 0;

      static bool firstCall = true;
      if( firstCall ) {
        firstCall = false;
      }

      int kh;
      double xh, yh, wh;

      int kH = 0;
      if( kHH < kHHMx ) {
	kH = kH5300;
        int khDisp = kOffHi + histN + kH;
        stringstream name0, title0;
        name0 << khDisp;
        title0 << "   ";
	CsHistograms::SetCurrentPath("/RICH");
        vRC5300.push_back( new CsHist2D( name0.str(), title0.str(),
//hh    -----------------------------------------------------------
        100, -10., 10., 60, 0., 360 ) );
        CsHistograms::SetCurrentPath("/");

	double thetaReco = ring->the();
        list<CsRCPhoton*> lPhotons = ring->lPhotons();
        list<CsRCPhoton*>::iterator ih;
        for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
          double thePhot = (*ih)->the();
          double phiPhot = (*ih)->phi();
	  double thetaB = (*ih)->theB();
	  double thetaM = (*ih)->theM();
	  double thetaBM = fabs( thetaB - thetaM );
	  double xx1 = thetaB - thetaReco;
	  double xx2 = thetaM - thetaReco;
	  //cout << xx1 << "  " << xx2 << endl;
	  int npt = int( thetaBM/0.15 + 1.);
	  double dBM = thetaBM / float( npt );
	  //cout << npt << "  " << dBM << endl;
	  for( int k=0; k<=npt; k++ ) {
	    xh = xx1 + k*dBM;
	    yh = phiPhot;
	    wh = 1.;
	    vRC5300[kHH]->Fill( xh, yh, wh );
//hh        --------------------------------
	  }
	  xh = thePhot - thetaReco;
	  yh = phiPhot;
          wh = 20.;
          vRC5300[kHH]->Fill( xh, yh, wh );
//hh      --------------------------------
	  xh = -9.;
	  yh = phiPhot;
	  wh = 10.;
	  vRC5300[kHH]->Fill( xh, yh, wh );
//hh      --------------------------------
	}
      }
      kHH++;
      kH5300++;
  }


//===========================================================================
  void CsRCEventAnalysis::moniCFRefInf( CsRCRing* ring ) {
//--------------------------------------------------------


//- Paolo  -  August 27, 2002
//  rev.      August     2008

    CsRCHistos& hist = CsRCHistos::Ref();
    float xh, yh;

    CsOpt* opt = CsOpt::Instance();
    CsRCRecConst *cons = CsRCRecConst::Instance();
    int kOffHi = cons->kOffHi();

    static std::vector<CsHist2D*> vRC6660;
    static std::vector<CsHist2D*> vRC6670;

//- defaults
//  --------
    static float momLimit  = 50.;
    static float tgxyLimit = 20.;
    static float nPhoLimit = 10.;

//- from rich1.options :
//  --------------------
    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;
      bool boo;
      vector<float> vPar;
      boo = opt->CsOpt::getOpt( "RICHONE", "n-1Limits", vPar );
      if( boo ) {
        momLimit  = vPar[0];
        tgxyLimit = vPar[1];
	nPhoLimit = vPar[2];
      }
//- Added 080822 ---
      for( int kh=0; kh<8; kh++ ) vRC6660.push_back( NULL );
      for( int kh=0; kh<5; kh++ ) vRC6670.push_back( NULL );
      if( CsRCHistos::Ref().bookHis() &&
	  CsRCHistos::Ref().levelBk() >= 2 ) {
	CsHistograms::SetCurrentPath("/RICH");
	vRC6660.clear();
        int kHist = 0;
        string hTitle = "n-1 from Ring";
	stringstream hN6661;
        kHist = kOffHi + 6661;
        hN6661 << kHist;
        vRC6660.push_back( new CsHist2D( hN6661.str(), hTitle,
					 200, 0.0002, 0.0022, 100, 0., 50.) );
	stringstream hN6662;
        kHist = kOffHi + 6662;
        hN6662 << kHist;
        vRC6660.push_back( new CsHist2D( hN6662.str(), hTitle,
					 200, 0.0002, 0.0022,
					 100, 0., 400.));
	stringstream hN6663;
        kHist = kOffHi + 6663;
        hN6663 << kHist;
        vRC6660.push_back( new CsHist2D( hN6663.str(), hTitle,
					 200, 0.0002, 0.0022,
					 100, -400., 400.) );
	stringstream hN6664;
        kHist = kOffHi + 6664;
        hN6664 << kHist;
        vRC6660.push_back( new CsHist2D( hN6664.str(), hTitle,
					 200, 0.0002, 0.0022,
					 100, -400., 400.) );
	stringstream hN6665;
        kHist = kOffHi + 6665;
        hN6665 << kHist;
        vRC6660.push_back( new CsHist2D( hN6665.str(), hTitle,
					 200, 0.0002, 0.0022,
					 100, -400., 400.) );
	stringstream hN6666;
        kHist = kOffHi + 6666;
        hN6666 << kHist;
        vRC6660.push_back( new CsHist2D( hN6666.str(), hTitle,
					 200, 0.0002, 0.0022,
					 100, 0., 50.) );
	stringstream hN6667;
        kHist = kOffHi + 6667;
        hN6667 << kHist;
        vRC6660.push_back( new CsHist2D( hN6667.str(), hTitle,
					 200, 0.0002, 0.0022,
					 90, 100., 400.));
	stringstream hN6668;
        kHist = kOffHi + 6668;
        hN6668 << kHist;
        vRC6660.push_back( new CsHist2D( hN6668.str(), hTitle,
					 200, 0.0002, 0.0022,
					 90, -300., 300.) );
	vRC6670.clear();
        for( int kh=0; kh<3; kh++ ) {
	  string hTitle = "n-1 from Ring - ID";
	  stringstream hN6670;
	  kHist = kOffHi + 6670 + kh + 1;
	  hN6670 << kHist;
	  vRC6670.push_back( new CsHist2D( hN6670.str(), hTitle,
					   200, 0.0002, 0.0022,
					   100, 0., 50.) );
	}
	CsHistograms::SetCurrentPath("/");
      }
    }

    bool bSele = true;

//- particle momentum :
    CsRCParticle* part = ring->pPartPhot()->pPart();
    float momPart = part->mom();
    //if( !( momPart > momLimit ) ) bSele = false;
    if( !( momPart < momLimit ) ) bSele = false;
//@@-------------------------------------------

//- particle angle at entrance window :
    double ttxy;
    ttxy = sqrt( pow( part->vDirIn().x()/part->vDirIn().z(), 2 ) +
  		 pow( part->vDirIn().y()/part->vDirIn().z(), 2 ) );
    ttxy *= 1000.;
    if( !( ttxy > tgxyLimit ) ) bSele = false;
//@@-----------------------------------------

//- photons in the ring :
    int nPhoRing = ring->lPhotons().size();
    if( !( nPhoRing >= int( nPhoLimit ) ) ) bSele = false;
//@@-----------------------------------------------------

    if( bSele ) {

      float theta = ring->the()/1000.;
      //double hthe2 = 0.5* theta*theta;
      //xh = hthe2;
//--- assume pion mass !
      float massIpo = CsRCRecConst::Instance()->massPartv()[8];
      double betaIpo = momPart / sqrt( momPart*momPart + massIpo*massIpo );
      double nMinus1 = 1./(betaIpo * cos( theta )) - 1.;
      xh = nMinus1;
      if( hist.hRC3059 ) hist.hRC3059->Fill( xh );
//                       ------------------------
    }


//- Added 080822 ---
    CsRCExeKeys *key = CsRCExeKeys::Instance();
    //if( CsRichOne::Instance()->UpRICHJob()  &&  key->selPMTonly() ) {

      bSele = true;

//    WARNING : using FIXED limits!
//--- 'momLimit' = lower limit
      double momLLimit = 5.;      // GeV/c
      if( momPart < momLLimit ) bSele = false;
//@@-----------------------------------------
//--- 'momLimit' = upper limit
      double momULimit = 50.;      // GeV/c
      if( momPart > momULimit ) bSele = false;
//@@-----------------------------------------
//--- 'tgxyLimit' = lower limit
      double tgxyLLimit = 20.;      // ~mrad
      if( ttxy < tgxyLLimit ) bSele = false;
//@@---------------------------------------
//--- 'nPhoLimit' = lower limit
      int nPhoLLimit = 10;
      //int nPhoLLimit =  5;
      if( int( ring->lPhotons().size() ) < nPhoLLimit ) bSele = false;
//@@-----------------------------------------------------------------

      if( bSele ) {
	CsRCPartPhotons* papho = ring->pPartPhot();
        int kDetPart = papho->kDetPart();
//----- particle angle to photon detector :
	double thPade = part->thPade()[kDetPart];
	thPade *= 1000.;
	double thPadeX = atan( papho->vDcPaReflW()[kDetPart].x() /
			       papho->vDcPaReflW()[kDetPart].z() );
	thPadeX *= 1000.;
	double tgx = part->vDirIn().x()/part->vDirIn().z();
	tgx *= 1000.;

	list<CsRCPhoton*> lPhotons = ring->lPhotons();
	list<CsRCPhoton*>::iterator ih;
	double theRing = 0.;
	for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	  theRing += (*ih)->the0();
	}
	if( lPhotons.size() > 0 ) theRing /= lPhotons.size();
	theRing /= 1000.;
//----- assume pion mass !
	float massIpo = CsRCRecConst::Instance()->massPartv()[8];
	double betaIpo = momPart / sqrt( momPart*momPart + massIpo*massIpo );
	double nMinus1 = 1./(betaIpo * cos( theRing )) - 1.;
	xh = nMinus1;
	if( hist.hRC3060 ) hist.hRC3060->Fill( xh );
//                         ------------------------

	int kHH;
	bool bHH;
	kHH = 0;
	bHH = kHH < int( vRC6660.size() );
	xh = nMinus1;
	yh = momPart;
	if( bHH  &&  vRC6660[kHH] ) vRC6660[kHH]->Fill( xh, yh );
//hh                                ----------------------------
	kHH = 1;
	bHH = kHH < int( vRC6660.size() );
	xh = nMinus1;
	yh = ttxy;
	if( bHH  &&  vRC6660[kHH] ) vRC6660[kHH]->Fill( xh, yh );
//hh                                ----------------------------
	kHH = 2;
	bHH = kHH < int( vRC6660.size() );
	xh = nMinus1;
	yh = tgx;
	if( bHH  &&  vRC6660[kHH] ) vRC6660[kHH]->Fill( xh, yh );
//hh                                ----------------------------
	kHH = 6;
	bHH = kHH < int( vRC6660.size() );
	xh = nMinus1;
	yh = thPade;
	if( bHH  &&  vRC6660[kHH] ) vRC6660[kHH]->Fill( xh, yh );
//hh                                ----------------------------
	kHH = 7;
	bHH = kHH < int( vRC6660.size() );
	xh = nMinus1;
	yh = thPadeX;
	if( bHH  &&  vRC6660[kHH] ) vRC6660[kHH]->Fill( xh, yh );
//hh                                ----------------------------

	double radDeg = cons->RadDeg();
	CsRCMirrors *mirr = CsRCMirrors::Instance();
	//CsRCPartPhotons* papho = ring->pPartPhot();
	//lPhotons = papho->lPhotons();
        //int kDetPart = papho->kDetPart();
	Hep3Vector vDcPartW = papho->vDcPartW()[kDetPart];
	for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	  int kDetClu = (*ih)->kDetClu();
	  float RR = mirr->RRv( kDetClu );
//------- assume photon from average point on part. traj. (vPoPhoEmW) :
	  Hep3Vector vPoPhoEmW = papho->vPoPhotW()[kDetClu];
//------- photon direction rel. to particle dir. :
          double al = sin((*ih)->the0()/1000.) * cos((*ih)->phiA()*radDeg);
          double am = sin((*ih)->the0()/1000.) * sin((*ih)->phiA()*radDeg);
          double an = cos((*ih)->the0()/1000.);
          Hep3Vector vDcPhoEmW( al, am, an );
          vDcPhoEmW = vDcPhoEmW.unit();
//------- photon direction to MWR :
          double theP, phiP;
          Hep3Vector vDcPhotW = CsRCUty::Instance()->
            rotfbcc( -1., vDcPartW, vDcPhoEmW, theP, phiP );
//          -----------------------------------------------
//------- photon impact on mirror (MWR) :
          Hep3Vector vPoC( 0., 0., 0. );
          Hep3Vector vPoPhoMir =
            mirr->vImpMir( vPoPhoEmW, vDcPhoEmW, vPoC, RR );
//          -----------------------------------------------
//------- normal to mirror at photon impact :
          Hep3Vector vDcNoPhoMir = (1./RR) * vPoPhoMir;
//------- photon reflected direction :
          double cosPhoMir = vDcNoPhoMir * vDcPhoEmW;
          Hep3Vector vDcPhoRefl = 2.*cosPhoMir * vDcNoPhoMir - vDcPhoEmW;
//------- photon incidence angles on detector :
          double tgxPho = vDcPhoRefl.x()/vDcPhoRefl.z();
          double tgyPho = vDcPhoRefl.y()/vDcPhoRefl.z();
	  double nMinus1 = 1./(betaIpo * cos( (*ih)->the0()/1000.) ) - 1.;
	  tgxPho *= 1000.;
	  kHH = 3;
	  bHH = kHH < int( vRC6660.size() );
	  xh = nMinus1;
	  yh = tgxPho;
	  if( bHH  &&  vRC6660[kHH] ) vRC6660[kHH]->Fill( xh, yh );
//hh                                  ----------------------------
	  tgyPho *= 1000.;
	  kHH = 4;
	  bHH = kHH < int( vRC6660.size() );
	  xh = nMinus1;
	  yh = tgyPho;
	  if( bHH  &&  vRC6660[kHH] ) vRC6660[kHH]->Fill( xh, yh );
//hh                                  ----------------------------
	  kHH = 5;
	  bHH = kHH < int( vRC6660.size() );
	  xh = nMinus1;
	  yh = momPart;
	  if( bHH  &&  vRC6660[kHH] ) vRC6660[kHH]->Fill( xh, yh );
//hh                                  ----------------------------
	}

//----- Added 081125 ---
	double likeDV = cons->likeDefVa();
        double* probaLK;
        probaLK = papho->probaLKAll();
//----- Production Likelihoods
        double likeBkg = papho->probaLKBgAll();
        double likeElec = probaLK[ 2];
        double likeMuon = probaLK[ 5];
        double likePion = probaLK[ 8];
        double likeKaon = probaLK[11];
        double likeProton = probaLK[14];
        bool bIDpart[4];
	bool bID = true;
        if( likeElec == likeBkg ) bID = false;      // NO photons from signal
//----- Select ID pions (analysis way)
        bool bIDpion = false;
        if( likePion > likeDV  &&  likePion > likeBkg  &&
            likePion > likeKaon  &&  likePion > likeProton )
	  bIDpion = true;
        if( momPart < 8.  &&  likePion < likeElec ) bIDpion = false;
	if( !bID ) bIDpion = false;
        bIDpart[0] = bIDpion;
//----- Select ID kaons (analysis way)
        bool bIDkaon = false;
        if( likeKaon > likeDV  &&  likeKaon > likeBkg  &&
            likeKaon > likePion  &&  likeKaon > likeProton )
	  bIDkaon = true;
        if( momPart < 8.  &&  likeKaon < likeElec ) bIDkaon = false;
	if( !bID ) bIDkaon = false;
        bIDpart[1] = bIDkaon;
//----- Select ID protons (analysis way) ?
        bool bIDproton = false;
        if( likeProton > likeDV  &&  likeProton > likeBkg  &&
            likeProton > likePion  &&  likeProton > likeKaon )
	  bIDproton = true;
        if( momPart < 8.  &&  likeProton < likeElec ) bIDproton = false;
	if( !bID ) bIDproton = false;
        bIDpart[2] = bIDproton;

	int kPart = 0;
	for( int kPaTy=8; kPaTy<=14; kPaTy+=3 ) {
	  kPart++;
	  if( !bIDpart[kPart-1] ) continue;
	  double massIpo = CsRCRecConst::Instance()->massPartv()[kPaTy];
	  double betaIpo = momPart/sqrt( momPart*momPart + massIpo*massIpo );
	  double nMinus1 = 0.;
	  if( betaIpo < 1. ) {
	    double den = betaIpo * cos( ring->the()/1000. );
	    if( den > 0. ) nMinus1 = 1./den - 1.;
	  }
	  kHH = kPart - 1;
	  bHH = kHH < int( vRC6670.size() );
	  xh = nMinus1;
	  yh = momPart;
	  if( bHH  &&  vRC6670[kHH] ) vRC6670[kHH]->Fill( xh, yh );
//hh                                  ----------------------------
	  break;
	}

      }
    //}

    return;
  }


  extern "C" {
    float hafitga_( const int, float&, float&,float&, float&, int&,
		    float* sig );
  }

//===========================================================================
  bool CsRCEventAnalysis::fitCFRefInf( double& index, double& sigIndex ) {
//------------------------------------------------------------------------


//- Paolo  -  August 27, 2002

    CsRCHistos& hist = CsRCHistos::Ref();

    float XMin = 0.00020;
    float XMax = 0.00245;
    size_t nBin = 100;
    if( hist.hRC3059 ) nBin = hist.hRC3059->GetDim(1).GetNBins();
    if( hist.hRC3059 ) XMin = hist.hRC3059->GetDim(1).GetMin();
    if( hist.hRC3059 ) XMax = hist.hRC3059->GetDim(1).GetMax();
    //cout << nBin << "  " << XMax << "  " << XMin << endl;
    float dX = (XMax - XMin)/float(nBin);
    size_t kBin;
    int kPoints = 0;
    double XFit[nBin];
    double binCto = 0.;
    double XAv = 0.;
    double XXAv = 0.;
    double sigX = 0.;
    double binMx = 0.;
    size_t kbinMx = 0;
//- fortran indexes for HBOOK...
    for( size_t kBin=1; kBin<=nBin; kBin++) {
      double binC = 0;
      if( hist.hRC3059 ) binC = hist.hRC3059->GetBinContent( kBin );
      if( binC > binMx ) {
	binMx = binC;
	kbinMx = kBin;
      }
    }
    size_t kb1 = kbinMx - 4;
    if( kb1 < 1 ) kb1 = 1;
    size_t kb2 = kbinMx + 4;
    if( kb2 > nBin ) kb2 = nBin;
    for( size_t kBin=kb1; kBin<=kb2; kBin++) {
      double binC = 0;
      if( hist.hRC3059 ) binC = hist.hRC3059->GetBinContent( kBin );
      if( binC > 0. ) {
        double X = XMin + kBin * dX - 0.5* dX;
        XAv += binC * X;
        XXAv += binC * X*X;
        binCto += binC;
      }
    }
    if( binCto > 0. ) {
      XAv /= binCto;
      sigX = XXAv/binCto - XAv*XAv;
      if( sigX < 0. ) sigX = 0.;
      sigX = sqrt( sigX/binCto );
    }
    index = XAv;
    sigIndex = sigX;

    return true;

  }



//===========================================================================
  bool CsRCEventAnalysis::fitGCFRefInd( double& index, double& sigIndex ) {
//-------------------------------------------------------------------------


//- Paolo  -  December 2003
//   rev.     April 2004


    CsRCExeKeys* key = CsRCExeKeys::Instance();

    static int binR[200];

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;
      std::cout << std::endl;
      key->acknoMethod( "CsRCEventAnalysis::fitGCFRefInd" );

      const int lineSize = 10;    //   test !!!
      char line[lineSize];
      ifstream f( "histo3559", ios::in );
      //ifstream f( "histo30580-09_11", ios::in );
      if( !f ) {}
      for( int kb=0; kb<200; kb++ ) {
        f.getline( line, lineSize, '\n' );
        istringstream s(line);
        s >> binR[kb];
        //std::cout << binR[kb] << std::endl;
      }    //   test !!!

    }

    CsRCHistos& hist = CsRCHistos::Ref();

    float XMin = 0.00050;
    float XMax = 0.00250;
    size_t nBin = 200;
    if( hist.hRC3559 ) nBin = hist.hRC3559->GetDim(1).GetNBins();
    if( hist.hRC3559 ) XMin = hist.hRC3559->GetDim(1).GetMin();
    if( hist.hRC3559 ) XMax = hist.hRC3559->GetDim(1).GetMax();
    //cout << nBin << "  " << XMax << "  " << XMin << endl;
    float dX = (XMax - XMin)/float(nBin);
    double binTot = 0.;
    size_t kBin;
    double binMx = 0.;
    size_t kbinMx = 0;
//- fortran indexes for HBOOK...
    for( size_t kBin=1; kBin<=nBin; kBin++) {
      double binC = 0;
      if( hist.hRC3559 ) binC = hist.hRC3559->GetBinContent( kBin );
      //binC = binR[kBin-1];    //   test !!!
      if( binC > binMx ) {
	binMx = binC;
	kbinMx = kBin;
      }
      binTot += binC;
    }
    //std::cout << binTot << "  " << kbinMx << std::endl;

    //if( binTot < 5000000. ) {
    if( binTot < 1000000. ) {
      std::cout << "  fitGCFRefInd : NOT enough statistics!" << std::endl;
      return  false;
    }
//@@-------------------------------------

    int nBinU = 7;
//----------------
    size_t kb1 = kbinMx - nBinU;
    if( kb1 < 1 ) kb1 = 1;
    size_t kb2 = kbinMx + nBinU;
    if( kb2 > nBin ) kb2 = nBin;
    int nBinFt = int( kb2 ) - int( kb1 ) + 1;
    double xCha[nBinFt];
    double yCont[nBinFt];
    int jBin = 0;
    for( size_t kBin=kb1; kBin<=kb2; kBin++) {
      //xCha[jBin] = XMin + kBin * dX + dX/2.;
      xCha[jBin] = XMin + kBin * dX - dX/2.;         //   040422
      double binC = 0;
      if( hist.hRC3559 ) binC = hist.hRC3559->GetBinContent( kBin );
      //binC = binR[kBin-1];    //   test !!!
      //std::cout << binC<< std::endl;
      yCont[jBin] = binC;
      jBin++;
    }
    double xBinMx = XMin + kbinMx * dX + dX/2.;

    int nParam = 10;
    double param[nParam];
    //param[0] = binTot/float( 4.* nBinFt );
    param[0] = binMx;
    //param[1] = 0.0014;
    param[1] = xBinMx;
    param[2] = 0.0001;
    param[3] = 0.;
    param[4] = 0.;
    param[5] = 1.;
    param[6] = 0.;
    param[7] = 0.;
    param[8] = 0.;
    param[9] = 0.;
    int iPaFit[nParam];
    iPaFit[0] = 0;
    iPaFit[1] = 0;
    iPaFit[2] = 0;
    iPaFit[3] = 1;
    iPaFit[4] = 1;
    iPaFit[5] = 1;
    iPaFit[6] = 1;
    iPaFit[7] = 1;
    iPaFit[8] = 1;
    iPaFit[9] = 1;

    CsRCGauPolFit oIndexFt( nBinFt, xCha, yCont, nParam, param, iPaFit );
//  --------------------------------------------------------------------

    if( oIndexFt.doChiFit() );
//      ---------------------

    if( !oIndexFt.flag() ) {
      std::cout << "  fitGCFRefInd : fit failure!" << std::endl;
      return  false;
    }


    //oIndexFt.print();
    //oIndexFt.doHist( vHist );
//  -------------------------

    index = oIndexFt.para()[1];
    //sigIndex = oIndexFt.para()[2];
    sigIndex = sqrt( oIndexFt.covPf()[2][2] );
    //std::cout << index << "  " << sigIndex << std::endl;

    return true;

  }


//===========================================================================
  void CsRCEventAnalysis::setProbsToTrack() {
//-------------------------------------------


//--- set RICH output variables to CsTrack
//    ------------------------------------
//--- Paolo  -  December 2002
//    Rev.      April 2007


//--- To set the RICH output buffer (to CsTrack) size
//    set in .../coral/src/track/CsTrack.h the line
//    #define CSTRACK_RICHDATASIZE (size)   with size = nProb
//    -------------------------------------------------------
//--- partProbs content list :
//----------------------------
// nProb = 15
//-   partProbs[  0 ]   background Like
//-   partProbs[  1 ]   pion Like
//-   partProbs[  2 ]   kaon Like
//-   partProbs[  3 ]   proton Like
//-   partProbs[  4 ]   pion Like derivative (or APV photons only)
//-   partProbs[  5 ]   kaon Like derivative (or APV photons only)
//-   partProbs[  6 ]   proton Like derivative (or APV photons only)
//-   partProbs[  7 ]   maximum likelihood angle
//-   partProbs[  8 ]   reconstr. ring angle
//-   partProbs[  9 ]   number of photons per ring
//-   partProbs[ 10 ]   fitted ring angle
//-   partProbs[ 11 ]   ring chisquare
//-   partProbs[ 12 ]   pion chisquare
//-   partProbs[ 13 ]   kaon chisquare
//-   partProbs[ 14 ]   proton chisquare
// nProb = 21 (new)
//-   partProbs[ 15 ]   electron Like
//-   partProbs[ 16 ]   muon Like
//-   partProbs[ 17 ]   electron Like derivative (or APV photons only)
//-   partProbs[ 18 ]   muon Like derivative (or APV photons only)
//-   partProbs[ 19 ]   electron chisquare
//-   partProbs[ 20 ]   muon chisquare

// nProb = 23 (proposed, but never implemented)!
//-   implemented 100305
//-   partProbs[ 21 ]   ring photon mean time (PMT only)
//-   partProbs[ 22 ]   number of ring PMT photons

// nProb = 30 (new new) --- (080506)
//-   partProbs[ 21 ]   ring photon mean time (PMT only)
//-   partProbs[ 22 ]   number of ring PMT photons
//-   partProbs[ 23 ]   pion Like derivative PMT photons only
//-   partProbs[ 24 ]   kaon Like derivative PMT photons only
//-   partProbs[ 25 ]   proton Like derivative PMT photons only
//-   partProbs[ 26 ]   electron Like derivative PMT photons only
//-   partProbs[ 27 ]   muon Like derivative PMT photons only
//-   partProbs[ 28 ]   number of Like APV photons
//-   partProbs[ 29 ]   number of Like PMT photons
//# nProb = 35   useless ???
//-   partProbs[ 30 ]   pion Like PMT photons only
//-   partProbs[ 31 ]   kaon Like PMT photons only
//-   partProbs[ 32 ]   proton Like PMT photons only
//-   partProbs[ 33 ]   electron Like PMT photons only
//-   partProbs[ 34 ]   muon Like PMT photons only


      CsRCExeKeys* key = CsRCExeKeys::Instance();
      CsRCRecConst *cons = CsRCRecConst::Instance();
      CsRCHistos& hist = CsRCHistos::Ref();
      int kh;
      double xh, wh;

      static bool firstCall = true;
      if( firstCall ) { 
        firstCall = false;

        key->acknoMethod( "CsRCEventAnalysis::setProbsToTrack" );
      }

      xh = 90.5;
      if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                       ------------------------

      static const int nProb = cons->outBufferSize();


//--- LOOP over PART-PHOTONS :
//    ------------------------
      CsRCEventPartPhotons* evepapho = CsRCEventPartPhotons::Instance();
      list<CsRCPartPhotons*> lPaPhotons = evepapho->lPartPhotons();
      list<CsRCPartPhotons*>::iterator ipf;

      for( ipf=lPaPhotons.begin(); ipf!=lPaPhotons.end(); ipf++ ) {
	if( !(*ipf) ) continue;
	CsRCPartPhotons* papho = (*ipf);

//      have an output for all processed particles (Yann)    //   070430
//      -------------------------------------------------
	//if( !papho->flag() ) continue;                     //   070430

	xh = 91.5;
	if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                         ------------------------


	CsRCRing* ring = papho->pRing();
	bool RING = false;
	if( ring != NULL  &&  ring->flag() ) RING = true;
//@@    ------------------------------------------------

	double likeDV = cons->likeDefVa();
	double* partProbs;
	double* probaLK;
	double* derivLK;
	double* derivLKPMT;
	double* qSquare;
	double qSquaR = 0.;
	int nPhoRing = 0;
	double pMinus1[nProb];
	double pZero[nProb];
	for( int k=0; k<nProb; k++ ) pMinus1[k] = -1.;
	for( int k=0; k<nProb; k++ ) pZero[k] = 0.;
	//partProbs = pMinus1;
	partProbs = pZero;
	double lZero[31];
	for( int k=0; k<31; k++ ) lZero[k] = 0.;
	probaLK = lZero;
	derivLK = lZero;
	derivLKPMT = lZero;
	qSquare = lZero;
	double theRecoLkR = 0.;
	//for( int k=0; k<31; k++ ) std::cout << qSquare[k] << "  ";
	//std::cout << std::endl;
	bool PIDset = false;

//----- set defaults to partProbs   (070509) :
//      --------------------------------------
//      Likelihoods (4) and DLike/dn (3) :
	for( int k=0; k<7; k++ ) partProbs[k] = likeDV;
//      Angles (3), nr. of photons (1) :
        for( int k=7; k<11; k++ ) partProbs[k] = 0.;
//      Chisqs (4) :
        for( int k=11; k<15; k++ ) partProbs[k] = 0.;
        if( nProb > 15 ) {
//      Likelihoods (2) and DLike/dn (2) :
          for( int k=15; k<19; k++ ) partProbs[k] = likeDV;
//      Chisqs (2) :
          for( int k=19; k<21; k++ ) partProbs[k] = 0.;
        }
//----- (080506) :
        if( nProb > 21 ) {
//      Ring mTime, ring PMT photons (2)
          for( int k=21; k<23; k++ ) partProbs[k] = 0.;
	}
//----- (100305) :
        if( nProb > 23 ) {
//      DLike/dn PMT only (5) :
	  for( int k=23; k<28; k++ ) partProbs[k] = likeDV;
//      Like photons APV and PMT (2) :
	  for( int k=28; k<30; k++ ) partProbs[k] = 0.;
//      Like PMT only (5) :
	  //# for( int k=30; k<35; k++ ) partProbs[k] = likeDV;
	}


	if( RING ) { 
//------- reconstructed Cerenkov angle :
//        ------------------------------
	  partProbs[  8 ] = ring->the();
	  partProbs[ 10 ] = ring->thetaRFit();
//@@      -----------------------------------

//------- from QSQUARE :
//        --------------
	  partProbs[  9 ] = ring->nPhotQsQ();
	  partProbs[ 11 ] = ring->ringQsQ();
//@@      ---------------------------------
	  qSquare = ring->qSquareRing();
	  partProbs[ 12 ] = qSquare[ 8];
	  partProbs[ 13 ] = qSquare[11];
	  partProbs[ 14 ] = qSquare[14];
//@@      -----------------------------
	  if( nProb > 15 ) {
    	    partProbs[ 19 ] = qSquare[ 2];
	    partProbs[ 20 ] = qSquare[ 5];
//@@        -----------------------------
	  }
	  //partProbs[ 12 ] = ring->partProb( 12 );    //   !!! WARNING
	  //partProbs[ 13 ] = ring->partProb( 13 );    //   !!! WARNING
	  //partProbs[ 14 ] = ring->partProb( 14 );    //   !!! WARNING
//@@      --------------------------------------
	  if( nProb > 21 ) {
	    double mT0 = ring->mT0();
	    double mTime = ring->mTime();
	    //double mTime = ring->mTime() + mT0;
    	    if( ring->mTimeSet() ) partProbs[ 21 ] = mTime;
	    partProbs[ 22 ] = ring->nPhotPMT();
//@@        ----------------------------------
	    //          std::cout << std::setprecision(4)<<ring->mTime() << "  " << ring->mT0() << "  " <<
            // ring->nPhotPMT() << std::endl;
	    // 	    if( ring->nPhotPMT() != ring->lPhotons().size() )
	    // 	      std::cout << std::precision(4)<<"@@@@@@@@@@@@@@  " << ring->nPhotPMT() << "  "
	    // 	         	  << ring->lPhotons().size() << std::endl;
	  }
	}

	if( cons->likeType() == "ALL"  &&  !RING ) {     //   041116
//----- 070509
	  //partProbs[  9 ] = papho->nPhotAll();         //   041116
//@@      -----------------------------------
	  //qSquare = papho->probaChiAll();
	  //partProbs[ 12 ] = qSquare[ 8];               //   041116
	  //partProbs[ 13 ] = qSquare[11];               //   041116
	  //partProbs[ 14 ] = qSquare[14];               //   041116
//@@      -----------------------------
	  //if( nProb > 15 ) {
	  //partProbs[ 19 ] = qSquare[ 2];
	  //partProbs[ 20 ] = qSquare[ 5];
//@@        -----------------------------
	  //}
//----- 070509
	}                                              //   041116

//----- from PART-PHOTONS  particle identification :
//      --------------------------------------------
	if( cons->likeType() == "ALL" ) {

	  partProbs[ 0 ] = papho->probaLKBgAll();
//@@      --------------------------------------
	  probaLK = papho->probaLKAll();
	  derivLK = papho->derivLKAll();
	  partProbs[ 1 ] = probaLK[ 8];
	  partProbs[ 2 ] = probaLK[11];
	  partProbs[ 3 ] = probaLK[14];
	  partProbs[ 4 ] = derivLK[ 8];
	  partProbs[ 5 ] = derivLK[11];
	  partProbs[ 6 ] = derivLK[14];
//@@      ----------------------------
	  if( nProb > 15 ) {
    	    partProbs[ 15 ] = probaLK[ 2];
	    partProbs[ 16 ] = probaLK[ 5];
	    partProbs[ 17 ] = derivLK[ 2];
	    partProbs[ 18 ] = derivLK[ 5];
//@@        -----------------------------
	  }
	  if( papho->thetaLikeSet() ) partProbs[ 7 ] = papho->thetaLike();
//@@                                  -----------------------------------
//------- (080506)
	  derivLKPMT = papho->derivLKAllPMT();
	  ///if( nProb > 21 ) {
//------- (100305)
	  if( nProb > 23 ) {
    	    partProbs[ 23 ] = derivLKPMT[ 8];
	    partProbs[ 24 ] = derivLKPMT[11];
	    partProbs[ 25 ] = derivLKPMT[14];
	    partProbs[ 26 ] = derivLKPMT[ 2];
	    partProbs[ 27 ] = derivLKPMT[ 5];
	    partProbs[ 28 ] = papho->nPhotAllAPV(); 
	    partProbs[ 29 ] = papho->nPhotAllPMT();
    	    //# partProbs[ 30 ] = probaLKPMT[ 8];
	    //# partProbs[ 31 ] = probaLKPMT[11];
	    //# partProbs[ 32 ] = probaLKPMT[14];
	    //# partProbs[ 33 ] = probaLKPMT[ 2];
	    //# partProbs[ 34 ] = probaLKPMT[ 5];
//@@        --------------------------------------
	  }
//------- (080506)

	  PIDset = true;
	}

//----- RING  particle identification :
//      -------------------------------
	else if( cons->likeType() == "RING"  &&  RING ) {

	  partProbs[ 0 ] = ring->probaLKBgRing();
//@@      --------------------------------------
	  probaLK = ring->probaLKRing();
	  derivLK = ring->derivLKRing();
	  partProbs[ 1 ] = probaLK[ 8];
	  partProbs[ 2 ] = probaLK[11];
	  partProbs[ 3 ] = probaLK[14];
	  partProbs[ 4 ] = derivLK[ 8];
	  partProbs[ 5 ] = derivLK[11];
	  partProbs[ 6 ] = derivLK[14];
//@@      ----------------------------
	  if( nProb > 15 ) {
	    partProbs[ 15 ] = probaLK[ 2];
	    partProbs[ 16 ] = probaLK[ 5];
	    partProbs[ 17 ] = derivLK[ 2];
	    partProbs[ 18 ] = derivLK[ 5];
//@@        ----------------------------
	  }
	  if( ring->thetaLikeSet() ) partProbs[ 7 ] = ring->thetaLike();
//                                   ----------------------------------
	  PIDset = true;
	}

	//else  continue;                 //   out (for CHISQ ) 030303


//      protect against out of range values :
//      -------------------------------------
	//cout << endl;
	//for( int k=0; k<nProb; k++ ) cout << partProbs[k] << "  ";
	//cout << endl;

//----- Likelihoods (4) :
        //for( int k=0; k<4; k++ ) {
	  //if( partProbs[k] > 1. ) {
	    //cout << "CsRCEventAnalysis::setProbsToTrack() : WARNING -- "
	    //    << "partProbs " << k << " > 1.  " << partProbs[k] << endl;
	  //}
	//}
//----- DLike/dn (3) :
        //for( int k=4; k<7; k++ ) {
	  //if( fabs( partProbs[k] ) > 10000. ) {
	    //cout << "CsRCEventAnalysis::setProbsToTrack(): WARNING -- "
	    //     << "partProbs " << k << " > 1000.  " << partProbs[k]
	    //     << endl;
	  //}
	//}
//----- Theta Like and Theta Rec (2) :
        for( int k=7; k<9; k++ ) {
	  if( partProbs[k] < 0. ) partProbs[k] = 0.;
	  if( partProbs[k] > 100. ) {
	    //cout << "CsRCEventAnalysis::setProbsToTrack() : WARNING -- "
	    //    << "partProbs " << k << " > 100.  " << partProbs[k] << endl;
	    partProbs[k] = 0.;
	  }
	}
//----- N photons (1) :
        for( int k=9; k<10; k++ ) {
	  if( partProbs[k] < 0. ) partProbs[k] = 0.;
	  //if( partProbs[k] > 200. ) {
	    //cout << "CsRCEventAnalysis::setProbsToTrack() : WARNING -- "
	    //     << "partProbs " << k << " > 200.  " << partProbs[k]
	    //     << endl;
	  //}
	}
//----- Theta Fit (1) :
        for( int k=10; k<11; k++ ) {
	  if( partProbs[k] < 0. ) partProbs[k] = 0.;
	  if( partProbs[k] > 100. ) {
	    //cout << "CsRCEventAnalysis::setProbsToTrack() : WARNING -- "
	    //    << "partProbs " << k << " > 100.  " << partProbs[k] << endl;
	    partProbs[k] = 0.;
	  }
	}
//----- Ring Chisquare (1) :
        for( int k=11; k<12; k++ ) {
	  if( partProbs[k] < 0. ) partProbs[k] = 0.;
	  if( partProbs[k] > 1000. ) {
	    //cout << "CsRCEventAnalysis::setProbsToTrack() : WARNING -- "
	    //     << "partProbs " << k << " > 1000.  " << partProbs[k]
	    //     << endl;
	    partProbs[k] = 1000.;
	  }
	}
//----- Chisquares (3) :
        for( int k=12; k<15; k++ ) {
	  if( partProbs[k] < 0. ) partProbs[k] = 0.;
	  if( partProbs[k] > 1000. ) {
	    //cout << "CsRCEventAnalysis::setProbsToTrack() : WARNING - "
	    //     << "partProbs " << k << " > 1000.  " << partProbs[k]
	    //     << endl;
	    partProbs[k] = 1000.;
	  }
	}
	if( nProb > 15 ) {
//------- Likelihoods (2) :
          //for( int k=15; k<17; k++ ) {
	  //}
//------- DLike/dn (2) :
          //for( int k=17; k<19; k++ ) {
	    //if( fabs( partProbs[k] ) > 10000. ) {
	    //}
	  //}
//------- Chisquares (2) :
          for( int k=19; k<21; k++ ) {
	    if( partProbs[k] < 0. ) partProbs[k] = 0.;
	    if( partProbs[k] > 1000. ) {
	      partProbs[k] = 1000.;
	    }
	  }
	}
//----- (080506)
	if( nProb > 21 ) {
//------- ring mTime, ring PMT photons (2) :
          //for( int k=21; k<23; k++ ) {
	  //}
	}
//----- (100305)
	if( nProb > 23 ) {
//------- DLike/dn PMT (5) :
          //for( int k=23; k<28; k++ ) {
	    //if( fabs( partProbs[k] ) > 10000. ) {
	      //cout << "CsRCEventAnalysis::setProbsToTrack(): WARNING -- "
	      //     << "partProbs " << k << " > 1000.  " << partProbs[k]
	      //     << endl;
	    //}
	  //}
//------- N photons (2) :
          for( int k=28; k<30; k++ ) {
	    if( partProbs[k] < 0. ) partProbs[k] = 0.;
	  }
//------- Like PMT (5) :
          for( int k=30; k<35; k++ ) {
	    //#
	  }
	}

//----- 070430 :
	if( !papho->flag()  ||  !PIDset ) {
//        Likelihoods (4) and DLike/dn (3) :
	  for( int k=0; k<7; k++ ) partProbs[k] = likeDV;
//        Angles (3), nr. of photons (1) :
	  for( int k=7; k<11; k++ ) partProbs[k] = 0.;
	  //use a code ?
	  partProbs[ 9 ] = -1.;
//@@--------------------------
//        Chisqs(4) :
	  for( int k=11; k<15; k++ ) partProbs[k] = 0.;
  	  if( nProb > 15 ) {
//          Likelihoods (2) and DLike/dn (2) :
	    for( int k=15; k<19; k++ ) partProbs[k] = likeDV;
//          Chisqs (2) :
	    for( int k=19; k<21; k++ ) partProbs[k] = 0.;
	  }
//------- (080506)
          if( nProb > 21 ) {
//          Ring mTime, ring PMT photons
            for( int k=21; k<23; k++ ) partProbs[k] = 0.;
	  }
//------- (100305)
          if( nProb > 23 ) {
//          DLike/dnPMT (5) :
	    for( int k=23; k<28; k++ ) partProbs[k] = likeDV;
//          Nr. of photons APV & PMT (2) :
	    for( int k=28; k<30; k++ ) partProbs[k] = 0.;
//          Like PMT (5) :
	    //# for( int k=30; k<35; k++ ) partProbs[k] = likeDV;
	  }


	  xh = 92.5;
	  if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                           ------------------------
	}
//----- 070430


//----- tests and prints
	//for( int k=0; k<nProb; k++ ) cout << partProbs[k] << "  ";
	//cout << endl;
	//if( PIDset ) {
	//  for( int k=0; k<nProb; k++ ) cout << partProbs[k] << "  ";
	//  cout << "PIDset" << endl;
	//}
	bool tout = false;
	if( tout  &&  PIDset ) {
	  //for( int k=0; k<nProb; k++ ) cout << partProbs[k] << "  ";
	  //cout << endl;
          bool dmp = true;
	  for( int k=0; k<4; k++ ) if( partProbs[k] != 0. ) dmp = false;
	  if( dmp ) {
	    CsRCParticle* part = ring->pPartPhot()->pPart();
	    float momPart = part->mom();
	    bool zero = false;
	    //bool zero = true;
	    //for( int k=4; k<nProb; k++ ) if( partProbs[k] != 0. ) 
	    //zero = true;
	    if( zero ) {
	      float momPart = ring->pPartPhot()->pPart()->mom();
	      for( int k=0; k<nProb; k++ ) cout << partProbs[k] << "  ";
	      cout << ring->lPhotons().size() << "  " << momPart << endl;
	    }
	  }
	}

//----- partProbs dump table:
//      ---------------------
	if( key->kPrintPartProbs() == 1 ) {
	  //if( partProbs[11] == 0.) {
	  //if( partProbs[16] > 100000.) {
	  //if( partProbs[2] > 0.) {

	  std::cout << std::endl;
	  std::cout << std::endl;
	  std::cout << " PartProbs";
	  //std::cout << " PartProbs (" << nProb << ")";
	  std::cout << "  -  Ev. " << CsRichOne::Instance()->kEvent();
	  std::cout << "  -  Part. " << papho->pPart()->kPart();
	  std::cout << "  -  mom " << papho->pPart()->mom();
	  std::cout << "  -  Buffer " << nProb;
	  std::cout << std::endl;

	  //std::cout << " " << RING << "  " << true << std::endl;
	  //std::cout << setprecision( 9 );
	  //std::cout << setiosflags(ios::scientific);

	  std::cout << "   ";
	  std::cout << "Like bk  " << "Like pi  " << "Like ka  "
		    << "Like pr  ";
	  ///if( nProb <= 21 ) {
	  if( nProb <= 23 ) {
	    std::cout << " DLpi  " << " DLka  " << " DLpr  ";
	  } else {
	    std::cout << "DLpiA  " <<" DLkaA  " << "DLprA  ";
	  }
	  std::cout << " ThML  " << " ThRg  " << "nPhRg  " << " ThFt  "
		    << "  ChRg  " << "  Chpi  " << "  Chka  " << "  Chpr  ";
	  std::cout << std::endl;

	  std::cout << "   ";
	  for( int k=0; k<4; k++ ) { 
	    std::string varStr = getVarStr( partProbs[k], 7, 5 );
	    std::cout << varStr << "  ";
	  }
	  for( int k=4; k<7; k++ ) {
	    std::string varStr = getVarStr( partProbs[k], 5, -2 );
	    std::cout << varStr << "  ";
	  }
	  for( int k=7; k<11; k++ ) {
	    std::string varStr = getVarStr( partProbs[k], 5, 1 );
	    std::cout << varStr << "  ";
	  }
	  for( int k=11; k<15; k++ ) {
	    std::string varStr = getVarStr( partProbs[k], 6, -2 );
	    std::cout << varStr << "  ";
	  }
	  std::cout << std::endl;
	  //std::cout << setprecision( 5 );
	  //std::cout << "   ";
	  //for( int k=0; k<4; k++ ) std::cout << partProbs[k] << "  ";
	  //std::cout << setprecision( 3 );
	  //for( int k=4; k<7; k++ ) std::cout << partProbs[k] << "  ";
	  //std::cout << setprecision( 2 );
	  //for( int k=7; k<15; k++ ) std::cout << partProbs[k] << "  ";
	  //std::cout << std::endl;

	  if( nProb > 15 ) {
	    std::cout << std::endl;
	    std::cout << "   ";
	    std::cout << "Like el  " << "Like mu  ";
	    ///if( nProb <= 21 ) {
	    if( nProb <= 23 ) {
	      std::cout << " DLel  " << " DLmu  ";
	    } else {
	      std::cout << "DLelA  " << "DLmuA  ";
	    }
	    std::cout << "  Chel  " << "  Chmu";
	    if( nProb == 21 ) std::cout << std::endl;
//--------- 100305
	    if( nProb > 21 ) {
	      std::cout << "  ";
	      std::cout << " RgTime  " << "nPhTm  ";
	    }
	    if( nProb == 23 ) std::cout << std::endl;

	    std::cout << "   ";
	    for( int k=15; k<17; k++ ) {
	      std::string varStr = getVarStr( partProbs[k], 7, 5 );
	      std::cout << varStr << "  ";
	    }
	    for( int k=17; k<19; k++ ) {
	      std::string varStr = getVarStr( partProbs[k], 5, -2 );
	      std::cout << varStr << "  ";
	    }
	    for( int k=19; k<21; k++ ) {
	      std::string varStr = getVarStr( partProbs[k], 6, -2 );
	      std::cout << varStr << "  ";
	    }
	    if( nProb == 21 ) std::cout << std::endl;
	    //std::cout << setprecision( 5 );
	    //std::cout << "   ";
	    //for( int k=15; k<17; k++ ) std::cout << partProbs[k] << "  ";
	    //std::cout << setprecision( 3 );
	    //for( int k=17; k<19; k++ ) std::cout << partProbs[k] << "  ";
	    //std::cout << setprecision( 2 );
	    //for( int k=19; k<21; k++ ) std::cout << partProbs[k] << "  ";
	    //std::cout << std::endl;
	    //if( nProb == 21 ) std::cout << std::endl;
	  }
//------- 100305
	  if( nProb > 21 ) {
	    std::string varStr = getVarStr( partProbs[21], 7, 0 );
	    std::cout << varStr << "  ";
	    varStr = getVarStr( partProbs[22], 5, 0 );
	    std::cout << varStr << "  ";
	    //for( int k=21; k<23; k++ ) std::cout << partProbs[k] << "  ";
	    if( nProb == 23 ) std::cout << std::endl;
	  }
//------- 080506
	  //if( nProb > 21 ) {
	  //  std::cout << std::endl;
	  //  std::cout << "   ";
	  //  std::cout << "RgTime    " << "nPhTm  ";
	  //  if( nProb == 23 ) std::cout << std::endl;
	  //}

//------- 100305
	  if( nProb > 23 ) {
	    std::cout << "DLpiP   " << "DLKaP  " << "DLprP  "
		      << "DLelP   " << "DLmuP   ";
	    std::cout << "nPhA  " << "nPhP ";
	    std::cout << std::endl;
	    std::cout << "   ";
	  }
	  //if( nProb > 21 ) {
	  //  for( int k=21; k<23; k++ ) std::cout << partProbs[k] << "  ";
	  //  if( nProb == 23 ) std::cout << std::endl;
	  //}
	  if( nProb > 23 ) {
	    std::cout << setprecision( 3 );
	    for( int k=23; k<28; k++ ) std::cout << partProbs[k] << "  ";
	    std::cout << setprecision( 2 );
	    for( int k=28; k<nProb; k++ ) std::cout << partProbs[k] << "  ";
	    //# for( int k=28; k<30; k++ ) std::cout << partProbs[k] << "  ";
	    std::cout << setprecision( 5 );
	    //# for( int k=30; k<35; k++ ) std::cout << partProbs[k] << "  ";
	    std::cout << setprecision( 2 );
	    std::cout << std::endl;
	  }
	}

//----- reset buffers :
//      ---------------
	if( cons->likeType() == "ALL" ) {
	  papho->setPartProbs( partProbs );
//@@--------------------------------------
	}
	else if( cons->likeType() == "RING"  &&  RING ) {
	  ring->setPartProbs( partProbs );
//@@-------------------------------------
	}
      	else if( cons->likeType() == "CHISQ"  &&  RING ) {
//        not implemented!
//@@-------------------------------------
	}

//--------------------------------------------------------
//----- put to track the RICH out buffer :
//      ----------------------------------
//      have an output for all particles (Yann)         //   070430
	//if( PIDset ) {                                //   070430
	  CsRCParticle* part = papho->pPart();
	  CsTrack* track = part->pTrack();
          //if( track ) {                               //   080729
	  if( track  &&  !key->readMyFile() ) {         //   080729
	    track->setRich1Probs( &partProbs[0] );
//@@        -------------------------------------
	    xh = 93.5;
	    if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                             ------------------------
	  }
	//}

      }   // end loop on part-photons

      xh = 99.5;
      if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                       ------------------------

      return;

  }


//===========================================================================
  void CsRCEventAnalysis::doMirrDetFit() {
//----------------------------------------


//- fits the ring on photon detectors for each mirror element
//  ---------------------------------------------------------
//- Paolo  -  January 2003


    CsRCRecConst *cons = CsRCRecConst::Instance();
    CsRCExeKeys *key = CsRCExeKeys::Instance();

    static int kOffHi = cons->kOffHi();
    CsRCHistos& hist = CsRCHistos::Ref();
    double xh, yh, wh;
    static vector<CsHist2D*> vRC2100, vRC2300, vRC2500, vRC2700;

    list<CsRCMirrorElem*> lMirrEleAlg =
      CsRCMirrors::Instance()->lMirrEleAlg();

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;
      key->acknoMethod( "CsRCEventAnalysis::doMirrDetFit" );

      CsHistograms::SetCurrentPath("/RICH");
//@@---------------------------------------
      int kmirr = 0;
      list<CsRCMirrorElem*>::iterator amirr;
      for( amirr=lMirrEleAlg.begin(); amirr!=lMirrEleAlg.end(); amirr++ ) {

	stringstream hTitle;
	stringstream hN2100;
	hN2100 << kOffHi + 2100 + kmirr;
	hTitle << (*amirr)->name();
	//hTitle << (*amirr)->name() <<
	//  "  C-offset y vs x - fit on det";
	vRC2100.push_back( new CsHist2D( hN2100.str(), hTitle.str(),
					 100, -10., 10., 100, -10., 10. ) );
	//	cout << hTitle.str() << "  " << kmirr << endl;
	hTitle.clear();

	stringstream hN2300;
	hTitle << (*amirr)->name() <<
	  "  R offset vs nPoints - fit on det";
	hN2300 << kOffHi + 2300 + kmirr;
	vRC2300.push_back( new CsHist2D( hN2300.str(), hTitle.str(),
					 100, -10., 10., 48, 0., 48. ) );
	hTitle.clear();

	stringstream hN2500;
	hTitle << (*amirr)->name() <<
	  "  chisq/nu vs n-iter - fit on det";
	hN2500 << kOffHi + 2500 + kmirr;
	vRC2500.push_back( new CsHist2D( hN2500.str(), hTitle.str(),
					 100, 0., 10., 20, 0., 20. ) );
	hTitle.clear();

	stringstream hN2700;
	hTitle << (*amirr)->name() <<
	  "  useful particle impacts";
	hN2700 << kOffHi + 2700 + kmirr;
	vRC2700.push_back( new CsHist2D( hN2700.str(), hTitle.str(),
			   100, -1500., 1500., 100, -1500., 1500. ) );
	hTitle.clear();

	kmirr++;
      }
      CsHistograms::SetCurrentPath("/");
//@@-----------------------------------
    }


//- loop over rings :
//  -----------------
    list<CsRCRing*> lRings = CsRCEventRings::Instance()->lRings();
    list<CsRCRing*>::iterator ir;
    for( ir=lRings.begin(); ir!=lRings.end(); ir++ ) {
      if( !(*ir)->flag() ) continue;

      //float CFRefInd = cons->CFRefInd();

      CsRCPartPhotons* ipf = (*ir)->pPartPhot();
      if( !ipf ) break;
      Hep3Vector vPoPaMir;
      vPoPaMir = ipf->pPart()->vPoPaMir0()[ipf->pPart()->kDetPart()];

      double rrElmQ = cons->rrAliMir() * cons->rrAliMir();
      string mirEve = "@@@";
      list<CsRCMirrorElem*>::iterator amirr;
      for( amirr=lMirrEleAlg.begin(); amirr!=lMirrEleAlg.end(); amirr++ ) {
	//Hep3Vector vCePos = (*amirr)->vpos() + (*amirr)->mirNo()->vC0();
	Hep3Vector vCePos = (*amirr)->vpos();   //   030122
	if( (vPoPaMir - vCePos).mag2() < rrElmQ ) {
	  mirEve = (*amirr)->name();  break;
	}
      }

//--- mirror element circle fit, use ring points :
//    --------------------------------------------
      double theReco = (*ir)->the();
      if( theReco > cons->theAliMir() ) {
	int kElmAli = 0;
	list<CsRCMirrorElem*>::iterator amirr;
	for( amirr=lMirrEleAlg.begin(); amirr!=lMirrEleAlg.end(); amirr++ ) {
	  if( (*amirr)->name() == mirEve ) {
	    vector<CsHist2D*> vDummy;
	    for( int k=0; k<7; k++ ) vDummy.push_back( NULL );
	    CsRCCircleFit oCircle;
	    if( (*ir)->getDetRFit( vDummy, oCircle ) ) {
//              ------------------------------------

	      CsRCCircleFit* pCircle = &oCircle;
	      if( pCircle->flag()  &&  pCircle->nIter() < 10 ) {

		xh = pCircle->para()[1];
		yh = pCircle->para()[2];
		if( vRC2100[kElmAli] ) vRC2100[kElmAli]->Fill( xh, yh );
//                                     --------------------------------
		xh = pCircle->para()[0] - (*ir)->the()* 3.36;
		yh = pCircle->degFree() + 3.;
		if( vRC2300[kElmAli] ) vRC2300[kElmAli]->Fill( xh, yh );
//                                     --------------------------------
		xh = pCircle->chiSquare();
		yh = pCircle->nIter();
		if( vRC2500[kElmAli] ) vRC2500[kElmAli]->Fill( xh, yh );
//                                     --------------------------------
		xh = vPoPaMir.x();
		yh = vPoPaMir.y();
		if( vRC2700[kElmAli] ) vRC2700[kElmAli]->Fill( xh, yh );
//                                     --------------------------------
	      }
	    }
	  }
	  kElmAli++;
	}
      }
    }

  }


//===========================================================================
  void CsRCEventAnalysis::checkLikeDistr() {
//------------------------------------------


//- Paolo  - January  2006


    if( !CsRCHistos::Ref().bookHis()  ||
        CsRCHistos::Ref().levelBk() < 2 ) return;
//  --------------------------------------------

    CsRCExeKeys* key = CsRCExeKeys::Instance();
    CsRCRecConst *cons = CsRCRecConst::Instance();
    static int kOffHi = cons->kOffHi();
    static double radDeg = cons->RadDeg();
    CsRCMirrors *mirr = CsRCMirrors::Instance();
    CsRCHistos& hist = CsRCHistos::Ref();
    int kh;
    double xh, yh, wh;

    static std::vector<CsHist2D*> hRCLikeDis;
    static std::vector<CsHist2D*> hRCLikeMDis;
    static std::vector<CsHist2D*> hRCLikeRMDis;
    static std::vector<CsHist2D*> hRCLikeMMDis;
    static std::vector<CsHist2D*> hRCRdist;
    
    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      key->acknoMethod( "CsRCEventAnalysis::checkLikeDistr" );

      std::string hTitle = "Like distr";
      int kName = 0;
      for( int kH=0; kH<4; kH++ ) {
	double xlm, xlx;
	double rlm, rlx;
	if( cons->likeLog() ) {
	  xlm = -800.;
	  xlx =  200.;
	  rlm =    0.;
	  rlx =  500.;
	} else {
	  xlm = 0.;
	  xlx = 2.;
	  rlm = 1.;
	  //rlx = 3.;
	  rlx = 10.;
	}
	hTitle = "Like distr";
	kName = kOffHi + 3450 + kH;
	stringstream hNLikeDis;
	hNLikeDis << kName;
	hRCLikeDis.push_back( new CsHist2D( hNLikeDis.str(), hTitle,
		  	    100, xlm, xlx, 6, 0., 600. ) );
	hTitle = "LikeM distr";
	kName = kOffHi + 3460 + kH;
	stringstream hNLikeMDis;
	hNLikeMDis << kName;
	hRCLikeMDis.push_back( new CsHist2D( hNLikeMDis.str(), hTitle,
		  	    100, xlm, xlx, 6, 0., 600. ) );
	hTitle = "Like RM distr";
	kName = kOffHi + 3470 + kH;
	stringstream hNLikeRMDis;
	hNLikeRMDis << kName;
	hRCLikeRMDis.push_back( new CsHist2D( hNLikeRMDis.str(), hTitle,
		  	    100, rlm, rlx, 6, 0., 600. ) );
	hTitle = "Like MM distr";
	kName = kOffHi + 3480 + kH;
	stringstream hNLikeMMDis;
	hNLikeMMDis << kName;
	hRCLikeMMDis.push_back( new CsHist2D( hNLikeMMDis.str(), hTitle,
		  	    100, rlm, rlx, 6, 0., 600. ) );
      }
      for( int kH=0; kH<6; kH++ ) {
	hTitle = "Rad dist.";
	kName = kOffHi + 3490 + kH;
	stringstream hNRdist;
	hNRdist << kName;
	hRCRdist.push_back( new CsHist2D( hNRdist.str(), hTitle,
			    100, -1200, 1200., 100, -1200., 1200. ) );
      }
    }

    static const int nProb = cons->outBufferSize();

    double likeDV = cons->likeDefVa();
//@@---------------------------------

//- LOOP over PART-PHOTONS :
//  ------------------------
    CsRCEventPartPhotons* evepapho = CsRCEventPartPhotons::Instance();
    list<CsRCPartPhotons*> lPaPhotons = evepapho->lPartPhotons();
    list<CsRCPartPhotons*>::iterator ipf;

    for( ipf=lPaPhotons.begin(); ipf!=lPaPhotons.end(); ipf++ ) {
      if( !(*ipf) ) continue;
      CsRCPartPhotons* papho = (*ipf);
      if( !papho->flag() ) continue;

      CsRCRing* ring = papho->pRing();
      bool RING = false;
      if( ring != NULL  &&  ring->flag() ) RING = true;
//@@  ------------------------------------------------

      double likeV[6];
      double* probaLK;
      double likeBkg = likeDV;
      double likeElec = likeDV;
      double likeMuon = likeDV;
      double likePion = likeDV;
      double likeKaon = likeDV;
      double likeProton = likeDV;

//--- from PART-PHOTONS  particle identification :
//    --------------------------------------------
      if( cons->likeType() == "ALL" ) {
	likeBkg = papho->probaLKBgAll();
	probaLK = papho->probaLKAll();
      }
//--- RING  particle identification :
//    -------------------------------
      else if( cons->likeType() == "RING"  &&  RING ) {
	likeBkg = ring->probaLKBgRing();
        probaLK = ring->probaLKRing();
      }
      else {
	probaLK = papho->probaLKAll();
      }
      likePion = probaLK[ 8];
      likeKaon = probaLK[11];
      likeProton = probaLK[14];
      if( nProb > 15 ) {
        likeElec = probaLK[ 2];
        likeMuon = probaLK[ 5];
      }
      //for( int k=0; k<15; k++ ) std::cout << probaLK[k] << "  ";
      //std::cout << std::endl;
      int nLike = 4;
      likeV[0] = likeBkg;
      likeV[1] = likePion;
      likeV[2] = likeKaon;
      likeV[3] = likeProton;
      likeV[4] = likeDV;
      likeV[5] = likeDV;
      if( nProb > 15 ) {
       	nLike = 6;
//@@-------------
        likeV[4] = likeElec;
	likeV[5] = likeMuon;
      }
      //for( int k=0; k<6; k++ ) std::cout << likeV[k] << "  ";
      //std::cout << std::endl;

      double* qSquare;
      double qsqElec = 0.;
      double qsqMuon = 0.;
      double qsqPion = 0.;
      double qsqKaon = 0.;
      double qsqProton = 0.;

//--- from QSQUARE :
//    --------------
      if( RING ) { 
        qSquare = ring->qSquareRing();
      }
      else if( cons->likeType() == "ALL"  &&  !RING ) {
	qSquare = papho->probaChiAll();
      }       
      else {
        qSquare = papho->probaChiAll();
      }
      qsqPion = qSquare[ 8];
      qsqKaon = qSquare[11];
      qsqProton = qSquare[14];
      if( nProb > 15 ) {
        qsqElec = qSquare[ 2];
        qsqMuon = qSquare[ 5];
      }

      CsRCParticle* part = papho->pPart();

      if( part->mom() < 2.5) continue;
      if( part->mom() > 60.) continue;
//@@---------------------------------

      double xPade = part->vPade()[part->kDetPart()].x();
      double yPade = part->vPade()[part->kDetPart()].y();
      double yPadeD = fabs(yPade) - 470.;
      double rPade = sqrt( xPade*xPade + yPadeD*yPadeD );

      int iPade = int( rPade/100.);
      if( iPade > 5 ) iPade = 5;
      xh = xPade;
      yh = yPade;
      if( hRCRdist[iPade] ) hRCRdist[iPade]->Fill( xh, yh );
//hh                        -------------------------------

      for( int kH=0; kH<4; kH++ ) {
	if( (cons->likeLog() && likeV[kH] != likeDV)  ||  likeV[kH] > 0.) {
	  xh = likeV[kH];
	  yh = rPade;
	  if( hRCLikeDis[kH] ) hRCLikeDis[kH]->Fill( xh, yh );
//hh                           ------------------------------
	}
      }

      double likeBack = likeV[0];
      if( likeBack == likeDV ) return;
//    -------------------------------

      double likeMax = -100000.;
      int kMax = -1;
      for( int kL=0; kL<nLike; kL++ ) {
	if( likeV[kL] == likeDV ) continue;
	if( likeV[kL] > likeMax ) {
	  likeMax = likeV[kL];
	  kMax = kL;
	}
      }
      double likeMax2 = -100000.;
      int kMax2 = -1;
      for( int kL=0; kL<nLike; kL++ ) {
	if( likeV[kL] == likeDV ) continue;
	if( kL == kMax ) continue;
	if( likeV[kL] > likeMax2 ) {
	  likeMax2 = likeV[kL];
	  kMax2 = kL;
	}
      }
      if( kMax >= 0  &&  kMax < 4 ) {
	if( cons->likeLog()  ||  likeMax > 0.) {
	  xh = likeMax;
	  yh = rPade;
	  if( hRCLikeMDis[kMax] ) hRCLikeMDis[kMax]->Fill( xh, yh );
//hh                              ---------------------------------
        }
	if( cons->likeLog()  ||  (likeMax > 0. && likeBack > 0.) ) {
	  if( cons->likeLog() ) {
	    xh = likeMax - likeBack;
	  } else {
	    xh = likeMax / likeBack;
	  }
          yh = rPade;
          if( hRCLikeRMDis[kMax] ) hRCLikeRMDis[kMax]->Fill( xh, yh );
//hh                               ----------------------------------
        }
      }
      if( kMax >= 0  &&  kMax < 4  &&  kMax2 >= 0 ) {
	if( cons->likeLog()  ||  (likeMax > 0. && likeMax2 > 0.) ) {
	  if( cons->likeLog() ) {
	    xh = likeMax - likeMax2;
	  } else {
	    xh = likeMax / likeMax2;
	  }
          yh = rPade;
          if( hRCLikeMMDis[kMax] ) hRCLikeMMDis[kMax]->Fill( xh, yh );
//hh                               ----------------------------------
        }
      }
    }

    return;
  }


//===========================================================================
  void CsRCEventAnalysis::checkThetaLikeMax() {
//---------------------------------------------


//- Paolo  - March  2007

//- Checks PMT ONLY!
//  Approximate analytic computation of Theta-Like-Max


    if( !CsRCHistos::Ref().bookHis()  ||
        CsRCHistos::Ref().levelBk() < 2 ) return;
//  --------------------------------------------

    CsRCExeKeys* key = CsRCExeKeys::Instance();
    CsRCRecConst *cons = CsRCRecConst::Instance();
    static int kOffHi = cons->kOffHi();
    static double radDeg = cons->RadDeg();
    CsRCMirrors *mirr = CsRCMirrors::Instance();
    CsRCHistos& hist = CsRCHistos::Ref();
    int kh;
    double xh, yh, wh;

    static std::vector<CsHist1D*> hRCTheDiff;
    
    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      key->acknoMethod( "CsRCEventAnalysis::checkThetaLikeMax" );

      std::string hTitle = "TheDiff distr";
      int kName = 0;
      for( int kH=0; kH<5; kH++ ) {
	kName = kOffHi + 3440 + kH;
	stringstream hNTheDiff;
	hNTheDiff << kName;
	hRCTheDiff.push_back( new CsHist1D( hNTheDiff.str(), hTitle,
		  	    100, -5., 5. ) );
      }
    }

    static const int nProb = cons->outBufferSize();
    static double likeDV = cons->likeDefVa();

//- LOOP over PART-PHOTONS :
//  ------------------------
    CsRCEventPartPhotons* evepapho = CsRCEventPartPhotons::Instance();
    list<CsRCPartPhotons*> lPaPhotons = evepapho->lPartPhotons();
    list<CsRCPartPhotons*>::iterator ipf;
    bool isPFFirst = true;
    for( ipf=lPaPhotons.begin(); ipf!=lPaPhotons.end(); ipf++ ) {
      if( !(*ipf) ) continue;
      CsRCPartPhotons* papho = (*ipf);
      if( !papho->flag() ) continue;
      if( isPFFirst ) papho->setLikeFirst( true );
      isPFFirst = false;

      CsRCRing* ring = papho->pRing();
      bool RING = false;
      if( ring != NULL  &&  ring->flag() ) RING = true;
//@@  ------------------------------------------------

      static CsRCLikeAll* likeAll = NULL;
      if( cons->backgrType() == "05" ) likeAll = new CsRCLikeAll05( papho );
      if( !likeAll ) return;
//    ---------------------

      double likeV[6];
      double* probaLK;
      double likeBkg = likeDV;

//--- from PART-PHOTONS  particle identification :
//    --------------------------------------------
      if( cons->likeType() == "ALL" ) {
	likeBkg = papho->probaLKBgAll();
	probaLK = papho->probaLKAll();
      }
//--- RING  particle identification :
//    -------------------------------
      else if( cons->likeType() == "RING"  &&  RING ) {
	likeBkg = ring->probaLKBgRing();
        probaLK = ring->probaLKRing();
      }
      else {
	probaLK = papho->probaLKAll();
      }
      //for( int k=0; k<15; k++ ) std::cout << probaLK[k] << "  ";
      //std::cout << std::endl;

      int nLike = 4;
      likeV[0] = likeBkg;
      if( likeBkg == likeDV ) continue;
//    --------------------------------
      likeV[1] = probaLK[ 8];
      likeV[2] = probaLK[11];
      likeV[3] = probaLK[14];
      likeV[4] = likeDV;
      likeV[5] = likeDV;
      if( nProb > 15 ) {
       	nLike = 6;
        likeV[4] = probaLK[ 2];
	likeV[5] = probaLK[ 5];
      }
      //for( int k=0; k<6; k++ ) std::cout << likeV[k] << "  ";
      //std::cout << std::endl;

      CsRCParticle* part = papho->pPart();
      if( part->mom() < 2.5 ) continue;
      if( part->mom() > 60. ) continue;
//@@----------------------------------

      double likeMax = -100000.;
      int kMax = -1;
      for( int kL=0; kL<nLike; kL++ ) {
	if( likeV[kL] == likeDV ) continue;
	if( likeV[kL] > likeMax ) {
	  likeMax = likeV[kL];
	  kMax = kL;
	}
      }
//    select KAONS only :
//    -------------------
      if( kMax == 2 ) {

	int kIpo = 11;
	double theIpo = papho->thetaIpoVS( kIpo );

	if( !papho->thetaLikeSet() ) continue;
	double theMaxApp = papho->thetaLike();

	double normS = likeAll->normSignal( theMaxApp );
//                     --------------------------------
	double normB = likeAll->normBackgr( theMaxApp );
//                     --------------------------------
	list<CsRCPhoton*> lPhotons = papho->lPhotons();
	list<CsRCPhoton*>::iterator ih;
	long double theNum = 0.;
	long double theDen = 0.;
	double theWNum = 0.;
	double theWDen = 0.;
	bool isPFirst = true;
	int kPhot = 0;
	for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	  if( !(*ih)->flag() ) continue;
	  if( !(*ih)->isPMT() ) continue;

	  if( isPFirst ) (*ih)->setLikeFirst( true );
	  isPFirst = false;

	  double signW = likeAll->likeSignal( (*ih), theMaxApp );
//                       ---------------------------------------
	  double backW = likeAll->likeBackgr( (*ih), theMaxApp );
//                       ---------------------------------------
	  double weight = signW / ( signW + backW );
	  double sigPhot = (*ih)->sigmaPhoPid( papho );
	  theNum += weight * (*ih)->the() / ( sigPhot*sigPhot );
	  theDen += weight / ( sigPhot*sigPhot );
	  theWNum += (*ih)->the() / ( sigPhot*sigPhot );
	  theWDen += 1./ ( sigPhot*sigPhot );

	  kPhot++;
	}
	double theMax = 0.;
	double theWgt = 0.;
	if( theDen > 0. ) {
	  //std::cout << setiosflags(ios::scientific);
	  //std::cout << theNum << "  " << theDen << "  " << kPhot
	  //          << std::endl;
	  //std::cout << resetiosflags(ios::scientific);
  	  theMax = theNum / theDen;
	  if( theWDen > 0. ) theWgt = theWNum / theWDen;
	  //std::cout << CsRichOne::Instance()->kEvent() << "  "
	  //	  << papho->kPaPhot() << "  " << theMax << "  "
	  //	  << theMaxApp << "  " << theWgt << "  " << theIpo
	  //      << std::endl;
	  int kH = 0;
	  xh = theMax - theIpo;
	  if( hRCTheDiff[kH] ) hRCTheDiff[kH]->Fill( xh );
//hh                           --------------------------
	  kH = 1;
	  xh = theMaxApp - theIpo;
	  if( hRCTheDiff[kH] ) hRCTheDiff[kH]->Fill( xh );
//hh                           --------------------------
	  kH = 2;
	  xh = theWgt - theIpo;
	  if( hRCTheDiff[kH] ) hRCTheDiff[kH]->Fill( xh );
//hh                           --------------------------
	  kH = 3;
	  xh = theMax - theMaxApp;
	  if( hRCTheDiff[kH] ) hRCTheDiff[kH]->Fill( xh );
//hh                           --------------------------
	}
      }
    }

    return;
  }


//===========================================================================
  void CsRCEventAnalysis::moniCFRefIndR() {
//-----------------------------------------


//- Paolo  -  December  2008


    if( !CsRCHistos::Ref().bookHis() ||
        CsRCHistos::Ref().levelBk() < 2 ) return;
//  --------------------------------------------

    CsRCHistos& hist = CsRCHistos::Ref();
    float xh, yh;

    CsOpt* opt = CsOpt::Instance();
    CsRCRecConst *cons = CsRCRecConst::Instance();
    CsRCDetectors *dets = CsRCDetectors::Instance();
    int kOffHi = cons->kOffHi();

    static std::vector<CsHist2D*> vRC6680;
    static std::vector<CsHist2D*> vRC6900;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      for( int kh=0; kh<20; kh++ ) vRC6680.push_back( NULL );
      for( int kh=0; kh<38; kh++ ) vRC6900.push_back( NULL );
      if( CsRCHistos::Ref().bookHis() && CsRCHistos::Ref().levelBk() >= 2 ) {
	CsHistograms::SetCurrentPath("/RICH");
	vRC6680.clear();
        for( int kh=0; kh<5; kh++ ) {
	  string hTitle = "n-1 from Photons";
	  stringstream hN6680;
	  int kHist = kOffHi + 6680 + kh + 1;
	  hN6680 << kHist;
	  vRC6680.push_back( new CsHist2D( hN6680.str(), hTitle,
					   200, 0.0002, 0.0022,
					   100, 0., 50.) );
	}
	for( int kh=5; kh<12; kh++ ) {
	  string hTitle = "n-1 from Rings";
	  stringstream hN6680;
	  int kHist = kOffHi + 6680 + kh + 1;
	  hN6680 << kHist;
	  vRC6680.push_back( new CsHist2D( hN6680.str(), hTitle,
					   400, 0.0002, 0.0022,
					   100, 0., 50.) );
	}
	for( int kh=12; kh<19; kh++ ) {
	  string hTitle = "n-1 from Photons";
	  stringstream hN6680;
	  int kHist = kOffHi + 6680 + kh + 1;
	  hN6680 << kHist;
	  vRC6680.push_back( new CsHist2D( hN6680.str(), hTitle,
					   200, 0.0002, 0.0022,
					   100, 0., 50.) );
	}
	vRC6900.clear();
        for( int kh=0; kh<12; kh++ ) {
	  string hTitle = "n-1 from Photons";
	  stringstream hN6900;
	  int kHist = kOffHi + 6900 + kh + 1;
	  hN6900 << kHist;
	  vRC6900.push_back( new CsHist2D( hN6900.str(), hTitle,
					   200, 0.0002, 0.0022,
					   100, 0., 50.) );
	}
        for( int kh=12; kh<38; kh++ ) {
	  string hTitle = "n-1 from Rings";
	  stringstream hN6900;
	  int kHist = kOffHi + 6900 + kh + 1;
	  hN6900 << kHist;
	  vRC6900.push_back( new CsHist2D( hN6900.str(), hTitle,
					   400, 0.0002, 0.0022,
					   100, 0., 50.) );
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

//--- particle momentum :
      CsRCParticle* part = papho->pPart();
      float momPart = part->mom();
//--- particle angle at RICH entrance window :
      double ttxy;
      ttxy = sqrt( pow( part->vDirIn().x()/part->vDirIn().z(), 2 ) +
  	  	   pow( part->vDirIn().y()/part->vDirIn().z(), 2 ) );
      ttxy *= 1000.;

//    SELECTION using FIXED limits
//    ----------------------------
//--- mom lower limit
      double momLLimit = 5.;       // GeV/c
      if( momPart < momLLimit ) continue;
//@@------------------------------------
//--- mom upper limit
      double momULimit = 50.;      // GeV/c
      if( momPart > momULimit ) continue;
//@@------------------------------------
//--- tgxy lower limit
      double tgxyLLimit = 20.;      // ~mrad
      if( ttxy < tgxyLLimit ) continue;
//@@----------------------------------

//--- assume PION mass
      float massIpo = CsRCRecConst::Instance()->massPartv()[8];

//--- n-1 from ALL PHOTONs angle
//    --------------------------
      list<CsRCPhoton*> lPhotonsA = papho->lPhotons();
      list<CsRCPhoton*>::iterator ih;
      for( ih=lPhotonsA.begin(); ih!=lPhotonsA.end(); ih++ ) {
        if( !(*ih)->flag() ) continue;
        double thePho = (*ih)->the0();
//      assume PION mass
	double betaIpo = momPart / sqrt( momPart*momPart + massIpo*massIpo );
	double nMinus1 = 1./(betaIpo * cos( thePho/1000. )) - 1.;
	xh = nMinus1;
	yh = momPart;
	if( (*ih)->isPMT() ) {
	  int kQua = -1;
	  int cClu = (*ih)->ptToClu()->ic();
//        determine photon Cathode
	  for( int kc=0; kc<4; kc++ ) {
	    if( cClu == dets->nCatPMT()[kc] ) kQua = kc;
	  }
	  int kHH = -1;
	  if( kQua >= 0 ) kHH = kQua;
	  if( kHH >= 0  &&  kHH < int( vRC6680.size() ) ) {
	    if( vRC6680[kHH] ) vRC6680[kHH]->Fill( xh, yh );
//hh                           ----------------------------
	  }
	}
	if( (*ih)->isAPV() ) {
	  int kHH = 4;
	  if( kHH >= 0  &&  kHH < int( vRC6680.size() ) ) {
	    if( vRC6680[kHH] ) vRC6680[kHH]->Fill( xh, yh );
//hh                           ----------------------------
	  }
	  int cClu = (*ih)->ptToClu()->ic();
//        determine photon Cathode
	  int kGG = cClu;
	  if( cClu > 3 ) kGG -= 1;
	  if( cClu > 5 ) kGG -= 1;
	  if( cClu > 10 ) kGG -= 1;
	  if( cClu > 12 ) kGG -= 1;
	  //std::cout << cClu << "  " << kGG << std::endl;
	  if( kGG >= 0  &&  kGG < int( vRC6900.size() ) ) {
	    if( vRC6900[kGG] ) vRC6900[kGG]->Fill( xh, yh );
//hh                           ----------------------------
	  }
	}
      }


//--- n-1 from RING angle
//    -------------------
      list<CsRCPhoton*> lPhotonsR = ring->lPhotons();

      double theRing = ring->the();
//--- assume PION mass
      double betaIpo = momPart / sqrt( momPart*momPart + massIpo*massIpo );
      double nMinus1 = 1./(betaIpo * cos( theRing/1000. )) - 1.;
      xh = nMinus1;
      yh = momPart;
      int cClu = -1;
      bool PMTonly = true;
      for( ih=lPhotonsR.begin(); ih!=lPhotonsR.end(); ih++ ) {
        if( !(*ih)->flag() ) continue;
	cClu = (*ih)->ptToClu()->ic();
	if( (*ih)->isPMT() ) continue;
	PMTonly = false;
	break;
      }
      if( cClu < 0 ) continue;
      bool APVonly = true;
      cClu = -1;
      for( ih=lPhotonsR.begin(); ih!=lPhotonsR.end(); ih++ ) {
        if( !(*ih)->flag() ) continue;
	cClu = (*ih)->ptToClu()->ic();
	if( (*ih)->isAPV() ) continue;
	APVonly = false;
	break;
      }
      if( cClu < 0 ) continue;

//--- n-Photon/Ring lower limit
      int nPhoLLimit = 0;
      if( PMTonly ) nPhoLLimit = 10;
      if( APVonly ) nPhoLLimit =  4;
      if( int( lPhotonsR.size() ) < nPhoLLimit ) continue;
//@@-----------------------------------------------------

      if( PMTonly ) {
	int kHH = 0 + 5;
        if( kHH >= 0  &&  kHH < int( vRC6680.size() ) ) {
	  if( vRC6680[kHH] ) vRC6680[kHH]->Fill( xh, yh );
//hh                         ----------------------------
	}
      }
      if( APVonly ) {
	int kHH = 1 + 5;
	if( kHH >= 0  &&  kHH < int( vRC6680.size() ) ) {
	  if( vRC6680[kHH] ) vRC6680[kHH]->Fill( xh, yh );
//hh                         ----------------------------
	}
      }
      if( !PMTonly  &&  !APVonly ) {
	int kHH = 2 + 5;
	if( kHH >= 0  &&  kHH < int( vRC6680.size() ) ) {
	  if( vRC6680[kHH] ) vRC6680[kHH]->Fill( xh, yh );
//hh                         ----------------------------
	}
      }

//    determine if RING photons are one Cathode only
      int kQua = -1;
      if( PMTonly ) {
        for( int kc=0; kc<4; kc++ ) {
	  for( ih=lPhotonsR.begin(); ih!=lPhotonsR.end(); ih++ ) {
	    if( !(*ih)->flag() ) continue;
	    cClu = (*ih)->ptToClu()->ic();
	    kQua = kc;
	    if( cClu == dets->nCatPMT()[kc] ) continue;
	    kQua = -1;
	    break;
	  }
	  if( kQua >= 0 ) break;
        }
        if( kQua >= 0 ) {
	  int kHH = kQua + 8;
	  if( kHH >= 0  &&  kHH < int( vRC6680.size() ) ) {
	    if( vRC6680[kHH] ) vRC6680[kHH]->Fill( xh, yh );
//hh                           ----------------------------
	  }
        }
      }
      int kCat = -1;
      if( APVonly ) {
        for( int kc=0; kc<16; kc++ ) {
	  for( ih=lPhotonsR.begin(); ih!=lPhotonsR.end(); ih++ ) {
	    if( !(*ih)->flag() ) continue;
	    cClu = (*ih)->ptToClu()->ic();
	    kCat = kc;
	    if( cClu == kCat ) continue;
	    kCat = -1;
	    break;
	  }
	  if( kCat >= 0 ) break;
        }
        if( kCat >= 0 ) {
	  int kGG = kCat + 12;
	  if( cClu > 3 ) kGG -= 1;
	  if( cClu > 5 ) kGG -= 1;
	  if( cClu > 10 ) kGG -= 1;
	  if( cClu > 12 ) kGG -= 1;
	  if( kGG >= 0  &&  kGG < int( vRC6900.size() ) ) {
	    if( vRC6900[kGG] ) vRC6900[kGG]->Fill( xh, yh );
//hh                           ----------------------------
	  }
        }
      }

//--- n-1 from RING PHOTONs angle
//    ---------------------------
      for( ih=lPhotonsR.begin(); ih!=lPhotonsR.end(); ih++ ) {
        if( !(*ih)->flag() ) continue;
        double thePho = (*ih)->the0();
//      assume PION mass
	double betaIpo = momPart / sqrt( momPart*momPart + massIpo*massIpo );
	double nMinus1 = 1./(betaIpo * cos( thePho/1000. )) - 1.;
	xh = nMinus1;
	yh = momPart;
	if( PMTonly ) {
	  int kHH = 0 + 12;
	  if( kHH >= 0  &&  kHH < int( vRC6680.size() ) ) {
	    if( vRC6680[kHH] ) vRC6680[kHH]->Fill( xh, yh );
//hh                           ----------------------------
	  }
	}
        if( APVonly ) {
	  int kHH = 1 + 12;
	  if( kHH >= 0  &&  kHH < int( vRC6680.size() ) ) {
	    if( vRC6680[kHH] ) vRC6680[kHH]->Fill( xh, yh );
//hh                           ----------------------------
	  }
        }
        if( !PMTonly  &&  !APVonly ) {
	  int kHH = 2 + 12;
	  if( kHH >= 0  &&  kHH < int( vRC6680.size() ) ) {
	    if( vRC6680[kHH] ) vRC6680[kHH]->Fill( xh, yh );
//hh                           ----------------------------
	  }
        }
	if( PMTonly ) {
          if( kQua >= 0 ) {
	    int kHH = kQua + 15;
	    if( kHH >= 0  &&  kHH < int( vRC6680.size() ) ) {
	      if( vRC6680[kHH] ) vRC6680[kHH]->Fill( xh, yh );
//hh                             ----------------------------
	    }
          }
	}
	if( APVonly ) {
          if( kCat >= 0 ) {
	    int kGG = kCat + 24;
	    if( cClu > 3 ) kGG -= 1;
	    if( cClu > 5 ) kGG -= 1;
	    if( cClu > 10 ) kGG -= 1;
	    if( cClu > 12 ) kGG -= 1;
	    if( kGG >= 0  &&  kGG < int( vRC6900.size() ) ) {
	      if( vRC6900[kGG] ) vRC6900[kGG]->Fill( xh, yh );
//hh                             ----------------------------
	    }
          }
	}
      }

    }

    return;
  }



//===========================================================================
  void CsRCEventAnalysis::checkRingEff() {
//----------------------------------------


//- Paolo  -  August  2009


    if( !CsRCHistos::Ref().bookHis() ||
        CsRCHistos::Ref().levelBk() < 2 ) return;
//  --------------------------------------------

    CsRCHistos& hist = CsRCHistos::Ref();
    float xh, yh;

    CsOpt* opt = CsOpt::Instance();
    CsRCRecConst *cons = CsRCRecConst::Instance();
    CsRCDetectors *dets = CsRCDetectors::Instance();
    int kOffHi = cons->kOffHi();

    static std::vector<CsHist1D*> vRC7000;
    static std::vector<CsHist1D*> vRC7050;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      for( int kh=0; kh<50; kh++ ) vRC7000.push_back( NULL );
      for( int kh=0; kh<50; kh++ ) vRC7050.push_back( NULL );
      if( CsRCHistos::Ref().bookHis() && CsRCHistos::Ref().levelBk() >= 2 ) {
	CsHistograms::SetCurrentPath("/RICH");
	vRC7000.clear();
        for( int kh=0; kh<50; kh++ ) {
	  string hTitle = " ";
	  stringstream hN7000;
	  int kHist = kOffHi + 7000 + kh + 1;
	  hN7000 << kHist;
	  vRC7000.push_back( new CsHist1D( hN7000.str(), hTitle,
					   100, 0., 50.) );
	}
	vRC7050.clear();
        for( int kh=0; kh<50; kh++ ) {
	  string hTitle = " ";
	  stringstream hN7050;
	  int kHist = kOffHi + 7050 + kh + 1;
	  hN7050 << kHist;
	  vRC7050.push_back( new CsHist1D( hN7050.str(), hTitle,
					   100, 0., 50.) );
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
      bool RING = false;
      if( ring  &&  ring->flag() ) RING = true;

//--- particle momentum :
      CsRCParticle* part = papho->pPart();
      float momPart = part->mom();
//--- particle angle at RICH entrance window :
      double ttxy;
      ttxy = sqrt( pow( part->vDirIn().x()/part->vDirIn().z(), 2 ) +
  	  	   pow( part->vDirIn().y()/part->vDirIn().z(), 2 ) );
      ttxy *= 1000.;

//    SELECTION using FIXED limits
//    ----------------------------
//--- mom lower limit
      ///double momLLimit = 5.;       // GeV/c
      double momLLimit = 2.5;       // GeV/c
      if( momPart < momLLimit ) continue;
//@@------------------------------------
//--- mom upper limit
      double momULimit = 50.;      // GeV/c
      if( momPart > momULimit ) continue;
//@@------------------------------------
//--- tgxy lower limit
      double tgxyLLimit = 20.;      // ~mrad
      if( ttxy < tgxyLLimit ) continue;
//@@----------------------------------

      xh = momPart;

      int kHH = 0;
      if( kHH >= 0  &&  kHH < int( vRC7000.size() ) ) {
	if( vRC7000[kHH] ) vRC7000[kHH]->Fill( xh );
//hh                       ------------------------
      }
      if( kHH >= 0  &&  kHH < int( vRC7050.size() ) ) {
        if( RING ) if( vRC7050[kHH] ) vRC7050[kHH]->Fill( xh );
//hh                                  ------------------------
      }

      int iCaPa = papho->iCaPa();
      //std::cout << iCaPa << std::endl;
      if( iCaPa < 0 ) continue;
      CsRCCathode* cat = dets->ptrToCat( iCaPa );
      double padx = cat->padx();
      double pady = cat->pady();
      double hCatx = cat->hCatx();
      double hCaty = cat->hCaty();
      double RingR = 55./1000. * 3300.;
      int kDetPart = papho->kDetPart();
      Hep3Vector vPoPaDetW = papho->vPoPaDetW()[kDetPart];
      Hep3Vector vOffCatWw = dets->vOffCatW( iCaPa );
      bool RingCo = false;
      if( vPoPaDetW.x()-RingR >= vOffCatWw.x()-hCatx  &&
          vPoPaDetW.x()+RingR <= vOffCatWw.x()+hCatx  &&
          vPoPaDetW.y()-RingR >= vOffCatWw.y()-hCaty  &&
          vPoPaDetW.y()+RingR <= vOffCatWw.y()+hCaty ) RingCo = true;

      kHH = iCaPa + 10;
      if( kHH >= 0  &&  kHH < int( vRC7000.size() ) ) {
        if( vRC7000[kHH] ) vRC7000[kHH]->Fill( xh );
//hh                       ------------------------
      }
      if( kHH >= 0  &&  kHH < int( vRC7050.size() ) ) {
        if( RING ) if( vRC7050[kHH] ) vRC7050[kHH]->Fill( xh );
//hh                                  ------------------------
      }
      if( RingCo ) {
	//std::cout << iCaPa << std::endl;
	kHH = iCaPa + 30;
	if( kHH >= 0  &&  kHH < int( vRC7000.size() ) ) {
	  if( vRC7000[kHH] ) vRC7000[kHH]->Fill( xh );
//hh                         ------------------------
	}
	if( kHH >= 0  &&  kHH < int( vRC7050.size() ) ) {
	  if( RING ) if( vRC7050[kHH] ) vRC7050[kHH]->Fill( xh );
//hh                                    ------------------------
	}
      }

      if( !RING ) continue;
      list<CsRCPhoton*> lPhotons = ring->lPhotons();
      list<CsRCPhoton*>::iterator ih;
      int cClu = -1;
      bool PMTonly = true;
      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
        if( !(*ih)->flag() ) continue;
	cClu = (*ih)->ptToClu()->ic();
	if( (*ih)->isPMT() ) continue;
	PMTonly = false;
	break;
      }
      if( cClu < 0 ) continue;
      bool APVonly = true;
      cClu = -1;
      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
        if( !(*ih)->flag() ) continue;
	cClu = (*ih)->ptToClu()->ic();
	if( (*ih)->isAPV() ) continue;
	APVonly = false;
	break;
      }
      if( cClu < 0 ) continue;

//--- n-Photon/Ring lower limit
      int nPhoLLimit = 0;
      if( PMTonly ) nPhoLLimit = 10;
      if( APVonly ) nPhoLLimit =  4;
      if( int( lPhotons.size() ) < nPhoLLimit ) continue;
//@@----------------------------------------------------

      xh = momPart;
      if( PMTonly ) {
	kHH = 1;
	//	if( kHH >= 0  &&  kHH < int( vRC7000.size() ) ) {
	//if( vRC7000[kHH] ) vRC7000[kHH]->Fill( xh );
//hh                         ------------------------
	//}
        //if( kHH >= 0  &&  kHH < int( vRC7050.size() ) ) {
	//if( RING ) if( vRC7050[kHH] ) vRC7050[kHH]->Fill( xh );
//hh                                    ------------------------
	//}
      }
      if( APVonly ) {
      }
      if( !PMTonly  &&  !APVonly ) {
      }

//    determine if RING photons are one Cathode only
      int kQua = -1;
      if( PMTonly ) {
        for( int kc=0; kc<4; kc++ ) {
	  for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	    if( !(*ih)->flag() ) continue;
	    cClu = (*ih)->ptToClu()->ic();
	    kQua = kc;
	    if( cClu == dets->nCatPMT()[kc] ) continue;
	    kQua = -1;
	    break;
	  }
	  if( kQua >= 0 ) break;
        }
        if( kQua >= 0 ) {
	  int kHH = kQua + 8;
        }
      }
      int kCat = -1;
      if( APVonly ) {
        for( int kc=0; kc<16; kc++ ) {
	  for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	    if( !(*ih)->flag() ) continue;
	    cClu = (*ih)->ptToClu()->ic();
	    kCat = kc;
	    if( cClu == kCat ) continue;
	    kCat = -1;
	    break;
	  }
	  if( kCat >= 0 ) break;
        }
        if( kCat >= 0 ) {
	  int kGG = kCat + 12;
	  if( cClu > 3 ) kGG -= 1;
	  if( cClu > 5 ) kGG -= 1;
	  if( cClu > 10 ) kGG -= 1;
	  if( cClu > 12 ) kGG -= 1;
	  //	  if( kGG >= 0  &&  kGG < int( vRC6900.size() ) ) {
	  //if( vRC6900[kGG] ) vRC6900[kGG]->Fill( xh, yh );
//hh                           ----------------------------
	  //}
        }
      }


    }

    return;
  }


//===========================================================================
  void CsRCEventAnalysis::checkMCPMTs() {
//---------------------------------------


//- Paolo  -  February  2009

    if( !CsRCHistos::Ref().bookHis() ||
        CsRCHistos::Ref().levelBk() < 2 ) return;
//  --------------------------------------------

    CsRCDetectors *dets = CsRCDetectors::Instance();
    CsRCMirrors *mirr = CsRCMirrors::Instance();
    CsRCHistos& hist = CsRCHistos::Ref();
    float xh, yh;

    int kEvent = CsRichOne::Instance()->kEvent();
    std::cout << std::endl;
    std::cout << "checkMCPMTs - Event " << kEvent << std::endl;

    CsRICH1Detector* rich = CsGeom::Instance()->getRich1Detector();
    list<CsMCHit*> lMCHits = rich->getMyMCHits();
    list<CsMCHit*>::iterator ii;
    //for( ii=lMCHits.begin(); ii!=lMCHits.end(); ii++ ) {
      //CsMCHit* hit = (*ii);
      //float xH = hit->getX();
      //float yH = hit->getY();
      //float tH = hit->getDTime();
      //std::cout << "hits " << tH << " " << xH << " " << yH << std::endl;
    //}
    //list<CsDigit*> lMyDigits = rich->getMyDigits();
    //list<CsDigit*>::iterator id;
    list<CsRCPhotonDet*> lPhoDet = dets->lPhoDet();
    list<CsRCMirrorNom*> lMirrNom = mirr->lMirrNom();
    CsRCEventPartPhotons* evepapho = CsRCEventPartPhotons::Instance();
    list<CsRCPartPhotons*> lPaPhotons = evepapho->lPartPhotons();
    list<CsRCPartPhotons*>::iterator ipf;
    for( ipf=lPaPhotons.begin(); ipf!=lPaPhotons.end(); ipf++ ) {
      if( !(*ipf) ) continue;
      CsRCPartPhotons* papho = (*ipf);

      CsRCParticle* part = papho->pPart();
      std::cout << "checkMCPMTs - Part " << part->kPart() << std::endl;

      list<CsRCPhoton*> lPhotons = papho->lPhotons();
      list<CsRCPhoton*>::iterator ih;
      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	if( !(*ih) ) continue;
	CsRCPhoton* pho = (*ih);

	CsRCCluster* clu = pho->ptToClu();
	if( !clu ) continue;
	int iCatc = clu->ic();
        double xc = clu->xc();
	double yc = clu->yc();
	int kClu = clu->kClu();
	list<CsRCPad*> lPads = clu->lPads();
	CsRCPad* pad = lPads.front();
	if( !pad ) continue;
	int kPad = pad->kPad();
	int ixp = pad->ix();
	int iyp = pad->iy();
	CsRCCathode* cat = dets->ptrToCat( iCatc );
	if( !cat ) continue;
	double padx = cat->padx();
	double pady = cat->pady();
	double hCatx = cat->hCatx();
	double hCaty = cat->hCaty();
	double xp = ((ixp+0.5)*padx - hCatx) + dets->vOffCatW( iCatc ).x();
	double yp = ((iyp+0.5)*pady - hCaty) + dets->vOffCatW( iCatc ).y();
	CsDigit* dig = pad->pDigit();
	if( !dig ) continue;
	int addr = dig->getAddress();
	int iCatd = rich->getCathode(addr)-1;
	//CsRCCathode* cath = dets->ptrToCat( iCatd );
	int ixd = rich->getPadX(addr)-1;
	int iyd = rich->getPadY(addr)-1;
	int nData = dig->getDataSize();
	double* data = dig->getData();
	bool print = false;
	//print = iCatd == 3 || iCatd == 5 || iCatd == 10 || iCatd == 12;
	print = !( iCatd == 3 || iCatd == 5 || iCatd == 10 || iCatd == 12 );
	list<CsMCHit*> lMCDgHits = (dynamic_cast<CsMCDigit*>(dig))->getHits();
	if( print ) std::cout << "h-size " << lMCDgHits.size() << std::endl;
        if( lMCDgHits.size() > 0 ) {
	  int kHit = 0;
	  for( ii=lMCDgHits.begin(); ii!=lMCDgHits.end(); ii++ ) {
	    CsMCRICH1Hit* hit = dynamic_cast<CsMCRICH1Hit*>( (*ii) );
	    if( !hit ) continue;
	    double th = hit->getDTime();
	    double xh = hit->getX();
	    double yh = hit->getY();
	    double zh = hit->getZ();
	    int iCath = hit->getCathode() - 1;
	    kHit++;
	    CsRCPhotonDet* phodet = lPhoDet.front();
	    if( yh < 0. ) phodet = lPhoDet.back();
	    if( !phodet ) continue;
	    CsRCMirrorNom* phomirr = lMirrNom.front();
	    if( yh < 0. ) phomirr = lMirrNom.back();
	    if( !phomirr ) continue;
	    Hep3Vector vPoC0( phomirr->vC0() );
	    //std::cout << "vPoC0  " << vPoC0 << std::endl;
	    //(0.00, 1600.00, 2535.00)
	    HepVector vrot( 3, 0 );
	    vrot[0] = xh;
	    vrot[1] = yh;
	    vrot[2] = zh;
	    for( int j=0; j<3; j++ ) vrot[j] -= vPoC0[j];
	    HepVector vMWR = phodet->rotMatrix() * vrot;
	    if( print ) {
	      //if( yh < 0. ) {
	      std::cout << "checkMCPMTs : cath " << iCath << "  ";
	      std::cout << "hit " << kHit << "  "
			<< xh << " " << yh << "  " << zh << " - ";
	      std::cout << "pad " << kPad
			<< "  " << xp << "  "  << yp << " - ";
	      std::cout << "clu " << kClu
			<< "  " << xc << "  "  << yc << "  -  ";
	      std::cout << "rot "
			<< vMWR[0] << "  " << vMWR[1] << "  " << vMWR[2];
	      std::cout << std::endl;
	      //}
	    }
	  }
	} else {
	//std::cout << "NO hits " << kDig << std::endl;
	}

      }
      //std::cout << "checkMCPMTs  " << lMCHits.size() << "  "
      //          << lMyDigits.size() << std::endl;
    }

    return;
  }


//===========================================================================
  void CsRCEventAnalysis::checkDataPMTs() {
//-----------------------------------------


//- Paolo  -  February  2009

    if( !CsRCHistos::Ref().bookHis() ||
        CsRCHistos::Ref().levelBk() < 2 ) return;
//  --------------------------------------------

    CsRCDetectors *dets = CsRCDetectors::Instance();
    CsRCMirrors *mirr = CsRCMirrors::Instance();
    CsRCHistos& hist = CsRCHistos::Ref();
    float xh, yh;

    CsRICH1Detector* rich = CsGeom::Instance()->getRich1Detector();
    //list<CsDigit*> lMyDigits = rich->getMyDigits();
    //list<CsDigit*>::iterator id;
    list<CsRCPhotonDet*> lPhoDet = dets->lPhoDet();
    list<CsRCMirrorNom*> lMirrNom = mirr->lMirrNom();
    CsRCEventPartPhotons* evepapho = CsRCEventPartPhotons::Instance();
    list<CsRCPartPhotons*> lPaPhotons = evepapho->lPartPhotons();
    list<CsRCPartPhotons*>::iterator ipf;
    for( ipf=lPaPhotons.begin(); ipf!=lPaPhotons.end(); ipf++ ) {
      if( !(*ipf) ) continue;
      CsRCPartPhotons* papho = (*ipf);

      CsRCParticle* part = papho->pPart();
      list<CsRCPhoton*> lPhotons = papho->lPhotons();
      list<CsRCPhoton*>::iterator ih;
      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	if( !(*ih) ) continue;
	CsRCPhoton* pho = (*ih);
	CsRCCluster* clu = pho->ptToClu();
	if( !clu ) continue;
	int iCat = clu->ic();
        double xc = clu->xc();
	double yc = clu->yc();
	int kClu = clu->kClu();
	list<CsRCPad*> lPads = clu->lPads();
	CsRCPad* pad = lPads.front();
	if( !pad ) continue;
	int kPad = pad->kPad();
	int ixp = pad->ix();
	int iyp = pad->iy();
	CsRCCathode* cat = dets->ptrToCat( iCat );
	if( !cat ) continue;
	double padx = cat->padx();
	double pady = cat->pady();
	double hCatx = cat->hCatx();
	double hCaty = cat->hCaty();
	double xp = ((ixp+0.5)*padx - hCatx) + dets->vOffCatW( iCat ).x();
	double yp = ((iyp+0.5)*pady - hCaty) + dets->vOffCatW( iCat ).y();
        //if( iCat == 3  ||  iCat == 5  ||  iCat == 10  ||  iCat == 12 ) {
	if( !( iCat == 3  ||  iCat == 5  ||  iCat == 10  ||  iCat == 12 ) ) {
	  std::cout << "checkDataPMTs : cath " << iCat << "  ";
	  std::cout << "pad " << kPad << "  " << xp << "  "  << yp << " - ";
	  std::cout << "clu " << kClu << "  " << xc << "  "  << yc << " - ";
	  std::cout << "dd "  << xp-xc << "  "  << yp-yc << " - ";
	//std::cout << "size "  << lPads.size() << " - ";
	//std::cout << padx << "  " << pady << "  " << hCatx << "  " << hCaty;
	  std::cout << std::endl;
	}
	CsDigit* dig = pad->pDigit();
	if( !dig ) continue;
	int addr = dig->getAddress();
	int iCath = rich->getCathode(addr)-1;
	CsRCCathode* cath = dets->ptrToCat( iCath );
	int ix = rich->getPadX(addr)-1;
	int iy = rich->getPadY(addr)-1;
	int nData = dig->getDataSize();
	double* data = dig->getData();
      }

    }

    return;
  }


//===========================================================================
  void CsRCEventAnalysis::checkOptCorr() {
//----------------------------------------


//- Paolo  -  December  2009


    if( !CsRichOne::Instance()->UpRICHJob() ) return;

    if( !CsRCHistos::Ref().bookHis() ||
        CsRCHistos::Ref().levelBk() < 2 ) return;
//  --------------------------------------------

    CsRCRecConst *cons = CsRCRecConst::Instance();
    int kOffHi = cons->kOffHi();   
    CsRCDetectors *dets = CsRCDetectors::Instance();
    CsRCMirrors *mirr = CsRCMirrors::Instance();

    CsRCHistos& hist = CsRCHistos::Ref();
    float xh, yh;
    int khh;

    static std::vector<std::string> mEName;
    static std::vector<Hep3Vector> vCePos0;

    static std::vector<CsHist2D*> vRC8500;
    static std::vector<CsHist2D*> vRC8600;
    static std::vector<CsHist2D*> vRC8700;

    static double xLmn = 0;
    static double xLmx = 0;
    static double yLmn = 0.;
    static double yLmx = 0.;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      mEName.push_back( "MB14" );
      mEName.push_back( "MB43" );
      mEName.push_back( "MT26" );
      mEName.push_back( "MT54" );
      std::list<CsRCMirrorElem*> lMirrEle = mirr->lMirrEle();
      std::list<CsRCMirrorElem*>::iterator emirr;
      for( int ke=0; ke<4; ke++ ) {
	for( emirr=lMirrEle.begin(); emirr!=lMirrEle.end(); emirr++ ) {
	  if( (*emirr)->name() == mEName[ke] ) {
	    vCePos0.push_back( (*emirr)->vpos() );
	    break;
	  }
	}
      }

      for( int kh=0; kh<100; kh++ ) vRC8500.push_back( NULL );
      for( int kh=0; kh<100; kh++ ) vRC8600.push_back( NULL );
      for( int kh=0; kh<100; kh++ ) vRC8700.push_back( NULL );
      if( CsRCHistos::Ref().bookHis() && CsRCHistos::Ref().levelBk() >= 2 ) {
	CsHistograms::SetCurrentPath("/RICH");
//----- 8500
	vRC8500.clear();
        for( int kh=0; kh<50; kh++ ) {
	  string hTitle = " ";
	  stringstream hN8500;
	  int kHist = kOffHi + 8500 + kh + 1;
	  hN8500 << kHist;
	  vRC8500.push_back( new CsHist2D( hN8500.str(), hTitle,
					   120, 30., 60., 180, 0., 360.) );
	}
	for( int kh=50; kh<75; kh++ ) {
	  string hTitle = " ";
	  stringstream hN8500;
	  int kHist = kOffHi + 8500 + kh + 1;
	  hN8500 << kHist;
	  vRC8500.push_back( new CsHist2D( hN8500.str(), hTitle,
					   120, -6., 6., 180, 0., 360.) );
	}
	for( int kh=75; kh<91; kh++ ) {
	  string hTitle = " ";
	  stringstream hN8500;
	  int kHist = kOffHi + 8500 + kh + 1;
	  hN8500 << kHist;
	  vRC8500.push_back( new CsHist2D( hN8500.str(), hTitle,
					   120, 30., 60., 180, 0., 360.) );
	}
	string hTitle = " ";
	stringstream hN8502;      //   dummy
	int kHist;
	kHist = kOffHi + 8500 + 91 + 1;
	hN8502 << kHist;
	vRC8500.push_back( new CsHist2D( hN8502.str(), hTitle,
					 120, 0.0, 0.1, 120, 0.0, 0.1 ) );
	stringstream hN8503;
	kHist = kOffHi + 8500 + 92 + 1;
	hN8503 << kHist;
	vRC8500.push_back( new CsHist2D( hN8503.str(), hTitle,
					 120, -0.06,  0.24,
					 120, -0.40, -0.10 ) );
	stringstream hN8504;
	kHist = kOffHi + 8500 + 93 + 1;
	hN8504 << kHist;
	vRC8500.push_back( new CsHist2D( hN8504.str(), hTitle,
					 120, -0.24,  0.06,
					 120, -0.40, -0.10 ) );
	stringstream hN8505;
	kHist = kOffHi + 8500 + 94 + 1;
	hN8505 << kHist;
	vRC8500.push_back( new CsHist2D( hN8505.str(), hTitle,
					 120, -0.06,  0.24,
					 120,  0.10,  0.40 ) );
	stringstream hN8506;
	kHist = kOffHi + 8500 + 95 + 1;
	hN8506 << kHist;
	vRC8500.push_back( new CsHist2D( hN8506.str(), hTitle,
					 120, -0.24,  0.06,
					 120,  0.10,  0.40 ) );
	for( int kh=96; kh<100; kh++ ) {
	  string hTitle = " ";
	  stringstream hN8500;
	  int kHist = kOffHi + 8500 + kh + 1;
	  hN8500 << kHist;
	  vRC8500.push_back( new CsHist2D( hN8500.str(), hTitle,
					   120, -0.06, 0.24,
					   120, -0.40, -0.10 ) );
	}
//----- 8600
	vRC8600.clear();
	int nCath = dets->nCathode();
	CsRCCathode* cat = dets->ptrToCat( 3 );      //   !
	double hCatx = cat->hCatx();
	double hCaty = cat->hCaty();
	Hep3Vector vOffCatWw;
	vOffCatWw = dets->vOffCatW( dets->nCatPMT()[1] );      //   !
	xLmn = vOffCatWw.x() - hCatx;
	yLmx = vOffCatWw.y() + hCaty;
	vOffCatWw = dets->vOffCatW( dets->nCatPMT()[2] );      //   !
	xLmx = vOffCatWw.x() + hCatx;
	yLmn = vOffCatWw.y() - hCaty;
	hTitle = " ";
	stringstream hN8601;
	kHist = kOffHi + 8601;
	hN8601 << kHist;
	vRC8600.push_back( new CsHist2D( hN8601.str(), hTitle,
					 360, -900., 900.,
					 360, -900., 900.) );
	stringstream hN8602;
	kHist = kOffHi + 8602;
	hN8602 << kHist;
	vRC8600.push_back( new CsHist2D( hN8602.str(), hTitle,
					 160, -200., 200.,
					 160, -200., 200.) );
	stringstream hN8603;
	kHist = kOffHi + 8603;
	hN8603 << kHist;
	vRC8600.push_back( new CsHist2D( hN8603.str(), hTitle,
					 180, 0., 360.,
					 100, 0., 100.) );
	stringstream hN8604;
	kHist = kOffHi + 8604;
	hN8604 << kHist;
	vRC8600.push_back( new CsHist2D( hN8604.str(), hTitle,
					 200, 100., 300.,
					 200, 400., 600.) );
	stringstream hN8605;
	kHist = kOffHi + 8605;
	hN8605 << kHist;
	vRC8600.push_back( new CsHist2D( hN8605.str(), hTitle,
					 160, -200., 200.,
					 160, -200., 200.) );
	for( int kh=5; kh<7; kh++ ) {
	  string hTitle = " ";
	  stringstream hN8600;
	  int kHist = kOffHi + 8600 + kh + 1;
	  hN8600 << kHist;
	  vRC8600.push_back( new CsHist2D( hN8600.str(), hTitle,
					   100, -25., 25.,
					   100, -25., 25.) );
	}
	for( int kh=0; kh<2; kh++ ) {
	  int khh = 7 + kh;
	  string hTitle = " ";
	  stringstream hN8600;
	  int kHist = kOffHi + 8600 + khh + 1;
	  hN8600 << kHist;
	  vRC8600.push_back( new CsHist2D( hN8600.str(), hTitle,
					   200, 0., 200.,
					   1, 0., 1.) );
	}
	stringstream hN8610;
	kHist = kOffHi + 8610;
	hN8610 << kHist;
	vRC8600.push_back( new CsHist2D( hN8610.str(), hTitle,
					 180, -360., 360.,
					 180, 0., 360.) );
	for( int kh=0; kh<8; kh++ ) {
	  int khh = 10 + kh;
	  string hTitle = " ";
	  stringstream hN8600;
	  int kHist = kOffHi + 8600 + khh + 1;
	  hN8600 << kHist;
	  vRC8600.push_back( new CsHist2D( hN8600.str(), hTitle,
					   100, -25., 25.,
					   100, -25., 25.) );
					   //100, -5., 5., 100, -5., 5.) );
	}
	for( int kh=0; kh<4; kh++ ) {
	  int khh = 18 + kh;
	  string hTitle = " ";
	  stringstream hN8600;
	  int kHist = kOffHi + 8600 + khh + 1;
	  hN8600 << kHist;
	  vRC8600.push_back( new CsHist2D( hN8600.str(), hTitle,
					   100, -5., 5., 100, -5., 5.) );
	}
	for( int kh=0; kh<4; kh++ ) {
	  int khh = 22 + kh;
	  string hTitle = " ";
	  stringstream hN8600;
	  int kHist = kOffHi + 8600 + khh + 1;
	  hN8600 << kHist;
	  vRC8600.push_back( new CsHist2D( hN8600.str(), hTitle,
					   //180, -360., 360.,
					   180, 0., 360.,
					   180, 0., 360.) );
	}
	for( int kh=0; kh<4; kh++ ) {
	  int khh = 26 + kh;
	  string hTitle = " ";
	  stringstream hN8600;
	  int kHist = kOffHi + 8600 + khh + 1;
	  hN8600 << kHist;
	  vRC8600.push_back( new CsHist2D( hN8600.str(), hTitle,
					   180, -360., 360.,
					   180, 0., 360.) );
	}
	for( int kh=0; kh<4; kh++ ) {
	  int khh = 30 + kh;
	  string hTitle = " ";
	  stringstream hN8600;
	  int kHist = kOffHi + 8600 + khh + 1;
	  hN8600 << kHist;
	  vRC8600.push_back( new CsHist2D( hN8600.str(), hTitle,
					   100, 0., 300.,
					   100, 0., 360.) );
	}
	for( int kh=0; kh<4; kh++ ) {
	  int khh = 34 + kh;
	  string hTitle = " ";
	  stringstream hN8600;
	  int kHist = kOffHi + 8600 + khh + 1;
	  hN8600 << kHist;
	  if( kh == 0 ) {
	    vRC8600.push_back( new CsHist2D( hN8600.str(), hTitle,
					     150, - 200., 1300.,
					     100, - 800.,  200.) );
	  }
	  if( kh == 1 ) {
	    vRC8600.push_back( new CsHist2D( hN8600.str(), hTitle,
					     150, -1300.,  200.,
					     100, - 800.,  200.) );
	  }
	  if( kh == 2 ) {
	    vRC8600.push_back( new CsHist2D( hN8600.str(), hTitle,
					     150, - 200., 1300.,
					     100, - 200.,  800.) );
	  }
	  if( kh == 3 ) {
	    vRC8600.push_back( new CsHist2D( hN8600.str(), hTitle,
					     150, -1300.,  200.,
					     100, - 200.,  800.) );
	  }
	}
	for( int kh=0; kh<4; kh++ ) {
	  int khh = 38 + kh;
	  string hTitle = " ";
	  stringstream hN8600;
	  int kHist = kOffHi + 8600 + khh + 1;
	  hN8600 << kHist;
	  vRC8600.push_back( new CsHist2D( hN8600.str(), hTitle,
					   180, -90., 270.,
					   1, 0., 1.) );
	}
	for( int kh=0; kh<4; kh++ ) {
	  int khh = 42 + kh;
	  string hTitle = " ";
	  stringstream hN8600;
	  int kHist = kOffHi + 8600 + khh + 1;
	  hN8600 << kHist;
	  if( kh == 0 ) {
	    vRC8600.push_back( new CsHist2D( hN8600.str(), hTitle,
					     100,     0., 1000.,
					     100,     0., 1000.) );
	  }
	  if( kh == 1 ) {
	    vRC8600.push_back( new CsHist2D( hN8600.str(), hTitle,
					     100, -1000.,    0.,
					     100,     0., 1000.) );
	  }
	  if( kh == 2 ) {
	    vRC8600.push_back( new CsHist2D( hN8600.str(), hTitle,
					     100,     0., 1000.,
					     100, -1000.,    0.) );
	  }
	  if( kh == 3 ) {
	    vRC8600.push_back( new CsHist2D( hN8600.str(), hTitle,
					     100, -1000.,    0.,
					     100, -1000.,    0.) );
	  }
	}
	for( int kh=0; kh<4; kh++ ) {
	  int khh = 46 + kh;
	  string hTitle = " ";
	  stringstream hN8600;
	  int kHist = kOffHi + 8600 + khh + 1;
	  hN8600 << kHist;
	  vRC8600.push_back( new CsHist2D( hN8600.str(), hTitle,
					   100, 0., 100.,
					   180, 0., 360.) );
	}

	for( int kh=0; kh<25; kh++ ) {
	  int khh = 50 + kh;
	  string hTitle = " ";
	  stringstream hN8600;
	  int kHist = kOffHi + 8600 + khh + 1;
	  hN8600 << kHist;
	  vRC8600.push_back( new CsHist2D( hN8600.str(), hTitle,
					   100, -25., 25.,
					   100, -25., 25.) );
	}
	for( int kh=0; kh<16; kh++ ) {
	  int khh = 75 + kh;
	  string hTitle = " ";
	  stringstream hN8600;
	  int kHist = kOffHi + 8600 + khh + 1;
	  hN8600 << kHist;
	  vRC8600.push_back( new CsHist2D( hN8600.str(), hTitle,
					   100, -25., 25.,
					   100, -25., 25.) );
	}
	for( int kh=0; kh<2; kh++ ) {
	  int khh = 91 + kh;
	  string hTitle = " ";
	  stringstream hN8600;
	  int kHist = kOffHi + 8600 + khh + 1;
	  hN8600 << kHist;
	  vRC8600.push_back( new CsHist2D( hN8600.str(), hTitle,
					   100, -5., 5.,
					   1, 0., 1.) );
	}
	for( int kh=0; kh<5; kh++ ) {
	  int khh = 93 + kh;
	  string hTitle = " ";
	  stringstream hN8600;
	  int kHist = kOffHi + 8600 + khh + 1;
	  hN8600 << kHist;
	  vRC8600.push_back( new CsHist2D( hN8600.str(), hTitle,
					   220, -2200., 2200.,
					   220, -2200., 2200. ) );
	}
//----- 8700
	vRC8700.clear();
	for( int kh=0; kh<16; kh++ ) {
	  int khh = 0 + kh;
	  string hTitle = " ";
	  stringstream hN8700;
	  int kHist = kOffHi + 8700 + khh + 1;
	  hN8700 << kHist;
	  vRC8700.push_back( new CsHist2D( hN8700.str(), hTitle,
					   100, -5., 5.,
                                           180, 0., 360.) );
	}
	for( int kh=0; kh<16; kh++ ) {
	  int khh = 16 + kh;
	  string hTitle = " ";
	  stringstream hN8700;
	  int kHist = kOffHi + 8700 + khh + 1;
	  hN8700 << kHist;
	  vRC8700.push_back( new CsHist2D( hN8700.str(), hTitle,
					   100, -5., 5.,
                                           180, 0., 360.) );
	}
//----- 8733-37:   NEW CLUSTER MAP ?
	for( int kh=0; kh<5; kh++ ) {
	  int khh = 32 + kh;
	  string hTitle = " ";
	  stringstream hN8700;
	  int kHist = kOffHi + 8700 + khh + 1;
	  hN8700 << kHist;
	  vRC8700.push_back( new CsHist2D( hN8700.str(), hTitle,
					   48, 0., 48.,
					   48, 0., 48. ) );
	}

	CsHistograms::SetCurrentPath("/");
      }
      std::cout << " RICHONE, CsRCEventAnalysis::checkOptCorr : histo buffer "
		<< vRC8500.size() << " + " << vRC8600.size() << " + "
		<< vRC8700.size() << std::endl;

//    8694
      for( emirr=lMirrEle.begin(); emirr!=lMirrEle.end(); emirr++ ) {
	int khh = 93;
	double xCePos = (*emirr)->vpos().x();
	double yCePos = (*emirr)->vpos().y();
	if( khh >= 0  &&  khh < int( vRC8600.size() ) ) {
	  //std::cout << xCePos << "  " << yCePos << std::endl;
	  if( vRC8600[khh] ) vRC8600[khh]->Fill( xCePos, yCePos );
//hh                         ------------------------------------
	}
      }
    }

//- NEW CLUSTER MAP (?)
    list<CsRCCathode*> lCathodes = dets->lCathodes();
    list<CsRCCathode*>::iterator ict;
    CsRCEventClusters* clus = CsRCEventClusters::Instance();
    list<CsRCCluster*> lClusters = clus->lClusters();
    list<CsRCCluster*>::iterator icl;
    for( icl=lClusters.begin(); icl!=lClusters.end(); icl++ ) {
      if( !(*icl)->flag() ) continue;
      CsRCCluster* clu = (*icl);

      list<CsRCPad*> lPads;
      CsRCPad* pad = clu->lPads().front();
      int ipx = pad->ix();
      int ipy = pad->iy();
      int iCat = pad->ic();
      int kQua = -1;
      if( iCat >=  0  &&  iCat <=  3 ) kQua = 0;
      if( iCat >=  4  &&  iCat <=  7 ) kQua = 1;
      if( iCat >=  8  &&  iCat <= 11 ) kQua = 2;
      if( iCat >= 12  &&  iCat <= 15 ) kQua = 3;
      if( kQua < 0 ) continue;
      CsRCCathode* cat = dets->ptrToCat( iCat );
      int kx = 0;
      int ky = 0;
      if( cat->isPMT() ) { kx = ipx/2; ky = ipy/2; }
      else { kx = ipx/3; ky = ipy/3; }
      int nChx = 24;
      int nChy = 24;
      int khh = 32 + kQua;
      int kxx = kx;
      int kyy = ky;
      if( iCat%4 ==  0 ) { kxx += nChx; kyy += nChy; }
      if( iCat%4 ==  1 ) { kxx += nChx; }
      if( iCat%4 ==  2 ) { kyy += nChy; }
      if( iCat%4 ==  3 ) { }
      xh = double( kxx );
      yh = double( kyy );
//--- hist # 8733-36
      if( khh >= 0  &&  khh < int( vRC8700.size() ) ) {
	if( vRC8700[khh] ) vRC8700[khh]->Fill( xh, yh );
//hh                       ----------------------------
      }
      /*
      int iCat = clu->ic();
      CsRCCathode* cat = dets->ptrToCat( iCat );
      //int kClu = clu->kClu();
      double xc = clu->xc();
      double yc = clu->yc();
      double hCatx = cat->hCatx();
      double xlCatn = dets->vOffCatW( iCat ).x() - hCatx;
      double xlCatx = dets->vOffCatW( iCat ).x() + hCatx;
      //double padx = cat->padx();
      int nChx = 0;
      if( cat->isPMT() ) nChx = cat->nPadx()/2;
      else  nChx = cat->nPadx()/3;
      double ddx = (xlCatx - xlCatn)/nChx;
      double hCaty = cat->hCaty();
      double ylCatn = dets->vOffCatW( iCat ).y() - hCaty;
      double ylCatx = dets->vOffCatW( iCat ).y() + hCaty;
      //double pady = cat->pady();
      int nChy = 0;
      if( cat->isPMT() ) nChy = cat->nPady()/2;
      else  nChy = cat->nPady()/3;
      double ddy = (ylCatx - ylCatn)/nChy;
      int kQua = -1;
      if( iCat >=  0  &&  iCat <=  3 ) kQua = 0;
      if( iCat >=  4  &&  iCat <=  7 ) kQua = 1;
      if( iCat >=  8  &&  iCat <= 11 ) kQua = 2;
      if( iCat >= 12  &&  iCat <= 15 ) kQua = 3;
      if( kQua < 0 ) continue;
      bool fill = false;
      for( int kx=0; kx<nChx; kx++ ) {
	double xln = xlCatn + kx*ddx;
	double xlx = xlCatn + (kx+1)*ddx;
	for( int ky=0; ky<nChy; ky++ ) {
	  double yln = ylCatn + ky*ddy;
	  double ylx = ylCatn + (ky+1)*ddy;
	  if( xc >= xln  &&  xc < xlx  &&  yc >= yln  &&  yc < ylx ) {
	    int khh = 32 + kQua;
	    int kxx = kx;
	    int kyy = ky;
	    if( iCat%4 ==  0 ) { kxx += nChx; kyy += nChy; }
	    if( iCat%4 ==  1 ) { kxx += nChx; }
	    if( iCat%4 ==  2 ) { kyy += nChy; }
	    if( iCat%4 ==  3 ) { }
	    xh = double( kxx );
	    yh = double( kyy );
	    if( khh >= 0  &&  khh < int( vRC8700.size() ) ) {
	      if( vRC8700[khh] ) vRC8700[khh]->Fill( xh, yh );
//hh                             ----------------------------
	    }
	    fill = true;
	  }
	  if( fill ) break;
	}
	if( fill ) break;
      }
      */

    }
    khh = 32 + 4;
    xh = double( 1 );
    yh = double( 1 );
    if( khh >= 0  &&  khh < int( vRC8700.size() ) ) {
      if( vRC8700[khh] ) vRC8700[khh]->Fill( xh, yh );
//hh                     ----------------------------
    }


    static const double DegRad = 0.0174533;
    static const int ncTxy = 7;
    static const double dDirx = 0.05;
    static const double xDirLmn = -0.06 - dDirx/2.;
    static const double xDirLmx =  0.24 + dDirx/2.;
    static const double dDiry = 0.05;
    static const double yDirLmn = -0.40 - dDiry/2.;
    static const double yDirLmx = -0.10 + dDiry/2.;

    CsRCEventPartPhotons* evepapho = CsRCEventPartPhotons::Instance();
    list<CsRCPartPhotons*> lPaPhotons = evepapho->lPartPhotons();
    list<CsRCPartPhotons*>::iterator ipf;
    for( ipf=lPaPhotons.begin(); ipf!=lPaPhotons.end(); ipf++ ) {
      if( !(*ipf) ) continue;
      CsRCPartPhotons* papho = (*ipf);
      if( !papho->flag() ) continue;

      CsRCRing* ring = papho->pRing();
      bool RING = false;
      if( ring  &&  ring->flag() ) RING = true;
      if( !RING ) continue;
//    --------------------

//--- Count Mirror Elements
      Hep3Vector vPoPaMir0 =
	papho->pPart()->vPoPaMir0()[papho->pPart()->kDetPart()];
      double rrElmQ = cons->rrAliMir() * cons->rrAliMir();
      int eMirrNb = 0;
      std::string eName = " ";
      std::list<CsRCMirrorElem*> lMirrEle = mirr->lMirrEle();
      std::list<CsRCMirrorElem*>::iterator emirr;
      std::list<CsRCMirrorElem*>::iterator emirrP;
      for( emirr=lMirrEle.begin(); emirr!=lMirrEle.end(); emirr++ ) {
	if ( (vPoPaMir0 - (*emirr)->vpos()).mag2() < rrElmQ ) {
	  eName = (*emirr)->name();
	  std::string nb = eName.substr( 2, 2 );
	  eMirrNb = atoi( nb.c_str() );
	  if( eName[1] == 'T' ) eMirrNb += 100;
	  emirrP = emirr;
	  break;
	}
      }
//    8608
      khh = 7;
      if( eMirrNb > 0 ) {
	xh = double( eMirrNb ) + 0.5;
	yh = 0.5;
	if( khh >= 0  &&  khh < int( vRC8600.size() ) ) {
	  if( vRC8600[khh] ) vRC8600[khh]->Fill( xh, yh );
//hh                         ----------------------------
	}
      }
//    8695
      khh = 94;
      if( eMirrNb > 0 ) {
	xh = (*emirrP)->vpos().x();
	yh = (*emirrP)->vpos().y();
	if( khh >= 0  &&  khh < int( vRC8600.size() ) ) {
	  if( vRC8600[khh] ) vRC8600[khh]->Fill( xh, yh );
//hh                         ----------------------------
	}
      }
//    8697
      khh = 96;
      if( papho->pPart()->kDetPart() == 1 ) khh = 97;
      xh = vPoPaMir0.x();
      yh = vPoPaMir0.y();
      if( khh >= 0  &&  khh < int( vRC8600.size() ) ) {
	if( vRC8600[khh] ) vRC8600[khh]->Fill( xh, yh );
//hh                       ----------------------------
      }

//--- Select rings fully contained in the 4 PMT cathodes :
//    ----------------------------------------------------
      int kDetPart = papho->kDetPart();
      double xPade = papho->vPoPaDetW()[kDetPart].x();
      double yPade = papho->vPoPaDetW()[kDetPart].y();
      double RCing =  ring->theReco() * 3.3;
      bool bCRing4 = false;
      //std::cout << xLmn << "  " << xPade << "  " << xLmx << " - "
      //	  << yLmn << "  " << yPade << "  " << yLmx << std::endl;
      if( (xPade - xLmn) > RCing  &&  (xLmx - xPade) > RCing  &&
	  (yLmx - yPade) > RCing  &&  (yPade - yLmn) > RCing ) bCRing4 = true;
      //std::cout << RCing << "  " << bCRing4 << "  " << true << std::endl;
      //if( !bCRing4 ) continue;
//    -----------------------
//--- Select rings fully contained in each PMT cathode :
//    --------------------------------------------------
      CsRCCathode *cat = dets->ptrToCat( dets->nCatPMT()[0] );      //   !
      double hCatx = cat->hCatx();
      double hCaty = cat->hCaty();
      bool bCRing1 = false;
      if( (xPade - xLmn) > RCing  &&  ((xLmn+2.*hCatx) - xPade) > RCing  &&
	  (yLmx - yPade) > RCing  &&  (yPade - (yLmx-2.*hCaty)) > RCing )
	bCRing1 = true;
      //std::cout << xLmn << "  " << xLmn+2.*hCatx << " - "
      //	  << yLmx << "  " << yLmx-2.*hCaty << std::endl;
      if( (xPade - (xLmx-2.*hCatx)) > RCing  &&  (xLmx - xPade) > RCing  &&
	  (yLmx - yPade) > RCing  &&  (yPade - (yLmx-2.*hCaty)) > RCing )
	bCRing1 = true;
      if( (xPade - xLmn) > RCing  &&  ((xLmn+2.*hCatx) - xPade) > RCing  &&
	  (yPade - yLmn) > RCing  &&  ((yLmn+2.*hCaty) - yPade) > RCing )
	bCRing1 = true;
      if( (xPade - (xLmx-2.*hCatx)) > RCing  &&  (xLmx - xPade) > RCing  &&
	  (yPade - yLmn) > RCing  &&  ((yLmn+2.*hCaty) - yPade) > RCing )
	bCRing1 = true;
      if( !bCRing1 ) continue;      ////
//    -----------------------

//--- Cut on low particle momentum :
      CsRCParticle* part = papho->pPart();
      double momPart = part->mom();
      if( momPart < 3. ) continue;
//    ---------------------------

//--- Select ID pions (analysis way)
      bool bIDpion = false;
      double likeDV = cons->likeDefVa();
      double* probaLK;
      probaLK = papho->probaLKAll();
      double likeBkg = papho->probaLKBgAll();
      double likeElec = probaLK[ 2];
      double likeMuon = probaLK[ 5];
      double likePion = probaLK[ 8];
      double likeKaon = probaLK[11];
      double likeProton = probaLK[14];
      if( likeElec == likeBkg ) continue;   // NO photons from signal
      if( likePion > likeDV  &&  likePion > likeBkg  &&
	  likePion > likeKaon  &&  likePion > likeProton ) bIDpion = true;
      if( momPart < 8.  &&  likePion < likeElec ) bIDpion = false;
      if( !bIDpion ) continue;
//    -----------------------

/*
//--- Count Mirror Elements
      Hep3Vector vPoPaMir0 =
	papho->pPart()->vPoPaMir0()[papho->pPart()->kDetPart()];
      double rrElmQ = cons->rrAliMir() * cons->rrAliMir();
      int eMirrNb = 0;
      std::string eName = " ";
      std::list<CsRCMirrorElem*> lMirrEle = mirr->lMirrEle();
      std::list<CsRCMirrorElem*>::iterator emirr;
      std::list<CsRCMirrorElem*>::iterator emirrP;
      for( emirr=lMirrEle.begin(); emirr!=lMirrEle.end(); emirr++ ) {
	if ( (vPoPaMir0 - (*emirr)->vpos()).mag2() < rrElmQ ) {
	  eName = (*emirr)->name();
	  std::string nb = eName.substr( 2, 2 );
	  eMirrNb = atoi( nb.c_str() );
	  if( eName[1] == 'T' ) eMirrNb += 100;
	  emirrP = emirr;
	  break;
	}
      }
//    8608
      khh = 7;
      if( eMirrNb > 0 ) {
	xh = double( eMirrNb ) + 0.5;
	yh = 0.5;
	if( khh >= 0  &&  khh < int( vRC8600.size() ) ) {
	  if( vRC8600[khh] ) vRC8600[khh]->Fill( xh, yh );
//hh                         ----------------------------
	}
      }
//    8695
      khh = 94;
      if( eMirrNb > 0 ) {
	xh = (*emirrP)->vpos().x();
	yh = (*emirrP)->vpos().y();
	if( khh >= 0  &&  khh < int( vRC8600.size() ) ) {
	  if( vRC8600[khh] ) vRC8600[khh]->Fill( xh, yh );
//hh                         ----------------------------
	}
      }
*/

//--- Select only 'RINGs' fully contained inside one mirror element
//--- Select only 4 mirror elements (one in each quadrant) (MRS)
      bool oneMirr = false;
      int kev = -1;
      for( int ke=0; ke<4; ke++ ) {
	if( (vPoPaMir0 - vCePos0[ke]).mag2() < rrElmQ ) {
	  //std::cout << vCePos0[ke] << std::endl;
	  oneMirr = true;
	  kev = ke;
	  break;
	}
      }
      if( !oneMirr ) continue;      ////
//    -----------------------
//    8609
      khh = 8;
      if( kev >= 0 ) {
	std::string nb = mEName[kev].substr( 2, 2 );
	int eMirrNb = atoi( nb.c_str() );
	if( mEName[kev][1] == 'T' ) eMirrNb += 100;
	xh = double( eMirrNb ) + 0.5;
	yh = 0.5;
	if( khh >= 0  &&  khh < int( vRC8600.size() ) ) {
	  if( vRC8600[khh] ) vRC8600[khh]->Fill( xh, yh );
//hh                         ----------------------------
	}
      }
//    8696
      khh = 95;
      if( eMirrNb > 0 ) {
	xh = (*emirrP)->vpos().x();
	yh = (*emirrP)->vpos().y();
	if( khh >= 0  &&  khh < int( vRC8600.size() ) ) {
	  if( vRC8600[khh] ) vRC8600[khh]->Fill( xh, yh );
//hh                         ----------------------------
	}
      }


//--- PARTICLE geometry
//    -----------------
      double TwoPI = cons->TwoPI();
      double RadDeg = cons->RadDeg();
      Hep3Vector vDcPartW = papho->vDcPartW()[kDetPart];
      double thePart = acos( vDcPartW.z() );
      double phiPart = 
	CsRCUty::Instance()->psatg( vDcPartW.y(), vDcPartW.x() );
//      --------------------------------------------------------
      thePart *= 1000.;          //   mrad
      phiPart *= RadDeg;         //   deg
      Hep3Vector vPoPaMirW = papho->vPoPaMirW()[kDetPart];
      Hep3Vector vcrossPaMirW = vDcPartW.cross( vPoPaMirW );
//    Particle plane
      double phiPaPl = 
	CsRCUty::Instance()->psatg( vcrossPaMirW.y(), vcrossPaMirW.x() );
//      ----------------------------------------------------------------
      phiPaPl *= RadDeg;         //   deg
      if( kDetPart == 1 ) phiPaPl -= 180.;
//    -----------------------------------
      if( phiPaPl < 0. ) phiPaPl += 360.;
      Hep3Vector vPade = papho->pPart()->vPade()[kDetPart];

      for(int kh=0; kh<4; kh++ ) {
	if( papho->iCaPa() == dets->nCatPMT()[kh] ) {
//        8631-34
	  khh = 30 + kh;
	  xh = thePart;
	  yh = phiPart;
	  if( khh >= 0  &&  khh < int( vRC8600.size() ) ) {
	    if( vRC8600[khh] ) vRC8600[khh]->Fill( xh, yh );
//hh                           ----------------------------
	  }
//        8635-38
	  khh = 34 + kh;
	  xh = vPoPaMirW.x();
	  yh = vPoPaMirW.y();
	  if( khh >= 0  &&  khh < int( vRC8600.size() ) ) {
	    if( vRC8600[khh] ) vRC8600[khh]->Fill( xh, yh );
//hh                           ----------------------------
	  }
//        8639-42
	  khh = 38 + kh;
	  xh = phiPaPl;
	  yh = 0.5;
	  if( khh >= 0  &&  khh < int( vRC8600.size() ) ) {
	    if( vRC8600[khh] ) vRC8600[khh]->Fill( xh, yh );
//hh                           ----------------------------
	  }
//        8643-46
	  khh = 42 + kh;
	  xh = vPade.x();
	  yh = vPade.y();
	  if( khh >= 0  &&  khh < int( vRC8600.size() ) ) {
	    if( vRC8600[khh] ) vRC8600[khh]->Fill( xh, yh );
//hh                           ----------------------------
	  }
	}
      }


//--- PART-PHOT photons
//    -----------------
      //list<CsRCPhoton*> lPhotons = papho->lPhotons();

//--- RING photons
//    ------------
      list<CsRCPhoton*> lPhotons = ring->lPhotons();
      int nPhotons = lPhotons.size();
//--- Select only RINGs with many photons
      if( nPhotons < 20 ) continue;      ////
//    ----------------------------
      list<CsRCPhoton*>::iterator ih;
      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	if( !(*ih) ) continue;
	CsRCPhoton* pho = (*ih);

	if( !pho->isPMT() ) continue;
//      ----------------------------
	CsRCCluster* clu = pho->ptToClu();
	if( !clu ) continue;
	int iCat = clu->ic();
        double xClu = clu->xc();
	double yClu = clu->yc();
	int kClu = clu->kClu();
	list<CsRCPad*> lPads = clu->lPads();
	CsRCPad* pad = lPads.front();
	if( !pad ) continue;
	int kDetClu = pho->kDetClu();
	double xw = clu->xc();
	double yw = clu->yc();
	double zw = dets->vOffCatW( iCat ).z();
        Hep3Vector vPoCluDet( xw, yw, zw );
	//std::cout << "zDetW = " << zw << "  " << iCat << std::endl;
	//- zw = 3461.54 ->cat 3, 5 - zw = 3452.72->cat 10, 12   ???
	//zw = zw + 10.;
	double zDetW = zw;

//----- Photon phi redefinition
//      -----------------------
//----- phi to particle-plane
	Hep3Vector vPoPhotW = papho->vPoPhotW()[kDetClu];
	double RR = mirr->RRv( kDetClu );
	double cosPaMir = vDcPartW * vPoPaMirW / RR;
	if( cosPaMir >  1.) cosPaMir =  1.;
	if( cosPaMir < -1.) cosPaMir = -1.;
	double thePaMir = acos( cosPaMir );
	bool flag = true;
        Hep3Vector vDcPhoEmW = 
	  papho->getCluAngle( vPoPhotW, vPoCluDet, RR, flag );
//        ---------------------------------------------------
	Hep3Vector vcrossPhoPaW = vDcPhoEmW.cross( vDcPartW );
	double cosPhi = vcrossPhoPaW * vcrossPaMirW;
	///double theClu = pho->theNorm()/1000.;
	double theClu = acos( vDcPhoEmW.dot( vDcPartW ) );
	double cosThClu = cos( theClu );
	double norm = (1.-cosThClu*cosThClu) * (1.-cosPaMir*cosPaMir);
	double phiClu = 0.;
	if( norm > 0. ) {
	  norm = sqrt( norm ) * RR;
	  cosPhi /= norm;
	  if( cosPhi >  1.) cosPhi =  1.;
	  if( cosPhi < -1.) cosPhi = -1.;
	  phiClu = acos( cosPhi );
	  if( vcrossPhoPaW.y() < 0. ) phiClu = TwoPI - phiClu;
	}  else  {
	  phiClu = 0.;
	}
	theClu *= 1000.;           //   mrad
	phiClu *= RadDeg;          //   deg
//----- phi to YZ-plane
	double theCluA = 0.;
	double phiCluA = 0.;
	Hep3Vector vDcPhoEmPa = CsRCUty::Instance()->
	  rotfbcc( +1., vDcPartW, vDcPhoEmW, theCluA, phiCluA );
//        ------------------------------------------------------
	theCluA *= 1000.;          //   mrad
	phiCluA *= RadDeg;         //   deg
	//if( kDetClu == 1 ) phiCluA -= 180.;
//      ---------------------------------
	if( phiCluA < 0. ) phiCluA += 360.;

//----- Check emiss. point
	bool checkp = true;
	std::vector<double> thep;
	if( checkp ) {
	  for( int ke=0; ke<4; ke++ ) {
	    double fr = -1200. + ke * 800.;
	    Hep3Vector vPoW = vPoPhotW + fr * vDcPartW;
	    Hep3Vector vDcW = papho->getCluAngle( vPoW, vPoCluDet, RR, flag );
	    thep.push_back( 1000. * acos( vDcW.dot( vDcPartW ) ) );
	  }
	  //std::cout << "checkp =  ";
	  //for( int ke=0; ke<4; ke++ ) std::cout << thep[ke] << "  ";
	  //std::cout << theClu << std::endl;
	}


//----- Photon angles from ABOVE = only ONE FULL mirror UP or DOWN used
//      ------------------------
	//double thePho = theCluA;        //   mrad
	//double phiPho = phiCluA;        //   deg
//----- Photon angles from RICHONE = all single mirror ELEMENTS used
//      --------------------------
	//double thePho = pho->theNorm();      //   mrad
	//double thePho = pho->the();          //   mrad
	double thePho = pho->the0();         //   mrad
	double phiPho = pho->phiA();         //   deg
	//if( kDetClu == 1 ) phiPho -= 180.;
//      ---------------------------------
	//std::cout << theCluA << "  " << phiCluA << " - " << pho->the0()
	//          << "  " << pho->phiA() << std::endl;
	if( phiPho < 0. ) phiPho += 360.;
	double phiPhow = pho->phi();
	if( kDetClu == 1 ) phiPhow -= 180.;
//      ----------------------------------
	if( phiPhow < 0. ) phiPhow += 360.;
	if( kDetClu == 1 ) phiPhow = 360. - phiPhow;
//      -------------------------------------------

	double dPhi = phiPho - phiPhow;
	if( dPhi < 0. ) dPhi += 360.;

//----- Cherenkov angle at saturation (assume pion mass)
//      ------------------------------------------------
	double massPi = CsRCRecConst::Instance()->massPartv()[8];
	double beta = momPart / sqrt( massPi*massPi + momPart*momPart );
	double thePhoSa = 1000.* acos( beta * cos( thePho/1000. ) );
	//std::cout << thePho << "  " << thePhoSa << std::endl;
    	double thePhoSaPi = 1000.* acos( 1./cons->CFRefIndVS() );
	//std::cout << thePhoSa << "  " << thePhoSaPi << std::endl;

//----- Check refr. index  (lam -> 200-700 um)
	bool checkn = true;
	std::vector<double> then;
	if( checkn ) {
	  for( int kn=0; kn<4; kn++ ) {
	    double lam = 240. + kn * 140.;
	    double ri = 1. + (0.2375 / (0.0001845 - 1./(lam*lam))) / 1000000.;
	    then.push_back( 1000.* acos( 1./ri ) );
	  }
	  //std::cout << "checkn =  ";
	  //for( int ke=0; ke<4; ke++ ) std::cout << then[ke] << "  ";
	  //std::cout << thePhoSaPi << std::endl;
	}

//----- Photon recostruction (MWR) (100420)
//      --------------------------
	double ll = sin( thePho/1000. ) * cos( phiPho*DegRad );
	double mm = sin( thePho/1000. ) * sin( phiPho*DegRad );
	double nn = cos( thePho/1000. );
	Hep3Vector vDcPhoEmP( ll, mm, nn );
	vDcPhoEmP = vDcPhoEmP.unit();
	double the = 0.;
	double phi = 0.;
	//Hep3Vector vDcPhoEmW;
	vDcPhoEmW = CsRCUty::Instance()->
	  rotfbcc( -1., vDcPartW, vDcPhoEmP, the, phi );
//        ---------------------------------------------
	vDcPhoEmW = vDcPhoEmW.unit();
	//std::cout << vDcPartW << "  " << vDcPhoEmP << "  " << vDcPhoEmW
	//	    << std::endl;
//----- Photon impact on mirror (only ONE FULL mirror UP or DOWN used)
//                               ------------------------------------
	Hep3Vector vPoC( 0., 0., 0. );
	Hep3Vector vPoPhoMir = mirr->vImpMir( vPoPhotW, vDcPhoEmW, vPoC, RR);
//                             ---------------------------------------------
//----- Normal to mirror at photon impact
	Hep3Vector vDcNoPhoMir = (1./RR) * vPoPhoMir;
//----- Direction of reflected photon
	double cosPhoMir = vDcNoPhoMir * vDcPhoEmW;
	Hep3Vector vDcPhoRefl = 2.*cosPhoMir * vDcNoPhoMir - vDcPhoEmW;
	vDcPhoRefl = vDcPhoRefl.unit();
//----- Photon incidence angle on detector
	double cosPhoDet = vDcPhoRefl.z();
	double sinPhoDet = sqrt(1.-cosPhoDet*cosPhoDet);
//----- Photon direction at detector (no QZW correction)
	Hep3Vector vDcPhoDetW = vDcPhoRefl;
//----- Photon impact on detector plane (no QZW correction)
	norm = (zDetW - vPoPhoMir.z()) / vDcPhoRefl.z();
	Hep3Vector vImpDet = vPoPhoMir + norm * vDcPhoRefl;
	double xPho = vImpDet.x();
	double yPho = vImpDet.y();

//----- Quartz-Window correction
	float CFRefInd = cons->CFRefInd();
	float qzRefInd = cons->qzRefInd();
	double CHRefInd = 1.000444;
	double rRatio = CFRefInd / qzRefInd;
	double rRatioQM = qzRefInd / CHRefInd;
	float ddQzW = cat->ddQzW();
	float ddGap = cat->ddGap();
	double zQzWOutside = zDetW + ddGap + ddQzW;
	double zQzWinside  = zDetW + ddGap;
//----- Photon impact on quartz win. (outside, C4F10 side)
	norm = (zQzWOutside - vPoPhoMir.z()) / vDcPhoRefl.z();
	Hep3Vector vPhoWinOut = vPoPhoMir + norm * vDcPhoRefl;
//----- Photon projection on quartz win. (inside, CH4 side)
	///norm = (zQzWinside - vPoPhoMir.z()) / vDcPhoRefl.z();
	///Hep3Vector vPhoWinIn = vPoPhoMir + norm * vDcPhoRefl;
//----- Photon direction inside quartz window
	double sinPhoQzR = sinPhoDet * rRatio;
	double cosPhoQzR = sqrt(1.-sinPhoQzR*sinPhoQzR);
	Hep3Vector vDcPhoQzR;
	vDcPhoQzR.setX( rRatio * vDcPhoRefl.x() );
	vDcPhoQzR.setY( rRatio * vDcPhoRefl.y() );
	vDcPhoQzR.setZ( cosPhoQzR );
	vDcPhoQzR = vDcPhoQzR.unit();
//----- Photon impact on quartz win. (inside, CH4 side)
	norm = - ddQzW / vDcPhoQzR.z();
	Hep3Vector vPhoWinQzR = vPhoWinOut + norm * vDcPhoQzR;
//----- Photon direction outside quartz window (outside, CH4 side)
	double sinPhoGap = sinPhoQzR * rRatioQM;
	double cosPhoGap = sqrt(1.-sinPhoGap*sinPhoGap);
	Hep3Vector vDcPhoGap;
	vDcPhoGap.setX( rRatioQM * vDcPhoQzR.x() );
	vDcPhoGap.setY( rRatioQM * vDcPhoQzR.y() );
	vDcPhoGap.setZ( cosPhoGap );
	vDcPhoGap = vDcPhoGap.unit();
//----- Photon impact on CsI plane
	norm = - ddGap / vDcPhoGap.z();
	Hep3Vector vPhoCsI = vPhoWinQzR + norm * vDcPhoGap;
	double xPhoQz = vPhoCsI.x();
	double yPhoQz = vPhoCsI.y();

	Hep3Vector vDcPhoDetWw = vDcPhoDetW;
	int kQua = -1;
//----- TopJura
	if( iCat == dets->nCatPMT()[0] ) {
	  kQua = 0;
	}
//----- TopSaleve
	else if( iCat == dets->nCatPMT()[1] ) {
	  vDcPhoDetWw.setX( -vDcPhoDetW.x() );
	  kQua = 1;
	}
//----- DownJura
	else if( iCat == dets->nCatPMT()[2] ) {
	  vDcPhoDetWw.setY( -vDcPhoDetW.y() );
	  kQua = 2;
	}
//----- DownSaleve
	else if( iCat == dets->nCatPMT()[3] ) {
	  vDcPhoDetWw.setX( -vDcPhoDetW.x() );
	  vDcPhoDetWw.setY( -vDcPhoDetW.y() );
	  kQua = 3;
	}
	int PMTnum = pad->PMTnum();
	int PMTcha = pad->PMTcha();

//      8601
	khh = 0;
	xh = xClu;
	yh = yClu;
	if( khh >= 0  &&  khh < int( vRC8600.size() ) ) {
	  if( vRC8600[khh] ) vRC8600[khh]->Fill( xh, yh );
//hh                         ----------------------------
	}
	double xPade = papho->pPart()->vPade()[pho->kDetClu()].x();
	double yPade = papho->pPart()->vPade()[pho->kDetClu()].y();
//      8602
	khh = 1;
	xh = xClu - xPade;
	yh = yClu - yPade;
	if( khh >= 0  &&  khh < int( vRC8600.size() ) ) {
	  if( vRC8600[khh] ) vRC8600[khh]->Fill( xh, yh );
//hh                         ----------------------------
	}
//      8603
	khh = 2;
	xh = CsRCUty::Instance()->psatg( yClu-yPade, xClu-xPade )
//           ----------------------------------------------------
	  * cons->RadDeg();
	yh = double( nPhotons );
	if( khh >= 0  &&  khh < int( vRC8600.size() ) ) {
	  if( vRC8600[khh] ) vRC8600[khh]->Fill( xh, yh );
//hh                         ----------------------------
	}
//      8604
	khh = 3;
	xh = xClu;
	yh = yClu;
	if( khh >= 0  &&  khh < int( vRC8600.size() ) ) {
	  if( vRC8600[khh] ) vRC8600[khh]->Fill( xh, yh );
//hh                         ----------------------------
	}
//      8605
	khh = 4;
	xh = xPho - xPade;
	yh = yPho - yPade;
	if( khh >= 0  &&  khh < int( vRC8600.size() ) ) {
	  if( vRC8600[khh] ) vRC8600[khh]->Fill( xh, yh );
//hh                         ----------------------------
	}
//      8606
	khh = 5;
	//std::cout << setprecision( 6 );
	//std::cout << xPho << "  " << xClu << std::endl;
	//std::cout << yPho << "  " << yClu << std::endl;
	xh = xPho - xClu;
	yh = yPho - yClu;
	if( khh >= 0  &&  khh < int( vRC8600.size() ) ) {
	  if( vRC8600[khh] ) vRC8600[khh]->Fill( xh, yh );
//hh                         ----------------------------
	}
//      8607
	khh = 6;
	xh = xPhoQz - xClu;
	yh = yPhoQz - yClu;
	if( khh >= 0  &&  khh < int( vRC8600.size() ) ) {
	  if( vRC8600[khh] ) vRC8600[khh]->Fill( xh, yh );
//hh                         ----------------------------
	}
//      8610
	khh = 9;
	xh = dPhi;
	yh = phiPho;
	if( khh >= 0  &&  khh < int( vRC8600.size() ) ) {
	  if( vRC8600[khh] ) vRC8600[khh]->Fill( xh, yh );
//hh                         ----------------------------
	}
	for(int kh=0; kh<4; kh++ ) {
	  if( iCat == dets->nCatPMT()[kh] ) {
//          8611-14
	    khh = 10 + kh;
	    xh = xPho - xClu;
	    yh = yPho - yClu;
	    if( khh >= 0  &&  khh < int( vRC8600.size() ) ) {
	      if( vRC8600[khh] ) vRC8600[khh]->Fill( xh, yh );
//hh                             ----------------------------
	    }
//          8615-18
    	    khh = 14 + kh;
	    xh = xPhoQz - xClu;
	    yh = yPhoQz - yClu;
	    if( khh >= 0  &&  khh < int( vRC8600.size() ) ) {
	      if( vRC8600[khh] ) vRC8600[khh]->Fill( xh, yh );
//hh                             ----------------------------
	    }
//          8619-22
    	    khh = 18 + kh;
	    xh = xPhoQz - xPho;
	    yh = yPhoQz - yPho;
	    if( khh >= 0  &&  khh < int( vRC8600.size() ) ) {
	      if( vRC8600[khh] ) vRC8600[khh]->Fill( xh, yh );
//hh                             ----------------------------
	    }
//          8623-26
    	    khh = 22 + kh;
	    xh = phiPho;
	    yh = phiPhow;
	    if( khh >= 0  &&  khh < int( vRC8600.size() ) ) {
	      if( vRC8600[khh] ) vRC8600[khh]->Fill( xh, yh );
//hh                             ----------------------------
	    }
//          8627-30
    	    khh = 26 + kh;
	    xh = dPhi;
	    yh = phiPho;
	    if( khh >= 0  &&  khh < int( vRC8600.size() ) ) {
	      if( vRC8600[khh] ) vRC8600[khh]->Fill( xh, yh );
//hh                             ----------------------------
	    }
//          8647-50
	    khh = 46 + kh;
	    xh = thePho;
	    yh = phiPho;
	    if( khh >= 0  &&  khh < int( vRC8600.size() ) ) {
	      if( vRC8600[khh] ) vRC8600[khh]->Fill( xh, yh );
//hh                             ----------------------------
	    }
	  }
	}
//      8692
	khh = 91;
	for( int kh=0; kh<4; kh++ ) {
	  xh = thep[kh] - theClu;
	  yh = 0.5;
	  if( khh >= 0  &&  khh < int( vRC8600.size() ) ) {
	    if( vRC8600[khh] ) vRC8600[khh]->Fill( xh, yh );
//hh                           ----------------------------
	  }
	}
//      8693
	khh = 92;
	for( int kh=0; kh<4; kh++ ) {
	  xh = thePhoSa - then[kh];
	  yh = 0.5;
	  if( khh >= 0  &&  khh < int( vRC8600.size() ) ) {
	    if( vRC8600[khh] ) vRC8600[khh]->Fill( xh, yh );
//hh                           ----------------------------
	  }
	}

	double tgx = vDcPhoDetWw.x()/vDcPhoDetWw.z();
	double tgy = vDcPhoDetWw.y()/vDcPhoDetWw.z();

	khh = 92 + kQua;
	xh = vDcPhoDetW.x()/vDcPhoDetW.z();
	yh = vDcPhoDetW.y()/vDcPhoDetW.z();
	if( khh >= 0  &&  khh < int( vRC8500.size() ) ) {
	  if( vRC8500[khh] ) vRC8500[khh]->Fill( xh, yh );
//hh                         ----------------------------
	}
	khh = 96 + kQua;
	xh = tgx;
	yh = tgy;
	if( khh >= 0  &&  khh < int( vRC8500.size() ) ) {
	  if( vRC8500[khh] ) vRC8500[khh]->Fill( xh, yh );
//hh                         ----------------------------
	}

	int ktx, kty;
	ktx = -1;
	kty = -1;
	for( int ky=0; ky<ncTxy; ky++ ) {
	  double yDirn = yDirLmn + ky * dDiry;
	  double yDirx = yDirn + dDiry;
	  if( tgy >= yDirn  &&  tgy < yDirx ) {
	    kty = ky;
	    ktx = -1;
	    for( int kx=0; kx<ncTxy; kx++ ) {
	      double xDirn = xDirLmn + kx * dDirx;
	      double xDirx = xDirn + dDirx;
	      if( tgx >= xDirn  &&  tgx < xDirx ) {
		ktx = kx;
		break;
	      }
	    }
	    break;
	  }
	}
	if( tgx < xDirLmn ) ktx = 0;
	if( tgx > xDirLmx ) ktx = ncTxy + 1;
	if( tgy < yDirLmn ) kty = 0;
	if( tgy > yDirLmx ) kty = ncTxy + 1;
	if( ktx <= 1 ) ktx = 1;
	if( ktx >= 5 ) ktx = 5;
	if( kty <= 1 ) kty = 1;
	if( kty >= 5 ) kty = 5;
	int kxy = (kty-1) * 5 + ktx - 1;

	int off = 0;
	off =  0;         //   all quadrants together
//@@------------
	khh = off + kxy;
	//xh = thePho;
	xh = thePhoSa;
	yh = phiPho;
	if( khh >= 0  &&  khh < int( vRC8500.size() ) ) {
	  if( vRC8500[khh] ) vRC8500[khh]->Fill( xh, yh );
//hh                         ----------------------------
	}

	//double dTheta = papho->getPMTOptCorr( clu, thePhoSa, phiPho );
	double dTheta = pho->getPMTOptCorr();
//      ------------------------------------
	//if( fabs( dTheta ) > 2. ) std::cout << dTheta << std::endl;

	off = 25;         //   all quadrants together
//@@------------
	khh = off + kxy;
	xh = thePhoSa + dTheta;
	//xh = thePhoSa - dTheta;
	yh = phiPho;
	if( khh >= 0  &&  khh < int( vRC8500.size() ) ) {
	  if( vRC8500[khh] ) vRC8500[khh]->Fill( xh, yh );
//hh                         ----------------------------
	}

	off = 50;
//@@------------
	khh =  off + kxy;
	xh = dTheta;
	yh = phiPho;
	//yh = phiPhow;
	if( khh >= 0  &&  khh < int( vRC8500.size() ) ) {
	  if( vRC8500[khh] ) vRC8500[khh]->Fill( xh, yh );
//hh                         ----------------------------
	}

	off = 75;
//@@------------
	khh =  off + PMTcha;
	xh = thePhoSa + dTheta;
	yh = phiPho;
	if( khh >= 0  &&  khh < int( vRC8500.size() ) ) {
	  if( vRC8500[khh] ) vRC8500[khh]->Fill( xh, yh );
//hh                         ----------------------------
	}


	if( iCat == dets->nCatPMT()[0] ) {      //  only TJ cathode
//        8651-75
	  khh = 50 + kxy;
	  xh = xPho - xClu;
	  yh = yPho - yClu;
	  if( khh >= 0  &&  khh < int( vRC8600.size() ) ) {
	    if( vRC8600[khh] ) vRC8600[khh]->Fill( xh, yh );
//hh                           ----------------------------
	  }
//        8676-91
	  khh = 75 + PMTcha;
	  xh = xPho - xClu;
	  yh = yPho - yClu;
	  if( khh >= 0  &&  khh < int( vRC8600.size() ) ) {
	    if( vRC8600[khh] ) vRC8600[khh]->Fill( xh, yh );
//hh                           ----------------------------
	  }
	}

	off = 0;
//@@------------
//      8701-16
	khh =  off + PMTcha;
	xh = thePhoSa - thePhoSaPi;
	yh = phiPho;
	if( khh >= 0  &&  khh < int( vRC8700.size() ) ) {
	  if( vRC8700[khh] ) vRC8700[khh]->Fill( xh, yh );
//hh                         ----------------------------
	}
	off = 16;
//@@------------
//      8717-32
	//std::cout << thePhoSa << "  " << thePhoSaPi << "  " << dTheta
	//          << std::endl;
	khh =  off + PMTcha;
	xh = thePhoSa+dTheta - thePhoSaPi;
	//xh = thePhoSa-dTheta - thePhoSaPi;
	yh = phiPho;
	if( khh >= 0  &&  khh < int( vRC8700.size() ) ) {
	  if( vRC8700[khh] ) vRC8700[khh]->Fill( xh, yh );
//hh                         ----------------------------
	}

      }

    }

    return;
  }



//===========================================================================
  std::string CsRCEventAnalysis::getVarStr( double var, int totLen,
//-----------------------------------------------------------------
					    int decLen0 ) {

//- Paolo - March  2010
//- From    std::string Cinema::getVarStr(...)

//- Fixed total string lenght, variable decimal digits if   decLen0 < 0
//- Fixed total string lenght, fixed decimal digits if      decLen0 >= 0

    std::string sVar = "";

    double pTot = pow( 10., totLen-1 );
    if( fabs( var ) > pTot ) {
      for( int ks=0; ks<totLen; ks++ ) sVar.append( "*" );
      return  sVar;
    }

    int varInt = int( var );
    int intLen = 0;
    int avarInt = abs( varInt );
    for( int k=0; k<10; k++ ) {
      if( avarInt <= 0 ) {
	intLen = k;
	break;
      }
      avarInt /= 10;
    }

    int decLen = decLen0;
    ///int decLen = totLen - intLen - 1;
    if( decLen0 < 0 ) decLen = totLen - intLen - 1;
    int pInt = int( pow( 10., decLen ) );
    if( var >= 0. ) var = int( var*pInt + 0.5 );
    if( var < 0. ) var = int( var*pInt - 0.5 );
    var /= pInt;
    bool pri = false;
    if( pri ) std::cout << var << "  " << varInt << "  " << intLen << "  ";

    varInt = int( var );
    if( pri ) std::cout << varInt << "  ";

    std::stringstream ssInt;
    ssInt << varInt;
    std::string sInt = ssInt.str();
    const char* ccInt = sInt.c_str();
    if( pri ) std::cout << ccInt << "  ";
    int nInt = strlen( ccInt );

    if( varInt == 0  &&  var < 0. ) totLen--;
    ///decLen = totLen - nInt - 1;
    if( decLen0 < 0 ) decLen = totLen - nInt - 1;
    if( pri ) std::cout << decLen << "  ";
    int pDec = int( pow( 10., decLen ) );
    if( pri ) std::cout << pDec << "  ";
    int varDec = int( fabs( var*pDec ) );
    if( pri ) std::cout << varDec << "  ";
    if( pri ) std::cout << varInt * pDec << "  ";
    varDec -= abs( varInt * pDec );
    if( pri ) std::cout << varDec << "  ";
    int nDecx = decLen;
    std::stringstream ssDec;
    ssDec << varDec;
    std::string sDec = ssDec.str();
    const char* ccDec = sDec.c_str();
    if( pri ) std::cout << ccDec << "  ";
    int nDec = strlen( ccDec );
    int nA = nDecx - nDec;
    if( pri ) std::cout << std::endl;

    sVar = "";
    //if( nInt > nIntx ) { sVar = "*"; nB--; sInt = sInt.substr( 2, nInt-2 );}
    if( varInt == 0 ) {
      if( var < 0. ) sVar.append( "-0." );
      else  sVar.append( "0." );
    } else  {
      sVar.append( sInt );
      sVar.append( "." );
    }
    if( varDec == 0 ) nA = decLen - 1;
    if( nA > 0 ) for( int k=0; k<nA; k++ ) sVar.append( "0" );
    sVar.append( sDec );

    int nB = totLen - sVar.size();
    if( decLen0 >= 0 ) {
      //int nB = totLen - sVar.size();
      std::string ss = "";
      if( nB > 0 ) for( int k=0; k<nB; k++ ) ss.append( " " );
      ss.append( sVar );
      sVar = ss;
    }
//- 100310
    if( nB < 0 ) {
      std::string ss = "*";
      int kB = -nB;
      //cout << nB << "  " << kB << endl;
      sVar = sVar.substr( kB, totLen );
      ss.append( sVar );
      sVar = ss;
    }

    return  sVar;
  }


//===========================================================================
  void CsRCEventAnalysis::singlePhotonCAT() {
//-------------------------------------------


//- Paolo  -  August  2010


    ///if( !CsRichOne::Instance()->UpRICHJob() ) return;

    if( !CsRCHistos::Ref().bookHis() ||
        CsRCHistos::Ref().levelBk() < 2 ) return;
//  --------------------------------------------

    CsRCRecConst *cons = CsRCRecConst::Instance();
    int kOffHi = cons->kOffHi();   
    CsRCDetectors *dets = CsRCDetectors::Instance();
    CsRCMirrors *mirr = CsRCMirrors::Instance();

    CsRCHistos& hist = CsRCHistos::Ref();
    float xh, yh;
    int khh;

    static std::vector<CsHist2D*> vRC8800;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      for( int kh=0; kh<100; kh++ ) vRC8800.push_back( NULL );
      if( CsRCHistos::Ref().bookHis() && CsRCHistos::Ref().levelBk() >= 2 ) {
	CsHistograms::SetCurrentPath("/RICH");

//----- 8800
	vRC8800.clear();
        for( int kh=0; kh<32; kh++ ) {
	  string hTitle = " ";
	  stringstream hN8800;
	  int kHist = kOffHi + 8800 + kh + 1;
	  hN8800 << kHist;
	  vRC8800.push_back( new CsHist2D( hN8800.str(), hTitle,
					   100, -10., 10., 90, 0., 360.) );
	}
        for( int kh=0; kh<16; kh++ ) {
	  string hTitle = " ";
	  stringstream hN8800;
	  int kHist = kOffHi + 8800 + 32 + kh + 1;
	  hN8800 << kHist;
	  vRC8800.push_back( new CsHist2D( hN8800.str(), hTitle,
					   100, -10., 10., 60, 0., 60.) );
	}
        for( int kh=0; kh<1; kh++ ) {
	  string hTitle = " ";
	  stringstream hN8800;
	  int kHist = kOffHi + 8800 + 48 + kh + 1;
	  hN8800 << kHist;
	  vRC8800.push_back( new CsHist2D( hN8800.str(), hTitle,
					   200, -1600., 1600.,
					   200, -1600., 1600.) );
	}
        for( int kh=0; kh<1; kh++ ) {
	  string hTitle = " ";
	  stringstream hN8800;
	  int kHist = kOffHi + 8800 + 49 + kh + 1;
	  hN8800 << kHist;
	  vRC8800.push_back( new CsHist2D( hN8800.str(), hTitle,
					   100, 0., 100., 16, 0., 16.) );
	}
        for( int kh=0; kh<16; kh++ ) {
	  string hTitle = " ";
	  stringstream hN8800;
	  int kHist = kOffHi + 8800 + 50 + kh + 1;
	  hN8800 << kHist;
	  vRC8800.push_back( new CsHist2D( hN8800.str(), hTitle,
					   100, -20., 20., 100, -20., 20.) );
	}
        for( int kh=0; kh<1; kh++ ) {
	  string hTitle = " ";
	  stringstream hN8800;
	  int kHist = kOffHi + 8800 + 66 + kh + 1;
	  hN8800 << kHist;
	  vRC8800.push_back( new CsHist2D( hN8800.str(), hTitle,
					   10, 0., 10., 3, 0., 3.) );
	}
	CsHistograms::SetCurrentPath("/");
      }
      std::cout << " RICHONE, CsRCEventAnalysis::singlePhotonCAT : " <<
	           "histo buffer " << vRC8800.size() << std::endl;

    }


    CsRCEventPartPhotons* evepapho = CsRCEventPartPhotons::Instance();
    list<CsRCPartPhotons*> lPaPhotons = evepapho->lPartPhotons();
    list<CsRCPartPhotons*>::iterator ipf;
    for( ipf=lPaPhotons.begin(); ipf!=lPaPhotons.end(); ipf++ ) {
      if( !(*ipf) ) continue;
      CsRCPartPhotons* papho = (*ipf);
      if( !papho->flag() ) continue;

      CsRCRing* ring = papho->pRing();
      bool RING = false;
      if( ring  &&  ring->flag() ) RING = true;
      if( !RING ) continue;
//    --------------------

//--- Select rings fully contained in each cathode :
//    ----------------------------------------------
      int kDetPart = papho->kDetPart();
      double xPade = papho->vPoPaDetW()[kDetPart].x();
      double yPade = papho->vPoPaDetW()[kDetPart].y();
      double RCing =  ring->theReco() * 3.3;
      int iCaPa = papho->iCaPa();
      if( iCaPa < 0 ) continue;      //***
//    ------------------------
      CsRCCathode* cat = dets->ptrToCat( iCaPa );
      double hCatx = cat->hCatx();
      double hCaty = cat->hCaty();
      bool bCRing1 = false;

      static int nPhoMinPMT = 20;
      static int nPhoMinAPV =  4;
//@@  --------------------------
      int nPhoMinRing = nPhoMinAPV;
      if( cat->isPMT() ) nPhoMinRing = nPhoMinPMT;

      Hep3Vector vOffCatW;
      vOffCatW = dets->vOffCatW( iCaPa );
      double xLmn = vOffCatW.x() - hCatx;
      double xLmx = vOffCatW.x() + hCatx;
      double yLmn = vOffCatW.y() - hCaty;
      double yLmx = vOffCatW.y() + hCaty;
      if( (xPade - xLmn) > RCing  &&  (xLmx - xPade) > RCing  &&
	  (yLmx - yPade) > RCing  &&  (yPade - yLmn) > RCing )
	bCRing1 = true;
      //std::cout << xLmn << "  " << xLmn+2.*hCatx << " - "
      //	  << yLmx << "  " << yLmx-2.*hCaty << std::endl;
      if( !bCRing1 ) continue;      //***
//    -----------------------

//--- Cut on low particle momentum :
      static double momMin = 3.;
//@@  -------------------------
      CsRCParticle* part = papho->pPart();
      double momPart = part->mom();
      if( momPart < momMin ) continue;
//    -------------------------------

//--- Select ID pions (analysis way)
      bool bIDpion = false;
      double likeDV = cons->likeDefVa();
      double* probaLK;
      probaLK = papho->probaLKAll();
      double likeBkg = papho->probaLKBgAll();
      double likeElec = probaLK[ 2];
      double likeMuon = probaLK[ 5];
      double likePion = probaLK[ 8];
      double likeKaon = probaLK[11];
      double likeProton = probaLK[14];
      if( likeElec == likeBkg ) continue;   // NO photons from signal
      if( likePion > likeDV  &&  likePion > likeBkg  &&
	  likePion > likeKaon  &&  likePion > likeProton ) bIDpion = true;
      if( momPart < 8.  &&  likePion < likeElec ) bIDpion = false;
      if( !bIDpion ) continue;      //***
//    -----------------------


//--- PART-PHOT photons
//    -----------------
      list<CsRCPhoton*> lPhotons = papho->lPhotons();
      int nPhotons = lPhotons.size();
//--- RING photons
//    ------------
      list<CsRCPhoton*> lPhotonsR = ring->lPhotons();
      int nPhotonsR = lPhotonsR.size();

//--- Select RINGs with many photons
      if( nPhotonsR < nPhoMinRing ) continue;      //***
//    --------------------------------------

//    8850
      khh = 49;
      xh = nPhotonsR;
      yh = iCaPa;
      if( khh >= 0  &&  khh < int( vRC8800.size() ) ) {
	if( vRC8800[khh] ) vRC8800[khh]->Fill( xh, yh );
//hh                       ----------------------------
      }

//--- assume pion mass
//    ----------------
      double massIpo = CsRCRecConst::Instance()->massPartv()[8];
      double betaIpo = momPart / sqrt( momPart*momPart + massIpo*massIpo );
      double cosThe = 0.;
      cosThe = 1./ ( betaIpo * cons->CFRefIndUV() );
      if( cosThe > 1.) cosThe = 1.;
      double thetaUV = acos( cosThe ) * 1000.;        //   mrad
      cosThe = 1./ ( betaIpo * cons->CFRefIndVS() );
      if( cosThe > 1.) cosThe = 1.;
      double thetaVS = acos( cosThe ) * 1000.;        //   mrad
      //std::cout << thetaUV << "  " << thetaVS << std::endl;
      double thePion = thetaUV;
      double thePionUV = thetaUV;
      if( cat->isPMT() ) thePion = thetaVS;
      double theRing = ring->the();        //   mrad

//--- PART-PHOT photons
      list<CsRCPhoton*>::iterator ih;
      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	if( !(*ih) ) continue;
	CsRCPhoton* pho = (*ih);

	double thePho = pho->the0();         //   mrad
	double phiPho = pho->phiA();         //   deg
	CsRCCluster* clu = pho->ptToClu();
	if( !clu ) continue;
	int iCat = clu->ic();
        double xc = clu->xc();
	double yc = clu->yc();
	int kClu = clu->kClu();
	list<CsRCPad*> lPads = clu->lPads();

//      8801-8816
	khh = iCat;
	xh = thePho - thePion;
	yh = phiPho;
	if( khh >= 0  &&  khh < int( vRC8800.size() ) ) {
	  if( vRC8800[khh] ) vRC8800[khh]->Fill( xh, yh );
//hh                         ----------------------------
	}

	double thePhoR = pho->the();         //   mrad
//      8817-8832
	khh = iCat + 16;
	xh = thePhoR - theRing;
	yh = phiPho;
	if( khh >= 0  &&  khh < int( vRC8800.size() ) ) {
	  if( vRC8800[khh] ) vRC8800[khh]->Fill( xh, yh );
//hh                         ----------------------------
	}
//      8849
	khh = 48;
	xh = xc;
	yh = yc;
	if( khh >= 0  &&  khh < int( vRC8800.size() ) ) {
	  if( vRC8800[khh] ) vRC8800[khh]->Fill( xh, yh );
//hh                         ----------------------------
	}
//      8867
	khh = 66;
	xh = lPads.size();
	yh = 0.5;
	if( khh >= 0  &&  khh < int( vRC8800.size() ) ) {
	  if( vRC8800[khh] ) vRC8800[khh]->Fill( xh, yh );
//hh                         ----------------------------
	}
	if( clu->isAPV() ) yh = 1.5;
	if( clu->isPMT() ) yh = 2.5;
	if( khh >= 0  &&  khh < int( vRC8800.size() ) ) {
	  if( vRC8800[khh] ) vRC8800[khh]->Fill( xh, yh );
//hh                         ----------------------------
	}
      }

//    8833-8848
      khh = iCaPa + 32;
      xh = theRing - thePionUV;
      yh = momPart;
      if( khh >= 0  &&  khh < int( vRC8800.size() ) ) {
	if( vRC8800[khh] ) vRC8800[khh]->Fill( xh, yh );
//hh                       ----------------------------
      }

//--- FITs :
//    ------
      //bool exeFitE = true;
      bool exeFitE = false;
      bool exeFitC = true;
      //bool exeFitC = false;
      static vector<CsHist2D*> vRC8890;
      static int nFitToE = 0;
      static int nFitFailE = 0;
      static int nFitToC = 0;
      static int nFitFailC = 0;
      if( abs( ring->kRingLoc() ) != 1 ) exeFitE = false;      // (?)
//@@----------------------------------------------------
      if( exeFitE || exeFitC ) {

	static bool fitDef = true;
	if( fitDef ) {
	  fitDef = false;
	  nFitToE = 0;
	  nFitFailE = 0;
	  nFitToC = 0;
	  nFitFailC = 0;
	  for( int j=0; j<10; j++) vRC8890.push_back( NULL );
	  if( CsRCHistos::Ref().bookHis() &&
	      CsRCHistos::Ref().levelBk() >= 2 ) {
	    CsHistograms::SetCurrentPath("/RICH");
	    if( exeFitE ) {
	      string hTitle;
	      vRC8890.clear();
	      stringstream hN8891;
	      hN8891 << kOffHi + 8891;
	      hTitle = "c0fit - c0in";
	      vRC8890.push_back( new CsHist2D( hN8891.str(), hTitle,
					       100, -10., 10.,
					       100, -10., 10. ) );
	      stringstream hN8892;
	      hN8892 << kOffHi + 8892;
	      hTitle = "AAfit - AAin vs ddet";
	      vRC8890.push_back( new CsHist2D( hN8892.str(), hTitle,
					       100, -10., 10.,
					       100, 0., 1000. ) );
	      stringstream hN8893;
	      hN8893 << kOffHi + 8893;
	      hTitle = "epsi vs ddet";
	      vRC8890.push_back( new CsHist2D( hN8893.str(), hTitle,
					       80, 0.8, 1.2,
					       100, 0., 1000. ) );
	      stringstream hN8894;
	      hN8894 << kOffHi + 8894;
	      hTitle = "chisq / nu";
	      vRC8890.push_back( new CsHist2D( hN8894.str(), hTitle,
					       100, 0., 10.,
					       100, 0., 100. ) );
	      stringstream hN8895;
	      hN8895 << kOffHi + 8895;
	      hTitle = "n iter";
	      vRC8890.push_back( new CsHist2D( hN8895.str(), hTitle,
					       20, 0., 20., 100, 0., 100. ) );
	      stringstream hN8896;
	      hN8896 << kOffHi + 8896;
	      hTitle = "fit fails";
	      vRC8890.push_back( new CsHist2D( hN8896.str(), hTitle,
					       10, 0., 10., 20, 0., 20. ) );
	      stringstream hN8897;
	      hN8897 << kOffHi + 8897;
	      hTitle = "pulls";
	      vRC8890.push_back( new CsHist2D( hN8897.str(), hTitle,
					       100, -5., 5.,
					       100, 0., 100. ) );
	      stringstream hN8898;
	      hN8898 << kOffHi + 8898;
	      hTitle = "phi vs phiPa";
	      vRC8890.push_back( new CsHist2D( hN8898.str(), hTitle,
					       90, -45., 45.,
					       90, 0., 360. ) );
	    }
	    if( exeFitC ) {
	      string hTitle;
	      vRC8890.clear();
	      stringstream hN8891;
	      hN8891 << kOffHi + 8891;
	      hTitle = "c0fit - c0in";
	      vRC8890.push_back( new CsHist2D( hN8891.str(), hTitle,
					       100, -10., 10.,
					       100, -10., 10. ) );
	      stringstream hN8892;
	      hN8892 << kOffHi + 8892;
	      hTitle = "AAfit - AAin vs ddet";
	      vRC8890.push_back( new CsHist2D( hN8892.str(), hTitle,
					       100, -10., 10.,
					       100, 0., 1000. ) );
	      stringstream hN8893;
	      hN8893 << kOffHi + 8893;
	      hTitle = "chisq / nu";
	      vRC8890.push_back( new CsHist2D( hN8893.str(), hTitle,
					       100, 0., 10.,
					       100, 0., 100. ) );
	      stringstream hN8894;
	      hN8894 << kOffHi + 8894;
	      hTitle = "n iter";
	      vRC8890.push_back( new CsHist2D( hN8894.str(), hTitle,
					       20, 0., 20., 100, 0., 100. ) );
	      stringstream hN8895;
	      hN8895 << kOffHi + 8895;
	      hTitle = "pulls";
	      vRC8890.push_back( new CsHist2D( hN8895.str(), hTitle,
					       100, -5., 5.,
					       100, 0., 100. ) );
	    }
	    CsHistograms::SetCurrentPath("/");
	  }
	}

	double xClu[nPhotonsR];
	double yClu[nPhotonsR];
	double errClu[nPhotonsR];
	double theR[nPhotonsR];
	double phiR[nPhotonsR];
	double errR[nPhotonsR];
	list<CsRCPhoton*>::iterator ih;
	double sigma = 5.;
//@@    -----------------
	double radDeg = cons->RadDeg();
	double sigPhoRec = papho->sigmaPhoRec( momPart );
//                         -----------------------------
	int nSig = 3;
//@@    ------------
	double dSigC = nSig * sigPhoRec;
	double theLoL = theRing - dSigC;
	double theUpL = theRing + dSigC;
	int kPhot = 0;
	for( ih=lPhotonsR.begin(); ih!=lPhotonsR.end(); ih++ ) {
	  if( !(*ih) ) continue;
	  CsRCPhoton* pho = (*ih);
	  double theW = pho->the();          //   mrad
	  double phiW = pho->phiA();         //   deg
	  if( theW < theLoL  ||  theW > theUpL ) continue;
	  CsRCCluster* clu = pho->ptToClu();
	  if( !clu ) continue;
	  xClu[kPhot] = clu->xc() - xPade;
	  yClu[kPhot] = clu->yc() - yPade;
	  errClu[kPhot] = sigma*sigma;
	  double sigW = pho->sigmaPhoPid( papho );
//@@                    -------------------------
	  theR[kPhot] = theW * cos( phiW/radDeg );
	  phiR[kPhot] = theW * sin( phiW/radDeg );
	  errR[kPhot] = sigW*sigW;
	  kPhot++;
	}
	int nClus = kPhot;
	int nPhot = kPhot;

	if( exeFitE ) {
//----- NOT WORKING !!! ??? (... HepMatrix(... ,0 )??? )
	  std::cout << nPhotonsR << "  " << nClus << std::endl;
	  for( int k=0; k<nClus; k++ ) {
	    std::cout << xClu[k] << " " << yClu[k] << std::endl;
	  }
	  int nParamE = 7;
	  double paramE[nParamE];
	  paramE[0] = 0.;                //   x-center (mm)
	  paramE[1] = 0.;                //   y-center (mm)
	  paramE[2] = theRing *3.35;     //   main axis (mm)
	  if( cat->isPMT() ) paramE[2] = papho->thetaUVtoVS( theRing ) * 3.35;
	  paramE[3] = 1.0;               //   eccentricity
	  //paramE[4] = 0.               //   main axis angle (rad)
	  if( xPade != 0. ) {
	    paramE[4] = atan( yPade/xPade );
	  } else {
	    paramE[4] = 1.5708;
	  }
	  paramE[4] = 0.;
	  paramE[5] = xPade;
	  paramE[6] = yPade;
	  int iPaFitE[nParamE];
	  iPaFitE[0] = 0;
	  iPaFitE[1] = 0;
	  iPaFitE[2] = 0;
	  iPaFitE[3] = 1;//*
	  iPaFitE[4] = 1;//*
	  iPaFitE[5] = 1;
	  iPaFitE[6] = 1;
	  nFitToE++;
	  //nClus = 5;//*
	  double errClf = errClu[0];//*
	  //*CsRCEllipseFit oEllipse( nClus, xClu, yClu, errClu, errClu,
	  CsRCEllipseFitTest oEllipse( nClus, xClu, yClu, errClu, errClu,//*
//        ---------------------------------------------------------------
				       nParamE, paramE, iPaFitE );
	  
	  if( oEllipse.doChiFit() ) {
	    std::cout << oEllipse.flag() << std::endl;
	    oEllipse.print();
	    oEllipse.doHist( vRC8890 );
//          --------------------------
	    std::cout << "done1" << std::endl;
	  } else {
	    nFitFailE++;
	    std::cout << "Fit failure " << nFitFailE << "  " << nFitToE
	              << std::endl;
	  }
	}

	if( exeFitC ) {
	  //std::cout << nPhotonsR << "  " << nClus << std::endl;
	  //for( int k=0; k<nClus; k++ ) {
	    //std::cout << xClu[k] << " " << yClu[k] << std::endl;
	  //}
	  int nParamC = 3;
	  double paramC[nParamC];
	  //paramC[0] = theRing;           //   Theta Ring radius (mrad)
	  paramC[0] = theRing * 3.35;    //   Ring radius (mm)
	  if( cat->isPMT() ) paramC[0] = papho->thetaUVtoVS( theRing ) * 3.35;
	  paramC[1] = 0.;                //   Ring x-center (mrad)
	  paramC[2] = 0.;                //   Ring y-center (mrad)
	  int iPaFitC[nParamC];
	  iPaFitC[0] = 0;
	  iPaFitC[1] = 0;
	  iPaFitC[2] = 0;
	  nFitToC++;
	  
 	  //CsRCCircleFit oCircle( nPhot, theR, phiR, errR,
//        -----------------------------------------------
	  //		  nParamC, paramC, iPaFitC );
 	  CsRCCircleFit oCircle( nClus, xClu, yClu, errClu,
//        -------------------------------------------------
				 nParamC, paramC, iPaFitC );
	  
	  if( oCircle.doChiFit() ) {
	    //std::cout << oCircle.flag() << std::endl;
	    //oCircle.print();
	    oCircle.doHist( vRC8890 );
//          --------------------------
	    //std::cout << "done1" << std::endl;

//          8851-8866
	    khh = iCaPa + 50;
	    xh = oCircle.para()[1];
	    yh = oCircle.para()[2];
	    if( khh >= 0  &&  khh < int( vRC8800.size() ) ) {
	      if( vRC8800[khh] ) vRC8800[khh]->Fill( xh, yh );
//hh                             ----------------------------
	    }
	  } else {
	    nFitFailC++;
	    //std::cout << "Fit failure " << nFitFailC << "  " << nFitToC
	    //          << std::endl;
	  }
	}
      }
      //std::cout << "done2" << std::endl;


//  16/12/2010
//  ----------
      bool exeFitD = true;
      //bool exeFitD = false;
      static vector<CsHist2D*> vRC8950;
      static int nFitToD = 0;
      static int nFitFailD = 0;
      Hep3Vector vPoPaFitW;
      for( int j=0; j<10; j++) vRC8950.push_back( NULL );
      if( exeFitD ) {

	static bool fitDef = true;
	if( fitDef ) {
	  fitDef = false;
	  vRC8950.clear();
	  string hTitle;
	  stringstream hN8951;
	  hN8951 << kOffHi + 8951;
	  hTitle = "c0fit - c0in";
	  vRC8950.push_back( new CsHist2D( hN8951.str(), hTitle,
					   100, -10., 10.,
					   100, -10., 10. ) );
	  stringstream hN8952;
	  hN8952 << kOffHi + 8952;
	  hTitle = "AAfit - AAin vs ddet";
	  vRC8950.push_back( new CsHist2D( hN8952.str(), hTitle,
					   100, -10., 10.,
					   100, 0., 1000. ) );
	  stringstream hN8953;
	  hN8953 << kOffHi + 8953;
	  hTitle = "chisq / nu";
	  vRC8950.push_back( new CsHist2D( hN8953.str(), hTitle,
					   100, 0., 10.,
					   100, 0., 100. ) );
	  stringstream hN8954;
	  hN8954 << kOffHi + 8954;
	  hTitle = "n iter";
	  vRC8950.push_back( new CsHist2D( hN8954.str(), hTitle,
					   20, 0., 20., 100, 0., 100. ) );
	  stringstream hN8955;
	  hN8955 << kOffHi + 8955;
	  hTitle = "pulls";
	  vRC8950.push_back( new CsHist2D( hN8955.str(), hTitle,
					   100, -5., 5.,
					   100, 0., 100. ) );
	  stringstream hN8956;
	  hN8956 << kOffHi + 8956;
	  hTitle = "dTgy vs dTgx";
	  vRC8950.push_back( new CsHist2D( hN8956.str(), hTitle,
					   100, -10., 10.,
					   100, -10., 10. ) );
	  stringstream hN8957;
	  hN8957 << kOffHi + 8957;
	  hTitle = "dTgy vs dTgx";
	  vRC8950.push_back( new CsHist2D( hN8957.str(), hTitle,
					   100, -10., 10.,
					   100, -10., 10. ) );
	  stringstream hN8958;
	  hN8958 << kOffHi + 8958;
	  hTitle = "dTgy vs dTgx";
	  vRC8950.push_back( new CsHist2D( hN8958.str(), hTitle,
					   100, -10., 10.,
					   100, -10., 10. ) );
	  stringstream hN8959;
	  hN8959 << kOffHi + 8959;
	  hTitle = "dTgy vs dTgx";
	  vRC8950.push_back( new CsHist2D( hN8959.str(), hTitle,
					   100, -10., 10.,
					   100, -10., 10. ) );
	  CsHistograms::SetCurrentPath("/");
	}

	nFitToC++;
	CsRCCircleFit oCircle;
	if( ring->getDetRFit( vRC8950, oCircle ) );
//          ------------------------------------
	CsRCCircleFit* pCircle = &oCircle;
	if( pCircle->flag()  &&  pCircle->nIter() < 10 ) {

	//CsRCEllipseFit oEllipse;                       // NOT working!
	//if( ring->getDetEFit( vRC8950, oEllipse ) );
//          -------------------------------------
	//CsRCEllipseFit* pEllipse = &oEllipse;
	//if( pEllipse->flag()  &&  pEllipse->nIter() < 10 ) {

	  double xCfit = pCircle->para()[1] + xPade;
	  double yCfit = pCircle->para()[2] + yPade;
	  double zCfit = dets->vOffCatW( iCaPa ).z();
	  Hep3Vector vPoPaDet( xCfit, yCfit, zCfit );
	  vPoPaFitW = vPoPaDet;
//------- particle 'back impact' on mirror (MWR) :
          Hep3Vector vPoC( 0., 0., 0. );
	  int kDetPart = papho->kDetPart();
	  double RR = mirr->RRv( kDetPart );
//------- assume vPoPaMirW as valid 'back impact' of particle on mirror
	  //Hep3Vector vPoPaMirW = papho->vPoPaMirW()[kDetPart];
	  //Hep3Vector vDcPaMirBk = vPoPaMirW - vPoPaDet;
	  //vDcPaMirBk.setX( vDcPaMirBk.x() * 1.1 );
	  //Hep3Vector vDcPaMirBk( 0., 0., 1. );
//------- assume vDcPaReflW as valid 'back direction' of particle to mirror
	  Hep3Vector vDcPaMirBk = papho->vDcPaReflW()[kDetPart];
	  vDcPaMirBk = vDcPaMirBk.unit();
	  //std::cout << setprecision( 6 );
	  //std::cout << vDcPaMirBk << "  " << std::endl;
          Hep3Vector vPoPaMirBk;
	  vPoPaMirBk = mirr->vImpMir( vPoPaDet, vDcPaMirBk, vPoC, RR );
//                     -----------------------------------------------
//------- normal to mirror at particle impact :
	  Hep3Vector vDcNoPaMirW = (1./RR) * (vPoPaMirBk - vPoC);
//------- particle 'reflected' direction :
	  double cosPaMir = vDcNoPaMirW * vDcPaMirBk;
	  Hep3Vector vDcPaReflW = 2.*cosPaMir * vDcNoPaMirW - vDcPaMirBk;
	  Hep3Vector vDcPartW = papho->vDcPartW()[kDetPart];
	  //std::cout << setprecision( 6 );
	  //std::cout << vDcPartW << "  " << vDcPaReflW << "  "
	  //<< vDcPartW-vDcPaReflW << std::endl;
	  double dTgx, dTgy;
	  dTgx = vDcPaReflW.x()/vDcPaReflW.z() - vDcPartW.x()/vDcPartW.z();
	  dTgy = vDcPaReflW.y()/vDcPaReflW.z() - vDcPartW.y()/vDcPartW.z();
//------- 8956
	  if( iCaPa ==  3 ) khh = 5;
	  if( iCaPa ==  5 ) khh = 6;
	  if( iCaPa == 10 ) khh = 7;
	  if( iCaPa == 12 ) khh = 8;
	  xh = dTgx * 1000.;
	  yh = dTgy * 1000.;
	  if( khh >= 0  &&  khh < int( vRC8950.size() ) ) {
	    if( vRC8950[khh] ) vRC8950[khh]->Fill( xh, yh );
//hh                           ----------------------------
	  }

	} else {
	  nFitFailD++;
	}
      }

      double dAlpha = 0.;
      double dBeta  = 0.;
//@@------------------------
      double xPaWin = part->vPosIn().x();
      double yPaWin = part->vPosIn().y();
      double zPaWin = part->vPosIn().z();
      double xPaWinP = xPaWin - zPaWin * dAlpha;
      double yPaWinP = yPaWin - zPaWin * dBeta;
      Hep3Vector vPoWin0P( xPaWinP, yPaWinP, zPaWin );
      double tgA = part->vDirIn().x() / part->vDirIn().z();
      double tgB = part->vDirIn().y() / part->vDirIn().z();
      double tgAP = tgA + dAlpha;
      double tgBP = tgB + dBeta;
      Hep3Vector vDcWin0P( tgAP, tgBP, 1. );
      vDcWin0P = vDcWin0P.unit();

      Hep3Vector vPoC0 = part->pMirPart()->vC0();
      double RR = mirr->RRv( kDetPart );
      //std::cout << RR << "  " << vPoC0 << std::endl;
      //vPoC0 += mirr->vCorPoPa0();    // wrong at this point
      //RR += mirr->corPaRR();         // wrong at this point
      Hep3Vector vPoC0c = vPoC0 + part->vCorPoPa0();   // needs patch 101220
      double RRc = RR + part->corPaRR();               // needs patch 101220
      //std::cout << mirr->corPaRR() << "  " << mirr->vCorPoPa0() <<std::endl;
      Hep3Vector vPoPaMir0P = mirr->vImpMir( vPoWin0P, vDcWin0P, vPoC0c, RRc);
//                            -----------------------------------------------
      Hep3Vector vDcNoPaMir0P = (1./RRc) * (vPoPaMir0P - vPoC0c);
      double cosPaMirP = vDcNoPaMir0P * vDcWin0P;
      Hep3Vector vDcPaRefl0P = 2.*cosPaMirP * vDcNoPaMir0P - vDcWin0P;
      Hep3Vector vDcDet0 = dets->vDcDetv( kDetPart );
      Hep3Vector vDet0 = dets->vDet0v( kDetPart );
      std::cout << setprecision(6);
      //std::cout << vDet0 << "  " << vDcDet0 << std::endl;
      double norm = ( (vDet0 - vPoPaMir0P) * vDcDet0 ) /
	( vDcPaRefl0P * vDcDet0 );
      Hep3Vector vPoPaDet0P = vPoPaMir0P + norm * vDcPaRefl0P;
      double ddd = vDcDet0 * (vPoPaDet0P - vDet0);      // test
      HepVector vrot( 3, 0 );
      //for( int j=0; j<3; j++ ) vrot[j] = (vPoPaDet0P - vPoC0c)[j];
      for( int j=0; j<3; j++ ) vrot[j] = (vPoPaDet0P - vPoC0)[j];   // ???
      list<CsRCPhotonDet*> lPhoDet = dets->lPhoDet();
      list<CsRCPhotonDet*>::iterator id = lPhoDet.begin();
      if( kDetPart == 1 ) id++;
      HepVector ww = (*id)->rotMatrix() * vrot;
      Hep3Vector vPoPaDetWP( ww[0], ww[1], ww[2] );
      //double the = atan( vDcDet0.y() / vDcDet0.z() );
      //std::cout <<  vDcDet0.y() << std::endl;
      //Hep3Vector vRotc = vPoPaDet0P - vPoC0c;   // test
      //Hep3Vector vRot = vPoPaDet0P - vPoC0;   // ???
      //vRot = vRot.unit();         // test
      //vRotc = vRotc.unit();       // test
      //double ccc = vRot * vRotc;  // test
      //std::cout <<  ccc << std::endl;
      //std::cout <<  vRot << std::endl;
      //Hep3Vector vPoPaDetWP;
      //vPoPaDetWP.setX(   vRot.x() );
      //vPoPaDetWP.setY( - vRot.z() * sin(the) + vRot.y() * cos(the) );
      //vPoPaDetWP.setZ(   vRot.z() * cos(the) + vRot.y() * sin(the) );
      Hep3Vector vPade = papho->pPart()->vPade()[kDetPart];
      //std::cout << ddd << "  " << vPoPaDetWP.z()-vPoPaFitW.z() << std::endl;
      //std::cout << setprecision(1);
      //std::cout << vPoPaDetWP << "  " << vPade << "  " << vPoPaFitW
      //       	<< std::endl;
    }

    return;
  }




//===========================================================================
  void CsRCEventAnalysis::checkLikeDisTheta() {
//---------------------------------------------


//- Paolo  -  September  2010


    if( !CsRichOne::Instance()->UpRICHJob() ) return;

    if( !CsRCHistos::Ref().bookHis() ||
        CsRCHistos::Ref().levelBk() < 2 ) return;
//  --------------------------------------------

    CsRCRecConst *cons = CsRCRecConst::Instance();
    int kOffHi = cons->kOffHi();   
    CsRCDetectors *dets = CsRCDetectors::Instance();
    CsRCMirrors *mirr = CsRCMirrors::Instance();

    CsRCHistos& hist = CsRCHistos::Ref();
    float xh, yh;
    int khh;

    static std::vector<CsHist2D*> vRC8900;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      for( int kh=0; kh<20; kh++ ) vRC8900.push_back( NULL );
      if( CsRCHistos::Ref().bookHis() && CsRCHistos::Ref().levelBk() >= 2 ) {
	CsHistograms::SetCurrentPath("/RICH");

//----- 8900 (01-20)
	vRC8900.clear();
        for( int kh=0; kh<20; kh++ ) {
	  string hTitle = " ";
	  stringstream hN8900;
	  int kHist = kOffHi + 8900 + kh + 1;
	  hN8900 << kHist;
	  vRC8900.push_back( new CsHist2D( hN8900.str(), hTitle,
					   200, 0., 5., 60, 0., 60.) );
	}
	CsHistograms::SetCurrentPath("/");
      }
      std::cout << " RICHONE, CsRCEventAnalysis::checkLikeDisTheta : " <<
	           "histo buffer " << vRC8900.size() << std::endl;

    }


    CsRCEventPartPhotons* evepapho = CsRCEventPartPhotons::Instance();
    list<CsRCPartPhotons*> lPaPhotons = evepapho->lPartPhotons();
    list<CsRCPartPhotons*>::iterator ipf;
    for( ipf=lPaPhotons.begin(); ipf!=lPaPhotons.end(); ipf++ ) {
      if( !(*ipf) ) continue;
      CsRCPartPhotons* papho = (*ipf);
      if( !papho->flag() ) continue;

      CsRCRing* ring = papho->pRing();
      bool RING = false;
      if( ring  &&  ring->flag() ) RING = true;
      if( !RING ) continue;
//    --------------------

//--- Select rings fully contained in each cathode :
//    ----------------------------------------------
      int kDetPart = papho->kDetPart();
      double xPade = papho->vPoPaDetW()[kDetPart].x();
      double yPade = papho->vPoPaDetW()[kDetPart].y();
      double RCing =  ring->theReco() * 3.3;
      int iCaPa = papho->iCaPa();
      if( iCaPa < 0 ) continue;      //***
//    ------------------------
      CsRCCathode* cat = dets->ptrToCat( iCaPa );
      double hCatx = cat->hCatx();
      double hCaty = cat->hCaty();
      bool bCRing1 = false;

      static int nPhoMinPMT = 20;
      static int nPhoMinAPV =  4;
//@@  --------------------------
      int nPhoMinRing = nPhoMinAPV;
      if( cat->isPMT() ) nPhoMinRing = nPhoMinPMT;

      Hep3Vector vOffCatW;
      vOffCatW = dets->vOffCatW( iCaPa );
      double xLmn = vOffCatW.x() - hCatx;
      double xLmx = vOffCatW.x() + hCatx;
      double yLmn = vOffCatW.y() - hCaty;
      double yLmx = vOffCatW.y() + hCaty;
      if( (xPade - xLmn) > RCing  &&  (xLmx - xPade) > RCing  &&
	  (yLmx - yPade) > RCing  &&  (yPade - yLmn) > RCing )
	bCRing1 = true;
      //std::cout << xLmn << "  " << xLmn+2.*hCatx << " - "
      //	  << yLmx << "  " << yLmx-2.*hCaty << std::endl;
      if( !bCRing1 ) continue;      //***
//    -----------------------

//--- Cut on low particle momentum :
      static double momMin = 3.;
//@@  -------------------------
      CsRCParticle* part = papho->pPart();
      double momPart = part->mom();
      if( momPart < momMin ) continue;      //***
//    -------------------------------

//--- Select RINGs with many photons
      list<CsRCPhoton*> lPhotonsR = ring->lPhotons();
      int nPhotonsR = lPhotonsR.size();
      if( nPhotonsR < nPhoMinRing ) continue;      //***
//    --------------------------------------


//--- Select ID particles
      double likeDV = cons->likeDefVa();
      double* probaLK;
      probaLK = papho->probaLKAll();
      double likeBkg = papho->probaLKBgAll();
      double likeElec = probaLK[ 2];
      double likeMuon = probaLK[ 5];
      double likePion = probaLK[ 8];
      double likeKaon = probaLK[11];
      double likeProton = probaLK[14];
      double likeVar[10];
      likeVar[0] = likeBkg;
      likeVar[1] = likeElec;
      likeVar[2] = likeMuon;
      likeVar[3] = likePion;
      likeVar[4] = likeKaon;
      likeVar[5] = likeProton;
      //std::cout << setprecision( 5 );
      //for( int k=0; k<6; k++ ) std::cout << likeVar[k] << "  ";
      //std::cout << std::endl;
      bool bID = true;
      if( likeBkg == 0. ) bID = false;
      if( likeElec == likeBkg ) bID = false;   // NO photons from signal

//--- Select ID pions (analysis way)
      bool bIDpion = false;
      if( likePion > likeDV  &&  likePion > likeBkg  &&
	  likePion > likeKaon  &&  likePion > likeProton ) bIDpion = true;
      if( momPart < 8.  &&  likePion < likeElec ) bIDpion = false;

//--- Select ID kaons (analysis way)
      bool bIDkaon = false;
      if( likeKaon > likeDV  &&  likeKaon > likeBkg  &&
	  likeKaon > likePion  &&  likeKaon > likeProton ) bIDkaon = true;

      if( bID ) {
//    8901-06
      int kh0 = 0;
      for( int kl=0; kl<6; kl++ ) {
	khh = kh0 + kl;
	if( likeVar[kl] > 0. ) {
	  xh = likeVar[kl];
	  yh = momPart;
	  if( khh >= 0  &&  khh < int( vRC8900.size() ) ) {
	    if( vRC8900[khh] ) vRC8900[khh]->Fill( xh, yh );
//hh                           ----------------------------
	  }
	}
      }

//    8907-12
      kh0 = 6;
      if( likeBkg > 0. ) {
	for( int kl=0; kl<6; kl++ ) {
	  int kMax = -1;
	  double likeMax = 0.;
	  for( int kl=0; kl<6; kl++ ) {
	    if( likeVar[kl] > likeMax ) {
	      likeMax = likeVar[kl];
	      kMax = kl;
	    }
	  }
	  if( kMax >= 0  &&  likeMax > 0. ) {
	    khh = kh0 + kMax;
	    xh = likeMax;
	    yh = momPart;
	    if( khh >= 0  &&  khh < int( vRC8900.size() ) ) {
	      if( vRC8900[khh] ) vRC8900[khh]->Fill( xh, yh );
//hh                             ----------------------------
	    }
	  }
	}
      }

      double fact = 0.2;
//    -----------------
//    8913-15   Pion
      kh0 = 12;
      if( bIDpion  &&  likeBkg > 0. ) {
	khh = kh0;
	xh = likePion;
	yh = momPart;
	if( khh >= 0  &&  khh < int( vRC8900.size() ) ) {
	  if( vRC8900[khh] ) vRC8900[khh]->Fill( xh, yh );
//hh                         ----------------------------
	}
	khh = kh0 + 1;
	xh = fact * likePion / likeBkg;
	yh = momPart;
	if( khh >= 0  &&  khh < int( vRC8900.size() ) ) {
	  if( vRC8900[khh] ) vRC8900[khh]->Fill( xh, yh );
//hh                         ----------------------------
	}
	int kMax2 = -1;
	double likeMax2 = 0.;
	for( int kl=0; kl<6; kl++ ) {
	  if( kl == 1 ) continue;
	  if( kl == 2 ) continue;
	  if( kl == 3 ) continue;
	  if( likeVar[kl] > likeMax2 ) {
	    likeMax2 = likeVar[kl];
	    kMax2 = kl;
	  }
	}
	if( kMax2 >= 0 ) {
	  khh = kh0 + 2;
	  xh = fact * likePion / likeMax2;
	  yh = momPart;
	  if( khh >= 0  &&  khh < int( vRC8900.size() ) ) {
	    if( vRC8900[khh] ) vRC8900[khh]->Fill( xh, yh );
//hh                           ----------------------------
	  }
	}
      }

//    8916-18   Kaon
      kh0 = 15;
      if( bIDkaon  &&  likeBkg > 0. ) {
	khh = kh0;
	xh = likeKaon;
	yh = momPart;
	if( khh >= 0  &&  khh < int( vRC8900.size() ) ) {
	  if( vRC8900[khh] ) vRC8900[khh]->Fill( xh, yh );
//hh                         ----------------------------
	}
	khh = kh0 + 1;
	xh = fact * likeKaon / likeBkg;
	yh = momPart;
	if( khh >= 0  &&  khh < int( vRC8900.size() ) ) {
	  if( vRC8900[khh] ) vRC8900[khh]->Fill( xh, yh );
//hh                         ----------------------------
	}
	int kMax2 = -1;
	double likeMax2 = 0.;
	for( int kl=0; kl<6; kl++ ) {
	  if( kl == 4 ) continue;
	  if( likeVar[kl] > likeMax2 ) {
	    likeMax2 = likeVar[kl];
	    kMax2 = kl;
	  }
	}
	if( kMax2 >= 0 ) {
	  khh = kh0 + 2;
	  xh = fact * likeKaon / likeMax2;
	  yh = momPart;
	  if( khh >= 0  &&  khh < int( vRC8900.size() ) ) {
	    if( vRC8900[khh] ) vRC8900[khh]->Fill( xh, yh );
//hh                           ----------------------------
	  }
	}
      }
      }

    }

    return;
  }
