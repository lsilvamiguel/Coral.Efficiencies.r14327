/*!
   \file    CsRCEventDisplay.cc
   \---------------------------
   \brief   CsRCEventDisplay class implementation.
   \author  Paolo Schiavon
   \version 0.02,  rev. 20/6/00, rev. October 2000
   \date    23 November 1999
*/


  #include <iostream>
  #include <ostream>
  #include <cstdio>

  #include "CsGeom.h"
  #include "CsHistograms.h"

//------------------------------
  #include "CsRCEventDisplay.h"

  #include "CsRCDetectors.h"
  #include "CsRCMirrors.h"

  #include "CsRichOne.h"
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

  #include "CsRCEventAnalysis.h"

  #include "CsRCExeKeys.h"
  #include "CsRCRecConst.h"

  #include "CsRCHistos.h"
  #include "CsRCUty.h"
//------------------------------

  using namespace std;
  using namespace CLHEP;

//===========================================================================
  CsRCEventDisplay::CsRCEventDisplay() {
//---------------------------------------
      flag_ = false;
  }

//===========================================================================
  CsRCEventDisplay::CsRCEventDisplay( const CsRCEventDisplay &disp ) {
//--------------------------------------------------------------------
    cout << "RICHONE : CsRCEventDisplay CopyConstructor" << endl;
    flag_ = disp.flag_;
  }

//===========================================================================
  void CsRCEventDisplay::print() {
//--------------------------------
    cout << endl;
    cout << " Event Display flag : " << flag_ << endl;
  }

//===========================================================================
  CsRCEventDisplay::~CsRCEventDisplay() {
//---------------------------------------
  }


//===========================================================================
  void CsRCEventDisplay::doEveDisplay() {
//---------------------------------------


//--- Procedure for event display.
//    ----------------------------
//--- Paolo  -  November 1999,  rev. October 2000
//              Rev. January 2004
//              -----------------
//              rev. August 2005
//              rev. May    2007


      CsRCRecConst *cons = CsRCRecConst::Instance();
      static int kOffHi = cons->kOffHi();
      static float CFRefInd = cons->CFRefInd();
      static int mcaScan = cons->mcaScan();

      static int nEveDSkip = cons->nEveDSkip();
      static int nEveDisplay = cons->nEveDisplay();

      CsRCHistos& hist = CsRCHistos::Ref();

//--- max = 15! beware space in memory
      static int nEveHisMax = 15;
//@@----------------------------
      static int jHisto;
      static int nevHis1;
      static int nevHis2;
      static int nEveDisTo;
      static int nRingDisTo;

      //^vector<CsHist2D*> vRC4000;
      //^int kH4000 = -1;
      //^vector<CsHist2D*> vRC4300;
      //^int kH4300 = -1;
      //^vector<CsHist1D*> vRC4600;
      //^int kH4600 = -1;
      //^vector<CsHist1D*> vRC4900;
      //^int kH4900 = -1;
      static const int nHFullA = 4000;
      vector<CsHist2D*> vRCFullA;
      int kHFullA = -1;
      static const int nHFullM = 4100; 
      vector<CsHist2D*> vRCFullM;
      int kHFullM = -1;
      static const int nHDetlA = 4200;
      vector<CsHist2D*> vRCDetlA;
      int kHDetlA = -1;
      static const int nHDetlM = 4400;
      vector<CsHist2D*> vRCDetlM;
      int kHDetlM = -1;
      static const int nHPattR = 4600;
      vector<CsHist1D*> vRCPattR;
      int kHPattR = -1;
      static const int nHHData = 4900;
      vector<CsHist1D*> vRCHData;
      int kHHData = -1;

      CsRCExeKeys *key = CsRCExeKeys::Instance();
      CsRCDetectors *dets = CsRCDetectors::Instance();

      static bool firstCall = true;
      if( firstCall ) { 
        firstCall = false;
        key->acknoMethod( "CsRCEventDisplay" );

	list<CsRCCathode*> lCathodes = dets->lCathodes();

        jHisto = -1;
        nevHis1 = nEveDSkip;
        nevHis2 = nEveDSkip + nEveDisplay - 1;
	nEveDisTo = 0;
	nRingDisTo = 0;
      }


      int kh;
      double xh, yh, wh;

      int kEvent = CsRichOne::Instance()->kEvent();

      CsRCEventAnalysis* ana = CsRCEventAnalysis::Instance();

      bool MCarloEvent = key->MCarloEvent();

      CsHistograms::SetCurrentPath("/RICH/DISP");
//@@--------------------------------------------


      if( kEvent >= nevHis1  &&  kEvent <= nevHis2 ) {

      bool bHisto = false;
//@@---------------------

      int nPartMax = 15;
//@@-------------------
      CsRCEventParticles* evparts = CsRCEventParticles::Instance();
      list<CsRCParticle*> lParticles = evparts->lParticles();
      int nPartEv = lParticles.size();
      list<CsRCParticle*>::iterator ip;
      int nPartUse = 0;
      for( ip=lParticles.begin(); ip!=lParticles.end(); ip++ )
	if( (*ip)->flag() ) nPartUse++;
      if( nPartEv > nPartMax ) nPartEv = nPartMax;
      if( nPartUse > nPartMax ) nPartUse = nPartMax;
//@@-----------------------------------------------

      int nPaPhoMax = 10;
//@@--------------------
      CsRCEventPartPhotons* evpapho = CsRCEventPartPhotons::Instance();
      list<CsRCPartPhotons*> lPartPhotons = evpapho->lPartPhotons();
      int nPaPhoEv = lPartPhotons.size();
      list<CsRCPartPhotons*>::iterator it;
      int nPaPhoUse = 0;
      for( it=lPartPhotons.begin(); it!=lPartPhotons.end(); it++ )
	if( (*it)->flag() ) nPaPhoUse++;
      if( nPaPhoEv > nPaPhoMax ) nPaPhoEv = nPaPhoMax;
      if( nPaPhoUse > nPaPhoMax ) nPaPhoUse = nPaPhoMax;
//@@---------------------------------------------------

      int nRingMax = 10;
//@@-------------------
      CsRCEventRings* evrings = CsRCEventRings::Instance();
      list<CsRCRing*> lRings = evrings->lRings();
      int nRingEv = lRings.size();
      list<CsRCRing*>::iterator ir;
      int nRingUse = 0;
      for( ir=lRings.begin(); ir!=lRings.end(); ir++ )
	if( (*ir)->flag() ) nRingUse++;
      if( nRingEv > nRingMax ) nRingEv = nRingMax;
      if( nRingUse > nRingMax ) nRingUse = nRingMax;
      int nRingDis = nRingUse;
      bool flagRingDis[nRingMax];
//@@-----------------------------------------------


//--- SET DISPLAY CONDITIONS
//    ----------------------
      if( MCarloEvent ) {
//    -----------------

//----- MONTECARLO run
//      --------------
	nRingDis = 0;
	int kRing = -1;
	for( ir=lRings.begin(); ir!=lRings.end(); ir++ ) {
	  kRing++;
	  flagRingDis[kRing] = false;
	  if( !(*ir)->flag() ) continue;

//------- Set mixed PMT and APV condition
	  list<CsRCPhoton*> lPhotons = (*ir)->lPhotons();
	  list<CsRCPhoton*>::iterator ih;
	  int nPhoRing = lPhotons.size();
	  bool mixed = false;
	  for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	    if( (*ih)->isPMT() != lPhotons.front()->isPMT() ) {
	      mixed = true;
	      break;
	    }
	  }
	  if( !mixed ) continue;
//@@      ---------------------

	  if( ana->flagMCMon() ) {
//--------- MCMonitor executed :
//          --------------------
            if( (*ir)->flagOverThrs() ) {
              bool dispCond = false;
              if( cons->dispMode() == "ALL" ) {
		flagRingDis[kRing] = true;
		nRingDis++;
		dispCond = true;
	      }
              if( cons->dispMode() == "REC" ) {
                if( (*ir)->flagReco() ) {
		  flagRingDis[kRing] = true;
		  nRingDis++;
		  dispCond = true;
		}
	      }
              if( cons->dispMode() == "NOR" ) {
                if( !(*ir)->flagReco() ) {
		  flagRingDis[kRing] = true;
		  nRingDis++;
		  dispCond = true;
		}
	      }

              if( dispCond ) bHisto = true;

            }   /* end if on flagOverThrs */
	  } else {
//--------- MCMonitor not executed :
//          ------------------------
	    flagRingDis[kRing] = true;
	    nRingDis++;
            bHisto = true;          //   display all
	  }

        }   /* end loop on Rings */

      } else  {

//----- DATA run :
//      ----------
        bool dispCond = false;
        //if( cons->dispMode() == "ALL" ) dispCond = true;
        if( cons->dispMode() == "ALL" ) {
	  nRingDis = 0;
	  int kRing = -1;
	  for( ir=lRings.begin(); ir!=lRings.end(); ir++ ) {
	    kRing++;
	    flagRingDis[kRing] = true;
	    nRingDis++;
	    dispCond = true;
	  }
	  //^if( lRings.size() == 0 ) dispCond = true;      //^
	}
        if( cons->dispMode() == "REC" ) {
	  nRingDis = 0;
	  int kRing = -1;
	  for( ir=lRings.begin(); ir!=lRings.end(); ir++ ) {
	    kRing++;
	    flagRingDis[kRing] = false;
	    if( !(*ir)->flag() ) continue;
	    if( !(*ir)->flagReco() ) continue;

//--------- further selection :
//          -------------------
	    CsRCParticle* part = (*ir)->pPartPhot()->pPart();
	    if( !part ) continue;
	    double ddPaDet = part->ddPaDet()[part->kDetPart()];
	    //if( !(ddPaDet > 800.) ) continue;
	    int nPhoRing = (*ir)->lPhotons().size();
	    //if( !(nPhoRing >=  6) ) continue;
	    //if( !(nPhoRing >=  9) ) continue;
	    //if( !(nPhoRing >=  30) ) continue;
	    //if( !(abs( (*ir)->kRingLoc() ) == 1) ) continue;
	    //if( !((*ir)->kRingLoc() == 0) ) continue;
//@@        -----------------------------------------------

	    flagRingDis[kRing] = true;
	    nRingDis++;
	    dispCond = true;
	    //break;
	  }
	}
	if( cons->dispMode() == "NOR" ) {
	  nRingDis = 0;
	  int kRing = -1;
	  for( ir=lRings.begin(); ir!=lRings.end(); ir++ ) {
	    kRing++;
	    flagRingDis[kRing] = false;
	    if( !(*ir)->flag() ) continue;
	    if( (*ir)->flagReco() ) continue;
//------------------------------------------
	    flagRingDis[kRing] = true;
	    nRingDis++;
            dispCond = true;
            //break;
          }
        }
	if( nRingDis > nRingMax ) nRingDis = nRingMax;

        if( dispCond ) bHisto = true;
      }


      if( bHisto ) {

	jHisto++;
	if( jHisto < nEveHisMax ) {
//@@    -------------------------
	  nEveDisTo++;
	  nRingDisTo += nRingDis;

//------- histo of event data :
//        ---------------------
	  int nBuff = 40;
	  //int nData = 20;
	  int nData = 22;                             //   040920
	  int nChanDa = nBuff + nRingMax * nData;
          //^int khData = kOffHi + 4900 + jHisto;
          int khData = kOffHi + nHHData + jHisto;
          stringstream name9, title9;
          name9 << khData;
          title9 << "   ";
          vRCHData.push_back( new CsHist1D( name9.str(), title9.str(),
//hh      -----------------------------------------------------------
					    nChanDa, 0. ,float( nChanDa ) ) );
          kHHData++;


//------- histo of event pad pattern :
//        ----------------------------

	  //^int nChanHs = 320;
	  //^float padxw = 8.;
//@@----------------------
	  //^float limitHs = (nChanHs*padxw)/2.;
	  //^//^float wind = 320.;
 	  //^float wind = 336.;
	  //^float ddYY = 204;
          //^int khDisp = kOffHi + 4000 + jHisto*10;
          //^stringstream name0, title0;
          //^name0 << khDisp;
          //^title0 << "   ";
          //^vRC4000.push_back( new CsHist2D( name0.str(), title0.str(),
//hh      -----------------------------------------------------------
	  //^nChanHs, -limitHs, limitHs, nChanHs, -limitHs, limitHs ) );
          //^kH4000++;

          int nChanHs, khDisp;
	  float limitHs = 1280.;
	  float wind = 336.;
	  float ddYY = 204.;
//------- APV part :
//        ----------
	  nChanHs = 320;
	  //padxw = 8.;
	  //limitHs = (nChanHs*padxw)/2.;
	  //wind = 320.;
          khDisp = kOffHi + nHFullA + jHisto;
          stringstream name0A, title0A;
          name0A << khDisp;
          title0A << "   ";
          vRCFullA.push_back( new CsHist2D( name0A.str(), title0A.str(),
//hh      --------------------------------------------------------------
	  nChanHs, -limitHs, limitHs, nChanHs, -limitHs, limitHs ) );
          kHFullA++;
//------- MAPMT part :
//        ------------
	  nChanHs = 214;
	  //padxw = 12.;
	  //limitHs = (nChanHs*padxw)/2.;
 	  //wind = 336.;
          khDisp = kOffHi + nHFullM + jHisto;
          stringstream name0M, title0M;
          name0M << khDisp;
          title0M << "   ";
          vRCFullM.push_back( new CsHist2D( name0M.str(), title0M.str(),
//hh      --------------------------------------------------------------
	  nChanHs, -limitHs, limitHs, nChanHs, -limitHs, limitHs ) );
          kHFullM++;


//------- plot event pads :
//        -----------------
          float PHPad = 1.;
//@@----------------------
	  list<CsRCPad*> lPadAll = CsRCEventPads::Instance()->lPads();
	  list<CsRCPad*>::iterator ia;
	  for( ia=lPadAll.begin(); ia!=lPadAll.end(); ia++ ) {

            int iCat = (*ia)->ic();
            int ixPa = (*ia)->ix();
            int iyPa = (*ia)->iy();
	    //PHPad = (*ia)->PH();
	    CsRCCathode* cat = dets->ptrToCat( iCat );
	    double padx = cat->padx();
	    double pady = cat->pady();
	    //std::cout << "CsRCEventDisplay : " << iCat << "  " << padx
	    //          << "  " << pady << std::endl;
	    double hCatx = cat->hCatx();
	    double hCaty = cat->hCaty();
	    //std::cout << "CsRCEventDisplay : " << iCat << "  " << ixPa
	    //	  << "  " << iyPa << "  "
	    //	  << hCatx << "  " << hCaty
	    //	  << std::endl;
            double xde = (ixPa*padx - hCatx) + dets->vOffCatW( iCat ).x();
            double yde = (iyPa*pady - hCaty) + dets->vOffCatW( iCat ).y();
	    xde += padx/2.;
	    yde += pady/2.;
            xh = xde;
            yh = yde;
	    if( yh >= 0. ) yh -= ddYY;
	    if( yh <  0. ) yh += ddYY;
            wh = PHPad;
            //^vRC4000[kH4000]->Fill( xh, yh, wh );
//hh        -----------------------------------
            if( cat->isAPV() ) vRCFullA[kHFullA]->Fill( xh, yh, wh );
//hh                           -------------------------------------
            if( cat->isPMT() ) vRCFullM[kHFullM]->Fill( xh, yh, wh );
//hh                           -------------------------------------
          }   /* end loop on Pads */

//------- plot event clusters (OPTIONAL) :
//        --------------------------------
          float PHClu = 100.;
//@@------------------------
	  list<CsRCCluster*> lCluAll =
	    CsRCEventClusters::Instance()->lClusters();
	  list<CsRCCluster*>::iterator ic;
	  for( ic=lCluAll.begin(); ic!=lCluAll.end(); ic++ ) {

	    int iCat = (*ic)->ic();
	    CsRCCathode* cat = dets->ptrToCat( iCat );
	    float xClu = (*ic)->xc();
	    float yClu = (*ic)->yc();
	    double xde = xClu;
	    double yde = yClu;
	    xh = xde;
	    yh = yde;
	    if( yh >= 0. ) yh -= ddYY;
	    if( yh <  0. ) yh += ddYY;
	    wh = PHClu;
            //vRC4000[kH4000]->Fill( xh, yh, wh );      //  !!!
//hh        -----------------------------------
            //if( cat->isAPV() ) vRCFullA[kHFullA]->Fill( xh, yh, wh );
//hh                           -------------------------------------
            //if( cat->isPMT() ) vRCFullM[kHFullM]->Fill( xh, yh, wh );
//hh                           -------------------------------------
          }   /* end loop on Clusters */


//------- record event data :
//        -------------------
          xh = 0.5;
          wh = float( kEvent );
          vRCHData[kHHData]->Fill( xh, wh );
//hh      ---------------------------------
          xh = 1.5;
          wh = float( nPartEv );
          //wh = float( nPartUse );
          vRCHData[kHHData]->Fill( xh, wh );
//hh      ---------------------------------
          xh = 2.5;
          wh = float( nPaPhoUse );
          vRCHData[kHHData]->Fill( xh, wh );
//hh      ---------------------------------
          xh = 3.5;
          //wh = float( nRingEv );
          wh = float( nRingUse );
          vRCHData[kHHData]->Fill( xh, wh );
//hh      ---------------------------------
          xh = 4.5;
	  wh = float( nRingDis );
          vRCHData[kHHData]->Fill( xh, wh );
//hh      ---------------------------------
          xh = 5.5;
          wh = float( nChanHs );
          vRCHData[kHHData]->Fill( xh, wh );
//hh      ---------------------------------
          xh = 6.5;
          wh = limitHs;
          vRCHData[kHHData]->Fill( xh, wh );
//hh      ---------------------------------
          xh = 7.5;
          wh = wind;
          vRCHData[kHHData]->Fill( xh, wh );
//hh      ---------------------------------
          xh = 8.5;
          wh = ddYY;
          vRCHData[kHHData]->Fill( xh, wh );
//hh      ---------------------------------

          xh = 8.5;
          for( ip=lParticles.begin(); ip!=lParticles.end(); ip++ ) {
	    if( (*ip)->vPade().size() == 0 ) continue;
	    int kDetPart = (*ip)->kDetPart();
	    double xPaDet = (*ip)->vPade()[kDetPart].x();
	    xh += 1.;
	    wh = xPaDet;
	    vRCHData[kHHData]->Fill( xh, wh );
//hh        ---------------------------------
	    double yPaDet = (*ip)->vPade()[kDetPart].y();
	    xh += 1.;
	    wh = yPaDet;
	    if( wh >= 0. ) wh -= ddYY;
	    if( wh <  0. ) wh += ddYY;
	    vRCHData[kHHData]->Fill( xh, wh );
//hh        ---------------------------------
	  }

//------- plot event rings cluster pattern :
//        ----------------------------------
	  int kRingDis = -1;
	  int kRing = -1;
          for( ir=lRings.begin(); ir!=lRings.end(); ir++ ) {
	    kRing++;
	    if( !(*ir)->flag() ) continue;

	    if( cons->dispMode() == "REC" ) if( !(*ir)->flagReco() ) continue;
//          -----------------------------------------------------------------
	    if( cons->dispMode() == "NOR" ) if( (*ir)->flagReco() ) continue;
//          ----------------------------------------------------------------
	    if( !flagRingDis[kRing] ) continue;
//          ----------------------------------
	    kRingDis++;

	    CsRCParticle* part = (*ir)->pPartPhot()->pPart();
            double theRing = (*ir)->the();
	    //^^list<CsRCPhoton*> lPhotons = (*ir)->lPhotons();
	    list<CsRCPhoton*> lPhotons = (*ir)->pPartPhot()->lPhotons();
	    list<CsRCPhoton*>::iterator ih;
	    int nPhoRing = lPhotons.size();

	    float PHPho = 100.;
//@@--------------------------
              double raRing = 0.;
              double raRingPMT = 0.;
              int kPhoUp = 0;
              int kPhoDw = 0;
	      int nPho = 0;
	      int nPhoPMT = 0;
	      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {

		CsRCCluster* clu = (*ih)->ptToClu();
                int iCat = clu->ic();
                float xClu = clu->xc();
                float yClu = clu->yc();
		//std::cout << "CsRCEventDisplay : " << iCat << "  " << xClu
		//          << "  " << yClu << "  " << std::endl;
                double xde = xClu;
		double yde = yClu;
                xh = xde;
                yh = yde;
		if( yh >= 0. ) yh -= ddYY;
		if( yh <  0. ) yh += ddYY;
                wh = PHPho;
		//^vRC4000[kH4000]->Fill( xh, yh, wh );
//hh            -----------------------------------

                if( part->vPade().size() == 0 ) continue;
                double xPaDet;
                double yPaDet;
                if( yde >= 0. )  {
		  kPhoUp++;
                  xPaDet = part->vPade()[0].x();
                  yPaDet = part->vPade()[0].y();
                }
                else  {
		  kPhoDw++;
                  xPaDet = part->vPade()[1].x();
                  yPaDet = part->vPade()[1].y();
                }
		if( (*ih)->isPMT() ) {
		  nPhoPMT++;
                  raRingPMT += sqrt( (xde-xPaDet)*(xde-xPaDet) +
				     (yde-yPaDet)*(yde-yPaDet) );
		}
		if( !(*ih)->isPMT() ) {
		  nPho++;
                  raRing += sqrt( (xde-xPaDet)*(xde-xPaDet) +
				  (yde-yPaDet)*(yde-yPaDet) );
		}

              }   /* end loop on Photons */

//----------- reconstructed ring radius
	      if( nPho > 0 ) {
		raRing /= float( nPho );
	      }
	      if( nPhoPMT > 0 ) {
		raRingPMT /= float( nPhoPMT );
	      }
	      if( nPhoPMT > nPho ) raRing = raRingPMT;

              double theCalc = 0.;
              double raCalc = 0.;
	      float partMom = part->mom();
              if( MCarloEvent )  {
                float partMass = part->paMass();
                double betaCalc = partMom / 
                      sqrt( partMom*partMom + partMass*partMass );
                double cosTheW = 1. / ( betaCalc * CFRefInd );
                if( cosTheW > 1. ) cosTheW = 1.;
                theCalc = acos( cosTheW ) * 1000.;
//------------- MC calculated ring radius (approximate)
                if( theRing > 0. ) raCalc = raRing * theCalc / theRing;
              }

	      double raPion = 0.;
	      double raKaon = 0.;
	      double raProton = 0.;
	      double thePion = 0.;
	      double theKaon = 0.;
	      double theProton = 0.;
              double raIpo = 0.;
	      for( int kPaTy=8; kPaTy <= 14; kPaTy+=3 ) {
		double theIpo = (*ir)->pPartPhot()->thetaIpo( kPaTy );
		//^if( dets->ptrToCat( part->getClosestCat() )->isPMT() )
		//^  theIpo = (*ir)->pPartPhot()->thetaIpoVS( kPaTy );
		if( theIpo >= 0. ) {
//--------------- mass hypo calculated ring radius (approximate)
		  if( theRing > 0. ) raIpo = raRing * theIpo / theRing;
		  if( kPaTy ==  8 ) {
		    raPion = raIpo;
		    thePion = theIpo;
		  }
		  if( kPaTy == 11 ) {
		    raKaon = raIpo;
 		    theKaon = theIpo;
		  }
		  if( kPaTy == 14 ) {
		    raProton = raIpo;
		    theProton = theIpo;
		  }
		}
	      }
	      //double probPion = (*ir)->partProbs()[1];
	      //double probKaon = (*ir)->partProbs()[2];
	      //double probProton = (*ir)->partProbs()[3];
	      double probPion = (*ir)->pPartPhot()->partProbs()[1];
	      double probKaon = (*ir)->pPartPhot()->partProbs()[2];
	      double probProton = (*ir)->pPartPhot()->partProbs()[3];
	      //std::cout << "Display : " << probPion << "  " << probKaon
	      //	  << "  " << probProton << std::endl;

              double xPaDetW = 1000000.;                      //   040223
              double yPaDetW = 1000000.;                      //   040223
	      int kDetPart = part->kDetPart();                //   040223
	      Hep3Vector vCorPoPhoWw = 
                (*ir)->pPartPhot()->vCorPoPhoW()[kDetPart];   //   040223
	      bool bPade = part->vPade().size() > 0;
	      //int iChan = nBuff + (kRing+1) * nData;
	      int iChan = nBuff + (kRingDis+1) * nData;
	      if( iChan > nChanDa ) break;

//----------- record event rings data :
//            -------------------------
              //xh = nBuff + kRing * nData + 0.5;
              xh = nBuff + kRingDis * nData + 0.5;
//            -----------------------------------
              if( kPhoUp >= kPhoDw ) {  
                if( bPade ) xPaDetW = part->vPade()[0].x();
		//if( bPade ) xPaDetW = part->vPade()[kDetPart].x(); // 040223
                //xPaDetW -= vCorPoPhoWw.x();                        // 040223
                wh = xPaDetW;
                vRCHData[kHHData]->Fill( xh   , wh );
//hh            ------------------------------------
                if( bPade ) yPaDetW = part->vPade()[0].y();
		//if( bPade ) yPaDetW = part->vPade()[kDetPart].y(); // 040223
                //yPaDetW -= vCorPoPhoWw.y();                        // 040223
                wh = yPaDetW;
		if( wh >= 0. ) wh -= ddYY;
		if( wh <  0. ) wh += ddYY;
                vRCHData[kHHData]->Fill( xh+1., wh );
//hh            ------------------------------------
              }
              else {
                if( bPade ) xPaDetW = part->vPade()[1].x();
		//if( bPade ) xPaDetW = part->vPade()[kDetPart].x(); // 040223
                //xPaDetW -= vCorPoPhoWw.x();                        // 040223
                wh = xPaDetW;
                vRCHData[kHHData]->Fill( xh   , wh );
//hh            ------------------------------------
                if( bPade ) yPaDetW = part->vPade()[1].y();
		//if( bPade ) yPaDetW = part->vPade()[kDetPart].y(); // 040223
                //yPaDetW -= vCorPoPhoWw.y();                        // 040223
                wh = yPaDetW;
		if( wh >= 0. ) wh -= ddYY;
		if( wh <  0. ) wh += ddYY;
                vRCHData[kHHData]->Fill( xh+1., wh );
//hh            -------------------------------------
	      }
	      //cout << xPaDetW << "  " << yPaDetW << endl;

              wh = raRing;
              vRCHData[kHHData]->Fill( xh+2., wh );
//hh          ------------------------------------
              wh = raCalc;
              vRCHData[kHHData]->Fill( xh+3., wh );
//hh          ------------------------------------
              wh = theRing;
              vRCHData[kHHData]->Fill( xh+4., wh );
//hh          ------------------------------------
              wh = theCalc;
              vRCHData[kHHData]->Fill( xh+5., wh );
//hh          ------------------------------------

              wh = raPion;
              vRCHData[kHHData]->Fill( xh+6., wh );
//hh          ------------------------------------
              wh = raKaon;
              vRCHData[kHHData]->Fill( xh+7., wh );
//hh          ------------------------------------
              wh = raProton;
              vRCHData[kHHData]->Fill( xh+8., wh );
//hh          ------------------------------------
              wh = thePion;
              vRCHData[kHHData]->Fill( xh+9., wh );
//hh          ------------------------------------
              wh = theKaon;
              vRCHData[kHHData]->Fill( xh+10., wh );
//hh          -------------------------------------
              wh = theProton;
              vRCHData[kHHData]->Fill( xh+11., wh );
//hh          -------------------------------------
              wh = probPion;
              vRCHData[kHHData]->Fill( xh+12., wh );
//hh          -------------------------------------
              wh = probKaon;
              vRCHData[kHHData]->Fill( xh+13., wh );
//hh          -------------------------------------
              wh = probProton;
              vRCHData[kHHData]->Fill( xh+14., wh );
//hh          -------------------------------------
              wh = partMom;
              vRCHData[kHHData]->Fill( xh+15., wh );
//hh          -------------------------------------
              wh = float( kPhoUp );
              vRCHData[kHHData]->Fill( xh+16., wh );
//hh          -------------------------------------
              wh = float( kPhoDw );
              vRCHData[kHHData]->Fill( xh+17., wh );
//hh          -------------------------------------

              float xRec = 0.;
//   MC-Monitor, DATA-Monitor
              if( !(*ir)->flagReco() ) xRec = 1.;
//   MC-Monitor only
              if( MCarloEvent ) if( !(*ir)->flagOverThrs() ) xRec = 2.;
	      if( (*ir)->kRingLoc() == 0 ) xRec += 10;

              wh = xRec;
              vRCHData[kHHData]->Fill( xh+18., wh );
//hh          -------------------------------------

	      double yRgLimUp = 0.;
	      double yRgLimDw = 0.;
	      if( (*ir)->pPartPhot()->vySplitLLim().size() > 0 ) {
		yRgLimUp = (*ir)->pPartPhot()->vySplitLLim()[0] - ddYY;
		yRgLimDw = (*ir)->pPartPhot()->vySplitLLim()[1] + ddYY;
	      }
              wh = yRgLimUp;
              vRCHData[kHHData]->Fill( xh+19., wh );
//hh          -------------------------------------
              wh = yRgLimDw;
              vRCHData[kHHData]->Fill( xh+20., wh );
//hh          -------------------------------------
              //std::cout <<  yRgLimUp << "  " << yRgLimDw << std::endl;

//----------- histo of ring cluster pattern :
//            -------------------------------
              float hLim = xPaDetW;
              float xhLim = hLim - wind;
              float xhLix = hLim + wind;
              if( xhLim < -limitHs)  {
                xhLim = -limitHs;
                xhLix = -limitHs + 2.*wind;
              }
              if( xhLix > limitHs )  {
                xhLim =  limitHs - 2.*wind;
                xhLix =  limitHs;
	      }
              hLim = yPaDetW;
	      if( hLim >= 0. ) hLim -= ddYY;
	      if( hLim <  0. ) hLim += ddYY;
              float yhLim = hLim - wind;
              float yhLix = hLim + wind;
              if( yhLim < -limitHs )  {
                yhLim = -limitHs;
                yhLix = -limitHs + 2.*wind;
      	      }
              if( yhLix >  limitHs )  {
                yhLim =  limitHs - 2.*wind;
                yhLix =  limitHs;
              }
	      //^//^int nChw = int( 2.* wind / 8.);      //^
	      //^int nChw = int( 2.* wind / 12.);      //^
              //^//int khDisZ = kOffHi + 4300 + jHisto*10 + kRing;
              //^int khDisZ = kOffHi + 4300 + jHisto*10 + kRingDis;
              //^stringstream name3, title3;
              //^name3 << khDisZ;
              //^title3 << "  ";
              //^vRC4300.push_back( new CsHist2D( name3.str(), title3.str(),
//hh          -----------------------------------------------------------
	      //^//^100, xhLim, xhLix,100, yhLim, yhLix ) );
              //^nChw, xhLim, xhLix,nChw, yhLim, yhLix ) );
              //^kH4300++;
	      int nChw, khDisZ;
//----------- APV part :
//            ----------
	      nChw = int( 2.* wind / 8.);
              khDisZ = kOffHi + nHDetlA + jHisto*10 + kRingDis;
              stringstream name3A, title3A;
              name3A << khDisZ;
              title3A << "  ";
              vRCDetlA.push_back( new CsHist2D( name3A.str(), title3A.str(),
//hh          --------------------------------------------------------------
	      nChw, xhLim, xhLix, nChw, yhLim, yhLix ) );
              kHDetlA++;
//----------- MAPMT part :
//            ------------
	      nChw = int( 2.* wind / 12.);
              khDisZ = kOffHi + nHDetlM + jHisto*10 + kRingDis;
              stringstream name3M, title3M;
              name3M << khDisZ;
              title3M << "  ";
              vRCDetlM.push_back( new CsHist2D( name3M.str(), title3M.str(),
//hh          --------------------------------------------------------------
              nChw, xhLim, xhLix, nChw, yhLim, yhLix ) );
              kHDetlM++;

//----------- plot ring pad pattern :
//            -----------------------
              float PHPad = 1.;
//@@--------------------------
	      for( ia=lPadAll.begin(); ia!=lPadAll.end(); ia++ ) {
                int iCat = (*ia)->ic();
                int ixPa = (*ia)->ix();
                int iyPa = (*ia)->iy();
		CsRCCathode* cat = dets->ptrToCat( iCat );
		double padx = cat->padx();
		double pady = cat->pady();
		double hCatx = cat->hCatx();
		double hCaty = cat->hCaty();
                double xde = (ixPa*padx-hCatx) + dets->vOffCatW( iCat ).x();
                double yde = (iyPa*pady-hCaty) + dets->vOffCatW( iCat ).y();
		xde += padx/2.;
		yde += pady/2.;
                xh = xde;
                yh = yde;
	        if( yh >= 0. ) yh -= ddYY;
	        if( yh <  0. ) yh += ddYY;
                wh = PHPad;
                //^vRC4300[kH4300]->Fill( xh, yh, wh );
//hh            -----------------------------------

              }   /* end loop on Pads */

//----------- plot ring cluster pattern :
//            ---------------------------
//            (or ALL clusters)
              PHPho = 100.;
//@@----------------------
   	      list<CsRCCluster*> lCluAll =
	        CsRCEventClusters::Instance()->lClusters();
	      list<CsRCCluster*>::iterator ic;
	      //^for( ic=lCluAll.begin(); ic!=lCluAll.end(); ic++ ) {
	        //^CsRCCluster* clu = (*ic);
	      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	        CsRCCluster* clu = (*ih)->ptToClu();
                int iCat = clu->ic();
		CsRCCathode* cat = dets->ptrToCat( iCat );
                float xClu = clu->xc();
                float yClu = clu->yc();
                double xde = xClu;
                double yde = yClu;
                xh = xde;
                yh = yde;
		if( yh >= 0. ) yh -= ddYY;
	        if( yh <  0. ) yh += ddYY;
                wh = PHPho;
                if( cat->isAPV() ) vRCDetlA[kHDetlA]->Fill( xh, yh, wh );
//hh                               -------------------------------------
                if( cat->isPMT() ) vRCDetlM[kHDetlM]->Fill( xh, yh, wh );
//hh                               -------------------------------------

              }   /* end loop on Ring Photons */


//----------- histo of ring pattern recognition :
//            -----------------------------------
              float *binCont;
	      binCont = (*ir)->binCont();
              //^int khBin = kOffHi + 4600 + jHisto*10 + kRingDis;
	      int khBin = kOffHi + nHPattR + jHisto*10 + kRingDis;
              stringstream name6, title6;
              name6 << khBin;
              title6 << "  ";
              vRCPattR.push_back( new CsHist1D( name6.str(), title6.str(),
//hh          ------------------------------------------------------------
						60, 0. ,60. ) );
     	      kHPattR++;

//----------- plot ring pattern recognition :
//            -------------------------------
              for ( int kBin=0; kBin < mcaScan; kBin++ )  {
                xh = float( kBin ) + 0.5;
                wh = binCont[kBin];
                vRCPattR[kHPattR]->Fill( xh, wh );
//hh            ---------------------------------
	      }

          }   /* end loop on Rings */


	  //sideDisplay( jHisto );
//        ---------------------

	  std::cout << std::endl;
  	  std::cout << "--- CsRCEventDisplay::doEveDisplay : "
		    << "Events Displayed  " << nEveDisTo 
		    << "  Rings Displayed  " << nRingDisTo << std::endl;

	}   /* end if on jHisto */

      }   /* end if on bHisto */

      }   /* end if on kEvent */

      CsHistograms::SetCurrentPath("/");
//@@-----------------------------------

      flag_ = true;

//--- conditional prints :
//    --------------------
      int kPrintEventDisplay = key->kPrintEventDisplay();
      if( kPrintEventDisplay == 1 ) print();

  }


//===========================================================================
  void CsRCEventDisplay::sideDisplay( int jHisto ) {
//--------------------------------------------------


//--- Procedure for event display.
//    ----------------------------
//--- Paolo  -  June 2001


      CsRCRecConst *cons = CsRCRecConst::Instance();
      static int kOffHi = cons->kOffHi();
      CsRCDetectors* dets = CsRCDetectors::Instance();
      CsRCMirrors* mirr = CsRCMirrors::Instance();

      vector<CsHist1D*> vRC5300;
      int kH5300 = -1;

      static int nHistoMax = 20;
      //int nData = 26;
      int nData = 40;                            //   top
      int nChanDa = nHistoMax * nData;

      int khData = kOffHi + 5300 + jHisto;
      stringstream name, title;
      name << khData;
      title << "   ";
      vRC5300.push_back( new CsHist1D( name.str(), title.str(),
//hh  ---------------------------------------------------------
      nChanDa, 0. ,float( nChanDa ) ) );
      kH5300++;

      double xh, yh, wh;

      CsRCEventRings* evrings = CsRCEventRings::Instance();
      list<CsRCRing*> lRings = evrings->lRings();
      list<CsRCRing*>::iterator ir;
      int kRing = -1;
      for( ir=lRings.begin(); ir!=lRings.end(); ir++ ) {
	//kRing++;
	if( !(*ir)->flag() ) continue;
        kRing++;                                  //   020906

	CsRCParticle* part = (*ir)->pPartPhot()->pPart();
	Hep3Vector vPoPart0 = part->vPosIn();
	Hep3Vector vDcPart0 = part->vDirIn();
	xh = kRing * nData + 0.5;
	wh = vPoPart0.z();
	vRC5300[kH5300]->Fill( xh   , wh );
//hh    ----------------------------------
	wh = vPoPart0.y();
	vRC5300[kH5300]->Fill( xh+1., wh );
//hh    ----------------------------------
	wh = vPoPart0.x();                           //   top
	vRC5300[kH5300]->Fill( xh+25., wh );         //   top
//hh    -----------------------------------
	double tanA = vDcPart0.x()/vDcPart0.z();
	double tanB = vDcPart0.y()/vDcPart0.z();

        list<CsRCMirrorNom*> lMirrNom = mirr->lMirrNom();
        list<CsRCMirrorNom*>::iterator in;
	int kMir = 0;
        for( in=lMirrNom.begin(); in!=lMirrNom.end(); in++ ) {
	  Hep3Vector vPoC0 = (*in)->vC0();
	  double RR = (*in)->RR();
//------- particle impact on mirrors (MRS) :
//        ----------------------------------
          Hep3Vector vPoPaMir0 =
	    mirr->vImpMir( vPoPart0, vDcPart0, vPoC0, RR );
//          ----------------------------------------------
//------- normal to mirror at particle impact :
//        ---------------------------------------
          Hep3Vector vDcNoPaMir0 = (1./RR) * (vPoPaMir0 - vPoC0);
//------- particle 'reflected' direction :
//        --------------------------------
          double cosPaMir = vDcNoPaMir0 * vDcPart0;
          Hep3Vector vDcPaRefl0 = 2.*cosPaMir * vDcNoPaMir0 - vDcPart0;
//------- particle 'impact' on detectors :
//        --------------------------------
	  Hep3Vector vDcDetvw = dets->vDcDetv( kMir );
          double norm = ( (dets->vDet0v( kMir ) - vPoPaMir0 ) * vDcDetvw ) /
                        ( vDcPaRefl0 * vDcDetvw );
          Hep3Vector vPoPaDet0 = vPoPaMir0 + norm * vDcPaRefl0;

	  wh = vPoPaMir0.z();
	  vRC5300[kH5300]->Fill( xh+2.+kMir, wh );
//hh      ---------------------------------------
	  wh = vPoPaMir0.y();
	  vRC5300[kH5300]->Fill( xh+4.+kMir, wh );
//hh      ---------------------------------------
	  wh = vPoPaMir0.x();                           //   top
	  vRC5300[kH5300]->Fill( xh+26.+kMir, wh );     //   top
//hh      ----------------------------------------
	  wh = vPoPaDet0.z();
	  vRC5300[kH5300]->Fill( xh+6.+kMir, wh );
//hh      ---------------------------------------
	  wh = vPoPaDet0.y();
	  vRC5300[kH5300]->Fill( xh+8.+kMir, wh );
//hh      ---------------------------------------
	  wh = vPoPaDet0.x();                           //   top
	  vRC5300[kH5300]->Fill( xh+28.+kMir, wh );     //   top
//hh      ----------------------------------------

//------- 'photon' impact on mirrors (MRS) and detectors :
//        ------------------------------------------------
	  double tanRing = tan( (*ir)->the()/1000. );
	  double tanPho = (tanB + tanRing) / (1. - tanB*tanRing);
          Hep3Vector vDcPho0u( tanA, tanPho, 1. );
	  vDcPho0u = vDcPho0u.unit();
          Hep3Vector vPhoMir0u =
	    mirr->vImpMir( vPoPart0, vDcPho0u, vPoC0, RR );
//          ---------------------------------------------
          Hep3Vector vDcNoPhoMir0u = (1./RR) * (vPhoMir0u - vPoC0);
          double cosPhoMir = vDcNoPhoMir0u * vDcPho0u;
          Hep3Vector vDcPhoRefl0u = 2.*cosPhoMir * vDcNoPhoMir0u - vDcPho0u;
          norm = ( (dets->vDet0v( kMir ) - vPhoMir0u ) * vDcDetvw ) /
	         ( vDcPhoRefl0u * vDcDetvw );
          Hep3Vector vPoPhoDet0u = vPhoMir0u + norm * vDcPhoRefl0u;

	  wh = vPhoMir0u.z();
	  vRC5300[kH5300]->Fill( xh+10.+kMir, wh );
//hh      ----------------------------------------
	  wh = vPhoMir0u.y();
	  vRC5300[kH5300]->Fill( xh+12.+kMir, wh );
//hh      ----------------------------------------
	  wh = vPhoMir0u.x();                           //   top
	  vRC5300[kH5300]->Fill( xh+30.+kMir, wh );     //   top
//hh      ----------------------------------------
	  wh = vPoPhoDet0u.z();
	  vRC5300[kH5300]->Fill( xh+14.+kMir, wh );
//hh      ----------------------------------------
	  wh = vPoPhoDet0u.y();
	  vRC5300[kH5300]->Fill( xh+16.+kMir, wh );
//hh      ----------------------------------------
	  wh = vPoPhoDet0u.x();                         //   top
	  vRC5300[kH5300]->Fill( xh+32.+kMir, wh );     //   top
//hh      ----------------------------------------

	  tanPho = (tanB - tanRing) / (1. + tanB*tanRing);
	  Hep3Vector vDcPho0d( tanA, tanPho, 1. );
	  vDcPho0d = vDcPho0d.unit();
          Hep3Vector vPhoMir0d =
	    mirr->vImpMir( vPoPart0, vDcPho0d, vPoC0, RR );
//          ---------------------------------------------
          Hep3Vector vDcNoPhoMir0d = (1./RR) * (vPhoMir0d - vPoC0);
          cosPhoMir = vDcNoPhoMir0d * vDcPho0d;
          Hep3Vector vDcPhoRefl0d = 2.*cosPhoMir * vDcNoPhoMir0d - vDcPho0d;
          norm = ( (dets->vDet0v( kMir ) - vPhoMir0d ) * vDcDetvw ) /
	         ( vDcPhoRefl0d * vDcDetvw );
          Hep3Vector vPoPhoDet0d = vPhoMir0d + norm * vDcPhoRefl0d;

	  wh = vPhoMir0d.z();
	  vRC5300[kH5300]->Fill( xh+18.+kMir, wh );
//hh      ----------------------------------------
	  wh = vPhoMir0d.y();
	  vRC5300[kH5300]->Fill( xh+20.+kMir, wh );
//hh      ----------------------------------------
	  wh = vPhoMir0d.x();                           //   top
	  vRC5300[kH5300]->Fill( xh+34.+kMir, wh );     //   top
//hh      ----------------------------------------
	  wh = vPoPhoDet0d.z();
	  vRC5300[kH5300]->Fill( xh+22.+kMir, wh );
//hh      ----------------------------------------
	  wh = vPoPhoDet0d.y();
	  vRC5300[kH5300]->Fill( xh+24.+kMir, wh );
//hh      ----------------------------------------
	  wh = vPoPhoDet0d.x();                         //   top
	  vRC5300[kH5300]->Fill( xh+36.+kMir, wh );     //   top
//hh      ----------------------------------------

	  kMir++;
	}

      }

  }
