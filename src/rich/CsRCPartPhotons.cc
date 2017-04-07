
/*!
   \file    CsRCPartPhotons.cc
   \--------------------------
   \brief   CsRCPartPhotons class implementation.
   \author  Paolo Schiavon
   \version 0.01
   \date    1 October 2000
*/


  #include <iostream>
  #include <ostream>
  #include <sstream>
  #include <cstdio>
  #include <cmath>

  #include <CLHEP/Vector/ThreeVector.h>
//-------------------------------------

  #include "CsErrLog.h"
  #include "CsMCHit.h"
  #include "CsMCRICH1Hit.h"
  #include "CsDigit.h"
  #include "CsMCDigit.h"
  #include "CsHelix.h"

  #include "CsStopwatch.h"

// ---------------------------
  #include "CsRCPartPhotons.h"

  #include "CsRichOne.h"
  #include "CsRCParticle.h"

  #include "CsRCPad.h"
  #include "CsRCEventClusters.h"
  #include "CsRCCluster.h"

  #include "CsRCDetectors.h"
  #include "CsRCPhotonDet.h"
  #include "CsRCCathode.h"
  #include "CsRCMirrors.h"

  #include "CsRCEventPartPhotons.h"
  #include "CsRCPhoton.h"
  #include "CsRCLikeAll.h"
  #include "CsRCLikeAll02.h"
  #include "CsRCLikeAll03.h"
  #include "CsRCLikeAll04.h"
  #include "CsRCLikeAll05.h"
  #include "CsRCLikeAllMap.h"
  #include "CsRCChiSqFit.h"
  #include "CsRCGauPolFit.h"

  #include "CsRCRing.h"

  #include "CsRCExeKeys.h"
  #include "CsRCRecConst.h"
  #include "CsRCHistos.h"
  #include "CsRCUty.h"

  #include "CsRCOptCorr.h"
// ---------------------------

  using namespace std;
  using namespace CLHEP;

//===========================================================================
  CsRCPartPhotons::CsRCPartPhotons() {
//------------------------------------
    clearPartPhotons();
  }

//===========================================================================
  CsRCPartPhotons::CsRCPartPhotons( const CsRCPartPhotons &partphot ) {
//---------------------------------------------------------------------
    //cout << "RICHONE : CsRCParPhotons CopyConstructor" << endl;
    kPaPhot_ = partphot.kPaPhot_;
    pPart_ = partphot.pPart_;
    pRing_ = partphot.pRing_;
    lPhotons_ = partphot.lPhotons_;
    lCluSignal_ = partphot.lCluSignal_;

    vPoPartW_ = partphot.vPoPartW_;
    vDcPartW_ = partphot.vDcPartW_;
    vPoPhotW_ = partphot.vPoPhotW_;
    vPoPaMirW_ = partphot.vPoPaMirW_;
    vCorPoPhoW_ = partphot.vCorPoPhoW_;
    vDcPaReflW_ = partphot.vDcPaReflW_;
    vPoPaDetW_ = partphot.vPoPaDetW_;
    vySplitLLim_ = partphot.vySplitLLim_;
    bSplitLLimSet_ = partphot.bSplitLLimSet_;

    //for( int j=1; j<=31; j++ ) thetaIpo_[j] = partphot.thetaIpo_[j];
    //for( int j=1; j<=31; j++ ) thetaIpoUV_[j] = partphot.thetaIpoUV_[j];
    //for( int j=1; j<=31; j++ ) thetaIpoVS_[j] = partphot.thetaIpoVS_[j];
//- corrected   100201   cppcheck
    for( int j=0; j<31; j++ ) thetaIpo_[j] = partphot.thetaIpo_[j];
    for( int j=0; j<31; j++ ) thetaIpoUV_[j] = partphot.thetaIpoUV_[j];
    for( int j=0; j<31; j++ ) thetaIpoVS_[j] = partphot.thetaIpoVS_[j];

    kDetPart_ = partphot.kDetPart_;
    pDetPart_ = partphot.pDetPart_;
    iCaPa_ = partphot.iCaPa_;
    iXpPa_ = partphot.iXpPa_;
    iYpPa_ = partphot.iYpPa_;
    lSelCats_ = partphot.lSelCats_;

    thetaLike_ = partphot.thetaLike_;
    pLikeMax_ = partphot.pLikeMax_;
    thetaLikeSet_ = partphot.thetaLikeSet_;
    fracUse_ = partphot.fracUse_;
    normSg_ = partphot.normSg_;
    normBk_ = partphot.normBk_;
    CsRCRecConst *cons = CsRCRecConst::Instance();
    int nProb = cons->outBufferSize();
    for( int j=0; j<nProb; j++ ) partProbs_[j] = partphot.partProbs_[j];
    partProbsSet_ = partphot.partProbsSet_;
    for( int j=0; j<31; j++ ) probaLKAll_[j] = partphot.probaLKAll_[j];
    for( int j=0; j<31; j++ ) derivLKAll_[j] = partphot.derivLKAll_[j];
    for( int j=0; j<31; j++ ) probaLKAllAPV_[j] = partphot.probaLKAllAPV_[j];
    for( int j=0; j<31; j++ ) probaLKAllPMT_[j] = partphot.probaLKAllPMT_[j];
    for( int j=0; j<31; j++ ) derivLKAllAPV_[j] = partphot.derivLKAllAPV_[j];
    for( int j=0; j<31; j++ ) derivLKAllPMT_[j] = partphot.derivLKAllPMT_[j];
    probaLKBgAll_ = partphot.probaLKBgAll_;
    nPhotAll_ = partphot.nPhotAll_;
    nPhotAllAPV_ = partphot.nPhotAllAPV_;
    nPhotAllPMT_ = partphot.nPhotAllPMT_;
    for( int j=0; j<31; j++ ) probaChiAll_[j] = partphot.probaChiAll_[j];

    nLikePro_ = partphot.nLikePro_;
    theLikeProMn_ = partphot.theLikeProMn_;
    theLikeProMx_ = partphot.theLikeProMx_;
    for( int j=0; j<nLikePro_; j++ ) theLikePro_[j] = partphot.theLikePro_[j];
    for( int j=0; j<nLikePro_; j++ ) pLikePro_[j] = partphot.pLikePro_[j];
    likeProSet_ = partphot.likeProSet_;

    likeONLY_ = partphot.likeONLY_;
    likeFirst_ = partphot.likeFirst_;
    pionONLY_ = partphot.pionONLY_;
    kaonONLY_ = partphot.kaonONLY_;

    likeAPVONLY_ = partphot.likeAPVONLY_;
    likePMTONLY_ = partphot.likePMTONLY_;

    flag_ = partphot.flag_;
  }

//===========================================================================
  void CsRCPartPhotons::clearPartPhotons() {
//------------------------------------------
    kPaPhot_ = -1;
    pPart_ = NULL;
    pRing_ = NULL;
    lPhotons_.clear();
    lCluSignal_.clear();

    vPoPartW_.clear();
    vDcPartW_.clear();
    vPoPhotW_.clear();
    vPoPaMirW_.clear();
    vCorPoPhoW_.clear();
    vDcPaReflW_.clear();
    vPoPaDetW_.clear();
    vySplitLLim_.clear();
    bSplitLLimSet_ = false;

    //for( int j=1; j<=31; j++ ) thetaIpo_[j] = -1.;
    //for( int j=1; j<=31; j++ ) thetaIpoUV_[j] = -1.;
    //for( int j=1; j<=31; j++ ) thetaIpoVS_[j] = -1.;
//- corrected   100201   cppcheck
    for( int j=0; j<31; j++ ) thetaIpo_[j] = -1.;
    for( int j=0; j<31; j++ ) thetaIpoUV_[j] = -1.;
    for( int j=0; j<31; j++ ) thetaIpoVS_[j] = -1.;

    kDetPart_ = -1;
    pDetPart_ = NULL;
    iCaPa_ = 0;
    iXpPa_ = 0;
    iYpPa_ = 0;
    lSelCats_.clear();

    thetaLike_ = 0.;
    pLikeMax_ = -100000.;
    thetaLikeSet_ = false;
    fracUse_ = 1.;
    normSg_ = 0.;
    normBk_ = 0.;
    CsRCRecConst *cons = CsRCRecConst::Instance();
    int nProb = cons->outBufferSize();
    for( int j=0; j<nProb; j++ ) partProbs_[j] = 0.;
    partProbsSet_ = false;
    for( int j=0; j<31; j++ ) probaLKAll_[j] = 0.;
    for( int j=0; j<31; j++ ) derivLKAll_[j] = 0.;
    for( int j=0; j<31; j++ ) probaLKAllAPV_[j] = 0.;
    for( int j=0; j<31; j++ ) probaLKAllPMT_[j] = 0.;
    for( int j=0; j<31; j++ ) derivLKAllAPV_[j] = 0.;
    for( int j=0; j<31; j++ ) derivLKAllPMT_[j] = 0.;
    probaLKBgAll_ = 0.;
    nPhotAll_ = 0;
    nPhotAllAPV_ = 0;
    nPhotAllPMT_ = 0;
    for( int j=0; j<31; j++ ) probaChiAll_[j] = 0.;

    nLikePro_ = 0;
    theLikeProMn_ = 0.;
    theLikeProMx_ = 0.;
    likeProSet_ = false;

    likeONLY_ = false;
    likeFirst_ = false;
    pionONLY_ = false;
    kaonONLY_ = false;
    likeAPVONLY_ = false;
    likePMTONLY_ = false;

    flag_ = false;
  }

//===========================================================================
  void CsRCPartPhotons::setCluSignal( CsRCCluster* ic ) {
//-------------------------------------------------------
    lCluSignal_.push_back( ic );
  }

//===========================================================================
  void CsRCPartPhotons::print() const {
//-------------------------------------
    cout << endl;
    cout << " PartPhotons  " << kPaPhot_ << " , Particle  "
	 << pPart_->kPart()
	 << " , from  photons : " << endl;
    list<CsRCPhoton*>::const_iterator ih;
    for( ih=lPhotons_.begin(); ih!=lPhotons_.end(); ih++ ) (*ih)->print();
  }

//===========================================================================
  CsRCPartPhotons::~CsRCPartPhotons() {
//-------------------------------------
    list<CsRCPhoton*>::iterator ih;
    for( ih=lPhotons_.begin(); ih!=lPhotons_.end(); ih++ ) delete *ih;
    lPhotons_.clear();
    //    for( is=lCluSignal_.begin(); ... ) delete *is;   //   NO
    lCluSignal_.clear();

    vDcPartW_.clear();
    vPoPhotW_.clear();
    vPoPaMirW_.clear();
    vCorPoPhoW_.clear();
    vDcPaReflW_.clear();
    vPoPaDetW_.clear();

    lSelCats_.clear();
  }


//===========================================================================
  CsRCPartPhotons::CsRCPartPhotons( int kPaPhot, CsRCParticle* part ) {
//---------------------------------------------------------------------


//--- Paolo - October 2000, rev. Sacha - March 2001


    clearPartPhotons();
//--------------------
    kPaPhot_ = kPaPhot;
    pPart_ = part;
    //pRing_ = NULL;
    //pDetPart_ = NULL;
    flag_ = true;

    doPartToMWR( part );
//---------------------

    doSelCathodes( part );
//-----------------------

    getPartPhotons( part );
//------------------------

    compThetaIpo( part );
//----------------------

  }


//===========================================================================
  void CsRCPartPhotons::getPartPhotons( CsRCParticle* part ) {
//------------------------------------------------------------


//--- 'Photon' reconstruction code.
//    -----------------------------
//--- Paolo  -  June 1999
//    from : RingRec-2.cc
//--- rev. : November 1999  from  ringrec-4.f
//--- rev. : October  2001
//--- rev. : April  2007


//--- from "CsRCExecKeys"
      CsRCExeKeys *key = CsRCExeKeys::Instance();
      bool CorrQzW = key->CorrQzW();
      bool CorrPMTOpt = key->CorrPMTOpt();

//--- from "ReconConst.h"
      CsRCRecConst *cons = CsRCRecConst::Instance();
      int kOffHi = cons->kOffHi();
      double TwoPI = cons->TwoPI();
      double RadDeg = cons->RadDeg();

      CsRCDetectors *dets = CsRCDetectors::Instance();
      CsRCMirrors *mirr = CsRCMirrors::Instance();

      CsRCHistos& hist = CsRCHistos::Ref();
      int kh;
      double xh, yh, wh;

      static int nPhoDet;
      static int nCathode;
      static int nHlfCat;

      static std::vector<CsHist2D*> vRC1100;
      static CsHist2D* hRC1110;

      static bool firstCall = true;
      if( firstCall ) {
        firstCall = false;

	nPhoDet = dets->lPhoDet().size();
        nCathode = dets->nCathode();
        nHlfCat = nCathode / 2;

	list<CsRCCathode*> lCathodes = dets->lCathodes();

//----- Monitoring histo :
	for( int kh=0; kh<5; kh++ ) vRC1100.push_back( NULL );
	hRC1110 = NULL;
        if( CsRCHistos::Ref().bookHis() &&
	    CsRCHistos::Ref().levelBk() >= 2 ) {
	  vRC1100.clear();
	  CsHistograms::SetCurrentPath("/RICH");
	  for( int kh=0; kh<4; kh++ ) {
	    stringstream hN1100;
	    int kHist = kOffHi + 1100 + kh + 1;
	    hN1100 << kHist;
	    string hTitle = " ";
	    vRC1100.push_back( new CsHist2D( hN1100.str(), hTitle,
//hh        ------------------------------------------------------
					     120, -3., 3., 70, 0., 70. ) );
	  }
	  stringstream hN1110;
	  int kHist = kOffHi + 1111;
	  hN1110 << kHist;
	  string hTitle = " ";
	  hRC1110 = new CsHist2D( hN1110.str(), hTitle,
//hh      ---------------------------------------------
				  120, -3., 3., 120, -3., 3. );
          CsHistograms::SetCurrentPath("/");
	}
      }

//--- tight cut on cathode background for split rings (added 01/10/18) :
//    ------------------------------------------------------------------
      vector<double> vySplitLLimit = getSplitLimits( part );
//                                   ----------------------

      lPhotons_.clear();

//--- SCAN selected cathodes for CLUSTERS :
//    -------------------------------------
      list<int>::iterator isc;
      CsRCEventClusters* clus = CsRCEventClusters::Instance();
      list<CsRCCluster*> lClusters = clus->lClusters();
      list<CsRCCluster*>::iterator icl;
      int kPhot = -1;
      for( icl=lClusters.begin(); icl!=lClusters.end(); icl++ ) {
	if( !(*icl)->flag() ) continue;

        int iCaw = (*icl)->ic();
        bool bProClu = false;
	for( isc=lSelCats_.begin(); isc!=lSelCats_.end(); isc++ ) {
          if( iCaw == (*isc) ) { bProClu = true;  break; }
	}
//----- scan ALL cathodes :
//      -------------------
	//bProClu = true;

        int kDetClu = -1;
	kDetClu = dets->cathodePos( iCaw );

	double xw = (*icl)->xc();
	double yw = (*icl)->yc();
	double zw = dets->vOffCatW( iCaw ).z();
	double zDetW = zw;
        Hep3Vector vPoCluDet( xw, yw, zw );

	int kDetPart = part->kDetPart();
	double ddPaDet = part->ddPaDet()[kDetPart];

//----- tight cut on cathode background for split rings (added 01/10/18) :
//      ------------------------------------------------------------------
	if( kDetClu != kDetPart ) {
  	  if( fabs( yw ) < fabs( vySplitLLimit[kDetClu] ) ) bProClu = false;
//@@                                                        ---------------
	}

//----- reject MIP clusters :
//      ---------------------
	if( MIPReject( (*icl) ) ) bProClu = false;
//@@                              ---------------


//----- kill 'halo' pads on detector :
//      ------------------------------
	//if( key->KillHalo() ) {
	//  if( killHalo( (*icl) ) ) bProClu = false;      //   dummy!
//@@                                 ---------------
	//}

//----- use or not the less important part of a split ring :
//      (central detector region only)                040203
//      ----------------------------------------------------
	if( !key->UseSplitRings() ) {
	  if( ddPaDet < 300.) {
            if( !(kDetClu == kDetPart) ) bProClu = false;
//@@                                     ---------------
	  }
	}

//----- full ring only for split rings (test 040924) :
//      ----------------------------------------------
	//if( !key->UseSplitRings() ) {
	//  if( ddPaDet < 300.) {
	//    if( kDetClu == kDetPart ) {
	//      if( fabs( yw ) < fabs( vySplitLLimit[kDetClu] ) ) 
	//	  bProClu = false;
//@@    //        ---------------
	//    }
	//  }
	//}


        if( bProClu ) {

          xh = xw;
          yh = yw;
          if( hist.hRC3506 ) hist.hRC3506->Fill( xh, yh );
//hh                         ----------------------------

//------- particle-detector angle :                   //   moved
//        -------------------------
          //Hep3Vector vDcPaReflWw = vDcPaReflW[kDetClu];
          //double thePaDet = acos( vDcPaReflWw.z() );

//------- photon positions :
//        ------------------
          float RR = mirr->RRv( kDetClu );
          Hep3Vector vPoCluDetw = vPoCluDet;
          Hep3Vector vPoPhotWw = vPoPhotW()[kDetClu];

/*-----------------------------------------------------------------------
//------- photon emission positions from MC (complicated...) :
//        ----------------------------------------------------
	  //cout << "---" << endl;
	  CsMCTrack* paMCTrack = part->pMCTrack();
	  Hep3Vector vPoCluDetMC( 0. );
          Hep3Vector vPosEmMC( 0. );
	  bool done = false;
          list<CsRCPad*> lPads = (*icl)->lPads();
          list<CsRCPad*>::iterator ip;
	  list<CsMCHit*> cluHits;
	  for( ip=lPads.begin(); ip!=lPads.end(); ip++ ) {
	    CsDigit* digi = (*ip)->pDigit();
	    cluHits = dynamic_cast<CsMCDigit*>(digi)->CsMCDigit::getHits();
	    list<CsMCHit*>::iterator ch;
	    //- take the last one...
	    for( ch=cluHits.begin(); ch!=cluHits.end(); ch++ ) {
	      CsMCTrack* phMCTrack = (*ch)->getMCTrack();
	      if( phMCTrack == paMCTrack ) {
	        CsMCRICH1Hit* hit = dynamic_cast<CsMCRICH1Hit*>( (*ch) );
		vPoCluDetMC.setX( hit->getX() );
		vPoCluDetMC.setY( hit->getY() );
		vPoCluDetMC.setZ( hit->getZ() );
	        vPosEmMC.setX( hit->getXprod() );
	        vPosEmMC.setY( hit->getYprod() );
	        vPosEmMC.setZ( hit->getZprod() );
	        done = true;
		//cout << vPosEmMC << endl;
	      }
	    }
	  }
	  if( done ) {
            list<CsRCPhotonDet*> lPhoDet = dets->lPhoDet();
            list<CsRCPhotonDet*>::iterator id = lPhoDet.begin();
	    for( int j=0; j<kDetClu; j++ ) id++;
	    Hep3Vector vddr;
            double xw, yw, zw;
	    HepVector vrot( 3, 0 );
	    for( int j=0; j<3; j++ ) vrot[j] =
	      (vPoCluDetMC - mirr->vC0v( kDetClu ))[j];
	    HepVector vPoCluDetMCW = (*id)->rotMatrix() * vrot;
            Hep3Vector vPoCluDetMCWw( vPoCluDetMCW[0], vPoCluDetMCW[1],
				      vPoCluDetMCW[2] );
	    //	    vPoCluDetw.setX( vPoCluDetMCWw.x() );       // true MC hit
	    //	    vPoCluDetw.setY( vPoCluDetMCWw.y() );       // true MC hit
	    //	    vPoCluDetw.setZ( vPoCluDetMCWw.z() );       // true MC hit
	    for( int j=0; j<3; j++ ) vrot[j] =
	      (vPosEmMC - mirr->vC0v( kDetClu ))[j];
	    HepVector vPosEmMCW = (*id)->rotMatrix() * vrot;
            Hep3Vector vPosEmMCWw( vPosEmMCW[0], vPosEmMCW[1], vPosEmMCW[2] );
	    vPoPhotWw.setX( vPosEmMCWw.x() );
	    vPoPhotWw.setY( vPosEmMCWw.y() );
	    vPoPhotWw.setZ( vPosEmMCWw.z() );
	    //cout << vPoPhotWw << endl;
	  }
	  else {
	    continue;                             // true em pt only
	    //	    if( cluHits.size() > 0 ) continue;    // true em pt + bkgr
	  }
//---------------------------------------------------------------------*/

	  Hep3Vector vPoCluDetNC = vPoCluDetw;
	  Hep3Vector vPoPhotWNC = vPoPhotWw;
	  float RRNC = RR;
	  double theCluPNC = 0.;

//------- correction for mirror centre misalignment (photon reflection) :
//        ---------------------------------------------------------------
	  std::string sMirrElePa = " ";
	  if( mirr->pMirrElePa() ) sMirrElePa = mirr->pMirrElePa()->name();
	  //cout << "getPartPhotons " << kDetClu << "  " << part->kDetPart()
	  //     << "  " << sMirrElePa
	  //     << "  " << mirr->corPhoRR()
	  //     << "  " << vCorPoPhoW()[kDetClu] 
	  //<< endl;
	  //cout << RR << "  " << vPoCluDetw << "  " << vPoPhotWw << endl;

          if( kDetClu == part->kDetPart() ) {          //   030120
	    RR += mirr->corPhoRR();
	    vPoCluDetw += vCorPoPhoW()[kDetClu];
	    //vPoPhotWw += vCorPoPhoW()[kDetClu];      //   040227 ???
	  }
	  //cout << RR << "  " << vPoCluDetw << "  " << vPoPhotWw << endl;

//------- detector (only) misalignment to trackers (and mirror) MWR ! :
//        -------------------------------------------------------------
          //vPoCluDetw += ???;                         //   pro-memoria
//------- mirror and detector misalignment to trackers (= particle) MWR ! :
//        -----------------------------------------------------------------
          //vPoPhotWw += ???;                          //   pro-memoria
//------- mirror (only) misalignment to trackers (and detector) MWR ! :
//        -------------------------------------------------------------
//        (both corrections equal)

//------- plot approximate photon hit on mirrors :
//        ----------------------------------------
          //phoToMirr( ... );                          //   pro-memoria     
//        ----------------

//------- particle at entrance (MWR) :
//        ----------------------------
          Hep3Vector vPoPartWw = vPoPartW()[kDetClu];
          Hep3Vector vDcPartWw = vDcPartW()[kDetClu];
          //cout << vPoPartWw << "  " << vPoPhotWw << "  " 
	  //     << vPoCluDetw << endl;

//------- compute emission angle (mrad) :
//        -------------------------------
	  bool aFlag = true;
	  Hep3Vector vDcPhoEmW;
          vDcPhoEmW = getCluAngle( vPoPhotWw, vPoCluDetw, RR, aFlag );
//                    -----------------------------------------------
	  if( !aFlag ) continue;
          double cosThClu = vDcPartWw * vDcPhoEmW;
          if( cosThClu >  1.) cosThClu =  1.;
          if( cosThClu < -1.) cosThClu = -1.;
          double theCluP = acos( cosThClu );
	  theCluP *= 1000.;

//------- use fiducial region around particle direction :
//        -----------------------------------------------
          if( theCluP < 70. ) {               /*   theta < 70. mrad   */
	  //if( cosThClu > 0.99755 ) {        /*   theta < 70. mrad   */
//@@--------------------------------

//--------------------------------------------------------------------------
//--------- monitor mirror and ... corrections :
//          ------------------------------------
	    bool aFlagNC = false;
            if( CsRCHistos::Ref().bookHis() &&
		CsRCHistos::Ref().levelBk() >= 2 ) aFlagNC = true;
	    if( aFlagNC ) {
              Hep3Vector vDcPhoEmWNC;
	      bool flg = true;
	      vDcPhoEmWNC = getCluAngle( vPoPhotWNC, vPoCluDetNC, RRNC, flg );
//            ---------------------------------------------------------------
	      if( flg ) {
  	        //std::cout << vDcPhoEmW << "  " << vDcPhoEmWNC << std::endl;
		double cosThCluNC = vDcPartWw * vDcPhoEmWNC;
		if( cosThCluNC >  1.) cosThCluNC =  1.;
		if( cosThCluNC < -1.) cosThCluNC = -1.;
		theCluPNC = acos( cosThCluNC );
		theCluPNC *= 1000.;
		xh = theCluPNC - theCluP;
		yh = theCluP;
		int kch = 0;
		if( kch >= 0  &&  kch < int( vRC1100.size() ) ) {
		  if( vRC1100[kch] ) vRC1100[kch]->Fill( xh, yh );
//                                   ----------------------------
		}
	      } else {
		aFlagNC = false;
	      }
	    }
//--------------------------------------------------------------------------

//--------- approximate correction for PMT OPTICS :
//          ---------------------------------------
	    bool oFlag = true;
	    Hep3Vector vPoCluDetO = vPoCluDetw;
	    Hep3Vector vDcPhoEmWO = vDcPhoEmW;
	    double theCluPO = theCluP;
	    double theCluPOt = -100.;
	    static int nEnt = 0;
	    static int nRej = 0;
	    if( CorrPMTOpt  &&  (*icl)->isPMT() ) {

	      vDcPhoEmWO = doPMTOptCorr( (*icl), vPoPhotWw, vDcPhoEmWO,
//                         --------------------------------------------
					 vPoCluDetO, RR, aFlag, oFlag );
	      //std::cout << vPoCluDetw << "  " << vPoCluDetO << std::endl;

//----------- discard uncorrected clu.s (TEST!)
	      //if( !oFlag ) continue;

//----------- //^compute 'corrected' emission angle :
//            //^------------------------------------
	      //^aFlag = true;
              //^vDcPhoEmW = getCluAngle( vPoPhotWw, vPoCluDetw, RR, aFlag );
//                        -----------------------------------------------
	      if( !aFlag ) continue;
              cosThClu = vDcPhoEmWO * vDcPartWw;
              if( cosThClu >  1.) cosThClu =  1.;
              if( cosThClu < -1.) cosThClu = -1.;
	      theCluPO = acos( cosThClu )*1000.;
	      if( oFlag ) theCluPOt = theCluPO;
	      nEnt++;
	      if( !oFlag ) nRej++;
	      //if( nEnt%1000 == 0 )
	      //  std::cout << nEnt << "  " << nRej << std::endl;

	      vPoCluDetw = vPoCluDetO;
	      vDcPhoEmW = vDcPhoEmWO;

	      int kch = 1;
	      xh = theCluPO - theCluP;
	      yh = theCluP;
	      if( kch >= 0  &&  kch < int( vRC1100.size() ) ) {
		if( vRC1100[kch] ) vRC1100[kch]->Fill( xh, yh );
//                                 ----------------------------
	      }

	    }

//--------- approximate correction for detector QUARTZ WINDOW :
//          ---------------------------------------------------
	    bool qFlag = true;
	    Hep3Vector vPoCluDetQ(0., 0., 0.);
	    Hep3Vector vDcPhoEmWQ(0., 0., 0.);
	    double theCluPQ = theCluP;

	    Hep3Vector vPoCluDetQx = vPoCluDetw;
	    Hep3Vector vDcPhoEmWQx = vDcPhoEmW;
	    double theCluPQx = theCluP;
            if( CorrQzW ) {

              vPoCluDetQ = doQzWCorr( vPoPhotWw, vDcPhoEmW, vPoCluDetw,
//                         --------------------------------------------
				      RR, zDetW );

//----------- compute 'corrected' emission angle :
//            ------------------------------------
	      aFlag = true;
              vDcPhoEmWQ = getCluAngle( vPoPhotWw, vPoCluDetQ, RR, aFlag );
//                         -----------------------------------------------
	      if( !aFlag ) continue;
	      //std::cout << setprecision(5) << std::endl;
	      //std::cout << " 1 " << vDcPhoEmWQ << "  "
	      //          << vPoCluDetQ << std::endl;
              cosThClu = vDcPhoEmWQ * vDcPartWw;
              if( cosThClu >  1.) cosThClu =  1.;
              if( cosThClu < -1.) cosThClu = -1.;
	      theCluPQ = acos( cosThClu )*1000.;

	      vPoCluDetw = vPoCluDetQ;
	      vDcPhoEmW = vDcPhoEmWQ;

//----------- Use instead : on test
              //vDcPhoEmWQ = doQzWCorr( vPoPhotWw, vDcPhoEmW, vPoCluDetQ,
//                         --------------------------------------------
	      //   		         RR, zDetW, aFlag );
	      //if( !aFlag ) continue;

//----------- Under TEST : DO NOT USE!
/*
              vDcPhoEmWQx = doQzWCorr( (*icl), vPoPhotWw, vDcPhoEmWQx,
//            --------------------------------------------------------
	        	    	       vPoCluDetQx, RR, zDetW, qFlag );
	      if( !qFlag ) continue;
              cosThClu = vDcPhoEmWQx * vDcPartWw;
              if( cosThClu >  1.) cosThClu =  1.;
              if( cosThClu < -1.) cosThClu = -1.;
	      theCluPQx = acos( cosThClu )*1000.;

	      vPoCluDetw = vPoCluDetQx;
	      vDcPhoEmW = vDcPhoEmWQx;
*/

	      int kch = 2;
	      xh = theCluPQ - theCluPO;
	      yh = theCluP;
	      if( kch >= 0  &&  kch < int( vRC1100.size() ) ) {
		if( vRC1100[kch] ) vRC1100[kch]->Fill( xh, yh );
//                                 ----------------------------
	      }
	      kch = 3;
	      xh = theCluPQx - theCluPO;
	      yh = theCluP;
	      if( kch >= 0  &&  kch < int( vRC1100.size() ) ) {
		if( vRC1100[kch] ) vRC1100[kch]->Fill( xh, yh );
//                                 ----------------------------
	      }
	      xh = theCluPQ - theCluPO;
	      yh = theCluPQx - theCluPO;
	      if( hRC1110 ) hRC1110->Fill( xh, yh );
//                          -----------------------

	    }


//--------- photon theta :
//          --------------
            double theClu, phiClu;
            theClu = acos( cosThClu );

//--------- photon phi :
//          ------------
            double alpa = vDcPartWw.x();
            double ampa = vDcPartWw.y();
            double anpa = vDcPartWw.z();
            double anorm = sqrt( (1.-ampa*ampa)*(1.-cosThClu*cosThClu) );
            if( anorm > 0.) {
              double cosPhClu = (vDcPhoEmW.y() - ampa*cosThClu) / anorm;
              if( cosPhClu >  1.) cosPhClu =  1.;
              if( cosPhClu < -1.) cosPhClu = -1.;
              phiClu = acos( cosPhClu );
              if( (vDcPhoEmW.z()*alpa - vDcPhoEmW.x()*anpa) < 0. )
                phiClu = TwoPI - phiClu;
            }  else  phiClu = 0.;
	    double phiCluS = phiClu;
	    phiCluS *= RadDeg;

//--------- photon phi in MWR (to vertical plane) :
//          ---------------------------------------
//          (for mirror alignments)
//          -----------------------
            double theCluA;
            double phiCluA;
	    Hep3Vector vDcPhoEmP = CsRCUty::Instance()->
              rotfbcc( +1., vDcPartWw, vDcPhoEmW, theCluA, phiCluA );
//            ------------------------------------------------------
	    //cout << "theClu, phiClu = " << theCluA
	    //     << "  " << phiCluA << endl;
            theCluA = theCluA * 1000. ;          /*   --> mrad  */
            phiCluA = phiCluA * RadDeg;          /*   --> deg   */

//--------- photon phi with respect to normal of 'particle plane' :
//          -------------------------------------------------------
//          (for sigma photon analysis, 22/11/99)
//          -------------------------------------
            Hep3Vector vPoPaMirWw = vPoPaMirW()[kDetClu];
	    //cout << "------" << endl;
	    //cout << kDetClu << "  " << vPoPaMirWw << "   "
	    //     << vDcPartWw << endl;
	    //vPoPaMirWw.setX( 0. );
	    //vDcPartWw.setX( 0. );
	    //vDcPartWw = vDcPartWw.unit();
	    //cout << vPoPaMirWw << "   " << vDcPartWw << endl;
	    //cout << vDcPartWw.cross( vPoPaMirWw ) << endl;
            double cosPaMir = vDcPartWw * vPoPaMirWw / RR;
            if( cosPaMir >  1.) cosPaMir =  1.;
            if( cosPaMir < -1.) cosPaMir = -1.;
            double thePaMir = acos( cosPaMir );
	    part->setThPamir( thePaMir );                        // ???
            Hep3Vector vcrossPhoPaw = vDcPhoEmW.cross( vDcPartWw );
            double cosPhi = vcrossPhoPaw * vDcPartWw.cross( vPoPaMirWw );
            double norm = (1.-cosThClu*cosThClu) * (1.-cosPaMir*cosPaMir);
            //norm = sqrt( norm ) * RR;
            if( norm > 0. ) {
	      norm = sqrt( norm ) * RR;
              cosPhi /= norm;
	      //cout << acos( cosPhi )*RadDeg << endl;
              if( cosPhi >  1.) cosPhi =  1.;
              if( cosPhi < -1.) cosPhi = -1.;
              phiClu = acos( cosPhi );
              if( vcrossPhoPaw.y() < 0. ) phiClu = TwoPI - phiClu;
            }  else  {
              phiClu = 0.;
	    }
            double theCluR = theClu;
            double phiCluR = phiClu;
            theCluR = theCluR * 1000. ;          /*   --> mrad  */
            phiCluR = phiCluR * RadDeg;          /*   --> deg   */
	    //cout << theCluR << "   " << phiCluR << endl;
	    //std::cout << theCluPNC << "  " << theCluP << "  " << theCluR
	    //          << std::endl;

	    xh = theCluR - theCluP;
	    yh = phiCluR;
	    if( hist.hRC3504 ) hist.hRC3504->Fill( xh, yh );
//hh                           ----------------------------
	    xh = phiCluR;
	    yh = phiCluA;
	    if( hist.hRC3502 ) hist.hRC3502->Fill( xh, yh );
//hh                           ----------------------------
	    xh = theCluPNC - theCluR;
	    yh = phiCluR;
	    if( aFlagNC ) if( hist.hRC3510 ) hist.hRC3510->Fill( xh, yh );
//hh                                         ----------------------------

/*--------------------------------------------------------------------------
	    Hep3Vector ss = vPoPaMirWw - vPoPartWw;
	    Hep3Vector ddss = ss * 0.1;
	    Hep3Vector vPoPartWX = vPoPartWw;
	    double theCluX[20];
	    for( int k=0; k<=10; k++ ) {
	      bool aFlagX = true;
              Hep3Vector vDcPhoEmWX = getCluAngle( vPoPartWX, vPoCluDetw, RR,
		  				   aFlagX );
	      if( aFlagX ) {
                double cosThClu = vDcPartWw * vDcPhoEmWX;
                if( cosThClu >  1.) cosThClu =  1.;
                if( cosThClu < -1.) cosThClu = -1.;
                theCluX[k] = acos( cosThClu );
	        theCluX[k] *= 1000.;
	      }
	      vPoPartWX += ddss;
	    }
	    if( theCluR < theCluX[0] && theCluR < theCluX[10] ||
		theCluR > theCluX[0] && theCluR > theCluX[10] ) {
	      for( int j=0; j<=10; j++ ) cout << theCluX[j] << "  ";
	      cout << "  " << theCluR << endl;
	    }
//------------------------------------------------------------------------*/

//--------- OBSOLETE ------
	    Hep3Vector vDcPhoEmWB, vDcPhoEmWM;
	    bool aFlagB = false;
	    //bool aFlagB = true;
            //vDcPhoEmWB = getCluAngle( vPoPartWw, vPoCluDetw, RR, aFlagB );
//                       ------------------------------------------------
	    bool aFlagM = false;
	    //bool aFlagM = true;
            //vDcPhoEmWM = getCluAngle( vPoPaMirWw, vPoCluDetw, RR, aFlagM );
//                       -------------------------------------------------
	    double theCluB = 1000000.;
	    double theCluM = 1000000.;
	    double phiCluB, phiCluM;
	    if( aFlagB && aFlagM ) {
              double cosThClu = vDcPartWw * vDcPhoEmWB;
              if( cosThClu >  1.) cosThClu =  1.;
              if( cosThClu < -1.) cosThClu = -1.;
              theCluB = acos( cosThClu );
	      theCluB *= 1000.;
              //double anorm = sqrt( (1.-ampa*ampa)*(1.-cosThClu*cosThClu) );
              //if( anorm > 0.) {
              //  double cosPhClu = (vDcPhoEmWB.y() - ampa*cosThClu) / anorm;
              //  if( cosPhClu >  1.) cosPhClu =  1.;
              //  if( cosPhClu < -1.) cosPhClu = -1.;
              //  phiCluB = acos( cosPhClu );
              //  if( (vDcPhoEmWB.z()*alpa - vDcPhoEmWB.x()*anpa) < 0. )
              //    phiCluB = TwoPI - phiCluB;
              //}  else  phiCluB = 0.;
	      //phiCluB *= RadDeg;
              cosThClu = vDcPartWw * vDcPhoEmWM;
              if( cosThClu >  1.) cosThClu =  1.;
              if( cosThClu < -1.) cosThClu = -1.;
              theCluM = acos( cosThClu );
	      theCluM *= 1000.;
              //anorm = sqrt( (1.-ampa*ampa)*(1.-cosThClu*cosThClu) );
              //if( anorm > 0.) {
              //  double cosPhClu = (vDcPhoEmWM.y() - ampa*cosThClu) / anorm;
              //  if( cosPhClu >  1.) cosPhClu =  1.;
              //  if( cosPhClu < -1.) cosPhClu = -1.;
              //  phiCluM = acos( cosPhClu );
              //  if( (vDcPhoEmWM.z()*alpa - vDcPhoEmWM.x()*anpa) < 0. )
              //    phiCluM = TwoPI - phiCluM;
              //}  else  phiCluM = 0.;
	      //phiCluM *= RadDeg;
	      //cout << theCluB << "  " << theCluR << "  "
	      //   << theCluM << "     "
	      //   << phiCluB << "  " << phiCluS << "  "
	      //   << phiCluM << endl;
//----------- theta entr - theta mirr vs phi photon :
//            -------------------------------------
	      //double ddBM = fabs(theCluB-theCluR) - fabs(theCluR-theCluM);
              //double thetaBM = theCluB - theCluM;
	      //cout << ddBM << "  " << fabs( ddBM / thetaBM ) << endl;
	    }
//------------------------------------------------------------------------

//--------- 080918 ---
//          theCluB, theCluM currently NOT USED
//          provisionally used for Opt. Corr. tests :
	    theCluB = theCluP;
	    theCluM = theCluPOt;
//--------- 090918 ------------

//--------- store useful clusters as 'photons' :
//          ------------------------------------
	    kPhot++;
            double PHw = (*icl)->PH();
	    list<CsRCPhotonDet*> lPhoDet = dets->lPhoDet();
	    list<CsRCPhotonDet*>::iterator id;
	    id = lPhoDet.end();
	    if( kDetClu == 0 ) { id = lPhoDet.begin(); }
	    if( kDetClu == 1 ) { id = lPhoDet.begin(); id++; }
            if( id != lPhoDet.end() ) lPhotons_.push_back( new CsRCPhoton(
                                         kPhot, theCluR, phiCluR,
				         phiCluA, theCluB, theCluM, PHw,
				         (*icl), (*id), this, kDetClu ) );

//cout << kDetClu << "  " << (*icl)->xc() << "  " << (*icl)->yc() << endl;

/*
// tests 01/03/19 -----------------------------------------------------------
	   Hep3Vector vPoCluDeB = photImpDet( vPoPartWw, vDcPhoEmW, zw, RR );
	   Hep3Vector vPoCluDeC = photImpDet( vPoPhotWw, vDcPhoEmW, zw, RR );
	   Hep3Vector vPoCluDeM = photImpDet( vPoPaMirWw, vDcPhoEmW, zw, RR );
//---------------------------------------------------------------------------
*/
          }   /* end if on cosThClu */

        }   /* end if on bProClu */

      }   /* end loop on clusters: kClu */

      //--- protect for zero photons
      if( lPhotons_.empty() ) flag_ = false;
      if( !flag_ ) {
	if( key->kPrintRejections() == 1 ) {
          std::cout << "RICHONE - getPartPhotons : Ev "
                    << CsRichOne::Instance()->kEvent() << "  PaPhoton "
		    << kPaPhot_ << "  mom "
		    << pPart_->mom() << " NO photons" <<  std::endl;
	}
      }


//--- monitoring histograms :
//    -----------------------
      int nPhotonPaPho = lPhotons_.size();
      xh = nPhotonPaPho;
      if( hist.hRC1131 ) hist.hRC1131->Fill( xh );
//hh                     ------------------------

  }



//===========================================================================
  Hep3Vector CsRCPartPhotons::getCluAngle( const Hep3Vector &vPoPho,
//--------------------------------------------------------------------------
                                           const Hep3Vector &vPoClu,
					   const float &RR, bool &flag ) {

//--- Paolo  -  June 1999. 
//--- from SegYpsC and segyps_.
//--- rev.  October  2000


//--- on 'photon' 'impact plane' (YS-plane) in the MWR :
//    --------------------------------------------------

      double rrPho = vPoPho.mag();
      double rrClu = vPoClu.mag();

//--- cos(omega-d) :
//    --------------
      double cosOmegaDe = (vPoPho * vPoClu) / (rrPho * rrClu);
      double sinw = 1.- cosOmegaDe*cosOmegaDe;
      if( sinw < 0. ) {
        string str = 
	  "RICHONE, CsRCPartPhotons::getCluAngle() : cosOmegaDe error";
        CsErrLog::Instance()->mes( elError, str );
      sinw = 0.;
      flag = false;
      return Hep3Vector(0., 0., 0.);
      }
      double sinOmegaDe = sqrt( sinw );

//--- cos(omega-e) (iterative procedure) :
//    ------------------------------------
      double sinOw = sinOmegaDe;
      //int nIter = 3;
      int nIter = 5;                        //   030721
//@@----------------
      double raPhoClu = rrPho / rrClu;
      double raPhoRR  = rrPho / RR;

      //std::cout << " " << sinOw << std::endl;
      //sinOw = 0.;
      //std::cout << " " << sinOw << std::endl;
      bool perr = true;
      for( int kk=0; kk < nIter; kk++ ) {
        double aswDe = raPhoClu * sinOw;
        if( aswDe > 1. && perr ) {
          string str = 
	    "RICHONE, CsRCPartPhotons::getCluAngle() : aswDe error";
          CsErrLog::Instance()->mes( elWarning, str );
          //            cout << "getCluAngle : aswDe error" << endl;
          aswDe = 1.;
          perr = false;
	  flag = false;
	  return Hep3Vector(0., 0., 0.);
        }
        double aswRR = raPhoRR * sinOw;
        if( aswRR > 1. && perr ) {
          string str =
	    "RICHONE, CsRCPartPhotons::getCluAngle() : aswRR error";
          CsErrLog::Instance()->mes( elWarning, str );
	  //            cout << "getCluAngle : aswRR error" << endl;
          aswRR = 1.;
          perr = false;
	  flag = false;
	  return Hep3Vector(0., 0., 0.);
        }
        double aswRR2 = aswRR*aswRR;
        double cosRR = sqrt(1.- aswRR2);
        double cosDe = sqrt(1.- aswDe*aswDe);
        double co2RR2 = 1.- 2.*aswRR2;
        double sinCorr = co2RR2*aswDe - 2.*aswRR*cosRR*cosDe;
        double cosCorr = co2RR2*cosDe + 2.*aswRR*cosRR*aswDe;
        sinOw = sinOmegaDe * cosCorr - cosOmegaDe * sinCorr;
	//std::cout << kk+1 << "  " << sinOw << std::endl;
      }

      double sinOmegaEm = sinOw;
      double cosw = 1.- sinOmegaEm*sinOmegaEm;
      if( cosw < 0. ) {
        string str =
	  "RICHONE, CsRCPartPhotons::getCluAngle() : cosOmegaEm error";
        CsErrLog::Instance()->mes( elError, str );
	//        cout << "getCluAngle : cosOmegaEm error" << endl;
        cosw = 0.;
	flag = false;
	return Hep3Vector(0., 0., 0.);
      }
      double cosOmegaEm = sqrt( cosw );

//--- photon polar angles :
//    ---------------------
      sinw = sinOmegaEm;
      double sinww = (sinOmegaEm*cosOmegaDe - sinOmegaDe*cosOmegaEm);
      double aa = sinw / rrClu;
      double bb = sinww / rrPho;
      Hep3Vector vvw = aa * vPoClu - bb * vPoPho;

      return  vvw.unit();

  }



//======================================================================
  Hep3Vector CsRCPartPhotons::doQzWCorr( 
//----------------------------------------------------------------------
                        Hep3Vector vPoPhotW, Hep3Vector vDcPhoEmW,
                        Hep3Vector vPoCluDet, float RR, float zDetW )  {


//--- approximate correction for detector quartz window :
//    ---------------------------------------------------

      CsRCExeKeys *key = CsRCExeKeys::Instance();


      static bool firstCall = true;
      if( firstCall ) {
        firstCall = false;

        key->acknoMethod( "CsRCPartPhotons::QuartzWindowCorr" );
      }

//--- from "ReconConst.h"
      CsRCRecConst *cons = CsRCRecConst::Instance();
      float CFRefInd = cons->CFRefInd();
      float qzRefInd = cons->qzRefInd();

      CsRCDetectors *dets = CsRCDetectors::Instance();
      list<CsRCCathode*> lCathodes = dets->lCathodes();
      CsRCCathode *cat  = lCathodes.front();
      static float ddQzW = cat->ddQzW();
      static float ddGap = cat->ddGap();

      static double zQzWOutside = zDetW + ddGap + ddQzW;
      static double zQzWinside  = zDetW + ddGap;

      static double rRatio = CFRefInd / qzRefInd;

//--- assume photon emitted from average point on part. traj.
//    -------------------------------------------------------
//--- 'photon' impact on mirror (MWR) :
//    ---------------------------------
      CsRCMirrors *mirr = CsRCMirrors::Instance();
      Hep3Vector vPoC( 0., 0., 0. );
      Hep3Vector vPoPhoMir = mirr->vImpMir( vPoPhotW, vDcPhoEmW, vPoC, RR );
//                           ----------------------------------------------

//--- normal to mirror at 'photon' impact :
//    -------------------------------------
      Hep3Vector vDcNoPhoMir = (1./RR) * vPoPhoMir;

//--- 'photon' reflected direction :
//    ------------------------------
      double cosPhoMir = vDcNoPhoMir * vDcPhoEmW;
      Hep3Vector vDcPhoRefl = 2.*cosPhoMir * vDcNoPhoMir - vDcPhoEmW;
//pp      cout << kMir << "  vDcPhoRefl  " << vDcPhoRefl << endl;

//--- 'photon' incidence angle on detector :
//    --------------------------------------
      double cosPhoDet = vDcPhoRefl.z();
      double sinPhoDet = sqrt(1.-cosPhoDet*cosPhoDet);

      double norm;
//--- 'photon' impact on quartz win. (outside) :
//    ------------------------------------------
      norm = (zQzWOutside - vPoPhoMir.z()) / vDcPhoRefl.z();
      Hep3Vector vPhoWinOut = vPoPhoMir + norm * vDcPhoRefl;

//--- 'photon' projection on quartz win. (inside) :
//    ---------------------------------------------
      norm = (zQzWinside - vPoPhoMir.z()) / vDcPhoRefl.z();
      Hep3Vector vPhoWinIn = vPoPhoMir + norm * vDcPhoRefl;

//--- 'photon' direction inside quartz win. (DRS) :
//    ---------------------------------------------
      double sinPhoQzR = sinPhoDet * rRatio;
      double cosPhoQzR = sqrt(1.-sinPhoQzR*sinPhoQzR);
      Hep3Vector vDcPhoQzR;
      vDcPhoQzR.setX( rRatio * vDcPhoRefl.x() );
      vDcPhoQzR.setY( rRatio * vDcPhoRefl.y() );
      vDcPhoQzR.setZ( cosPhoQzR );

//--- refracted 'photon' impact on quartz win. (inside) :
//    ---------------------------------------------------
      norm = - ddQzW / vDcPhoQzR.z();
      Hep3Vector vPhoWinQzR = vPhoWinOut + norm * vDcPhoQzR;

//--- corrections to x and y :
//    ------------------------
      Hep3Vector vcorr = vPhoWinQzR - vPhoWinIn;
      vcorr.setZ( 0.);
//pp      cout << "vcorr = " << vcorr << endl;

      Hep3Vector vout = vPoCluDet - vcorr;

      return vout;
  }


//===========================================================================
  void CsRCPartPhotons::doPartToMWR( CsRCParticle* part ) {
//---------------------------------------------------------


//--- from "CsRCRecConst"
      CsRCRecConst *cons = CsRCRecConst::Instance();
      float partPathFr = cons->partPathFr();

      CsRCDetectors *dets = CsRCDetectors::Instance();
      double zEntrWind = dets->zEntrWind();
      double zExitWind = dets->zExitWind();
      CsRCMirrors *mirr = CsRCMirrors::Instance();

      CsRCHistos& hist = CsRCHistos::Ref();
      double xh, yh, wh;

//--- particle variables :
//    --------------------
      Hep3Vector vPoPart0 = part->vPosIn();
      Hep3Vector vDcPart0 = part->vDirIn();

      vector<double> thPadew;
      vector<double> ddPaDet;

//--- compute variables in both MWR :
//    -------------------------------
//--- LOOP over up and down MIRROR :
//    ------------------------------
      list<CsRCMirrorNom*> lMirrNom = mirr->lMirrNom();
      list<CsRCMirrorNom*>::iterator in = lMirrNom.begin();
      list<CsRCPhotonDet*> lPhoDet = dets->lPhoDet();
      list<CsRCPhotonDet*>::iterator id = lPhoDet.begin();
      int kMir = 0;
      for( in=lMirrNom.begin(); in!=lMirrNom.end(); in++ ) {

        Hep3Vector vPoC0( (*in)->vC0() );
        double RR = (*in)->RR();
        //cout << "vPoPart0   " << vPoPart0 << endl;
        //cout << "vDcPart0   " << vDcPart0 << endl;
        //cout << "vPoC0   " << vPoC0 << endl;
	//cout << "doPartToMWR " << RR << "  " << vPoC0 << endl;

//----- particle impact on mirrors (MWR) :
//      ----------------------------------
	Hep3Vector vPoPaMir0 = part->vPoPaMir0()[kMir];
	//cout << endl << kMir << "  vPoPaMir0  " << vPoPaMir0 << endl;


//----- correction for mirror centre misalignment (particle 'reflection') :
//      -------------------------------------------------------------------
//      (i.e. mirror ELEMENT hit by the particle)
//      -----------------------------------------
	if( kMir == part->kDetPart() ) {           //   0301023
	  vPoC0 += mirr->vCorPoPa0();
//@@--------------------------------
	  RR += mirr->corPaRR();
//@@---------------------------
	}
//----- WARNING : vPoPaMir0 should be ricalculated and set again (201220)
        ////vPoPaMir0 = mirr->vImpMir( vPoPart0, vDcPart0, vPoC0, RR );
        ////            ----------------------------------------------
	////part->setPoPaMir0( ... );
	//cout << "doPartToMWR - corr " << RR << "  " << vPoC0 << endl;

	Hep3Vector vDcDetvw = dets->vDcDetv( kMir );

//----- normal to mirror at particle impact :
//      -------------------------------------
        Hep3Vector vDcNoPaMir0 = (1./RR) * (vPoPaMir0 - vPoC0);

//----- particle 'reflected' direction :
//      --------------------------------
        double cosPaMir = vDcNoPaMir0 * vDcPart0;
        Hep3Vector vDcPaRefl0 = 2.*cosPaMir * vDcNoPaMir0 - vDcPart0;
	//cout << kMir << "  vDcPaRefl0  " << vDcPaRefl0 << endl;

//----- particle 'impact' on detectors :
//      --------------------------------
        double norm = ( ( dets->vDet0v( kMir ) - vPoPaMir0 ) * vDcDetvw ) /
                      ( vDcPaRefl0 * vDcDetvw );
        Hep3Vector vPoPaDet0 = vPoPaMir0 + norm * vDcPaRefl0;
	//cout << "vPoPaDet0   " << vPoPaDet0 << endl;

//----- monitor of geometry :
//      ---------------------
	if( kMir == part->kDetPart() ) {
          xh = vPoPaMir0.z();
          yh = fabs( vPoPaMir0.y() );
          if( hist.hRC3503 ) hist.hRC3503->Fill( xh, yh );
//                           ----------------------------
          xh = vPoPaDet0.z();
          yh = fabs( vPoPaDet0.y() );
          if( hist.hRC3503 ) hist.hRC3503->Fill( xh, yh );
//                           ----------------------------
	}

//----- photon 'emission' point (assume particle path 'middle' point) :
//      ---------------------------------------------------------------
	Hep3Vector vPoPhot0 = vPoPart0 + partPathFr * (vPoPaMir0 - vPoPart0);
	//cout << "vPoPhot0   " << vPoPhot0 << endl;

//----- particle DIRECTION at photon 'emission' point (extrapolated) :
//      --------------------------------------------------------------
//      corrected to 'medium' point <entr.-wind.--impact-on-mirror>
	Hep3Vector vDcPartEm0 = part->vDirEmP();
	double dZem = vPoPhot0.z() - part->vPosEmP().z();
	double dXem = vDcPartEm0.x()/vDcPartEm0.z() * dZem;
	double dYem = vDcPartEm0.y()/vDcPartEm0.z() * dZem;
	if( part->vPosEmP().z() < 1000000. ) {
	  vPoPhot0.setX( part->vPosEmP().x() + dXem );
	  vPoPhot0.setY( part->vPosEmP().y() + dYem );
	}
	//cout << "      " << vPoPhot0 << endl;

//----- CORRECTIONS for average mcs through RICH :
//      ------------------------------------------
	bool PartCorr = CsRCExeKeys::Instance()->PartCorr();
	if( PartCorr ) {
	  Hep3Vector vDcPartCorrPo(0., 0., 0.);
	  Hep3Vector vDcPartCorrDc(0., 0., 0.);
	  if( mcsPartCorr( part, vDcPartCorrPo, vDcPartCorrDc ) ) {
//@@-------------------------------------------------------------
	    vPoPhot0 += vDcPartCorrPo;
	    double tana = vDcPartEm0.x()/vDcPartEm0.z() + vDcPartCorrDc.x();
	    double tanb = vDcPartEm0.y()/vDcPartEm0.z() + vDcPartCorrDc.y();
	    vDcPartEm0.setX( tana );
	    vDcPartEm0.setY( tanb );
	    vDcPartEm0.setZ( 1. );
	    vDcPartEm0 = vDcPartEm0.unit();
	  }
	}


//----- TRANSFORM to up or down MWR system :
//      ------------------------------------
        double xw, yw, zw;
	HepVector vrot( 3, 0 );

//----- ... particle position at RICH entrance :
//      ----------------------------------------
	for( int j=0; j<3; j++ ) vrot[j] = (part->vPosIn() - vPoC0)[j];
	HepVector vPoPartW = (*id)->rotMatrix() * vrot;
        Hep3Vector vPoPartWw( vPoPartW[0], vPoPartW[1], vPoPartW[2] );
        vPoPartW_.push_back( vPoPartWw );
	//cout << "vPoPartW   " << vPoPartWw << endl;

//----- ... photon 'emission' point :
//      -----------------------------
	for( int j=0; j<3; j++ ) vrot[j] = (vPoPhot0 - vPoC0)[j];
	HepVector vPoPhotW = (*id)->rotMatrix() * vrot;
	Hep3Vector vPoPhotWw( vPoPhotW[0], vPoPhotW[1], vPoPhotW[2] );
        vPoPhotW_.push_back( vPoPhotWw );
	//cout << "vPoPhotWw   " << vPoPhotWw << endl;

//----- ... particle direction at photon 'emission' point :
//      ---------------------------------------------------
	for( int j=0; j<3; j++ ) vrot[j] = vDcPartEm0[j];
	HepVector vDcPartW = (*id)->rotMatrix() * vrot;
        Hep3Vector vDcPartWw( vDcPartW[0], vDcPartW[1], vDcPartW[2] );
        vDcPartW_.push_back( vDcPartWw );
	//cout << "vDcPartW   " << vDcPartWw << endl;

//----- ... particle impact on mirror :
//      -------------------------------
	for( int j=0; j<3; j++ ) vrot[j] = (vPoPaMir0 - vPoC0)[j];
	HepVector vPoPaMirW = (*id)->rotMatrix() * vrot;
        Hep3Vector vPoPaMirWw( vPoPaMirW[0], vPoPaMirW[1], vPoPaMirW[2] );
        vPoPaMirW_.push_back( vPoPaMirWw );

//----- ... correction to mirror centre ( photon reflection ) :
//      -------------------------------------------------------
	for( int j=0; j<3; j++ ) vrot[j] = mirr->vCorPoPho0()[j];
	HepVector vCorPoPhoW = (*id)->rotMatrix() * vrot;
        Hep3Vector vCorPoPhoWw( vCorPoPhoW[0], vCorPoPhoW[1], vCorPoPhoW[2] );
        vCorPoPhoW_.push_back( vCorPoPhoWw );
	//cout << "vCorPoPhoW   " << vCorPoPhoWw << endl;

//----- ... 'reflected' particle direction :
//      ------------------------------------
	for( int j=0; j<3; j++ ) vrot[j] = vDcPaRefl0[j];
	HepVector vDcPaReflW = (*id)->rotMatrix() * vrot;
        Hep3Vector vDcPaReflWw( vDcPaReflW[0],  vDcPaReflW[1], vDcPaReflW[2] );
        vDcPaReflW_.push_back( vDcPaReflWw );
	//cout << "vDcPaReflW   " << vDcPaReflWw << endl;

//----- angle particle detector :
//      -------------------------
	thPadew.push_back( acos( vDcPaReflWw.z() ) );

//----- ... particle 'impact' on detector :
//      -----------------------------------
	for( int j=0; j<3; j++ ) vrot[j] = (vPoPaDet0 - vPoC0)[j];
	HepVector vPoPaDetW = (*id)->rotMatrix() * vrot;
        Hep3Vector vPoPaDetWw( vPoPaDetW[0], vPoPaDetW[1], vPoPaDetW[2] );
        vPoPaDetW_.push_back( vPoPaDetWw );
	//cout << "vPoPaDetW   " << vPoPaDetWw << endl;

//----- distance particle 'impact' to 'beam spot'on detectors :   021212
//      -------------------------------------------------------
	double ddPaDetw = 0.;
	double xPade0 =   0.;
	double yPade0 = 470.;
//@@------------------------
	double xPade = vPoPaDetWw.x();
	double yPade = vPoPaDetWw.y();
        if( yPade >= 0. )
	  ddPaDetw = sqrt( pow( xPade-xPade0, 2 ) + pow( yPade-yPade0, 2 ) );
        if( yPade < 0. )
	  ddPaDetw = sqrt( pow( xPade-xPade0, 2 ) + pow( yPade+yPade0, 2 ) );
	ddPaDet.push_back( ddPaDetw );

        kMir++;
	id++;
      }   /* end of loop on mirrors: kMir */

      part->setThPade( thPadew );
//    --------------------------
      part->setvPade( vPoPaDetW_ );
//    ----------------------------
      part->setDdPaDet( ddPaDet );
//    ---------------------------
      if( !(ddPaDet.size() == 2) )
	cout << "RICHONE, doPartToMWR() : "
             << "wrong variable computation ! " << ddPaDet.size() << endl;

  }


//==========================================================================
  void CsRCPartPhotons::doSelCathodes( CsRCParticle* part ) {
//-----------------------------------------------------------


//--- Paolo - October 2000


      CsRCDetectors *dets = CsRCDetectors::Instance();
      CsRCMirrors *mirr = CsRCMirrors::Instance();
      CsRCHistos& hist = CsRCHistos::Ref();
      double xh, yh, wh;

      static int nMirror;
      static int nCathode;
      static int nHlfCat;

      static float ddCluUsex, ddCluUsey;

      int k;

      static bool firstCall = true;
      if( firstCall ) {
	firstCall = false;

        nMirror = mirr->nMirror();
        nCathode = dets->nCathode();
        nHlfCat = nCathode / 2;

        list<CsRCCathode*> lCathodes = dets->lCathodes();
        CsRCCathode *cat = lCathodes.front();

        CsRCRecConst *cons = CsRCRecConst::Instance();
        ddCluUsex = cons->ddCluUsex();
        ddCluUsey = cons->ddCluUsey();

      }

      int kDetPart = part->kDetPart();
//@@---------------------------------
      lSelCats_.clear();

//--- determine particle path length inside the radiator :
//    ----------------------------------------------------
      double pathLen = ( vPoPaMirW_[kDetPart] - vPoPartW_[kDetPart] ).mag();
      part->setPathLen( pathLen );
//    ---------------------------
      //cout << pathLen << endl;
      xh = pathLen;
      yh = part->mom();
      if( hist.hRC3501 ) hist.hRC3501->Fill( xh, yh );
//                       ----------------------------

      Hep3Vector vPoPaDetW = vPoPaDetW_[kDetPart];

//--- determine the cathode 'hit' by 'reflected' particle :
//    -----------------------------------------------------
//    modif. 01/10/23
//--- (iCaPa = -1 : 'hit' out of cathodes).
//--- note : iCaPa from 0 to 15,  iXpPa,iYpPa from 0 to 71.
//    -----------------------------------------------------
      int iCaPa = -1;
      int iXpPa = 1000000;
      int iYpPa = 1000000;
      for( int kc=0; kc<nCathode; kc++ ) {
	CsRCCathode* cat = dets->ptrToCat( kc );
	double padx = cat->padx();
	double pady = cat->pady();
	double hCatx = cat->hCatx();
	double hCaty = cat->hCaty();
        Hep3Vector vOffCatWw = dets->vOffCatW( kc );
        if( vPoPaDetW.x() >= vOffCatWw.x()-hCatx  &&
            vPoPaDetW.x() <= vOffCatWw.x()+hCatx  &&
            vPoPaDetW.y() >= vOffCatWw.y()-hCaty  &&
            vPoPaDetW.y() <= vOffCatWw.y()+hCaty ) {
          iCaPa = kc;
	  if( iCaPa >= 0  &&  iCaPa < nCathode ) {
            iXpPa = int( (vPoPaDetW.x() - dets->vOffCatW( iCaPa ).x()
      		    + hCatx - padx/2.) / padx );
            iYpPa = int( (vPoPaDetW.y() - dets->vOffCatW( iCaPa ).y()
      		    + hCaty - pady/2.) / pady );
          }
          break;
        }
      }
      //cout << vPoPaDetW << "  ";
      //cout << iCaPa << "  " << iXpPa << "  " << iYpPa << endl;


//--- sort 'useful' cathodes :
//    ------------------------
      float ddCluLix = ddCluUsex;
      float ddCluLiy = ddCluUsey;

/*----------------------------------------------------------------------
//--- modif. 01/10/23  (version 3)
//    ----------------------------
      float xPaLim[4], yPaLim[4];
      float yyOffP = yDetLLimit( kDetPart );
      xPaLim[0] = vPoPaDetW.x() - ddCluLix;
      xPaLim[1] = vPoPaDetW.x() + ddCluLix;
      xPaLim[2] = vPoPaDetW.x() + ddCluLix;
      xPaLim[3] = vPoPaDetW.x() - ddCluLix;
      yPaLim[0] = vPoPaDetW.y() + ddCluLiy - yyOffP;
      yPaLim[1] = vPoPaDetW.y() + ddCluLiy - yyOffP;
      yPaLim[2] = vPoPaDetW.y() - ddCluLiy - yyOffP;
      yPaLim[3] = vPoPaDetW.y() - ddCluLiy - yyOffP;
      for( int kc=0; kc<nCathode; kc++ ) {
	float yyOffC = yDetLLimit( dets->cathodePos( kc ) );
        Hep3Vector vOffCatWw = dets->vOffCatW( kc );
        for( int kp=0; kp < 4; kp++ ) {
          if( xPaLim[kp] >= vOffCatWw.x()-hCatx  &&
              xPaLim[kp] <= vOffCatWw.x()+hCatx  &&
              yPaLim[kp] >= vOffCatWw.y()-hCaty-yyOffC  &&
              yPaLim[kp] <= vOffCatWw.y()+hCaty-yyOffC ) {
	    lSelCats_.push_back( kc );
	    break;
	  }
        }
      }
      //list<int>::iterator ii;
      //for( ii=lSelCats_.begin(); ii!=lSelCats_.end(); ii++ )
      //  cout << (*ii) << " ";
      //cout << endl;
//--------------------------------------------------------------------*/

//--- modif. 02/06/25 (version 4, see version 2)
//    ------------------------------------------
//    after some comments of WK

      bool split = getSplitCondition();
//                 -------------------
      int iCat1 = 0;
      int iCat2 = nHlfCat;
      int kMir = 0;
      if( !split && kDetPart == 1 ) {
	iCat1 = nHlfCat;
        iCat2 = nCathode;
	kMir = 1;
      }
      list<CsRCMirrorNom*> lMirrNom = mirr->lMirrNom();
      list<CsRCMirrorNom*>::iterator in;
      for( in=lMirrNom.begin(); in!=lMirrNom.end(); in++ ) {

        Hep3Vector vPoPaDetW = vPoPaDetW_[kMir];
	float xPaLimMn, xPaLimMx, yPaLimMn, yPaLimMx;
	xPaLimMn = vPoPaDetW.x() - ddCluLix;
	xPaLimMx = vPoPaDetW.x() + ddCluLix;
	yPaLimMn = vPoPaDetW.y() - ddCluLiy;
	yPaLimMx = vPoPaDetW.y() + ddCluLiy;
	for( int kc=iCat1; kc<iCat2; kc++ ) {
	  CsRCCathode* cat = dets->ptrToCat( kc );
	  double hCatx = cat->hCatx();
	  double hCaty = cat->hCaty();
	  Hep3Vector vOffCatWw = dets->vOffCatW( kc );
	  if( vOffCatWw.x()+hCatx >= xPaLimMn &&
	      vOffCatWw.x()-hCatx <= xPaLimMx &&
	      vOffCatWw.y()+hCaty >= yPaLimMn &&
	      vOffCatWw.y()-hCaty <= yPaLimMx ) lSelCats_.push_back( kc );
	}
	if( !split ) break;
        iCat1 = nHlfCat;
        iCat2 = nCathode;
        kMir++;
      }
      //list<int>::iterator ii;
      //for( ii=lSelCats_.begin(); ii!=lSelCats_.end(); ii++ )
      //  cout << (*ii) << " ";
      //cout << endl;
//--------------------------------------------------------------------

      if( kDetPart == -1 ) {
        string str =
      "RICHONE, CsRCPartPhotons::doSelCathodes() : wrong particle reflection";
        CsErrLog::Instance()->mes( elFatal, str );
      }

      kDetPart_ = kDetPart;
      iCaPa_ = iCaPa;
      iXpPa_ = iXpPa;
      iYpPa_ = iYpPa;

      if( CsRCExeKeys::Instance()->kPrintPartPhotons() == 3 ) {
	list<int>::iterator ii;
	cout << endl;
	cout << " Particle  " << part->kPart() << ", Selected Cathodes :  ";
	for( ii=lSelCats_.begin(); ii!=lSelCats_.end(); ii++ ) {
	  cout << (*ii) << "  ";
	}
	cout << endl;
      }

  }


//===========================================================================
  bool CsRCPartPhotons::getSplitCondition( CsRCParticle* part ) {
//---------------------------------------------------------------


//--- Paolo  -  June 2002

      CsRCDetectors *dets = CsRCDetectors::Instance();
      CsRCMirrors *mirr = CsRCMirrors::Instance();
      CsRCHistos& hist = CsRCHistos::Ref();
      double xh, yh, wh;

      static double thCmax;

      static bool firstCall = true;
      if( firstCall ) {
        firstCall = false;
        CsRCRecConst *cons = CsRCRecConst::Instance();
        int kOffHi = cons->kOffHi();
        thCmax = acos( 1./cons->CFRefInd() );
        thCmax = 0.055;
//@@    --------------
      }

      bool split = false;

      Hep3Vector vPosIn0 = part->vPosIn();
      Hep3Vector vDirIn0 = part->vDirIn();
      //cout << vPosIn0 << "  " << vDirIn0 
      //     << "  " << part->vPoPaMir0()[part->kDetPart()] << endl;
      double tga = vDirIn0.x()/vDirIn0.z();
      Hep3Vector vDcPaPro( tga, 0., 1. );
      vDcPaPro = vDcPaPro.unit();
      Hep3Vector vPoPaPro( vPosIn0.x(), 0., vPosIn0.z() );
      //cout << vPoPaPro << "  " << vDcPaPro << endl;
      Hep3Vector vPoC0 = part->pMirPart()->vC0();
      double RR = part->pMirPart()->RR();
      Hep3Vector vPoMirPro = mirr->vImpMir( vPoPaPro, vDcPaPro, vPoC0, RR );
      //cout << vPoMirPro << endl;
      double ddSplit = sqrt( (vPoMirPro - vPosIn0).mag2() -
			     pow( (vPoMirPro - vPosIn0)*vDirIn0, 2 ) );
      double rCMaxMirr = part->pathLen() * tan( thCmax );
      //cout << ddSplit << "  " << rCMaxMirr << endl;

      //double hhSplit = part->vPoPaMir0()[part->kDetPart()].y();
      //cout << ddSplit << "  " << hhSplit << "  " << rCMaxMirr << endl;

      if( ddSplit < rCMaxMirr ) split = true;

      return split;

  }


//==========================================================================
  bool CsRCPartPhotons::mcsPartCorr( CsRCParticle* part, 
//---------------------------------------------------------
				     Hep3Vector &corrPo,
				     Hep3Vector &corrDc ) {

//--- Paolo - February 2001


      CsRCHistos& hist = CsRCHistos::Ref();
      int kh;
      double xh, yh, wh;

      static double coPo = 0;
      static double coDc = 0;

      static bool firstCall = true;
      if( firstCall ) {
        firstCall = false;

	CsRCExeKeys* key = CsRCExeKeys::Instance();
        key->acknoMethod( "CsRCPartPhotons::CMSPartCorr" );

        float partPathFr = CsRCRecConst::Instance()->partPathFr();
        CsRCDetectors *dets = CsRCDetectors::Instance();
        double zEntrWind = dets->zEntrWind();
        double zExitWind = dets->zExitWind();

	coPo = partPathFr*sqrt( partPathFr );
	coDc = sqrt( 3.*partPathFr ) / (zExitWind - zEntrWind);
	//coDc = sqrt( 3.)* coPo / (zExitWind - zEntrWind);      //   ???
      }
      bool flag = false;

//--- use extrapolation forward to the RICH exit window
//    (dummy (unit) cov matrix for MC tracks)
//    and measured point downstream.
//    ------------------------------
      corrPo.set(0., 0., 0.);
      corrDc.set(0., 0., 0.);
      if( part->zo() < 1.e+06 && part->vPosExW().z() < 1.e+06 ) {
        double ddXX = part->xo() - part->vPosExW().x();
        double ddYY = part->yo() - part->vPosExW().y();
        corrPo.setX( coPo*ddXX );
        corrPo.setY( coPo*ddYY );
        corrPo.setZ( 0. );
        corrDc.setX( coDc*ddXX );
        corrDc.setY( coDc*ddYY );
        corrDc.setZ( 0. );
        flag = true;
        //cout << corrPo << "  " << corrDc << "  " << flag << endl;
        xh = corrPo.x();
        yh = corrPo.y();
        if( hist.hRC3613 ) hist.hRC3613->Fill( xh, yh );
//                         ----------------------------
        xh = corrDc.x() * 1000.;
        yh = corrDc.y() * 1000.;
        if( hist.hRC3614 ) hist.hRC3614->Fill( xh, yh );
//                         ----------------------------
      }

      return flag;

  }



//==========================================================================
  Hep3Vector CsRCPartPhotons::photImpDet( Hep3Vector vPoPhotW,
//--------------------------------------------------------------------------
                                          Hep3Vector vDcPhoEmW,
					  const double zDetW,
					  const float RR ) {


//- Paolo  -  March 2001.


//- 'photon' impact on mirror (MWR) :
//  ---------------------------------
    CsRCMirrors *mirr = CsRCMirrors::Instance();
    Hep3Vector vPoC( 0., 0., 0. );
    Hep3Vector vPoPhoMir = mirr->vImpMir( vPoPhotW, vDcPhoEmW, vPoC, RR );
//                         ---------------------------------------------

//- normal to mirror at 'photon' impact :
//  -------------------------------------
    Hep3Vector vDcNoPhoMir = (1./RR) * vPoPhoMir;

//- 'photon' reflected direction :
//  ------------------------------
    double cosPhoMir = vDcNoPhoMir * vDcPhoEmW;
    Hep3Vector vDcPhoRefl = 2.*cosPhoMir * vDcNoPhoMir - vDcPhoEmW;

//- 'photon' impact on detector plane :
//  -----------------------------------
    double norm = (zDetW - vPoPhoMir.z()) / vDcPhoRefl.z();
    Hep3Vector vImpDet = vPoPhoMir + norm * vDcPhoRefl;

    return vImpDet;

  }


//===========================================================================
  bool CsRCPartPhotons::MIPReject( CsRCCluster* clu ) {
//-----------------------------------------------------


//- Paolo  -  April 2001.


//- rejection of MIP clusters on PH basis :
//  ---------------------------------------
    bool MIPReject = false;

    static float PHmipCut = CsRCRecConst::Instance()->PHmipCut();

    if( clu->PH() > PHmipCut ) MIPReject = true;

    return MIPReject;

  }


//===========================================================================
  bool CsRCPartPhotons::killHalo( CsRCCluster* clu ) {
//----------------------------------------------------


//- Paolo  -  October 2001.


//- rejection of 'beam halo' clusters :
//  -----------------------------------
//  pro-memoria

    return false;

  }


//===========================================================================
  vector<double> CsRCPartPhotons::getSplitLimits( CsRCParticle* part ) {
//----------------------------------------------------------------------


//--- Paolo  -  October 2001.

      CsRCDetectors *dets = CsRCDetectors::Instance();
      CsRCMirrors *mirr = CsRCMirrors::Instance();
      CsRCHistos& hist = CsRCHistos::Ref();
      double xh, yh, wh;

      static double tgCmax;

      static bool firstCall = true;
      if( firstCall ) {
        firstCall = false;
        CsRCRecConst *cons = CsRCRecConst::Instance();
        int kOffHi = cons->kOffHi();
        tgCmax = tan( acos( 1./cons->CFRefInd() ) );
        tgCmax = 0.055;
//@@    --------------
      }

      Hep3Vector vPosIn0 = part->vPosIn();
      Hep3Vector vDirIn0 = part->vDirIn();
      double tgb = vDirIn0.y()/vDirIn0.z();
      double tgbl = 0.;
      if( tgb > 0. ) tgbl = (tgb-tgCmax)/(1+tgb*tgCmax);
      if( tgb < 0. ) tgbl = (tgb+tgCmax)/(1-tgb*tgCmax);
      double tga = vDirIn0.x()/vDirIn0.z();
      double norm = 1.;
      norm = sqrt( tga*tga + tgbl*tgbl + 1. );
      Hep3Vector vDcPaL( tga/norm, tgbl/norm, 1./norm );
      Hep3Vector vPoC0 = part->pMirPart()->vC0();
      double RR = part->pMirPart()->RR();
      Hep3Vector vPoMirL = mirr->vImpMir( vPosIn0, vDcPaL, vPoC0, RR );

      vector<double> vySplitLLimit;
      vySplitLLimit.push_back( 0. );
      vySplitLLimit.push_back( 0. );
      bool exe = false;
      if( kDetPart_ == 0  &&  vPoMirL.y() < 0. ) exe = true;
      if( kDetPart_ == 1  &&  vPoMirL.y() > 0. ) exe = true;
      if( exe ) {
	vySplitLLimit.clear();

        double tgap = (tga+tgCmax)/(1-tga*tgCmax);
        double tgam = (tga-tgCmax)/(1+tga*tgCmax);
        Hep3Vector vPoPaC( vPosIn0.x(), 0., vPosIn0.z() );
        Hep3Vector vDcPaC( vDirIn0.x(), 0., vDirIn0.z() );
	Hep3Vector vPoMirC = mirr->vImpMir( vPoPaC, vDcPaC, vPoC0, RR );

	int kDetPaw = 0;
	list<CsRCPhotonDet*> lPhoDet = dets->lPhoDet();
	list<CsRCPhotonDet*>::iterator id;
	std::list<CsRCMirrorNom*> lMirrNom = mirr->lMirrNom();
	std::list<CsRCMirrorNom*>::iterator im;
	im = lMirrNom.begin();
	for( id=lPhoDet.begin(); id!=lPhoDet.end(); id++ ) {

	  Hep3Vector vPoC0 = (*im)->vC0();
	  double RR = (*im)->RR();
	  im++;

          Hep3Vector vDcNor(0., 0., 0.);
          double cosMir = 0.;
	  Hep3Vector vDcRefl(0., 0., 0.);
  	  Hep3Vector vDcDet = dets->vDcDetv( kDetPaw );
	  Hep3Vector vPoDet = dets->vDet0v( kDetPaw );

	  norm = sqrt( tgap*tgap + 1. );
          Hep3Vector vDcPaP( tgap/norm, 0., 1./norm );
	  Hep3Vector vPoMirP = mirr->vImpMir( vPoPaC, vDcPaP, vPoC0, RR );
	  Hep3Vector vDcMirP = (vPoMirP - vPosIn0).unit();
          vDcNor = (1./RR) * (vPoMirP - vPoC0);
          cosMir = vDcNor * vDcMirP;
          vDcRefl = 2.*cosMir * vDcNor - vDcMirP;
          norm = ( ( vPoDet - vPoMirP ) * vDcDet ) / ( vDcRefl * vDcDet );
          Hep3Vector vPoDetP = vPoMirP + norm * vDcRefl;

	  norm = sqrt( tgam*tgam + 1. );
          Hep3Vector vDcPaM( tgam/norm, 0., 1./norm );
	  Hep3Vector vPoMirM = mirr->vImpMir( vPoPaC, vDcPaM, vPoC0, RR );
	  Hep3Vector vDcMirM = (vPoMirM - vPosIn0).unit();
          vDcNor = (1./RR) * (vPoMirM - vPoC0);
          cosMir = vDcNor * vDcMirM;
          vDcRefl = 2.*cosMir * vDcNor - vDcMirM;
          norm = ( ( vPoDet - vPoMirM ) * vDcDet ) / ( vDcRefl * vDcDet );
          Hep3Vector vPoDetM = vPoMirM + norm * vDcRefl;
	  //cout << vPoMirP << "  " << vPoMirC << "  " << vPoMirM << endl;
	  //cout << vPoDetP << "  " << vPoDetM << endl;

  	  HepVector vrot( 3, 0 );
	  for( int j=0; j<3; j++ ) vrot[j] = (vPoDetP - vPoC0)[j];
	  HepVector vPoDetPWw = (*id)->rotMatrix() * vrot;
	  Hep3Vector vPoDetPW( vPoDetPWw[0], vPoDetPWw[1], vPoDetPWw[2] );
	  for( int j=0; j<3; j++ ) vrot[j] = (vPoDetM - vPoC0)[j];
          HepVector vPoDetMWw = (*id)->rotMatrix() * vrot;
	  Hep3Vector vPoDetMW( vPoDetMWw[0], vPoDetMWw[1], vPoDetMWw[2] );
	  float yDetLL = yDetLLimit( kDetPaw );
	  double splitLmin = 0.;
	  if( kDetPaw == 0 ) {
	    splitLmin = vPoDetPW.y();
	    if( splitLmin > vPoDetMW.y() ) splitLmin = vPoDetMW.y();
	    if( splitLmin < yDetLL ) splitLmin = yDetLL;
	  }
	  if( kDetPaw == 1 ) {
	    splitLmin = vPoDetPW.y();
	    if( splitLmin < vPoDetMW.y() ) splitLmin = vPoDetMW.y();
	    if( splitLmin > yDetLL ) splitLmin = yDetLL;

	  }
	  vySplitLLimit.push_back( splitLmin );

	  //cout << vPoMirL.y() << "  " << kDetPaw << " - " 
	  //     << vPoPaDetW_[kDetPaw] << " P "
	  //     << vPoDetPW << " M " << vPoDetMW << endl;
          kDetPaw++;
	}
	vySplitLLim_.push_back( vySplitLLimit[0] );
	vySplitLLim_.push_back( vySplitLLimit[1] );
	//std::cout << kDetPart_ << "  " << vPoPaDetW_[kDetPart_].y()
	//	  << "  " << vySplitLLimit[0] << "  " << vySplitLLimit[1] 
	//	  << std::endl;
      }

//--- monitor the limits :
//    --------------------
      xh = vySplitLLimit[0];
      yh = vySplitLLimit[1];
      if( hist.hRC3508 ) hist.hRC3508->Fill( xh, yh );
//hh                     ----------------------------
      xh = vPoPaDetW_[kDetPart_].y();
      yh = vPoPaDetW_[1-kDetPart_].y();
      if( hist.hRC3509 ) hist.hRC3509->Fill( xh, yh );
//hh                     ----------------------------

      return vySplitLLimit;

    }


//===========================================================================
  float CsRCPartPhotons::yDetLLimit( const int kDet ) {
//-----------------------------------------------------


//--- Paolo  -  October 2001.

//--- rev. December 2006


      static float yyLimDet[2];

      static bool firstCall = true;
      if( firstCall ) {
        firstCall = false;

	CsRCDetectors *dets = CsRCDetectors::Instance();
        list<CsRCCathode*> lCathodes = dets->lCathodes();
	double hCaty = 0.;
        float yyLimMin =  10000.;
        float yyLimMax = -10000.;
	int kCat = -1;
        list<CsRCCathode*>::iterator it;
        for( it=lCathodes.begin(); it!=lCathodes.end(); it++ ) {
          float yyOff = (*it)->vOffCatW().y();
          if( yyOff > 0. ) {
	    if( yyOff < yyLimMin ) {
	      yyLimMin = yyOff;
	      hCaty = (*it)->hCaty();
	    }
	  }
          if( yyOff < 0. ) {
	    if( yyOff > yyLimMax ) {
	      yyLimMax = yyOff;
	      hCaty = (*it)->hCaty();
	    }
	  }
	}
	yyLimDet[0] = yyLimMin - hCaty;
	yyLimDet[1] = yyLimMax + hCaty;

      }
      return yyLimDet[kDet];

  }


//========================================================================
  double CsRCPartPhotons::sigmaPhoRec( const double mom ) const {
//---------------------------------------------------------------

//- interface function :
//  --------------------
//  December  1999,  rev. December  2002


    CsRCHistos& hist = CsRCHistos::Ref();
    double xh, yh, wh;

    double sigmaPhot = 1000000.;
    float sigma = 0.;

    //sigma = sigmaPhot5( beta );
//@@--------------------------

    sigma = sigmaPhot7( mom );                 //   02/12/10
//@@-------------------------
    sigma *= getCorrFactor();
//@@------------------------

    if( pPart() ) {
      //  float ddpdet = pPart()->ddPaDet()[pPart()->kDetPart()];
      //  sigma = sigmaPhot8( ddpdet );            //   04/02/05
//@@------------------------------
    }

    if( sigma > 0. ) sigmaPhot = sigma;

//- monitor error values
    xh = mom;
    yh = sigmaPhot;
    if( hist.hRC1701 ) hist.hRC1701->Fill( xh, yh );
//hh                   ----------------------------

    //std::cout << "sigmaPhoRec  " << sigmaPhot << std::endl;

    return  sigmaPhot;

  }


//========================================================================
  double CsRCPartPhotons::sigmaPhot5( const double beta ) const {
//---------------------------------------------------------------


//    December  1999,  rev. April  2001
//    Revised  August 2002

//--- sigma single-photon from MC :
//    -----------------------------

      static float par[10];

      static bool firstCall = true;
      if( firstCall ) {
	firstCall = false;

        CsOpt* opt = CsOpt::Instance();
        CsRCRecConst *cons = CsRCRecConst::Instance();

//----- defaults
//      ( from Dima's 48k file ),
//      ( Paolo file  rc5-48k-20sg-cor.hist = open cuts ) :
//      ---------------------------------------------------
        par[0] = 62.2;                              //   10/11/99

//      from Vadim's 6*3k ref files/lepto_full ( co-18k-15t.hist )
//      ----------------------------------------------------------
	//par[0] =   47.8;                             //   4/4/01

//----- from rich1.options :
//      --------------------
	vector<float> vPar;
        bool boo = opt->CsOpt::getOpt( "RICHONE", "sigmaPhoRec", vPar );
        if( boo ) {
	  for( unsigned int k=0; k<vPar.size(); k++ ) par[k] = vPar[k];

	  if( cons->printConsts() ) {
	    cout << " RICHONE, CsRCPartPhotons::sigmaPhot5 :   ";
	    for( unsigned int k=0; k<vPar.size(); k++ ) cout << par[k] <<"  ";
	    cout << endl;
	  }
	} else {
	  //cout << " RICHONE, CsRCPartPhotons::sigmaPhot5 :   ";
	  //cout << "NO parameters read, default used!" << endl;
          string mess = "RICHONE, CsRCPartPhotons::sigmaPhot5 : ";
	  string err = "NO parameters read, default used!";
	  mess.append( err );
	  CsErrLog::Instance()->mes( elError, mess );
	}
      }

//--- MC : from sigphorec.kumac/sigmapho5.f - histo.s : 3521 (+off)
//    WARNING : histo 3521 does not exist anymore !!!
//    --------------------------------------------------------
//--- with mcs correction
//    -------------------
//--- f = f0 + p0 * (1-beta**2)
//    -------------------------

      double betaw = beta;

      //float sigmaPhot = 0.795 + par[0] * (1.- betaw*betaw);
//--------------------------------------------------------
      float sigmaPhot = par[0] + par[1] * (1.- betaw*betaw);  //   020826
//---------------------------------------------------------

      //cout << "sigmaPhot5  " << sigmaPhot << endl;

      return  sigmaPhot;

  }


//========================================================================
  double CsRCPartPhotons::sigmaPhot7( const double mom ) const {
//--------------------------------------------------------------


//- Paolo   -   December 2002

//- sigma single-photon from Data :
//  -------------------------------
//  n-phomom.kumac & histo 3533 (3535/3540)

    static float par[10];

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      CsOpt* opt = CsOpt::Instance();
      CsRCRecConst *cons = CsRCRecConst::Instance();

//--- defaults ( run 20330, 2002 ) :   //   02/12/10
//    -----------------------------
      par[0] =  1.34;
      par[1] = -0.019;
      par[2] =  0.00026;
//--- defaults ( run 51908, 2006 ) :   //   06/11/28
//    -----------------------------
      par[0] =  1.31;
      par[1] = -0.0031;
      par[2] = -0.00007;

//--- from rich1.options :
//    --------------------
      vector<float> vPar;
      bool boo = opt->CsOpt::getOpt( "RICHONE", "sigmaPhoRec", vPar );
      if( boo ) {
	for( unsigned int k=0; k<vPar.size(); k++ ) par[k] = vPar[k];

	if( cons->printConsts() ) {
	  cout << " RICHONE, CsRCPartPhotons::sigmaPhot7 :   ";
	  for( unsigned int k=0; k<vPar.size(); k++ ) cout << par[k] <<"  ";
	    cout << endl;
	}
      } else {
	//cout << " RICHONE, CsRCPartPhotons::sigmaPhot7 :   ";
	//cout << "NO parameters read, default used!" << endl;
        string mess = "RICHONE, CsRCPartPhotons::sigmaPhot7 : ";
	string err = "NO parameters read, default used!";
	mess.append( err );
	CsErrLog::Instance()->mes( elError, mess );
      }
    }

//- f = p0 + p1* mom + p2* mom*mom
//  ------------------------------

    float sigmaPhot = par[0] + par[1]* mom + par[2]* mom*mom;
//----------------------------------------------------------

    //cout << "   sigmaPhot7  " << sigmaPhot << "  " << mom << endl;

    return  sigmaPhot;

  }


//========================================================================
  double CsRCPartPhotons::sigmaPhot8( const double ddpdet ) const {
//-----------------------------------------------------------------


//- Paolo   -   February 2004

//- sigma single-photon from Data :
//  -------------------------------
//  n-phoddet.kumac & histo 3545

    static float par[10];

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      CsOpt* opt = CsOpt::Instance();
      CsRCRecConst *cons = CsRCRecConst::Instance();

//--- defaults ( run 31815, 2003 ) :   //   13/02/05
//    -----------------------------
      par[0] =  1.29;
      par[1] =  0.;
      par[2] =  0.0000025;

//--- from rich1.options :
//    --------------------
      vector<float> vPar;
      //bool boo = opt->CsOpt::getOpt( "RICHONE", "sigmaPhoRec", vPar );
      bool boo = opt->CsOpt::getOpt( "RICHONE", "sigmaPhoPatt", vPar );
      if( boo ) {
	for( unsigned int k=0; k<vPar.size(); k++ ) par[k] = vPar[k];

	if( cons->printConsts() ) {
	  cout << " RICHONE, CsRCPartPhotons::sigmaPhot8 :   ";
	  for( unsigned int k=0; k<vPar.size(); k++ ) cout << par[k] <<"  ";
	    cout << endl;
	}
      } else {
	//cout << " RICHONE, CsRCPartPhotons::sigmaPhot8 :   ";
	//cout << "NO parameters read, default used!" << endl;
        string mess = "RICHONE, CsRCPartPhotons::sigmaPhot8 : ";
	string err = "NO parameters read, default used!";
	mess.append( err );
	CsErrLog::Instance()->mes( elError, mess );
      }
    }

//- f = p0 + p1* ddpdet + p2* ddpdet*ddpdet
//  ---------------------------------------
    float ddp = ddpdet;
    if( ddp > 700. ) ddp = 700.;
//@@---------------------------

    float sigmaPhot = par[0] + par[1]* ddp + par[2]* ddp*ddp;
//----------------------------------------------------------

    //cout << "   sigmaPhot8  " << sigmaPhot << "  " << mom << endl;

    return  sigmaPhot;

  }


//===========================================================================
  double CsRCPartPhotons::getCorrMom( const double sigPhoRec ) const {
//--------------------------------------------------------------------

//- Paolo  -  December 2002

    double corrMom = sigPhoRec;             // default

    corrMom = sigPhoRec / 1.05;             // 2002

    if( CsRichOne::Instance()->UpRICHJob() ) corrMom = 1.;
//  -----------------------------------------------------

    return  corrMom;

  }


//===========================================================================
  double CsRCPartPhotons::getCorrFactor() const {
//-----------------------------------------------

//- Paolo  -  December 2002

    double corrFact = 1.;                  // default

    corrFact = 1.15;                       // 2002

    if( CsRichOne::Instance()->UpRICHJob() ) corrFact = 1.;
//  ------------------------------------------------------

    return  corrFact;

  }


//===========================================================================
  long double CsRCPartPhotons::getLikelihood( const double theIpo,
//----------------------------------------------------------------
					      int &nPhotons ) {
//- interface method
//  ----------------
//- Paolo  -  December 2002
//- rev.      January  2005

    CsRCRecConst* cons = CsRCRecConst::Instance();

//- getLikeAll timing :
//  -------------------
    //CsStopwatch chronos( 1 );
    //static int chroLk = 0;
    //static double likeTime = 0.;
    //static int calls = 0;
    //static const int printF = CsRCExeKeys::Instance()->printFrequency();
    //static int kEventP = 0;
    //if( calls == 0 ) {
    //  chroLk = chronos.start();
    //  likeTime = 0.;
    //}
//  -----------
    double likeDV = cons->likeDefVa();
    //long double like = 0.;
    long double like = likeDV;

    if( cons->backgrType() == "BKGMAP" ) {
      like = getLikeAllMap( theIpo, nPhotons );      //   OBSOLETE!
//           ---------------------------------
    }
    else {
      like = getLikeAll( theIpo, nPhotons );
//           ------------------------------
    }

//  -----------
    //likeTime += chronos.inter( chroLk );
    //calls++;
    //int kEvent = CsRichOne::Instance()->kEvent();
    //if( kEvent % printF == 0  &&  kEvent != kEventP ) {
    //  std::cout << "  getLikelihood time = " << likeTime
    //            << "  calls " << calls << std::endl;
    //}
    //kEventP = kEvent;
//  -----------

    return  like;

  }


//===========================================================================
  long double CsRCPartPhotons::getLikelihoodPMT( const double theIpo,
//----------------------------------------------------------------
					         int &nPhotons ) {
//- interface method
//  ----------------
//- Paolo  -  May 2008

    CsRCRecConst* cons = CsRCRecConst::Instance();
    double likeDV = cons->likeDefVa();
    long double like = likeDV;

    setLikePMTONLY( true );
//  ----------------------
    if( cons->backgrType() == "BKGMAP" ) {
      like = 0.;                           //   OBSOLETE!
//    ---------
    }
    else {
      like = getLikeAll( theIpo, nPhotons );
//           ------------------------------
    }
    setLikePMTONLY( false );
//  -----------------------

    return  like;
  }


//===========================================================================
  long double CsRCPartPhotons::getLikelihoodAPV( const double theIpo,
//----------------------------------------------------------------
					         int &nPhotons ) {
//- interface method
//  ----------------
//- Paolo  -  May 2008

    CsRCRecConst* cons = CsRCRecConst::Instance();
    double likeDV = cons->likeDefVa();
    long double like = likeDV;

    setLikeAPVONLY( true );
//  ----------------------
    if( cons->backgrType() == "BKGMAP" ) {
      like = 0.;                           //   OBSOLETE!
//    ---------
    }
    else {
      like = getLikeAll( theIpo, nPhotons );
//           ------------------------------
    }
    setLikeAPVONLY( false );
//  -----------------------

    return  like;
  }


  extern "C" {
    float prob_( const float&, const int& );
    float erf_( const float& );
  }

//===========================================================================
  long double CsRCPartPhotons::getLikeAll( const double theIpo,
//-------------------------------------------------------------
				           int &nPhotons ) {

//- compute  Likelihood value for a given hypothesis (All photons)
//  --------------------------------------------------------------
//- Paolo  -  December 2002
//-           rev. May 2004
//-           rev. May 2006


    CsRCExeKeys* key = CsRCExeKeys::Instance();
    CsRCRecConst* cons = CsRCRecConst::Instance();

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      key->acknoMethod( "CsRCPartPhotons::getLikelihood(All)" );
      std::cout << "               Type  " << cons->backgrType() << std::endl;
    }

    static CsRCLikeAll* likeAll = NULL;
    if( cons->backgrType() == "02" ) likeAll = new CsRCLikeAll02( this );
    if( cons->backgrType() == "03" ) likeAll = new CsRCLikeAll03( this );
    if( cons->backgrType() == "04" ) likeAll = new CsRCLikeAll04( this );
    static int kEventC = -1;
    int kEvent = CsRichOne::Instance()->kEvent();
    static const CsRCPartPhotons* pPartPhotC = NULL;
    if( cons->backgrType() == "05" ) {
      if( kEvent != kEventC  ||  pPartPhotC != this ) {
	delete likeAll;
	likeAll = new CsRCLikeAll05( this );
      }
      kEventC = kEvent;
      pPartPhotC = this;
    }
//@@--------------------------------------------------------------------
    if( !likeAll ) {
      cout << " RICHONE, CsRCPartPhotons::getLikeAll() :   ";
      cout << "reference to a NOT existing background parametrization!";
      cout << endl;
      string mess = "RICHONE, CsRCPartPhotons::getLikeAll() : ";
      string err = "Reference to a NOT existing background parametrization!";
      mess.append( err );
      CsErrLog::Instance()->mes( elFatal, mess );
    }
    bool likeRatio = cons->likeRatio();
    if( cons->backgrType() == "04" ) likeRatio = true;

    list<CsRCPhoton*> lPhotons = lPhotons_;
    list<CsRCPhoton*>::iterator ih;

//- Provisional: to be improved!
//  ----------------------------
    //^double theIpow = 0.;
    //^if( lPhotons.size() > 0 )
    //^  theIpow = lPhotons.front()->getThetaIpo( theIpo );
    CsRCDetectors* dets = CsRCDetectors::Instance();
    double theIpow = theIpo;
    if( dets->ptrToCat( this->pPart()->getClosestCat() )->isPMT() )
      theIpow = thetaUVtoVS( theIpo );

    double normS = likeAll->normSignal( theIpow );
//                 ------------------------------
    double normB = likeAll->normBackgr( theIpow );
//                 ------------------------------
    double likeNorm = normS + normB;

    if( pionONLY_ ) {
      normSg_ = normS;
      normBk_ = normB;
      //std::cout << normSg_ << "  " << normBk_ << std::endl;
    }

    CsRCHistos& hist = CsRCHistos::Ref();
    int level = hist.levelBk();
    bool doHist = false;
    if ( level >= 2 ) doHist = true;
    if( doHist  &&  likeONLY_ ) {
      double xh, yh, wh;
      xh = lPhotons.size() - likeNorm;
      if( hist.hRC1571 ) hist.hRC1571->Fill( xh );
//    -------------------------------------------
    }
    //if( likeONLY_ ) {
    //std::cout << theIpow << "  -  " << normS << "  " << normB << "  "
    //	  << likeNorm << "  " << lPhotons.size() << std::endl;
    //}

    double likeDV = cons->likeDefVa();
//@@---------------------------------
    long double pLike = -1.;
    long double pLikeBk = -1.;
    long double logLike = 0.;
    long double logLikeBk = 0.;

    bool isFirst = true;
    int kPhot = 0;
    for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
      if( !(*ih)->flag() ) continue;

      if( isFirst ) (*ih)->setLikeFirst( true );
      isFirst = false;

//(080506)
      if( likeAPVONLY_  &&   (*ih)->isPMT() ) continue;
      if( likePMTONLY_  &&   (*ih)->isAPV() ) continue;
//@@--------------------------------------------------

      double theIpow = (*ih)->getThetaIpo( theIpo );

      double signW = likeAll->likeSignal( (*ih), theIpow );
//                   -------------------------------------
      double backW = likeAll->likeBackgr( (*ih), theIpow );
//                   -------------------------------------

      //if( fabs(signW) > 1000. || fabs(backW) > 1000.) {
      //if( backW == 0.) {
      //if( likeONLY_ ) {
      //std::cout << setprecision( 8 ) << "PPhot  ";
      //std::cout << kPhot << "  " << (*ih)->kPhot() << "  "
      //<< signW << "  " << backW << std::endl;
      //}
      //}

      //checkLikeRatios( (*ih), signW, backW );
//    --------------------------------------

      if( theIpow <= 0. ) signW = 0.;
      double likeW = signW + backW;

      pLike *= likeW;
      pLike = fabs( pLike );
      pLikeBk *= backW;
      pLikeBk = fabs( pLikeBk );
      if( likeW > 0. ) logLike += log( likeW );
      if( backW > 0. ) logLikeBk += log( backW );

      (*ih)->setLikeFirst( false );
      kPhot++;
    }   //   end for on Photons


    nPhotons = kPhot;
    if( nPhotons > 0 ) {
      long double pLike0 = pLike;
      long double power = 1./float( nPhotons );
      if( pLike > 0. ) pLike = pow( pLike, power );
      if( likeNorm > 0. ) {
        pLike /= likeNorm;
      } else {
	pLike = likeDV;
      }

      if( cons->backgrType() == "05" ) {
	//@@pLike = pLike;
//----- Standard 2003-04 production Like -----
        pLike *= exp( -normS * power ) * likeNorm;

        //0905/pLike = pLike0;
	//0905/pLike /= pow( likeNorm, double( nPhotons ) );
        //0905/pLike *= exp( -normS );
	//0905/pLike *= likeNorm * exp( -normS );
	//@pLike *= exp( -likeNorm * power ) * likeNorm;
	//-pLike = pLike0;
	//-pLike = pLike0 / pow( likeNorm, nPhotons );
        //long double nFact = 1.;
	//for( int k=1; k<=nPhotons; k++ ) nFact *= double( k );
	//-pLike /= nFact;
	if( cons->likeLog() ) {
	  //@@pLike = logLike - likeNorm;
	  pLike = logLike;
	  //@pLike = logLike - nPhotons*log( likeNorm );
	  //-pLike -= log( nFact );
	  //0905/pLike = logLike - nPhotons*log( likeNorm );
	  //0905/pLike = logLike/nPhotons - log( likeNorm );
	  //0905/pLike = logLike/nPhotons - normS/nPhotons;
	}
      }
      //std::cout << setprecision( 10 );
      //if( likeONLY_ ) cout << pLike << "  " << logLike << "  "
      //   << likeNorm << endl; 
      //cout << pLike << "  " << normB << "  " << normS << "  " 
      //     << likeNorm << "  " << theIpow << "  " << power << endl;

      if( pLikeBk > 0. ) pLikeBk = pow( pLikeBk, power );
      if( normB > 0. ) {
        pLikeBk /= normB;
      } else {
        pLikeBk = likeDV;
      }

      logLike *= power;
      if( likeNorm > 0. ) logLike -= log( likeNorm );
      logLikeBk *= power;
      if( normB > 0. ) logLikeBk -= log( normB );

      if( likeRatio  &&  theIpow > 0. ) {
        if( pLikeBk > 0. ) {
          pLike /= pLikeBk;
        } else {
          pLike = likeDV;
        }
	logLike -= logLikeBk;
      }

    } else {
      pLike = likeDV;
      logLike = likeDV;
    }
    //cout << pLike << "  " << pLikeBk << endl;

    if( !(cons->backgrType() == "05") ) delete likeAll;

    //if( isnan( pLike ) ) std::cout << "getLikeAll " << theIpo << "  "
    //	               <<  theIpow << "  " << likeONLY_ << std::endl;
    return pLike;

  }


//===========================================================================
  void CsRCPartPhotons::getLikeDeriv( const int kPaTy, const double theIpo,
//-------------------------------------------------------------------------
				      const double pLikeTOT,
				      double& derivLKAPV, double& derivLKPMT,
				      double& pLikeLKAPV, double& pLikeLKPMT,
				      int& nPhotAPV, int& nPhotPMT ) {

//- Paolo  -  May  2008

    CsRCExeKeys* key = CsRCExeKeys::Instance();
    CsRCRecConst* cons = CsRCRecConst::Instance();

    double likeDV = cons->likeDefVa();
    float CFRefInd = cons->CFRefInd();
    float CFRefIndUV = cons->CFRefIndUV();
    float CFRefIndVS = cons->CFRefIndVS();

    int nProb = cons->outBufferSize();
    float momPart = pPart_->mom();
    float massPart = pPart_->mass( kPaTy );

    double theIpoD = 0.;
    double theIpoDP = 0.;
    double theIpoDM = 0.;
    double likeDerAPV = 0.;
    double likeDerPMT = 0.;
    double likeDer = 0.;
    long double pLikeAPV = likeDV;
    long double pLikePMT = likeDV;
    long double pLike = likeDV;
    long double pLikeD = likeDV;
    long double pLikeDP = likeDV;
    long double pLikeDM = likeDV;
    int nPhot = 0;
    double CFRefIndD = 0.;
    double CFRefIndDP = 0.;
    double CFRefIndDM = 0.;

    bool twoPtDer = false;
//@@---------------------

    bool print = false;

//- WARNING : take into account the specific form of the Likelihood,
//  i.e. NOT log, extended Like, (N)root( Like )
  
//- compute dLike/dn using dn as DCFRefInd :
    //double DCFRefInd = 0.000025;
    double DCFRefInd = 0.000050;      // ~ 3.3 % for UV (3.7 % for VS)
//@@---------------------------

    if( CsRichOne::Instance()->UpRICHJob()  &&  nProb > 21 ) {

      twoPtDer = true;      // only little more accurate...
//@@-----------------
      if( twoPtDer ) DCFRefInd /= 2.;

      CFRefIndDP = CFRefInd + DCFRefInd;
      theIpoDP = getThetaIpo( CFRefIndDP, momPart, massPart );
//               --------------------------------------------
      if( twoPtDer ) {
        CFRefIndDM = CFRefInd - DCFRefInd;
        theIpoDM = getThetaIpo( CFRefIndDM, momPart, massPart );
//                 --------------------------------------------
      }
//--- APV's
      pLikeAPV = getLikelihoodAPV( theIpo, nPhotAPV );
//               ------------------------------------
      if( nPhotAPV > 0 ) {
	if( theIpoDP > 0. ) {
	  pLikeDP = getLikelihoodAPV( theIpoDP, nPhotAPV );
//                  --------------------------------------
	  likeDerAPV = (pLikeDP - pLikeAPV) / DCFRefInd;
	  //std::cout << "likeDerAPV  " << likeDerAPV << "  ";
	}
	if( twoPtDer ) {
	  if( theIpoDP > 0.  &&  theIpoDM > 0. ) {
	    pLikeDM = getLikelihoodAPV( theIpoDM, nPhotAPV );
//                    --------------------------------------
	    likeDerAPV = (pLikeDP - pLikeDM) / (2.*DCFRefInd);
	    //std::cout << likeDerAPV << "  ";
	  }
	}
      }
//--- PMT's
      pLikePMT = getLikelihoodPMT( theIpo, nPhotPMT );
//               ------------------------------------
      if( nPhotPMT > 0 ) {
	if( theIpoDP > 0. ) {
	  pLikeDP = getLikelihoodPMT( theIpoDP, nPhotPMT );
//                  --------------------------------------
	  likeDerPMT = (pLikeDP - pLikePMT) / DCFRefInd;
	  //std::cout << "likeDerPMT  " << likeDerPMT << "  ";
	}
	if( twoPtDer ) {
	  if( theIpoDP > 0.  &&  theIpoDM > 0. ) {
	    pLikeDM = getLikelihoodPMT( theIpoDM, nPhotPMT );
//                    --------------------------------------
	    likeDerPMT = (pLikeDP - pLikeDM) / (2.*DCFRefInd);
	    //std::cout << likeDerPMT << std::endl;
	  }
	}
      }
      pLikeLKAPV = pLikeAPV;
      pLikeLKPMT = pLikePMT;
      //derivLKAPV = likeDerAPV * pLikePMT;      // obso
      //derivLKPMT = likeDerPMT * pLikeAPV;      // obso
//--- added 081001
      int nPhotTOT = nPhotAPV + nPhotPMT;
      double factAPV = 0.;
      double factPMT = 0.;
      if( nPhotTOT > 0 ) {
	if( pLikeLKAPV > likeDV ) {
	  factAPV = pLikeTOT / pLikeLKAPV;
	  factAPV *= double( nPhotAPV )/double( nPhotTOT );
	}
	if( pLikeLKPMT > likeDV ) {
	  factPMT = pLikeTOT / pLikeLKPMT;
	  factPMT *= double( nPhotPMT )/double( nPhotTOT );
	}
	derivLKAPV = factAPV * likeDerAPV;
	derivLKPMT = factPMT * likeDerPMT;
      } else {
	derivLKAPV = 0.;
	derivLKPMT = 0.;
      }
//--- TOT for TEST
      long double pLikeTOT = likeDV;
      double likeDerTOT = 0.;
      //bool test = true;
      bool test = false;
      if( test ) {
	pLikeTOT = getLikelihood( theIpo, nPhotTOT );
//                 ---------------------------------
        if( theIpoDP > 0. ) {
	  pLikeDP = getLikelihood( theIpoDP, nPhotTOT );
//                  -----------------------------------
	  likeDerTOT = (pLikeDP - pLikeTOT) / DCFRefInd;
	  if( twoPtDer ) {
	    if( theIpoDP > 0.  &&  theIpoDM > 0. ) {
	      pLikeDM = getLikelihood( theIpoDM, nPhotTOT );
//                      -----------------------------------
	      likeDerTOT = (pLikeDP - pLikeDM) / (2.*DCFRefInd);
	      //std::cout << likeDerTOT << "  ";
	    }
	  }
	}
      }
//--- 081001
      print = false;
      if( print ) {
	if( nPhotAPV > 0  &&  nPhotPMT > 0 ) {
	  CsRCLikeAll* likeAll = new CsRCLikeAll05( this );
	  double normS = likeAll->normSignal( theIpo );
	  long double xAPV = pLikeAPV * exp( normS/float( nPhotAPV ) );
	  xAPV = pow( xAPV, nPhotAPV );
	  long double xPMT = pLikePMT * exp( normS/float( nPhotPMT ) );
	  xPMT = pow( xPMT, nPhotPMT );
	  long double yTOT = 1./float( nPhotTOT );
	  long double xTOT = exp( -normS*yTOT ) * pow( xAPV*xPMT, yTOT );
          std::cout << std::endl;
          std::cout << nPhotTOT << " = " << nPhotAPV << " + " << nPhotPMT
		    << std::endl;
          std::cout << setprecision(12);
          std::cout << pLikeTOT << "  =  " << pLikeAPV << "  &  " << pLikePMT
		    << "  ?  " << xTOT << std::endl;
          std::cout << setprecision(2);
          std::cout << derivLKAPV << " = " << factAPV << " * " << likeDerAPV
		    << "   /   ";
          std::cout << derivLKPMT << " = " << factPMT << " * " << likeDerPMT
		    << std::endl;
          std::cout << derivLKAPV+derivLKPMT << " = " << likeDerTOT
		    << std::endl;
	}
      }

      if( cons->likeLog() ) {
//----- WORNING : NOT YET IMPLEMENTED!
	derivLKAPV = likeDerAPV;
	derivLKPMT = likeDerPMT;
      }

    } else {

      twoPtDer = true;      // only little more accurate...
//@@-----------------
      if( twoPtDer ) DCFRefInd /= 2.;

      CFRefIndDP = CFRefInd + DCFRefInd;
      theIpoDP = getThetaIpo( CFRefIndDP, momPart, massPart );
//               --------------------------------------------
      if( twoPtDer ) {
        CFRefIndDM = CFRefInd - DCFRefInd;
        theIpoDM = getThetaIpo( CFRefIndDM, momPart, massPart );
//                 --------------------------------------------
      }
      pLike = getLikelihood( theIpo, nPhot );
//            ------------------------------
      if( nPhot > 0 ) {
        if( theIpoDP > 0. ) {
	  pLikeDP = getLikelihood( theIpoDP, nPhot );
//                  --------------------------------
	  likeDer = (pLikeDP - pLike) / DCFRefInd;
        }
	if( twoPtDer ) {
	  if( theIpoDP > 0.  &&  theIpoDM > 0. ) {
	    pLikeDM = getLikelihood( theIpoDM, nPhot );
//                    --------------------------------
	    likeDer = (pLikeDP - pLikeDM) / (2.*DCFRefInd);
	  }
	}
      }
      nPhotAPV = nPhot;
      pLikeLKAPV = pLike;
      pLikeLKPMT = 1.;
      derivLKAPV = likeDer;
      derivLKPMT = 0.;

      print = false;
      if( print ) {
	std::cout << kPaTy << "  " << nPhot << "  " << setprecision(12)
		  << pLike << "  " << setprecision(2)
		  << likeDer << std::endl;
      }

      //std::cout << kPaTy << "  " << pLike << "  " << pLikeDP << "  "
      //     << likeDer << std::endl;
    }

    return;
  }


//===========================================================================
  long double CsRCPartPhotons::getLikeAllMap( const double theIpo,
//----------------------------------------------------------------
				              int &nPhotons ) {

//- compute  Likelihood value for a given hypothesis (All photons)
//  --------------------------------------------------------------
//- Stefano  -  November 2004
//-
//- OBSOLETE!

    std::cout << "getLikeAllMap : method NOT YET UPGRADED to RICH with PMT.s!"
	      << std::endl;
    string mess =
      "RICHONE, CsRCPartPhotons::getLikeAllMap() : ";
    string err = "method NOT YET UPGRADED to RICH with PMT.s!";
    mess.append( err );
    CsErrLog::Instance()->mes( elFatal, mess );
//  ------------------------------------------
    return  0.;

  }


//===========================================================================
  bool CsRCPartPhotons::getThetaLikeMax( long double& likeMax,
//------------------------------------------------------------
					 double& theMax ) {


//- theta ring from max. Likelihood
//  -------------------------------
//- Paolo  -  September 2002
//    rev.    December  2002 
//    rev.    February  2003
//    rev.    October   2004
//    rev.    January   2005

//  WARNING : NOT likeLog compatible!
//  ---------------------------------

    CsRCExeKeys* key = CsRCExeKeys::Instance();
    static const bool tTrue = key->thetaLikeMaxTrue();

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      key->acknoMethod( "CsRCPartPhotons::getThetaLikeMax" );

      if( !tTrue ) {
        CsRCEventPartPhotons::Instance()->getBackgrParam04();   // 050208
      }
    }

    CsRCRecConst* cons = CsRCRecConst::Instance();

    const std::string backgrType0 = cons->backgrType();         // 050208
    if( !tTrue ) {
      //std::cout << "---b0 " << backgrType0 << std::endl;
      const std::string backgrTypeU = "04";                     // 050208
//@@--------------------------------------
      cons->setBackgrType( backgrTypeU );                       // 050208
      //const std::string bt = cons->backgrType();
      //std::cout << "bt " << bt << std::endl;
    }

    bool bRet = false;
    double delta, den, theMs, theCt, thePs, likeMs, likeCt, likePs;

//- Coarse search ( from checkLikelihood - 02/11/14 )
//  -------------------------------------------------
    theMax = 0.;
    double likeMaxW = -1.;
    double theMaxW = 0.;
    long double likeM = 0.;
    //double theLiMin = 5.;
    //double theLiMax = 55.;
    double theLiMin = 6.;
    double theLiMax = 56.;
    //double theLiMax = 60.;
    int nIpo = 25;
    double dThe = (theLiMax - theLiMin) / float( nIpo );
    long double likeV[nIpo];
    int kIpoV = 0;
    int nPhot = 0;
    for( int kIpo=0; kIpo<nIpo; kIpo++ ) {
      double theIpo = theLiMin + kIpo * dThe + dThe/2.;
      likeM = getLikelihood( theIpo, nPhot );
//            ------------------------------
      likeV[kIpo] = likeM;
      if( likeM < 0. ) continue;
      if( likeM > likeMaxW ) {
	likeMaxW = likeM;
	theMaxW = theIpo;
	kIpoV = kIpo;
      }
    }
    //cout << likeMax << "  " << theMax << "    " << likeMaxW << "  "
    //     << theMaxW << endl;

    likeMax = likeMaxW;
    theMax = theMaxW;
//  ----------------
    if( likeMax < 0. ) { bRet = false; goto pRet; }
    if( theMax <= 0. ) { bRet = false; goto pRet; }
    if( theMax >= theLiMax ) { 
      theMax = theLiMax;
      bRet = true; goto pRet;
    }


//- Fine search
//  -----------
    if(kIpoV <= 0 ) kIpoV = 1;
    if(kIpoV >= nIpo-1 ) kIpoV = nIpo - 2;
    delta = dThe;
    theMs = theMax - delta;
    theCt = theMax;
    thePs = theMax + delta;
    likeMs = likeV[kIpoV-1];
    if( likeMs < 0. ) { bRet = true; goto pRet; }
    likeCt = likeMax;
    likePs = likeV[kIpoV+1];
    if( likePs < 0. ) { bRet = true; goto pRet; }
    den = (likePs - 2.* likeCt + likeMs);
    if( den == 0. ) { bRet = true; goto pRet; }
    theMax = theCt - delta/2. * (likePs - likeMs) / den;
//  ---------------------------------------------------
    if( theMax <= 0. ) { bRet = false; goto pRet; }
    if( theMax < theMs ) theMax = theMs;
    if( theMax > thePs ) theMax = thePs;
    if( theMax >= theLiMax ) {
      likeMax = getLikelihood( theMax, nPhot );
      bRet = true; goto pRet;
    }

    delta /= 3.;
    //cout << "------ " << theMax << "  " << delta << endl;
    for( int kRep=0; kRep<3; kRep++ ) {
      int nPhot = 0;
      theMs = theMax - delta;
      theCt = theMax;
      thePs = theMax + delta;
      likeMs = getLikelihood( theMs, nPhot );
      if( likeMs < 0. ) break;
      likeCt = getLikelihood( theCt, nPhot );
      if( likeCt < 0. ) break;
      likePs = getLikelihood( thePs, nPhot );
      if( likePs < 0. ) break;
      //cout << likeMs << "  " << likeCt << "  " << likePs << endl;
      den = (likePs - 2.* likeCt + likeMs);
      if( den == 0. ) break;
      theMax = theCt - delta/2. * (likePs - likeMs) / den;
//    ---------------------------------------------------
      if( theMax <= 0. ) { bRet = false; goto pRet; }
      if( theMax < theMs ) theMax = theMs;
      if( theMax > thePs ) theMax = thePs;
      if( theMax >= theLiMax ) break;
      delta /= 3.;
    }

    likeMax = getLikelihood( theMax, nPhot );
    if( likeMax < 0. ) { bRet = false; goto pRet; }
    //cout << "=== " << likeMax << "  " << theMax << "  " << delta << endl;
    bRet = true;

    pRet : 
    if( !tTrue ) {
      cons->setBackgrType( backgrType0 );                   // 050208
      //const std::string b0 = cons->backgrType();
      //std::cout << " b0 " << b0 << std::endl;
    }

    return  bRet;

  }


//===========================================================================
  void CsRCPartPhotons::checkLikelihood( const double theReco ) {
//---------------------------------------------------------------


//- Paolo  -  July 2001
//    rev.    December 2002


    CsRCRecConst *cons = CsRCRecConst::Instance();
    int kOffHi = cons->kOffHi();

    static vector<CsHist1D*> vRC3080;
    static int kH3080 = 0;
    static int kHH = 0;
    static int kHHMx = 20;
    if( !CsRCHistos::Ref().bookHis() ||
	CsRCHistos::Ref().levelBk() < 2 ) kHHMx = 0;
    static CsHist1D* hRC3070;
    static double dTheRe = 4.;
    static int nIpoRe = 80;
    static int nPhoRing = 0;

    static int nPhoMin = 20;
//@@-----------------------

    int kTyLL = 2;
//@@-------------

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      CsRCExeKeys* key = CsRCExeKeys::Instance();
      key->acknoMethod( "CsRCPartPhotons::checkLikelihood" );

      hRC3070 = NULL;
      if( CsRCHistos::Ref().bookHis() &&
	  CsRCHistos::Ref().levelBk() >= 2 ) {
        int kHist = kOffHi + 3070;
	stringstream hN3070;
        hN3070 << kHist;
        string hTitle = "likelihood shape";
        CsHistograms::SetCurrentPath("/RICH");
        hRC3070 = 
	  new CsHist1D( hN3070.str(), hTitle, nIpoRe, -dTheRe, dTheRe );
//hh    ---------------------------------------------------------------
        CsHistograms::SetCurrentPath("/");
      }
    }

    if( kHH > kHHMx ) return;

    nPhoRing = lPhotons_.size();
    int iPartTy = 0;
    if( CsRCExeKeys::Instance()->MCarloEvent() ) {
      iPartTy = pPart_->iPartT();
      if( iPartTy > 200 ) iPartTy -= 200;
      if( iPartTy > 30 ) iPartTy = 0;   //   wrong myFile?
//@@--------------------------------
    }
    if( nPhoRing > nPhoMin ) {
      //if( iPartTy ==  8 || iPartTy ==  9 ) {     // pions only
      //if( iPartTy == 11 || iPartTy == 12 ) {     // kaons only
      //if( iPartTy == 14 || iPartTy == 15 ) {     // protons only
      double theIpo = 0.;
      int nPhot = 0;
      double theMin = 0.;
      double theMax = 80.;
      int nIpo = 321;                            // corrected 020716 !
      double dThe = (theMax - theMin)/float( nIpo );

      if( kHH < kHHMx ) {
	int kHist = kOffHi + 3080 + kHH;
        stringstream hN3080;
	hN3080 << kHist;
	string hTitle = "likelihood shape";
	CsHistograms::SetCurrentPath("/RICH");
        vRC3080.push_back( new CsHist1D( hN3080.str(), hTitle,
//hh    ------------------------------------------------------
					 nIpo, theMin, theMax ) );
	CsHistograms::SetCurrentPath("/");
      }
      //cout << "checkL---------------  " << theReco << endl;

      double pLike0 = -100000.;
      double theIpo0 = -1.;
      double pLike[nIpo];
      for( int kIpo=0; kIpo<nIpo; kIpo++ ) {
	//cout << kIpo << "  ";
        theIpo = theMin + kIpo * dThe + dThe/2.;
        pLike[kIpo] = getLikelihood( theIpo, nPhot );
//                    ------------------------------
	//cout << pLike[kIpo] << endl;
        if( pLike[kIpo] > pLike0 ) {
	  pLike0 = pLike[kIpo];
	  theIpo0 = theIpo;
	}
      }

      if( kHH < kHHMx ) {
	//int nIpo4 = nIpo - 4;
	int nIpo4 = nIpo - 6;
        for( int kIpo=0; kIpo<nIpo4; kIpo++ ) {
          theIpo = theMin + kIpo * dThe + dThe/2.;
	  double xh = theIpo;
	  //double wh = pLike[kIpo] / pLike0;
          double wh = pLike[kIpo];
	  if( vRC3080[kHH] ) vRC3080[kHH]->Fill( xh, wh );
//hh                         ----------------------------
	}
        //int iIpo = 4;
        int iIpo = 6;
	float momPart = pPart_->mom();
	double xh = theMin + (nIpo-iIpo) * dThe + dThe/2.;
	double wh = momPart;
	if( vRC3080[kHH] ) vRC3080[kHH]->Fill( xh, wh );
//hh                       ----------------------------
        //iIpo = 3;
        iIpo = 5;
        for( int kPaTy=kTyLL; kPaTy<=14; kPaTy+=3 ) {
//      -------------------------------------------
	  theIpo = thetaIpo_[kPaTy];
	  double xh = theMin + (nIpo-iIpo) * dThe + dThe/2.;
          double wh = theIpo;
          if( vRC3080[kHH] ) vRC3080[kHH]->Fill( xh, wh );
//hh                         ----------------------------
	  //cout << kHH << "  " << momPart << endl;
          iIpo--;
	}
      }
      kHH++;

      theMin = theReco - dTheRe;
      theMax = theReco + dTheRe;
      dThe = (theMax - theMin)/float( nIpoRe );
      for( int kIpo=0; kIpo<=nIpoRe; kIpo++ ) {
        theIpo = theMin + kIpo * dThe;
        double pLike = getLikelihood( theIpo, nPhot );
//                     ------------------------------
	//pLike /= pLike0;
	double xh = - dTheRe + kIpo * dThe;       // !
	xh -= 0.05;
        double wh = pLike;
        if( hRC3070 ) hRC3070->Fill( xh, wh );
//hh                  -----------------------
      }
    }

  }


//===========================================================================
  void CsRCPartPhotons::checkLikeRatios( const CsRCPhoton* pPhot,
//---------------------------------------------------------------
			   	         const double signal,
					 const double backgr ) {


//- compute some values for Likelihood
//  ----------------------------------
//- Paolo  -  October 2004


    CsRCExeKeys* key = CsRCExeKeys::Instance();
    CsRCRecConst* cons = CsRCRecConst::Instance();

    double xh, yh, wh;
    static CsHist1D* hRC1301 = NULL;
    static CsHist1D* hRC1302 = NULL;
    static CsHist1D* hRC1303 = NULL;
    static CsHist1D* hRC1305 = NULL;
    static CsHist1D* hRC1309 = NULL;
    static CsHist1D* hRC1310 = NULL;

    static const long double epsi = 1.e-320;
    static const long double lepsi = log( epsi );
//@@--------------------------------------------

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      key->acknoMethod( "CsRCPartPhotons::checkLikeRatios" );

      if( CsRCHistos::Ref().bookHis() ) {
	int kOffHi = cons->kOffHi();
	stringstream hN1301;
	stringstream hN1302;
	stringstream hN1303;
	stringstream hN1305;
	stringstream hN1309;
	stringstream hN1310;
	string hTitle = "   ";
	hN1301 << kOffHi + 1301;
	hRC1301 = new CsHist1D( hN1301.str(), hTitle, 100, -140., 10. );
	hN1302 << kOffHi + 1302;
	hRC1302 = new CsHist1D( hN1302.str(), hTitle, 100, -140., 10. );
	hN1303 << kOffHi + 1303;
	hRC1303 = new CsHist1D( hN1303.str(), hTitle, 100, -140., 10. );
	hN1305 << kOffHi + 1305;
	hRC1305 = new CsHist1D( hN1305.str(), hTitle, 100, -140., 10. );
	hN1309 << kOffHi + 1309;
	hRC1309 = new CsHist1D( hN1309.str(), hTitle, 100, -140., 10. );
	hN1310 << kOffHi + 1310;
	hRC1310 = new CsHist1D( hN1310.str(), hTitle, 100, -140., 10. );
      }
    }

    static std::vector<double> sigVec;
    static std::vector<double> bakVec;

    static CsRCPhoton* pPhotF = NULL;
    pPhotF = lPhotons_.front();
    static CsRCPhoton* pPhotB = NULL;
    pPhotB = lPhotons_.back();

    static bool make3 = false;
    static bool makek = false;
    static int count3 = 0;
    static int countk = 0;
    if( pPhot == pPhotF ) {
      sigVec.clear();
      bakVec.clear();
      count3++;
      make3 = false;
      if( count3 == 100 ) {
//@@---------------------
	count3 = 0;
	make3 = true;
      }
      countk++;
      makek = false;
      if( countk == 100 ) {
//@@---------------------
	countk = 0;
	makek = true;
      }
    }

    sigVec.push_back( signal );
    bakVec.push_back( backgr );

    if( pPhot == pPhotB ) {
      int nPhot = sigVec.size();

//--- Prod 1
      for( int kf=0; kf<nPhot; kf++ ) {
	double ratio1 = 0.;
	if( bakVec[kf] > 0. ) {
	  ratio1 = sigVec[kf]/bakVec[kf];
	  //std::cout << ratio1 << "  " 
	  //<< sigVec[kf] << "  " << bakVec[kf] << std::endl;
	  xh = 0.;
          if( ratio1 > 0. ) xh = log( ratio1 );
          if( hRC1301 ) hRC1301->Fill( xh );
//hh                    -------------------
	}
      }

//--- Prod 2
      for( int kf=0; kf<nPhot; kf++ ) {
	for( int jf=0; jf<nPhot; jf++ ) {
	  if( jf <= kf ) continue;
	  double ratio2 = 0.;
	  if( bakVec[kf] > 0.  &&  bakVec[jf] > 0. ) {
	    ratio2 = sigVec[kf]/bakVec[kf] * sigVec[jf]/bakVec[jf];
	    xh = 0.;
	    if( ratio2 > 0. ) xh = log ( ratio2 ) / 2.;
	    if( hRC1302 ) hRC1302->Fill( xh );
//hh                       ------------------
	  }
	}
      }

//--- Prod 3
      if( make3 ) {
        for( int kf=0; kf<nPhot; kf++ ) {
	  for( int jf=0; jf<nPhot; jf++ ) {
	    if( jf <= kf ) continue;
	    for( int lf=0; lf<nPhot; lf++ ) {
	      if( lf <= jf ) continue;
	      double ratio3 = 0.;
	      if( bakVec[kf] > 0. && bakVec[jf] > 0. && bakVec[lf] > 0. ) {
	        ratio3 = sigVec[kf]/bakVec[kf] * sigVec[jf]/bakVec[jf] *
		         sigVec[lf]/bakVec[lf];
                xh = 0.;
	        if( ratio3 > 0. ) xh = log ( ratio3 ) / 3.;
	        if( hRC1303 ) hRC1303->Fill( xh );
//hh                          -------------------
	      }
	    }
	  }
        }
      }

/*
//------------------------------------------------------
      double sigPro = 0.;
      double bakPro = 0.;
      bool doPro = false;
      static int nUse = 50;
      nUse = nPhot;
      //----------
      if( nPhot <= nUse ) doPro = true;
      static int nPro = 4;
      if( nPro > nUse ) doPro = false;
      if( doPro ) {
	int ix[nUse];
	ix[0] = 0;
	for( int kx=1; kx<nPro; kx++ ) ix[kx] = ix[kx-1] + 1;
	for( int kLoop=nPro-1; kLoop>=0; kLoop-- ) {
	  if( kLoop < nPro-1 ) ix[kLoop]++;
	  for( int kx=kLoop+1; kx<nPro; kx++ ) ix[kx] = ix[kx-1] + 1;
	  for( int kLp=nPro-1; kLp>=kLoop; kLp-- ) {
	    if( kLp < nPro-1 ) ix[kLp]++;
	    //for( int kx=kLp+1; kx<nPro; kx++ ) ix[kx] = ix[kx-1] + 1;
	    for( int kl=0; kl<nUse; kl++ ) {
	      if( kl > 0 ) ix[kLp]++;
	      if( ix[kLp] < nUse ) {
		//std::cout << nUse << "  " << kLoop << "  " << kLp
		//	  << "  " << kl << std::endl;
                //std::cout << ix[0] << "  " << ix[1] << "  " << ix[2]
                //	  << "  " << ix[3] << std::endl;
		sigPro = 0.;
		bakPro = 0.;
		for( int kx=0; kx<nPro; kx++ ) {
	          sigPro *= sigVec[ix[kx]];
	          bakPro *= bakVec[ix[kx]];
	        }
	        //
	      } else  break;
	    }
	  }
	}
      }
//------------------------------------------------------
*/

      long double sigPro = 0.;
      long double bakPro = 0.;

//--- Prod k
      bool doPro = false;
      static int nUseMx = 70;
      //nUseMx = nPhot;
      //nPhot = nUseMx;
      //--------------
      if( nPhot <= nUseMx ) doPro = true;
      int nUse = nPhot;
      if( !makek ) doPro = false;

      static int nPro = 3;
      //static int nPro = nPhot - 1;
//@@---------------------
      if( nPro > nUse ) doPro = false;
      int countPro = 0;
      if( doPro ) {
	int ix[nUse];
	ix[0] = 0;
	for( int kx=1; kx<nPro; kx++ ) ix[kx] = ix[kx-1] + 1;
	bool pro = true;
	while ( pro ) {
	  int kil = nPro - 1;
	  bool good = true;
	  do {
	    good = ix[kil] < nUse;
	    if( good ) {
	      //std::cout << ix[0] << "  " << ix[1] << "  " << ix[2]
	      //	  << "  " << ix[3] << std::endl;
	      sigPro = 0.;
	      bakPro = 0.;
	      for( int kx=0; kx<nPro; kx++ ) {
		double sigw = sigVec[ix[kx]];
		sigw = log( sigw );
		sigw /= float( nPro );
		if( sigw <= lepsi ) sigw = lepsi;
		sigPro += sigw;
		if( sigPro <= lepsi ) sigPro = lepsi;
		double bakw = bakVec[ix[kx]];
		bakw = log( bakw );
		bakw /= float( nPro );
		if( bakw <= lepsi ) bakw = lepsi;
		bakPro += bakw;
		if( bakPro <= lepsi ) bakPro = lepsi;
	      }
	      //std::cout << sigPro << "  " << bakPro << std::endl;
	      if( sigPro != 0.  &&  bakPro != 0. ) {
		double ratiok = sigPro - bakPro;
		//ratiok /= float( nPro );
		xh = ratiok;
		if( hRC1305 ) hRC1305->Fill( xh );
//hh                          -------------------
	      }
	      countPro++;
	    }
	    ix[kil]++;
	  } while ( good );
	  good = true;
	  do {
	    kil--;
	    if( kil < 0 ) { pro = false; break; }
	    ix[kil]++;
	    if( ix[kil] >= nUse ) continue;
	    for( int kx=kil+1; kx<nPro; kx++ ) {
	      ix[kx] = ix[kx-1] + 1;
	      if( ix[kx] >= nUse ) { good = true; break; }
	      good = false;
	    }
	    //std::cout << "1 " << ix[0] << "  " << ix[1] << "  " << ix[2]
	    //	  << "  " << ix[3] << std::endl;
	  } while ( good );
	}
	//std::cout << countPro <<std::endl;
      }

//--- Prod N-1
      for( int kf=0; kf<nPhot; kf++ ) {
	sigPro = 0.;
	bakPro = 0.;
	for( int jf=0; jf<nPhot; jf++ ) {
	  if( jf == kf ) continue;
	  double sigw = sigVec[jf];
	  sigw = log( sigw );
	  sigw /= float( nPhot-1 );
	  if( sigw <= lepsi ) sigw = lepsi;
	  sigPro += sigw;
	  if( sigPro <= lepsi ) sigPro = lepsi;
	  double bakw = bakVec[jf];
	  bakw = log( bakw );
	  bakw /= float( nPhot-1 );
	  if( bakw <= lepsi ) bakw = lepsi; 
	  bakPro += bakw;
	  if( bakPro <= lepsi ) bakPro = lepsi;
	}
        //std::cout << sigPro << "  " << bakPro << std::endl;
	if( sigPro != 0.  &&  bakPro != 0. ) {
	  double ration = sigPro - bakPro;
	  //ration /= float( nPhot-1 );
	  //std::cout << ration << std::endl;
	  xh = ration;
	  if( hRC1309 ) hRC1309->Fill( xh );
//hh                    -------------------
	}
      }

//--- Prod N
      sigPro = 0.;
      bakPro = 0.;
      for( int kf=0; kf<nPhot; kf++ ) {
	double sigw = sigVec[kf];
	sigw = log( sigw );
	sigw /= float( nPhot );
	if( sigw <= lepsi ) sigw = lepsi;
	sigPro += sigw;
	if( sigPro <= lepsi ) sigPro = lepsi;
	double bakw = bakVec[kf];
	bakw = log( bakw );
	bakw /= float( nPhot );
	if( bakw <= lepsi ) bakw = lepsi; 
	bakPro += bakw;
	if( bakPro <= lepsi ) bakPro = lepsi;
      }
      //std::cout << sigPro << "  " << bakPro << std::endl;
      double ration = sigPro - bakPro;
      //ration /= float( nPhot );
      xh = ration;
      if( hRC1310 ) hRC1310->Fill( xh );
//hh                -------------------

      //exit(0);
    }


    return;
  }


//===========================================================================
  double CsRCPartPhotons::getQsquare( const double theIpo, int &nPhoRing ) {
//--------------------------------------------------------------------------


//- compute  qSquare value for a given hypothesis
//  ---------------------------------------------
//- Paolo  -  October 2004


    CsRCExeKeys *key = CsRCExeKeys::Instance();

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      key->acknoMethod( "CsRCPartPhotons::getQsquare" );
    }


    list<CsRCPhoton*> lPhotons = lPhotons_;
    list<CsRCPhoton*>::const_iterator ih;

    long double qSq = 0.;
    int kPhot = 0;
    for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
      if( !(*ih)->flag() ) continue;

      double theIpow = (*ih)->getThetaIpo( theIpo );

      double thePhot = (*ih)->the();
      double sigPhot = (*ih)->sigmaPhoPid( this );
//                     --------------------------
      double qSw = (thePhot - theIpow) / sigPhot;
      qSq += qSw*qSw;

      kPhot++;
    }   /* end for on Photons */

    nPhoRing = kPhot;
    double qSquare = 0.;
    if( nPhoRing > 0) {
      qSquare = qSq / float( nPhoRing );
    }
    //cout << theIpow << "  " << nPhoRing << "  " << qSq
    //     << "  " << qSquare << endl;

    return  qSquare;

  }


//===========================================================================
  void CsRCPartPhotons::compThetaIpo( CsRCParticle* part ) {
//----------------------------------------------------------


//- compute the Ch angle for the given mass hypothesis
//  --------------------------------------------------
//  for UV and VS ref. Index
//- Paolo  -  August 2005

    CsRCRecConst* cons = CsRCRecConst::Instance();
    bool UpRICHJob = CsRichOne::Instance()->UpRICHJob();
    float CFRefInd = cons->CFRefInd();
    float CFRefIndVS = 1.;
    if( UpRICHJob ) CFRefIndVS = cons->CFRefIndVS();
    double* massPartv = cons->massPartv();

    for( int kPaTy=2; kPaTy<=14; kPaTy+=3 ) {

      thetaIpo_[kPaTy] = -1.;
      thetaIpo_[kPaTy+1] = -1.;
      thetaIpoUV_[kPaTy] = -1.;
      thetaIpoUV_[kPaTy+1] = -1.;
      float massIpo = massPartv[kPaTy];
      float momPart = part->mom();
      double betaIpo = momPart / sqrt( momPart*momPart + massIpo*massIpo );
      double cosTheW = 1./ (betaIpo * CFRefInd);
      if( cosTheW <= 1.) {
	thetaIpo_[kPaTy] = acos( cosTheW ) * 1000.;
	thetaIpo_[kPaTy+1] = thetaIpo_[kPaTy];
	thetaIpoUV_[kPaTy] = thetaIpo_[kPaTy];
	thetaIpoUV_[kPaTy+1] = thetaIpo_[kPaTy];
      }
      thetaIpoVS_[kPaTy] = -1.;
      thetaIpoVS_[kPaTy+1] = -1.;
      if( UpRICHJob ) {
        cosTheW = 1./ (betaIpo * CFRefIndVS);
        if( cosTheW <= 1.) {
	  thetaIpoVS_[kPaTy] = acos( cosTheW ) * 1000.;
	  thetaIpoVS_[kPaTy+1] = thetaIpoVS_[kPaTy];
        }
      }
    }

  }


//===========================================================================
  void CsRCPartPhotons::setPhoTheNorm() {
//---------------------------------------


//- Paolo  -  August 2005


    list<CsRCPhoton*>::iterator ih;
    for( ih=lPhotons_.begin(); ih!=lPhotons_.end(); ih++ ) {

//--- change photon the() angle definition (for RING analysis) :
//    ----------------------------------------------------------
      (*ih)->setTheToNorm();
//@@-----------------------

    }

  }


//===========================================================================
  double CsRCPartPhotons::howPMT() {
//----------------------------------


//- Paolo  -  May 2006


    int nPhoPMT = 0;
    list<CsRCPhoton*>::iterator ih;
    for( ih=lPhotons_.begin(); ih!=lPhotons_.end(); ih++ ) {
      if( (*ih)->isPMT() ) nPhoPMT++;
    }
    double ratio = 0.;
    if( lPhotons_.size() > 0 ) ratio = float(nPhoPMT)/float(lPhotons_.size());

    return  ratio;

  }


//===========================================================================
  void CsRCPartPhotons::checkLikeReso() {
//---------------------------------------


//- Paolo  -  March 2006
//  rev.      March 2007


    CsRCExeKeys* key = CsRCExeKeys::Instance();
    CsRCRecConst *cons = CsRCRecConst::Instance();
    int kOffHi = cons->kOffHi();

    static std::vector<CsHist2D*> vRC3270;
    static std::vector<CsHist1D*> vRC8000;
    static std::vector<CsHist1D*> vRC8020;
    static std::vector<CsHist1D*> vRC8050;
    static int nFGood = 0;
    static int nHGood = 0;
    static int nFFail = 0;
//@@---------------------

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      key->acknoMethod( "CsRCPartPhotons::checkLikeReso" );

      for( int kh=0; kh<7; kh++ ) vRC3270.push_back( NULL );
      if( CsRCHistos::Ref().bookHis() &&
	  CsRCHistos::Ref().levelBk() >= 2 ) {
	CsHistograms::SetCurrentPath("/RICH");
	vRC3270.clear();
        string hTitle = "likelihood reso";
        int kHist = 0;
	stringstream hN3270;
	stringstream hN3271;
	stringstream hN3272;
	stringstream hN3273;
	stringstream hN3274;
	stringstream hN3275;
	stringstream hN3276;
	stringstream hN3277;
        kHist = kOffHi + 3270;
        hN3270 << kHist;
        vRC3270.push_back( new CsHist2D( hN3270.str(), hTitle,
					 //120, 0., 60., 100, 0., 5.) );
					 //120, 0., 60., 200, 0., 20.) );
					 60, 0., 60., 200, 0., 5.) );
	                                 //60, 0., 60., 200, 1., 3.) );
        kHist = kOffHi + 3271;
        hN3271 << kHist;
        vRC3270.push_back( new CsHist2D( hN3271.str(), hTitle,
					 60, 0., 60., 100, 0., 3.) );
        kHist = kOffHi + 3272;
        hN3272 << kHist;
        vRC3270.push_back( new CsHist2D( hN3272.str(), hTitle,
					 60, 0., 60., 100, -5., 5.) );
        kHist = kOffHi + 3273;
        hN3273 << kHist;
        vRC3270.push_back( new CsHist2D( hN3273.str(), hTitle,
					 //120, 0., 60., 100, 0., 3.) );
					 //120, 0., 60., 100, -5., 5.) );
					 60, 0., 60., 100, -1., 1.) );
        kHist = kOffHi + 3274;
        hN3274 << kHist;
        vRC3270.push_back( new CsHist2D( hN3274.str(), hTitle,
					 //120, 0., 60., 100, -5., 5.) );
					 60, 0., 60., 100, -0.25, 0.25) );
        kHist = kOffHi + 3275;
        hN3275 << kHist;
        vRC3270.push_back( new CsHist2D( hN3275.str(), hTitle,
					 60, 0., 60., 100, 0., 1.) );
        kHist = kOffHi + 3276;
        hN3276 << kHist;
        vRC3270.push_back( new CsHist2D( hN3276.str(), hTitle,
					 60, 0., 60., 100, 0., 1.) );
        kHist = kOffHi + 3277;
        hN3277 << kHist;
        vRC3270.push_back( new CsHist2D( hN3277.str(), hTitle,
					 60, 0., 60., 100, -5., 5.) );
        CsHistograms::SetCurrentPath("/");
      }
    }

    double xh, yh, wh;


    int iPartTy = 0;
    if( key->MCarloEvent() ) {
      iPartTy = pPart_->iPartT();
      if( iPartTy > 200 ) iPartTy -= 200;
      if( iPartTy > 30 ) iPartTy = 0;   //   wrong myFile?
//@@--------------------------------
    }


    float momPart = pPart_->mom();

//- Particle momentum selection :
//  -----------------------------
    if( momPart < 2.5 ) return;
//@@--------------------------

    double phiMass  = CsRichOne::Instance()->phiMass();

//- PHI MASS SELECTION! (for PHI Sample ONLY!!!)
//---------------------
    if( key->myFileType() == 3 ) {
      if( fabs( phiMass - 1.020) > 0.012 ) return;
//@@---------------------------------------------
    }

    double thePart = acos( pPart_->nd() );

//- Particle angle selection :
//  --------------------------
    //if( thePart < 50./1000.) return;
//@@-------------------------------

    double likeDV = cons->likeDefVa();
//@@---------------------------------


//- Background Likelihood :
//  -----------------------
    int nPhot = 0;
    double pLikeBk = likeDV;
    if( partProbsSet_ ) pLikeBk = probaLKBgAll_;
    else  pLikeBk = getLikelihood( 0., nPhot );
//                  --------------------------

//- Fine step Likelihood and ChiSquare :
//  ------------------------------------
    double theMin = 0.;
    double theMax = 0.;
    int nIpo = 0;
    if( !likeProSet_ ) if( !getLikeProfile() ) return;
//                          ----------------
    theMin = theLikeProMn_;
    theMax = theLikeProMx_;
    nIpo = nLikePro_;
    double dThe = (theMax - theMin)/float( nIpo-1 );
    double pLikek[nIpo];
    double qSquak[nIpo];
    double thek[nIpo];
    for( int kIpo=0; kIpo<nIpo; kIpo++ ) {
      thek[kIpo] = theLikePro_[kIpo];
      pLikek[kIpo] = pLikePro_[kIpo];
      qSquak[kIpo] = getQsquare( thek[kIpo], nPhot );
//                   -------------------------------
    }

    double pLikeMx = -100000.;
    double theIpoMx = -1.;
    int kIpoMx = -1;
    for( int kIpo=0; kIpo<nIpo; kIpo++ ) {
      if( pLikek[kIpo] > pLikeMx ) {
        pLikeMx = pLikek[kIpo];
        theIpoMx = thek[kIpo];
	kIpoMx = kIpo;
      }
      //std::cout << setprecision( 7 );
      //std::cout << kIpo << "  " << thek[kIpo] << "  " << pLikek[kIpo] 
      //<< std::endl;
      //std::cout << kIpo << "  " << thek[kIpo] << "  " << qSquak[kIpo] 
      //<< std::endl;
    }
    //std::cout << kIpoMx << "  " << theIpoMx << "  " << pLikeMx << std::endl;

    //getRelMinima( nIpo, pLikek );      //   pro-memoria
//@@----------------------------

    if( theIpoMx <  5. ) return;
    if( theIpoMx > 60. ) return;
//@@---------------------------


//- Theta(likemax) and Sigma calculation :
//  --------------------------------------
    double sigCal = 0.;
    double meaCal = 0.;

    bool exeCal = true;
//@@------------------
    if( exeCal ) {
      int npc = 17;
      int nof = 0;
      int nst = 1;
      int lpc = npc/2;
      double pLikec[npc];
      double thec[npc];
      int nc = 0;
      for( int kpc=0; kpc<npc; kpc++ ) {
	int jpc = kpc - lpc;
	int kpcc = kIpoMx - nof + jpc*nst;
	if( kpcc < 0 ) continue;
	if( kpcc >= nIpo ) continue;
	pLikec[kpc] = pLikek[kpcc] - pLikeBk;
//    ---------------------------------------
        if( cons->likeLog() ) {
        } else {
	  if( pLikec[kpc] <= 0.)  pLikec[kpc] = 0.0001;
        }
        thec[kpc] = thek[kpcc];
        nc++;
      }
      lpc = nc/2;
      double ach = dThe;
    nc = 0;
    for( int kc=4; kc<=lpc; kc++ ) {
      double loLikeC = 0.;
      double loLikeL = 0.;
      double loLikeR = 0.;
      if( cons->likeLog() ) {
	loLikeC = pLikec[lpc];
	loLikeL = pLikec[lpc-kc];
	loLikeR = pLikec[lpc+kc];
      } else {
	loLikeC = log( pLikec[lpc] );
	loLikeL = log( pLikec[lpc-kc] );
	loLikeR = log( pLikec[lpc+kc] );
      }
      //std::cout << setprecision(6);
      //std::cout << pLikec[lpc-kc] << "  " << pLikek[lpc] << "  "
      //          << pLikeMx << "  " << pLikek[lpc+kc] << std::endl;
      //std::cout << setprecision(2);
      //std::cout << loLikeL << "  " << loLikeC << "  " << log(pLikeMx)
      //          << "  " << loLikeR << std::endl;
      double ddt = ach * kc*nst;
      double ddlo = 2.*loLikeC - loLikeL - loLikeR;
      if( ddlo <= 0.) continue; 
      sigCal += ddt / sqrt( ddlo );
      meaCal += thec[lpc] + ddt* (loLikeR - loLikeL)/ (2.* ddlo);
      //std::cout << theIpoMx << "  " << meaCal << " " << sigCal
      //          << std::endl;
      nc++;
    }
    if( nc > 0 ) {
      sigCal /= nc;
      meaCal /= nc;
    }
    //std::cout << theIpoMx << "  " << meaCal << " " << sigCal << std::endl;
    //std::cout << std::endl;
    //std::cout << setprecision( 6 );
    //std::cout << CsRCRecConst::Instance()->CFRefInd() << std::endl;
    }


//- Particle ID :
//  -------------
    int kTyLL = 2;
//@@-------------
    double theIpo[5];
    double pLikem[5];
    int kIpo = 0;
    for( int kPaTy=kTyLL; kPaTy<=14; kPaTy+=3 ) {
//  -------------------------------------------
      theIpo[kIpo] = thetaIpo( kPaTy );
      if( partProbsSet_ )  pLikem[kIpo] = probaLKAll_[kPaTy];
      else  pLikem[kIpo] = getLikelihood( theIpo[kIpo], nPhot );
//                         ------------------------------------
      kIpo++;
    }
    //int selIpoMx = 2;      //      Pions
    int selIpoMx = 3;      //      Kaons
//@@----------------
    int nMassIpo = kIpo;
    double pLikemMx = -100000.;
    double theIpomMx = -1.;
    int kIpomMx = -1;
    for( int kIpo=0; kIpo<5; kIpo++ ) {
      if( pLikem[kIpo] == likeDV ) continue;
      if( pLikem[kIpo] > pLikemMx ) {
	pLikemMx = pLikem[kIpo];
	theIpomMx = theIpo[kIpo];
	kIpomMx = kIpo;
      }
    }
    //for( int k=0; k<5; k++ ) std::cout << pLikem[k] << "  " << theIpo[k]
    //                               << "  ";
    //if( momPart > 10.  &&  momPart < 18. ) {
    //for( int k=0; k<5; k++ ) std::cout << pLikem[k] << "  ";
    //std::cout << pLikeBk << "  " << kIpomMx << " " << pLikemMx
    //	  << "  " << momPart << std::endl;
    //}

//- Theta(likemax) and Sigma FIT :
//  ------------------------------
    double norm   = -1.;
    double meaFit = -1.;
    double sigFit = -1.;
    double po0    = -1.;
    double po1    = -1.;
    double po2    = -1.;

    bool exeFit = true;
//@@------------------
    if( kIpomMx == selIpoMx  &&  exeFit ) {

      //int nPf = 24;
      int nPf = 12;
//@@--------------
      double thef[nPf];
      double pLikef[nPf];
      double eLikef[nPf];
      int kf = 0;
      for( int kk=0; kk<nPf; kk++ ) {
        int ku = kIpoMx - nPf/2 + kk;
        if( ku < 0 ) continue;
        if( ku >= nIpo ) continue;
        thef[kf] = thek[ku];
	//std::cout << log(pLikek[ku]) << "  " << log(pLikeBk) << std::endl;
	//pLikek[ku] = pow( pLikek[ku], nPhotAll() );   // !!!
	//pLikeBk    = pow( pLikeBk, nPhotAll() );      // !!!
	//std::cout << log(pLikek[ku]) << "  " << log(pLikeBk) << "  "
	//  	    << normSg() << "  " << nPhotAll() << std::endl;
        pLikef[kf] = pLikek[ku] - pLikeBk;
	//std::cout << thef[kf] << "  " << pLikef[kf]
        //          << "  " << log(pLikef[kf]) << std::endl;
//@@-------------------------------------
	if( cons->likeLog() ) {
	} else {
          if( pLikef[kf] < 0. ) pLikef[kf] = 0.;
	}
        eLikef[kf] = 0.01* pLikek[ku];
        if( cons->likeLog() ) eLikef[kf] = fabs( eLikef[kf] );
//@@---------------------------------------------------------
        //std::cout << kf << "  " << thef[kf] << "  " << pLikef[kf] 
	//          << std::end l;
        kf++;
      }

      int nPoint = kf;
      int nParam = 10;
      double param[nParam];
      double toll[nParam];
      int iPaFit[nParam];
      for( int kp=0; kp<nParam; kp++ ) {
        param[kp] = 0.;
        toll[kp] = 0.;
        iPaFit[kp] = 1;
      }
      if( cons->likeLog() ) {
	param[6] = 0.5*theIpoMx*theIpoMx;
	param[7] = -theIpoMx;
	param[8] = 0.5;
	toll[6] = 0.01;
	toll[7] = 0.01;
	toll[8] = 0.01;
	iPaFit[6] = 0;
        iPaFit[7] = 0;
	iPaFit[8] = 0;
      } else {
	param[0] = pLikeMx - pLikeBk;
	param[1] = theIpoMx;
	param[2] = 2.;
        //param[6] = pLikeBk;
	param[6] = 0.;
	param[7] = 0.;
	toll[0] = 0.01;
	toll[1] = 0.0001;
	toll[2] = 0.0001;
	toll[6] = 0.01;
	toll[7] = 0.01;
	iPaFit[0] = 0;
	iPaFit[1] = 0;
	iPaFit[2] = 0;
	if( pLikeBk > 0.) {
	  iPaFit[6] = 0;
	  iPaFit[7] = 0;
	}
      }


      CsRCGauPolFit oLikeFt( nPoint, thef, pLikef, eLikef,
//    ----------------------------------------------------
		  	     nParam, param, toll, iPaFit );
      if( oLikeFt.doChiFit() );
//        ------------------

      if( oLikeFt.flag() ) {
        nFGood++;
	if( cons->likeLog() ) {
	  po0 = oLikeFt.para()[6];
	  po1 = oLikeFt.para()[7];
	  po2 = oLikeFt.para()[8];
	  meaFit = - po1 / (2.*po2);
	  sigFit = 1./ sqrt( 2.*fabs( po2 ) );
	} else {
	  norm   = oLikeFt.para()[0];
	  meaFit = oLikeFt.para()[1];
	  sigFit = oLikeFt.para()[2];
	  po0    = oLikeFt.para()[6];
	  po1    = oLikeFt.para()[7];
	}
	//std::cout << "LikeReso fit " << norm << "  " << meaFit << "  ";
	//std::cout << sigFit << "  " << po0 << "  " << po1 << std::endl;
      } else {
	bool print = false;
//@@----------------------
	if( print ) {
          std::cout << "  checkLikeReso : fit failure!  " << nFFail << "  ";
          std::cout << oLikeFt.nIter() << "  ";
          std::cout << theIpoMx << "  " << pLikeMx << " " << pLikeBk
		    << std::endl;
	}
        nFFail++;

	if( !histLikeProfile( vRC8000, 8000, 20 ) ) {}
//           ------------------------------------

      }
      //std::cout << setprecision( 6 );
      //std::cout << norm << "  " << meaFit << "  " << sigFit << "  "
      //          << po0 << "  " << po1 << "  "
      //  	  << theIpoMx << std::endl;
      //std::cout << po0 << "  " << po1 << "  " << po2 << std::endl;
      //std::cout << setprecision( 2 );
      //std::cout << nFGood << "  " << nFFail << std::endl;
    }
    //std::cout << theIpoMx << "  " << meaCal << " " << sigCal << "  ";
    //std::cout << theIpoMx << "  " << meaFit << "  " << sigFit << std::endl;

    double theWgt = getThetaWgav();
    //std::cout << theIpoMx << "  " << meaFit << "  " << theWgt << std::endl;

//- HISTOGRAMMING :
//  ---------------

    double dTheIpo = fabs(theIpo[kIpomMx]-theIpo[kIpomMx-1]);
    if( dTheIpo == 0. ) return;

    int kHH = 0;
    bool bHH = kHH < int( vRC3270.size() );
    xh = momPart;
    ///yh = sigFit;
    yh = sigFit / sqrt( nPhotAll() );
    if( bHH  &&  vRC3270[kHH] ) vRC3270[kHH]->Fill( xh, yh );
//hh                            ----------------------------

    //double sigmaLikeMax = 0.65;
    double sigmaLikeMax = 0.51;
//@@--------------------------
    double sigmaDTheta = 29.*exp( -0.31*momPart ) + 0.61;
//@@----------------------------------------------------
    //std::cout << setprecision( 6 );
    //std::cout << pLikemMx << "  " << pLikeBk << "  ";
    //std::cout << kIpomMx << std::endl;
    if( pLikemMx > pLikeBk ) {
      ///if( kIpomMx == selIpoMx ) {
	//std::cout << pLikemMx << "  " << kIpomMx << "  " << std::endl;

        kHH = 1;
	bHH = kHH < int( vRC3270.size() );
        xh = momPart;
        //yh = sigFit / fabs(theIpo[kIpomMx]-theIpo[kIpomMx-1]);
	yh = sigFit / dTheIpo;
        if( bHH  &&  vRC3270[kHH] ) vRC3270[kHH]->Fill( xh, yh );
//hh                                ----------------------------
        kHH = 2;
	bHH = kHH < int( vRC3270.size() );
        xh = momPart;
        //yh = theIpoMx - theIpo[kIpomMx];
        //yh = theWgt - theIpo[kIpomMx];
	if( !thetaLikeSet_ ) if( !getThetaLikeMax() ) return;
//      ----------------------------------------------------
	yh = thetaLike_ - theIpo[kIpomMx];
        if( bHH  &&  vRC3270[kHH] ) vRC3270[kHH]->Fill( xh, yh );
//hh                                ----------------------------
	kHH = 5;
	bHH = kHH < int( vRC3270.size() );
        xh = momPart;
        //yh = sigmaLikeMax / fabs(theIpo[kIpomMx]-theIpo[kIpomMx-1]);
	yh = sigmaLikeMax / dTheIpo;
        if( bHH  &&  vRC3270[kHH] ) vRC3270[kHH]->Fill( xh, yh );
//hh                                ----------------------------
	kHH = 6;
	bHH = kHH < int( vRC3270.size() );
        xh = momPart;
        //yh = sigmaDTheta / fabs(theIpo[kIpomMx]-theIpo[kIpomMx-1]);
	yh = sigmaDTheta / dTheIpo;
        if( bHH  &&  vRC3270[kHH] ) vRC3270[kHH]->Fill( xh, yh );
//hh                                ----------------------------
        kHH = 7;
	bHH = kHH < int( vRC3270.size() );
        xh = momPart;
        yh = theWgt - theIpo[kIpomMx];
        if( bHH  &&  vRC3270[kHH] ) vRC3270[kHH]->Fill( xh, yh );
//hh                                ----------------------------

	kHH = 3;
	bHH = kHH < int( vRC3270.size() );
        xh = momPart;
        yh = meaCal - meaFit;
        if( bHH  &&  vRC3270[kHH] ) vRC3270[kHH]->Fill( xh, yh );
//hh                                ----------------------------
        kHH = 4;
	bHH = kHH < int( vRC3270.size() );
        xh = momPart;
        yh = sigCal - sigFit;
        if( bHH  &&  vRC3270[kHH] ) vRC3270[kHH]->Fill( xh, yh );
//hh                                ----------------------------

	nHGood++;

	if( !histLikeProfile( vRC8020, 8020, 20 ) ) {}
//           ------------------------------------

        if( !histQsqProfile( vRC8050, 8050, 20 ) ) {}
//           -----------------------------------

      ///}
    }

    return;
  }


//===========================================================================
  void CsRCPartPhotons::getLikeProb() {
//-------------------------------------


//- Paolo  -  April 2006


    CsRCExeKeys* key = CsRCExeKeys::Instance();
    CsRCRecConst *cons = CsRCRecConst::Instance();
    int kOffHi = cons->kOffHi();

    static std::vector<CsHist2D*> vRC3290;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      key->acknoMethod( "CsRCPartPhotons::getLikeProb" );

      if( CsRCHistos::Ref().bookHis() &&
	  CsRCHistos::Ref().levelBk() >= 2 ) {
	CsHistograms::SetCurrentPath("/RICH");
        string hTitle = "likelihood prob";
	int kHist = 0;
	stringstream hN3290;
	stringstream hN3291;
	stringstream hN3292;
	stringstream hN3293;
	stringstream hN3294;
	stringstream hN3295;
	stringstream hN3296;
	stringstream hN3297;
        kHist = kOffHi + 3290;
        hN3290 << kHist;
        vRC3290.push_back( new CsHist2D( hN3290.str(), hTitle,
					 100, 0., 50., 100, 0., 1.) );
        kHist = kOffHi + 3291;
        hN3291 << kHist;
        vRC3290.push_back( new CsHist2D( hN3291.str(), hTitle,
					 100, 0., 50., 100, -1., 1.) );
        kHist = kOffHi + 3292;
        hN3292 << kHist;
        vRC3290.push_back( new CsHist2D( hN3292.str(), hTitle,
					 100, 0., 50., 100, -1., 1.) );
        kHist = kOffHi + 3293;
        hN3293 << kHist;
        vRC3290.push_back( new CsHist2D( hN3293.str(), hTitle,
					 100, 0., 50., 100, -2.5, 2.5) );
        kHist = kOffHi + 3294;
        hN3294 << kHist;
        vRC3290.push_back( new CsHist2D( hN3294.str(), hTitle,
					 100, 0., 50., 100, -2.5, 2.5) );
        kHist = kOffHi + 3295;
        hN3295 << kHist;
        vRC3290.push_back( new CsHist2D( hN3295.str(), hTitle,
					 100, -5., 5., 100, 0., 1.) );
        kHist = kOffHi + 3296;
        hN3296 << kHist;
        vRC3290.push_back( new CsHist2D( hN3296.str(), hTitle,
					 100, -5., 5., 100, 0., 1.) );
        kHist = kOffHi + 3297;
        hN3297 << kHist;
        vRC3290.push_back( new CsHist2D( hN3297.str(), hTitle,
					 100, 0., 50., 100, -1., 1.) );
        CsHistograms::SetCurrentPath("/");
      }
    }

    double xh, yh, wh;

    float momPart = pPart_->mom();

//- Particle momentum selection :
//  -----------------------------
    if( momPart < 2.5 ) return;
//@@--------------------------

    double phiMass  = CsRichOne::Instance()->phiMass();

//- PHI MASS SELECTION! (for PHI Sample ONLY!!!)
//---------------------
    if( key->myFileType() == 3 ) {
      if( fabs( phiMass - 1.020) > 0.012 ) return;
//@@---------------------------------------------
    }

    double thePart = acos( pPart_->nd() );

//- Particle angle selection :
//  --------------------------
    //if( thePart < 50./1000.) return;
//@@-------------------------------

    double likeDV = cons->likeDefVa();
//@@---------------------------------


//- Background Likelihood :
//  -----------------------
    int nPhot = 0;
    double pLikeBk = likeDV;
    if( partProbsSet_ ) pLikeBk = probaLKBgAll_;
    else  pLikeBk = getLikelihood( 0., nPhot );
//                  --------------------------


//- Fine step Likelihood :
//  ----------------------
    double theMinc = 0.;
    double theMaxc = 0.;
    int nIpoc = 0;
    if( !likeProSet_ ) if( !getLikeProfile() ) return;
//  -------------------------------------------------
    theMinc = theLikeProMn_;
    theMaxc = theLikeProMx_;
    nIpoc   = nLikePro_;
    double dThec = (theMaxc - theMinc)/float( nIpoc-1 );
    double pLikekc[nIpoc];
    double thekc[nIpoc];
    for( int kIpo=0; kIpo<nIpoc; kIpo++ ) {
      thekc[kIpo] = theLikePro_[kIpo];
      pLikekc[kIpo] = pLikePro_[kIpo];
    }
    double pLikeMxc = -100000.;
    double theIpoMxc = -1.;
    int kIpoMxc = -1;
    for( int kIpo=0; kIpo<nIpoc; kIpo++ ) {
      if( pLikekc[kIpo] > pLikeMxc ) {
        pLikeMxc = pLikekc[kIpo];
        theIpoMxc = thekc[kIpo];
	kIpoMxc = kIpo;
      }
      //std::cout << setprecision( 6 );
      //std::cout << kIpo << "  " << thekc[kIpo] << "  " << pLikekc[kIpo] 
      //          << std::endl;
    }
    //std::cout << kIpoMxc << "  " << theIpoMxc << "  " << pLikeMxc
    //          << std::endl;
    if( theIpoMxc <  5. ) return;
    if( theIpoMxc > 60. ) return;
//@@----------------------------


//- Particle ID :
//  -------------
    int kTyLL = 2;
//@@-------------
    double theIpo[5];
    double pLikem[5];
    int kIpo = 0;
    for( int kPaTy=kTyLL; kPaTy<=14; kPaTy+=3 ) {
//  -------------------------------------------
      pLikem[kIpo] = likeDV;
      theIpo[kIpo] = thetaIpo( kPaTy );
      if( theIpo[kIpo] > 0. ) {
	if( partProbsSet_ ) pLikem[kIpo] = probaLKAll_[kPaTy];
	else  pLikem[kIpo] = getLikelihood( theIpo[kIpo], nPhot );
//                           ------------------------------------
      }
      kIpo++;
    }
    int nMassIpo = kIpo;
    double pLikemMx = -100000.;
    double theIpomMx = -1.;
    int kIpomMx = -1;
    for( int kIpo=0; kIpo<nMassIpo; kIpo++ ) {
      if( pLikem[kIpo] > pLikemMx ) {
	pLikemMx = pLikem[kIpo];
	theIpomMx = theIpo[kIpo];
	kIpomMx = kIpo;
      }
    }

    if( pLikemMx > pLikeBk ) {
      if( kIpomMx == 3 ) {
//@@----------------------

/*
//----- Optimize theta Max. :
//      ---------------------
	double theMax = 0.;
	int kMax = 0;
	if( kIpoMxc > 0  &&  kIpoMxc < nIpoc-1 ) {
	  int kmx = kIpoMxc;
  	  double delta = dThec;
	  double theMs = thekc[kmx] - delta;
	  double theCt = thekc[kmx];
	  double thePs = thekc[kmx] + delta;
	  double likeMs = pLikekc[kmx-1];
	  double likeCt = pLikekc[kmx];
	  double likePs = pLikekc[kmx+1];
	  double den = (likePs - 2.* likeCt + likeMs);
	  if( den != 0. ) {
	    theMax = theCt - delta/2.* (likePs - likeMs) / den;
	    kMax = kIpoMxc;
	    if( theMax < thekc[kmx] ) kMax--;
	  }
	}
	//std::cout << theMax << "  " << thekc[kIpoMxc] << std::endl;
	//std::cout << kMax << "  " << kIpoMxc << std::endl;
*/

	if( !thetaLikeSet_ ) if( !getThetaLikeMax() ) return;
//      ----------------------------------------------------
	double theMax = thetaLike_;
	int kMax = kIpoMxc;
	//if( theMax < thekc[kmx] ) kMax--;
//@@------------------------------------

//----- Compute like integral :
//      -----------------------
	double pLiOff = pLikeBk;
	double intTot = 0.;
	double intTots = 0.;
     	for( int kC=0; kC<nIpoc-1; kC++ ) {
          double bin = thekc[kC+1] - thekc[kC];
          double pLiS = (pLikekc[kC]+pLikekc[kC+1])/2.;
          intTot += pLiS * bin;
          pLiS -= pLiOff;
          if( pLiS < 0. ) pLiS = 0.;
          intTots += pLiS * bin;
        }
        //std::cout << intTot << "  " << intTots << " - ";

	double deThe[5];
	double lkProb[5];
	double lkProbs[5];
//----- For Pion and Kaon only :
//      ------------------------
	for( int kIw=2; kIw<4; kIw++ ) {
	  lkProb[kIw] = 0.;
	  lkProbs[kIw] = 0.;
	  if( theMax == 0.  ||  theMax == theMaxc ) continue;
	  if( theIpo[kIw] < 0. ) continue;

	  double theIw = theIpo[kIw];
	  double pLikeIw = pLikem[kIw];
	  //std::cout << theIw << "  " << theMax << " - ";
	  //if( theIw-theMax == -2. ) std::cout << "!!!!!";
	  double intProb = 0.;
	  double intProbs = 0.;
//------- theIpo > theMax :
//        -----------------
	  if( theIw > theMax ) {
	    //std::cout << " ( + ) ";
//--------- forward :
//          ---------
	    int kC1 = 0;
  	    for( int kC=kMax-1; kC<nIpoc-1; kC++ ) {
	      if( theIw > thekc[kC]  &&  theIw <= thekc[kC+1] ) {
		double bin = thekc[kC+1] - theIw;
	        double pLiS = (pLikeIw+pLikekc[kC+1])/2.;
		intProb = bin * pLiS;
	        pLiS -= pLiOff;
	        if( pLiS < 0. ) pLiS = 0.;
	        intProbs = bin * pLiS;
	        kC1 = kC + 1;
	        break;
	      }
	    }
	    //std::cout << intProb << "  " << intProbs << "  " << kC1 << "* ";
  	    for( int kC=kC1; kC<nIpoc-1; kC++ ) {
	      double bin = thekc[kC+1] - thekc[kC];
	      double pLiS = (pLikekc[kC]+pLikekc[kC+1])/2.;
	      intProb += pLiS * bin;
	      pLiS -= pLiOff;
	      if( pLiS < 0. ) pLiS = 0.;
	      intProbs += pLiS * bin;
	    }
	    //std::cout << intProb << "  " << intProbs << "  ";
//--------- backward :
//          ----------
	    double theIm = theMax - (theIw-theMax);
	    //std::cout << theIm << "*  ";
	    if( theIm < 0. ) theIm = 0.;
	    kC1 = 0;
	    for( int kC=kMax+1; kC>0; kC-- ) {
	      if( theIm < thekc[kC]  &&  theIm >= thekc[kC-1] ) {
		double bin = theIm - thekc[kC-1];
	        double pLiS = (pLikekc[kC]+pLikekc[kC-1])/2.;
		intProb += bin * pLiS;
	        pLiS -= pLiOff;
	        if( pLiS < 0. ) pLiS = 0.;
	        intProbs += bin * pLiS;
	        kC1 = kC - 1;
	        break;
	      }
	    }
	    //std::cout << intProb << "  " << intProbs << "  " << kC1 << "* ";
	    for( int kC=kC1; kC>0; kC-- ) {
	      double bin = thekc[kC] - thekc[kC-1];
	      double pLiS = (pLikekc[kC]+pLikekc[kC-1])/2.;
	      intProb += pLiS * bin;
	      pLiS -= pLiOff;
	      if( pLiS < 0. ) pLiS = 0.;
	      intProbs += pLiS * bin;
	    }
	    //std::cout << intProb << "  " << intProbs << " - ";

//------- theIpo < theMax :
//        -----------------
	  } else {
	    //std::cout << " ( - ) ";
//--------- backward :
//          ----------
	    int kC1 = 0;
  	    for( int kC=kMax+1; kC>0; kC-- ) {
	      if( theIw < thekc[kC]  &&  theIw >= thekc[kC-1] ) {
		double bin = theIw - thekc[kC-1];
	        double pLiS = (pLikeIw+pLikekc[kC-1])/2.;
		intProb = bin * pLiS;
	        pLiS -= pLiOff;
	        if( pLiS < 0. ) pLiS = 0.;
	        intProbs = bin * pLiS;
	        kC1 = kC - 1;
	        break;
	      }
	    }
	    //std::cout << intProb << "  " << intProbs << "  ";
	    for( int kC=kC1; kC>0; kC-- ) {
	      double bin = thekc[kC] - thekc[kC-1];
	      double pLiS = (pLikekc[kC]+pLikekc[kC-1])/2.;
	      intProb += pLiS * bin;
	      pLiS -= pLiOff;
	      if( pLiS < 0. ) pLiS = 0.;
	      intProbs += pLiS * bin;
	    }
	    //std::cout << intProb << "  " << intProbs << "  ";
//--------- forward :
//          ---------
	    double theIm = theMax + (theMax-theIw);
	    if( theIm > theMaxc ) theIm = theMaxc;
	    kC1 = nIpoc;;
	    for( int kC=kMax-1; kC<nIpoc-1; kC++ ) {
	      if( theIm > thekc[kC]  &&  theIm <= thekc[kC+1] ) {
		double bin = thekc[kC+1] - theIm;
	        double pLiS = (pLikekc[kC]+pLikekc[kC+1])/2.;
		intProb += bin * pLiS;
	        pLiS -= pLiOff;
	        if( pLiS < 0. ) pLiS = 0.;
	        intProbs += bin * pLiS;
	        kC1 = kC + 1;
	        break;
	      }
	    }
	    //std::cout << intProb << "  " << intProbs << "  ";
	    for( int kC=kC1; kC<nIpoc-1; kC++ ) {
	      double bin = thekc[kC+1] - thekc[kC];
	      double pLiS = (pLikekc[kC]+pLikekc[kC+1])/2.;
	      intProb += pLiS * bin;
	      pLiS -= pLiOff;
	      if( pLiS < 0. ) pLiS = 0.;
	      intProbs += pLiS * bin;
	    }
	    //std::cout << intProb << "  " << intProbs << " - ";
	  }
	  deThe[kIw] = theIw - theMax;
	  lkProb[kIw] = 0.;
	  if( intProb > intTot ) intProb = intTot;
	  if( intTot > 0. ) lkProb[kIw] = intProb / intTot;
	  if( intProbs > intTots ) intProbs = intTots;
	  if( intTots > 0. ) lkProbs[kIw] = intProbs / intTots;
	  //std::cout << lkProb[kIw] << "  " << lkProbs[kIw] << std::endl;

	}

        double meaCal = 0.;
	double sigCal = 0.;
	double sigmaP = sigCal;
        double thetaP = meaCal;
        double eK = (thetaP - theIpo[kIpomMx])/sigmaP;
        double probK = 1. - erf_( fabs( eK/sqrt(2.) ) );
        double ep = (thetaP - theIpo[kIpomMx-1])/sigmaP;
        double probp = 1. - erf_( fabs( ep/sqrt(2.) ) );
        //std::cout << eK << " probK = " << probK << "  ";
        //std::cout << ep << " probp = " << probp << std::endl;

//----- HISTOGRAMMING :
//      ---------------
        //std::cout << setprecision( 6 );
        //std::cout << pLikemMx << "  " << kIpomMx << "  " << std::endl;

	bool bHH = false;
	int kHH = 100;
	if( probK > 0.) {
  	  kHH = 0;
	  bHH = kHH < int( vRC3290.size() );
          xh = momPart;
	  yh = probK;
          if( bHH  &&  vRC3290[kHH] ) vRC3290[kHH]->Fill( xh, yh );
//hh                                  ----------------------------
 	  kHH = 1;
	  bHH = kHH < int( vRC3290.size() );
          xh = momPart;
	  yh = probK - probp;
          if( bHH  &&  vRC3290[kHH] ) vRC3290[kHH]->Fill( xh, yh );
//hh                                  ----------------------------
  	  kHH = 2;
	  bHH = kHH < int( vRC3290.size() );
          xh = momPart;
	  yh = (probK - probp)/probK;
          if( bHH  &&  vRC3290[kHH] ) vRC3290[kHH]->Fill( xh, yh );
//hh                                  ----------------------------
	}
	if( thetaLikeSet_ ) {
	  int kHH = 3;
	  bHH = kHH < int( vRC3290.size() );
	  xh = momPart;
	  yh = thetaLike_ - meaCal;
          if( bHH  &&  vRC3290[kHH] ) vRC3290[kHH]->Fill( xh, yh );
//hh                                  ----------------------------
	  kHH = 4;
	  bHH = kHH < int( vRC3290.size() );
	  xh = momPart;
	  yh = thetaLike_ - theIpoMxc;
          if( bHH  &&  vRC3290[kHH] ) vRC3290[kHH]->Fill( xh, yh );
//hh                                  ----------------------------
	}
	kHH = 5;
	bHH = kHH < int( vRC3290.size() );
        xh = deThe[3];
        yh = lkProb[3];
        if( bHH  &&  vRC3290[kHH] ) vRC3290[kHH]->Fill( xh, yh );
//hh                                ----------------------------
        kHH = 6;
	bHH = kHH < int( vRC3290.size() );
        xh = deThe[3];
        yh = lkProbs[3];
        if( bHH  &&  vRC3290[kHH] ) vRC3290[kHH]->Fill( xh, yh );
//hh                                ----------------------------
	kHH = 7;
	bHH = kHH < int( vRC3290.size() );
        xh = momPart;
        yh = lkProbs[3] - lkProbs[2];
        if( bHH  &&  vRC3290[kHH] ) vRC3290[kHH]->Fill( xh, yh );
//hh                                ----------------------------
      }
    }

    return;
  }


//===========================================================================
  int CsRCPartPhotons::getRelMaxima( int npo, double* array,
//----------------------------------------------------------
				     int& nRelMax, int* iRelMax ) {

//- Paolo  -  April 2006


    CsRCExeKeys* key = CsRCExeKeys::Instance();
    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      key->acknoMethod( "CsRCPartPhotons::getRelMaxima" );
    }

    //std::cout << "getRelMaxima  " << setprecision(5);
    //for( int k=0; k<npo; k++ ) std::cout << k << " " << array[k] << "  ";
    //std::cout << std::endl;

    int retCode = -1;

    static const int nDim = 20;
    bool dump = false;
//  -----------------

    int nMinMx = 10;
//@@---------------
    double relMax[nDim];
    int kRelMax[nDim];
    int kRelMinDw[nDim];
    int kRelMinUp[nDim];
    int nPass = 10;

    int pLimDw[nDim];
    int pLimUp[nDim];
    int pLimDwO[nDim];
    int pLimUpO[nDim];

    int npu = npo - 1;
    pLimDw[0] = 0;
    pLimUp[0] = npu;

    int nPoUse = 0;
    int kMin = 0;
    int nZone = 1;
    if( nZone > nDim ) nZone = nDim;
    bool bBreak = false;
    for( int kPs=0; kPs<nPass; kPs++ ) {
      //std::cout << "kPs " << kPs << std::endl;

      int kZone = 0;
      for( int kZo=0; kZo<nZone; kZo++ ) {
	//std::cout << "   kZo " << kZo << std::endl;
        int np1 = pLimDw[kZo];
        int np2 = pLimUp[kZo];
	if( np2 <= np1 ) continue;

        double maxArr = 0.;
        int maxPo = -1;
        for( int kp=np1; kp<np2; kp++ ) {
          if( array[kp] > maxArr ) {
	    maxArr = array[kp];
	    maxPo = kp;
          }
        }
        int minPoDw = np1;
        for( int kp=maxPo-1; kp>=np1; kp-- ) {
          if( array[kp] < array[kp+1] ) continue;
          else {
	    minPoDw = kp;
	    break;
          }
        }
	if( minPoDw < 0 ) {
	  retCode = -1;
	  if( dump ) std::cout << "getRelMaxima : L-limit ovfl" << std::endl;
	  return  retCode;
	}
        int minPoUp = np2;
        for( int kp=maxPo+1; kp<np2; kp++ ) {
          if( array[kp] < array[kp-1] ) continue;
          else {
	    minPoUp = kp;
	    break;
          }
        }
	if( minPoUp > npu ) {
	  retCode = -1;
	  if( dump ) std::cout << "getRelMaxima : U-limit ovfl" << std::endl;
	  return  retCode;
	}
	if( minPoDw > minPoUp ) {
	  retCode = -1;
	  if( dump ) std::cout << "getRelMaxima : limits ovfl" << std::endl;
	  return  retCode;
	}

        relMax[kMin] = maxArr; 
        kRelMax[kMin] = maxPo;
        kRelMinDw[kMin] = minPoDw;
        kRelMinUp[kMin] = minPoUp;
	nPoUse += minPoUp - minPoDw + 1;
	//std::cout << "nPoUse = " << nPoUse << std::endl;
	kMin++;
	if( kMin >= nDim ) {
	  retCode = 1;
	  bBreak = true;
	  if( dump ) std::cout << "getRelMaxima : --- n-Min ovfl - D"
			       << std::endl;
	  break;
	}
	if( kMin >= nMinMx ) {
	  retCode = 1;
	  bBreak = true;
	  if( dump ) std::cout << "getRelMaxima : --- n-Min ovfl - R"
			       << std::endl;
	  break;
	}
	//std::cout << "      kMin " << kMin << std::endl;
	//std::cout << maxArr << "  " << maxPo << "  " << minPoDw
	//  	    << "  " << minPoUp << std::endl;
	if( nPoUse == npo ) {
	  retCode = 0;
	  bBreak = true;
	  if( dump ) std::cout << "getRelMaxima : ------ reg end - Use"
			       << std::endl;
	  break;
	}
	if( minPoDw == 0  &&  minPoUp == npu ) {
	  retCode = 0;
	  bBreak = true;
	  if( dump ) std::cout << "getRelMaxima : ------ reg end - 0n"
			       << std::endl;
	  break;
	}

	if( minPoDw-1 >  np1 ) {
  	  pLimDwO[kZone] = np1;
	  pLimUpO[kZone] = minPoDw-1;
	  kZone++;
	}
	if( np2 > minPoUp+1 ) {
  	  pLimDwO[kZone] = minPoUp+1;
	  pLimUpO[kZone] = np2;
	  kZone++;
	}
        //std::cout << np1 << "  " << minPoDw-1 << "  " << minPoUp+1
	//	    << "  " << np2 << std::endl;
      }
      if( bBreak ) break;
      if( kZone > nDim ) {
	  retCode = 1;
	  if( dump ) std::cout << "getRelMaxima : --- n-Zon ovfl"
			       << std::endl;
	  break;
      }
      nZone = kZone;
      for( int kp=0; kp<nZone; kp++ ) {
	pLimDw[kp] = pLimDwO[kp];
	pLimUp[kp] = pLimUpO[kp];
      }
      //std::cout << "Pass = " << kPs << std::endl;
      //for( int kp=0; kp<nZone; kp++ ) std::cout << pLimDw[kp] << "-"
      //					  << pLimUp[kp] << "  ";
      //std::cout << std::endl;
    }

    int kLikeEq = 0;
    for(int kp=0; kp<kMin; kp++ ) {
      int kEq = 0;
      for(int kq=kp+1; kq<kMin; kq++ ) {
	if( array[kRelMax[kp]] == array[kRelMax[kq]] ) kEq++;
      }
      if( kEq > 0 ) kLikeEq++;
    }
    if( kLikeEq >= kMin-2 ) {
      retCode = 2;
	  if( dump ) std::cout << "getRelMaxima : ! const like" << std::endl;
      return  retCode;
    }
    nRelMax = kMin;
    for(int kp=0; kp<nRelMax; kp++ ) iRelMax[kp] = kRelMax[kp];

    return  retCode;

  }


//===========================================================================
  bool CsRCPartPhotons::getLikeProfile() {
//----------------------------------------

//- Paolo  -  April 2006


    CsRCExeKeys* key = CsRCExeKeys::Instance();

    if( likeProSet_ ) return  true;
//  ------------------------------            

    //theLikeProMn_ =   0.;
    theLikeProMn_ =  10.;
    theLikeProMx_ =  60.;
    //nLikePro_ =  241;
    //nLikePro_ =   51;
    nLikePro_ =   101;
//@@--------------------

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;
      key->acknoMethod( "CsRCPartPhotons::getLikeProfile" );

      std::cout << " RICHONE, CsRCPartPhotons::getLikeProfile on   "
		<< nLikePro_ << "   points" << std::endl;
//  -----------------------------------------------------------------
    }

    double dThe = (theLikeProMx_ - theLikeProMn_)/float( nLikePro_-1 );
    int nPhot = 0;
    for( int kIpo=0; kIpo<nLikePro_; kIpo++ ) {
      theLikePro_[kIpo] = theLikeProMn_ + kIpo * dThe;
      pLikePro_[kIpo] = getLikelihood( theLikePro_[kIpo], nPhot );
//                      -----------------------------------------
    }
    likeProSet_ = true;
    //for( int kIpo=0; kIpo<nLikePro_; kIpo++ ) std::cout <<
    //					pLikePro_[kIpo] << "  ";
    //std::cout << std::endl;
    //std::cout << "getLikeProfile: like profile calculated, ev  "
    //<< CsRichOne::Instance()->kEvent() << "  " << this << std::endl;

    return  true;
  }


//===========================================================================
  bool CsRCPartPhotons::getThetaLikeMax() {
//-----------------------------------------


//- theta ring from Maximum Likelihood
//  ----------------------------------
//- Paolo  -  September 2002
//    rev.    December  2002 
//    rev.    February  2003
//    rev.    October   2004
//    rev.    January   2005
//    ----------------------
//    rev.    April     2006


    CsRCExeKeys* key = CsRCExeKeys::Instance();
    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      key->acknoMethod( "CsRCPartPhotons::getThetaLikeMax2" );
    }

    CsRCRecConst* cons = CsRCRecConst::Instance();

    double delta, den, theMs, theCt, thePs, likeMs, likeCt, likePs;


    double likeDV = cons->likeDefVa();
//@@---------------------------------

    bool exeRM = true;
    //bool exeRM = false;
//@@------------------
    static int nGood = 0;
    bool dump = false;
//  -----------------

    double theKinLim = acos( 1./cons->CFRefInd() ) * 1000.;
    //std::cout << setprecision(6);
    //std::cout << theKinLim << "  " << cons->CFRefInd() << std::endl;
//- include approx. experimental resolution (3sig, mrad)
    theKinLim += 1.2;

//- Fine scan :
//  -----------
    double theLiMin = 0.;
    double theLiMax = 0.;
    int nIpo = 0;
    if( !likeProSet_ ) if( !getLikeProfile() ) return  false;
//                          ----------------   -------------

    theLiMin = theLikeProMn_;
    theLiMax = theLikeProMx_;
    nIpo = nLikePro_;
    double dThe = (theLiMax - theLiMin) / float( nIpo-1 );
    double theV[nIpo];
    double likeV[nIpo];
    int nPhot = 0;
    for( int kIpo=0; kIpo<nIpo; kIpo++ ) {
      theV[kIpo] = theLikePro_[kIpo];
      likeV[kIpo] = pLikePro_[kIpo];
    }

    int nRelMax = 0;
    int iRelMax[10];

    double likeM = likeDV;
    double likeMax = -100000.;
    double theMax = -1.;
    int kMax = 0;
    for( int kIpo=0; kIpo<nIpo; kIpo++ ) {
      likeM = likeV[kIpo];
      if( !cons->likeLog() ) if( likeM < 0. ) goto pRet;
//                                            ---------
      if( likeM > likeMax ) {
	likeMax = likeM;
	theMax = theV[kIpo];
	kMax = kIpo;
      }
    }
//- Avoid constant like
    if( kMax <= 1 ) goto pRet;
//                  ---------

    if( dump ) {
      std::cout << "getThetaLikeMax  ";
      std::cout << likeMax << "  " << theMax << "  " << kMax << std::endl;
    }
    if( theMax <= theLiMin ) goto pRet;
//                           ---------

    if( exeRM ) {
      int retCode = getRelMaxima( nIpo, &likeV[0], nRelMax, &iRelMax[0] );
//                  -----------------------------------------------------
      //if( retCode == -1 ) goto pRet;
//                        ---------
//--- Const like
      if( retCode == 2 ) goto pRet;
//                       ---------
      if( retCode == 0  ||  retCode == 1 ) {
	double theMaxR = 0.;
	int kmR = -1;
        for(int km=0; km<nRelMax; km++ ) {
          theMaxR = theLikeProMn_ +  dThe* iRelMax[km];
	  //if( theMaxR >= theLiMax ) continue;
	  if( theMaxR > theKinLim ) continue;
  	  kmR = km;
	  break;
        }
	if( kmR >= 0 ) {
	  kMax = iRelMax[kmR];
	  theMax = theMaxR;
	}
	if( dump ) {
	  std::cout << "                 " << likeV[kMax] << "  ";
	  std::cout << theMax << "  " << kMax << "  " << kmR << std::endl;
	  for(int k=2; k<5; k++) std::cout << thetaIpo_[3*k+2] << "  ";
	  std::cout << std::endl;
	}
      }
    }
    //if( theMax >= theLiMax ) goto pRet;
    if( theMax >= theKinLim ) goto pRet;
//                            ---------

//- Interpolation :
//  ---------------
    if( kMax <= 0 ) kMax = 1;
    if( kMax >= nIpo-1 ) kMax = nIpo - 2;
    delta = dThe;
    theMs = theMax - delta;
    theCt = theMax;
    thePs = theMax + delta;
    likeMs = likeV[kMax-1];
    likeCt = likeMax;
    likePs = likeV[kMax+1];
    den = (likePs - 2.* likeCt + likeMs);
    if( den == 0. ) goto pRet;
//                  ---------
    theMax = theCt - delta/2. * (likePs - likeMs) / den;
//  ---------------------------------------------------
    if( theMax <= 0. ) goto pRet;
//                     ---------
    if( theMax < theMs ) theMax = theMs;
    if( theMax > thePs ) theMax = thePs;
    if( theMax >= theLiMax ) goto pRet;
//                           ---------
//- Interpolation 2 :
//  -----------------
    delta /= 2.;
    theMs = theMax - delta;
    theCt = theMax;
    thePs = theMax + delta;
    likeMs = getLikelihood( theMs, nPhot );
//           -----------------------------
    if( !cons->likeLog() ) if( likeMs < 0. ) goto pRet;
//                                           ---------
    likeCt = getLikelihood( theCt, nPhot );
//           -----------------------------
    if( !cons->likeLog() ) if( likeCt < 0. ) goto pRet;
//                                           ---------
    likePs = getLikelihood( thePs, nPhot );
//           -----------------------------
    if( !cons->likeLog() ) if( likePs < 0. ) goto pRet;
//                                           ---------
    //cout << likeMs << "  " << likeCt << "  " << likePs << endl;
    den = (likePs - 2.* likeCt + likeMs);
    if( den == 0. ) goto pRet;
//                  ---------
    theMax = theCt - delta/2. * (likePs - likeMs) / den;
//  ---------------------------------------------------
    if( theMax <= 0. ) goto pRet;
//                     ---------
    if( theMax < theMs ) theMax = theMs;
    if( theMax > thePs ) theMax = thePs;
    if( theMax >= theLiMax ) goto pRet;
//                           ---------
    likeMax = getLikelihood( theMax, nPhot );
//            ------------------------------
    if( !cons->likeLog() ) if( likeMax < 0. ) goto pRet;
//                                            ---------
    if( dump ) {
      std::cout << "                 ";
      std::cout << likeMax << "  " << theMax << std::endl;
    }

    thetaLike_ = theMax;
    pLikeMax_ = likeMax;
    thetaLikeSet_ = true;
//  --------------------

    nGood++;
    //std::cout << "nGood  " << nGood << std::endl;
    return  true;


    pRet :
      return  false;

  }


//===========================================================================
  bool CsRCPartPhotons::histLikeProfile( std::vector<CsHist1D*> &vRCh,
//--------------------------------------------------------------------
					 const int nRCh, const int nHHMx ) {

//- Paolo  -  May 2006


    CsRCRecConst *cons = CsRCRecConst::Instance();
    int kOffHi = cons->kOffHi();

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      CsRCExeKeys* key = CsRCExeKeys::Instance();
      key->acknoMethod( "CsRCPartPhotons::histLikeProfile" );
    }

    if( int( vRCh.size() ) >= nHHMx ) return  false;
//                                    -------------

    double likeDV = cons->likeDefVa();
//@@---------------------------------

//- get Profile :
//  -------------
    double theMin = 0.;
    double theMax = 0.;
    int nIpo = 0;
    if( !likeProSet_ ) if( !getLikeProfile() ) return  false;
//                          ----------------   -------------
    theMin = theLikeProMn_;
    theMax = theLikeProMx_;
    nIpo = nLikePro_;
    double dThe = (theMax - theMin)/float( nIpo-1 );
    double thek[nIpo];
    double pLikek[nIpo];
    for( int kIpo=0; kIpo<nIpo; kIpo++ ) {
      thek[kIpo] = theLikePro_[kIpo];
      pLikek[kIpo] = pLikePro_[kIpo];
    }
    int nPhot = 0;
    double pLikeBk = likeDV;
    if( partProbsSet_ ) pLikeBk = probaLKBgAll_;
    else  pLikeBk = getLikelihood( 0., nPhot );
//                  --------------------------

    int kHist = kOffHi + nRCh + vRCh.size();
    stringstream hNh;
    hNh << kHist;
    string hTitle = "likelihood shape";
    CsHistograms::SetCurrentPath("/RICH");
    CsHist1D* hist = new CsHist1D( hNh.str(), hTitle,
//hh------------------------------------------------
				   nIpo+10, theMin, theMax+10*dThe );
    vRCh.push_back( hist );
    CsHistograms::SetCurrentPath("/");
    double xh, yh, wh;

    for( int kIpo=0; kIpo<nIpo; kIpo++ ) {
      xh = thek[kIpo];
      wh = pLikek[kIpo];
      if( hist ) hist->Fill( xh, wh );
//hh             --------------------
    }
    xh = theMax + dThe;
    wh = pPart_->mom();
    if( hist ) hist->Fill( xh, wh );
//hh           --------------------
    int kTyLL = 2;
//@@-------------
    int kIpo = 0;
    for( int kPaTy=kTyLL; kPaTy<=14; kPaTy+=3 ) {
      xh = theMax + (kIpo+2)*dThe;
      wh = thetaIpo( kPaTy );
      if( wh < 0. ) wh = 0.;
      if( hist ) hist->Fill( xh, wh );
//hh             --------------------
      kIpo++;
    }
    int nMassIpo = kIpo;
    xh = theMax + 7*dThe;
    wh = pLikeBk;
    if( hist ) hist->Fill( xh, wh );
//hh           --------------------
    xh = theMax + 8*dThe;
    wh = lPhotons_.size();
    if( hist ) hist->Fill( xh, wh );
//hh           --------------------

    std::cout << "histLikeProfile : Profile nr. " << int( vRCh.size() )
	      << "  produced (" << hist->GetName() << ")" << std::endl;
//  -------------------------------------------------------------------

    return  true;

  }


//===========================================================================
  bool CsRCPartPhotons::histQsqProfile( std::vector<CsHist1D*> &vRCh,
//-------------------------------------------------------------------
					const int nRCh, const int nHHMx ) {

//- Paolo  -  May 2006


    CsRCRecConst *cons = CsRCRecConst::Instance();
    int kOffHi = cons->kOffHi();

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      CsRCExeKeys* key = CsRCExeKeys::Instance();
      key->acknoMethod( "CsRCPartPhotons::histQsqProfile" );
    }

    if( int( vRCh.size() ) >= nHHMx ) return  false;
//                                    -------------

//- get Profile :
//  -------------
    double theMin = 10.;
    double theMax = 60.;
    int nIpo = 51;
    if( likeProSet_ ) {
      theMin = theLikeProMn_;
      theMax = theLikeProMx_;
      nIpo = nLikePro_;
    }
    double dThe = (theMax - theMin)/float( nIpo-1 );
    double thek[nIpo];
    double qSquak[nIpo];
    int nPhot = 0;
    for( int kIpo=0; kIpo<nIpo; kIpo++ ) {
      if( likeProSet_ ) thek[kIpo] = theLikePro_[kIpo];
      else  thek[kIpo] = theMin + kIpo * dThe;
      qSquak[kIpo] = getQsquare( thek[kIpo], nPhot );
//                   -------------------------------
    }

    int kHist = kOffHi + nRCh + vRCh.size();
    stringstream hNh;
    hNh << kHist;
    string hTitle = "likelihood shape";
    CsHistograms::SetCurrentPath("/RICH");
    CsHist1D* hist = new CsHist1D( hNh.str(), hTitle,
//hh------------------------------------------------
				   nIpo+10, theMin, theMax+10*dThe );
    vRCh.push_back( hist );
    CsHistograms::SetCurrentPath("/");
    double xh, yh, wh;

    for( int kIpo=0; kIpo<nIpo; kIpo++ ) {
      xh = thek[kIpo];
      wh = qSquak[kIpo];
      if( hist ) hist->Fill( xh, wh );
//hh             --------------------
    }
    xh = theMax + dThe;
    wh = pPart_->mom();
    if( hist ) hist->Fill( xh, wh );
//hh           --------------------
    int kTyLL = 2;
//@@-------------
    int kIpo = 0;
    for( int kPaTy=kTyLL; kPaTy<=14; kPaTy+=3 ) {
      xh = theMax + (kIpo+2)*dThe;
      wh = thetaIpo( kPaTy );
      if( wh < 0. ) wh = 0.;
      if( hist ) hist->Fill( xh, wh );
//hh             --------------------
      kIpo++;
    }

    return  true;

  }


//===========================================================================
  double CsRCPartPhotons::nPhotExpct( double theIpo ) {
//-----------------------------------------------------


//- Paolo  -  May 2006


    CsRCRecConst *cons = CsRCRecConst::Instance();
    static const int kOffHi = cons->kOffHi();
    static const double TwoPi = cons->TwoPI();
    CsRCDetectors* dets = CsRCDetectors::Instance();
    int nCathode = dets->nCathode();

    static const float focLen = pPart_->pMirPart()->RR() / 2.;

    static int iCaUp =  3;      //   pos Y
    static int iCaDw = 10;      //   neg Y
    static int iCaLf =  5;      //   neg X
    static int iCaRg =  3;      //   pos X
    static double xLimLf = 0.;
    static double xLimRg = 0.;
    static double yLimUp = 0.;
    static double yLimDw = 0.;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      CsRCExeKeys* key = CsRCExeKeys::Instance();
      key->acknoMethod( "CsRCPartPhotons::nPhotExpct" );

      CsRCCathode* cat = NULL;
      double hCatx, hCaty;
      Hep3Vector vOffCatWw;
      cat = dets->ptrToCat( iCaUp );
      hCaty = cat->hCaty();
      vOffCatWw = dets->vOffCatW( iCaUp );
      yLimUp = vOffCatWw.y() + hCaty;
      cat = dets->ptrToCat( iCaDw );
      hCaty = cat->hCaty();
      vOffCatWw = dets->vOffCatW( iCaDw );
      yLimDw = vOffCatWw.y() - hCaty;
      cat = dets->ptrToCat( iCaLf );
      hCatx = cat->hCatx();
      vOffCatWw = dets->vOffCatW( iCaLf );
      xLimLf = vOffCatWw.x() - hCatx;
      cat = dets->ptrToCat( iCaRg );
      hCatx = cat->hCatx();
      vOffCatWw = dets->vOffCatW( iCaRg );
      xLimRg = vOffCatWw.x() + hCatx;
    }

    double nPhotEx = 0;

    double theIrad = theIpo / 1000.;
    float pathLen = pPart_->pathLen();
    nPhotEx = cons->nZero() * pathLen/10. * theIrad*theIrad;
//---------------------------------------------------------

    if( !CsRichOne::Instance()->UpRICHJob() ) return  nPhotEx;
//  ---------------------------------------------------------

    double nPhotExUV = nPhotEx;
    double nPhotExVS = cons->nZeroVS() * pathLen/10. * theIrad*theIrad;
//--------------------------------------------------------------------
    double rTheIpo = theIrad * focLen;
    int kDetPart = pPart_->kDetPart();
    Hep3Vector vPade = pPart_->vPade()[kDetPart];

    int iCaPa = iCaPa_;
    bool catPMT = false;
    double HPi = TwoPi/4;
    if( iCaPa >= 0 ) catPMT = dets->ptrToCat( iCaPa )->isPMT();
    else {
      if( vPade.y() < yLimUp  &&  vPade.y() > yLimDw  &&
          vPade.x() < xLimRg  &&  vPade.x() > xLimLf ) catPMT = true;
      else  catPMT = false;
    }
    double frcUV = 0.;
    double frcVS = 0.;
    if( catPMT ) {
//--- Internal cathodes (PMT)
      double ddx = 0.;
      double ddy = 0.;
      if( vPade.y() >= 0.) ddy = yLimUp - vPade.y();
      if( vPade.y() <  0.) ddy = vPade.y() - yLimDw;
      if( vPade.x() >= 0.) ddx = xLimRg - vPade.x();
      if( vPade.x() <  0.) ddx = vPade.x() - xLimLf;
      frcUV = 0.;
      if( rTheIpo > 0. ) {
        if( ddx >= rTheIpo  &&  ddy >= rTheIpo ) frcUV = 0.;
        else if( ddx >= rTheIpo  &&  ddy < rTheIpo ) {
	  frcUV = 2.* acos( ddy/rTheIpo );
        }
        else if( ddx < rTheIpo  &&  ddy >= rTheIpo ) {
	  frcUV = 2.* acos( ddx/rTheIpo );
        }
        else if( ddx < rTheIpo  &&  ddy < rTheIpo ) {
	  frcUV = acos( ddy/rTheIpo ) + acos( ddx/rTheIpo ) + HPi;
        }
	frcUV /= TwoPi;
	if( frcUV > 1. ) frcUV = 1.;
      }
      nPhotEx = nPhotExVS * (1.- frcUV) + nPhotExUV * frcUV;
      if( nPhotEx < 0 ) nPhotEx = 0;
      //if( isnan( nPhotEx ) )
      //std::cout << " UV  " << frcUV << "  " << rTheIpo << "  ";
    } else {
//--- External cathodes (CsI)
      double ddx = 0.;
      double ddy = 0.;
      if( vPade.y() >= 0.) ddy = vPade.y() - yLimUp;
      if( vPade.y() <  0.) ddy = yLimDw - vPade.y();
      if( vPade.x() >= 0.) ddx = vPade.x() - xLimRg;
      if( vPade.x() <  0.) ddx = xLimLf - vPade.x();
      frcVS = 0.;
      if( rTheIpo > 0. ) {
        if( ddx >= rTheIpo  ||  ddy >= rTheIpo ) frcVS = 0.;
        else if( ddx < 0.  &&  ddy > 0. ) {
	  if( fabs( ddx ) >= rTheIpo ) {
	    frcVS = 2.* rTheIpo * acos( ddy/rTheIpo );
	  } else {
	    frcVS = acos( ddy/rTheIpo ) - acos( fabs(ddx)/rTheIpo ) + HPi;
	  }
	}
	else if( ddx > 0.  &&  ddy < 0. ) {
	  if( fabs( ddy ) >= rTheIpo ) {
  	    frcVS = 2.* rTheIpo * acos( ddx/rTheIpo );
	  } else {
	    frcVS = acos( ddx/rTheIpo ) - acos( fabs(ddy)/rTheIpo ) + HPi;
	  }
        }
        else if( ddx > 0.  &&  ddy > 0. ) {
	  frcVS = acos( ddy/rTheIpo ) + acos( ddx/rTheIpo ) - HPi;
        }
	frcVS /= TwoPi;
	if( frcVS > 1. ) frcVS = 1.;
      }
      nPhotEx = nPhotExUV * (1.- frcVS) + nPhotExVS * frcVS;
      if( nPhotEx < 0 ) nPhotEx = 0;
    }


//- Histogramms :
//  -------------
    CsRCHistos& hist = CsRCHistos::Ref();
    int level = hist.levelBk();
    bool doHist = false;
    if ( level >= 2 ) doHist = true;
    if( doHist ) {
      if( likeONLY_ ) {
        double xh, yh, wh;
        xh = vPade.x();
        yh = vPade.y();
        wh = double( nPhotEx );
        if( hist.hRC1505 ) hist.hRC1505->Fill( xh, yh, wh );
//      ---------------------------------------------------
        xh = vPade.x();
        yh = vPade.y();
        if( hist.hRC1506 ) hist.hRC1506->Fill( xh, yh );
//      -----------------------------------------------
	if( catPMT ) {
	  xh = frcUV;
	  if( hist.hRC1507 ) hist.hRC1507->Fill( xh );
//        -------------------------------------------
	  xh = nPhotEx;
	  if( hist.hRC1572 ) hist.hRC1572->Fill( xh );
//        -------------------------------------------
	} else {
	  xh = frcVS;
	  if( hist.hRC1508 ) hist.hRC1508->Fill( xh );
//        -------------------------------------------
	  xh = nPhotEx;
	  if( hist.hRC1573 ) hist.hRC1573->Fill( xh );
//        -------------------------------------------
	}
      }
    }

    return  nPhotEx;

  }


extern "C" { double rndm_( double& ); }

//===========================================================================
  Hep3Vector CsRCPartPhotons::doPMTOptCorr( CsRCCluster* clu,
//-----------------------------------------------------------
					    const Hep3Vector vPoPhotW,
					    const Hep3Vector vDcPhoEmW,
					    Hep3Vector& vPoCluDet,
					    const float RR,
					    bool& aflag, bool& oflag ) {

//- Paolo  -  December 2006


    Hep3Vector vPoCluDetO = vPoCluDet;
    Hep3Vector vDcPhoCorO = vDcPhoEmW;

    static bool dump = false;
    //static bool dump = true;

//- AMPS ---
//- Default choice between AM or PS Optical Corrections :
//  -----------------------------------------------------
    static bool bAMCorr = true;      // Use AM Corrections ( Default ) 
    //static bool bAMCorr = false;
//@@---------------------------
    //static bool bPSCorr = true;      // Use PS Corrections
    static bool bPSCorr = false;
//@@-----------------------------
    if( !bAMCorr  &&  !bPSCorr ) bAMCorr = true;
    //if( !bAMCorr  &&  !bPSCorr ) return  vDcPhoEmW;
//                               -----------------
    if( bAMCorr  &&  bPSCorr ) bAMCorr = false;      // For production
//- can be commented out for comparison tests...

    static int nTot = 0;
    static int nCor = 0;
    static int nOutx = 0;
    static int nOuty = 0;
    static int n12 = 0;
    static int nQx[4];
    static int nQy[4];
    if( nTot == 0 ) {
      for( int k=0; k<4; k++ ) nQx[k] = 0;
      for( int k=0; k<4; k++ ) nQy[k] = 0;
    }
    static double xDiffM = 0.;
    static double yDiffM = 0.;

    nTot++;
    if( nTot%100 == 0 ) {
      //std::cout << "- doPMTOptCorr : " << nCor << "  " << nTot << std::endl;
    }

    if( !clu->isPMT() ) {
      if( dump ) std::cout << "doPMTOptCorr - NO PMT" << std::endl;
      oflag = false;
      return  vDcPhoEmW;
//    -----------------
    }

    CsRCDetectors *dets = CsRCDetectors::Instance();
    CsRCRecConst *cons = CsRCRecConst::Instance();

    CsRCHistos& hist = CsRCHistos::Ref();
    int kOffHi = cons->kOffHi();
    static std::vector<CsHist2D*> vRC3320;
    static std::vector<CsHist2D*> vRC6620;
    double xh, yh, wh;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      CsRCExeKeys* key = CsRCExeKeys::Instance();
      key->acknoMethod( "CsRCPartPhotons::doPMTOptCorr" );

      CsOpt* opt = CsOpt::Instance();
      static const std::string AMCorr = "AMCO";
      static const std::string PSCorr = "PSCO";
      bool boo;
      std::string sOpw;
      boo = opt->CsOpt::getOpt( "RICHONE", "PMTOptCorrType", sOpw );
      if( boo ) {
	if ( sOpw == AMCorr ) { bAMCorr = true; bPSCorr = false; }
	if ( sOpw == PSCorr ) { bPSCorr = true; bAMCorr = false; }
        if( bAMCorr ) std::cout << " RICHONE, "
	  << "CsRCPartPhotons::doPMTOptCorr : AM corrections used (set)"
		  	        << std::endl;
        if( bPSCorr ) std::cout << " RICHONE, "
          << "CsRCPartPhotons::doPMTOptCorr : PS corrections used (set)"
			        << std::endl;
      } else {
        if( bAMCorr ) std::cout << " RICHONE, "
	  << "CsRCPartPhotons::doPMTOptCorr : AM corrections used (default)"
			        << std::endl;
        if( bPSCorr ) std::cout << " RICHONE, "
          << "CsRCPartPhotons::doPMTOptCorr : PS corrections used (default)"
			        << std::endl;
      }

      for( int kh=0; kh<70; kh++ ) vRC3320.push_back( NULL );
      for( int kh=0; kh<10; kh++ ) vRC6620.push_back( NULL );
      if( CsRCHistos::Ref().bookHis() &&
	  CsRCHistos::Ref().levelBk() >= 2 ) {
	CsHistograms::SetCurrentPath("/RICH");
	vRC3320.clear();
	int nch = 20;
        int kHist = 0;
	for( int kPMTch=0; kPMTch<16; kPMTch++ ) {
	  string hTitle = "xCorrAM-PS";
	  stringstream hN3321;
	  kHist = kOffHi + 3321 + 4*kPMTch;
	  hN3321 << kHist;
	  vRC3320.push_back( new CsHist2D( hN3321.str(), hTitle,
					   nch, -0.1, 0.4, nch, -0.5, 0. ) );
	  hTitle = "yCorrAM-PS";
	  stringstream hN3322;
	  kHist = kOffHi + 3322 + 4*kPMTch;
	  hN3322 << kHist;
	  vRC3320.push_back( new CsHist2D( hN3322.str(), hTitle,
					   nch, -0.1, 0.4, nch, -0.5, 0. ) );
	  hTitle = "c-xCorrAM-PS";
	  stringstream hN3323;
	  kHist = kOffHi + 3323 + 4*kPMTch;
	  hN3323 << kHist;
	  vRC3320.push_back( new CsHist2D( hN3323.str(), hTitle,
					   nch, -0.1, 0.4, nch, -0.5, 0. ) );
	  hTitle = "c-yCorrAM-PS";
	  stringstream hN3324;
	  kHist = kOffHi + 3324 + 4*kPMTch;
	  hN3324 << kHist;
	  vRC3320.push_back( new CsHist2D( hN3324.str(), hTitle,
					   nch, -0.1, 0.4, nch, -0.5, 0. ) );
	}
	vRC6620.clear();
	for( int kQua=0; kQua<4; kQua++ ) {
	  string hTitle = "tgY vs tgX";
	  stringstream hN6620;
	  kHist = kOffHi + 6620 + kQua + 1;
	  hN6620 << kHist;
	  vRC6620.push_back( new CsHist2D( hN6620.str(), hTitle,
					   50, -0.1, 0.4, 50, -0.5, 0. ) );
	}
	for( int kQua=0; kQua<4; kQua++ ) {
	  string hTitle = "PMT occ";
	  stringstream hN6620;
	  kHist = kOffHi + 6620 + kQua + 6;
	  hN6620 << kHist;
	  vRC6620.push_back( new CsHist2D( hN6620.str(), hTitle,
					   12, 0., 12., 12, 0., 12. ) );
	}
        CsHistograms::SetCurrentPath("/");
      }
    }
    //std::cout << " RICHONE, " << bAMCorr << "  " << bPSCorr << std::endl;

    int cClu = clu->ic();
    std::list<CsRCPad*> lPads = clu->lPads();
    if( lPads.empty() ) {
      if( dump ) std::cout << "doPMTOptCorr - lPads empty" << std::endl;
      oflag = false;
      return  vDcPhoEmW;
//    -----------------
    }
//- For PMTs only one pad per cluster!
    CsRCPad* pad = lPads.front();

//- Compute the 'photon' incidence angles on detector :
//  ---------------------------------------------------
//- Assume photon emitted from average point on part. traj.
//- 'photon' impact on mirror (MWR) :
    CsRCMirrors *mirr = CsRCMirrors::Instance();
    Hep3Vector vPoC( 0., 0., 0. );
    Hep3Vector vPoPhoMir = mirr->vImpMir( vPoPhotW, vDcPhoEmW, vPoC, RR );
//                         ----------------------------------------------
//- Normal to mirror at 'photon' impact :
    Hep3Vector vDcNoPhoMir = (1./RR) * vPoPhoMir;
//- 'photon' reflected direction :
    double cosPhoMir = vDcNoPhoMir * vDcPhoEmW;
    Hep3Vector vDcPhoRefl = 2.*cosPhoMir * vDcNoPhoMir - vDcPhoEmW;
    Hep3Vector vDcPho = vDcPhoRefl.unit();
    Hep3Vector vDcPhow = vDcPho;
    if( bAMCorr ) {
//--- Opposite sign to cope Andreas def.
//    Uses impinging (reflected) photon direction
      vDcPhow.setX( -vDcPhow.x() );
      vDcPhow.setY( -vDcPhow.y() );
      vDcPhow.setZ( -vDcPhow.z() );
    }
    //std::cout << "doPMTOptCorr " << cClu << "  " << vDcPhow.y() <<std::endl;

    int PMTnum = pad->PMTnum();
    int PMTcha = pad->PMTcha();
    //std::cout << "doPMTOptCorr " << cClu << "  " << PMTnum
    //          << "  " << PMTcha << std::endl;

    int chaw = PMTcha;
    int iPadx = PMTcha%4;
    int iPady = int( PMTcha/4 );
    int kQua = -1;
    //std::cout << cClu << "  " << chaw << std::endl;
    if( dets->nCatPMT().size() <= 0 ) return  vDcPhoEmW;
//                                    -----------------
//- Mirroring for cathodes different from #03 (TopJura)
//  Optical Correction Table is for TopJura ONLY
//- TopJura
    if( cClu == dets->nCatPMT()[0] ) {
//--- Why iPady seems inverted (AM) ??? 
      if( bAMCorr ) chaw = (4-iPady-1)*4 + iPadx;            //AM
      //0chaw = iPady*4 + iPadx;
      if( bPSCorr ) chaw = (4-iPady-1)*4 + (4-iPadx-1);      //PS
      kQua = 0;
    }
//- TopSaleve
    else if( cClu == dets->nCatPMT()[1] ) {
      vDcPhow.setX( -vDcPhow.x() );
      if( bAMCorr ) chaw = (4-iPady-1)*4 + 4-iPadx-1;        //AM
      //0chaw = iPady*4 + 4-iPadx-1;
      if( bPSCorr ) chaw = (4-iPady-1)*4 + iPadx;            //PS
      kQua = 1;
    }
//- DownJura
    else if( cClu == dets->nCatPMT()[2] ) {
      vDcPhow.setY( -vDcPhow.y() );
      if( bAMCorr ) chaw = iPady*4 + iPadx;                  //AM
      //0chaw = (4-iPady-1)*4 + iPadx;
      if( bPSCorr ) chaw = iPady*4 + (4-iPadx-1);            //PS
      kQua = 2;
    }
//- DownSaleve
    else if( cClu == dets->nCatPMT()[3] ) {
      vDcPhow.setX( -vDcPhow.x() );
      vDcPhow.setY( -vDcPhow.y() );
      if( bAMCorr ) chaw = iPady*4 + 4-iPadx-1;              //AM
      //0chaw = (4-iPady-1)*4 + 4-iPadx-1;
      if( bPSCorr ) chaw = iPady*4 + iPadx;                  //PS
      kQua = 3;
    }
    //std::cout << cClu << "  " << vDcPhow.y() << std::endl;
    //std::cout << cClu << "  " << chaw << std::endl;

//test---obsolete----------------------------------------------
    //if( fabs( vDcPhow.x()-0.02 ) > 0.01 ) return  vDcPhoEmW;
    //if( fabs( vDcPhow.y()-0.35 ) > 0.02 ) return  vDcPhoEmW;
    //if( fabs( vDcPhow.x()+0.025 ) > 0.01 ) return  vDcPhoEmW;
    //if( fabs( vDcPhow.y()-0.30  ) > 0.01 ) return  vDcPhoEmW;
    //vDcPhow.setX( -0.025 );      // ave
    //vDcPhow.setY(  0.30  );      // ave
    //vDcPhow.setX( -0.087 );      // 5.   deg
    //vDcPhow.setY(  0.270 );      // 15.5 deg
    //vDcPhow.setX(  0.0   );      //
    //vDcPhow.setY(  0.350 );      //
    //vDcPhow.setY(  0.200 );      //
    //vDcPhow.setY(  0.0 );        //
//test---------------------------------------------------------

//- Monitoring histograms :
//  -----------------------
    if( kQua >= 0  &&  kQua < int( hist.vRC2120.size() ) ) {
      xh = float( pad->PMTcha() ) + 0.5;
      if( hist.vRC2120[kQua] ) hist.vRC2120[kQua]->Fill( xh );
//    -------------------------------------------------------
    }
    int kcj = kQua + 4;
    if( kcj >= 0  &&  kcj < int( hist.vRC2120.size() ) ) {
      xh = float( chaw ) + 0.5;
      if( hist.vRC2120[kcj] ) hist.vRC2120[kcj]->Fill( xh );
//    -----------------------------------------------------
    }

    double tgx = vDcPhow.x()/vDcPhow.z();
    double tgy = vDcPhow.y()/vDcPhow.z();
    int kck = kQua;
    if( kck >= 0  &&  kck < int( vRC6620.size() ) ) {
      xh = tgx;
      yh = tgy;
      if( vRC6620[kck] ) vRC6620[kck]->Fill( xh, yh );
//    -----------------------------------------------
    }
    kck = kQua + 4;
    if( kck >= 0  &&  kck < int( vRC6620.size() ) ) {
      xh = PMTnum%12 + 0.5;
      yh = int( PMTnum/12 ) + 0.5;
      if( vRC6620[kck] ) vRC6620[kck]->Fill( xh, yh );
//    -----------------------------------------------
    }

//test---obsolete------------------------
    //if( chaw !=  9 ) return  vDcPhoEmW;
    //if( chaw != 10 ) return  vDcPhoEmW;
    //if( chaw !=  6 ) return  vDcPhoEmW;
    //if( chaw !=  5 ) return  vDcPhoEmW;
    //chaw = 3;
    //double ww;
    //chaw = int( rndm_( ww )*15.999 );
//test-----------------------------------
    //std::cout << setprecision(7);
    //std::cout << "= " << chaw << "  " << vDcPhow.x() << "  " << vDcPhow.y();
    //std::cout << std::endl;

    double xCluC = 0.;
    double yCluC = 0.;
    double xCluCAM = 0.;
    double yCluCAM = 0.;
    if( bAMCorr ) {
      //CsRCOptCorr OptCorr;
      //if( !OptCorr.GetPhotonPosition( chaw, vDcPhow.x(), vDcPhow.y(),
      CsRCOptCorr* optcorr = CsRCOptCorr::Instance();
      if( !optcorr->GetPhotonPosition( chaw, vDcPhow.x(), vDcPhow.y(),
//    ----------------------------------------------------------------
          		 	       xCluCAM, yCluCAM ) ) {
	if( dump ) std::cout << "doPMTOptCorr - NO AM Corr  "
			     << nTot << std::endl;
        if( dump ) {
          //std::cout << std::endl;
          //std::cout << cClu << "  doPMTOptCorr - AMCorr - ";
          //std::cout << chaw << "  " << vDcPhow.x() << "  " << vDcPhow.y();
          //std::cout << std::endl;
        }

        oflag = false;
        return  vDcPhoEmW;
//      -----------------
      }

      //xCluCAM = - xCluCAM;      // 080930
//--- From cm -> mm !
      xCluCAM *= 10.;
      yCluCAM *= 10.;
      xCluC = xCluCAM;
      yCluC = yCluCAM;
      nCor++;

    }
    //if( nTot%100 == 0 ) {
    //std::cout << "- " << nCor << "  " << nTot << std::endl;
    //}
    //std::cout << vDcPhow.x() << "  " << vDcPhow.y();
    //std::cout << " - AM " << chaw << "  " << xCluC << "  " << yCluC;
    //std::cout << std::endl;

    double xCluCPS = 0.;
    double yCluCPS = 0.;
    if( bPSCorr ) {
      if( !GetPMTPhotonPosition( kQua, chaw, tgx, tgy, xCluCPS, yCluCPS ) ) {
//    -----------------------------------------------------------------------
	if( dump ) std::cout << "doPMTOptCorr - NO PS Corr  "
			     << nTot << std::endl;
        if( dump ) {
          //std::cout << std::endl;
          //std::cout << cClu << "  doPMTOptCorr - PSCorr - ";
          //std::cout << chaw << "  " << vDcPhow.x() << "  " << vDcPhow.y();
          //std::cout << std::endl;
        }
        oflag = false;
        return  vDcPhoEmW;
//      -----------------
      }
      //xCluCPS = - xCluCPS;      // 080930
      //yCluCPS = - yCluCPS;      // 080930
      xCluC = xCluCPS;
      yCluC = yCluCPS;
      nCor++;
      //std::cout << "  - PS   " << xCluC << "  " << yCluC << std::endl;
      //std::cout << "AM - PS   " << xh << "  " << yh << std::endl;
    }

    kcj = kQua + 8;
    if( kcj >= 0  &&  kcj < int( hist.vRC2120.size() ) ) {
      xh = float( chaw ) + 0.5;
      if( hist.vRC2120[kcj] ) hist.vRC2120[kcj]->Fill( xh );
//    -----------------------------------------------------
    }

//- Mirroring for cathodes different from #03 (TopJura)
//- TopJura
    if( cClu == dets->nCatPMT()[0] ) {
    }
//- TopSaleve
    else if( cClu == dets->nCatPMT()[1] ) {
      xCluC = -xCluC;
      xCluCAM = -xCluCAM;
      xCluCPS = -xCluCPS;
      //vDcPhow.setX( -vDcPhow.x() );
    }
//- DownJura
    else if( cClu == dets->nCatPMT()[2] ) {
      yCluC = -yCluC;
      yCluCAM = -yCluCAM;
      yCluCPS = -yCluCPS;
      //vDcPhow.setY( -vDcPhow.y() );
    }
//- DownSaleve
    else if( cClu == dets->nCatPMT()[3] ) {
      xCluC = -xCluC;
      xCluCAM = -xCluCAM;
      xCluCPS = -xCluCPS;
      yCluC = -yCluC;
      yCluCAM = -yCluCAM;
      yCluCPS = -yCluCPS;
      //vDcPhow.setX( -vDcPhow.x() );
      //vDcPhow.setY( -vDcPhow.y() );
    }

    CsRCCathode* cat = dets->ptrToCat( cClu );
    double xcePMT = 0.;
    double ycePMT = 0.;
    if( !cat->ccePMTWC( PMTnum, xcePMT, ycePMT ) ) {
//       ---------------------------------------
      if( dump ) std::cout << "doPMTOptCorr - NO ccePMT" << std::endl;
      oflag = false;
      return  vDcPhoEmW;
//    -----------------
    }
    double xCluCM = vPoCluDet.x() - xcePMT - dets->vOffCatW( cClu ).x();
    double yCluCM = vPoCluDet.y() - ycePMT - dets->vOffCatW( cClu ).y();

    bool print = false;
    //bool print = true;
    if( print ) {
      std::cout << vPoCluDet.x() << "  " << vPoCluDet.y() << "  "
	        << cClu << "  " << PMTnum << "  "
                << xcePMT << "  " << ycePMT << "  "
                << dets->vOffCatW( cClu ).x() << "  "
                << dets->vOffCatW( cClu ).y() << "  "
                << xCluCM << "  " << yCluCM << "  "
                << xCluCAM << "  " << yCluCAM << "  "
                << xCluCPS << "  " << yCluCPS << std:: endl;
    }

    double xDiff = 0.;
    double yDiff = 0.;
    if( bAMCorr ) {
      xDiff = xCluCAM - xCluCM;
      yDiff = yCluCAM - yCluCM;
    }
    if( bPSCorr ) {
      xDiff = xCluCPS - xCluCM;
      yDiff = yCluCPS - yCluCM;
      //if( fabs( xDiff ) > 12. ) std::cout << "X  " << xDiff << std::endl;
      if( fabs( xDiff ) > 12. ) nOutx++;
      //if( fabs( yDiff ) > 12. ) std::cout << "Y  " << yDiff << std::endl;
      if( fabs( yDiff ) > 12. ) nOuty++;
    }
    if( nTot%10000 == 0 ) {
      //std::cout << "doPMTOptCorr : " << nOutx << "  " << nOuty << "  "
      //	  << nTot << std::endl;
    }
    if( bAMCorr  &&   bPSCorr ) {
      xDiff = xCluCPS - xCluCAM;
      yDiff = yCluCPS - yCluCAM;
    }
    xh = xDiff;
    yh = yDiff;
    int kHH = 0;
    kHH = 0 + chaw*4;
    if( vRC3320[kHH] ) vRC3320[kHH]->Fill( tgx, tgy, xh );
//hh                   ----------------------------------
    kHH = 1 + chaw*4;
    if( vRC3320[kHH] ) vRC3320[kHH]->Fill( tgx, tgy, yh );
//hh                   ----------------------------------
    kHH = 2 + chaw*4;
    if( vRC3320[kHH] ) vRC3320[kHH]->Fill( tgx, tgy, 1. );
//hh                   ----------------------------------
    kHH = 3 + chaw*4;
    if( vRC3320[kHH] ) vRC3320[kHH]->Fill( tgx, tgy, 1. );
//hh                   ----------------------------------
    print = false;
    //print = true;
    if( print ) {
      double limit = 12.;
      if( fabs( xDiff ) > limit ) {
	if( fabs( xDiff ) > xDiffM ) xDiffM = fabs( xDiff );
        n12++;
        double fra = float( n12 )/float( nTot );
        std::cout << "X-Diff  " << xDiff << " - " << kQua << "  " << chaw
		  << setprecision(4)
                  << "  " << tgx << "  " << tgy << "  "
		  << fra << setprecision(2);
	if( tgx <= -0.06 || tgx >=  0.24 ) {
	  std::cout << " -- TGX / ";
	  nQx[kQua]++;
	}
	if( tgy <= -0.40 || tgy >= -0.10 ) {
	  std::cout << " -- TGY / ";
	  nQy[kQua]++;
	}
      }
      if( fabs( yDiff ) > limit ) {
	if( fabs( yDiff ) > yDiffM ) yDiffM = fabs( yDiff );
        n12++;
        double fra = float( n12 )/float( nTot );
        std::cout << "Y-Diff  " << yDiff << " - " << kQua << "  " << chaw
		  << setprecision(4)
                  << "  " << tgx << "  " << tgy << "  "
		  << fra << setprecision(2);
	if( tgx <= -0.06 || tgx >=  0.24 ) {
	  std::cout << " -- TGX / ";
	  nQx[kQua]++;
	}
	if( tgy <= -0.40 || tgy >= -0.10 ) {
	  std::cout << " -- TGY / ";
	  nQy[kQua]++;
	}
      }
      if( fabs( xDiff ) > limit || fabs( yDiff ) > limit ) {
        std::cout << "- " << nTot << " - ";
        for( int k=0; k<4; k++ ) std::cout << nQx[k] << " ";
        for( int k=0; k<4; k++ ) std::cout << nQy[k] << " ";
        std::cout << " - " << xDiffM << "  " << yDiffM;
        std::cout << std::endl;
        std::cout << vPoCluDet.x() << "  " << vPoCluDet.y() << "  "
	          << cClu << "  " << PMTnum << "  "
                  << xcePMT << "  " << ycePMT << "  "
                  << dets->vOffCatW( cClu ).x() << "  "
                  << dets->vOffCatW( cClu ).y() << "  "
                  << xCluCM << "  " << yCluCM << "  ";
	if( bAMCorr ) std::cout << xCluCAM << "  " << yCluCAM << "  ";
	if( bPSCorr ) std::cout << xCluCPS << "  " << yCluCPS;
	std::cout << std:: endl;
      }
    }

    double xCluCWC = xCluC + xcePMT + dets->vOffCatW( cClu ).x();
    double yCluCWC = yCluC + ycePMT + dets->vOffCatW( cClu ).y();
    Hep3Vector vPoCluCorr( xCluCWC, yCluCWC, vPoCluDet.z() );
    Hep3Vector CorrClu = vPoCluCorr - vPoCluDet;
    //if( fabs( CorrClu.x() ) > 12. ) {
    //if( fabs( CorrClu.y() ) > 12. ) {
    //std::cout << vPoCluDet << "  " << vPoCluCorr << "  "
    //          << cClu << "  " << nCor << std::endl;

    vPoCluDetO = vPoCluCorr;

    aflag = true;
    vDcPhoCorO = getCluAngle( vPoPhotW, vPoCluDetO, RR, aflag );
//               ----------------------------------------------
    if( !aflag ) {
      if( dump ) std::cout << "doPMTOptCorr - NO Clu Ang" << std::endl;
      return  vDcPhoEmW;
//    -----------------
    }

//test----------------------------------------------------------
/*
    int iPx = chaw%4;
    int iPy = int( chaw/4 );
    int chac = 0;
    if( cClu == dets->nCatPMT()[0] ) chac = (4-iPy-1)*4 + iPx;
    else if( cClu == dets->nCatPMT()[1] ) chac = (4-iPy-1)*4 + 4-iPx-1;
    else if( cClu == dets->nCatPMT()[2] ) chac = iPy*4 + iPx;
    else if( cClu == dets->nCatPMT()[3] ) chac = iPy*4 + 4-iPx-1;
    double xcePPad = cat->xcePPadPMT()[chac];
    double ycePPad = cat->ycePPadPMT()[chac];
    //std::cout << xCluC << "  " << yCluC << "  "
    //          << xcePPad << "  " << ycePPad << std::endl;
    CorrClu.setX( xCluC - xcePPad );
    CorrClu.setY( yCluC - ycePPad );
*/
//test----------------------------------------------------------

//- Monitoring histograms :
//  -----------------------
    double radDeg = cons->RadDeg();
    double da = 2.;
    double tgxm, tgxx, tgym, tgyx;
//  tgx > tg(5-2) and < 5+2 deg  and  tgy > -15.5-2 and < -15.5+2 deg / 0
    tgxm = tan( (5.-da)/radDeg );
    tgxx = tan( (5.+da)/radDeg );
    tgym = tan( (-15.5-da)/radDeg );
    tgyx = tan( (-15.5+da)/radDeg );
//  tgx > tg(-1-2) and < -1+2 deg  and  tgy > -20-2 and < -20+2 deg   / 1
    //tgxm = tan( (-1.-da)/radDeg );
    //tgxx = tan( (-1.+da)/radDeg );
    //tgym = tan( (-20.-da)/radDeg );
    //tgyx = tan( (-20.+da)/radDeg );
//  tgx > tg(11-2) and < 11+2 deg  and  tgy > -10-2 and < -10+2 deg   / 1
    //tgxm = tan( (11.-da)/radDeg );
    //tgxx = tan( (11.+da)/radDeg );
    //tgym = tan( (-10.-da)/radDeg );
    //tgyx = tan( (-10.+da)/radDeg );
    if( tgx > tgxm  &&  tgx < tgxx  &&  tgy > tgym  &&  tgy < tgyx ) {
    //if( PMTnum == 0 ) {
    //if( PMTnum == 6 ) {
    //if( PMTnum == 72 ) {
    //if( PMTnum == 78 ) {
    //if( PMTcha == 5 ) {
    //if( PMTcha == 9 ) {
    if( kQua >= 0  &&  kQua < int( hist.vRC2100.size() ) ) {
      //xh = vDcPho.x();
      xh = vDcPho.x()/vDcPho.z();
      //yh = CorrClu.x();
      yh = xDiff;
      if( hist.vRC2100[kQua] ) hist.vRC2100[kQua]->Fill( xh, yh );
//    -----------------------------------------------------------
    }
    if( kQua >= 0  &&  kQua < int( hist.vRC2105.size() ) ) {
      //xh = CorrClu.x();
      //yh = CorrClu.y();
      xh = xCluC;
      yh = yCluC;
      if( hist.vRC2105[kQua] ) hist.vRC2105[kQua]->Fill( xh, yh );
//    -----------------------------------------------------------
    }
    if( kQua >= 0  &&  kQua < int( hist.vRC2110.size() ) ) {
      //xh = vDcPho.y();
      xh = vDcPho.y()/vDcPho.z();
      //yh = CorrClu.y();
      yh = yDiff;
      if( hist.vRC2110[kQua] ) hist.vRC2110[kQua]->Fill( xh, yh );
//    -----------------------------------------------------------
    }
    }

    oflag = true;
//- return corrected photon position :
//  ----------------------------------
    vPoCluDet = vPoCluDetO;

//- return corrected photon direction :
//  -----------------------------------
    return  vDcPhoCorO;

  }


//===========================================================================
  bool CsRCPartPhotons::getHitCathodes( const double theIpo,
//----------------------------------------------------------
					int* kCInt, int* nCInt,
					double* aCInt ) {

//- Paolo  -  december 2006

    CsRCDetectors *dets = CsRCDetectors::Instance();
    CsRCRecConst *cons = CsRCRecConst::Instance();
    CsRCHistos& hist = CsRCHistos::Ref();
    double xh, yh, wh;

    static int nCathode = dets->nCathode();
    double xLimL[nCathode], xLimR[nCathode];
    double yLimD[nCathode], yLimU[nCathode];

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      CsRCExeKeys* key = CsRCExeKeys::Instance();
      key->acknoMethod( "CsRCPartPhotons::getHitCathodes" );

    }

    list<CsRCCathode*> lCathodes = dets->lCathodes();
    list<CsRCCathode*>::iterator ic;
    for( ic=lCathodes.begin(); ic!=lCathodes.end(); ic++ ) {
      CsRCCathode* cat = (*ic);
      int kc = cat->kCat();
      xLimL[kc] = cat->xLimMnW();
      xLimR[kc] = cat->xLimMxW();
      yLimD[kc] = cat->yLimMnW();
      yLimU[kc] = cat->yLimMxW();
    }

    CsRCParticle* part = pPart_;
    int kDetPart = part->kDetPart();
    double partX = part->vPade()[kDetPart].x();
    double partY = part->vPade()[kDetPart].y();

/*
//- get cathode 'hit' by the 'reflected' particle
    int kCat = -1;
    CsRCCathode* cat = NULL;
    for( ic=lCathodes.begin(); ic!=lCathodes.end(); ic++ ) {
      int kc = (*ic)->kCat();
      if( partX >= xLimL[kc]  &&  partX <= xLimR[kc]  &&
          partY >= yLimD[kc]  &&  partY <= yLimU[kc] ) {
	kCat = kc;
	cat = (*ic);
	break;
      }
    }
    if( kCat < 0 ) {
//--- discard 'reflected' particles outside PD.s
      if( partX < xDetLimits( 0 )[0]  ||  partX > xDetLimits( 0 )[1] ||
	  partX < xDetLimits( 1 )[0]  ||  partX > xDetLimits( 1 )[1] )
	return  false;
//      -------------
      if( partY < yDetLimits( 0 )[0]  ||  partY > yDetLimits( 0 )[1] ||
	  partY < yDetLimits( 1 )[0]  ||  partY > yDetLimits( 1 )[1] )
	return  false;
//      -------------
//--- assign cath to 'reflected' particles 'hits' between adj. cathodes
      double disttMn = 1000000.;
      for( ic=lCathodes.begin(); ic!=lCathodes.end(); ic++ ) {
	int kc = (*ic)->kCat();
	double distt = pow( partX - (*ic)->vOffCatW().x(), 2 ) +
	  pow( partY - (*ic)->vOffCatW().y(), 2 );
	if( distt < disttMn ) {
	  disttMn = distt;
	  kCat = kc;
	}
      }
    }
*/

//- get cathode closest to the 'reflected' particle 'hit'
    int kCat = part->getClosestCat();
    //if( iCaPa_ != kCat ) std::cout << iCaPa_ << "  " << kCat << std::endl;
    if( kCat < 0 ) {
      std::cout << "getHitCathodes ret 0  " << std::endl;
      return  false;
//    -------------
    }

    double xxc[4];
    double yyc[4];

    double theIpow = theIpo;
    double Rad = theIpow/1000.* part->pMirPart()->RR()/2.;   //   approximate!
    double Rad2 = Rad * Rad;
    double delta = 0.;
    double xInt, yInt;
    for( ic=lCathodes.begin(); ic!=lCathodes.end(); ic++ ) {
      int kc = (*ic)->kCat();
      kCInt[kc] = -1;
      nCInt[kc] = 0;
      aCInt[kc] = 0.;
    }
    int kCc = 0;
    int kCx = 0;
    double yLim = yLimD[kCat];
    for( int kk=0; kk<2; kk++ ) {
      delta = Rad2 - (yLim-partY)*(yLim-partY);
      if( delta > 0. ) {
	delta = sqrt( delta );
	for( int ks=-1; ks<=+1; ks+=2 ) {
	  xInt = partX + ks * delta;
	  if( xInt >= xLimL[kCat]  &&  xInt <= xLimR[kCat] ) {
	    xxc[kCc] = xInt;
	    yyc[kCc] = yLim;
	    kCc++;
	    kCx++;
	  }
	}
      }
      yLim = yLimU[kCat];
    }
    int kCy = 0;
    double xLim = xLimL[kCat];
    for( int kk=0; kk<2; kk++ ) {
      delta = Rad2 - (xLim-partX)*(xLim-partX);
      if( delta > 0. ) {
	delta = sqrt( delta );
	for( int ks=-1; ks<=+1; ks+=2 ) {
	  yInt = partY + ks * delta;
	  if( yInt >= yLimD[kCat]  &&  yInt <= yLimU[kCat] ) {
	    xxc[kCc] = xLim;
	    yyc[kCc] = yInt;
	    kCc++;
	    kCy++;
	  }
	}
      }
      xLim = xLimR[kCat];
    }
    //std::cout << "getHitCathodes 0c  " << kCx << "  " << kCy << std::endl;

    int kInt = kCx + kCy;
    if( kInt%2 != 0 ) {
      std::cout << "getHitCathodes ret 1a  " << kCx << " " << kCy <<std::endl;
      return  false;
//    -------------
    }
    if( kInt > 4 )  {
      std::cout << "getHitCathodes ret 1b  " << kInt << std::endl;
      return  false;
//    -------------
    }

    int jCat = 0;
    double delTo = 0.;
    kCInt[kCat] = jCat;
    nCInt[kCat] = kInt;
    aCInt[kCat] = 0.;
    if( kInt == 0 ) {
      double del = cons->TwoPI();
      del -= 0.008;
      delTo  = del;
      aCInt[kCat] = del;
      if( likeONLY_ ) {
        xh = del * cons->RadDeg();
	yh = 0.5;
	if( hist.hRC1623 ) hist.hRC1623->Fill( xh, yh );
//      -----------------------------------------------
	xh = delTo * cons->RadDeg();
	yh = 4.5;
	if( hist.hRC1623 ) hist.hRC1623->Fill( xh, yh );
//      -----------------------------------------------
      }
      return  true;
//    ------------
    }
    else if( kInt == 2 ) {
      double del = acos( ((xxc[0]-partX)*(xxc[1]-partX) +
			  (yyc[0]-partY)*(yyc[1]-partY) ) / Rad2 );
      del = cons->TwoPI() - del;
      delTo += del;
      aCInt[kCat] = del;
      if( likeONLY_ ) {
	xh = del * cons->RadDeg();
	yh = 1.5;
	if( hist.hRC1623 ) hist.hRC1623->Fill( xh, yh );
//      -----------------------------------------------
      }
    }
    else if( kInt == 4 ) {
      double del1 = acos( ((xxc[0]-partX)*(xxc[1]-partX) + 
			   (yyc[0]-partY)*(yyc[1]-partY) ) / Rad2 );
      double del2 = acos( ((xxc[2]-partX)*(xxc[3]-partX) + 
		 	   (yyc[2]-partY)*(yyc[3]-partY) ) / Rad2 );
      double del = cons->TwoPI() - del1 - del2;
      delTo += del;
      aCInt[kCat] = del;
      if( likeONLY_ ) {
	xh = del * cons->RadDeg();
	yh = 2.5;
	if( hist.hRC1623 ) hist.hRC1623->Fill( xh, yh );
//      -----------------------------------------------
      }
    }
    int kInt0 = kInt;

    int kIntTo = kInt0;
    for( ic=lCathodes.begin(); ic!=lCathodes.end(); ic++ ) {
      int kc = (*ic)->kCat();
      if( kc == kCat ) continue;
      int kCc = 0;
      int kCx = 0;
      double yLim = yLimD[kc];
      for( int kk=0; kk<2; kk++ ) {
	delta = Rad2 - (yLim-partY)*(yLim-partY);
        if( delta > 0. ) {
	  delta = sqrt( delta );
	  for( int ks=-1; ks<=+1; ks+=2 ) {
	    xInt = partX + ks * delta;
	    if( xInt >= xLimL[kc]  &&  xInt <= xLimR[kc] ) {
	      xxc[kCc] = xInt;
	      yyc[kCc] = yLim;
	      kCc++;
	      kCx++;
	    }
	  }
	}
	yLim = yLimU[kc];
      }
      int kCy = 0;
      double xLim = xLimL[kc];
      for( int kk=0; kk<2; kk++ ) {
        delta = Rad2 - (xLim-partX)*(xLim-partX);
        if( delta > 0. ) {
	  delta = sqrt( delta );
	  for( int ks=-1; ks<=+1; ks+=2 ) {
	    yInt = partY + ks * delta;
	    if( yInt >= yLimD[kc]  &&  yInt <= yLimU[kc] ) {
	      xxc[kCc] = xLim;
	      yyc[kCc] = yInt;
	      kCc++;
	      kCy++;
	    }
	  }
	}
	xLim = xLimR[kc];
      }

      int kInt = kCx + kCy;
      if( kInt%2 != 0 ) {
	std::cout << "getHitCathodes ret 2a  " << std::endl;
	return  false;
//      -------------
      }
      if( kInt > 2 ) {
	std::cout << "getHitCathodes ret 2b  " << std::endl;
	return  false;
//      -------------
      }

      if( kInt == 2 ) {
	jCat++;
	kCInt[kc] = jCat;
	nCInt[kc] = kInt;
        double del = acos( ((xxc[0]-partX)*(xxc[1]-partX) +
			    (yyc[0]-partY)*(yyc[1]-partY) ) / Rad2 );
	delTo += del;
	aCInt[kc] = del;
	if( likeONLY_ ) {
	  xh = del * cons->RadDeg();
	  yh = 3.5;
	  if( hist.hRC1623 ) hist.hRC1623->Fill( xh, yh );
//        -----------------------------------------------
	  xh = delTo * cons->RadDeg();
	  yh = 4.5;
	  if( hist.hRC1623 ) hist.hRC1623->Fill( xh, yh );
//        -----------------------------------------------
	}
      }
      kIntTo += kInt;
      //if( likeONLY_ ) {
      //  xh = delTo * cons->RadDeg();
      //  yh = 4.5;
      //  if( hist.hRC1623 ) hist.hRC1623->Fill( xh, yh );
//      -----------------------------------------------
      //}
    }
    if( kIntTo == 0 ) {
      std::cout << "getHitCathodes ret 3  " << std::endl;
      return  false;
//    -------------
    }

    return  true;
  }


//===========================================================================
  std::vector<double> CsRCPartPhotons::getSplitLimits() {
//-------------------------------------------------------


//- Paolo  -  January 2007
//  from existing getSplitLimits( p )


    CsRCDetectors *dets = CsRCDetectors::Instance();
    CsRCMirrors *mirr = CsRCMirrors::Instance();
    CsRCHistos& hist = CsRCHistos::Ref();
    double xh, yh, wh;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;
    }

    CsRCRecConst *cons = CsRCRecConst::Instance();
    double tgCmax = tan( acos( 1./cons->CFRefIndUV() ) );   //   ???

    CsRCParticle* part = pPart_;
    Hep3Vector vPosIn0 = part->vPosIn();
    Hep3Vector vDirIn0 = part->vDirIn();
    double tgb = vDirIn0.y()/vDirIn0.z();
    double tgbl = 0.;
    if( tgb > 0. ) tgbl = (tgb-tgCmax) / (1+tgb*tgCmax);
    if( tgb < 0. ) tgbl = (tgb+tgCmax) / (1-tgb*tgCmax);
    double tga = vDirIn0.x()/vDirIn0.z();
    Hep3Vector vDcPaL( tga, tgbl, 1.);
    vDcPaL = vDcPaL.unit();
    Hep3Vector vPoC0 = part->pMirPart()->vC0();
    double RR = part->pMirPart()->RR();
    Hep3Vector vPoMirL = mirr->vImpMir( vPosIn0, vDcPaL, vPoC0, RR );

    std::vector<double> vySplitLLimit;
    vySplitLLimit.push_back( yDetLimits( 0 )[1] );
    vySplitLLimit.push_back( yDetLimits( 1 )[0] );
    bool exe = false;
    if( kDetPart_ == 0  &&  vPoMirL.y() < 0. ) exe = true;
    if( kDetPart_ == 1  &&  vPoMirL.y() > 0. ) exe = true;
    if( bSplitLLimSet_ ) exe = false;
    if( exe ) {
      vySplitLLimit.clear();

        double tgap = (tga+tgCmax)/(1-tga*tgCmax);
        double tgam = (tga-tgCmax)/(1+tga*tgCmax);
        Hep3Vector vPoPaC( vPosIn0.x(), 0., vPosIn0.z() );
        Hep3Vector vDcPaC( vDirIn0.x(), 0., vDirIn0.z() );
	vDcPaC = vDcPaC.unit();
	Hep3Vector vPoMirC = mirr->vImpMir( vPoPaC, vDcPaC, vPoC0, RR );

	int kDetPaw = 0;
	list<CsRCPhotonDet*> lPhoDet = dets->lPhoDet();
	list<CsRCPhotonDet*>::iterator id;
	std::list<CsRCMirrorNom*> lMirrNom = mirr->lMirrNom();
	std::list<CsRCMirrorNom*>::iterator im;
	im = lMirrNom.begin();
	for( id=lPhoDet.begin(); id!=lPhoDet.end(); id++ ) {

	  Hep3Vector vPoC0 = (*im)->vC0();
	  double RR = (*im)->RR();
	  im++;

          Hep3Vector vDcNor(0., 0., 0.);
          double cosMir = 0.;
	  Hep3Vector vDcRefl(0., 0., 0.);
	  double norm = 0.;
  	  Hep3Vector vDcDet = dets->vDcDetv( kDetPaw );
	  Hep3Vector vPoDet = dets->vDet0v( kDetPaw );

	  Hep3Vector vDcPaP( tgap, 0., 1.);
	  vDcPaP = vDcPaP.unit();
	  Hep3Vector vPoMirP = mirr->vImpMir( vPoPaC, vDcPaP, vPoC0, RR );
	  Hep3Vector vDcMirP = (vPoMirP - vPosIn0).unit();
          vDcNor = (1./RR) * (vPoMirP - vPoC0);
          cosMir = vDcNor * vDcMirP;
          vDcRefl = 2.*cosMir * vDcNor - vDcMirP;
          norm = ( ( vPoDet - vPoMirP ) * vDcDet ) / ( vDcRefl * vDcDet );
          Hep3Vector vPoDetP = vPoMirP + norm * vDcRefl;

	  Hep3Vector vDcPaM( tgam, 0., 1.);
	  vDcPaM = vDcPaM.unit();
	  Hep3Vector vPoMirM = mirr->vImpMir( vPoPaC, vDcPaM, vPoC0, RR );
	  Hep3Vector vDcMirM = (vPoMirM - vPosIn0).unit();
          vDcNor = (1./RR) * (vPoMirM - vPoC0);
          cosMir = vDcNor * vDcMirM;
          vDcRefl = 2.*cosMir * vDcNor - vDcMirM;
          norm = ( ( vPoDet - vPoMirM ) * vDcDet ) / ( vDcRefl * vDcDet );
          Hep3Vector vPoDetM = vPoMirM + norm * vDcRefl;
	  //cout << vPoMirP << "  " << vPoMirC << "  " << vPoMirM << endl;
	  //cout << vPoDetP << "  " << vPoDetM << endl;

  	  HepVector vrot( 3, 0 );
	  for( int j=0; j<3; j++ ) vrot[j] = (vPoDetP - vPoC0)[j];
	  HepVector vPoDetPWw = (*id)->rotMatrix() * vrot;
	  Hep3Vector vPoDetPW( vPoDetPWw[0], vPoDetPWw[1], vPoDetPWw[2] );
	  for( int j=0; j<3; j++ ) vrot[j] = (vPoDetM - vPoC0)[j];
          HepVector vPoDetMWw = (*id)->rotMatrix() * vrot;
	  Hep3Vector vPoDetMW( vPoDetMWw[0], vPoDetMWw[1], vPoDetMWw[2] );
	  double yDetLL = 0.;
	  double splitLmin = 0.;
	  if( kDetPaw == 0 ) {
	    splitLmin = vPoDetPW.y();
	    if( splitLmin > vPoDetMW.y() ) splitLmin = vPoDetMW.y();
	    yDetLL = yDetLimits( kDetPaw )[1];
	    if( splitLmin < yDetLL ) splitLmin = yDetLL;
	  }
	  if( kDetPaw == 1 ) {
	    splitLmin = vPoDetPW.y();
	    if( splitLmin < vPoDetMW.y() ) splitLmin = vPoDetMW.y();
	    yDetLL = yDetLimits( kDetPaw )[0];
	    if( splitLmin > yDetLL ) splitLmin = yDetLL;

	  }
	  vySplitLLimit.push_back( splitLmin );

	  //cout << vPoMirL.y() << "  " << kDetPaw << " - " 
	  //     << vPoPaDetW_[kDetPaw] << " P "
	  //     << vPoDetPW << " M " << vPoDetMW << endl;
          kDetPaw++;
	}
	vySplitLLim_.push_back( vySplitLLimit[0] );
	vySplitLLim_.push_back( vySplitLLimit[1] );
	bSplitLLimSet_ = true;
	//std::cout << kDetPart_ << "  " << vPoPaDetW_[kDetPart_].y()
	//	  << "  " << vySplitLLimit[0] << "  " << vySplitLLimit[1] 
	//	  << std::endl;
      }

//--- monitor the limits :
//    --------------------
      xh = vySplitLLimit[0];
      yh = vySplitLLimit[1];
      if( hist.hRC3508 ) hist.hRC3508->Fill( xh, yh );
//hh                     ----------------------------
      xh = vPoPaDetW_[kDetPart_].y();
      yh = vPoPaDetW_[1-kDetPart_].y();
      if( hist.hRC3509 ) hist.hRC3509->Fill( xh, yh );
//hh                     ----------------------------

      return vySplitLLimit;

    }


//===========================================================================
  std::vector<double> CsRCPartPhotons::xDetLimits( const int kDet ) {
//-------------------------------------------------------------------


//- Paolo  -  January 2007


    static std::vector<double> xxLimDetUp;
    static std::vector<double> xxLimDetDw;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      CsRCDetectors *dets = CsRCDetectors::Instance();
      list<CsRCCathode*> lCathodes = dets->lCathodes();
      double hCatx = 0.;
      float xxLimMinU =  1000000.;
      float xxLimMaxU = -1000000.;
      float xxLimMinD =  1000000.;
      float xxLimMaxD = -1000000.;
      int kCat = -1;
      list<CsRCCathode*>::iterator ic;
      for( ic=lCathodes.begin(); ic!=lCathodes.end(); ic++ ) {
	float xxOff = (*ic)->vOffCatW().x();
	if( xxOff > 0. ) {
	  if( xxOff < xxLimMinU ) xxLimMinU = xxOff - (*ic)->hCatx();
	  if( xxOff > xxLimMaxU ) xxLimMaxU = xxOff + (*ic)->hCatx();
	}
	if( xxOff < 0. ) {
	  if( xxOff < xxLimMinD ) xxLimMinD = xxOff - (*ic)->hCatx();
	  if( xxOff > xxLimMaxD ) xxLimMaxD = xxOff + (*ic)->hCatx();
	}
      }
      xxLimDetUp.push_back( xxLimMaxU );
      xxLimDetUp.push_back( xxLimMinU );
      xxLimDetDw.push_back( xxLimMaxD );
      xxLimDetDw.push_back( xxLimMinD );

    }

    if( kDet == 0 ) return  xxLimDetUp;
    else  return  xxLimDetDw;

  }


//===========================================================================
  std::vector<double> CsRCPartPhotons::yDetLimits( const int kDet ) {
//-------------------------------------------------------------------


//- Paolo  -  January 2007


    static std::vector<double> yyLimDetUp;
    static std::vector<double> yyLimDetDw;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      CsRCDetectors *dets = CsRCDetectors::Instance();
      list<CsRCCathode*> lCathodes = dets->lCathodes();
      double hCaty = 0.;
      float yyLimMinU =  1000000.;
      float yyLimMaxU = -1000000.;
      float yyLimMinD =  1000000.;
      float yyLimMaxD = -1000000.;
      int kCat = -1;
      list<CsRCCathode*>::iterator ic;
      for( ic=lCathodes.begin(); ic!=lCathodes.end(); ic++ ) {
	float yyOff = (*ic)->vOffCatW().y();
	if( yyOff > 0. ) {
	  if( yyOff < yyLimMinU ) yyLimMinU = yyOff - (*ic)->hCaty();
	  if( yyOff > yyLimMaxU ) yyLimMaxU = yyOff + (*ic)->hCaty();
	}
	if( yyOff < 0. ) {
	  if( yyOff < yyLimMinD ) yyLimMinD = yyOff - (*ic)->hCaty();
	  if( yyOff > yyLimMaxD ) yyLimMaxD = yyOff + (*ic)->hCaty();
	}
      }
      yyLimDetUp.push_back( yyLimMaxU );
      yyLimDetUp.push_back( yyLimMinU );
      yyLimDetDw.push_back( yyLimMaxD );
      yyLimDetDw.push_back( yyLimMinD );

    }

    if( kDet == 0 ) return  yyLimDetUp;
    else  return  yyLimDetDw;

  }


//===========================================================================
  bool CsRCPartPhotons::getSplitCondition() {
//-------------------------------------------

//- Paolo  -  January 2007

    bool split = false;

    double tanCmax = tan( acos( 1./CsRCRecConst::Instance()->CFRefIndUV() ) );
    double rCMaxMirr = pPart_->pathLen() * tanCmax;
    double hhSplit = pPart_->vPoPaMir0()[pPart_->kDetPart()].y();
    //std::cout << hhSplit << "  " << rCMaxMirr << std::endl;

    if( fabs( hhSplit ) < rCMaxMirr ) split = true;

    return split;
  }



//========================================================================
  double CsRCPartPhotons::thetaVStoUV( const double the ) const {
//---------------------------------------------------------------

//- Paolo
//- January  2007

//- angles in mrad

    CsRCRecConst *cons = CsRCRecConst::Instance();

    double beta1 = cons->CFRefIndVS() * cos( the/1000. );
    double cos = 1.;
    if( cons->CFRefIndUV() > 0.) cos = beta1 / cons->CFRefIndUV();
    if( cos > 1.) cos = 1.;
    double theta = 1000. * acos( cos );
    //std::cout << "thetaVStoUV " << the << "  " << theta << std::endl;

    return  theta;

  }


//========================================================================
  double CsRCPartPhotons::thetaUVtoVS( const double the ) const {
//---------------------------------------------------------------

//- Paolo
//- January  2007

//- angles in mrad

    CsRCRecConst *cons = CsRCRecConst::Instance();

    double beta1 = cons->CFRefIndUV() * cos( the/1000. );
    double cos = 1.;
    if( cons->CFRefIndVS() > 0.) cos = beta1 / cons->CFRefIndVS();
    if( cos > 1.) cos = 1.;
    double theta = 1000. * acos( cos );
    //std::cout << "thetaUVtoVS " << the << "  " << theta << std::endl;

    return  theta;

  }


//===========================================================================
  double CsRCPartPhotons::nPhotExpctCat( double theIpo ) {
//--------------------------------------------------------


//- Paolo  -  January 2007

//  WARNING : SPLIT RING part still MISSING!


    CsRCRecConst *cons = CsRCRecConst::Instance();
    static const int kOffHi = cons->kOffHi();
    static const double TwoPi = cons->TwoPI();
    CsRCDetectors* dets = CsRCDetectors::Instance();
    static int nCathode = dets->nCathode();

    static const float focLen = pPart_->pMirPart()->RR() / 2.;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      CsRCExeKeys* key = CsRCExeKeys::Instance();
      key->acknoMethod( "CsRCPartPhotons::nPhotExpctCat" );

    }

    double nPhotEx = 0;

    double sTheIrad2 = pow( sin( theIpo / 1000.), 2 );
    float pathLen = pPart_->pathLen();
    double nPhotExUV = cons->nZero() * pathLen/10. * sTheIrad2;
//------------------------------------------------------------
    double nPhotExVS = cons->nZeroVS() * pathLen/10. * sTheIrad2;
//--------------------------------------------------------------
    //double rTheIpo = theIrad * focLen;
    int kDetPart = pPart_->kDetPart();
    Hep3Vector vPade = pPart_->vPade()[kDetPart];

    int kCatHit[nCathode];
    int nCatHit[nCathode];
    double aCatHit[nCathode];
    bool catHit = false;
    if( theIpo > 0. ) catHit = getHitCathodes( theIpo, &kCatHit[0],
//                    ---------------------------------------------
			       &nCatHit[0], &aCatHit[0] );
    if( !catHit ) return  0.;
//                ----------
/*
std::cout << likeONLY_ << std::endl;
if( likeONLY_ ) {
std::cout << "nPhotExpctCat  ";
double aTo = 0.;
for( int k=0; k<16; k++ ) {
  if( kCatHit[k] >= 0. ) {
std::cout << kCatHit[k] << "  " << aCatHit[k]*180./3.14159 << "   ";
    aTo += aCatHit[k];
  }
}
std::cout << " -   " << aTo*180./3.14159 << std::endl;
}
*/

    list<CsRCCathode*> lCathodes = dets->lCathodes();
    list<CsRCCathode*>::iterator ic;
    int kCat = -1;
    for( int kc=0; kc<nCathode; kc++ ) {
      if( kCatHit[kc] == 0 ) {
	kCat = kc;
	break;
      }
    }
    if( kCat < 0 ) return  0.;
//                 ----------
    double nPhotExCat[nCathode];
    double angTo = 0.;
    double nPhotTo = 0;
    double nPhotIn  = 0;
    double nPhotOut = 0;
//- thetaIpo is UV or VS...
    double sTheSatUV2 = pow( sin( acos( 1./cons->CFRefIndUV() ) ), 2 );
    double sTheSatVS2 = pow( sin( acos( 1./cons->CFRefIndVS() ) ), 2 );
    double frcUV = 0.;
    double frcVS = 0.;
    for( ic=lCathodes.begin(); ic!=lCathodes.end(); ic++ ) {
      int kc = (*ic)->kCat();
      nPhotExCat[kc] = 0;
      if( kCatHit[kc] >= 0. ) {
	double nPhotExw = nPhotExUV;
	if( (*ic)->isPMT() ) nPhotExw = nPhotExVS;
	if( cons->bPhotAtSat() ) {
	  double sTheSatw2 = sTheSatUV2;
	  if( (*ic)->isPMT() ) sTheSatw2 = sTheSatVS2;
	  if( sTheIrad2 > sTheSatw2 ) sTheIrad2 = sTheSatw2;
          nPhotExw = cons->nPhotAtSat()[kc] * sTheIrad2 / sTheSatw2;
	}
	double nPhot = nPhotExw * aCatHit[kc] / cons->TwoPI();
	nPhotExCat[kc] = nPhot;
	nPhotTo += nPhot;
	if( kc == kCat) nPhotIn   = nPhot;
	if( kc != kCat) nPhotOut += nPhot;
	angTo += aCatHit[kc];
	if( (*ic)->isAPV() ) frcUV += aCatHit[kc] / cons->TwoPI();
	if( (*ic)->isPMT() ) frcVS += aCatHit[kc] / cons->TwoPI();
      }
    }
    nPhotEx = nPhotTo;


//- Histogramms :
//  -------------
    CsRCHistos& hist = CsRCHistos::Ref();
    int level = hist.levelBk();
    bool doHist = false;
    if ( level >= 2 ) doHist = true;
    if( doHist ) {
      if( likeONLY_  &&  catHit ) {
        double xh, yh, wh;
        xh = vPade.x();
        yh = vPade.y();
        wh = double( nPhotEx );
        if( hist.hRC1505 ) hist.hRC1505->Fill( xh, yh, wh );
//      ---------------------------------------------------
        xh = vPade.x();
        yh = vPade.y();
        if( hist.hRC1506 ) hist.hRC1506->Fill( xh, yh );
//      -----------------------------------------------
	if( dets->ptrToCat( kCat )->isPMT() ) {
	  xh = frcUV;
	  if( hist.hRC1507 ) hist.hRC1507->Fill( xh );
//        -------------------------------------------
	  xh = nPhotEx;
	  if( hist.hRC1572 ) hist.hRC1572->Fill( xh );
//        -------------------------------------------
	} else {
	  xh = frcVS;
	  if( hist.hRC1508 ) hist.hRC1508->Fill( xh );
//        -------------------------------------------
	  xh = nPhotEx;
	  if( hist.hRC1573 ) hist.hRC1573->Fill( xh );
//        -------------------------------------------
	}
	if( pionONLY_ ) {
  	  for( ic=lCathodes.begin(); ic!=lCathodes.end(); ic++ ) {
	    int kc = (*ic)->kCat();
	    if( kCatHit[kc] < 0 ) continue;
	    if( kc < 0  ||  kc >= int( hist.vRC1350.size() ) ) continue;
	    yh = -0.5;
	    if( kc == kCat ) {
	      //if( kCatHit[kc] == 0 ) {
	      xh = nPhotExCat[kc];
  	      if( nCatHit[kc] == 0 ) yh = 0.5;
  	      if( nCatHit[kc] == 2 ) yh = 1.5;
  	      if( nCatHit[kc] == 4 ) yh = 2.5;
	      if( hist.vRC1350[kc] ) hist.vRC1350[kc]->Fill( xh, yh );
//            -------------------------------------------------------
	      xh = nPhotOut;
	      yh = 4.5;
	      if( hist.vRC1350[kc] ) hist.vRC1350[kc]->Fill( xh, yh );
//            -------------------------------------------------------
	      xh = nPhotEx;
	      yh = 5.5;
	      if( hist.vRC1350[kc] ) hist.vRC1350[kc]->Fill( xh, yh );
//            -------------------------------------------------------
	      xh = nPhotIn;
	      yh = nPhotOut;
	      if( hist.vRC1370[kc] ) hist.vRC1370[kc]->Fill( xh, yh );
//            -------------------------------------------------------
	    } else {
	      xh = nPhotExCat[kc];
  	      if( nCatHit[kc] == 2 ) yh = 3.5;
	      if( hist.vRC1350[kc] ) hist.vRC1350[kc]->Fill( xh, yh );
//            -------------------------------------------------------
	    }
	  }
	}
      }
    }

    return  nPhotEx;
  }


//===========================================================================
  double CsRCPartPhotons::getThetaWgav() {
//----------------------------------------


//- theta from weighted average
//  ---------------------------
//- Paolo  -  March 2007


    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;
      CsRCExeKeys::Instance()->acknoMethod( "CsRCPartPhotons::getThetaWgav" );
    }

    CsRCHistos& hist = CsRCHistos::Ref();
    double xh, yh, wh;

    CsRCRecConst * cons = CsRCRecConst::Instance();

    list<CsRCPhoton*> lPhotons = lPhotons_;
    list<CsRCPhoton*>::iterator ih;
    double coSig = 0.;
    double coBkg = 0.;
    for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
      double thePhot = (*ih)->the();
      double sigPhot = (*ih)->sigmaPhoPid( this );
//@@                   --------------------------
      coSig += thePhot/(sigPhot*sigPhot);
      coBkg += 1. /(sigPhot*sigPhot);
    }
    double thetaWgav = 0.;
    if( coBkg > 0. ) thetaWgav = coSig / coBkg;

    return  thetaWgav;

  }


//======================================================================
  Hep3Vector CsRCPartPhotons::doQzWCorr( const Hep3Vector vPoPhotW,
//----------------------------------------------------------------------
					 const Hep3Vector vDcPhoEmW,
					 Hep3Vector& vPoCluDet,
					 const float RR, const float zDetW,
					 bool& aflag ) {


//--- approximate correction for detector quartz window :
//    ---------------------------------------------------
//    new form 15/06/07

      CsRCExeKeys *key = CsRCExeKeys::Instance();


      static bool firstCall = true;
      if( firstCall ) {
        firstCall = false;

        key->acknoMethod( "CsRCPartPhotons::QuartzWindowCorr" );
      }

//--- from "ReconConst.h"
      CsRCRecConst *cons = CsRCRecConst::Instance();
      float CFRefInd = cons->CFRefInd();
      float qzRefInd = cons->qzRefInd();

      CsRCDetectors *dets = CsRCDetectors::Instance();
      list<CsRCCathode*> lCathodes = dets->lCathodes();
      CsRCCathode *cat  = lCathodes.front();
      static float ddQzW = cat->ddQzW();
      static float ddGap = cat->ddGap();

      static double zQzWOutside = zDetW + ddGap + ddQzW;
      static double zQzWinside  = zDetW + ddGap;

      static double rRatio = CFRefInd / qzRefInd;

//--- assume photon emitted from average point on part. traj.
//    -------------------------------------------------------
//--- 'photon' impact on mirror (MWR) :
//    ---------------------------------
      CsRCMirrors *mirr = CsRCMirrors::Instance();
      Hep3Vector vPoC( 0., 0., 0. );
      Hep3Vector vPoPhoMir = mirr->vImpMir( vPoPhotW, vDcPhoEmW, vPoC, RR );
//                           ----------------------------------------------

//--- normal to mirror at 'photon' impact :
//    -------------------------------------
      Hep3Vector vDcNoPhoMir = (1./RR) * vPoPhoMir;

//--- 'photon' reflected direction :
//    ------------------------------
      double cosPhoMir = vDcNoPhoMir * vDcPhoEmW;
      Hep3Vector vDcPhoRefl = 2.*cosPhoMir * vDcNoPhoMir - vDcPhoEmW;
//pp      cout << kMir << "  vDcPhoRefl  " << vDcPhoRefl << endl;

//--- 'photon' incidence angle on detector :
//    --------------------------------------
      double cosPhoDet = vDcPhoRefl.z();
      double sinPhoDet = sqrt(1.-cosPhoDet*cosPhoDet);

      double norm;
//--- 'photon' impact on quartz win. (outside) :
//    ------------------------------------------
      norm = (zQzWOutside - vPoPhoMir.z()) / vDcPhoRefl.z();
      Hep3Vector vPhoWinOut = vPoPhoMir + norm * vDcPhoRefl;

//--- 'photon' projection on quartz win. (inside) :
//    ---------------------------------------------
      norm = (zQzWinside - vPoPhoMir.z()) / vDcPhoRefl.z();
      Hep3Vector vPhoWinIn = vPoPhoMir + norm * vDcPhoRefl;

//--- 'photon' direction inside quartz win. (DRS) :
//    ---------------------------------------------
      double sinPhoQzR = sinPhoDet * rRatio;
      double cosPhoQzR = sqrt(1.-sinPhoQzR*sinPhoQzR);
      Hep3Vector vDcPhoQzR;
      vDcPhoQzR.setX( rRatio * vDcPhoRefl.x() );
      vDcPhoQzR.setY( rRatio * vDcPhoRefl.y() );
      vDcPhoQzR.setZ( cosPhoQzR );

//--- refracted 'photon' impact on quartz win. (inside) :
//    ---------------------------------------------------
      norm = - ddQzW / vDcPhoQzR.z();
      Hep3Vector vPhoWinQzR = vPhoWinOut + norm * vDcPhoQzR;

//--- corrections to x and y :
//    ------------------------
      Hep3Vector vcorr = vPhoWinQzR - vPhoWinIn;
      vcorr.setZ( 0.);
//pp      cout << "vcorr = " << vcorr << endl;

      Hep3Vector vPoCluQzw = vPoCluDet - vcorr;

      aflag = true;
      Hep3Vector vDcPhoCorr = getCluAngle( vPoPhotW, vPoCluQzw, RR, aflag );
//                            ---------------------------------------------
      if( !aflag ) return  vDcPhoEmW;

//--- return corrected photon pos. on quartz outer wind. :
//    ----------------------------------------------------
      //^return  vPoCluQzw;
      vPoCluDet = vPoCluQzw;

//--- return corrected photon direction :
//    -----------------------------------
      return  vDcPhoCorr;

  }


//======================================================================
  Hep3Vector CsRCPartPhotons::doQzWCorr( CsRCCluster* clu,
//----------------------------------------------------------------------
					 const Hep3Vector vPoPhotW,
					 const Hep3Vector vDcPhoEmW,
					 Hep3Vector& vPoCluDet,
					 const float RR, const float zDetW,
					 bool& aflag ) {


//--- approximate correction for photon detector quartz window :
//    ----------------------------------------------------------
//    new version, upgraded to PMT standard 5/07 

      CsRCExeKeys *key = CsRCExeKeys::Instance();
      CsRCRecConst *cons = CsRCRecConst::Instance();

      static bool firstCall = true;
      if( firstCall ) {
        firstCall = false;

        key->acknoMethod( "CsRCPartPhotons::QuartzWindowCorr" );
      }

      float CFRefInd = cons->CFRefInd();
      if( clu->isPMT() ) CFRefInd = cons->CFRefIndVS();
      if( clu->isAPV() ) CFRefInd = cons->CFRefIndUV();
      static float qzRefInd = cons->qzRefInd();

      CsRCDetectors *dets = CsRCDetectors::Instance();
      CsRCCathode *cat = dets->ptrToCat( clu->ic() );
      static float ddQzW = cat->ddQzW();
      static float ddGap = cat->ddGap();

      static double zQzWOutside = zDetW + ddGap + ddQzW;
      static double zQzWInside  = zDetW + ddGap;

      static double rRatioFQ = CFRefInd / qzRefInd;
      static double CHRefInd = 1.000444;
      static double rRatioQM = qzRefInd / CHRefInd;

//--- assume photon emitted from average point on part. traj.
//    -------------------------------------------------------
//--- 'photon' impact on mirror (MWR) :
//    ---------------------------------
      CsRCMirrors *mirr = CsRCMirrors::Instance();
      Hep3Vector vPoC( 0., 0., 0. );
      Hep3Vector vPoPhoMir = mirr->vImpMir( vPoPhotW, vDcPhoEmW, vPoC, RR );
//                           ----------------------------------------------
//--- normal to mirror at 'photon' impact :
//    -------------------------------------
      Hep3Vector vDcNoPhoMir = vPoPhoMir.unit();

//--- 'photon' reflected direction :
//    ------------------------------
      double cosPhoMir = vDcNoPhoMir * vDcPhoEmW;
      Hep3Vector vDcPhoRefl = 2.*cosPhoMir * vDcNoPhoMir - vDcPhoEmW;

//--- 'photon' incidence angle on detector :
//    --------------------------------------
      double cosPhoDet = vDcPhoRefl.z();
      double sinPhoDet = sqrt(1.-cosPhoDet*cosPhoDet);

      double norm;
//--- 'photon' impact on quartz wind. (outside=C4F10 side) :
//    ------------------------------------------------------
      norm = (zQzWOutside - vPoPhoMir.z()) / vDcPhoRefl.z();
      Hep3Vector vPhoWinOut = vPoPhoMir + norm * vDcPhoRefl;
      //std::cout << setprecision(5) << std::endl;
      //std::cout << " pro-qou " << vPhoWinOut << "  " << sinPhoDet << std::endl;

//--- 'photon' projection on quartz wind. (inside=CH4 side) :
//    -------------------------------------------------------
      norm = (zQzWInside - vPoPhoMir.z()) / vDcPhoRefl.z();
      Hep3Vector vPhoWinIn = vPoPhoMir + norm * vDcPhoRefl;
      //std::cout << " pro-qin " << vPhoWinIn << "  " << std::endl;

//--- 'photon' direction inside quartz wind. (DRS) :
//    ----------------------------------------------
      double sinPhoQzR = sinPhoDet * rRatioFQ;
      double cosPhoQzR = sqrt(1.-sinPhoQzR*sinPhoQzR);
      Hep3Vector vDcPhoQzR;
      vDcPhoQzR.setX( rRatioFQ * vDcPhoRefl.x() );
      vDcPhoQzR.setY( rRatioFQ * vDcPhoRefl.y() );
      vDcPhoQzR.setZ( cosPhoQzR );
//--- refracted 'photon' impact on quartz wind. (inside) :
//    ----------------------------------------------------
      norm = - ddQzW / vDcPhoQzR.z();
      Hep3Vector vPhoWinQzR = vPhoWinOut + norm * vDcPhoQzR;
      //std::cout << " rfr-qin " << vPhoWinQzR << "  " << sinPhoQzR << std::endl;

//--- 'photon' projection on CsI plane :
//    ----------------------------------
      norm = (zDetW - vPoPhoMir.z()) / vDcPhoRefl.z();
      Hep3Vector vPhoCsIP = vPoPhoMir + norm * vDcPhoRefl;
      //std::cout << " pro-csi " << vPhoCsIP << "  " << std::endl;
      //std::cout << " clu-csi " << vPoCluDet << "  " << std::endl;

//--- 'photon' direction outside quartz wind. (CH4 side) :
//    ----------------------------------------------------
      double sinPhoGap = sinPhoQzR * rRatioQM;
      double cosPhoGap = sqrt(1.-sinPhoGap*sinPhoGap);
      Hep3Vector vDcPhoGap;
      vDcPhoGap.setX( rRatioQM * vDcPhoQzR.x() );
      vDcPhoGap.setY( rRatioQM * vDcPhoQzR.y() );
      vDcPhoGap.setZ( cosPhoGap );
//--- refracted 'photon' impact on CsI plane :
//    ----------------------------------------
      norm = - ddGap / vDcPhoGap.z();
      Hep3Vector vPhoCsI = vPhoWinQzR + norm * vDcPhoGap;
      //std::cout << " rfr-csi " << vPhoCsI << "  " << sinPhoGap << std::endl;

      Hep3Vector vPoCluQzw;
      vPoCluQzw.setX( vPhoWinOut.x() + vPoCluDet.x() - vPhoCsI.x() );
      vPoCluQzw.setY( vPhoWinOut.y() + vPoCluDet.y() - vPhoCsI.y() );
      vPoCluQzw.setZ( zQzWOutside );
      //std::cout << " crr-qou " << vPoCluQzw << std::endl;

      aflag = true;
      Hep3Vector vDcPhoCorr = getCluAngle( vPoPhotW, vPoCluQzw, RR, aflag );
//                            ---------------------------------------------
      if( !aflag ) return  vDcPhoEmW;

//--- return corrected photon pos. on quartz outer wind. :
//    ----------------------------------------------------
      //^return  vPoCluQzw;
      vPoCluDet = vPoCluQzw;

//--- return corrected photon direction :
//    -----------------------------------
      return  vDcPhoCorr;

  }


//===========================================================================
  bool CsRCPartPhotons::GetPMTPhotonPosition( int kQua, int PMTchan,
//------------------------------------------------------------------
					      double xTgIn, double yTgIn,
					      double& xPosC, double& yPosC ) {

//- Paolo - June  2008

    static double parCorrX[16][6];
    static double parCorrY[16][6];

    CsRCExeKeys* key = CsRCExeKeys::Instance();
    CsRCRecConst *cons = CsRCRecConst::Instance();
    int kOffHi = cons->kOffHi();
    static std::vector<CsHist2D*> vRC3310;

    static int nTot = 0;
    static int nTgx = 0;
    static int nTgy = 0;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      key->acknoMethod( "CsRCPartPhotons::GetPMTPhotonPosition" );

      for( int kh=0; kh<4; kh++ ) vRC3310.push_back( NULL );
      if( CsRCHistos::Ref().bookHis() &&
	  CsRCHistos::Ref().levelBk() >= 2 ) {
	CsHistograms::SetCurrentPath("/RICH");
	vRC3310.clear();
        string hTitle = "tgPhoy vs x";
        int kHist = 0;
	stringstream hN3311;
        kHist = kOffHi + 3311;
        hN3311 << kHist;
        vRC3310.push_back( new CsHist2D( hN3311.str(), hTitle,
					 100, -0.5, 0.5, 100, -1., 1. ) );
	stringstream hN3312;
        kHist = kOffHi + 3312;
        hN3312 << kHist;
        vRC3310.push_back( new CsHist2D( hN3312.str(), hTitle,
					 100, -0.5, 0.5, 100, -1., 1. ) );
	stringstream hN3313;
        kHist = kOffHi + 3313;
        hN3313 << kHist;
        vRC3310.push_back( new CsHist2D( hN3313.str(), hTitle,
					 100, -0.5, 0.5, 100, -1., 1. ) );
	stringstream hN3314;
        kHist = kOffHi + 3314;
        hN3314 << kHist;
        vRC3310.push_back( new CsHist2D( hN3314.str(), hTitle,
					 100, -0.5, 0.5, 100, -1., 1. ) );
        CsHistograms::SetCurrentPath("/");
      }
    }

    double xTgA = xTgIn;
    double yTgA = yTgIn;

    double xh, yh, wh;
    int kHH = kQua;
    xh = xTgA;
    yh = yTgA;
    if( vRC3310[kHH] ) vRC3310[kHH]->Fill( xh, yh );
//hh                   ----------------------------

//- TopJura ONLY!

//- Fit   7 x 7
/*
    parCorrX[0][0] = -0.0151314;
    parCorrX[0][1] = 19.507;
    parCorrX[0][2] = -121.888;
    parCorrX[0][3] = 1.90156;
    parCorrX[0][4] = -7.0509;
    parCorrX[0][5] = 37.3838;

    parCorrY[0][0] = 2.35244;
    parCorrY[0][1] = -9.82821;
    parCorrY[0][2] = 39.5807;
    parCorrY[0][3] = 0.0106238;
    parCorrY[0][4] = 18.8585;
    parCorrY[0][5] = -117.598;

    parCorrX[1][0] = 2.74944;
    parCorrX[1][1] = -8.25521;
    parCorrX[1][2] = -46.9514;
    parCorrX[1][3] = 1.80271;
    parCorrX[1][4] = -4.17334;
    parCorrX[1][5] = 30.6329;

    parCorrY[1][0] = 1.99756;
    parCorrY[1][1] = -2.82464;
    parCorrY[1][2] = 16.3711;
    parCorrY[1][3] = 0.2456;
    parCorrY[1][4] = 19.7401;
    parCorrY[1][5] = -121.95;

    parCorrX[2][0] = 4.91152;
    parCorrX[2][1] = -34.0045;
    parCorrX[2][2] = 27.198;
    parCorrX[2][3] = 2.1765;
    parCorrX[2][4] = -7.61338;
    parCorrX[2][5] = 38.4747;

    parCorrY[2][0] = 3.02788;
    parCorrY[2][1] = -10.7496;
    parCorrY[2][2] = 22.8511;
    parCorrY[2][3] = 0.415853;
    parCorrY[2][4] = 20.9447;
    parCorrY[2][5] = -130.981;

    parCorrX[3][0] = 6.08624;
    parCorrX[3][1] = -51.573;
    parCorrX[3][2] = 96.2699;
    parCorrX[3][3] = 3.10862;
    parCorrX[3][4] = -18.4759;
    parCorrX[3][5] = 68.6659;

    parCorrY[3][0] = 3.18516;
    parCorrY[3][1] = -11.2442;
    parCorrY[3][2] = 17.7578;
    parCorrY[3][3] = 1.08019;
    parCorrY[3][4] = 14.7786;
    parCorrY[3][5] = -116.843;

    parCorrX[4][0] = -0.460823;
    parCorrX[4][1] = 25.146;
    parCorrX[4][2] = -135.375;
    parCorrX[4][3] = 3.1947;
    parCorrX[4][4] = -18.6698;
    parCorrX[4][5] = 57.6507;

    parCorrY[4][0] = 1.26225;
    parCorrY[4][1] = 0.221739;
    parCorrY[4][2] = 19.9642;
    parCorrY[4][3] = 2.95381;
    parCorrY[4][4] = -11.7015;
    parCorrY[4][5] = -40.0735;

    parCorrX[5][0] = 2.34297;
    parCorrX[5][1] = -4.48576;
    parCorrX[5][2] = -53.4881;
    parCorrX[5][3] = 3.01822;
    parCorrX[5][4] = -15.0128;
    parCorrX[5][5] = 45.9208;

    parCorrY[5][0] = 2.51276;
    parCorrY[5][1] = -10.116;
    parCorrY[5][2] = 36.2157;
    parCorrY[5][3] = 2.42412;
    parCorrY[5][4] = -5.25947;
    parCorrY[5][5] = -53.1776;

    parCorrX[6][0] = 4.32919;
    parCorrX[6][1] = -29.9489;
    parCorrX[6][2] = 27.2356;
    parCorrX[6][3] = 2.8645;
    parCorrX[6][4] = -12.3364;
    parCorrX[6][5] = 43.5658;

    parCorrY[6][0] = 8.87311;
    parCorrY[6][1] = -37.1742;
    parCorrY[6][2] = 82.9674;
    parCorrY[6][3] = 1.06625;
    parCorrY[6][4] = -3.26463;
    parCorrY[6][5] = -16.8219;

    parCorrX[7][0] = 4.88217;
    parCorrX[7][1] = -40.4011;
    parCorrX[7][2] = 74.1556;
    parCorrX[7][3] = 3.48731;
    parCorrX[7][4] = -14.5178;
    parCorrX[7][5] = 52.8625;

    parCorrY[7][0] = 3.31578;
    parCorrY[7][1] = -12.2616;
    parCorrY[7][2] = 18.6933;
    parCorrY[7][3] = 4.20204;
    parCorrY[7][4] = -19.1828;
    parCorrY[7][5] = -29.9072;

    parCorrX[8][0] = -0.247924;
    parCorrX[8][1] = 22.1887;
    parCorrX[8][2] = -127.629;
    parCorrX[8][3] = 3.1615;
    parCorrX[8][4] = -16.565;
    parCorrX[8][5] = 43.9261;

    parCorrY[8][0] = 1.89973;
    parCorrY[8][1] = -9.13054;
    parCorrY[8][2] = 42.285;
    parCorrY[8][3] = 4.48164;
    parCorrY[8][4] = -33.4172;
    parCorrY[8][5] = 26.1717;

    parCorrX[9][0] = 2.54054;
    parCorrX[9][1] = -5.41096;
    parCorrX[9][2] = -54.8382;
    parCorrX[9][3] = 2.85997;
    parCorrX[9][4] = -13.0768;
    parCorrX[9][5] = 31.2261;

    parCorrY[9][0] = 2.29114;
    parCorrY[9][1] = -12.7876;
    parCorrY[9][2] = 43.6088;
    parCorrY[9][3] = 4.63304;
    parCorrY[9][4] = -33.3815;
    parCorrY[9][5] = 24.6482;

    parCorrX[10][0] = 4.26206;
    parCorrX[10][1] = -26.414;
    parCorrX[10][2] = 15.1397;
    parCorrX[10][3] = 3.22276;
    parCorrX[10][4] = -14.1999;
    parCorrX[10][5] = 34.9172;

    parCorrY[10][0] = 3.26009;
    parCorrY[10][1] = -19.0385;
    parCorrY[10][2] = 50.0844;
    parCorrY[10][3] = 4.42616;
    parCorrY[10][4] = -30.6398;
    parCorrY[10][5] = 19.5769;

    parCorrX[11][0] = 3.71774;
    parCorrX[11][1] = -30.2784;
    parCorrX[11][2] = 56.3024;
    parCorrX[11][3] = 5.34971;
    parCorrX[11][4] = -25.1384;
    parCorrX[11][5] = 74.5304;

    parCorrY[11][0] = 0.00456999;
    parCorrY[11][1] = -0.0197194;
    parCorrY[11][2] = 0.0321399;
    parCorrY[11][3] = 3409.44;
    parCorrY[11][4] = -20578.5;
    parCorrY[11][5] = 3144.83;

    parCorrX[12][0] = 0.0127803;
    parCorrX[12][1] = 21.2016;
    parCorrX[12][2] = -131.435;
    parCorrX[12][3] = 3.06532;
    parCorrX[12][4] = -13.8497;
    parCorrX[12][5] = 24.7518;

    parCorrY[12][0] = 2.5566;
    parCorrY[12][1] = -24.5402;
    parCorrY[12][2] = 97.248;
    parCorrY[12][3] = 3.77414;
    parCorrY[12][4] = -31.4924;
    parCorrY[12][5] = 53.0026;

    parCorrX[13][0] = 2.93164;
    parCorrX[13][1] = -9.32024;
    parCorrX[13][2] = -43.3619;
    parCorrX[13][3] = 2.51425;
    parCorrX[13][4] = -4.03899;
    parCorrX[13][5] = -6.8395;

    parCorrY[13][0] = 2.53525;
    parCorrY[13][1] = -24.2726;
    parCorrY[13][2] = 83.1439;
    parCorrY[13][3] = 4.64471;
    parCorrY[13][4] = -38.0315;
    parCorrY[13][5] = 64.063;

    parCorrX[14][0] = 3.97919;
    parCorrX[14][1] = -22.4648;
    parCorrX[14][2] = 7.76546;
    parCorrX[14][3] = 3.20061;
    parCorrX[14][4] = -7.16681;
    parCorrX[14][5] = 0.610435;

    parCorrY[14][0] = 3.48531;
    parCorrY[14][1] = -27.7353;
    parCorrY[14][2] = 81.4562;
    parCorrY[14][3] = 4.39783;
    parCorrY[14][4] = -36.0281;
    parCorrY[14][5] = 63.7598;

    parCorrX[15][0] = 3.3969;
    parCorrX[15][1] = -25.4168;
    parCorrX[15][2] = 44.6132;
    parCorrX[15][3] = 4.86843;
    parCorrX[15][4] = -8.21762;
    parCorrX[15][5] = 3.33895;

    parCorrY[15][0] = 3.36804;
    parCorrY[15][1] = -24.0525;
    parCorrY[15][2] = 64.4272;
    parCorrY[15][3] = 5.07836;
    parCorrY[15][4] = -38.8856;
    parCorrY[15][5] = 62.7314;
*/


//- Fit   5 x 5
//- 080805
    parCorrX[0][0] = 0.754319;
    parCorrX[0][1] = 2.57719;
    parCorrX[0][2] = -95.8263;
    parCorrX[0][3] = 1.90035;
    parCorrX[0][4] = -7.79282;
    parCorrX[0][5] = 42.8449;

    parCorrY[0][0] = 1.71826;
    parCorrY[0][1] = -6.19579;
    parCorrY[0][2] = 38.7801;
    parCorrY[0][3] = 0.922878;
    parCorrY[0][4] = 2.26862;
    parCorrY[0][5] = -95.5399;

    parCorrX[1][0] = 2.4687;
    parCorrX[1][1] = -19.527;
    parCorrX[1][2] = -19.6916;
    parCorrX[1][3] = 1.95275;
    parCorrX[1][4] = -9.18531;
    parCorrX[1][5] = 49.9597;

    parCorrY[1][0] = 2.16099;
    parCorrY[1][1] = -13.5417;
    parCorrY[1][2] = 59.0578;
    parCorrY[1][3] = 1.09876;
    parCorrY[1][4] = 2.83054;
    parCorrY[1][5] = -98.7784;

    parCorrX[2][0] = 2.93441;
    parCorrX[2][1] = -24.7734;
    parCorrX[2][2] = -2.20294;
    parCorrX[2][3] = 1.71859;
    parCorrX[2][4] = -5.03808;
    parCorrX[2][5] = 34.7057;

    parCorrY[2][0] = 2.43687;
    parCorrY[2][1] = -16.6557;
    parCorrY[2][2] = 62.709;
    parCorrY[2][3] = 1.2496;
    parCorrY[2][4] = 4.19256;
    parCorrY[2][5] = -107.828;

    parCorrX[3][0] = 2.6491;
    parCorrX[3][1] = -28.2081;
    parCorrX[3][2] = 58.4731;
    parCorrX[3][3] = 2.4993;
    parCorrX[3][4] = -6.28056;
    parCorrX[3][5] = 44.7864;

    parCorrY[3][0] = 2.35319;
    parCorrY[3][1] = -9.21603;
    parCorrY[3][2] = 21.2339;
    parCorrY[3][3] = 1.63544;
    parCorrY[3][4] = -0.092624;
    parCorrY[3][5] = -94.6966;

    parCorrX[4][0] = 0.64539;
    parCorrX[4][1] = 3.57493;
    parCorrX[4][2] = -95.8803;
    parCorrX[4][3] = 2.24562;
    parCorrX[4][4] = -14.7217;
    parCorrX[4][5] = 62.9569;

    parCorrY[4][0] = 1.69446;
    parCorrY[4][1] = -5.64783;
    parCorrY[4][2] = 38.6924;
    parCorrY[4][3] = 2.30762;
    parCorrY[4][4] = -19.9848;
    parCorrY[4][5] = -16.4366;

    parCorrX[5][0] = 2.1385;
    parCorrX[5][1] = -14.9952;
    parCorrX[5][2] = -29.2448;
    parCorrX[5][3] = 2.10256;
    parCorrX[5][4] = -13.8358;
    parCorrX[5][5] = 58.8106;

    parCorrY[5][0] = 2.03192;
    parCorrY[5][1] = -12.18;
    parCorrY[5][2] = 57.0572;
    parCorrY[5][3] = 2.24421;
    parCorrY[5][4] = -17.4766;
    parCorrY[5][5] = -20.612;

    parCorrX[6][0] = 2.82123;
    parCorrX[6][1] = -23.2772;
    parCorrX[6][2] = 4.89029;
    parCorrX[6][3] = 1.94017;
    parCorrX[6][4] = -10.5007;
    parCorrX[6][5] = 48.4605;

    parCorrY[6][0] = 2.34229;
    parCorrY[6][1] = -13.5065;
    parCorrY[6][2] = 50.3732;
    parCorrY[6][3] = 2.3726;
    parCorrY[6][4] = -17.0289;
    parCorrY[6][5] = -22.7841;

    parCorrX[7][0] = 2.87798;
    parCorrX[7][1] = -28.5379;
    parCorrX[7][2] = 56.6256;
    parCorrX[7][3] = 2.58323;
    parCorrX[7][4] = -11.1772;
    parCorrX[7][5] = 54.5515;

    parCorrY[7][0] = 0.251752;
    parCorrY[7][1] = -0.853093;
    parCorrY[7][2] = 2.05168;
    parCorrY[7][3] = 26.4132;
    parCorrY[7][4] = -190.137;
    parCorrY[7][5] = -178.983;

    parCorrX[8][0] = 0.627406;
    parCorrX[8][1] = 1.93647;
    parCorrX[8][2] = -83.5079;
    parCorrX[8][3] = 2.2564;
    parCorrX[8][4] = -12.8881;
    parCorrX[8][5] = 47.6527;

    parCorrY[8][0] = 1.36406;
    parCorrY[8][1] = 1.33206;
    parCorrY[8][2] = 12.9431;
    parCorrY[8][3] = 2.07908;
    parCorrY[8][4] = -23.3674;
    parCorrY[8][5] = 10.3376;

    parCorrX[9][0] = 2.01306;
    parCorrX[9][1] = -12.8855;
    parCorrX[9][2] = -30.9624;
    parCorrX[9][3] = 2.04699;
    parCorrX[9][4] = -11.2424;
    parCorrX[9][5] = 39.7165;

    parCorrY[9][0] = 1.72596;
    parCorrY[9][1] = -4.07223;
    parCorrY[9][2] = 25.2703;
    parCorrY[9][3] = 1.96609;
    parCorrY[9][4] = -20.9369;
    parCorrY[9][5] = 9.59968;

    parCorrX[10][0] = 2.55565;
    parCorrX[10][1] = -17.6448;
    parCorrX[10][2] = -6.56169;
    parCorrX[10][3] = 2.01472;
    parCorrX[10][4] = -8.21276;
    parCorrX[10][5] = 29.7954;

    parCorrY[10][0] = 2.00288;
    parCorrY[10][1] = -6.6613;
    parCorrY[10][2] = 25.8683;
    parCorrY[10][3] = 2.16763;
    parCorrY[10][4] = -20.9568;
    parCorrY[10][5] = 6.6894;

    parCorrX[11][0] = 3.04396;
    parCorrX[11][1] = -27.4402;
    parCorrX[11][2] = 51.7735;
    parCorrX[11][3] = 2.39997;
    parCorrX[11][4] = -6.48283;
    parCorrX[11][5] = 26.1237;

    parCorrY[11][0] = 2.12523;
    parCorrY[11][1] = -3.74309;
    parCorrY[11][2] = 2.41474;
    parCorrY[11][3] = 2.67892;
    parCorrY[11][4] = -24.2405;
    parCorrY[11][5] = 6.16572;

    parCorrX[12][0] = 0.507961;
    parCorrX[12][1] = 3.35701;
    parCorrX[12][2] = -87.278;
    parCorrX[12][3] = 2.20817;
    parCorrX[12][4] = -8.20246;
    parCorrX[12][5] = 19.4872;

    parCorrY[12][0] = 1.22922;
    parCorrY[12][1] = 4.69464;
    parCorrY[12][2] = 3.98457;
    parCorrY[12][3] = 2.31371;
    parCorrY[12][4] = -29.4359;
    parCorrY[12][5] = 65.097;

    parCorrX[13][0] = 2.18544;
    parCorrX[13][1] = -14.6463;
    parCorrX[13][2] = -22.2393;
    parCorrX[13][3] = 1.98991;
    parCorrX[13][4] = -4.01999;
    parCorrX[13][5] = 1.70169;

    parCorrY[13][0] = 1.77142;
    parCorrY[13][1] = -3.31407;
    parCorrY[13][2] = 26.7721;
    parCorrY[13][3] = 2.20029;
    parCorrY[13][4] = -26.5827;
    parCorrY[13][5] = 59.8189;

    parCorrX[14][0] = 2.37619;
    parCorrX[14][1] = -15.4882;
    parCorrX[14][2] = -4.51194;
    parCorrX[14][3] = 2.28717;
    parCorrX[14][4] = -2.06541;
    parCorrX[14][5] = -7.24181;

    parCorrY[14][0] = 2.29913;
    parCorrY[14][1] = -5.98765;
    parCorrY[14][2] = 28.3199;
    parCorrY[14][3] = 2.15243;
    parCorrY[14][4] = -24.5195;
    parCorrY[14][5] = 55.3347;

    parCorrX[15][0] = 2.9255;
    parCorrX[15][1] = -24.2422;
    parCorrX[15][2] = 44.2135;
    parCorrX[15][3] = 2.63699;
    parCorrX[15][4] = -2.31273;
    parCorrX[15][5] = -1.02296;

    parCorrY[15][0] = 2.48373;
    parCorrY[15][1] = 0.135761;
    parCorrY[15][2] = -11.4553;
    parCorrY[15][3] = 2.37885;
    parCorrY[15][4] = -25.7152;
    parCorrY[15][5] = 56.1766;
//- 080805


    static double xDirec[7], yDirec[7];
    static double xDirLm, xDirLx, yDirLm, yDirLx;
    static double dDirx, dDiry;
    xDirec[0] = -0.24;
    xDirec[1] = -0.19;
    xDirec[2] = -0.14;
    xDirec[3] = -0.09;
    xDirec[4] = -0.04;
    xDirec[5] =  0.01;
    xDirec[6] =  0.06;
    dDirx = 0.05;
    //xDirLm = -xDirec[6] - dDirx/2.;
    //xDirLx = -xDirec[0] + dDirx/2.;
    xDirLm = -xDirec[6] - dDirx;
    xDirLx = -xDirec[0] + dDirx;
//@@---------------------------
    yDirec[0] =  0.1;
    yDirec[1] =  0.15;
    yDirec[2] =  0.2;
    yDirec[3] =  0.25;
    yDirec[4] =  0.3;
    yDirec[5] =  0.35;
    yDirec[6] =  0.4;
    dDiry = 0.05;
    //yDirLm = -yDirec[6] - dDiry/2.;
    //yDirLx = -yDirec[0] + dDiry/2.;
    yDirLm = -yDirec[6] - dDiry;
    yDirLx = -yDirec[0] + dDiry;
//@@---------------------------

//- Angular limits for correct function evaluation
    static double xDirLmw, xDirLxw, yDirLmw, yDirLxw;
    xDirLmw = -xDirec[6];
    xDirLxw = -xDirec[0];
    yDirLmw = -yDirec[6];
    yDirLxw = -yDirec[0];
//@@-------------------------------
    //std::cout << xDirLmw << "  " << xDirLxw << "  " << yDirLmw << "  "
    //          << yDirLxw << std::endl;;

    xPosC = 0.;
    yPosC = 0.;

//- check and count angular range outliers
    nTot++;
    if( xTgA < xDirLm  ||  xTgA > xDirLx ) {
      nTgx++;
      //std::cout << "GetPMTPhotonPosition : tgX out of lims  " << xTgA
      //          << "  " << xDirLm << "  " << xDirLx << std::endl;
      //return  false;
    }
    if( yTgA < yDirLm  ||  yTgA > yDirLx ) {
      nTgy++;
      //std::cout << "GetPMTPhotonPosition : tgY out of lims  " << yTgA
      //          << "  " << yDirLm << "  " << yDirLx << std::endl;
      //return  false;
    }
    if( nTot%10000 == 0 ) {
      //std::cout << "GetPMTPhotonPosition : tot " << nTot << "  tgx out "
      // 	  << nTgx << "  tgy out " << nTgy << std::endl; 
    }

    if( xTgA < xDirLmw ) xTgA = xDirLmw;
    if( xTgA > xDirLxw ) xTgA = xDirLxw;
    if( yTgA < yDirLmw ) yTgA = yDirLmw;
    if( yTgA > yDirLxw ) yTgA = yDirLxw;
//@@-----------------------------------

//- 'PMTchan' (as in/from TestPMT.cc, the corr.s comp. code) is 'wrong' or
//  inverted, i.e. = 16-kPPad-1 instead of = kPPad; this explain the inverted
//  order in 'chaw' comp. from iPadx iPady in doPMTOptCorr

//- xTgAoff, yTgAoff come from the Fit procedure in TestAMC.cc
    double xTgAoff = -0.035;
    double yTgAoff = -0.375;
    xTgA -= xTgAoff;  
    yTgA -= yTgAoff;
//@@---------------

    int kch = PMTchan;
    double xPoX = parCorrX[kch][0] + parCorrX[kch][1]*xTgA +
                  parCorrX[kch][2]*xTgA*xTgA;
    double yPoX = parCorrX[kch][3] + parCorrX[kch][4]*yTgA +
                  parCorrX[kch][5]*yTgA*yTgA;
    xPosC = xPoX * yPoX;

    double xPoY = parCorrY[kch][0] + parCorrY[kch][1]*xTgA +
                  parCorrY[kch][2]*xTgA*xTgA;
    double yPoY = parCorrY[kch][3] + parCorrY[kch][4]*yTgA +
                  parCorrY[kch][5]*yTgA*yTgA;
    yPosC = xPoY * yPoY;

    //std::cout << setprecision(6);
    //std::cout << "GetPMTPhotonPosition  " << xTgA << "  " << yTgA << "  "
    //          << xPosC << "  " << yPosC << "  ";
    //std::cout << std::endl;

//- See comment on PMTchan above
    int kPad = 16 - kch - 1;
    double xPC = 0.;
    double yPC = 0.;
    GetPadCentre( kPad, xPC, yPC );
//--------------------------------
    xPosC += xPC;
    yPosC += yPC;
    //std::cout << kPad << "  " << xPC << "  " << yPC << " ";
    //std::cout << xPosC << "  " << yPosC;
    //std::cout << std::endl;

    return  true;
  }


//===========================================================================
  bool CsRCPartPhotons::GetPadCentre( int ppad, double& xc, double& yc ) {
//------------------------------------------------------------------------

//  Paolo   April 2008
//  Provisional!

//- TopJura ONLY!

    int k, j;
    double dx, dy1, dy2;
    double xx1, xx2, yy11, yy12, yy21, yy22;

    k = ppad%4;
    j = int( ppad/4 );
    k -= 2;
    j -= 2;
    dx  = 11.97;
    dy1 =  0.29;
    dy2 = 11.62;
    xx1 = 0. + k * dx;
    xx2 = xx1 + dx;
    yy11 = 0. + j * dy2 + k * dy1;
    yy12 = yy11 + dy2;
    yy21 = yy11 - dy1;
    yy22 = yy12 - dy1;
    //std::cout << xx1 << "  " << yy11 << "  "
    //      << xx1 << "  " << yy12 << "  "
    //      << xx2 << "  " << yy21 << "  "
    //      << xx2 << "  " << yy22 << std::endl;
    xc = (xx1 + xx2) / 2.;
    yc = (yy11 + yy12 + yy21 + yy22) / 4.;
    //std::cout << xc << "  " << yc << std::endl;

    return  true;
  }


//===========================================================================
  double CsRCPartPhotons::getThetaIpo( double refrInd, double momPart, 
//--------------------------------------------------------------------
				       double massPart ) {

//  Paolo   October 2008

    double thetaIpo = 0.;
    double betaPart = momPart / sqrt( momPart*momPart + massPart*massPart );
    double cosTheta = 1000.;
    if( betaPart > 0. && refrInd > 0. ) cosTheta = 1./( betaPart * refrInd );
    if( cosTheta <= 1. ) thetaIpo = acos( cosTheta ) * 1000.;

    return thetaIpo;
  }



// WARNING : TOO MANY LINES (>8000) -> NO emacs STYLE !!!

