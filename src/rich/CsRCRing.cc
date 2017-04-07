/*!
   \file    CsRCRing.cc
   \-------------------
   \brief   CsRCRing class implementation.
   \author  Paolo Schiavon
   \version 0.01
   \date    October 2000
*/

//- revised  August 2005


  #include <iostream>
  #include <ostream>
  #include <cstdio>
  #include <cmath>
  #include <cstdlib>

  #include <CLHEP/Vector/ThreeVector.h>
  #include "CLHEP/Matrix/Matrix.h"
  #include "CLHEP/Random/RandFlat.h"
  #include "CLHEP/Random/RandGauss.h"
//-------------------------------------

  #include "CsErrLog.h"
  #include "CsTrack.h"
  #include "CsGeom.h"

// ----------------------------
  #include "CsRCRing.h"

  #include "CsRichOne.h"

  #include "CsRCDetectors.h"
  #include "CsRCPhotonDet.h"
  #include "CsRCMirrors.h"

  #include "CsRCParticle.h"
  #include "CsRCPad.h"
  #include "CsRCCluster.h"

  #include "CsRCEventPartPhotons.h"
  #include "CsRCPartPhotons.h"
  #include "CsRCPhoton.h"
  #include "CsRCLikeRing02.h"
  #include "CsRCLikeRing03.h"
  #include "CsRCLikeRing04.h"
  #include "CsRCLikeRing05.h"

  #include "CsRCEventAnalysis.h"

  #include "CsRCChiSqFit.h"
  #include "CsRCCircleFit.h"
  #include "CsRCEllipseFit.h"
  #include "CsRCGauPolFit.h"

  #include "CsRCExeKeys.h"
  #include "CsRCRecConst.h"
  #include "CsRCHistos.h"
// ----------------------------

  using namespace std;
  using namespace CLHEP;

//=========================================================================== 
  CsRCRing::CsRCRing() {
//----------------------
    kRing_ = -1;
    pPartPhot_ = NULL;
    the_ = 0.;
    lPhotons_.clear();
    kRingLoc_ = 0;
    theLoL_ = 0.;
    theUpL_ = 0.;
    theReco_ = 0.;
    thetaWgav_ = 0.;
    thetaRFit_ = 0.;
    thetaLike_ = 0.;
    thetaLikeSet_ = false;
    mTime_ = 0.;
    mT0_ = 0.;
    bool mTimeSet_ = false;
    nPhotPMT_ = 0;
    CsRCRecConst *cons = CsRCRecConst::Instance();
    int  mcaScan = cons->mcaScan();
    for( int j=0; j < mcaScan; j++ ) binCont_[j] = 0.;
    //for( int j=0; j<15; j++ ) partProbs_[j] = 0.;
    int nProb = cons->outBufferSize();
    for( int j=0; j<nProb; j++ ) partProbs_[j] = 0.;
    partProbsSet_ = false;
    for( int j=0; j<31; j++ ) probaLKRing_[j] = 0.;
    for( int j=0; j<31; j++ ) derivLKRing_[j] = 0.;
    for( int j=0; j<31; j++ ) qSquareRing_[j] = 0.;
    probaLKBgRing_ = 0.;
    ringQsQ_ = 0.;
    nPhotQsQ_ = 0;
    flagReco_ = false;
    flagOverThrs_ = false;
    flag_ = false;
    flagBack_ = false;
    ringBack_ = 0.;
  }

//=========================================================================== 
  CsRCRing::CsRCRing( const CsRCRing &ring ) {
//--------------------------------------------
    //cout << "RICHONE : CsRCRing CopyConstructor" << endl;
    kRing_ = ring.kRing_;
    pPartPhot_ = ring.pPartPhot_;
    the_ = ring.the_;
    lPhotons_ = ring.lPhotons_;
    kRingLoc_ = ring.kRingLoc_;
    theLoL_ = ring.theLoL_;
    theUpL_ = ring.theUpL_;
    theReco_ = ring.theReco_;
    thetaWgav_ = ring.thetaWgav_;
    thetaRFit_ = ring.thetaRFit_;
    thetaLike_ = ring.thetaLike_;
    thetaLikeSet_ = ring.thetaLikeSet_;
    mTime_ = ring.mTime_;
    mT0_ = ring.mT0_;
    mTimeSet_ = ring.mTimeSet_;
    nPhotPMT_ = ring.nPhotPMT_;
    CsRCRecConst *cons = CsRCRecConst::Instance();
    int  mcaScan = cons->mcaScan();
    for( int j=0; j < mcaScan; j++ ) binCont_[j] = ring.binCont_[j];
    //for( int j=0; j<15; j++ ) partProbs_[j] = ring.partProbs_[j];
    int nProb = cons->outBufferSize();
    for( int j=0; j<nProb; j++ ) partProbs_[j] = ring.partProbs_[j];
    partProbsSet_ = ring.partProbsSet_;
    for( int j=0; j<31; j++ ) probaLKRing_[j] = ring.probaLKRing_[j];
    for( int j=0; j<31; j++ ) derivLKRing_[j] = ring.derivLKRing_[j];
    for( int j=0; j<31; j++ ) qSquareRing_[j] = ring.qSquareRing_[j];
    probaLKBgRing_ = ring.probaLKBgRing_;
    ringQsQ_ = ring.ringQsQ_;
    nPhotQsQ_ = ring.nPhotQsQ_;
    flagReco_ = ring.flagReco_;
    flagOverThrs_ = ring.flagOverThrs_;
    flag_ = ring.flag_;
    flagBack_ = ring.flagBack_;
    ringBack_ = ring.ringBack_;
  }

//=========================================================================== 
  void CsRCRing::print() const {
//------------------------------
    cout << endl;
    cout << " Ring  " << kRing()
	 << " , PartPhotons  " << pPartPhot()->kPaPhot()
         << " , Particle  " << pPartPhot()->pPart()->kPart() << endl;
    cout << "   theta = " << the() << ", photons  " << lPhotons_.size()
	 << " , n. :   ";
    list<CsRCPhoton*>::const_iterator ih;
    for( ih=lPhotons_.begin();  ih!=lPhotons_.end(); ih++ ) {
      cout << (*ih)->kPhot() << "  ";
    }
    cout << endl;
    if( CsRCExeKeys::Instance()->MCarloEvent() ) {
      cout << "   thetaMC = " << pPartPhot()->pPart()->thetaCer()
	   << ", photons  " << pPartPhot()->pPart()->nPhoCer() << endl;
    }
    CsRCParticle* part = pPartPhot()->pPart();
    cout << "   part. on PD (vPade) (MWR)  "
      //<< pPartPhot()->vCorPoPhoW()[part->kDetPart()] << "  "
         << part->vPade()[part->kDetPart()] << endl;
  }

//=========================================================================== 
  void CsRCRing::printPhotons() const {
//-------------------------------------
    list<CsRCPhoton*>::const_iterator ih;
    for( ih=lPhotons_.begin();  ih!=lPhotons_.end(); ih++ ) {
      (*ih)->print();
      if( CsRCExeKeys::Instance()->kPrintEventRings() == 4 ) {
	(*ih)->ptToClu()->print();
	list<CsRCPad*> lPads = (*ih)->ptToClu()->lPads();
	list<CsRCPad*>::iterator ia;
	for( ia=lPads.begin(); ia!=lPads.end(); ia++ ) (*ia)->print();
      }
    }
  }

//=========================================================================== 
  void CsRCRing::printPhoFlag() const {
//-------------------------------------
    cout << " Ring  " << kRing() << "  flags : ";
    list<CsRCPhoton*>::const_iterator ih;
    for( ih=lPhotons_.begin();  ih!=lPhotons_.end(); ih++ )
      cout << (*ih)->flag() << "  ";
    cout << endl;
  }

//=========================================================================== 
  CsRCRing::~CsRCRing() {
//-----------------------
    //for( it=lPhotons_.begin(); ... ) delete *it;   //   NO
    lPhotons_.clear();
  }


//=========================================================================== 
  void CsRCRing::getRingPk( const int kRing, CsRCPartPhotons* papho ) {
//---------------------------------------------------------------------

//--- Paolo - October 2000
//---         revised September 2002


//--- 'Peak' ring constructor.
//    ------------------------


//--- from "CsRCExecKeys"
      CsRCExeKeys *key = CsRCExeKeys::Instance();
      bool DoClu = key->DoClu();
      //bool DoCluCo = key->DoCluCo();
      bool UseCluPH = key->UseCluPH();
      int kDoClu = 0;
      if( DoClu ) kDoClu = 0;
      if( !DoClu ) kDoClu = 1;

//--- from "CsRCRecConst"
      CsRCRecConst *cons = CsRCRecConst::Instance();
      int *mcanWind = cons->mcanWind();
      int mcaScan = cons->mcaScan();
      float xlScan = cons->xlScan();
      float xuScan = cons->xuScan();
      float binScan = cons->binScan();
      int nMore = cons->nMore();                          //   040213
      //int nMore = getVarNMore( papho, cons->nMore() );      //   040213
      int nScanMx = mcaScan;

      static bool firstCall = true;
      if( firstCall ) {
        firstCall = false;
        key->acknoMethod( "CsRCRing::getRingPk" );
      }

      list<CsRCPhoton*> lPhotons = papho->lPhotons();
      list<CsRCPhoton*>::iterator ih;

      int mcan = mcanWind[kDoClu];                        //   040213
      //int mcan = getVarMcan( papho, mcanWind[kDoClu] );     //   040213
      int mcanmx = 0;
      int ncan1 = 0;
      int ncan2 = 0;
      int kcanmx = -1;

      int kCanMx = 0;
      string peakSrcMode = cons->peakSrcMode();
      if( peakSrcMode == "COUNT" ) {
	if( peakCountSearch( papho, kCanMx ) ) kcanmx = kCanMx;
//        ------------------------------------
      }
      else if( peakSrcMode == "CHI" ) {
	//^if( peakChiSearch( papho, kCanMx ) ) kcanmx = kCanMx;
	//^if( peakCountSearchMC( papho, kCanMx ) ) kcanmx = kCanMx;
	//^if( peakCountSearchVB( papho, kCanMx ) ) kcanmx = kCanMx;
	//^if( peakCountSearchTEST( papho, kCanMx ) ) kcanmx = kCanMx;
//        ----------------------------------------
	string mess = "RICHONE, getRingPk() : ";
	string err = "peak search mode requested not VALID";
	mess.append( err );
	CsErrLog::Instance()->mes( elError, mess );
      }
      else if( peakSrcMode == "MASS" ) {
	//^if( peakMassSearch( papho, kCanMx ) ) kcanmx = kCanMx;
//        -----------------------------------
	string mess = "RICHONE, getRingPk() : ";
	string err = "peak search mode requested not VALID";
	mess.append( err );
	CsErrLog::Instance()->mes( elError, mess );
      }
      else if( peakSrcMode == "LIKE" ) {
	//^if( peakLikeSearch( papho, kCanMx ) ) kcanmx = kCanMx;
//        -----------------------------------
	string mess = "RICHONE, getRingPk() : ";
	string err = "peak search mode requested not VALID";
	mess.append( err );
	CsErrLog::Instance()->mes( elError, mess );
      }
      else {
	cout << "RICHONE, getRingPk() : WRONG peak search mode requested" 
	     << endl;
	cout << "-------------------------------------------------------"
	     << endl;
	string mess = "RICHONE, getRingPk() : ";
	string err = "WRONG peak search mode requested";
	mess.append( err );
	CsErrLog::Instance()->mes( elError, mess );
      }
      if( kcanmx < 0 ) { 
	flag_ = false;
	//std::cout << kcanmx << std::endl;
	return;
      }

      float theMax = (float( kcanmx ) + float( mcan )/2. + 0.5)* binScan;
      //cout << theMax << endl;

//--- monitor pattern rec. windows :
//    ------------------------------
      CsRCHistos& hist = CsRCHistos::Ref();
      double xh, yh, wh;
      xh = float( mcan );
      yh = float( nMore );
      if( hist.hRC1708 ) hist.hRC1708->Fill( xh, yh );

      histCntRatio();
//    --------------

//--- compute average around peak (width = mcan +/- nMore) :
//    ------------------------------------------------------
      ncan1 = kcanmx - nMore;
      ncan2 = kcanmx + mcan - 1 + nMore;
      //ncan1 = kcanmx;              //   030730
      //ncan2 = kcanmx + mcan - 1;   //   030730
//@@---------------------------------------
      float theMN = xlScan + ncan1 * binScan;
      float theMX = xlScan + ncan2 * binScan + binScan;
      //cout << theMN << "  " << theMX << "  ";

      double thePhotPk = 0.;
      float totCntPk = 0.;
      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
        if( !(*ih)->flag() ) continue;
        double thew = (*ih)->the();
        if( thew >= theMN  &&  thew <= theMX) {
          float PHw = 1.;
//@@--------------------
          if( UseCluPH ) PHw = (*ih)->PH() + 1.;
//        -------------------------------------
	  if( (*ih)->isPMT() ) PHw = 1.;
//        -----------------------------
          thePhotPk += thew * PHw;
          totCntPk += PHw;
	}
      }
      if( totCntPk > 0. ) thePhotPk /= totCntPk;

//    sort ring clusters inside window (+/- nMore) around average :
//    -------------------------------------------------------------
      float hWind = float( mcan + 2*nMore )* binScan /2.;
      theMN = thePhotPk - hWind;
      theMX = thePhotPk + hWind;
      //cout << theMN << "  " << theMX << endl;
      nPhotPMT_ = 0;
      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
        if( !(*ih)->flag() ) continue;
	double thew = (*ih)->the();
        if( thew >= theMN  &&  thew <= theMX) {
	  lPhotons_.push_back( (*ih) );
//        ----------------------------
	  if( (*ih)->isPMT() ) nPhotPMT_++;
	}
      }

//--- 03/05/14
/*
//    WARNING: phiA wrong! others ambiguous...
      //std::cout << lPhotons_.size() << std::endl;
      lDoublePhotons_.clear();
      for( ih=lPhotons_.begin(); ih!=lPhotons_.end(); ih++ ) {
	float phiD = (*ih)->phi() + 180.;
	if( phiD > 360.) phiD -= 360.;
	//std::cout << (*ih)->phi() << "  " << phiD << std::endl;
	lDoublePhotons_.push_back( new CsRCPhoton( (*ih)->kPhot(),
	  (*ih)->the(), phiD, (*ih)->phiA(), (*ih)->theB(),
          (*ih)->theM(), (*ih)->PH(), (*ih)->ptToClu(), (*ih)->ptToDet(),
          (*ih)->pPartPhot(), (*ih)->kDetClu() ) );
      }
      for( ih=lDoublePhotons_.begin(); ih!=lDoublePhotons_.end(); ih++ ) {
	lPhotons_.push_back( (*ih) );
//      ----------------------------
      }
      //std::cout << lPhotons_.size() << std::endl;
*/
//--- 03/05/14

      if( !lPhotons_.empty() ) {
//----- ring location ( kRingLoc_: up = +1, split = 0, down = -1 ) :
//      ------------------------------------------------------------
        CsRCPhotonDet* detFirst = lPhotons_.front()->ptToDet();
        kRingLoc_ = detFirst->kDet();
        if( kRingLoc_ == 1 ) kRingLoc_ = -1;
        if( kRingLoc_ == 0 ) kRingLoc_ = +1;
        for( ih=lPhotons_.begin(); ih!=lPhotons_.end(); ih++ ) {
	  if( (*ih)->ptToDet() != detFirst ) { kRingLoc_ = 0;   break; }
        }

//----- fill Ring-Pk object :
//      ---------------------
        kRing_ = kRing;
        pPartPhot_ = papho;
        the_ = thePhotPk;
	partProbsSet_ = false;
        flagReco_ = false;
	flagOverThrs_ = false;
        flag_ = true;
	flagBack_ = false;
	ringBack_ = 1.;

        papho->setpRing( this );
        papho->pPart()->setpRing( this );

//----- monitoring photons in ring :
//      ----------------------------
        int nPhotPk = lPhotons_.size();
        xh = nPhotPk;
	if( hist.hRC1545 ) hist.hRC1545->Fill( xh );
//hh                       ------------------------

      } else { flag_ = false; }   // end if not empty
      //if( flag_ == false ) {
      //std::cout << "getRingPk 3  " << CsRichOne::Instance()->kEvent()
      //          << "  " << the_ << std::endl;
      //}
      //if( papho->thetaLike() == 0. ) {
      //std::cout << "getRingPk 3  " << CsRichOne::Instance()->kEvent()
      //	  << "  " << the_ << std::endl;
      //}

      if( flag_ ) peakScan();
//                ----------

      if( key->RejWroPeak() ) {
	if( rejectWrongPeak( theMax ) ) flag_ = false;
//          -------------------------
      }

  }


//=========================================================================== 
  void CsRCRing::getRingMaxLike( const int kRing, CsRCPartPhotons* papho ) {
//--------------------------------------------------------------------------

//--- Paolo - September 2004


//--- 'Maximum Likelihood' ring constructor.
//    -------------------------------------

      CsRCExeKeys *key = CsRCExeKeys::Instance();
      bool DoClu = key->DoClu();
      bool UseCluPH = key->UseCluPH();
      int kDoClu = 0;
      if( DoClu ) kDoClu = 0;
      if( !DoClu ) kDoClu = 1;

      CsRCRecConst *cons = CsRCRecConst::Instance();
      int *mcanWind = cons->mcanWind();
      int mcaScan = cons->mcaScan();
      float xlScan = cons->xlScan();
      float xuScan = cons->xuScan();
      float binScan = cons->binScan();
      int nMore = cons->nMore();  ;

      static bool firstCall = true;
      if( firstCall ) {
        firstCall = false;
        key->acknoMethod( "CsRCEventRings::getRingMaxLike" );
      }

      if( !papho->thetaLikeSet() ) {
	flag_ = false;
	//std::cout << "getRingMaxLike 1  " << kRing << std::endl;
	return;
      }
      double thetaLike = papho->thetaLike();
      if( thetaLike <= 0. ) {
	flag_ = false;
	//std::cout << "getRingMaxLike 2  " << kRing << std::endl;
	return;
      }

      list<CsRCPhoton*> lPhotons = papho->lPhotons();
      list<CsRCPhoton*>::iterator ih;

//    sort ring clusters inside window around thetaLike :
//    ---------------------------------------------------
      double thePhotLk = 0.;
      float totCntLk = 0.;
      int mcan = mcanWind[kDoClu];                          //   ???
      float hWind = float( mcan + 2*nMore )* binScan /2.;
      float theMN = thetaLike - hWind;
      float theMX = thetaLike + hWind;
      //cout << theMN << "  " << theMX << endl;
      nPhotPMT_ = 0;
      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
        if( !(*ih)->flag() ) continue;
	double thew = (*ih)->the();
        if( thew >= theMN  &&  thew <= theMX) {
	  lPhotons_.push_back( (*ih) );
//        ----------------------------
          float PHw = 1.;
//@@--------------------
          if( UseCluPH ) PHw = (*ih)->PH() + 1.;
//        -------------------------------------
	  if( (*ih)->isPMT() ) PHw = 1.;
//        -----------------------------
          thePhotLk += thew * PHw;
          totCntLk += PHw;
	  if( (*ih)->isPMT() ) nPhotPMT_++;
	}
      }
      if( totCntLk > 0. ) thePhotLk /= totCntLk;

      if( !lPhotons_.empty() ) {
//----- ring location ( kRingLoc_: up = +1, split = 0, down = -1 ) :
//      ------------------------------------------------------------
        CsRCPhotonDet* detFirst = lPhotons_.front()->ptToDet();
        kRingLoc_ = detFirst->kDet();
        if( kRingLoc_ == 1 ) kRingLoc_ = -1;
        if( kRingLoc_ == 0 ) kRingLoc_ = +1;
        for( ih=lPhotons_.begin(); ih!=lPhotons_.end(); ih++ ) {
	  if( (*ih)->ptToDet() != detFirst ) { kRingLoc_ = 0; break; }
        }

//----- fill RingMaxLike object :
//      -------------------------
        kRing_ = kRing;
        pPartPhot_ = papho;
        the_ = thePhotLk;
	partProbsSet_ = false;
        flagReco_ = false;
	flagOverThrs_ = false;
        flag_ = true;
	flagBack_ = false;
	ringBack_ = 1.;

        papho->setpRing( this );
        papho->pPart()->setpRing( this );

//----- monitoring photons in ring :
//      ----------------------------
        CsRCHistos& hist = CsRCHistos::Ref();
        double xh, wh;
        int nPhotLk = lPhotons_.size();
        xh = nPhotLk;
	//if( hist.hRC???? ) hist.hRC????->Fill( xh );
//hh                       ------------------------

      } else { flag_ = false; }   // end if not empty
      //if( flag_ == false ) {
      //std::cout << "getRingMaxLike 3  " << CsRichOne::Instance()->kEvent()
      //  << "  " << thetaLike << "  " << hWind << std::endl;
      //}
      
      return;
  }


//===========================================================================
  void CsRCRing::histCntRatio() {
//-------------------------------


//--- Paolo  -  May 1999, rev.  October 2000


//--- from "CsRCExecKeys"
      CsRCExeKeys *key = CsRCExeKeys::Instance();
      bool DoClu = key->DoClu();
      int kDoClu = 0;
      if( DoClu ) kDoClu = 0;
      if( !DoClu ) kDoClu = 1;

//--- from "CsRCRecConst"
      CsRCRecConst *cons = CsRCRecConst::Instance();
      int kOffHi = cons->kOffHi();
      int *mcanWind = cons->mcanWind();
      int mcaScan = cons->mcaScan();
      float xuScan = cons->xuScan();
      float binScan = cons->binScan();
      int nScanMx = mcaScan;


//--- monitoring histo.s for the scan for peak in theta :
//    ---------------------------------------------------
      CsRCHistos& hist = CsRCHistos::Ref();
      int kh;
      double xh, wh;

      int mcan = mcanWind[kDoClu];
//@@-----------------------------
      //int mcanmx = nScanMx - mcan + 1;          //   scan up to xuScan mrad
      int mcanmx = nScanMx - mcan;   //   030723
//@@---------------------------------
      if( mcanmx > mcaScan ) mcanmx = mcaScan;

      int ncan1 = 0;
      int ncan2 = mcanmx;
//@@---------------------
      float cntPeakv[mcaScan];
      float cntPeak;
      float cntPeakmx = 0.;
      int kcanmx = -1;
      int kc;
      //for( kc=ncan1; kc<ncan2; kc++ ) {
      for( kc=ncan1; kc<=ncan2; kc++ ) {
        cntPeak = 0.;
        for( int kp=kc; kp<kc+mcan; kp++ ) cntPeak += binCont_[kp];
        cntPeakv[kc] = cntPeak;
        if( cntPeak > cntPeakmx ) { cntPeakmx = cntPeak;  kcanmx = kc; }

        xh = float( kc ) + binScan/2.;
        wh = cntPeak;
	if( hist.hRC1585 ) hist.hRC1585->Fill( xh, wh );
//hh                       ----------------------------
      }

      xh = cntPeakmx;
      if( hist.hRC1581 ) hist.hRC1581->Fill( xh );
//hh                     ------------------------

      float cntPeakav = 0.;
      ncan1 = 0;
      //ncan2 = kcanmx - 1;
      ncan2 = kcanmx;
      for( kc=ncan1; kc<ncan2; kc++ ) cntPeakav += cntPeakv[kc];
      //ncan1 = kcanmx + mcan - 1;
      ncan1 = kcanmx + mcan;
      ncan2 = mcanmx;
      //for( kc=ncan1; kc<ncan2; kc++ ) cntPeakav += cntPeakv[kc];
      for( kc=ncan1; kc<=ncan2; kc++ ) cntPeakav += cntPeakv[kc];
      //cntPeakav = mcan * cntPeakav / (mcanmx - mcan);
      cntPeakav = cntPeakav / (mcanmx - mcan);

      xh = cntPeakav;
      if( hist.hRC1582 ) hist.hRC1582->Fill( xh );
//hh                     ------------------------

      if( cntPeakmx > 0.)  {
        float cntRatio = cntPeakav / cntPeakmx;
        xh = cntRatio;
	if( hist.hRC1583 ) hist.hRC1583->Fill( xh );
//hh                       ------------------------
      }

//    plot theta inside window (+/- nPlot) :
//    --------------------------------------
      int nPlot = 30;
      ncan1 = kcanmx - nPlot;
      if( ncan1 < 0 ) ncan1 = 0;
      ncan2 = kcanmx + mcan + nPlot;
      if( ncan2 > nScanMx ) ncan2 = nScanMx;
      kh = 0;
      for( kc=ncan1; kc<ncan2; kc++ ) {
	xh = kh + binScan/2.;
	wh = binCont_[kc];
	if( hist.hRC1584 ) hist.hRC1584->Fill( xh, wh );
//hh                       ----------------------------
	kh++;
      }

  }


//===========================================================================
  void CsRCRing::peakScan() {
//---------------------------


//--- Paolo  -  September, 2002


//--- from "CsRCExecKeys"
      CsRCExeKeys *key = CsRCExeKeys::Instance();
      bool DoClu = key->DoClu();
      int kDoClu = 0;
      if( DoClu ) kDoClu = 0;
      if( !DoClu ) kDoClu = 1;

//--- from "CsRCRecConst"
      CsRCRecConst *cons = CsRCRecConst::Instance();
      int kOffHi = cons->kOffHi();
      int *mcanWind = cons->mcanWind();
      int mcaScan = cons->mcaScan();
      float xuScan = cons->xuScan();
      float binScan = cons->binScan();
      int nScanMx = mcaScan;


//--- monitoring histo.s for the scan for peak in theta :
//    ---------------------------------------------------
      CsRCHistos& hist = CsRCHistos::Ref();
      int kh;
      double xh, wh;

      int mcan = mcanWind[kDoClu];
      //int mcan = getVarMcan( papho, mcanWind[kDoClu] );     //   040213
//@@-----------------------------
      //int mcanmx = nScanMx - mcan + 1;          //   scan up to xuScan mrad
      int mcanmx = nScanMx - mcan;   //   030723
//@@---------------------------------
      if( mcanmx > mcaScan ) mcanmx = mcaScan;

      int ncan1 = 0;
      int ncan2 = mcanmx;
//@@---------------------
      float cntPeakv[mcaScan];
      float cntPeak;
      float cntPeakmx = 0.;
      int kcanmx = -1;
      int kc;
      for( kc=ncan1; kc<=ncan2; kc++ ) {
        cntPeak = 0.;
        for( int kp=kc; kp<kc+mcan; kp++ ) cntPeak += binCont_[kp];
        cntPeakv[kc] = cntPeak;
        if( cntPeak > cntPeakmx ) { cntPeakmx = cntPeak;  kcanmx = kc; }
      }

      float theMax = (float( kcanmx ) + float( mcan )/2. + 0.5)* binScan;
      double theIpo[5];
      double *massPartv = cons->massPartv();
      float CFRefInd = cons->CFRefInd();
      if( !pPartPhot_ ) return;
      CsRCParticle* part = pPartPhot_->pPart();
      float momPart = part->mom();
      int kPart = 0;
      //cout << CsRichOne::Instance()->kEvent() << endl;
      for( int kPaTy=8; kPaTy <= 14; kPaTy+=3 ) {
        theIpo[kPart] = 0.;
        float massIpo = massPartv[kPaTy];
        double betaIpo = momPart / sqrt( momPart*momPart + massIpo*massIpo );
        double cosTheW = 1./ (betaIpo * CFRefInd);
        if( cosTheW > 1.) cosTheW = 1.;
        theIpo[kPart] = acos( cosTheW ) * 1000.;
        //cout << momPart << "  " <<  massIpo << "  "
        //     << theIpo[kPart] << endl;
        kPart++;
      }

      static int kHHMx = 20;
      if( !CsRCHistos::Ref().bookHis() ||
	  CsRCHistos::Ref().levelBk() < 2 ) kHHMx = 0;
      static vector<CsHist1D*> vRC6000;
      static int kPlotB = 0;
      static bool bPlotB = true;
//@@---------------------------
      if( bPlotB  &&  kPlotB < kHHMx ) {
	CsHistograms::SetCurrentPath("/RICH");
	int khPlots = kOffHi + 6000 + kPlotB;
	stringstream name0, title0;
	name0 << khPlots;
	title0 << "   ";
        vRC6000.push_back( new CsHist1D( name0.str(), title0.str(),
//hh    -----------------------------------------------------------
					 //mcanmx, 0., mcanmx*binScan ) );
					 mcaScan, 0., mcaScan*binScan ) );
	int ncan1 = 0;
	int ncan2 = mcaScan;
//@@-----------------------
	xh = (float( ncan1 ) + 0.5)* binScan;
	//wh = theMax;
	wh = the_;
	if( vRC6000[kPlotB] ) vRC6000[kPlotB]->Fill( xh, wh );
//hh                          -------------------------------
        for( int kPa=0; kPa<3; kPa++ ) {
	  xh = (ncan1+kPa+1 + 0.5)* binScan;
	  wh = theIpo[kPa];
	  if( vRC6000[kPlotB] ) vRC6000[kPlotB]->Fill( xh, wh );
//hh                            -------------------------------
	}
	for( int kh=ncan1+4; kh<ncan2; kh++ ) {
	  xh = (kh + 0.5)* binScan;
	  wh = binCont_[kh];
	  if( vRC6000[kPlotB] ) vRC6000[kPlotB]->Fill( xh, wh );
//hh                            -------------------------------
	}
	CsHistograms::SetCurrentPath("");
	kPlotB++;
      }

      if( !CsRCHistos::Ref().bookHis() ||
	  CsRCHistos::Ref().levelBk() < 2 ) kHHMx = 0;
      static vector<CsHist1D*> vRC6100;
      static int kPlotP = 0;
      static bool bPlotP = true;
//@@---------------------------
      if( bPlotP  &&  kPlotP < kHHMx ) {
	CsHistograms::SetCurrentPath("/RICH");
	int khPlots = kOffHi + 6100 + kPlotP;
	stringstream name0, title0;
	name0 << khPlots;
	title0 << "   ";
        vRC6100.push_back( new CsHist1D( name0.str(), title0.str(),
//hh    -----------------------------------------------------------
					 mcanmx, 0., mcanmx*binScan ) );
	int ncan1 = 0;
	int ncan2 = mcanmx;
//@@-----------------------
	xh = (float( ncan1 ) + 0.5)* binScan;
	//wh = theMax;
	wh = the_;
	if( vRC6100[kPlotP] ) vRC6100[kPlotP]->Fill( xh, wh );
//hh                          -------------------------------
        for( int kPa=0; kPa<3; kPa++ ) {
	  xh = (ncan1+kPa+1 + 0.5)* binScan;
	  wh = theIpo[kPa];
	  if( vRC6100[kPlotP] ) vRC6100[kPlotP]->Fill( xh, wh );
//hh                            -------------------------------
	}
	for( int kh=ncan1+4; kh<ncan2; kh++ ) {
	  xh = (kh + 0.5)* binScan;
	  wh = cntPeakv[kh];
	  if( vRC6100[kPlotP] ) vRC6100[kPlotP]->Fill( xh, wh );
//hh                            -------------------------------
	}
	CsHistograms::SetCurrentPath("");
	kPlotP++;
      }

      //bool ringFlag = false;
//@@--//---------------------
      //float hWind = 4.;
      //for( int kPa=0; kPa<3; kPa++ ) {
	//cout << theMax << "  " << theIpo[kPa]-hWind << "  "
	//     << theIpo[kPa]+hWind << endl;
      //if( theMax >= theIpo[kPa]-hWind  &&  theMax < theIpo[kPa]+hWind )
        //ringFlag = true;
      //}
      //flag_ = ringFlag;

  }


//===========================================================================
  bool CsRCRing::rejectWrongPeak( const float theMax ) {
//------------------------------------------------------


//--- Paolo  -  September, 2002


      CsRCRecConst *cons = CsRCRecConst::Instance();

      int kTyLL = 2;
//@@---------------

      bool ringFlag = true;      //   --> reject Ring

      float hWind = 4.;
//@@------------------
      double *massPartv = cons->massPartv();
      float CFRefInd = cons->CFRefInd();
      if( !pPartPhot_ ) return ringFlag;
      CsRCParticle* part = pPartPhot_->pPart();
      float momPart = part->mom();
      for( int kPaTy=kTyLL; kPaTy <= 14; kPaTy+=3 ) {
//    ---------------------------------------------
        float massIpo = massPartv[kPaTy];
        double betaIpo = momPart / sqrt( momPart*momPart + massIpo*massIpo );
        double cosTheW = 1./ (betaIpo * CFRefInd);
        if( cosTheW > 1.) cosTheW = 1.;
        double theIpo = acos( cosTheW ) * 1000.;
        //cout << momPart << "  " <<  massIpo << "  " << theIpo << endl;
	if( theMax >= theIpo-hWind  &&  theMax < theIpo+hWind ) 
	  ringFlag = false;
      }

      return ringFlag;

  }


//=========================================================================== 
  bool CsRCRing::peakCountSearch( const CsRCPartPhotons* papho, int &kcMax )
//--------------------------------------------------------------------------
  {


//--- Paolo - November 2002

//--- standard version


      CsRCExeKeys *key = CsRCExeKeys::Instance();
      static bool firstCall = true;
      if( firstCall ) {
	firstCall = false;

        key->acknoMethod( "CsRCRing::peakCountSearch" );
      }

      bool DoClu = key->DoClu();
      bool UseCluPH = key->UseCluPH();
      int kDoClu = 0;
      if( DoClu ) kDoClu = 0;
      if( !DoClu ) kDoClu = 1;

      CsRCRecConst *cons = CsRCRecConst::Instance();
      int *mcanWind = cons->mcanWind();
      int mcaScan = cons->mcaScan();
      float xlScan = cons->xlScan();
      float xuScan = cons->xuScan();
      float binScan = cons->binScan();
      int nScanMx = mcaScan;

      for( int k=0; k<mcaScan; k++ ) binCont_[k] = 0.;
      //memset( binCont_, 0, sizeof( binCont_ ) );
      //if( mcaScan > sizeof( binCont_ ) )
      //std::cout << mcaScan << "  " << sizeof( binCont_ ) << std::endl;

//--- bin the theta distribution (counts only or PH) :
//    ------------------------------------------------
      list<CsRCPhoton*> lPhotons = papho->lPhotons();
      list<CsRCPhoton*>::iterator ih;
      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
        double theCluR = (*ih)->the();
        float PHw = 1.; 
//@@------------------
	if( UseCluPH ) PHw = (*ih)->PH() + 1.; // 020903
//      -------------------------------------
	if( (*ih)->isPMT() ) PHw = 1.;
//      -----------------------------
        int kCan = int( (theCluR - xlScan) / binScan );
        if( kCan >= 0  &&  kCan < mcaScan ) binCont_[kCan] += PHw; // 020909
      }

//--- scan for peak in theta :
//    ------------------------
      int mcan = mcanWind[kDoClu];                        //   040213
      //int mcan = getVarMcan( papho, mcanWind[kDoClu] );     //   040213
//@@------------------------------
      //int mcanmx = nScanMx - mcan + 1;   //   scan up to xuScan mrad
      int mcanmx = nScanMx - mcan;   //   030723
//@@-----------------------------
      if( mcanmx > mcaScan ) mcanmx = mcaScan;
      int ncan1 = 0;
      int ncan2 = mcanmx;
//@@---------------------
      float cntPeak;
      float cntPeakmx = 0.;
      int kcanmx = -1;
      for( int kc=ncan1; kc<=ncan2; kc++ ) {
        cntPeak = 0.;
        for( int kp=kc; kp<kc+mcan; kp++ ) cntPeak += binCont_[kp];
        if( cntPeak > cntPeakmx ) { cntPeakmx = cntPeak;  kcanmx = kc; }
        //test for V.Frolov - 030723
	for( int kp=kc; kp<kc+mcan; kp++ ) if( kp >= mcaScan )
	  std::cout << mcaScan << "  " << kp << "  " <<
	  binCont_[kp] << std::endl;
      }
      kcMax = kcanmx;

      return true;

  }


//=========================================================================== 
  void CsRCRing::histPeakSearch( const CsRCPartPhotons* papho ) {
//---------------------------------------------------------------


//- Paolo - January  2003


    CsRCRecConst *cons = CsRCRecConst::Instance();
    int mcaScan = cons->mcaScan();
    float binScan = cons->binScan();
    int nScanMx = mcaScan;
    int kOffHi = cons->kOffHi();

    CsRCExeKeys *key = CsRCExeKeys::Instance();
    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      key->acknoMethod( "CsRCRing::histPeakSearch" );
    }

    CsRCHistos& hist = CsRCHistos::Ref();
    int kh;
    double xh, yh, wh;
    static int kHHMx = 20;
    if( !CsRCHistos::Ref().bookHis() ||
	CsRCHistos::Ref().levelBk() < 2 ) kHHMx = 0;
    static vector<CsHist2D*> vRC6300;
    static int kPlotP = 0;
    static bool bPlotP = true;
//@@-------------------------
    if( bPlotP  &&  kPlotP < kHHMx ) {
      CsHistograms::SetCurrentPath("/RICH");
      int khPlots = kOffHi + 6300 + kPlotP;
      stringstream name0, title0;
      name0 << khPlots;
      title0 << "   ";
      vRC6300.push_back( new CsHist2D( name0.str(), title0.str(),
//hh  -----------------------------------------------------------
				       nScanMx, 0., nScanMx*binScan,
				       90, 0., 360. ) );

      CsHistograms::SetCurrentPath("");

      list<CsRCPhoton*> lPhotons = papho->lPhotons();
      list<CsRCPhoton*>::iterator ih;
      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {

	xh = (*ih)->the();
	yh = (*ih)->phi();
	//wh = (*ih)->PH();
	wh = 1.;
	if( vRC6300[kPlotP] ) vRC6300[kPlotP]->Fill( xh, yh, wh );
//hh                          -----------------------------------
      }

    }
    kPlotP++;
  }


//===========================================================================
  void CsRCRing::setThetaType() {
//-------------------------------


//- Paolo  - December  2002


    CsRCRecConst* cons = CsRCRecConst::Instance();

//- set Ring limits :
//  -----------------
//  needed by the angle computations
    double momPart = pPartPhot_->pPart()->mom();
    double sigPhoRec = pPartPhot_->sigmaPhoRec( momPart );
//@@                               ----------------------
    float dSigC = cons->sigCut() * sigPhoRec;
//@@              --------------
    theLoL_ = the_ - dSigC;
    theUpL_ = the_ + dSigC;

//- set reconstructed Ring Cerenkov angle :
//  ---------------------------------------
    double theReco = the_;                             //   default

    if( cons->thetaType() == "REC" ) theReco = the_;
//                                   --------------
    setPartProb( 8, the_ );
//@@----------------------

    thetaWgav_ = getThetaWgav();
//               --------------
    if( thetaWgav_ == 0. ) thetaWgav_ = the_;
    if( cons->thetaType() == "WAVE" ) theReco = thetaWgav_;
//                                    --------------------

    thetaRFit_ = getThetaRFit();
//               --------------
    if( thetaRFit_ == 0. ) thetaRFit_ = the_;
    if( cons->thetaType() == "FIT" ) theReco = thetaRFit_;
//                                   --------------------
    setPartProb( 10, thetaRFit_ );
//@@-----------------------------

    theReco_ = theReco;

//- Redefine Ring limits :
//  ----------------------
//  according to the angle type used
    //theLoL_ = theReco - dSigC;        //   out 030303
    //theUpL_ = theReco + dSigC;        //   out 030303

    //cout << the_ << "  " << thetaWgav_ << "  " << thetaRFit_
    // << "  " << theReco << endl;

  }



//===========================================================================
  double CsRCRing::getThetaWgav() {
//---------------------------------


//- theta ring from weighted average
//  -----------------------------------
//- Paolo  -  January 2001
//    rev.    December 2002


    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;
      CsRCExeKeys::Instance()->acknoMethod( "CsRCRing::getThetaWgav" );
    }

    CsRCHistos& hist = CsRCHistos::Ref();
    double xh, yh, wh;

    CsRCRecConst * cons = CsRCRecConst::Instance();
    static float CFRefInd = cons->CFRefInd();
    double theReco = the_;

    list<CsRCPhoton*> lPhotons = lPhotons_;
    list<CsRCPhoton*>::iterator ih;
    double coSig = 0.;
    double coBkg = 0.;
    for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
      double thePhot = (*ih)->the();
      if( thePhot >= theLoL_  &&  thePhot <= theUpL_ )  {
	double sigPhot = (*ih)->sigmaPhoPid( pPartPhot_ );
//@@                     --------------------------------
	coSig += thePhot/(sigPhot*sigPhot);
	coBkg += 1. /(sigPhot*sigPhot);
      }
    }   /* end for on Photons */

    //double thetaCer = pPartPhot_->pPart()->thetaCer();

    double thetaWgav = 0.;
    if( coBkg > 0. ) thetaWgav = coSig / coBkg;

//- 'pulls' :
//  ---------
    for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
      double thePhot = (*ih)->the();
      if( thePhot >= theLoL_  &&  thePhot <= theUpL_ )  {
	double sigPhot = (*ih)->sigmaPhoPid( pPartPhot_ );
//@@                     --------------------------------
	double den = 0.;
	if( coBkg != 0. ) den = sigPhot*sigPhot - 1./coBkg;
	double pull = 1000000.;
	if( den > 0. ) {
	  pull = (thePhot - thetaWgav) / sqrt( den );
	  xh = pull;
	  if( hist.hRC1632 ) hist.hRC1632-> Fill( xh );
//hh                         -------------------------
	}
      }
    }   /* end for on Photons */

    return  thetaWgav;

  }


//===========================================================================
  double CsRCRing::getThetaRFit() {
//---------------------------------

//- theta ring from fit to the ring
//  -------------------------------
//- Paolo  -  April  2001
//    rev.    December 2002


    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      CsRCExeKeys::Instance()->acknoMethod( "CsRCRing::getThetaRFit" );
    }

    CsRCHistos& hist = CsRCHistos::Ref();
    double xh, yh, wh;

    CsRCRecConst* cons = CsRCRecConst::Instance();
    double radDeg = cons->RadDeg();
    int nPhotMinRing = cons->nPhotMinRing();

    double thetaFit = 0.;

    list<CsRCPhoton*> lPhotons = lPhotons_;
    int nPhotons = lPhotons.size();
    if( nPhotons > 0 ) {

      int kPhot = 0;
      double thePhot[nPhotons];
      double phiPhot[nPhotons];
      double sigPhot[nPhotons];

      double theReco = the_;

      list<CsRCPhoton*>::iterator ih;
      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {

	double theW = (*ih)->the();
        if( theW >= theLoL_  &&  theW <= theUpL_ ) {
	  double phiW = (*ih)->phi();
	  double sigW = (*ih)->sigmaPhoPid( pPartPhot_ );
//@@                    --------------------------------
	  thePhot[kPhot] = theW * cos( phiW/radDeg );
	  phiPhot[kPhot] = theW * sin( phiW/radDeg );
	  sigPhot[kPhot] = sigW*sigW;
	  kPhot++;
	}
      }   /* end for on Photons */
      int nPoint = kPhot;

//--- fit ring in S.Y. plane :
//    ------------------------
      bool exeFit = true;
      if( nPoint < nPhotMinRing ) exeFit = false;
//@@--------------------------------------------
      if( exeFit ) {
	int nParam = 3;
	double param[nParam];
	param[0] = the_;
	param[1] = 0.;
	param[2] = 0.;
	int iPaFit[nParam];
	iPaFit[0] = 0;
	iPaFit[1] = 0;
	iPaFit[2] = 0;

        CsRCCircleFit oCirclev( nPoint, thePhot, phiPhot, sigPhot,
//      ----------------------------------------------------------
                                nParam, param, iPaFit );

        if( oCirclev.doChiFit() ) {
//          -------------------

	  thetaFit = oCirclev.para()[0];
          if( thetaFit > cons->xuScan() ) thetaFit = 0.;  //   021010
	  if( thetaFit < theLoL_  ||  thetaFit > theUpL_ ) thetaFit = 0.;
//@@      --------------------------------------------------------------

	  //oCirclev.print();
	  oCirclev.doHist( hist.vRC3700 );
//        -------------------------------

        } else {
	}

      }      /*   end if exeFit   */

    }

    return  thetaFit;

  }


//===========================================================================
  long double CsRCRing::getLikelihood( const double theReco,
//----------------------------------------------------------
				       const double theIpo,
				       int &nPhotons ) {

//- interface method
//  ----------------
//- Paolo  -  December 2002


    long double like = getLikeRing( theReco, theIpo, nPhotons );
//                     ----------------------------------------

    return  like;

  }


  extern "C" { float erf_( const float& ); }

//===========================================================================
  long double CsRCRing::getLikeRing( const double theReco,
//--------------------------------------------------------
				     const double theIpo,
				     int &nPhoRing ) {


//- compute  Likelihood value for a given hypothesis (ring photons only)
//  --------------------------------------------------------------------
//- Paolo  -  December 2002


    CsRCExeKeys *key = CsRCExeKeys::Instance();
    CsRCRecConst * cons = CsRCRecConst::Instance();

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      key->acknoMethod( "CsRCRing::getLikelihood(Ring)" );
    }

    CsRCLikeRing* likeRing = NULL;
    if( cons->backgrType() == "02" ) likeRing = new CsRCLikeRing02( this );
    if( cons->backgrType() == "03" ) likeRing = new CsRCLikeRing03( this );
    if( cons->backgrType() == "04" ) likeRing = new CsRCLikeRing04( this );
    if( cons->backgrType() == "05" ) likeRing = new CsRCLikeRing05( this );
//@@----------------------------------------------------------------------
    if( !likeRing ) {
      cout << " RICHONE, CsRCRing::getLikeRing() :   ";
      cout << "reference to a NOT existing background parametrization!";
      cout << endl;
      string mess = "RICHONE, CsRCRing::getLikeRing() : ";
      string err = "Reference to a NOT existing background parametrization!";
      mess.append( err );
      //CsErrLog::Instance()->mes( elError, mess );
      CsErrLog::Instance()->mes( elFatal, mess );
    }
    bool likeRatio = cons->likeRatio();
    if( cons->backgrType() == "04" ) likeRatio = true;
    if( cons->backgrType() == "05" ) likeRatio = true;

    double normS = likeRing->normSignal( theIpo );
//                 ------------------------------
    double normB = likeRing->normBackgr( theReco, theIpo );
//                 ---------------------------------------
    double likeNorm = normS + normB;

    //ringBack_ = likeRing->getRingBackground( theReco );
//@@--------------------------------------------------
    //flagBack_ = true;
//@@----------------

    list<CsRCPhoton*> lPhotons = lPhotons_;
    list<CsRCPhoton*>::const_iterator ih;

    long double pLike = -1.;
    long double pLikeBk = -1.;                     //   030609
    int kPhot = 0;
    for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
      double thePhot = (*ih)->the();
      if( thePhot >= theLoL_  &&  thePhot <= theUpL_ ) {
	double sigPhot = (*ih)->sigmaPhoPid( pPartPhot_ );
//                       --------------------------------

	double signW = likeRing->likeSignal( thePhot, theIpo, sigPhot );
//                     ------------------------------------------------
	//double backW = likeRing->likeBackgr( thePhot, theIpo );
        double backW = likeRing->likeBackgr( (*ih), theIpo );    //   040511
//                     -------------------------------------

	//backW = 0.;                                //   040505
	//double backWgt = 0.;                       //   040505
	//if( (*ih)->ptToClu()->getBackWgt( backWgt ) ) {
//@@    //---------------------------------------------
	//  backW = 0.095 * backWgt * thePhot;       //   040505
	//}                                          //   040505

	double likeW = signW + backW;

	pLike *= likeW;
	pLike = fabs( pLike );
	pLikeBk *= backW;                          //   030609
	pLikeBk = fabs( pLikeBk );                 //   030609
	kPhot++;
      }
    }   /* end for on Photons */
    //cout << kPhot << "  " << pLike << endl;

    nPhoRing = kPhot;
    if( nPhoRing > 0 ) {
      long double power = 1./float( nPhoRing );
      if( pLike > 0. ) pLike = pow( pLike, power );
      if( likeNorm > 0. ) {
	pLike /= likeNorm;
      } else {
	pLike = -1.;
      }
      if( pLikeBk > 0. ) pLikeBk = pow( pLikeBk, power );     //   030609
      //if( normB > 0. ) {                                      //
      //  pLikeBk /= normB;                                     //
      if( likeNorm > 0. ) {                                   //   040510
        pLikeBk /= likeNorm;                                  //   040510
      } else {                                                //
	pLikeBk = -1.;                                        //
      }                                                       //
      if( likeRatio ) {                                       //
	if( pLikeBk > 0. ) {                                  //
	  pLike /= pLikeBk;                                   //
	} else {                                              //
	  pLike = -1.;                                        //   030609
	}
      }
      //if( pLike < 1. ) cout << "  " << pLike << endl;
    } else {
      pLike = -1.;
    }
    //std::cout << pLike << std::endl;

    delete likeRing;

    return  pLike;

  }


//===========================================================================
  bool CsRCRing::getThetaLikeMax( long double& likeMax, double& theMax ) {
//------------------------------------------------------------------------


//- theta ring from max. Likelihood
//  -------------------------------
//- Paolo  -  September 2002
//    rev.    December  2002 
//    rev.    February  2003
//    rev.    October   2004


    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      CsRCExeKeys::Instance()->acknoMethod( "CsRCRing::getThetaLikeMax" );
    }

    CsRCRecConst * cons = CsRCRecConst::Instance();

    double theReco = the_;

//- Coarse search ( from checkLikelihood - 02/11/14 )
//  -------------------------------------------------
    theMax = 0.;
    double likeMaxW = -1.;
    double theMaxW = 0.;
    long double likeM = 0.;
    double theLiMin = 5.;
    double theLiMax = 55.;
    int nIpo = 25;
    double dThe = (theLiMax - theLiMin) / float( nIpo );
    long double likeV[nIpo];
    int kIpoV = 0;
    int nPhot = 0;
    for( int kIpo=0; kIpo<nIpo; kIpo++ ) {
      double theIpo = theLiMin + kIpo * dThe + dThe/2.;
      likeM = getLikelihood( theReco, theIpo, nPhot );
//            ---------------------------------------
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
    if( likeMax < 0. ) return  false;
    if( theMax <= 0. ) return  false;
    if( theMax >= theLiMax ) { theMax = theLiMax; return  true; }


//- Fine search
//  -----------
    if(kIpoV <= 0 ) kIpoV = 1;
    if(kIpoV >= nIpo-1 ) kIpoV = nIpo - 2;
    //cout << endl;
    //cout << "*** " << theMax << "  " << likeMaxW << "  "
    //     << likeV[kIpoV-1] << "  " << likeV[kIpoV+1] << endl;
    double delta = dThe;
    double theMs = theMax - delta;
    double theCt = theMax;
    double thePs = theMax + delta;
    double likeMs = likeV[kIpoV-1];
    if( likeMs < 0. ) return  true;
    double likeCt = likeMax;
    double likePs = likeV[kIpoV+1];
    if( likePs < 0. ) return  true;
    double den = (likePs - 2.* likeCt + likeMs);                 //   041007
    if( den == 0. ) return  true;                                //   041007
    theMax = theCt - delta/2. * (likePs - likeMs) / den;         //   041007
//  ---------------------------------------------------
    if( theMax <= 0. ) return  false;
    if( theMax < theMs ) theMax = theMs;
    if( theMax > thePs ) theMax = thePs;
    if( theMax >= theLiMax ) {
      likeMax = getLikelihood( theMax, theMax, nPhot );
      return  true;
    }

    delta /= 3.;
    //cout << "------ " << theMax << "  " << delta << endl;
    for( int kRep=0; kRep<3; kRep++ ) {
      int nPhot = 0;
      double theMs = theMax - delta;
      double theCt = theMax;
      double thePs = theMax + delta;
      double likeMs = getLikelihood( theReco, theMs, nPhot );
      if( likeMs < 0. ) break;
      double likeCt = getLikelihood( theReco, theCt, nPhot );
      if( likeCt < 0. ) break;
      double likePs = getLikelihood( theReco, thePs, nPhot );
      if( likePs < 0. ) break;
      double den = (likePs - 2.* likeCt + likeMs);               //   040513
      if( den == 0. ) break;                                     //   040513
      theMax = theCt - delta/2. * (likePs - likeMs) / den;       //   040513
//    ---------------------------------------------------
      if( theMax <= 0. ) return  false;
      if( theMax < theMs ) theMax = theMs;
      if( theMax > thePs ) theMax = thePs;
      if( theMax >= theLiMax ) break;
      delta /= 3.;
    }
    likeMax = getLikelihood( theMax, theMax, nPhot );
    if( likeMax < 0. ) return  false;
    //cout << "=== " << likeMax << "  " << theMax << "  " << delta << endl;

    return  true;

  }


//===========================================================================
  void CsRCRing::checkChiSquare( const double theReco ) {
//-------------------------------------------------------


//- Paolo  -  August 2003


    CsRCRecConst *cons = CsRCRecConst::Instance();
    int kOffHi = cons->kOffHi();

    static vector<CsHist1D*> vRC3120;
    static int kH3120 = 0;
    static int kHH = 0;
    static int kHHMx = 20;
    if( !CsRCHistos::Ref().bookHis() ||
	CsRCHistos::Ref().levelBk() < 2 ) kHHMx = 0;
    //static double dTheRe = 4.;
    //static int nIpoRe = 80;

    static int nPhoMin = 15;
//@@-----------------------

    int kTyLL = 2;
//@@-------------

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      CsRCExeKeys::Instance()->acknoMethod( "CsRCRing::checkChiSquare" );
    }

    if( kHH > kHHMx ) return;

    int nPhoRing = lPhotons_.size();
    int iPartTy = 0;
    if( CsRCExeKeys::Instance()->MCarloEvent() ) {
      iPartTy = pPartPhot_->pPart()->iPartT();
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
      int nIpo = 321;
      double dThe = (theMax - theMin)/float( nIpo );

      if( kHH < kHHMx ) {
	int kHist = kOffHi + 3120 + kHH;
        stringstream hN3120;
	hN3120 << kHist;
	string hTitle = "likelihood shape";
	CsHistograms::SetCurrentPath("/RICH");
        vRC3120.push_back( new CsHist1D( hN3120.str(), hTitle,
//hh    ------------------------------------------------------
					 nIpo, theMin, theMax ) );
	CsHistograms::SetCurrentPath("/");
      }

      double chiSq[nIpo];
      for( int kIpo=0; kIpo<nIpo; kIpo++ ) {

        theIpo = theMin + kIpo * dThe + dThe/2.;
        double chiSquare = getQsquare( theIpo, nPhot );
//                         ---------------------------
	chiSq[kIpo] = exp( - 0.5* chiSquare );
	//cout << chiSq[kIpo] << endl;
      }

      if( kHH < kHHMx ) {
	int nIpo6 = nIpo - 6;
        for( int kIpo=0; kIpo<nIpo6; kIpo++ ) {
          theIpo = theMin + kIpo * dThe + dThe/2.;
	  double xh = theIpo;
          double wh = chiSq[kIpo];
	  if( vRC3120[kHH] ) vRC3120[kHH]->Fill( xh, wh );
//hh                         ----------------------------
	}
        int iIpo = 6;
	float momPart = pPartPhot_->pPart()->mom();
	double xh = theMin + (nIpo-iIpo) * dThe + dThe/2.;
	double wh = momPart;
	if( vRC3120[kHH] ) vRC3120[kHH]->Fill( xh, wh );
//hh                       ----------------------------
        iIpo = 5;
        for( int kPaTy=kTyLL; kPaTy<=14; kPaTy+=3 ) {
//      -------------------------------------------
	  float massIpo = CsRCRecConst::Instance()->massPartv()[kPaTy];
          double betaIpo = momPart /
            sqrt( momPart*momPart + massIpo*massIpo );
          double cosTheW = 1./(betaIpo* CsRCRecConst::Instance()->CFRefInd());
	  theIpo = -1.;
          if( cosTheW <= 1.) theIpo = acos( cosTheW ) * 1000.;
	  double xh = theMin + (nIpo-iIpo) * dThe + dThe/2.;
          double wh = theIpo;
          if( vRC3120[kHH] ) vRC3120[kHH]->Fill( xh, wh );
//hh                         ----------------------------
	  //cout << kHH << "  " << momPart << endl;
          iIpo--;
	}
      }
      kHH++;

    }

  }


//===========================================================================
  void CsRCRing::checkLikelihood( const double theReco ) {
//--------------------------------------------------------


//- Paolo  -  July 2001
//    rev.    December 2002
//    rev.    January  2004


    CsRCRecConst *cons = CsRCRecConst::Instance();
    int kOffHi = cons->kOffHi();

    static vector<CsHist1D*> vRC3080;
    static int kH3080 = 0;
    static int kHH = 0;
    static int kHHMx = 20;
    static int kHHt = 0;
    static int kHHtMx = 1000;
    if( !CsRCHistos::Ref().bookHis() ||
	CsRCHistos::Ref().levelBk() < 2 ) kHHMx = 0;
    static CsHist1D* hRC3070;
    static double dTheRe = 4.;
    static int nIpoRe = 80;
    static int nPhoRing = 0;

    static int nPhoMin = 15;
//@@-----------------------

    int kTyLL = 2;
//@@-------------

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      CsRCExeKeys::Instance()->acknoMethod( "CsRCRing::checkLikelihood" );

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

    //if( kHH > kHHMx ) return;
    if( kHH < kHHMx ) {

      nPhoRing = lPhotons_.size();
      int iPartTy = 0;
      if( CsRCExeKeys::Instance()->MCarloEvent() ) {
        iPartTy = pPartPhot_->pPart()->iPartT();
        if( iPartTy > 200 ) iPartTy -= 200;
        if( iPartTy > 30 ) iPartTy = 0;   //   wrong myFile?
//@@----------------------------------
      }
      if( nPhoRing > nPhoMin ) {
        //if( iPartTy ==  8 || iPartTy ==  9 ) {     // MC pions only
        //if( iPartTy == 11 || iPartTy == 12 ) {     // MC kaons only
        //if( iPartTy == 14 || iPartTy == 15 ) {     // MC protons only
        double theIpo = 0.;
        int nPhot = 0;
        double theMin = 0.;
        double theMax = 80.;
        int nIpo = 321;                            // corrected 020716 !
        double dThe = (theMax - theMin)/float( nIpo );

	//if( kHH < kHHMx ) {
	int kHist = kOffHi + 3080 + kHH;
        stringstream hN3080;
	hN3080 << kHist;
	string hTitle = "likelihood shape";
	CsHistograms::SetCurrentPath("/RICH");
        vRC3080.push_back( new CsHist1D( hN3080.str(), hTitle,
//hh    ------------------------------------------------------
					 nIpo, theMin, theMax ) );
	CsHistograms::SetCurrentPath("/");
	//}
        double pLike0 = -1.;
        double theIpo0 = -1.;
        double pLike[nIpo];
        for( int kIpo=0; kIpo<nIpo; kIpo++ ) {
	  //cout << kIpo << "  ";
          theIpo = theMin + kIpo * dThe + dThe/2.;
          pLike[kIpo] = getLikelihood( theReco, theIpo, nPhot );
//                      ---------------------------------------
  	  //cout << pLike[kIpo] << endl;
          if( pLike[kIpo] > pLike0 ) {
	    pLike0 = pLike[kIpo];
	    theIpo0 = theIpo;
	  }
	}
	//if( kHH < kHHMx ) {
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
	float momPart = pPartPhot_->pPart()->mom();
	double xh = theMin + (nIpo-iIpo) * dThe + dThe/2.;
	double wh = momPart;
	if( vRC3080[kHH] ) vRC3080[kHH]->Fill( xh, wh );
//hh                       ----------------------------
        //iIpo = 3;
        iIpo = 5;
        for( int kPaTy=kTyLL; kPaTy<=14; kPaTy+=3 ) {
//      -------------------------------------------
	  float massIpo = CsRCRecConst::Instance()->massPartv()[kPaTy];
          double betaIpo = momPart /
            sqrt( momPart*momPart + massIpo*massIpo );
          double cosTheW = 1./(betaIpo* CsRCRecConst::Instance()->CFRefInd());
	  theIpo = -1.;
          if( cosTheW <= 1.) theIpo = acos( cosTheW ) * 1000.;
	  double xh = theMin + (nIpo-iIpo) * dThe + dThe/2.;
          double wh = theIpo;
          if( vRC3080[kHH] ) vRC3080[kHH]->Fill( xh, wh );
//hh                         ----------------------------
	  //cout << kHH << "  " << momPart << endl;
          iIpo--;
	}
	//}
	kHH++;
      }
    }
    //kHH++;

    if( kHHt < kHHtMx ) {
      nPhoRing = lPhotons_.size();
      if( nPhoRing > nPhoMin ) {
        double theIpo = 0.;
        int nPhot = 0;
        double theMin = theReco - dTheRe;
        double theMax = theReco + dTheRe;
        double dThe = (theMax - theMin)/float( nIpoRe );
        for( int kIpo=0; kIpo<=nIpoRe; kIpo++ ) {
          theIpo = theMin + kIpo * dThe;
          long double pLike = -1.;
          pLike = getLikelihood( theReco, theIpo, nPhot );
//                ---------------------------------------
	  //pLike /= pLike0;
       	  double xh = - dTheRe + kIpo * dThe;       // !
	  xh -= 0.05;
          double wh = pLike;
          if( hRC3070 ) hRC3070->Fill( xh, wh );
//hh                    -----------------------
	}
	kHHt++;
      }
    }

  }


//===========================================================================
  void CsRCRing::checkLikeVsIndex( const double theReco ) {
//---------------------------------------------------------


//- Paolo  -  August 2003


    CsRCRecConst *cons = CsRCRecConst::Instance();
    int kOffHi = cons->kOffHi();
    float CFRefInd = cons->CFRefInd();

    static vector<CsHist1D*> vRC3160;
    static int kH3160 = 0;
    static int kHH = 0;
    static int kHHMx = 20;
    if( !CsRCHistos::Ref().bookHis() ||
 	 CsRCHistos::Ref().levelBk() < 2 ) kHHMx = 0;

    double xh, yh, wh;

    static int nPhoMin = 15;
//@@-----------------------

    int kTyLL =  8;
    int kTyMM =  8;
//@@--------------

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      CsRCExeKeys::Instance()->acknoMethod( "CsRCRing::checkLikeVsIndex" );
    }

    if( kHH > kHHMx ) return;

    int nPhoRing = lPhotons_.size();
    int iPartTy = 0;
    if( CsRCExeKeys::Instance()->MCarloEvent() ) {
      iPartTy = pPartPhot_->pPart()->iPartT();
      if( iPartTy > 200 ) iPartTy -= 200;
      if( iPartTy > 30 ) iPartTy = 0;   //   wrong myFile?
//@@--------------------------------
    }
    if( nPhoRing > nPhoMin ) {

      //if( iPartTy ==  8 || iPartTy ==  9 ) {     // pions only
      //if( iPartTy == 11 || iPartTy == 12 ) {     // kaons only
      //if( iPartTy == 14 || iPartTy == 15 ) {     // protons only

      float refIndM1 = CFRefInd - 1.;
      float fracIndM1 = 0.15;
//@@------------------------
      float CFRefIMin = 1. + refIndM1 * ( 1.- fracIndM1 );
      float CFRefIMax = 1. + refIndM1 * ( 1.+ fracIndM1 );
      //std::cout << CFRefIMin << "  " << CFRefIMax << std::endl;
      double theIpo = 0.;
      int nPhot = 0;
      double theMin = 0.;
      double theMax = 0.;
      int lIpo =  7;
      int nIpo = 160;
//@@---------------
      int nIpo6 = nIpo + lIpo;
      float dRefI = (CFRefIMax - CFRefIMin)/float( nIpo );
      float refIMn6 = CFRefIMin;
      float refIMx6 = CFRefIMax + float( lIpo ) * dRefI;
      //std::cout << dRefI << "  " << refIMx6 << std::endl;

      float momPart = pPartPhot_->pPart()->mom();
      for( int kPaTy=kTyLL; kPaTy<=kTyMM; kPaTy+=3 ) {
//    ----------------------------------------------
	float massIpo = CsRCRecConst::Instance()->massPartv()[kPaTy];
        double betaIpo = momPart /
          sqrt( momPart*momPart + massIpo*massIpo );

//---------------------------------------------------
        CFRefIMin = 1./ (betaIpo*1.00000);
	CFRefIMax = 1./ (betaIpo*0.99849);             // cos( 55 mard )
	dRefI = (CFRefIMax - CFRefIMin)/float( nIpo );
	refIMn6 = CFRefIMin;
	refIMx6 = CFRefIMax + float( lIpo ) * dRefI;
//---------------------------------------------------

	if( kHH < kHHMx ) {
	  int kHist = kOffHi + 3160 + kHH;
	  stringstream hN3160;
	  hN3160 << kHist;
	  string hTitle = "likelihood vs index";
	  CsHistograms::SetCurrentPath("/RICH");
          vRC3160.push_back( new CsHist1D( hN3160.str(), hTitle,
//hh      ------------------------------------------------------
					   nIpo6, refIMn6, refIMx6 ) );

	  CsHistograms::SetCurrentPath("/");
	}

        double pLike0 = -1.;
	double theIpo0 = -1.;
	double pLike[nIpo];
	for( int kIpo=0; kIpo<nIpo; kIpo++ ) {
	  float refInd = CFRefIMin + kIpo * dRefI + dRefI/2.;
	  double cosThew = 1./ (betaIpo*refInd);
	  theIpo = 0.;
	  if( cosThew <= 1.) theIpo = acos( cosThew ) * 1000.;
	  pLike[kIpo] = getLikelihood( theReco, theIpo, nPhot );
//                      ---------------------------------------
	  if( pLike[kIpo] > pLike0 ) {
	    pLike0 = pLike[kIpo];
	    theIpo0 = theIpo;
	  }

	  if( kHH < kHHMx ) {
	    xh = refInd;
	    wh = pLike[kIpo];
	    if( vRC3160[kHH] ) vRC3160[kHH]->Fill( xh, wh );
//hh                           ----------------------------
	  }
	}

	int iIpo = 0;
	if( kHH < kHHMx ) {
          iIpo = 1;
	  xh = CFRefIMin + (nIpo+iIpo) * dRefI - dRefI/2.;
	  wh = momPart;
	  if( vRC3160[kHH] ) vRC3160[kHH]->Fill( xh, wh );
//hh                         ----------------------------
	  iIpo = 2;
	  xh = CFRefIMin + (nIpo+iIpo) * dRefI - dRefI/2.;
	  wh = CFRefInd;
	  if( vRC3160[kHH] ) vRC3160[kHH]->Fill( xh, wh );
//hh                         ---------------------------
	  iIpo = 3;
	  //double cosThew = 1./ (betaIpo*CFRefInd);
	  for( int kPa=2; kPa<=14; kPa+=3 ) {
  	    float massPa = CsRCRecConst::Instance()->massPartv()[kPa];
            double betaPa = momPart /
              sqrt( momPart*momPart + massPa*massPa );
	    double cosThew = 1./ (betaPa*CFRefInd);
	    double refInd = 1./ (cosThew*betaIpo);
	    xh = CFRefIMin + (nIpo+iIpo) * dRefI - dRefI/2.;
	    wh = refInd;
	    if( vRC3160[kHH] ) vRC3160[kHH]->Fill( xh, wh );
//hh                           ----------------------------
	    iIpo++;
	  }
	}
	kHH++;
      }

    }

  }


//===========================================================================
  double CsRCRing::getQsquare( const double theIpo, int &nPhoRing ) {
//-------------------------------------------------------------------


//- compute  qSquare value for a given hypothesis
//  ---------------------------------------------
//- Paolo  -  October 2002
//    rev.    December 2002


    CsRCExeKeys *key = CsRCExeKeys::Instance();

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      key->acknoMethod( "CsRCRing::getQsquare" );
    }

    double theReco = the_;

    list<CsRCPhoton*> lPhotons = lPhotons_;
    list<CsRCPhoton*>::const_iterator ih;

    double qSquare = 0.;
    int kPhot = 0;
    for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
      double thePhot = (*ih)->the();
      if( thePhot >= theLoL_  &&  thePhot <= theUpL_ ) {

	double sigPhot = (*ih)->sigmaPhoPid( pPartPhot_ );
//                       --------------------------------
	double qSw = (thePhot - theIpo) / sigPhot;
	qSquare += qSw*qSw;

	kPhot++;
      }
    }   /* end for on Photons */
    //cout << kPhot << "  " << qSquare << endl;

    nPhoRing = kPhot;
    if( nPhoRing > 0) {
      qSquare /= float( nPhoRing );
      //cout << theIpo << "  " << nPhoRing << "  " << qSquare << endl;
    }

    if( nPhoRing == 0 ) {
      //cout << theReco << "  " << theLoL << "  " << theUpL << endl;
      //cout << the_ << endl;
      //for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) 
      //cout << (*ih)->the() << "  ";
      //cout << endl;
    }

    return  qSquare;

  }


//===========================================================================
  double CsRCRing::getRingQsq( const double theReco, int &nPhoRing ) {
//--------------------------------------------------------------------


//- compute reconstructed ring chi-square
//  -------------------------------------
//- Paolo  -  July 2001
//    rev.    December 2002 


    list<CsRCPhoton*> lPhotons = lPhotons_;
    list<CsRCPhoton*>::const_iterator ih;

    int kPhot = 0;
    double coSig = 0.;
    double coBkg = 0.;
    for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
      double thePhot = (*ih)->the();
      if( thePhot >= theLoL_  &&  thePhot <= theUpL_ ) {
	double sigPhot = (*ih)->sigmaPhoPid( pPartPhot_ );
//                       --------------------------------
	coSig += thePhot/(sigPhot*sigPhot);
	coBkg += 1. /(sigPhot*sigPhot);
	kPhot++;
      }
    }   // end for on Photons
    nPhoRing = kPhot;
    double theRecAv = 0.;
    if( nPhoRing > 0 ) theRecAv = coSig / coBkg;
//                     ------------------------


    double qSquare = 0.;
    for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
      double thePhot = (*ih)->the();
      if( thePhot >= theLoL_  &&  thePhot <= theUpL_ ) {
	double sigPhot = (*ih)->sigmaPhoPid( pPartPhot_ );
//@@                     --------------------------------
	double qSw = (thePhot - theRecAv) / sigPhot;
	qSquare += qSw*qSw;
      }
    }   // end for on Photons

    if( nPhoRing > 0 ) qSquare /= float( nPhoRing );

    return qSquare;

  }


//===========================================================================
  bool CsRCRing::getDetRFit( const vector<CsHist2D*>& vHist,
//----------------------------------------------------------
			     CsRCCircleFit& oCircle ) {

//- Paolo  -  November  2002

//- fit to the ring radius on detectors :
//  -------------------------------------

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;
      CsRCExeKeys::Instance()->acknoMethod( "CsRCRing::getDetRFit" );
    }

    CsRCHistos& hist = CsRCHistos::Ref();
    double xh, yh, wh;

    CsRCRecConst* cons = CsRCRecConst::Instance();
    double radDeg = cons->RadDeg();
    static float CFRefInd = cons->CFRefInd();
    //int nPhotMinRing = cons->nPhotMinRing();
    int nPhotMinRing = int( cons->nPhoAliMir() );        //   040331

    int kDetPart = pPartPhot_->kDetPart();
    double xPade = pPartPhot_->vPoPaDetW()[kDetPart].x();
    double yPade = pPartPhot_->vPoPaDetW()[kDetPart].y();

    list<CsRCPhoton*> lPhotons = lPhotons_;
    list<CsRCPhoton*>::iterator ih;

    int nPhotons = lPhotons.size();
    if( nPhotons == 0 ) return false;

    double xClu[nPhotons];
    double yClu[nPhotons];
    double errClu[nPhotons];

    double theReco = the_;
    double momPart = pPartPhot_->pPart()->mom();
    double sigPhoRec = pPartPhot_->sigmaPhoRec( momPart );
//                                 ----------------------
    //float sigC = cons->sigCut();
    float sigC = cons->nSigAliMir();                     //   040331
//@@-------------------------------
    float dSigC = sigC * sigPhoRec;
    double theLoL = theReco - dSigC;
    double theUpL = theReco + dSigC;

    double sigma = 5.;
//@@-----------------
    int kPhot = 0;
    for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {

      if( (*ih)->kDetClu() != kDetPart ) continue;       //   040329
      double theW = (*ih)->the();
      if( theW >= theLoL  &&  theW <= theUpL ) {

	xClu[kPhot] = (*ih)->ptToClu()->xc() - xPade;
	yClu[kPhot] = (*ih)->ptToClu()->yc() - yPade;
	errClu[kPhot] = sigma*sigma;

	kPhot++;
      }
    }   /* end for on Photons */
    int nPoint = kPhot;

//- fit :
//  -----
    bool exeFit = true;
    if( nPoint < nPhotMinRing ) exeFit = false;
//@@------------------------------------------
    //if( abs( kRingLoc_ ) != 1 ) exeFit = false;         //   021212//040329
//@@------------------------------------------
    if( exeFit ) {

      static bool circleFitDef = true;
      if( circleFitDef ) {
	circleFitDef = false;
	CsRCExeKeys::Instance()->acknoMethod( "CsRCRing::getDetRFit(Fit)" );
      }

      int nParam = 3;
      double param[nParam];
      //^param[0] = the_ * 3.36;
      param[0] = the_ * 3.20;
      param[1] = 0.;
      param[2] = 0.;
      int iPaFit[nParam];
      iPaFit[0] = 0;
      iPaFit[1] = 0;
      iPaFit[2] = 0;

      CsRCCircleFit oCirclev( nPoint, xClu, yClu, errClu,
//    ---------------------------------------------------
			      nParam, param, iPaFit );

      if( oCirclev.doChiFit() );
//        -------------------
      oCircle = oCirclev;
      if( !oCirclev.flag() ) {
	oCircle.setFlag( false );
	return  false;
      }

      //oCirclev.print();
      oCirclev.doHist( vHist );
//    ------------------------

      if( yPade > 0.) {
	xh = oCirclev.para()[1];
	yh = oCirclev.para()[2];
	if( vHist[5] ) vHist[5]->Fill( xh, yh );
//                     ------------------------
      }
      if( yPade < 0.) {
	xh = oCirclev.para()[1];
	yh = oCirclev.para()[2];
	if( vHist[6] ) vHist[6]->Fill( xh, yh );
//                     ------------------------
      }
    }

    if( !exeFit ) {
      oCircle.setFlag( false );
      return false;
    }

    return  true;

  }


//===========================================================================
  bool CsRCRing::getDetEFit( const vector<CsHist2D*>& vHist,
//----------------------------------------------------------
			     CsRCEllipseFit& oEllipse ) {

//- Paolo  -  February  2003


//- fit of an ellipse to the ring on detectors :
//  --------------------------------------------

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;
      CsRCExeKeys::Instance()->acknoMethod( "CsRCRing::getDetEFit" );
    }

    CsRCHistos& hist = CsRCHistos::Ref();
    double xh, yh, wh;

    CsRCRecConst* cons = CsRCRecConst::Instance();
    double radDeg = cons->RadDeg();
    static float CFRefInd = cons->CFRefInd();
    int nPhotMinRing = cons->nPhotMinRing();


    int kDetPart = pPartPhot_->kDetPart();
    double xPade = pPartPhot_->vPoPaDetW()[kDetPart].x();
    double yPade = pPartPhot_->vPoPaDetW()[kDetPart].y();

    list<CsRCPhoton*> lPhotons = lPhotons_;
    list<CsRCPhoton*>::iterator ih;

    int nPhotons = lPhotons.size();
    if( nPhotons == 0 ) return false;
//  --------------------------------

    double xClu[nPhotons];
    double yClu[nPhotons];
    double errClu[nPhotons];

    double theReco = the_;
    double momPart = pPartPhot_->pPart()->mom();
    double sigPhoRec = pPartPhot_->sigmaPhoRec( momPart );
//                                 ----------------------
    float sigC = cons->sigCut();
    float dSigC = sigC * sigPhoRec;
    double theLoL = theReco - dSigC;
    double theUpL = theReco + dSigC;

    double sigma = 5.;
    int kPhot = 0;
    for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {

      double theW = (*ih)->the();
      if( theW >= theLoL  &&  theW <= theUpL ) {

	xClu[kPhot] = (*ih)->ptToClu()->xc() - xPade;
	yClu[kPhot] = (*ih)->ptToClu()->yc() - yPade;
	errClu[kPhot] = sigma*sigma;

	kPhot++;
      }
    }   /* end for on Photons */
    int nPoint = kPhot;

//- fit :
//  -----
    bool exeFitE = true;
    if( nPoint < nPhotMinRing ) exeFitE = false;
    if( nPoint < 7 ) exeFitE = false;             //^   no crash(?)
//@@-------------------------------------------
    if( abs( kRingLoc_ ) != 1 ) exeFitE = false;
//@@-------------------------------------------
    if( exeFitE ) {

      static int kOffHi = cons->kOffHi();
      static vector<CsHist2D*> vRC3740;
      static bool ellipseFitDef = true;
      if( ellipseFitDef ) {
	ellipseFitDef = false;
	CsRCExeKeys::Instance()->acknoMethod( "CsRCEllipseFit from CsRCRing");

	for( int j=0; j<10; j++) vRC3740.push_back( NULL );
	int level = 1;
        bool bk = CsRCHistos::Ref().bookHis()  &&
	  level <= CsRCHistos::Ref().levelBk();
        if( bk ) {
	  string hTitle;
	  vRC3740.clear();

          stringstream hN3741;
          hN3741 << kOffHi + 3741;
          hTitle = "c0fit - c0in";
          vRC3740.push_back( new CsHist2D( hN3741.str(), hTitle,
					   100, -10., 10., 100, -10., 10. ) );
          stringstream hN3742; 
          hN3742 << kOffHi + 3742;
          hTitle = "AAfit - AAin vs ddet";
          vRC3740.push_back( new CsHist2D( hN3742.str(), hTitle,
					   100, -10., 10., 100, 0., 1000. ) );
          stringstream hN3743;
	  hN3743 << kOffHi + 3743;
	  hTitle = "epsi vs ddet";
          vRC3740.push_back( new CsHist2D( hN3743.str(), hTitle,
					   80, 0.8, 1.2, 100, 0., 1000. ) );
	  stringstream hN3744;
	  hN3744 << kOffHi + 3744;
	  hTitle = "chisq / nu";
	  vRC3740.push_back( new CsHist2D( hN3744.str(), hTitle,
					   100, 0., 10., 100, 0., 100. ) );
	  stringstream hN3745;
	  hN3745 << kOffHi + 3745;
	  hTitle = "n iter";
	  vRC3740.push_back( new CsHist2D( hN3745.str(), hTitle,
					   20, 0., 20., 100, 0., 100. ) );
	  stringstream hN3746;
	  hN3746 << kOffHi + 3746;
	  hTitle = "fit fails";
          vRC3740.push_back( new CsHist2D( hN3746.str(), hTitle,
					   10, 0., 10., 20, 0., 20. ) );
	  stringstream hN3747;
	  hN3747 << kOffHi + 3747;
	  hTitle = "pulls";
          vRC3740.push_back( new CsHist2D( hN3747.str(), hTitle,
					   100, -5., 5., 100, 0., 100. ) );
	  stringstream hN3748;
	  hN3748 << kOffHi + 3748;
	  hTitle = "phi vs phiPa";
          vRC3740.push_back( new CsHist2D( hN3748.str(), hTitle,
					   90, -45., 45., 90, 0., 360. ) );
	}
      }
      int nParamE = 7;
      double paramE[nParamE];
      paramE[0] = 0.;
      paramE[1] = 0.;
      paramE[2] = the_ * 3.36;
      paramE[3] = 1.0;
      if( xPade != 0. ) {
	paramE[4] = atan( yPade/xPade );
      } else {
	paramE[4] = 1.5708;
      }
      //paramE[4] = 0.;
      paramE[5] = xPade;
      paramE[6] = yPade;
      int iPaFitE[nParamE];
      iPaFitE[0] = 0;
      iPaFitE[1] = 0;
      iPaFitE[2] = 0;
      iPaFitE[3] = 0;
      iPaFitE[4] = 0;
      iPaFitE[5] = 1;
      iPaFitE[6] = 1;

      CsRCEllipseFit oEllipsev( nPoint, xClu, yClu, errClu, errClu,
//    -------------------------------------------------------------
				nParamE, paramE, iPaFitE );

      if( oEllipsev.doChiFit() );
      oEllipse = oEllipsev;
      if( !oEllipsev.flag() ) {
	oEllipse.setFlag( false );
	return  false;
      }
      //oEllipsev.print();
      oEllipsev.doHist( vRC3740 );
//    ---------------------------

    }

    if( !exeFitE ) {
      oEllipse.setFlag( false );
      return false;
    }

    return  true;

  }


//===========================================================================
  long double CsRCRing::ringPoissProb( const double theIpo ) {
//------------------------------------------------------------


//- Paolo  -  August 2001


    long double poissPro = 0.;

    double nPhoExp = getNPhotExpected( theIpo );
//                   --------------------------
    //float corrBack = 1.;
    //if( flagBack_ ) corrBack = ringBack_;
    //int nPhoRingCo = int( lPhotons_.size() * corrBack );
    double nPhoRingCo = getNPhotCorrected( theIpo );
//                      ---------------------------
    if( nPhoRingCo == 0. ) return  poissPro;

    //long double poissPro = 1.;
    for( int kPhot=1; kPhot<=nPhoRingCo; kPhot++ )
      poissPro *= nPhoExp / kPhot;
    poissPro *= exp( -nPhoExp );
    //cout << nPhoRingCo << "  " << nPhoExp << "  " << poissPro << endl;

    return  poissPro;

  }


//===========================================================================
  double CsRCRing::ringGaussProb( const double theIpo ) {
//-------------------------------------------------------


//- Paolo  -  February 2003


    double gaussPro = 0.;

    double nPhoExp = getNPhotExpected( theIpo );
//                   --------------------------
    double nPhoRingCo = getNPhotCorrected( theIpo );
//                      ---------------------------
    if( nPhoRingCo == 0. ) return  gaussPro;

    double sigma = 5.0;
//@@------------------
    double expw = (nPhoRingCo - nPhoExp) / sigma;
    //double gaussPro = exp( -0.5* expw*expw ) / (sigma*2.51);
    gaussPro = 0.5* ( 1.+ erf_( expw/1.41421 ) );

    //cout << nPhoRingCo << "  " << nPhoExp << "  " << gaussPro << endl;

    return  gaussPro;

  }


//===========================================================================
  double CsRCRing::getNPhotExpected( const double theIpo ) {
//----------------------------------------------------------


//- Paolo  -  August 2001
//    rev.    December 2002


    float nZero = CsRCRecConst::Instance()->nZero();
    double pathLen = pPartPhot_->pPart()->pathLen();
    double nPhoExp = nZero * pathLen/10.* (theIpo/1000.)*(theIpo/1000.); 

    //cout << theIpo << "  " << nZero << "  " << pathLen
    //     << "  " << nPhoExp << endl;

    return  nPhoExp;

  }


//===========================================================================
  double CsRCRing::getNPhotCorrected( const double theIpo ) {
//-----------------------------------------------------------


//- Paolo  -  March 2003


    float corrBack = 0.;
    if( flagBack_ ) corrBack = ringBack_;                  //   !!!

    double nPhoRingCo = lPhotons_.size() * corrBack;
    //cout << lPhotons_.size() << "  " << ringBack_ << endl;

    return  nPhoRingCo;

  }


//===========================================================================
  bool CsRCRing::getRingness( const list<CsRCPhoton*>& lPhotons,
//--------------------------------------------------------------
			      double& ringNess ) {


//- Paolo  -  November  2002


    //list<CsRCPhoton*> lPhotons = lPhotons_;
    list<CsRCPhoton*>::const_iterator ih;

//- order phi of photons in ascending order :
//  -----------------------------------------
    int nPhotons = lPhotons.size();
    if( nPhotons == 0 ) return false;

    double phiw[nPhotons];
    int kPhot = 0;
    for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
      phiw[kPhot] = (*ih)->phi();
      kPhot++;
    }
    //for( int k=0; k<nPhotons; k++ ) cout << phiw[k] << "  ";
    //cout << endl;
    double phio[nPhotons];
    for( int kPhot=0; kPhot<nPhotons; kPhot++ ) {
      double phiMin = 10000;
      int kMin = 0;
      for( int jPhot=0; jPhot<nPhotons; jPhot++ ) {
	double phi = phiw[jPhot];
	if( phi < phiMin  &&  phi < 10000 ) {
	  phiMin = phi;
	  kMin = jPhot;
	}
      }
      phio[kPhot] = phiMin;
      phiw[kMin] = 10000;
    }
    //for( int k=0; k<nPhotons; k++ ) cout << phio[k] << "  ";
    //cout << endl;

//- compute ringness :
//  ------------------
    double dPhi = 360. / nPhotons;
    double term = ( fabs( 360.+phio[0]-phio[nPhotons-1] ) - dPhi ) / dPhi;
    ringNess = term*term;
    for( int kPhot=1; kPhot<nPhotons; kPhot++ ) {
      term = ( fabs( phio[kPhot]-phio[kPhot-1] ) - dPhi ) / dPhi;
      ringNess += term*term;
    }
    ringNess /= nPhotons;
    //cout << "ringNess = " << ringNess << endl;

    return true;
 
  }


//===========================================================================
  float CsRCRing::recoCut( const double beta ) {
//----------------------------------------------

//--- interface function :
//    --------------------
//    December 1999,   rev. November 2002


      //float recoCut = 3. * sigmaRing( beta );
      float recoCut = 2. * sigmaRing( beta );
//@@----------------------------------------

      static bool firstCall = true;
      if( firstCall ) {
	firstCall = false;

        if( CsRCRecConst::Instance()->printConsts() ) {
	  //cout << endl;
	  cout << " RICHONE, CsRCRing::recoCut :   " 
	       << recoCut << " (1rst call)" << endl;
	}
      }

      return recoCut;

  }


//===========================================================================
  float CsRCRing::sigmaRing( const double beta ) {
//------------------------------------------------

//--- interface function :
//    --------------------
//    March  2001, rev, December 2002

      CsRCHistos& hist = CsRCHistos::Ref();
      double xh, yh, wh;

      float sigmaRing = 1000000.;

      float sigma = sigmaRing4( beta );
//@@----------------------------------

      if( sigma > 0. ) sigmaRing = sigma;

//--- monitor error values
      xh = beta - 1.;;
      yh = sigmaRing;
      if( hist.hRC1705 ) hist.hRC1705->Fill( xh, yh );
//hh                     ----------------------------

      return sigmaRing;

  }


//===========================================================================
  float CsRCRing::sigmaRing3( const double beta ) {
//-------------------------------------------------


//    December  1999,  rev. April  2001


//--- sigma ring from MC :
//    --------------------

      static float par[10];
      static float a1;
      static float a9;
      static float x9;

      static bool firstCall = true;
      if( firstCall ) {
	firstCall = false;

        CsOpt* opt = CsOpt::Instance();
        CsRCRecConst *cons = CsRCRecConst::Instance();

        a1 = 0.178;
        a9 = 1.70;
        x9 = 0.9984;

//----- defaults
//      ( from Dima's 48k file ) :
//      --------------------------
	par[0] =   0.551;                       //   10/11/99

//      from Vadim's 6*3k ref files/lepto_full ( co-18k-15t.hist )
//      ----------------------------------------------------------
	//par[0] =   15.62;                      //(?)   02/04/01

//----- from rich1.options :
//      --------------------
	bool boo;
        float pPar;
        boo = opt->CsOpt::getOpt( "RICHONE", "sigmaRing", pPar );
        if( boo ) {
	  par[0] = pPar;
        }
        if( cons->printConsts() ) {
	  //cout << endl;
	  cout << " RICHONE, CsRCRing::sigmaRing3 :   " 
	       << par[0] << endl;
	}

      }

//--- from sigring.kumac/sigmarin3.f - histo.s : 3010 (+off)
//    ------------------------------------------------------
//    x = beta
//    function : [ a9 * (1-x9)**p0 ] / [ (a9/a1)**1/p0 * (x-x9) - (x-1) ]**p0
//    f = a1  per x = 1.
//    f = a9  per x = 0.9984 = x9
//    ---------------------------

      float betaw = beta;
      if( betaw < 0.9984 ) betaw = 0.9984;

      double sigmaRing3 = a9 * pow( double(1.-x9), double(par[0]) ) / 
          pow( ( pow( double(a9/a1), double(1./par[0]) ) * (betaw-x9) - (betaw-1.) ),
          double(par[0]) );

      //cout << "sigmaRing3  " << sigmaRing3 << endl;

      return sigmaRing3;

  }


//========================================================================
  float CsRCRing::sigmaRing4( const double beta ) {
//-------------------------------------------------


//- Paolo   -   December 2202


//- sigma-ring from Data :
//  ----------------------

//- from sigring-data.kumac & histo 3053 (+off)

    static float par[10];

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

//--- defaults
//    ( from run 20330 & 22642, 2002 )
//    --------------------------------
      par[0] =     0.47;    //   02/12/10
      par[1] = - 855.;      //   02/12/10

//    ( from run 51908, 2006 )
//    --------------------------------
      par[0] =     0.19;    //   06/11/28
      par[1] = - 836.;      //   06/11/28

//--- from rich1.options :
//    --------------------
      CsOpt* opt = CsOpt::Instance();
      vector<float> vPar;
      bool boo = opt->CsOpt::getOpt( "RICHONE", "sigmaRing", vPar );
      if( boo ) {
	for( unsigned int k=0; k<vPar.size(); k++ ) par[k] = vPar[k];

	if( CsRCRecConst::Instance()->printConsts() ) {
	  //cout << endl;
	  cout << " RICHONE, CsRCRing::sigmaRing4 :   ";
	  for( unsigned int k=0; k<vPar.size(); k++ ) cout << par[k] << "  ";
	  cout << endl;
	}
      } else {
	//cout << " RICHONE, CsRCRing::sigmaRing4 :   ";
	//cout << "NO parameters read, default used!" << endl;
        string mess = "RICHONE, CsRCRing::sigmaRing4 : ";
	string err = "NO parameters read, default used!";
	mess.append( err );
	CsErrLog::Instance()->mes( elError, mess );
      }
    }

    float betaw = beta;
    if( betaw < 0.9984 ) betaw = 0.9984;

    float sigmaRing = par[0] + par[1]* (betaw - 1.);
//@@-----------------------------------------------
    //cout << "      sigmaRing4 " << sigmaRing << endl;

    return sigmaRing;

  }


//========================================================================
  int CsRCRing::getVarMcan( const CsRCPartPhotons* papho,
//-------------------------------------------------------
			    const int mCan ) {


//- Paolo   -   February 2004

    int mCanVar = 0;

    CsRCPartPhotons* pt = const_cast<CsRCPartPhotons*>(papho);
    float fact = 1.;
    CsRCParticle* part = pt->pPart();
    float ddpdet = part->ddPaDet()[part->kDetPart()];
    float norm = pt->sigmaPhot8( 0.);
    if( norm > 0. ) 
      mCanVar = int( mCan* fact* pt->sigmaPhot8( ddpdet ) / norm );
    if( mCanVar < mCan ) mCanVar = mCan;

    //if( mCanVar > mCan )
    //  std::cout << "getVarMcan  " << mCanVar << std::endl;

    mCanVar = mCan;
//@@--------------

    return mCanVar;

  }


//========================================================================
  int CsRCRing::getVarNMore( const CsRCPartPhotons* papho, 
//--------------------------------------------------------
			     const int nMore ) {


//- Paolo   -   February 2004
 
    int nMoreVar = 0;

    CsRCPartPhotons* pt = const_cast<CsRCPartPhotons*>(papho);
    float fact = 1.;
    CsRCParticle* part = pt->pPart();
    float ddpdet = part->ddPaDet()[part->kDetPart()];
    float norm = pt->sigmaPhot8( 0.);

    if( norm > 0. ) 
      nMoreVar = int( nMore* fact* pt->sigmaPhot8( ddpdet ) / norm );
    if( nMoreVar < nMore )  nMoreVar = nMore;

    //if( nMoreVar > nMore )
    //  std::cout << "getVarNMore  " << nMoreVar << std::endl;

    nMoreVar = nMore;
//@@----------------

    return nMoreVar;
  
  }


//========================================================================
  void CsRCRing::checkChiReso() {
//-------------------------------

//- Paolo  -  April 2006


    CsRCExeKeys* key = CsRCExeKeys::Instance();
    CsRCRecConst *cons = CsRCRecConst::Instance();
    int kOffHi = cons->kOffHi();

    static std::vector<CsHist2D*> vRC3280;
    static std::vector<CsHist1D*> vRC8070;
    static int nPhoRing = 0;
    static int nPhoMin = 20;
    static int nFGood = 0;
    static int nHGood = 0;
    static int nFFail = 0;
//@@-----------------------

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      key->acknoMethod( "CsRCRing::checkChiReso" );

      for( int kh=0; kh<7; kh++ ) vRC3280.push_back( NULL );
      if( CsRCHistos::Ref().bookHis() &&
	  CsRCHistos::Ref().levelBk() >= 2 ) {
	CsHistograms::SetCurrentPath("/RICH");
	vRC3280.clear();
        string hTitle = "chisq reso";
        int kHist = 0;
	stringstream hN3280;
        kHist = kOffHi + 3280;
        hN3280 << kHist;
        vRC3280.push_back( new CsHist2D( hN3280.str(), hTitle,
					 100, 0., 50., 100, 0., 5.) );
	stringstream hN3281;
        kHist = kOffHi + 3281;
        hN3281 << kHist;
        vRC3280.push_back( new CsHist2D( hN3281.str(), hTitle,
					 100, 0., 50., 100, 0., 3.) );
	stringstream hN3282;
        kHist = kOffHi + 3282;
        hN3282 << kHist;
        vRC3280.push_back( new CsHist2D( hN3282.str(), hTitle,
					 100, 0., 50., 100, -5., 5.) );
	stringstream hN3283;
        kHist = kOffHi + 3283;
        hN3283 << kHist;
        vRC3280.push_back( new CsHist2D( hN3283.str(), hTitle,
					 100, 0., 50., 100, 0., 3.) );
	stringstream hN3284;
        kHist = kOffHi + 3284;
        hN3284 << kHist;
        vRC3280.push_back( new CsHist2D( hN3284.str(), hTitle,
					 100, 0., 50., 100, -5., 5.) );
	stringstream hN3285;
        kHist = kOffHi + 3285;
        hN3285 << kHist;
        vRC3280.push_back( new CsHist2D( hN3285.str(), hTitle,
					 100, 0., 50., 100, -5., 5.) );
	stringstream hN3286;
        kHist = kOffHi + 3286;
        hN3286 << kHist;
        vRC3280.push_back( new CsHist2D( hN3286.str(), hTitle,
					 100, 0., 50., 100, -1., 1.) );
        CsHistograms::SetCurrentPath("/");
      }
    }

    double xh, yh, wh;
    //std::cout << cons->CFRefInd() << std::endl;

    nPhoRing = lPhotons_.size();
    int iPartTy = 0;
    if( key->MCarloEvent() ) {
      iPartTy = pPartPhot_->pPart()->iPartT();
      if( iPartTy > 200 ) iPartTy -= 200;
      if( iPartTy > 30 ) iPartTy = 0;   //   wrong myFile?
//@@--------------------------------
    }
    //if( nPhoRing < nPhoMin ) return;
//@@-------------------------------

    float momPart = pPartPhot_->pPart()->mom();

    if( momPart < 2.5 ) return;
//@@--------------------------

    double phiMass  = CsRichOne::Instance()->phiMass();
    //std::cout << setprecision( 4 );
    //std::cout << "=======   " << phiMass << std::endl;
    //std::cout << setprecision( 2 );

//- PHI MASS SELECTION! (PHI Sample)
//---------------------
    if( key->myFileType() == 3 ) {
      if( fabs( phiMass - 1.020) > 0.012 ) return;
//@@---------------------------------------------
    }

    double thePart = acos( pPartPhot_->pPart()->nd() );
    //if( thePart < 50./1000.) return;
//@@-------------------------------

    list<CsRCPhoton*> lPhotons = lPhotons_;
    list<CsRCPhoton*>::const_iterator ih;
    double theWgAv = 0.;
    double sigWgAv = 0.;
    for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
      double thePhot = (*ih)->the();
      if( thePhot >= theLoL_  &&  thePhot <= theUpL_ ) {
	double sigPhot = (*ih)->sigmaPhoPid( pPartPhot_ );
//                       --------------------------------
	double sig2 = sigPhot*sigPhot;
	theWgAv += thePhot / sig2;
	sigWgAv += 1./ sig2;
      }
    }
    theWgAv /= sigWgAv;
    sigWgAv = 1./ sqrt( sigWgAv );
    //std::cout << "checkChiReso " << theWgAv << "  " << sigWgAv << std::endl;

    int nPhot = 0;
    double theMin =  0.;
    double theMax = 80.;
    int nIpo = 321;
    double dThe = (theMax - theMin)/float( nIpo );

    double qSquaMn = 1000000.;
    double theIpoMn = -1.;
    double qSquak[nIpo];
    double thek[nIpo];
    int kIpoMn = -1;
    for( int kIpo=0; kIpo<nIpo; kIpo++ ) {
      thek[kIpo] = theMin + kIpo * dThe + dThe/2.;
      qSquak[kIpo] = getQsquare( thek[kIpo], nPhot );
//                   -------------------------------
      if( qSquak[kIpo] < qSquaMn ) {
        qSquaMn = qSquak[kIpo];
        theIpoMn = thek[kIpo];
	kIpoMn = kIpo;
      }
      //std::cout << setprecision( 7 );
      //std::cout << kIpo << "  " << thek[kIpo] << "  " << qSquak[kIpo] 
      //<< std::endl;
    }
    //std::cout << kIpoMn << "  " << theIpoMn << "  " << qSquaMn << std::endl;

    if( theIpoMn <  5. ) return;
    if( theIpoMn > 60. ) return;
//@@---------------------------

    int kTyLL = 2;
//@@-------------
    double theIpo[5];
    double qSquam[5];
    int kIpo = 0;
    for( int kPaTy=kTyLL; kPaTy<=14; kPaTy+=3 ) {
//  -------------------------------------------
      float massIpo = CsRCRecConst::Instance()->massPartv()[kPaTy];
      double betaIpo = momPart / sqrt( momPart*momPart + massIpo*massIpo );
      double cosTheW = 1./(betaIpo* CsRCRecConst::Instance()->CFRefInd());
      theIpo[kIpo] = -1.;
      if( cosTheW <= 1.) theIpo[kIpo] = acos( cosTheW ) * 1000.;
      qSquam[kIpo] = getQsquare( theIpo[kIpo], nPhot );
//                   ---------------------------------
      kIpo++;
    }
    int nMassIpo = kIpo;

    double qSquamMn = 1000000.;
    double theIpomMn = -1.;
    int kIpomMn = -1;
    for( int kIpo=0; kIpo<5; kIpo++ ) {
      if( qSquam[kIpo] < qSquamMn ) {
	qSquamMn = qSquam[kIpo];
	theIpomMn = theIpo[kIpo];
	kIpomMn = kIpo;
      }
    }

    double norm  = -1.;
    double mean  = -1.;
    double sigma = -1.;

    bool exeFit = true;
//@@------------------
    if( exeFit ) {
    if( kIpomMn == 3 ) {

    int nPf = 24;
//@@------------
    double thef[nPf];
    double qSquaf[nPf];
    double eSquaf[nPf];
    int kf = 0;
    for( int kk=0; kk<nPf; kk++ ) {
      int ku = kIpoMn - nPf/2 + kk;
      if( ku < 0 ) continue;
      if( ku >= nIpo ) continue;
      thef[kf] = thek[ku];
      qSquaf[kf] = exp( -0.5* qSquak[ku] );
      eSquaf[kf] = 0.01* qSquak[ku];
//@@-------------------------------
      //std::cout << kf << " " << thef[kf] << "  " << qSquaf[kf] << std::endl;
      kf++;
    }

    int nPoint = kf;
    int nParam = 10;
    double param[nParam];
    double toll[nParam];
    int iPaFit[nParam];
    for( int kp=0; kp<nParam; kp++ ) {
      param[kp] = 0.;
      param[kp] = 0.;
      iPaFit[kp] = 1;
    }
    //param[0] = exp( -0.5* qSquaMn )/2.;
    param[0] = 0.5;
    param[1] = theIpoMn;
    param[2] = 1.;
    toll[0] = 0.01;
    toll[1] = 0.0001;
    toll[2] = 0.0001;
    iPaFit[0] = 0;
    iPaFit[1] = 0;
    iPaFit[2] = 0;

    CsRCGauPolFit oSquaFt( nPoint, thef, qSquaf, eSquaf,
//  ----------------------------------------------------
			   nParam, param, toll, iPaFit );
    if( oSquaFt.doChiFit() );
//      ------------------

    if( oSquaFt.flag() ) {
      nFGood++;
      norm  = oSquaFt.para()[0];
      mean  = oSquaFt.para()[1];
      sigma = oSquaFt.para()[2];
    } else {
      std::cout << "  checkChiReso : fit failure!  " << nFFail << "  ";
      std::cout << oSquaFt.nIter() << "  ";
      std::cout << theIpoMn << "  " << qSquaMn << std::endl;
      nFFail++;
      if(nFFail <= 20 ) {
        //int kFF = nFFail - 1;
        //int kHist = kOffHi + 8070 + kFF;
        //stringstream hN8070;
        //hN8070 << kHist;
        //string hTitle = "qSqua shape FFail";
        //CsHistograms::SetCurrentPath("/RICH");
        //vRC8070.push_back( new CsHist1D( hN8070.str(), hTitle,
//hh    ------------------------------------------------------
	//  				 nIpo, theMin, theMax ) );
        //CsHistograms::SetCurrentPath("/");
        //for( int kIpo=0; kIpo<nIpo; kIpo++ ) {
	//xh = thek[kIpo];
	//wh = qSquak[kIpo];
	//if( vRC8070[kFF] ) vRC8070[kFF]->Fill( xh, wh );
//hh                         ----------------------------
	//}
      }
      //return;
    }
    //std::cout << "Fit  " << norm << "  " << mean << "  " 
    //          << sigma << std::endl;
    //std::cout << nFGood << "  " << nFFail << std::endl;
    }
    }


    int kHH = 0;
    xh = momPart;
    yh = sigma;
    if( vRC3280[kHH] ) vRC3280[kHH]->Fill( xh, yh );
//hh                   ----------------------------

    /*int kTyLL = 2;
//@@-------------
    double theIpo[5];
    double qSquam[5];
    int kIpo = 0;
    for( int kPaTy=kTyLL; kPaTy<=14; kPaTy+=3 ) {
//  -------------------------------------------
      float massIpo = CsRCRecConst::Instance()->massPartv()[kPaTy];
      double betaIpo = momPart / sqrt( momPart*momPart + massIpo*massIpo );
      double cosTheW = 1./(betaIpo* CsRCRecConst::Instance()->CFRefInd());
      theIpo[kIpo] = -1.;
      if( cosTheW <= 1.) theIpo[kIpo] = acos( cosTheW ) * 1000.;
      qSquam[kIpo] = getQsquare( theIpo[kIpo], nPhot );
//                   ---------------------------------
      kIpo++;
    }
    int nMassIpo = kIpo;

    double qSquamMn = 1000000.;
    double theIpomMn = -1.;
    int kIpomMn = -1;
    for( int kIpo=0; kIpo<5; kIpo++ ) {
      if( qSquam[kIpo] < qSquamMn ) {
	qSquamMn = qSquam[kIpo];
	theIpomMn = theIpo[kIpo];
	kIpomMn = kIpo;
      }
      }*/

    //double sigmaChiMax = 0.64;
//@@-------------------------
    //double sigmaDTheta = 84.*exp( -0.41*momPart ) + 0.65;
//@@----------------------------------------------------
      if( kIpomMn == 3 ) {
        kHH = 1;
        xh = momPart;
        yh = sigma / fabs(theIpo[kIpomMn]-theIpo[kIpomMn-1]);
        if( vRC3280[kHH] ) vRC3280[kHH]->Fill( xh, yh );
//hh                       ----------------------------
        kHH = 2;
        xh = momPart;
        yh = theIpoMn - theIpo[kIpomMn];
        if( vRC3280[kHH] ) vRC3280[kHH]->Fill( xh, yh );
//hh                       ----------------------------
	kHH = 5;
        xh = momPart;
        //yh = sigmaChiMax / fabs(theIpo[kIpomMn]-theIpo[kIpomMn-1]);
	yh = sigma - sigWgAv;
        if( vRC3280[kHH] ) vRC3280[kHH]->Fill( xh, yh );
//hh                       ----------------------------
	kHH = 6;
        xh = momPart;
        //yh = sigmaDTheta / fabs(theIpo[kIpomMn]-theIpo[kIpomMn-1]);
	yh = mean - theWgAv;
        if( vRC3280[kHH] ) vRC3280[kHH]->Fill( xh, yh );
//hh                       ----------------------------
      }
      if( kIpomMn == 2 ) {
        kHH = 3;
        xh = momPart;
        yh = sigma / fabs(theIpo[kIpomMn]-theIpo[kIpomMn-1]);
        if( vRC3280[kHH] ) vRC3280[kHH]->Fill( xh, yh );
//hh                       ----------------------------
        kHH = 4;
        xh = momPart;
        yh = theIpoMn - theIpo[kIpomMn];
        if( vRC3280[kHH] ) vRC3280[kHH]->Fill( xh, yh );
//hh                       ----------------------------
      }
      if( kIpomMn == 3 ) {
	nHGood++;
	if( nHGood < 20 ) {
	  int kFF = nHGood - 1;
	  int kHist = kOffHi + 8070 + kFF;
	  stringstream hN8070;
	  hN8070 << kHist;
	  string hTitle = "qSquare shape Good";
	  CsHistograms::SetCurrentPath("/RICH");
	  vRC8070.push_back( new CsHist1D( hN8070.str(), hTitle,
//hh      ------------------------------------------------------
	  			   nIpo, theMin, theMax ) );
	  CsHistograms::SetCurrentPath("/");
	  for( int kIpo=0; kIpo<nIpo-7; kIpo++ ) {
	    xh = thek[kIpo];
	    wh = qSquak[kIpo];
	    if( vRC8070[kFF] ) vRC8070[kFF]->Fill( xh, wh );
//hh                           ----------------------------
	  }
	  xh = thek[nIpo-1];
	  wh = momPart;
	  if( vRC8070[kFF] ) vRC8070[kFF]->Fill( xh, wh );
//hh                         ----------------------------
	  for( int kIpo=0; kIpo<nMassIpo; kIpo++ ) {
	    xh = thek[nIpo-kIpo-2];
	    wh = theIpo[kIpo];
	    if( vRC8070[kFF] ) vRC8070[kFF]->Fill( xh, wh );
//hh                           ----------------------------
	  }

	}
      }

    return;
  }


//========================================================================
  bool CsRCRing::ringPMTonly() {
//------------------------------

//- Paolo  -  October 2006

    CsRCRecConst *cons = CsRCRecConst::Instance();
    static const int kOffHi = cons->kOffHi();
    static const double TwoPi = cons->TwoPI();
    CsRCDetectors* dets = CsRCDetectors::Instance();
    int nCathode = dets->nCathode();

    CsRCParticle* part =  pPartPhot_->pPart();
    static const float focLen = part->pMirPart()->RR() / 2.;

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
      key->acknoMethod( "CsRCRing::ringPMTonly" );

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

    double rRing = the_/1000. * focLen;
    rRing *= 1.1;
    int kDetPart = part->kDetPart();
    Hep3Vector vPade = part->vPade()[kDetPart];

    int iCaPa = pPartPhot_->iCaPa();
    bool catPMT = false;
    if( iCaPa >= 0 ) catPMT = dets->ptrToCat( iCaPa )->isPMT();
    else {
      if( vPade.y() < yLimUp  &&  vPade.y() > yLimDw  &&
          vPade.x() < xLimRg  &&  vPade.x() > xLimLf ) catPMT = true;
      else  catPMT = false;
    }
    if( catPMT ) {
//--- Internal cathodes (PMT)
      double ddx = 0.;
      double ddy = 0.;
      if( vPade.y() >= 0.) ddy = yLimUp - vPade.y();
      if( vPade.y() <  0.) ddy = vPade.y() - yLimDw;
      if( vPade.x() >= 0.) ddx = xLimRg - vPade.x();
      if( vPade.x() <  0.) ddx = vPade.x() - xLimLf;
      if( rRing > 0. ) {
        if( ddx >= rRing  &&  ddy >= rRing ) return  true;
      }
    }

    return  false;
  }


//========================================================================
  bool CsRCRing::ringAPVonly() {
//------------------------------

//- Paolo  -  October 2006

    CsRCRecConst *cons = CsRCRecConst::Instance();
    static const int kOffHi = cons->kOffHi();
    static const double TwoPi = cons->TwoPI();
    CsRCDetectors* dets = CsRCDetectors::Instance();
    int nCathode = dets->nCathode();

    CsRCParticle* part =  pPartPhot_->pPart();
    static const float focLen = part->pMirPart()->RR() / 2.;

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
      key->acknoMethod( "CsRCRing::ringAPVonly" );

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

    double rRing = the_/1000. * focLen;
    rRing *= 1.1;
    int kDetPart = part->kDetPart();
    Hep3Vector vPade = part->vPade()[kDetPart];

    int iCaPa = pPartPhot_->iCaPa();
    bool catAPV = false;
    if( iCaPa >= 0 ) catAPV = dets->ptrToCat( iCaPa )->isAPV();
    else {
      if( vPade.y() < yLimUp  &&  vPade.y() > yLimDw  &&
          vPade.x() < xLimRg  &&  vPade.x() > xLimLf ) catAPV = true;
      else  catAPV = false;
    }
    if( catAPV ) {
//--- External cathodes (CsI)
      double ddx = 0.;
      double ddy = 0.;
      if( vPade.y() >= 0.) ddy = vPade.y() - yLimUp;
      if( vPade.y() <  0.) ddy = yLimDw - vPade.y();
      if( vPade.x() >= 0.) ddx = vPade.x() - xLimRg;
      if( vPade.x() <  0.) ddx = xLimLf - vPade.x();
      if( rRing > 0. ) {
        if( ddx >= rRing  ||  ddy >= rRing ) return  true;
      }
    }

    return  false;
  }


//========================================================================
  void CsRCRing::setTime() {
//--------------------------

//- Paolo   -   November 2006

//- Average Photon (=Pad) time (PMT only)
//  -------------------------------------

    CsRCExeKeys* key = CsRCExeKeys::Instance();
    CsRCRecConst* cons = CsRCRecConst::Instance();
    CsRICH1Detector* rich = CsGeom::Instance()->getRich1Detector();

    static float offset;
    static std::vector<CsHist1D*> vRC8140;
    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;
 
      key->acknoMethod( "CsRCRing::getTime" );

      offset = -1681.;
      CsOpt* opt = CsOpt::Instance();
      bool boo;
      vector<float> vflo;
      boo = opt->CsOpt::getOpt( "RICHONE", "PMTTimeCuts", vflo );
      if( boo ) offset  = vflo[0];

      if( CsRCHistos::Ref().bookHis() && CsRCHistos::Ref().levelBk() >= 2 ) {
	int kOffHi = cons->kOffHi();
        CsHistograms::SetCurrentPath("/RICH");
        int kHist = 0;
        string hTitle;
        kHist = kOffHi + 8141;
        stringstream hN8141;
	hN8141 << kHist;
        hTitle = "Ring Time";
        vRC8140.push_back( new CsHist1D( hN8141.str(), hTitle,
//hh    ------------------------------------------------------
	  	  	  	         200, -50., 50. ) );
        kHist = kOffHi + 8142;
        stringstream hN8142;
	hN8142 << kHist;
        hTitle = "Ring-Track Time";
        vRC8140.push_back( new CsHist1D( hN8142.str(), hTitle,
//hh    ------------------------------------------------------
	  	  	  	         200, -50., 50. ) );
        kHist = kOffHi + 8143;
        stringstream hN8143;
	hN8143 << kHist;
        hTitle = "Ring-Track Time";
        vRC8140.push_back( new CsHist1D( hN8143.str(), hTitle,
//hh    ------------------------------------------------------
	  	  	  	         200, -50., 50. ) );
	CsHistograms::SetCurrentPath("/");
      }
    }

    double xh, yh, wh;
    int kh;

    // TRACK TIME:- NOTA BENE:
    // i) Most tracks have no good timing. Exception are tracks w/ scifi or
    //   or hodo hits.
    // ii) Some tracks have no timing at all.
    //     As a substitute, the best possible estimate is:
    //   - Track associated to the primary vertex: event's time.
    //   - Pile-up                               : 0.
    //   Would remain to be able to distinguish pile-ups from secondaries. Which
    //   is possible if it's far halo: parallel, far away from the axis.
    //     For the time being: assign them a 0 value.
    // => Time cut:
    // i) In any case, mitigate it w/ uncertainty on track's time.
    // ii) For tracks associated to vertex, better use event's time.
    double trackTime = 0.; if (!key->readMyFile()) {
      CsTrack *track = pPartPhot_->pPart()->pTrack();
      if (track->hasMeanTime()) trackTime = track->getMeanTime();
    }

    double mTime = 0.;
    double mT0 = 0.;
    mTimeSet_ = false;
    list<CsRCPhoton*> lPhotons = lPhotons_;
    list<CsRCPhoton*>::const_iterator ih;
    int nPhot = 0;
    for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
      if( !(*ih)->flag() ) continue;
      if( !(*ih)->isPMT() ) continue;
      mTime += (*ih)->PH();
      mT0 += (*ih)->ptToClu()->lPads().front()->PMTT0();
      nPhot++;
    }
    if( nPhot > 0 ) {
      mTime /= float( nPhot );
      if( rich->T0_flag()  ||  key->readMyFile() ) offset = 0.;
      mTime -= offset;
      mT0 /= float( nPhot );
    }
    //std:: cout << mTime << "  " << trackTime << std::endl;

    kh = 0;
    if( int( vRC8140.size() ) > kh ) {
      xh = mTime;
      if( vRC8140[kh] ) vRC8140[kh]->Fill( xh );
//hh  ---------------------------------------
    }
    kh = 1;
    if( int( vRC8140.size() ) > kh ) {
      xh = mTime - trackTime;
      if( vRC8140[kh] ) vRC8140[kh]->Fill( xh );
//hh  ---------------------------------------
    }
    kh = 2;
    if( int( vRC8140.size() ) > kh ) {
      for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
	if( !(*ih)->flag() ) continue;
	if( !(*ih)->isPMT() ) continue;
	xh = (*ih)->PH() - offset - mTime;
	if( vRC8140[kh] ) vRC8140[kh]->Fill( xh );
//hh  -------------------------------------------
      }
    }


    mTime_ = mTime;
    if( !key->readMyFile() ) mT0_ = mT0;
    //std:: cout << mTime << "  " << mT0 << std::endl;
    mTimeSet_ = true;

    return;
  }
