/*!
   \file    CsRichOne.cc
   \---------------------
   \brief   CsRichOne.cc class implementation.
   \author  Paolo Schiavon
   \version 0.02B,  rev. 20/6/2000
   \date    December 1999, rev.  October 2000
*/


  #include <iostream>
  #include <ostream>

  #include <cstdio>
  #include <cstdlib>

  #include <unistd.h>
  #include <fcntl.h>
  #include <sys/types.h>
  #include <sys/stat.h>

//---------------------------
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
//----------------------------

  #include "Coral.h"
  #include "CsHistograms.h"

  #include "CsStopwatch.h"

  #include "CsMCDigit.h"
  #include "CsTrack.h"
  #include "CsMCTrack.h"

  #include "CsVertex.h"
  #include "CsParticle.h"

//-----------------------------
  #include "CsRICH1Detector.h"

  #include "CsRichOne.h"

  #include "CsRCMirrors.h"
  #include "CsRCDetectors.h"

  #include "CsRCEventParticles.h"
  #include "CsRCParticle.h"

  #include "CsRCEventPads.h"
  #include "CsRCPad.h"
  #include "CsRCEventClusters.h"
  #include "CsRCCluster.h"

  #include "CsRCEventPartPhotons.h"
  #include "CsRCPartPhotons.h"

  #include "CsRCEventRings.h"
  #include "CsRCRing.h"

  #include "CsRCEventAnalysis.h"
  #include "CsRCEventDisplay.h"
  #include "CsRichOneDisplay.h"

  #include "CsRCExeKeys.h"
  #include "CsRCRecConst.h"
  #include "CsRCHistos.h"
  #include "CsRCnTup.h"
//------------------------------

  #include <CLHEP/Vector/ThreeVector.h>
//-------------------------------------

  using namespace std;
  using namespace CLHEP;

  CsRichOne* CsRichOne::instance_ = 0;

//===========================================================================
  CsRichOne* CsRichOne::Instance() {
//----------------------------------
    if( instance_ == 0 ) instance_ = new CsRichOne();
    return instance_;
  }

//===========================================================================
  CsRichOne::CsRichOne() {
//------------------------
    CsRegistry reg;
    reg.EOJRegistration( this );
    flag_ = 0;
    UpRICHJob_ = false;
    endMyFile_ = false;
  }

//===========================================================================
  CsRichOne::CsRichOne( const CsRichOne &rich ) {
//-----------------------------------------------

    cout << "RICHONE : CsRichOne CopyConstructor" << endl;
    instance_ = instance_;

    kEvent_ = rich.kEvent_;
    skipMyEvent_ = rich.skipMyEvent_;
    kMyWrEvent_ = rich.kMyWrEvent_;

    probPart_ = rich.probPart_;
    flag_ = rich.flag_;
    stopMyFileJob_ = rich.stopMyFileJob_;
    UpRICHJob_ = rich.UpRICHJob_;
    endMyFile_ = rich.endMyFile_;

    richTime_ = rich.richTime_;
    partsTime_ = rich.partsTime_;
    padsTime_ = rich.padsTime_;
    clusTime_ = rich.clusTime_;
    paphsTime_ = rich.paphsTime_;
    ringsTime_ = rich.ringsTime_;
    anasyTime_ = rich.anasyTime_;
    dispTime_ = rich.dispTime_;

    evGood_ = rich.evGood_;
    kMuonPart_ = rich.kMuonPart_;
    zVertex_ = rich.zVertex_;
    phiMass_ = rich.phiMass_;
  }

//===========================================================================
  void CsRichOne::print() const {
//-------------------------------
    cout << "CsRichOne flag = " << flag_ << endl;
  }

//===========================================================================
  CsRichOne::~CsRichOne() { }
//---------------------------


  extern "C" {
    void histoput_();
  }

//===========================================================================
  void CsRichOne::doRichOne() {
//-----------------------------

//--- December 1999

      //cout << "enter RICHONE" << endl;

      int iret;
      int k = 0;

//--- from "CsRCExecKeys"
      CsRCExeKeys *key = CsRCExeKeys::Instance();

      bool DoClu = key->DoClu();
      bool EventDisp = key->EventDisp();
      bool EventROOTDisp = key->EventROOTDisp();
      bool EventAnasy = key->EventAnasy();

//--- from "CsRCRecConst"
      CsRCRecConst *cons = CsRCRecConst::Instance();
      int kOffHi = cons->kOffHi();

      static int fh;
      int pSw;
      static int nPartTot;
      static int kEvent;

//--- from "CsRCHistos"
      CsRCHistos& hist = CsRCHistos::Ref();
      double xh, yh, wh;

      int chrTot;
      int chrPrt, chrPad, chrClu, chrPph, chrRng, chrAna, chrDsp;
      CsStopwatch chronos( 10 );

      static bool readMyFile;
      static bool proMyFile = true;
      static int  nProMyEvs = 1;
      static bool writeNtup;;

      static bool firstCall = true;
      if( firstCall ) {
	firstCall = false;

	richTime_ = 0.;
	partsTime_ = 0.;
	padsTime_ = 0.;
	clusTime_ = 0.;
	paphsTime_ = 0.;
	ringsTime_ = 0.;
	anasyTime_ = 0.;
	dispTime_ = 0.;

        nPartTot = 0;
        kEvent = 0;

	kMyWrEvent_ = 0;

        readMyFile = key->readMyFile();
	if( readMyFile ) nProMyEvs = CsRCExeKeys::Instance()->nProMyEvs();

        writeNtup = key->writeNtup();

	if( cons->printConsts() ) cons->print();
//---------------------------------------------
	CsRCDetectors::Instance()->setLocalGeo();
//@@--------------------------------------------
      }


//--  Find the RICH1 pointer
//    ----------------------
      string richName = "RIC1";
      CsRICH1Detector* rich = CsGeom::Instance()->getRich1Detector();
      if( rich == NULL ) {
	CsErrLog::mes( elFatal,
	  "RICHONE, CsRichOne::doRichOne() : RICH1 detector not found." );
      }

      stopMyFileJob_ = false;

      for( int kMyEve=0; kMyEve<nProMyEvs; kMyEve++ ) {
//----------------------------------------------------------------------------
      if( !proMyFile ) break;

      chrTot = chronos.start();

      static int nEvStartPrint = 1000000000;
//    -------------------------------------
      if( kEvent >= nEvStartPrint ) setFullPrintOn( 1 );
//                                  -------------------
      int kPrintRichRec = key->kPrintRichRec();
      int printFrequency = key->printFrequency();


//--- get EVENTS :
//    ------------
      kEvent_ = kEvent;

//--- skip myFile events ------------------------------------------------
//----------------------
//    WARNING : not always working! Works on regular events...
      skipMyEvent_ = false;
      static int nMySkip = 0;
      if( readMyFile ) nMySkip = CsRCExeKeys::Instance()->nSkipMyEvs();
      if( readMyFile  &&  nMySkip > 0 ) {
	if( kMyEve < nMySkip ) skipMyEvent_ = true;
	if( kMyEve == nMySkip ) {
  	  std::cout << "--- RICHONE :  CsRichOne::doRichOne : "
		    << nMySkip << " events skipped on MyFile" << std::endl;
	}
      }
      if( readMyFile ) {
//----- Warning readMyEvent needs the myFile initialization!
        //if( kEvent > 1  &&  kEvent <= 100 ) {
	  //readMyEvent();
	  //kEvent++;
	  //continue;
	//}
      }
      //if( readMyFile  &&  kMyEve > nMySkip ) setFullPrintOn( 1 );
//                                           -------------------

      if( skipMyEvent_ ) kPrintRichRec = 0;
      if( kPrintRichRec == 1 ) {
	cout << endl;
        cout << " RichOne ======  Event  " << kEvent
	     << "  ==============================================" << endl;
      }


//------- EVENT PARTICLES ---------
//---------------------------------
      chrPrt = chronos.start();

      int nPartEv = 0;
      Coral* coral = Coral::Instance();

//--- init the EventParticles object :
//    --------------------------------
      CsRCEventParticles* evparts = CsRCEventParticles::Instance();
//    ------------------------------------------------------------
      evparts->clearEventParticles();
      evparts->getEventParticles();

      if( key->endWriteMyFile() ) break;

      nPartEv = evparts->lParticles().size();
      nPartTot += nPartEv;

      partsTime_ += chronos.stop( chrPrt ) * 1000.;


//------- EVENT PADS ---------
//----------------------------
      chrPad = chronos.start();

//--- init the EventPads object :
//    ---------------------------
      CsRCEventPads* evpads = CsRCEventPads::Instance();
//    -------------------------------------------------
      evpads->clearEventPads();
      int nPadEv = 0;
      if( evpads->getEventPads() ) {
	nPadEv = evpads->lPads().size();
      }
      if( endMyFile_ ) continue;

      padsTime_ += chronos.stop( chrPad ) * 1000.;

      if( readMyFile ) physAppend();                //   040518
//                     ------------

//--- skip myFile events --- MOVED UP ------------- //   060816
//----------------------
      if( readMyFile ) {
	if( skipMyEvent_ ) {
	  if( stopMyFileJob_ ) break;
//----------------------------------
	  else {
	    kEvent++;
	    continue;
	  }
	}
      }


//--- write myFile for faster MC or Data tests --------------------------
//--------------------------------------------
      if( key->writeMyFile() ) {
        if( nPartEv > 0  &&  nPadEv > 0 ) writeMyFile();      // (011121)
//                                        -------------
      }
//-----------------------------------------------------------------------


//--- PARTICLE selection and analysis :
//    ---------------------------------
      if( nPartEv > 0 ) evparts->partAnalysis();
//    -----------------------------------------


//--- process current event particles with previous event pads (test!) --
//    'SCRAMBLED'
//--------------------------------------------------------------------
      bool padPrev = key->padScrambled();
//@@------------------------------------
      if( padPrev ) {
	static list<CsRCPad*> lPadsPrev;
	//cout << " Event  " << kEvent << endl;
        //cout << endl << "Previous :" << endl;
        //list<CsRCPad*>::iterator ia;
        //for( ia=lPadsPrev.begin(); ia!=lPadsPrev.end(); ia++) (*ia)->print();
        list<CsRCPad*> lPadsTemp = lPadsPrev;
	//cout << endl << "Temp :" << endl;
	//for( ia=lPadsTemp.begin(); ia!=lPadsTemp.end(); ia++) (*ia)->print();
        lPadsPrev = evpads->lPads();
	//cout << endl << "New Prev (=Current) :" << endl;
	//for( ia=lPadsPrev.begin(); ia!=lPadsPrev.end(); ia++) (*ia)->print();
	//cout << endl << "Current :" << endl;
	//evpads->print();
	evpads->clearPadList();
	evpads->copyPadList( lPadsTemp );
	//cout << endl << "Actual :" << endl;
	//evpads->print();
	static bool firstEve = true;
	if( firstEve ) { 
	  firstEve = false;
	  kEvent++;
	  if( readMyFile ) continue;
	  else return;
	}
      }
//-----------------------------------------------------------------------


//--- test on MC track hits ---------------------------------------------
//-------------------------
      bool trkTest = false;
      if( trkTest ) testMCTracks();
//-----------------------------------------------------------------------


      //     ******************** PHOTONS <-> TRACKS ********************
      // - This block performs, as I(Y.B.) understands, PHOTONS <-> TRACKS
      //  association.
      // - Was, in CsRichOne.cc,v.161, conditioned by condition (I), viz.:

      //if( nPadEv > 0 ) {

      //  which triggers a segmentation violation,
      //   - when condition (I) "nPadEv>0" is not fulfilled,
      //   - while condition (II) "nPartEv>0" is fulfilled,
      //  in "CsRCEventAnalysis::setProbsToTrack", which is called upon
      //  condition (II), but not necessarily (I), being fulfilled, because it
      //  bypasses "CsRCEventPartPhotons::clearEventPartPhotons".
      // - In addition, condition (I) is per se suspicious, at the present
      //  stage but also maybe at other stages, for it contradicts our
      //  objective (which is elsewhere in the code labelled "070430") to
      //  compute, and output, likelihoods even if no single photon can be
      //  associated to the track.
      // => Therefore I(Y.B.) comments out the condition (I) supra.
      {

//----- EVENT CLUSTERS ---------
//------------------------------
        chrClu = chronos.start();

//----- get the EventClusters object pt :
//      ---------------------------------
        CsRCEventClusters* evclus = CsRCEventClusters::Instance();
//      ---------------------------------------------------------
	evclus->clearEventClusters();
	evclus->getEventClusters();

        clusTime_ += chronos.stop( chrClu ) * 1000.;



//----- EVENT PARTICLE-PHOTONS ---------
//--------------------------------------
        chrPph = chronos.start();

//----- get the PartPhotons object pt :
//      -------------------------------
        CsRCEventPartPhotons* evphos = CsRCEventPartPhotons::Instance();
//      ---------------------------------------------------------------
        evphos->clearEventPartPhotons();
        if( nPartEv > 0 ) evphos->getEventPartPhotons();

	int nPaPhot = evphos->lPartPhotons().size();

        paphsTime_ += chronos.stop( chrPph ) * 1000.;



//----- EVENT RINGS ---------
//-----------------------------
        chrRng = chronos.start();

//----- get the EventRings object pt :
//      ------------------------------
        CsRCEventRings* evrings = CsRCEventRings::Instance();
//     ----------------------------------------------------
        evrings->clearEventRings();

	if( key->SearchForRings() ) {
	  if( nPaPhot > 0 ) evrings->getEventRings();
	}

        int nRingEv = evrings->lRings().size();

        ringsTime_ += chronos.stop( chrRng ) * 1000.;



//----- EVENT ANALYSIS ---------
//--------------------------------
//----- perform the EventAnalysis and PID :
//      -----------------------------------
        chrAna = chronos.start();
        if( EventAnasy ) {

	  CsRCEventAnalysis* ana = CsRCEventAnalysis::Instance();
	  //if( nRingEv > 0 ) ana->doEveAnalysis();   //   041116
//                            --------------------
	  ana->doEveAnalysis();                       //   041116
//        --------------------

        }
        anasyTime_ += chronos.stop( chrAna ) * 1000.;



//----- EVENT DISPLAY ---------
//-----------------------------
//----- do the EventDisplay :
//      ---------------------
        chrDsp = chronos.start();

        if( EventDisp ) {

          CsRCEventDisplay oEveDispv;
          oEveDispv.doEveDisplay();
//        ------------------------

	}
        dispTime_ += chronos.stop( chrDsp ) * 1000.;

//----- do the interactive (Andrea's) EventDisplay :
//      --------------------------------------------
        if( EventROOTDisp ) {

	  CsRichOneDisplay::Instance()->doRichOneDisplay();
//        ------------------------------------------------
	}



//----- FILL RICH NTUPLE ---------
//--------------------------------   MOVED down   070626 !
	//^if( writeNtup ) {
  	  //if( nRingEv > 0 ) CsRCnTup::Instance()->fillOutput();
  	  //^CsRCnTup::Instance()->fillOutput();         //   041116
//        //^----------------------------------
	//^}


      }       /* end if on nPadEv */


//--- OUTPUT BUFFER TO TRACKS  ---------
//-------------------------------------- added   070430
      if( nPartEv > 0 ) {
	CsRCEventAnalysis* ana = CsRCEventAnalysis::Instance();
	ana->setProbsToTrack();
//      ----------------------
      }


//--- FILL RICH NTUPLE ---------
//------------------------------ moved here   070626
      if(  nPartEv > 0  &&  nPadEv > 0 ) {
	if( writeNtup ) {
	  CsRCnTup::Instance()->fillOutput();
//        ----------------------------------
	}
      }


//pro-memoria
	//pp        evparts->printMem();
	//pp        evpads->printMem();
	//pp        evclus->printMem();
	//pp        evrings->printMem();
	//pp	    cout << endl;

	//      }       /* end if on nPadEv,nPartEv */

      richTime_ += chronos.stop( chrTot ) * 1000.;


      kEvent++;
      if( kEvent%printFrequency == 0 )  {
        cout << endl;
        iret = printf("events = %i   particles = %i \n",kEvent,nPartTot );
	std::cout << setprecision( 7 );
	std::cout << "refractive index " << cons->CFRefInd();
	if( UpRICHJob_ ) std::cout << " ( " << cons->CFRefIndVS() << " ) ";
	std::cout << std::endl;
	std::cout << setprecision( 0 );
	cout << "rich time = " << richTime_ << endl;
        cout << "  parts time = " << partsTime_;
        cout << "  pads time = " << padsTime_;
        cout << "  clus time = " << clusTime_;
        cout << "  part-phos time = " << paphsTime_;
        cout << "  rings time = " << ringsTime_;
        cout << "  analysis time = " << anasyTime_;
        cout << "  display time = " << dispTime_;
	std::cout << setprecision( 2 );
	cout << endl;
      }
      //printMemory( kEvent );      //   crash fixed? (24/1/01)

      //if( !readMyFile ) checkOutput();
//                      -------------

//----------------------------------------------------------------------------
      if( !readMyFile ) break;
//---------------------------
      if( stopMyFileJob_ ) break;
//------------------------------
      }
      if( readMyFile ) proMyFile = false;

      if( readMyFile ) {
	end();
//      -----
	std::cout << endl;
        std::string str = 
	  "RICHONE, CsRichOne::doRichOne() : end of reading myFile";
	std::cout << str << std::endl;
	str = "-------------------------------------------------------";
	std::cout << str << std::endl;
	std::cout << endl;
	//CsErrLog::Instance()->mes( elFatal, str );
	//exit(0);                       //   out 090211, see main.cc...
      }

      //cout << "exit RICHONE" << endl;

  }


//=========================================================================
  bool CsRichOne::end() {
//-----------------------

//--- June  2000


    CsRCExeKeys *key = CsRCExeKeys::Instance();

//----------------------------------------------------------------
    cout<< endl << endl;
    cout << "RICHONE end of processing :" << endl;
    cout << "---------------------------" << endl;
//----------------------------------------------------------------

    if( key->writeMyFile() ) {
      cout << "  " << kMyWrEvent_ << " events written to myFile" << endl;
      closeMyFile();
//-----------------
    }

    if( key->MCarloEvent() ) MCRecEffic();
//                           ------------

    dataRecEffic();
//----------------

//- print refractive index :
    double CFRefInd = 0.;
    double sigma = 0.;
    //if( CsRCEventAnalysis::Instance()->fitCFRefInf( CFRefInd, sigma ) ) {
    if( CsRCEventAnalysis::Instance()->fitGCFRefInd( CFRefInd, sigma ) ) {
//                                     -------------------------------
      int ksf = 0;
      double sigmaw = sigma;
      for( int k=1; k<=20; k++ ) {
	sigmaw *= 10.;
	if( sigmaw > 1. ) {
	  ksf = k;
	  break;
	}
      }

      
      /*
      double norm = pow( 10., double( ksf+1 ) );
      //cout << CFRefInd << "  " << norm << endl;
      long index = long( CFRefInd *norm ) + int( norm );
      //cout << index << endl;
      long double dindex = index;
      dindex /= norm;
      //CFRefInd = double( index ) / norm;
      //cout <<  dindex << "  " << CFRefInd << endl;
      cout << endl;
      cout << "  Average  C4F10  Refr.Index = ";
      cout << setprecision( ksf+1 );
      cout << dindex << endl;
      //cout << setprecision( ksf+1 ) << dindex << endl;
      stringstream n0;
      for( int k=1; k<ksf; k++ ) n0 << '0';
      cout << "                           +/- ";
      cout << "0."<< n0.str() << int( sigmaw ) << endl;
      */

      int pindex = int( CFRefInd * 1000000. );
      std::cout << "  Average  C4F10  Refr.Index = "
		<< " 1.00" << pindex << std::endl;

//--- print momentum thresholds :
      CFRefInd += 1.;
      cout << setprecision( 3 );
      cout << endl;
      cout << "  Momentum thresholds for: elec., muon,  pion,  kaon,  proton:"
  	   << endl;
      cout << "                           ";
      for( int kPaTy=2; kPaTy<=14; kPaTy+=3 ) {
        double massIpo = CsRCRecConst::Instance()->massPartv()[kPaTy];
        double momThresh = massIpo / sqrt( CFRefInd*CFRefInd - 1.);
        cout << momThresh << "   ";
      }
      cout << "  (GeV/c)" << endl;
    }

    execTime();
//------------

    bool AliMirror = key->AliMirror();
    CsRCMirrors *mirr = CsRCMirrors::Instance();
    //if( AliMirror ) mirr->CsRCMirrors::fitAliMirrors();
//----------------------------------------------------
    //if( AliMirror ) mirr->CsRCMirrors::fitAliMirrAll();
//----------------------------------------------------

//----------------------------------------------------------------
    cout<< endl << endl;
    cout << "RICHONE end of job :" << endl;
    cout << "--------------------" << endl;
//----------------------------------------------------------------

    return ( true );

  }



//=========================================================================
  void CsRichOne::MCRecEffic() {
//------------------------------

//--- June  2000


      CsRCHistos& hist = CsRCHistos::Ref();
      CsRCEventAnalysis* ana = CsRCEventAnalysis::Instance();
      CsRCExeKeys *key = CsRCExeKeys::Instance();
      bool MCMoni = key->MCMoni();
      bool PartIdent = key->PartIdent();

      if( MCMoni ) {
        float effic = 0.;
        float num = ana->nMCRecPart();
        float den = ana->nMCProcPart();
        if( den != 0 ) {
          effic = num / den;
	  effic = int(effic*1000)/1000.;
	}
//----------------------------------------------------------------
        cout << "  MC reconstruction efficiency = " << effic
	     << "  (" << num << "/" << den << ")" << endl;
//----------------------------------------------------------------
      }

      if( PartIdent && hist.hRC1655 ) {
      double binC[36];
      size_t nBinX = 4;
      size_t nBinY = 4;
      size_t kBin;
      size_t kBin0;
      size_t cBin = 1;
      for( size_t kBiny=1; kBiny <= nBinY; kBiny++ ) {
	kBin0 = kBiny * (nBinX + 2);         // un.flw & ov.flw !
        for( size_t kBinx=1; kBinx <= nBinX; kBinx++ ) {
	  kBin = kBin0 + kBinx;
          binC[cBin] = 0;
          if( hist.hRC1655 ) binC[cBin] = hist.hRC1655->GetBinContent( kBin );
	  //	  cout << kBin << "/ " << binC[cBin] << "   ";
          cBin++;
        }
	//	cout << endl;
      }
      double SbinC;
      double pionID = 0.;
      double kaonID = 0.;
      double protonID = 0.;
      SbinC = binC[2] + binC[6] + binC[10] + binC[14];
      //      cout << binC[6] << "  " << SbinC << endl;
      if( SbinC != 0. ) pionID = binC[6]/SbinC;
      pionID = int(pionID*1000)/1000.;
      SbinC = binC[3] + binC[7] + binC[11] + binC[15];
      //      cout << binC[11] << "  " << SbinC << endl;
      if( SbinC != 0. ) kaonID = binC[11]/SbinC;
      kaonID = int(kaonID*1000)/1000.;
      SbinC = binC[4] + binC[8] + binC[12] + binC[16];
      //      cout << binC[16] << "  " << SbinC << endl;
      if( SbinC != 0. ) protonID = binC[16]/SbinC;
      protonID = int(protonID*1000)/1000.;
      double pionBK = 0.;
      double kaonBK = 0.;
      double protonBK = 0.;
      SbinC = binC[6] + binC[7] + binC[8];
      //      cout << binC[7]+binC[8] << "  " << SbinC << endl;
      if( SbinC != 0. ) pionBK = (binC[7]+binC[8])/SbinC;
      pionBK = int(pionBK*1000)/1000.;
      SbinC = binC[10] + binC[11] + binC[12];
      //      cout << binC[10]+binC[12] << "  " << SbinC << endl;
      if( SbinC != 0. ) kaonBK = (binC[10]+binC[12])/SbinC;
      kaonBK = int(kaonBK*1000)/1000.;
      SbinC = binC[14] + binC[15] + binC[16];
      //      cout << binC[14]+binC[15] << "  " << SbinC << endl;
      if( SbinC != 0. ) protonBK = (binC[14]+binC[15])/SbinC;
      protonBK = int(protonBK*1000)/1000.;

      //      cout << endl;
//----------------------------------------------------------------
      cout << "  MC pion ID efficiency   = " << pionID
	   << ",     ID bkgr = " << pionBK << endl;
      cout << "  MC kaon ID efficiency   = " << kaonID
	   << ",     ID bkgr = " << kaonBK << endl;
      cout << "  MC proton ID efficiency = " << protonID
	   << ",     ID bkgr = " << protonBK << endl;
//----------------------------------------------------------------
      }

  }


//=========================================================================
  void CsRichOne::dataRecEffic() {
//--------------------------------

//--- June  2000


      CsRCEventAnalysis* ana = CsRCEventAnalysis::Instance();
      CsRCExeKeys *key = CsRCExeKeys::Instance();
      bool DataMoni = key->DataMoni();

      if( DataMoni ) {
        float effic = 0.;
        float num = ana->nDaRecPart();
        float den = ana->nDaProcPart();
        if( den != 0 ) {
          effic = num / den;
	  effic = int(effic*1000)/1000.;
	}

//------------------------------------------------------------------------
        //      cout << endl;
        cout << "  Data monitored reconstruction efficiency = " << effic
	     << "  (" << num << "/" << den << ")" << endl;
//-----------------------------------------------------------------------
      }

  }


//=========================================================================
  void CsRichOne::execTime() {
//----------------------------

//--- April  2001


    if( kEvent_ > 0 ) {
      double event = kEvent_;
      double richTime = richTime_ / event;
      double partsTime = partsTime_ / event;
      double padsTime = padsTime_ / event;
      double clusTime = clusTime_ / event;
      double paphsTime = paphsTime_ / event;
      double ringsTime = ringsTime_ / event;
      double anasyTime = anasyTime_ / event;
      double dispTime = dispTime_ / event;
      richTime = int(richTime*10)/10.;
      partsTime = int(partsTime*10)/10.;
      padsTime = int(padsTime*10)/10.;
      clusTime = int(clusTime*10)/10.;
      paphsTime = int(paphsTime*10)/10.;
      ringsTime = int(ringsTime*10)/10.;
      anasyTime = int(anasyTime*10)/10.;
      dispTime = int(dispTime*10)/10.;

//----------------------------------------------------------------
      cout << endl;
      cout << "RICHONE times :" << endl;
      cout << "---------------" << endl;
      cout << "  input of particles  " << partsTime << "  ms/ev" << endl;
      cout << "  input of pads       " << padsTime << "  ms/ev" << endl;
      cout << "  pad clustering      " << clusTime << "  ms/ev" << endl;
      cout << "  photons to particle " << paphsTime << "  ms/ev" << endl;
      cout << "  ring recognition    " << ringsTime << "  ms/ev" << endl;
      cout << "  monitor and PID     " << anasyTime << "  ms/ev" << endl;
      cout << "  ring display        " << dispTime << "  ms/ev" << endl;
      cout << "  TOTAL RICHONE ----- " << richTime << "  ms/ev" << endl;
    } else {
      cout << "RICHONE times :" << endl;
      cout << "---------------" << endl;
      cout << "No statistics  " << endl;
    }
//----------------------------------------------------------------

  }


//=========================================================================
  bool CsRichOne::getPartProbs( CsTrack* track, double* partProbs ) {
//-------------------------------------------------------------------

//- Paolo - October  2000
//  rev. March  2003
   

    double *partPro;
    bool calc = false;

    CsRCRecConst *cons = CsRCRecConst::Instance();

    //static int nProb = 15;
    static const int nProb = cons->outBufferSize();
    for( int k=0; k<nProb; k++ ) partProbs[k] = 0.;

      if( cons->likeType() == "ALL" ) {
      CsRCEventPartPhotons* evepapho = CsRCEventPartPhotons::Instance();
      list<CsRCPartPhotons*> lPaPhotons = evepapho->lPartPhotons();
      list<CsRCPartPhotons*>::iterator ipf;
      for( ipf=lPaPhotons.begin(); ipf!=lPaPhotons.end(); ipf++ ) {
	if( !(*ipf) ) continue;
	//^if( !(*ipf)->flag() ) continue;                    //   070613
	if( !(track == (*ipf)->pPart()->pTrack()) ) continue;
	partPro = (*ipf)->partProbs();
	for( int k=0; k<nProb; k++ ) partProbs[k] = partPro[k];
	//for( int k=0; k<8; k++ ) partProbs[k] = partPro[k];
	//CsRCRing* ring = (*ipf)->pRing();
	//if( ring  &&  ring->flag() ) {
	//  partPro = ring->partProbs();
	//  for( int k=8; k<nProb; k++ ) partProbs[k] = partPro[k];
	//}
	//else {
	//  for( int k=8; k<nProb; k++ ) partProbs[k] = 0.;
	//}
        //for( int k=0; k<nProb; k++ ) cout <<  partPro[k] << "  " ;
        //cout << endl;
	calc = true;
      }
    }

    if( cons->likeType() == "RING" ) {
      CsRCEventRings* evrings = CsRCEventRings::Instance();
      list<CsRCRing*> lRings = evrings->lRings();
      list<CsRCRing*>::iterator ir;
      for( ir=lRings.begin(); ir!=lRings.end(); ir++ ) {
	if( !(*ir) ) continue;
	if( !(*ir)->flag() ) continue;
	if( !(track == (*ir)->pPartPhot()->pPart()->pTrack()) ) continue;
	partPro = (*ir)->partProbs();
	for( int k=0; k<nProb; k++ ) partProbs[k] = partPro[k];
        //for( int k=0; k<nProb; k++ ) cout <<  partPro[k] << "  " ;
        //cout << endl;
	calc = true;
      }
    }

    if( cons->likeType() == "CHISQ" ) {
      CsRCEventRings* evrings = CsRCEventRings::Instance();
      list<CsRCRing*> lRings = evrings->lRings();
      list<CsRCRing*>::iterator ir;
      for( ir=lRings.begin(); ir!=lRings.end(); ir++ ) {
	if( !(*ir) ) continue;
	if( !(*ir)->flag() ) continue;
	if( !(track == (*ir)->pPartPhot()->pPart()->pTrack()) ) continue;
	partPro = (*ir)->partProbs();
	for( int k=0; k<8; k++ ) partProbs[k] = 0.;
	for( int k=8; k<nProb; k++ ) partProbs[k] = partPro[k];
        //for( int k=0; k<nProb; k++ ) cout <<  partPro[k] << "  " ;
        //cout << endl;
	calc = true;
      }
    }

    return calc;

  }


//===========================================================================
  bool CsRichOne::getPartProbs( CsMCTrack* track, double* partProbs ) {
//---------------------------------------------------------------------

//--- Paolo - October  2000


    double *partPro;
    bool calc = false;
    CsRCRecConst *cons = CsRCRecConst::Instance();
    static const int nProb = cons->outBufferSize();

    CsRCEventRings* evrings = CsRCEventRings::Instance();
    list<CsRCRing*> lRings = evrings->lRings();
    list<CsRCRing*>::iterator ir;
    for( ir=lRings.begin(); ir!=lRings.end(); ir++ ) {
      if( track == (*ir)->pPartPhot()->pPart()->pMCTrack() ) {
	if( (*ir)->flag() ) {
	  partPro = (*ir)->partProbs();
	  int nPhotons = (*ir)->lPhotons().size();
	  partPro[9] = nPhotons;
	  //for( int k=0; k<15; k++ ) partProbs[k] = partPro[k];
	  for( int k=0; k<nProb; k++ ) partProbs[k] = partPro[k];
	  calc = true;
	}
      }
    }
    return calc;

  }


//===========================================================================
  void CsRichOne::printMemory( int kEvent ) const {
//-------------------------------------------------

//--- Paolo - November  2000


      CsRCHistos& hist = CsRCHistos::Ref();
      double xh, yh, wh;

      static int count = 0;

      if( kEvent % CsRCExeKeys::Instance()->printFrequency() == 0 ) {

	CsRCEventParticles::Instance()->printMem();
	cout << " ";
	CsRCEventPads::Instance()->printMem();
	cout << "      ";
	CsRCEventClusters::Instance()->printMem();
	cout << "  ";
	cout << "                                ";
	CsRCEventRings::Instance()->printMem();
        cout << endl;
	list<CsRCParticle*> lParticles = 
	  CsRCEventParticles::Instance()->lParticles();
	if( lParticles.size() != 0 ) {
	  cout << " " << &lParticles << " " << lParticles.front()
	       << " " << lParticles.back();
	}
	list<CsRCPad*> lPads = CsRCEventPads::Instance()->lPads();
	if( lPads.size() != 0 ) {
	  cout << "  " << &lPads << " " << &lPads.front()
	       << " " << &lPads.back();
	}
        list<CsRCCluster*> lClusters =
	  CsRCEventClusters::Instance()->lClusters();
	if( lClusters.size() != 0 ) {
	  cout << "  " << &lClusters << " " << lClusters.front()
	       << " " << lClusters.back();
	}
	list<CsRCPartPhotons*> lPartPhotons = 
	  CsRCEventPartPhotons::Instance()->lPartPhotons();
        list<CsRCRing*> lRings = CsRCEventRings::Instance()->lRings();
	if( lPartPhotons.size() != 0 ) {
  	  cout << "  " << &lPartPhotons << " " << lPartPhotons.front()
	       << " " << lPartPhotons.back();
	  cout << "  " << &lRings << " " << lRings.front()
	       << " " << lRings.back();
	}
        cout << endl;
	//        cout << endl;

        if( hist.hRC1031 ) {
        count++;
	xh = count + 0.5;
	yh = 0.5;
	size_t ih = 0;
	if( lParticles.size() != 0 ) {
	  ih = size_t(lParticles.front());
	  wh = ih/10000. - 22000.;
	  hist.hRC1031->Fill( xh, yh, wh );
//hh      --------------------------------
	  ih = size_t(lParticles.back());
	  yh++;
	  wh = ih/10000. - 22000.;
	  hist.hRC1031->Fill( xh, yh, wh );
//hh      --------------------------------
	}
	if( lPads.size() != 0 ) {
	  ih = size_t(&lPads.back());
	  yh++;
	  wh = ih/10000. - 22000.;
	  hist.hRC1031->Fill( xh, yh, wh );
//hh      --------------------------------
	}
	if( lClusters.size() != 0 ) {
	  ih = size_t(lClusters.back());
	  yh++;
	  wh = ih/10000. - 22000.;
	  hist.hRC1031->Fill( xh, yh, wh );
//hh      --------------------------------
	}
	if( lPartPhotons.size() != 0 ) {
  	  ih = size_t(lPartPhotons.back());
	  yh++;
	  wh = ih/10000. - 22000.;
	  hist.hRC1031->Fill( xh, yh, wh );
//hh      --------------------------------
	  list<CsRCPhoton*> lPhotonP = lPartPhotons.back()->lPhotons();
	  ih = size_t(lPhotonP.back());
	  yh++;
	  wh = ih/10000. - 22000.;
	  hist.hRC1031->Fill( xh, yh, wh );
//hh      --------------------------------
	  ih = size_t(lRings.back());
	  yh++;
	  wh = ih/10000. - 22000.;
	  hist.hRC1031->Fill( xh, yh, wh );
//hh      --------------------------------
	  list<CsRCPhoton*> lPhotonR = lRings.back()->lPhotons();
	  ih = size_t(lPhotonR.back());
	  yh++;
	  wh = ih/10000. - 22000.;
	  hist.hRC1031->Fill( xh, yh, wh );
//hh      --------------------------------
	}
	}

      }

  }


//===========================================================================
  void CsRichOne::testMCTracks() {
//--------------------------------

//--- Paolo - February  2001


	cout << endl;
        cout << "----------------------" << endl;

        Coral* coral = Coral::Instance();
        CsRICH1Detector* rich = CsGeom::Instance()->getRich1Detector();
        list<CsMCTrack*> lTrackTk;
        CsRCEventParticles* evparts = CsRCEventParticles::Instance();
        list<CsRCParticle*> lParticles = evparts->lParticles();
        list<CsRCParticle*>::iterator ipk;
        for( ipk=lParticles.begin(); ipk!=lParticles.end(); ipk++ ) {
	  lTrackTk.push_back( (*ipk)->pMCTrack() );
        }
        list<CsMCTrack*>::iterator it;
	for( it=lTrackTk.begin(); it!=lTrackTk.end(); it++ ) {
          list<CsMCHit*> mchits = (*it)->getMCHits();
          list<CsMCHit*>::iterator ih;
	  int nCerPhot = 0;
          for( ih=mchits.begin(); ih!=mchits.end(); ih++ ) {
	    CsDetector* det =  dynamic_cast<CsDetector*>((*ih)->getDet());
	    if( det == dynamic_cast<CsDetector*>( rich ) ) {
	      nCerPhot++;
	    }
	  }
	  if( nCerPhot > 0 ) {
	    cout << "track  " << (*it)->getGnum() << "   "
		 << nCerPhot << "  photons" << endl;
	  }
	}
        list<CsMCTrack*> lTrackHt;
        list<CsDigit*> mydigits = rich->getMyDigits();
        list<CsDigit*>::iterator id;
        int kDig = 0;
        for( id=mydigits.begin(); id!=mydigits.end(); id++ ) {
          //-- tests MC tracks to tracks from MC RICH hits 
          list<CsMCHit*> lHits = 
              dynamic_cast<CsMCDigit*>(*id)->CsMCDigit::getHits();
          list<CsMCHit*>::iterator ih;
          if( !lHits.empty() ) {
            for( ih=lHits.begin(); ih!=lHits.end(); ih++ ) {
              CsMCTrack* it = (*ih)->getMCTrack();
  	      int iTrack = it->getGnum();
	      //	      if( (*ih)->getOrigin() > 0 ) {
	      //cout << kDig << "  iTrack Hit  " << iTrack << endl;
	      //	      }
	      //cout << (*ih)->getP() << endl;
              list<CsMCTrack*>::iterator ikh;
       	      bool neq = true;
              for( ikh=lTrackHt.begin(); ikh!=lTrackHt.end(); ikh++ ) {
                if( iTrack == (*ikh)->getGnum() ) { neq = false;   break; }
              }
              if( neq ) { lTrackHt.push_back( it ); }
	    }
	  }
	  kDig++;
	}
        //cout << "lTrackHt     " << lTrackHt.size() << endl;

        list<CsMCTrack*>::iterator ikh;
	for( ikh=lTrackHt.begin(); ikh!=lTrackHt.end(); ikh++ ) {
          //  cout << (*ikh)->getGnum() << "  ";
	}
        //cout << endl;

        //-- checks on input tracks
        list<CsMCTrack*> lTrackEq;
        list<CsMCTrack*> lTrackTkS;
        list<CsMCTrack*> lTrackHtS;
        bool Exe = true;
        if( Exe ) {
	  list<CsMCTrack*>::iterator ikk;
	  list<CsMCTrack*>::iterator ikh;
	  for( ikk=lTrackTk.begin(); ikk!=lTrackTk.end(); ikk++ ) {
	    bool equ = false;
	    for( ikh=lTrackHt.begin(); ikh!=lTrackHt.end(); ikh++ ) {
	      if( (*ikh)->getGnum() == (*ikk)->getGnum() ) {
                equ = true;   break;
	      }
	    }
   	    if( equ ) { lTrackEq.push_back( (*ikh) ); }
	    if( !equ ) {
	      lTrackTkS.push_back( (*ikk) );
	      //	    lTrackTk.remove( (*ikk) );   // NO
	      //	    lTrackHt.erase( ikh );       // NO
	    }
	  }
	  for( ikh=lTrackHt.begin(); ikh!=lTrackHt.end(); ikh++ ) {
	    bool equ = false;
	    for( ikk=lTrackTk.begin(); ikk!=lTrackTk.end(); ikk++ ) {
	      if( (*ikh)->getGnum() == (*ikk)->getGnum() ) {
                equ = true;  break;
	      }
	    }
	    if( !equ ) { lTrackHtS.push_back( (*ikh) ); }
	  }
	  bool pri = true;
	  if( pri ) {
            list<CsMCTrack*>::iterator ik;
	    cout << endl;
	    //            cout << "-----------------" << endl;
	    list<CsMCTrack*> lTrackAll = coral->getMCTracks();
	    cout << "MC Tracks all         " << lTrackAll.size() << ":  ";
	    for( ik=lTrackAll.begin(); ik!=lTrackAll.end(); ik++ ) {
	      cout << (*ik)->getGnum() << "  ";
	    }
      	    cout << endl;
	    cout << "MC Tracks from Parts   " << lTrackTk.size() << ":  ";
            for( ik=lTrackTk.begin(); ik!=lTrackTk.end(); ik++ ) {
	      cout << (*ik)->getGnum() << "  ";
	    }
      	    cout << endl;
	    cout << "MC Tracks from Hits    " << lTrackHt.size() << ":  ";
            for( ik=lTrackHt.begin(); ik!=lTrackHt.end(); ik++ ) {
	      cout << (*ik)->getGnum() << "  ";
	    }
            cout << endl;
	    //	    if( !lTrackHtS.empty() ) {
	      //              cout << "-----------------" << endl; 
              cout << "common Tracks        " << lTrackEq.size() << ":  ";
              for( ik=lTrackEq.begin(); ik!=lTrackEq.end(); ik++ ) {
	        cout << (*ik)->getGnum() << "  ";
	      }
	      cout << endl;
  	      cout << "Tracks left (Parts)  " << lTrackTkS.size() << endl;;
	      if( lTrackTkS.size() > 0 ) {
		cout << "=====================================" << endl;
	      }
	      for( ik=lTrackTkS.begin(); ik!=lTrackTkS.end(); ik++ ) {
		//cout << (*ik)->getParticle()->getGeantNumber() << "  ";
	        cout << (*ik)->getGnum() << "  ";
	        cout << (*ik)->getP() << "  ";
	        cout << (*ik)->getM() << "  ";
	        const CsMCVertex* vertex = (*ik)->getInVertex();
	        if( vertex != 0 ) {
		  cout << "  Vertex. x: " << vertex->getX()
		       << " y: " << vertex->getY()
		       << " z: " << vertex->getZ() << endl;
		}
	        list<CsMCHit*> lHitsS = (*ik)->getMCHits();
	        list<CsMCHit*>::iterator ihS;
                for( ihS=lHitsS.begin(); ihS!=lHitsS.end(); ihS++ ) {
		  cout << (*ihS)->getZ() << "  ";
		}
		cout << endl;
	      }
	      cout << "Tracks left (Hits)   " << lTrackHtS.size() << endl;
	      for( ik=lTrackHtS.begin(); ik!=lTrackHtS.end(); ik++ ) {
		//cout << (*ik)->getParticle()->getGeantNumber() << "  ";
	        cout << (*ik)->getGnum() << "  ";
	        cout << (*ik)->getP() << "  ";
	        cout << (*ik)->getM() << "  ";
	        const CsMCVertex* vertex = (*ik)->getInVertex();
	        if( vertex != 0 ) {
		  cout << "  Vertex. x: " << vertex->getX()
		       << " y: " << vertex->getY()
		       << " z: " << vertex->getZ() << endl;
		}
	        list<CsMCHit*> lHitsS = (*ik)->getMCHits();
	        list<CsMCHit*>::iterator ihS;
                for( ihS=lHitsS.begin(); ihS!=lHitsS.end(); ihS++ ) {
		  //                  cout << (*ihS)->getZ() << "  ";
		}
		//	        cout << endl;
	      }
	      cout << endl;
	      //	    }
	  }
	}

  }


//===========================================================================
  extern "C" {
    void writebuff_( const char*,
		     int&,
                     float*, float*, float*, float*, float*, float*, 
                     float*, float*, float*,
                     float*, float*, float*, float*, float*, float*,
                     float*, float*, float*, float*, float*, float*,
                     float*,
		     int*, float*,
		     float*, float*, int*, float*,         //   020322
		     int*, int*, int*,
		     int&, 
		     int*, int*, int*, float*, int );
    void writeapp_( const char*, int&, int&, float&, float&, int );
    void closewrtfile_();
  }

//===========================================================================
  void CsRichOne::writeMyFile() {
//-------------------------------

//--- Paolo - February  2001, rev. March 2001


    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      kMyWrEvent_ = 0;
    }

    CsRCEventParticles* evparts = CsRCEventParticles::Instance();
    list<CsRCParticle*> lParts = evparts->lParticles();
    int nPartEv = evparts->lParticles().size();
    list<CsRCParticle*>::iterator ip;
    int NPAEV = 30;
    float xpa[NPAEV], ypa[NPAEV], zpa[NPAEV];
    float lpa[NPAEV], mpa[NPAEV], npa[NPAEV];
    float xpo[NPAEV], ypo[NPAEV], zpo[NPAEV];
    float xpem[NPAEV], ypem[NPAEV], zpem[NPAEV];
    float lpem[NPAEV], mpem[NPAEV], npem[NPAEV];
    float xpew[NPAEV], ypew[NPAEV], zpew[NPAEV];
    float lpew[NPAEV], mpew[NPAEV], npew[NPAEV];
    float ppa[NPAEV];
    int iTrack[NPAEV], iPartT[NPAEV], ioExWd[NPAEV];
    int nPhoCer[NPAEV];
    float thetaCer[NPAEV];
    float zHelix0[NPAEV], zHelix1[NPAEV], chiSq[NPAEV];    //   020322
    int nClus[NPAEV];                                      //   020322
    int kPart = 0;
    for( ip=lParts.begin(); ip!=lParts.end(); ip++ ) {
      //^if( !(*ip)->flag() ) continue;
      if( kPart >= NPAEV ) break;                          //   020722
      xpa[kPart] = (*ip)->xa();
      ypa[kPart] = (*ip)->ya();
      zpa[kPart] = (*ip)->za();
      lpa[kPart] = (*ip)->ld();
      mpa[kPart] = (*ip)->md();
      npa[kPart] = (*ip)->nd();
      xpo[kPart] = (*ip)->xo();
      ypo[kPart] = (*ip)->yo();
      zpo[kPart] = (*ip)->zo();
      xpem[kPart] = (*ip)->vPosEmP().x();
      ypem[kPart] = (*ip)->vPosEmP().y();
      zpem[kPart] = (*ip)->vPosEmP().z();
      lpem[kPart] = (*ip)->vDirEmP().x();
      mpem[kPart] = (*ip)->vDirEmP().y();
      npem[kPart] = (*ip)->vDirEmP().z();
      xpew[kPart] = (*ip)->vPosExW().x();
      ypew[kPart] = (*ip)->vPosExW().y();
      zpew[kPart] = (*ip)->vPosExW().z();
      lpew[kPart] = (*ip)->vDirExW().x();
      mpew[kPart] = (*ip)->vDirExW().y();
      npew[kPart] = (*ip)->vDirExW().z();
      ppa[kPart] = (*ip)->mom() * (*ip)->charge();         //   011121
      iTrack[kPart] = (*ip)->iTrack();
      iPartT[kPart] = (*ip)->iPartT();
      ioExWd[kPart] = (*ip)->ioExWd();
      nPhoCer[kPart] = (*ip)->nPhoCer();
      thetaCer[kPart] = (*ip)->thetaCer();
//--- as from 020322
      zHelix0[kPart] = (*ip)->zHelix0();
      zHelix1[kPart] = (*ip)->zHelix1();
      nClus[kPart] = (*ip)->nClus();
      chiSq[kPart] = (*ip)->chiSq();
      kPart++;
    }
    int nPartW = kPart;

//- as from 040511 - provisional?
    float index = CsRCRecConst::Instance()->CFRefInd();
    int iindex = int( (index-1.) * 1000000. );
    if( nPartW > 0 ) ioExWd[0] += iindex*1000;
//- as from 070130 - provisional?
    float indexVS = CsRCRecConst::Instance()->CFRefIndVS();
    if( UpRICHJob_  &&  !CsRCExeKeys::Instance()->MCarloEvent() ) {
      int iindex = int( (index-1.) * 10000000. + 0.5 );
      int iindexVS = int( (indexVS-1.) * 10000000. + 0.5 );
      if( nPartW > 0 ) ioExWd[0] = iindex + iindexVS*100000;
    }
    //cout << setprecision(7) << index << "  " << indexVS << "  "
    // << ioExWd[0] << "  " << CsRCRecConst::Instance()->CFRefInd()
    // << "  " << CsRCRecConst::Instance()->CFRefIndVS() << endl;
//@@-----------------------------

    CsRCEventPads* evpads = CsRCEventPads::Instance();
    list<CsRCPad*> lPads = evpads->lPads();
    int nPadEv = evpads->lPads().size();
    list<CsRCPad*>::iterator id;
    int MAXPAD = 5000;
    int ict[MAXPAD],ixt[MAXPAD],iyt[MAXPAD];
    float pht[MAXPAD];
    int kPad = 0;
    for( id=lPads.begin(); id!=lPads.end(); id++ ) {
      //^if( !(*id)->flag() ) continue;
      if( kPad >= MAXPAD ) break;                          //   020722
      ict[kPad] = (*id)->ic() + 1;
      ixt[kPad] = (*id)->ix() + 1;
      iyt[kPad] = (*id)->iy() + 1;
      //^pht[kPad] = (*id)->PH();
      pht[kPad] = (*id)->PHPack();
      kPad++;
    }
    int nPadW = kPad;

    const char *fileName;
    string myFileName = CsRCExeKeys::Instance()->WmyFileName();
    fileName = myFileName.c_str();

    writebuff_( fileName,
		nPartW,
                xpa, ypa, zpa, lpa, mpa, npa, xpo, ypo, zpo,
		xpem, ypem, zpem, lpem, mpem, npem,
		xpew, ypew, zpew, lpew, mpew, npew,
                ppa,
		nPhoCer, thetaCer,
		zHelix0, zHelix1, nClus, chiSq,            //   020322
		iTrack, iPartT, ioExWd,
		nPadW,
                ict, ixt, iyt, pht, strlen( fileName ) );

    if( CsRCExeKeys::Instance()->readMyFile() ) {
//  -------------------------------------------            //   040616
      ///if( CsRCExeKeys::Instance()->myFileType() >= 3 ) {
      if( CsRCExeKeys::Instance()->myFileType() == 3 ) {
	writeapp_( fileName, evGood_, kMuonPart_, zVertex_, phiMass_,
		   strlen( fileName ) );
      }
    }

    kMyWrEvent_++;

  }


//===========================================================================
  void CsRichOne::closeMyFile() {
//-------------------------------

//--- Paolo - March 2002

    closewrtfile_();

  }


//===========================================================================
  void CsRichOne::setFullPrintOn( const int level ) {
//---------------------------------------------------

//--- Paolo - May 2001


      CsRCExeKeys *key = CsRCExeKeys::Instance();

      int lev = 0;
      lev = level; if( lev > 1 ) lev = 1;
      key->setPrintRichRec( lev );
      lev = level; if( lev > 2 ) lev = 2;
      key->setPrintEventParts( lev );
      lev = level; if( lev > 1 ) lev = 1;
      key->setPrintEventSimul( lev );
      lev = level; if( lev > 2 ) lev = 2;
      key->setPrintEventPads( lev );
      lev = level; if( lev > 2 ) lev = 2;
      key->setPrintEventClus( lev );
      lev = level; if( lev > 3 ) lev = 3;
      key->setPrintPartPhotons( lev );
      lev = level; if( lev > 4 ) lev = 4;
      key->setPrintEventRings( lev );
      lev = level; if( lev > 1 ) lev = 1;
      key->setPrintEventAnalysis( lev );
      lev = level; if( lev > 1 ) lev = 1;
      key->setPrintPartProbs( lev );
      lev = level; if( lev > 1 ) lev = 1;
      key->setPrintEventDisplay( lev );

  }


//===========================================================================
  void CsRichOne::checkOutput() {
//-------------------------------

//- Paolo - January 2003


//- loop on tracks :
//  ----------------
    Coral* coral = Coral::Instance();
    CsEvent* event = coral->getEvent();
    list<CsTrack*> tracks = event->getTracks();
    list<CsTrack*>::iterator it;
 
    //static int nProb = 15;
    CsRCRecConst *cons = CsRCRecConst::Instance();
    static const int nProb = cons->outBufferSize();
    double partProbs[nProb];
    const double* partProbsT;
    int kTrack = 0;
    for( it=tracks.begin(); it!=tracks.end(); it++ ) {
      if( !(*it) ) continue;

      if( getPartProbs( (*it), partProbs ) ) {
	cout << "RichOne : " << kTrack << " - ";
	for( int k=0; k<nProb; k++ ) cout << partProbs[k] << "  ";
	cout << endl;
      }
 
      if( (*it)->hasRich1Probs() ) {
	partProbsT = (*it)->getRich1Probs();
	cout << "CsTrack : " << kTrack << " - ";
	for( int k=0; k<nProb; k++ ) cout << partProbsT[k] << "  ";
	cout << endl;
      }

      kTrack++;
    }

  }


//===========================================================================
  void CsRichOne::physSelec() {
//-----------------------------

//- Paolo - May 2004


    CsRCExeKeys *key = CsRCExeKeys::Instance();
    CsRCRecConst *cons = CsRCRecConst::Instance();

    if( !key->physSelection() ) return;
//  ----------------------------------
    if( key->readMyFile() ) return;
//  ------------------------------

    double xh = 0.;
    static CsHist1D* hRC401 = NULL;
    static CsHist1D* hRC402 = NULL;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      key->acknoMethod( "CsRichOne::physSelec" );

      CsHistograms::SetCurrentPath("/RICH");
//@@---------------------------------------
      static int kOffHi = cons->kOffHi();
      string hTitle;
      stringstream hN401;
      hN401 << kOffHi + 401;
      hTitle = "z target";
      hRC401 = new CsHist1D( hN401.str(), hTitle, 200, -1400., 600. );
      stringstream hN402;
      hN402 << kOffHi + 402;
      hTitle = "phi mass";
      hRC402 = new CsHist1D( hN402.str(), hTitle, 100, 1.0, 1.1 );
      CsHistograms::SetCurrentPath("/");
//@@-----------------------------------
    }

    //std::cout << "CsRichOne::physSelec" << std::endl;

    const char *fileName;
    string myFileName = key->WmyFileName();
    fileName = myFileName.c_str();

    const list<CsVertex*> &vrts = CsEvent::Instance()->getVertices();
    const list<CsTrack* > &Trks = CsEvent::Instance()->getTracks();
    const vector<CsParticle* > &parts = CsEvent::Instance()->getParticles();

    int nPartEv = CsRCEventParticles::Instance()->lParticles().size();
    int nTracks = Trks.size();
    //std::cout << nPartEv << "  " << nTracks << std::endl;

    bool bGood = true;
    CsVertex *vPr(0);                   // pointer to primary vertex
    const CsTrack *pToMu(0);            // pointer to beam
    const CsTrack *pToMuP(0);           // pointer to mu'
    vector<CsParticle*>::const_iterator ip;
    for( ip=parts.begin(); ip!=parts.end(); ip++ ) {
      if( (*ip)->getType() != CsParticle::SPECIAL ) continue;
      const CsTrack *trk = (*ip)->getTrack();
      if( trk == 0 ) {
        cout<<"ERROR: particle is special but no track!!!"<<endl;
        continue;
      }
      vector<CsHelix> v = trk->getHelices();
      if( v[0].getZ() < 0 ) {          // beam
	pToMu = trk;
	continue;
      }
      CsVertex* vrt = trk->getFirstVertex();
      if( vrt==0 || !vrt->isPrimary() ) continue;  // no primary
      vPr = vrt;                       // OK there is primary with mu/mu'
      pToMuP = trk;                    // mu'
      vector<CsHelix> vHelix = trk->getHelices();
      double cop = vHelix[0].getCop();
      double mom = 1./fabs( cop );
      //std::cout << "muP mom = " << mom << std::endl;
      break;
    }
    if( vPr == 0 ) bGood = false;

    float zVertex = 0.;
    float PhiMass = 0.;
    double Xv(0),Yv(0),Zv(0);          // coordinates of primary vertex
    CsHelix h1, h2;
    CsTrack* pToTrk[10];
    if( bGood  &&  vPr != 0 ) {
      Xv = vPr->getX();
      Yv = vPr->getY();
      Zv = vPr->getZ();
      //std::cout << Xv << "  " << Yv << "  " << Zv << std::endl;
      xh = Zv;
      if( hRC401 ) hRC401->Fill( xh );
//hh               ------------------
      zVertex = Zv;
      const list<CsTrack*> &trks = vPr->getTracks();
      list<CsTrack*>::const_iterator it = trks.begin();
      int kTrk = 0;
      //CsTrack* pToTrk[10];
      for( it=trks.begin(); it!=trks.end(); it++) {
	if( (*it) == pToMu ) continue;
	if( (*it) == pToMuP ) continue;
	pToTrk[kTrk] = (*it);
	kTrk++;
	if( kTrk >= 9 ) break;
      }
      if( kTrk == 2 ) {
	//std::cout << kTrk << std::endl;
        h1 = pToTrk[0]->getHelices()[0];
        h2 = pToTrk[1]->getHelices()[0];
	//std::cout << kTrk << std::endl;
      } else {
	bGood = false;
      }
      //if( bGood &&  !(kTrk == 2) ) std::cout << "--- " << kTrk << std::endl;
      if( bGood ) {
	//std::cout << "----" << std::endl;
	for( it=trks.begin(); it!=trks.end(); it++) {
	  vector<CsHelix> vHelix = (*it)->getHelices();
	  double cop = vHelix[0].getCop();
	  double mom = 1./fabs( cop );
	  //std::cout << mom << "  ";
	}
	//std::cout << std::endl;
      }
    }

    //std::cout << zVertex << std::endl;
    if( bGood  &&  h1(5)/fabs(h1(5)) != h2(5)/fabs(h2(5)) ) {
      //std::cout << h1(5) << "  " << h2(5) << std::endl;
      CsHelix ht1,ht2;
      h1.Extrapolate(Zv,ht1);
      h2.Extrapolate(Zv,ht2);
      TVector3 v1,v2;
      double Cop, dXdZ, dYdZ;
      Cop = ht1.getCop();
      dXdZ= ht1.getDXDZ();
      dYdZ= ht1.getDYDZ();
      v1[2] = 1/fabs(Cop)/sqrt(1+dXdZ*dXdZ+dYdZ*dYdZ);
      v1[0] = v1[2] * dXdZ;
      v1[1] = v1[2] * dYdZ;
      Cop = ht2.getCop();
      dXdZ= ht2.getDXDZ();
      dYdZ= ht2.getDYDZ();
      v2[2] = 1/fabs(Cop)/sqrt(1+dXdZ*dXdZ+dYdZ*dYdZ);
      v2[0] = v2[2] * dXdZ;
      v2[1] = v2[2] * dYdZ;
      TVector3 v0(v1+v2);
//Hep3Vector vv1( v1[0], v1[1], v1[2] );
//Hep3Vector vv2( v2[0], v2[1], v2[2] );
//Hep3Vector vv0( v0[0], v0[1], v0[2] );
//std::cout << std::endl << "physSelec-1  " << std::endl;
//std::cout << setprecision( 5 );
//std::cout << vv1 << "  " << vv2 << "  " << vv0 << std::endl;
//std::cout << 1./h1.getCop() << "  " << 1./ht1.getCop() << "  "
//  	    << 1./h2.getCop() << "  " << 1./ht2.getCop() << std::endl;
      double KMass = 0.49368;
      double E1 = sqrt( KMass*KMass + v1.Mag2() );
      double E2 = sqrt( KMass*KMass + v2.Mag2() );
      double E0 = E1 + E2;
      PhiMass = sqrt(E0*E0-v0.Mag2());
//std::cout << E1 << "  " << E2 << "  " << E0 << std::endl;
//std::cout << PhiMass << std::endl;
      xh = PhiMass;
      if( hRC402 ) hRC402->Fill( xh );
//hh               ------------------
    } else {
      bGood = false;
    }


    if( !key->writeMyFile() ) return;

    int evGood = 0;
    if( !bGood ) evGood = 1;

    int kMuonPart = 0;
    int kK1 = 0;
    int kK2 = 0;
    int kPart = 1;
    int nRP = 0;
    int kRP = 0;
    if ( bGood ) {
      CsRCEventParticles* evparts = CsRCEventParticles::Instance();
      list<CsRCParticle*> lParticles = evparts->lParticles();
      list<CsRCParticle*>::iterator ipk;
      for( ipk=lParticles.begin(); ipk!=lParticles.end(); ipk++ ) {
	if( (*ipk)->pTrack() == pToMuP ) {
	  kMuonPart = kPart;
        }
	if( (*ipk)->pTrack() == pToTrk[0] ) {
	  kK1 = kPart;
	  if( !(*ipk)->pTrack()->hasRich1Probs() ) 
	    nRP += int( pow( 10.0, kRP ) );
	  kRP++;
	}
	if( (*ipk)->pTrack() == pToTrk[1] ) {
	  kK2 = kPart;
	  if( !(*ipk)->pTrack()->hasRich1Probs() ) 
	    nRP += int( pow( 10.0, kRP ) );
	  kRP++;
	}
	kPart++;
      }
      kMuonPart = kMuonPart + kK1*100 + kK2*10000;
      for( ipk=lParticles.begin(); ipk!=lParticles.end(); ipk++ ) {
	//std::cout <<  (*ipk)->mom() << "  ";
      }
      //std::cout << std::endl;
    }
    if( evGood == 0 ) evGood -= nRP;
    //if( evGood < 0 ) std::cout << evGood << "  " << kPart-1 << "  " 
    //     		         << kMuonPart << std::endl;
    //std::cout << evGood << "  " << kPart-1 << "  " 
    //    << kMuonPart << "  " << zVertex << "  " << PhiMass << std::endl;

    //int evGood_ = evGood;
    //int kMuonPart_ = kMuonPart;
    //float zVertex_ = zVertex;
    //float phiMass_ = PhiMass;
    evGood_ = evGood;
    kMuonPart_ = kMuonPart;
    zVertex_ = zVertex;
    phiMass_ = PhiMass;
//std::cout << "physSelec-2  " << std::endl;
//std::cout << evGood_ << "  " << kMuonPart_ << "  " << zVertex_
//          << "  " << phiMass_ << std::endl;

    writeapp_( fileName, evGood, kMuonPart, zVertex, PhiMass, 
	       strlen( fileName )  );

//^^    CsRCEventRings::Instance()->testPhiKzero();
//^^------------------------------------------

    return;

  }



#include "coral_config.h"
#if USE_RFIO
#  include <shift.h>
#endif


//===========================================================================
  void CsRichOne::physAppend() {
//------------------------------

//- Paolo - May 2004


    CsRCExeKeys *key = CsRCExeKeys::Instance();

    ///if( !(key->myFileType() >= 3) ) return;
    if( !(key->myFileType() == 3) ) return;
//  --------------------------------------

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      key->acknoMethod( "CsRichOne::physAppend" );
    }

    int bytePerInt = sizeof(int);
    int bytePerFloat = sizeof(float);
    static int fh;
    int iret;
    int pSw;

    int ibuff[10];
    float fbuff[2000];

    //std::cout << "CsRichOne::physAppend" << std::endl;

    fh = CsRCEventParticles::Instance()->myFile();

    const char *fileName;
    string myFileName = key->WmyFileName();
    fileName = myFileName.c_str();

// WARNING: to be checked and updated for g-files written under SLC5 !!!
// 100309 --------------------------------------------------------------
    pSw = 4 * bytePerInt;
    iret = read( fh, ibuff, pSw );
// -----------------------------------
    evGood_ = ibuff[2];
    kMuonPart_ = ibuff[3];

    pSw = 2 * bytePerFloat;
    iret = read( fh, fbuff, pSw );
// -----------------------------------
    zVertex_ = fbuff[0];
    phiMass_ = fbuff[1];
    //std::cout << evGood_ << "  " << kMuonPart_ << "  "
    //<< zVertex_ << "  " << phiMass_ << std::endl;

    return;

  }


//===========================================================================
  bool CsRichOne::readMyEvent() {
//-------------------------------

//- Paolo - August  2006

    CsRCExeKeys *key = CsRCExeKeys::Instance();
    CsRCRecConst *cons = CsRCRecConst::Instance();

    int bytePerInt = sizeof(int);
    int bytePerFloat = sizeof(float);
    static int fh;
    int iret;
    int pSw;

    int ibuff[10];
    float fbuff[2000];

    bool dump = true;
    bool dumpt = true;
    bool dumpd = true;
//@@------------------

    fh = CsRCEventParticles::Instance()->myFile();
    if( fh == 0 ) return  false;
//  ---------------------------

    pSw = 3 * bytePerInt;
    iret = read( fh, ibuff, pSw );
// -----------------------------------
    if( iret < 12 ) return  false;
//  -----------------------------
    int nPartEv = ibuff[2];
    if( dump ) std::cout << "nPartEv = " << nPartEv << std::endl;
    if( nPartEv > 30 ) nPartEv = 30;
//@@-------------------------------

    for( int kPart=0; kPart<nPartEv; kPart++ ) {

      pSw = 9 * bytePerFloat;
      iret = read( fh, fbuff, pSw );
// -------------------------------------
      if( dumpt ) {
        Hep3Vector vPoInPart( fbuff[0],fbuff[1],fbuff[2] );
        Hep3Vector vDcInPart( fbuff[3],fbuff[4],fbuff[5] );
        Hep3Vector vPoExPart( fbuff[6],fbuff[7],fbuff[8] );
	std::cout << "vPoInPart, vDcInPart = "
		  << vPoInPart << " "
		  << setprecision(5) << vDcInPart;
	std::cout << "    " << setprecision(2);
	std::cout << "vPoExPart = " << vPoExPart << std::endl;
      }
      pSw = 12 * bytePerFloat;
      iret = read( fh, fbuff, pSw );
// -------------------------------------
      if( dumpt ) {
	Hep3Vector vPosEmP( fbuff[0],fbuff[1],fbuff[2] );
	Hep3Vector vDirEmP( fbuff[3],fbuff[4],fbuff[5] );
	Hep3Vector vPosExW( fbuff[6],fbuff[7],fbuff[8] );
	Hep3Vector vDirExW( fbuff[9],fbuff[10],fbuff[11] );
	std::cout << "vPosEmP, vDirEmP = " << vPosEmP << " "
		  << setprecision(5) << vDirEmP;
	std::cout << "   " << setprecision(2);
	std::cout << "vPosExW, vDirExW = " << vPosExW << " "
		  << setprecision(5) << vDirExW;
	std::cout << setprecision(2) << std::endl;
      }
      pSw = 1 * bytePerFloat;
      iret = read( fh, fbuff, pSw );
// -------------------------------------
      if( dumpt ) {
        float momPart = fabs( fbuff[0] );
	std::cout << "momPart = " << momPart;
	std::cout << std::endl;
      }

      int nPhoCer = 0;
      double thetaCer = 0.;
      if( key->myFileType() >= 1 ) {
        pSw = bytePerInt;
        iret = read( fh, ibuff, pSw );
// ---------------------------------------
        if( dumpt ) nPhoCer = ibuff[0];
        pSw = bytePerFloat;
        iret = read( fh, fbuff, pSw );
// ---------------------------------------
	if( dumpt ) {
          thetaCer = fbuff[0];
	  std::cout << "nPhoCer, thetaCer = " << nPhoCer << " " << thetaCer;
	  std::cout << std::endl;
	}
      }

      float zHelix0 = 0.;
      float zHelix1 = 0.;
      int nClus = 0;
      float chiSq = 0.;
      if( key->myFileType() >= 2 ) {
        pSw = 2 * bytePerFloat;
        iret = read( fh, fbuff, pSw );
// ---------------------------------------
	if( dumpt ) {
          zHelix0 = fbuff[0];
          zHelix1 = fbuff[1];
	}
        pSw = 1 * bytePerInt;
        iret = read( fh, ibuff, pSw );
// ---------------------------------------
	if( dumpt ) nClus = ibuff[0];
        pSw = 1 * bytePerFloat;
        iret = read( fh, fbuff, pSw );
// -----------------------------------------
        if( dumpt ) {
	  chiSq = fbuff[0];
  	  std::cout << "zHelix0, zHelix1 = " << zHelix0 << " " << zHelix1;
	  std::cout << std::endl;
	}
      }

      pSw = 3 * bytePerInt;
      iret = read( fh, ibuff, pSw );
// -------------------------------------
      if( dumpt ) {
        int iTrack = ibuff[0];
        int iPartTy = ibuff[1];
        int ioExWd = ibuff[2];
        std::cout << "iTrack, iPartTy, ioExWd = "
		  << iTrack << " " <<  iPartTy << " " << ioExWd;
        std::cout << std::endl;
      }
    }
    std::cout << std::endl;


    pSw = bytePerInt;
    iret = read( fh, ibuff, pSw );
//  ------------------------------------
    int nPadEv = ibuff[0];
    if( dump ) std::cout << "nPadEv = " << nPadEv << std::endl;
    if( nPadEv > 5000 ) nPadEv = 5000;
//@@---------------------------------

    for( int kPad=0; kPad<nPadEv; kPad++ ) {

      pSw = 3 * bytePerInt;
      iret = read( fh, ibuff, pSw );
//    ----------------------------------
      if( dumpd ) {
	int iCat = ibuff[0];
	int iXPad = ibuff[1];
	int iYPad = ibuff[2];
	iCat--;
	iXPad--;   iYPad--;
	std::cout << "iCat, iXPad, iYPad = "
		  << iCat << " " << iXPad << " " << iYPad;
	//std::cout << std::endl;
      }
      pSw = bytePerFloat;
      iret = read( fh, fbuff, pSw );
//--------------------------------------
      if( dumpd ) {
	float PHPad = fbuff[0];
	std::cout << "  ";
	std::cout << "PHPad = " << PHPad << "   ";
	//std::cout <<std::endl;
      }

    }

    return  true;
  }
