/*!
   \file    CsRCEventParticles.cc
   \-----------------------------
   \brief   CsRCEventParticles class implementation.
   \author  Paolo Schiavon
   \version 0.01
   \date    1 October 2000
*/

    #include <ostream>
    #include <cstdio>

    #include <unistd.h>
    #include <fcntl.h>
    #include <sys/types.h>
    #include <sys/stat.h>

    #include "CLHEP/Matrix/Matrix.h"
    #include "CLHEP/Matrix/SymMatrix.h"
    #include "CLHEP/Matrix/DiagMatrix.h"
    #include "CLHEP/Matrix/Vector.h"
//--------------------------------------

    #include "Coral.h"
    #include "CsTrack.h"
    #include "CsMCTrack.h"
    #include "CsMCRICH1Hit.h"
    #include "CsMCUtils.h"
    #include "CsOpt.h"

//---------------------------------
    #include "CsRCEventParticles.h"

    #include "CsRCMirrors.h"
    #include "CsRCMirrorNom.h"
    #include "CsRCDetectors.h"
    #include "CsRichOne.h"

    #include "CsRCEventAnalysis.h"

    #include "CsRCExeKeys.h"
    #include "CsRCRecConst.h"
    #include "CsRCHistos.h"
//---------------------------------

    #include "TMtx.h"
    #include "CsRCTrackMomFit.h"

    #include <CLHEP/Vector/ThreeVector.h>
//---------------------------------------

    using namespace std;
    using namespace CLHEP;

    CsRCEventParticles* CsRCEventParticles::instance_ = 0;

//==========================================================================
    CsRCEventParticles* CsRCEventParticles::Instance() {
//------------------------------------------------------
      if( instance_ == 0 ) instance_ = new CsRCEventParticles();
      return instance_; 
    }

//==========================================================================
    CsRCEventParticles::CsRCEventParticles() {
//--------------------------------------------
      lParticles_.clear();
      flag_ = true;
      myFile_ = 0;
      myFileRun_ = 0;
    }

//==========================================================================
    CsRCEventParticles::CsRCEventParticles( 
//-----------------------------------------
                        const CsRCEventParticles &evparts ) {
      cout << " RICHONE : CsRCEventParticles CopyConstructor" << endl;
      instance_ = evparts.instance_;
      lParticles_ = evparts.lParticles_;
      flag_ = evparts.flag_;
      flagSimul_ = evparts.flagSimul_;
      myFile_ = evparts.myFile_;
      myFileRun_ = evparts.myFileRun_;
    }

//==========================================================================
    CsRCEventParticles& CsRCEventParticles::operator=( 
//----------------------------------------------------
                        const CsRCEventParticles &evparts ) {
      if( this != &evparts ) {
	instance_ = evparts.instance_;
	lParticles_ = evparts.lParticles_;
	flag_ = evparts.flag_;
        flagSimul_ = evparts.flagSimul_;
        myFile_ = evparts.myFile_;
	myFileRun_ = evparts.myFileRun_;
      }
      return ( *this );
    }

//==========================================================================
    void CsRCEventParticles::clearEventParticles() {
//--------------------------------------------------
      list<CsRCParticle*>::iterator ip;
      for( ip=lParticles_.begin(); ip!=lParticles_.end(); ip++ ) delete *ip;
      lParticles_.clear();
    }

//==========================================================================
    void CsRCEventParticles::print() const {
//------------------------------------------
      list<CsRCParticle*>::const_iterator ip;
      for( ip=lParticles_.begin(); ip!=lParticles_.end(); ip++ ) {
        (*ip)->print();
      }
    }

//==========================================================================
    CsRCEventParticles::~CsRCEventParticles() {
//---------------------------------------------
      clearEventParticles();
    }


//==========================================================================
  void CsRCEventParticles::getEventParticles() {
//----------------------------------------------


//--- Paolo  -  October 2000;  Rev. November 2001


      CsRCExeKeys *key = CsRCExeKeys::Instance();
      CsRCRecConst *cons = CsRCRecConst::Instance();

//--- from "CsRCHistos"
      CsRCHistos& hist = CsRCHistos::Ref();
      double xh, yh, wh;

      Coral* coral = Coral::Instance();

      static bool firstCall = true;
      static bool readMyFile = false;
      static bool tracking = false;
      if( firstCall ) {
	firstCall = false;
        key->acknoMethod( "CsRCEventParticles::getEventParticles" );

        readMyFile = key->readMyFile();
//@@----------------------------------
	tracking = CsOpt::Instance()->getOpt( "", "make tracking" );
      }

//--- construct the EventParticles object :
//    -------------------------------------

      if( key->MCarloEvent() ) {

//----- myfile MC event (old and new type) :
//      ------------------------------------
        if( readMyFile ) setMyParticle();
//                       ---------------

        else {

//------- COMGeant event :
//        ----------------
	  if( tracking ) setMCRecParticle();
//                       ------------------
	  else  setMCParticle();
//              ---------------
        }
      }

      if( key->DataEvent() ) {

//----- myfile Data event (new type) :
//      ------------------------------
        if( readMyFile ) setMyParticle();
//                       ---------------

        else {

//------- Raw DATA event :
//        ----------------
	  setDataParticle();
//        -----------------
	  //	setMCRecParticle();        //   DATA (Display) tests
	  //	lParticles_.clear();       //   DATA (Display) tests
	}
      }

//--- monitoring histograms :
//    -----------------------
      int nPartEv = lParticles_.size();
      xh = nPartEv;
      if( hist.hRC1531 ) hist.hRC1531->Fill( xh );
//hh                     ------------------------

//--- conditional prints :
//    --------------------
      int kPrintEventParts = key->kPrintEventParts();
      if( kPrintEventParts == 1 ) {
	cout << endl;
        cout << " Particles : " << nPartEv << endl;
      }

  }



#include "coral_config.h"
#if USE_RFIO
#  include <shift.h>
#endif


//==========================================================================
  void CsRCEventParticles::setMyParticle() {
//------------------------------------------


//--- Paolo  -  November 2000;  
//    Rev.      November 2001


      CsRCExeKeys *key = CsRCExeKeys::Instance();
      CsRCRecConst *cons = CsRCRecConst::Instance();

      int bytePerInt = sizeof(int);
      int bytePerFloat = sizeof(float);
      static int fh;
      int iret = 0;
      int pSw;

      //bool dumpBf = true;
      bool dumpBf = false;

      int ibuff[10];
      float fbuff[2000];

      static bool firstCall = true;
      if( firstCall ) {
	firstCall = false;
        key->acknoMethod( "CsRCEventParticles::setMyParticle" );

	int nn;
	if( !initMyFile( nn ) ) return;
//          -----------------

	fh = myFile_;
	CsRichOne::Instance()->setEndMyFile( false );

      }   /* end if firstCall */

      int nPartEv = 0;

      if( CsRichOne::Instance()->endMyFile() ) {
	if( !initMyFile( nPartEv ) ) return;
//          ----------------------
	fh = myFile_;
	CsRichOne::Instance()->setEndMyFile( false );
      }
      else {
//----- skip myFile record header
        ///pSw = 3 * bytePerInt;
        ///iret = read( fh, ibuff, pSw );
        ///nPartEv = ibuff[2];
	int nSkip = 1;
	if( key->readGf5() ) nSkip = 2;
        pSw = nSkip * bytePerInt;
        iret = read( fh, ibuff, pSw );
// ----------------------------------
	if( dumpBf ) std::cout << "myFile record length " << ibuff[0]
			       << std::endl;
        ///if( iret < 12 ) {
        ///if( iret < pSw ) {
	///  if( !initMyFile( nPartEv ) ) return;
        ///              ----------------------
  	///  fh = myFile_;
	///  CsRichOne::Instance()->setEndMyFile( false );
        ///}
//----- read nPartEv
        pSw = 1 * bytePerInt;
        iret = read( fh, ibuff, pSw );
// ----------------------------------
	nPartEv = ibuff[0];
      }
      if( dumpBf ) std::cout << "CsRCEventParticles::setMyParticle"
			     << std::endl;
      if( dumpBf ) std::cout << "nPartEv " << nPartEv << std::endl;
      //if( nPartEv > 30 ) {
      while( nPartEv > 30 ) {
	std::cout << std::endl;
        std::cout << " RICHONE, CsRCEventParticles::setMyParticle() :"
	  ///<< " stop reading myFile, nPartEv > 30  "
	  ///<< nPartEv << std::endl;
        	  << " stop reading myFile, EOF reached (nPartEv > 30)"
		  << std::endl;
	if( !initMyFile( nPartEv ) ) return;
//          ----------------------
	fh = myFile_;
	CsRichOne::Instance()->setEndMyFile( false );
      }

      for( int kPart=0; kPart < nPartEv; kPart++ )  {

	pSw = 9 * bytePerFloat;
	iret = read( fh, fbuff, pSw );
// ----------------------------------
	Hep3Vector vPoInPart( fbuff[0],fbuff[1],fbuff[2] );
	Hep3Vector vDcInPart( fbuff[3],fbuff[4],fbuff[5] );
	Hep3Vector vPoExPart( fbuff[6],fbuff[7],fbuff[8] );
	if( dumpBf ) {
	  //std::cout << "vPoInPart " << vPoInPart << std::endl;
	}

	pSw = 12 * bytePerFloat;
	iret = read( fh, fbuff, pSw );
// ----------------------------------
	Hep3Vector vPosEmP( fbuff[0],fbuff[1],fbuff[2] );
	Hep3Vector vDirEmP( fbuff[3],fbuff[4],fbuff[5] );
	Hep3Vector vPosExW( fbuff[6],fbuff[7],fbuff[8] );
	Hep3Vector vDirExW( fbuff[9],fbuff[10],fbuff[11] );

	pSw = 1 * bytePerFloat;
	iret = read( fh, fbuff, pSw );
// ----------------------------------
        float momPart = fabs( fbuff[0] );             //   011121
	int charge = 0;
	if( fbuff[0] > 0. ) charge = +1;
	if( fbuff[0] < 0. ) charge = -1;              //   011121
	//cout << vPoInPart << endl;
	//cout << vDcInPart << endl;
	//cout << vPoExPart << endl;
	//cout << momPart << endl;

	int nPhoCer = 0;
	double thetaCer = 0.;
        if( key->myFileType() >= 1 ) {
          pSw = bytePerInt;
	  iret = read( fh, ibuff, pSw );
// ------------------------------------
  	  nPhoCer = ibuff[0];
          pSw = bytePerFloat;
	  iret = read( fh, fbuff, pSw );
// ------------------------------------
  	  thetaCer = fbuff[0];
	}

//----- as from 020322
	float zHelix0 = 0.;
	float zHelix1 = 0.;
	int nClus = 0;
	float chiSq = 0.;
	if( key->myFileType() >= 2 ) {
	  pSw = 2 * bytePerFloat;
	  iret = read( fh, fbuff, pSw );
// ------------------------------------
	  zHelix0 = fbuff[0];
	  zHelix1 = fbuff[1];
	  pSw = 1 * bytePerInt;
	  iret = read( fh, ibuff, pSw );
// ------------------------------------
	  nClus = ibuff[0];
	  pSw = 1 * bytePerFloat;
	  iret = read( fh, fbuff, pSw );
// ------------------------------------
	  chiSq = fbuff[0];
	  //cout << zHelix0 << "  " << zHelix1 << endl;
	}

        pSw = 3 * bytePerInt;
        iret = read( fh, ibuff, pSw );
// ----------------------------------
  	int iTrack = ibuff[0];
        int iPartTy = ibuff[1];
        int ioExWd = ibuff[2];

//----- check skip event
//      ----------------
	if( CsRichOne::Instance()->skipMyEvent() ) continue;
//      ---------------------------------------------------

//----- as from 040511 - provisional?
	if( !key->MCarloEvent() ) {
//---------------------------------
	  if( myFileRun_ == 0 ) {
	    if( kPart == 0 ) {
  	      //int iindex = ioExWd / 1000;
	      //float index = 1. + float( iindex ) / 1000000.;
	      //cons->setCFRefInd( index );
//            //--------------------------
	      //ioExWd -= iindex * 1000;
	      float index = 0.;
	      float indexVS = 0.;
//----------- as from 070130 - provisional?
	      if( CsRichOne::Instance()->UpRICHJob() ) {
	        int iindexVS = ioExWd / 100000;
	        indexVS = 1. + float( iindexVS ) / 10000000.;
	        cons->setCFRefIndVS( indexVS );
//              ------------------------------
	        int iindex = ioExWd - iindexVS * 100000;
	        index = 1. + float( iindex ) / 10000000.;
	        cons->setCFRefInd( index );
//              --------------------------
	        if( key->selPMTonly() ) cons->setCFRefInd( indexVS );
	        if( key->selAPVonly() ) cons->setCFRefIndVS( index );
	      }
	      else {
		int iindex = ioExWd / 1000;
		index = 1. + float( iindex ) / 1000000.;
		cons->setCFRefInd( index );
//              --------------------------
		ioExWd -= iindex * 1000;
		indexVS = 0.;
	      }
	      //cout << setprecision(6);
	      //cout << ioExWd << " " << index << " " << cons->CFRefInd()
	      //<< "  " << indexVS << "  " << cons->CFRefIndVS() << endl;
	    }
	  }
	}
        int jPart = lParticles_.size();
        if( key->MCarloEvent() ) {
          if( key->myFileType() == 0 ) {
            lParticles_.push_back( new CsRCParticle( jPart,
		  		   vPoInPart, vDcInPart, vPoExPart, momPart,
		                   iTrack, iPartTy, ioExWd ) );
	  }
          if( key->myFileType() >= 1 ) {
            lParticles_.push_back( new CsRCParticle( jPart,
				   vPoInPart, vDcInPart, vPoExPart,
				   vPosEmP, vDirEmP, vPosExW, vDirExW,
				   momPart, nPhoCer, thetaCer,
		                   iTrack, iPartTy, ioExWd ) );
	  }
          flagSimul_ = true;
	}
        if( key->DataEvent() ) {

//------- set geometrical corrections --- 021113
	  bool bCorr = true;
	  if( bCorr) {
	    //cout << vDcInPart << "  ";
	    setGeoCorrections( vPoInPart, vDcInPart );
	    setGeoCorrections( vPosEmP, vDirEmP );
//          -----------------------------------------
	    //cout << vDcInPart << endl;
	  }

          if( key->myFileType() >= 1 ) {
            lParticles_.push_back( new CsRCParticle( jPart, NULL,
		  	           vPoInPart, vDcInPart, vPoExPart,
				   vPosEmP, vDirEmP, vPosExW, vDirExW,
				   momPart, charge ) );
	  }
          flagSimul_ = false;
	}

	lParticles_.back()->setQualityPars( zHelix0, zHelix1, nClus,
					    chiSq );
        if( CsRCExeKeys::Instance()->kPrintEventParts() == 2 ) {
	  lParticles_.back()->print();
	}

      }   /* end loop on particles: kPart */

  }


//==========================================================================
  void CsRCEventParticles::setMCParticle() {
//------------------------------------------

//--- Paolo  -  November 2000


      CsRCHistos& hist = CsRCHistos::Ref();
      int kh;
      double xh, yh, wh;

      CsRCExeKeys* key = CsRCExeKeys::Instance();
      CsRCRecConst *cons = CsRCRecConst::Instance();
      CsRICH1Detector* rich = CsGeom::Instance()->getRich1Detector();
      CsRCDetectors *dets = CsRCDetectors::Instance();

      static float partPathFr = 0.;
      static double zEntrWind = 0.;
      static double zExitWind = 0.;
      static double zPoPhot0 = 0.;
      list<CsRCMirrorNom*> lMirrNom = CsRCMirrors::Instance()->lMirrNom();
      list<CsRCPhotonDet*> lPhoDet = CsRCDetectors::Instance()->lPhoDet();
      static double cov[15];

      static bool firstCall = true;
      if( firstCall ) {
        firstCall = false;
	key->acknoMethod( "CsRCEventParticles::setMCParticle" );

	partPathFr = cons->partPathFr();
        zEntrWind = dets->zEntrWind();
        zExitWind = dets->zExitWind();
	double zMirrMin = lMirrNom.front()->vC0().z() +
	  lMirrNom.front()->RR() * cos( lPhoDet.front()->detAng() );
        zPoPhot0 = zEntrWind + partPathFr *(zMirrMin - zEntrWind);
	for( int k=0; k<15; k++ ) cov[k] = 0.;
	int j = 1;
	for( int k=0; k<15; k+=j ) { cov[k] = 1.; j++; }
      }


      //-- MC Tracks
      Coral* coral = Coral::Instance();
      list<CsMCTrack*> mctracks = coral->getMCTracks();
      list<CsMCTrack*>::iterator it;
      //-- loop on Monte Carlo true tracks
      list<CsMCHit*> lHits;
      for( it=mctracks.begin(); it!=mctracks.end(); it++ ) {

        //-- get the particle associated to the track
        CsMCParticle* part = (*it)->getParticle();
        //-- only charged particles
        if( part->getCharge() == 0 ) continue;

        list<CsMCHit*> mchits = (*it)->getMCHits();
        list<CsMCHit*>::iterator ih;
        list<CsMCHit*>::iterator ihMin = mchits.end();
        double zDiff = 0.;
        double zDiffMin = 1000000.;
        int nPhoCer = 0;
        double thetaCer = 0.;
        bool flagOrig = true;
        lHits.clear();
        for( ih=mchits.begin(); ih!=mchits.end(); ih++ ) {

	  CsDetector* det =  dynamic_cast<CsDetector*>((*ih)->getDet());
	  if( det != NULL ) {
	    zDiff = zEntrWind - (*ih)->getZ();
	    if( zDiff > 0. ) {
	      if( zDiff < zDiffMin ) {
		zDiffMin = zDiff;
		ihMin = ih;
	      }
	    }
	  }
	  if( det == dynamic_cast<CsDetector*>( rich ) ) {
	    nPhoCer++;
	    lHits.push_back( (*ih) );
	    if( (*ih)->getOrigin() != 0 ) flagOrig = false;
//cout << "rich  " << (*ih)->getZ() << "  ";
	    CsMCRICH1Hit* hit = dynamic_cast<CsMCRICH1Hit*>( (*ih) );
	    thetaCer += hit->getCherAngle();
//cout << hit->getX() << "  " << hit->getY() << "  " << hit->getCathode();
//cout << "  " << hit->getYdetDRS() << "  " << hit->getZdetDRS() << "  ";
//cout << hit->getCherAngle()*1000. << "  " << (*ih)->getOrigin() << endl;
//cout << hit->getPhotEnergy() << "    ";
//cout << hit->getZprod() << endl;
	    Hep3Vector vPosEmMC( hit->getXprod(), hit->getYprod(),
				 hit->getZprod() );
	  }
	  else { continue; }
	}
	if( nPhoCer != 0 ) {
          thetaCer /= float( nPhoCer );
          thetaCer *= 1000.;

          //- monitoring histograms :
	  double massPart = part->getMass();
          double momPart = (*mchits.begin())->getP().mag();
          double betaPart = momPart /
            sqrt( momPart*momPart + massPart*massPart );
          for( ih=mchits.begin(); ih!=mchits.end(); ih++ ) {
	    if( (*ih)->getOrigin() != 0 ) continue;
            CsDetector* det =  dynamic_cast<CsDetector*>((*ih)->getDet());
	    if( det == dynamic_cast<CsDetector*>( rich ) ) {
	      CsMCRICH1Hit* hit = dynamic_cast<CsMCRICH1Hit*>( (*ih) );
	      double theHit = hit->getCherAngle()*1000.;
	      double ePhot = hit->getPhotEnergy();

	      xh = theHit - thetaCer;
	      yh = thetaCer;
	      if( hist.hRC1061 ) hist.hRC1061->Fill( xh, yh );
//                               ----------------------------
	      xh = thetaCer;
	      yh = ePhot;
	      if( hist.hRC1062 ) hist.hRC1062->Fill( xh, yh );
//                               ----------------------------
	      xh = theHit - thetaCer;
	      yh = ePhot;
	      if( hist.hRC1063 ) hist.hRC1063->Fill( xh, yh );
//                               ----------------------------
              xh = nPhoCer;
	      yh = ePhot;
	      if( hist.hRC1066 ) hist.hRC1066->Fill( xh, yh );
//                               ----------------------------
	    }
	  }
	  xh = nPhoCer;
          yh = betaPart;
          if( hist.hRC1068 ) hist.hRC1068->Fill( xh, yh );
//                           ----------------------------
	}

        if( !flagOrig ) continue;
//@@----------------------------

	//-- skip particle under C threshold
        //if( nPhoCer == 0 ) continue;       // MC ONLY   out 030325   !!!
//@@-------------------------------
	//-- skip particle with too few photons
	//!!!	if( nPhoCer <= 20 ) continue;     // MC ONLY, test
//@@-----------------------------------------
        if( ihMin == mchits.end() ) continue;
        if( zDiffMin > 1000. ) continue;      // provisional
	//cout << (*ihMin)->getDet()->getId() << endl;

	float momPart = (*ihMin)->getP().mag();
        float momMinAcc = cons->momMinAcc();
	//-- low momentum particle cut
	if( momPart < momMinAcc ) continue;
//@@--------------------------------------
	//	CsRCEventAnalysis::Instance()->hitDisplay( 5000, lHits );

        double xDetUS = (*ihMin)->getX();
        double yDetUS = (*ihMin)->getY();
        double zDetUS = (*ihMin)->getZ();
//----- project forward to the RICH entrance window
        double tana = (*ihMin)->getP().x()/(*ihMin)->getP().z();
        double xEntWin = xDetUS + tana * ( dets->zEntrWind() - zDetUS );
        double tanb = (*ihMin)->getP().y()/(*ihMin)->getP().z();
        double yEntWin = yDetUS + tanb * ( dets->zEntrWind() - zDetUS );
        double zEntWin = dets->zEntrWind();
        Hep3Vector vDcInPart( (*ihMin)->getP() );

//----- extrapolate forward to the RICH entrance window
//      dummy (unit) cov matrix
        CsHelix helix = CsHelix( xDetUS, yDetUS, zDetUS,
			 	 tana, tanb, part->getCharge()/momPart, cov );
	CsHelix hxEntWin;
	if( helix.Extrapolate( zEntrWind, hxEntWin ) ) {
	  //cout << momPart << endl;
	  //cout << xDetUS << "  " << yDetUS << "  " << zDetUS << endl;
	  //cout << xEntWin << "  " << yEntWin << "  " << zEntWin << endl;
	  xEntWin = hxEntWin.getX();
	  yEntWin = hxEntWin.getY();
          //cout << xEntWin << "  " << yEntWin << "  " << zEntWin << endl;
          vDcInPart.setX( hxEntWin.getDXDZ() );
          vDcInPart.setY( hxEntWin.getDYDZ() );
          vDcInPart.setZ( 1. );
          vDcInPart = vDcInPart.unit();
	}
        Hep3Vector vPoInPart( xEntWin, yEntWin, zEntWin );

        xh = vPoInPart.x();
        yh = vPoInPart.y();
        if( hist.hRC1201 ) hist.hRC1201->Fill( xh, yh );
//                         ----------------------------

        //-- project forward to the RICH exit window (test)
	//        double xExtWinP = xDetUS + tana * ( zExitWind - zDetUS );
	//        double yExtWinP = yDetUS + tanb * ( zExitWind - zDetUS );
	//        double zExtWinP = zExitWind;

        ihMin = mchits.end();
        zDiff = 0.;
        zDiffMin = 100000.;
        for( ih=mchits.begin(); ih!=mchits.end(); ih++ ) {
          CsDetector* det = dynamic_cast<CsDetector*>((*ih)->getDet());
          if( det == NULL ) continue;
          zDiff = (*ih)->getZ() - zExitWind;
          if( zDiff > 0. ) {
            if( zDiff < zDiffMin ) {
	      zDiffMin = zDiff;
	      ihMin = ih;
	    }
	  }
	}
        Hep3Vector vPoExPart( 1000000., 1000000., 1000000. );
        if( ihMin != mchits.end() ) {
          if( zDiffMin < 1000. ) {                // provisional !!!
            double xDetDS = (*ihMin)->getX();
            double yDetDS = (*ihMin)->getY();
            double zDetDS = (*ihMin)->getZ();

//--------- project backward to the RICH exit window
//          assume no mag field
            double tana = (*ihMin)->getP().x()/(*ihMin)->getP().z();
            double xExtWin = xDetDS + tana * (zExitWind - zDetDS);
            double tanb = (*ihMin)->getP().y()/(*ihMin)->getP().z();
            double yExtWin = yDetDS + tanb * (zExitWind - zDetDS);
            double zExtWin = zExitWind;
	    vPoExPart.setX( xExtWin );
	    vPoExPart.setY( yExtWin );
	    vPoExPart.setZ( zExtWin );

            xh = vPoExPart.x();
            yh = vPoExPart.y();
            if( hist.hRC1202 ) hist.hRC1202->Fill( xh, yh );
//                             ----------------------------
	    //cout << vPoExPart << endl;
	    //cout << momPart << endl;
	  }
	}

//----- particle forward extrapolation to prov. photon 'emission point'
        Hep3Vector vPoEmssPt( 1000000., 1000000., 1000000. );
        Hep3Vector vDcEmssPt( 1000000., 1000000., 1000000. );
        CsHelix hxEmssPt;
        double zPoPhot0 = zEntrWind + partPathFr *(zExitWind - zEntrWind);
       	if( helix.Extrapolate( zPoPhot0, hxEmssPt ) ) {
          vPoEmssPt.setX( hxEmssPt.getX() );
          vPoEmssPt.setY( hxEmssPt.getY() );
          vPoEmssPt.setZ( hxEmssPt.getZ() );
          vDcEmssPt.setX( hxEmssPt.getDXDZ() );
          vDcEmssPt.setY( hxEmssPt.getDYDZ() );
          vDcEmssPt.setZ( 1. );
          vDcEmssPt = vDcEmssPt.unit();
        }

//----- extrapolate forward to the RICH exit window
//      dummy (unit) cov matrix
        Hep3Vector vPoExPartEx( 1000000., 1000000., 1000000. );
        Hep3Vector vDcExPartEx( 1000000., 1000000., 1000000. );
        CsHelix hxExtWin;
        if( helix.Extrapolate( zExitWind, hxExtWin ) ) {
          vPoExPartEx.setX( hxExtWin.getX() );
          vPoExPartEx.setY( hxExtWin.getY() );
          vPoExPartEx.setZ( hxExtWin.getZ() );
          vDcExPartEx.setX( hxExtWin.getDXDZ() );
          vDcExPartEx.setY( hxExtWin.getDYDZ() );
          vDcExPartEx.setZ( 1. );
          vDcExPartEx = vDcExPartEx.unit();
        }

  	int iTrack = (*it)->getGnum();

        int kPart = lParticles_.size();
        lParticles_.push_back( new CsRCParticle( kPart, (*it),
			       vPoInPart, vDcInPart, vPoExPart,
			       vPoEmssPt, vDcEmssPt,
			       vPoExPartEx, vDcExPartEx,
                               momPart, nPhoCer, thetaCer ) );
        flagSimul_ = true;

        if( CsRCExeKeys::Instance()->kPrintEventParts() == 2 ) {
	  lParticles_.back()->print();
    	  lParticles_.back()->printMC();
	}
      }

    }


//==========================================================================
  void CsRCEventParticles::setMCRecParticle() {
//---------------------------------------------

//--- Paolo  -  November 2000
//    WARNING : to be upgraded to setDataPaticle standard ! (5/2007)


      CsRCExeKeys* key = CsRCExeKeys::Instance();
      CsRCRecConst *cons = CsRCRecConst::Instance();
      CsRCDetectors *dets = CsRCDetectors::Instance();
      CsRCHistos& hist = CsRCHistos::Ref();
      int kh;
      double xh, yh, wh;

      static float partPathFr = 0.;
      static double zEntrWind = 0.;
      static double zExitWind = 0.;
      static double zPoPhot0 = 0.;
      list<CsRCMirrorNom*> lMirrNom = CsRCMirrors::Instance()->lMirrNom();
      list<CsRCPhotonDet*> lPhoDet = dets->lPhoDet();

      static float momMinAcc = 0.;
      static double cov[15];

      static bool firstCall = true;
      if( firstCall ) {
        firstCall = false;
	key->acknoMethod( "CsRCEventParticles::setMCRecParticle" );

	partPathFr = cons->partPathFr();
        zEntrWind = dets->zEntrWind();
        zExitWind = dets->zExitWind();
	double zMirrMin = lMirrNom.front()->vC0().z() +
	  lMirrNom.front()->RR() * cos( lPhoDet.front()->detAng() );
        zPoPhot0 = zEntrWind + partPathFr *(zMirrMin - zEntrWind);

	momMinAcc = cons->momMinAcc();

	for( int k=0; k<15; k++ ) cov[k] = 0.;
	int j = 1;
	for( int k=0; k<15; k+=j ) { cov[k] = 1.; j++; }
      }

      Coral* coral = Coral::Instance();
      list<CsMCTrack*> mctracks = coral->getMCTracks();
      CsRICH1Detector* rich = CsGeom::Instance()->getRich1Detector();
      CsEvent* event = coral->getEvent();
      list<CsTrack*> tracks = event->getTracks();
      list<CsTrack*>::iterator it;

//--- loop on tracks :
//    ----------------
      //cout << endl;
      //cout << "tracks : " << tracks.size() << endl;
      //cout << "MCtracks : " << mctracks.size() << endl;
      for( it=tracks.begin(); it!=tracks.end(); it++ ) {

        vector<CsHelix> vHelix = (*it)->getHelices();

	if( vHelix.size() == 0 ) continue;
//@@-------------------------------------
	double zHelix0 = vHelix[0].getZ();
	double zHelix1 = 0.;
	//if( vHelix.size() == 2 ) zHelix1 = vHelix[1].getZ();
        //---  030409
	//if( vHelix.size() >= 2 ) zHelix1 = vHelix[1].getZ();
	//---  030619
	if( vHelix.size() >= 2 ) zHelix1 = vHelix[vHelix.size()-1].getZ();

	//if( !(zHelix0 < 3500. && zHelix1 > 3500.) ) continue;  // MC01
        //if( !(zHelix0 > 3500. && zHelix1 > 18000.) ) continue;

        if( !(zHelix0 > 0. && zHelix1 > 3500.) ) continue;       // MC67
//@@-----------------------------------------------------

//----- low momentum particle cut
	double cop = vHelix[0].getCop();
	double momPart = 1.e+09;
	if( cop != 0. ) momPart = 1./fabs( cop );
        if( momPart > 300. ) continue;
        if( momPart < momMinAcc ) continue;
//@@--------------------------------------
	//cout << "***  " << zHelix0 << "  " << zHelix1 << "  " 
	//     << momPart << endl;
        //if( !(zHelix0 > 0. && zHelix1 > 3500.) ) {
	//if( vHelix.size() <= 2 ) {
	//  for( unsigned int k=0; k<vHelix.size(); k++ ) 
	//cout << vHelix[k].getZ() << "  * ";
	//cout << momPart << endl;
	//}

//----- get MC theta-Cherenkov and n-photons
        int nPhoCer = 0;
        double thetaCer = 0.;

	list<CsCluster*> lClusters = (*it)->getClusters();
	int nClus = lClusters.size();
	double chiSq = (*it)->getChi2();

	int nHits = 0;
//----- get MC associate track
	CsMCTrack* mctrack = CsMCUtils::getAssociatedMCTrack( (*it), nHits );
        if( !mctrack ) continue;
//@@---------------------------

        CsMCParticle* part = mctrack->getParticle();
	double charge = part->getCharge();
	if( charge == 0. ) continue;
//@@-------------------------------
	bool bOrig = true;
        list<CsMCHit*> mchits = mctrack->getMCHits();
        list<CsMCHit*>::iterator ih;
        double momMC = (*mchits.begin())->getP().mag();
        //cout << endl;
        //cout << "MCpart = " << mctrack->getGnum() << "  " << momMC << "  "
        //     << charge << "  " << nHits << endl;
        list<CsMCHit*>::iterator ihMin = mchits.end();
        double zDiff = 0.;
        double zDiffMin = 100000.;
        for( ih=mchits.begin(); ih!=mchits.end(); ih++ ) {
          //cout << (*ih)->getZ() << "  ";
	  CsDetector* det =  dynamic_cast<CsDetector*>((*ih)->getDet());
	  if( det == dynamic_cast<CsDetector*>( rich ) ) {
	    nPhoCer++;
	    CsMCRICH1Hit* hit = dynamic_cast<CsMCRICH1Hit*>( (*ih) );
	    thetaCer += hit->getCherAngle();
	    if( hit->getOrigin() != 0 ) {
	      bOrig = false;
	      //cout << "===  " << hit->getOrigin() << endl;
	    }
	  }
	  else {
            if( det == NULL ) continue;
            zDiff = (*ih)->getZ() - zExitWind;
            if( zDiff > 0. ) {
              if( zDiff < zDiffMin ) {
	        zDiffMin = zDiff;
	        ihMin = ih;
	      }
	    }
	  }
	}
        //cout << endl;
	if( nPhoCer != 0 ) {
          thetaCer /= float( nPhoCer );
          thetaCer *= 1000.;
	}
	//cout << "mom = " << momPart << "   momMC = " << momMC << endl;
	//cout << "thC = " << thetaCer << "   nphoC = " << nPhoCer << endl;

	//-- skip particle under C threshold
        //if( nPhoCer == 0 ) continue;        //   030325   !!!
//@@-------------------------------
	if( !bOrig ) continue;
	//if( bOrig ) continue;               //   useful test!
//@@-------------------------

	float trkTime = - 9999.;                             //   030409
	if( (*it)->hasMeanTime() ) trkTime = (*it)->getMeanTime();
	//cout << trkTime << endl;
	xh = trkTime;
	//if( hRC101 ) hRC101->Fill( xh );
//                   ------------------
	if( trkTime == 9999 ) {
	  //cout << "***  " << zHelix0 << "  " << zHelix1 << "  "
	  //     << momPart << endl;
	}
	if( fabs( trkTime ) > 10.) continue;                 //   030409
//@@---------------------------------------

//----- particle helix extrapolation to 0 (Vadim cut)        //   030409
	CsHelix hxPTg;
	double xPTg = 0.;
	double yPTg = 0.;
	if( vHelix[0].Extrapolate( 0., hxPTg ) ) {
	  xPTg = hxPTg.getX();
	  yPTg = hxPTg.getY();
  	  //if( !(zHelix0 > 0. && zHelix1 > 3500.) ) {
	  //cout << zHelix0 << "  " << zHelix1 << "  " << momPart << "  "
	  //     << xPTg << "  " << yPTg << endl;
	  //}
	  if( fabs( xPTg ) > 500.  ||  fabs( yPTg ) > 300. ) continue;
//@@-----------------------------------------------------------------
	}

//----- use MC associate track
//      ----------------------
/*	bool entwMC = false;
        Hep3Vector vDcInPartMC( 0. );
	if( nHits >= 50 && ihMin != mchits.end() && zDiffMin < 1000. ) {
//------- project forward to the RICH entrance window
          double xDetUS = (*ihMin)->getX();
          double yDetUS = (*ihMin)->getY();
          double zDetUS = (*ihMin)->getZ();
          double tana = (*ihMin)->getP().x()/(*ihMin)->getP().z();
          double xEntWin = xDetUS + tana * ( dets->zEntrWind() - zDetUS );
          double tanb = (*ihMin)->getP().y()/(*ihMin)->getP().z();
          double yEntWin = yDetUS + tanb * ( dets->zEntrWind() - zDetUS );
          double zEntWin = zEntrWind;
	  vDcInPartMC.setX( (*ihMin)->getP().x() );
	  vDcInPartMC.setY( (*ihMin)->getP().y() );
	  vDcInPartMC.setZ( (*ihMin)->getP().z() );
          //-- extrapolate forward to the RICH entrance window
          //   dummy (unit) cov matrix
          CsHelix helix = CsHelix( xDetUS, yDetUS, zDetUS,
			  	 tana, tanb, part->getCharge()/momPart, cov );
	  CsHelix hxEntWin;
	  if( helix.Extrapolate( zEntrWind, hxEntWin ) ) {
	    //cout << momPart << endl;
	    //cout << xDetUS << "  " << yDetUS << "  " << zDetUS << endl;
	    //cout << xEntWin << "  " << yEntWin << "  " << zEntWin << endl;
	    //cout << hxEntWin.getX() << "  " << hxEntWin.getY() << "  "
	    //     << hxEntWin.getZ() << endl;
	    xEntWin = hxEntWin.getX();
	    yEntWin = hxEntWin.getY();
            vDcInPartMC.setX( hxEntWin.getDXDZ() );
            vDcInPartMC.setY( hxEntWin.getDYDZ() );
            vDcInPartMC.setZ( 1. );
            vDcInPartMC = vDcInPartMC.unit();
	  }
          Hep3Vector vPoInPartMC( xEntWin, yEntWin, zEntWin );
	  entwMC = true;
	  //cout << nHits << " ===== vDcInPartMC = " << vDcInPartMC << endl;
	}
*/

//----- particle at RICH entrance window
//      helix extrapolation (use mag.field)
	bool entw = false;
	CsHelix hxEntWin;
	int kex = 0;
	if( zHelix1 < zEntrWind ) kex = 1;

	if( vHelix[kex].Extrapolate( zEntrWind, hxEntWin ) ) {

          Hep3Vector vPoInPart( hxEntWin.getX(), hxEntWin.getY(),
	  		        hxEntWin.getZ() );

          xh = vPoInPart.x();
          yh = vPoInPart.y();
          if( hist.hRC1201 ) hist.hRC1201->Fill( xh, yh );
//                           ----------------------------

	  //!!!	  Hep3Vector vDcInPart( vDcInPartMC );
	  Hep3Vector vDcInPart( hxEntWin.getDXDZ(), hxEntWin.getDYDZ(), 1. );
          vDcInPart = vDcInPart.unit();

	  //cout << " ===== vPoInPart = " << vPoInPart << endl;
	  //cout << " ===== vDcInPart = " << vDcInPart << endl;
	  entw = true;

	  bool monHist = true;
//@@-------------------------
//------- resolution monitoring histo.s
	  if( monHist ) {
            list<CsMCHit*>::iterator ihXMin = mchits.end();
            double zXDiff = 0.;
            double zXDiffMin = 100000.;
            for( ih=mchits.begin(); ih!=mchits.end(); ih++ ) {
	      CsDetector* det =  dynamic_cast<CsDetector*>((*ih)->getDet());
	      if( det != dynamic_cast<CsDetector*>( rich ) ) {
                zXDiff = (*ih)->getZ() - zHelix0;
                if( fabs( zXDiff ) < zXDiffMin ) {
	          zXDiffMin = fabs( zXDiff );
	          ihXMin = ih;
		}
	      }
	    }
    	    //cout << (*ihXMin)->getZ() << endl;
  	    if( ihXMin != mchits.end() ) {
	      CsHelix hxHit0;
	      int kex = 0;
	      if( vHelix[kex].Extrapolate( (*ihXMin)->getZ(), hxHit0 ) ) {
	        double *covMx = hxHit0.getCov();
	        //for( int j=0; j<15; j++ ) cout << covMx[j] << "  ";
	        //cout << endl;
       	        xh = hxHit0.getX() - (*ihXMin)->getX();
	        xh /= sqrt( covMx[0] );
	        yh = hxHit0.getY() - (*ihXMin)->getY();
	        yh /= sqrt( covMx[2] );
                if( hist.hRC1204 ) hist.hRC1204->Fill( xh, yh );
//                                 ----------------------------
  	        xh = hxHit0.getDXDZ() -
	          (*ihXMin)->getP().x()/(*ihXMin)->getP().z();
	        xh /= sqrt( covMx[5] );
	        yh = hxHit0.getDYDZ() -
	          (*ihXMin)->getP().y()/(*ihXMin)->getP().z();
	        yh /= sqrt( covMx[9] );
                if( hist.hRC1205 ) hist.hRC1205->Fill( xh, yh );
//                                 ----------------------------
	        xh = hxHit0.getDXDZ() -
	          (*ihXMin)->getP().x()/(*ihXMin)->getP().z();
	        xh *= 1000.;
	        yh = hxHit0.getDYDZ() -
	          (*ihXMin)->getP().y()/(*ihXMin)->getP().z();
	        yh *= 1000.;
                if( hist.hRC1208 ) hist.hRC1208->Fill( xh, yh );
//                                 ----------------------------
	      }
	    }
            zXDiff = 0.;
            zXDiffMin = 100000.;
            for( ih=mchits.begin(); ih!=mchits.end(); ih++ ) {
	      CsDetector* det =  dynamic_cast<CsDetector*>((*ih)->getDet());
	      if( det != dynamic_cast<CsDetector*>( rich ) ) {
	        //cout << (*ih)->getZ() << "  ";
                zXDiff = (*ih)->getZ() - zEntrWind;
                if( fabs( zXDiff ) < zXDiffMin ) {
	          zXDiffMin = fabs( zXDiff );
	          ihXMin = ih;
		}
	      }
	    }
	    //cout << "   " << (*ihXMin)->getZ() << endl;
	    if( ihXMin != mchits.end() ) {
  	      CsHelix hxHitEw;
	      int kex = 0;
	      if( zHelix1 < zEntrWind ) kex = 1;
	      if( vHelix[kex].Extrapolate( (*ihXMin)->getZ(), hxHitEw ) ) {
	        double *covMx = hxHitEw.getCov();
	        xh = hxHitEw.getX() - (*ihXMin)->getX();
	        xh /= sqrt( covMx[0] );
	        yh = hxHitEw.getY() - (*ihXMin)->getY();
	        yh /= sqrt( covMx[2] );
                if( hist.hRC1206 ) hist.hRC1206->Fill( xh, yh );
//                                 ----------------------------
	        xh = hxHitEw.getDXDZ() -
	          (*ihXMin)->getP().x()/(*ihXMin)->getP().z();
	        xh /= sqrt( covMx[5] );
	        yh = hxHitEw.getDYDZ() -
	          (*ihXMin)->getP().y()/(*ihXMin)->getP().z();
	        yh /= sqrt( covMx[9] );
                if( hist.hRC1207 ) hist.hRC1207->Fill( xh, yh );
//                                 ----------------------------
	        xh = hxHitEw.getDXDZ() -
	          (*ihXMin)->getP().x()/(*ihXMin)->getP().z();
	        xh *= 1000.;
	        yh = hxHitEw.getDYDZ() -
	          (*ihXMin)->getP().y()/(*ihXMin)->getP().z();
	        yh *= 1000.;
                if( hist.hRC1209 ) hist.hRC1209->Fill( xh, yh );
//                                 ----------------------------
	      }
	    }
	  }

	  bool doExtrap = true;
//@@--------------------------
//------- particle to provisional photon 'emission point'
          Hep3Vector vPoEmssPt( 1000000., 1000000., 1000000. );
          Hep3Vector vDcEmssPt( 1000000., 1000000., 1000000. );
//------- particle straight projection from entrance window (no field)
          double xEmss = vPoInPart.x() +
	    vDcInPart.x()/vDcInPart.z() * (zPoPhot0 - zEntrWind);
          double yEmss = vPoInPart.y() +
	    vDcInPart.y()/vDcInPart.z() * (zPoPhot0 - zEntrWind);
          double zEmss = vPoInPart.z();
	  vPoEmssPt.setX( xEmss );
	  vPoEmssPt.setY( yEmss );
          vPoEmssPt.setZ( zEmss );
          vDcEmssPt = vDcInPart;
          CsHelix hxEmssPt;
	  if( doExtrap ) {
//------- particle forward helix extrapolation (use mag.field)
	    int kox = 0;
	    if( zHelix1 < zEntrWind ) kox = 1;
       	    if( vHelix[kox].Extrapolate( zPoPhot0, hxEmssPt ) ) {
	      vPoEmssPt.setX( hxEmssPt.getX() );
	      vPoEmssPt.setY( hxEmssPt.getY() );
	      vPoEmssPt.setZ( hxEmssPt.getZ() );
	      vDcEmssPt.setX( hxEmssPt.getDXDZ() );
              vDcEmssPt.setY( hxEmssPt.getDYDZ() );
              vDcEmssPt.setZ( 1. );
	      vDcEmssPt = vDcEmssPt.unit();
	    }
	  }

//------- particle to RICH exit window
          Hep3Vector vPoExPartEx( 1000000., 1000000., 1000000. );
          Hep3Vector vDcExPartEx( 1000000., 1000000., 1000000. );
//------- particle straight projection from entrance window (no field)
          double xExit = vPoInPart.x() +
	    vDcInPart.x()/vDcInPart.z() * (zExitWind - zEntrWind);
          double yExit = vPoInPart.y() +
	    vDcInPart.y()/vDcInPart.z() * (zExitWind - zEntrWind);
	  vPoExPartEx.setX( xExit );
	  vPoExPartEx.setY( yExit );
          vPoExPartEx.setZ( zExitWind );
          vDcExPartEx = vDcInPart;
          CsHelix hxExtWin;
	  if( doExtrap ) {
//------- particle for/backward helix extrapolation at RICH exit window
	    int kox = 0;
	    if( zHelix1 > 3500. && zHelix1 < 15000. ) kox = 1;
       	    if( vHelix[kox].Extrapolate( zExitWind, hxExtWin ) ) {
	      vPoExPartEx.setX( hxExtWin.getX() );
              vPoExPartEx.setY( hxExtWin.getY() );
              vPoExPartEx.setZ( hxExtWin.getZ() );
	      vDcExPartEx.setX( hxExtWin.getDXDZ() );
              vDcExPartEx.setY( hxExtWin.getDYDZ() );
              vDcExPartEx.setZ( 1. );
	      vDcExPartEx = vDcExPartEx.unit();
	    }
	  }
          xh = vPoExPartEx.x();
          yh = vPoExPartEx.y();
          if( hist.hRC1202 ) hist.hRC1202->Fill( xh, yh );
//                           ----------------------------
	  //cout << " ===== vPoExPartEx = "<< vPoExPartEX << endl;

//------- particle backward projection to point at RICH exit window
	  vector<string> detNamePo;
	  //detNamePo.push_back( "AAS1" );
	  //detNamePo.push_back( "AAS2" );
	  //detNamePo.push_back( "AAS3" );
	  //detNamePo.push_back( "AAS4" );
	  detNamePo.push_back( "STV1" );
	  detNamePo.push_back( "STY1" );
	  detNamePo.push_back( "STX1" );
	  detNamePo.push_back( "STV2" );
	  detNamePo.push_back( "STY2" );
	  detNamePo.push_back( "STX2" );
	  detNamePo.push_back( "STV3" );
	  detNamePo.push_back( "STY3" );
	  detNamePo.push_back( "STX3" );
	  detNamePo.push_back( "STV4" );
	  detNamePo.push_back( "STY4" );
	  detNamePo.push_back( "STX4" );
	  detNamePo.push_back( "STV5" );
	  detNamePo.push_back( "STY5" );
	  detNamePo.push_back( "STX5" );
	  detNamePo.push_back( "STV6" );
	  detNamePo.push_back( "STY6" );
	  detNamePo.push_back( "STX6" );
	  //int unitPo = 1;
	  int unitPo = 7;
          Hep3Vector vPoExPart( 1000000., 1000000., 1000000. );
	  const char *detName;
  	  list<CsCluster*> lClusters = (*it)->getClusters();
	  list<CsCluster*>::iterator ic;
	  std::vector<double> cClu, zClu;
	  std::vector<double> cosDet, sinDet;
	  for( ic=lClusters.begin(); ic!=lClusters.end(); ic++ ) {
	    list<CsDetector*> detes = (*ic)->getDetsList();
	    if( detes.size() > 1 ) continue;
	    detName = detes.front()->getName();
	    int unit = detes.front()->getUnit();
            //for( unsigned int kDet=0; kDet<4; kDet++ ) {
	    for( unsigned int kDet=0; kDet<detNamePo.size(); kDet++ ) {
	      if( detName != detNamePo[kDet] ) continue;
	      if( unit != unitPo ) continue;
	      //cout << detName << "  " << unit << "  " << (*ic)->getW()
	      //     << "  " << (*ic)->getU() << endl;
	      //if( !((*ic)->getW() >  9200. && (*ic)->getW() <  9300.) )
	      //continue;
	      //if( !((*ic)->getW() > 10100. && (*ic)->getW() < 10300.) )
	      //continue;
	      cosDet.push_back(detes.front()->getCosAng());
	      sinDet.push_back(detes.front()->getSinAng());
	      //cout << detName << "  " << unit << "  "
	      //     << cosDet[kClu] << "  " << sinDet[kClu] << "  " << endl;
	      cClu.push_back((*ic)->getU());
	      zClu.push_back((*ic)->getW() - zExitWind);
	    }
	  }
	  int nClu = cClu.size();
	  //if( nClu > 0 ) {
	  //cout << endl;
	  //for( int k=0; k<nClu; k++ )
	  //  cout << cClu[k] << "  " << zClu[k] << "  " << cosDet[k] << endl;
	  //cout << endl;
	  //}
//------- use PCPOINT algo
	  double tana = hxExtWin.getDXDZ();
	  double tanb = hxExtWin.getDYDZ();
	  HepVector XX( 2, 0 );
	  if( nClu >= 3 ) {
	    //if( nClu >= 3 && nClu <=4 ) {
	    HepMatrix GX( 2, 2, 0 ), EX( 2, 2, 0 );
	    HepVector YM( nClu, 0 );
	    HepMatrix AA( nClu, 2, 0 );
	    for( int k=0; k<nClu; k++) {
	      YM[k] = cClu[k] - 
		(tana*cosDet[k] + tanb*sinDet[k]) * zClu[k];
	      AA[k][0] = cosDet[k];
	      AA[k][1] = sinDet[k];
	    }
//--------- assume unit (virtual) error matrix.
	    GX = AA.T() * AA;
	    int ifl = 0;
	    EX = GX.inverse( ifl );
	    if( ifl == 0 ) XX = EX * AA.T() * YM;
	    vPoExPart.setX( XX[0] );
	    vPoExPart.setY( XX[1] );
	    vPoExPart.setZ( zExitWind );
	    //cout << "  vPoExPartBkw = "<< vPoExPart << endl;
	  }
	  //cout << "  vPoExPartBkw = "<< vPoExPart << endl;

//------- construct particle object

          int kPart = lParticles_.size();
	  lParticles_.push_back( new CsRCParticle( kPart, (*it), mctrack,
		  	         vPoInPart, vDcInPart, vPoExPart,
				 vPoEmssPt, vDcEmssPt,
				 vPoExPartEx, vDcExPartEx,
                                 momPart, nPhoCer, thetaCer ) );
          flagSimul_ = true;

	  lParticles_.back()->setQualityPars( zHelix0, zHelix1, nClus,
					      chiSq );

	  lParticles_.back()->setMCmTime( trkTime );      //   040127
	  //std::cout << lParticles_.back()->MCmTime() << std::endl;

          if( CsRCExeKeys::Instance()->kPrintEventParts() == 2 ) {
	    lParticles_.back()->print();
	    lParticles_.back()->printMC();
	  }

	}

      }

    }


//==========================================================================
  void CsRCEventParticles::setDataParticle() {
//--------------------------------------------

//--- Paolo  -  November 2000


      CsRCExeKeys* key = CsRCExeKeys::Instance();
      CsRCRecConst *cons = CsRCRecConst::Instance();
      CsRCDetectors *dets = CsRCDetectors::Instance();
      CsRCHistos& hist = CsRCHistos::Ref();
      int kh;
      double xh, yh, wh;

      static float partPathFr = 0.;
      static double zEntrWind = 0.;
      static double zExitWind = 0.;
      static double zPoPhot0 = 0.;
      list<CsRCMirrorNom*> lMirrNom = CsRCMirrors::Instance()->lMirrNom();
      list<CsRCPhotonDet*> lPhoDet = dets->lPhoDet();

      static float momMinAcc = 0.;
      static float momMaxAcc = 0.;

      static bool firstCall = true;
      if( firstCall ) {
        firstCall = false;
	key->acknoMethod( "CsRCEventParticles::setDataParticle" );

	partPathFr = cons->partPathFr();
        zEntrWind = dets->zEntrWind();
        zExitWind = dets->zExitWind();
	double zMirrMin = lMirrNom.front()->vC0().z() +
	  lMirrNom.front()->RR() * cos( lPhoDet.front()->detAng() );
        zPoPhot0 = zEntrWind + partPathFr *(zMirrMin - zEntrWind);
	momMinAcc = cons->momMinAcc();
	momMaxAcc = cons->momMaxAcc();
      }

      xh = 0.5;
      if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                       ------------------------

      Coral* coral = Coral::Instance();
      CsEvent* event = coral->getEvent();
      list<CsTrack*> tracks = event->getTracks();
      list<CsTrack*>::iterator it;

//--- loop on tracks :
//    ----------------
      int kTrack = -1;
      for( it=tracks.begin(); it!=tracks.end(); it++ ) {

	kTrack++;

	xh = 1.5;
	if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                         ------------------------

	//-- track chi-square
	double chiSq = (*it)->getChi2();

	//-- helices
        vector<CsHelix> vHelix = (*it)->getHelices();
	if( vHelix.size() == 0 ) {
	  if( key->kPrintRejections() == 1 ) {
            std::cout << "RICHONE - setDataParticle : Ev "
                      << CsRichOne::Instance()->kEvent() << "  Track "
	       	      << kTrack << ",  NO Helices"
		      << std::endl;
	  }
	  xh = 2.5;
	  if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                           ------------------------
	  continue;
//@@--------------
	}
	double zHelix0 = vHelix[0].getZ();
	double zHelix1 = 0.;
	//if( vHelix.size() == 2 ) zHelix1 = vHelix[1].getZ();
	//---  020909
        //if( vHelix.size() >= 2 ) zHelix1 = vHelix[1].getZ();
        //---  030619
	if( vHelix.size() >= 2 ) zHelix1 = vHelix[vHelix.size()-1].getZ();

	//cout << zHelix0 << "  " << zHelix1 << endl;
	//for(int k=0; k<vHelix.size(); k++ ) cout << vHelix[k].getZ()
	//                                         << "  ";
	//cout << endl;

        //if( !(zHelix0 < 3500. && zHelix1 > 3500.) ) continue;
        //if( !(zHelix0 >  9000. && zHelix0 < 18000.) ) continue;
        //if( !(zHelix1 > 18000. && zHelix1 < 35000.) ) continue;
//@@------------------------------------------------------------

	double zTarget = CsGeom::Instance()->getTargetCenter();
        //-- for DST processing (Vinicio)
	if( CsInit::Instance()->getDataType() == "dst" ) {    // 060519
	  //if( !(zHelix0 > 0.) ) {                           // 090204
	  if( !(zHelix0 > zTarget ) ) {
	    if( key->kPrintRejections() == 1 ) {
              std::cout << "RICHONE - setDataParticle : Ev "
                        << CsRichOne::Instance()->kEvent() << "  Track "
	       	        << kTrack << ",  bad z-helix "
			<< zHelix0 << std::endl;
	    }
	    xh = 3.5;
	    if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                             ------------------------
            continue;                                         // 060519
//@@----------------
	  }
	} else {                                              // 060519
	  //if( !(zHelix0 > 0.  &&  zHelix1 > 3500.) ) {      // 090204
	  if( !(zHelix0 > zTarget  &&  zHelix1 > 3500.) ) {
	    //-- do not count beam tracks
	    //if( zHelix0 > 0.  &&  zHelix1 > 0. ) {          // 090204
	    if( zHelix0 > zTarget  &&  zHelix1 > zTarget ) {
	      if( key->kPrintRejections() == 1 ) {
                std::cout << "RICHONE - setDataParticle : Ev "
                          << CsRichOne::Instance()->kEvent() << "  Track "
	         	  << kTrack << ",  bad z-helix "
			  << zHelix0 << "  " << zHelix1 << "  "
			  << 1./vHelix[0].getCop() << std::endl;
	      } 
	      xh = 4.5;
	      if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                               ------------------------
	    } else {
	      xh = 5.5;
	      if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                               ------------------------
	    }
            continue;                                         // 060519
//@@----------------
	  }
	}                                                     // 060519

	int SM12 = -1;
        if( zHelix0 < 3500.  &&  zHelix1 < 18000. ) SM12 = 1;
        if( zHelix0 > 3500.  &&  zHelix1 > 18000. ) SM12 = 2;
        if( zHelix0 < 3500.  &&  zHelix1 > 18000. ) SM12 = 3;

//----- particle momentum :
//      -------------------
        double cop = vHelix[0].getCop();
        double momPart = 1.e+09;
	if( cop != 0. ) {
	  momPart = 1./fabs( cop );
	} else {
	  //^if( momPart > 300. ) {
	  if( key->kPrintRejections() == 1 ) {
            std::cout << "RICHONE - setDataParticle : Ev "
                      << CsRichOne::Instance()->kEvent() << "  Track "
	              << kTrack << "  NO mom ("
		      << momPart << ")" << std::endl;
	  }
	  xh = 6.5;
	  if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                           ------------------------
	  continue;
//@@--------------
	}
        if( momPart < momMinAcc ) {
	  if( key->kPrintRejections() == 1 ) {
            std::cout << "RICHONE - setDataParticle : Ev "
                      << CsRichOne::Instance()->kEvent() << "  Track "
	              << kTrack << "  mom "
		      << momPart << ", low momentum"
		      << std::endl;
	  }
	  xh = 7.5;
	  if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                           ------------------------
	  continue;
//@@--------------
	}
        if( momPart > momMaxAcc ) {
	  if( key->kPrintRejections() == 1 ) {
            std::cout << "RICHONE - setDataParticle : Ev "
                      << CsRichOne::Instance()->kEvent() << "  Track "
	              << kTrack << "  mom "
		      << momPart << ", high momentum"
		      << std::endl;
	  } 
	  xh = 8.5;
	  if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                           ------------------------
	  continue;
//@@--------------
	}
//cout << "***  " << zHelix0 << "  " << zHelix1 << "  " << momPart << endl;

	if (hist.hRC1220)                             // HISTO MOMENTUM VS. TIME
	  hist.hRC1220->Fill((*it)->hasMeanTime()?(*it)->getMeanTime():1000,
			     momPart );

	vector<CsHelix>::iterator ix;
	//for( ix=vHelix.begin(); ix!=vHelix.end(); ix++ ) {
	//cout << "zHelix = " << (*ix).getZ() << endl;
	//}
	list<CsCluster*> lClusters = (*it)->getClusters();
	//list<CsCluster*>::iterator ic;
	//for( ic=lClusters.begin(); ic!=lClusters.end(); ic++ ) {
	//  cout << (*ic)->getW() << "  ";
	//}
	//cout << endl;
	int nClus = lClusters.size();


//----- particle helix extrapolation to PT :
//      ------------------------------------
	CsHelix hxPTg;
	double xPTg = 0.;
	double yPTg = 0.;
	//if( vHelix[0].Extrapolate( -350., hxPTg ) ) {         // 090204
	if( vHelix[0].Extrapolate( zTarget, hxPTg ) ) {
	  xPTg = hxPTg.getX();
	  yPTg = hxPTg.getY();
	  //cout << hxPTg.getX() << "  " << hxPTg.getY() << endl; 
	}

//----- particle at RICH entrance window :
//      ----------------------------------
//      helix extrapolation (use mag.field)
	CsHelix hxEntWin;
	int kex = 0;
        if( vHelix.size() > 1 ) {                          // 060519
	  if( zHelix1 < zEntrWind ) kex = 1;
	}                                                  // 060519
	if( vHelix[kex].Extrapolate( zEntrWind, hxEntWin ) ) {

          Hep3Vector vPoInPart( hxEntWin.getX(), hxEntWin.getY(),
	  		        hxEntWin.getZ() );

//------- entrance window size - as from Silvia  (021107)
          if( fabs( vPoInPart.x() ) >  1660. ) {
	    if( key->kPrintRejections() == 1 ) {
              std::cout << "RICHONE - setDataParticle : Ev "
                        << CsRichOne::Instance()->kEvent() << "  Track "
	         	<< kTrack << "  mom "
			<< momPart << ",  Out of Entr. Wind (X), "
			<< vPoInPart.x() << std::endl;
	    }
	    xh = 9.5;
	    if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                             ------------------------
	    continue;
//@@----------------
	  }
          if( fabs( vPoInPart.y() ) >  1205. ) {
	    if( key->kPrintRejections() == 1 ) {
              std::cout << "RICHONE - setDataParticle : Ev "
                        << CsRichOne::Instance()->kEvent() << "  Track "
	         	<< kTrack << "  mom "
			<< momPart << ",  Out of Entr. Wind. (Y), "
			<< vPoInPart.y() << std::endl;
	    }
	    xh = 10.5;
	    if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                             ------------------------
	    continue;
//@@----------------
	  }

	  int charge = 0;
	  if( hxEntWin.getCop() > 0. )  charge =  1;
	  if( hxEntWin.getCop() < 0. )  charge = -1;
	  if( charge == 0 ) {
	    if( key->kPrintRejections() == 1 ) {
              std::cout << "RICHONE - setDataParticle : Ev "
                        << CsRichOne::Instance()->kEvent() << "  Track "
	         	<< kTrack << "  mom "
			<< momPart << ",  charge = 0" << std::endl;
	    }
	    xh = 11.5;
	    if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                             ------------------------
	    continue;
//@@----------------
	  }
          Hep3Vector vDcInPart( hxEntWin.getDXDZ(), hxEntWin.getDYDZ(), 1.);
          vDcInPart = vDcInPart.unit();

//------- set particle direction parallel to z-axis (test):
	  //vDcInPart.setX( 0.);  vDcInPart.setY( 0.);  vDcInPart.setZ( 1.);
	  //cout << " ===== vPoInPart = " << vPoInPart << endl;
	  //cout << " ===== vDcInPart = " << vDcInPart << endl;

//------- geometrical cuts (NOT ACTIVE!)
	  bool bPT = true;
          double xyPT = ( pow( xPTg/ 100., 2 ) +
			  pow( yPTg/ 100., 2 ) );
	  bPT = xyPT > 1.;
	  //if( !bPT ) continue;
	  bool ppoo = true;
          double ppxy = ( pow( vPoInPart.x()/ 200., 2 ) +
			  pow( vPoInPart.y()/ 200., 2 ) );
	  ppoo = ppxy > 1.;
	  //if( !ppoo ) continue;
	  bool ttgg = true;
	  double ttxy = ( pow( vDcInPart.x()/vDcInPart.z()/ 0.010, 2 ) +
			  pow( vDcInPart.y()/vDcInPart.z()/ 0.010, 2 ) );
	  ttgg = ttxy > 1.;
	  //if( !ttgg ) continue;

          xh = xPTg;
          yh = yPTg;
          if( hist.hRC1200 ) hist.hRC1200->Fill( xh, yh );
//                           ----------------------------
          xh = vPoInPart.x();
          yh = vPoInPart.y();
          if( hist.hRC1201 ) hist.hRC1201->Fill( xh, yh );
//                           ----------------------------
          xh = chiSq;
          //xh = chiSq / float( nClus );
          yh = float( SM12 ) - 0.5;
          if( hist.hRC1203 ) hist.hRC1203->Fill( xh, yh );
//                           ----------------------------


//------- particle (provisional) photon 'emission point' :
//        ------------------------------------------------
          Hep3Vector vPoEmssPt( 1000000., 1000000., 1000000. );
          Hep3Vector vDcEmssPt( 1000000., 1000000., 1000000. );
//------- use trajectory smoothing at emission point (zPoPhot0) :
//        -------------------------------------------------------
//        if helix available
	  CsHelix hxEmssPt;
	  bool smoo = false;
	  if( vHelix.size() > 2 ) {                     //  020904
	    vector<CsHelix>::iterator ixp = vHelix.begin();
	    ixp++;
	    vector<CsHelix>::iterator ixm = vHelix.end();
	    ixm--;
	    for( ix=ixp; ix!=ixm; ix++ ) {
	      if( fabs( (*ix).getZ() - zPoPhot0 ) < 500. ) {
		//cout << "zHelix em = " << (*ix).getZ() << endl;
                hxEmssPt = (*ix);
		smoo = true;
		break;
	      }
	    }
	  }
	  if( smoo ) {
	    vPoEmssPt.setX( hxEmssPt.getX() );
	    vPoEmssPt.setY( hxEmssPt.getY() );
	    vPoEmssPt.setZ( hxEmssPt.getZ() );
            vDcEmssPt.setX( hxEmssPt.getDXDZ() );
            vDcEmssPt.setY( hxEmssPt.getDYDZ() );
            vDcEmssPt.setZ( 1. );
            vDcEmssPt = vDcEmssPt.unit();
	  }                                             //  020904
	  else {
//------- particle straight projection from entrance window (no field)
            double xEmss = vPoInPart.x() +
	      vDcInPart.x()/vDcInPart.z() * (zPoPhot0 - zEntrWind);
            double yEmss = vPoInPart.y() +
	      vDcInPart.y()/vDcInPart.z() * (zPoPhot0 - zEntrWind);
            double zEmss = vPoInPart.z();
	    vPoEmssPt.setX( xEmss );
	    vPoEmssPt.setY( yEmss );
            vPoEmssPt.setZ( zEmss );
            vDcEmssPt = vDcInPart;
 	    bool doEnExtrap = true;
//          ----------------------
//--------- particle forward extrapolation (use mag. field)
	    if( doEnExtrap ) {
	      int kox = 0;
	      if( vHelix.size() > 1 ) {                  // 061018
	        if( zHelix1 < zEntrWind ) kox = 1;
	      }
       	      if( vHelix[kox].Extrapolate( zPoPhot0, hxEmssPt ) ) {
	        vPoEmssPt.setX( hxEmssPt.getX() );
	        vPoEmssPt.setY( hxEmssPt.getY() );
   	        vPoEmssPt.setZ( hxEmssPt.getZ() );
	        vDcEmssPt.setX( hxEmssPt.getDXDZ() );
                vDcEmssPt.setY( hxEmssPt.getDYDZ() );
                vDcEmssPt.setZ( 1. );
	        vDcEmssPt = vDcEmssPt.unit();
	      }
	    }
	    //cout << "Extrap  " << zPoPhot0 << "  " << vDcEmssPt << endl;
	  }

	  /* //^
//------- use trajectory smoothing at emission point (zPoPhot0)
//        if helix available
	  if( vHelix.size() > 2 ) {                     //  020904
	    vector<CsHelix>::iterator ixp = vHelix.begin();
	    ixp++;
	    vector<CsHelix>::iterator ixm = vHelix.end();
	    ixm--;
	    for( ix=ixp; ix!=ixm; ix++ ) {
	      if( fabs( (*ix).getZ() - zPoPhot0 ) < 500. ) {
		//cout << "zHelix em = " << (*ix).getZ() << endl;
                hxEmssPt = (*ix);
		break;
	      }
	    }
	    vPoEmssPt.setX( hxEmssPt.getX() );
	    vPoEmssPt.setY( hxEmssPt.getY() );
	    vPoEmssPt.setZ( hxEmssPt.getZ() );
            vDcEmssPt.setX( hxEmssPt.getDXDZ() );
            vDcEmssPt.setY( hxEmssPt.getDYDZ() );
            vDcEmssPt.setZ( 1. );
            vDcEmssPt = vDcEmssPt.unit();
	  }                                              //  020904
	  //cout << "Smooth  " << vDcEmssPt << endl;
	  */ //^


//------- particle to RICH exit window :
//        ------------------------------
          Hep3Vector vPoExPartEx( 1000000., 1000000., 1000000. );
          Hep3Vector vDcExPartEx( 1000000., 1000000., 1000000. );
//------- particle straight projection from entrance window (no field)
          double xExit = vPoInPart.x() +
	    vDcInPart.x()/vDcInPart.z() * (zExitWind - zEntrWind);
          double yExit = vPoInPart.y() +
	    vDcInPart.y()/vDcInPart.z() * (zExitWind - zEntrWind);
	  vPoExPartEx.setX( xExit );
	  vPoExPartEx.setY( yExit );
          vPoExPartEx.setZ( zExitWind );
          vDcExPartEx = vDcInPart;
 	  bool doExExtrap = true;
//        -----------------------
//------- particle forw/backward extrapolation to RICH exit window
          CsHelix hxExtWin;
	  if( doExExtrap ) {
	    int kox = 0;
	    if( vHelix.size() > 1 ) {                     //  061018
  	      if( zHelix1 > 3500. && zHelix1 < 15000. ) kox = 1;
	    }
       	    if( vHelix[kox].Extrapolate( zExitWind, hxExtWin ) ) {
	      vPoExPartEx.setX( hxExtWin.getX() );
              vPoExPartEx.setY( hxExtWin.getY() );
              vPoExPartEx.setZ( hxExtWin.getZ() );
	      vDcExPartEx.setX( hxExtWin.getDXDZ() );
              vDcExPartEx.setY( hxExtWin.getDYDZ() );
              vDcExPartEx.setZ( 1. );
	      vDcExPartEx = vDcExPartEx.unit();
	    }
	  }
          xh = vPoExPartEx.x();
          yh = vPoExPartEx.y();
          if( hist.hRC1202 ) hist.hRC1202->Fill( xh, yh );
//                           ----------------------------
	  //cout << " ===== vPoExPartEx = " << vPoExPartEx << endl;
	  //cout << "-  " << vPoInPart << "  " << vPoEmssPt
	  //     << "  " << vPoExPartEx << endl;
	  //cout << "   " << vDcInPart << "  " << vDcEmssPt
	  //     << "  " << vDcExPartEx << endl;


//------- particle backward projection to a point at RICH exit window :
//        -------------------------------------------------------------
          Hep3Vector vPoExPart( 1000000., 1000000., 1000000. );
	  bool doBkExtrap = true;
//        ----------------------
	  if( doBkExtrap ) {
  	    vector<string> detNamePo;
//          Use MWPCs
  	    //detNamePo.push_back( "AAS1" );
	    //detNamePo.push_back( "AAS2" );
	    //detNamePo.push_back( "AAS3" );
	    //detNamePo.push_back( "AAS4" );
	    //int unitPo = 1;
//          Use STRAWs
	    //detNamePo.push_back( "STV1" );
	    //detNamePo.push_back( "STY1" );
	    //detNamePo.push_back( "STX1" );
	    //detNamePo.push_back( "STV2" );
	    //detNamePo.push_back( "STY2" );
	    //detNamePo.push_back( "STX2" );
	    //detNamePo.push_back( "STV3" );
	    //detNamePo.push_back( "STY3" );
	    //detNamePo.push_back( "STX3" );
	    //detNamePo.push_back( "STV4" );
	    //detNamePo.push_back( "STY4" );
	    //detNamePo.push_back( "STX4" );
	    //detNamePo.push_back( "STV5" );
	    //detNamePo.push_back( "STY5" );
	    //detNamePo.push_back( "STX5" );
	    //detNamePo.push_back( "STV6" );
	    //detNamePo.push_back( "STY6" );
	    //detNamePo.push_back( "STX6" );
	    //int unitPo = 7;
//          Use RICHWALL (2006>)
	    detNamePo.push_back( "RWX1" );
	    detNamePo.push_back( "RWY1" );
	    detNamePo.push_back( "RWX2" );
	    detNamePo.push_back( "RWY2" );
            //int unitPo = 1;
	    if( !detNamePo.empty() ) {
	      const char *detName;
      	      list<CsCluster*> lClusters = (*it)->getClusters();
	      list<CsCluster*>::iterator ic;
	      std::vector<double> cClu, zClu;
	      std::vector<double> cosDet, sinDet;
	      for( ic=lClusters.begin(); ic!=lClusters.end(); ic++ ) {
	        list<CsDetector*> lDetectors = (*ic)->getDetsList();
	        if( lDetectors.size() > 1 ) continue;
	        detName = lDetectors.front()->getName();
	        int unit = lDetectors.front()->getUnit();
	        for( unsigned int kDet=0; kDet<detNamePo.size(); kDet++ ) {
	          if( detName != detNamePo[kDet] ) continue;
	          //^if( unit != unitPo ) continue;
	          //cout << detName << "  " << unit << "  " << (*ic)->getW()
	          //     << "  " << (*ic)->getU() << endl;
	          cosDet.push_back(lDetectors.front()->getCosAng());
	          sinDet.push_back(lDetectors.front()->getSinAng());
	          //cout << detName << "  " << unit << "  "
	          //<< cosDet[kClu] << "  " << sinDet[kClu] << "  " << endl;
	          cClu.push_back((*ic)->getU());
	          zClu.push_back((*ic)->getW() - zExitWind);
	        }
    	      }
 	      int nClu = cClu.size();
	      //if( nClu > 0 ) {
              //cout << endl;
	      //cout << nClu << endl;
	      //for( int k=0; k<nClu; k++ )
	      //  cout << cClu[k] << "  " << zClu[k]
	      //<< "  " << cosDet[k] << endl;
	      //}
//----------- use PCPOINT algo
	      double tana = hxExtWin.getDXDZ();
	      double tanb = hxExtWin.getDYDZ();
	      HepVector XX( 2, 0 );
	      if( nClu >= 3 ) {
	        HepMatrix GX( 2, 2, 0 ), EX( 2, 2, 0 );
	        HepVector YM( nClu, 0 );
	        HepMatrix AA( nClu, 2, 0 );
	        for( int k=0; k<nClu; k++) {
	          YM[k] = cClu[k] - 
	    	    (tana*cosDet[k] + tanb*sinDet[k]) * zClu[k];
	          AA[k][0] = cosDet[k];
	          AA[k][1] = sinDet[k];
	        }
//------------- assume unit (virtual) error matrix.
	        GX = AA.T() * AA;
	        int ifl = 0;
	        EX = GX.inverse( ifl );
	        if( ifl == 0 ) XX = EX * AA.T() * YM;
	        vPoExPart.setX( XX[0] );
	        vPoExPart.setY( XX[1] );
	        vPoExPart.setZ( zExitWind );
		//cout << "  vPoExPartBkw = " << vPoExPart 
		//     << "  " << vPoExPartEx << endl;
	      }
	    }
	  }

//------- set geometrical corrections
	  bool bCorr = true;
	  if( bCorr) {
	    setGeoCorrections( vPoInPart, vDcInPart );
	    setGeoCorrections( vPoEmssPt, vDcEmssPt );
//          -----------------------------------------
	  }

/*
//%%%%---------------------------------------
	  int kox = 0;
	  std::cout << setprecision( 8 );
          if( vHelix.size() > 1 ) {
	    if( zHelix1 > 3500. && zHelix1 < 15000. ) kox = 1;
          }
	  std::cout << std::endl;
	  std::cout << "tarck-in     = " << vHelix[kox].getX() << "  "
	            << vHelix[kox].getY() << "  " << vHelix[kox].getZ()
	            << "  " << vHelix[kox].getDXDZ() << "  "
	            << vHelix[kox].getDYDZ() << "  "
	            << 1./vHelix[kox].getCop() << std::endl;
	  //set infinite momentum
	  double mx;
	  TMtx Mv(5);
	  TMtx Mm(5,5);
	  const CsHelix cHx = vHelix[kox];
	  THlx tHx; 
	  tHx.ImportHelix( cHx );
	  tHx.Get( mx, Mv, Mm );
	  //for(int k=1;k<=5;k++) std::cout << Mv(k) << " ";
	  //std::cout << std::endl;
	  Mv(5) = 0.;
	  Mv(5) = 1./0.01;
	  tHx.Set( mx, Mv, Mm );
	  //for(int k=1;k<=5;k++) std::cout << Mv(k) << " ";
	  //std::cout << std::endl;
	  //vHelix[kox] = tHx.ExportHelix();
	  std::cout << std::endl;
	  std::cout << "mom set " << 1./vHelix[kox].getCop() << std::endl;
	  //
//@@-----------------------------------
	  zExitWind = 6000.;
	  Hep3Vector vPoCoo( 0.), vDcCoo( 0.);
	  if( vHelix[kox].Extrapolate( zExitWind, hxExtWin ) ) {
	    std::cout << std::endl;
	    std::cout << "track-extrap = " << hxExtWin.getX() << "  "
		      << hxExtWin.getY() << "  " << hxExtWin.getZ()
		      << "  " << hxExtWin.getDXDZ() << "  "
		      << hxExtWin.getDYDZ() << "  "
		      << 1./hxExtWin.getCop() << std::endl;
	  }
	  std::cout << std::endl;
	  vPoCoo.setX( vHelix[kox].getX() );
	  vPoCoo.setY( vHelix[kox].getY() );
	  vPoCoo.setZ( vHelix[kox].getZ() );
	  vDcCoo.setX( vHelix[kox].getDXDZ() );
	  vDcCoo.setY( vHelix[kox].getDYDZ() );
	  vDcCoo.setZ( 1. );
	  vDcCoo = vDcCoo.unit();
	  double chgOvMom = 0.0003 * vHelix[kox].getCop();
	  if( CsRCTrackMomFit::Instance()->test( vPoCoo, vDcCoo, chgOvMom,
						 zExitWind ) ) {
  	    std::cout << std::endl;
	    std::cout << "track-trackd = " << vPoCoo.x() << "  "
	              << vPoCoo.y() << "  " << vPoCoo.z() << "  "
	              << vDcCoo.x()/vDcCoo.z() << "  "
	              << vDcCoo.y()/vDcCoo.z() << std::endl;
	  }
	  exit(0);
//%%%%-------------------------------------------------------------------
*/

//------- construct particle object

          int kPart = lParticles_.size();
          lParticles_.push_back( new CsRCParticle( kPart, (*it),
		  	         vPoInPart, vDcInPart, vPoExPart,
				 vPoEmssPt, vDcEmssPt,
				 vPoExPartEx, vDcExPartEx,
				 momPart, charge ) );
          flagSimul_ = false;

	  lParticles_.back()->setQualityPars( zHelix0, zHelix1, nClus,
					      chiSq );

	  xh = 13.5;
	  if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                           ------------------------

          if( CsRCExeKeys::Instance()->kPrintEventParts() == 2 ) {
	    lParticles_.back()->print();
	  }

	} else {

	  xh = 12.5;
	  if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                           ------------------------

	  if( key->kPrintRejections() == 1 ) {
            std::cout << "RICHONE - setDataParticle : Ev "
                      << CsRichOne::Instance()->kEvent() << "  Track "
	              << kTrack << "  mom "
		      << momPart << ",  NO extrap to Entr Wind"
		      << std::endl;
	  }
	}

      }

      if( lParticles_.empty() ) {
	xh = 14.5;
	if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                         ------------------------
      }

      xh = 19.5;
      if( hist.hRC1000 ) hist.hRC1000->Fill( xh );
//                       ------------------------

      return;

    }


//==========================================================================
  void CsRCEventParticles::setGeoCorrections( Hep3Vector& vPo, 
//---------------------------------------------------------------
					      Hep3Vector& vDc ) {

//--- corrections due to RICH positioning
//    -----------------------------------
//--- Paolo  -  November 2002


      CsOpt* opt = CsOpt::Instance();
      bool boo;
      vector<float> vflo;

      static Hep3Vector vCorrPo( 0.);
      static Hep3Vector vCorrDc( 0.);

      static bool firstCall = true;
      if( firstCall ) {
        firstCall = false;

	CsRCExeKeys* key = CsRCExeKeys::Instance();
        key->acknoMethod( "CsRCEventParticles::setGeoCorrections" );

//----- read corrections, provisional!
	boo = opt->CsOpt::getOpt( "RICHONE", "partPosCorrs", vflo );
	if( boo ) {
	  vCorrPo.setX( vflo[0] );
	  vCorrPo.setY( vflo[1] );
	  vCorrPo.setZ( vflo[2] );
	  cout << " RICHONE, CsRCEventParticles::setGeoCorrections() :";
	  cout << " corrections to particle position  " << vCorrPo << endl;
          cout << "--------------------------------------------------------"
	       << "----------------------------" << endl;
	}
	boo = opt->CsOpt::getOpt( "RICHONE", "partDirCorrs", vflo );
	if( boo ) {
	  vCorrDc.setX( vflo[0] );
	  vCorrDc.setY( vflo[1] );
	  vCorrDc.setZ( vflo[2] );
	  cout << " RICHONE, CsRCEventParticles::setGeoCorrections() :";
	  cout << " corrections to particle direction  " << vCorrDc << endl;
          cout << "--------------------------------------------------------"
	       << "----------------------------" << endl;
	}

      }

//--- set corrections, provisional!
      vPo += vCorrPo;
      vDc += vCorrDc;
      vDc = vDc.unit();

  }


//==========================================================================
  void CsRCEventParticles::partAnalysis() {
//-----------------------------------------


//--- particle selection in input to ring recognition
//    -----------------------------------------------
//--- Paolo  -  December 1999, rev. October 2000


//--- from "CsRCExecKeys"
      CsRCExeKeys *key = CsRCExeKeys::Instance();
      bool PartFilt = key->PartFilt();
      bool ExitWind = key->ExitWind();
      bool PartSele = key->PartSele();
      bool PartCorr = key->PartCorr();
      bool PartReso = key->PartReso();

//--- from "CsRCRecConst"
      CsRCRecConst *cons = CsRCRecConst::Instance();

      CsRCHistos& hist = CsRCHistos::Ref();
      int kh;
      double xh, yh, wh;

      bool MCarloEvent = key->MCarloEvent();
      bool DataEvent = key->DataEvent();

//--- loop over particles :
//    ---------------------
      //list<CsRCParticle*> lParticles = lParticles();
      list<CsRCParticle*>::iterator ip;
      for( ip=lParticles_.begin(); ip!=lParticles_.end(); ip++ ) {

//----- monitor input particles (MC) :
//      ------------------------------
        if( MCarloEvent ) if( (*ip)->flag() ) {
	  int iPartT = (*ip)->iPartT();
	  if( iPartT > 200) iPartT -= 200;
	  if( hist.hRC3615 ) hist.hRC3615->Fill( double( iPartT ), 0.5 );
//hh                         -------------------------------------------
	}

//----- filter particle in :   ( not anymore useful... )
//      --------------------
        if( PartFilt ) {
	  if( !partFilter( (*ip) ) ) (*ip)->setFlag( false );
//            --------------------
        }

//----- select particle type ( MC ) :
//      -----------------------------
//      WARNING : referenced also from MCMonitor !
        if( PartSele  &&  MCarloEvent ) {
          if( !partSelect( (*ip) ) ) (*ip)->setFlag( false );
//            --------------------
	}

//----- check particle deflection at exit window :
//      ------------------------------------------
        if( ExitWind ) {
	  if( !exitWindow( (*ip) ) ) (*ip)->setFlag( false );
//            --------------------
        }

//----- correct particle direction for m.c.s. :
//      ---------------------------------------
        if( PartCorr ) {
// !!!    partCorr( (*ip) );   //  moved to CsRCPartPhotons::mcsPartCorr()
//        -----------------
	}

//-==-- simulate particle experimental resolution :
//      -------------------------------------------
        if( PartReso  &&  MCarloEvent )  partResol( (*ip) );
//                                       ------------------
      }


//-------------------------------------------------------------------
      bool exeTest = false;
      //bool exeTest = true;
//@@----------------------
      if( DataEvent  &&  exeTest ) {
//----- checks for 'correlated' background : (050225)
//      ------------------------------------
        //cout << lParticles_.size() << endl;
	std::list<CsRCParticle*> lParticles_cp;
	for( ip=lParticles_cp.begin(); ip!=lParticles_cp.end(); ip++ )
	  delete (*ip);
	lParticles_cp.clear();
	for( ip=lParticles_.begin(); ip!=lParticles_.end(); ip++ ) {
	  CsRCParticle* part = (*ip);
	  if( key->myFileType() >= 1 ) {
	    int kPart = part->kPart();
	    kPart += 100;
	    CsTrack* pTrack = part->pTrack();
	    Hep3Vector vPosIn = part->vPosIn();
	    vPosIn.setX( -vPosIn.x() );
	    vPosIn.setY( -vPosIn.y() );
	    Hep3Vector vDirIn = part->vDirIn();
	    vDirIn.setX( -vDirIn.x() );
	    vDirIn.setY( -vDirIn.y() );
	    Hep3Vector vPoExPart( 0. );
	    Hep3Vector vPosEmP = part->vPosEmP();
	    vPosEmP.setX( -vPosEmP.x() );
	    vPosEmP.setY( -vPosEmP.y() );
	    Hep3Vector vDirEmP = part->vDirEmP();
	    vDirEmP.setX( -vDirEmP.x() ) ;
	    vDirEmP.setY( -vDirEmP.y() ) ;
	    Hep3Vector vPosExW( 0. );
	    Hep3Vector vDirExW( 0. );
	    double momPart = part->mom();
	    int charge = part->charge();
	    lParticles_cp.push_back( new CsRCParticle( kPart, NULL,
				     vPosIn, vDirIn, vPoExPart,
			  	     vPosEmP, vDirEmP, vPosExW, vDirExW,
				     momPart, charge ) );
	    lParticles_cp.back()->setQualityPars( 0., 0., 0, 0. );
	  }
        }
	for( ip=lParticles_cp.begin(); ip!=lParticles_cp.end(); ip++ ) {
	  lParticles_.push_back( (*ip) );
	}
      }
      //cout << lParticles_.size() << endl;
//-------------------------------------------------------------------

      if( MCarloEvent ) {
        for( ip=lParticles_.begin(); ip!=lParticles_.end(); ip++ ) {
          if( !(*ip)->flag() ) {
            xh = (*ip)->iPartT();
            if( hist.hRC3617 ) hist.hRC3617->Fill( xh );
//hh                           ------------------------
          }
        }
      }

      int nPartEv = 0;
      for( ip=lParticles_.begin(); ip!=lParticles_.end(); ip++ ) {
	if( (*ip)->flag() ) nPartEv++;
      }
      xh = nPartEv;
      if( hist.hRC1533 ) hist.hRC1533->Fill( xh );
//hh                     ------------------------

    }


//==========================================================================
  bool CsRCEventParticles::partFilter( CsRCParticle* part ) {
//-----------------------------------------------------------


//--- December  1999, rev. October 2000
//--- Revised   September 2002


      CsRCExeKeys *key = CsRCExeKeys::Instance();
      CsRCRecConst *cons = CsRCRecConst::Instance();

      CsRCMirrors *mirr = CsRCMirrors::Instance();
      static std::list<CsRCMirrorElem*> lMirrEleProc;
      std::list<CsRCMirrorElem*>::iterator im;

      static bool firstCall = true;
      if( firstCall ) {
        firstCall = false;
        key->acknoMethod( "CsRCEventParticles::ParticleFilter" );

	std::list<CsRCMirrorElem*>lMirrEle = mirr->lMirrEle();
	std::list<string> lMirr;
	CsOpt* opt = CsOpt::Instance();
        bool boo;
        boo = opt->CsOpt::getOpt( "RICHONE", "MirrorEleProcName", lMirr );
        if( !boo ) lMirr.clear();
	//lMirr.push_back( "MB42" );
	if( !lMirr.empty() ) {
          std::list<string>::iterator is;
          for( is=lMirr.begin(); is!=lMirr.end(); is++ ) {
	    for( im=lMirrEle.begin(); im!=lMirrEle.end(); im++ ) {
	      if( (*is) == (*im)->name() ) {
	        lMirrEleProc.push_back( (*im) );
	        //std::cout << lMirrEleProc.back()->name() << "  ";
	        break;
	      }
	    }
	  }
	} else {
	  lMirrEleProc.clear();
	}
	if( !lMirrEleProc.empty() ) {
	  std::cout << " RICHONE :  CsRCEventParticles::ParticleFilter "
		    << "process mirror elements  ";
	  for( im=lMirrEleProc.begin(); im!=lMirrEleProc.end(); im++ ) {
	    //std::cout << lMirrEleProc.back()->name() << "  ";
	    std::cout << (*im)->name() << "  ";
	  }
	  std::cout << "only (MIRROR CENTRE!)" << std::endl;
	}

      }

      bool flag = true;

      list<CsRCParticle*>::iterator ip;


//--- check data flag :
//    -----------------
//--- data particle selection
      if( key->DataEvent() ) {

	if( part->mom() < cons->momMinAcc() ) flag = false;

	if( part->mom() > cons->momMaxAcc() ) flag = false;

	if( !flag ) return  flag;
//------------------------------

//----- process selected mirror elements (centre) only :
//      ------------------------------------------------
	float rrElmQ = 10000.;                   //   100.**2
//@@---------------------------------------------------------
	Hep3Vector vPoPaMir(0., 0., 0.);
	flag = false;
        for( im=lMirrEleProc.begin(); im!=lMirrEleProc.end(); im++ ) {
	  Hep3Vector vCePos = (*im)->vpos();
//------- approx. particle impact on mirrors :
	  Hep3Vector vPosIn0 = part->vPosIn();
	  Hep3Vector vDirIn0 = part->vDirIn();
	  float zMir = 8700.;
          float yMir = vPosIn0.y() + vDirIn0.y() * (zMir - vPosIn0.z());
	  std::list<CsRCMirrorNom*>::iterator in;
	  std::list<CsRCMirrorNom*> lMirrNom = mirr->lMirrNom();
	  if( yMir >= 0.) in = lMirrNom.begin();
	  if( yMir <  0.) { in = lMirrNom.begin(); in++; }
	  Hep3Vector vPoC0 = (*in)->vC0();
	  double RR = (*in)->RR();
	  Hep3Vector vPoPaMir0 = mirr->vImpMir( vPosIn0, vDirIn0, vPoC0, RR );
//                               --------------------------------------------
	  if ( (vPoPaMir0 - vCePos).mag2() < rrElmQ ) {
	    flag = true;
	    break;
	  }
	}
	if( lMirrEleProc.empty() ) flag = true;

        if( !flag ) return  flag;
//------------------------------
      }

//--- check simulation flag :
//    -----------------------
      if( key->MCarloEvent() ) {

//----- check two particle momenta 'equal'(!) (18/11/99) :
//      --------------------------------------------------
//----- obsolete ( 020415 )
        double momC = part->mom();
        for( ip=lParticles_.begin(); ip!=lParticles_.end(); ip++ ) {
          if( (*ip) == part ) continue;
          if( momC == (*ip)->mom() ) {
            flag = false;  ;
            break;
          }
        }   /* end loop on particles: ip */
        if( !flag ) return  flag;
//------------------------------

        for( ip=lParticles_.begin(); ip!=lParticles_.end(); ip++ ) {

//$$          if( (*ip)->ioEx() != 2 )  {
//$$            bool flag = false;
//$$            (*ip)->setFlag( flag );
//$$          }
//@@----------------------------------
        }   /* end loop on particles: ip */

        for( ip=lParticles_.begin(); ip!=lParticles_.end(); ip++ ) {

//$$          bool flag = false;
//$$          (*ip)->setFlag( flag );
//  //@@----------------------------
//$$          if( (*ip)->iPartT() > 200 )  {
//$$            if( (*ip)->ioEx() == 0 )  {
//$$              bool flag = true;
//$$              (*ip)->setFlag( flag );
//  //@@--------------------------------
//$$            }
//$$          }

        }   /* end loop on particles: ip  */
      }   /* end if MCarloEvent */

      return flag;

  }


//========================================================================
  bool CsRCEventParticles::partSelect( CsRCParticle* part ) {
//-----------------------------------------------------------


//    December  1999,  rev.  October 2000


      CsRCExeKeys *key = CsRCExeKeys::Instance();
 
      static bool firstCall = true;
      if( firstCall ) {
        firstCall = false;

	key->acknoMethod( "CsRCEventParticles::ParticleSelection" );
      }

      bool bAccp = false;
//@@--------------------

      if( part->flag() ) {

        int iPartT = part->iPartT();
        if( iPartT > 200) iPartT -= 200;
//@@------------------------------------
	int charge = part->charge();

//----- select particle interacting in RICH :
//      -------------------------------------
//$$        if( iPartT > 200) bAccp = true;
//$$//@@----------------------------------

//----- select all :
//      ------------
//$$        if( iPartT >= 1 && iPartT <= 15) bAccp = true;
//$$//@@-------------------------------------------------

//----- select electrons :
//      ------------------
        //if( iPartT ==  2 || iPartT ==  3) bAccp = true;
//@@--------------------------------------------------

//----- select muons :
//      --------------
	//if( iPartT ==  5 || iPartT ==  6) bAccp = true;
//@@--------------------------------------------------

//----- Select pions :
//      --------------
	//if( iPartT ==  8 || iPartT ==  9) bAccp = true;
//@@--------------------------------------------------

//----- select kaons :
//      --------------
	//if( iPartT == 11 || iPartT == 12) bAccp = true;
//@@--------------------------------------------------

//----- select protons :
//      ----------------
	//if( iPartT == 14 || iPartT == 15) bAccp = true;
//@@--------------------------------------------------


//----- select pions under thresh :
//      ---------------------------
	if( iPartT ==  8 || iPartT ==  9 ) {
	  double mass = CsRCRecConst::Instance()->massPartv()[iPartT];
	  double CFRefInd = CsRCRecConst::Instance()->CFRefInd();
	  double momThresh = mass / sqrt( CFRefInd*CFRefInd - 1.);
	  //if( part->mom() < momThresh ) bAccp = true;
//@@------------------------------------------------
	}
//----- select kaons under thresh :
//      ---------------------------
	if( iPartT == 11 || iPartT == 12 ) {
	  double mass = CsRCRecConst::Instance()->massPartv()[iPartT];
	  double CFRefInd = -CsRCRecConst::Instance()->CFRefInd();
	  double momThresh = mass / sqrt( CFRefInd*CFRefInd - 1.);
	  //if( part->mom() < momThresh ) bAccp = true;
//@@------------------------------------------------
	}
//----- select protons under thresh :
//      -----------------------------
	if( iPartT == 14 || iPartT == 15 ) {
	  double mass = CsRCRecConst::Instance()->massPartv()[iPartT];
	  double CFRefInd = CsRCRecConst::Instance()->CFRefInd();
	  double momThresh = mass / sqrt( CFRefInd*CFRefInd - 1.);
	  //if( part->mom() < momThresh ) bAccp = true;
//@@------------------------------------------------
	}
	//if( bAccp ) cout << iPartT << "  " << part->mom() << endl;

//----- select particle charge :
//      ------------------------
//$$        if( charge != + 1) bAccp = false;
//$$        if( charge != - 1) bAccp = false;
//$$//@@-                      -------------

      }

//--- OVERWRITE SELECTION !
//    ---------------------
      bAccp = true;

      return bAccp;

  }


//===========================================================================
  bool CsRCEventParticles::exitWindow( CsRCParticle* part ) {
//-----------------------------------------------------------


//--- check particle deflection at exit window :
//    ------------------------------------------
//--- Paolo  -  December 1999,  rev.  October 2000


      CsRCExeKeys *key = CsRCExeKeys::Instance();

      CsRCHistos& hist = CsRCHistos::Ref();
      int kh;
      double xh, yh, wh;
 
      static double cov[15];
      static bool firstCall = true;
      if( firstCall ) {
        firstCall = false;
        key->acknoMethod( "CsRCEventParticles::ExitWindowTest" );

	for( int k=0; k<15; k++ ) cov[k] = 0.;
	int j = 1;
	for( int k=0; k<15; k+=j ) { cov[k] = 1.; j++; }
      }

//--- from "CsRCRecConst"
      CsRCRecConst *cons = CsRCRecConst::Instance();

      bool flag = true;

      if( !part->flag() ) return  false;

      if( part->zo() < 1000000.) {

        double mom = part->mom();

//----- project particle to the exit window from entrance window :
//      ----------------------------------------------------------
        double tanA = part->ld() / part->nd();
        double tanB = part->md() / part->nd();
        double ddZZ = part->zo() - part->za();
        double xpaWind = part->xa() + tanA * ddZZ;
        double ypaWind = part->ya() + tanB * ddZZ;
//----- or, better
//----- use extrapolation forward to the RICH exit window
//      (dummy (unit) cov matrix)
//      -------------------------
	if( part->vPosExW().z() < 1000000.) {
	  xpaWind = part->vPosExW().x();
	  ypaWind = part->vPosExW().y();
	}
        double ddXX = part->xo() - xpaWind;
        double ddYY = part->yo() - ypaWind;

        double ddRRq = ddXX*ddXX + ddYY*ddYY;
	//double cirCutq = cirCut( mom );
	//cirCutq = cirCutq*cirCutq;
	//if( ddRRq > cirCutq ) flag = false;

	double cutXq = getCutXq( mom );
	double cutYq = getCutYq( mom );
	double ddRRqW = 1000000.;
	if( cutXq != 0.  &&  cutYq != 0. ) {
	  ddRRqW = ddXX*ddXX/cutXq + ddYY*ddYY/cutYq;
	  if( ddRRqW > 2.) flag = false;
//@@-----------------------------------
	}

//----- test (BUMP?) :
//      --------------
	//if( xpaWind < -450.  ||  xpaWind > -350. ) flag = false;
	//if( ypaWind < -100.  ||  ypaWind >  100. ) flag = false;
//@@-----------------------------------------------------------

	if( flag ) {
          xh = mom;
          yh = ddRRq;
          if( hist.hRC3605 ) hist.hRC3605->Fill( xh, yh );
//hh                         ----------------------------
          xh = xpaWind;
          yh = ypaWind;
          if( hist.hRC3636 ) hist.hRC3636->Fill( xh, yh );
//hh                         ----------------------------
	}
        xh = ddXX;
        yh = ddYY;
        if( hist.hRC3606 ) hist.hRC3606->Fill( xh, yh );
//hh                       ----------------------------
        xh = mom;
        yh = ddRRqW;
        if( hist.hRC3607 ) hist.hRC3607->Fill( xh, yh );
//hh                       ----------------------------
        xh = mom;
        yh = ddXX;
        if( hist.hRC3608 ) hist.hRC3608->Fill( xh, yh );
//hh                       ----------------------------
        xh = mom;
        yh = ddYY;
        if( hist.hRC3609 ) hist.hRC3609->Fill( xh, yh );
//hh                       ----------------------------
        xh = xpaWind;
        yh = ypaWind;
        if( hist.hRC3635 ) hist.hRC3635->Fill( xh, yh );
//hh                       ----------------------------

      } else {
//--- select particles with hit on tracker downsteam only :
//    -----------------------------------------------------
	//flag = false;
//@@----------------
      }

//--- swap function action :
//    ----------------------
      //bool Xflag = flag;             //   !!!
      //if( Xflag ) flag = false;      //   !!!
      //if( !Xflag ) flag = true;      //   !!!

      return flag;

  }


//===========================================================================
  double CsRCEventParticles::getCutXq( const double mom ) {
//---------------------------------------------------------


//    February 2001

//--- 3 sigma X-cut at the exit window, square
//    ----------------------------------------

      static float par0;
      static float par1;
      static float dd;

      static bool firstCall = true;
      if( firstCall ) {
	firstCall = false;

        CsOpt* opt = CsOpt::Instance();
        CsRCRecConst *cons = CsRCRecConst::Instance();

//----- defaults
//      from Vadim's 6*3k ref files/lepto_full ( co-18k.hist )
//      ------------------------------------------------------
//----- parameters from fit :
//      ---------------------
        par0 = 36.7;             //   8/02/01
        par1 =  0.89;            //   8/02/01
//----- ad hoc correction :
//      -------------------
        dd = 1.;                 //  15/02/01
//@@-----------

//----- from rich1.options :
//      --------------------
	bool boo;
        vector<float> vPar;
        boo = opt->CsOpt::getOpt( "RICHONE", "getCutXq", vPar );
        if( boo ) {
	  par0 = vPar[0];
	  par1 = vPar[1];
        }
        if( cons->printConsts() ) {
	  cout << " RICHONE, CsRCEventParticles::getCutXq :   " 
	       << par0 << "  " << par1 << endl;
	}

      }

//--- from exwind.kumac/circu2.f - histo.s 3608 & 3609 (+off):
//    --------------------------------------------------------

      double cutXq = 0.6*dd + par0 / pow( double(mom-dd), double(par1) );
//--------------------------------------------------------

      //cout << "cutXq  " << cutXq << endl;

      return  cutXq*cutXq;

  }

//===========================================================================
  double CsRCEventParticles::getCutYq( const double mom ) {
//---------------------------------------------------------


//    February 2001

//--- 3 sigma Y-cut at the exit window, square
//    ----------------------------------------

      static float par0;
      static float par1;
      static float dd;

      static bool firstCall = true;
      if( firstCall ) {
	firstCall = false;

        CsOpt* opt = CsOpt::Instance();
        CsRCRecConst *cons = CsRCRecConst::Instance();

//----- defaults
//      from Vadim's 6*3k ref files/lepto_full ( co-18k.hist )
//      ------------------------------------------------------
//----- parameters from fit :
//      ---------------------
        par0 = 27.9;             //   8/02/01
        par1 =  0.93;            //   8/02/01
//----- ad hoc correction :
//      -------------------
        dd = 1.;                 //  15/02/01
//@@-----------

//----- from rich1.options :
//      --------------------
	bool boo;
        vector<float> vPar;
        boo = opt->CsOpt::getOpt( "RICHONE", "getCutYq", vPar );
        if( boo ) {
	  par0 = vPar[0];
	  par1 = vPar[1];
        }
        if( cons->printConsts() ) {
	  cout << " RICHONE, CsRCEventParticles::getCutYq :   " 
	       << par0 << "  " << par1 << endl;
	}

      }

//--- from exwind.kumac/circu2.f - histo.s 3608 & 3609 (+off):
//    --------------------------------------------------------

      double cutYq = 0.6*dd + par0 / pow( double(mom-dd), double(par1) );
//--------------------------------------------------------

      //cout << "cutYq  " << cutYq << endl;

      return  cutYq*cutYq;

  }


//===========================================================================
  double CsRCEventParticles::cirCut( double mom )   {
//---------------------------------------------------


//--- interface function :
//    --------------------
//    December  1999


      double cirCut = 3.* 1.4142 * cirCu1( mom );
//----------------------------------------------
      return cirCut;

  }

//===========================================================================
  double CsRCEventParticles::cirCu1( double mom )   {
//---------------------------------------------------


//--- radius (1 sigma) of circular cut at the exit window
//    ---------------------------------------------------
//    from MC (Dima file rich1_50000_all.ntp) :
//    -----------------------------------------
//--- from circu.kumac :
//    ------------------
//    December  1999,  rev.  October 2000

 
      static float par[10];
//@@----------------------
      par[1] = 111.3;                    //   26/10/99
      par[2] =   0.112;                  //   26/10/99
//@@-------------------

      double cirCu1 = par[1] / sqrt( (mom+5.)*(mom+5.)*(mom+5.) ) +
//-----------------------------------------------------------------
                      par[2];
      return cirCu1;

  }


//===========================================================================
  void CsRCEventParticles::partCorr( CsRCParticle* part )   {
//-----------------------------------------------------------


//--- correct particle direction for m.c.s. :
//    ---------------------------------------
//--- Paolo  -  December 1999, rev.  October 2000

//--- OBSOLETE as from 19/02/01
//    moved to CsRCPartPhotons::mcsPartCorr()


      CsRCExeKeys *key = CsRCExeKeys::Instance();

      CsRCHistos& hist = CsRCHistos::Ref();
      int kh;
      double xh, yh, wh;
 
      static bool firstCall = true;
      if( firstCall ) {
        firstCall = false;
        key->acknoMethod( "CsRCEventParticles::ParticleCorrection" );
      }

//--- from "CsRCRecConst"
      CsRCRecConst *cons = CsRCRecConst::Instance();

      if( part->flag() ) {
	if( part->zo() < 1.e+06 ) {            // provisional!!!

//----- particle at the exit window vs entrance window :
//      ------------------------------------------------
        double tanA = part->ld() / part->nd();
        double tanB = part->md() / part->nd();
        double ddXX = part->xo() - part->xa();
        double ddYY = part->yo() - part->ya();
        double ddZZ = part->zo() - part->za();
        double tanAw = ddXX / ddZZ;
        double tanBw = ddYY / ddZZ;
        double norm = sqrt( 1.+ tanAw*tanAw + tanBw*tanBw );
        double ldw = tanAw / norm;
        double mdw = tanBw / norm;
        double ndw =    1. / norm;
        Hep3Vector vDirw( ldw, mdw, ndw );
	//!!!	part->setDirIn( vDirw );         //!!!

        ddXX = tanAw - tanA;
        ddYY = tanBw - tanB;
        xh = ddXX;
        yh = ddYY;
        if( hist.hRC3637 ) hist.hRC3637->Fill( xh, yh );
//hh                       ----------------------------

	}
      }

  }


//==========================================================================
  void CsRCEventParticles::partResol( CsRCParticle* part ) {
//----------------------------------------------------------


//--- simulate particle experimental resolution :
//    -------------------------------------------
//--- Paolo  -  December 1999


//--- from "CsRCRecConst"
      CsRCRecConst *cons = CsRCRecConst::Instance();


//pro-memoria .............................


  }


//==========================================================================
  bool CsRCEventParticles::initMyFile( int& nPartEv ) {
//-----------------------------------------------------


//--- Paolo  -  November 2001
//    rev.      March 2003


      CsRCExeKeys *key = CsRCExeKeys::Instance();

      static int kMyFile;
      static int nMyFiles;
      static list<string> RmyFileName;
      static list<string>::iterator is;

      int bytePerInt = sizeof(int);
      static int fh;
      int iret;
      int pSw;

      int ibuff[10];
      char *fileName;

      static bool firstCall = true;
      if( firstCall ) {
	firstCall = false;

	kMyFile = 0;
	RmyFileName = key->RmyFileName();
	nMyFiles = RmyFileName.size();

	is = RmyFileName.begin();

	do {

	  if( kMyFile >= nMyFiles ) {
	    key->setEndWrMyFile();
//          ---------------------
	    return  false;
	  }

	  fileName = const_cast<char*>((*is).c_str());
	  fh = openMyFile( fileName );
//        ---------------------------
	  kMyFile++;
	  is++;

	} while( fh < 0 );

	myFile_ = fh;

	setRunIndex( fileName );
//      -----------------------

	return  true;

      }   /* end if firstCall */

      close( myFile_ );               //   030318
//    ----------------

      do {

	if( kMyFile >= nMyFiles ) {
	  key->setEndWrMyFile();
//        ---------------------
	  return false;
	}

	fileName = const_cast<char*>((*is).c_str());
	fh = openMyFile( fileName );
//      ---------------------------
	kMyFile++;
	is++;

      } while( fh < 0 );

      myFile_ = fh;

      setRunIndex( fileName );
//    -----------------------

//--- skip myFile record header
      ///pSw = 3 * bytePerInt;
      ///iret = read( fh, ibuff, pSw );
      ///nPartEv = ibuff[2];
      int nSkip = 1;
      if( key->readGf5() ) nSkip = 2;
      pSw = nSkip * bytePerInt;
      iret = read( fh, ibuff, pSw );
// --------------------------------
      //std::cout << "myFile record length " << ibuff[0]
      //	  << std::endl;
//--- read nPartEv
      pSw = 1 * bytePerInt;
      iret = read( fh, ibuff, pSw );
// ----------------------------------
      nPartEv = ibuff[0];
      //std::cout << "nPartEv " << ibuff[0] << std::endl;

      return true;

  }


//==========================================================================
  int CsRCEventParticles::openMyFile( char *fileName ) {
//------------------------------------------------------


//--- Paolo  -  November 2001
//    rev.      March 2003


      CsRCExeKeys *key = CsRCExeKeys::Instance();

      int bytePerInt = sizeof(int);
      static int fh;
      int iret;
      int pSw;

      int ibuff[10];

      fh = open( fileName, O_RDONLY, 0 );
//--------------------------------------
      if( fh < 0 ) {
	cout << endl;
        cout << " RICHONE, CsRCEventParticles::setMyParticle() : file ";
	cout << fileName << " not found or error in opening ";
	cout << "( " << fh << " )" << endl;
	return fh;
      }

      cout << endl;
      cout << " RICHONE, CsRCEventParticles::setMyParticle() : file ";
      cout << fileName  << " opened" << endl;

//--- added   090211
      CsOpt* opt = CsOpt::Instance();
      bool boo = false;
      string outFile = "   ";
      boo = opt->CsOpt::getOpt( "histograms", "home", outFile );
      if( boo ) {
	cout << endl;
	cout << " RICHONE, CsRCEventParticles::setMyParticle() : file ";
	cout << outFile  << " in use as output file" << endl;
      }

//--- skip myFile header
      ///pSw = 2 * bytePerInt;
      ///iret = read( fh, ibuff, pSw );
      int nSkip = 3;
      if( key->readGf5() ) nSkip = 5;
      pSw = nSkip * bytePerInt;
      iret = read( fh, ibuff, pSw );
//---------------------------------
      int fhl = ibuff[0];
      int fhd = ibuff[1];
      if( key->readGf5() ) fhd = ibuff[2];
      int ftr = ibuff[2];
      if( key->readGf5() ) ftr = ibuff[3];
      std::cout << " RICHONE, CsRCEventParticles::openMyFile() : ";
      std::cout << "File header length " << fhl << ", file header "
		<< fhd << ", file header trailer " << ftr << std::endl;
      if( iret == -1 ) {
        cout << " RICHONE, CsRCEventParticles::openMyFile() : ";
	cout << "error ( " << iret << " ) reading file "<< fileName << endl;
	cout << "File skipped" << endl;
	return  iret;
        //exit(1);
      }

      return fh;
  }


//==========================================================================
  void CsRCEventParticles::setRunIndex( char *fileName ) {
//--------------------------------------------------------


//- Paolo  -  November 2002

//- intended for myFile only
//  ------------------------

    CsOpt* opt = CsOpt::Instance();
    CsRCRecConst *cons = CsRCRecConst::Instance();

    static vector<double> vRunNo;
    static vector<double> vRunIndex;
    static vector<double> vRunIndexVS;
    static bool twoind = false;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;

      vRunNo.clear();
      vRunIndex.clear();
      vRunIndexVS.clear();
      vector<float> vflo;
      while( opt->CsOpt::getOptRec( "RICHONE", "runRefrIndex", vflo ) ) {
        vRunNo.push_back( vflo[0] );
        vRunIndex.push_back( vflo[1] );
	if( vflo.size() == 3 ) twoind = true;
	vRunIndexVS.push_back( vflo[2] );
      }
      //std::cout << setprecision(6);
      //for( size_t kRun=0; kRun<vRunNo.size(); kRun++ )
      //  std::cout << vRunNo[kRun] << "  " << vRunIndex[kRun] << "  "
      //     << vRunIndexVS[kRun] << std::endl;
      //std::cout << setprecision(2);
    }

    myFileRun_ = 0;
    if( vRunNo.size() == 0 ) return;

    string sFileName;
    sFileName.assign( fileName );
    //cout << sFileName.size() << endl;
    for( size_t kRun=0; kRun<vRunNo.size(); kRun++ ) {
      stringstream runNumber;
      //runNumber << vRunNo[kRun];
      runNumber << vRunNo[kRun];
      size_t pos = sFileName.find( runNumber.str(), 0 );
      //cout << sFileName << "  ";
      //cout << runNumber.str() << "  " << pos << endl;
      if( pos > sFileName.size() ) continue;
      //cout << runNumber.str() << "  xx  " << pos << endl;
      cons->setCFRefInd( vRunIndex[kRun] );
      if( twoind ) cons->setCFRefIndVS( vRunIndexVS[kRun] );
      myFileRun_ = int( vRunNo[kRun] );

      //int iindex = int( (cons->CFRefInd()-1.)*1000000. );
      double index = cons->CFRefInd();
      cout << " RICHONE, CsRCEventParticles::setRunIndex (myFile) : run "
      //     << vRunNo[kRun];
           << int( vRunNo[kRun] );
      CsRICH1Detector* rich = CsGeom::Instance()->getRich1Detector();
      cout << setprecision( 6 );
      if( rich->getCalibFlag() ) {
	//cout << ", Refractive Index from calibrations : 1.00"
	//     << iindex << endl;
	cout << ", Refractive Index from calibrations : " << index << endl;
        cout << " --------------------------------------------------------"
	     << "----------------------------------------------------"
             << endl;
      } else {
	//cout << ", Refractive Index set to 1.00" << iindex << endl;
	cout << ", Refractive Index set to " << index << endl;
        cout << " --------------------------------------------------------"
	     << "---------------------------------------" << endl;
      }
      cout << setprecision( 2 );
      break;
    }

  }
