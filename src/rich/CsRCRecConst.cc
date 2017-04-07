/*!
   \file    CsRCRecConst.cc
   \brief   Constants for reconstruction for CsRichOne.
   \author  Paolo Schiavon
   \version 0.02,  rev. 20/6/00
   \date    $Date: 2010/03/10 15:10:28 $
*/


#include <iostream>
#include <ostream>
#include <cstdio>
#include <cstdlib>

#include <cmath>

#include "CsRCRecConst.h"

#include "CsErrLog.h"
#include "CsOpt.h"
#include "CsGeom.h"
#include "CsEvent.h"
#include "CsRICH1Detector.h"

#include "CsRCDetectors.h"
#include "CsRCExeKeys.h"

  using namespace std;



  CsRCRecConst* CsRCRecConst::instance_ = 0;

//===========================================================================
  CsRCRecConst* CsRCRecConst::Instance() {
//----------------------------------------
    if( instance_ == 0 ) instance_ = new CsRCRecConst();
    return instance_; 
  }

//===========================================================================
  void CsRCRecConst::print() {
//----------------------------
    cout << setprecision( 3 );
    cout << "   " << endl;
    cout << "CsRCRecConst::print() : Reconstruction Constants" << endl;
    cout << "------------------------------------------------" << endl;
    cout << "TwoPi = " << TwoPI_ << ",  rad>deg = " << RadDeg_ << endl;
    cout << "hOffHi = " << kOffHi_ << endl;
    cout << "momMinAcc = " << momMinAcc_ << endl;
    cout << "momMaxAcc = " << momMaxAcc_ << endl;
    cout << "outBufferSize = " << outBufferSize_ << endl;
    cout << "momMinProc = " << momMinProc_ << endl;
    cout << "momMaxProc = " << momMaxProc_ << endl;
    cout << "PHmipCut = " << PHmipCut_ << endl;
    cout << "maxCluSize = " << maxCluSize_ << endl;
    cout << "ddCluUse = " << ddCluUsex_ << ", " <<  ddCluUsey_ << endl;
    cout << setprecision( 7 );
    cout << "CFRefInd = " << CFRefInd_ << endl;
    cout << "CFRefIndUV = " << CFRefIndUV_ << " - ";
    cout << "CFRefIndVS = " << CFRefIndVS_ << endl;
    if( MCCFRefIndSet_ ) cout << "MCCFRefInd = " << MCCFRefInd_ << endl;
    if( MCCFRefIndSetVS_ ) cout << "MCCFRefIndVS = " << MCCFRefIndVS_ << endl;
    cout << "qzRefInd = " << qzRefInd_ << endl;
    cout << setprecision( 3 );
    cout << "partPathFr = " << partPathFr_ << endl;
    cout << "ringDefMode = " << ringDefMode_ << endl;
    cout << "peakSrcMode = " << peakSrcMode_ << endl;
    cout << "mcanWind = " << mcanWind_[0] << ", " << mcanWind_[1] << endl;
    cout << "nMore = " << nMore_ << endl;
    cout << "mcaScan = " << mcaScan_ << endl;
    cout << "xl/uScan = " << xlScan_ << ", " << xuScan_ << endl;
    cout << "binScan = " << binScan_ << endl;
    cout << "nPhotMin = " << nPhotMin_ << endl;
    cout << "sigBfCut = " << sigBfCut_ << endl;
    cout << "nSigmaBf = " << nSigmaBf_[0] << ", " << nSigmaBf_[1] << endl;
    cout << "nBinPhi = " << nBinPhi_ << endl;
    cout << "phiBin = " << phiBin_ << endl;
    cout << "binThr = " << binThr_ << endl;
    cout << "weight = " << weight_ << endl;
    cout << "sig3SCut = " << sig3SCut_ << endl;
    cout << "theAliMir = " << theAliMir_ << endl;
    cout << "rrAliMir = " << rrAliMir_ << endl;
    cout << "nPhoAliMir = " << nPhoAliMir_ << endl;
    cout << "nSigAliMir = " << nSigAliMir_ << endl;
    cout << "massPartv = ";
    for( int k=1; k<16; k++ ) cout << massPartv_[k] << ", ";
    cout << endl;
    cout << "nPhotMinRing = " << nPhotMinRing_ << endl;
    cout << "sigCut = " << sigCut_ << endl;
    cout << "theMaxRgQ = " << theMaxRgQ_ << endl;
    cout << "nPhotAtSat = ";
    for( int k=1; k<16; k++ ) cout << nPhotAtSat_[k] << ", ";
    cout << endl;
    cout << "nZero = " << nZero_ << endl;
    cout << "nZeroUV = " << nZeroUV_ << endl;
    cout << "nZeroVS = " << nZeroVS_ << endl;
    cout << "likeType = " << likeType_ << endl;
    cout << "thetaType = " << thetaType_ << endl;
    cout << "backgrType = " << backgrType_ << endl;
    string action = "NO";
    if( likeRatio_ ) action = "YES";
    cout << "likeRatio = " << action << endl;
    action = "NO";
    if( likeLog_ ) action = "YES";
    cout << "likeLog = " << action << endl;
    cout << "likeDefVa = " << likeDefVa_ << endl;
    cout << "dispMode = " << dispMode_ << "  "
	 << "nEveDSkip = " << nEveDSkip_ << "  "
	 << "nEveDisplay = " << nEveDisplay_ << endl;

    cout << setprecision( 1 );
    cout << "cathode pad PH cath thresholds = ";
    for( int k=0; k<16; k++ ) cout << "(" << k << ")" << cathodeThre_[k]
				   << ",  ";
    cout << endl;
    cout << "cathode pad PH pad thresholds = ";
    for( int k=0; k<16; k++ ) cout << "(" << k << ")" << cathPadThre_[k]
				   << ",  ";
    cout << endl;
    cout << setprecision( 3 );
    cout << "cathode N0 (bkgmap) = ";
    for( int k=0; k<16; k++ ) cout << "(" << k << ")" << cathN0_[k]
				   << ",  ";
    cout << endl;
    cout << endl;
    //    cout << "CsRCRecConst pt = " << this << "  " << massPartv_ << endl;
    //    cout << endl;
    cout << setprecision( 2 );
  }

//===========================================================================
  CsRCRecConst::~CsRCRecConst() { }
//-------------------------------


//===========================================================================
  CsRCRecConst::CsRCRecConst() {
//------------------------------


//--- June  2000, rev. 23/8/00


//--- numerical constants :
//    ---------------------
      //      TwoPI_ = 2. * M_PI;
      //      RadDeg_ = 180. / M_PI;
      TwoPI_ = 2. * 3.1415926;
      RadDeg_ = 180. / 3.1415926;


      CsOpt* opt = CsOpt::Instance();
      bool boo = false;
      int kin = 0;
      float flo = 0.;
      string sflo = "   ";
      vector<int> vin;
      vector<float> vflo;
      string YES = "YES";

      bool lprint = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "PrintConsts", kin );
      if( boo ) { if( kin > 0 ) { lprint = true; } }
      printConsts_ = lprint;

//--- RICH output buffer (to CsTrack) size :
//    fixed in coral/src/track/CsTrack.h as
//    #define CSTRACK_RICHDATASIZE (size)
//    --------------------------------------
      outBufferSize_ = CSTRACK_RICHDATASIZE;
//--- Force LOCAL VALUE
      //outBufferSize_ = 21;
      //outBufferSize_ = 30;
//@@---------------------
      std::cout << std::endl;
      std::cout << "RICHONE, CsRCRecConst::CsRCRecConst() :"
		<< " Using a RICH output buffer of " << outBufferSize_
		<< " words" << std::endl;

//--- monitoring HISTOGRAMS number offset :
//    -------------------------------------
      kOffHi_ = 2000;
      boo = opt->CsOpt::getOpt( "RICHONE", "HistoOffset", kin );
      if( boo ) { kOffHi_ = kin; }

//--- In CsRCEventPads::getEventPads :
//    --------------------------------
      bool UsePadThresh = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "UsePadThresh", sflo );
      if( boo ) { if ( sflo == YES ) UsePadThresh = true; }

      for( int k=0; k<16; k++ ) cathPadThre_[k] = 2.5;
      boo = opt->CsOpt::getOpt( "RICHONE", "cathPadThresh", vflo );
      if( boo ) {
	for( int k=0; k<16; k++ ) cathPadThre_[k] = vflo[k];
      }

      //if( UsePadThresh ) {
      //  for( int k=0; k<16; k++ ) cathodeThre_[k] = 2.5;
      //} else {
      //  for( int k=0; k<16; k++ ) cathodeThre_[k] = 0.0;
      //}

      for( int k=0; k<16; k++ ) cathodeThre_[k] = 0.0;
      boo = opt->CsOpt::getOpt( "RICHONE", "cathodeThresh", vflo );
      if( boo ) {
	for( int k=0; k<16; k++ ) cathodeThre_[k] = vflo[k];
      }

//--- map 050204 (Stefano) OBSOLETE
//---------------------------------
//--- In CsRCLikeAllMap::normSignal :
//    -------------------------------
      for( int k=0; k<16; k++ ) cathN0_[k] = 0.;             // default!
      boo = opt->CsOpt::getOpt( "RICHONE", "cathN0", vflo );
      if( boo ) {
	for( int k=0; k<16; k++ ) cathN0_[k] = vflo[k];
      }

//--- In CsRCEventParticles::getEventParticles :
//    ------------------------------------------
      momMinAcc_ = 0.;
      boo = opt->CsOpt::getOpt( "RICHONE", "momMinAcc", flo );
      if( boo ) { momMinAcc_ = flo; }
//--- In CsRCEventParticles::getEventParticles :
//    ------------------------------------------
      momMaxAcc_ = 180.;
      boo = opt->CsOpt::getOpt( "RICHONE", "momMaxAcc", flo );
      if( boo ) { momMaxAcc_ = flo; }

//--- In CsRCPartPhotons::getPartPhotons :
//    ------------------------------------
      PHmipCut_ = 2000000.;
      boo = opt->CsOpt::getOpt( "RICHONE", "PHmipCut", flo );
      if( boo ) {
	PHmipCut_ = flo;
      }

//--- In CsRCPartPhotons::doSelCathodes :
//    -----------------------------------
      ddCluUsex_ = 200.;
      ddCluUsey_ = 200.;
      boo = opt->CsOpt::getOpt( "RICHONE", "CathodeUseReg", vflo );
      if( boo ) {
	ddCluUsex_ = vflo[0];
	ddCluUsey_ = vflo[1];
      }

//--- In CsRCEventClusters.cc :
//    -------------------------
//--- Maximun cluster size :
//    ----------------------
      maxCluSize_ = 1000.;
      boo = opt->CsOpt::getOpt( "RICHONE", "MaxCluSize", flo );
      if( boo ) {
	maxCluSize_ = flo;
      }

//--- In CsRCMirrors.cc, CsRCPartPhotons.cc, CsRCEventRings.cc,
//    CsRCRing.cc, CsRCEventAnalysis.cc, CsRCEventDisplay.cc :
//    --------------------------------------------------------
//--- C4F10 nominal refraction index :
//    --------------------------------
      CFRefInd_ = 1.00153;
      boo = opt->CsOpt::getOpt( "RICHONE", "C4F10RefrIndex", flo );
      if( boo ) { CFRefInd_ = flo; }
      CFRefIndUV_ = CFRefInd_;
      boo = opt->CsOpt::getOpt( "RICHONE", "C4F10RefrIndexUV", flo );
      if( boo ) { CFRefIndUV_ = flo; }
      CFRefIndVS_ = CFRefInd_ - 0.0002;
      boo = opt->CsOpt::getOpt( "RICHONE", "C4F10RefrIndexVS", flo );
      if( boo ) { CFRefIndVS_ = flo; }

      MCCFRefInd_ = CFRefInd_;
      //MCCFRefIndUV_ = CFRefInd_;           //   090213
      //MCCFRefIndVS_ = CFRefInd_;           //   090213
      MCCFRefIndUV_ = CFRefIndUV_;
      MCCFRefIndVS_ = CFRefIndVS_;
      MCCFRefIndSet_ = false;
      MCCFRefIndSetUV_ = false;
      MCCFRefIndSetVS_ = false;

//--- In CsRCPartPhotons.cc :
//    -----------------------
//--- quartz refraction index :
//    -------------------------
      qzRefInd_ = 1.6;
      boo = opt->CsOpt::getOpt( "RICHONE", "QuartzRefrIndex", flo );
      if( boo ) { qzRefInd_ = flo; }

//--- fraction of particle path for PHOTON 'EMISSION' point
//    (assume particle path middle point) :
//    -----------------------------------------------------
      partPathFr_ = 0.5;
      boo = opt->CsOpt::getOpt( "RICHONE", "ParticlePathFrac", flo );
      if( boo ) { partPathFr_ = flo; }

//--- In CsRCEventRings.cc :
//    ----------------------
      ringDefMode_ = "PEAK";
      boo = opt->CsOpt::getOpt( "RICHONE", "RingDefMode", sflo );
      if( boo ) {
	if( sflo == "MAXLIKE" ) ringDefMode_ = "MAXLIKE";
      }

//--- In CsRCRing.cc :
//    ----------------
//--- peak search mode :
//    ------------------
      peakSrcMode_ = "COUNT";
      boo = opt->CsOpt::getOpt( "RICHONE", "PeakSearchMode", sflo );
      if( boo ) {
	if( sflo == "COUNT" ) peakSrcMode_ = "COUNT";
	else if( sflo == "MASS" ) peakSrcMode_ = "MASS";
	else if( sflo == "LIKE" ) peakSrcMode_ = "LIKE";
	else if( sflo == "CHI" ) peakSrcMode_ = "CHI";
	else {
          string str = 
	    "RICHONE, CsRCRecConst(): UNAVAILABLE Peak Search Mode";
          CsErrLog::Instance()->mes( elError, str );
	}
      }
//--- peak scan window :
//    ------------------
      mcanWind_[0] = 2;             // for kDoClu = 0 and kDoClu = 1
      mcanWind_[1] = 5;
      //mcanWind_[1] = 7;           // for kDoClu = 0 and kDoClu = 1 ???
      boo = opt->CsOpt::getOpt( "RICHONE", "PeakScanWindow", vin );
      if( boo ) {
	mcanWind_[0] = vin[0];
	mcanWind_[1] = vin[1];
      }
//--- enlarged peak window :
//    ----------------------
      nMore_ = 3;
      boo = opt->CsOpt::getOpt( "RICHONE", "nMore", kin );
      if( boo ) {
	nMore_ = kin;
      }

//--- binning of theta peak :
//    -----------------------
      mcaScan_ = 60;
      xlScan_ =   0.;
      xuScan_ = 60.;
      boo = opt->CsOpt::getOpt( "RICHONE", "ThetaPeakBins", kin );
      if( boo ) { mcaScan_ = kin; }
      boo = opt->CsOpt::getOpt( "RICHONE", "ThetaPeakRange", vflo );
      if( boo ) {
        xlScan_ = vflo[0];
        xuScan_ = vflo[1];
      }
      if( mcaScan_ > 200 ) mcaScan_ = 200;
//@@-------------------------------------
      binScan_ = (xuScan_ - xlScan_) / mcaScan_;


//--- In CsRCEventRings.cc ( background FILTER ) :
//    --------------------------------------------
      nPhotMin_ = 0;         //  minimum number of photons per ring
      boo = opt->CsOpt::getOpt( "RICHONE", "nMinPhotons", kin );
      if( boo ) { nPhotMin_ = kin; }

      sigBfCut_ = 1.0;       //  average sigma single photon for Bf cut, mrad
      boo = opt->CsOpt::getOpt( "RICHONE", "SigmaSingPhot", flo );
      if( boo ) { sigBfCut_ = flo; }

      nSigmaBf_[0] = 1.;     //  cut at nSigmaBf sigma's
      nSigmaBf_[1] = 3.;     //  for kDoClu = 0 and kDoClu = 1
      boo = opt->CsOpt::getOpt( "RICHONE", "nSigmaCut", vflo );
      if( boo ) {
        nSigmaBf_[0] = vflo[0];
        nSigmaBf_[1] = vflo[1];
      }

//---                      ( PHI Peak FILTER - Weighted average ) :
//    -------------------------------------------------------------
      nBinPhi_ = 120.;       //  number of bins in phi
      boo = opt->CsOpt::getOpt( "RICHONE", "nBinsInPhi", flo );
      if( boo ) { nBinPhi_ = flo; }
      phiBin_ = 360. / nBinPhi_;     //  bin size in phi
//@@---------------------------

      binThr_ = 1;           //  threshold bin cont.
      boo = opt->CsOpt::getOpt( "RICHONE", "nBinThreshold", kin );
      if( boo ) { binThr_ = kin; }

      weight_ = 0.5;         //  weight of the bin
      boo = opt->CsOpt::getOpt( "RICHONE", "BinWeight", flo );
      if( boo ) { weight_ = flo; }

//---                      ( THREE sigma (single photon) CUT ) :
//    ----------------------------------------------------------
      sig3SCut_ = 1.000;      //  average sigma single photon, mrad, redefined
      boo = opt->CsOpt::getOpt( "RICHONE", "Sigma3SCut", flo );
      if( boo ) { sig3SCut_ = flo; }

//--- mirror alignment :
//    ------------------
      theAliMir_ = 53.;
      boo = opt->CsOpt::getOpt( "RICHONE", "MinThetaAlign", flo );
      if( boo ) { theAliMir_ = flo; }

      rrAliMir_ = 70.;
      boo = opt->CsOpt::getOpt( "RICHONE", "MaxInternalRad", flo );
      if( boo ) { rrAliMir_ = flo; }

      nPhoAliMir_ = 5.;
      boo = opt->CsOpt::getOpt( "RICHONE", "nMinPhotons", flo );
      if( boo ) { nPhoAliMir_ = flo; }

      nSigAliMir_ = 3.;
      boo = opt->CsOpt::getOpt( "RICHONE", "nSigmaPhoCut", flo );
      if( boo ) { nSigAliMir_ = flo; }


//--- In CsRCEventAnalysis.cc :
//    -------------------------
//--- particle masses (GEANT coding) :
//    --------------------------------
      for ( int km=0; km < 31; km++ ) { massPartv_[km] = -1.; }

      massPartv_[ 1] =  0.;           //  1   gamma 
      massPartv_[ 2] =  0.00051;      //  2   e+ mass, GeV/c**2
      massPartv_[ 3] =  0.00051;      //  3   e- mass, GeV/c**2

      massPartv_[ 5] =  0.1057;       //  5   mu+ mass, GeV/c**2
      massPartv_[ 6] =  0.1057;       //  6   mu- mass, GeV/c**2

      massPartv_[ 8] =  0.13957;      //  8   pi+ mass, GeV/c**2
      massPartv_[ 9] =  0.13957;      //  9   pi- mass, GeV/c**2

      massPartv_[11] =  0.49368;      // 11   k+ mass, GeV/c**2
      massPartv_[12] =  0.49368;      // 12   k- mass, GeV/c**2

      massPartv_[14] =  0.93827;      // 14   p mass, GeV/c**2
      massPartv_[15] =  0.93827;      // 15   antip mass, GeV/c**2

      nPhotMinRing_ = 0;
      boo = opt->CsOpt::getOpt( "RICHONE", "nPhotMinRing", kin );
      if( boo ) { nPhotMinRing_ = kin; }

      momMinProc_ = 0.;
      boo = opt->CsOpt::getOpt( "RICHONE", "momMinProc", flo );
      if( boo ) { momMinProc_ = flo; }

      momMaxProc_ = 180.;
      boo = opt->CsOpt::getOpt( "RICHONE", "momMaxProc", flo );
      if( boo ) { momMaxProc_ = flo; }

      sigCut_ = 3.;
      boo = opt->CsOpt::getOpt( "RICHONE", "nSigmaCut", flo );
      if( boo ) { sigCut_ = flo; }

      float thetaMaxLk = 70.;
      boo = opt->CsOpt::getOpt( "RICHONE", "thetaMaxLk", flo );
      if( boo ) { thetaMaxLk = flo; }
      theMaxRgQ_ = thetaMaxLk*thetaMaxLk;

//--- n-Photon-at-satutation :
//    ------------------------
      for( int k=0; k<16; k++ ) nPhotAtSat_[k] = 30.;
      bPhotAtSat_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "nPhotAtSat", vflo );
      if( boo ) {
	for( int k=0; k<16; k++ ) nPhotAtSat_[k] = vflo[k];
	bPhotAtSat_ = true;
      }

      nZero_ = 30.;
      boo = opt->CsOpt::getOpt( "RICHONE", "nZero", flo );
      if( boo ) { nZero_ = flo; }
      nZeroUV_ = 30.;
      boo = opt->CsOpt::getOpt( "RICHONE", "nZeroUV", flo );
      if( boo ) { nZeroUV_ = flo; }
      nZeroVS_ = 30.;
      boo = opt->CsOpt::getOpt( "RICHONE", "nZeroVS", flo );
      if( boo ) { nZeroVS_ = flo; }

      thetaType_ = "REC";
      boo = opt->CsOpt::getOpt( "RICHONE", "ThetaType", sflo );
      if( boo ) {
	if( sflo == "REC" ) thetaType_ = "REC";
	else if( sflo == "LIKE" ) thetaType_ = "LIKE";
	else if( sflo == "FIT" ) thetaType_ = "FIT";
	else if( sflo == "WAVE" ) thetaType_ = "WAVE";
	else {
          string str = 
	    "RICHONE, CsRCRecConst(): UNAVAILABLE Cherenkov angle calc. type";
          CsErrLog::Instance()->mes( elError, str );
	}
      }
      likeType_ = "ALL";
      boo = opt->CsOpt::getOpt( "RICHONE", "LikeType", sflo );
      if( boo ) { 
	if( sflo == "RING" ) likeType_ = "RING";
	else if( sflo == "ALL" ) likeType_ = "ALL";
	else if( sflo == "MIX" ) likeType_ = "MIX";
	else if( sflo == "CHISQ" ) likeType_ = "CHISQ";
	else {
          string str = 
	    "RICHONE, CsRCRecConst(): UNAVAILABLE Likelihood type";
          CsErrLog::Instance()->mes( elError, str );
	}
      }
      backgrType_ = "05";
      boo = opt->CsOpt::getOpt( "RICHONE", "BackgrType", sflo );
      if( boo ) { 
	if( sflo == "01" ) backgrType_ = "01";
	else if( sflo == "02" ) backgrType_ = "02";
	else if( sflo == "03" ) backgrType_ = "03";
	else if( sflo == "04" ) backgrType_ = "04";
	else if( sflo == "05" ) backgrType_ = "05";
//map 050204 (Stefano) OBSOLETE!
	else if( sflo == "BKGMAP" ) {
	  backgrType_ = "BKGMAP";
	  boo = opt->CsOpt::getOpt( "RICHONE", "BackgrFile", sflo );
	  FILE *file;
	  const char *filename = sflo.c_str();
	  file = fopen(filename, "r");
	  if( !file ) {
	    std::string str =
	      " RICHONE, CsRCRecConst(): UNAVAILABLE Background Map";
            CsErrLog::Instance()->mes( elFatal, str );
	  }
	  else {
	    for(int j=0; j<400; j++) {
	      for(int i=0; i<400; i++) {
		fscanf (file,"%lf;", &backgrmap_[i][j]);
	      }
	    }
	    fclose( file );
	    std::cout << std::endl;
	    std::cout << "RICHONE, CsRCRecConst() : ";
	    std::cout << "Background map for BKGMAP likelihood in use : ";
	    std::cout << filename << std::endl;
	    std::cout << "--------------------------";
	    std::cout << "---------------------------------------------";
	    std::cout << std::endl;
	  }
	}
//map
	else {
          string str = 
	    "RICHONE, CsRCRecConst(): UNAVAILABLE Background parametrization";
          CsErrLog::Instance()->mes( elError, str );
	}
      }

      likeRatio_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "LikeRatio", sflo );
      if( boo ) { if ( sflo == YES ) likeRatio_ = true; }

      likeLog_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "LikeLog", sflo );
      if( boo ) { if ( sflo == YES ) likeLog_ = true; }

//--- Liklelihood Default Value :
//    ---------------------------
      likeDefVa_ = 0.;
      if( likeLog_ ) likeDefVa_ = -9999.;
//    ----------------------------------

//--- In CsRCEventDisplay.cc :
//    ------------------------
      nEveDSkip_ = 0;
      boo = opt->CsOpt::getOpt( "RICHONE", "nEveDSkip", kin );
      if( boo ) { nEveDSkip_ = kin; }
      nEveDisplay_ = 0;
      boo = opt->CsOpt::getOpt( "RICHONE", "nEveDisplay", kin );
      if( boo ) { nEveDisplay_ = kin; }
      string strd;
      dispMode_ = "ALL";
      boo = opt->CsOpt::getOpt( "RICHONE", "DisplayMode", strd );
      if( boo ) {
	if( strd == "ALL" ) dispMode_ = "ALL";
	else if( strd == "REC" ) dispMode_ = "REC";
	else if( strd == "NOREC" ) dispMode_ = "NOR";
	else {
          string str = 
	    "RICHONE, CsRCRecConst(): UNAVAILABLE Display mode request";
          CsErrLog::Instance()->mes( elError, str );
	}
      }

      //if( lprint ) print();          //   moved to CsRichOne::doRichOne()
                                       //   090213
  }


//===========================================================================
  double CsRCRecConst::CFRefInd() {
//---------------------------------


//- Paolo
//- October  2002
//- rev. May 2006


    CsRCExeKeys *key = CsRCExeKeys::Instance();

    static bool firstCall = true;

    double index = CFRefInd_;

    if( key->MCarloEvent() ) {
      if( MCCFRefIndSet_ ) {
	index = MCCFRefInd_;
	if( firstCall ) {
	  firstCall = false;
	  std::cout << setprecision( 6 ) << endl;
	  if( CsRichOne::Instance()->UpRICHJob() ) {
	    std::cout << " RICHONE, CsRCRecConst::CFRefIndUV(), ";
	    std::cout << " Refractive IndexUV from MC detectors.dat : "
		      << index;
	  } else {
  	    std::cout << " RICHONE, CsRCRecConst::CFRefInd(), ";
	    std::cout << " Refractive Index from MC detectors.dat : "
		      << index;
	  }
	  std::cout << setprecision( 2 ) << std::endl;
	}
      } else {
	if( firstCall ) {
	  firstCall = false;
	  std::cout << setprecision( 6 ) << endl;
	  if( CsRichOne::Instance()->UpRICHJob() ) {
	    std::cout << " RICHONE, CsRCRecConst::CFRefIndUV(), ";
	    std::cout << " Refractive IndexUV from MC rich1.opt : " << index;
	  } else {
	    std::cout << " RICHONE, CsRCRecConst::CFRefInd(), ";
	    std::cout << " Refractive Index from MC rich1.opt : " << index;
	  }
	  std::cout << setprecision( 2 ) << std::endl;
	}
      }
      return  index;
    }

    if( key->readMyFile() ) {
      if( firstCall ) {
	firstCall = false;
	std::cout << setprecision( 6 ) << endl;
        std::cout << " RICHONE, CsRCRecConst::CFRefInd(), ";
        std::cout << " Refractive Index from MyFile rich1.opt : " << index;
        std::cout << setprecision( 2 ) << std::endl;
      }
      return  index;
    }

    CsRICH1Detector* rich = CsGeom::Instance()->getRich1Detector();
    if( rich->getCalibFlag() ) index = rich->getIndex();

    static int oldRun = 0;
    int newRun = CsEvent::Instance()->getRunNumber();
    if( newRun != oldRun ) {
      oldRun = newRun;
      std::cout << std::endl;
      std::cout << " RICHONE, CsRCRecConst::CFRefInd(), ";
      std::cout << " Refractive Index";
//--- from Andrea 050615
      //<< std::setiosflags( std::ios::fixed ) << std::setfill(' ')
      //<< std::setw(8) << std::setprecision(6) << index;
      if( rich->getCalibFlag() ) {
	std::cout << " for Run " << newRun;
	std::cout << " from DB : ";
      } else {
	std::cout << "  NOT from DB : ";
      }
      std::cout << setprecision( 6 ) << index;
      std::cout << setprecision( 2 ) << std::endl;
    }

    return  index;

  }


//===========================================================================
  double CsRCRecConst::CFRefIndVS() {
//----------------------------------


//- Paolo
//- August 2005
//- rev. May 2006

    CsRCExeKeys *key = CsRCExeKeys::Instance();

    static bool firstCall = true;

    double index = CFRefIndVS_;

    if( key->MCarloEvent() ) {
      if( MCCFRefIndSetVS_ ) {
	index = MCCFRefIndVS_;
	if( firstCall ) {
	  firstCall = false;
	  std::cout << setprecision( 6 ) << endl;
  	  std::cout << " RICHONE, CsRCRecConst::CFRefIndVS(), ";
	  std::cout << " Refractive IndexVS from MC detectors.dat : "<< index;
	  std::cout << setprecision( 2 ) << std::endl;
	}
      } else {
	if( firstCall ) {
	  firstCall = false;
	  std::cout << setprecision( 6 ) << endl;
	  std::cout << " RICHONE, CsRCRecConst::CFRefIndVS(), ";
	  std::cout << " Refractive IndexVS from MC rich1.opt : " << index;
	  std::cout << setprecision( 2 ) << std::endl;
	}
      }
      return  index;
    }

    if( key->readMyFile() ) {
      if( firstCall ) {
	firstCall = false;
	std::cout << setprecision( 6 ) << endl;
        std::cout << " RICHONE, CsRCRecConst::CFRefIndVS(), ";
        std::cout << " Refractive IndexVS from MyFile rich1.opt : " << index;
        std::cout << setprecision( 2 ) << std::endl;
      }
      return  index;
    }

    CsRICH1Detector* rich = CsGeom::Instance()->getRich1Detector();
    if( rich->getCalibFlag() ) index = rich->getIndexVS();

    static int oldRun = 0;
    int newRun = CsEvent::Instance()->getRunNumber();
    if( newRun != oldRun ) {
      oldRun = newRun;
      std::cout << endl;
      std::cout << " RICHONE, CsRCRecConst::CFRefIndVS(), ";
      std::cout << " Refractive IndexVS";
//--- from Andrea 050615
      //<< std::setiosflags( std::ios::fixed ) << std::setfill(' ')
      //<< std::setw(8) << std::setprecision(6) << index;
      if( rich->getCalibFlag() ) {
	std::cout << "  for Run  " << newRun;
	std::cout << "  from DB : ";
      } else {
	std::cout << "  NOT from DB : ";
      }
      std::cout << setprecision( 6 ) << index;
      std::cout << setprecision( 2 ) << std::endl;
    }

    return  index;
  }


//===========================================================================
  float* CsRCRecConst::cathodeThre() {
//------------------------------------


//- October  2002


    int nCathode = CsRCDetectors::Instance()->nCathode();
    float* thresh = cathodeThre_;

    CsRICH1Detector* rich = CsGeom::Instance()->getRich1Detector();
    if( rich->getCalibFlag() ) {
      for( int k=0; k<nCathode; k++ ) thresh[k] = rich->getThresh();
    }

    return  thresh;

  }
