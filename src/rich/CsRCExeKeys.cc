/*!
   \file    CsRCExeKeys.cc
   \brief   Execution Keys for CsRichOne.
   \author  Paolo Schiavon
   \version 0.02,  rev. 20/6/00
   \date    June  2000
*/


  #include <iostream>
  #include <ostream>
  #include <cstdio>
  #include <cstdlib>

  #include <cmath>

  #include "CsRCExeKeys.h"

  #include "Coral.h"
  #include "CsOpt.h"

  using namespace std;

  CsRCExeKeys* CsRCExeKeys::instance_ = 0;

//===========================================================================
  CsRCExeKeys* CsRCExeKeys::Instance() {
//--------------------------------------
    if( instance_ == 0 ) instance_ = new CsRCExeKeys();
    return instance_; 
  }

//===========================================================================
  CsRCExeKeys::~CsRCExeKeys() { }
//-------------------------------


//===========================================================================
  CsRCExeKeys::CsRCExeKeys() {
//----------------------------


//--- June  2000,  rev. July  2000


      CsOpt* opt = CsOpt::Instance();

      std::cout << setiosflags(ios::boolalpha);
//@@  ----------------------------------------

      bool boo;
      int bOpw;
      string sOpw;
      string YES = "YES";
      string NO = "NO";
      int kOpw;

      bool lprint = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "PrintKeys", kOpw );
      if( boo ) { if( kOpw > 0 ) { lprint = true; } }


//--- use MONTECARLO or DATA :
//    ------------------------
      Coral* coral = Coral::Instance();
      MCarloEvent_ = false;
      if( coral->isAMonteCarloEvent() ) MCarloEvent_ = true;
      DataEvent_ = false;
      if( coral->isADataEvent() ) DataEvent_ = true;
      //      DataEvent_ = true;      // !!!

//--- myFile :
//    --------
      readMyFile_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "readMyFile", sOpw );
      if( boo ) { if( sOpw == YES ) readMyFile_ = true; }
      RmyFileName_.clear();
      while( opt->CsOpt::getOptRec( "RICHONE", "readMyFileName", sOpw ) ) {
        RmyFileName_.push_back( sOpw );
      }

      int myFileChunks = 1;
      int nRobs = 9;
      boo = opt->CsOpt::getOpt( "RICHONE", "myFileChunks", kOpw );
      if( boo ) myFileChunks = kOpw;
      boo = opt->CsOpt::getOpt( "RICHONE", "myFileRobs", kOpw );
      if( boo ) nRobs = kOpw;
      //cout << myFileChunks << "  " << nRobs << endl;
      if( RmyFileName_.size() == 1  &&  myFileChunks > 1 )
	fillMyFileList( myFileChunks, nRobs );
//      -------------------------------------

      myFileType_ = 2;
      boo = opt->CsOpt::getOpt( "RICHONE", "myFileType", kOpw );
      if( boo ) myFileType_ = kOpw;

      nProMyEvs_ = 0;
      boo = opt->CsOpt::getOpt( "RICHONE", "nProMyEvs", bOpw );
      if( boo ) nProMyEvs_ = bOpw;

      nSkipMyEvs_ = 0;
      boo = opt->CsOpt::getOpt( "RICHONE", "nSkipMyEvs", bOpw );
      if( boo ) nSkipMyEvs_ = bOpw;

      writeMyFile_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "writeMyFile", sOpw );
      if( boo ) { if( sOpw == YES ) writeMyFile_ = true; }
      //      if( readMyFile_ ) writeMyFile_ = false;
//@@----------------------------------------

      WmyFileName_ = " ";
      boo = opt->CsOpt::getOpt( "RICHONE", "writeMyFileName", sOpw );
      if( boo ) WmyFileName_ = sOpw;

      endWriteMyFile_ = false;

      physSelection_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "physSelection", sOpw );
      if( boo ) { if( sOpw == YES ) physSelection_ = true; }

//--- PARTICLE SELECTION in input :
//    -----------------------------

//--- filter-out unwanted particles :
//    -------------------------------
      PartFilt_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "ParticleFilter", sOpw );
      if( boo ) { if ( sOpw == YES ) PartFilt_ = true; }

//--- particle type selection ( MC only ) :
//    -------------------------------------
      PartSele_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "ParticleSelection", sOpw );
      if( boo ) { if ( sOpw == YES ) PartSele_ = true; }

//--- check on exit window selection :
//    --------------------------------
      ExitWind_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "ExitWindowCheck", sOpw );
      if( boo ) { if ( sOpw == YES ) ExitWind_ = true; }

//--- mcs correction to particle direction :
//    --------------------------------------
      PartCorr_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "ParticleCorrection", sOpw );
      if( boo ) { if ( sOpw == YES ) PartCorr_ = true; }

//--- measurement errors on particle ( MC only ) :
//    --------------------------------------------
      PartReso_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "ParticleResolution", sOpw );
      if( boo ) { if ( sOpw == YES ) PartReso_ = true; }

//--- PAD THRESHOLD :
//    ---------------
      UsePadThresh_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "UsePadThresh", sOpw );
      if( boo ) { if ( sOpw == YES ) UsePadThresh_ = true; }

//--- PAD TYPE SELECTION :
//    --------------------
      selPMTonly_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "selectPMTonly", sOpw );
      if( boo ) { if ( sOpw == YES ) selPMTonly_ = true; }
      selAPVonly_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "selectAPVonly", sOpw );
      if( boo ) { if ( sOpw == YES ) selAPVonly_ = true; }
      if( selPMTonly_  &&  selAPVonly_ ) {
	selPMTonly_ = false;
	selAPVonly_ = false;
      }

//--- PAD HALO :
//    ----------
      KillHaloPads_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "killHaloPads", sOpw );
      if( boo ) { if ( sOpw == YES ) KillHaloPads_ = true; }

//--- PROCESS PADS SCRAMBLED :
//    ------------------------
      padScrambled_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "padScrambled", sOpw );
      if( boo ) { if ( sOpw == YES ) padScrambled_ = true; }

//--- CLUSTERING :
//    ------------
      DoClu_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "doClustering", sOpw );
      if( boo ) { if ( sOpw == YES ) DoClu_ = true; }

//--- CLUSTER HALO :
//    --------------
      KillHaloClus_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "killHaloClus", sOpw );
      if( boo ) { if ( sOpw == YES ) KillHaloClus_ = true; }

//--- No Ring Splitting :
//    -------------------
      //NoRingSplitting_ = false;
      //boo = opt->CsOpt::getOpt( "RICHONE", "NoRingSplitting", sOpw );
      //if( boo ) { if ( sOpw == YES ) NoRingSplitting_ = true; }
      UseSplitRings_ = true;
      boo = opt->CsOpt::getOpt( "RICHONE", "UseSplitRings", sOpw );
      if( boo ) { if ( sOpw == NO ) UseSplitRings_ = false; }

//--- CORRECTION for detector QUARTZ WINDOW :
//    ---------------------------------------
      CorrQzW_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "doQuartzWcorr", sOpw );
      if( boo ) { if ( sOpw == YES ) CorrQzW_ = true; }

//--- CORRECTION for PMT Optics :
//    ---------------------------
      CorrPMTOpt_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "doPMTOptCorr", sOpw );
      if( boo ) { if ( sOpw == YES ) CorrPMTOpt_ = true; }

//--- KEYS for LIKELIHOOD TESTS FLOW :
//    -------------------------------
      likeONLY_ = false;
//@@-------------------
      //likeONLY_ = true;      //@@@
//@@--------------------
      doCheckLike_ = true;
      boo = opt->CsOpt::getOpt( "RICHONE", "DoCheckLikelihood", sOpw );
      if( boo ) { if ( sOpw == NO ) doCheckLike_ = false; }

      doThetaLikeMax_ = true;
      boo = opt->CsOpt::getOpt( "RICHONE", "DoThetaLikeMax", sOpw );
      if( boo ) { if ( sOpw == NO ) doThetaLikeMax_ = false; }
      thetaLikeMaxTrue_ = true;
      boo = opt->CsOpt::getOpt( "RICHONE", "ThetaLikeMaxTrue", sOpw );
      if( boo ) { if ( sOpw == NO ) thetaLikeMaxTrue_ = false; }

//--- SEARCH FOR RINGS :
//    ------------------
      SearchForRings_ = true;
      boo = opt->CsOpt::getOpt( "RICHONE", "SearchForRings", sOpw );
      if( boo ) { if ( sOpw == NO ) SearchForRings_ = false; }

//--- use cluster PH in PEAK search :
//    -------------------------------
      UseCluPH_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "useClusterPH", sOpw );
      if( boo ) { if ( sOpw == YES ) UseCluPH_ = true; }

//--- reject PEAK :
//    -------------
      RejWroPeak_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "rejWrongPeak", sOpw );
      if( boo ) { if ( sOpw == YES ) RejWroPeak_ = true; }

//--- background FILTER :
//    -------------------
      DoBkgFil_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "doBkgFilter", sOpw );
      if( boo ) { if ( sOpw == YES ) DoBkgFil_ = true; }

//--- PHI Peak FILTER (Weighted average) :
//    ------------------------------------
      DoWgtAv_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "doPhiPeakFilter", sOpw );
      if( boo ) { if ( sOpw == YES ) DoWgtAv_ = true; }
      //      kDoWgtAv_ = 1;

//--- THREE sigma (single photon) CUT :
//    ---------------------------------
      SigmaCut_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "doThreeSigmaCut", sOpw );
      if( boo ) { if ( sOpw == YES ) SigmaCut_ = true; }


//--- Event DISPLAY :
//    ---------------
      EventDisp_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "EventDisplay", sOpw );
      if( boo ) { if ( sOpw == YES ) EventDisp_ = true; }

//--- Interactiv Event DISPLAY (Andrea) :
//    -----------------------------------
      EventROOTDisp_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "EventROOTDisplay", sOpw );
      if( boo ) { if ( sOpw == YES ) EventROOTDisp_ = true; }

//--- Event ANALYSIS :
//    ----------------
      EventAnasy_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "EventAnalysis", sOpw );
      if( boo ) { if ( sOpw == YES ) EventAnasy_ = true; }

//--- Event selection (re-analysis) :
//    -------------------------------
      PartPhoSelect_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "PartPhoSelection", sOpw );
      if( boo ) { if ( sOpw == YES ) PartPhoSelect_ = true; }

      RingSelect_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "RingSelection", sOpw );
      if( boo ) { if ( sOpw == YES ) RingSelect_ = true; }

//--- MC monitoring :
//    ---------------
      MCMoni_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "MCMonitor", sOpw );
      if( boo ) { if ( sOpw == YES ) MCMoni_ = true; }

//--- data monitoring :
//    -----------------
      DataMoni_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "DataMonitor", sOpw );
      if( boo ) { if ( sOpw == YES ) DataMoni_ = true; }

//--- particle identification :
//    -------------------------
      PartIdent_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "PartIdent", sOpw );
      if( boo ) { if ( sOpw == YES ) PartIdent_ = true; }

//--- mirror data correction procedure :
//    ----------------------------------
      CorrMirror_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "doMirrorCorr", sOpw );
      if( boo ) { if ( sOpw == YES ) CorrMirror_ = true; }
//--- mirror alignment procedure :
//    ----------------------------
      AliMirror_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "doMirrorAlign", sOpw );
      if( boo ) { if ( sOpw == YES ) AliMirror_ = true; }


//--- PRINTING ( debugging ) keys :
//    -----------------------------
//    ( kPrint... = 0 : NO print,  ---  note different meaning.
//      kPrint... = k : print level )

      AcknoMethods_ = 0;
      boo = opt->CsOpt::getOpt( "RICHONE", "AcknoMethods", kOpw );
      if( boo ) AcknoMethods_ = kOpw;

      kPrintRichRec_ = 0;
      boo = opt->CsOpt::getOpt( "RICHONE", "kPrintRichRec", kOpw );
      if( boo ) kPrintRichRec_ = kOpw;

      kPrintEventParts_ = 0;
      boo = opt->CsOpt::getOpt( "RICHONE", "kPrintEventParts", kOpw );
      if( boo ) kPrintEventParts_ = kOpw;

      kPrintEventSimul_ = 0;
      boo = opt->CsOpt::getOpt( "RICHONE", "kPrintEventSimul", kOpw );
      if( boo ) kPrintEventSimul_ = kOpw;

      kPrintEventPads_ = 0;
      boo = opt->CsOpt::getOpt( "RICHONE", "kPrintEventPads", kOpw );
      if( boo ) kPrintEventPads_ = kOpw;

      kPrintEventClus_ = 0;
      boo = opt->CsOpt::getOpt( "RICHONE", "kPrintEventClus", kOpw );
      if( boo ) kPrintEventClus_ = kOpw;

      kPrintPartPhotons_ = 0;
      boo = opt->CsOpt::getOpt( "RICHONE", "kPrintPartPhotons", kOpw );
      if( boo ) kPrintPartPhotons_ = kOpw;

      kPrintEventRings_ = 0;
      boo = opt->CsOpt::getOpt( "RICHONE", "kPrintEventRings", kOpw );
      if( boo ) kPrintEventRings_ = kOpw;

      kPrintEventDisplay_ = 0;
      boo = opt->CsOpt::getOpt( "RICHONE", "kPrintEventDisplay", kOpw );
      if( boo ) kPrintEventDisplay_ = kOpw;

      kPrintEventAnalysis_ = 0;
      boo = opt->CsOpt::getOpt( "RICHONE", "kPrintEventAnalysis", kOpw );
      if( boo ) kPrintEventAnalysis_ = kOpw;

      kPrintPartProbs_ = 0;
      boo = opt->CsOpt::getOpt( "RICHONE", "kPrintPartProbs", kOpw );
      if( boo ) kPrintPartProbs_ = kOpw;

      kPrintRejections_ = 0;
      boo = opt->CsOpt::getOpt( "RICHONE", "kPrintRejections", kOpw );
      if( boo ) kPrintRejections_ = kOpw;

      printFrequency_ = 10000000;
      boo = opt->CsOpt::getOpt( "RICHONE", "printFrequency", kOpw );
      if( boo ) printFrequency_ = kOpw;
//--- avoid zero value
      if( printFrequency_ == 0 ) printFrequency_ = 1;

//--- Event NTUPLE :
//    --------------
      writeNtup_ = false;
      boo = opt->CsOpt::getOpt( "RICHONE", "writeNtuple", sOpw );
      if( boo ) { if ( sOpw == YES ) writeNtup_ = true; }

      readGf5_ = false;
      if( myFileType_ == 4 ) readGf5_ = true;


      if( lprint ) print();

  }


//===========================================================================
  void CsRCExeKeys::acknoMethod( string sMethod ) {
//-------------------------------------------------

//--- July  2000

    if( AcknoMethods_ > 0 ) {
      cout << "--- RICHONE :  " << sMethod << " in use" << endl;
    }

  }


//===========================================================================
  void CsRCExeKeys::print() {
//---------------------------


    cout << "   " << endl;
    cout << "CsRCExeKeys::print() : Execution Keys" << endl;
    cout << "-------------------------------------" << endl;

    if( MCarloEvent_ ) {
      cout << "RICHONE MonteCarlo Job" << endl;
      cout << "----------------------" << endl;
    }
    if( DataEvent_ ) {
      cout << "RICHONE Data Job" << endl;
      cout << "----------------" << endl;
    }

    string action = "NO";
    if( readMyFile_ ) action = "YES";
    cout << "read pre-processed g-files    " << action << endl;
    if( readMyFile_ && myFileType_ == 0 ) {
      cout << "OLD type MC files    " << endl;
    }
    if( readMyFile_ && myFileType_ >= 1 ) {
      cout << "NEW type g-files, type " << myFileType_ << endl;
    }
    action = "NO";
    if( writeMyFile_ ) action = "YES";
    cout << "write pre-processed g-files   " << action << endl;
    cout << endl;
    action = "NO";
    if( physSelection_ ) action = "YES";
    cout << "make phys.selection to g-files   " << action << endl;
    cout << endl;

    action = "NO";
    if( SearchForRings_ ) action = "YES";
    cout << "make SearchForRings            " << action << endl;

    action = "NO";
    if( EventDisp_ ) action = "YES";
    cout << "make EventDisplay              " << action << endl;
    action = "NO";
    if( EventROOTDisp_ ) action = "YES";
    cout << "make EventROOTDisplay          " << action << endl;

    action = "NO";
    if( EventAnasy_ ) action = "YES";
    cout << "make EventAnalysis             " << action << endl;
    action = "NO";
    if( PartPhoSelect_ ) action = "YES";
    cout << "make PartPhotonSelection       " << action << endl;
    action = "NO";
    if( RingSelect_ ) action = "YES";
    cout << "make RingSelection             " << action << endl;

    action = "NO";
    if( MCMoni_ ) action = "YES";
    cout << "make MCMonitor                 " << action << endl;
    action = "NO";
    if( DataMoni_ ) action = "YES";
    cout << "make DataMonitor               " << action << endl;
    action = "NO";
    if( PartIdent_ ) action = "YES";
    cout << "make PartIdent                 " << action << endl;
    cout << "   " << endl;
    action = "NO";
    if( PartFilt_ ) action = "YES";
    cout << "make ParticleFilter            " << action << endl;
    action = "NO";
    if( PartSele_ ) action = "YES";
    cout << "make ParticleSelection         " << action << endl;
    action = "NO";
    if( ExitWind_ ) action = "YES";
    cout << "make ExitWindowCheck           " << action << endl;
    action = "NO";
    if( PartCorr_ ) action = "YES";
    cout << "make ParticleCorrection        " << action << endl;
    action = "NO";
    if( PartReso_ ) action = "YES";
    cout << "make ParticleResolution        " << action << endl;

    cout << "   " << endl;
    action = "NO";
    if( UsePadThresh_ ) action = "YES";
    cout << "set Pad Thresholds             " << action << endl;
    action = "NO";
    if( KillHaloPads_ ) action = "YES";
    cout << "make killHaloPads              " << action << endl;
    cout << "   " << endl;
    action = "NO";
    if( selPMTonly_ ) action = "YES";
    cout << "select PMT Pads only           " << action << endl;
    cout << "   " << endl;
    action = "NO";
    if( selAPVonly_ ) action = "YES";
    cout << "select APV Pads only           " << action << endl;
    cout << "   " << endl;

    cout << "   " << endl;
    action = "NO";
    if( DoClu_ ) action = "YES";
    cout << "make Clustering                " << action << endl;
    action = "NO";
    //if( DoCluCo_ ) action = "YES";
    //cout << "make ClusteringCo              " << action << endl;
    if( UseCluPH_ ) action = "YES";
    cout << "use Cluster PH                 " << action << endl;
    action = "NO";
    if( KillHaloClus_ ) action = "YES";
    cout << "make killHaloClus              " << action << endl;
    cout << "   " << endl;

    action = "NO";
    //if( NoRingSplitting_ ) action = "YES";
    if( UseSplitRings_ ) action = "YES";
    cout << "Use Rings Splitted             " << action << endl;
    cout << "   " << endl;
    action = "NO";
    if( CorrQzW_ ) action = "YES";
    cout << "make QuartzWinCorrection       " << action << endl;
    cout << "   " << endl;
    action = "NO";
    if( CorrPMTOpt_ ) action = "YES";
    cout << "make PMTOpticalCorrection      " << action << endl;
    cout << "   " << endl;
    action = "NO";
    if( doCheckLike_ ) action = "YES";
    cout << "use checkLikelihood            " << action << endl;
    cout << "   " << endl;
    action = "NO";
    if( doThetaLikeMax_ ) action = "YES";
    cout << "use getThetaLikeMax            " << action << endl;
    cout << "   " << endl;
    action = "NO";
    if( thetaLikeMaxTrue_ ) action = "YES";
    cout << "use thetaLikeMaxTrue           " << action << endl;
    cout << "   " << endl;
    action = "NO";
    if( RejWroPeak_ ) action = "YES";
    cout << "reject 'wrong' peak            " << action << endl;
    action = "NO";
    if( DoBkgFil_ ) action = "YES";
    cout << "make BackgroundFilter          " << action << endl;
    action = "NO";
    if( DoWgtAv_ ) action = "YES";
    cout << "make PhiPeakFilter             " << action << endl;
    action = "NO";
    if( SigmaCut_ ) action = "YES";
    cout << "make ThreeSigmaCut             " << action << endl;
    cout << "   " << endl;
    action = "NO";
    if( CorrMirror_ ) action = "YES";
    cout << "make MirrorCorrections         " << action << endl;
    action = "NO";
    if( AliMirror_ ) action = "YES";
    cout << "make MirrorAlignment           " << action << endl;
    cout << "   " << endl;
    cout << "Print RichRec                  " << kPrintRichRec_ << endl;
    cout << "Print EventParts               " << kPrintEventParts_ << endl;
    cout << "Print EventPads                " << kPrintEventPads_ << endl;
    cout << "Print EventClus                " << kPrintEventClus_ << endl;
    cout << "Print EventPartPhotons         " << kPrintPartPhotons_ << endl;
    cout << "Print EventRings               " << kPrintEventRings_ << endl;
    cout << "Print EventDisplay             " << kPrintEventDisplay_ << endl;
    cout << "Print EventAnalysis            " << kPrintEventAnalysis_ << endl;
    cout << "Print PartProbs                " << kPrintPartProbs_ << endl;
    cout << "Print Rejections               " << kPrintRejections_ << endl;
    cout << "   " << endl;
    cout << "Print Frequency              " << printFrequency_ << endl;
    cout << "   " << endl;
    action = "NO";
    if( writeNtup_ ) action = "YES";
    cout << "write ntuple                  " << action << endl;
    cout << "   " << endl;
  }


//===========================================================================
  void CsRCExeKeys::fillMyFileList( int nChunks, int nRobs ) {
//------------------------------------------------------------

//- Paolo - March  2003


    if( RmyFileName_.size() > 1 ) return;

    string sFileName = RmyFileName_.front();
    size_t pos = sFileName.find( ".gfile", 0 );
    pos -= 11;
    string path;
    path.assign( sFileName, 0 , pos );
    string snumb1, snumb2;
    snumb1.assign( sFileName, pos , 2 );
    //pos += 4;
    //snumb2.assign( sFileName, pos , 1 );
    pos += 3;
    snumb2.assign( sFileName, pos , 2 );
    stringstream ssnumb1, ssnumb2;
    ssnumb1 << snumb1;
    ssnumb2 << snumb2;
    int numb1, numb2;
    ssnumb1 >> numb1;
    ssnumb2 >> numb2;
    numb2++;
    //pos += 1;
    pos += 2;
    string trail;
    trail.assign( sFileName, pos , 12 );
    int kChunk = 1;
    int numbt = numb1 + nChunks/nRobs;
    //for( int knu1=numb1; knu1<=20; knu1++ ) {
    for( int knu1=numb1; knu1<=numbt; knu1++ ) {
      if( knu1 > numb1 ) numb2 = 1;
      for( int knu2=numb2; knu2<=nRobs; knu2++ ) {
	stringstream ssFileName;
	ssFileName.clear();
	//if( knu1 == 9 ) {
	//ssFileName << path << "0" << knu1 << "00" << knu2 << trail;
	if( knu1 <= 9 ) {
	  if( knu2 <= 9 ) {
	    ssFileName << path << "0" << knu1 << "00" << knu2 << trail;
	  } else {
	    ssFileName << path << "0" << knu1 << "0" << knu2 << trail;
	  }
	} else {
	  //ssFileName << path << knu1 << "00" << knu2 << trail;
	  if( knu2 <= 9 ) {
	    ssFileName << path << knu1 << "00" << knu2 << trail;
	  } else {
	    ssFileName << path << knu1 << "0" << knu2 << trail;
	  }
	}
        RmyFileName_.push_back( ssFileName.str() );
	kChunk++;
	if( kChunk >= nChunks ) break;
      }
      if( kChunk >= nChunks ) break;
    }

    std::cout << "List of the g-files requested :" << std::endl;
    std::cout << "-------------------------------" << std::endl;
    std::cout << std::endl;
    list<string>::iterator is;
    for( is=RmyFileName_.begin(); is!=RmyFileName_.end(); is++ )
    std::cout << (*is) << std::endl;
    std::cout << std::endl;

  }
