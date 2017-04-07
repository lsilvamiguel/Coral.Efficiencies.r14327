// $Id: CsInit.cc 14082 2015-10-29 15:58:42Z lsilva $

/*!
   \file    CsInit.cc
   \brief   Compass Initialization Class.
   \author  Benigno Gobbo
   \version $Revision: 14082 $
   \date    $Date: 2015-10-29 16:58:42 +0100 (Thu, 29 Oct 2015) $
*/

/// to enable strptime in time.h
#define _XOPEN_SOURCE 1

#include <unistd.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include "coral_config.h"
#include "CsInit.h"
#include "CsOpt.h"
#include "CsRandom.h"
#include "CsHistograms.h"
#include "CsEvent.h"
#include "CsZebraProto.h"
#include "CsHbookProto.h"
#if USE_NewEDIS
#  include "CsEvdis.h"
#endif
#include "CsRegistrySing.h"
#include "CDB.h"

#include "FileDB.h"
#if USE_MySQL
#  include "MySQLDB.h"
#  define DBSERVER "lxfs1657.cern.ch"
#  define DBUSER "anonymous"
#  define DBNAME "runlb"
#endif

#include "CoralRelease.h"
#include "Reco/DataBase.h"
#include "CsStopwatch.h"

using namespace std;

using Reco::DataBase;
DataBase *data_base=NULL;

extern "C" int inithbook_( int& );

CsInit* CsInit::instance_ = NULL;

int CsInit::zebra_[] = {0};

CsInit* CsInit::Instance() {
 if( instance_ != 0 )
   return instance_;
 else {
   cerr << "CsInit FATAL: wrong singleton usage." << endl;
   exit(1);
 }
}


CsInit* CsInit::Instance( int argc, char **argv ) {
 if( instance_ == 0 ) {
   instance_ = new CsInit( argc, argv );
   // Initialization of other singletons:
   CsEvent* event = CsEvent::Instance();
 }
 return instance_;
}

CsInit::CsInit( int argc, char **argv )  : 
  startOfRun_(0), endOfRun_(0), // avoid those in any case!
  fromDB_(false),  // avoid uninitialised variable for MC jobs
  updateSolenoidField_(false) // do not update solenoid field by default
{
  CsStopwatch chronos; int chrono = chronos.start();

  string tag, key, str; int n; list<int> ln; list<string> ls;
  CsOpt *opt = CsOpt::Instance(argc,argv); // Instantiate the Option Interpreter
  CsRegistrySing::Instance();              // Instantiate the Registry Singleton
  CsRandom::Instance();   // Instantiate the Random number Generator
  CsErrLog::Instance();   // Instantiate the Error Logger

  // Clear...
  mcJob_ = false; dtJob_ = false; hadronJob_ = false; BCSJob_ = false;
  detTable_.erase();  pitchTable_.erase();


  //       ******************** DATA RUN ********************
  tag = "Data"; key = "job";
  if( opt->getOpt( tag, key ) ) dtJob_ = true;
  key = "year";
  if( opt->getOpt( tag, key, n ) ) year_ = n;
  key = "period";
  if( opt->getOpt( tag, key, str ) ) period_ = str;
  key = "type";
  if( opt->getOpt( tag, key, str ) ) dataType_ = str;
  key = "run select";
  if( opt->getOpt( tag, key, ln ) ) {
    list<int>::iterator li;
    for( li=ln.begin(); li!=ln.end(); li++ ) {
      runs_.push_back( (*li) );
    }
  }
  key = "run skip";
  if( opt->getOpt( tag, key, ln ) ) {
    list<int>::iterator li;
    for( li=ln.begin(); li!=ln.end(); li++ ) {
      runs_.remove( (*li) );
    }
  }
  key = "container";
  if( opt->getOpt( tag, key, str ) ) {
    container_ = str;
  }
  else {
    container_.erase();
  }
  key = "file";
  while( opt->getOptRec( tag, key, str ) ) {
    dtFiles_.clear(); // to save only last "Data file" option.
    dtFiles_.push_back( str );
  }
  key = "files";
  while( opt->getOptRec( tag, key, str ) ) {
    dtFiles_.push_back( str ); // Read in several files recursively
  }

  key = "parallel reconstruction";
  if( opt->getOpt( tag, key ) ) {
    parallel_ = true;
  }
  else {
    parallel_ = false;
  }

  key = "store reconstructed events on DB";
  if( opt->getOpt( tag, key ) ) {
    saveRecoEvents_ = true;
  }
  else {
    saveRecoEvents_ = false;
  }


  //     ******************** MONTE CARLO JOB ********************
  tag = "Monte Carlo";
  key = "job";
  if( opt->getOpt( tag, key ) ) mcJob_ = true;
  key = "file";
  while( opt->getOptRec( tag, key, str ) ) {
    mcFiles_.push_back( str );
  }
  key = "ntfile";
  while( opt->getOptRec( tag, key, str ) ) {
    mcNtFiles_.push_back( str );
  }
  if (mcJob_) {  // ***** RESET RANDOM SEED EVERY NEW EVENT *****
    if (opt->getOpt("","reset random seed every new event"))
      // It can be useful, for debugging purposes, to reset the random at the
      // beginning of each event
      _resetRandomSeed = true;
    else
      _resetRandomSeed = false;
  }

  // Mickey Mouse MC
  tag = "Mickey";
  key = "hits";
  if( opt->getOpt( tag, key ) ) {
    _mickeyhits = true;
  }
  else {
    _mickeyhits = false;
  }
  key = "all";
  if( opt->getOpt( tag, key ) ) {
    _mickeyall = true;
  }
  else {
    _mickeyall = false;
  }
  key = "onreco";
  if( opt->getOpt( tag, key ) ) {
    _mickeyonreco = true;
  }
  else {
    _mickeyonreco = false;
  }
  key = "number of events";
  if( opt->getOpt( tag, key, n ) ) {
    _mickeynevents = n;
  }
  else {
    _mickeynevents = 100;   // default....
  }

  //   ******************** EVENTS TO READ, SKIP,... ********************
  tag = "";  key = "events to read";
  if( opt->getOpt( tag, key, n ) ) {
    maxEvents_ = (unsigned int) n;
  }
  else {
    maxEvents_ = 0;
  }
  key = "events to skip";
  if( opt->getOpt( tag, key, n ) ) {
    skipEvents_ = (unsigned int) n;
  }
  else {
    skipEvents_ = 0;
  }
  key = "tolerated spate of errors";
  if (opt->getOpt(tag,key,n)) {
    maxConsecutiveSkips_ = (unsigned int) n;
  }
  else maxConsecutiveSkips_ = 0;

  // Skipping first part of the spill
  if (opt->getOpt("events","BOS_skip",minTimeInSpill_)) {
    CsErrLog::msg(elWarning,__FILE__,__LINE__,
		  "Skipping first %f s in spill",minTimeInSpill_);
  }
  else minTimeInSpill_ = -1;

  //   ******************** DETECTOR TABLE ********************
  tag = ""; key = "detector table";
  if( opt->getOpt( tag, key, str ) ) detTable_ = str;

  if( detTable_.empty() ) {
    cerr << "CsInit FATAL: no 'detector table' pointer available." << endl
	 << "You must supply this option." << endl;
    exit(1);
  }

  //           ********** VARIABLE-SIZED PITCH TABLE **********
  tag = ""; key = "VarPitch table";
  if( opt->getOpt( tag, key, str ) ) pitchTable_ = str;

  if( pitchTable_.empty() ) {
    cerr << "CsInit WARNING: no variable-sized pitch table file defined." << endl;
    varPitchTable_ = false;
  }
  else varPitchTable_ = true;

  //                 ********** TRIGGERS **********
  // - Read their description from option. It's typically described in:
  //       "$CORAL/src/pkopt/trigger[.<extension>].opt".
  // - Expected are entries like:
  //			name  mask
  //Trigger mask	I	1
  //Trigger mask	M	2 // M covering both MT and InclMT.
  //etc...
  // - They are otherwise, for real, as opposed to MC, data, described in the
  //  mapping. The two sources of info could be X-checked, once the run# is
  //  known: not yet done, as of 2010/07. 
  list<string> names;
  while (opt->getOptRec("Trigger","mask",names)) { // Recursive option...
    // ...can be specified several times: earlier entries being not overriden.
    if (names.size()!=2) CsErrLog::msg(elFatal,__FILE__,__LINE__,
      "Bad \"Trigger mask\" entry, w/ # of arguments = %d",names.size());
    unsigned int mask; istringstream(names.back().c_str())>>mask;
    string name = names.front(); TCSMasks_[mask] = name[0];
  }

  //   ******************** EVERYTHING'S OK FOR MC JOB? ********************
  bool useTGEANT =  false;
#if USE_TGEANT
  std::string path;
  if (CsOpt::Instance()->getOpt( "CsTGEANTFile", "file", path )) {
    useTGEANT = true;
  }
#endif

  if (mcJob_  && !useTGEANT) {
    if( mcFiles_.empty() && mcNtFiles_.empty() ) {
      cerr << "CsInit FATAL: 'Monte Carlo job' set but no MC files specified.\n";
      exit(1);
    }
    if( !mcFiles_.empty() ) {
      list <string>::iterator Is;
      for (Is=mcFiles_.begin(); Is!=mcFiles_.end(); Is++ ) {
	mcFilesPtr_.push_back( &(*Is) );
      }
    }
    else if( !mcNtFiles_.empty() ) {
      list <string>::iterator Is;
      for (Is=mcNtFiles_.begin(); Is!=mcNtFiles_.end(); Is++ ) {
	mcNtFilesPtr_.push_back( &(*Is) );
      }
    }
  }
  else{ // if not mcJob_ clear lists
    mcFiles_.clear();
    mcFilesPtr_.clear();
    mcNtFiles_.clear();
    mcNtFilesPtr_.clear();
  }

  if( !mcJob_ ) {
    if( _mickeyall && !_mickeyonreco ) {
      mcJob_ = true;
    }
  }
  else {
    if( _mickeyall ) {
      cerr << "CsInit FATAL: 'Monte Carlo job and Mickey Mouse job both set."
	   << endl;
      exit(1);
    }
  }

  //   ******************** EVERYTHING'S OK FOR DATA JOB? ********************
  if( dtJob_ ) { // if dtJob_ ask for Det. tab, at least a Run, etc...
    if( ( runs_.empty() && container_.empty() )
	|| year_ == 0 ) {
      if( dtFiles_.empty() ) {
	cerr << "CsInit FATAL: 'data job' set but not enough options set."
	     << endl;
	exit(1);
      }
      else {
	fromDB_ = false;
	list <string>::iterator Is;
	for (Is=dtFiles_.begin(); Is!=dtFiles_.end(); Is++ ) {
	  dtFilesPtr_.push_back( &(*Is) );
	}
      }
    }
    else {
      if( dtFiles_.empty() ) {
	fromDB_ = true;
	runs_.sort();
	runs_.unique();
	if( !container_.empty() && runs_.size() > 1 ) {
	  cerr <<
	    "CsInit FATAL: A container and more than one run specified."
	       << endl;
	  exit(1);
	}
	if( !container_.empty() && parallel_ ) {
	  cerr <<
	    "CsInit FATAL: A container and Parallel processing set."
	       << endl;
	  exit(1);
	}
	if( runs_.size() != 1 && parallel_ ) {
	  cerr <<
	    "CsInit FATAL: More than a Run and Parallel processing set."
	       << endl;
	  exit(1);
	}
      }
      else {
	cerr << "CsInit FATAL: both DB and DATE readout mode set." << endl;
	exit(1);
      }
    }
  }
  else{ // if not dtJob_ clear lists etc.
    dtFiles_.clear();
    runs_.clear();
    period_ = "";
    year_ = 0;
  }


  if( dtJob_ && mcJob_ ) { // do not allow mixed access!!!
    cerr << "CsInit FATAL: both 'Data job' and 'Monte Carlo job' set." << endl;
    exit(1);
  }

  if( opt->getOpt( "", "DataBase", str ) )
  {
    cerr<<"\n\nERROR !! Old calibration database system is outdated !!\n\n";
    cerr<<"  Please use CDB instead of DataBase in option file\n";
    cerr<<"    see opt.rd.2002 for example\n\n"<<endl;
    exit(1);
//     assert(data_base==NULL);
//     data_base = new DataBase(str,DataBase::WRITE|DataBase::READ|DataBase::CREATE);
  }

  // Inizialize zebra if MCJob or USE_HBOOK and HBOOK histo package set
  tag = "";
  key = "histograms package";
  string histoPackage;
  opt->getOpt( tag, key, histoPackage );

  initzebradone_ = false;
  if( mcJob_ ) {
    initZebra();
    int nzebra = nZebra_;
    if( ! inithbook_( nzebra ) ) exit(1);
    initzebradone_ = true;
  }
  else {
#ifdef USE_HBOOK
    if( histoPackage=="HBOOK" ) {
      int zebra = 0;
      if( ! inithbook_( zebra ) ) exit(1);
      initzebradone_ = true;
    }
#endif
  }

  //   ******************** HISTOGRAMS PACKAGE ********************
  CsHistograms::Init();

  //   ******************** EVENT DISPLAY ********************
#if USE_NewEDIS
  // use event display package?
  tag = "";
  key = "event display";
  bool useEvDisplay(false);
  if (opt->getOpt(tag, key))
    useEvDisplay = true;

  // Initialize Event Display package
  if( useEvDisplay ) {
    if (!CsOpt::Instance()->getOpt("CsROOTGeometry", "file"))
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
                    "The ROOT based event display can only be used if ROOT geometries are used, not with material maps.");
    CsEvdis::Init();
  }
#endif

  //          ***** SOURCE of CALIBRATION DB *****

  {
    // we need to set the cdbSwitch_ even when calibration is not used to
    // avoid invalid memory access in CsDet::CsDet()
    string cdbSwitch;
    opt->getOpt("","CDB use",cdbSwitch);
    if      (cdbSwitch=="ConditionsDB") cdbSwitch_ = USE_ConditionsDB;
    else if (cdbSwitch== "FileDB")      cdbSwitch_ = USE_FileDB;
    else if (cdbSwitch== "MySQLDB")     cdbSwitch_ = USE_MySQLDB;
    else                                cdbSwitch_ = USE_NoCDB;
  }

  use_calibration = false;

  if(opt->getOpt( "use", "calibration")) {
    string cdbLocation, cdbUseTimeStr;
    opt->getOpt("","CDB location", cdbLocation);
    cdbUseTime = 0;
    if (opt->getOpt("","CDB usetime", cdbUseTimeStr)) {
      struct tm stm;
      strptime(cdbUseTimeStr.c_str(), "%Y-%m-%d-%H:%M:%S", &stm);
      cdbUseTime = mktime(&stm);
      if (cdbUseTime<0) cdbUseTime = 0;
    }

    if      (cdbSwitch_==USE_ConditionsDB) {
      // cdbLocation is boot file
    }
    else if (cdbSwitch_==USE_FileDB) {
      _cdb = new FileDB(cdbLocation);
      use_calibration = true;
    }
    else if (cdbSwitch_==USE_MySQLDB) {
#if USE_MySQL
      string cdbserver, cdbusername, cdbuserpasswd, cdbdbname, entrytime;
      string cdbspecialplace;
      int cdbportnumber;
      if (!opt->getOpt("", "CDB server", cdbserver)) {
        cdbserver = DBSERVER;
      }
      if (!opt->getOpt("", "CDB username", cdbusername)) {
        cdbusername = DBUSER;
      }
      if (!opt->getOpt("", "CDB userpasswd", cdbuserpasswd)) {
        cdbuserpasswd = "";
      }
      if (!opt->getOpt("", "CDB dbname", cdbdbname)) {
        cdbdbname = DBNAME;
      }

      MySQLDB* _mysqldb = new MySQLDB(cdbserver.c_str(), cdbusername.c_str(),
                                      cdbuserpasswd.c_str(), cdbdbname.c_str());
      _cdb = _mysqldb;
      if (opt->getOpt("", "CDB entrytime", entrytime)) {
        string entrymysqltime = MySQLInterface::toMySQLtime(entrytime.c_str());
        _mysqldb->setEntryTime(entrymysqltime);
        cout<<"CsInit::CsInit: MySQLDB will use entries before "<<_mysqldb->getEntryTime()<<endl;
      }
      if (opt->getOpt("", "CDB specialplace", cdbspecialplace)) {
        _mysqldb->setSpecialPlace(cdbspecialplace);
        cout<<"CsInit::CsInit: MySQLDB will take calib files in special places ";
        cout<<"keyed as "<<cdbspecialplace<<" in database"<<endl;
      }
      if (opt->getOpt("", "CDB portnumber", cdbportnumber)) {
        _mysqldb->setNumPort(cdbportnumber);
	CsErrLog::
	  msg(elWarning,__FILE__,__LINE__,
	      "MySQLDB will use port number \"%d\" to connect to the server",
	      cdbportnumber);
      }
      use_calibration = true;
#else
      cerr<<"CsInit::CsInit => Error: can't use MySQLDB calibration DB without MySQL\n";
      cerr<<"  Please reconfigure Coral with --with-MySQL option to use it"<<endl;
      exit(1);
#endif
    }
  }

  cout << endl
       <<coralLogo[0] << endl
       <<coralLogo[1] << endl
       <<coralLogo[2] << endl
       <<coralLogo[3] << endl
       <<coralLogo[4] << endl
       <<coralLogo[5] << endl
       << endl
       << "Release: " << CORAL_MAJOR_RELEASE
       << "."         << CORAL_MINOR_RELEASE
       << "."         << CORAL_BUILD
       << endl << endl;

  cout << "Total initialization time: "
       << chronos.stop( chrono )
       << " s."
       << endl;

  if (dtJob_) {

    // ******************** RD: SELECT "DETECTORS.DAT" ********************

    // - "detectors.dat" can be specified in options files via
    //   i) either a definite file,
    //  ii) or a directory.
    // - In case (ii), the most appropriate file is selected on view of run#,
    //  and some other characteristics (target polarity only considered so far).
    // - Run# taken from "runs_" or extracted from "dtFiles_" names.

    if (!runs_.empty() && !dtFiles_.empty())
      CsErrLog::mes(elFatal,
		    "Both \"Data file\" and \"Data runs select\" options "
		    "are specified. Choose one of them!");
    if (runs_.empty() && dtFiles_.empty())
      CsErrLog::mes(elFatal,
		    "None of \"Data file\" and \"Data runs select\" options "
		    "is specified. Choose one of them!");

    int runNum =                                              // ********** RUN#
      !runs_.empty() ? runs_.front() : -1;       // FROM LIST of RUNS, if any...
    if (!dtFiles_.empty()) {       // ...else (overriding) BY <run#> in FILENAME
      int ExtractRunNumber(string fn); // Prototype
      runNum = ExtractRunNumber(dtFiles_.front());
    }

    //                                                 ***** GET "DETECTORS.DAT"
    //  Depending possibly on target field polarity...
    //  (N.B.: The target's polarity is retrieved a 2nd time, later on, when it
    //  comes to re-scale the target's field map, w/ the same algorithm, cf.
    //  "CsMagInfo::readPolarization", for no other reason than bad programming.
    // ...or upon beam charge polarity.
    if      (CsOpt::Instance()->getOpt("hadron","run")) {
      //            ********** HADRON JOB
      hadronJob_ = true; // Flag "hadronJob_"
      detTable_ = GetDetectorsDat(detTable_,runNum,1);
    }
    else if (CsOpt::Instance()->getOpt("BCS","run")) {
      //            ********** HADRON JOB
      BCSJob_ = true;    // Flag "BCSJob_"
      detTable_ = GetDetectorsDat(detTable_,runNum,2);
    }
    else
      detTable_ = GetDetectorsDat(detTable_,runNum);

    // Updating solenoid field enabled? This is only useful when processing data
    // taken during field rotation. And we  would like it to be enabled only
    // then. Unfortunately there is so far (as of 05/05) no method (asking,
    //  e.g., MySQL and) returning the info, on a run by run basis => In order
    // not to risk wasting resources, let's condition the updating by option.
    if (CsOpt::Instance()->getOpt("CsMagInfo","Update solenoid")) {
#if USE_MySQL
      updateSolenoidField_ = true;
#else
      CsErrLog::mes(elFatal,"Requesting \"Upadate solenoid\" while MySQL off!");
#endif
    }
    updateSolenoidField_ = false;

#if USE_MySQL
    // This is an example how to read SM2 NMR values measured by E. Weise
    // The getSM2NMR(runnumber) method calculates an average of the measurements
    //  over the run
//     double sm2nmr = 0;
//     if (useMySQLDB()) {
//       MySQLDB* _mysqldb = (MySQLDB*) getDB();
//       bool connectfg = _mysqldb->isConnected();
//       if (!connectfg) _mysqldb->ConnectDB();
//       sm2nmr = _mysqldb->getSM2NMR(runNum);
//       if (!connectfg) _mysqldb->DisconnectDB();
// 
//       if (sm2nmr > 0) cout << "SM2 NMR average value "<<sm2nmr<<endl;
//       if (sm2nmr == 0) cout << "SM2 NMR average value not found"<<endl;
//       if (sm2nmr < 0) cout << "Error when reading SM2 NMR average value"<<endl;
//     }

    // This is an example how to get the shift path to the 1st event file.
    // This file is not present for all runs.
//     string firstevtpath("");
//     if (useMySQLDB()) {
//       char* pathstr;
//       MySQLDB* _mysqldb = (MySQLDB*) getDB();
//       bool connectfg = _mysqldb->isConnected();
//       if (!connectfg) _mysqldb->ConnectDB();
//       pathstr = _mysqldb->get1stEvtPath(runNum);
//       if (pathstr) firstevtpath = pathstr;
//       if (!connectfg) _mysqldb->DisconnectDB();
// 
//       if (firstevtpath != "") cout << "First event file path: "<<firstevtpath<<endl;
//       else cout << "The first event file for run "<<runNum<<" does not yet exist..."<<endl;
//     }

#endif


  }// end of if(dtJob_)

  if (mcJob_) { // In MC case, check MC updating solenoid field is not requested
    if (CsOpt::Instance()->getOpt("CsMagInfo","Update solenoid"))
      CsErrLog::mes(elFatal,"Requesting \"Upadate solenoid\" for MC data!");
  }

} // End of CsInit constructor

void CsInit::initZebra() {

  // init zebra
  mzebra(-1);
  int ixstore = 0;
  mzstor( ixstore, "z", " ", &zebra_[5], &zebra_[10],
	  &zebra_[11], &zebra_[11], &zebra_[2000],
	  &zebra_[nZebra_-1], 1, 1 );
}

int CsInit::getCoralMajorRelease() {
  return( CORAL_MAJOR_RELEASE );
}

int CsInit::getCoralMinorRelease() {
  return( CORAL_MINOR_RELEASE );
}

int CsInit::getCoralBuild() {
  return( CORAL_BUILD );
}

inline void dumpOptions( int dumplevel ) {
  CsOpt::Instance()->dump( dumplevel);
}

//  **********************************************************************
//  ***************   EXTRACT RUN NUMBER FROM FILE NAME    ***************
// Expecting name to be
// - Either <anything>-<run#>.raw, e.g. data file on castor =
//      cdr<eb#><chunk#>-<run#>.raw
// - Or calib/random file =
//      cdrpccoeb<eb#>-<run#>-<mask>.<chunk#>.raw
// Returning 0 if run# could not be found.
//  **********************************************************************
int ExtractRunNumber(string fnam)
{
  int pos = 0, npos = 0; // Character pos./#pos delimiting <run#> field
  int dot, dot2, dash, dash2; // Character pos. of intermediate fields
  if ((dot  = fnam.rfind('.'))>0 &&           // Dot found (not @ beginning)
      (dash = fnam.rfind('-',dot-1))>=0) {    // Dash found upstream of dot
    if ((dot2 = fnam.rfind('.',dot-1))<dash) {// No intervening dot: data file..
      pos = dash+1; npos = dot-dash-1;
    }
    else if ((dash  = fnam.rfind('-',dot2-1))>0 && // ...else look for "-...-"
	     (dash2 = fnam.rfind('-',dash-1))>=0) {
      pos = dash2+1; npos = dash-dash2-1;
    }
  }
  if (!pos) {
    printf("ExtractRunNumber\a ==> No run# field in input data file name:\n"
	   "     \"%s\"\n",fnam.c_str());
    printf("   Expected syntax: either \"<anything>-<run#>.raw\" or \"<anything>-<run#>-<mask>.<chunk#>.raw\".\n");
    printf("   Returning run# = 0!\n");
    return 0; // 0 being a non physical run#
  }
  string runNumField(fnam,pos,npos);
  char *end, **endptr = &end;       // Check run# field houses a numerical value
  int num = strtol(runNumField.c_str(),endptr,0);
  if (**endptr!='\0') {
    printf("ExtractRunNumber\a ==> Non-numerical value found in run# field of input data file name:\n"
	   "     \"%s\"\n",fnam.c_str());
    printf("   Assuming syntax is either \"<anything>-<run#>.raw\" or \"<anything>-<run#>-<mask>.<chunk#>.raw\".\n");
    printf("   Returning run# = 0!\n");
    return 0; // 0 being a non physical run#
  }
  else return num;
}

//  **************************************************************************
//  *****  FIND MOST RELEVANT detectors.dat FOR GIVEN RUN and JOB TYPE  ******
// - Among all detectors.dat in given path, w/ <run#> and <type> extensions.
// - MOST RELEVANT is CLOSEST EARLIER <run#> 
//                 w/ MAGNETS SETTING RELEVANT for <type> 
// - JOB TYPE is:
//   0: spin asymmetry:   relevant depends on target (solenoid,dipole) field
//                   - (+/-1,0) <-> <type> = "plus" or "minus" 
//                   - (  0 ,1) <-> <type> = "transv"
//   1: hadron                      <type> = "hadron"
//   2: charge asymmetry: relevant depends on SM1/2 polarity
//                   - +/-1     <-> <type> = "mu+" or "mu-"
//  **********************************************************************
#include <dirent.h>
string CsInit::GetDetectorsDat(string path, // Path to either file or directory
			       int run,
			       int jobType)
{
  cout<<endl;
  //                                          ***** IS "path" DIRECTORY OR FILE?
  if (path[path.length()-1]!='/') { // Path not ending w/ "/"
    cout<<"GetDetectorsDat ==> \""<<path<<"\" will be used.\n\n";
    return path;                           // ...=> FILE RETURN INPUT "path" ARG
  }
  else {                                   // ...ELSE DIRECTORY
    if (run<0) {
      cout<<"GetDetectorsDat ==> Run number can't be extracted from input data's file name\n";
      cout<<"Please check your options file (cf. tips for helping yourself out in the file's header)\n";
      exit(1);
    }
  }

  printf  ("GetDetectorsDat ==> Run# to be processed "
	   "(as was 'guessed' by parsing input file name) is %d\n\n",run);
  // Now that the case of FILE has been dealt w/, get info about magnets.
  int solSign = 0, dipole = 0, smSign;
  if      (jobType==2) GetMagnetInfo(run,0x1,solSign,dipole,smSign);
  else if (jobType==0) GetMagnetInfo(run,0x2,solSign,dipole,smSign);
  if      (jobType==0) {
    if (solSign)
      printf("GetDetectorsDat ==> This run had solenoid field sign : \"%c\"\n",
	     solSign==+1 ? '+' : '-');
    else if (dipole)
      printf("GetDetectorsDat ==> This run had dipole field ON\n");
    else
      printf("GetDetectorsDat ==> This run had target field OFF\n");
  }
  else if (jobType==2)
    printf("GetDetectorsDat ==> This run had SM1/2 polarity %c0\n",
	   smSign>0?'>':'<');

  string targetedType;
  if      (jobType==0) {
    if      (solSign==+1) targetedType = "plus";
    else if (solSign==-1) targetedType = "minus";
    else if (dipole)      targetedType = "transv";
    else                  targetedType = "zero";
  }
  else if (jobType==1)    targetedType = "hadron";
  else if (jobType==2) {
    if      (smSign==+1)  targetedType = "mu+";
    else if (smSign==-1)  targetedType = "mu-";
  }

  DIR *dp; if (!(dp = opendir(path.c_str())))         // ***** OPENING DIRECTORY
	     CsErrLog::msg(elFatal,__FILE__,__LINE__,
			   "Can't open directory \"%s\"",path.c_str());
  struct dirent *ep; map<int,string> mRun2File;      // ***** PARSING FILE NAMES
  while ((ep = readdir(dp))) { // Loop over file names
    string fnam(ep->d_name);
    //                  ***** REQUIRE "detectors[.<anything>].<run#>.<type>.dat"
    int dot1; if (fnam.find("detectors.")!=0 ||
		  (dot1 = fnam.find(".dat"))!=(int)fnam.length()-4) continue;
    int dot2 = fnam.rfind('.', dot1-1); if (dot2<0) continue;
    int dot3 = fnam.rfind('.', dot2-1); if (dot3<0) continue;

    string type(fnam,dot2+1,dot1-dot2-1);    // ***** FILTER TYPE w/ TARGET INFO
    if (type!=targetedType) continue;

    string runNumField(fnam,dot3+1,dot2-dot3-1);           // ***** EXTRACT RUN#
    char *end, **endptr = &end;     // Check run# field houses a numerical value
    int num = strtol(runNumField.c_str(),endptr,0);
    if (**endptr!='\0') continue;

    mRun2File[num] = fnam;                       // ***** MAP RUN# <-> FILE NAME
  } // End of loop over file names
  closedir (dp);
  if (mRun2File.empty())
    CsErrLog::msg(elFatal,__FILE__,__LINE__,
      "No files w/ syntax \"detectors[.<anything>].<run#>.%s.dat\""
		  "in directory \"%s\"",targetedType.c_str(),path.c_str());

  // Files w/ well formed syntax are expected to be ordered by increasing run#. 
  int bestRun; map<int,string>::reverse_iterator irun;
  for (irun = mRun2File.rbegin(), bestRun = -1; irun!=mRun2File.rend();
       irun++) {
    if ((*irun).first>run) continue;
    bestRun = (*irun).first; break;
  }

  if (bestRun!=-1)
    printf("GetDetectorsDat ==> \"%s%s\" will be used to process it.\n",
	   path.c_str(),mRun2File[bestRun].c_str());
  else {
    printf("GetDetectorsDat ==> Files w/ syntax\n"
	   "   \"detectors[.<anything>].<run#>.%s.dat\" in \"%s\" are:\n",
	   targetedType.c_str(),path.c_str());
    for (irun = mRun2File.rbegin(); irun!=mRun2File.rend(); irun++)
      printf("\"%s\" ",(*irun).second.c_str());
    cout<<endl;
    CsErrLog::mes(elFatal,"No appropriate \"detectors.dat\" could be found.\n"
      "=> Check your options file. In last resort, "
      "you can always specify an explicit \"detectors.dat\" file name.");
  }
  cout<<endl;
  return path+mRun2File[bestRun];
}

//------------------------------------------------------------------------------

void CsInit::GetMagnetInfo(int Run,
			   unsigned int ignore, // 0x1: Target, 0x2: SM1/2
			   int &solSign,        // +/-1, 0 or -2: not available
			   int &dipole,         // 1 or 0 or -2: not available
			   int &smSign)         // +/-1 or -2: not available
{
  // - Get sign of solenoid and SM1/2 fields, and ON/OFF status of dipole.
  // - From either MySQLDB or file, depending upon option.
  // - The routine exits w/ fatal error if:
  //   - Info requested cannot be accessed (be it from MySQLDB or File).
  //   - MySQLDB is used, as opposed to File (the idea being that File can be
  //    used as a work-around method in case MySQL triggers a fatal error and
  //    cannot be cured), and...
  //     ...Both solenoid and dipole =0, because it's a combination not expected
  //     in polarized muon data taking. (=> Avoid "GetMagnetInfo" or book
  //     "ignore==1" while MySQLDB if your data are of a hadron, or otherwise
  //     unpolarized, type.)
  //     ...Or solenoid =0 alone in 2002,2003 where dipole is not in DB. (=>
  //     Avoid "GetMagnetInfo" while MySQLDB for 2002,2003 transversity.
  //   - In all cases, if inconsistencies: both solenoid and dipole, SM1 and
  //    SM2 of opposite polarity, off, w/ unexpected values of current.

  // - As to target, the raw info to be retrieved is the value of current(s).
  // - If MySQLDB, the sign of the solenoid is:
  //   - If < 2004: from a query to "tb_offlinepolar".
  //   - If >=2004: from a query to "tb_tgtCurrents".
  //  (Cf. "../geom/CsField.cc:readPolarization" for further details.)
  // - In any  case, the absolute value of the current is checked against a cut,
  //  cf. "SolCurrentCut" infra: all values below cut are considered = 0. The
  //  idea being simply to single out the cases where the current compatible w/
  //  0 taking into account the finite precision of the measurement of the
  //  current.
  //   This cut was originally set = 300(A). Which turned out to be too high for
  //  the 2006 case, where for a 1T setting, |current| ~= 260 A.  I(Y.B.) decide
  //  to decrease it to some lower value => Question is: what lower value?
  //    I) It has to be < 260 A. Therefore something like 200 A would do, since
  //      there is in principle no data taken @ intermediate values: it's
  //      either all or nothing (in 2002-2004), or 1T or 2.5T or nothing (2006).
  //   II) Then there is the case of field rotation runs. To cover that case, it
  //      has been decided, cf. "CsMagInfo::scaleSolenoidMap", that, upon
  //      option ("CsMagInfo Update solenoid"), the target field is scaled
  //      every so often, w/ the, algebraic, ratio of the value of the current,
  //      retrieved using same query as the one used here but via a different
  //      class, over that associated to the map file. => So that it would not
  //      really matter what polarity is returned here, except that it is
  //      used to select between the ".plus" and ".minus" versions of the
  //      "detectors.dat", cf. "CsInit::CsInit". The bottom line then is: do the
  //      ".plus" and ".minus" differ in any other way than the solenoid scale
  //      they contain (entry "mag 3", 6th numerical value)? The answer is:
  //         - (a) In generality, "maybe",
  //         - In practice:
  //           - (b) 2003 on: "no" , since Marcin introduced the idea of tuning
  //            the alignment simultaneously on data from both polarities,
  //           - (c) 2002: "yes", since alignment does depend on polarity. But
  //            then so does it depend on the exact instantenous value of the
  //            field, and only an interpolation between ".plus" and ".minus"
  //            would be accurate, which our present scheme does not provide
  //            for. In the absence of interpolation, the best choice is still
  //            to guess the correct polarity. Except for low enough field where
  //            it does not really matter.
  //   => Conclusion:
  //      i) The numerical value of the cut only matters in very special cases.
  //     ii) It has to be low enough in regard to case (II.c).
  //    iii) It has to be larger than the precision.
  //         => I opt for "SolCurrentCut = 50".
  // (Note: the same "SolCurrentCut" is applied to the dipole case.)

  const float SolCurrentCut = 50;

  // - Last, what we're interested in is the sign of the field, not that of the
  //  current. Conversion:
  //  2002-2004, i.e Run < 45064 = 1st 06T01 run: field has OPPOSITE sign,
  //  2006...? :                                  field has SAME sign as current

  int sign, current2Field = Run<45064 ? -1 : +1;

  // Uncertainty on current. For SMC, I(Y.B.) determined it from the
  // deviation from nominal observed in DB = .08A. Assuming same precision
  // for the COMPASS/OD case and rounding up =>
  float soldCur_ = .1;

  float solCur = 0, targCurStatus = -1, dipCur = 0, sm1Cur = 0, sm2Cur = 0;
  unsigned int error = 0;
  bool fromMySQL = CsOpt::Instance()->getOpt("CsMagInfo","MySQLDB");
  if (fromMySQL) {

    //       ******************** MySQL ********************
#if USE_MySQL
    if (!useMySQLDB()) CsErrLog::mes(elFatal,
      "Can't read MySQLDB offline polarization without MySQLDB on!\n"
      "  => Book \"CDB use MySQLDB\" in your options file.");
    pair<double,double> *tgtCurrents = 0;
    MySQLDB *_mysqldb = (MySQLDB*)getDB();
    bool connectfg = _mysqldb->isConnected();
    if (!connectfg) {
      _mysqldb->ConnectDB();
      if (!_mysqldb->isConnected())  CsErrLog::mes(elFatal,
	"No connecting to MySQLDB while trying to retrieve magnet info.");
    }

    if (!(ignore&0x1)) {                       // ***** QUERY TARGET CURRENTS...
      // ...both from "tb_offlinePolar" and "tb_tgtcurrents" if 2004 on.
      map<string,double> offpol = _mysqldb->getTgtOfflinePolar(Run);
      if (Run>=33000) {  // I.e. if later or equal 2004
	tgtCurrents = new pair<double,double>;
	*tgtCurrents = _mysqldb->getTgtCurrents(Run);
      }
      //                         ***** PARSE TARGET CURRENTS for INCONSISTENCIES
      solCur = offpol["solenoid"];
      targCurStatus = offpol["status"]; // Meaning of status: 0 no value, 1 offline value, 2 online value, 3 old online value
      bool targCurError =
	targCurStatus==-1 || targCurStatus==0 || targCurStatus==3;

      if (Run<33000) {// If earlier than 2004: "tb_offlinepolar" only available
	if (targCurError) {
	  printf("\nGetMagnetInfo ==> run %d:"
  "No reliable target solenoid current found in MySQLDB's \"tb_offlinepolar\".",
		 Run);
	  error |= 0x1; // 0x1: error affecting target info
	}
      }
      else {          // If later or equal 2004: use "tb_tgtCurrents"...
	// Compare the 2 values, still. And warn if disagreement > 2*precision.
	if (fabs(tgtCurrents->first-solCur)>2*soldCur_)
	  CsErrLog::msg(elError,__FILE__,__LINE__,
	    "Run %d: Solenoid current in MySQL's \"tb_offlinepolar\" (=%.2f) "
	    "and \"tb_tgtcurrents\" (=%.2f) differ. The latter is retained.",
			Run,solCur,tgtCurrents->first);
	solCur = tgtCurrents->first;
	dipCur = tgtCurrents->second;
	delete tgtCurrents;
      }
    }
    if (!(ignore&0x2)) {                           // ***** QUERY SM1/2 CURRENTS
      pair<double,double> smCurrents = _mysqldb->getSMcurrents(Run);
      sm1Cur = smCurrents.first; sm2Cur = smCurrents.second;
    }
    if (!connectfg) _mysqldb->DisconnectDB(); // If initially connected: keep so

    if (!error)
      printf("\nGetMagnetInfo from MySQLDB\n  ==> run %d :",Run);
#else
    CsErrLog::mes(elFatal,
      "Can't read magnet info from MySQLDB without MySQL!"
      "  => Reconfigure Coral with --with-MySQL option "
      "or use option \"CsMagInfo File\".");
#endif
  }
  else {     // ******************** File ********************

    string path; if (!CsOpt::Instance()->getOpt("CsMagInfo","File",path))
      CsErrLog::mes(elFatal,
	"Option \"CsMagInfo\" nowhere specified in options file\n"
        "   => No guessing most appropriate \"detectors.dat\".\n"
        "   => Either supply this info via option \"CsMagInfo File <file>\"\n"
        "      Or specify explicit \"detectors.dat\" file in options file.");

    fstream in; in.open( path.c_str(), ios::in);
    if (!in ) CsErrLog::msg(elFatal,__FILE__,__LINE__,
			    "Error opening CsMagInfo file \"%s\"",path.c_str());

    // File is expected to have been created by:
    //    "$CORAL/scripts/magInfo.pl"
    // Let's check this is indeed the case.
    string ctmp, label; in>>ctmp>>label;
    if (ctmp!="/*" || (int)label.find("magInfo.pl")<0)
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
	"CsMagInfo file \"%s\" NOT of a \"magInfo.pl\" type!",path.c_str());

    int run; float soleno, dipole, dummy, sm1, sm2;
    bool found = false, comment = true;
    while (in.good() && !in.eof()) {
      in>>ctmp;
      if (comment) {
	if (ctmp=="*/") comment = false;
      }
      else if (ctmp=="/*") comment = true;
      else {
	run = atoi(ctmp.c_str());
	// One gets successively:
	//   type  soleno  dipole    up    +/-   down    +/-  target#  nmr   nmr#   sm1  sm2
	in>>dummy>>soleno>>dipole>>dummy>>dummy>>dummy>>dummy>>dummy>>dummy>>dummy>>sm1>>sm2;
	if (Run==run) {
	  solCur = soleno; dipCur = dipole; sm1Cur = sm1; sm2Cur = sm2;
	  found = true; break;
	}
      }
    }
    in.close();
    if (!found) CsErrLog::msg(elFatal,__FILE__,__LINE__,
      "No magnet info for run #%d in \"CsMagInfo File\" file \"%s\"",
			      Run,path.c_str());

    printf("\nGetMagnetInfo from File \"%s\"\n  ==> run %d: ",path.c_str(),Run);
  }

  //     *************** CHECK RETRIEVED VALUES ***************

  //                 ***** ARE THEY MEANUNGFUL ? *****
  // Currents: Large values (1e6,1e9) signal missing info, in MySQLDB (at least
  // in its "tb_beamvalues"), and in FileDB's produced by the "magInfo.pl" (at
  // least for target currents). => Let's checked for these large values.
  if (!(ignore&0x1) && solCur>1e5 && dipCur>1e5) {
    printf(" Target magnet info not available.");
    error |= 0x1;
  }
  if (!(ignore&0x2) && sm1Cur>1e5 || sm2Cur>1e5) {// Note: one OK while other...
    // ... KO, in MySQLDB, this can't be, cf. "MySQLDBInterface::getSMcurrents".
    printf(" SM1/2 magnets info not available.");
    error |= 0x2;
  }
  //           ***** IF MEANINGFULL: PRINT THEM to STDOUT *****
  if (error!=0x3) {
    if (!(ignore&0x1)) printf(" soleno = %.2fA, dipole = %.2fA",solCur,dipCur);
    if (!(ignore&0x2)) printf(" SM1 = %.2fA, SM2 = %.2fA",sm1Cur,sm2Cur);
  }
  printf("\n");

  if (!(ignore&0x1) && !(error&0x1)) {
    //                  ***** TARGET CONSISTENCY CHECKS *****
    //    ***** DETERMINE SIGN of SOLENOID FIELD, STATUS od DIPOLE *****
    if (solCur==0 && dipCur==0) {// Both values = 0 probably means error
      printf("Target solenoid and dipole currents = 0 "
	     "in MySQLDB's \"tb_tgtcurrent\": Inconsistency!\n");
      error |= 0x1;
    }
    if (!(error&0x1)) {
      if      (solCur<-SolCurrentCut) solSign = -1*current2Field;
      else if (solCur< SolCurrentCut) solSign =  0;
      else                            solSign = +1*current2Field;
      if (fabs(dipCur)>SolCurrentCut) dipole = 1;
      else                            dipole = 0;

      if (solSign && dipole) {
	printf("Both solenoid and dipole are ON. Suspicious!\n");
	error |= 0x1;
      }
    }
  }
  if (!(ignore&0x2) && !(error&0x2)) {
    //                 ***** SM1/SM2 CONSISTENCY CHECKS *****
    //             ***** DETERMINE SPECTROMETER POLARITY *****
    // - Polarity (i.e. "smSign"): used to select detectors.dat of "mu+/mu-"
    //  type in GPD data taking.
    // - Consistency checks. Are considered inconsistent:
    //   - Opposite polarities (in SM1 and SM2) 
    //   - SM1 XOR SM2
    //   - |SM1| < 2495, |SM2| < 3995
    //  Keeping in mind that what one wants here is to guess the most relevant
    //  det.dat and that one always has the possibility to supply an explicitly
    //  det.dat in order to get along w/ special cases.
    if (sm1Cur*sm2Cur<0) {     // SM1/2 of opposite sign.
      printf("SM1(=%.2fA), SM2(=%.2fA) w/ opposite sign: Suspicious!\n",
	     sm1Cur,sm2Cur);
      error |= 0x2;
    }
    else if (sm1Cur*sm2Cur==0) {// SM1 OR SM2 off, but not both: cf. supra.
      printf("SM1(=%.2fA) XOR SM2(=%.2fA): Suspicious!\n",
	     sm1Cur,sm2Cur);
      error |= 0x2;
    }
    else if (fabs(fabs(sm1Cur)-2500)>20 || // Non nominal values.
	     fabs(fabs(sm2Cur)-4000)>20 &&
	     fabs(fabs(sm2Cur)-5000)>20) {
      printf("SM1(=%.2fA), SM2(=%.2fA): Not nominal => Suspicious!\n",  
		  sm1Cur,sm2Cur);
      error |= 0x2;
    }
    else smSign = sm1Cur>0 ? +1 : -1;
  }

  if (error) {           // ***** EXIT UPON ERROR *****
    if (fromMySQL) CsErrLog::mes(elFatal,
        "  => No guessing most appropriate \"detectors.dat\".\n"
        "  => Either specify explicit \"detectors.dat\" file in options file.\n"
        "  => Or supply relevant info via option \"CsMagInfo File <file>\".");
    else {
      string path; CsOpt::Instance()->getOpt("CsMagInfo","File",path);
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
        "Bad \"CsMagInfo File\" file \"%s\"\n"
	"  => No guessing most appropriate \"detectors.dat\".",path.c_str());
    }
  }
}
