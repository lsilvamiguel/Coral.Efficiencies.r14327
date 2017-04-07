// $Id: CsEvent.cc,v 1.338 2011/03/01 01:20:34 ybedfer Exp $

/*!
   \file    CsEvent.cc
   \brief   Compass Event Class.
   \author  Benigno Gobbo
   \version $Revision: 1.338 $
   \date    $Date: 2011/03/01 01:20:34 $
*/

#include <algorithm>
#include <wordexp.h>
#include <stdlib.h>

#include "coral_config.h"
#include "CsInit.h"
#if USE_FileStore
#  include "CsFileStore.h"
#endif
//#if USE_COMPASS_Date
//#  include "CsDateStore.h"
//#endif
#  include "CsDummyStore.h"
#if USE_ORACLE
#  include "CsOraStore.h"
#  include "PatchEventHeaderDB.h"
#endif
#if USE_MySQL
#  include "MySQLDB.h"
#  include "CsPP.h"
#endif
#include "CsDDDStore.h"

#include "CsRandom.h"
#include "CsOpt.h"
#include "CsGeom.h"
#include "CDB.h"
#include "FileDB.h" 
#include "CsEvent.h"
#include "CsErrLog.h"
#include "CsRegistrySing.h"
#include <CLHEP/Matrix/Matrix.h>
#include "CsGeant3.h"

// TRAFFIC
#include "CsTrafficPrepattern.h"
#include "CsTrafficBridging.h"
#include "CsTrafficFitting.h"

// VERTEX
#include "CsAverPattern.h"
#include "CsRolandPattern.h"
#include "CsKalmanFitting.h"

#include "CsMCDigit.h"
#include "CsGauss.h"
#include "DaqDataDecoding/DaqEvent.h"
#include "DaqDataDecoding/DaqEventsManager.h"
#include "DaqDataDecoding/Chip.h"
#include "DaqDataDecoding/ChipF1.h"
#include "DaqDataDecoding/Scaler.h"
#include "DaqDataDecoding/TriggerTime.h"
#include "CsEventUtils.h"
#include "CsCalorimeter.h"
#include "CsCEDARDetector.h"
#include "Reco/CalorimeterParticle.h"

#include "CsHist.h"
#include "CsRichOne.h"
#include "CsParticle.h"
#include "CsBeamRecons.h"
#include "CsBeamReconstruction.h"
#include "CsRwRecons.h"
#include "CsRwChargeRecons.h"
#include "CsBuildParticles.h"

#include "CsStopwatch.h"
#include "CsRegistry.h"

#include "CsDetector.h"
#include "CsMiscDetector.h"


using namespace std;
using namespace CLHEP;


// DST data
#define DST_VERSION 4
#define DST_FATNESS 0

struct HodoDataStruct_type { 
  string tbname;
  int lower;
  int upper; 
} HodoDataStruct[24] = {
  { "HI04X1_d", -2500, -2430 },
  { "HI04X1_u", -2500, -2430 },
  { "HI05X1_d", -2550, -2450 },
  { "HI05X1_u", -2550, -2450 },
  { "HM04X1_d", -2550, -2450 },
  { "HM04X1_u", -2550, -2450 },
  { "HM05X1_d", -2600, -2500 },
  { "HM05X1_u", -2600, -2500 },
  { "HM05Y1_d", -2500, -2400 },
  { "HM05Y1_u", -2500, -2400 },
  { "HL04X1_m", -2250, -2050 },
  { "HL05X1_m", -2250, -2050 },
  { "HO04Y1_m", -2200, -1900 },
  { "HO04Y1_m", -2200, -1900 },
  { "HO04Y2_m", -2200, -1900 },
  { "HMSC1",    -2500,  -500 },
  { "VT01P1_I", -7800, -7000 },
  { "VT02P1_I", -7800, -7000 },
  { "VT01P1B1", -7800, -6900 },
  { "VT01X1_O", -7800, -6700 },
  { "VT01X1mO", -7800, -6400 },
  { "VTsum",    -7400, -6400 },
  { "VT01P1sf", -7600, -7400 },
  { "VT02P1sf", -7500, -7300 }
};

extern void PID_doMuonID( vector<CsParticle*>& parts );
extern void PID_doBeamID( vector<CsParticle*>& parts );
extern void muIDinMW1( vector<CsParticle*>& parts );
//extern void PID_doCalID( vector<CsParticle*>& parts );

// ------------To check time usages, Benigno 20030618---------------
#define CHECK_TIMEUSAGE 1

#if CHECK_TIMEUSAGE

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>

/* it's written in pure C to be used elsewhere... */

void timeusage( const char* str, int diff ) {

  struct rusage usage;
  double ctime, utime, stime;
  static double prevCtime, prevUtime = 0, prevStime = 0;
  struct timeval clock;
  char date[40], zone[10];
  char format1[] = "%a, %d/%b/%Y %T", format2[] = "(%Z)";

  printf( "+--------------------------------------------------------------------------+\n" );
  if( getrusage( RUSAGE_SELF, &usage ) == 0 ) {
    gettimeofday( &clock, NULL );
    strftime( date, 40, format1, gmtime( (time_t*) &(clock.tv_sec) ) );
    strftime( zone, 10, format2, gmtime( (time_t*) &(clock.tv_sec) ) );

    ctime = clock.tv_sec + clock.tv_usec/1000000.;
    utime = usage.ru_utime.tv_sec + usage.ru_utime.tv_usec/1000000.;
    stime = usage.ru_stime.tv_sec + usage.ru_stime.tv_usec/1000000.;
    printf( "%s (%s %s)\n", str, date, zone );
    printf( "+--------------------------------------------------------------------------+\n" );
    if( diff == 0 ) {
      printf( "Time since the Epoch : %16.4f s\n", ctime );
      printf( "CPU Time (user)      :     %12.4f s\n", utime );
      printf( "CPU Time (sys)       :     %12.4f s\n", stime );
      prevCtime = ctime;
      prevUtime = utime;
      prevStime = stime;
    }
    else if( diff == 1 ) {
      printf( "Time since the Epoch : %16.4f s, Difference:   %12.4f s\n", ctime, ctime-prevCtime );
      printf( "CPU Time (user)      :     %12.4f s, Difference:   %12.4f s\n", utime, utime-prevUtime );
      printf( "CPU Time (sys)       :     %12.4f s, Difference:   %12.4f s\n", stime, stime-prevStime );
    }
    else {
      printf( "ERROR: diff=%d, it can be only 0 or 1", diff );
    }
  }
  else {
    printf( "ERROR: getrusage failed.\n" );
  }
  
  printf( "+--------------------------------------------------------------------------+\n" );  
  
}
#endif // CHECK_TIMEUSAGE
// ---------------------------------------------------------------------


CsEvent* CsEvent::_instance = 0;

#ifdef TIME_CHECK
extern bool TiClust, TiDecod;

extern double clu_times[8]; //clustering time for given detector
extern int clu_nb[8];

extern double dec_times, dec1, dec2, dec3, dec4, dec5, dec6, dec7;  //decoding (digitalization) for event
extern int n_dig;

#endif

// timing stuff (local quantities)

static float _loadCalibTime = 0;
static float _eventLoadTime = 0;
static int   _neventLoadTime = 0;
static float _digitizTime = 0;
static int   _ndigitizTime = 0;
static float _clusterizeTime = 0;
static int   _nclusterizeTime = 0;
static float _beamTime = 0;
static int   _nbeamTime = 0;
static float _trackingTime = 0;
static int   _nTrackingTime = 0;
static float _caloTime = 0;
static int   _ncaloTime = 0;
static float _rich1Time = 0;
static int   _nrich1Time = 0;
static float _cedarTime = 0;
static int   _ncedarTime = 0;
static float _vertexTime = 0;
static int   _nvertexTime = 0;
static float _pidTime = 0;
static int   _npidTime = 0;
static float _partTime = 0;
static int   _npartTime = 0;
static float _dstDownTime = 0;
static int   _ndstDownTime = 0;
static float _dstUpTime = 0;
static int   _ndstUpTime = 0;
static float _finishTime = 0;
static int   _nfinishTime = 0;
static float _totalTime = 0;
static int   _ntotalTime = 0;

static int   _nSkippedEvents = 0;

int   CsEvent::_nspillmax = 200;

static CsStopwatch _chronometer(1);
static int _chrono = -1; 

CsEvent* CsEvent::Instance() {
  if( _instance == 0 ) {
    _instance = new CsEvent();
    CsRegistry reg;
    reg.EOJRegistration(_instance);
  }
  return( _instance );
}

// ***************************************************************************
// *************************  CsEvent INSTANTIATION  *************************
// ***************************************************************************
CsEvent::CsEvent() :
  _triggerTimeCorr(0.)  // set default for trigger time correction to 0, as in
                        // case of MC this value is never touched, but should
                        // still be initialized. For real data it is set for
                        // every event.
{

  CsInit *init = CsInit::Instance(); CsOpt *opt = CsOpt::Instance();

  // Oracle changing:
  //*******************
  const int ORACLE_STORAGE = 2;
  const int FILE_STORAGE   = 3;

  int storageSource = 0;
  string src = "";

  if(CsOpt::Instance()->getOpt("Database", "select", src) || 
     CsOpt::Instance()->getOpt("Data", "CsStore", src) ) // just for compatibility with previous version
    {
      std::transform(src.begin(), src.end(), src.begin(), ::tolower );
      if(src == "oracle")
	{
	  storageSource = ORACLE_STORAGE;
	  cout << "ORACLE DB reading has been selected by user." << endl;
	}
      else if(src == "file")
	{
	  storageSource = FILE_STORAGE;
	  cout << "Reading DST from the file has been selected by user." << endl;
	}
      else
	{
	  cerr << "CsEvent::CsEvent: incorrect storage source: " << src
	       << "\nThis option will be ignored." << endl;
	}
    }
  //*********************

  // clear stuffs
  _thisRun = 0;
  _previousRun = 0;
  _nOutStreams = 0;
  _outStreamsMask = 0;
  _selTrigMask = 0xffffffff; // DEFAULT: all bits on. 
  _selTrigStrict = false;
  _selTrigMaskZero = false;
  _nmuinspillsofar = 0;

  // clear the event counter and set eventual limits 
  _nEvents = 0;
  _nSelEvents = 0;
  _maxEvents = init->getMaxEvents();
  _skipEvents = init->getSkipEvents();
  if( _maxEvents != 0 ) _maxEvents += _skipEvents;
  _maxConsecutiveSkips = init->getMaxConsecutiveSkips(); // These are skips upon decoding error, as opposed to ``a priori skips'' concerned by the "_skipEvents" supra

  // Monte Carlo Event.
  if( init->IsAMonteCarloJob() ) {
    _MCEvent = true;
    _DTEvent = false; // just to be sure... 
    // Instantiate the CsGeom singleton class
    CsGeom* geom = CsGeom::Instance(init->getDetectorTable());
    // Instantiate the CsGeant3 singleton class
    _geant3MC = CsGeant3::Instance();

  }
  else if( init->IsADataJob() ) {
    _DTEvent = true; 
    _MCEvent = false; // just to be sure... 
    // Instantiate the CsGeom singleton class (To fill CsDetector object)
    if ( init->IsVarPitch() )
      CsGeom* geom = CsGeom::Instance(init->getDetectorTable(),init->getPitchTable());
    else
      CsGeom* geom = CsGeom::Instance(init->getDetectorTable());
    if( init->isFromDB() ) {
      if(storageSource == ORACLE_STORAGE)
	{
#if USE_ORACLE
	  _store = CsOraStore::Instance();
#else
	  _store = CsDummyStore::Instance();
#endif
	}
#if USE_FileStore
      else if(storageSource == FILE_STORAGE)
	{
	  _store = CsFileStore::Instance();
	}
#endif
      else
	{
	  _store = CsDummyStore::Instance();
	}
    }
    else {
//#     if USE_COMPASS_Date
      _store = CsDDDStore::Instance();
//#     else
//      _store = CsDummyStore::Instance();
//#     endif
    }
	
#   if CHECK_TIMEUSAGE
    timeusage( "DATABASE INIT started", 0 );
#   endif
    if( !_store->init() ) {
      CsErrLog::Instance()->mes(elFatal, "CsStore initialization failed." );
    }
#   if CHECK_TIMEUSAGE
    timeusage( "DATABASE INIT ended", 1 );
#   endif

#   if CHECK_TIMEUSAGE
    timeusage( "DATABASE SCAN started", 0 );
#   endif
    if( !_store->scan() ) {
      CsErrLog::Instance()->mes(elFatal, "CsStore scan failed." );
    }
#   if CHECK_TIMEUSAGE
    timeusage( "DATABASE SCAN ended", 1 );
#   endif
  }

  // ***** OK NOW WE ARE READY TO ACCESS INPUT DATA. WHAT SHOULD I DO? *****

  //                                        ***** RECONSTRUCTION SCHEMA?
  if (!opt->getOpt("","reconstruction schema",recoSchema_)) {
    CsErrLog::mes(elError,"No reconstruction schema set => Default one used.");
    recoSchema_ = 1; // Set it to 1 if none selected.
  }

  //                                                     ***** DECODING?
  _decoding = opt->getOpt("","make decoding");
  _MCDecodingExact = false; // MC decoding mode: Trivial Digit-Hit association?
  if (_MCEvent) {
    string mode;
    if (opt->getOpt("","make decoding",mode) && mode=="MCExact")
      _MCDecodingExact = true;
  }

  _clustering = opt->getOpt("","make clustering");  // ***** CLUSTERING?
  if (_clustering && !_decoding)   // Check: decoding -> clustering
    CsErrLog::Instance()->mes(elFatal,
      "Wrong options: `make clustering' but not `make decoding'");

  if (_MCEvent) {
    // Check consistency betweek decoding and clustering options
    if( _clustering ) {
      CsOpt::Instance()->getOpt( "", "make clustering", _MCClusteringMode );
      if( _MCDecodingExact ) {
	if( _MCClusteringMode!="MCExact" && 
	    _MCClusteringMode!="MCSmeared" && 
	    _MCClusteringMode!="MCQuantized" ) {
	  CsErrLog::mes( elFatal, 
	    "Wrong association of MC decoding and clustering modes" );
	}
      }
      else {
	if( _MCClusteringMode=="MCExact" || 
	    _MCClusteringMode=="MCSmeared" || 
	    _MCClusteringMode=="MCQuantized" ) {
	  CsErrLog::mes( elFatal, 
	    "Wrong association of MC decoding and clustering modes" );
	}
      }  
    }
  }

  _tracking = opt->getOpt( "", "make tracking" );     // ***** TRACKING?
  // Check: decoding -> tracking
  if (_tracking && (!_decoding || !_clustering))
    CsErrLog::mes(elFatal,
      "Wrong options: `make tracking' but not decoding or clustering" );
  if (_tracking) _setTrackingPackages();   // Set tracking packages...

  //                                                      ***** LR MODE? 
  _doClustersAssociation = opt->getOpt("LR","cluster association") ||
    opt->getOpt("LR","make associations");
  if (_doClustersAssociation) CsErrLog::mes(elInfo,
      "No raising LR ambiguities by cluster associations in CsEvent" );

  // Calorimeters Reconstruction?
  _calorimeters = 
    CsOpt::Instance()->getOpt( "", "make calorimeters reconstruction" );

  // CEDAR Reconstruction?
  _cedars =
    CsOpt::Instance()->getOpt( "", "make cedar reconstruction" );

  // RW Reconstruction?
  _RWcalorimeters = CsOpt::Instance()->getOpt( "", "make RW reconstruction" );
  if(_RWcalorimeters){ CsRwRecons::Instance();
                cout << "CsRwRecons::  RICH Wall reconstruction is ON" << endl;}
  else cout << "CsRwRecons::  RICH Wall reconstruction is OFF!!" << endl;

  // RW Charge Reconstruction
  _RWchargeEcal1 = CsOpt::Instance()->getOpt( "", "make RW charge reconstruction" );
  if(_RWchargeEcal1){ CsRwChargeRecons::Instance();
                cout << "CsRwChargeRecons::  RICH Wall charge reconstruction is ON" << endl;}
  else cout << "CsRwChargeRecons::  RICH Wall charge reconstruction is OFF!!" << endl;

  // Rich1?
  _rich1 = CsOpt::Instance()->getOpt( "", "make rich1 reconstruction" );

  // Beam?
  _beam = CsOpt::Instance()->getOpt( "", "make beam reconstruction" );
  if(CsOpt::Instance()->getOpt( "BeamRecons", "useTRAFFIC",_beamTR ));
  else _beamTR=0;
  if( _beam ){
     if(CsOpt::Instance()->getOpt( "beam", "method", _beamMethod )){
        if(_beamMethod == 0){
           CsBeamRecons::Instance();
        }
        else if(_beamMethod == 1){
           CsBeamReconstruction::Instance();
        }
        else{
           CsErrLog::Instance()->mes( elFatal, "Wrong value of 'beam method'");
        }
     }else{
        CsErrLog::Instance()->mes( elFatal, "Beam reconstruction method not set" );
     }
  }

  // Vertex?
  _vertex = CsOpt::Instance()->getOpt( "", "make vertex reconstruction" );
  // Check: tracking, beam
  if( _vertex && !_tracking && init->getDataType() != "dst" ) {
    CsErrLog::Instance()->mes( elFatal, 
      "Wrong options: `make vertex reconstruction' but there is not tracking reconstruction" );
  }
  // Set vertex packages...
  if( _vertex ) {
    _setVertexPackages();
  }

  // DST Production ?
  _dst = CsOpt::Instance()->getOpt( "", "make DST" );
  if( _dst ) {
    list<int> dstopt;
    if( CsOpt::Instance()->getOpt( "", "make DST", dstopt ) ) {
      if( dstopt.size() > 0 ) {
	list<int>::iterator ii = dstopt.begin();
	_dstversion = *ii;
	if( dstopt.size() > 1 ) {
	  ii++;
	  _dstfatness = *ii;
	}
	else {
	  _dstfatness = (DST_FATNESS);
	}
      }
      else{
	_dstversion = (DST_VERSION);
	_dstfatness = (DST_FATNESS);
      }
    }
  }

  // These are the two key strings to find TBNames string 
  const string startTBNStr = ":Start of TBNames string:";
  const string endTBNStr   = ":End of TBNames string:";
  // The same for DST version
  const string startDSTVStr = ":Start of DST version string:";
  const string endDSTVStr   = ":End of DST version string:";
  // Fine, what should I need for this run?
  if( _dst ) {  // I've to reconstruct and produce DSTs... 
    // I've to do three things: 
    // - build a string of detector TBnames to be stored on Run header
    // - build a map of detectors
    // - store on Run HEader DST version
    list<CsDetector*> dets = CsGeom::Instance()->getDetectors();
    int i = 0;
    _detvec.clear();
    _detvec.resize( dets.size() ); 
    _detmap.clear();
    string allTBNames = startTBNStr;
    for( list<CsDetector*>::iterator id=dets.begin(); id!=dets.end(); id++ ) {
      allTBNames += (*id)->GetTBName();
      // Now _detvec is used in miniDST production too...
      _detvec[i] = (*id);
      _detmap[*id] = i++;
    }
    allTBNames += endTBNStr;

    string DSTversion = startDSTVStr;
    ostringstream o;
    o << _dstversion;
    DSTversion += o.str();
    DSTversion += endDSTVStr;

    // IF REQUIRED: upload allTBNames to DTS Run Header (only in DST prod.)
    if( _dst && CsOpt::Instance()->getOpt( "", "store TBNames on DST" ) ) {
      if( _store != NULL ) {
        bool status;
        status = _store->uploadTBNamesString( allTBNames );	
        if( !status ) {
	  CsErrLog::mes( elFatal, "Something wrong during TBnames string upload" );
        } 
	// Upload also DST version...
        status = _store->uploadTBNamesString( DSTversion );
        if( !status ) {
	  CsErrLog::mes( elFatal, "Something wrong during DST version string upload" );
        }
	// If just upload strings, everything was done...
	if( CsOpt::Instance()->getOpt( "", "store TBNames on DST only" ) ) {
	  cout << "CsEvent: INFO, Just wrote the TBnames and DST version strings" << endl; 
	  _store->saveAndExit();
	  exit(0);
	}
      }
      else {
	CsErrLog::mes( elFatal, "Cannot access CsStore object" );
      }
    }
  }
  else if( init->getDataType() == "dst" ) { // I've to read DSTs...
    // I need three things:
    // - get the vector of detector TBnames from Run header
    // - build a vector of detectors.
    // - know DST version
    string str;
    string allTBNames;
    // download allTBNames from DTS Run Header...
    if( _store != NULL ) {
      bool status;
      status = _store->downloadTBNamesString( str );
      if( !status ) {
	CsErrLog::mes( elFatal, "Something wrong during TBnames string upload" );
      }
    }

    if( str.empty() ) {
      CsErrLog::mes( elFatal, "TBnames string from DST is empty" );
    }
    const unsigned int TBsize = 8; // Size of TB detector name string
    size_t startPos = str.rfind( startTBNStr );
    size_t endPos   = str.rfind( endTBNStr );
    if( startPos == string::npos || endPos == string::npos ) {
      CsErrLog::mes( elFatal, "TBName string from DBs has no limit key strings" );
    }
    startPos += startTBNStr.size();
    allTBNames = str.substr( startPos, endPos - startPos ); 
    _detvec.resize( allTBNames.size()/TBsize ); 
    for( unsigned int i=0; i<allTBNames.size()/TBsize; i++ ) {
      bool found = false;
      list<CsDetector*> dets = CsGeom::Instance()->getDetectors();
      list<CsDetector*>::iterator id;
      for( id=dets.begin(); id!=dets.end(); id++ ) {
	if( allTBNames.substr( i*TBsize, TBsize ) == (*id)->GetTBName() ) {
	  _detvec[i] = (*id);
	  found = true;
	}
      }
      if( !found ) {
	string str = "No detector found with ";
	str.append( allTBNames.substr( i*TBsize, TBsize ) );
	str.append( " TBname, check used geometry database." );
	CsErrLog::mes( elFatal, str );
      }
    }
#if USE_ORACLE
    startPos = str.rfind( startDSTVStr );
    endPos   = str.rfind( endDSTVStr );
    if( startPos != string::npos ) {
		 if ( endPos == string::npos ) {
			 CsErrLog::mes( elFatal, "DST version string from DBs has no limit key strings" );
		 }
	 }
    startPos += startDSTVStr.size();
    _dstversion = atoi( str.substr( startPos, endPos ).c_str() );
    if( _dstversion == 0 ) {
		 _dstversion = 4;
		 //CsErrLog::mes( elFatal, "DST version is 0" );
    }
#endif
  } 
  else { // Here we're NOT building NOR reading DST 
    // So let's build detector vector&map needed for tracking...
    list<CsDetector*> dets = CsGeom::Instance()->getDetectors();
    int i=0;
    _detvec.clear();
    _detvec.resize( dets.size() ); 
    _detmap.clear();
    for( list<CsDetector*>::iterator id=dets.begin();
	 id!=dets.end(); id++ ) {
      _detvec[i] = (*id);
      _detmap[*id] = i++;
    }
  }
  _actualRICHDataSize=CSTRACK_RICHDATASIZE;
  CsOpt::Instance()->getOpt("Rich", "datasize", _actualRICHDataSize);

  // Any specified trigger mask?
  int triggermask = 0;
  CsOpt::Instance()->setIntMode( "hex" ); // set hex
  if( CsOpt::Instance()->getOpt( "","selection trigger mask",triggermask )) {
    _selTrigMask = triggermask;
    if( CsOpt::Instance()->getOpt( "","selection trigger strict" )) {
      _selTrigStrict = true;
    }
    cout << "-------------------------------------------------------------------" 
	      << endl
	      << "Selection Trigger Mask ";
    if( _selTrigStrict ) cout << "(strict)";
    cout << ": ";
    for( int ii=31; ii>=0; ii-- ) {
      cout << ((triggermask>>ii)&1);
      if( ii==24 || ii==16 || ii==8 ) cout << " ";
    }
    cout << endl;
    cout << "-------------------------------------------------------------------" 
	 << endl;
  }
  if( CsOpt::Instance()->getOpt( "","selection zero trigger mask" )) {
     _selTrigMaskZero = true;
     cout <<" Events with zero trigger mask will be selected " << endl;
  }
  CsOpt::Instance()->setIntMode( "dec" ); // reset to dec

  //  Extra time to be added to detector time windows in order to account for
  // event's trigger jitter being larger than that of the reference trigger
  // for which time windows were specified (In practice, this reference is
  // a mixture of I, L and M).
  //  In case of triggers superposition worse contributor takes precedence.
  // Except for Calorimeter, which intervenes only when alone.
  list<string> strings; // List to house time entries as a function trigger type
  list<string>::iterator istr; double time;
  CsOpt::Instance()->getOpt("Trigger","ExtraTimeWidth",strings);
  for (istr = strings.begin(); istr!=strings.end(); istr++) {
    istringstream(*istr) >> time;
    _extraTimeWidths.push_back(time);
  }
  // Trigger jitter. Solely for MC. Where the jitter for a given event is
  // generated randomly by sampling a Gauss distribution w/ width = MCJitter of
  // the event's trigger bit. If there are several bits in the trigger pattern,
  // the earliest is chosen. However, some bits get a handicap of 5ns, cf.
  // "Trigger MCDelayed" option infra (These are typically Calo and HighQ2).
  strings.clear(); CsOpt::Instance()->getOpt("Trigger","MCJitter",strings);
  for (istr = strings.begin(); istr!=strings.end(); istr++) {
    istringstream(*istr) >> time; _triggerMCJitters.push_back(time);
  }
  _triggerMCDelayed = 0; // Pattern of triggers which are delayed by 5ns
  CsOpt::Instance()->setIntMode("hex"); // Get integer in hexadecimal numbering
  if (CsOpt::Instance()->getOpt("Trigger","MCDelayed",triggermask))
    _triggerMCDelayed = (unsigned int)triggermask;
  else if (init->IsAMonteCarloJob() && _extraTimeWidths.size()>0)
    // Require "MCDelayed" when any "ExtraTimeWidth" has been specified, so as
    // to avoid opening unduly large time windows for events w/ a trigger
    // pattern combining a badly timed trigger (which is typically delayed) and
    // a precisely timed one.
    CsErrLog::mes(elFatal,"\"Trigger\" options: With MC data, \"MCDelayed\" is mandatory if any \"ExtraTimeWidth\" is specified => Check your options file!");
  CsOpt::Instance()->setIntMode("dec"); // Back to decimal numbering
  // Trigger prescaling. Can be useful in the alignment procedure, in order to
  // distribute more uniformly the tracks over phase space, by e.g. prescaling
  // the BeamTrigger.
  strings.clear(); CsOpt::Instance()->getOpt("Trigger","PreScale",strings);
  for (istr = strings.begin(); istr!=strings.end(); istr++) {
    int f; istringstream(*istr)>>f; _triggerPreScales.push_back(f);
  }


  //     *************** BOOK SOME SPILL-RELATED HISTOGRAMS ****************
  CsHistograms::SetCurrentPath("/EVENT");
  _hnmuperspill = new CsHist1D("hnmuperspill","Number of Muons per Spill",
			       CsEvent::_nspillmax,1,1+_nspillmax);
  _hntimeinspill = new CsHist1D("hntimeinspill","event time in spill",100,0,20);
}

const CS::DaqEvent& CsEvent::getDaqEvent   (void) const
{
    return CsInit::Instance()->getDaqEventsManager().GetEvent();
}

const CS::Chip::Digits &CsEvent::getChipDigits(void) const
{
    return CsInit::Instance()->getDaqEventsManager().GetEventDigits();
}


bool CsEvent::_getNextDataEvent() {

# if CHECK_TIMEUSAGE
  static bool firstevent = true;
  if( firstevent ) {
    timeusage( "1ST EVENT DOWNLOAD started", 0 );
  }
# endif

  // Reset the mask for raw event output to streams.
  _outStreamsMask = 0;
  
  // now let see if next event exists...
  
  bool status = false;
  bool dstok  = false;
  bool trigmaskok = false;
  const unsigned int allTrigs = 0xffff; // Although TCS can handle 0x7fffff,
					// it was decided to limit
					// max. #triggers to 16. And take
					// advantage of the fact to make use
					// of the other bits of the trigger
					// TDC module.
  CsInit* init = CsInit::Instance();



  if( ( init->isFromDB() && init->getDataType() == "raw" ) || 
      !(init->isFromDB()) ) {
    do {
      status = _store->next();
      if( status ) {
	unsigned int triggermask = getTriggerMask();
	if (_selTrigStrict)  // If strict: care to strip away cinderella bits
	  trigmaskok = triggermask && 
	    (triggermask&allTrigs)==(triggermask&_selTrigMask);
	else {
	  trigmaskok = (triggermask&_selTrigMask)!=0;
          if( _selTrigMaskZero && triggermask==0 )
            trigmaskok = true;
        }  
      }
    } while( status && !trigmaskok );
  } 
  else{
    if( _store != NULL ) {
      do {
	status = _store->next();
	if( status ) dstok  = _store->dst(_store->getSelectedSlot());
	if( dstok ) {
	  unsigned int triggermask = getTriggerMask();
	  if (_selTrigStrict)  // If strict: care to strip away cinderella bits
	    trigmaskok = (triggermask&allTrigs)==_selTrigMask;
	  else 
	    trigmaskok = (triggermask&_selTrigMask)!=0;
	}
      } while( status && !dstok && !trigmaskok );
    }
  }

# if CHECK_TIMEUSAGE
  if( firstevent ) {
    firstevent = false;
    timeusage( "1ST EVENT DOWNLOAD ended", 1 );
  }
# endif

  return( status );

}


bool CsEvent::getNextEvent()
{

 try_again:
  clear();
  
  if (_chrono==-1) _chrono = _chronometer.start();
  float interT = _chronometer.inter(_chrono), ignoreT = 0;
  bool gotIt = true;
  CsInit *init  = CsInit::Instance();// Instantiate CsInit to access user options.

  if( _maxEvents!=0 && _nEvents>=_maxEvents)   // ***** CUT on #EVENTS PROCESSED
    return false;
  while (_nEvents<_skipEvents && gotIt) {                   // ***** SKIP EVENTS
    _nEvents++;  // Increase the event counter
    if (_MCEvent) gotIt = _geant3MC->getNextEvent();
    else          gotIt = _getNextDataEvent();
    if      (_nEvents==_skipEvents) printf("\rSkipping events %d\n",_nEvents);
    else if (_nEvents%10==0) {
      printf("\rSkipping events %d...",_nEvents); cout<<flush;
    }
  }
  if (gotIt) {
    _nEvents++;                              // ***** INCREASE THE EVENT COUNTER
    int rejected; do {
      if (_MCEvent) gotIt = _geant3MC->getNextEvent( _selTrigMask ); // MC event
      else          gotIt = _getNextDataEvent();                     // RD event
      rejected = 0; if (gotIt) {                             // ***** PRESCALING
	unsigned int triggerMask = getTriggerMask();
	int bit; for (bit = 0; bit<(int)_triggerPreScales.size(); bit++) {
	  if ((1<<bit&triggerMask)) continue;
	  int preScale = _triggerPreScales[bit];
	  if (!preScale) rejected = 1;
	  else if (preScale>1)
	    rejected |= CsRandom::flat()>1./_triggerPreScales[bit];
	}
      }
    }
    while (rejected);
  }

  static bool reloadCalibrations = true;
  if (gotIt) {
    _thisRun = getRunNumber();
    if (_previousRun!=_thisRun) {             // ***** NEW RUN? BOOK CALIBRATION
      _previousRun = _thisRun;
      CsRegistrySing::Instance()->callSorMethods();
      reloadCalibrations = true;
    }
  }
  if (reloadCalibrations) {     // ***** IF BOOKED: GET CALIBS FOR ALL DETECTORS
    reloadCalibrations = false;
    if (!_MCEvent) {
      static string path;
      if( path.empty() )
          if( !CsOpt::Instance()->getOpt( "", "decoding map", path ) )
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
		      "CsEvent::getNextEvent() Cannot can not find \"decoding map\" option.  ");
      CS::Chip::Maps &daq_maps = init->getDaqMaps();
      CS::DaqOption &daq_options = init->getDaqOptions();
      vector<string> dets;
      daq_maps.clear(); daq_options.DaqMap::Clear(); daq_options.Clear();
      CS::Chip::ReadMaps(getRunNumber(),path,daq_maps,daq_options,dets);
      for( map<string,CsDet*>::iterator it=CsDet::GetAllDetectors().begin();
           it!=CsDet::GetAllDetectors().end(); it++ ){
          it->second->SetDaqDataDecodingInfoGeneral(daq_maps);
      }
    }

    _loadCalibTime -= _chronometer.inter(_chrono);

    if (_MCEvent) {               // ***** READ calibrations in MC case... *****
      int time = init->getStartOfRun();
      for( map<string,CsDet*>::iterator it=CsDet::GetAllDetectors().begin();
           it!=CsDet::GetAllDetectors().end(); it++ ){
        it->second->readMCCalibration(time);
	//LS Eff                                                                         
        //to be called here ...
        it->second->readMCEffMaps(time);
      }
    }
    
    if (gotIt && init->useCalibration()) {
      if (init->useCDB()) {
	int time = init->getCDBUseTime();
	if (time==0) time = init->getStartOfRun();
	if (time==0) time = getEventTime().secFrEpoch(); //- 2*60*60;//Should be improve!
	cout << "reading calibration" << endl;
        if (!init->getDB()->ConnectDB()) 
          throw CS::Exception("CsEvent::getNextEvent: can't connect to CDB database");
	for( map<string,CsDet*>::iterator it=CsDet::GetAllDetectors().begin();
	     it!=CsDet::GetAllDetectors().end(); it++ ) {
	  try {
	    it->second->readCalibration(time);
	  } catch( const std::exception &e ) {
	    cerr << "Exception in readCalibration(): "
		 << it->second->GetTBName() << ": " << e.what() << endl;
	    exit(1);
	  }
	}
//read calibrations for tcs phase corrections
	_readTcsCalibration(time);
        //init->getDB()->DisconnectDB(); // Connection is closed infra...
	cout << "end of reading calibration" << endl;
      }
    }
    
    if (!_MCEvent) {                 // ***** READ VARIOUS CsMagInfo *****
      CsField *field = CsGeom::Instance()->getCsField();
      if (field->getNumOfMags()>0) {
        CsMagInfo *magp = field->getMagInfo();
	if (magp[0].fsc!=0)      // Polar, target (solenoid and dipole) currents
	  magp[0].readPolarization(field,_store->getRun());
	if (field->getNumOfMags()>2) // SM2 NMR
	  magp[2].scaleSM2withNMR(field,_store->getRun());
      }
    }
    
    if (init->useCDB() && init->getDB()->isConnected())
      init->getDB()->DisconnectDB();

    if( gotIt ) {
      _loadCalibTime += _chronometer.inter(_chrono);
      ignoreT         = _loadCalibTime;
    }

  }

  static bool first_time = true; static int iStream = -1; if (first_time) {
    first_time = false;
    // ********** FIRST EVENT (processed) => INITIALISATIONS **********

    //                               ***** BOOK DETECTORS' HISTOS *****
    cout << "CsEvent::getNextEvent - INFO: booking histograms" << endl;
    list<CsDetector*> dets = CsGeom::Instance()->getDetectors();
    for( list<CsDetector*>::iterator Id = dets.begin(); Id != dets.end(); Id++ )
      (*Id)->BookHistograms();
    cout << "CsEvent::getNextEvent - INFO: end of histogram booking" << endl;

    list<string> streams;        // ********** OUTPUT STREAMS **********
    // Syntax: CsEvent WriteBack <file_name> Random
    // Note: The same "OutputStream" can be used to output events interesting
    // from a physics point of view, cf. "CsKalmanFitting::doFitting"
    if (CsOpt::Instance()->getOpt("CsEvent","WriteBack",streams)) {
      list<string>::iterator is = streams.begin();
      //      ***** FIRST FIELD = FILE NAME *****
      if ((iStream = CsEvent::Instance()->openRawEventOutputStream(*is))<0)
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
		      "Cannot open WriteBack stream \"%s\"",(*is).c_str());
      is++;// ***** FOLLOWING FIELDS: SPECIF. KINDS of EVENT TO WRITE-BACK *****
      int doStream; for (doStream = 0; is!=streams.end(); is++) {
	if (*is=="Randoms") doStream |= 0x1;
	else
	  CsErrLog::msg(elFatal,__FILE__,__LINE__,
			"Unknown WriteBack selection \"%s\"",(*is).c_str());
      }
      if (doStream==0)
	CsErrLog::mes(elFatal,"WriteBack: No selection specified!");
    }
  }

  if( gotIt ) {
    _eventLoadTime += (_chronometer.inter(_chrono) - interT -ignoreT);
    _neventLoadTime++;
  }

  if (gotIt) { //    ***** WRITE BACK RANDOMS (if requested) *****
    // Random triggers pattern HARD-CODED, taken from "2006.xml/DAQ.xml":
    //    SlowRandomTrigger   = "10"
    //    RandomTrigger       = "11"
    // !!WARNING!! This trigger bit assignment is not always right: case of
    // 2009 Primakoff.
    unsigned int triggerMask = getTriggerMask();
    const unsigned int randomTrigs = 0xc00;
    const unsigned int calibTrigs = 0x1ffef000;
    if (iStream>=0 && ((triggerMask&randomTrigs) || (triggerMask&calibTrigs)))
      outputRawEventToStream(iStream);
  }

  // Call the decoding and reconstruction schema handler...
 
  // MC? Raw from DB? Raw from file? Then make event reconstruction
  // DST? Then just download reconstructed event
  if( ( ( ( init->isFromDB() && init->getDataType() == "raw" ) || 
	  !(init->isFromDB()) ) && _DTEvent ) || _MCEvent ) { 

    // The data have to be reconstructed...

    // Decode no more in reco schema. BG 2003/06/26
    if( gotIt ) {
      _digitizTime -= _chronometer.inter(_chrono);
      gotIt = _decode();
      _digitizTime += _chronometer.inter(_chrono);
      _ndigitizTime++;

      static unsigned int consecutiveSkips = 0;
      if (!gotIt) {// ***** SKIP EVENT IF DECODING TROUBLES. BG 2003/06/26 *****
	CsErrLog::msg(elError, __FILE__, __LINE__,
		      "+----------------------------------------------------------+ \n"
		      "Event skipped due to decoding troubles. \n"
		      "Run: %d, Event: %d \n"
		      "+----------------------------------------------------------+",
		      getRunNumber(), getEventNumberInRun() );
	_nSkippedEvents++; ;
	if (_maxConsecutiveSkips && ++consecutiveSkips>_maxConsecutiveSkips) {
	  CsErrLog::msg(elError, __FILE__, __LINE__,
			"Consecutive events skipped = %d, > tolerated spate of errors(=%d) => give up",
			consecutiveSkips,_maxConsecutiveSkips);
	  return false ;
	}
	goto try_again;
      }
      else consecutiveSkips = 0;
    }
      
    if( gotIt ) {
#if USE_MySQL
      if (init->updateSolenoidField()) { // This implies !MC, cf. "CsInit.h"
	// Rescale solenoid field if required, e.g. in case of rotation data.
	CsField *field = CsGeom::Instance()->getCsField();
	if (field->getNumOfMags()>0) {
	  time_t evTime = getEventTime().secFrEpoch();
	  MySQLDB *_mysqldb = (MySQLDB*) init->getDB();
	  // ***** GET SOLENOID CURRENT AS A F(TIME) *****
	  //  Note (from Damien): "getTgtCurrentsTime" keeps in memory the
	  // results of the db accesses (times and values of the currents for 2
	  // measurements around given time), so that, if next call concerns the
	  // same time interval, there will be no access to db. In 2004 there
	  // is a measurement every minute. If we say there are 10
	  // chunks/EB/run, a chunk covers ~ 6 minutes, so let say there will be
	  // 7 or 8 accesses for each job. In my mind it is not a lot.
	  //  No need to _mysqldb->connectDB(): it's handled internally.
	  pair<double,double> curs = _mysqldb->getTgtCurrentsTime(evTime);
	  field->getMagInfo()[0].scaleSolenoidMap(field,(float)curs.first);
	}
	else CsErrLog::mes(elFatal,
			   "Solenoid update required while no target magnet!");
      }
#endif
      gotIt = _reconstructionSchemaHandler();
    }

    // DST production?
    if( gotIt && _dst ) {

      _dstUpTime -= _chronometer.inter(_chrono);

      _upload();

      _dstUpTime += _chronometer.inter(_chrono);
      _ndstUpTime++;
    }
  }
  else if( init->getDataType() == "dst" && gotIt ) {  

    _dstDownTime -= _chronometer.inter(_chrono);

    // The data have to be downloaded from DST
    while( ! _download() && gotIt ) {
      cerr << "Error found in DST data, the event will be skipped" << endl;
      // read next event...
      gotIt = _getNextDataEvent();

      // Verify run switch...
      if( gotIt ) {
	_thisRun = getRunNumber();
	if( _previousRun != _thisRun ) {
	  _previousRun = _thisRun;
	  CsRegistrySing::Instance()->callSorMethods();
	}
      }
    }

    _dstDownTime += _chronometer.inter(_chrono);
    _ndstDownTime--;
  }


  if( gotIt ) {
    _finishTime -= _chronometer.inter(_chrono);
  }

  if( gotIt ) {
    _finishTime += _chronometer.inter(_chrono);
    _nfinishTime++;
    
    _totalTime += (_chronometer.inter(_chrono) - interT -ignoreT );
    _ntotalTime++;
  }

  if( gotIt) { _nSelEvents++; }

  return( gotIt );
}

void CsEvent::clear(void)
{
  if( CsInit::Instance()->getDataType() == "dst" ) {
    _recoEvent.clearBeamTracks();
  }
  _recoEvent.clear();


  list<CsDetector*> dets = CsGeom::Instance()->getDetectors();
  for( list<CsDetector*>::iterator Id=dets.begin(); Id!=dets.end(); Id++ ) {
    (*Id)->clearDigitsList();
  }

  CsRICH1Detector* rich = CsGeom::Instance()->getRich1Detector();
  if( rich != NULL ) {
    rich->clearDigitsList();
  }

  for( map<string,CsDet*>::iterator it=CsDet::GetAllDetectors().begin();
       it!=CsDet::GetAllDetectors().end(); it++ )
    it->second->Clear();
}

const list <CsMCTrack*> &CsEvent::getMCTracks(void) const
{
  if( _MCEvent )
    return _geant3MC->getMCTracks();
  else
  {
    cerr << "CsEvent::getMCTracks(): this is not MC event!\n";
    exit(1);
  }
}

const list <CsMCVertex*> &CsEvent::getMCVertices(void) const
{
  if( _MCEvent )
    return _geant3MC->getMCVertices();
  else
  {
    cerr << "CsEvent::getMCVertices():  this is not MC event!\n";
    exit(1);
  }
}

const list <CsMCHit*> &CsEvent::getMCHits(void) const
{
  if( _MCEvent )
    return( _geant3MC->getMCHits() );
  else
  {
    cerr << "CsEvent::getMCHits():  this is not MC event!\n";
    exit(1);
  }
}

uint32 CsEvent::getEventSize(void) const
{
  if( _DTEvent )
    return( *(_recoEvent.getEventHeader()) );
  else
  {
    cerr << "CsEvent::getEventSize():  Not available on MCs!\n";
    exit(1);
  }
}

uint32 CsEvent::getRunNumber(void) const
{
  if( _DTEvent ) {
    if( CsInit::Instance()->getDataType() == "dst" && _dstversion > 3 ) {
      return( *(_recoEvent.getEventHeader()+1) );
    }
    else {
      return _store->getRun();
    }
  }
  else if( _MCEvent) {
    return _geant3MC->getRunNb();
  }
  else
  {
    cerr << "CsEvent::getRunNumber():  Unknown event type!\n";
    exit(1);
  }
}

uint32 CsEvent::getEventNumberInRun(void) const
{
  if( _DTEvent ) {
    if( CsInit::Instance()->getDataType() == "dst" && _dstversion > 3 ) {
      return( *(_recoEvent.getEventHeader()+2) );
    }
    else {
      return _store->getEventInRun();
    }
  }
  else if( _MCEvent) {
    return _geant3MC->getEventNb();
  }
  else
  {
    cerr << "CsEvent::getEventNumberInRun():  Unknown event type!\n";
    exit(1);
  }
}

uint32 CsEvent::getBurstNumber(void) const
{
  if( _DTEvent ) {
    if( CsInit::Instance()->getDataType() == "dst" && _dstversion > 3 ) {
      return( *(_recoEvent.getEventHeader()+3) );
    }
    else {
      return _store->getBurst();
    }
  }
  else {
    return (1+int(this->getEventNumberInRun()/1000));
  }
}

uint32 CsEvent::getEventNumberInBurst(void) const
{
  if( _DTEvent ) {
    if( CsInit::Instance()->getDataType() == "dst" && _dstversion > 3 ) {
      return( *(_recoEvent.getEventHeader()+4) );
    }
    else {
      return( _store->getEventInBurst() );
    }
  }
  else {
    return (1+this->getEventNumberInRun()%1000);
  }
}

uint32 CsEvent::getTriggerNumber(void) const
{
  if( _DTEvent )
    return( *(_recoEvent.getEventHeader()+5) );
  else
    return 0;
}

pair<uint32,uint32> CsEvent::getTime(void) const
{
  if( _DTEvent ) {
    pair<uint32,uint32> p(*(_recoEvent.getEventHeader()+6),
			  *(_recoEvent.getEventHeader()+7));
    return( p ) ;
  }
  else {
    pair<uint32,uint32> p(0,0);
    return( p );
  }
}

uint32 CsEvent::getErrorCode(void) const
{
  if( _DTEvent )
    return( *(_recoEvent.getEventHeader()+8) );
  else
    return 0;
}

uint32 CsEvent::getTriggerMask(void) const
{
    if( _DTEvent )
		 if ( CsInit::Instance()->getDaqEventsManager().IsEventAvailable() ) {
			 return getDaqEvent().GetTrigger(); // like rather daq info if any...
		 } else {
			 return _store->getTriggerMask();   // ...to those stored in DB
		 }
    else if( _MCEvent )
        return _geant3MC->getTriggerMask();
    else
        throw "CsEvent::getTriggerMask(): Internal problem: unknown event type!";
}

CsTime CsEvent::getEventTime() const
{
  if( _DTEvent )
    return _store->getTime();
  else
  if( _MCEvent )
    return _geant3MC->getEventTime();
  else
  {
    cerr << "CsEvent::getEventTime(): unknown event type."<<endl;
    exit(1);
  }
}

bool CsEvent::_reconstructionSchemaHandler()
{
  bool status = false;

  // Simple at present. Something more nice in future...
  switch( recoSchema_ ) {
  case 0:
    status = _reconstructionSchemaPP();
    break;
  case 1:
    status = _reconstructionSchema001();
    break;
  case 2:
    status = _reconstructionSchema002();
    break;
  default:
    CsErrLog::Instance()->mes( elFatal, 
      "No reconstruction schema available. Stop processing" );
    break;
  }

  return status;

}


void CsEvent::_decodeMC(void) {
  
  //    ************************************************************
  //    ******************** MONTE-CARLO EVENT *********************
  //    ************************************************************
  
  if (CsInit::Instance()->resetRandomSeed()) {
    // ***** optionally: RESET THE RANDOM NUMBER GENERATOR *****
    // Can be useful for debugging purposes, so as to get simillar events
    // despite diverging random generation histories.
    CsRandom *random = CsRandom::Instance();
    long oldSeed = random->getSeed();
    // In order no to have the same sequence of random numbers for all
    // successive events, one should reset a (possibly only slightly)
    // different seed, by, e.g., taking into account the event number. This
    // is done in "../evmc/CsGeant3.cc". Here, we simply re-instate the
    // CsGeant3's seed.
    random->setSeed(oldSeed);
    // Reset the, pair-wise, generation of gaussian numbers.
    random->flagGauss();
    CsErrLog::msg(elInfo,__FILE__,__LINE__,"Resetting CsGeant3 random seed: 0x%lx\n",oldSeed);
  }

  //                   ***** TRIGGER JITTER *****
  // - Determine the trigger jitter of the current event. All triggers w/in
  // the trigger pattern are independently randomly timed and made to compete
  // for the definition of the overall trigger time (as I(Y.B.) understand is
  // done in reality).
  // - Taking into account the 5ns delaying the worst timed triggers. Which
  // is achieved by first delaying and letting delayed triggers to compete w/
  // the rest. Then, if all component triggers are of a delayed type,
  // subtracting out the 5ns delay. If not, and even if the winner is of a
  // delayed type => no subtraction. This is what must have been enforced by
  // the decoding library in the latest productions of [2002,2006] data. But
  // the cases of 2007, and 2008 w/ its early accidental triggering, may not
  // fit this description...
  unsigned int triggerMask = getTriggerMask();
  int bit, status;
  for (bit = status = 0, _triggerMCOffset = 0;
       bit<(int)_triggerMCJitters.size(); bit++)
    if (1<<bit&triggerMask) {
      double time = _triggerMCJitters[bit] * CsRandom::gauss();
      bool delayed = 1<<bit&_triggerMCDelayed; if (delayed) time += 5;
      if (time<_triggerMCOffset || status==0) {
	_triggerMCOffset = time; status = delayed ? 2 : 1;
      }
    }
  if (status==2 && (triggerMask&_triggerMCDelayed)==triggerMask)
    // If all trigger bits are of "delayed" type, and the 5ns delay has been
    // effective...
    _triggerMCOffset -= 5;  // ...subtract it out.
  
  _extraTimeWidth = 0; // ***** EXTRA TIME WIDTH for CURRENT EVENT *****
  // (This is meant to make for the large trigger jitter associated w/ some
  // of the triggers. In MC, contrary to the RD case, the expectation for the
  // jitter is defined once and for all, derived from the trigger pattern.
  // Therefore the determination of "_extraTimeWidth" is simplified: the
  // largest width is retained, except if it corresponds to a delayed trigger
  // and there are besides non-delayed trigger in the trigger pattern.
  unsigned int itrig; int extraTimeSet;
  for (itrig = 0, extraTimeSet = 0; itrig<_extraTimeWidths.size(); itrig++) {
    if ((triggerMask&1<<itrig) &&
	// Consider only non-delayed trigger in a first step
	!(_triggerMCDelayed&1<<itrig) &&
	_extraTimeWidth<_extraTimeWidths[itrig]) {
      // ...Retain worst such trigger...
      _extraTimeWidth = _extraTimeWidths[itrig]; extraTimeSet = 1;
    }
  }
  if (!extraTimeSet) {  // ...If extra time not yet set non-delayed trigger...
    for (itrig = 0; itrig<_extraTimeWidths.size(); itrig++) {
      if ((triggerMask&1<<itrig) &&
	  (_triggerMCDelayed&1<<itrig) && // ...Consider the delayed triggers
	  _extraTimeWidth<_extraTimeWidths[itrig])
	_extraTimeWidth = _extraTimeWidths[itrig];
    }
  }

  // store TCS phase and Time in Spill into Calorimeters classes
  setTCSPhaseTime ( 40. );        // ...Dummy TCS Phase for MC
  {
    const double tcsphase = getTCSPhaseTime();
    const double tis      = -1;   // disable time in spill energy corrections for MC
    for( vector<CsCalorimeter*>::iterator it = CsGeom::Instance()->getCalorimeters().begin(); 
         it!=CsGeom::Instance()->getCalorimeters().end(); it++ ) {
      (*it)->SetTCSPhase   ( tcsphase );
      (*it)->SetTimeInSpill( tis );
      (*it)->SetEventIDInfo();
    }
  }

  if (_MCDecodingExact) {
    // At the moment do all digits, regardless what's set...
    list<CsMCHit*> hits = getMCHits();
    for (list<CsMCHit*>::iterator Ih = hits.begin(); Ih!=hits.end(); Ih++) {
      CsMCParticle* particle = ((*Ih)->getMCTrack())->getParticle();
      CsDet* det = (*Ih)->getDet();
      // only charged particles in active detector area
      CsDetector* trkdet = dynamic_cast<CsDetector*>(det);
      
      if(  particle->getCharge()!=0 && trkdet!=0 ) {
	// WARNING: the digit IS dummy!!!!
	CsMCDigit* mcdigit = new CsMCDigit( *det, 0 );
	CsMCHit* hit = (*Ih);
	mcdigit->addHit( *hit );
	// add this digit to its detector list
	det->addDigit( *mcdigit );
	// add this digit to Reconstructed Event list
	addDigit( *mcdigit );
      }
    }
    // set all detector decoding flags as "done"
    list<CsDetector*> dets = CsGeom::Instance()->getDetectors();
    for( list<CsDetector*>::iterator Id=dets.begin(); Id!=dets.end(); Id++ )
      (*Id)->setDecodingDone();
    CsRICH1Detector* rich = CsGeom::Instance()->getRich1Detector();
    if( rich != NULL ) {
      rich ->setDecodingDone();
    }
    
    vector<CsCalorimeter*> calos = CsGeom::Instance()->getCalorimeters();
    for( vector<CsCalorimeter*>::iterator it = calos.begin(); it!=calos.end(); it++ ) {
      (*it)->makeMCDecoding();
    } 
  }
  else {
    for (map<string,CsDet*>::iterator it = CsDet::GetAllDetectors().begin();
	 it!=CsDet::GetAllDetectors().end(); it++)
      it->second->makeMCDecoding();
  }
}


bool CsEvent::_decodeRD(void) {

  //   ************************************************************
  //   ********************* REAL DATA EVENT **********************
  //   ************************************************************
  try {
    
#   ifdef TIME_CHECK
    int stardet=stopwatch.start();
#   endif    

//       static string path;
//       if( path.empty() )
//           if( !CsOpt::Instance()->getOpt( "", "decoding map", path ) )
//               throw CS::Exception("CsEvent::_decode(): can not find \"decoding map\" option.");
//       static unsigned run_old = unsigned(-1);
    CS::Chip::Maps &daq_maps = CsInit::Instance()->getDaqMaps();
//       CS::DaqOption &daq_options = CsInit::Instance()->getDaqOptions();
//       if (run_old!=getRunNumber()) {  // ********** NEW RUN NUMBER **********
// 	vector<string> dets;
// 	daq_maps.clear(); daq_options.DaqMap::Clear(); daq_options.Clear();
// 	CS::Chip::ReadMaps(getRunNumber(),path,daq_maps,daq_options,dets);
// 	run_old = getRunNumber();
//         for( map<string,CsDet*>::iterator it=CsDet::GetAllDetectors().begin();
//            it!=CsDet::GetAllDetectors().end(); it++ ){
//           it->second->SetDaqDataDecodingInfoGeneral(daq_maps);
//         }
//       }

    //  *************** TIME GATE for CURRENT EVENT ***************
    
    _extraTimeWidth = 0; // ***** EXTRA TIME WIDTH for CURRENT EVENT *****
    // (This is meant to make for the large trigger jitter associated w/ some
    // of the triggers. The expectation for the jitter will be refined later
    // on, after all data including the trigger pattern TDC have been decoded:
    // the present "getExtraTimeWidth" is an upper bound and is to be used to
    // cut away hits that can't in any case be retained, in order to speed
    // up the processing.)
    unsigned int triggerMask = getTriggerMask();
    for (unsigned int itrig = 0; itrig<_extraTimeWidths.size(); itrig++) {
      if ((triggerMask&1<<itrig) &&
	  _extraTimeWidth<_extraTimeWidths[itrig]) // Let worst trigger set...
	_extraTimeWidth = _extraTimeWidths[itrig]; // ...extra time width
    }
    
#   ifdef TIME_CHECK
    if(TiDecod){ dec1+=stopwatch.stop(stardet); stardet=stopwatch.start(); }
    if(TiDecod){ dec2+=stopwatch.stop(stardet); }          
#   endif

    if (daq_maps.size()==0)
      throw CS::Exception("CsEvent::_decode():WW: daq maps was not found for run %d",getRunNumber());

    // Data Events...
#   ifdef TIME_CHECK 
    if(TiDecod){ stardet=stopwatch.start(); }   
#   endif  

    // When reading from db, DaqEvent buffer must be 
    if ( CsInit::Instance()->isFromDB() ) {
      CsInit::Instance()->getDaqEventsManager().SetDaqEvent(_store->rawBuffer());
      //duic: Patch to restore trigger mask during production
#if USE_ORACLE
      if ( getRunNumber() > 45063 && getRunNumber() < 54994 ) { // only for 2006 data
	PatchEventHeaderDB::Instance()->DBPatchTriggerMaskForEvent(CsEvent::Instance(), CsEvent::Instance()->GetStore(), const_cast<CS::DaqEvent*>(&(getDaqEvent())));
      }
#endif
    }

    CsInit::Instance()->getDaqEventsManager().DecodeEvent();
    
#   ifdef TIME_CHECK
    if(TiDecod){ stardet=stopwatch.start(); }  
    if(TiDecod){ dec5+=stopwatch.stop(stardet); stardet=stopwatch.start(); }
#   endif
    
  }
  catch(std::exception &e) {
    cerr << "CsEvent::_decode(): exception:\n" << e.what() << "\n";
    return( false ); 
  } 
  catch(std::string &e) {
    cerr << "CsEvent::_decode(): exception:\n" << e << "\n";
    return( false ); 
  } 
  catch(char *e) {
    cerr << "CsEvent::_decode(): exception:\n" << e << "\n";
    return( false ); 
  } 
  catch(...) {
    cerr << "CsEvent::_decode(): unknown exception\n";
    return( false );
  }
  
  // Make TCS Phase available to the decoding step of detector classes
	if( _tcscorr.find(getBurstNumber()) == _tcscorr.end() )
		setTCSPhaseTime( CsEvent::Instance()->getDaqEvent().GetTT().GetPhaseTCS() );
	else
		setTCSPhaseTime( CsEvent::Instance()->getDaqEvent().GetTT().GetPhaseTCS() - _tcscorr[getBurstNumber()]);

  // get time in spill from libDDD
  _recoEvent.setTimeInSpill( CsEvent::Instance()->getDaqEvent().GetTT().GetTimeInSpill() );

  // store TCS phase and Time in Spill into Calorimeters classes
  {
    vector<CsCalorimeter*> &calos = CsGeom::Instance()->getCalorimeters();
    const double tcsphase = getTCSPhaseTime();
    const double tis      = getTimeInSpill();
    for (vector<CsCalorimeter*>::iterator it = calos.begin(); it!=calos.end(); it++) {
      (*it)->SetTCSPhase   ( tcsphase );
      (*it)->SetTimeInSpill( tis );
      (*it)->SetEventIDInfo();
    }
  }
  
  //     *************** DECODE ALL DETECTORS ***************
  
  for( map<string,CsDet*>::iterator det=CsDet::GetAllDetectors().begin(); det!=CsDet::GetAllDetectors().end(); det++ ) {
    // in case decoding for one detector fails, skip to the next detector
    try {
      det->second->DecodeChipDigits( getChipDigits() );
    }
    catch(std::exception &e) {
      CsErrLog::msg(elError, __FILE__, __LINE__,
		    "CsEvent::_decode(): exception:\n%s\nDecoding error for %s, detector may be empty during reconstruction.",
		    e.what(), det->second->GetTBName().c_str());
    } 
    catch(std::string &e) {
      CsErrLog::msg(elError, __FILE__, __LINE__,
		    "CsEvent::_decode(): exception:\n%s\nDecoding error for %s, detector may be empty during reconstruction.",
		    e.c_str(), det->second->GetTBName().c_str());
    } 
    catch(char *e) {
      CsErrLog::msg(elError, __FILE__, __LINE__,
		    "CsEvent::_decode(): exception:\n%s\nDecoding error for %s, detector may be empty during reconstruction.",
		    e, det->second->GetTBName().c_str());
    }
    catch(...) {
      CsErrLog::msg(elError, __FILE__, __LINE__,
		    "CsEvent::_decode(): unknown exception caught:\nDecoding error for %s, detector may be empty during reconstruction.",
		    det->second->GetTBName().c_str());
    }
  }
  
# ifdef TIME_CHECK
  if(TiDecod){ dec6+=stopwatch.stop(stardet); }
  if(TiDecod){ stardet=stopwatch.start(); dec7+=stopwatch.stop(stardet);}
# endif

  //             *************** SCALERS ***************
  typedef multimap<CS::DetID,CS::Chip::Digit*>::const_iterator mIt;
  pair<mIt,mIt> mRange;
  if (getRunNumber()<24885) { // If in 2002..
    // (Note: The run# serving here as an upper bound corresponds to what's in
    // in 2002's SCALER maps, which turn out to be consistent w/ 2003's.
    string scaler = "SC99P2__";
    mRange = getChipDigits().equal_range(scaler);
    for( mIt dIt = mRange.first; dIt != mRange.second; dIt++ ) {
      CS::Scaler::Digit* sd = dynamic_cast<CS::Scaler::Digit*>(dIt->second);
      if( sd != NULL ) {
	switch ( sd->GetChannel() ) {
	  // The counts below are used to determine the flux in 2002 data.
	case 16: _recoEvent.setSC99P2_16( sd->GetValue() ); break;
	case 17: _recoEvent.setSC99P2_17( sd->GetValue() ); break;
	case 18: _recoEvent.setSC99P2_18( sd->GetValue() ); break;
	case 19: _recoEvent.setSC99P2_19( sd->GetValue() ); break;
	case 20: _recoEvent.setSC99P2_20( sd->GetValue() ); break;
	case 21: _recoEvent.setSC99P2_21( sd->GetValue() ); break;	  
	}
      }
      else {
	throw CS::Exception("CsEvent::_decode(): scaler digit not casted as scaler object");
      }
    }
    scaler = "SC01P1__";  // This does not exist any longer beyond 2002!
    mRange = getChipDigits().equal_range(scaler);
    for( mIt dIt = mRange.first; dIt != mRange.second; dIt++ ) {
      CS::Scaler::Digit* sd = dynamic_cast<CS::Scaler::Digit*>(dIt->second);
      if( sd != NULL ) {
	if( sd->GetChannel() == 1 ) _recoEvent.setSC01P1_01( sd->GetValue() );
	if( sd->GetChannel() == 9 ) _recoEvent.setSC01P1_09( sd->GetValue() );
      }
      else {
	throw CS::Exception("CsEvent::_decode(): scaler digit not casted as scaler object");
      }
    }
  }
  //else if (getRunNumber()<45000) // In 2003, 2004
  // No flux from scaler count in 2003/4 yet...
  else if (getRunNumber()>=45000) { // In 2006, 2007...

    CsOpt *opt = CsOpt::Instance();
    list<string> FluxScalers;
    if ( !opt->getOpt("", "FluxScalers", FluxScalers) )
      CsErrLog::mes(elFatal, "FluxScalers option missing!  It must contain all "
		    "tbnames of scalers used for flux measurement.");

    for ( list<string>::const_iterator it = FluxScalers.begin(); it != FluxScalers.end(); it++) {
      static unsigned good_count = 0, bad_count  = 0;
      unsigned int mcount = 0;
      const char proj = it->at(6);
      assert ( proj == 'X' || proj == 'Y' );
      mRange = getChipDigits().equal_range(*it);
      for (mIt dIt = mRange.first; dIt!=mRange.second; dIt++) {
	CS::Scaler::Digit *sd = dynamic_cast<CS::Scaler::Digit*>(dIt->second);
	if (sd==NULL) {
	  throw CS::Exception("CsEvent::_decode(): scaler digit not casted as scaler object");
	}
	if (sd->GetChannel()<32)
	  mcount += (unsigned int)sd->GetValue();
      }

      if (mcount == 0)
	bad_count++;
      else
	good_count++;

      if ( bad_count > 500 && good_count/7 < bad_count ) {
	stringstream s;
	s << "Too many bad flux scaler readings: " << bad_count
	  << " bad readings and only " << good_count << " good readings. "
	  << "Probably your FluxScaler option is set incorrectly. "
	  << "Present bad reading is from " << *it << ".";
	static int prvEvt = -1, warning_count = -1; Severity severity;
	int evt = getEventNumberInRun(); if (evt!=prvEvt) {
	  prvEvt = evt; warning_count++;
	}
#define CsE_SKIP_FSC_ERROR 50
	if (!warning_count)
	  severity = elError;	
	else if (!(warning_count%CsE_SKIP_FSC_ERROR)) {
	  severity = elError;
	  if (it==FluxScalers.begin())
	    s << "\n(Note: " << CsE_SKIP_FSC_ERROR <<
	      " bunches of such errors have been downgraded to warning level.)";
	}
	else severity = elWarning;
	CsErrLog::mes(severity,s.str());
      }

      if (proj=='X') _recoEvent.addFluxX(mcount);
      else           _recoEvent.addFluxY(mcount);
    }
  }

  for( int i=0; i<24; i++ ) {
    mRange = getChipDigits().equal_range( HodoDataStruct[i].tbname );
    for( mIt dIt = mRange.first; dIt != mRange.second; dIt++ ) {
      CS::ChipF1::Digit* fd = dynamic_cast<CS::ChipF1::Digit*>(dIt->second);
      if( fd != NULL ) {
	int channel = fd->GetChannel();
	//int time   = fd->GetTime() - CS::ChipF1::GetTT().GetTimeNorm();
	double time=fd->GetTimeDecoded()/fd->GetTimeUnit();
	
	if( HodoDataStruct[i].lower<time&&HodoDataStruct[i].upper>time ) {
	  unsigned int datum = int( - time * 100.0 + 0.5 )&0x1fffff;
	  datum = datum | ((channel&0x3f)<<21);
	  datum = datum | ((i&0x1f)<<27);
	  _recoEvent.addHodoDatum( datum );
	}
      }
      else {
	throw CS::Exception( "CsEvent::_decode(): Hodo digit not casted as F1" );
      }
    }
  }

  return true;
}


bool CsEvent::_decode(void)
{
  if (!_decoding) return(true);
#ifdef TIME_CHECK
  CsStopwatch stopwatch;
  int star = stopwatch.start();
#endif    

  // ******* DECODING ******** //
  if (_MCEvent) {
    _decodeMC();
  } else {
    if ( ! _decodeRD() )
      return false;
  }


  //        *************** STORE ALL DIGITS IN _recoEvent ***************

  list<CsDetector*> dets = CsGeom::Instance()->getDetectors();
  for( list<CsDetector*>::iterator idt=dets.begin(); idt!=dets.end(); idt++ ) {
    list<CsDigit*> digs = (*idt)->getMyDigits();
    for( list<CsDigit*>::iterator idg=digs.begin(); idg!=digs.end(); idg++ ) { 
      _recoEvent.addDigit( *(*idg) );
    }
  }
  CsRICH1Detector*  rich = CsGeom::Instance()->getRich1Detector();
  if( rich != NULL ) {
    list<CsDigit*> digs = rich->getMyDigits();
    for( list<CsDigit*>::iterator idg=digs.begin(); idg!=digs.end(); idg++ ) { 
      _recoEvent.addDigit( *(*idg) );
    }
  }

  vector<CsCalorimeter*> calos = CsGeom::Instance()->getCalorimeters();
  for( vector<CsCalorimeter*>::iterator it = calos.begin(); it!=calos.end(); it++ ) {
      list<CsDigit*> digs = (*it)->getMyDigits();
      for( list<CsDigit*>::iterator idg=digs.begin(); idg!=digs.end(); idg++ ) { 
	_recoEvent.addDigit( *(*idg) );
      }  
  } 

  list<CsMiscDetector*> others = CsGeom::Instance()->getMiscDetectors();
  for( list<CsMiscDetector*>::iterator idt=others.begin(); idt!=others.end(); idt++ ) {
    list<CsDigit*> digs = (*idt)->getMyDigits();
    for( list<CsDigit*>::iterator idg=digs.begin(); idg!=digs.end(); idg++ ) { 
      _recoEvent.addDigit( *(*idg) );
    }
  }

  _recoEvent.sortDigits(); // *************** SORT DIGITS ***************


  if( _DTEvent ) {
    // TCS setting is at the begining of decoding  
    try {                                // Get Trigger Time
      _recoEvent.setTriggerTime( getDaqEvent().GetTT().GetTime(0) );
    }
    catch (std::exception &e) {
      cerr << "CsEvent::_decode(): exception:" << endl << e.what() << endl;
      return( false );
    }
    catch(...) {
      cerr << "CsEvent::_decode(): unknown exception" << endl;
      return false;
    }
  } 
  else               // ***** MONTECARLO...
    _recoEvent.setTriggerTime( 0 ); // ...Is this useful?

  if (_DTEvent) {// *************** STORE Date EVENT HEADER DATA ***************
#warning TODO: Implement proper fix to missing fields in DATE 5 header.
    // Two lines commented out w/ the advent of Date,v5:
    // _daqEvent->GetTriggerNumber() and _daqEvent->GetErrorCode() throw exceptions when processing data recorded with Date,v5
    // Apparently those two information are not avalaible in new header format.
    // We should check if those informations are needed, and implement a proper fix.
    _recoEvent.setEventHeader( getDaqEvent().GetLength(),                  0 );
    _recoEvent.setEventHeader( getDaqEvent().GetRunNumber(),               1 );
    _recoEvent.setEventHeader( getDaqEvent().GetEventNumberInRun(),        2 );
    _recoEvent.setEventHeader( getDaqEvent().GetBurstNumber(),             3 );
    _recoEvent.setEventHeader( getDaqEvent().GetEventNumberInBurst(),      4 );
	 if ( getDaqEvent().GetHeader().GetVersion()<0xffff ) {
		 _recoEvent.setEventHeader( getDaqEvent().GetTriggerNumber(),           5 );
	 } else {
		 _recoEvent.setEventHeader(0, 5); // UNAVAILABLE for Date,v5
	 }
    _recoEvent.setEventHeader( getDaqEvent().GetTime().first,              6 );
    _recoEvent.setEventHeader( getDaqEvent().GetTime().second,             7 );
    //_recoEvent.setEventHeader( _daqEvent->GetErrorCode(),               8 ); // Commented out in order to cope w/ Date,v5 data
    _recoEvent.setEventHeader( getDaqEvent().GetHeader().GetTrigger(),     9 );

  }


# ifdef TIME_CHECK
  if(TiDecod){ dec_times=stopwatch.stop(star); n_dig=getChipDigits().size(); }
# endif  

  return( true );

}

bool CsEvent::_clusterize() {

  // Make clusters if needed
  if( ! _clustering ) return( true );

  int buildMode;
  if( _MCClusteringMode == "MCExact" ) {
    // Cluster correspond exactly to hit
    buildMode = 1;
  }
  else if( _MCClusteringMode == "MCSmeared" ) {
    // Cluster correspond exactly to smeared hit
    buildMode = 2;
    }
  else if( _MCClusteringMode == "MCQuantized" ) {
    // Cluster correspond exactly to "quantized" (by det. wire) hit
    buildMode = 3;
  }
  else {
    // Cluster built from digit.
    _MCClusteringMode = "Standard";
    buildMode = 0;
  }

  // An alternative master (i.e. defining reference T0) trigger may be defined
  // infra in the case of real data ("standard clustering"). Init it to -1 (
  // meaning no alternative) in any case.
  _alterMasterTrigger = -1;

  if (buildMode==0) {  // *************** STANDARD CLUSTERING ***************
    if( CsInit::Instance()->IsADataJob()) {
    
      //      *************** TIME GATE for CURRENT EVENT ***************
      // The expectation for the jitter has been refined, now that all data
      // including the trigger pattern TDC have been decoded: the time gate can
      // possibly be decreased, prior to the clusterisation step.

      // Calo based triggers have bad timing => Single them out (cf. new run#
      // block infra).
      static unsigned int caloBasedTrigs = 0;
      static unsigned run_old = unsigned(-1);
      if (run_old!=getRunNumber()) {  // New run number
	const CS::DaqOption &opts = CsInit::Instance()->getDaqEventsManager().GetDaqOptions();
	int bit = -1;// High-Q2 trigger: Exists? Get its bit.
	try { bit = opts.GetTriggerBit("CalorimeterTrigger"); } catch (...) {}
	if (bit>=0) caloBasedTrigs |= 1<<bit;
	bit = -1;// High-Q2 trigger: Exists? Get its bit.
	try { bit = opts.GetTriggerBit("LargeQ2Trigger"); } catch (...) {}
	if (bit>=0) caloBasedTrigs |= 1<<bit;
	run_old = getRunNumber();
      }
      unsigned int triggerMask = getTriggerMask(); // Init w/ trigger pattern
      const unsigned int allTrigs = 0xffff;
      triggerMask &= allTrigs; // Strip away on-line filter
      // Case of combination of caloBased & hodoscope triggers
      bool hasCaloAndHodo = (triggerMask&caloBasedTrigs) &&
	triggerMask!=(triggerMask&caloBasedTrigs);

      //                                                ***** GET MASTER TRIGGER
      // (I.e. trigger used to set master time by the decoding library on view
      // of the trigger pattern TDC.)
      const CS::Trigger *masterTrig = getDaqEvent().GetTT().GetTriggerMT();
      int masterTrigBit = -1;
      if (!masterTrig)
	CsErrLog::msg(elError,__FILE__,__LINE__,
	  "Master trigger empty while trigger pattern = 0x%x",triggerMask);
      else {
	masterTrigBit = (int)masterTrig->GetBit();
	if (masterTrigBit<0 || !(1<<masterTrigBit&triggerMask)) {
	  CsErrLog::msg(elError,__FILE__,__LINE__,"Master trigger (=%d) not in trigger pattern (=0x%x)",masterTrigBit,triggerMask);
	  masterTrigBit = -1;
	}
      }
      //#define CsEvt_DEBUG_MASTER_TIME 2
#ifdef CsEvt_DEBUG_MASTER_TIME
      if (masterTrigBit>=0) printf("Master Trigger: %s bit=%d  precision=%gns\n",masterTrig->GetName().c_str(),masterTrigBit,masterTrig->GetPrecision());
#endif

      //                                         ***** CHECK TRIGGER PATTERN TDC
      // - Singling out cases where single hodo-based TDC (single TDC hit on
      //  single hodo-trigger bit). The corresponding hodo trigger may be
      //  selected as an alternative master trigger, cf. infra.
      // - Retrieving the time diff of that TDC w.r.t. master trigger.
      //  => Setting "_triggerTimeCorr" accordingly.
      // - Checking fired TDC channels correspond to trigger mask.
      const CS::DaqOption &opts = CsInit::Instance()->getDaqEventsManager().GetDaqOptions();
      const CS::DetID *trigPatTDCId = &opts.GetTTConfig().trigger_mask_DetID;
      if (!trigPatTDCId) CsErrLog::mes(elFatal,
	"CsEvt_DEBUG_MASTER_TIME: No CsDet corresponding to triggerPatternTDC");
      typedef CS::Chip::Digits::const_iterator m_it; // Iterator type
      pair<m_it,m_it> m_range = getChipDigits().equal_range(*trigPatTDCId);
      m_it d_it; unsigned int bitsPat; int hodoTrigBit;
#if defined CsEvt_DEBUG_MASTER_TIME && CsEvt_DEBUG_MASTER_TIME > 1
      printf("Trigger pattern TDC:");
#endif
      int iter; for (iter = 0, _triggerTimeCorr = 0; iter<2; iter++) {
	// 1st iter: Get TDC of master trigger: store it.
	// 2nd iter: Get all TDCs
	static double masterTDC;
	for (d_it = m_range.first, bitsPat = 0, hodoTrigBit = -1;
	     d_it!=m_range.second; d_it++) {
	  const CS::ChipF1::Digit *f1 =
	    dynamic_cast<const CS::ChipF1::Digit*>(d_it->second);
	  int bit = f1->GetChannel(); double tdc = f1->GetTimeDecoded();
	  if (iter) {
	    unsigned int pat = 1<<bit; bitsPat |= pat;
#if defined CsEvt_DEBUG_MASTER_TIME && CsEvt_DEBUG_MASTER_TIME > 1
	    printf(" 0x%x %.2f %.2f",pat,tdc,tdc-masterTDC);
#endif
	    if (fabs(tdc-masterTDC)<8 && !(pat&caloBasedTrigs)) {
	      if (hodoTrigBit==-1) hodoTrigBit = bit;
	      else                 hodoTrigBit = -2;
	      _triggerTimeCorr += tdc;
	    }
	  }
	  else {
	    if (bit!=masterTrigBit) continue;
	    masterTDC = tdc; _triggerTimeCorr -= tdc; break;
	  }
	}
      }
#if defined CsEvt_DEBUG_MASTER_TIME && CsEvt_DEBUG_MASTER_TIME > 1
      printf("\n");
#endif
      if (bitsPat!=triggerMask) CsErrLog::msg(elError,__FILE__,__LINE__,
 "Trigger pattern TDCs 0x%x != trigger pattern 0x%x\n",bitsPat,triggerMask);
      if (hodoTrigBit<0) // Too complex a case: defining no hodo based time...
	_triggerTimeCorr = 0; // ... => Reset "_triggerTimeCorr"

      if (hasCaloAndHodo &&    // ***** COMBINATION caloBased & hodo TRIGGERS...
	  masterTrigBit>=0) {        // ...while master trigger is defined
	//                                        ***** ...REDEFINE TIME GATE...
	if (1<<masterTrigBit&caloBasedTrigs) {// I) Master trigger is calo-based
	  // If unique hodo trigger can be defined, let's use it instead
	  if (0<hodoTrigBit) {
	    if (hodoTrigBit<(int)_extraTimeWidths.size())
	      _extraTimeWidth =_extraTimeWidths[hodoTrigBit];
	    _alterMasterTrigger = hodoTrigBit;
	    // And "_triggerTimeCorr" wil be used to correct the offset of the
	    // hodo trigger w.r.t. master trigger.
	  }
	}
	else {                               // II) Master trigger is hodo-based
	  _extraTimeWidth = masterTrigBit<(int)_extraTimeWidths.size() ? 
	    _extraTimeWidths[masterTrigBit] : 0;
	  _triggerTimeCorr = 0; // ... => Reset "_triggerTimeCorr"
	}
      }
      else _triggerTimeCorr = 0; // ... => Reset "_triggerTimeCorr"
    }

    list<CsDetector*> dets = CsGeom::Instance()->getDetectors();
    list<CsDetector*>::iterator i;
#ifdef TIME_CHECK  
   CsStopwatch stopwatch;
#endif  
    for( i=dets.begin(); i!=dets.end(); i++ ) {
#ifdef TIME_CHECK
     string Det_tmp = (*i)->GetTBName();
     string Detekt;
     Detekt.append(Det_tmp,0,2); 	  
     int star=stopwatch.start(); 
#warning "Histograms for 'Standard' clustering (data only?) mode"                              
#endif
      (*i)->clusterize();
#ifdef TIME_CHECK    
if(TiClust){
     if(!Detekt.compare("BM")){
        clu_times[0] +=stopwatch.stop(star);
//	clu_nb[0] +=((*i)->getMyClusters()).size();
       }
     else if(!Detekt.compare("MM")){
        clu_times[1] +=stopwatch.stop(star);
	clu_nb[1] +=((*i)->getMyClusters()).size();
	}
     else if(!Detekt.compare("DC")){
        clu_times[2] +=stopwatch.stop(star);
	clu_nb[2] +=((*i)->getMyClusters()).size();
	}
     else if(!Detekt.compare("ST")){
        clu_times[3] +=stopwatch.stop(star);
	clu_nb[3] +=((*i)->getMyClusters()).size();
        }
     else if(!Detekt.compare("MB")){
        clu_times[4] +=stopwatch.stop(star);
	clu_nb[4] +=((*i)->getMyClusters()).size();
	}
     else if(Detekt.compare("GM")==0 || Detekt.compare("GP")==0){
        clu_times[5] +=stopwatch.stop(star);
	clu_nb[5] +=((*i)->getMyClusters()).size();
	}
     else if(!Detekt.compare("FI")){ 
        clu_times[6] +=stopwatch.stop(star);
	clu_nb[6] +=((*i)->getMyClusters()).size();
	}         
     else if(Detekt.compare("PS")==0 || Detekt.compare("PA")==0 || Detekt.compare("PB")==0){     
        clu_times[7] +=stopwatch.stop(star);
	clu_nb[7] +=((*i)->getMyClusters()).size();
	}
     else 
     CsErrLog::Instance()->mes(elFatal, 
       "CsEvent::_clusterize: Time debugging do not recognize detector type" );
}
#endif
    }
    if( _doClustersAssociation ) CsEventUtils::associateClusters(); 
  }
  else {
    _mkClustersFromHits( buildMode );
  }

  //store all clusters in _recoEvent

  list<CsDetector*> dets = CsGeom::Instance()->getDetectors();
  for( list<CsDetector*>::iterator idt=dets.begin(); idt!=dets.end(); idt++ ) {
    list<CsCluster*> clus = (*idt)->getMyClusters();
    for(list<CsCluster*>::iterator idc=clus.begin();idc!=clus.end();idc++) { 
      _recoEvent.addCluster( *(*idc) );
    }
  }

  // Sort Clusters
  sortClusters();

  return( true );

}

bool CsEvent::_reconstructionSchemaPP() {

#if USE_MySQL
  CsPP::Instance()->ProcessEvent();

  return true;
#else
  CsErrLog::Instance()->mes(elError,"Calling preprocessor, but CsPPI/MySQL DB not enabled");

  return false;
#endif

};


bool CsEvent::_reconstructionSchema001()
{
  bool status = false;
  
  // Decoding is out of reco. schema. BG 2003/06/26

  //   *************** CLUSTERIZE TRACKING DETECTORS ***************
  _clusterizeTime -= _chronometer.inter(_chrono);
  status = _clusterize();
  _clusterizeTime += _chronometer.inter(_chrono); _nclusterizeTime++;

  if (!isAMonteCarloEvent()) {
    static bool skippingBoS = false;   // *************** BoS ***************
    if (getTimeInSpill()<CsInit::Instance()->getMinTimeInSpill()) {
      if (!skippingBoS) cout<<"Skipping BoS..."<<std::flush;
      skippingBoS = true; return true;
    }
    if (skippingBoS) { skippingBoS = false; cout << endl; }
    _hntimeinspill->Fill( getTimeInSpill() );
  }

  if (_DTEvent) {
    //    *************** FILL FLUX/SPILL HISTOGRAM ***************     
    int spill = getDaqEvent().GetBurstNumber(); unsigned int nsf1;
    if (getRunNumber()<24885) // If in 2002..
      nsf1 = getSC99P2_16() + getSC99P2_17() + getSC99P2_18() +
	getSC99P2_19() + getSC99P2_20() + getSC99P2_21();     
    else
      nsf1 = _recoEvent.getFluxX();
    int nsincelast = nsf1<_nmuinspillsofar ? nsf1 /* Scaler was reset */ :
      nsf1 - _nmuinspillsofar;
    _hnmuperspill->Fill(spill,nsincelast);
    _nmuinspillsofar = nsf1;
  }

  int iter; for (iter = 0, _reTrackingON = false; iter<2; iter++) {
    // *********************************************************************
    // ******************** 2 ITERATIONS ON RECO SCHEMA ********************
    //  I) T0 = 0
    // II) T0 = beam track's time of (I), if CsVrtPattern::doPattern decides so,
    //    i.e., typically, if beam track's time large
    if (iter) {
      if (_reTrackingON) {      // Upon re-tracking...
	_recoEvent.reset();     // ...erase all products of reconstruction
	if (_MCEvent) {         // ...erase X-references to MC tracks
	  list<CsMCTrack*> MCTracks = CsGeant3::Instance()->getMCTracks();
	  list<CsMCTrack*>::iterator mct = MCTracks.begin();
	  while (mct!=MCTracks.end()) {
	    (*mct)->clearAssociatedTrackID(); mct++;
	  }
	}
      }
      else break;
    }

    if (_tracking) {  // *************** TRACK reconstruction ***************
      _trackingTime -= _chronometer.inter(_chrono);

      list<CsZone*>    zones = CsGeom::Instance()->getZones();
      list<CsCluster*> unusedClusters;
      list<CsTrack*>   tracks;
      list<CsCluster*> clusters;// Build the list of clusters to be used...
      clusters = getClusters();             // ...At the moment use all clusters

      if (!clusters.empty()) {
	if (_prepattern1->doPrepattern(clusters,zones)) {      // Do prepattern
	  status  = _prepattern1->getPatterns( tracks );
	  if (status) setTracks(tracks);     // Fill the event track list
	  unusedClusters = _prepattern1->getUnusedClusters();
	  if (_bridging->doBridging(tracks,unusedClusters)) {    // Do bridging
	    setTracks(tracks);                 // Update the event tracks list
	    if (_fitting->doFitting(tracks,unusedClusters)) {      // Do fitting
	      if (status) setTracks(tracks);     // Update the event tracks list
	    }
	    else CsErrLog::Instance()->mes(elError,"Track fitting problems");
	  }
	  else CsErrLog::Instance()->mes(elError,"Track bridging problems");
	}
	else CsErrLog::Instance()->mes(elError,"Track prepattern problems");
	// IMPORTANT! clean local tracks: CsRecoEvent has a copy of them
	list<CsTrack*>::iterator It;
	for (It=tracks.begin(); It!=tracks.end(); It++) {
	  delete *It;
	}
	tracks.clear();
      }
      else {
	CsErrLog::Instance()->mes(elError,"No clusters"); status = true;
      }
      _trackingTime += _chronometer.inter(_chrono); if (!iter) _nTrackingTime++;
    }
    else status = true;

    if (_beam) {      // *************** BEAM reconstruction ***************
      _beamTime -= _chronometer.inter(_chrono);
    
      if(_beamMethod == 0){
         CsBeamRecons *beamReco = CsBeamRecons::Instance();
         if (_beamTR) beamReco->bmreconsTR(_recoEvent.getTracks());  // Use TRAFFIC..
         else         beamReco->bmrecons();                          // ...or dedicated beam package
         _recoEvent.setBeamTracksList(beamReco->getBeam());          // store the pointers to the reconstructed beam tracks
      }

      if(_beamMethod == 1){
         CsBeamReconstruction *beamReco = CsBeamReconstruction::Instance();
         beamReco->make_beam_reconstruction();
         _recoEvent.setBeamTracksList(beamReco->getBeam());
      }

      _beamTime += _chronometer.inter(_chrono); if (!iter) _nbeamTime++;    
    }

    if (_rich1) {      // *************** RICH reconstruction ***************
      _rich1Time -= _chronometer.inter(_chrono);

      CsRichOne::Instance()->doRichOne(); // Expecting direct output to tracks

      _rich1Time += _chronometer.inter(_chrono); if (!iter) _nrich1Time++;
    }

    // *************** CEDAR reconstruction ***************
    if (_cedars) {
      _cedarTime -= _chronometer.inter(_chrono);
      list<CsTrack*>::iterator It;
      for (It=_recoEvent.getTracks().begin(); It!=_recoEvent.getTracks().end(); ++It)
        (*It)->addCEDARInfo();
      _cedarTime += _chronometer.inter(_chrono); if (!iter) _ncedarTime ++;
    }

    if (_calorimeters) {  // *************** CALO reconstruction ***************
      _caloTime -= _chronometer.inter(_chrono);

      vector<Reco::CalorimeterParticle> particles_from_calorimeters;
      vector<CsCalorimeter*> calos = CsGeom::Instance()->getCalorimeters();
      for (vector<CsCalorimeter*>::iterator itc = calos.begin();
	   itc!=calos.end(); itc++) {
	vector<Reco::CalorimeterParticle> particles_rec;
	(*itc)->Reconstruction();
	particles_rec = (*itc)->GetCalorimeterParticles();
	particles_from_calorimeters.insert(particles_from_calorimeters.end(),
					   particles_rec.begin(), 
					   particles_rec.end() );
      } 

      _recoEvent.setCalObjsVector( particles_from_calorimeters );

      _caloTime += _chronometer.inter(_chrono); if (!iter) _ncaloTime ++;
    }

    //         ******************** BUILD FINAL PARTICLES ********************

    vector<CsParticle*> particles;

    // Track + Calo Cluster association
    _partTime -= _chronometer.inter(_chrono);
    CsBuildParticles::Instance()->Build(_recoEvent, particles);
    //    CsBuildParticles::Instance()->SetClusterZ( particles );
    _partTime += _chronometer.inter(_chrono); if (!iter) _npartTime++;

    _pidTime -= _chronometer.inter(_chrono);

    PID_doMuonID(particles);       // (so-called) scattered muon ID
    PID_doBeamID(particles);       // Beam ID

    setParticles( particles );     // Update of event particles vector 

    if (_RWcalorimeters) {         // Combined Calorimeter+RICHwall reco
      CsRwRecons::Instance()->RwRecons();
      int len = CsRwRecons::Instance()->dataout.size();
      if (len>0) {
	for (int i = 0; i<len; i++) {
	  int j = CsRwRecons::Instance()->dataout[i]->npart;
	  particles[j]->addCalobj(CsRwRecons::Instance()->dataout[i]->calobj);
	}
      }
    }
    if (_RWchargeEcal1){
      CsRwChargeRecons::Instance()->RwChargeRecons();
    }

    _pidTime += _chronometer.inter(_chrono); if (!iter) _npidTime++;

    //          ******************** VERTEX reconstruction ********************
    if (_vertex) {
      _vertexTime -= _chronometer.inter(_chrono);

      if (!particles.empty()) {
	list<CsVertex*> vrts; map<CsTrack*,bool> specials;
    
	muIDinMW1(particles);  // MW1 muID (should be moved to "PIDdoMuonID"?)

	double *reTrackT0; // Request for a re-tracking. "doPattern" is...
	// ...expected to store there the T0 to be used in a 2nd pass of
	// tracking and reconstruction (typically the beam track time of the
	// best vertex, when it turns out to differ significantly from 0)
	// provided it is not a null pointer. Therefore...
	reTrackT0 = iter?0:&_reTrackT0;// ...set = null, if already in 2nd pass.

	bool vPatternOK =             // ***** DO VERTEX P(attern) R(ecognition)
	  _vpattern->doPattern(particles,reTrackT0);

	// get list of vertices in any case to be able to clear it later
	vPatternOK &= _vpattern->getPatterns(vrts,specials);

	if (reTrackT0 && *reTrackT0) {// If re-tracking requested and granted...
	  // ...validate "_reTrackT0" w/ boolean (more robust)
	  _reTrackingON = true;
	  vPatternOK = false;    // ...give up current reconstruction.
	}
	if (vPatternOK) {
	  // The T0 used in the final tracks fit, be it via the original
	  // fitting or via re-fitting) was decided from an analysis on the
	  // vertices from PR, before vertex fitting. => Have to check that
	  // this decision is compatible w/ fitted vertices. Let us prepare
	  // the ground for this check.
	  bool t0FromPR = false;
	  double trackFittingT0 = 0; // T0 eventually used when fitting tracks
	  const CsVertex *v; if ((v = _vpattern->getT0SettingVertex())) {
	    if (_reTrackingON) // When in the 2nd pass of event reco...
	      // ..."_vpattern::doPattern" should not try to reset tracking T0
	      // (and hence request a refitting the tracks).
	      CsErrLog::mes(elFatal, "Refitting of tracks was requested "
			    "while in 2nd pass of event reconstruction");
	    const std::list<CsTrack*> tracks = v->getTracks();
	    if (!tracks.empty() && tracks.front()->hasMeanTime()) {
	      trackFittingT0 = tracks.front()->getMeanTime(); t0FromPR = true;
	    }
	  }
	  if (_vfitting) {                          // ***** DO VERTEX FITTING
	    double T0 = // On input to "doFitting", argument T0 must hold the
	      // value of T0 eventually used in tracking.
	      _reTrackingON ? _reTrackT0 : trackFittingT0;
	    if (!_vfitting->doFitting(vrts,specials,
				     _reTrackingON,&T0)) {
	      _reTrackingON = true; _reTrackT0 = T0;
	    }
	  }
	  else
	    CsErrLog::mes(elWarning,"Vertex fitting problems");
	}
	// set list of vertices in any case, otherwise the vertices might never
	// be deleted
	setVertices(vrts);
      }
      _vertexTime += _chronometer.inter(_chrono); if (!iter) _nvertexTime++;
    }
  }
  CsRichOne::Instance()->physSelec();

  return status;
}


//----------------------------------------------------------------

bool CsEvent::_reconstructionSchema002()
{
  bool status = false;

  // decode out of reco. schema. BG 2003/06/26

  // Clusterize (find hits of) tracking detectors

  _clusterizeTime -= _chronometer.inter(_chrono);

  status =_clusterize();

  _clusterizeTime += _chronometer.inter(_chrono);
  _nclusterizeTime++;


  // Tracks Reconstruction
  if( _tracking ) {

    _trackingTime -= _chronometer.inter(_chrono);

    list<CsZone*>    zones = CsGeom::Instance()->getZones();
    list<CsCluster*> unusedClusters;
    list<CsTrack*>   tracks;
   
    // Build the list of clusters to be used
    //list<CsCluster*> clusters = getClusters();

    list<CsCluster*> clusters;

    // At the moment use all clusters...

    if( false ) {

      list<CsDetector*> dets = CsGeom::Instance()->getDetectors();
      for(list<CsDetector*>::iterator id=dets.begin();id!=dets.end();id++ ) {
	list<CsCluster*> mycls = (*id)->getMyClusters();
	for(list<CsCluster*>::iterator ic=mycls.begin();ic!=mycls.end();ic++) { 
	  double time;
	  if( (*ic)->getTime( time ) ) { // This cluser has time...
	    // select cluster 
	    if( (*id)->GetTBName().substr( 0, 2 )=="FI" && 
		fabs(time)<9999999.){
	      clusters.push_back( *ic );
	    }
	    else {
	      clusters.push_back( *ic );
	    }
	  }
	  else { // This cluster has no time value set...
	    clusters.push_back( *ic );
	  }
	}
      }
    }
    else {
      clusters = getClusters();
    }

    if( !clusters.empty() ) {
      // do prepattern
      if( _prepattern1->doPrepattern( clusters, zones ) ) {
	status  = _prepattern1->getPatterns( tracks );
	// fill the event track list
	if( status ) setTracks( tracks );
	unusedClusters = _prepattern1->getUnusedClusters();
	// do bridging
	if( _bridging->doBridging( tracks, unusedClusters ) ) {
	  // update the event tracks list 
	  if( status ) setTracks( tracks );
	  // do fitting
	  if( _fitting->doFitting( tracks, unusedClusters ) ) {
	    // update the event tracks list 
	    if( status ) setTracks( tracks );
	  }
	  else {
	    CsErrLog::Instance()->mes( elError, "Track fitting problems" );
	  }
	}
	else {
	  CsErrLog::Instance()->mes( elError, "Track bridging problems" );
	}
      }  
      else {
	CsErrLog::Instance()->mes( elError, "Track prepattern problems" );
      }

      // do prepattern for Recon
      if( _prepattern2->doPrepattern( clusters, zones ) ) {
	status  = _prepattern2->getPatterns( tracks );
	// fill the event track list
	if( status ) setTracks( tracks );
      }
      

      // IMPORTANT! clean local tracks: CsRecoEvent has a copy of them
      list<CsTrack*>::iterator It;
      for( It=tracks.begin(); It!=tracks.end(); It++ ) {
	delete *It;
      }
      tracks.clear();
    }
    else {
      CsErrLog::Instance()->mes( elError, "No clusters" );
      status = true;
    }

    _trackingTime += _chronometer.inter(_chrono);
    _nTrackingTime++;

  }
  else {
    status = true;
  }
  
  // Beam reconstruction
  if( _beam ) {
    
    _beamTime -= _chronometer.inter(_chrono);
    
    if(_beamMethod == 0){
      CsBeamRecons* beamReco = CsBeamRecons::Instance();
      if(_beamTR!=0) beamReco->bmreconsTR(_recoEvent.getTracks());
      else beamReco->bmrecons();
      _recoEvent.setBeamTracksList( beamReco->getBeam() );
    }
    
    if(_beamMethod == 1){
       CsBeamReconstruction *beamReco = CsBeamReconstruction::Instance();
       beamReco->make_beam_reconstruction();
       _recoEvent.setBeamTracksList(beamReco->getBeam());
    }

    _beamTime += _chronometer.inter(_chrono);
    _nbeamTime++;
    
  }

  // RICH1 reconstruction
  if( _rich1 ) {

    _rich1Time -= _chronometer.inter(_chrono);

    // should store directly data to tracks
    CsRichOne::Instance()->doRichOne();

    _rich1Time += _chronometer.inter(_chrono);
    _nrich1Time++;

  }

  // Calorimeters reconstruction
  if( _calorimeters ) {

    _caloTime -= _chronometer.inter(_chrono);

    vector<Reco::CalorimeterParticle> particles_from_calorimeters;
    vector<CsCalorimeter*> &calos = CsGeom::Instance()->getCalorimeters();
    for( vector<CsCalorimeter*>::iterator itc = calos.begin(); itc!=calos.end(); itc++ ) 
    {
	vector<Reco::CalorimeterParticle> particles_rec;
	(*itc)->Reconstruction();
	particles_rec = (*itc)->GetCalorimeterParticles();
	particles_from_calorimeters.insert(particles_from_calorimeters.end(),
					   particles_rec.begin(), 
					   particles_rec.end() );
    } 
    _recoEvent.setCalObjsVector( particles_from_calorimeters );

    _caloTime += _chronometer.inter(_chrono);
    _ncaloTime ++;

  }


  // Build final particles
  
  _pidTime -= _chronometer.inter(_chrono);

  vector<CsParticle*> particles;

  // Track + Calo Cluster association
  _partTime -= _chronometer.inter(_chrono);
  CsBuildParticles::Instance()->Build(_recoEvent, particles);
  //  CsBuildParticles::Instance()->SetClusterZ( particles );
  _partTime += _chronometer.inter(_chrono);
  _npartTime++;


  // Muon wall information is used for muons selection
  PID_doMuonID( particles );

  // mu' information is used for beam selection
  PID_doBeamID( particles );

  // merging of cal clusters with tracks
  // ==> Now done directly in _buildparticles. BG 2002/08/07
  //if( _calorimeters ) PID_doCalID( particles );

  // update of the event particles vector 
  setParticles( particles );

  // Combined Calorimeter+RICHwall reco
  if(_RWcalorimeters) {
    CsRwRecons::Instance()->RwRecons();
    int len = CsRwRecons::Instance()->dataout.size();
    if(len>0){
      for(int i=0;i<len; i++){
	int j=CsRwRecons::Instance()->dataout[i]->npart;
	particles[j]->addCalobj(CsRwRecons::Instance()->dataout[i]->calobj);
      }
    }
  }

  if (_RWchargeEcal1){
    CsRwChargeRecons::Instance()->RwChargeRecons();
  }

  _pidTime += _chronometer.inter(_chrono);
  _npidTime++;

  // vertex reconstruction
  if( _vertex ) {

    _vertexTime -= _chronometer.inter(_chrono);
    
    list<CsVertex*> vrts;
    map<CsTrack*,bool> specials;
    
    if( !particles.empty() ) {
    
      // do prepattern
      if( _vpattern->doPattern( particles ) ) {
	
	if( _vpattern->getPatterns( vrts, specials ) ) {
	  
	  if( _vfitting != 0 ) {
	    // do fitting
	    if( _vfitting->doFitting( vrts, specials ) ) {
	      
	      // update the event vertices list 
	      setVertices( vrts );
	      
	    } else
	      CsErrLog::Instance()->mes( elWarning, "Vertex fitting problems" );  
	  }
	}

      } else
	CsErrLog::Instance()->mes( elWarning, "Vertex pattern problems" );
    }

    _vertexTime += _chronometer.inter(_chrono);
    _nvertexTime++;

  }

  return status;
}


//----------------------------------------------------------------

void CsEvent::_setTrackingPackages() {

  // PREPATTERN
  string tag = "", key = "track prepattern method", pkg;
  if (CsOpt::Instance()->getOpt(tag,key,pkg)) {
    if (pkg=="traffic") _prepattern1 = new CsTrafficPrepattern;
    else
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
	"Track prepattern package \"%s\" not available.",pkg.c_str());
  }
  else
    CsErrLog::mes(elFatal,"No track prepattern package specified!" );

  // BRIDGING
  key = "track bridging method";
  if (CsOpt::Instance()->getOpt(tag,key,pkg)) {
    if (pkg=="traffic") _bridging = new CsTrafficBridging;
    else
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
	"Track bridging package \"%s\" not available.",pkg.c_str());
  }
  else
    CsErrLog::mes(elFatal,"No track bridging package specified!" );

  // FITTING
  key = "track fitting method";
  if (CsOpt::Instance()->getOpt(tag,key,pkg)) {
    if (pkg=="traffic") _fitting = new CsTrafficFitting;
    else
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
	"Track fitting package \"%s\" not available.",pkg.c_str());
  }
  else
    CsErrLog::mes(elFatal,"No track fitting package specified!" );
};


//---------------------------------------------------------

void CsEvent::_setVertexPackages() {

  // PATTERN
  string tag = "";
  string key = "vertex pattern method";
  string pkg;
  if( CsOpt::Instance()->getOpt( tag, key, pkg ) ) {
    if( pkg == "averaging" ) {
      _vpattern = new CsAverPattern();
    }
    else if( pkg == "roland" ) {
      _vpattern = new CsRolandPattern();
    }
    else {
      CsErrLog::Instance()->mes( elFatal, 
				 "Vertex pattern package not available" );
    }
  }
  else {
    CsErrLog::Instance()->mes( elFatal, "No vertex pattern package set" );
  }

  // FITTING
  key = "vertex fitting method";
  if( CsOpt::Instance()->getOpt( tag, key, pkg ) ) {
    if( pkg == "kalman" ) {
      _vfitting = new CsKalmanFitting();
    }
    else {
      CsErrLog::Instance()->mes( elFatal, 
				 "Vertex fitting package not available" );
    }
  }
  else {
    _vfitting = 0;
    CsErrLog::Instance()->mes( elWarning, "No vertex fitting package set" );
  }

}



void CsEvent::_mkClustersFromHits( int mode ) {

  // At the moment do all clusters, regardless what's set...

  // set random number generator (seed event dependent)
  CsGauss c; 
  c.setSeed( getEventNumberInRun() | 69 );

  list<CsDigit*>::iterator Id;
  list<CsDigit*> digits = getDigits();

  // a protection
  if( digits.empty() ) return;

  for( Id=digits.begin(); Id!=digits.end(); Id++ ) {
    // At present, make a different cluster for each digit...
    CsDet* dtc = (*Id)->getDet();
    CsDetector* det = dynamic_cast<CsDetector*>(dtc);

    if( det == 0 ) continue; // only tracking devices

    double      ang  = det->getAng();
    double      wirD = det->getWirD();
    double      wirP = det->getWirP();
    int         nWir = det->getNWir();
    double      zcm  = det->getZcm();
    bool        skip = false;

    CsDigit* dig = *Id;
    CsMCDigit* mcdig = dynamic_cast<CsMCDigit*>(dig);
    if( mcdig == 0 ) continue;               // should never happen
    CsMCHit* hit=(mcdig->getHits()).front(); // There MUST be a single Hit
    
    double x = hit->getX(); // hit coordinates
    double y = hit->getY();
    double z = hit->getZ();
    double t = hit->getDTime(); // hit delay time

    HepMatrix rotM(3,3); rotM = det->getRotWRS();        

    // Rotate back ?????
    int err;
    HepMatrix irotM(3,3); irotM = rotM.inverse( err );
    double u = irotM(1,1)*x+irotM(1,2)*y+irotM(1,3)*z; // WRS
    double v = irotM(2,1)*x+irotM(2,2)*y+irotM(2,3)*z;
    double w = irotM(3,1)*x+irotM(3,2)*y+irotM(3,3)*z;

    // Set errors: Unitary matrix...
    HepMatrix cov(3,3,0); 

    // Detector Resolution
    double res;
    if( det->hasDrift() ) {
      res = (det->getVel())*(det->getSpSli()); // detector resolution
    }
    else {
      res = det->getWirP() / sqrt(12.);  // wire pitch / sqrt(12)
    }
    cov(1,1) = res*res;
    if( ang == 90 ) {           // wire length/2
      cov(2,2) = pow( (det->getXsiz())/2., 2 );
    }
    else {
      cov(2,2) = pow( (det->getYsiz())/2./cos(ang/180.*(M_PI)), 2 );
    }
    cov(3,3) = 1.;               // assume 1 mm resolution in Z

    double uDrift = 0;
    // if buildMode == 1 Keep exact hits coordinates
    if( mode == 2 ) { // Make "Smeared" clusters
      u = u + res * c.random();
    }
    else if( mode == 3 ) { // Make "Quantized" clusters
      double uExact = u;
      int wire = int( (u-wirD) / wirP + 0.5 ); // wire for this hit
      if( wire < 0 || wire >= nWir ) {
	skip = true;
	ostringstream ost;
	ost << "Unreliable wire number: "<<wire<<" (0,"<<nWir<<").";
	CsErrLog::Instance()->mes( elError, ost.str() );
      }
      else {
	skip = false;
	u = wirD + wire * wirP;
	v = 0;
	w = zcm;
	if( det->hasDrift() ) { // drift detectors
	  double drift  = fabs( uExact - u );
	  do {
	    double random = c.random();
	    uDrift = drift + res * random;
	  } while( ( uDrift<0 || uDrift>(wirP/2) ) 
		   && wire != 0 && wire != nWir  ); 
	}
      }
    }
	
    // Save the cluster(s):
    if( !skip ) {
      if( mode != 3 || !det->hasDrift() ) { 
	// single cluster for non drift 
	CsCluster* cluster = new CsCluster( u, v, w, cov );
	cluster->addDigit( *(*Id) );
	cluster->addDet( *det );
	cluster->setTime( t );
	addCluster( *cluster );
      }
      else if( mode == 3 && det->hasDrift() ) { 
	// 2 clusters for drift
	CsCluster* cluster = new CsCluster( u+uDrift, v, w, cov );
	cluster->addDigit( *(*Id) );
	cluster->addDet( *det );
	addCluster( *cluster );
	
	cluster = new CsCluster( u-uDrift, v, w, cov );
	cluster->addDigit( *(*Id) );
	cluster->addDet( *det );
	addCluster( *cluster );
      }
      else {
	// should NEVER happen, but...
	string str = "A strange ERROR happens. Please check code...";
	CsErrLog::Instance()->mes( elFatal, str );
      }
    }
  }
}


void CsEvent::_upload(void) {

# if CHECK_TIMEUSAGE
  static bool firstevent = true;
  if( firstevent ) {
    timeusage( "1ST EVENT UPLOAD started", 0 );
  }
# endif

  // simple instantiation of a persistent event
  if( _store  ) {
    if( !_store->uploadDST( &_recoEvent, _dstfatness, _dstversion ) ) {
      CsErrLog::Instance()->mes(elFatal, "Error in DST upload." );
    }
  }

# if CHECK_TIMEUSAGE
  if( firstevent ) {
    firstevent = false;
    timeusage( "1ST EVENT UPLOAD ended", 1 );
  }
# endif

}


bool CsEvent::_download(void) {

  if( _store  ) {
    _recoEvent.clearBeamTracks();
    _recoEvent.clear();
    if( !_store->downloadDST( &_recoEvent, _dstversion ) ) {
      CsErrLog::Instance()->mes(elError, "Error in DST download." );
      return false;
    }
  }
  return true;

}


void CsEvent::testDump( bool alsoClus ) {

  cout << "-------------------------------------------" << endl
       << "    EVENT " << getEventNumberInRun() << " DUMP" << endl
       << "-------------------------------------------" << endl;
  cout << "Scalers Data" << endl
       << "SC99P2_16: " << getSC99P2_16() 
       << ", SC99P2_17: " << getSC99P2_17()
       << ", SC99P2_18: " << getSC99P2_18() << endl
       << "SC99P2_19: " << getSC99P2_19()
       << ", SC99P2_20: " << getSC99P2_20()
       << ", SC99P2_21: " << getSC99P2_21() << endl
       << "SC01P1_01: " << getSC01P1_01()
       << ", SC01P1_09: " << getSC01P1_09() << endl;
  cout << "Event Header" << endl
       << "Size: " << getEventSize() << endl
       << "Run #: " << getRunNumber() 
       << ", # in Run: " << getEventNumberInRun() << endl
       << "Burst #: " << getBurstNumber()
       << ", # in Burst: " << getEventNumberInBurst() << endl
       << "Trigger #: " << getTriggerNumber() 
       << ", Time: " << getTime().first 
       << "." << getTime().second << endl
       << "Error Code: " << getErrorCode() 
       << ", Trigger Mask: " << getTriggerMask() << endl;
  cout << "Trigger Time: " << getTriggerTime() << endl;
  cout << "Hodoscopes Data" << endl;
  string tbname; int channel; double time;
  for( unsigned int i=0; i<=getHodoDataSize(); i++ ) {
    getHodoDatum( i, tbname, channel, time );
    cout << tbname << ", ch.: " << channel << ", time: " << time << endl;
  }
  cout << "-------------------------------------------" << endl;

  int npart = 0;
  vector<CsParticle*> part = _recoEvent.getParticles();
  for( unsigned int i=0; i<part.size(); i++ ) {
    cout << "Particle " << ++npart 
	 << ", charge: " << part[i]->getCharge() 
	 << ", PID: " << part[i]->getName()
	 << ", Type: " << part[i]->getType()
	 << endl;
    const CsTrack* trk = part[i]->getTrack();
    if( trk != NULL ) {
      CsTrack track = *trk;
      cout << "Track ID: " << track.getId();
      if( dynamic_cast<const CsBeam*>(trk) ) {
	cout << ", This is a BEAM track" << endl;
      }
      else {
	cout << endl;
      }
      vector<CsHelix> hlx = track.getHelices(); 
      for( unsigned int j=0; j<hlx.size(); j++ ) {
	cout << "  Track helix " << j << ", X=("
	     << hlx[j].getX() << ", "
	     << hlx[j].getY() << ", "
	     << hlx[j].getZ() << ")"
	     << ", dx/dz=" << hlx[j].getDXDZ() 
	     << ", dy/dz=" << hlx[j].getDYDZ() 
	     << ", c/p="   << hlx[j].getCop() 
	     << endl;
      } 
      if( track.hasRich1Probs() ) {
	const double* richp = track.getRich1Probs();
	cout << " Track RICH1 probs : ";
	unsigned int jmax = 10;
	if( _dstversion > 3 || _dstversion == 0 ) jmax = _actualRICHDataSize; 
	for( unsigned j=0; j<jmax; j++ ) {
	  cout << *(richp+j) << " ";
	}
	cout << endl;
      }
      if( track.hasMeanTime() ) {
	cout << "Track time: " << track.getMeanTime() 
	     << " +/- " << track.getMeanTimeError() << endl;
      }
      cout << "Track xed rad. len.: " << track.getXX0() << endl;
      cout << "Track Chi2: " << track.getChi2() << endl;

      const unsigned int* exp = track.getExpectedDetsBitmap();
      const unsigned int* fir = track.getFiredDetsBitmap();
      cout << "Expected Detectors BitMap" << endl;
      for( int jj=CSTRACK_MAPSIZE-1; jj>=0; jj-- ) {
	for( int kk=31; kk>=0; kk-- ) {
	  cout << (((exp[jj])>>kk)&1);
	  if( kk==24 || kk==16 || kk==8 ) cout << " ";
	}
	cout << " ";
	if( jj%2 == 0 ) cout << endl;
      }
      cout << "Fired Detectors BitMap" << endl;
      for( int jj=CSTRACK_MAPSIZE-1; jj>=0; jj-- ) {
	for( int kk=31; kk>=0; kk-- ) {
	  cout << (((fir[jj])>>kk)&1);
	  if( kk==24 || kk==16 || kk==8 ) cout << " ";
	}
	cout << " ";
	if( jj%2 == 0 ) cout << endl;
      }	
    }
    vector<Reco::CalorimeterParticle*> calob = part[i]->getCalObjects(); 
    if( ! calob.empty() ) {
      for( unsigned int j=0; j<calob.size(); j++ ) {
	cout << " CalObj: X=("
	     << calob[j]->GetX() << ", "
	     << calob[j]->GetY() << ", "
	     << calob[j]->GetZ() << ")"
	     << ", E=" << calob[j]->GetE() 
	     << endl;
      }
    }
  }

  if( alsoClus ) {
    // dump also clusters 
    list<CsCluster*> clusters = _recoEvent.getClusters();
    list<CsCluster*>::iterator Ic;
    for( Ic=clusters.begin(); Ic!=clusters.end(); Ic++ ) {
      cout << " Cluster: X=("
	   << (*Ic)->getU() << ", "
	   << (*Ic)->getV() << ", "
	   << (*Ic)->getW() << ")";
      double time, analog;
      if( (*Ic)->getTime( time ) ) {
	cout << ", T=" << time;
      }
      if( (*Ic)->getAnalog( analog ) ) {
	cout << ", A=" << analog;
      }
      list<CsDetector*> dets = (*Ic)->getDetsList();
      if( !dets.empty() ) {
	cout << ", Det=" << dets.front()->GetTBName();
      }
      else {
	cout << ", Det=Unknown";
      }
      cout << endl;
    }
  }

  list<CsVertex*> vertices = _recoEvent.getVertices();
  for( list<CsVertex*>::iterator i=vertices.begin(); i!=vertices.end(); i++ ) {
    cout << "Vertex: x=" << (*i)->getX() 
	 << ", y=" << (*i)->getY() 
	 << ", z=" << (*i)->getZ() 
	 << "; Chi2: " << (*i)->getChi2() << endl;
    HepMatrix m(3,3); 
    vector<HepMatrix*> cov = (*i)->getCov();
    for( unsigned int j=0; j<cov.size(); j++ ) {
      m = *(cov[j]);
      cout << " cov[" << j << "]: "
	   << m(1,1) << " " << m(1,2) << " " << m(1,3) << endl
	   << "         " << m(2,1) << " " << m(2,2) << " " << m(2,3) << endl
	   << "         " << m(3,1) << " " << m(3,2) << " " << m(3,3) << endl;
    }
    list<CsTrack*> vtxTracks = (*i)->getTracks();
    for( list<CsTrack*>::iterator vt=vtxTracks.begin();
	 vt!=vtxTracks.end(); vt++ ) {
      Cs3Vector vec;
      (*i)->getPar( (*vt), vec );
            cout << "Track ID: " << (*vt)->getId() 
		 << ", dx/dz: " << vec.dXdZ
		 << ", dy/dz: " << vec.dYdZ
		 << ", c/p: " << vec.Cop
		 << endl;
    }
  }
 
}


bool CsEvent::end() {

#if USE_MySQL
  if(recoSchema_==0) {
    CsPP::Instance()->end();
  }
#endif

    cout << endl 
	 << "+-----------------------------------------------------+" << endl
	 << "|  Number of Processed Events       :  " 
	 << std::setw(10) << std::setfill(' ') << _nSelEvents << "     |" <<endl 
	 << "|  Events skipped (Decoding errors) :  " 
	 << std::setw(10) << std::setfill(' ') << _nSkippedEvents << "     |" <<endl 
	      << "+-----------------------------------------------------+" << endl << endl;

  if( _ntotalTime > 0 ) {

    cout << endl 
	 << "+-----------------------------------------------------+" << endl
	 << "| getNextEvent Mean Elapsed Time (" 
	 << std::setw(10) << std::setfill(' ') << _ntotalTime << " events). |" <<endl 
	 << "+-----------------------------------------------------+" << endl;
    if( _loadCalibTime != 0. ) {
      cout << std::setprecision(4)
	   << std::setfill( ' ' )
	   << std::setiosflags(std::ios::fixed)
	   << "Calibration loading:  "
	   << std::setw(8) << fabs(_loadCalibTime) << " s." << endl
	   << "+-----------------------------------------------------+" << endl;
    }    
    if( _neventLoadTime != 0 ) {
      _eventLoadTime /= _neventLoadTime;
      cout << std::setprecision(4)
	   << std::setfill( ' ' )
	   << std::setiosflags(std::ios::fixed)
	   << "Event loading:        "
	   << std::setw(8) << fabs(_eventLoadTime) << " s." << endl;
    }
    if( _ndigitizTime != 0 ) {
      _digitizTime /= _ndigitizTime;
      cout << std::setprecision(4)
	   << std::setfill( ' ' )
	   << std::setiosflags(std::ios::fixed)
	   << "Digitization:         "
	   << std::setw(8) << fabs(_digitizTime) << " s." << endl;
    }
    if( _nclusterizeTime != 0 ) {
      _clusterizeTime /= _nclusterizeTime;
      cout << std::setprecision(4)
	   << std::setfill( ' ' )
	   << std::setiosflags(std::ios::fixed)
	   << "Clusterization:       " 
	   << std::setw(8) << fabs(_clusterizeTime) << " s." << endl;
    }
    if( _nbeamTime != 0 ) {
      _beamTime /= _nbeamTime; 
      cout << std::setprecision(4)
	   << std::setfill( ' ' )
	   << std::setiosflags(std::ios::fixed)
	   << "Beam Reconstruction:  " 
	   << std::setw(8) << fabs(_beamTime) << " s." << endl;
    }
    if( _nTrackingTime != 0 ) {
      _trackingTime /= _nTrackingTime;
      cout << std::setprecision(4)
	   << std::setfill( ' ' )
	   << std::setiosflags(std::ios::fixed)
	   << "Tracking:             " 
	   << std::setw(8) << fabs(_trackingTime) << " s." << endl;
    }
    
    if( _ncaloTime != 0 ) {
      _caloTime /= _ncaloTime;
      cout << std::setprecision(4)
	   << std::setfill( ' ' )
	   << std::setiosflags(std::ios::fixed)
	   << "Calo Reconstruction:  " 
	   << std::setw(8) << fabs(_caloTime) << " s." << endl;
    }
    if( _nrich1Time != 0 ) {
      _rich1Time /= _nrich1Time;
      cout << std::setprecision(4)
	   << std::setfill( ' ' )
	   << std::setiosflags(std::ios::fixed)
	   << "Rich1 Reconstruction: " 
	   << std::setw(8) << fabs(_rich1Time) << " s." << endl;
    }
    if( _ncedarTime != 0 ) {
      _cedarTime /= _ncedarTime;
      cout << std::setprecision(4)
	   << std::setfill( ' ' )
	   << std::setiosflags(std::ios::fixed)
	   << "CEDAR Reconstruction:  "
	   << std::setw(8) << fabs(_cedarTime) << " s." << endl;
    }

    if( _npartTime != 0 ) {
      _partTime /= _npartTime;
      cout << std::setprecision(4)
	   << std::setfill( ' ' )
	   << std::setiosflags(std::ios::fixed)
	   << "Particle Production:  " 
	   << std::setw(8) << fabs(_partTime) << " s." << endl;
    }
    if( _npidTime != 0 ) {
      _pidTime /= _npidTime;
      cout << std::setprecision(4)
	   << std::setfill( ' ' )
	   << std::setiosflags(std::ios::fixed)
	   << "Particle Id:          " 
	   << std::setw(8) << fabs(_pidTime) << " s." << endl;
    }
    if( _nvertexTime != 0 ) {
      _vertexTime /= _nvertexTime;
      cout << std::setprecision(4)
	   << std::setfill( ' ' )
	   << std::setiosflags(std::ios::fixed)
	   << "Vertexing:            "
	   << std::setw(8) << fabs(_vertexTime) << " s." << endl;
    }
    
    if( _ndstDownTime != 0 ) {
      _dstDownTime /= _ndstDownTime;
      cout << std::setprecision(4)
	   << std::setfill( ' ' )
	   << std::setiosflags(std::ios::fixed)
	   << "DST Downloading:      " 
	   << std::setw(8) << fabs(_dstDownTime) << " s." << endl;
    }
    if( _ndstUpTime != 0 ) {
      _dstUpTime /= _ndstUpTime;
      cout << std::setprecision(4)
	   << std::setfill( ' ' )
	   << std::setiosflags(std::ios::fixed)
	   << "DST Uploading:        " 
	   << std::setw(8) << fabs(_dstUpTime) << " s." << endl;
    }
    if( _nfinishTime != 0 ) {
      _finishTime /= _nfinishTime;
      cout << std::setprecision(4)
	   << std::setfill( ' ' )
	   << std::setiosflags(std::ios::fixed)
	   << "Finish Time:          " 
	   << std::setw(8) << fabs(_finishTime) << " s." << endl;
    }
    if( _ntotalTime != 0 ) {
      _totalTime /= _ntotalTime;
      cout << std::setprecision(4)
	   << std::setfill( ' ' )
	   << std::setiosflags(std::ios::fixed)
	   << "+-----------------------------------------------------+" << endl
	   << "Mean Time Spent on Event: " 
	   << std::setw(8) << fabs(_totalTime) << " s." << endl
	   << "+-----------------------------------------------------+" << endl;
    }
    cout << endl << endl;
  }

  return true;
}

int CsEvent::openRawEventOutputStream( const string name ) {
  
  if( _nOutStreams > 31 ) {
    CsErrLog::mes( elError, "Too many output streams opened: max: 32" );
    return( -1 );
  }

  _nOutStreams++;

  wordexp_t* exp = (wordexp_t*) malloc(1024);
  wordexp( name.c_str(), exp, 0 );
  _outStreams[_nOutStreams-1].open( exp->we_wordv[0] );
  if( _outStreams[_nOutStreams-1].fail() ) {
    string str = "Cannot open ";
    str.append( exp->we_wordv[0] );
    str.append("." ); 
    CsErrLog::mes( elFatal, str );
    wordfree( exp );
    return( -1 );
  }
  else {
    wordfree( exp );
    return( _nOutStreams - 1 );
  }
}

void CsEvent::closeRawEventOutputStream( const int stream ) {
  if( _nOutStreams > stream ) {
    _outStreams[stream].close();
	 _nOutStreams--;
  }
}

void CsEvent::closeAllRawEventOutputStream( void ) {
  for( int i=0; i<_nOutStreams; i++ ) {
    _outStreams[i].close();
	 _nOutStreams = 0;
  }
}

bool CsEvent::outputRawEventToStream( const int stream ) {
  if( _nOutStreams <= stream ) {
    CsErrLog::mes( elError, "Stream number too big." );
    return false;
  }
  else {
    if( ((_outStreamsMask>>stream)&1) == 0 ) { // not already written
      const char* buff = (const char*) _store->rawBuffer();
      int buffsize = *(int*)buff;
      _outStreams[stream].write( buff, buffsize );  // write it...
      _outStreamsMask |= (1<<stream);               // set as written
    }
  }
  return true;
}


bool CsEvent::getPolarization( pair<double,double> &pol ) const
{ 
  
  CsMagInfo* magp = CsGeom::Instance()->getCsField()->getMagInfo();

  if( _DTEvent )
    return magp[0].getPolarization( _store->getRun(), pol );

  return false;
}

bool CsEvent::getHodoDatum( const unsigned int i, string& tbname, 
			    int& channel, double& time ) {
		
  vector<unsigned int> hodoData = _recoEvent.getHodoData();
  if( i >= hodoData.size() ) return false;
  tbname = HodoDataStruct[((hodoData[i]>>27)&0x1f)].tbname;
  channel = ((hodoData[i]>>21)&0x3f);
  time = -float(hodoData[i]&0x1fffff)/100.0;
  return true;

}


bool CsEvent::rebuildDigitsAndHits() {

  bool status = false;

  // Check: accessing DST in read mode?
  if( CsInit::Instance()->getDataType() != "dst" ) {
    CsErrLog::mes( elError, "This function is dummied accessing raw data." );
    return false;
  }
  
  // Check: empty digit list
  if( ! _recoEvent.getDigits().empty() ) {
    CsErrLog::mes( elFatal, "DST with non empty digit list." );
    return false;
  }
    
  // decode raw buffer
  status = _decode();

  if( status ) {
    // Check: empty digit list
    if( ! _recoEvent.getClusters().empty() ) {
      CsErrLog::mes( elInfo, "DST with non empty digit list. Nothing to clusterize." );
    }

    // hit reconstruction
    status = _clusterize();
  }

  return( status );

}

void CsEvent::_readTcsCalibration(time_t timePoint){
	CDB *cdb_;
	std::string filedbloc;
	if (CsOpt::Instance()->getOpt("tcsphase","FileDB",filedbloc)) {
		cerr<<"Warning on "<<"tcsphase"<<": using private FileDB for calib, location: "<<filedbloc<<endl;
		cdb_ = new FileDB(filedbloc);
	} else {
		cdb_ = CsInit::Instance()->getDB();
	}
#if USE_MySQL
	if (CsInit::Instance()->useMySQLDB()) { 
		std::string entrytime;
		if (CsOpt::Instance()->getOpt("tcsphase","CDBentrytime",entrytime) ){
			MySQLDB* mysqldb_ = dynamic_cast<MySQLDB*>(cdb_);
			if (mysqldb_) {
				std::string entrymysqltime = MySQLInterface::toMySQLtime(entrytime.c_str());
				mysqldb_->setSpecificEntryTime("tcsphase", entrymysqltime.c_str());
			} else {
				std::cerr<<"Warning in CsEvent::_readTcsCalibration: CDBentrytime entry for  tcsphase  works only with MySQLDB"<<endl;
			}
		}
	}
#endif

	CDB::Time tp(timePoint,0);
	tm *t = localtime(&tp.first);
	try {
		string s("");
		cdb_->read("tcsphase", s, tp, "_jumps");
		if (s == "") throw CS::Exception("empty string from calib file");
		std::istringstream is(s);
		std::string str;
		while( getline(is,str) ) {
			unsigned int run,spill;
			float jump;
			if(str=="" || str.c_str()[0]=='#') continue;
			int ret = sscanf(str.c_str(),"%d %d %f",&run,&spill,&jump);
			assert( ret == 3 );
			if( run == getRunNumber()){
				if(_tcscorr.count(spill)  == 0){
					_tcscorr.insert( std::pair<int,double>(spill,jump) );
				}
			}
		}

	}
	catch( const std::exception &e ) { 
		std::cerr << "TCSphase jump correction, no calibration for local time ";
		std::cerr << t << ": " << e.what() << endl;
	}


}

