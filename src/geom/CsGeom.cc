// $Id: CsGeom.cc 14094 2015-11-06 15:28:48Z lsilva $

/*!
   \file    CsGeom.cc
   \brief   Compass Geometry Interface Class.
   \author  Benigno Gobbo
   \version $Revision: 14094 $
   \date    $Date: 2015-11-06 16:28:48 +0100 (Fri, 06 Nov 2015) $
*/

#include <cmath>
#include "CsInit.h"
#include "CsGeom.h"
#include "CsErrLog.h"
#include "CsOpt.h"
#include <cstring>
#include "CsStrawTubesDetector.h"
#include "CsDriftChamberDetector.h"
#include "CsMWPCDetector.h"
#include "CsGEMDetector.h"
#include "CsPixelGEMDetector.h"
#include "CsMicroMegaDetector.h"
#include "CsPixelMumegaDetector.h"
#include "CsSiTrackerDetector.h"
#include "CsDFiberHodoDetector.h"
#include "CsJFiberHodoDetector.h"
#include "CsTriggerHodoDetector.h"
#include "CsBMSDetector.h"
#include "CsRICH1Detector.h"
#include "CsRICH1UpGrade.h"
#include "CsDriftTubeDetector.h"
#include "CsRichWallDetector.h"
#include "CsDWDetector.h"
#include "CsRPDetector.h"
#include "CsCalorimeter.h"
#include "CsHCAL1.h"
#include "CsHCAL2.h"

#include "CsECAL0.h"
#include "CsECAL1.h"
#include "CsECAL2.h"

#include "CsMW1Detector.h"
#include <functional>
#include <cstdlib>
#include "CsRCDetectors.h"
#include "CsRCMirrors.h"
#include "CLHEP/Matrix/Vector.h"
#include "CsStopwatch.h"
#include <sstream>
#include <map>
//LS Eff
//#define DEBUG_EFF


// ------------To check resource usages, Benigno 20020226---------------
#define CHECK_RESOURCES 0

#if CHECK_RESOURCES

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>

/*
  it's written in pure C to be used elsewhere...
*/

void usage( char* str ) {

  struct rusage usage;
  double utime, stime;
  char filename[80];
  FILE *fid;
  long size=0, resident=0, share=0, trs=0, drs=0, lrs=0, drt=0;
  char *buffer, *cursor;
  int status, bufsize;
  long vmsize=0, vmlck=0, vmrss=0, vmdata=0, vmstk=0, vmexe=0, vmlib=0;

  /* A mess: on Linux no memory info available from getrusage... */

  if( getrusage( RUSAGE_SELF, &usage ) == 0 ) {
    utime = usage.ru_utime.tv_sec + usage.ru_utime.tv_usec/1000000.;
    stime = usage.ru_stime.tv_sec + usage.ru_stime.tv_usec/1000000.;
    printf( "+----------------------------------------------------+\n" );
    printf( "%s\n", str );
    printf( "+----------------------------------------------------+\n" );
    printf( "TIME (user):   %11.2f s\n", utime );
    printf( "TIME (sys):    %11.2f s\n", stime );
    printf( "SWAP num.:     %8d\n", usage.ru_nswap );
  }
  else {
    printf( "getrusage failed.\n" );
  }

  sprintf( filename, "/proc/%d/statm", getpid() );

  fid = fopen( filename, "r" );

  if( fid != NULL ) {
    fscanf( fid, "%d %d %d %d %d %d %d", &size, &resident,
      &share, &trs, &drs, &lrs, &drt );
    fclose( fid );
    printf( "STM. total:      %6d pages.\n", size );
    printf( "STM. resident:   %6d pages.\n", resident );
    printf( "STM. shared:     %6d pages.\n", share );
    printf( "STM. code:       %6d pages.\n", trs );
    printf( "STM. data/stack: %6d pages.\n", drs );
    printf( "STM. library:    %6d pages.\n", lrs );
    printf( "STM. dirty:      %6d pages.\n", drt );
  }
  else {
    printf( "%s opening failed\n", filename );
  }

  sprintf( filename, "/proc/%d/status", getpid() );

  fid = fopen( filename, "r" );

  if( fid != NULL ) {
/*  this does not work, why???
    status = fseek( fid, 0, SEEK_END );
    size = ftell( fid );
    rewind( fid );
*/
    bufsize = 1024;
    buffer = (char*) malloc( bufsize );
    status = fread( buffer, 1, bufsize, fid );
    fclose( fid );
    if( status > 0 ) {
      cursor = strstr( buffer, "VmSize:" );
      if( cursor ) {
	sscanf( cursor,
		"VmSize: %lu kB\n"
		"VmLck: %lu kB\n"
		"VmRSS: %lu kB\n"
		"VmData: %lu kB\n"
		"VmStk: %lu kB\n"
		"VmExe: %lu kB\n"
		"VmLib: %lu kB\n",
		&vmsize, &vmlck, &vmrss, &vmdata,
		&vmstk, &vmexe, &vmlib );
	printf( "MEM. total:      %6d kB.\n", vmsize );
	printf( "MEM. locked:     %6d kB.\n", vmlck );
	printf( "MEM. resident:   %6d kB.\n", vmrss );
	printf( "MEM. data:       %6d kB.\n", vmdata );
	printf( "MEM. stack:      %6d kB.\n", vmstk );
	printf( "MEM. execut.:    %6d kB.\n", vmexe );
	printf( "MEM. library:    %6d kB.\n", vmlib );
      }
    }
    free( buffer );
  }
  else {
    printf( "%s opening failed\n", filename );
  }

  printf( "+----------------------------------------------------+\n" );

}
#endif // CHECK_RESOURCES
// ---------------------------------------------------------------------

using namespace std;
using namespace CLHEP;

CsGeom* CsGeom::_instance = 0;

CsGeom* CsGeom::Instance( const string* detTableName ) {
  if( _instance == 0 ) {
    _instance = new CsGeom( detTableName );
  }
  return( _instance );
}

CsGeom* CsGeom::Instance( const string* detTableName, const string* pitchTableName ) {
  if( _instance == 0 ) {
    _instance = new CsGeom( detTableName, pitchTableName );
  }
  return( _instance );
}

CsGeom* CsGeom::Instance() {
  if( _instance != 0 )
    return( _instance );
  else {
    CsErrLog::mes( elFatal, "wrong CsGeom singleton instance." );
    return( NULL );
  }
}


CsGeom::CsGeom( const string* detTableName ) {

  CsStopwatch chronos;
  int chrono = chronos.start();

  // clearage...
  _matrices.clear();
  _dets.clear();
  _others.clear();
  _rich1 = NULL;
  _zones.clear();
  _muonSetup = false;
  _hadronSetup = false;
  _targetCenter = 0;     // Default value for _targetCenter

  if( !readDetTable( detTableName ) ) {
    string str = "Geant detector table file ";
    str.append( *(detTableName) );
    str.append(" not found." );
    CsErrLog::Instance()->mes( elFatal, str );
  }

  if( !makeZones() ) {
    string str = "Problems in making zones ";
    CsErrLog::Instance()->mes( elFatal, str );
  }

//=== This as been moved to CsEvent.cc, so that histograms are booked after
//=== Calibration DB is read.
//   // book histograms for all CsDetectors
//   cout << "CsGeom::CsGeom - INFO: booking histograms" << endl;
//   for(list<CsDetector*>::iterator idet=_dets.begin(); idet!=_dets.end();idet++)
//   (*idet)->BookHistograms();
//   cout << "CsGeom::CsGeom - INFO: end of histogram booking" << endl;

  setDetsZones();

  //        ********** DETECTOR ASSOCIATION **********
  CsOpt* opt = CsOpt::Instance();

  //       ***** DRIFT-LIKE DETECTOR ASSOCIATION *****
  string tag = "LR"; string key = "make associations";
  if (!opt->getOpt(tag,key) && !opt->getOpt(tag,"detector association"))
    CsErrLog::mes(elInfo,
      "Detector association for LR ambiguity raising not requested");
  else if (!associateDets())
    CsErrLog::mes(elInfo,
      "Detector association for LR ambiguity raising cancelled.");

  //         ***** CsGEM AND CsPixelGEM DETECTORS ASSOCIATION *****
  // The "tag" used here is "GEM" and not "GM" or "GP", although the option
  // "make associations" typically appears next to other options, e.g.
  // "ampCorrelationMC" or "Master", introduced by "GM" or "GP", for:
  //  - It then can be made to apply to both "GM" and "GP" GEMs.
  //  - The use of the 1st two characters of a TB name indicates elswhere in
  //   coral that the option can be enabled for particular instances, as opposed
  //   to the whole class, by specifying their full TB name instead. Which
  //   is not the case here.
  int GEMAssociation = 0;
  tag = "GEM"; if (opt->getOpt(tag,key)) {
    CsErrLog::mes(elInfo,"GEM detectors association.");     GEMAssociation |= 1;
  }
  else if (opt->getOpt("GM","ampCorrelationMC")) CsErrLog::mes(elFatal,
 "\"GM\" \"ampCorrelationMC\" requested while \"GEM make associations\" isn't");
  else if (opt->getOpt("GP","ampCorrelationMC")) CsErrLog::mes(elFatal,
 "\"GP\" \"ampCorrelationMC\" requested while \"GEM make associations\" isn't");
  //         ***** CsPixelGEM PIECES PIXELISED<->STRIPPED ASSOCIATION *****
  // - This association is used to get a single COMGeant ID, encompassing both
  //  the pixelised central piece of a pixelGEM and one of its 2 orthogonal
  //  stripped external pieces, feed MC hits to both kind of pieces. It has had
  //  to be introduced in order to provide for COMGeant setup version 2007.05_1.
  //  May not be usefull in the long term.
  tag = "pixelGEM"; if (opt->getOpt(tag,key)) {
    CsErrLog::mes(elInfo,"CsPixelGEM->CsGEM association."); GEMAssociation |= 2;
  }
  if (GEMAssociation) associateGEMs(GEMAssociation);

  cout << endl << endl
       << "---------------------------------------------------" << endl
       << " Geometry related objects setup ended in "
       << chronos.stop( chrono )
       << " s." << endl
       << "---------------------------------------------------" << endl
       << endl;

}


CsGeom::CsGeom( const string* detTableName, const string* pitchTableName ) {

  CsStopwatch chronos;
  int chrono = chronos.start();

  // clearage...
  _matrices.clear();
  _dets.clear();
  _others.clear();
  _rich1 = NULL;
  _zones.clear();
  _muonSetup = false;
  _hadronSetup = false;
  _targetCenter = 0;     // Default value for _targetCenter

  if( !readDetTable( detTableName ) ) {
    string str = "Geant detector table file ";
    str.append( *(detTableName) );
    str.append(" not found." );
    CsErrLog::mes( elFatal, str );
  }

//   string pitchfile = "/afs/cern.ch/user/a/alehmann/w0/Coral-28Dec05/coral/src/hardcoral/PitchTable.dat";
//   const string* pitchTableName = &pitchfile;
  if( !readPitchTable( pitchTableName ) ) {
    string str = "Pitch table file ";
    str.append( *(pitchTableName) );
    str.append(" not found." );
    CsErrLog::mes( elFatal, str );
  }

  if (!makeZones())    //  ********** TRACK RECONSTRUCTION ZONES **********
    CsErrLog::mes(elFatal,"Problems in making zones");
  setDetsZones();

  //        ********** DETECTOR ASSOCIATION **********
  CsOpt* opt = CsOpt::Instance();

  //       ***** DRIFT-LIKE DETECTOR ASSOCIATION *****
  string tag = "LR"; string key = "make associations";
  if (!opt->getOpt(tag,key) && !opt->getOpt(tag,"detector association"))
    CsErrLog::mes(elInfo,
      "Detector association for LR ambiguity raising not requested");
  else if (!associateDets())
    CsErrLog::mes(elInfo,
      "Detector association for LR ambiguity raising cancelled.");

  //         ***** CsGEM AND CsPixelGEM DETECTORS ASSOCIATION *****
  // The "tag" used here is "GEM" and not "GM" or "GP", although the option
  // "make associations" typically appears next to other options, e.g.
  // "ampCorrelationMC" or "Master", introduced by "GM" or "GP", for:
  //  - It then can be made to apply to both "GM" and "GP" GEMs.
  //  - The use of the 1st two characters of a TB name indicates elswhere in
  //   coral that the option can be enabled for particular instances, as opposed
  //   to the whole class, by specifying their full TB name instead. Which
  //   is not the case here.
  int GEMAssociation = 0;
  tag = "GEM"; if (opt->getOpt(tag,key)) {
    CsErrLog::mes(elInfo,"GEM detectors association.");     GEMAssociation |= 1;
  }
  else if (opt->getOpt("GM","ampCorrelationMC")) CsErrLog::mes(elFatal,
 "\"GM\" \"ampCorrelationMC\" requested while \"GEM make associations\" isn't");
  else if (opt->getOpt("GP","ampCorrelationMC")) CsErrLog::mes(elFatal,
 "\"GP\" \"ampCorrelationMC\" requested while \"GEM make associations\" isn't");
  //         ***** CsPixelGEM PIECES PIXELISED<->STRIPPED ASSOCIATION *****
  // - This association is used to get a single COMGeant ID, encompassing both
  //  the pixelised central piece of a pixelGEM and one of its 2 orthogonal
  //  stripped external pieces, feed MC hits to both kind of pieces. It has had
  //  to be introduced in order to provide for COMGeant setup version 2007.05_1.
  //  May not be usefull in the long term.
  tag = "pixelGEM"; if (opt->getOpt(tag,key)) {
    CsErrLog::mes(elInfo,"CsPixelGEM->CsGEM association."); GEMAssociation |= 2;
  }
  if (GEMAssociation) associateGEMs(GEMAssociation);


  cout << endl << endl
       << "---------------------------------------------------" << endl
       << " Geometry related objects setup ended in "
       << chronos.stop( chrono )
       << " s." << endl
       << "---------------------------------------------------" << endl
       << endl;

}


struct CsGeom::sortDetectors_ :
  public binary_function<CsDetector*, CsDetector*, bool> {
  bool operator() ( CsDetector* d1, CsDetector* d2 ) {
    double z1 = d1->getZcm(), z2 = d2->getZcm();

    if (z1==z2) {
      // Special cases: pixel GEM and pixel MM
      const string &name1 = d1->GetTBName(), &name2 = d2->GetTBName();
      if ((name1.find("GP")==0 && name2.find("GP")==0) ||
	  (name1.find("GM")==0 && name2.find("GM")==0)) {
	// [Pixel]GEMs: Enforcing ordering (P1<U<V < Y<X<P2) needed by trafdic
	char c1 = name1[4], c2 = name2[4];
	string corrSeq("YXPUV");
	short t1 = corrSeq.find(c1), t2 = corrSeq.find(c2);
	return t1<=t2 ;
      }
      if (name1.find("MP")==0 && name2.find("MP")==0) {
	// PixelMMs: Enforcing ordering (M1<X|V or Y|U<M2) needed by trafdic
	char c1 = name1[4]; if (c1=='P') c1='M';// Older MPs: 'P' instead of 'M'
	char c2 = name2[4]; if (c2=='P') c2='M';
	char *cXYUV = c1=='M' ? &c2 : (c2=='M' ? &c1 : 0);
	if (cXYUV) {
	  if      (*cXYUV=='X' || *cXYUV=='V') return cXYUV==&c2;
	  else if (*cXYUV=='Y' || *cXYUV=='U') return cXYUV==&c1;
	  else cXYUV = 0;
	}
	if (!cXYUV)
	  CsErrLog::msg(elFatal,__FILE__,__LINE__,
			"PixelMMs \"%s\",\"%s\" have same abscissa(=%.5f) while"
			" they're not a (X|Y|U|V,M|P) pair of coordinates",
			name1.c_str(),name2.c_str(),z1);
      }
      // Other cases like ST or MB ('a','b','c' or 'r','l') slices
      double x1 = d1->getXcm(), x2 = d2->getXcm();
      if (x1==x2) {
	double y1 = d1->getYcm(), y2 = d2->getYcm();
	return y1<y2;
      }
      return x1<x2;
    }
    else return z1<z2;
  }
};

bool CsGeom::makeZones() {
  /*
  \brief   Build the list of CsZones = track reconstruction zones

  - From options "define zone", w/ arguments: Zmin, Zmax, Name.
  - Ordering is according to the sucession of option entries. 
  - By default: single zone "All Apparatus".
  - Option is recursive. Yet, re-specifying a zone w/ same name but different
  range, overwrites previous range.
  */

  list<CsDetector*> dets = _dets;
  CsOpt* opt = CsOpt::Instance();
  bool   found = false;

  vector<float> z;

  string tag = "define";
  string key = "zone";
  list<string> limits;
  while( opt->getOptRec( tag, key, limits ) ) {
    found = true;
    float zmin, zmax;
    string name;
    list<string>::iterator Is;
    int i = 0;
    for( Is=limits.begin(); Is!=limits.end(); Is++, i++ ) {
      if( i == 0 ) {
        istringstream( *Is ) >> zmin;
      }
      else if( i == 1 ) {
        istringstream( *Is ) >> zmax;
      }
      else if( i == 2 ) {
	name = (*Is);
      }
      else {
	name += " ";
	name += (*Is);
      }
    }
    dets.clear();
    list<CsDetector*>::iterator idet;
    for( idet=_dets.begin(); idet!=_dets.end(); idet++ ) {
      double z = (*idet)->getZcm();
      if( zmin <= z && z < zmax ) {
	dets.push_back( (*idet) );
      }
    }

    CsZone *zone = new CsZone(zmin,zmax,name,dets);

    // ***** CURRENT ZONE IS A REDEFINITION of A PREVIOUSLY ENTERED ONE?.. *****
    list<CsZone*>::iterator iZ; bool overWrite;
    for (iZ = _zones.begin(), overWrite = false; iZ!=_zones.end(); iZ++) {
      const CsZone *zi = *iZ; if (zi->getName()==name) {
	CsErrLog::msg(elBasicInfo,__FILE__,__LINE__,
		      "Redefining zone \"%s\": [%.2f,%.2f] -> [%.2f,%.2f]",
		      name.c_str(),zi->getZMin(),zi->getZMax(),zmin,zmax);
	overWrite = true;            // ... => OVERWRITE PREVIOUS ONE
	// Note that this still allows to define overlapping zones: suffices to
	// declare them w/ different names.
	_zones.insert(iZ,zone); _zones.erase(iZ); break;
      }
    }
    if (!overWrite) { //   ***** ...ELSE APPEND NEW CsZone *****
      _zones.push_back(zone);
    }
  }

  if (!found) {
    CsZone* zone = new CsZone( 0., 99999., "All Apparatus", dets );
    _zones.push_back( zone );
  }

  return true;
}

void CsGeom::setDetsZones() {

  list<CsDetector*>::iterator Id;
  list<CsZone*>::iterator Iz;
  for( Id=_dets.begin(); Id!=_dets.end(); Id++ ) {
    float z = (*Id)->getZcm();
    for( Iz=_zones.begin(); Iz!=_zones.end(); Iz++ ) {
      float z1 = (*Iz)->getZMin();
      float z2 = (*Iz)->getZMax();
      if( z1<z && z<z2 ) {
	(*Id)->addZone( *(*Iz) );
      }
    }
  }
}


bool CsGeom::readPitchTable( const string* pitchTableName ) {
  // read the table containing the variable pitch corrections

  const int lineSize = 16384;        //  needed for detectors with 1024 wires
  char   line[lineSize];
  string opt;
  double pcorr;
  bool failure = false;

  // open input file
  int len = pitchTableName->size();
  ifstream f( pitchTableName->c_str(), ios::in );
  if( !f.good() ) {
    _pitchTableName.erase();
    return( false );
  }


  CsOpt* options = CsOpt::Instance();
  list<string> _noVarPitch;
  if (options->getOpt("VarPitch","DetOffVarP")) options->getOpt( "VarPitch","DetOffVarP",_noVarPitch);
  typedef list<string>::iterator NoPI;
  for(NoPI nopi=_noVarPitch.begin(); nopi!=_noVarPitch.end(); nopi++) {
    cout << "Detectors " << (*nopi) << " are excluded from the variabled-sized pitch list!" << endl;
  }


//   _PitchCorr.clear();

  while ( f.good() ) {
    f.getline( line, lineSize, '\n' );
    if( line[0] == '\0' || line[0] != ' ' ) continue;

    // something good in this line...
    istringstream s(line);
    s >> opt; if( ! s ) continue;

    char TBname[9]="";

    _PitchCorr.clear();

    // detector?
    if( opt == "det" ) {
      s >> TBname;       // det. TB name
      int  Nwire;         s >> Nwire;        // number of wires
      for ( int i = 0; i < Nwire; i++ ) {
        s >> pcorr;
        _PitchCorr[i] = pcorr;
      }
    }
    // Store eventual failures
    failure = s.fail();
    s.clear();


    // bypass detectors excluded from variable-sized pitch list
    string detname = TBname;
//     typedef list<string>::iterator NoPI;
    for(NoPI nopi=_noVarPitch.begin(); nopi!=_noVarPitch.end(); nopi++) {
      if ( detname.find(*nopi) != string::npos ) goto next_line;
    }

    // setup the variable-sized pitches
    typedef list<CsDetector*>::iterator DI;
    for(DI di=_dets.begin(); di!=_dets.end(); di++) {
      if((*di)->GetTBName()==TBname) {
	(*di)->SetVarPitch(&_PitchCorr);
	cout << TBname << " has a variable-sized pitch" << endl;
      }
    }

//     cout << TBname;
//     for (std::map<int,double>::iterator ipitch=_PitchCorr.begin(); ipitch!=_PitchCorr.end(); ipitch++ ) {
//       cout << "  " << ipitch->second;
//     }
//     cout << endl;

  next_line:;

  }

  // Check if f readout ended "somehow" correctly:
  if( failure ) {
    CsErrLog::mes( elFatal, "It seems there's something wrong in Pitchtable.dat file, please have a look to it." );
    exit(1);
  }

  return( true );
}


bool CsGeom::readDetTable( const string* detTableName ) {

# if CHECK_RESOURCES
  //+---------- Check Resource Usages, Benigno 20020226
  usage( "Entering CsGeom::readDetTable" );
# endif

  //  const int lineSize = 256;
  const int lineSize = 512;        //  needed for RICH 1 lines
  char   line[lineSize];
  string opt;
  int i, j;
  int magCount = 0;
  int ndets = 0;

  int CGEAversion = 0;
  int CGEArelease = 0;
  int CGGeomVersion = 0;
  int CGGeomRelease = 0;

  bool failure = false;

  // open input file
  int len = detTableName->size();
  ifstream f( detTableName->c_str(), ios::in );
  if ( !f.good() ) {
    _detTableName.erase();
    return( false );
  }

  while ( f.good() ) {
    f.getline( line, lineSize, '\n' );
    if( line[0] == '\0' || line[0] != ' ' ) continue;
    if( line[1] == ' ' || line[1] == '-' || line[1] == '=' ) continue;

    // something good in this line...
    istringstream s(line);
    s >> opt; if( ! s ) continue;

    // Comgeant version?
    if( opt == "ver" ) {
      s >> _comgeantVersion;
      CGEAversion = atoi( _comgeantVersion.c_str() + 1 );
      CGEArelease = atoi( _comgeantVersion.c_str() + 7 );
      ostringstream ostr;
      ostr << "COMGEANT version " << CGEAversion
	   << ", release " << CGEArelease;
      CsErrLog::mes( elInfo, ostr.str() );
    }
    else if( opt == "mgeo" ) {
      _muonSetup = true;
      s >> _geomVersion;
      CGGeomVersion = atoi( _geomVersion.c_str() + 1 );
      CGGeomRelease = atoi( _geomVersion.c_str() + 7 );
      ostringstream ostr;
      ostr << "COMGEANT Muon Setup, geometry version " << CGEAversion
	   << ", release " << CGEArelease;
      CsErrLog::mes( elInfo, ostr.str() );
    }
    else if( opt == "hgeo" ) {
      _hadronSetup = true;
      s >> _geomVersion;
      CGGeomVersion = atoi( _geomVersion.c_str() + 1 );
      CGGeomRelease = atoi( _geomVersion.c_str() + 7 );
      ostringstream ostr;
      ostr << "COMGEANT Hadron Setup, geometry version " << CGEAversion
	   << ", release " << CGEArelease;
      CsErrLog::mes( elInfo, ostr.str() );
    }

    // orientation?
    else if( opt == "ori" ) {
      for( i=0; i<3; i++ ) s >> _orient[i];  // not used at present...
    }

    // magnetic field?
    else if( opt == "mag" ) {  // magnet parameters
      magCount++;
      CsMagInfo magp;
      s >> magp.mag;   // magnet
      s >> magp.xcm; magp.xcm *= 10.;  // centre
      s >> magp.ycm; magp.ycm *= 10.;
      s >> magp.zcm; magp.zcm *= 10.;
      s >> magp.rot;   // rotation
      s >> magp.fsc;   // field scale
      s >> magp.flg1;  // flags
      s >> magp.flg2;
      if( magCount == 3 ) s >> magp.curr;  // current. Applies only for SM2
      else magp.curr = 0;
      CsGeom::GeaRef2CsRefVec( magp.xcm,
				 magp.ycm,
				 magp.zcm );

      if( magp.mag == 3 ) {       // If solenoid magnet...
	if (CGEAversion<7 || ( CGEAversion==7 && CGEArelease<3 ) ) {
	  //=== Early releases: no target specific info (cf. "targ" keyword)
	  // => Equate target center position w/ target magnet's
	  _targetCenter = magp.zcm; // (in mm)
	  CsErrLog::msg(elInfo,__FILE__,__LINE__,
	    "Target default assignment = target magnet's center = %f",
				    _targetCenter);
	  // N.B.: - In the hadron setup, the target magnet does not exist as a
	  // physical entity, it's defined still and used (in the early geometry
	  // releases) to indicate the target position.
	  //       - The detailed distribution of material in the target is
	  //        specified elsewhere and independently.
	}
      }

      _field.addMagInfo( magp );
    }

    else if (opt=="targ") {     // ***** TARGET DESCRIPTION *****
      // - Rudimentary. A detailed description can (and should) be supplied via
      //  material map.
      // - The description given in the "detectors.dat" allows a variety of
      //  orientations and shapes: here we assume that it's a cylinder aligned
      //  w/ beam axis.
      CsTargetCell cell; int iCell, iData; string name; double xData;
      s >> iCell; if (iCell!=(int)_targetCells.size()+1)
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
	  "Error reading target info: %dth cell has # = %d",
				  _targetCells.size()+1,iCell);
      s >> name; s >> iData; // Info disregarded
      s >> iData; if (iData!=5)
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
	  "Error reading target info: shape = %d != 5",iData);
      s >> xData;            // Info disregarded.
      s >> cell.radius; cell.radius *= 5; // "detectors.dat" gives the diameter
      s >> cell.length; cell.length *= 10;
      s >> cell.x[0]; cell.x[0] *= 10.;
      s >> cell.x[1]; cell.x[1] *= 10.;
      s >> cell.x[2]; cell.x[2] *= 10.;
      _targetCells.push_back(cell);
      // Derive "_targetCenter" from set of cells
      double sWeights; for (iCell = 0, _targetCenter=sWeights = 0;
			    iCell<(int)_targetCells.size(); iCell++) {
	cell = _targetCells[iCell]; double w = cell.length;
	_targetCenter += cell.x[0]*w; sWeights += w;
      }
      _targetCenter /= sWeights;
    }

    else if( opt == "rot" ) {   // rotation matrices
      HepMatrix rotM(3,3);
      int dummy; s >> dummy;
      for( i=0; i<3; i++ ) {
	for( j=0; j<3; j++ ) {
          s >> rotM(i+1,j+1);
	}
      }
      CsGeom::GeaRef2CsRefMat( rotM );
      _matrices.push_back( rotM );
    }
    // beam back propagation 
    else if( opt == "bbp" ) {

      //read line
      char     bbp_plane_name[9]="";   s >> bbp_plane_name;
      double   bbp_cor;                s >> bbp_cor;
      double   bbp_res;                s >> bbp_res;

      //check if resolution not negative
      if( bbp_res <= 0. ){
         CsErrLog::msg(elFatal,__FILE__,__LINE__,
	    "Forbidden value of beam back propagation coefficient, resolution = %f for %s", bbp_res, bbp_plane_name);
      }

      //create object
      CsBeamBackPropagationCoeff bbp(bbp_plane_name, bbp_cor, bbp_res);

      //add
      _bbp.push_back( bbp);

    }
    // detector?
    else if( opt == "det" ) {

      ndets++;

      int    id;      s >> id;          // detector number
      char TBname[9]=""; s>>TBname; // det. TB name
      char   name[5]; s >> name;        // detector name
      int    unit;    s >> unit;        // detector number in station
      int    type;    s >> type;        // detector type
      double rdLen;   s >> rdLen;       // radiation length
      double xsiz;    s >> xsiz; xsiz *= 10.; // detector size (mm)
      double ysiz;    s >> ysiz; ysiz *= 10.;
      double zsiz;    s >> zsiz; zsiz *= 10.;
      CsGeom::GeaRef2CsRefVec( xsiz, ysiz, zsiz );
      double xcm;     s >> xcm; xcm *= 10.; // detector centre (MRS) (mm)
      double ycm;     s >> ycm; ycm *= 10.;
      double zcm;     s >> zcm; zcm *= 10.;
      CsGeom::GeaRef2CsRefVec( xcm, ycm, zcm );
      int    rotMNum; s >> rotMNum;     // rotation matrix number
      double wirD;    s >> wirD; wirD *= 10.;  // 1st wire offset (mm)
      double ang;     s >> ang;         // angle of wires in DRS (degrees)
      int    nWir;    s >> nWir;        // number of wires
      double wirP;    s >> wirP; wirP *= 10.;  // wires pitch
      double eff;     // detector efficiency
      //LS Eff
      string sTBname = TBname, sTBn = sTBname.substr(0,2), key = "mcEffMapsEnabled"; 
      if (CsOpt::Instance()->getOpt(sTBn,key)){
	  //CsOpt::Instance()->getOpt("ALL",key)) { ... working on a smart tag
        eff = 1.0; //if Eff Maps enabled then set TBname detector global eff to 1. 
#ifdef DEBUG_EFF
        cout<<"**** LS *** : CsGeom::readDetTable : "<< sTBn <<" global Eff set to "<< eff <<" due to "<<key<<endl;
#endif 
      }
      else      s >> eff;         // detector efficiency
      double bkg;     s >> bkg;         // detector background
      double tGate;   s >> tGate;       // detector time gate

      // Put needed quantities in the WRS (wire reference system):
      // X' orthogonal to the wires direction
      // Y' parallel to the wires direction
      // Z' == Z MRS
      // Origin == Origin of the MRS

      // Rotation angle around Z of the DRS w.r.t. the Y axis of the MRS:
      double angDRS = atan2(_matrices[rotMNum-1](2,1),_matrices[rotMNum-1](1,1));
      ang = ang * (M_PI) / 180.;   // angle in radians...

      // Rotation matrice of wires w.r.t. DRS
      HepMatrix rotWire(3,3,1);
      rotWire(1,1) = cos(ang);
      rotWire(1,2) = -sin(ang);
      rotWire(2,1) = sin(ang);
      rotWire(2,2) = cos(ang);

      // Rotation matrice of WRS w.r.t. MRS:
      HepMatrix rotWRS(3,3,1);
      rotWRS = _matrices[rotMNum-1] * rotWire;

      // Calculate the angle of wires w.r.t. the Y axis of the MRS:
      double angWRS = atan2( rotWRS(2,1), rotWRS(1,1) );

      // Recalculate the first wire offset w.r.t. the WRS:
      double wirDWRS = xcm * cos(angWRS) + ycm * sin(angWRS) + wirD;

      angWRS = angWRS * 180. / (M_PI); // angle in degrees...

      // Detectors with drift measurement.
      double vel  =0.;        // velocity in mm/time_slices
      double t0   =0.;        // T0 in ns
      double thRes=0.;        // 2 hits resolution (time slices)
      double spSli=0.;        // space resolution (time slices)
      double tiSli=0.;        // time slices in ns

      // Store eventual failures
      failure = s.fail();

      if( s >> vel ) {
	vel *= 10.; // Vel drom COMGEANT is in cm/s
	s >> t0;
	s >> thRes;
	s >> spSli;
	s >> tiSli;
	if( s.fail() ) {
	  tiSli = 1.;     // if not available assume tsSli = 1 ns
          ostringstream ostr;
	  ostr << "Time Slice size information not available for "
	       << name << " detector. I set it to "
	       << tiSli
	       << " ns.";
	  CsErrLog::mes( elAnomaly, ostr.str() );
	}
      }
      s.clear();

      // Detector TDC resolution (if available). By hand ...
      double tRes = 0.;

      if( CGEAversion>6 || ( CGEAversion==6 && CGEArelease>2 ) ) {
	tRes = spSli;
      }
      else {
	// Here we still use names as TBnames were not available on those
	// Comgeant releases...
	// Sci. Fi.
	if( strncmp( name, "VSA", 3 ) == 0 ) tRes =  0.35;
	if( strncmp( name, "VSB", 3 ) == 0 ) tRes =  0.35;
	if( strncmp( name, "VSD", 3 ) == 0 ) tRes =  0.35;
	// Hodoscopes (no exact info yet...), use this (pessimistic) value
	if( strncmp( name, "H4V", 3 ) == 0 ) tRes =  0.2;
	if( strncmp( name, "H5V", 3 ) == 0 ) tRes =  0.2;
	if( strncmp( name, "HD4", 3 ) == 0 ) tRes =  0.2;
	if( strncmp( name, "HD5", 3 ) == 0 ) tRes =  0.2;
      }

      // Fill The correct detector type
      bool exist=false;

      //  Check whether any detector w/ same TBname as current "det" entry
      // has already been instantiated.
      //  If so, complete its construction w/ the info from current "det" entry.
      // It can be an adding an extra slice of wires to a variable picth
      // detector. Or an extra dimension to a pixel GEM.
      typedef list<CsDetector*>::iterator DI;
      for(DI di=_dets.begin(); di!=_dets.end(); di++) {
	if((*di)->GetTBName()==TBname) {
	  (*di)->AddSubDetector(ndets, id, name, TBname,unit,  type,
				rdLen, xsiz,  ysiz, zsiz,
				xcm,   ycm,   zcm,  _matrices[rotMNum-1],
				rotWRS, wirDWRS,  angWRS,   nWir, wirP,
				eff,   bkg, tGate) ;
	  exist=true;
	}
      }

      if(exist) continue;

      // Drift Chamber like detectors
      if( strncmp( TBname, "DC", 2 ) == 0 ) {  // Drift Chamber

#       if CHECK_RESOURCES
	//+---------- Check Resource Usages, Benigno 20020226
	usage( "Drift Chamber" );
#       endif

	CsDriftChamberDetector* det =
	  new CsDriftChamberDetector( ndets, id,    name, TBname,unit,  type,
				      rdLen, xsiz,  ysiz, zsiz,
				      xcm,   ycm,   zcm,  _matrices[rotMNum-1],
				      rotWRS, wirDWRS,  angWRS,   nWir, wirP,
				      eff,   bkg, tGate,
				      vel,  t0, thRes, spSli, tiSli );
	_dets.push_back( det );

#       if CHECK_RESOURCES
	//+---------- Check Resource Usages, Benigno 20020226
	usage( "End Drift Chamber" );
#       endif

      }

      if( strncmp( TBname, "DW", 2 ) == 0 ) { // DW = w45

#       if CHECK_RESOURCES
	//+---------- Check Resource Usages, Benigno 20020226
	usage( "DW" );
#       endif

	CsDWDetector* det =
	  new CsDWDetector( ndets, id,    name, TBname,unit,  type,
				      rdLen, xsiz,  ysiz, zsiz,
				      xcm,   ycm,   zcm,  _matrices[rotMNum-1],
				      rotWRS, wirDWRS,  angWRS,   nWir, wirP,
				      eff,   bkg, tGate,
				      vel,  t0, thRes, spSli, tiSli );
	_dets.push_back( det );

#       if CHECK_RESOURCES
	//+---------- Check Resource Usages, Benigno 20020226
	usage( "End DW" );
#       endif

      }


      if( strncmp( TBname, "MB", 2 ) == 0 ) {  // Drift Tubes (near mu wall)

#       if CHECK_RESOURCES
	//+---------- Check Resource Usages, Benigno 20020226
	usage( "Drift Tubes" );
#       endif

	CsDriftTubeDetector* det =
	  new CsDriftTubeDetector( ndets, id,    name, TBname,unit,  type,
				   rdLen, xsiz,  ysiz, zsiz,
				   xcm,   ycm,   zcm,  _matrices[rotMNum-1],
				   rotWRS, wirDWRS,  angWRS,   nWir, wirP,
				   eff,   bkg, tGate,
				   vel,  t0, thRes, spSli, tiSli );
	_dets.push_back( det );

#       if CHECK_RESOURCES
	//+---------- Check Resource Usages, Benigno 20020226
	usage( "End Drift Tubes" );
#       endif

      }


      if( strncmp( name, "ST", 2 ) == 0 ) {  // Straw Tubes

#       if CHECK_RESOURCES
	//+---------- Check Resource Usages, Benigno 20020226
	usage( "Straw Tubes" );
#       endif

	CsStrawTubesDetector* det =
	  new CsStrawTubesDetector( ndets, id,    name, TBname,unit,  type,
				    rdLen, xsiz,  ysiz, zsiz,
				    xcm,   ycm,   zcm,  _matrices[rotMNum-1],
				    rotWRS, wirDWRS,  angWRS,   nWir, wirP,
				    eff,   bkg, tGate,
				    vel,  t0, thRes, spSli, tiSli );
	_dets.push_back( det );

#       if CHECK_RESOURCES
	//+---------- Check Resource Usages, Benigno 20020226
	usage( "End Straw Tubes" );
#       endif

      }

      if( strncmp( TBname, "DR", 2 ) == 0 ||  // Rich Wall official TB name
	  strncmp( TBname, "WD", 2 ) == 0) {  // Alternative name used in earlier MC

#       if CHECK_RESOURCES
	//+---------- Check Resource Usages, Benigno 20020226
	usage( "Rich Wall" );
#       endif

	CsRichWallDetector* det =
	  new CsRichWallDetector( ndets, id,    name, TBname,unit,  type,
		        	  rdLen, xsiz,  ysiz, zsiz,
				  xcm,   ycm,   zcm,  _matrices[rotMNum-1],
				  rotWRS, wirDWRS,  angWRS,   nWir, wirP,
				  eff,   bkg, tGate,
				  vel,  t0, thRes, spSli, tiSli );
	_dets.push_back( det );

#       if CHECK_RESOURCES
	//+---------- Check Resource Usages, Benigno 20020226
	usage( "End Rich Wall" );
#       endif

      }


      // MWPC like detectors

      if( strncmp( TBname, "PA", 2 ) == 0 ||
	  strncmp( TBname, "PB", 2 ) == 0 ||
	  strncmp( TBname, "PS", 2 ) == 0 ) {  // MWPC

#       if CHECK_RESOURCES
	//+---------- Check Resource Usages, Benigno 20020226
	usage( "MWPC" );
#       endif

	CsMWPCDetector* det =
	  new CsMWPCDetector( ndets, id,    name, TBname,unit,  type,
			      rdLen, xsiz,  ysiz, zsiz,
			      xcm,   ycm,   zcm,  _matrices[rotMNum-1],
			      rotWRS, wirDWRS,  angWRS,   nWir, wirP,
			      eff,   bkg, tGate, vel, t0, thRes,
                              spSli, tiSli );
	_dets.push_back( det );

#       if CHECK_RESOURCES
	//+---------- Check Resource Usages, Benigno 20020226
	usage( "End MWPC" );
#       endif

      }

      // Muon Wall 1 detectors

      if( strncmp( TBname, "MA", 2 ) == 0 ) {  // Muon Wall 1

#       if CHECK_RESOURCES
	//+---------- Check Resource Usages, Benigno 20020226
	usage( "MW1" );
#       endif

	CsMW1Detector* det =
	  new CsMW1Detector( ndets, id,    name, TBname,unit,  type,
			     rdLen, xsiz,  ysiz, zsiz,
			     xcm,   ycm,   zcm,  _matrices[rotMNum-1],
			     rotWRS, wirDWRS,  angWRS,   nWir, wirP,
			     eff,   bkg, tGate );
	_dets.push_back( det );

#       if CHECK_RESOURCES
	//+---------- Check Resource Usages, Benigno 20020226
	usage( "End MW1" );
#       endif

      }

      if( strncmp( TBname, "GM", 2 ) == 0 ) {  // GEMs

#       if CHECK_RESOURCES
	//+---------- Check Resource Usages, Benigno 20020226
	usage( "GEM" );
#       endif

        double spSig    = vel;     // space resolution, mm
        double eGain    = t0;      // effective gain a.u.
        double eGSig    = thRes;   // gain sigma, a.u.
	double sWidth   = spSli;   // signal width, mm
	double tmRes    = tiSli;   // time resolution, ns

	CsGEMDetector* det =
	  new CsGEMDetector( ndets, id,    name, TBname,unit,  type,
			     rdLen, xsiz,  ysiz, zsiz,
			     xcm,   ycm,   zcm,  _matrices[rotMNum-1],
			     rotWRS, wirDWRS,  angWRS,   nWir, wirP,
			     eff,   bkg, tGate, spSig, eGain, eGSig,
                             sWidth, tmRes );
	_dets.push_back( det );

#       if CHECK_RESOURCES
	//+---------- Check Resource Usages, Benigno 20020226
	usage( "End GEM" );
#       endif

      }

      if( strncmp( TBname, "GP", 2 ) == 0 ) {  // pixelGEMs

#       if CHECK_RESOURCES
	usage( "PixelGEM" );
#       endif

        double spSig    = vel;     // space resolution, mm
        double eGain    = t0;      // effective gain a.u.
        double eGSig    = thRes;   // gain sigma, a.u.
	double sWidth   = spSli;   // signal width, mm
	double tmRes    = tiSli;   // time resolution, ns  (Note: as of 06/12, overwritten by built-in value)

	CsPixelGEMDetector* det =
	  new CsPixelGEMDetector( ndets, id,    name, TBname,unit,  type,
				  rdLen, xsiz,  ysiz, zsiz,
				  xcm,   ycm,   zcm,  _matrices[rotMNum-1],
				  rotWRS, wirDWRS,  angWRS,   nWir, wirP,
				  eff,   bkg, tGate, spSig, eGain, eGSig,
				  sWidth, tmRes );
	_dets.push_back( det );

#       if CHECK_RESOURCES
	usage( "End PixelGEM" );
#       endif

      }

      if( strncmp( TBname, "MM", 2 ) == 0 ) {  // Micro-Mega

#       if CHECK_RESOURCES
	//+---------- Check Resource Usages, Benigno 20020226
	usage( "MicroMega" );
#       endif

	CsMicroMegaDetector* det =
	  new CsMicroMegaDetector( ndets, id,    name, TBname,unit,  type,
				   rdLen, xsiz,  ysiz, zsiz,
				   xcm,   ycm,   zcm,  _matrices[rotMNum-1],
				   rotWRS, wirDWRS,  angWRS,   nWir, wirP,
				   eff,   bkg, tGate );
	_dets.push_back( det );

#       if CHECK_RESOURCES
	//+---------- Check Resource Usages, Benigno 20020226
	usage( "End MicroMega" );
#       endif

      }

      if( strncmp( TBname, "MP", 2 ) == 0 ) {  // pixelMumegas

#       if CHECK_RESOURCES
	usage( "PixelMumega" );
#       endif

        double spSig    = vel;     // space resolution, mm
        double eGain    = t0;      // effective gain a.u.
        double eGSig    = thRes;   // gain sigma, a.u.
	double sWidth   = spSli;   // signal width, mm
	double tmRes    = tiSli;   // time resolution, ns 

	CsPixelMumegaDetector* det =
	  new CsPixelMumegaDetector( ndets, id,    name, TBname,unit,  type,
				  rdLen, xsiz,  ysiz, zsiz,
				  xcm,   ycm,   zcm,  _matrices[rotMNum-1],
				  rotWRS, wirDWRS,  angWRS,   nWir, wirP,
				  eff,   bkg, tGate, spSig, eGain, eGSig,
				  sWidth, tmRes );
	_dets.push_back( det );

#       if CHECK_RESOURCES
	usage( "End PixelMumega" );
#       endif

      }

      if( strncmp( TBname, "RP", 2 ) == 0 ) {  // RPD

#       if CHECK_RESOURCES
        usage( "RPD" );
#       endif

        CsRPDetector* det =
          new CsRPDetector( ndets, id,    name, TBname,unit,  type,
                            rdLen, xsiz,  ysiz, zsiz,
                            xcm,   ycm,   zcm,  _matrices[rotMNum-1],
                            rotWRS, wirDWRS,  angWRS,   nWir, wirP,
                            eff,   bkg, tGate );
        _dets.push_back( det );

#       if CHECK_RESOURCES
        usage( "End RPD" );
#       endif

      }


      if( strncmp( TBname, "SI", 2 ) == 0 ) {  // Silicon Tracker

#       if CHECK_RESOURCES
	//+---------- Check Resource Usages, Benigno 20020226
	usage( "Silicon" );
#       endif

	CsSiTrackerDetector* det =
	  new CsSiTrackerDetector( ndets, id,    name, TBname,unit,  type,
				   rdLen, xsiz,  ysiz, zsiz,
				   xcm,   ycm,   zcm,  _matrices[rotMNum-1],
				   rotWRS, wirDWRS,  angWRS,   nWir, wirP,
				   eff,   bkg, tGate );
	_dets.push_back( det );

#       if CHECK_RESOURCES
	//+---------- Check Resource Usages, Benigno 20020226
	usage( "End Silicon" );
#       endif

      }


      // Detectors with TDC read out

      if( strncmp( TBname, "FI", 2 ) == 0 ) {  // Fiber Hodoscopes

#       if CHECK_RESOURCES
	//+---------- Check Resource Usages, Benigno 20020226
	usage( "Fiber Hodo" );
#       endif

	int fibnumber = atoi( TBname+2 );

	if( fibnumber > 0 && fibnumber < 5 || fibnumber == 35) {  // Japanese Fibers
	  CsJFiberHodoDetector* det =
	    new CsJFiberHodoDetector( ndets, id,    name, TBname,unit,  type,
				      rdLen, xsiz,  ysiz, zsiz,
				      xcm,   ycm,   zcm,  _matrices[rotMNum-1],
				      rotWRS, wirDWRS,  angWRS,   nWir, wirP,
				      eff,   bkg, tGate, tRes );
	  _dets.push_back( det );
	}

	else if( ( fibnumber > 4 && fibnumber < 9 )
		 || fibnumber == 15 || fibnumber == 55  // German and Warsaw Fibers
		 || fibnumber == 35) /* DY Vertex detector */ {
	  CsDFiberHodoDetector* det =
	    new CsDFiberHodoDetector( ndets, id,    name, TBname,unit,  type,
				      rdLen, xsiz,  ysiz, zsiz,
				      xcm,   ycm,   zcm,  _matrices[rotMNum-1],
				      rotWRS, wirDWRS,  angWRS,   nWir, wirP,
				      eff,   bkg, tGate, tRes );
	  _dets.push_back( det );
	}

	else
	  CsErrLog::msg(elFatal,__FILE__,__LINE__,
			"\"%s \": wrong SciFi TBname",TBname);

#       if CHECK_RESOURCES
	//+---------- Check Resource Usages, Benigno 20020226
	usage( "End Fiber Hodo" );
#       endif

      }

      if(  strncmp( TBname, "H",  1 ) == 0 ||
	   strncmp( TBname, "V",  1 ) == 0) {  // Trigger Hodoscopes

#       if CHECK_RESOURCES
	//+---------- Check Resource Usages, Benigno 20020226
	usage( "Trigger Hodo" );
#       endif

	CsTriggerHodoDetector* det =
	  new CsTriggerHodoDetector( ndets, id,    name, TBname,unit,  type,
				     rdLen, xsiz,  ysiz, zsiz,
				     xcm,   ycm,   zcm,  _matrices[rotMNum-1],
				     rotWRS, wirDWRS,  angWRS,   nWir, wirP,
				     eff,   bkg, tGate, tRes );
	_dets.push_back( det );

#       if CHECK_RESOURCES
	//+---------- Check Resource Usages, Benigno 20020226
	usage( "End Trigger Hodo" );
#       endif

      }

      if( strncmp( TBname, "BM", 2 ) == 0 ) {  // BMS

#       if CHECK_RESOURCES
	//+---------- Check Resource Usages, Benigno 20020226
	usage( "BMS" );
#       endif

	CsBMSDetector* det =
	  new CsBMSDetector( ndets, id,    name, TBname,unit,  type,
			     rdLen, xsiz,  ysiz, zsiz,
			     xcm,   ycm,   zcm,  _matrices[rotMNum-1],
			     rotWRS, wirDWRS,  angWRS,   nWir, wirP,
			     eff,   bkg, tGate, tRes );
	_dets.push_back( det );

#       if CHECK_RESOURCES
	//+---------- Check Resource Usages, Benigno 20020226
	usage( "End BMS" );
#       endif

      }


      // Check if there are some shifts on this detector for alignment tests
      if( CsInit::Instance()->IsAMonteCarloJob() ) { // just a protection...
	CsDetector* det = _dets.back();
	stringstream mykey;
	mykey << det->getName() << " " << det->getUnit();
	double shift;
	string mytag = "detector X shift";
	if( CsOpt::Instance()->getOpt( mytag, mykey.str(), shift ) ||
	    CsOpt::Instance()->getOpt( mytag, det->GetTBName(), shift ) ) {
	  cout << "--------------------------------------------------" << endl
	       << "WARNING!!! You are adding a "
	       << shift
	       << " mm shift on "
	       << det->GetTBName() << " (" << mykey.str() << ") "
	       << " X position!" << endl
	       << "--------------------------------------------------" << endl;
	  det->setShiftOnX( shift );
	}
	mytag = "detector Y shift";
	if( CsOpt::Instance()->getOpt( mytag, mykey.str(), shift ) ||
	    CsOpt::Instance()->getOpt( mytag, det->GetTBName(), shift ) ) {
	  cout << "--------------------------------------------------" << endl
	       << "WARNING!!! You are adding a "
	       << shift
	       << " mm shift on "
	       << det->GetTBName() << " (" << mykey.str() << ") "
	       << " Y position!" << endl
	       << "--------------------------------------------------" << endl;
	  det->setShiftOnY( shift );
	}
      }

    }
    else if( opt == "cedar" ) {
      string TBname; s>>TBname;
      if( strncmp( TBname.c_str(), "CE", 2 ) == 0 ) {  // CEDAR

        #       if CHECK_RESOURCES
        //+---------- Check Resource Usages, Benigno 20020226
        usage( "CEDAR" );
        #       endif

        CsCEDARDetector* det =
          new CsCEDARDetector( TBname,  *detTableName ) ;
        _cedar_dets.push_back(det);
        std::sort(_cedar_dets.begin(), _cedar_dets.end(), det->compById);

        #       if CHECK_RESOURCES
        //+---------- Check Resource Usages, Benigno 20020226
        usage( "End CEDAR" );
        #       endif

      }
    }
    // calo
    else if( opt == "calo" )
    {

      int id;
      string TBname;

      s >> TBname;
      std::cout << "###### TBname= " << TBname << endl;
      // It is not possible to add calorimeters to _dets list, because
      // calorimeters are not derived from CsDetector.  So all calorimeters
      // will be added only to internal static list of CsDet class.
      
      if ( TBname == "EC01P1__" )
	_calorimeters.push_back( new CsECAL1(TBname, *detTableName ) );
      else if ( TBname == "EC02P1__" )
	_calorimeters.push_back( new CsECAL2(TBname, *detTableName ) );
      else if ( TBname == "HC01P1__" )
	_calorimeters.push_back( new CsHCAL1(TBname, *detTableName ) );
      else if ( TBname == "HC02P1__" )
	_calorimeters.push_back( new CsHCAL2(TBname, *detTableName ) );
      else if ( TBname == "EC00P1__" )
      	_calorimeters.push_back( new CsECAL0(TBname, *detTableName ) );
      else
	throw CS::Exception( "Unknown calorimeter %s in %s!", TBname.c_str(), detTableName->c_str() );
    }
    // detector dead zone?
    else if( opt == "dead" ) {

      int    id;      s >> id;          // detector number
      char TBname[9] = " ";
      if ((CGEAversion==6 && CGEArelease>=7) || CGEAversion>6 )
	s>>TBname;                      // det. TB name
      char   name[5]; s >> name;        // detector name
      int    unit;    s >> unit;        // detector number (in list of likes)
      int    shape;   s >> shape;       // shape; rectangular or circular
      double xsiz;    s >> xsiz; xsiz *= 10.; // dead zone size (mm)
      double ysiz;    s >> ysiz; ysiz *= 10.;
      double zsiz;    s >> zsiz; zsiz *= 10.;
      CsGeom::GeaRef2CsRefVec( xsiz, ysiz, zsiz );
      double xcm;     s >> xcm; xcm *= 10.;   // dead zone centre (MRS) (mm)
      double ycm;     s >> ycm; ycm *= 10.;
      double zcm;     s >> zcm; zcm *= 10.;
      CsGeom::GeaRef2CsRefVec( xcm, ycm, zcm );
      int    rotMNum; s >> rotMNum;     // rotation matrix number

      // Crude check: When don't pass through, would otherwise pose problem
      if (!(0<id&&id<99999 && 0<unit&&unit<99 && 0<shape&&shape<99 &&
	    0<rotMNum&&rotMNum<99))
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
		      "Dead Zone - Corrupted file @ %s!",TBname);

      // Loop on CsDetector's list for currently referenced det
      list<CsDetector*>::iterator idet;	// Iterator on CsDetector's list
      for (idet = _dets.begin(); idet!=_dets.end(); idet++) {
	if (
	  // Same name, and unit and id (for backward compatibility)...
	  ( strncmp((*idet)->GetID().GetName().c_str(),name,4)==0 &&
            (*idet)->getUnit()==unit &&
            (*idet)->GetID().GetNumber()==(unsigned int)id )
	   // ... or same TBname
	   || strncmp((*idet)->GetTBName().c_str(),TBname,8)==0 )
	  break;
      }
      if (idet==_dets.end()) {
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
		      "Dead Zone - Non existing det id=%d %s,%d",id,name,unit);
      }

      // We want to implement in "coral" all the possibilities permitted by the
      // very general description of the dead zones in "DetTable". (Still
      // we do not provide for Z offset: dead zone zcm = detector's zcm.)
      // BUT we  also want to have, besides this general description, a
      // simplified one which allows to speed up "inActiveArea" method in some
      // simple but prevalent cases. We retain such 2...
      //  i) A rectangular dead zone with sides parallel to the detector's
      //    own sides.
      // ii) A circular dead zone.
      // ..and flag them
      // => Check that specification just read agrees w/ either of them
      // => Fill CsDetector members:
      //    - dZtype = Shape + 16*FAST_ENABLED
      //    - dZdim = X 1/2size or radius squared
      //    - dZysiz = Y 1/2size (if rectangular shape)
      //    - dZxdrs, dZydrs
      //    - rotDZRS

      // ********** CHECK DZ's ORIENTATION wrt. DETECTOR's **********

      HepMatrix mDZ  = _matrices[rotMNum-1];
      // Get detector's orientation
      // Looks like "(*idet)->getRotDRS()" does not fit => instead
      CsDetector *det = *idet;

      if ( det->getDZDim() != 0. || det->getDZYdim() != 0. ) {
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
		      "Duplicate Dead Zone @ id=%d %s,%d",id,name,unit);
      }
      
      HepMatrix mdet = det->getRotDRS();
      int iD2DZ; HepMatrix D2DZ = mDZ.inverse(iD2DZ)*mdet;
      if (iD2DZ) {
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
		      "Dead Zone - Bad Orientation @ id=%d %s,%d",id,name,unit);
      }
      if (shape==1) {                            // Rectangular shape...
	// Require matrix to be equal to detector's
	// "if (mdet!=mDZ)"  does not work => intead
	bool equal = true;
	for (int i = 0; i<3; i++) for (int j = 0; j<3; j++)
	  if (mdet[i][j]!=mDZ[i][j]) { equal = false; break; }
	if (equal) shape |= 0x10;                    // Set ENABLE_FAST flag
      }
      else if (shape==5) {                       // Circular shape..
	// Require DRS to DZRS to be a rotation about Y axis times
	// the circular permutation of axes bringing Y along Z
	int i, diff;
	for (i=diff = 0; i<3; i++) for (int j = 0; j<3; j++) if ((i==1)^(j==2))
	  if (fabs(D2DZ[i][j])>1e-6) diff = 1;
	if (fabs(D2DZ[1][2]-1)>1e-6) diff = 1;
	if (fabs(D2DZ[2][0]*D2DZ[0][1]-D2DZ[0][0]*D2DZ[2][1]-1)>
	    5e-5) // Given the bad precision of rotation matrices in "detectors.dat"
	  diff = 1;
	if (!diff) shape |= 0x10;                    // Set ENABLE_FAST flag
      }
      else {
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
	  "Dead Zone - Unexpected Shape (=%d) @ id=%d %s,%d",shape,id,name,unit);
      }
      if (!(shape&0x10)) {
	CsErrLog::msg(elWarning,__FILE__,__LINE__,
	  "Dead Zone: Skewed Orientation @ id=%d %s,%d => No fast inActiveArea",
		      id,name,unit);
      }
      else if (shape==0x11 || shape==0x15) {
	if (zcm!=(*idet)->getZcm()) {
	  shape &= 0xf;                              // Cancel ENABLE_FAST flag
	  CsErrLog::msg(elWarning,__FILE__,__LINE__,
	    "Dead Zone - Skewed Centering @ id=%d %s,%d => No fast inActiveArea",
			id,name,unit);
	}
      }
      // Set dead zone members
      double xcd, ycd, zcd;             // Dead zone center...
      det->rotatePointMRS2DRS(xcm,ycm,0,xcd,ycd,zcd);
      if ((shape&0xf)==1)                     // Rectangular shape
	det->setDeadZone(shape,xsiz/2,ysiz/2,xcd,ycd,D2DZ);
      else                                    // Circular shape...
	det->setDeadZone(shape,xsiz*xsiz/4,   // ...= radius squared
			 0,xcd,ycd,D2DZ);
    }

    // rich1 / added by Paolo : 23/8/00
    else if( opt == "rich" ) {

#     if CHECK_RESOURCES
      //+---------- Check Resource Usages, Benigno 20020226
      usage( "RICH1" );
#     endif

      int id;            s >> id;            // detector number

      bool rich_up_grade = CsOpt::Instance()->getOpt("CsRICH1UpGrade","CONFIG");
      if( rich_up_grade )
        _rich1 = new CsRICH1UpGrade( id, "RI01P___" );
      else
        _rich1 = new CsRICH1Detector( id, "RI01P___" );

#     if CHECK_RESOURCES
      //+---------- Check Resource Usages, Benigno 20020226
      usage( "End RICH1" );
#     endif

    }
//- rev. by Paolo 090213
    else if( opt == "rich1" ) {
      string rindex;     s >> rindex;
      if( rindex == "rindex" ) {
	double index;    s >> index;
	if( _rich1 != NULL ) {
	  _rich1->setMCCFRefInd( index );
	}
      }
    //}
    //else if( opt == "rich1" ) {
      //string rindex;     s >> rindex;
      if( rindex == "rindexUV" ) {
	double index;    s >> index;
	if( _rich1 != NULL ) {
	  _rich1->setMCCFRefInd( index );
	}
      }
    //}
    //else if( opt == "rich1" ) {
      //string rindex;     s >> rindex;
      if( rindex == "rindexVS" ) {
	double index;    s >> index;
	if( _rich1 != NULL ) {
	  _rich1->setMCCFRefIndVS( index );
	}
      }
    }
    else if( opt == "phot0" ) {
      string name;       s >> name;
      double xDet0;      s >> xDet0; xDet0 *= 10.; // data are in COMGEANT ref
      double yDet0;      s >> yDet0; yDet0 *= 10.;
      double zDet0;      s >> zDet0; zDet0 *= 10.;
      CsGeom::GeaRef2CsRefVec( xDet0, yDet0, zDet0 );

      HepMatrix rotMRC(3,3);
      int dummy; s >> dummy;
      for( i=0; i<3; i++ ) {
	for( j=0; j<3; j++ ) {
          s >> rotMRC(i+1,j+1);
	}
      }
      CsGeom::GeaRef2CsRefMat( rotMRC );

      if( _rich1 != NULL ) {
	_rich1->setPhotonDet( name, xDet0, yDet0, zDet0, rotMRC );
      }
    }
    else if( opt == "cath" ) {
      int id;            s >> id;
      string TBname;     s >> TBname;
      string name;       s >> name;
      double xOffCat0;   s >> xOffCat0;	xOffCat0 *= 10.;
      double yOffCat0;   s >> yOffCat0;	yOffCat0 *= 10.;
      double zOffCat0;   s >> zOffCat0;	zOffCat0 *= 10.;
      CsGeom::GeaRef2CsRefVec( xOffCat0, yOffCat0, zOffCat0 );

      int nPadx;         s >> nPadx;
      int nPady;         s >> nPady;
      double padx;       s >> padx; padx *= 10.;
      double pady;       s >> pady; pady *= 10.;
      double ddQzW;      s >> ddQzW; ddQzW *= 10.;
      double ddGap;      s >> ddGap; ddGap *= 10.;

      HepMatrix rotMRC(3,3);
      int dummy; s >> dummy;
      for( i=0; i<3; i++ ) {
	for( j=0; j<3; j++ ) {
          s >> rotMRC(i+1,j+1);
	}
      }
      CsGeom::GeaRef2CsRefMat( rotMRC );

      if( _rich1 != NULL ) {
	_rich1->setCathode( id, TBname, name,
			    xOffCat0, yOffCat0, zOffCat0,
			    nPadx, nPady, padx, pady, ddQzW, ddGap,
			    rotMRC );
      }
    }
    else if( opt == "mirr0" ) {
      string name;       s >> name;
      double xC0;      s >> xC0; xC0 *= 10.;
      double yC0;      s >> yC0; yC0 *= 10.;
      double zC0;      s >> zC0; zC0 *= 10.;
      CsGeom::GeaRef2CsRefVec( xC0, yC0, zC0 );
      double RR;       s >> RR; RR *= 10;

      if( _rich1 != NULL ) {
	_rich1->setMirrNom( name, xC0, yC0, zC0, RR );
      }
    }
    else if( opt == "mirr" ) {
      string name;       s >> name;
      double theta;      s >> theta;
      double phi;        s >> phi;
      double RR;         s >> RR; RR *= 10.;
      double deTheta;    s >> deTheta;
      double dePhi;      s >> dePhi;
      double delta;      s >> delta;
      double qfact;      s >> qfact;
      int    align;      s >> align;
      if( _rich1 != NULL ) {
	_rich1->setMirrEle( name, theta, phi, RR, deTheta, dePhi,
			    delta, qfact, align );
      }
    }

    else if( opt == "misc" ) {
      char TBname[9]=""; s>>TBname; // det. TB name
      // Here two detectors that just need decoding
      if(  strncmp( TBname, "HMSC1", 5 ) == 0 ||
	   strncmp( TBname, "HVETO", 5 ) == 0) {
	CsMiscDetector* det = new CsMiscDetector( TBname );
	_others.push_back( det );
      }
    }


    // RPD detectror.
    // for time being it's just like normal tracking "det".
    // To be modified for more specific geometry description
    //

    else if( opt == "rpd" ) {
      double ddummy = 0; // Init, in order to avoid "might be used uninitialized ..."
      HepMatrix rotDummy(3,3,1);

      int    id;      s >> id;          // detector number
      char TBname[9]=""; s>>TBname; // det. TB name
      char   name[5]; s >> name;        // detector name
      int    unit;    s >> unit;        // detector number in station
      int    type;    s >> type;        // detector type
      double rdLen;   s >> rdLen;       // radiation length
      double xsiz;    s >> xsiz; xsiz *= 10.; // detector size (mm)
      double ysiz;    s >> ysiz; ysiz *= 10.;
      double zsiz;    s >> zsiz; zsiz *= 10.;
      CsGeom::GeaRef2CsRefVec( xsiz, ysiz, zsiz );
      double xcm;     s >> xcm; xcm *= 10.; // detector centre (MRS) (mm)
      double ycm;     s >> ycm; ycm *= 10.;
      double zcm;     s >> zcm; zcm *= 10.;
      CsGeom::GeaRef2CsRefVec( xcm, ycm, zcm );
      int    rotMNum; s >> rotMNum;     // rotation matrix number
      double wirD;    s >> wirD; wirD *= 10.;  // 1st wire offset (mm)
      double ang;     s >> ang;         // angle of wires in DRS (degrees)
      int    nWir;    s >> nWir;        // number of wires
      double wirP;    s >> wirP; wirP *= 10.;  // wires pitch
      double eff;     s >> eff;         // detector efficiency
      double bkg;     s >> bkg;         // detector background
      double tGate;   s >> tGate;       // detector time gate

      CsRPDetector* det =
	new CsRPDetector( ndets, id,    name, TBname,unit,  type,
			  rdLen, xsiz,  ysiz, zsiz,
			  xcm,   ycm,   zcm,  _matrices[rotMNum-1],
			  rotDummy, ddummy,  ddummy,   nWir, wirP,
			  eff,   bkg, tGate);
      _dets.push_back( det );

    } // end of RPD block

  } // End of loop on reading detector table

  // Check if f readout ended "somehow" correctly:
  if( failure ) {
    CsErrLog::mes( elFatal, "It seems there's something wrong in detectors.dat file, please have a look to it." );
    exit(1);
  }

  if (CGEAversion>=7 && (CGEAversion>7 || CGEArelease>=3)) {
    //      ***** LATER RELEASES of COMGeant: >= 7.3 *****
    // => CHECK TARGET INFO ("targ" entries in the detector table) has been read
    if (_targetCells.size()==0)
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
	"No target info (\"targ\" entry) found, in this COMGeant,v%d.%d table",
		    CGEAversion,CGEArelease);
    // => REPORT on TARGET CENTER eventually derived from the info
    CsErrLog::msg(elInfo,__FILE__,__LINE__,
      "Target center assigned to weighted average of %d cell(s) = %f",
		  _targetCells.size(),_targetCenter);
  }

  // Calorimeters initialization
  for( vector<CsCalorimeter*>::iterator it = getCalorimeters().begin(); it!=getCalorimeters().end(); it++ ) {
      (*it)->Initialize();
  } 
  
  _dets.sort( CsGeom::sortDetectors_() );

  // Read the Magnetic Field maps
  _field.ReadMaps();

  // Read the Material maps
  _matmap.ReadMaps( _geomVersion );

  _detTableName = *detTableName;

# if CHECK_RESOURCES
  //+---------- Check Resource Usages, Benigno 20020226
  usage( "Exiting CsGeom::readDetTable" );
# endif

  return( true );

}

void CsGeom::GeaRef2CsRefVec( double& x, double& y, double& z ) {
  // From GeantRS : X // beam, Z vertical
  // To MRS :       Z // beam, Y verical
  double tmp = x;
  x          = y;
  y          = z;
  z          = tmp;
}

void CsGeom::GeaRef2CsRefMat( HepMatrix& a ) {
  // From GeantRS : X // beam, Z vertical
  // To MRS :       Z // beam, Y verical
  HepMatrix r(3,3);
  double set[] = { 0, 1, 0, 0, 0, 1, 1, 0, 0 };
  for( int i=0; i<9; i++ ) r( i/3+1, i%3+1 ) = set[i];
  a = r * a;
  a = a * r.T();
}

//____________________________
bool CsGeom::associateDets() {

  CsErrLog::mes(elInfo,"======= Detector association for LR ======");

  CsOpt* opt = CsOpt::Instance();

  list<string> assocList;    // must have three items
  list<string>::iterator Is; // corresponding iterator
  CsDetector* det0 = NULL; // detectors to be associated
  CsDetector* det1 = NULL;
  list<CsDetector*>::iterator idet;	// iterator on CsDetector lists

  string TBName0, TBName1; double Acut = -1, LRcut = 0.5;

  //             ********** TraFFiC TURNED-OFF DETECTORS **********
  // (We may or may not want to perform LR ambiguity raising on detectors
  // that are turned-off in tracking.
  //  I) The tricky case is when only one plane in a given associated pair is
  //   turned off. This can only be done purposefully in order to monitor the
  //   response (in efficiency, resolution,...) of the turned-off plane. And
  //   performing LR ambiguity raising on it in these conditions introduces a
  //   bias, by facilitating the reco of those tracks where the response of the
  //   plane under exam has been OK.
  //  => The best way to avoid the problem is to forbid the turning off single
  //    planes. This is done by default.
  // II) Then remains the case where a pair of detectors are turned off. In that
  //    case, we are sure LR raising will have no impact on the reconstruction.
  //    And it can be useful for the analysis the residuals.
  //  => We decide to perform then the ambiguity raising.
  //  But bypassing ambiguity raising can be enabled by a special option. If
  //  indeed enabled, the possibility to turn off a single plane is consented.
  list<string> detOff; opt->getOpt( "TraF","DetNameOff",detOff);
  bool doTurnOff = opt->getOpt("LR","DetNameOff");

  int nMWPCs = 0;
  while (opt->getOptRec("LR","associate",assocList)) {

    // *************** LOOP on ALL "LR associate" ENTRIES ***************

    det0 = 0; det1 = 0;
    // Parameters #3 and #4 are optionnal.
    // par #3 is LRCut;  Default is 0.5, or last given value.
    // par #4 is LRMode; Default is 0 (target pointing).
    if (assocList.size()<3)
      CsErrLog::mes(elFatal,"Wrong format in association list!");
    int i = 0, LRMode = 0;

    // ********** GET ASSOCIATION TBNAMES, CUT and OPTIONAL PAR's **********
    for (Is = assocList.begin(); Is!=assocList.end(); Is++, i++) {
      istringstream s( *Is );
      switch (i) {
      case 0: s >> TBName0;     break;
      case 1: s >> TBName1;     break;
      case 2: s >> Acut;        break;
      case 3: s >> LRcut; 	break;
      case 4: s >> LRMode; 	break;
      default: break;
      }
    }

    //      ********** FIND DETECTORS MATCHING TBName0 and TBName1 **********
    for (idet = _dets.begin(); idet!=_dets.end(); idet++) {
      if ((*idet)->GetTBName()==TBName0) det0 = *idet;
      if ((*idet)->GetTBName()==TBName1) det1 = *idet;
    }
    CsErrLog::msg(elInfo,__FILE__,__LINE__,"%s <-> %s.",
		  TBName0.c_str(),TBName1.c_str());
    if (det0==0) {
      CsErrLog::msg(elWarning,__FILE__,__LINE__,"%s not found.",
		    TBName0.c_str()); continue;
    }
    if (det1==0) {
      CsErrLog::msg(elWarning,__FILE__,__LINE__,"%s not found.",
		    TBName1.c_str()); continue;
    }

    //      ********** CHECK IF DETECTOR ARE OFF FOR TraFFiC **********
    int isOff; list<string>::iterator IDO;
    for (IDO = detOff.begin(), isOff = 0; IDO!=detOff.end(); IDO++) {
      if (TBName0.find(*IDO)==0) isOff |= 0x1;
      if (TBName1.find(*IDO)==0) isOff |= 0x2;
    }
    if (isOff) {
      if (isOff!=0x3 && !doTurnOff) {
	// One only plane is off AND special option not asked for => Fatal
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
"%s %s: One plane is OFF while counterpart is ON. This is not recommended. Enter \"LR DetNameOff\" if you still want to carry on...",TBName0.c_str(),TBName1.c_str());
      }
      if (doTurnOff) {
	CsErrLog::msg(elError,__FILE__,__LINE__,
"=> Bypassing LR ambiguity raising on (%s,%s)",TBName0.c_str(),TBName1.c_str());
	continue;
      }
    }

    if (LRMode<0 || 2<LRMode) // ********** WRONG OPTION ARG's **********
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
	"LR mode for (\"%s\",\"%s\") = %d non valid",
		    TBName0.c_str(),TBName1.c_str(),LRMode);

    bool accept = false;

    // DRIFT CHAMBERS
    CsDriftChamberDetector* driftdet0 = dynamic_cast<CsDriftChamberDetector*>(det0);
    CsDriftChamberDetector* driftdet1 = dynamic_cast<CsDriftChamberDetector*>(det1);
    if( driftdet0 != NULL
      && driftdet1 != NULL
      && driftdet0->isMWPC() == false
      && driftdet1->isMWPC() == false ) {
      accept = true;
      // associate detectors
      if( ! ( driftdet0->hasAssociateDet() || driftdet1->hasAssociateDet() ) ) {
	driftdet0->setAssociateDet( *det1 );
	driftdet0->setAssociationCut( Acut );
	driftdet0->setLRProbCut( LRcut );
	driftdet0->setLRMode( LRMode );

	driftdet1->setAssociateDet( *det0 );
	driftdet1->setAssociationCut( Acut );
	driftdet1->setLRProbCut( LRcut );
	driftdet1->setLRMode( LRMode );
      }
    }

    // STRAW TUBES
    CsStrawTubesDetector* strawdet0 = dynamic_cast<CsStrawTubesDetector*>(det0);
    CsStrawTubesDetector* strawdet1 = dynamic_cast<CsStrawTubesDetector*>(det1);
    if( strawdet0 != NULL
      && strawdet1 != NULL
      && strawdet0->isMWPC() == false
      && strawdet1->isMWPC() == false ) {
      accept = true;
      // associate detectors
      if( ! ( strawdet0->hasAssociateDet() || strawdet1->hasAssociateDet() ) ) {
	strawdet0->setAssociateDet( *det1 );
	strawdet0->setAssociationCut( Acut );
	strawdet0->setLRProbCut( LRcut );
	strawdet0->setLRMode( LRMode );

	strawdet1->setAssociateDet( *det0 );
	strawdet1->setAssociationCut( Acut );
	strawdet1->setLRProbCut( LRcut );
	strawdet1->setLRMode( LRMode );
      }
    }

    // DRIFT TUBES
    CsDriftTubeDetector* dtubedet0 = dynamic_cast<CsDriftTubeDetector*>(det0);
    CsDriftTubeDetector* dtubedet1 = dynamic_cast<CsDriftTubeDetector*>(det1);
    if( dtubedet0 != NULL
      && dtubedet1 != NULL
      && dtubedet0->isMWPC() == false
      && dtubedet1->isMWPC() == false ) {
      accept = true;
      // associate detectors
      if( ! ( dtubedet0->hasAssociateDet() || dtubedet1->hasAssociateDet() ) ) {
	dtubedet0->setAssociateDet( *det1 );
	dtubedet0->setAssociationCut( Acut );
	dtubedet0->setLRProbCut( LRcut );
	dtubedet0->setLRMode( LRMode );

	dtubedet1->setAssociateDet( *det0 );
	dtubedet1->setAssociationCut( Acut );
	dtubedet1->setLRProbCut( LRcut );
	dtubedet1->setLRMode( LRMode );
      }
    }

    // RICH WALL
    CsRichWallDetector* richwdet0 = dynamic_cast<CsRichWallDetector*>(det0);
    CsRichWallDetector* richwdet1 = dynamic_cast<CsRichWallDetector*>(det1);
    if( richwdet0 != NULL
      && richwdet1 != NULL
      && richwdet0->isMWPC() == false
      && richwdet1->isMWPC() == false ) {
      accept = true;
      // associate detectors
      if( ! ( richwdet0->hasAssociateDet() || richwdet1->hasAssociateDet() ) ) {
	richwdet0->setAssociateDet( *det1 );
	richwdet0->setAssociationCut( Acut );
	richwdet0->setLRProbCut( LRcut );
	richwdet0->setLRMode( LRMode );

	richwdet1->setAssociateDet( *det0 );
	richwdet1->setAssociationCut( Acut );
	richwdet1->setLRProbCut( LRcut );
	richwdet1->setLRMode( LRMode );
      }
    }

    // W45
    CsDWDetector* dwdet0 = dynamic_cast<CsDWDetector*>(det0);
    CsDWDetector* dwdet1 = dynamic_cast<CsDWDetector*>(det1);
    if( dwdet0 != NULL
      && dwdet1 != NULL
      && dwdet0->isMWPC() == false
      && dwdet1->isMWPC() == false ) {
      accept = true;
      // associate detectors
      if( ! ( dwdet0->hasAssociateDet() || dwdet1->hasAssociateDet() ) ) {
	dwdet0->setAssociateDet( *det1 );
	dwdet0->setAssociationCut( Acut );
	dwdet0->setLRProbCut( LRcut );
	dwdet0->setLRMode( LRMode );

	dwdet1->setAssociateDet( *det0 );
	dwdet1->setAssociationCut( Acut );
	dwdet1->setLRProbCut( LRcut );
	dwdet1->setLRMode( LRMode );
      }
    }

    // some outputs
    if( accept )
      CsErrLog::msg(elInfo,__FILE__,__LINE__,"Acut: %f mm, LRcut %s",
		    Acut,LRMode ? "COvlap" : "TPoint");
    else {
      cout<<"CsGeom::associateDets: "<<TBName0<<","<<TBName1<<" !associated\n";
      nMWPCs++;
    }
  }

  if (nMWPCs)
    CsErrLog::mes(elError,"At least one det treated as MWPC.");

  return( true );
}



CsMaterialMap* CsGeom::getCsMaterialMap()
{

  if( !_matmap.usingROOTGeometry() && !_matmap.getNofMaps() )
    CsErrLog::mes(elFatal,"No one Material Map was read.");

  return( &_matmap );
}


string CsGeom::getGeomVers() {
  int version = atoi( _comgeantVersion.c_str() + 1 );
  int release = atoi( _comgeantVersion.c_str() + 7 );
  if( version > 6 || ( version == 6 && release > 2 ) ) {
    return ( _geomVersion );
  }
  else {
    CsErrLog::mes(elError,"Information not available for this COMGEANT release");
    string dummy = "";
    return( dummy );
  }
}


bool CsGeom::isHadronSetup() {
  int version = atoi( _comgeantVersion.c_str() + 1 );
  int release = atoi( _comgeantVersion.c_str() + 7 );
  if( version > 6 || ( version == 6 && release > 2 ) ) {
    return ( _hadronSetup );
  }
  else {
    CsErrLog::mes(elError,"Information not available for this COMGEANT release");
    return( false );
  }
}

bool CsGeom::isMuonSetup() {
  int version = atoi( _comgeantVersion.c_str() + 1 );
  int release = atoi( _comgeantVersion.c_str() + 7 );
  if( version > 6 || ( version == 6 && release > 2 ) ) {
    return ( _muonSetup );
  }
  else {
    CsErrLog::mes(elError,"Information not available for this COMGEANT release");
    return( false );
  }
}

//____________________________
void CsGeom::associateGEMs(int mode) {

  // ******************** mode&0x1: ASSOCIATE GEM PLANES ********************
  // - I.e. CsGEMDetector or CsPixelGEMDetector
  // - W/ same TB name, apart from the orientation (viz. name[4]).
  // - X<->Y and U<->V.
  // -
  // *************** mode&0x2: ASSOCIATE pixelGEM's SUB-PIECES ***************
  // - I.e. a "GP..P", the logical entity that describes the inner pixelised
  //  core of the physical pixelGEMs entities, w/ a "GP..XYUV" stripped piece.
  // - W/ same TB name, apart from the orientation (viz. name[4]) and, possibly,
  //  the plane number (viz. name[5], e.g. GP02P2->GP02X1).
  // - W/ same abscissa, given +/-0.1mm.
  // - This association is used to get a single COMGeant ID, encompassing both
  //  the pixelised central piece of a pixelGEM and one of its 2 orthogonal
  //  stripped external pieces, feed MC hits to both kind of pieces. It has had
  //  to be introduced in order to provide for COMGeant setup version 2007.05_1.
  //  May not be usefull in the long term.

  CsOpt *opt = CsOpt::Instance(); CsInit *init = CsInit::Instance();

  //             ********** TraFFiC TURNED-OFF DETECTORS **********
  // (At present GEM association is only used for MC simulation of amplitude
  // correlation. And there's no point in bypassing this simulation when any
  // one of the associated detector planes is turned off in tracking. It can
  // even be harmful, for it will have the random generation depend upon ON/OFF.
  //  But in the future this association (the CsGEM member data hereafter
  // defined, that is) may be used explicitly in the reco. Then we may want not
  // to associate turned-off detectors. A special option is available for
  // achieving this.)
  list<string> detOff;
  if (opt->getOpt("GEM","DetNameOff")) opt->getOpt( "TraF","DetNameOff",detOff);

  list<CsDetector*>::iterator idet0;
  for (idet0 = _dets.begin(); idet0!=_dets.end(); idet0++) {

    //  *************** LOOP ON CsDetector's ***************
    const string &TBName0 = (*idet0)->GetTBName();

    CsPixelGEMDetector *pixGEM0 = dynamic_cast<CsPixelGEMDetector*>(*idet0);
    if (pixGEM0 && TBName0[4]=='P') {
      if (!(mode&0x2)) continue;

      //     ********** CsDetector IS A PIXEL PIECE OF A CsPixelGEM **********
      // ***** ASSOCIATE PIXEL PIECE -> NEIGHBOUR [MASTER] STRIPPED PIECE *****

      list<CsDetector*>::iterator idet1;
      for (idet1 = _dets.begin(); idet1!=_dets.end(); idet1++) {
	if (idet1==idet0) continue;
	CsPixelGEMDetector *pixGEM1 = dynamic_cast<CsPixelGEMDetector*>(*idet1);
	if (!pixGEM1) continue;
	if (TBName0.compare(0,4,pixGEM1->GetTBName(),0,4)) continue;
	if (fabs(pixGEM1->getZcm()-pixGEM0->getZcm())>.11) continue;
	if ((opt->getOpt(pixGEM1->GetTBName(),"ampCorrelationMC") ||
	     // Amplitude correlation requested? require associate to be master
	     opt->getOpt("GP","ampCorrelationMC")) && !pixGEM1->isMaster())
	  continue;
	pixGEM0->setAssociateDet(*pixGEM1);
      }
      if (!pixGEM0->getAssociateDet())
	CsErrLog::msg(elFatal,__FILE__,__LINE__,"CsPG \"%s\" -> CsGEM association failed",TBName0.c_str());
      continue;
    }

    CsGEMDetector *GEM0 = dynamic_cast<CsGEMDetector*>(*idet0);
    if ( ( !GEM0 && !pixGEM0 ) || !(mode&0x1)) continue;

    //       ********** CsDetector IS A Cs(Pixel)GEMDetector **********
    //         ********** ASSOCIATE X<->Y and V<->U **********
    bool doAmpCorrMC = init->IsAMonteCarloJob() &&
      (opt->getOpt(TBName0,"ampCorrelationMC") ||
       opt->getOpt(TBName0.substr(0,2),"ampCorrelationMC"));

    string TBName1(TBName0); char let[] = "XYUV", let1[] = "YXVU";
    for (int i = 0; i<4; i++) if (TBName0[4]==let[i]) TBName1[4]=let1[i];
    CsDet *pdet1 = CsDet::FindDetector(TBName1); if (pdet1) {
      CsGEMDetector *GEM1 = dynamic_cast<CsGEMDetector*>(pdet1);
      CsPixelGEMDetector *pixGEM1 = dynamic_cast<CsPixelGEMDetector*>(pdet1);
      if (GEM1 || pixGEM1) {
	list<string>::iterator IDO;        // ***** CHECK DETECTORS ARE NOT OFF
	for (IDO = detOff.begin(); IDO!=detOff.end(); IDO++) {
	  if (TBName0.find(*IDO)==0 || TBName1.find(*IDO)==0) continue;
	}

	if (GEM0) {
	  if (!GEM1)
	    CsErrLog::msg(elFatal,__FILE__,__LINE__,"CsGEM \"%s\" -> CsGEM association failed",TBName0.c_str());
	  if (doAmpCorrMC && // ***** CASE MONTE CARLO  AMPLITUDE CORRELATION...
	      !(GEM0->isMaster()^GEM1->isMaster()))//...REQUIRE MASTER XOR SLAVE
	    CsErrLog::msg(elFatal,__FILE__,__LINE__,"GEM coupled planes %s and %s both declared as \"s\", check your options",TBName0.c_str(),TBName1.c_str(),GEM0->isMaster()?"Master":"Slave");
	  GEM0->setAssociateDet(*GEM1);
	}
	else {
	  if (!pixGEM1)
	    CsErrLog::msg(elFatal,__FILE__,__LINE__,"CsPixelGEM \"%s\" -> CsPixelGEM association failed",TBName0.c_str());
	  if (doAmpCorrMC && // ***** CASE MONTE CARLO  AMPLITUDE CORRELATION...
	      !(pixGEM0->isMaster()^pixGEM1->isMaster()))//... MASTER XOR SLAVE
	    CsErrLog::msg(elFatal,__FILE__,__LINE__,"PixelGEM coupled planes %s and %s both declared as \"s\", check your options",TBName0.c_str(),TBName1.c_str(),pixGEM0->isMaster()?"Master":"Slave");
	  pixGEM0->setAssociateDet(*pixGEM1);
	}
	CsErrLog::msg(elInfo,__FILE__,__LINE__,"GEM association %s <- %s.",
		      TBName0.c_str(),TBName1.c_str());

      }
    }
    if ( (doAmpCorrMC || pdet1)
         && !( (GEM0 && GEM0->getAssociateDet() )
               || ( pixGEM0 && pixGEM0->getAssociateDet() ) ) )
      CsErrLog::msg(elFatal,__FILE__,__LINE__,"GEM association failed %s <- %s",
                    TBName0.c_str(),TBName1.c_str());

  }
  return;

}
