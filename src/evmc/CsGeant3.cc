// $Id: CsGeant3.cc,v 1.90 2010/12/21 15:54:01 schluter Exp $
 
/*!
   \file    CsGeant3.cc
   \brief   Geant Interface Class.
   \author  Benigno Gobbo
   \version $Revision: 1.90 $
   \date    $Date: 2010/12/21 15:54:01 $
*/

#include <cmath>

#include "CsInit.h"
#include "CsGeant3.h"
#include "CsGeom.h"
#include "CsDet.h"
#include "CsOpt.h"
#include "CsComgNtCommons.h"
#include "CsErrLog.h"
#include "CsZebraProto.h"
#include "CsHbookProto.h"
#include <algorithm>
#include <cstdlib>

#include "CsGeom.h"
#include "CsRandom.h"
#include "CsMCTrkHit.h"
#include "CsMCRICH1Hit.h"
#include "CsHelix.h"
#include "CsTriggerHodoDetector.h"
#include "CsGEMDetector.h" 
#include "CsPixelGEMDetector.h"
#include "CsCalorimeter.h"

using namespace std;
using namespace CLHEP;
using CS::DetID;

CsGeant3* CsGeant3::_instance = 0;

// This are the COMGEANT ntuples
QbeaType Qbea;
QheaType Qhea;
QkinType Qkin;
QhitType Qhit;

//      COMMON/GCLINK/JDIGI ,JDRAW ,JHEAD ,JHITS ,JKINE ,JMATE ,JPART
//     +      ,JROTM ,JRUNG ,JSET  ,JSTAK ,JGSTAT,JTMED ,JTRACK,JVERTX
//     +      ,JVOLUM,JXYZ  ,JGPAR ,JGPAR2,JSKLT
struct gclink_struct {
  int jdigi,  jdraw, jhead, jhits, jkine, jmate, jpart,
    jrotm, jrung, jset, jstak, jgstat, jtmed, jtrack, jvertx,
    jvolum, jxyx, jgpar, jgpar2, jsklt; } gclink; 

//      COMMON/COLINK/JOHIT, JNULL
struct colink_struct { int johit, jnull; } colink;

// Use old good a-la fortran syntax
static union{ int  *iq; float* q; };
int  *lq;

// ZEBRA file logical unit
int lunfz;

// HBOOK file logical unit
int lunhb = 100;

int levent;
int lentries;

extern "C" {
  extern struct { int iquest[100]; } quest_;
}

CsGeant3* CsGeant3::Instance() {
  if( _instance == 0 ) {
    _instance = new CsGeant3();
  }
  return( _instance );
}

CsGeant3::CsGeant3() {

  _TGeantInUse = false;
#if USE_TGEANT
    // ###### TGEANT CODE START ######
    _TGeantInUse = false;
    std::string path;
    if (CsOpt::Instance()->getOpt( "CsTGEANTFile", "file", path ))
    {
      _TGeantInUse = true;
        openTgeantFile(path);
    }
    
    // ###### TGEANT CODE END ######
#endif

  int* zebra = CsInit::Instance()->getZebra();
  lq = &zebra[9];
  iq = &zebra[17];

  int ixdiv = 0;
  mzlink(ixdiv,"/GCLINK/",&(gclink.jdigi),&(gclink.jsklt),&(gclink.jdigi),8);
  mzlink(ixdiv,"/COLINK/",&(colink.johit),&(colink.jnull),&(colink.johit),8);
  //mzwipe(0);

  CsInit* init = CsInit::Instance();

  bool mickeyall = init->mickeyAll();
  bool mickeyonreco = init->mickeyOnReco();

  if(!_TGeantInUse) {
    if( !mickeyall && !mickeyonreco ) {
      _MCFiles = init->getMCFilesList();     // get the list of MC files
      _MCNtFiles = init->getMCNtFilesList(); // get the list of MC Nt files
      if( !_MCNtFiles.empty() ) {  // There are only Ntuple Files...
        _currentMCNtFile = _MCNtFiles.begin();   // get the fist file in list
        // open the fist file
        int status = openGeantNTFile( *_currentMCNtFile );  
        _clearMCstructs();
        if( status != 0 ) {
          string str = "Geant input file "; 
          str.append(*(*_currentMCNtFile));
          str.append(" not opened." );
          CsErrLog::mes( elFatal, str );
        }
        else {
          setNtuple();
        }
        _NTFile = true;
      }
      else if( !_MCFiles.empty() ) { 
        _currentMCFile = _MCFiles.begin();     // get the fist file in list
        // open the fist file
        int status = openGeantFile( *_currentMCFile );  
        _clearMCstructs();
        if( status != 0 ) {
          string str = "Geant input file "; 
          str.append( *(*_currentMCFile) );
          str.append(" not opened." );
          CsErrLog::mes( elFatal, str );
        }
        _NTFile = false;
      }
      else {
        // should NEVER enter here, but in any case...
        CsErrLog::mes( elFatal, "No MC Files to process." );
      }
    } 
  }
  // some preset...
  _run   = 0;
  _event = 0;
  _mickeyfirst = true;

  _CGVersion  = 0.0;
  _GeoVersion = 0.0;

  // and some cleaning...
  _clearMCstructs();

  //-am--------------------------- 0
  // Read a trigger map
  if( !_trig.ReadTriggerMap() ) {
    CsErrLog::mes( elInfo, " No temporary trigger file in options." );
  }
  //-am--------------------------- 0

}

bool CsGeant3::getNextEvent(unsigned long selTrigMask) {

  clear(); // Clear my private lists

  bool mickeyall    = CsInit::Instance()->mickeyAll();
  bool mickeyhits   = CsInit::Instance()->mickeyHits();
  bool mickeyonreco = CsInit::Instance()->mickeyOnReco();

  if (CsInit::Instance()->resetRandomSeed()) {
    // ***** optionally: RESET THE RANDOM NUMBER GENERATOR *****
    // Can be useful for debugging purposes, so as to get simillar events in
    // various versions of coral, despite diverging random generation histories.
    CsRandom *random = CsRandom::Instance();
    static long seed; static bool first = true; if (first) {
      seed = random->getSeed(); first = false;
    }
    long newSeed = getEventNb()|seed; random->setSeed(newSeed);
    random->flagGauss();// Reset the, pair-wise, generation of gaussian numbers.
    CsErrLog::msg(elInfo,__FILE__,__LINE__,"Resetting random seed: 0x%lx^0x%x = 0x%lx\n",seed,getEventNb(),newSeed);
  }
  
#if USE_TGEANT
    // ###### TGEANT CODE START ######
  if (_TGeantInUse) {
      if (_outputBackEnd->streamGetNextEvent()) {
          readTgeantEvent();
          return true;
      } else
          return false;
  }
    // ###### TGEANT CODE END ######
#endif
  
  bool gotIt = true, trigmaskok = false; do {  // Until trigger OK...
    if (mickeyall) {
      int nevents = CsInit::Instance()->mickeyNumberOfEvents();
      if (_event<nevents) {
	readMickeyMouseHead(); readMickeyMouseKine(); readMickeyMouseHits();
	gotIt = true;
      }
      else
	gotIt = false;
    }
    else if (mickeyonreco) {
      readMickeyMouseHead(); readMickeyMouseKine(); readMickeyMouseHits();
      gotIt = true;
    } 
    else if (!_NTFile) {        // ********** ZEBRA FILE **********
      if (readGeantEvent()) {                            // ***** NEXT EVENT?...
	readGeantHead(); readGeantKine();
	if (mickeyhits) readMickeyMouseHits();
	else            readGeantHits();
	gotIt = true;
	//-am--------------------------- 1
	if (!_trig.CheckTmpTrigger()) {
	  CsErrLog::mes( elInfo, " Event not triggered (skipped)." );
	  getNextEvent();
	}
	//-am--------------------------- 1
      }
      else {                                       // ***** ...Else NEXT FILE...
	_currentMCFile++; if (_currentMCFile!=_MCFiles.end()) {
	  int status = openGeantFile(*_currentMCFile); 
	  if (status!=0) CsErrLog::msg(elFatal,__FILE__, __LINE__,
  "Error opening Geant input file \"%s\"",(*_currentMCFile)->c_str());
	  _clearMCstructs(); _NTFile = false;
	  if (readGeantEvent()) {
	    readGeantHead(); readGeantKine();
	    if (mickeyhits) readMickeyMouseHits();
	    else            readGeantHits();
	    gotIt = true;
	    //-am--------------------------- 2
	    if (!_trig.CheckTmpTrigger()) {
	      CsErrLog::mes( elInfo, " Event not triggered (skipped)." );
	      getNextEvent();
	    }
	    //-am--------------------------- 2
	  }
	  else {                           // ***** ...Else NEXT FILE is NTuple?
	    _currentMCNtFile++; if (_currentMCNtFile!=_MCNtFiles.end()) {
	      int status = openGeantNTFile( *_currentMCNtFile ); 
	      if (status!=0) CsErrLog::msg(elFatal,__FILE__, __LINE__,
  "Error opening Geant input file \"%s\"",(*_currentMCNtFile)->c_str());
	      _clearMCstructs(); _NTFile = true;
	      if (readGeantNTEvent()) {
		readGeantNTHead(); readGeantNTKine();
		if (mickeyhits) readMickeyMouseHits();
		else            readGeantHits();
		gotIt = true;
	      }
	      else
		gotIt = false;
	    }
	    else
	      gotIt = false;
	  }
	}
	else
	  gotIt = false;
      }
    }
    else {                      // ********** NTuple FILE **********
      if (readGeantNTEvent()) {
	readGeantNTHead(); readGeantNTKine();
	if (mickeyhits) readMickeyMouseHits();
	else            readGeantHits();
	gotIt = true;
	//-am--------------------------- 3
	if (!_trig.CheckTmpTrigger()) {
	  CsErrLog::mes( elInfo, " Event not triggered (skipped)." );
	  getNextEvent();
	}
	//-am--------------------------- 3
      }
      else {
	_currentMCNtFile++; if (_currentMCNtFile!=_MCNtFiles.end()) {
	  int status = openGeantNTFile( *_currentMCNtFile ); 
	  if (status!=0) CsErrLog::msg(elFatal,__FILE__, __LINE__,
  "Error opening Geant input file \"%s\"",(*_currentMCNtFile)->c_str());
	  _clearMCstructs(); _NTFile = true;
	  if (readGeantNTEvent()) {
	    readGeantNTHead(); readGeantNTKine();
	    if( mickeyhits ) readMickeyMouseHits();
	    else             readGeantHits();
	    gotIt = true;
	    //-am--------------------------- 4
	    if (!_trig.CheckTmpTrigger()) {
	      CsErrLog::mes( elInfo, " Event not triggered (skipped)." );
	      getNextEvent();
	    }
	    //-am--------------------------- 4
	  }
	  else
	    gotIt = false;
	}
	else
	  gotIt = false;
      }
    }
    if (gotIt ) {               // ********** TRIGGER SLECTION **********
      setTriggerMask(); if ((_TrigMask&selTrigMask)!=0) trigmaskok = true;
    }    
  } while (gotIt&&!trigmaskok && 
	   !mickeyall && !mickeyhits && 
	   selTrigMask!=0xffffffff);
  return gotIt;
}

string CsGeant3::getMCFileName() { 
  if( _NTFile )
    return( *(*_currentMCNtFile) ); 
  else
    return( *(*_currentMCFile) );
} 

/*! \fn int CsGeant3::openGeantFile(const string* fname, int lrecl, const char *option)
    \author R. Brun (modified by B.Gobbo)
    \date 24 April 1999
    \version 0.0
    \brief This function assumes a Geant file in Zebra/FZ format
    created with the C I/O option "L"
*/

int CsGeant3::openGeantFile(const string* fname, int lrecl, 
			     const char *option) {

  _Gname = *fname;
  int   ier;
  int   len  = fname->size();
  char* name = new char[len+1];
  strcpy( name, fname->c_str() );
  cfopen( lunfz, 0, 0, "r ", 0, name, ier, 2, len );
  if( ier == 0 ) {
    quest_.iquest[0] = lunfz;
    int lopt = strlen(option);
    fzfile( lunfz, lrecl, option, lopt );
  }
  delete [] name;
  return( ier );
}

/*! \fn int CsGeant3::openGeantNTFile(const string* fname, const int lrecl )
    \author B.Gobbo
    \date 04 August 1999
    \version 0.0
    \brief Opens an NT file
*/

int CsGeant3::openGeantNTFile(const string* fname, const int lrecl ) {

  lentries = 0;
  levent = 0;
  _Gname = *fname;
  int   status = -1;
  int   len  = fname->size();
  char* name = new char[len+1];
  strcpy( name, fname->c_str() );
  lunhb++;
  hropen( lunhb, "NTUPLE", name, " ", lrecl, status);
  if( status == 0 ) {
    hldir( " ", " ");
    hrin( 1, 999, 0 );
    hrin( 2, 999, 0 );
    hnoent( 2, lentries );
  }
  delete [] name;
  return( status );
}

void CsGeant3::clear() {

  // Remove all MC Tracks and clear the _tracks list
  if( !_tracks.empty() ) {
    list<CsMCTrack*>::iterator i;
    for( i=_tracks.begin(); i!=_tracks.end(); i++ ) {
      delete *i;
    }
    _tracks.clear();
  }

  // Remove all MC Vertices and clear the _vertices list
  if( !_vertices.empty() ) {
    list<CsMCVertex*>::iterator i;
    for( i=_vertices.begin(); i!=_vertices.end(); i++ ) {
      delete *i;
    }
    _vertices.clear();
  }

  // Remove all MC Hits and clear the _hits list
  if( !_hits.empty() ) {
    list<CsMCHit*>::iterator i;
    for( i=_hits.begin(); i!=_hits.end(); i++ ) {
      CsMCTrkHit* trkhit = dynamic_cast<CsMCTrkHit*>(*i);
      CsMCRICH1Hit* rich1hit = dynamic_cast<CsMCRICH1Hit*>(*i);
      if( trkhit != NULL ) { delete trkhit; }
      else if( rich1hit != NULL ) { delete rich1hit; }
      else { delete *i; }
    }
    _hits.clear();
  }

  // Clear all detector MC hit lists
  list<CsDetector*>::iterator id;
  list<CsDetector*> det = CsGeom::Instance()->getDetectors();
  for( id=det.begin(); id!=det.end(); id++ ) {
    (*id)->clearMCHitList();
  }
  CsRICH1Detector* rich = CsGeom::Instance()->getRich1Detector();
  if( rich != NULL ) {
    rich->clearMCHitList();
  }
}

struct CsGeant3::_sortMCHits : 
  public binary_function<CsMCHit*, CsMCHit*, bool> {
  bool operator() (CsMCHit* hh1, CsMCHit* hh2) { 
    
    CsMCTrkHit* h1 = dynamic_cast<CsMCTrkHit*>(hh1);
    CsMCTrkHit* h2 = dynamic_cast<CsMCTrkHit*>(hh2);
    if(h1&&h2) {
          if(*h1<*h2) return true;
          return false;
    }
    if((h1!=0)&&(h2==0)) return true;
    if((h1==0)&&(h2!=0)) return false;
    CsMCRICH1Hit* rh1 = dynamic_cast<CsMCRICH1Hit*>(hh1);
    CsMCRICH1Hit* rh2 = dynamic_cast<CsMCRICH1Hit*>(hh2);
    if(*rh1<*rh2) return true;
    return false;
  }
};
 
/*! \fn bool CsGeant3::readGeantEvent(const char *option)
    \author R. Brun (modified by B.Gobbo)
    \date  22 April 1999
    \brief Read the next event from the Geant FZ file.
    Search for the next HEAD record.
    If found, read VERT,KINE,HITS and DIGI records
      returns 0 if nothing read or EOF reached.
      returns 1 otherwise.
*/

bool CsGeant3::readGeantEvent(const char *option) {

  // This needs BANK name in Zebra Output. To be implemented in COMGEANT 
  const int maxhead = 100;
  int iuhead[maxhead];
  const char *H = strstr(option,"H");
  const char *D = strstr(option,"D");
  // BG 990422 Add OHIT option
  const char *O = strstr(option,"O");
  int nuhead;

  // Clear my private lists
  clear();
  
  mzwipe(0);

  // 990429. Some problems: endianess and char swap. OHIT are
  // correclty swaped. HEAD, VERT, KINE and DIGI are not. 
  // Simple solution: accept both swapped or not names...

  //       Get next event header
  while(1) {
    nuhead = maxhead;
    fzin(lunfz,0,0,0,"S",nuhead,iuhead,1);
    if (quest_.iquest[0] >= 2) return( false );
    if (nuhead >= 3) {
      if (strncmp((char*)&iuhead[nuhead-1],"HEAD",4) == 0 ||
	  strncmp((char*)&iuhead[nuhead-1],"DAEH",4) == 0 ) break;
    }
  }
 
  //        header found, read all events data structures
  fzin(lunfz,0,&(gclink.jhead),1,"A",nuhead,iuhead,1);
  while(1) {
    nuhead = maxhead;
    fzin(lunfz,0,0,0,"S",nuhead,iuhead,1);
    if (quest_.iquest[0] >= 2) return( false );
    if (strncmp((char*)&iuhead[nuhead-1],"VERT",4)  == 0 ||
	strncmp((char*)&iuhead[nuhead-1],"TREV",4)  == 0 ) {
      fzin(lunfz,0,&(gclink.jvertx),1,"A",nuhead,iuhead,1);
    } else if (strncmp((char*)&iuhead[nuhead-1],"KINE",4) == 0 ||
	       strncmp((char*)&iuhead[nuhead-1],"ENIK",4)  == 0 ) {
      fzin(lunfz,0,&(gclink.jkine),1,"A",nuhead,iuhead,1);
      if (!H && !D && !O ) return( true );
    } else if (strncmp((char*)&iuhead[nuhead-1],"HITS",4) == 0 ||
	       strncmp((char*)&iuhead[nuhead-1],"STIH",4)  == 0 ) {
      fzin(lunfz,0,&(gclink.jhits),1,"A",nuhead,iuhead,1);
      if (!D && !O ) return( true );
    } else if (strncmp((char*)&iuhead[nuhead-1],"DIGI",4) == 0 ||
	       strncmp((char*)&iuhead[nuhead-1],"IGID",4)  == 0 ) {
      fzin(lunfz,0,&(gclink.jdigi),1,"A",nuhead,iuhead,1);
      // Start: BG 990422 Add OHIT readout
      if( !O ) return( true );      
    } else if (strncmp((char*)&iuhead[nuhead-1],"OHIT",4) == 0 ||
	       strncmp((char*)&iuhead[nuhead-1],"TIHO",4)  == 0 ) {
      fzin(lunfz,0,&(colink.johit),1,"A",nuhead,iuhead,1);
      // End
      return( true );
    }
  }
}

/*! \fn bool CsGeant3::readGeantNTEvent()
    \author B.Gobbo
    \date   04 August 1999
    \brief Read the next event from the Geant Ntuple file.
      returns 0 if nothing read or EOF reached.
      returns 1 otherwise.
*/
bool CsGeant3::readGeantNTEvent() {

  // Clear my private lists
  clear();

  int status;
  levent++;

  if( levent <= lentries ) {
    hgnt( 2, levent, status );
  }
  else {
    return( false );
  }
  if( status == 0 ) 
    return( true );
  else
    return( false );
}

/*! \fn int CsGeant3::readGeantHead()
    \author B.Gobbo
    \date 09 September 1999
    \brief  Fill the informations from Geant3 HEAD bank.
*/
int CsGeant3::readGeantHead() {

  int jhead   = gclink.jhead;
  if (jhead  <= 0) return( 0 );
  _run   = iq[jhead+1];
  _event = iq[jhead+2];

  // check comgeant versions 
  string CGverFromDet = CsGeom::Instance()->getComgeantVers();
  int ver = atoi( CGverFromDet.c_str() + 1 );
  int rel = atoi( CGverFromDet.c_str() + 7 );
  string CGgeoFromDet = "";
  if( ver > 6 || ( ver==6 && rel>2 ) ) {
    CGgeoFromDet = CsGeom::Instance()->getGeomVers();
  }

  static string oldMCFile = "";
  static bool first = true;
  string thisMCFile = getMCFileName();
  if( oldMCFile != thisMCFile ) {
    oldMCFile = thisMCFile;
    first = true;
  }

  if( first ) {
    char CGver[9];
    char CGgeo[9];
    strncpy( CGver, (char*)&iq[jhead+7], 8 );
    CGver[8] = '\0';
    strncpy( CGgeo, (char*)&iq[jhead+9], 8 );
    CGgeo[8] = '\0';
    first = false;

    // check correctness (just of
    if( CGver[0] == 'v'                    &&
	'0' <= CGver[1] && CGver[1] <= '9' &&
	'0' <= CGver[2] && CGver[2] <= '9' &&
	'0' <= CGver[3] && CGver[3] <= '9' &&
	strncmp( CGver+4, "rel", 3 ) == 0  &&
	'0' <= CGver[7] && CGver[7] <= '9' ) {

      cout << "COMGEANT version from Zebra file:    " << CGver << endl 
	   << "COMGEANT version from detectors.dat: " << CGverFromDet << endl
	   << "Geometry version from Zebra file:    " << CGgeo << endl
	   << "Geometry version from detectors.dat: " << CGgeoFromDet << endl;
    
      _CGVersion  = float(atoi(CGver+1))+float(atoi(CGver+7))/10.0;
      _GeoVersion = float(atoi(CGgeo+1))+float(atoi(CGgeo+6))/100.0;
	

      if( CGver != CGverFromDet ) {
	CsErrLog::mes( elFatal, "detectors.dat, Zebra file version missmatch" );
      }
    }
    else {
      cout << "No way to perform detectors.dat <-> zebra file consistency check on this release." << endl;
    }
  }

  return( 1 );
}

bool CsGeant3::readGeantLund() {  

  int jhead = gclink.jhead;
  if( jhead <= 0 ) return false;

  int jhaux = lq[jhead-1];
  if( jhaux <= 0 ) return false;

  int jtlnd = lq[jhaux-1];
  if( jtlnd <= 0 ) return false;
  int tlndsize = iq[jtlnd-1];

  // cleanage...
  _halfclearMCstructs();

  int np = iq[jtlnd-2];

  if( tlndsize == 27 ) {  // Old COMGEANT TLND banks size 

    _useludatanew = false;
    _tlndok = true;

    _ludata.x  = q[jtlnd+1];
    _ludata.y  = q[jtlnd+2];
    _ludata.w2 = q[jtlnd+3];
    _ludata.q2 = q[jtlnd+4];
    _ludata.u  = q[jtlnd+5];
    for(int i= 6; i<=17; i++ ) _ludata.lst[i+14] = int(q[jtlnd+i]);
    for(int i=18; i<=27; i++ ) _ludata.parl[i+2] = q[jtlnd+i];

    int jplnd;
    lujet part;
    for( int i=1; i<=np; i++ ) {
      jplnd = lq[jtlnd-i];
      for( int j=1; j<=5; j++ ) {
	part.k[j-1] = int(q[jplnd+j]);
      }
    
      part.p[0] = q[jplnd+7]; //
      part.p[1] = q[jplnd+8]; // NB: Rotation from COMG ref. to Coral ref.
      part.p[2] = q[jplnd+6]; //
      part.p[3] = q[jplnd+9];
      part.p[4] = q[jplnd+10];
    
      // at the moment there's no V vector on data 
      //for( int j=11; j<=15; j++ ) {
      //  part.v[j-11] = q[jplnd+j];
      //}
      for( int i=0; i<5; i++ ) part.v[i] = 0.0;
      //part.lu2kine = int(q[jplnd+16]);
      part.lu2kine = int(q[jplnd+11]);

      _lujets.push_back( part );
      
    }

    int jrlnd = lq[jtlnd];
    if( jrlnd != 0 ) { 
      int rlndsize = iq[jrlnd-1]; 

      if( rlndsize == 56 || rlndsize == 79 ) { // Check on RLND bank size 
	_rlndok = true;
	_ludata.genType = int(q[jrlnd+1]);
	for(int i= 2; i<=15; i++ ) _ludata.cut[i-2]   = q[jrlnd+i];
	for(int i=16; i<=35; i++ ) _ludata.lst[i-16]  = int(q[jrlnd+i]);
	for(int i=36; i<=37; i++ ) _ludata.lst[i-3]   = int(q[jrlnd+i]);
	for(int i=38; i<=46; i++ ) _ludata.parl[i-38] = q[jrlnd+i];
	for(int i=47; i<=56; i++ ) _ludata.parl[i-37] = q[jrlnd+i];
	if( rlndsize == 79 ) {
	  for(int i=57; i<=62; i++ ) _ludata.parhfl[i-57] = q[jrlnd+i];
	  for(int i=63; i<=70; i++ ) _ludata.cuthfl[i-63] = q[jrlnd+i];
	  for(int i=71; i<=74; i++ ) _ludata.lsthfl[i-71] = int(q[jrlnd+i]);
	  for(int i=75; i<=79; i++ ) _ludata.lsthfl[i-70] = int(q[jrlnd+i]);
	}
      }
      else {
	CsErrLog::msg( elFatal,__FILE__, __LINE__,
		       "Wrong RLND bank size: %d, expected: 56 or 79", 
		       rlndsize );
      }
    }
  }

  else if( tlndsize == 48 ) { // New COMGEANT TLND banks size (Lepto/Aroma)

    _useludatanew = true;
    _tlndok = true;

    _ludatanew.x  = q[jtlnd+2];
    _ludatanew.y  = q[jtlnd+3];
    _ludatanew.w2 = q[jtlnd+4];
    _ludatanew.q2 = q[jtlnd+5];
    _ludatanew.u  = q[jtlnd+6];
    for(int i=7; i<=26; i++ )  _ludatanew.uservar[i-7] = q[jtlnd+i];
    for(int i=27; i<=38; i++ ) _ludatanew.lst[i-7]     = int(q[jtlnd+i]);
    for(int i=39; i<=48; i++ ) _ludatanew.parl[i-19]   = q[jtlnd+i];

    int jplnd;
    lujet part;
    for( int i=1; i<=np; i++ ) {
      jplnd = lq[jtlnd-i];
      for( int j=1; j<=5; j++ ) {
	part.k[j-1] = int(q[jplnd+j]);
      }
      part.p[0] = q[jplnd+7]; //
      part.p[1] = q[jplnd+8]; // NB: Rotation from COMG ref. to Coral ref.
      part.p[2] = q[jplnd+6]; //
      part.p[3] = q[jplnd+9];
      part.p[4] = q[jplnd+10];
      part.lu2kine = int(q[jplnd+11]);
      _lujets.push_back( part );
    }

    int jrlnd = lq[jtlnd];
    if( jrlnd != 0 ) { 
      int rlndsize = iq[jrlnd-1]; 

      if( rlndsize == 63 || rlndsize == 91 ) { // Check on RLND bank size 
	_rlndok = true;
	_ludatanew.genType = int(q[jrlnd+1]);
	for(int i=2; i<=15; i++ )  _ludatanew.cut[i-2]   = q[jrlnd+i];
	for(int i=16; i<=35; i++ ) _ludatanew.lst[i-16]  = int(q[jrlnd+i]);
	for(int i=36; i<=43; i++ ) _ludatanew.lst[i-4]   = int(q[jrlnd+i]);
	for(int i=44; i<=63; i++ ) _ludatanew.parl[i-44] = q[jrlnd+i];
	if( rlndsize == 91 ) {
	  for( int i=64; i<=73; i++ ) _ludatanew.parhfl[i-64] = q[jrlnd+i]; 
	  for( int i=74; i<=81; i++ ) _ludatanew.cuthfl[i-74] = q[jrlnd+i]; 
	  for( int i=82; i<=91; i++ ) _ludatanew.lsthfl[i-82] = int(q[jrlnd+i]); 
	}
      }
      else {
	CsErrLog::msg( elFatal,__FILE__, __LINE__,
		       "Wrong RLND bank size: %d, expected: 63 or 91", 
		       rlndsize ); 
      }
    }
  }

  else if( tlndsize == 826 ) { // New COMGEANT TLND banks size (Pythia)

    _useludatanew = true;
    _tlndok = true;

    _ludatanew.x  = q[jtlnd+2];
    _ludatanew.y  = q[jtlnd+3];
    _ludatanew.w2 = q[jtlnd+4];
    _ludatanew.q2 = q[jtlnd+5];
    _ludatanew.u  = q[jtlnd+6];
    for(int i=  7; i<= 26; i++ ) _ludatanew.uservar[i-7] = q[jtlnd+i];
    for(int i= 27; i<=226; i++ ) _pypars.mstp[i-27]      = int(q[jtlnd+i]); 
    for(int i=227; i<=426; i++ ) _pypars.parp[i-227]     = q[jtlnd+i]; 
    for(int i=427; i<=626; i++ ) _pypars.msti[i-427]     = int(q[jtlnd+i]); 
    for(int i=627; i<=826; i++ ) _pypars.pari[i-627]     = q[jtlnd+i]; 

    int jplnd;
    lujet part;
    for( int i=1; i<=np; i++ ) {
      jplnd = lq[jtlnd-i];
      for( int j=1; j<=5; j++ ) {
	part.k[j-1] = int(q[jplnd+j]);
      }
      part.p[0] = q[jplnd+7]; //
      part.p[1] = q[jplnd+8]; // NB: Rotation from COMG ref. to Coral ref.
      part.p[2] = q[jplnd+6]; //
      part.p[3] = q[jplnd+9];
      part.p[4] = q[jplnd+10];
      part.lu2kine = int(q[jplnd+11]);
      _lujets.push_back( part );
    }

    int jrlnd = lq[jtlnd];
    if( jrlnd != 0 ) { 
      int rlndsize = iq[jrlnd-1]; 

      if( rlndsize == 863 ) { // Check on RLND bank size 
	_rlndok = true;
	_ludatanew.genType = int(q[jrlnd+1]);
	_pysubs.msel       = int(q[jrlnd+2]);
	_pysubs.mselpd     = int(q[jrlnd+3]);
	for(int i=  4; i<=503; i++ ) _pysubs.msub[i-4]   = int(q[jrlnd+i]);
	for(int i=504; i<=663; i++ ) 
	  _pysubs.kfin[(i-504)/80][(i-504)%80] = int(q[jrlnd+i]); 
	for(int i=664; i<=863; i++ ) _pysubs.ckin[i-664] = q[jrlnd+i];
      }
      else {
	CsErrLog::msg( elFatal,__FILE__, __LINE__,
		       "Wrong RLND bank size: %d, expected: 863", 
		       rlndsize );    
      }
    }
  }
  else {
    CsErrLog::msg( elFatal,__FILE__, __LINE__,
		   "Wrong TLND bank size: %d, accepted: 27, 28, 826", 
		       tlndsize );      
  }
    
  return true;
}

/*! \fn int CsGeant3::readGeantNTHead()
    \author B.Gobbo
    \date 09 September 1999
    \brief  Fill the informations from Header block of Comgeant Ntuples.
*/
int CsGeant3::readGeantNTHead() {
  _run   = Qhea.irun;
  _event = Qhea.ieve;
  return( 1 );
}

/*! \fn int CsGeant3::readMickeyMouseHead()
    \author B.Gobbo
    \date 22 August 2000
    \brief  Fill the informations of Mickey Mouse MC.
*/
int CsGeant3::readMickeyMouseHead() {
  _run    = 69;
  _event += 1;
  return( 1 );
}

/*! \fn int CsGeant3::readGeantKine()
    \author B.Gobbo
    \date 20 April 1999
    \brief Fill MC Track and Vertex classes from Geant3 KINE and VERT banks
*/
int CsGeant3::readGeantKine() {

  readGeantLund();

  int jkine   = gclink.jkine;
  int jvertx  = gclink.jvertx;
  if (jkine  <= 0) return( 0 );
  if (jvertx <= 0) return( 0 );
 
  double x, y, z, t;
  int nvertx = iq[jvertx+1];   // number of MC vertices
  int iv, jv;
  //CsMCVertex** vtxPtr = new CsMCVertex*[nvertx]; 
  map< int, CsMCVertex*, less<int> > mv;
  // Loop over vertices...
  for( iv=1; iv<=nvertx; iv++ ) {
    jv = lq[jvertx-iv]; if( jv<=0 ) continue;
    x  = q[jv+1] * 10.;                   // vertex coordinates (mm)
    y  = q[jv+2] * 10.;
    z  = q[jv+3] * 10.;
    CsGeant3::_geaRef2CsRefVec( x, y, z ); // Rotate fr. GeantRS to CompassMRS
    t  = q[jv+4];                         // Time of flight (hope in us)
    // add a vertex (with no inTrack) to the list;
    // inTrack will be set later
    CsMCVertex* vertex = new CsMCVertex( iv, x, y, z, t ); 
    _vertices.push_back( vertex ); 
    //vtxPtr[iv-1] = vertex;
    mv[iv-1] = vertex;
  }

  double px, py, pz; int ip;
  int ntrack = iq[jkine+1];   // number of MC tracks
  int it, jk;
  //CsMCTrack** trkPtr = new CsMCTrack*[ntrack]; 
  map< int, CsMCTrack*, less<int> > mt;
  // Loop over tracks...
  for( it=1; it<=ntrack; it++ ) {
    jk = lq[jkine-it]; if( jk<=0 ) continue;
    ip = int(q[jk+5]);   // particle nr in JPART
    px = q[jk+1];        // Particle momentum (GeV) 
    py = q[jk+2];        //
    pz = q[jk+3];        //
    CsGeant3::_geaRef2CsRefVec( px, py, pz );

    int ovx = int(q[jk+6]); // origin vertex

    //CsMCTrack* track = new CsMCTrack( it, px, py, pz, CsMCParticle( ip ),
    //				      *vtxPtr[ovx-1] );
    CsMCTrack* track = new CsMCTrack( it, px, py, pz, CsMCParticle( ip ),
				      *mv[ovx-1] );
    _tracks.push_back( track );
    //trkPtr[it-1] = track;
    mt[it-1] = track;

    //vtxPtr[ovx-1]->addOutTrack( *track ); 
    mv[ovx-1]->addOutTrack( *track );  // add out track to vertex

  }
    
  // Loop again on tracks to set:
  //   - vertices inTrack
  //   - tracks outTrack list
  for( it=1; it<=ntrack; it++ ) {
    jk = lq[jkine-it]; if( jk<=0 ) continue;
    int nev = int(q[jk+7]);      // n end vertices
    if( nev > 0 ) {
      for( int i=1; i<=nev; i++ ) {
	iv = int(q[jk+7+i]); 
	//vtxPtr[iv-1]->setInTrack( *trkPtr[it-1] ); 
	mv[iv-1]->setInTrack( *mt[it-1] ); // Set inTracks in vertices
	jv = lq[jvertx-iv]; if( jv<=0 ) continue;
	int nok = int(q[jv+7]);      // n out tracks
	//trkPtr[it-1]->addOutVertex( *vtxPtr[iv-1] ); 
	mt[it-1]->addOutVertex( *mv[iv-1] ); //Add outVertex in track
	for( int j=1; j<=nok; j++ ) {
	  int ot = int(q[jv+7+j]);    
	  //trkPtr[it-1]->addOutTrack( *trkPtr[ot-1] ); 
	  mt[it-1]->addOutTrack( *mt[ot-1] ); 
	}
      }
    }
  }

  //delete [] vtxPtr;
  //delete [] trkPtr;
  return ntrack;
}

/*! \fn int CsGeant3::readGeantNTKine()
    \author B.Gobbo
    \date 04 August 1999
    \brief Fill MC Track and Vertex classes from Geant3 Ntuples
*/
int CsGeant3::readGeantNTKine() {

  double x, y, z, t; int ivp;
  int nvertx = Qkin.nver;   // number of MC vertices
  //CsMCVertex** vtxPtr = new CsMCVertex*[nvertx]; 
  map< int, CsMCVertex*, less<int> > mv;
  // Loop over vertices...
  for( int iv=0; iv<nvertx; iv++ ) {
    ivp = Qkin.igev[iv];                   // Geant Vertex Number
    x   = Qkin.vert[iv][0] * 10.;          // vertex coordinates (mm)
    y   = Qkin.vert[iv][1] * 10.;
    z   = Qkin.vert[iv][2] * 10.;
    CsGeant3::_geaRef2CsRefVec( x, y, z ); // Rotate fr. GeantRS to CompassMRS
    t  = double(Qkin.ltimv[iv]) / 10000.; // Time of flight (us)
    // add a vertex (with no inTrack) to the list;
    // inTrack will be set later
    CsMCVertex* vertex = new CsMCVertex( ivp, x, y, z, t ); 
    _vertices.push_back( vertex );
    //vtxPtr[iv] = vertex;
    mv[iv] = vertex;
  }

  double px, py, pz; int itp, ip;
  int ntrack = Qkin.ntra;   // number of MC tracks
  //CsMCTrack** trkPtr = new CsMCTrack*[ntrack]; 
  map< int, CsMCTrack*, less<int> > mt;
  // Loop over tracks...
  for( int it=0; it<ntrack; it++ ) {
    ip  = Qkin.iget[it];    // Geant Track Number
    itp = Qkin.itra[it];    // particle nr in JPART
    px  = Qkin.ptra[it][0]; // Particle momentum (GeV) 
    py  = Qkin.ptra[it][1];
    pz  = Qkin.ptra[it][2];
    CsGeant3::_geaRef2CsRefVec( px, py, pz );

    int ovx = Qkin.itvb[it]; // origin vertex

    //CsMCTrack* track = new CsMCTrack( ip, px, py, pz, CsMCParticle( itp ),
    //		  		        *vtxPtr[ovx-1] );
    CsMCTrack* track = new CsMCTrack( ip, px, py, pz, CsMCParticle( itp ),
    		  		        *mv[ovx-1] );
    _tracks.push_back( track );
    //trkPtr[it] = track;
    mt[it] = track;
  }

  for( int iv=0; iv<nvertx; iv++ ) {
    int intrk = Qkin.imov[iv];  // origin track
    int nok = Qkin.ntdv[iv];    // n out tracks
    int ot1 = Qkin.itdv[iv];    // 1st out track
    int otn = ot1 + nok - 1;    // last out track

    for( int i=ot1; i<=otn; i++ ) {
      //vtxPtr[iv]->addOutTrack( *trkPtr[i-1] ); 
      mv[iv]->addOutTrack( *mt[i-1] ); // Add outTrack in vertex
    }

    if( intrk > 0 ) {
      //vtxPtr[iv]->setInTrack( *trkPtr[intrk-1] ); 
      mv[iv]->setInTrack( *mt[intrk-1] );   // Set inTracks in vertices
      //trkPtr[intrk-1]->addOutVertex( *vtxPtr[iv] ); 
      mt[intrk-1]->addOutVertex( *mv[iv] ); //Add outVertex in track
      for( int i=ot1; i<=otn; i++ ) {
	//trkPtr[intrk-1]->addOutTrack( *trkPtr[i-1] ); 
	mt[intrk-1]->addOutTrack( *mt[i-1] ); // Add outTrack in track
      }
    }
  }


  //delete [] vtxPtr;
  //delete [] trkPtr;

  return ntrack;
}

/*! \fn int CsGeant3::readMickeyMouseKine()
    \author B.Gobbo
    \date 
    \brief generates MC Track and Vertex objects (this is a porting of
     Sergei Gerassimov's Core Simulator).
*/
int CsGeant3::readMickeyMouseKine() {

  int ntracks = 0;
 
  if( _mickeyfirst ) {

    _mickeyfirst = false;

    // Vertices and Tracks from Options...
    CsOpt* opt = CsOpt::Instance();

    string tag = "Mickey";
    string key = "vertex";
    vector<float> vtxdata;
    while( opt->getOptRec( tag, key, vtxdata ) ) { 
      _mickeyvertices.push_back( vtxdata );
    }
    key = "track";
    vector<float> trkdata;
    while( opt->getOptRec( tag, key, trkdata ) ) { 
      _mickeytracks.push_back( trkdata );
    }  
  }

  vector<int> nvxused;
  map< int, CsMCVertex*, less<int> > mv;
  for( unsigned int i=0; i<_mickeyvertices.size(); i++ ) {
    vector<float> vtxdata = _mickeyvertices[i];
    int   nv   = int(vtxdata[0]); 
    if( find( nvxused.begin(), nvxused.end(), nv ) != nvxused.end() ) {
      CsErrLog::mes( elFatal, "More vertices with same number." );
    }
    nvxused.push_back( nv );
    double x   = vtxdata[1];   // vertex coordinates (mm)
    double y   = vtxdata[2];
    double z   = vtxdata[3];
    double t   = vtxdata[4];   // Time of flight (hope in us)
    // add a vertex (with no inTrack) to the list;
    // inTrack will be set later
    CsMCVertex* vertex = new CsMCVertex( nv, x, y, z, t ); 
    _vertices.push_back( vertex ); 
    mv[nv-1] = vertex;
    vtxdata.clear();
  }
  
  if( _vertices.size() == 0 ) {
    CsErrLog::mes( elFatal, "No vertices found in option file." );
  }
    
  vector<int> ntkused;
  map< int, CsMCTrack*, less<int> > mt;
  for( unsigned int i=0; i<_mickeytracks.size(); i++ ) {
    vector<float> trkdata = _mickeytracks[i];
    ntracks ++;
    int    nt  = int(trkdata[0]);
    if( find( ntkused.begin(), ntkused.end(), nt ) != ntkused.end() ) {
      CsErrLog::mes( elFatal, "More tracks with same number." );
    }
    ntkused.push_back( nt );
    double px  = trkdata[1];        // Particle momentum (GeV) 
    double py  = trkdata[2];        
    double pz  = trkdata[3];        
    int ip     = int(trkdata[4]);   // particle type
    int invx   = int(trkdata[5]);   // origin vertex
    int outvx  = int(trkdata[6]);   // end vertex
    // check existence of in and out vertex:
    if( find( nvxused.begin(), nvxused.end(), invx ) == nvxused.end() ) {
      CsErrLog::mes( elFatal, "Track with associated unexisting in vertex." );
    }
    if( outvx != 0 ) {
      if( find( nvxused.begin(), nvxused.end(), outvx ) == nvxused.end() ) {
	CsErrLog::mes( elFatal, 
		       "Track with associated unexisting out vertex." );
      }
    }
    
    CsMCTrack* track = new CsMCTrack( nt, px, py, pz, CsMCParticle( ip ),
				      *mv[invx-1] );
    _tracks.push_back( track );
    mt[nt-1] = track;
    if( outvx > 0 ) {
      mv[outvx-1]->setInTrack( *mt[nt-1] );    //set inTracks in vertex
      mt[nt-1]->addOutVertex( *mv[outvx-1] );  //add outVertex in track
    }
    mv[invx-1]->addOutTrack( *mt[nt-1] );    //add outTrack to vertex 
  }
  
  if( _tracks.size() == 0 ) {
    CsErrLog::mes( elFatal, "No tracks found in option file." );
  }
    
  list<CsMCTrack*>::iterator It;
  for( It=_tracks.begin(); It!=_tracks.end(); It++ ) {
    list<CsMCVertex*> outvertices = (*It)->getOutVertices();
    if( !outvertices.empty() ) {
      CsMCVertex* outvtx = outvertices.front();
      list<CsMCTrack*> vtxTrks = outvtx->getOutTracks();
      list<CsMCTrack*>::iterator Itt;
      for( Itt=vtxTrks.begin(); Itt!=vtxTrks.end(); Itt++ ) {
	(*It)->addOutTrack( *(*Itt) ); 
      }
    }
  }

  return ntracks;

}

/*! \fn int CsGeant3::readGeantHits()
    \author B.Gobbo
    \date 21 April 1999
    \brief Fill MC Hit classe from Geant3 OHIT/HITS bank
*/

int  CsGeant3::readGeantHits() {

  // At present (990422) uses OHIT bank only
  // Implement JHIT in future? Will see...

  int johit = colink.johit;
  if( johit <= 0 ) return 0;

  int savehits = 0;

  int ntraject = iq[johit+1];     // number of trajectories

  for( int it=1; it<=ntraject; it++ ) {
    int johtj = lq[johit-it];
    int track = iq[johtj+1];      // Geant Track number
    int nhits = iq[johtj+2];      // Number of hits


    list<CsMCTrack*>::iterator It;  // find out my CsMCTrack...
    for( It=_tracks.begin(); 
	 It!=_tracks.end() && (*It)->getGnum()!=track; It++ ); 

    for( int ih=1; ih<=nhits; ih++ ) {
      //      savehits ++;
      int johth = lq[johtj-ih];

      int ndat=iq[johth-1]  ; 

      // number of words per hit 
      // 19 for 'det' and 'sla' - hodoscopes will appear soon
      // 4 for calorimeters
      // 20 for RICH

      if(ndat==19) {  // trackers (can be hodoscopes also - if 41 < detc < 50 )

	double xmm = q[johth+1] * 10.; // Centre of Tr. in sens vol (MRS) (mm)
	double ymm = q[johth+2] * 10.;
	double zmm = q[johth+3] * 10.;
	CsGeant3::_geaRef2CsRefVec( xmm, ymm, zmm );
	double uid = q[johth+4] * 10.; // Entrance point in sens vol (DRS) (mm)
	double vid = q[johth+5] * 10.;
	double wid = q[johth+6] * 10.;
	CsGeant3::_geaRef2CsRefVec( uid, vid, wid );
	double uod = q[johth+7] * 10.; // Exit point from sens vol (DRS) (mm)
	double vod = q[johth+8] * 10.;
	double wod = q[johth+9] * 10.;
	CsGeant3::_geaRef2CsRefVec( uod, vod, wod );
	double elos = q[johth+10];       // total energy lost 
	double eion = q[johth+11];       // energy lost by ionisation
	double p   = q[johth+12];        // momentum (GeV)
	int    tk0 = int(q[johth+13]);   // 0=this, N=Other Geant Track  
	double time = q[johth+14];       // DeltaT from time zero (ns?)
	int    detc = int(q[johth+15]);  // detector code
	int    detn = int(q[johth+16]);  // detector number
	double cxm = q[johth+17];        // cos of trk at entrance point (MRS)
	double cym = q[johth+18];
	double czm = q[johth+19];
	CsGeant3::_geaRef2CsRefVec( cxm, cym, czm );

	if (detc==626 || detc==627) {

	  //   ********** RECOIL PROTON DETECTOR **********
	  // Note: Geant hits from RPD have associated detector number == 0!
	  //      => Therefore they will skip the piece of code infra where
	  //        hits from all other tracking detectors are assigned to
	  //        their associated detector.
	  // This is most probably due to a bug in COMGeant. In order to work
	  // around this bug, the hit assignment is performed here in a distinct
	  // block.

	  if( tk0 != 0 ) tk0 = CsMCParticle( tk0 ).getNumber();
	  CsGeom* geom = CsGeom::Instance();
	  list<CsDetector*> dets = geom->getDetectors();
	  list<CsDetector*>::iterator Id; 
	  for( Id=dets.begin(); Id!=dets.end(); Id++ ){ // Find associated detector. 
	    DetID did = (*Id)->GetID();
	    if(int(did.GetNumber()) == detc){ // found
	      // cout<<"============= MC Hit for det "<<(*Id)->GetTBName()<<" will be added"<<endl;
	      break;
	    }
	  };
	  if (Id==dets.end()) {
	    cout<<"CsGeant3 ==> Corresponding RP detector do not exists"<<endl;
	    continue;
	  }
	  
	  double detEff = 0.;
	  detEff = (*Id)->getEff();
	  CsMCHit* hit = new CsMCTrkHit( xmm, ymm, zmm, 
					 uid, vid, wid, 
					 uod, vod, wod,
					 elos, eion, time,
					 Hep3Vector( p*cxm, p*cym, p*czm ),
					 *(*It), tk0, *(*Id), detc);
	  _hits.push_back( hit );
	  (*Id)->addMCHit( *hit );
	}//  end of "recoil proton detector" block

	if (detn!=0) {

	  // ********** GEANT HIT HAS ASSOCIATED TRACKING DETECTOR **********

	  // if tk0!=0 set its particle number using PDG convention
	  if( tk0 != 0 ) tk0 = CsMCParticle( tk0 ).getNumber();
	  
	  list<CsDetector*> dets = CsGeom::Instance()->getDetectors();
	  list<CsDetector*>::iterator Id; // Find my associated detector... 
	  for (Id = dets.begin(); Id!=dets.end() && !(*Id)->IsMyHit(detn); Id++);
	  if (Id==dets.end()) {
	    CsErrLog::msg (elError,__FILE__,__LINE__,
 "Hit associated to non existing detector %d, please check detectors.dat",detn);
	    continue;
	  }
	  CsDetector *csDet = *Id;
	  
	  //            ********** CsPixelGEM (CsPG) ->CsGEM  **********
          CsPixelGEMDetector *pixelGEM = dynamic_cast<CsPixelGEMDetector*>(csDet);
#define CsG3_DISPATCH_pixelGEM
#ifdef CsG3_DISPATCH_pixelGEM
	  // This patch is needed to process 2007.05_1 and 2007.05_2 data, where
	  // a single ID covers both the pixelised and stripped parts.
	  if (pixelGEM) {
	    if (pixelGEM->GetTBName()[4]=='P') { // ***** CsPG OF THE PIXEL KIND
	      if (pixelGEM->getAssociateDet()) {
		if (!csDet->inActiveArea(xmm,ymm)) { // If CsPG is associated...
		  // ...(to a CsPG of the strip kind): redirect its hits to the
		  // latter when they fall outside its sensitive area.
		  csDet = pixelGEM->getAssociateDet();
		  detn = csDet->GetID().GetNumber();
		  pixelGEM = dynamic_cast<CsPixelGEMDetector*>(csDet);
		}
		else // ...in order to prevent any further association 
		  pixelGEM = 0;
	      }
	    }
	    else {                               // ***** CsPG OF THE STRIP KIND
	      if (!pixelGEM->doAmpCorrelationMC()) {
		// If the amplitude correlation scheme is not applied, we have
		// to discard hits that fall into the pixel region
		if (!csDet->inActiveArea(xmm,ymm)) continue;
	      }
	      //else we must be dealing w/ a slave and hits wil be discarded
	      //infra
	    }
	  }
#endif

	  //              ********** INEFFICIENCY **********
	  double detEff = csDet->getEff(); bool fired;
	  if (detEff<.9999 && CsRandom::flat()>=detEff) fired = false;
	  else                                          fired = true;

	  CsDetector *csDet2 = 0; // ***** CORRELATION in GEMs *****
	  // To ensure that hits from GEM (CsGEM or CsPG) associated planes be
	  // identical which is necessary for correct amplitude correlation:
	  // - Save only the hits from "Master" (cf. CsGEM/PG's option "Master")
	  //  detectors, but feed them to both master and slave planes.
	  // - Skip those from "slave" ones.
	  // - Perform time smearing here (and not in the CsGEM/PG class).
	  // - As to the inefficiency: it should also be correlated. Not yet
	  //  done. (Could be done by assigning a relative inefficiency, to
	  //  the slave detector, which would be very low. Or by assuming
	  //  this relative inefficiency to be 0%: then would suffice to
	  //  cancel the condition infra checking for "random2<detEff2".)
	  // Maybe also usefull for SI (in future).
	  CsGEMDetector *GEMd1 = dynamic_cast<CsGEMDetector*>(csDet);
	  if (pixelGEM || GEMd1) {
	    bool doAmpCorr, isMaster; float tRes;
	    if (GEMd1) {
	      doAmpCorr = GEMd1->doAmpCorrelationMC();
	      isMaster = doAmpCorr && GEMd1->isMaster();
	      tRes = GEMd1->getTRes();    csDet2 = GEMd1->getAssociateDet();
	    }
	    else {
	      doAmpCorr = pixelGEM->doAmpCorrelationMC();
	      isMaster = doAmpCorr && pixelGEM->isMaster();
	      tRes = pixelGEM->getTRes(); csDet2 = pixelGEM->getAssociateDet();
	    }
	    if (doAmpCorr) {
	      if (isMaster) {
		time = time + tRes * CsRandom::gauss();
		if (csDet2) {
		  double detEff2 = csDet2->getEff();
		  double random2 = CsRandom::flat();
		  // The "if" condition below could be cancelled. On the ground
		  // that the main contribution to the inefficiency of the
		  // GEM is via the triggering of the avalanche, which affects
		  // equally the 2 coordinate planes...
		  if (random2>detEff2)
		    csDet2 = 0;
		}
		else
		  CsErrLog::msg(elFatal,__FILE__,__LINE__,
  "%s: has no associate while \"ampCorrelationMC\" requested",csDet->GetTBName().c_str());
	      }
	      else {
		fired = false; csDet2 = 0; // Skip hits from slave
	      }
	    }
	    else csDet2 = 0;
	  }

	  if (fired) {
	    savehits ++;
	    CsMCHit* hit = new CsMCTrkHit( xmm, ymm, zmm, 
					   uid, vid, wid, 
					   uod, vod, wod,
					   elos, eion, time,
					   Hep3Vector( p*cxm, p*cym, p*czm ),
					   *(*It), tk0, *csDet, detn);
	    _hits.push_back( hit );
	    csDet->addMCHit( *hit );
	  }

      	  if (csDet2) { 
	    // Case of GEM amplidtude correlation: fill associated plane
	    savehits ++;
	    int detn2 = csDet2->GetID().GetNumber();
	    CsMCHit* hit2  = new CsMCTrkHit( xmm, ymm, zmm, 
					     uid, vid, wid, 
					     uod, vod, wod,
					     elos, eion, time,
					     Hep3Vector( p*cxm, p*cym, p*czm ),
					     *(*It), tk0, *csDet2, detn2);
	    _hits.push_back( hit2 );
	     csDet2->addMCHit( *hit2 );
	  }
	}
      }
      else if (ndat==4) {

	// ******************** CALORIMETERS ********************

	CsCalorimeter::CalorimeterMCData d;
	d.dE  = q[johth+1];            // dE (GeV)

        d.track_id = int(q[johth+2]);  // 0=original track entered cell ,
	                               // N=Other Geant Track 
	                               //(N=particle type)

	//double  time_c  = q[johth+3];  // t-t0 (ns)
	d.dT       = q[johth+3];       // t-t0 (ns)
	d.cell_id = int(q[johth+4]);   //cell ID = matrix ID + module number
	                               // I can put these two numbers 
	                               // separetly, if needed
        CsDet::AddMCHitAll(ndat,&d);
      }
      else if(ndat==20){
//
// RICH data 
//
	double cher_xm=q[johth+1]*10.; // X of phot. detect. point MRS [mm]
	double cher_ym=q[johth+2]*10.; // Y
	double cher_zm=q[johth+3]*10.; // Z
	CsGeant3::_geaRef2CsRefVec(cher_xm, cher_ym, cher_zm );
	
	double cher_yd=q[johth+4]*10.; // Y photon detection point, DRS
	double cher_zd=q[johth+5]*10.; // Z 
	
	double cher_xp=q[johth+6]*10.; // X point of production MRS [mm]
	double cher_yp=q[johth+7]*10.; // Y
	double cher_zp=q[johth+8]*10.; // Z 
	CsGeant3::_geaRef2CsRefVec(cher_xp, cher_yp, cher_zp );
	
	double cher_mpx=q[johth+9];  // Px of mother particle
	double cher_mpy=q[johth+10]; // Py
	double cher_mpz=q[johth+11]; // Pz
	CsGeant3::_geaRef2CsRefVec(cher_mpx, cher_mpy, cher_mpz );
	
	double cher_eph=q[johth+12]; // photon energy (eV)
	
	int cher_im=int(q[johth+13]); // if generated by product = IPART
	
	double cher_tim=q[johth+14]; // t-t0 [ns], should be around zero for particles
	                             // from the main vertex, non-zero for pile-up
	
	int cher_ityp=int(q[johth+15]); // detector type * 1000
	int cher_id  =int(q[johth+16]); // detector ID = 900+cathode number
//	int detn = (cher_id/100)*100;   // detector ID
	int detn = 900;   // detector ID
	int cher_cathode  =int(q[johth+16])-detn; // cathode number
	
	double cher_xr=q[johth+17]*10.; // X point of reflection MRS [cm]
	double cher_yr=q[johth+18]*10.; // Y
	double cher_zr=q[johth+19]*10.; // Z
	CsGeant3::_geaRef2CsRefVec(cher_xr, cher_yr, cher_zr );
	
	double cher_ang=q[johth+20]; // cher angle, [rad]

	savehits ++;
	// set particle number using PDG convention
	if( cher_im != 0 ) cher_im = CsMCParticle( cher_im ).getNumber();

	CsRICH1Detector* rich = CsGeom::Instance()->getRich1Detector();
	if( rich == NULL ) {
	  CsErrLog::mes( elError, 
	   "Hit associated to a non esisting detector, please check detectors.dat content" );
	  continue;
	}
	
	CsMCHit* hit = new CsMCRICH1Hit( 
					cher_xm, cher_ym, cher_zm, 
					cher_yd, cher_zd, 
					cher_xp, cher_yp, cher_zp,
					cher_xr, cher_yr, cher_zr,
					cher_eph, cher_tim, cher_ang,
					Hep3Vector( cher_mpx, cher_mpy, 
						    cher_mpz ),
					*(*It),
					cher_im, cher_cathode, *rich );
	_hits.push_back( hit );
	rich->addMCHit( *hit ); 
      }
    }
  }

  // sort hits
  _hits.sort( _sortMCHits() ); 

  list<CsMCHit*>::iterator Ih;
  for( Ih=_hits.begin(); Ih!=_hits.end(); Ih++ ) {
    CsMCTrack* mytrack = (*Ih)->getMCTrack();
    mytrack->addMCHit( *(*Ih) );   // add the hit to its track	       
  }

  return savehits;
}

/*! \fn int CsGeant3::readGeantNTHits()
    \author B.Gobbo
    \date 04 August 1999
    \brief Fill MC Hit classe from Comgeant Ntuples
*/

int  CsGeant3::readGeantNTHits() {

  int CGvers = atoi( (CsGeom::Instance()->getGeomVers()).c_str() + 1 );

  int savehits = 0;

  int nhits = Qhit.nhit;      // Number of hits
  
  for( int ih=0; ih<nhits; ih++ ) {

    int track  = (Qhit.ip1hit[ih] & 0x0000ffff);  // Ntuple Track Number
    int gtrack = Qkin.iget[track-1];              // Geant Track Number

    list<CsMCTrack*>::iterator It;  // find out my CsMCTrack...
    for( It=_tracks.begin(); 
	 It!=_tracks.end() && (*It)->getGnum()!=gtrack; It++ ); 

    savehits ++;

    double xmm = 0.; // Will be set later on 
    double ymm = Qhit.hit[ih][0]*10.; // Centre of Tr. in sens. vol. (MRS) (mm)
    double zmm = Qhit.hit[ih][1]*10.;
    CsGeant3::_geaRef2CsRefVec( xmm, ymm, zmm );
    double uid = 0.; // Entrance point in sens. vol. (DRS) (mm). NOT AVAILABLE!
    double vid = 0.;
    double wid = 0.;
    // CsGeant3::_geaRef2CsRefVec( uid, vid, wid );
    double uod = 0.; // Exit point from sens. vol. (DRS) (mm). NOT AVAILABLE!
    double vod = 0.;
    double wod = 0.;
    //CsGeant3::_geaRef2CsRefVec( uod, vod, wod );
    double elos = 0.;     // total energy lost. NOT AVAILABLE!
    double eion = 0.;     // energy lost by ionisation. NOT AVAILABLE!
    double p    = 0.;     // momentum (GeV). NOT AVAILABLE!
    int    tk0 = (Qhit.ip1hit[ih] & 0xffff0000)>>16; // 0=this, N=Geant Trk  
    // DeltaT from T0 (ns)
    double time = double(((Qhit.ip2hit[ih]&0xffff0000)>>16)-32768)/10.;  
    int    detc = 0;  // detector code  NOT AVAILABLE
    int    detn = (Qhit.ip2hit[ih] & 0x0000ffff);  // detector number
    double cxm = 0.;  // cos of trk at entrance point (MRS). NOT AVAILABLE!
    double cym = 0.;
    double czm = 0.;
    //CsGeant3::_geaRef2CsRefVec( cxm, cym, czm );

    // if tk0!=0 set its particle number using PDG convention
    if( tk0 != 0  ) tk0 = CsMCParticle( tk0 ).getNumber();

    //if( detc > 1000 && detn != 0 ) { // if associated to a detector...
    if( detn != 0 ) { // if associated to a detector...
      CsGeom* geom = CsGeom::Instance();
      list<CsDetector*> dets = geom->getDetectors();
      list<CsDetector*>::iterator Id; // Find my associated detector... 
      if( CGvers < 5 ) {
	// This was COMGEANT 4:
	// A stupid thing: ntuple does not contain the detector id but
	// the row number in detector.dat... So, temporarly use this...
	for( Id=dets.begin(); Id!=dets.end() && (*Id)->getRow()!=detn; Id++ );
	if( Id == dets.end() ) {
	  CsErrLog::mes( elError, 
			 "Hit associated to a non esisting detector, please check detectors.dat content" );
	  continue;
	}
      }
      else {
	// On COMGEANT 5 it works correctly:
	for( Id=dets.begin(); Id!=dets.end() && (*Id)->GetID()!=detn; Id++ );
	if( Id == dets.end() ) {
	  CsErrLog::mes( elError, 
			 "Hit associated to a non esisting detector, please check detectors.dat content" );
	  continue;
	}
      }
      zmm = (*Id)->getZcm();  // Set 3rd coordinate from Detector Centre (mm)

      CsMCHit* hit = new CsMCTrkHit( xmm, ymm, zmm, 
				     uid, vid, wid, 
				     uod, vod, wod,
				     elos, eion, time,
				     Hep3Vector( p*cxm, p*cym, p*czm ),
				     *(*It), tk0, *(*Id), detn);  
      _hits.push_back( hit );
      (*Id)->addMCHit( *hit );
    }
  }

  // sort hits
  _hits.sort( _sortMCHits() ); 

  list<CsMCHit*>::iterator Ih;
  for( Ih=_hits.begin(); Ih!=_hits.end(); Ih++ ) {
    CsMCTrack* mytrack = (*Ih)->getMCTrack();
    mytrack->addMCHit( *(*Ih) );   // add the hit to its track	       
  }

  return savehits;
}

/*! \fn int CsGeant3::readMickeyMouseHits()
    \author B.Gobbo
    \date 21 August 2000
    \brief Generate MC Hit objects (this is a porting of
     Sergei Gerassimov's Core Simulator).
*/

int  CsGeant3::readMickeyMouseHits() {

  int savehits = 0;

  // No multiple scattering ???
  string tag = "Mickey";
  string key = "multiple scattering off";
  bool multiplescattering = ! CsOpt::Instance()->getOpt( tag, key ); 


  list<CsDetector*> dets = CsGeom::Instance()->getDetectors();

  list<CsMCTrack*>::iterator It;
  It=_tracks.begin(); It++;  // skip 1st track (beam)
  for( ; It!=_tracks.end(); It++ ) {
    double px = (*It)->getPX();
    double py = (*It)->getPY();
    double pz = (*It)->getPZ();
    double p  = sqrt( px*px + py*py + pz*pz );

    // Instance an helix...
    double x    = ((*It)->getInVertex())->getX();
    double y    = ((*It)->getInVertex())->getY();
    double z    = ((*It)->getInVertex())->getZ();
    double dxdz = px/pz;
    double dydz = py/pz;
    double cop  = double(((*It)->getParticle())->getCharge()) / p;
    double cov[15] = { 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1 };

    CsHelix hi( x, y, z, dxdz, dydz, cop, cov );

    list<CsDetector*>::iterator Id;
    for( Id=dets.begin(); Id!=dets.end(); Id++ ) {
      
      // temporary protection...
      if( (*Id)->GetID() == 900 ) continue;

      double detx = (*Id)->getXcm(); // Detector centre coordinates (MRS)
      double dety = (*Id)->getYcm();
      double detz = (*Id)->getZcm();
      if( z > detz ) continue;  // Track starts after this detector...

      if( !((*It)->getOutVertices()).empty() ) {
	double zend = ((*It)->getOutVertices()).front()->getZ();
	if( zend < detz ) continue; // Track died before this detector...
      }

      CsHelix he;

      if( !hi.Extrapolate(detz,he) ) continue; // Ignore if something wrong...

      x    = he.getX();
      y    = he.getY();
      z    = he.getZ();
      dxdz = he.getDXDZ();
      dydz = he.getDYDZ();
      cop  = he.getCop();
      p    = fabs( 1. / cop );
      double cosx   = dxdz / sqrt( 1. + dxdz*dxdz + dydz*dydz );
      double cosy   = dydz / sqrt( 1. + dxdz*dxdz + dydz*dydz );
      double cosz   = 1.   / sqrt( 1. + dxdz*dxdz + dydz*dydz );

      // Is the extrapolated helix inside the detector? 
      int err;
      HepMatrix irotd = ((*Id)->getRotDRS()).inverse( err );

      double x_drs = irotd(1,1)*x + irotd(1,2)*y + irotd(1,3)*z - detx;
      double y_drs = irotd(2,1)*x + irotd(2,2)*y + irotd(2,3)*z - dety;
      double z_drs = irotd(3,1)*x + irotd(3,2)*y + irotd(3,3)*z - detz;
      
      double detXsiz = (*Id)->getXsiz();
      double detYsiz = (*Id)->getYsiz();
      double detZsiz = (*Id)->getZsiz();

      // in any case, move from hi to he:
      hi = he;

      // Skip if out of detector volume
      if( x_drs > 0 && ( detXsiz/2 - x_drs ) < 0 ) continue;  
      if( x_drs < 0 && (-detXsiz/2 - x_drs ) > 0 ) continue;
      if( y_drs > 0 && ( detYsiz/2 - y_drs ) < 0 ) continue; 
      if( y_drs < 0 && (-detYsiz/2 - y_drs ) > 0 ) continue;

      double cosx_drs = irotd(1,1)*cosx + irotd(1,2)*cosy + irotd(1,3)*cosz;
      double cosy_drs = irotd(2,1)*cosx + irotd(2,2)*cosy + irotd(2,3)*cosz;
      double cosz_drs = irotd(3,1)*cosx + irotd(3,2)*cosy + irotd(3,3)*cosz;

      if( multiplescattering ) {
	// Multiple scattering (completely taken from Sergei's software)
	double path = detZsiz / cosz_drs;
	double len = path / (*Id)->getRdLen(); 
	double sigTheta = 0.0136 * fabs(cop) * sqrt(len) * (1.+0.038*log(len));

	double z1 = CsRandom::gauss();
	double z2 = CsRandom::gauss();
	double dr = z1*path*sigTheta/sqrt(12.0)+z2*path*sigTheta/2.;
	double da = z2*sigTheta;
	x = x + dr;
	dxdz = tan( atan(dxdz) + da );
      
	z1 = CsRandom::gauss();
	z2 = CsRandom::gauss();
	dr = z1*path*sigTheta/sqrt(12.0)+z2*path*sigTheta/2.;
	da = z2*sigTheta;
	y = y + dr;
	dydz = tan( atan(dydz) + da );

	hi = CsHelix( x, y, z, dxdz, dydz, cop, cov );
      }

      // OK: track hits detector...

      // Centre of Tr. in sens. vol. (MRS) (mm)
      double xmm = x; 
      double ymm = y;
      double zmm = z;
      // Entrance point in sens. vol (DRS) (mm)
      double uid = x_drs + cosx_drs / cosz_drs * ( -detZsiz/2 - z_drs );
      double vid = y_drs + cosy_drs / cosz_drs * ( -detZsiz/2 - z_drs );
      double wid = - detZsiz/2;
      // Exit point from sens. vol. (DRS) (mm)
      double uod = x_drs + cosx_drs / cosz_drs * (  detZsiz/2 - z_drs );
      double vod = y_drs + cosy_drs / cosz_drs * (  detZsiz/2 - z_drs );
      double wod =   detZsiz/2;

      double elos = 0;       // total energy lost 
      double eion = elos;    // energy lost by ionisation
      int    tk0  = 0;       // 0=this, N=Other Geant Track  
      double time = 0;       // DeltaT from time zero (ns?)

      savehits ++;
	  
      CsMCHit* hit = new CsMCTrkHit( xmm, ymm, zmm, 
				     uid, vid, wid, 
				     uod, vod, wod,
				     elos, eion, time,
				     Hep3Vector( p*cosx, p*cosy, p*cosz ),
				     *(*It), tk0, *(*Id)); 
      _hits.push_back( hit );
      (*Id)->addMCHit( *hit );
    }	
  }

  // sort hits
  _hits.sort( _sortMCHits() ); 

  list<CsMCHit*>::iterator Ih;
  for( Ih=_hits.begin(); Ih!=_hits.end(); Ih++ ) {
    CsMCTrack* mytrack = (*Ih)->getMCTrack();
    mytrack->addMCHit( *(*Ih) );   // add the hit to its track	       
  }

  return savehits;
}

void CsGeant3::setNtuple() {

  // Well, no comments needed...

  hbname( 2, " ", 0, "$CLEAR");

  hbname( 2, "RUN", &Qhea.ieve, "$SET:ieve" );
  hbname( 2, "RUN", &Qhea.irun, "$SET:irun" );
  hbname( 2, "RUN", &Qhea.iend, "$SET:iend" );

  hbname( 2, "BEAM", &Qbea.ibtyp, "$SET:ibtyp" );
  hbname( 2, "BEAM", &Qbea.ibfla, "$SET:ibfla" );
  hbname( 2, "BEAM", (int*)&Qbea.bpara, "$SET:bpara" );

  hbname( 2, "KINE", &Qkin.nver, "$SET:nver" );
  hbname( 2, "KINE", (int*)&Qkin.igev, "$SET:igev" );
  hbname( 2, "KINE", (int*)&Qkin.vert, "$SET:vert" );
  hbname( 2, "KINE", (int*)&Qkin.ltimv, "$SET:ltimv" );
  hbname( 2, "KINE", (int*)&Qkin.imov, "$SET:imov" );
  hbname( 2, "KINE", (int*)&Qkin.ntdv, "$SET:ntdv" );
  hbname( 2, "KINE", (int*)&Qkin.itdv, "$SET:itdv" );
  hbname( 2, "KINE", &Qkin.ntra, "$SET:ntra" );
  hbname( 2, "KINE", (int*)&Qkin.iget, "$SET:iget" );
  hbname( 2, "KINE", (int*)&Qkin.ptra, "$SET:ptra" );
  hbname( 2, "KINE", (int*)&Qkin.itra, "$SET:itra" );
  hbname( 2, "KINE", (int*)&Qkin.itvb, "$SET:itvb" );
  hbname( 2, "KINE", (int*)&Qkin.itve, "$SET:itve" );
  hbname( 2, "KINE", &Qkin.nkinc, "$SET:nkinc" );
  hbname( 2, "KINE", (int*)&Qkin.xkinc, "$SET:xkinc" );
  hbname( 2, "KINE", (int*)&Qkin.pkinc, "$SET:pkinc" );

  hbname( 2, "hit", &Qhit.nhit, "$SET:nhit" );
  hbname( 2, "hit", &Qhit.nhitall, "$SET:nhitall" );
  hbname( 2, "hit", (int*)&Qhit.ip1hit, "$SET:ip1hit" );
  hbname( 2, "hit", (int*)&Qhit.ip2hit, "$SET:ip2hit" );
  hbname( 2, "hit", (int*)&Qhit.hit, "$SET:hit" );

  hbname( 2, "dig", &Qhit.ndig, "$SET:ndig" );
  hbname( 2, "dig", &Qhit.ndigall, "$SET:ndigall" );
  hbname( 2, "dig", (int*)&Qhit.ip1dig, "$SET:ip1dig" );
  hbname( 2, "dig", (int*)&Qhit.ip2dig, "$SET:ip2dig" );
  hbname( 2, "dig", &Qhit.npdig, "$SET:npdig" );
  hbname( 2, "dig", (int*)&Qhit.jpdig, "$SET:jpdig" );
}
 
void CsGeant3::_geaRef2CsRefVec( double& x, double& y, double& z ) {
  // From GeantRS : X // beam, Z vertical
  // To MRS :       Z // beam, Y verical
  double tmp = x;
  x          = y;
  y          = z;
  z          = tmp;
}

void CsGeant3::_geaRef2CsRefMat( HepMatrix& a ) {
  // From GeantRS : X // beam, Z vertical
  // To MRS :       Z // beam, Y verical
  HepMatrix r(3,3);
  double set[] = { 0, 1, 0, 0, 0, 1, 1, 0, 0 };
  for( int i=0; i<9; i++ ) r( i/3+1, i%3+1 ) = set[i];
  a = r * a;
  a = a * r.T();
}


void CsGeant3::setTriggerMask( void ) {
  
  // themporary solution. Names of hodoscopes should be optional!
  
  _TrigMask = 0;
  bool FI4(1),FL4(1),FM4(1),FO3(1); // against double hits in one hodo
  int  NHI(0),NHL(0),NHM(0),NHO(0); // number of hits 
  
  // V.Alexakhin 05.10.03 trigger mask from GEANT (if available)
  // Would be also nice to add in proper place the checking of mu-prime candidate 
  // if it satisfies trigger matrices, like it is done for real data... 
  //  
  int jhead   = gclink.jhead;
  if (jhead  > 0) {    
    int tmaskMC = iq[jhead+5];     // trigger mask from GEANT 
    if( tmaskMC >= 0 ) {
        _TrigMask = tmaskMC ;
	//        cout<<" tmaskMC set "<<  _TrigMask<<endl;
       	return;                 
    }   
  }

  const list<CsMCTrack*> &mctracks = CsGeant3::Instance()->getMCTracks();
  list<CsMCTrack*>::const_iterator it;
  for(it=mctracks.begin();it!=mctracks.end();it++) {
    if( (*it)->getGnum() != 2 ) continue;
    const list<CsMCHit*> &hits = (*it)->getMCHits();
    list<CsMCHit*>::const_iterator ih;
    for(ih=hits.begin();ih!=hits.end();ih++) {
      const CsDet* det = (*ih)->getDet();
      const string& name = det->GetTBName();
      
           if( FI4 && name.find("HI04") == 0 ) { NHI++; FI4 = false; }
      else if( FL4 && name.find("HL04") == 0 ) { NHL++; FL4 = false; }
      else if( FM4 && name.find("HM04") == 0 ) { NHM++; FM4 = false; }
      else if( FO3 && name.find("HO03") == 0 ) { NHO++; FO3 = false; }
      
      else if( name.find("HI05") == 0 ) NHI++;
      else if( name.find("HL05") == 0 ) NHL++;
      else if( name.find("HM05") == 0 ) NHM++;
      else if( name.find("HO04") == 0 ) NHO++;

    }
  }
  
  if( NHI > 1 ) { _TrigMask = _TrigMask|1; }
  if( NHL > 1 ) { _TrigMask = _TrigMask|4; }
  if( NHM > 1 ) { _TrigMask = _TrigMask|2; }
  if( NHO > 1 ) { _TrigMask = _TrigMask|8; }
  
  //cout<<"TrigMask="<<_TrigMask<<"  NHI="<<NHI<<" NHL="<<NHL<<" NHM="<<NHM<<" NHO="<<NHO<<endl<<endl;
}

void CsGeant3::_clearMCstructs() {

  //flags
  _useludatanew = true;
  _tlndok       = false;
  _rlndok       = false;

  //old ludata
  _ludata.x       = 0.0;
  _ludata.y       = 0.0;
  _ludata.w2      = 0.0;
  _ludata.q2      = 0.0;
  _ludata.u       = 0.0;
  for( int i=0; i<35; i++ ) _ludata.lst[i]     =   0;
  for( int i=0; i<30; i++ ) _ludata.parl[i]    = 0.0;
  for( int i=0; i< 6; i++ ) _ludata.parhfl[i]  = 0.0;
  for( int i=0; i< 8; i++ ) _ludata.cuthfl[i]  = 0.0;
  for( int i=0; i<10; i++ ) _ludata.lsthfl[i]  =   0;
  _ludata.genType = 0;
  for( int i=0; i<14; i++ ) _ludata.cut[i]     = 0.0;

  //new ludata (ludatanew)
  _ludatanew.x       = 0.0;
  _ludatanew.y       = 0.0;
  _ludatanew.w2      = 0.0;
  _ludatanew.q2      = 0.0;
  _ludatanew.u       = 0.0;
  for( int i=0; i<20; i++ ) _ludatanew.uservar[i] = 0.0;
  for( int i=0; i<40; i++ ) _ludatanew.lst[i]     =   0;
  for( int i=0; i<30; i++ ) _ludatanew.parl[i]    = 0.0;
  for( int i=0; i<10; i++ ) _ludatanew.parhfl[i]  = 0.0;
  for( int i=0; i< 8; i++ ) _ludatanew.cuthfl[i]  = 0.0;
  for( int i=0; i<10; i++ ) _ludatanew.lsthfl[i]  =   0;
  _ludatanew.genType = 0;
  for( int i=0; i<14; i++ ) _ludatanew.cut[i]     = 0.0;

  // clear lujets vector
  _lujets.clear();
  
  // pysubs structure
  _pysubs.msel   = 0;
  _pysubs.mselpd = 0;
  for( int i=0; i<500; i++ ) _pysubs.msub[i] = 0;
  for( int j=0; j<2; j++ ) for( int i=0; i<80; i++ ) _pysubs.kfin[j][i] = 0;
  for( int i=0; i<200; i++ ) _pysubs.ckin[i] = 0.0;
  
  // pypars structure
  for( int i=0; i<200; i++ ) {
    _pypars.mstp[i] = 0; _pypars.parp[i] = 0.0; 
    _pypars.msti[i] = 0; _pypars.pari[i] = 0.0;
  }
}

void CsGeant3::_halfclearMCstructs() {

  //flags
  _useludatanew = true;
  _tlndok       = false;

  //old ludata
  _ludata.x       = 0.0;
  _ludata.y       = 0.0;
  _ludata.w2      = 0.0;
  _ludata.q2      = 0.0;
  _ludata.u       = 0.0;
  for( int i=20; i<32; i++ ) _ludata.lst[i]     =   0;
  for( int i=20; i<30; i++ ) _ludata.parl[i]    = 0.0;
 
  //new ludata (ludatanew)
  _ludatanew.x       = 0.0;
  _ludatanew.y       = 0.0;
  _ludatanew.w2      = 0.0;
  _ludatanew.q2      = 0.0;
  _ludatanew.u       = 0.0;
  for( int i= 0; i<20; i++ ) _ludatanew.uservar[i] = 0.0;
  for( int i=20; i<32; i++ ) _ludatanew.lst[i]     =   0;
  for( int i=20; i<30; i++ ) _ludatanew.parl[i]    = 0.0;

  // clear lujets vector
  _lujets.clear();
  
  // pypars structure
  for( int i=0; i<200; i++ ) {
    _pypars.mstp[i] = 0; _pypars.parp[i] = 0.0; 
    _pypars.msti[i] = 0; _pypars.pari[i] = 0.0;
  }
}


#if USE_TGEANT
// ###### TGEANT CODE START ######

void CsGeant3::openTgeantFile(std::string _fileName)
{
    // some preset...
    _run = 0;
    _event = 0;
    _mickeyfirst = true;
    _NTFile = false;

    _CGVersion = 42.0;
    _GeoVersion = 0.0;

    // and some cleaning...
    _clearMCstructs();

    if (!_trig.ReadTriggerMap()) {
        CsErrLog::mes( elInfo, " No temporary trigger file in options." );
    }

    // TGEANT code
    std::string tgeantFile = _fileName;
#ifdef TGEANT_ROOT
    if (tgeantFile.find(".root") != std::string::npos)
        _outputBackEnd = new T4OutputROOT();
    else 
#endif      
      if (tgeantFile.find(".tgeant") != std::string::npos)
        _outputBackEnd = new T4OutputASCII();
    else
        std::cerr
                << "Error in openTgentFile: Unknown TGEANT file! Use ROOT (.root) or ASCII (.tgeant) format."
                << std::endl;

    
    _outputBackEnd->streamLoad(tgeantFile);
    std::cout << "Eventsize of loaded TGEANT file: " << _outputBackEnd->streamGetEventNumber()
              << std::endl;

    int runstart = tgeantFile.find_last_of("_run");
    int runstop = tgeantFile.find_last_of(".");

    _run = 1; 
}

void CsGeant3::readTgeantEvent(void)
{
    currentEvent = _outputBackEnd->streamGetEventPointer();
    
    zBeamStart = currentEvent->beamData.trajectories.at(0).position[2];
    
    readTgeantLund();
    readTgeantKine();
    readTgeantHits();

    // at end of readTgeantEvent() prepare for next CORAL event loop
    _event++;
}

void CsGeant3::readTgeantLund(void)
{
    T4BeamData* beamData = &currentEvent->beamData;
    T4BeamParticle* currentParticle;

    _CGVersion = 42.0;
    
    lujet part;
    // cleanage...
    _halfclearMCstructs();
    
    // check what generator was used
    // resembled by this int
#define LEPTO_TGEANT   2
#define HEPGEN_TGEANT  3
#define PYTHIA_TGEANT  1
#define PRIMAKOFF_TGEANT  6

    _useludatanew = true;
    _tlndok = true;

    _ludatanew.x = beamData->x_bj;
    _ludatanew.y = beamData->y;
    _ludatanew.w2 = beamData->w2;
    _ludatanew.q2 = beamData->q2;
    _ludatanew.u = beamData->nu;

    if (sizeof(float) != sizeof(Float_t) || sizeof(int) != sizeof(Int_t))
        std::cout
                << "memcpy in next line is going to fail miserably!!! DO NOT USE RESULTS!!! (CsGeant3::readTgeantLund)"
                << std::endl;

    //iterate over all the particles
    for (unsigned int a = 0; a < beamData->nBeamParticle; a++) {
        currentParticle = &(beamData->beamParticles.at(a));
        memcpy(part.p, currentParticle->p, sizeof(Float_t) * 5);
        memcpy(part.k, currentParticle->k, sizeof(Float_t) * 5);
        memset(part.v, 0 , sizeof(Float_t) * 5);
        //this means we cant find a CsMCTrack corresponding to this
        //this is bad, but will hopefully works
        part.lu2kine = -1;
        _lujets.push_back(part);
    }

    if (beamData->generator == HEPGEN_TGEANT
            || beamData->generator == LEPTO_TGEANT) {
        memcpy(_ludatanew.cut, beamData->cut, sizeof(Float_t) * 14);
        memcpy(_ludatanew.lst, beamData->lst, sizeof(Int_t) * 40);
        memcpy(_ludatanew.parl, beamData->parl, sizeof(Float_t) * 30);
        memcpy(_ludatanew.uservar, beamData->uservar, sizeof(Float_t) *20);
    
          
        _ludatanew.genType = 5;
    }

    //this is for lepto-style generators, like lepto or hepgen
    else if (beamData->generator == PYTHIA_TGEANT){ //otherwise its pythia here we dont need to check, because there is no other possibility
        memcpy(_ludatanew.uservar, beamData->uservar, sizeof(Float_t) *20);
	memcpy(_pypars.mstp, beamData->pypars.mstp, sizeof(Int_t) *200);
        memcpy(_pypars.msti, beamData->pypars.msti, sizeof(Int_t) *200);
        memcpy(_pypars.parp, beamData->pypars.parp, sizeof(Float_t) *200);
        memcpy(_pypars.pari, beamData->pypars.pari, sizeof(Float_t) *200);
    }
    else if (beamData->generator == PRIMAKOFF_TGEANT) {
    }
    else
        //should never be here!
        std::cout << "WRONG GENERATOR CHOSEN IN TGEANT: "
                  << beamData->generator << std::endl;
}

int CsGeant3::getIndexOfVertex(double t, double x, double y, double z)
{
   for (unsigned int i = 0; i < _tVertices.size(); i++) {
    if (fabs(_tVertices.at(i)->getT() - t) < 0.0002)
        if (fabs(_tVertices.at(i)->getX() - x) < 0.0002)
            if (fabs(_tVertices.at(i)->getY() - y)< 0.0002)
                if (fabs(_tVertices.at(i)->getZ() - z) < 0.0002)
                     return i;
   }
   return -1;
}

void CsGeant3::readTgeantKine(void)
{
    vector<int> badShit;

    _tTracks.clear();
    _tVertices.clear();
    _tTrackId.clear();

    //base we add a vertex for each track-start if it does not exist yet
    //first we check if a vertex exists for the given kinematics
    //if so we use it, if not we create one
    T4BeamData* beamData = &currentEvent->beamData;

    for (unsigned int i = 0; i < beamData->nTrajectories; i++) { // cut out all tracks with (daughter-)ids that are large
		   // CORAL doesn't know particles with Id > 1000000000
		if (beamData->trajectories.at(i).particleId > 1000000000
			|| std::find(badShit.begin(), badShit.end(),
			beamData->trajectories.at(i).parentId) != badShit.end()) {
		  badShit.push_back(beamData->trajectories.at(i).trackId);
		  continue;
		}

		double x = beamData->trajectories.at(i).position[0]; /*mm*/
		double y = beamData->trajectories.at(i).position[1]; /*mm*/
		double z = beamData->trajectories.at(i).position[2]; /*mm*/
		double time = beamData->trajectories.at(i).time * 1e-9; /*ns*/

		int index = getIndexOfVertex(time, x, y, z);
		bool vertexnew = false;

		//create new vertex
		if (index == -1) { // make a new vertex
			CsMCVertex* vertex;
			if (beamData->trajectories.at(i).trackId != 1) {
				// make a new vertex for all incoming particles except primary ones
				vertex = new CsMCVertex(_tVertices.size() + 1, x, y, z, time);
//				cout << "new vert ("<<x<<", "<<y<<", "<<z<<") t="<<time<<" add="<<vertex<<endl;
			} else if (beamData->trajectories.at(i).trackId == 1) { // incoming muon vertex with vertexTime, cause incoming muon does not propagate backwards
				x = beamData->vertexPosition[0]; /*mm*/
				y = beamData->vertexPosition[1]; /*mm*/
				z = beamData->vertexPosition[2]; /*mm*/
				time = beamData->vertexTime * 1e-9; /*ns*/
				vertex = new CsMCVertex(_tVertices.size() + 1, x, y, z, time);
//				cout << "new vert for inc mu ("<<x<<", "<<y<<", "<<z<<") t="<<time<<" add="<<vertex<<endl;
			} else
				std::cout << " Should not see this in " << __LINE__ << " in file " << __FILE__ << std::endl;
  
			index = _tVertices.size(); // enumerate the index upwards
			_tVertices.push_back(vertex); // save vertex in the list
			vertexnew = true;
		}

		//skipping pre-vertex-delta electrons for now.
		//we in every case need a new track for this thing, so we make one and link it to either the existing or the newly created vertex
		CsMCTrack* track;
		if (beamData->trajectories.at(i).trackId != 1) { // do it for all tracks except bad-tracks and primary
			track = new CsMCTrack(_tTracks.size() + 1, // set the number equal to the ladder track-number
								  beamData->trajectories.at(i).momentum[0] / 1000./*GeV/c*/,
								  beamData->trajectories.at(i).momentum[1] / 1000./*GeV/c*/,
								  beamData->trajectories.at(i).momentum[2] / 1000./*GeV/c*/,
								  CsMCParticle((int) beamData->trajectories.at(i).particleId, 0),
								  *_tVertices.at(index)); // link track to vertex
		}
		else {
		  if (beamData->beamParticles.size() == 0)
		  {
		    return;
		  }
			track = new CsMCTrack(_tTracks.size() + 1, // track for incoming particle with negative momentum and anti-particle
								  -beamData->beamParticles.at(0).p[0] /*GeV/c*/,
								  -beamData->beamParticles.at(0).p[1] /*GeV/c*/,
								  -beamData->beamParticles.at(0).p[2] /*GeV/c*/,
								  CsMCParticle((int) beamData->beamParticles.at(0).k[1], 0),
								  *_tVertices.at(0)); // link this track to vertex number 0 - primary
		}
		_tTrackId[beamData->trajectories.at(i).trackId] = _tTracks.size(); // track map (tgeant Id - vector position coral)
		_tTracks.push_back(track); // track saved

		// until now:
		// all tracks looped, except the bad ones
		// vertices added if not added already (index = position in _tVertices vector)
		// alle tracks added, incoming muon with negative momentum
		// and track linked with _tVertices.at(index)
		// tgeant trackId linked with _tTracks vektor-position in map

		// we make sure everything gets linked together correctly
		_tVertices.at(index)->addOutTrack(*track); // add track to vertex
		//this is a bit of a hack for geant4 delivering the beam myon with parentId 0
		if (beamData->trajectories.at(i).parentId > 0) {
			if (index != 0 && vertexnew) {
				_tVertices.at(index)->setInTrack(
					*_tTracks.at(_tTrackId[beamData->trajectories.at(i).parentId]));
			}
			if (vertexnew) {
				_tTracks.at(_tTrackId[beamData->trajectories.at(i).parentId])
				->addOutVertex(*_tVertices.at(index));
			}
		}
	}

    //now we copy over our generated vectors into the lists of Coral
    for (unsigned int i = 0; i < _tVertices.size(); i++)
        _vertices.push_back(_tVertices.at(i));
    for (unsigned int i = 0; i < _tTracks.size(); i++)
        _tracks.push_back(_tTracks.at(i));  
}

void CsGeant3::readTgeantHits(void)
{
    // TGEANT units:
    // positions in mm
    // time in ns
    // energy in MeV
    T4Event* currentEvent = _outputBackEnd->streamGetEventPointer();
    vector<T4HitData>::iterator hIterator;

    // first add the tracking hits
    for (hIterator = currentEvent->tracking.begin(); hIterator != currentEvent->tracking.end(); hIterator++)
        createTgeantTrackingHit(*hIterator);

    // add the trigger hits
    for (hIterator = currentEvent->trigger.begin(); hIterator != currentEvent->trigger.end(); hIterator++)
        createTgeantTrackingHit(*hIterator, true);

    _hits.sort(_sortMCHits());

    // now add the calorimeter hits
    for (hIterator = currentEvent->calorimeter.begin(); hIterator != currentEvent->calorimeter.end(); hIterator++)
        createTgeantCaloHit(*hIterator);
    
    
    
//     vector<T4RichData>::iterator RICHerator;
    
//     for (RICHerator = currentEvent->rich.begin(); RICHerator != currentEvent->rich.end(); RICHerator++)
//         createTgeantRICHHit(*RICHerator);
    
    
    

    _TrigMask = currentEvent->trigMask;
}

void CsGeant3::createTgeantTrackingHit(T4HitData& hitData, bool isTrigger)
{
  double time = calcTgeantTime(hitData.time, hitData.primaryHitPosition);

  CsGeom* geom = CsGeom::Instance();
  list<CsDetector*> dets = geom->getDetectors(); // detector list of CORAL
  list<CsDetector*>::iterator Id;

  for (Id = dets.begin(); Id != dets.end(); Id++) // search detector using tbName (this is a unique identifier in CORAL)
    if ((*Id)->GetTBName() == hitData.detectorName) // found
      break;
  if (Id == dets.end()) {
    cout << "TGEANT Hit in detector " << hitData.detectorName
        << " and detector id " << hitData.detectorId << " not found." << endl;
    return;
  }
  
  double xin, yin, zin;
  double xout, yout, zout;

  double xbar = hitData.hitPosition[0]; /*mm*/
  double ybar = hitData.hitPosition[1]; /*mm*/
  double zbar = hitData.hitPosition[2]; /*mm*/

  (*Id)->rotatePointMRS2DRSOppanCOMGEANTStyle(hitData.primaryHitPosition[0],
       hitData.primaryHitPosition[1], hitData.primaryHitPosition[2], xin, yin,
       zin,hitData.detectorId);
  (*Id)->rotatePointMRS2DRSOppanCOMGEANTStyle(hitData.lastHitPosition[0],
	   hitData.lastHitPosition[1], hitData.lastHitPosition[2], xout, yout, zout,hitData.detectorId);
  
  // special case for CAMERA hits: we fill the primaryHitPosition in the MCHit
  if (hitData.detectorName[0] == 'C') {
    xbar = hitData.primaryHitPosition[0]; /*mm*/
    ybar = hitData.primaryHitPosition[1]; /*mm*/
    zbar = hitData.primaryHitPosition[2]; /*mm*/
    xout = xin;
    yout = yin;
    zout = zin;
  }

  // special case for hodoscopes:
  int additionalId = 0;
  if (isTrigger) {
    // we have to look at the TGEANT channel number for HO04 hits
    // because the different planes have different channel numbers
    if (hitData.detectorName == "HO04Y1_m"
        || hitData.detectorName == "HO04Y2_m") {
      int channelNo = hitData.channelNo;
      if (channelNo < 5) // ch 0-4
        additionalId = channelNo + 1;
      else if (channelNo < 7) // ch 5-6
        additionalId = channelNo - 4;
      else if (channelNo < 11) // ch 9-10
        additionalId = channelNo - 8;
      else if (channelNo < 16) // ch 11-15
        additionalId = channelNo - 10;
      else if (channelNo < 21) // ch 16-20
        additionalId = channelNo - 15;
      else if (channelNo < 23) // ch 21-22
        additionalId = channelNo - 20;
      else if (channelNo < 25) // ch 23-24
        additionalId = channelNo - 22;
      else if (channelNo < 27) // ch 25-26
        additionalId = channelNo - 24;
      else
        // ch 27-31
        additionalId = channelNo - 26;
      
    // we have to look at the TGEANT channel number for HO03 hits
    // because the different planes have different channel numbers
    } else if (hitData.detectorName == "HO03Y1_m") {
      int channelNo = hitData.channelNo;
      if (channelNo < 6) // ch 0-5
        additionalId = channelNo + 1;
      else if (channelNo < 8) // ch 6-7
        additionalId = channelNo - 5;
      else if (channelNo < 10) // ch 8-9
        additionalId = channelNo - 7;
      else // ch 10-15
        additionalId = channelNo - 9;
      
    // we have to look at the TGEANT channel number for H1 hits
    // because the different planes have different channel numbers
    } else if (hitData.detectorName == "HG01Y1__") {
      int channelNo = hitData.channelNo;
      if (channelNo < 7) // ch 0-6
        additionalId = channelNo + 1;
      else if (channelNo < 13) // ch 7-12
        additionalId = channelNo - 6;
      else if (channelNo < 19) // ch 13-18
        additionalId = channelNo - 12;
      else if (channelNo < 25) // ch 19-24
        additionalId = channelNo - 18;
      else // ch 25-31
        additionalId = channelNo - 24;
      
      // for these detectors we take the TGEANT channel number modulo 8
    } else if (hitData.detectorName == "HL04X1_m"
        || hitData.detectorName == "HL05X1_m"
        || hitData.detectorName == "HM04Y1_d"
        || hitData.detectorName == "HM04Y1_u"
        || hitData.detectorName == "HM05Y1_d"
        || hitData.detectorName == "HM05Y1_u")
      additionalId = hitData.channelNo % 8 + 1;
    
    // here we just add the TGEANT channel number
    else
      additionalId = hitData.channelNo + 1;
  }
  
  CsMCHit* newHit;
  int originID=0;
//   if (hitData.particleId == 11)
//     originID=11;
//   

  newHit = new CsMCTrkHit(xbar, ybar, zbar, xin, yin, zin, xout, yout,
      zout, hitData.energyDeposit / 1000 /*GeV*/,
      hitData.energyDeposit / 1000 /*GeV*/, time /*ns*/,
      Hep3Vector(hitData.momentum[0] / 1000 /*GeV/c*/,
          hitData.momentum[1] / 1000 /*GeV/c*/,
          hitData.momentum[2] / 1000 /*GeV/c*/),
      *(_tTracks.at(_tTrackId[hitData.trackId])), originID, *(*Id),
      hitData.detectorId + additionalId);

  
  _hits.push_back(newHit);
  (*Id)->addMCHit(*newHit);

  
  _tTracks.at(_tTrackId[hitData.trackId])->addMCHit(*newHit);
 
}

void CsGeant3::createTgeantRICHHit(T4RichData& richData)
{
  CsRICH1Detector* rich = CsGeom::Instance()->getRich1Detector();
  if( rich == NULL ) {
	  CsErrLog::mes( elError, 
	   "Hit associated to a non esisting detector, please check detectors.dat content" );
	  return;
	}
	
	//photon detection point MRS [mm]
	double cher_xm = richData.photonHitPosition[0];
	double cher_ym = richData.photonHitPosition[1];
	double cher_zm = richData.photonHitPosition[2];
	//photon production point[mm]
	double cher_xp = richData.photonProductionPosition[0];
	double cher_yp = richData.photonProductionPosition[1];
	double cher_zp = richData.photonProductionPosition[2];
	//TODO Check this pad stuff!!! (note the rotation CG->CORAL)
	double cher_yd = richData.xPadPosition;
	double cher_zd = richData.yPadPosition;
	
	//photon reflection point - set it to 0/0/0
	double cher_xr = 0.0;
	double cher_yr = 0.0;
	double cher_zr = 0.0;
	
	//photon energy
	double cher_eph = richData.photonEnergy*1e6;
	
	//cherenkov time
	double cher_tim = calcTgeantTime( richData.time , richData.photonHitPosition);
	
	//cherenkov angle
	double cher_ang = richData.cerenkovAngle;
	
	//momentum of motherparticle
	Hep3Vector mothersMomentum = Hep3Vector(richData.momentumMotherParticle[0],richData.momentumMotherParticle[1],richData.momentumMotherParticle[1]);
	
	//parent TrackID
	int parentTrackId = richData.parentTrackId;
	
	
	//particle-origin-indes (0 for from original track, PDG-ID otherwise)
	int cher_im = 0;
	
	//cathode-number is detector ID
	int cher_cathode = richData.detectorId-900;
	
	
	
	CsMCHit* hit = new CsMCRICH1Hit( 
					cher_xm, cher_ym, cher_zm, 
					cher_yd, cher_zd, 
					cher_xp, cher_yp, cher_zp,
					cher_xr, cher_yr, cher_zr,
					cher_eph, cher_tim, cher_ang,
					mothersMomentum,
					*(_tTracks.at(_tTrackId[parentTrackId])),
					cher_im, cher_cathode, *rich );
	_hits.push_back( hit );
	rich->addMCHit( *hit ); 
  
  
  
  
  
}

void CsGeant3::createTgeantCaloHit(T4HitData& hitData)
{
    CsCalorimeter::CalorimeterMCData d;
    d.dE = hitData.energyDeposit / 1000.0; // dE (MeV => GeV)
    int originID=0;
    //if (hitData.particleId == 11)
    //  originID=3;
    d.track_id = originID; // 0=original track entered cell ,
    
    double time = calcTgeantTime(hitData.time, hitData.primaryHitPosition);
 
    d.dT = time; // t-t0 (ns)
    d.cell_id = hitData.channelNo; //cell ID = matrix ID + module number
                                       // I can put these two numbers
                                       // separetly, if needed
        CsDet::AddMCHitAll(4, &d);
}

double CsGeant3::calcTgeantTime(double _time, double* _hitPosition)
{
    double lightspeed = 299.793; /*mm/ns*/
    double distance = sqrt(
      pow(_hitPosition[0] , 2)
          + pow(_hitPosition[1] , 2)
          + pow(_hitPosition[2] - zBeamStart, 2)); /*mm*/
    double time = _time - distance/lightspeed; /*ns*/
    return time;
}

// ###### TGEANT CODE END ######
#endif
