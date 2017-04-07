// $Id: CsDateStore.cc,v 1.25 2010/08/04 07:25:25 suhl Exp $

/*!
   \file    CsDateStore.cc
   \brief   Implementation of CsDateStore class member functions
   \author  Sergei.Gerassimov@cern.ch
   \version $Revision: 1.25 $
   \date    $Date: 2010/08/04 07:25:25 $   
*/

#include "coral_config.h"

#include <ctype.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "CsInit.h"

#if USE_COMPASS_Date
// long32 defined by the compiler
#define long64 long long
#include "monitor.h"
#include "event.h"
#else
struct eventStructV3 { unsigned int eventHeader[20]; unsigned short rawData[1]; };
struct eventStructV5 { unsigned int eventHeader[17]; unsigned short rawData[1]; };
#define ERR_BASE          0x80000000
#define MON_ERR_EOF       (ERR_BASE|0x8000)
#define MON_ERR_SYS_ERROR (ERR_BASE|0x1000)
#define MON_ERR_BAD_EVENT (ERR_BASE|0x0008)
#define MON_ERR_MALLOC    (ERR_BASE|0x0080)
#endif

#include "CsDateStore.h"
#include "CsOpt.h"

#ifdef COMPASS_USE_OSPACE_STD
# include <ospace/std/algorithm>
# include <ospace/std/iomanip>
#else
# include <algorithm>
# include <iomanip>
#endif

#if USE_RFIO
# include <shift.h>
#else
# include <sys/stat.h>
# include <sys/types.h>
# include <unistd.h>
# include <fcntl.h>
#endif

using namespace std;

CsDateStore* CsDateStore::instance_ = 0;


CsDateStore* CsDateStore::Instance() {
  if( instance_ ==  0 ) {
    instance_ = new CsDateStore();
  }
  return( instance_ );
}

CsDateStore::CsDateStore() {
  init_ = false;
  scan_ = false;
  // Added: 20050109, Benigno
  _fileid = 0;

  // Added: 20011001, Benigno
  _files = CsInit::Instance()->getDateFilesList();
  _currentFile = _files.begin();
  // Added 20030306, Duic

}

//--------------------------------------------------------------

CsDateStore::~CsDateStore() {
  instance_ = 0;
}

//--------------------------------------------------------------

bool CsDateStore::init() {

  if(init_) {
    cerr << "CsDateStore::init: init already performed " << endl;
    return(true);
  }

  CsInit* coralInit = CsInit::Instance();
  init_ = true;
  return(true);
}

//--------------------------------------------------------------

bool CsDateStore::scan() {

  // Modified: 20011001, Benigno

  if( !init_ ) {
    cerr << "CsDateStore::scan() ERROR: scan already performed" << endl;
    exit(1);
  }

  if( _files.empty() ) {
    cerr << "CsDateStore:scan() ERROR: No files to be read." << endl;
    exit(1);
  }

  // fill DATE monitoring policy table

  int status;
  bool allFine = true;

# if USE_COMPASS_Date
  static char* MonTab[3] = {strdup("ALL"), strdup("yes"), NULL};
  status = monitorDeclareTable(MonTab); // declare monitoring policy (DATE lib)
  if ( status != 0 ) {
    cerr<< "CsDateStore::scan() ERROR from monitorDeclareTable():" 
	<< monitorDecodeError( status ) << "."
	<< endl;
    return(false);
  } else {
    cout << "CsDateStore::scan() INFO: the monitoring policy table is " ;
    for(int i = 0 ; i < 2; i++) cout << MonTab[i] << " ";
    cout << endl;
  }
# else
  
# endif

  // file(s) stream or date stream?

  if( _files.size() == 1 && _files.front()->c_str()[0] == '@' ) {
# if USE_COMPASS_Date
    char* data_source = const_cast<char*>(_files.front()->c_str());
    status = monitorSetDataSource(data_source); 
    if( status == 0 ) { 
      allFine = true;
      cout<< "CsDateStore::scan() INFO: input is a data stream " << endl 
	  << "DATE stream " << data_source << " is set as input."
	  << endl;
    }
    else {
      allFine = false;
      cout<< "CsDateStore::scan() ERROR from monitorSetDataSource: " 
	  << endl;
    }
# else
    cerr<< "CsDateStore::scan() ERROR You need to build CORAL with DATE Libs to access date streams"
	<< endl;
    return(false);
# endif
  }
  else {
  
  //  It should be a list of file(s). Check if all files exist 

    list<string*>::iterator is;
  
    for( is=_files.begin(); is!=_files.end(); is++ ) {

      char* data_source = const_cast<char*>((*is)->c_str());

      // use rfio_stat or stat funct. to check if files exist (Benigno 20011106)
      struct stat statbuf;

#if USE_RFIO
      status = rfio_stat( data_source, &statbuf ); 
#else
      status = stat( data_source, &statbuf ); 
#endif  

      if ( status != 0 ) {
	cerr<< "CsDateStore::scan() ERROR from monitorSetDataSource: " 
	    << endl;
#if USE_RFIO
	rfio_perror( data_source );
#else
	perror( data_source );
#endif
	allFine = false;
      } 
      else {
	cout << data_source << " file found." << endl;
      }
    }

    // If all OK tell DATE to set on 1st file...
    if( allFine ) {
      char* data_source = const_cast<char*>((*_currentFile)->c_str());
#     if USE_COMPASS_Date
      status = monitorSetDataSource(data_source); 
      cout<< "CsDateStore::scan() INFO: the scan was fine." << endl 
	  << "File " << data_source << " is set as first input one."
	  << endl;
#     else
#     if USE_RFIO
      _fileid = rfio_open( data_source, O_RDONLY );
#     else
      _fileid = open( data_source, O_RDONLY );
#     endif
      if( _fileid == -1 ) {
	cerr<< "CsDateStore::scan() ERROR: the scan failed." << endl 
	    << "File " << data_source << " not opened."
	    << endl;
	allFine = false;
      }
      else {
	cout<< "CsDateStore::scan() INFO: the scan was fine." << endl 
	    << "File " << data_source << " ( id:" << _fileid
	    << ") is set as first input one."
	    << endl;
      }
#     endif
    }
  }

  scan_ = allFine;
  return(scan_);

}

//--------------------------------------------------------------

bool CsDateStore::next() {
  
  // Midified: 20011001, Benigno
  
  bool dump = false; // debug printout flag
  
  if( !scan_ ) {
    cerr << "CsDateStore::next() ERROR: no scan performed, exit!" << endl;
    exit(1);
  }  
  
  static void *ptr = NULL;                 // this is plane C style
  if(ptr != NULL) { free(ptr); ptr=NULL; } // never do like this in C++ !
  
  list<uint32>::iterator ityp;

  //
  // DATE input
  //
  while(true){ // "endless" loop (untill one get requested event type)

    // get next event (DATE lib)
#   if USE_COMPASS_Date
    unsigned int status = monitorGetEventDynamic( &ptr ); 
#   else
    unsigned int status = _getEvent( &ptr ); 
#   endif
    
    while( status == MON_ERR_EOF ) { // end of file reached. Move to another...
      _currentFile++;
      if( _currentFile == _files.end() ) {
	cout<< "CsDateStore::next() INFO: end of all files reached" << endl;
	return( false );
      }
      char* data_source = const_cast<char*>((*_currentFile)->c_str());
      cout <<  "CsDateStore::next() INFO: Switch to " 
	   << data_source << " input file."
	   << endl;
#     if USE_COMPASS_Date
      status = monitorLogout();
#     else
#     if USE_RFIO
      status = rfio_close( _fileid );
#     else
      status = close( _fileid );
#     endif
#     endif      
#     if USE_COMPASS_Date
      status = monitorSetDataSource(data_source); 
#     else
#     if USE_RFIO
      _fileid = rfio_open( data_source, O_RDONLY );
#     else
      _fileid = open( data_source, O_RDONLY );
#     endif
      if( _fileid == -1 ) {
	cerr<< "CsDateStore::next() : ERROR." << endl 
	    << "File " << data_source << " not opened."
	    << endl;
	return MON_ERR_SYS_ERROR;
      }
#     endif
      if(ptr != NULL) { free(ptr); ptr=NULL; } 
#     if USE_COMPASS_Date
      status = monitorGetEventDynamic( &ptr ); 
#     else
      status =_getEvent( &ptr ); 
#     endif
    }
      
    if ( status != 0 && status != MON_ERR_EOF ) {
#     if USE_COMPASS_Date      
      cout<< "CsDateStore::next() ERROR from monitorGetEventDynamic: "
	  << monitorDecodeError( status ) << "." 
	  << endl;
#     endif
      if( status == MON_ERR_SYS_ERROR ) {
#       if USE_RFIO
	rfio_perror( "CsDateStore::next() ERROR, system related: " );
#       else
	perror( "CsDateStore::next() ERROR, system related: " );
#       endif
      }  
      return(false); // seems, end of file. Terminate the job. 
    }
      
    // exit loop if physics event had been read-in.
#   if USE_COMPASS_Date    
    struct eventStruct* event = (struct eventStruct*) ptr; 
#   if defined (EVENT_H_ID) && EVENT_H_ID == 0x00020018
    if( (event->eventHeader.type & EVENT_TYPE_MASK) == physicsEvent ) break; // see eventType enumeration in event.h)
#   else
#   if defined (EVENT_MAJOR_VERSION_NUMBER) &&  EVENT_MAJOR_VERSION_NUMBER == 0x0003
#   if defined (EVENT_MINOR_VERSION_NUMBER) &&  EVENT_MINOR_VERSION_NUMBER == 0x0006
    if( event->eventHeader.eventType == physicsEvent ) break; // see eventType enumeration in event.h)
#   endif
#   endif
#   endif
#   else
    if( _dateVersionGuess == 3 ) {
      struct eventStructV3* event = (struct eventStructV3*) ptr; 
      if( (event->eventHeader[2] & 0x0000ffff) == 7 ) break; 
    }
    else if( _dateVersionGuess == 5 ) {
      struct eventStructV5* event = (struct eventStructV5*) ptr; 
      if( event->eventHeader[4] == 7 ) break;
    }
#   endif
    
  } // end of "endless" loop 
  
  // save pointer to DATE buffer
  this->date_buffer_ptr = (uint8*) ptr;
  
  return(true);
}

//--------------------------------------------------------------

uint8*  CsDateStore::rawBuffer() {
  return(date_buffer_ptr);
}

//--------------------------------------------------------------

uint32  CsDateStore::getEventInRun() const {
#if USE_COMPASS_Date
  struct eventStruct* event = (struct eventStruct*) this->date_buffer_ptr;
# if defined (EVENT_H_ID) && EVENT_H_ID == 0x00020018
  return(uint32(event->eventHeader.nbInRun));
# else
# if defined (EVENT_MAJOR_VERSION_NUMBER) &&  EVENT_MAJOR_VERSION_NUMBER == 0x0003
# if defined (EVENT_MINOR_VERSION_NUMBER) &&  EVENT_MINOR_VERSION_NUMBER == 0x0006
  return(uint32(event->eventHeader.eventId[0]));
# endif
# endif
# endif
#else
  if( _dateVersionGuess == 3 ) {
    struct eventStructV3* event = (struct eventStructV3*) this->date_buffer_ptr; 
    return(uint32(event->eventHeader[6])); 
  }
  else if( _dateVersionGuess == 5 ) {
    struct eventStructV5* event = (struct eventStructV5*) this->date_buffer_ptr; 
    return(event->eventHeader[6]);
  }
  else { //should never happen
    cerr << "CsDateStore::getEventInRun(): wrong dateVersion guess.\n";
    exit(1);
  }    
#endif
  
}

//--------------------------------------------------------------

uint32  CsDateStore::getRun()  const {
#if USE_COMPASS_Date
  struct eventStruct* event = (struct eventStruct*) this->date_buffer_ptr;
# if defined (EVENT_H_ID) && EVENT_H_ID == 0x00020018
  return(uint32(event->eventHeader.runNb));
# else
# if defined (EVENT_MAJOR_VERSION_NUMBER) &&  EVENT_MAJOR_VERSION_NUMBER == 0x0003
# if defined (EVENT_MINOR_VERSION_NUMBER) &&  EVENT_MINOR_VERSION_NUMBER == 0x0006
  return(uint32(event->eventHeader.eventRunNb));
# endif
# endif
# endif
#else
  if( _dateVersionGuess == 3 ) {
    struct eventStructV3* event = (struct eventStructV3*) this->date_buffer_ptr; 
    return(uint32(event->eventHeader[4])); 
  }
  else if( _dateVersionGuess == 5 ) {
    struct eventStructV5* event = (struct eventStructV5*) this->date_buffer_ptr; 
    return(event->eventHeader[5]);
  }
  else { //should never happen
    cerr << "CsDateStore::getRun(): wrong dateVersion guess.\n";
    exit(1);
  }    
#endif
} 

//--------------------------------------------------------------

uint32  CsDateStore::getEventInBurst()  const
{
#if USE_COMPASS_Date
  struct eventStruct* event = (struct eventStruct*) this->date_buffer_ptr;
# if defined (EVENT_H_ID) && EVENT_H_ID == 0x00020018
  return(uint32(event->eventHeader.nbInBurst));
# else
# if defined (EVENT_MAJOR_VERSION_NUMBER) &&  EVENT_MAJOR_VERSION_NUMBER == 0x0003
# if defined (EVENT_MINOR_VERSION_NUMBER) &&  EVENT_MINOR_VERSION_NUMBER == 0x0006
  return(uint32(event->eventHeader.eventId[1] & 0x000fffff));
# endif
# endif
# endif
#else
  if( _dateVersionGuess == 3 ) {
    struct eventStructV3* event = (struct eventStructV3*) this->date_buffer_ptr; 
    return(uint32(event->eventHeader[7])); 
  }
  else if( _dateVersionGuess == 5 ) {
    struct eventStructV5* event = (struct eventStructV5*) this->date_buffer_ptr; 
    return(event->eventHeader[7] & 0x000fffff);
  }
  else { //should never happen
    cerr << "CsDateStore::getEventInBurst(): wrong dateVersion guess.\n";
    exit(1);
  }    
#endif
}

//--------------------------------------------------------------

uint32  CsDateStore::getBurst() const
{
#if USE_COMPASS_Date
  struct eventStruct* event = (struct eventStruct*) this->date_buffer_ptr;
# if defined (EVENT_H_ID) && EVENT_H_ID == 0x00020018
  return(uint32(event->eventHeader.burstNb));
# else
# if defined (EVENT_MAJOR_VERSION_NUMBER) &&  EVENT_MAJOR_VERSION_NUMBER == 0x0003
# if defined (EVENT_MINOR_VERSION_NUMBER) &&  EVENT_MINOR_VERSION_NUMBER == 0x0006
  return(uint32(((event->eventHeader.eventId[1])>>20) & 0x00000fff ));
# endif
# endif
# endif
#else
  if( _dateVersionGuess == 3 ) {
    struct eventStructV3* event = (struct eventStructV3*) this->date_buffer_ptr; 
    return(uint32(event->eventHeader[5])); 
  }
  else if( _dateVersionGuess == 5 ) {
    struct eventStructV5* event = (struct eventStructV5*) this->date_buffer_ptr; 
    return( ((event->eventHeader[7])>>20) & 0x00000fff);
  }
  else { //should never happen
    cerr << "CsDateStore::getBurst(): wrong dateVersion guess.\n";
    exit(1);
  }    
#endif
}

//--------------------------------------------------------------

uint32  CsDateStore::getTriggerMask() const
{
#if USE_COMPASS_Date
  struct eventStruct* event = (struct eventStruct*) this->date_buffer_ptr;
# if defined (EVENT_H_ID) && EVENT_H_ID == 0x00020018
  return(uint32(event->eventHeader.typeAttribute[1]));
# else
# if defined (EVENT_MAJOR_VERSION_NUMBER) &&  EVENT_MAJOR_VERSION_NUMBER == 0x0003
# if defined (EVENT_MINOR_VERSION_NUMBER) &&  EVENT_MINOR_VERSION_NUMBER == 0x0006
  return(uint32(event->eventHeader.eventTriggerPattern[0]));
# endif
# endif
# endif
#else
  if( _dateVersionGuess == 3 ) {
    struct eventStructV3* event = (struct eventStructV3*) this->date_buffer_ptr; 
    return(uint32(event->eventHeader[19])); 
  }
  else if( _dateVersionGuess == 5 ) {
    struct eventStructV5* event = (struct eventStructV5*) this->date_buffer_ptr; 
    return(event->eventHeader[8]); // to be tested
  }
  else { //should never happen
    cerr << "CsDateStore::getTriggerMask(): wrong dateVersion guess.\n";
    exit(1);
  }    
#endif
}

//--------------------------------------------------------------

uint32  CsDateStore::getErrorCode()        const
{
#if USE_COMPASS_Date
  struct eventStruct* event = (struct eventStruct*) this->date_buffer_ptr;
# if defined (EVENT_H_ID) && EVENT_H_ID == 0x00020018
  return(uint32(event->eventHeader.errorCode));
# else
# if defined (EVENT_MAJOR_VERSION_NUMBER) &&  EVENT_MAJOR_VERSION_NUMBER == 0x0003
# if defined (EVENT_MINOR_VERSION_NUMBER) &&  EVENT_MINOR_VERSION_NUMBER == 0x0006
  cerr << "CsDateStore::getErrorCode(): To be implemented with Date V5\n";
  exit(1);
# endif
# endif
# endif
#else
  if( _dateVersionGuess == 3 ) {
    struct eventStructV3* event = (struct eventStructV3*) this->date_buffer_ptr; 
    return(uint32(event->eventHeader[15])); 
  }
  else if( _dateVersionGuess == 5 ) {
    struct eventStructV5* event = (struct eventStructV5*) this->date_buffer_ptr; 
    cerr << "CsDateStore::getErrorCode(): To be implemented with Date V5\n";
    exit(1);
  }
  else { //should never happen
    cerr << "CsDateStore::getErrorCode(): wrong dateVersion guess.\n";
    exit(1);
  }    
#endif
} 
 
//--------------------------------------------------------------

CsTime  CsDateStore::getTime()  const
{
#if USE_COMPASS_Date
  struct eventStruct* event = (struct eventStruct*) this->date_buffer_ptr;
# if defined (EVENT_H_ID) && EVENT_H_ID == 0x00020018
  CsTime t(event->eventHeader.time, event->eventHeader.usec);
  return(t);
# else
# if defined (EVENT_MAJOR_VERSION_NUMBER) &&  EVENT_MAJOR_VERSION_NUMBER == 0x0003
# if defined (EVENT_MINOR_VERSION_NUMBER) &&  EVENT_MINOR_VERSION_NUMBER == 0x0006
  CsTime t(event->eventHeader.eventTimestamp,0);
  return(t);
# endif
# endif
# endif
#else
  if( _dateVersionGuess == 3 ) {
    struct eventStructV3* event = (struct eventStructV3*) this->date_buffer_ptr; 
    CsTime t(event->eventHeader[13], event->eventHeader[14]);
    return(t);
  }
  else if( _dateVersionGuess == 5 ) {
    struct eventStructV5* event = (struct eventStructV5*) this->date_buffer_ptr; 
    CsTime t(event->eventHeader[16], 0);
    return(t);
  }
  else { //should never happen
    cerr << "CsDateStore::getTime(): wrong dateVersion guess.\n";
    exit(1);
  }    
#endif
}

//--------------------------------------------------------------

unsigned int CsDateStore::_getEvent( void *ptr ) {

  _dateVersionGuess = 0;

  int status = 0;
# define MY_BUFF_SIZE 16
  unsigned int buff[MY_BUFF_SIZE/4];
  uint8 *p;

  // read first 4 entries
# if USE_RFIO
  status = rfio_read( _fileid, buff, MY_BUFF_SIZE ); 
# else
  status = read( _fileid, buff, MY_BUFF_SIZE ); 
# endif  
  if( status != MY_BUFF_SIZE ) {
    return status == 0 ? MON_ERR_EOF : MON_ERR_SYS_ERROR;
  }

  // check "event magic"
  if( buff[1] != 0xDA1E5AFE ) return MON_ERR_BAD_EVENT;
      
  // guess date version: in V5.13 event.h version is 3.06 and header size is 68....
  if( buff[3] == 0x00030006 && buff[2] == 68 ) {
    _dateVersionGuess = 5;
  } 
  // in (compass) date version 3.7.1 header size is 80....
  else if( buff[3] == 80 ) {
    _dateVersionGuess = 3;
  }
  else {
    _dateVersionGuess = 0;
    return MON_ERR_SYS_ERROR;
  }

  // that's fine, read all the rest of the event.

  if( ( p = (uint8*)malloc( buff[0] ) ) == NULL ) {
    return MON_ERR_MALLOC;
  }
  memcpy( p, buff, MY_BUFF_SIZE );

# if USE_RFIO
  status = rfio_read( _fileid, p+MY_BUFF_SIZE, buff[0]-MY_BUFF_SIZE ); 
# else
  status = read( _fileid, p+MY_BUFF_SIZE, buff[0]-MY_BUFF_SIZE ); 
# endif  
  if( (unsigned int)status != (buff[0]-MY_BUFF_SIZE) ) {
    return MON_ERR_BAD_EVENT;
  }
  
  *(void**) ptr = (void*) p;

  return 0;

}







