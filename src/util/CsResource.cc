// $Id: CsResource.cc,v 1.4 2010/01/18 09:41:32 tnagel Exp $ 

/*!
  \file CsResource.cc
  \brief Debugging object to read and store coral voracity
  \author Hugo Pereira
  \version $Revision: 1.4 $
  \date $Date: 2010/01/18 09:41:32 $
*/


#include <iostream>
#include <cstdlib>   // malloc()

#include "CsResource.h"

CsResource* CsResource::instance_ = 0;

//_____________________________________________________________
CsResource* CsResource::Instance( void ) 
{
  if( !instance_ ) instance_ = new CsResource( ); 
  return instance_;
}

//_____________________________________________________________
CsResource::CsResource( void )
{
  //=== Set Log File Name
  logFName_ = "coral.resource.log";
  OldRC_ = new CsResource::CsRCData();
  NewRC_ = new CsResource::CsRCData();

  //=== Write Header to LogFile
  FILE *out;
  out = fopen( logFName_.c_str(), "w" );
  if( !out ) {
    std::cout << "CsResource::CsResource - ERROR: cannot open file \"" << logFName_ << "\".\n";
    return;
  } else {  
    std::cout << "CsResource::CsResource - INFO: RC log file is \"" << logFName_ << "\".\n";
    fprintf( out, "%20s | %6s | %6s | %6s | %6s | %6s | %6s | %6s \n",
      " ", "total", "locked", "resid", "data", "stack", "exec", "lib");
    fclose( out );
  }
  
}

//_____________________________________________________________
CsResource::~CsResource( void )
{
  delete OldRC_;
  delete NewRC_;
  return;
}

//_____________________________________________________________
void CsResource::dumpRCToLog( std::string tag )
{

  *OldRC_ = *NewRC_;  // Replace OldRC_ content 
  NewRC_->update();  // Update NewRC_ content
  
  // Open LogFile
  FILE *out;
  out = fopen( logFName_.c_str(), "a" );
  if( !out ) {
    std::cout << "CsResource::dumpRCToLog - ERROR: cannot open file \"" << logFName_ << "\".\n";
    return;
  } else {
    fprintf( out, "%20s | %6i | %6i | %6i | %6i | %6i | %6i | %6i\n",
      tag.substr(0,20).c_str(),
      NewRC_->vmSize_, 
      NewRC_->vmLock_, 
      NewRC_->vmResident_, 
      NewRC_->vmData_, 
      NewRC_->vmStack_, 
      NewRC_->vmExec_, 
      NewRC_->vmLib_);
    fclose( out );
  }

  return;
}

//_____________________________________________________________
void CsResource::dumpRCDiffToLog( std::string tag )
{

  *OldRC_ = *NewRC_;  // Replace OldRC_ content 
  NewRC_->update();  // Update NewRC_ content
  
  // Open LogFile
  FILE *out;
  out = fopen( logFName_.c_str(), "a" );
  if( !out ) {
    std::cout << "CsResource::dumpRCToLog - ERROR: cannot open file \"" << logFName_ << "\".\n";
    return;
  } else {
    fprintf( out, "%20s | %6i | %6i | %6i | %6i | %6i | %6i | %6i\n",
      tag.substr(0,20).c_str(),
      int(NewRC_->vmSize_)     - int(OldRC_->vmSize_),
      int(NewRC_->vmLock_)     - int(OldRC_->vmLock_),
      int(NewRC_->vmResident_) - int(OldRC_->vmResident_), 
      int(NewRC_->vmData_)     - int(OldRC_->vmData_),
      int(NewRC_->vmStack_)    - int(OldRC_->vmStack_), 
      int(NewRC_->vmExec_)     - int(OldRC_->vmExec_), 
      int(NewRC_->vmLib_)      - int(OldRC_->vmLib_)
    );
    fclose( out );
  }

  return;
}

//_____________________________________________________________

//_____________________________________________________________
CsResource::CsRCData::CsRCData( void )
{
  if( !( success_ = update() ) ) 
  std::cout << "CsResource::CsRCData::CsRCData - ERROR: Troubles reading resources\n";
  return;
} 
  
//_____________________________________________________________
CsResource::CsRCData& CsResource::CsRCData::operator = ( const CsResource::CsRCData& rc )
{
  if( this != &rc ) {
    success_    = rc.success_;
    vmSize_     = rc.vmSize_;
    vmLock_     = rc.vmLock_;
    vmResident_ = rc.vmResident_;
    vmData_     = rc.vmData_;
    vmStack_    = rc.vmStack_;
    vmExec_     = rc.vmExec_;
    vmLib_      = rc.vmLib_; 
  } 
  return *this;
}   
  
//_____________________________________________________________
bool CsResource::CsRCData::update( void )
{

  //=== Get proc fileName 
  char filename[80];
  sprintf( filename, "/proc/%d/status", getpid() );
  
  FILE *fid;
  fid = fopen( filename, "r" );
  bool success = true;

  if( fid ) {
    int   bufsize = 1024;
    char* buffer = (char*) malloc( bufsize );
    int   status = fread( buffer, 1, bufsize, fid );
    fclose( fid );
    
    if( status > 0 ) {
      char* cursor = strstr( buffer, "VmSize:" );
      if( cursor ) {
	      sscanf( cursor, 
      		"VmSize: %iu kB\n"
      		"VmLck:  %iu kB\n"
      		"VmRSS:  %iu kB\n"
      		"VmData: %iu kB\n"
      		"VmStk:  %iu kB\n"
      		"VmExe:  %iu kB\n"
      		"VmLib:  %iu kB\n",
      		&vmSize_,  &vmLock_, &vmResident_, &vmData_,
      		&vmStack_, &vmExec_, &vmLib_ );
      } else success = false;
    }   else success = false;
    free( buffer );
  }
  else {
    std::cout << "CsResource::CsRCData::Update - ERROR: cannot open \"" << filename << "\"\n";
    success = false;
  }

  return success;

}

  
