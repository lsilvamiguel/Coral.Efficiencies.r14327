// $Id: CsStopwatch.cc,v 1.6 2007/02/05 10:18:54 gobbo Exp $

/*!
  \file    CsStopwatch.cc
  \brief   Compass class for elapsed time measurements
  \author  Benigno Gobbo
  \version $Revision: 1.6 $
  \date    $Date: 2007/02/05 10:18:54 $
*/  

#include <sys/time.h>
//#include <sys/times.h>
#include <sys/resource.h>
#include <unistd.h>
#include "CsStopwatch.h"

CsStopwatch::CsStopwatch() {
  _ntim   = 0;
  _maxtim = 10;
  _startS = new int[_maxtim];
  _startU = new int[_maxtim];
  _used   = new bool[_maxtim];
  for( int i=0; i<_maxtim; i++ ) {
    _startS[i] = 0;
    _startU[i] = 0;
    _used[i]   = false;
  }
}

CsStopwatch::CsStopwatch( int maxtim ) {
  _ntim   = 0;
  _maxtim = maxtim;
  _startS = new int[_maxtim];
  _startU = new int[_maxtim];
  _used   = new bool[_maxtim];
  for( int i=0; i<_maxtim; i++ ) {
    _startS[i] = 0;
    _startU[i] = 0;
    _used[i]   = false;
  }
}

CsStopwatch::~CsStopwatch() {
  delete [] _startS;
  delete [] _startU;
  delete [] _used;
}

CsStopwatch::CsStopwatch( const CsStopwatch& sw ) {
  _maxtim  = sw._maxtim;
  _ntim    = sw._ntim;
  _startS  = new int[_maxtim];
  _startU  = new int[_maxtim];
  _used    = new bool[_maxtim];
  // Not elegant, but what's the size of bool?
  for( int i = 0; i<_maxtim; i++ ) {
    //memcpy( _startU, sw._startU, _maxtim*4 );
    _startS[i] = sw._startS[i];
    _startU[i] = sw._startU[i];
    _used[i]   = _used[i]; 
  }
}

CsStopwatch& CsStopwatch::operator=( const CsStopwatch& sw ) {
  if( this  != &sw ) {
    _maxtim  = sw._maxtim;
    _ntim    = sw._ntim;
    _startS  = new int[_maxtim];
    _startU  = new int[_maxtim];
    _used    = new bool[_maxtim];
    // Not elegant, but what's the size of bool?
    for( int i = 0; i<_maxtim; i++ ) {
      //memcpy( _startU, sw._startU, _maxtim*4 );
      _startS[i] = sw._startS[i];
      _startU[i] = sw._startU[i];
      _used[i]   = _used[i]; 
    }
  }
  return( *this );
}          

int CsStopwatch::start() {

  // Get the current time first.
  //struct timeval tp;
//# ifdef __hpux
  //void *tzp = 0;
  //gettimeofday( &tp, tzp );
//# else
  //struct timezone tzp;
  //gettimeofday( &tp, &tzp );
//# endif
  //int s = tp.tv_sec;
  //int u = tp.tv_usec;

  // Get current process time first.
  //struct tms procTime;
  //int tick = sysconf( _SC_CLK_TCK );
  //clock_t sysup = times( &procTime );
  //int s = int( procTime.tms_utime + procTime.tms_stime ) / tick; 
  //int u = int( procTime.tms_utime + procTime.tms_stime ) % tick;
  //u = u * 1000000 / tick;

  // Get current process time first.
  struct rusage usage;
  int status = getrusage( RUSAGE_SELF, &usage );
  struct timeval usrTime = usage.ru_utime;
  struct timeval sysTime = usage.ru_stime;
  int s = usrTime.tv_sec + sysTime.tv_sec;
  int u = usrTime.tv_usec + sysTime.tv_usec;

  // Look for a free chrono
  if( _ntim >= _maxtim ) {
    return( -1 );
  }
  else {
    for( int i=0; i<_maxtim; i++ ) {
      if( !_used[i] ) {
	_ntim++;        // add a stopwatch to used ones
	_startS[i] = s; // set the counters
	_startU[i] = u;
	_used[i]   = true;
	return( i );    // return the stopwatch number
      }
    }
  }
  return( -1 );   // well... Just to avoid boring warnings...
}

double CsStopwatch::inter( int chrono ) {

  // Get the current time first.
  //struct timeval tp;
//# ifdef __hpux
  //void *tzp = 0;
  //gettimeofday( &tp, tzp );
//# else
  //struct timezone tzp;
  //gettimeofday( &tp, &tzp );
//# endif
  //int s = tp.tv_sec;
  //int u = tp.tv_usec;

  // Get current process time first.
  //int tick = sysconf( _SC_CLK_TCK );
  //struct tms procTime;
  //clock_t sysup = times( &procTime );
  //int s = int( procTime.tms_utime + procTime.tms_stime ) / tick; 
  //int u = int( procTime.tms_utime + procTime.tms_stime ) % tick;
  //u = u * 1000000 / tick;

  // Get current process time first.
  struct rusage usage;
  int status = getrusage( RUSAGE_SELF, &usage );
  struct timeval usrTime = usage.ru_utime;
  struct timeval sysTime = usage.ru_stime;
  int s = usrTime.tv_sec + sysTime.tv_sec;
  int u = usrTime.tv_usec + sysTime.tv_usec;

  double t = double(s-_startS[chrono]);
  double v = double(u-_startU[chrono])/1000000.;
  t = t + v;

  return( t );   // return the elapsed time. Do not reset counters
}
  
double CsStopwatch::stop( int chrono ) {

  // Get the current time first.
  //struct timeval tp;
//# ifdef __hpux
  //void *tzp = 0;
  //gettimeofday( &tp, tzp );
//# else
  //struct timezone tzp;
  //gettimeofday( &tp, &tzp );
//# endif
  //int s = tp.tv_sec;
  //int u = tp.tv_usec;

  // Get current process time first.
  //int tick = sysconf( _SC_CLK_TCK );
  //struct tms procTime;
  //clock_t sysup = times( &procTime );
  //int s = int( procTime.tms_utime + procTime.tms_stime ) / tick; 
  //int u = int( procTime.tms_utime + procTime.tms_stime ) % tick;
  //u = u * 1000000 / tick;

  // Get current process time first.
  struct rusage usage;
  int status = getrusage( RUSAGE_SELF, &usage );
  struct timeval usrTime = usage.ru_utime;
  struct timeval sysTime = usage.ru_stime;
  int s = usrTime.tv_sec + sysTime.tv_sec;
  int u = usrTime.tv_usec + sysTime.tv_usec;

  double t = double(s-_startS[chrono]);
  double v = double(u-_startU[chrono])/1000000.;
  t = t + v;

  _startS[chrono] = 0; // reset the counters
  _startU[chrono] = 0;
  _used[chrono]   = false;
  _ntim--;            // free one stopwatch
  return( t );        // return the elapsed time
}
