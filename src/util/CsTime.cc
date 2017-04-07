// $Id: CsTime.cc,v 1.14 2010/01/28 12:51:26 tnagel Exp $

/*!
  \file    CsTime.cc
  \brief   CsTime class
  \author  Benigno Gobbo
  \version $Revision: 1.14 $
  \date    $Date: 2010/01/28 12:51:26 $
*/  

//#include <ctime>
#include <sys/time.h>
#include "CsTime.h"
# ifdef COMPASS_USE_OSPACE_STD
#  include <ospace/std/iomanip>
# else
#  ifdef __HP_aCC
#   include <iomanip.h>
#  else
#   include <iomanip>
#  endif // __HP_aCC
# endif // COMPASS_USE_OSPACE_STD

using namespace std;

const int CsTime::_usecInSec = 1000000; 

bool CsTime::_dmy = true;

CsTime::CsTime() {

  struct timeval tp;
# ifdef __hpux
  void *tzp = 0;
  gettimeofday( &tp, tzp );
# else
  //struct timezone tzp;
  //gettimeofday( &tp, &tzp );
  gettimeofday( &tp, NULL );
# endif

  _secFrEpoch = tp.tv_sec;
  _usec       = tp.tv_usec;

}

CsTime::CsTime( time_t secFrEpoch, int usec ) {
  int secs;

  if( usec >= 0 ) {
    secs = usec / _usecInSec;
    usec = usec % _usecInSec;
  }
  else {
    secs = usec / _usecInSec - 1;
    usec = _usecInSec + usec % _usecInSec;
  }
  _secFrEpoch = secFrEpoch + secs;
  _usec = usec; 
}

CsTime::CsTime( int year, int month, int day, 
		int hour, int min, 
		int sec, int usec ) {

  struct tm tmtime; 
  tmtime.tm_year = year - 1900; 
  tmtime.tm_mon  = month - 1;
  tmtime.tm_mday = day;
  tmtime.tm_hour = hour;
  tmtime.tm_min  = min;
  tmtime.tm_sec  = sec;
 
  _secFrEpoch = mktime( &tmtime );
  _usec       = usec;
}

CsTime::CsTime( const CsTime & ct ) :
  _secFrEpoch(ct._secFrEpoch), _usec(ct._usec) 
{}
        
CsTime& CsTime::operator=( const CsTime & ct ) {
  if( this != &ct ) {
    _secFrEpoch = ct._secFrEpoch;  
    _usec = ct._usec; 
  }
  return( *this );
}          

bool CsTime::operator==( const CsTime & ct ) const {
  if( _secFrEpoch == ct._secFrEpoch && _usec == ct._usec ) 
    return( true );
  else
    return( false );
}
    
bool CsTime::operator!=( const CsTime & ct ) const {    
  return !operator==( ct );
}

bool CsTime::operator>( const CsTime & ct ) const {    
  if( _secFrEpoch > ct._secFrEpoch ||
      ( _secFrEpoch == ct._secFrEpoch && _usec > ct._usec )) 
    return( true );
  else
    return( false );
}

bool CsTime::operator<=( const CsTime & ct ) const {
  return !operator>( ct );
}

bool CsTime::operator<( const CsTime & ct ) const {   
  if( _secFrEpoch < ct._secFrEpoch ||
      ( _secFrEpoch == ct._secFrEpoch && _usec < ct._usec )) 
    return( true );
  else
    return( false );
}

bool CsTime::operator>=( const CsTime & ct ) const {    
  return !operator<( ct );
}

ostream& operator<<( ostream & st, const CsTime & dt ) {
  char gDate[40]; 
  char gZone[10];
  char fmt1[] = "%a, %d/%b/%Y %T";
  char fmt2[] = "(%Z)";
  time_t sfe = dt._secFrEpoch;

  if( dt._dmy ) {
    strftime( gDate, 40, fmt1, gmtime( &sfe ) );
    strftime( gZone, 10, fmt2, gmtime( &sfe ) );
    char oldfill=st.fill('0');
    st << gDate << "." << setw(6) << dt._usec << " " << gZone;
    st.fill(oldfill);
  }
  else {
    char oldfill=st.fill('0');
    st << dt._secFrEpoch << "." << setw(6) << dt._usec;
    st.fill(oldfill);
  }

  return( st ); 
} 

ostream& dmy(ostream& os) {
  CsTime::_setDmy(true);
  return os;
}

ostream& sec(ostream& os) {
  CsTime::_setDmy(false);
  return os;
}
