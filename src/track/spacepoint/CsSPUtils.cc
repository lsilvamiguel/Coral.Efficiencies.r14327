	// $Id: CsSPUtils.cc,v 1.13 2010/01/28 12:51:26 tnagel Exp $

/*!
   \file    CsEventUtils.cc
   \brief   Compass Event Utilities Class.
   \author  Hugo Pereira
   \version $Revision: 1.13 $
   \date    $Date: 2010/01/28 12:51:26 $
*/

#include "CsSPUtils.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CsOpt.h"
#include "CsErrLog.h"

#include "CsDetector.h"
#include "CsDetFamily.h"
#include "CsSPMaker.h"
#include <time.h>

using namespace std;

//_____________________________________________________________________________
CsSPUtils* CsSPUtils::instance_ = 0;

CsSPUtils* CsSPUtils::Instance( void ) {
  if( instance_ == 0 ) instance_ = new CsSPUtils( );
  return( instance_ );
}

//_____________________________________________________________________________
CsSPUtils::CsSPUtils( void )
{}
	
//_____________________________________________________________________________
bool CsSPUtils::getOptForFamily( const CsDetFamily *pf, string tag, string key, double &par)
{
  CsOpt* opt = CsOpt::Instance();
  list<string> p;    			// parameter list << must contain 2 items
  list<string>::iterator Ip;    			// parameter list << must contain 2 items
	int id = pf->getID();
  int newid;
	bool found = false;
	double c;
	
	while( opt->getOptRec( tag, key, p ) ) {

	  if( p.size()!= 2 ) {
	    CsErrLog::Instance()->mes(elWarning,"Wrong format for family parameter");
	    continue;
	  }
	  
	  int i = 0;
	  for( Ip=p.begin(); Ip!=p.end(); Ip++, i++ ) {
	    switch ( i ) {
	    case 0: istringstream( (*Ip) ) >> newid; break;
	    case 1: istringstream( (*Ip) ) >> c;     break;
	    default: break;
	    }
	  }
	  if( newid == id ) { par = c; found = true; break;}
	}
	return found;
}	

//_____________________________________________________________________________
bool CsSPUtils::getOptForFamily( const CsDetFamily *pf, string tag, string key, string &word)
{
  CsOpt* opt = CsOpt::Instance();
  list<string> p;    			// parameter list << must contain 2 items
	int id = pf->getID();
	int newid;
	bool found = false;
	
	while( opt->getOptRec( tag, key, p ) ) {

	  if( p.size()!= 2 ) {
	    CsErrLog::Instance()->mes(elWarning,"Wrong format for family parameter");
	    continue;
	  }
	 
	  istringstream( p.front() ) >> newid;
	  if( newid == id ) { word = p.back(); found = true; break;}
	}
	return found;
}	
//_____________________________________________________________________________
bool CsSPUtils::getOptForFamily(  const CsDetFamily *pf, string tag, string key, list<string> &words)
{
  CsOpt* opt = CsOpt::Instance();
  list<string> p;    			// parameter list << must contain 2 items
  list<string>::iterator Ip;    			// parameter list << must contain 2 items
	int id = pf->getID();
	int newid;
	bool found = false;
	
  words.clear();
	while( opt->getOptRec( tag, key, p ) ) {

	  if( p.size()< 2 ) {
	    CsErrLog::Instance()->mes(elWarning,"Wrong format for family parameter");
	    continue;
	  }
	  
	  istringstream( p.front() ) >> newid;
	  if( newid == id ) { 
	    for( Ip = p.begin(); Ip != p.end(); Ip++ ) 
	      if( Ip != p.begin() ) words.push_back( *Ip ); 
	    found = true; 
	    break;
	  }
	}
	return found;
}	

//______________________________________________________________________
string CsSPUtils::getDate( void )
{
	struct tm *tPoint;
	time_t now=time( NULL );	
  
	tPoint=localtime( &now );
	if( tPoint == NULL ) {
    CsErrLog::Instance()->mes( elError, "CsSPUtils::getDate.  Can't eval Local time" );
		return " ";
  }
	
  char dayS[3];
	if ( tPoint->tm_mday<10 ) sprintf( dayS, "0%1i", tPoint->tm_mday ); 
	else sprintf( dayS, "%2i", tPoint->tm_mday );

  char monthS[3];
	if (tPoint->tm_mon+1<10) sprintf( monthS, "0%1i", tPoint->tm_mon+1); 
	else sprintf( monthS, "%2i", tPoint->tm_mon+1 );

  char yearS[5];
  sprintf( yearS,"%4i",1900+tPoint->tm_year );
  
	return string(dayS)+"/"
    +string(monthS)+"/"
    +string(yearS);
}

//______________________________________________________________________
string CsSPUtils::getTime( void )
{
	struct tm *tPoint;
	time_t now=time( NULL );	;

	tPoint=localtime( &now );
	if( tPoint == NULL ) {
    CsErrLog::Instance()->mes( elError, "CsSPUtils::getDate.  Can't eval Local time" );
		return " ";
  }
	 
	char hS[3];
	char mS[3];
	char sS[3];

	// hour
	if( tPoint->tm_hour< 10 ) sprintf(hS,"0%1i",tPoint->tm_hour);
	else sprintf(hS,"%2i",tPoint->tm_hour);
	
	// minute 
	if( tPoint->tm_min< 10 ) sprintf(mS,"0%1i",tPoint->tm_min);
	else sprintf(mS,"%2i",tPoint->tm_min);
	
	// second 
	if( tPoint->tm_sec< 10 ) sprintf(sS,"0%1i",tPoint->tm_sec);
	else sprintf(sS,"%2i",tPoint->tm_sec);

  return string(hS)+":"
    +string(mS)+":"
    +string(sS);
}
