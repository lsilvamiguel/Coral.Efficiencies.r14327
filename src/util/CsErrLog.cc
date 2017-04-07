// $Id: CsErrLog.cc,v 1.16 2010/01/28 12:51:26 tnagel Exp $

/*!
   \file    CsErrLog.cc
   \brief   Compass Error and Message Logger.
   \author  Benigno Gobbo 
   \version $Revision: 1.16 $
   \date    $Date: 2010/01/28 12:51:26 $
*/

#include <cstdio>
#include <cstdarg>
#include <cstdlib>   // for exit()
#include "CsOpt.h"
#include "CsErrLog.h"

CsErrLog*        CsErrLog::_instance         = NULL;

Severity         CsErrLog::_minLogSeverity   = elWarning;   
Severity         CsErrLog::_minStoreSeverity = elBasicInfo; 
int              CsErrLog::_verbosity        = 0;        
std::vector<logger>   CsErrLog::_log;
bool             CsErrLog::_exitonerror      = false;

CsErrLog* CsErrLog::Instance() {
 if( _instance == NULL ) {
   _instance = new CsErrLog();
 }
 return _instance; 
}

CsErrLog::CsErrLog( ) {

  // Instance the Option Interpreter
  CsOpt* opt = CsOpt::Instance( );

  std::string tag = "error logger";
  std::string key = "log level";
  std::string loglevel;
  if( opt->getOpt( tag, key, loglevel ) ) {
    if( loglevel      == "debug" )   _minLogSeverity = elDebugging;
    else if( loglevel == "verbose" ) _minLogSeverity = elVerbose;
    else if( loglevel == "info" )    _minLogSeverity = elInfo;
    else if( loglevel == "anomaly" ) _minLogSeverity = elAnomaly;
    else if( loglevel == "warning" ) _minLogSeverity = elWarning;
    else if( loglevel == "error" )   _minLogSeverity = elError;
    else if( loglevel == "fatal" )   _minLogSeverity = elFatal;
    else                             _minLogSeverity = elWarning;
  }
  else {
    _minLogSeverity = elWarning;
  }

  key = "verbosity";
  std::string verbosity;
  if( opt->getOpt( tag, key, verbosity ) ) {
    if( verbosity      == "low" )     _verbosity =  0;
    else if( verbosity == "normal" )  _verbosity =  2;
    else if( verbosity == "high" )    _verbosity =  3;    
    else                              _verbosity =  0;
  }
  else {
    _verbosity = 0;
  }

  key = "store level";
  std::string storelevel;
  if( opt->getOpt( tag, key, storelevel ) ) {
    if( storelevel      == "debug" )   _minStoreSeverity = elDebugging;
    else if( storelevel == "verbose" ) _minStoreSeverity = elVerbose;
    else if( storelevel == "info" )    _minStoreSeverity = elInfo;
    else if( storelevel == "anomaly" ) _minStoreSeverity = elAnomaly;
    else if( storelevel == "warning" ) _minStoreSeverity = elWarning;
    else if( storelevel == "error" )   _minStoreSeverity = elError;
    else if( storelevel == "fatal" )   _minStoreSeverity = elFatal;
    else if( storelevel == "none" )    _minStoreSeverity = elBasicInfo;
    else                               _minStoreSeverity = elBasicInfo;
  }
  else {
    _minStoreSeverity = elBasicInfo;
  }

  std::string a;
  if(      _minLogSeverity == elDebugging ) a = "DEBUG";
  else if( _minLogSeverity == elVerbose )   a = "VERBOSE";
  else if( _minLogSeverity == elInfo )      a = "INFO";
  else if( _minLogSeverity == elAnomaly )   a = "ANOMALY";
  else if( _minLogSeverity == elWarning )   a = "WARNING";
  else if( _minLogSeverity == elError )     a = "ERROR";
  else if( _minLogSeverity == elFatal )     a = "FATAL";
  else if( _minLogSeverity == elBasicInfo ) a = "BASIC INFO";
  std::string b;
  if(      _minStoreSeverity == elDebugging ) b = "DEBUG";
  else if( _minStoreSeverity == elVerbose )   b = "VERBOSE";
  else if( _minStoreSeverity == elInfo )      b = "INFO";
  else if( _minStoreSeverity == elAnomaly )   b = "ANOMALY";
  else if( _minStoreSeverity == elWarning )   b = "WARNING";
  else if( _minStoreSeverity == elError )     b = "ERROR";
  else if( _minStoreSeverity == elFatal )     b = "FATAL";
  else if( _minStoreSeverity == elBasicInfo ) b = "BASIC INFO";
  std::string c;
  if(      _verbosity == 0 )  c =  "LOW";
  else if( _verbosity == 2 )  c =  "MEDIUM";
  else if( _verbosity == 2 )  c =  "HIGH";    


  std::cout << std::endl 
       << " ---------------------------------------------------- " << std::endl
       << "                 Coral Error Logger " << std::endl
       << " Minimum Log Severity:   " << a << std::endl
       << " Minimum Store Severity: " << b << std::endl
       << " Verbosity:              " << c << std::endl
       << " ---------------------------------------------------- " << std::endl
       << std::endl;

}

void CsErrLog::msg( Severity severity, std::string message, 
		    std::string source, int code ) {
  logger log;
  log._severity = severity;
  log._source   = source;
  log._code     = code;
  log._date     = CsTime();
  log._msg      = message;

  if( severity >= _minStoreSeverity ) {
    _log.push_back( log );
  }

  // send to cout or cerr (depending on the severity)
  if( severity >= _minLogSeverity ) {
    if( log._severity < elWarning || log._severity == elBasicInfo ) 
      _omsg( std::cout, log );
    else 
      _omsg( std::cerr, log ); 
  }
  
  // if fatal dump warnings, errors and fatal and abort.
  // do the same if error and _exitonerror = true.
  if( severity == elFatal ||
      ( severity == elError && _exitonerror ) ) {
    std::cout << std::endl << std::endl 
	 << "A FATAL ERROR APPEARED" << std::endl
	 << "The program will exit after logger dump." << std::endl;
    dump( elWarning );
    exit(1);
  }
}

void CsErrLog::msg( Severity severity, const char *source, int code, const char *format,...)
{
  char message[10000];
  va_list ap;
  va_start(ap,format);
  vsprintf(message, format, ap);
  msg(severity,message,source,code);
  va_end(ap);
}

void CsErrLog::dump() {
  dump( _minLogSeverity );
}
  
void CsErrLog::dump( Severity minSeverity ) {

  std::string severity;
  if( minSeverity == elDebugging )      severity = "DEBUG";
  else if( minSeverity == elVerbose )   severity = "VERBOSE";
  else if( minSeverity == elInfo )      severity = "INFO";
  else if( minSeverity == elAnomaly )   severity = "ANOMALY";
  else if( minSeverity == elWarning )   severity = "WARNING";
  else if( minSeverity == elError )     severity = "ERROR";
  else if( minSeverity == elFatal )     severity = "FATAL";
  else if( minSeverity == elBasicInfo ) severity = "NOTA BENE";

  std::cout << std::endl << std::endl
       << "Logger Dump, severity above " << severity << std::endl
       << "---------------------------------------" << std::endl;
  std::vector<logger>::iterator i;
  for( i=_log.begin(); i!=_log.end(); i++ ) {
    if( (*i)._severity >= minSeverity ) {
      _omsg( std::cout, (*i), 1 );
    }
  }
  std::cout << "End of Logger Dump." << std::endl;
}

void CsErrLog::_omsg( std::ostream& str, logger log ) {
  _omsg( str, log, _verbosity );
}

void CsErrLog::_omsg( std::ostream& str, logger log, int verbosity ) {

  std::string severity;
  if( log._severity == elDebugging )      severity = "DEBUG";
  else if( log._severity == elVerbose )   severity = "VERBOSE";
  else if( log._severity == elInfo )      severity = "INFO";
  else if( log._severity == elAnomaly )   severity = "ANOMALY";
  else if( log._severity == elWarning )   severity = "WARNING";
  else if( log._severity == elError )     severity = "ERROR";
  else if( log._severity == elFatal )     severity = "FATAL";
  else if( log._severity == elBasicInfo ) severity = "NOTA BENE";
  
  if( verbosity == 0 ) {
    str << severity 
	<< ": "
	<< "`" << log._msg << "'" << std::endl;
  }
  else if( verbosity == 1 ) {
    str << severity
	<< ", " << dmy << log._date 
	<< ", " << log._source 
	<< ", " << log._code  
	<< ": " << log._msg << std::endl;
  }
  else if( verbosity == 2 ) {
    str << std::endl
	<< severity
	<< ", on " << dmy << log._date 
	<<std::endl 
	<< " from: " << log._source 
	<< "   " << log._code  
	<<std::endl 
	<< "`" << log._msg << "'" << std::endl;
  }
  else {
    str << std::endl
	<< " Severity level: " << severity << std::endl
	<< " Date:           " << dmy << log._date << std::endl 
	<< " File/Facility:  " << log._source << std::endl
	<< " Line/Code:      " << log._code <<std::endl 
	<< " Message:" << std::endl
	<< "`" << log._msg << "'" << std::endl;
  }
}
