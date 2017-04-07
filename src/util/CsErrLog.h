// $Id: CsErrLog.h,v 1.10 2003/03/13 16:57:31 benigno Exp $

/*!
   \file    CsErrLog.h
   \brief   Compass Error and Message Logger.
   \author  Benigno Gobbo 
   \version $Revision: 1.10 $
   \date    $Date: 2003/03/13 16:57:31 $
*/

#ifndef CsErrLog_h
#define CsErrLog_h

#include "CsSTD.h"
#include "CsTime.h"

//! \enum Severity
//! \brief defines the severity level of the message. 
enum Severity {
  elDebugging = -2, //!< Max information level.
  elVerbose   = -1, //!< More infos that the default.
  elInfo      = 0,  //!< Normal procedural messages.
  elAnomaly,        //!< Something anomalous. But understood... 
  elWarning,        //!< Some unexpected behaviour. Result could be wrong.
  elError,          //!< Requests or action not produced. Severe.
  elFatal,          //!< No way to continue. The program must end.
  elBasicInfo = 100 //!< Mandatory info
};

//! \struct logger
//! \brief defines all quantities to be logged.
typedef struct { 
  Severity _severity; //!< severity level (see Severity).
  std::string   _source;   //!< source of the message (defaults to file).
  int      _code;     //!< code (defaults to line).
  CsTime   _date;     //!< time.
  std::string   _msg;      //!< message.
} logger;

/*! \class CsErrLog
    \brief Compass Error and Message Logger. 
*/

class CsErrLog {

 public:

  /*! \fn static CsErrLog* Instance();
      \brief singleton instantiation.
  */
  static CsErrLog* Instance();

  /*! \fn void msg(Severity severity, string message,string source, int code); 
      \brief Store a message to the logger and output it to cout or cerr
      depending on the severity. It is preferable to access it using the
      macro (in this case source will be set as \c __FILE__  and code will 
      be set as \c __LINE__): 
      <p><b> mes( Severity severity, string message ); </b>
      \param severity severity level.
      \param message the error message.
      \param source of the message. 
      \param code message code. 
      \param message the related message.
  */

  static void msg( Severity severity, std::string message, std::string source, int code ); 

  /*! \brief Error message with arguments. See above function.
       Example: CsErrLog->Instance()->msg(elWarning,__FILE__,__LINE__,"Bad argument: %s",arg);
  */
  static void msg( Severity severity, const char *source, int code, const char *message,...);

# ifdef mess
# undef mess
# endif
# define mes( sev, message ) msg( sev, message, __FILE__, __LINE__ )
 
  /*! \fn void dump();
      \brief Dumps all messages with severity >= minLogLevel_.
   */
  static void dump( void );

  /*! \fn static void dump ( Severity minSeverity );
      \brief Dumps all messages with severity >= minSeverity
      \param minSeverity the value of Severity above which the message
      must be dumped.
   */
  static void dump( Severity minSeverity );

  //! clear the logged data
  static void clear( void ) { _log.clear(); }

  //! set exit in case of Error level message
  static void setExitOnError( void ) { _exitonerror = true; }

  //! unset exit in case of Error level message
  static void unsetExitOnError( void ) { _exitonerror = false; }

 protected:

  CsErrLog();                    //!< Default Constructor.

 private:

  static CsErrLog* _instance;         //!< The singleton static pointer.
  static Severity         _minLogSeverity;   //!< Minimum logging level.
  static Severity         _minStoreSeverity; //!< Minimum level to be stored.
  static int              _verbosity;        //!< Verbosity of messages.
  static std::vector<logger>   _log;              //!< List of logged messages.
  static bool             _exitonerror;      //!< Exit if level=elError

  /*! \fn static void _omsg( std::ostream& str, logger log, int verbosity );    
    \brief Write the message to ostream str.
    \param str reference to ostream (cout or cerr).
    \param log logger structure to be output.
    \param verbosity verbosity of messages.
  */
  static void _omsg( std::ostream& str, logger log, int verbosity );    

  /*! \fn static void _omsg( std::ostream& str, logger log );    
    \brief Write the message to ostream str.
    \param str reference to ostream (cout or cerr).
    \param log logger structure to be output.
  */
  static void _omsg( std::ostream& str, logger log );    

};

#endif // CsErrLog_h
