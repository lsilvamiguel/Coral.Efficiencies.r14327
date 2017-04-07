#ifndef CsTimer_h
#define CsTimer_h

/*!
   \file CsTimer.h
   \brief Engine to time events
   \version 0.0
   \author  Massimo Lamanna
   \date    18 May 1999

   ver. 0.1 10 May 1999 - Massimo Lamanna
*/

#include <CsPlatforms.h>
#include <CsTypes.h>

#include <CsSTD.h>
#include <CsTime.h>

#include <map>

/*! \class CsTimer
    \brief Engine to time events

    This is a singleton object which keep track
    of all the quantities to be timed in a program.

    These quantities are identified by string,
    which can have descriptive values like
   
    string elapsedTime("Elapsed time");

    Each of them has to be started at the
    appropriate place with the method
    start(string). This method accept up to
    5 strings (i.e. 5 independent quantities to be 
    with a common start).

    To stop, a stop(string) method is provided.

    The value of the time difference for each
    string is kept in the CsTimer object and
    it can be tested (method isValid(string)
    which is true only if a start and a stop have
    been issue for a given string. The value of
    each time difference can be retreeved
    with the getDeltaT(string) method which
    returns a float.

    Multiple starts (before a stop) reset the start time.
    Multiple stops are discarded.

    Problems: internally timeb is used; on Linux
    is gives a 1 s resolution.

*/

class CsTimer  {

public:

  static CsTimer* Instance();

  void start(std::string);
  void start(std::string,std::string);
  void start(std::string,std::string,std::string);
  void start(std::string,std::string,std::string,std::string);
  void start(std::string,std::string,std::string,std::string,std::string);
  void stop(std::string);
  bool isValid(std::string);
  float gimmeDeltaT(std::string);

 protected:

  CsTimer();

private:

  static CsTimer* instance_;

  static time_t         baseTime; // Time from epoch
  static unsigned short millisec; // Milliseconds

  std::map<std::string,float,std::less<std::string> > timerMap;

};

#endif // CsTimer_h
