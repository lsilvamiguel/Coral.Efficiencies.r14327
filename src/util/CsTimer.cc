// $Id: CsTimer.cc,v 1.6 2010/01/28 12:51:26 tnagel Exp $

/*!
  \file    CsTimer.cc
  \brief   CsTimer class
  \author  Massimo Lamanna
  \version $Revision: 1.6 $
  \date    $Date: 2010/01/28 12:51:26 $
*/  

#include "CsTimer.h"

#ifdef COMPASS_USE_OSPACE_STD
#  include <ospace/std/algorithm>
#  include <ospace/std/iostream>
#else
#  include <algorithm>
#  include <iostream>
#endif

#include <time.h>
#include <sys/timeb.h>

using namespace std;

CsTimer* CsTimer::instance_ = 0;
time_t         CsTimer::baseTime(0);
unsigned short CsTimer::millisec(0);

CsTimer::CsTimer() {

#ifdef _WIN32
  timeb_ time;
  _ftime( &time );
#else
  timeb time;
  ftime( &time );
#endif

  baseTime = time.time;
  millisec = time.millitm;

}

CsTimer* CsTimer::Instance() {
  if( instance_ ==  0 ) {
    instance_ = new CsTimer();
  }
  return( instance_ );
}

void CsTimer::start(string label) {

#ifdef CS_DUMP

  bool dump(false);

  if(getenv("CSTIMER_DUMP")!=NULL) {
    dump = true;
  }
  
  dump && cout << "CsTimer::start dump version" << endl;

#endif  

#ifdef _WIN32
  struct _timeb time;
  _ftime( &time );
#else
  struct timeb time;
  ftime( &time );
#endif

  timerMap[label] = -( time.time - baseTime + (time.millitm - millisec)/1000. ); // Negative

#ifdef CS_DUMP
  dump && cout << "CsTime::start " << baseTime << " s from epoch (baseTime)" << endl;
  dump && cout << "CsTime::start " << millisec << " ms correction (baseTime)" << endl;
  dump && cout << "CsTime::start " << time.time << " s from epoch" << endl;
  dump && cout << "CsTime::start " << time.millitm << " ms correction" << endl;
  dump && cout << "CsTime::start " << timerMap[label] << " s start time * (-1)" << endl;
  dump && cout << "CsTime::start " << label << " label" << endl;
#endif //CS_DUMP

}

void CsTimer::stop(string label) {

#ifdef CS_DUMP

  bool dump(false);

  if(getenv("CSTIMER_DUMP")!=NULL) {
    dump = true;
  }
  
  dump && cout << "CsTimer::stop dump version" << endl;

#endif //CS_DUMP

  if(timerMap.find(label) != timerMap.end()) {

#ifdef _WIN32
    struct _timeb time1,time2;
    _ftime( &time2 );
#else
    struct timeb time1,time2;
    ftime( &time2 );
#endif

    timerMap[label] = (time2.time -baseTime) + (time2.millitm-millisec)/1000. + timerMap[label];

#ifdef CS_DUMP
    dump && cout << "CsTime::stop " << baseTime << " s from epoch (baseTime)" << endl;
    dump && cout << "CsTime::stop " << millisec << " ms correction (baseTime)" << endl;
    dump && cout << "CsTime::stop " << time2.time << " s from epoch" << endl;
    dump && cout << "CsTime::stop " << time2.millitm << " ms correction" << endl;
    dump && cout << "CsTime::stop " << timerMap[label] << " s timer" << endl;
    dump && cout << "CsTime::stop " << label << " label" << endl;
#endif //CS_DUMP

  }
#ifdef CS_DUMP
  else {
    dump && cout << "CsTime::stop no valid quantity (label) " << label << endl;
  }
#endif //CS_DUMP
  
}

float CsTimer::gimmeDeltaT(string label) {

  if(timerMap.find(label) != timerMap.end()) {
    return(timerMap[label]);
  }
  else {
    return(-1);
  }
}

bool CsTimer::isValid(string label) {

#ifdef CS_DUMP

  bool dump(false);

  if(getenv("CSTIMER_DUMP")!=NULL) {
    dump = true;
  }
  
  dump && cout << "CsTimer::bool dump version" << endl;
#endif //CS_DUMP

  if(timerMap.find(label) != timerMap.end()) {

#ifdef CS_DUMP
    dump && cout << "CsTimer::bool no start for label " << label << endl;
#endif //CS_DUMP

    return (false);

  }

  if(timerMap[label]<0) {
#ifdef CS_DUMP
    dump && cout << "CsTimer::bool no stop for label " << label << endl;
#endif //CS_DUMP

    return (false);
  }

  return (true);

}


// multiple sting CsTimer::start (3,4,5 arguments)

void CsTimer::start(string label1,string label2) {

  start(label1);
  start(label2);

}

void CsTimer::start(string label1,string label2,string label3) {

  start(label1);
  start(label2);
  start(label3);

}

void CsTimer::start(string label1,string label2,string label3, string label4) {

  start(label1);
  start(label2);
  start(label3);
  start(label4);

}

void CsTimer::start(string label1,string label2,string label3, string label4, string label5) {

  start(label1);
  start(label2);
  start(label3);
  start(label4);
  start(label5);

}
