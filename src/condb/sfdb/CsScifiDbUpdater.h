#ifndef _CsScifiDBUpdater_h_
#define _CsScifiDBUpdater_h_

#include "CsCondDbUpdater.h"
#include "CsCalibScifi.h"

//---- CsScifiDBUpdater ----------------------------

class CsScifiDbUpdater : public CsCondDbUpdater
{
public:
  CsScifiDbUpdater(const string &dbName="sfd");
  virtual ~CsScifiDbUpdater();
  bool store(CsCalibScifi &scifi,CsTime &startTime,CsTime &endTime);	
};

#endif
