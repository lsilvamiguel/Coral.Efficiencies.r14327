#if __GNUG__ >= 2
#  pragma implementation
#endif
 
#include "CsSTD.h"
#include "CsOpt.h"
/*!
   \file    CsTrigHodDbUpdater.cc
   \brief   Compass Some CDB updater class implementation file.
   \author  Vassili Motchalov
   \version $Revision: 1.1 $
   \date    $  $
*/
#include "CsErrLog.h"
#include "CsOpt.h"
#include "CsTrigHodDbUpdater.h"
#include "CsCalibAbstract.h"
#include "CsTrigHodConstants.h"

CsTrigHodDbUpdater::CsTrigHodDbUpdater(string  dbName) :
  CsCondDbUpdater(dbName) {}


CsTrigHodDbUpdater::~CsTrigHodDbUpdater() {
  
  CsErrLog::Instance()->mes(elDebugging, "CsTrigHodDbUpdater is going to die.");
  
}


bool CsTrigHodDbUpdater::update(const string& Container_remark,
				      const CsTime& startTime, 
				      const CsTime& endTime)  {

  if(!findContainer(Container_remark)) return false; 
  CsTrigHodConstants TempConst;
  if (!TempConst.readCoeff(container,  datafile.c_str()) ||
      !updateContainer(container,
		       TempConst.getLength(container),
		       TempConst.getCoeff(container),
		       startTime,   endTime))
    return false;  
  return true; 
}
 
bool CsTrigHodDbUpdater::checkType(const CsCalibAbstract& TempConst) {
  return true;
}




