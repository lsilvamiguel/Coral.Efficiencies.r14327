/*!
   \file    CsBmsDbUpdater.h
   \brief   Compass BMS CDB updater class definition file.
   \author  Yoshiyuki Miyachi
   \version $Revision: 1.6 $
   \date    $Date: 2000/06/06 09:48:35 $
*/
//-*-Mode: C++;-*-
#ifndef _CsBmsDbUpdater_h_
#define _CsBmsDbUpdater_h_

#if __GNUG__ >= 2
#  pragma interface
#endif

#include "CsSTD.h"
#include "CsCondDbUpdater.h"
/*!
  \class	CsBmsDbUpdater
  \brief	Compass BMS CDB updater class.
*/
class CsBmsDbUpdater : public CsCondDbUpdater {
public:

/*!
  \fn		CsBmsDbUpdater(const string& dbName = "bms");
  \brief	default constructor
*/
    CsBmsDbUpdater(const string& dbName = "bms");

    ~CsBmsDbUpdater(); //!< destructor

/*!
  \fn bool updateAll();
  \return	a boolen
  \brief	store BMS calibration constants.
*/
  bool updateAll();


private:

};

#endif
