/*!
   \file    CsTrigHodDbUpdater.h
   \brief   Compass Some  CDB updater class definition file.
   \author  Vassili Motchalov
   \version $ $
   \date    $ $
*/
//-*-Mode: C++;-*-
#ifndef CsTrigHodDbUpdater_h
#define CsTrigHodDbUpdater_h

#if __GNUG__ >= 2
#  pragma interface
#endif

#include "CsSTD.h"
#include "CsTime.h"
#include "CsCondDbUpdater.h"
#include "CsTrigHodConstants.h"
/*!
	\class	CsTrigHodDbUpdater
	\brief	Compass Some CDB updater class.
*/
class CsTrigHodDbUpdater : public CsCondDbUpdater {
public:
  
  /*!
    \fn	CsTrigHodDbUpdater( string dbName ="trh"); 
    \brief	default constructor
    */
  /*    CsTrigHodDbUpdater(const string& dbName = "ec1"); OLD */
  CsTrigHodDbUpdater( string dbName ="trh"); 
  
  ~CsTrigHodDbUpdater(); 

  /*!
    \fn		bool update(const string& Container_name,  const CsTime& startTime, const CsTime& endTime );

    \param	startTime is a CsTime.
    \param	endTime is a CsTime.
    \return	a boolen
    \brief	store Some calibration constants information in BMS CDB with the given time interval.
    */ 

  /* bool update_all( const CsTime& startTime, const CsTime& endTime );
    bool update_all( CsCalibAbstract& TempConst,  const CsTime& startTime, const CsTime& endTime );
   bool update(const string& Container_name,CsCalibAbstract& TempConst,  const CsTime& startTime, const CsTime& endTime );
   */
  bool update(const string& Container_name,  const CsTime& startTime, const CsTime& endTime );
private:
  /*bool updateContainer(const string& Container_name, const int length, double * ptr,  const CsTime& startTime, const CsTime& endTime );
   */
  bool checkType(const CsCalibAbstract& TempConst);
};

#endif




