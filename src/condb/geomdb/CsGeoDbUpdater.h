#ifndef _CsGeoDbUpdater_h_
#define _CsGeoDbUpdater_h_
/*!
  \file		CsGeoDbUpdater.h
  \brief	Geometry Database ("geo") handler for updating definition file.
  \author	$Author: miyachi $
  \version	$Revision: 1.2 $
  \date		$Date: 2000/06/07 08:36:43 $
*/
#if __GNUG__ >= 2
#  pragma interface
#endif

#include "CsSTD.h"
#include "CsTime.h"
#include "CsCondDbUpdater.h"

//---- CsGeoDbUpdater -----------------------------------------------------------
/*!
  \class	CsGeoDbUpdater
  \brief	Geometry DB Updater Object DB.
*/
class CsGeoDbUpdater : public CsCondDbUpdater {
public:

/*!
  \fn		CsGeoDbUpdater(const string& dbName = "geo");
  \brief	default constructor
*/
  CsGeoDbUpdater(const string& dbName = "geo");

  ~CsGeoDbUpdater(); //!< destructor

/*!
  \fn		bool store(	const CsDetectorFile& detectorTable,
					  const CsTime& startTime, const CsTime& endTime);
  \param	detectorTable is a CsDetectorFile.
  \param	startTime is a CsTime.
  \param	endTime is a CsTime.
  \return	a boolen
  \brief	store geometry information in Geometry database using CsDetectorFile,
		  "detectors.dat", with the given time interval.
*/
  virtual bool update(const CsTime& startTime, const CsTime& endTime);

/*
  \fn		update
  \param	Container_name  Name of the Container to save    
  \param 	startTime  The beginning time of the Interval
  \param	endTime    The end time of the Interval
  \brief 	Update specific Container to  the DataBase
*/
  virtual bool update(
	const string& Container_name,  
	const CsTime& startTime, const CsTime& endTime ) {
	return true;
  }

/*!
  \fn		string makeContainerName( 
			const string detName, const int& i , const string gversion );
  \param	detName is a string
  \param	i is a int
  \param	gversion is a string
  \brief	define container name 
*/
  static string makeContainerName( 
	const string detName, const int& i , const string gversion );

private:
/*!
  \fn		bool checkType(const CsCalibAbstract& TempConst) 
  \param	Container_name  Name of the Container to save    
  \brief	Check the type of  CsCalibAbstrac object
  \warning	Not ready yet - specific for each class
*/
  virtual bool checkType(const CsCalibAbstract& TempConst) {
	return true;
  }

};

#endif
