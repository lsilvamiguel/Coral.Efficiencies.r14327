#ifndef _CsGeoDbReader_h_
#define _CsGeoDbReader_h_
/*!
  \file		CsGeoDbReader.h
  \brief	Geometry Database ("geo") handler for reading definition file.
  \author	$Author: miyachi $
  \version	$Revision: 1.2 $
  \date		$Date: 2000/06/07 08:36:43 $
*/

#if __GNUG__ >= 2
#  pragma interface
#endif

#include "CsCondDbReader.h"
#include "CsGeoDbDetector.h"
/*!
  \class	CsGeoDbReader
  \brief	Geometry DB Reader Object DB.
*/
class CsGeoDbReader : public CsCondDbReader {

public:
/*!
  \fn		CsGeoDbReader(const string& dbName = "geo");
  \param	dbName is a string, which has a default parameter, "geo".
  \brief	a default constructor. 
			Normally, one dose not need to specify DB name.
			Only when CDB developpers want to create own test DB, 
			the optional argument dose have meaning. 
			The given name must be equal to the name in configuration file.
  
*/
    CsGeoDbReader(const string& dbName = "geo");

    ~CsGeoDbReader(); //!< destructor

	void check(ostream& os);	//!< test routine, which will go sometime..

/*!
  \fn		list<CsGeoDbDetector*> getDetectorTable( 
			  const CsTime& runTime, const string& gversion );
  \param	runTime is a CsTime
  \param	gversion is a string
  \return	a list<CsGeoDbDetector*>
  \brief	Read Geometry data of gversion at runTime from CDB.
			And return detector table information found in CDB.
*/
	list<CsGeoDbDetector*> getDetectorTable( 
		const CsTime& runTime, const string& gversion );

//! dump information to ostream
	friend ostream& operator<<(ostream& os, CsGeoDbReader& reader);

private:
/*!
  \fn		void getAllConstants( const CsTime& runTime, const string& gversion );
  \param	runTime is a CsTime
  \param	gversion is a string
  \brief	Try to find Geometry data, which has the given version, 
			at the given time (runTime). 
			The found data are holded by this handler and 
			one can access by getDetectorTable() method of this object.
*/
	void getAllConstants( const CsTime& runTime, const string& gversion );

/*!
  \fn		bool checkType(const CsCalibAbstract& TempConst)
  \param 	Container_name  Name of the Container to save    
  \brief	Check the type of  CsCalibAbstrac object
  \warning	Not ready yet - should be defined for each class
*/
	virtual bool checkType(const CsCalibAbstract& TempConst){return true;}

};

#endif
