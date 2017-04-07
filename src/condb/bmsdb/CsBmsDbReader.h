/*!
   \file    CsBmsDbReader.h
   \brief   Compass BMS CDB reader class definition file.
   \author  Yoshiyuki Miyachi
   \version $Revision: 1.5 $
   \date    $Date: 2000/06/06 09:48:35 $
*/

//-*-Mode: C++;-*-
#ifndef _CsBmsDbReader_h_
#define _CsBmsDbReader_h_

#if __GNUG__ >= 2
#  pragma interface
#endif

#include "CsBMSconstants.h"
#include "CsCondDbReader.h"
#include "CsInterval.h"

//---- CsBmsDbReader -----------------------------------------------------------
/*!
  \class	CsBmsDbReader
  \brief	BMS calibration constants DB Reader Object.
*/
class CsBmsDbReader : public CsCondDbReader {

public:
/*!
  \fn		CsBmsDbReader(const string& dbName = "bms");
  \brief	default constructor. Default name of DB is defined as "bms".
*/
  CsBmsDbReader(const string& dbName = "bms");

  ~CsBmsDbReader(); //!< destructor

/*!
  \fn		CsBMSconstants getConstants(const CsTime& runTime); [public]
  \param	runTime is a CsTime
  \return	a CsBMSconstants
  \brief	get CsBMSconstants object from DB.
*/
  CsBMSconstants getConstants(const CsTime& runTime);

/*!
  \fn		list<CsBMSconstants> getAllConstants(const CsTime& runTime) [public]
  \param	runTime is a CsTime
  \return	a CsBMSconstants
  \brief	check all container which hold CsBMSconstants object and get them.
*/
  list<CsBMSconstants> getAllConstants(const CsTime& runTime);

/*!
  \fn		map<CsInterval, CsBMSconstants> getAllConstants(
			const CsTime& beginTime, const CsTime& endTime) [public]
  \param	beginTime is a CsTime
  \param	endTime is a CsTime
  \return	map<CsInterval, CsBMSconstants>
  \brief	check all container which hold CsBMSconstants object and get them.
*/
  map<CsInterval, CsBMSconstants> getAllConstants(
	  const CsTime& beginTime, const CsTime& endTime);

/*!
  \struct	ObjInfo
  \brief	This contains object identifieres (OID) 
			for container, interval, and persistent objects which 
			have been found by the reader.
*/
  struct ObjInfo {

	string containerName;	//!< container name
	string containerID;		//!< container id
	string intervalID;		//!< interval id
	string objectID;		//!< object id

  //!	equal to operator
	bool operator ==(const ObjInfo& objInfo){
	  return (	
		( containerName == objInfo.containerName )	&&
		( containerID	== objInfo.containerID )	&&
		( intervalID 	== objInfo.intervalID )		&&
		( objectID 		== objInfo.objectID )			); 
	}

  //! less than operator
	bool operator <(const ObjInfo& objInfo){
	  return containerName < objInfo.containerName ;
	}

  //! out put operator
	friend ostream& operator<<( ostream& os, const ObjInfo& objInfo );

  };

private:

/*!
  \var		list<ObjInfo> objInfos_ [private]
  \brief	list of Object Information for read object from DB.
*/
  list<ObjInfo> objInfos_;
  
/*!
  \fn		CsBMSconstants getFrom(	
			const CsCondDbConf::Container& container, 
			const CsTime& beginTime, const CsTime& endTime); [private]
  \param	container is a CsCondDbConf::Container
  \param	beginTime is a CsTime
  \param	endTime is a CsTime
  \return	CsBMSconstants
  \brief	Get CsBMSconstants from the given container.
*/
  map< CsInterval, CsBMSconstants > getFrom(	
	  const CsCondDbConf::Container& container, 
	  const CsTime& beginTime, const CsTime& endTime);

/*!
  \fn		CsBMSconstants getFrom(	
			const CsCondDbConf::Container& container, 
			const CsTime& runTime); [private]
  \param	container is a CsCondDbConf::Container
  \param	runTime is a CsTime
  \return	CsBMSconstants
  \brief	Get CsBMSconstants from the given container.
*/
  CsBMSconstants getFrom(	const CsCondDbConf::Container& container, 
						  const CsTime& runTime);

/*!
  \fn		bool getPobj(	HepRef( calibInterval ) pInterval, 
				  CsBMSconstants& bmsConst );
  \param	pInterval is a HepRef( calibInterval )
  \param	bmsConst is a CsBMSconstants
  \return	a bool
  \brief	get persistent object from the given interval and 
			set information of bmsConst using found persistent object.
			If everything are succeed, return true.
*/
  bool getPobj(	HepRef( calibInterval ) pInterval, 
				  CsBMSconstants& bmsConst );

/*!
  \fn		void dumpObjInfo(ostream& os);
  \param	os is a ostream
  \brief	dump contents in found object list to os.
*/
  void dumpObjInfo(ostream& os);

};

#endif
