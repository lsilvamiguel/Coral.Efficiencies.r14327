#if __GNUG__ >= 2
#  pragma implementation
#endif

/*!
   \file    CsBmsDbReader.cc
   \brief   Compass BMS CDB reader class implementation file.
   \author  Yoshiyuki Miyachi
   \version $Revision: 1.5 $
   \date    $Date: 2000/06/06 09:48:35 $
*/

#include "CsErrLog.h"
#include "CsOpt.h"
#include "CsBmsDbReader.h"
#include "CsBMSpConstants.h"

CsBmsDbReader::CsBmsDbReader(const string& dbName) : 
	CsCondDbReader(dbName) {
	this->init();
}


CsBmsDbReader::~CsBmsDbReader(){
//	CsErrLog::Instance()->mes(elDebugging, "Reader is going to die.");
	
// DUMP container information.
	ofstream objInfoFile("bmsdb.out");	
	this->dumpObjInfo(objInfoFile);

}

map<CsInterval, CsBMSconstants> CsBmsDbReader::getAllConstants(
		const CsTime& beginTime, const CsTime& endTime){
// get container list by tagging remark data
	list<CsCondDbConf::Container> 
		containerList = configuration_.findByRemark("BMSconstants");
	list<CsCondDbConf::Container>::iterator iCont;


//	open DB.
	this->open();
	
	
//	list<CsBMSconstants> lBmsConsts;
//	Fill list
//	for(iCont = containerList.begin(); iCont != containerList.end(); iCont++){
//
//		list<CsBMSconstants> foundConstant = 
//			this->getFrom( (*iCont), beginTime,  endTime);
//		lBmsConsts.insert(lBmsConsts.end(),  foundConstant.begin(), foundConstant.end());
//
//	}

	
	map<CsInterval, CsBMSconstants> lBmsConsts;
	if( containerList.size()>0 ){
		map<CsInterval, CsBMSconstants> foundConstant = 
			this->getFrom( (*containerList.begin()), beginTime,  endTime);

		lBmsConsts.insert( foundConstant.begin(), foundConstant.end());
	}

//	close DB												
	this->close();
	return lBmsConsts;

}


list<CsBMSconstants> CsBmsDbReader::getAllConstants(const CsTime& runTime){

// get container list by tagging remark data
	list<CsCondDbConf::Container> 
		containerList = configuration_.findByRemark("BMSconstants");

	list<CsBMSconstants> lBmsConsts;
	list<CsCondDbConf::Container>::iterator iCont;

//	open DB.
	this->open();

//	Fill list
	for(iCont = containerList.begin(); iCont != containerList.end(); iCont++){
		lBmsConsts.push_back(this->getFrom( (*iCont), runTime ));
	}

//	close DB												
	this->close();

	return lBmsConsts;
}

map<CsInterval, CsBMSconstants> CsBmsDbReader::getFrom(
	const CsCondDbConf::Container& container,
	const CsTime& beginTime, const CsTime& endTime){


//	list<CsBMSconstants> lBmsConsts;
	map<CsInterval, CsBMSconstants> lBmsConsts;
	
	HepTime bTime = this->convertTime(beginTime);
	HepTime eTime = this->convertTime(endTime);

// find interval in the container
	HepRef( calibInterval ) pInterval;
	
	int res = this->getDatabase()->findInterval(
				pInterval, container.name.c_str(), bTime );

	if( res ){

	// find object
		do{

			CsBMSconstants constant;
			if( this->getPobj(pInterval, constant) ){
//				lBmsConsts.push_back( constant );

				CsInterval newInterval(
					this->convertTime(pInterval->beginTime()),
					this->convertTime(pInterval->endTime()));
//				newInterval.beginTime	= 
//					this->convertTime(pInterval->beginTime());
//				newInterval.endTime		= 
//					this->convertTime(pInterval->endTime());
//
				lBmsConsts[newInterval] = constant;

			}

			pInterval = pInterval->next();

		} while ( (pInterval != 0) && (pInterval->beginTime() < eTime) );
		

	} else {	// in case interval is not found.
		CsErrLog::Instance()->mes(elError, 
			"CsBmsDbReader::getFrom\tcould not find interval in " 
			+ container.name);
	}
	
	return lBmsConsts;

}

bool CsBmsDbReader::getPobj(
  HepRef( calibInterval ) pInterval, CsBMSconstants& bmsConst ){
  bool status(false);


  HepRef(HepContObj) refCont = pInterval.containedIn();

  HepRef(CsBMSpConstants) pPerConst = 
	(HepRef(CsBMSpConstants)) pInterval->getObject();


  if( pPerConst ){

	// add objects information to list

	ObjInfo newObjInfo;
	newObjInfo.containerName	= refCont.name();
	newObjInfo.containerID		= refCont.sprint();
	newObjInfo.intervalID		= pInterval.sprint();
	newObjInfo.objectID		 	= pPerConst.sprint();

	objInfos_.push_back(newObjInfo);
	status = true;

	bmsConst = CsBMSconstants( *pPerConst );


  } else {	// in case object is not found.
  
	ostrstream out;
	out << "Persistent object was not found in ["
		<< refCont.name()
		<< "] while"
		<< pInterval 
		<< endl;
	CsErrLog::Instance()->mes(elError, out.str());

  }

  return status;
}



CsBMSconstants CsBmsDbReader::getFrom(
	const CsCondDbConf::Container& container,
	const CsTime& runTime){
	
//	Create CsBmsConstants object to return.
	CsBMSconstants bmsConst;
	HepTime rTime = this->convertTime(runTime);

// find interval in the container
	HepRef( calibInterval ) pInterval;
	
	int res = this->getDatabase()->findInterval(
				pInterval, container.name.c_str(), rTime );
	if( res ){
	// find object
		CsBMSconstants constant;
		if( this->getPobj(pInterval, constant) ){
			bmsConst = constant;
		}


	} else {	// in case interval is not found.
		CsErrLog::Instance()->mes(elError, 
			"CsBmsDbReader::getFrom\tcould not find interval in " 
			+ container.name);
	}

	return bmsConst;
}


CsBMSconstants CsBmsDbReader::getConstants(const CsTime& runTime){

// get container list by tagging remark data
	list<CsCondDbConf::Container> 
		containerList = configuration_.findByRemark("BMSconstants");

	if(containerList.size() == 0) {
		CsErrLog::Instance()->mes( elError, 
			"No container defined in configuration file." );
		return CsBMSconstants();
	}

	ostrstream out;
	out	<< configuration_ << endl;
	CsErrLog::Instance()->mes(elDebugging, out.str() );

	if(containerList.size() != 1){
		CsErrLog::Instance()->mes( elInfo, 
			"More than one BMSconstants containers are defined." );		
	}

//	open DB.
	this->open();

//	Create CsBmsConstants object to return.
	CsBMSconstants bmsConst = 
		this->getFrom(	*( containerList.begin() ), runTime );

//	close DB												
	this->close();

	return bmsConst;
}

void CsBmsDbReader::dumpObjInfo(ostream& os){

	os	<< "#############################################################\n"
		<< "#[cdb name]\t[container]\t[container ID]\t[inteval ID]\t[object ID]" 
		<< "#\n"
		<< endl;
	list<CsBmsDbReader::ObjInfo>::iterator iObjInfo;
	for( iObjInfo = objInfos_.begin(); iObjInfo != objInfos_.end() ; iObjInfo++ )
		os	<< this->getDbName() << '\t' 
			<< (*iObjInfo)
			<< endl;

}


ostream& operator<<( ostream& os, const CsBmsDbReader::ObjInfo& objInfo ){
	os	<< objInfo.containerName	<< '\t'
		<< objInfo.containerID		<< '\t'
		<< objInfo.intervalID		<< '\t'
		<< objInfo.objectID ;
	return os;
}
