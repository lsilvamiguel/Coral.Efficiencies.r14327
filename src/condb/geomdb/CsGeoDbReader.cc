/*!
  \file		CsGeoDbReader.cc
  \brief	Geometry Database ("geo") handler for reading implementation file.
  \author	$Author: miyachi $
  \version	$Revision: 1.2 $
  \date		$Date: 2000/06/07 08:36:42 $
*/
#if __GNUG__ >= 2
#  pragma implementation
#endif

# ifdef COMPASS_USE_OSPACE_STD
#  include <ospace/std/algorithm>
# else
#  include <algorithm>
# endif // COMPASS_USE_OSPACE_STD


#include "CsErrLog.h"
#include "CsGeoDbReader.h"
#include "GdbDBuilder.h"
#include "CsGeoDbUpdater.h"

CsGeoDbReader::CsGeoDbReader(const string& dbName) : 
	CsCondDbReader(dbName) {
	
	// initialize Transaction Manager Don't forget it
	this->init();

	string mes = "CsGeoDbReader object has been constructed in ";
	manager_->getMode() == READ ? mes += "READ " : mes += "UPDATE ";
	mes += "Mode";

	CsErrLog::Instance()->mes(elDebugging, mes);

}

CsGeoDbReader::~CsGeoDbReader(){
}

void CsGeoDbReader::getAllConstants(	
		const CsTime& runTime, const string& gversion ){

	GdbDBuilder* pBuild = GdbDBuilder::Instance( GdbDBuilder::GEODB );

	if( pBuild->isUpdated() ) return ;
	
//	Try to find Containers of gversion
	CsErrLog::Instance()->mes(elDebugging, "Start to read geometry DB");

	list<CsCondDbConf::Container> 
			containerList = configuration_.findBySource(gversion);

	list<CsCondDbConf::Container>::iterator iCont;

	map<string, int> detectorMap;

	for(	iCont = containerList.begin(); 
		iCont != containerList.end(); 
		iCont++){

		string detectorName = (*iCont).remark;
		
		if(	detectorMap.find(detectorName) == detectorMap.end() ){
			detectorMap[detectorName] = 0;
		}

		detectorMap[detectorName] += 1;

	}

	
//      open DB.
	this->open();
	HepTime rTime = this->convertTime(runTime);

	map<string, int>::iterator itrMap;
	int detId = 0;
	for(	itrMap =  detectorMap.begin();
		itrMap != detectorMap.end();
		itrMap++){
		// Detector loop start
		GdbDetector& newDetector = 
			pBuild->add( 
				GdbDetector( detId++, (*itrMap).first.c_str(), (*itrMap).second )
			);

		int nStation = (*itrMap).second;
		for(int i =0; i< nStation; i++){
			// station loop start

			string containerName = 
				CsGeoDbUpdater::makeContainerName( 
						(*itrMap).first, i, gversion );

			CsCondDbConf::Container foundContainer = 
				configuration_.findByName(containerName);

			if( foundContainer.name == "!" ){
				CsErrLog::Instance()->mes( elError,
					"Can not find container:" + containerName);
			} else {

//				cerr	<< "Found container["
//						<< foundContainer.name << "]"
//						<< endl;
//// find interval in the container

				HepRef( calibInterval ) pInterval;
				if( this->getDatabase()->findInterval(  
					pInterval, foundContainer.name.c_str(), rTime)){

				// find object
					HepRef( GdbPerStation ) pPerConst = 
							(HepRef(GdbPerStation)) pInterval->getObject();

					if( pPerConst ){

						newDetector.station(i) = GdbStation( *pPerConst );

					} else {        // in case object is not found.
						ostrstream out;
						out << "CsGeoDbReader::getFrom\tcould not find persistent object in "
							<< foundContainer.name 
							<< " while " << pInterval << endl;
						CsErrLog::Instance()->mes(elError, out.str());
					}

				} else {        // in case interval is not found.
					CsErrLog::Instance()->mes(elError, 
							"CsGeoDbReader::getFrom\tcould not find interval in " 
							+ foundContainer.name);
				}


			}
			
			// statiton loop end
		}
		

//		cerr	<< "Try to create detector table in "
//				<< newDetector.name()
//				<< " system."
//				<< endl;
//				
		newDetector.createDetectorTable();
		
		// detector loop end
	}

	//      close DB
	this->close();

	pBuild->isUpdated(true);

}

//GdbStation CsGeoDbReader::getFrom(
//		const CsCondDbConf::Container& container,
//        const CsTime& runTime){
//
//	GdbStation station;
//	HepTime rTime = this->convertTime(runTime);
//
//// find interval in the container
//	HepRef( calibInterval ) pInterval;
//	if( this->getDatabase()->findInterval(  
//		pInterval, container.name.c_str(), rTime)){
//
//	// find object
//		HepRef( GdbPerStation ) pPerConst = 
//				(HepRef(GdbPerStation)) pInterval->getObject();
//
//		if( pPerConst ){
//
//				station = GdbStation( *pPerConst );
//
//		} else {        // in case object is not found.
//				ostrstream out;
//				out << "CsGeoDbReader::getFrom\tcould not find persistent object in "
//					<< container.name << " while " << pInterval << endl;
//				CsErrLog::Instance()->mes(elError, out.str());
//		}
//
//	} else {        // in case interval is not found.
//			CsErrLog::Instance()->mes(elError, 
//					"CsGeoDbReader::getFrom\tcould not find interval in " 
//					+ container.name);
//	}
//
//	return station;
//}

void CsGeoDbReader::check(ostream& os){

	// create geoIndex.DB
	
	string mainDbName = this->getDbName();
	mainDbName += "Index.DB";
	
	os << "Try to find[" << mainDbName << "]" << endl;

	ooHandle(ooDBObj)& mainDb = 
		manager_->db( mainDbName.c_str() ) ;
	
	os << "DB name:\t["	<< mainDb.name()		<< "]"	<< endl;
	os << "DB path:\t["	<< mainDb.pathName()	<< "]"	<< endl;
	
	ooItr( ooContObj ) contI;
	
	mainDb.contains( contI );
	os	<< "Container list in " << mainDb.name()	<< "]\t";
	do{

		os	<< "name:["	<< contI.name()	<< "]\t";
		os	<< "Type Number:["	<< contI.typeN()	<< "]\t";
		os	<< "Type Name:["	<< contI.typeName()	<< "]\t";
		os	<< endl;
		
	} while( contI.next());

}

list<CsGeoDbDetector*> CsGeoDbReader::getDetectorTable( 
	const CsTime& runTime, const string& gversion ){

	list<CsGeoDbDetector*> newList;

	GdbDBuilder* pBuild = 
		GdbDBuilder::Instance( GdbDBuilder::GEODB );

	if( ! pBuild->isUpdated() ) this->getAllConstants(runTime, gversion);

	list<GdbDetector>::iterator itrDet;
	for(	itrDet  = pBuild->detectors().begin();
		itrDet != pBuild->detectors().end();
		itrDet++){
		// detector loop
		
//		cerr	<< "Detecot:"
//				<< (*itrDet).name()
//				<< "\thas "
//				<< (*itrDet).refTable().size()
//				<< " active area."
//				<< endl;

		list<GdbDetTableCont>::iterator itrDetTable;
		for(	itrDetTable =  (*itrDet).refTable().begin();
			itrDetTable != (*itrDet).refTable().end();
			itrDetTable++){
			// detector table loop
			
			newList.push_back( &(*itrDetTable) );
			
			// end detector table loop
		}
		// end detector loop
	}

	return newList;

}

ostream& operator<<(ostream& os, CsGeoDbReader& reader){
	GdbDBuilder* pBuild = GdbDBuilder::Instance( GdbDBuilder::GEODB );

	if( ! pBuild->isUpdated() ) {
		os	<< "This has not been ready..... You should call getAllConstants() method."
			<< endl;

		reader.check(os);


	} else {

		os	<< "================ DB structure dump =====================" << endl;
		os	<< (*pBuild) << endl;

	}

	return os;
}
