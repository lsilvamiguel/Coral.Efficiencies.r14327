/*!
  \file		CsGeoDbUpdater.cc
  \brief	Geometry Database ("geo") handler for updating implementation file.
  \author	$Author: miyachi $
  \version	$Revision: 1.3 $
  \date		$Date: 2000/06/21 00:10:52 $
*/
#if __GNUG__ >= 2
#  pragma implementation
#endif

#include "CsErrLog.h"
#include "CsGeoDbUpdater.h"
#include "GdbPerStation.h"
#include "GdbDBuilder.h"

CsGeoDbUpdater::CsGeoDbUpdater(const string& dbName) :
	CsCondDbUpdater(dbName) {

	// initialize Transaction Manager Don't forget it
	this->init();
	
	string mes = "CsGeoDbUpdater object has been constructed in ";
	manager_->getMode() == READ ? mes += "READ " : mes += "UPDATE ";
	mes += "Mode";
	
	CsErrLog::Instance()->mes(elDebugging, mes);

}


CsGeoDbUpdater::~CsGeoDbUpdater() {
	CsErrLog::Instance()->mes(elDebugging, "CsGeoDbUpdater is going to die.");
}

bool CsGeoDbUpdater::update(	
	const CsTime& startTime, const CsTime& endTime) {

//	create transient Geometry volumes 
//		using FFR cards, g3call file, and detectors.dat
	GdbDBuilder& builder = *( GdbDBuilder::Instance( GdbDBuilder::G3FILE ) );

	builder.build();


//	Convert CsTime to HepTime
	HepTime sTime = this->convertTime( startTime );
	HepTime eTime = this->convertTime( endTime );


//	Open DB
	this->open();

	list<GdbDetector>::iterator iDetList;
	
	for(	iDetList =  builder.detectors().begin();
		iDetList != builder.detectors().end();
		iDetList++){
		
		int nS = (*iDetList).nStation();
		
		for(	int i = 0; i < nS ; i++){
//		
//		
//			iStaList =  (*iDetList).station().begin();
//			iStaList !=  (*iDetList).station().end();
//			iStaList++){
  
			GdbPerStation newPerStation = 
			  GdbPerStation( (*iDetList).station(i) );

			HepRef(GdbPerStation) refDetector =
				new( this->getDatabase()->hint() ) 
					GdbPerStation(  newPerStation );

			string containerName = 
				CsGeoDbUpdater::makeContainerName( 
						(*iDetList).name(),
//						refDetector->volumeName(), 
						i, 
						refDetector->gVersion() );

			CsCondDbConf::Container newContainer = 
				configuration_.findByName(containerName);

			if( newContainer.name == "!" ){
				// update configuration file entry
				newContainer.name 		= containerName;
				newContainer.objectName	= "GdbPerStation";
				newContainer.dataSource	= refDetector->gVersion();
				newContainer.remark		= (*iDetList).name();
				configuration_.lContainer().push_back( newContainer );

			}

			// store to DB
			this->getDatabase()->store( 
				refDetector,
				newContainer.name.c_str(),
				sTime, eTime
			);
			
		}

	}

	ofstream outfile("GEODB.conf");
	outfile << configuration_ << endl;

// set status flag to HAVE
	status_=HAVE;

// close DB
	this->close();

	return true;
}

string CsGeoDbUpdater::makeContainerName( 
	const string volName, const int& i , const string gversion ){

	ostrstream strID;

	strID << volName;
	strID << "-";
	strID << int(i + 1);
	strID << "-";
	strID << gversion << '\0';

	return strID.str();
}
