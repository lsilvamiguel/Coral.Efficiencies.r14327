/*!
  \file		GdbDetector.cc
  \brief	Geometry DB detector object implementation file
  \author	$Author: miyachi $
  \version	$Revision: 1.3 $
  \date		$Date: 2000/06/21 00:10:52 $
*/
#if __GNUG__ >= 2
#  pragma implementation
#endif

#include "CsErrLog.h"
#include "CsOpt.h"
#include "GdbDetector.h"
#include "GdbDBuilder.h"
#include "CsCOMGEANT.h"

//---- GdbDetector ---------------------------------------------------------
//	char name_[16];				//!< detecotor name
//	GdbDetectorType type_;					//!< detector type id
//	list<GdbStation> station_;	//!< Station Information
//const int GdbDetector::maxNameLength = 8;

GdbDetector::GdbDetector() :
	id_(0), 
	type_(), 
	station_(), 
	table_(),
	idMap_()
	{
	this->name( "--------" );
}

GdbDetector::GdbDetector(const int& id, const char* pName, const int& unit) :
	id_(id),
	type_(), 
	station_(unit),
	table_(),
	idMap_()
	{
	this->name( pName );
}

GdbDetector::GdbDetector( const GdbDetector& detector ) :
	id_(			detector.id()),
	type_(			detector.type()), 
	station_(		detector.stationCopy()),
	table_(			detector.table()),
	idMap_(			detector.idMap())
	{
	this->name( detector.name() );
}


GdbDetector::~GdbDetector(){
}


bool GdbDetector::createIdMap(){
	string tag = ""; 
	string key = "detector table";
	string detectorsTable;
	if( ! CsOpt::Instance()->getOpt( tag, key, detectorsTable ) ){
		CsErrLog::Instance()->mes(elError, "No detectors.dat have been specified.");
		return false;
	}
	
	ifstream detTablStream( detectorsTable.c_str() );
	if(detTablStream.fail()){
		CsErrLog::Instance()->mes(elError, "Fail to open detectors.dat");
		return false;
	}

	const int lineSize = 2048;
	char line[lineSize];

	do{
		if( detTablStream.eof() ) continue;
		detTablStream.getline( line, lineSize, '\n' );
		istrstream s(line);
		string opt;
		s >> opt;
	
		if(opt == "det"){
			string name;
			int id, unit;
			s >> id >> name >> unit;
			
			PlaneId pi(name, unit); 
			
			idMap_[pi] = id;
			
		}	
	
	} while( !detTablStream.eof() );
	/*	
	ofstream detMapFile("detectorID.map");
	map<PlaneId, int>::iterator iM;
	detMapFile	<< "ID\t[NAME]\tunit\t"	<< endl;
	for( iM = idMap_.begin(); iM != idMap_.end(); iM++ )
	detMapFile	<< (*iM).second	<< '\t'
				<< "["	<< (*iM).first.name	<< "]\t"
				<< (*iM).first.unit	<< '\t'
				
				<< endl;
	*/
	return true;
}

int GdbDetector::findIdInMap(const PlaneId& pid){
	map<PlaneId, int>::iterator iM;
//
//	cerr << "Try to find :[" << pid.name << "] " << pid.unit << endl;
//
	for(	iM = idMap_.begin(); iM != idMap_.end(); iM++ ){
		string	name = (*iM).first.name;
		int		unit = (*iM).first.unit;
//	
//		cerr	<< "MAP ENTRY:[" << name << "] " << unit << "\t"
//				<< "ID: "	<< (*iM).second
//				<< endl;
//				
//				
		if( name == pid.name ){
//			cerr << "\tCheck unit" << unit << endl;
			if( unit == pid.unit ) return (*iM).second;
		}	

	}
	return 0;
}


const char* GdbDetector::name( const char* pName ){
	return strncpy(name_, pName, GdbDetector::maxNameLength);
}

bool GdbDetector::createDetectorTable(){

	if( this->refTable().size() > 0 ) return true;


	for(	int iSt = 0; iSt < this->nStation(); iSt++){
	
		list<GdbVolume::ActivePlane> activePlane;
		GdbStation* pStation = &( this->station()[iSt] );


//		ostrstream debugBuff;
//		debugBuff 	<< "NAME:["	<< pStation->volumeName() << "]\t"
//					<< "UNIT:["	<< iSt	<< "]"	
//					<< endl;
//
//		CsErrLog::Instance()->mes( elDebugging, debugBuff.str() );


//		cerr 	<< "NAME:["	<< pStation->volumeName() << "]\t"
//				<< "UNIT:["	<< iSt	<< "]"	
//				<< endl;


		pStation->findActivePlane( activePlane );

	
		list<GdbVolume::ActivePlane>::iterator iA;
		for(	iA  = activePlane.begin();
			iA != activePlane.end();
			iA++){

//			GdbVolume* pV = (*iA).pVolume;
			GdbVolume* pV;
			
			if( (*iA).index < 0 ) {
				pV = dynamic_cast<GdbVolume*>( pStation );
			} else {
//				pV = &( pStation->refVolumes()[ (*iA).index ] );
				pV = &( pStation->volume( (*iA).index ) );
			}

			if( pV->wire().number() >0 ){ // ingnore station.
				
				const char* vName = pV->volumeName();
				const int   unit  = iSt + 1;

				// Try to find detId
				if( idMap_.size() == 0) this->createIdMap();

				PlaneId pid(vName, unit);
				int detId;
				
//				detId = idMap_[ pid ];
//				
//				if(detId ==0) 
				detId = this->findIdInMap(pid);
				
//				cerr	<< "Found ID:"	<< detId
//						<< "\tFor ["	<< vName << "]\t"
//						<< unit
//						<< endl;

				if( pV->wire().id() != detId ) pV->wire().id(detId);


				// Get center point vector and rotation matrix in MRS
				GdbPoint3D center = 
					pStation->center() +
					pStation->rotMatrix() * (*iA).center;

				GdbRotation rotation =
					pStation->rotMatrix() * (*iA).rotation;

				// Get Solid parameter
				// Define Solid Type Information
				string solid = pV->solid().type();
				long nParam = pV->solid().nParam();
				vector<double> param(nParam, 0.0);
				for( long i = 0 ; i< nParam; i++) 
					param[i] = pV->solid().param(i);
			
//				HepPoint3D size(0.0, 0.0, 0.0);
//				if( solid == "BOX" ) {
//					size =	HepPoint3D(
//								pV->solid().param(0),
//								pV->solid().param(1),
//								pV->solid().param(2)
//							);
//
//				}


				// assign detecot type
				const GdbDetectorType& detType = pStation->detectorType();

				GdbDetTableCont newDet(
					detId,
					vName,
					unit,
					detType.type(),
					pV->material().radLen(),
//					size,
					COMGEANT::solidId( pV->solid().type() ),
					param,
					GdbDBuilder::convert( center ),
					HepMatrix( GdbDBuilder::convert( rotation ) ),
					pV->wire().firstPosition(),
					pV->wire().angle(),
					pV->wire().number(),
					pV->wire().pitch(),
					detType.eff(),
					detType.bkg(),
					detType.gate()
				);


				if(detType.hasDrift()){
					newDet.setDriftInfo(
						detType.driftV(), detType.T0(), 
						detType.dHitRes(), detType.sRes(), detType.tSlic()
					);
				}

		
//////////// Now only the case the acitve plane containes one dead volume is available.		
				if( pV->nSubVolume() == 1 ){

//					GdbVolume::VolumeContainer lSV= pV->subVolume();
//					GdbVolume dV = *(lSV.begin());
					long nextIndex = (*iA).index + 1;

//					GdbVolume& dV = pStation->refVolumes()[ nextIndex ];
					GdbVolume& dV = pStation->volume( nextIndex );
					
					GdbSolidType sT = dV.solid();
					long nP = sT.nParam();
				
					GdbPoint3D dCenter = center + rotation * dV.center();
					GdbRotation dRotation = rotation * dV.rotMatrix();

					vector<double> par(nP, 0.0);
					for( long i = 0 ; i< nP; i++) par[i] = sT.param(i);

					newDet.deadZone(
						COMGEANT::solidId( sT.type() ),
						par, 
						GdbDBuilder::convert( dCenter ), 
						HepMatrix( GdbDBuilder::convert( dRotation ))
					);

				}

				this->refTable().push_back( newDet );
			}
		}
	}
	return true;
}

void GdbDetector::dumpDetectorTable(ostream& os){

	list<GdbDetTableCont>::iterator iC;
	for(	iC =  this->refTable().begin();
		iC != this->refTable().end();
		iC++){

//		(*iC).dumpDetTableHeader(os);
		os	<< (*iC)	<< endl;
	}	

}

void GdbDetector::dumpVolumeTree(ostream& os){
	for(	int iSt = 0; iSt < this->nStation(); iSt++){
		os	<< "#----------- Station:\t"	<< "Unit:" << iSt + 1 << endl;
		os << this->station()[iSt] << endl;
	}
}

ostream& operator<<(ostream& os, GdbDetector& detector){
	os	<< "#################################################\n"
		<< "#\t"
		<< "GdbDetector:\t" 
		<< "[ID: "	<< detector.id()	<< "\t"
		<< "NAME: ("	<< detector.name()	<< ")\t"
		<< "Unit: "	<< detector.nStation() 	<< "\t"
		<< "]"		<< endl;

	os	<< "#---- Volume Tree" << endl; 
	detector.dumpVolumeTree(os);

	os	<< "#---- Detector Table" << endl; 
	detector.dumpDetectorTable(os);
	
	os	<< "\n#################################################" << endl;

	return os;
}

void GdbDetector::dump(ostream& os){

  os	<< "detector"	<< endl;
  os	<< "\t"	<< "name"		<< "\t"	<< this->name()	<< endl;
  os	<< "\t"	<< "type"		<< "\t"	<< this->type()	<< endl;

  int nSta = this->nStation();
  if( nSta > 0){
	  os	<< "\t"	<< "station" ;
	  for(int i = 0; i< nSta ; i++) 
		  os  << "\t"	<< this->station(i).volumeName();
	  os	<< endl;
  }

//	os	<< "\t"	<< "gate"		<< "\t"	<<
  os	<< "end detector" << endl;

  this->dumpVolumeTree(os);
  
}

void GdbDetector::dumpXML(ostream& os){

  const char* pTab = "\t";

  os	<< "<detector>"	<< endl;
  os	<< *pTab	<< "<name>"	<< this->name()	<< "</name>"	<< endl;
  os	<< *pTab	<< "<type>"	<< this->type()	<< "</type>"	<< endl;

  int nSta = this->nStation();
  for(int i = 0; i< nSta ; i++) this->station(i).dumpXML(os, pTab );

  os	<< "</detector>" << endl;

  this->dumpVolumeTree(os);
  
}


void GdbDetector::createDetectorFile(){
	string outputFile = this->name();
	outputFile += ".gdat";
	
	ofstream outFileStream(outputFile.c_str());
//	this->dump(outFileStream);
	this->dumpXML(outFileStream);
	
	outFileStream << "#-------------- end ------------------" << endl;
}
