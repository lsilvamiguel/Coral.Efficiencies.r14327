/*!
  \file		GdbDBuilder.cc
  \brief	Geometry Data builder object implementation file
  \author	$Author: miyachi $
  \version	$Revision: 1.2 $
  \date		$Date: 2000/06/07 08:36:44 $
*/
#if __GNUG__ >= 2
#  pragma implementation
#endif

#include "GdbDBuilder.h"
#include "CsSTD.h"
#include "CsRegistry.h"
#include "CsErrLog.h"
#include "CsGeom.h"
#include <cmath>
#include "FfrReader.h"
#include "CsG3CallFile.h"

#include "CsCOMGEANT.h"


//---- GdbDBuilder ---------------------------------------------------------

GdbDBuilder* GdbDBuilder::instance_ = NULL ;

GdbDBuilder* GdbDBuilder::Instance(){
	if(instance_ == NULL){
		CsErrLog::Instance()->mes(elFatal, 
			"Wrong GdbDBuilder instanciation, you must specify the mode.");;		
	}
	return instance_;
}

bool GdbDBuilder::end(){
	return true;
}

GdbDBuilder* GdbDBuilder::Instance( const GdbDBuilder::MODE& mode ){
	if( instance_ != NULL ){
		if( mode != instance_->mode() ){
			CsErrLog::Instance()->mes(elFatal, 
				"Wrong GdbDBuilder instanciation, you given wrong mode flag.");;		
			return instance_;
		}
		return instance_;
	}
	
	instance_ = new GdbDBuilder(mode);
	CsRegistry csReg_;
	if( csReg_.EOJRegistration( (CsEndOfJob*) instance_ ) ){
		CsErrLog::Instance()->mes(elDebugging, 
			"GdbDBuilder has been registerd successfully.");;
	} else {
		CsErrLog::Instance()->mes(elError, 
			"GdbDBuilder has not been registerd successfully.");;
	}
	
	return instance_;
}

GdbDBuilder::GdbDBuilder( const GdbDBuilder::MODE& mode ) :
	mode_(mode), detectors_(), updated_(false) {
}


GdbDBuilder::~GdbDBuilder(){
}

bool GdbDBuilder::build(){

//	CsErrLog::Instance()->mes( elDebugging, 
//		"Detector structure construction is going to start.");

	if(this->mode() == GdbDBuilder::G3FILE ){

		if( ! this->detector() ){
			CsErrLog::Instance()->mes(elError, "FAIL");
		}

		this->isUpdated(true);
		
	} else {

		


	}

//	CsErrLog::Instance()->mes( elDebugging, 
//		"Detector structure construction has been finished.");

	return true;

}

bool GdbDBuilder::detector(){
// Get Detector Information from FFR cards.
	ostrstream debugLog;

	for( int id = 4; id < 21; id++ ){
		if( this->getDetectorFromFFR(id) ){
			debugLog << "Detector(" << id << ") has been found in FFR." << endl;
		} else {
			debugLog << "Detector(" << id << ") has not been found in FFR." << endl;
		}
	}

// Fill Detector information using G3CallFile
	if(  this->detectors().size() == 0 ){
		CsErrLog::Instance()->mes(elError, "No detector information were found.");
		return false;
	}

	for(	itrDetector iDet =  this->detectors().begin();
		iDet != this->detectors().end();
		iDet++){
// Create volume tree in station
		this->fillDetectorByG3Call( (*iDet) );

// Create Detector table in detecot
		(*iDet).createDetectorTable();
		
	}


// Create RICH volumes
	this->rich();

// Muon dependent detectors
	if(FfrReader::Instance()->muon() != ""){
		this->muonFilter();
	}
	

// Create Passive Volumes as GdbDetector
	this->CreatePassiveVolumes();


	for(	itrDetector iDet =  this->detectors().begin();
		iDet != this->detectors().end();
		iDet++){
		debugLog	<< "Detector:["
					<< (*iDet).name()
					<< "] has been successfull created."
					<< endl;
	}
	
	CsErrLog::Instance()->mes(elDebugging, debugLog.str() );

	return true;
}

bool GdbDBuilder::CreatePassiveVolumes(){

	string dName, vTop;
//	Magnet
//	SM1
	dName = "SM1"; vTop = "SM1M";
	this->passiveVolume(dName, vTop);
	
//	SM2
	dName = "SM2"; vTop = "MAG2";
	this->passiveVolume(dName, vTop);
	
//	Walls
	dName = "WALLS"; vTop = "WALL";
	this->passiveVolume(dName, vTop);

//	Muon Wall 1
	dName = "MW1W"; vTop = "MW1W";
	this->passiveVolume(dName, vTop);

//	HCAL1
	dName = "HCAL1"; vTop = "HC1";
	this->passiveVolume(dName, vTop);

//	HCAL2
	dName = "HCAL2"; vTop = "HC2";
	this->passiveVolume(dName, vTop);

//	HODOSCOPS
	dName = "HOD1"; vTop = "HOD1";
	this->passiveVolume(dName, vTop);

//	HODOSCOPS
	dName = "HOD2"; vTop = "HOD2";
	this->passiveVolume(dName, vTop);

//	HODOSCOPS
	dName = "HOD3"; vTop = "HOD3";
	this->passiveVolume(dName, vTop);

//	HODOSCOPS
	dName = "HOD4"; vTop = "HOD4";
	this->passiveVolume(dName, vTop);

//	HODOSCOPS
	dName = "H4VS"; vTop = "H4VS";
	this->passiveVolume(dName, vTop);


//	ECAL1
	dName = "ECAL1"; vTop = "LG1";
	this->passiveVolume(dName, vTop);

//	ECAL1
	dName = "ECAL2"; vTop = "LG2";
	this->passiveVolume(dName, vTop);

//	Muon Filter 1
	dName = "MU1F"; vTop = "MU1F";
	this->passiveVolume(dName, vTop);

//	Target
	dName = "TAGT"; vTop = "PTSU";
	if(FfrReader::Instance()->muon() != ""){
		vTop = "PTSU";
	}
	this->passiveVolume(dName, vTop);

	return true;
}

bool GdbDBuilder::passiveVolume( 
	const string& detectorName, const string& name ){

	int stationUnit(1);
	int detectorId(0);

// create detector
	this->detectors().push_back( 
		GdbDetector( detectorId, detectorName.c_str(), stationUnit ) );

	itrDetector iNewDet = this->detectors().end(); iNewDet--;

// Create station
	// set serch starting point
	Gspos* pSerchStartPoiint = &( CsG3CallFile::Instance()->HALL() );
	bool findResult(false);

	Gspos& gspos = 
		pSerchStartPoiint->findPosition( name, 1 , findResult );

	if(! findResult ){
		CsErrLog::Instance()->mes(elError, 
			name + "was not found in volume tree.");
		return false;
	}

//	Set Version tag to station.
	(*iNewDet).station()[0].gVersion( FfrReader::Instance()->gVersion() );

// Get Geometry Information in MRS.
	(*iNewDet).station()[0].center(
		GdbDBuilder::convert(
			CsG3CallFile::Instance()->centerInHALL( gspos )
		)
	);

	(*iNewDet).station()[0].rotMatrix(
		GdbDBuilder::convert(
			CsG3CallFile::Instance()->rotationToHALL( gspos )
		)
	);

// create subvolumes
	this->setVolumeInformation( (*iNewDet).station()[0], gspos );
	this->createVolumesIn( (*iNewDet).station()[0], gspos );

	return true;
}



bool GdbDBuilder::rich(){
	string dName("RIC1GASM"), volumeTop("GASM");
// mirror ?
	dName = "RIC1GASM"; volumeTop = "GASM";
	this->passiveVolume(dName, volumeTop);
// upperside photon detector ?
	dName = "RIC1UPDT"; volumeTop = "UPDT";
	this->passiveVolume(dName, volumeTop);
// bottom side photon detector ?
	dName = "RIC1DNDT"; volumeTop = "DNDT";
	this->passiveVolume(dName, volumeTop);
	

	return true;
}


bool GdbDBuilder::muonFilter(){
// Muon Filter 2
	string detectorName("MU2F");
	int detectorId(40);

// Try to find station information
	int stationUnit(0);
	vector<string> stationBuff;
	string sbuff;
	
	while( ( sbuff = this->findMF2Station( stationUnit + 1 ) ) != "" ){
		FfrReader::Instance()->strip( sbuff, FfrReader::nDelim );
		stationBuff.push_back( sbuff );
		stationUnit++;		
	}
	
	ostrstream out;
	out	<< "Station found:\t" << stationUnit << endl;
	CsErrLog::Instance()->mes(elDebugging, out.str());

	if( stationUnit == 0 ) {
		CsErrLog::Instance()->mes(elError, 
			"No station for MUF2 could not be found.");
		return false;
	}

// create detector
	this->detectors().push_back( 
		GdbDetector( detectorId, detectorName.c_str(), stationUnit ) 
	);

	itrDetector iNewDet = this->detectors().end(); iNewDet--;

// Common Gate Length
	sbuff = FfrReader::Instance()->dataLine()[ "MUF2DEFI01" ];
	istrstream inputStream(sbuff.c_str());
	
	double	tGate(0.0), dTmp,  eff(1.0), bkg(0.0), 
			driftV(0.0), t0(0.0), dRes(0.0), sRes(0.0), tSli(0.0);
			
	inputStream >> tGate 
				>> dTmp 
				>> dTmp 
				>> eff >> bkg 
				>> driftV >> t0
				>> dRes >> sRes >> tSli ;

	GdbDetectorType 
		detType( (*iNewDet).name(), 11, tGate, eff, bkg);
	
	detType.setDriftInfo(driftV, t0, dRes, sRes, tSli);
	
	(*iNewDet).type(detType);

// create station
	for( int iSta =0; iSta < stationUnit; iSta++ ){
		istrstream iss( stationBuff[iSta].c_str() );

		Gspos* pSerchStartPoiint = 
			&( CsG3CallFile::Instance()->HALL() );
		
		bool findResult(false);
		string name; iss >> name ;
		Gspos& gspos = 
			pSerchStartPoiint->findPosition( name, 1 , findResult );

		if(! findResult ){
			CsErrLog::Instance()->mes(elError, 
				name + "was not found in volume tree.");
			return false;
		}

//	Set Version tag to station.
		(*iNewDet).station()[iSta].gVersion( 
			FfrReader::Instance()->gVersion() );

// Get Geometry Information in MRS.
		(*iNewDet).station()[iSta].center(
			GdbDBuilder::convert(
				CsG3CallFile::Instance()->centerInHALL( gspos )
			)
		);

		(*iNewDet).station()[iSta].rotMatrix(
			GdbDBuilder::convert(
				CsG3CallFile::Instance()->rotationToHALL( gspos )
			)
		);

// copy detector information to detector type

		(*iNewDet).station()[iSta].detectorType( (*iNewDet).type() );

		this->setVolumeInformation( (*iNewDet).station()[iSta], gspos );

		this->createVolumesIn( (*iNewDet).station()[iSta], gspos );
		
	}
	return true;	
}


string GdbDBuilder::findMF2Station( const int& id ){
	ostrstream strID;
	strID << setfill( '0' );
	strID << setw(2) << id;
	string tag("MUF2STAT");
	return FfrReader::Instance()->dataLine()[ tag + strID.str() ];
}


bool GdbDBuilder::getDetectorFromFFR( int const& id ){
	if( id < 1 || id > 99 ){
		ostrstream out;
		out << "ID must be within range (1-99) :\t" << id << endl;
		CsErrLog::Instance()->mes(elError, out.str() );
		return false;
	}


	ostrstream strID;
	strID << setfill( '0' );
	strID << setw(2) << id;

	const int ntag_opt=2;
	vector<string> tag_opt( ntag_opt, "" );
	tag_opt[0] = "GSET";
	tag_opt[1] = "SET";

//---------------- Get Detector Name  ---------------------
	string tag(""), name("");
	for(int i =0; i < ntag_opt; i++) {
//		tag  = tag_opt[i] + strID.str();
		string detectorID = strID.str();
		tag  = tag_opt[i] + detectorID[0] + detectorID[1];
		name = FfrReader::Instance()->dataLine()[ tag + "NAME" ];
		if( name != "" ) break;
	}

	if( name == "" ){
		ostrstream out;
		out << "FFR KEY(" << tag + "NAME" << ") was not used in ffr card." << endl;
		CsErrLog::Instance()->mes(elDebugging, out.str() );
		return false;
	}
	
	FfrReader::Instance()->strip(name, FfrReader::space);
	FfrReader::Instance()->strip(name, FfrReader::nDelim);

//---------------- Get UNIT Information ---------------------
	string strBuff = FfrReader::Instance()->dataLine()[ tag + "UNITS" ];
	int unit = FfrReader::Instance()->stoi( strBuff );
	
//---------------- Create Detector Object -------------------
	this->detectors().push_back( GdbDetector( id, name.c_str(), unit ) );

//---------------- Obtain iterator to new GdbDetector object
	itrDetector iDet = this->detectors().end();
	iDet--;	


//---------------- Get Detector Type Information ---------------------
	strBuff = FfrReader::Instance()->dataLine()[ tag + "TYPE" ];
	int type = FfrReader::Instance()->stoi( strBuff );

//---------------- Get Time Gate Information ---------------------
	strBuff = FfrReader::Instance()->dataLine()[ tag + "GATE" ];
	double tGate = FfrReader::Instance()->stod( strBuff );

//---------------- Get Efficiency Information ---------------------
	strBuff = FfrReader::Instance()->dataLine()[ tag + "DETEFFI" ];
	int nPara(0); double eff(1.0);
	FfrReader::Instance()->multi(strBuff, nPara, eff );
	
//---------------- Get Background Information ---------------------
	strBuff = FfrReader::Instance()->dataLine()[ tag + "DETBACK" ];
	double bkg(0.0);
	FfrReader::Instance()->multi(strBuff, nPara, bkg );


	GdbDetectorType detType( (*iDet).name(), type, tGate, eff, bkg );

//----------------- Drift Information -----------------------------


	if( type == 11 ) { // straw
		double	velo(0.0), t0(0.0),
				dHitRes(0.0), sRes(0.0), tSlic(0.0);

		strBuff = FfrReader::Instance()->dataLine()[ tag + "DETVELO" ];
		FfrReader::Instance()->multi(strBuff, nPara, velo );
		
		strBuff = FfrReader::Instance()->dataLine()[ tag + "DETTIM0" ];
		FfrReader::Instance()->multi(strBuff, nPara, t0 );
		
		strBuff = FfrReader::Instance()->dataLine()[ tag + "DETDRES" ];
		FfrReader::Instance()->multi(strBuff, nPara, dHitRes );
		
		strBuff = FfrReader::Instance()->dataLine()[ tag + "DETSRES" ];
		FfrReader::Instance()->multi(strBuff, nPara, sRes );
		
		strBuff = FfrReader::Instance()->dataLine()[ tag + "DETSLIC" ];
		FfrReader::Instance()->multi(strBuff, nPara, tSlic );
		
		detType.setDriftInfo(velo, t0, dHitRes, sRes, tSlic);
		
	}

	(*iDet).type(detType); 

	return true;
}

bool GdbDBuilder::fillDetectorByG3Call(GdbDetector& detector){

	unsigned int nSta = detector.nStation(); // number of station.
	Gspos* pSerchStartPoiint = &( CsG3CallFile::Instance()->HALL() );

	
//	ostrstream out;
//	out << "Detector Creation for:\t" << detector.name() << endl;
//	CsErrLog::Instance()->mes( elDebugging, out.str() );

	for( unsigned int iSta = 0; iSta < nSta ; iSta++){
		bool findResult;

		ostrstream out;
		out << "Try to find [" 
			<< detector.name() << ", Unit " << iSta + 1 << "]"<< "\tin\n"
			<< pSerchStartPoiint->name()
			<< endl;

		Gspos& station = 
			pSerchStartPoiint->findPosition( 
				detector.name(), (iSta + 1) , findResult );
			
		if( ! findResult ){
			out << " has not been found in tree." << endl;
			CsErrLog::Instance()->mes( elError, out.str() );
		}

//		out	<< " has been found in tree as" << endl;
//		out << station << endl;
//		
//		CsErrLog::Instance()->mes( elDebugging, out.str() );


//	Set Version tag to station.
		detector.station()[iSta].gVersion( FfrReader::Instance()->gVersion() );

//	Each Staion must be stay in HALL system
		detector.station()[iSta].center( 
			GdbDBuilder::convert(
				CsG3CallFile::Instance()->centerInHALL( station ) 
			)
		);
			
		detector.station()[iSta].rotMatrix(
			GdbDBuilder::convert(
				CsG3CallFile::Instance()->rotationToHALL( station )
			)
		);

// copy detector information to detector type
		detector.station()[iSta].detectorType( detector.type() );
		this->setVolumeInformation( detector.station()[iSta], station );
		this->createVolumesIn( detector.station()[iSta], station );

	}
	
	return true;
}

bool GdbDBuilder::createVolumesIn( GdbStation& station, Gspos& pos ){

//=============================================================================
//	Update Station Information.
//=============================================================================
//---- create subvolumes in volume
//---- 

	list<Gspos*> ptrPosList = pos.pDaughter();
	list<Gspos*>::iterator itrPtrPos;
	
	for(	itrPtrPos = ptrPosList.begin();
		itrPtrPos != ptrPosList.end();
		itrPtrPos++){

//---- Add subvolume to station
		station.addVolume(
			GdbVolume(	GdbDBuilder::convert( (*itrPtrPos)->center() ),
						GdbDBuilder::convert( (*itrPtrPos)->rM())    ) );
		
//		long nSub = volume.nSubVolume();
//		this->createVolumesIn( volume.refSubVolume( nSub -1 ), *(*itrPtrPos) );

//---- create daughter volumes in station
		long nSub = station.nVolume();
		this->setVolumeInformation( 
			station.volume( nSub -1 ), *(*itrPtrPos) );
//			station.refVolumes()[ nSub -1 ], *(*itrPtrPos) );
		this->createVolumesIn( station, *(*itrPtrPos) );

//		GdbVolume::VolumeContainer::iterator 
//			iNew = volume.refSubVolume().end(); iNew--;
//		this->createVolumesIn( (*iNew), *(*itrPtrPos) );

	}
		
	return true;
}

void GdbDBuilder::setVolumeInformation( GdbVolume& volume, Gspos& pos ){
//---- Get Volume Name information
  volume.volumeName( pos.name().c_str() );	// Name

//---- Get volume information
  Gsvolu& gsvol = *( pos.pVolume() );
  vector<double> param = gsvol.param();
	
  if(param.size() == 0){	// parametarized volume
	 param = pos.param();
  }

  GdbSolidType newSolid( gsvol.type().c_str(), param.size() );

  for( unsigned long iparam = 0; iparam < param.size(); iparam++ ){
	newSolid.addParam( iparam, param[iparam]);
  }

  volume.solid(newSolid);

//---- Get SubVolume Information
  volume.nSubVolume( pos.pDaughter().size() );


//---- Get Material information
  CsMaterial& material = *( gsvol.pMedium()->pMaterial() );

  GdbMaterial newMaterial( 
	material.name().c_str(), material.density(),
	material.radiationLength()
  );

  volume.material(newMaterial);

//---- Get Wire information.
  if( pos.pDetector() != NULL ) { // Yes this is active

	Gsdet& detInfo = *( pos.pDetector() );
	vector<double> detParam = detInfo.user();

	int wirN = 	int( detParam[Gsdet::WirN] );		
	double wirP = detParam[Gsdet::WirP];
	double firstPosition = detParam[Gsdet::WirD];

	double angle = 0.0;

	double cosAng = detParam[Gsdet::CosAng] ;
	double sinAng = detParam[Gsdet::SinAng] ;

	double factor = 180.0 / M_PI;
	if( cosAng > 0.0 ) angle =  asin( sinAng ) * factor ;
	else {
	  angle = acos( cosAng );
	  angle *=( sinAng > 0.0 ? factor  : -1.0*factor );
	}

	GdbWireInfo wireInfo( wirN, wirP, angle, firstPosition );

	volume.wire(wireInfo);
  }

}


ostream& operator<<( ostream& os, GdbDBuilder& builder ){

	if( ! builder.isUpdated() ){
		if(builder.mode() == GdbDBuilder::G3FILE ){
			CsErrLog::Instance()->mes(elInfo,
				"GdbDbuilder is not ready for use. You should call build() method.");
		} else {
			CsErrLog::Instance()->mes(elInfo,
				"GdbDbuilder is not ready for use. You should build detector before.");
		}		
		return os;
	}

	if(builder.mode() == GdbDBuilder::G3FILE )
		os	<< *( FfrReader::Instance() ) << endl;

	GdbDBuilder::itrDetector iDet;
	if(  builder.detectors().size() > 0 ){
		os	<< "#   -------------- Detector list -----------------" << endl;
		for(	iDet =  builder.detectors().begin();
			iDet != builder.detectors().end();
			iDet++){
			os	<< (*iDet) << endl;
		}
	} else {
		os	<< "Empty detector list." << endl;
	}
	return os;

}

void GdbDBuilder::makeDetectorFiles(){
	itrDetector iDet;
	for(	iDet = this->detectors().begin(); 
		iDet != this->detectors().end();
		iDet++)
		(*iDet).createDetectorFile();
}

GdbPoint3D GdbDBuilder::convert(const HepPoint3D& point) {
	return GdbPoint3D(point.x(), point.y(), point.z());
}

HepPoint3D GdbDBuilder::convert(const GdbPoint3D& point) {
	return HepPoint3D(point.x(), point.y(), point.z());
}

GdbRotation GdbDBuilder::convert(const CsRotation& rotM) {
	return GdbRotation(
		rotM.xx(), rotM.xy(), rotM.xz(),
		rotM.yx(), rotM.yy(), rotM.yz(),
		rotM.zx(), rotM.zy(), rotM.zz()
	);
}

CsRotation GdbDBuilder::convert(const GdbRotation& rotM) {
	return CsRotation(
		rotM.xx(), rotM.xy(), rotM.xz(),
		rotM.yx(), rotM.yy(), rotM.yz(),
		rotM.zx(), rotM.zy(), rotM.zz()
	);
}

GdbDetector& GdbDBuilder::add( const GdbDetector& detector ) {
	detectors_.push_back(detector) ;
	list<GdbDetector>::iterator iNewDet = detectors_.end();
	iNewDet--;
	return (*iNewDet);
}

