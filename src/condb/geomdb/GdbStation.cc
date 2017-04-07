/*!
  \file		GdbStation.cc
  \brief	Geometry DB transient station object implementation file
  \author	$Author: miyachi $
  \version	$Revision: 1.2 $
  \date		$Date: 2000/06/07 08:36:47 $
*/
#if __GNUG__ >= 2
#  pragma implementation
#endif

#include "CsSTD.h"
#include "GdbStation.h"
#include "GdbDBuilder.h"
#include "CsErrLog.h"


//---- GdbStation ---------------------------------------------------------
GdbStation::GdbStation() : 
	GdbVolume(), 
	gVersion_(""), 
	detectorType_(), 
	references_(),
	volumes_(),
	volumeMap_()
	{
}

GdbStation::GdbStation(const GdbStation& station) : 
	GdbVolume(	station.volumeName(),
				station.solid(),
				station.material(),
				station.center(),
				station.rotMatrix(),
				station.nSubVolume(),
				station.isActive(),
				station.wire() ),
	gVersion_(""),
	detectorType_(station.detectorType()), 
	references_(station.references()),
	volumes_(station.volumes()),
	volumeMap_(station.volumeMap())
	{
}


GdbStation::GdbStation(GdbPerStation& perStation) : 
	GdbVolume(	perStation.volumeName(),
				perStation.solid(),
				perStation.material(),
				perStation.center(),
				perStation.rotMatrix(),
				perStation.nSubVolume(),
				perStation.isActive(),
				perStation.wire() ),
	gVersion_(perStation.gVersion()),
	detectorType_(perStation.detectorType()), 
	references_(),
	volumes_(),
	volumeMap_()
	{	
	long nR = perStation.referenceNumber();
	for(	long i=0; i< nR; i++){
		this->addReference( 
			GdbDBuilder::convert( perStation.referencePoint(i) ) );
	}
	
	long nV = perStation.nVolume();
	for( long i=0; i<nV; i++){
		this->addVolume( perStation.volume(i) );
	}
	
}
		
GdbStation::~GdbStation(){
}

GdbStation& GdbStation::operator=(const GdbStation& station){
	if( this != &station){

		this->volumeName(	station.volumeName()	);
		this->solid(		station.solid()			);
		this->material(		station.material()		); 
		this->center(		station.center()		);
		this->rotMatrix(	station.rotMatrix()		); 
		this->nSubVolume(	station.nSubVolume()	);
		this->isActive(		station.isActive()		); 
		this->wire(			station.wire()			);
		this->detectorType(	station.detectorType()	);
		this->references(	station.references()	);
		this->gVersion(		station.gVersion()		);
		this->volumes(		station.volumes()		);

	}
	return *this;
}

void GdbStation::dumpDetectorsTable(ostream& os){
	
}

ostream& operator<<(ostream& os, GdbStation& station){
//	station.dump(os);
	station.dumpXML(os);
	return os;
}

void GdbStation::dumpVolume(ostream& os, const char* pTab){
	for( long i = 0; i < this->nVolume(); i++){
		this->dumpVolume(os, i, pTab);
	}
}

void GdbStation::dumpVolume(ostream& os, const long& index, const char* pTab){
	string top("");
	if( pTab != NULL ) top = pTab ;

	string TAB = top + '\t';

	GdbVolume& volume = this->volume(index);

	os	<< top	<< "volume"	<< endl;
	os	<< TAB	<< "name"		<< "\t"	<< volume.volumeName() << endl;
	os	<< TAB	<< volume.solid() << endl;
	os	<< TAB	<< "position"	<< "\t"	<< volume.center()	<< endl;
	os	<< TAB	<< volume.rotMatrix() << endl;
	os	<< TAB	<< volume.material() << endl;


// dump subvolume name
	vector<long>& 
		vSubVolumeIndex = volumeMap_[index];

	long nSub = vSubVolumeIndex.size();
	if(nSub > 0) os	<< TAB	<< "volumes";
	for(	long i = 0; i < nSub; i++) {
		os	<< "\t"	
			<< this->volume( vSubVolumeIndex[i] ).volumeName();
	}
	os << endl;

	if(this->isActive()){
		os	<< TAB	<< volume.wire() << endl;		
	}

	os << top	<< "end volume" << endl;

}

void GdbStation::dump(ostream& os, const char* pTab){
	string top("");
	if( pTab != NULL ) top = pTab ;

	string TAB = top + '\t';
	os	<< "#"	
		<< top	<< "VERSION"	<< "\t"	<< this->gVersion()	<< endl;
	os	<< top	<< "station"	<< endl;
	os	<< TAB	<< "name"		<< "\t"	<< this->volumeName() << endl;
	os	<< TAB	<< this->solid() << endl;
	os	<< TAB	<< "position"	<< "\t"	<< this->center()	<< endl;
	os	<< TAB	<< this->rotMatrix() << endl;
	os	<< TAB	<< this->material() << endl;

	if( this->volumeMap().size() == 0) {
		long index(-1);
		this->createVolumeMap(index);		
	}

// dump subvolume name
	vector<long>& 
		vSubVolumeIndex = volumeMap_[-1];

	long nSub = vSubVolumeIndex.size();
	if(nSub > 0) os	<< TAB	<< "volumes";
	for(	long i = 0; i < nSub; i++) {
		os	<< "\t"	
			<< this->volume( vSubVolumeIndex[i] ).volumeName();
	}
	os << endl;

	if(this->isActive()){
		os	<< TAB	<< this->wire() << endl;		
	}

	list<HepPoint3D> refs = this->references();
	for(	list<HepPoint3D>::iterator iRef = refs.begin();
		iRef != refs.end();
		iRef++){
		os	<< TAB	<< "reference"		<< '\t'	<< (*iRef)	<< endl;
	}	

	os << TAB	<< this->detectorType() << endl;

	os << top	<< "end station" << endl;

	this->dumpVolume(os, pTab);

}

void GdbStation::dumpXML(ostream& os, const char* pTab){

	string top("");

	if( pTab != NULL ) top = pTab ;

	string TAB = top + '\t';
	os	<< "<--"	
		<< top	<< "VERSION"	<< "\t"	<< this->gVersion()	
		<< "-->"
		<< endl;

	os	<< top	<< "<station>"	<< endl;
	os	<< TAB	<< "<name>"	<< this->volumeName() << "</name>"	<< endl;
	os	<< TAB	<< "<solid>"	<< this->solid()	<< "</solid>"	<< endl;
	os	<< TAB	<< "<position>"	<< this->center()	<< "</position>"<< endl;
	os	<< TAB	<< "<rotation>"	<< this->rotMatrix() << "</rotation>"	<< endl;
	os	<< TAB	<< "<material>"	<< this->material() << "</material>"	<< endl;

	if( this->volumeMap().size() == 0) {
		long index(-1);
		this->createVolumeMap(index);		
	}

// referene point
  list<HepPoint3D> refs = this->references();
  for(	list<HepPoint3D>::iterator iRef = refs.begin();
	  iRef != refs.end();
	  iRef++){
	  os	<< TAB	<< "<reference>"	<< (*iRef)	<< "</reference>"<< endl;
  }	

// active area
  if(this->isActive()){
	  os	<< TAB	<< "<wire>"	<< this->wire() << "</wire>"	<< endl;		
  }

// detector type
  os << TAB	<< "<type>"	<< this->detectorType() << "</type>"	<< endl;


// dump subvolume name
  vector<long>& vSubVolumeIndex = volumeMap_[-1];
  long nSub = vSubVolumeIndex.size();
  const char* pTAB = TAB.c_str();
  for(	long i = 0; i < nSub; i++) {
	this->volumeXML( os, vSubVolumeIndex[i] ,  pTAB);
  }

  os << top	<< "</station>" << endl;

}

void GdbStation::volumeXML(ostream& os, const long& index, const char* pTab){
  string top("");
  if( pTab != NULL ) top = pTab ;

  string TAB = top + '\t';

  GdbVolume& volume = this->volume(index);

// start

  os	<< top	<< "<volume>"	<< endl;
  os	<< TAB	<< "<name>"	<< volume.volumeName() << "</name>"	<< endl;
  os	<< TAB	<< "<solid>"	<< volume.solid()	<< "</solid>"	<< endl;
  os	<< TAB	<< "<position>"	<< volume.center()	<< "</position>"<< endl;
  os	<< TAB	<< "<rotation>"	<< volume.rotMatrix() << "</rotation>"	<< endl;
  os	<< TAB	<< "<material>"	<< volume.material() << "</material>"	<< endl;

// active
  if(this->isActive()){
	os	<< TAB	<< "<wire>"	<< volume.wire() << "</wire>"	<< endl;		
  }

// dump subvolume name
  vector<long>& vSubVolumeIndex = volumeMap_[index];
  long nSub = vSubVolumeIndex.size();
  const char* pTAB = TAB.c_str();
  for(	long i = 0; i < nSub; i++) {
	this->volumeXML( os, vSubVolumeIndex[i] ,  pTAB);
  }

// end tag
  os << top	<< "</volume>" << endl;
  
}

GdbStation::operator GdbPerStation (){

// create persitent station
	GdbPerStation pS( (*this), this->nReference() );
//	CsErrLog::Instance()->mes(elDebugging, 
//		"GdbPerStation(const GdbVolume&, const long&) has been finished.");

// set gversion
	pS.gVersion( this->gVersion().c_str() );

//	CsErrLog::Instance()->mes(elDebugging, 
//		"Set geometry version to persitent object. (DONE)");

// set detector type
	pS.detectorType( this->detectorType() );
//	CsErrLog::Instance()->mes(elDebugging, 
//		"Set detector type to persitent object. (DONE)");
	
	
// set reference point	
	long i = 0;
	list<HepPoint3D> ref = this->references();
	list<HepPoint3D>::iterator iRef;
	for(	iRef =  ref.begin();
		iRef != ref.end();
		iRef++ ){
		pS.addReference( GdbDBuilder::convert( (*iRef) ), i );
		i++;
	}

//	CsErrLog::Instance()->mes(elDebugging, 
//		"Set reference point to persitent object. (DONE)");


// create volume list
	long nV = this->nVolume();
	for( long i=0; i<nV; i++){
		pS.addVolume( this->volume(i) );
	}

//	ostrstream out;
//	out	<< "Reference Number:("	<< pS.referenceNumber() 
//		<< ")\tOriginal("		<< this->nReference()	<< ")\n"
//		<< "SolidType:"		<< pS.solid()	<< endl;
//
//	CsErrLog::Instance()->mes(elDebugging, out.str() );
//	
//	cerr	<< "Reference Number:("	<< pS.referenceNumber() 
//			<< ")\tOriginal("		<< this->nReference()	<< ")\n"
//			<< "SolidType:"		<< pS.solid()	<< endl;

	
	return pS;
}

void GdbStation::createVolumeMap( long& index ){

	if( !( index < this->nVolume() ) ){
		CsErrLog::Instance()->mes( elError,
			"The given index exceeded volume size.");
		return;
	}

	if( index < 0 ) {
		index = -1;

//		 << "Initialize volume map" << endl;
//
		long nV = this->nVolume();
		
		for(long i=-1; i< nV; i++){
			long nSub = 
				( i == -1 ? this->nSubVolume() : this->volume(i).nSubVolume() );		
			vector<long> newVector(nSub, 0);
			volumeMap_[i] = newVector;
		}

	}
	
	long nSub	= 
		(index < 0 ? this->nSubVolume() : this->volume(index).nSubVolume() );

	if( nSub != 0 ){

		vector<long>& vectSubV = volumeMap_[index];

		long point = index;

//		cerr	<< "Create volume Map at [" << index 
//				<< "] which consists " << nSub << " volumes" << endl;

		for( long i=0; i<nSub; i++){

			index++;
//			cerr	<< "Add index [" << index << "] "
//					<< " to [" << point << "]"
//					<< " as [" << i << "] contents."
//					<< endl;
			vectSubV[i] = index ;

			
			this->createVolumeMap( index );
					
		}

	}

}

bool GdbStation::findActivePlane( 
	long& index, list<ActivePlane>& foundPlane ){


//	ostrstream debugBuff;
//	debugBuff	<< "Try to look if ["
//				<< index << "] plane is active or not, in total ["
//				<< this->nVolume() << "] volumes." << endl;
//	CsErrLog::Instance()->mes(elDebugging, debugBuff.str());

//	cerr	<< "Try to look if ["
//			<< index << "] plane is active or not, in total ["
//			<< this->nVolume() << "] volumes." << endl;

	bool findStatus(false);

//	if( this->refVolumes()[index].isActive() ){
	if( this->volume(index).isActive() ){
	
		ActivePlane newPlane;
		newPlane.index		= index;
		newPlane.center		= GdbPoint3D();
		newPlane.rotation	= GdbRotation();

		foundPlane.push_back( newPlane );

		findStatus = true;

	}
	
	
	long nSub = this->volume(index).nSubVolume();
	index++;
	
	if(nSub == 0) return (findStatus);

	
	for( long i = 0; i< nSub; i++ )	{
			
		list<ActivePlane> newList;
		
//		GdbPoint3D	center		= this->refVolumes()[index].center();
//		GdbRotation rotation	= this->refVolumes()[index].rotMatrix();
		GdbPoint3D	center		= this->volume(index).center();
		GdbRotation rotation	= this->volume(index).rotMatrix();

		if( this->findActivePlane( index, newList ) ){

			list<ActivePlane>::iterator iA;
			for(	iA = newList.begin();
				iA != newList.end();
				iA++){
				
				(*iA).center	= center + rotation * (*iA).center ;
				(*iA).rotation	= rotation * (*iA).rotation;

			}
			
			foundPlane.insert(
				foundPlane.end(), newList.begin(),  newList.end());
			
			findStatus = true;

		}

	}
	
	
	return (findStatus);
}

bool GdbStation::findActivePlane( list<ActivePlane>& foundPlane ){

	bool findStatus(false);

	if( this->isActive() ){
	
		ActivePlane newPlane;
		newPlane.index		= -1;
		newPlane.center		= GdbPoint3D();
		newPlane.rotation	= GdbRotation();

		foundPlane.push_back( newPlane );

		findStatus = true;

	}

	
	long startPoint(0);

	long nSub = this->nSubVolume();
	for( long i = 0; i < nSub ; i++ ){
	
		list<ActivePlane> newList;
		if( this->findActivePlane( startPoint, newList ) ){

			GdbVolume& volume = this->volume(0);

			list<ActivePlane>::iterator iA;
			for(	iA = newList.begin();
				iA != newList.end();
				iA++){

				(*iA).center = 
					volume.center() +
					volume.rotMatrix() * (*iA).center ;

				(*iA).rotation =
					volume.rotMatrix() * (*iA).rotation;

			}

			foundPlane.insert(
				foundPlane.end(), newList.begin(),  newList.end());

			findStatus = true;
		}

	}
//
//
//
	return (findStatus);
}
