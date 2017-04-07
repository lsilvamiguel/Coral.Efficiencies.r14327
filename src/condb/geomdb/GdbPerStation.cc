/*!
  \file		GdbPerStation.cc
  \brief	Geometry DB persistent station object implementation file
  \author	$Author: miyachi $
  \version	$Revision: 1.2 $
  \date		$Date: 2000/06/07 08:36:46 $
*/
#if __GNUG__ >= 2
#  pragma implementation
#endif

#include "CsSTD.h"
#include "CsErrLog.h"
#include "GdbPerStation.h"

//---- GdbPerStation ---------------------------------------------------------

GdbPerStation::GdbPerStation(	const char* name,
								const GdbDetectorType& detectorType,
								const GdbSolidType& solid,
								const GdbMaterial& material,
								const long& nReferences ) : 
	GdbVolume(name, solid, material, GdbPoint3D(0.0, 0.0, 0.0)),
	nRef_(nReferences),
	detectorType_(detectorType),
	references_(nReferences),
	volumes_() {
	
	this->gVersion("----------------");

}

GdbPerStation::GdbPerStation(	const GdbVolume& volume,
								const long& nReferences  ) :
	GdbVolume(	volume.volumeName(),
				volume.solid(),
				volume.material(),
				volume.center(),
				volume.rotMatrix(),
				volume.nSubVolume(),
				volume.isActive(),
				volume.wire() ),
	nRef_(nReferences),	
	detectorType_(),
	references_(nReferences),
	volumes_() {
//
//	CsErrLog::Instance()->mes(elDebugging, "GdbPerStation(const GdbVolume&, const long&)");
//
	this->gVersion("----------------");

}

GdbPerStation::GdbPerStation( const GdbPerStation& station ) :
	GdbVolume(	station.volumeName(),
				station.solid(),
				station.material(),
				station.center(),
				station.rotMatrix(),
				station.nSubVolume(),
				station.isActive(),
				station.wire() ),
	nRef_(station.referenceNumber()),	
	detectorType_(station.detectorType()), 
	references_(station.references()),
	volumes_( station.volumes() ) {

	this->gVersion(station.gVersion());

}

GdbPerStation::~GdbPerStation(){
}

const char* GdbPerStation::gVersion() const {
	return geomVersion_;
}

const char* GdbPerStation::gVersion( const char* pgVersion ){
	return strncpy(geomVersion_, pgVersion, GdbPerStation::MAXCHAR );
}

void GdbPerStation::detectorType(const GdbDetectorType& type){
	detectorType_ = type;
}

GdbDetectorType GdbPerStation::detectorType() const {
	return detectorType_;
}

long GdbPerStation::referenceNumber() const {
	return nRef_;
}

bool GdbPerStation::addReference(const GdbPoint3D& point, const long& index){
	if(this->checkIndex(index)){
		this->references()[index] = point;
	}
	return false;
}

GdbPoint3D GdbPerStation::referencePoint(const long& index) {
	if(this->checkIndex(index)) {
		return this->references()[index];
	}
	CsErrLog::Instance()->mes(elError, "return NULL vector");
	return GdbPoint3D(0.0, 0.0, 0.0);
}


d_Varray<GdbPoint3D>& GdbPerStation::refReferences() {
	return references_;
}

d_Varray<GdbPoint3D> GdbPerStation::references() const {
	return references_;
}

bool GdbPerStation::checkIndex( const long& index ) const {
	if( index < this->referenceNumber() ) {
		return true;
	}

	ostrstream out;
	out << "the given index "
		<< index
		<< " exceedes the maximun length "
		<< this->referenceNumber()
		<< endl;
	
	CsErrLog::Instance()->mes(elError, out.str());
	
	return false;
}

d_Varray<GdbVolume>& GdbPerStation::refVolumes() {
	return volumes_;
}

d_Varray<GdbVolume> GdbPerStation::volumes() const {
	return volumes_;
}

void GdbPerStation::volumes(const d_Varray<GdbVolume> volumes) {
	volumes_ = volumes;
}

void GdbPerStation::addVolume( const GdbVolume& volume ){
	volumes_.extend(volume);
}

long GdbPerStation::nVolume() const {
	return volumes_.size();
}

GdbVolume& GdbPerStation::volume(const long& i ) const {
	return volumes_[i];
}


ostream& operator<<(ostream& os, GdbPerStation& station){

	station.dump(os);
	
	return os;
}

void GdbPerStation::dumpVolume(ostream& os, const char* pTab){
	for( long i = 0; i < this->nVolume(); i++){
		this->dumpVolume(os, i, pTab);
	}
}

void GdbPerStation::dumpVolume(ostream& os, const long& index, const char* pTab){
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


	if( volume.nSubVolume() > 0)
		os << TAB	<< "nvolume"	<< '\t' << volume.nSubVolume() << endl;

	if(this->isActive()){
		os	<< TAB	<< this->wire() << endl;		
	}

	os << top	<< "end volume" << endl;

}


void GdbPerStation::dump(ostream& os, const char* pTab){
	string top("");
	if( pTab != NULL ) top = pTab ;

	string TAB = top + '\t';

	os	<< "#"	
		<< top	<< "VERSION"	<< "\t"	<< this->gVersion()	<< endl;
	os	<< top	<< "station"	<< endl;
	os	<< TAB	<< "name"		<< "\t"	<< this->volumeName() << endl;
	os	<< TAB	<< this->solid() << endl;
	os	<< TAB	<< "position"	<< '\t'	<< this->center()	<< endl;
	os	<< TAB	<< this->rotMatrix() << endl;
	os	<< TAB	<< this->material() << endl;

	if(this->isActive()){
		os	<< TAB	<< this->wire() << endl;		
	}

	long nRef = this->referenceNumber();
	for( long i=0; i< nRef ; i++){
		GdbPoint3D ref = this->refReferences()[i];
		os	<< TAB	<< "reference"		<< '\t'	<< ref <<	endl;
	}

	if(this->nSubVolume() > 0)
		os << TAB	<< "nvolume"	<< '\t' << this->nSubVolume() << endl;

	os << TAB	<< this->detectorType() << endl;

	os << top	<< "end station" << endl;

	this->dumpVolume(os, pTab);

}
