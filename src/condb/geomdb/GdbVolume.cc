/*!
  \file		GdbVolume.cc
  \brief	Geometry DB volume object implementation file
  \author	$Author: miyachi $
  \version	$Revision: 1.2 $
  \date		$Date: 2000/06/07 08:36:48 $
*/
#if __GNUG__ >= 2
#  pragma implementation
#endif

#include "CsSTD.h"
#include "GdbVolume.h"

#  include <ospace/std/algorithm>
#  include <ospace/std/functional>



//---- GdbVolume ---------------------------------------------------------
GdbVolume::GdbVolume() :
	solid_(),
	material_(), 
	center_(), 
	rotMatrix_(), 
//	subVolume_(),
	nSubVolume_(0),
	isActive_(false), 
	wire_() 	
	{
	this->volumeName("--------");
}

GdbVolume::GdbVolume(	const GdbPoint3D& center,
						const GdbRotation& rotation	) :
	solid_(), 
	material_(), 
	center_(center), 
	rotMatrix_(rotation), 
//	subVolume_(),
	nSubVolume_(0),
	isActive_(false), 
	wire_() 
	{
	this->volumeName("--------");
}


GdbVolume::GdbVolume(
  const char* volumeName, const GdbSolidType solid, const GdbMaterial& material,
  const GdbPoint3D& center) :
	solid_(solid), 
	material_(material), 
	center_(center), 
	rotMatrix_(), 
	nSubVolume_(0),
	isActive_(false), 
	wire_() 
	{
	this->volumeName(volumeName);
}

GdbVolume::GdbVolume(	const char* volumeName, 
						const GdbSolidType solid,
						const GdbMaterial material,
						const GdbPoint3D& center,
						const GdbRotation& rotation,
//						const GdbVolume::VolumeContainer& subVolume,
						const long& nSubVolume,
						const bool& active,
						const GdbWireInfo& wire) : 
	solid_(solid), 
	material_(material), 
	center_(center), 
	rotMatrix_(rotation), 
//	subVolume_(subVolume),
	nSubVolume_(nSubVolume),
	isActive_(active), 
	wire_(wire)  {
	this->volumeName(volumeName);
}

GdbVolume::GdbVolume(const GdbVolume& gVolume) :
	solid_(gVolume.solid()),
	material_(gVolume.material()), 
	center_(gVolume.center()),
	rotMatrix_(gVolume.rotMatrix()), 
//	subVolume_(gVolume.subVolume()),
	nSubVolume_(gVolume.nSubVolume()),	
	isActive_(gVolume.isActive()), 
	wire_(gVolume.wire())   {
	this->volumeName(gVolume.volumeName());
}

GdbVolume::~GdbVolume(){
}

GdbVolume& GdbVolume::operator=(const GdbVolume& gVolume){
	if(this != &gVolume){
		this->volumeName(gVolume.volumeName());
		this->solid(gVolume.solid());
		this->material( gVolume.material()); 
		this->center( gVolume.center() );
		this->rotMatrix( gVolume.rotMatrix() ); 
//		this->subVolume( gVolume.subVolume() );
		this->nSubVolume(gVolume.nSubVolume()),	
		this->isActive( gVolume.isActive() ); 
		this->wire( gVolume.wire() );
	}
	return *this;
}

const char* GdbVolume::volumeName() const {
	return volumeName_;
}

GdbSolidType GdbVolume::solid() const {
	return solid_;
}

GdbMaterial GdbVolume::material() const {
	return material_;
}

GdbPoint3D GdbVolume::center() const {
	return center_;
}

GdbRotation GdbVolume::rotMatrix() const {
	return rotMatrix_;
}
	
GdbWireInfo GdbVolume::wire() const {
	return wire_;
}

bool GdbVolume::isActive() const {
	return isActive_;
}

void GdbVolume::isActive(const bool& active){
	isActive_ = active;
}

void GdbVolume::volumeName(const char* volumeName){
	strncpy( volumeName_, volumeName, GdbVolume::maxNameLength );		
}

void GdbVolume::solid(const GdbSolidType& solid){
	solid_ = solid;
}

void GdbVolume::material(const GdbMaterial& material){
	material_ = material;
} 
void GdbVolume::center(const GdbPoint3D& center){
	center_ = center;
}
void GdbVolume::rotMatrix(const GdbRotation& rotMartix){
	rotMatrix_ = rotMartix;
}


bool GdbVolume::operator<(const GdbVolume& gVolume) const {
	return true;
}

bool GdbVolume::operator==(const GdbVolume& gVolume) const {
	return (this == &gVolume);
}

bool GdbVolume::operator==(const char* name) const {
	if( strcmp( this->volumeName(), name ) == 0 ) return true;
	return false;
}

void GdbVolume::wire(const GdbWireInfo& wire){
	wire_ = wire;
	isActive_ = true;
}

ostream& operator<<(ostream& os, GdbVolume& volume){
	volume.dump(os);
	return os;
}

void GdbVolume::dump(ostream& os, const char* pTab){
	string top("");
	if( pTab != NULL ) top = pTab ;

	string TAB = top + '\t';

	os	<< top	<< "volume"	<< endl;
	os	<< TAB	<< "name"		<< "\t"	<< this->volumeName() << endl;
	os	<< TAB	<< this->solid() << endl;
	os	<< TAB	<< "position"	<< "\t"	<< this->center()	<< endl;
	os	<< TAB	<< this->rotMatrix() << endl;
	os	<< TAB	<< this->material() << endl;

	if(this->isActive()){
		os	<< TAB	<< this->wire() << endl;		
	}
	os	<< top	<< "end volume" << endl;

}


long GdbVolume::nSubVolume() const {
	return nSubVolume_;
}

void GdbVolume::nSubVolume(const long& nSub) {
	nSubVolume_ = nSub;
}
