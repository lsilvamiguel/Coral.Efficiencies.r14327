/*!
  \file		GdbMaterial.cc
  \brief	Geometry DB material description object implementation file
  \author	$Author: miyachi $
  \version	$Revision: 1.2 $
  \date		$Date: 2000/06/07 08:36:45 $
*/
#if __GNUG__ >= 2
#  pragma implementation
#endif

#include "CsSTD.h"
#include "CsErrLog.h"
#include "GdbMaterial.h"

//---- GdbMaterial ---------------------------------------------------------
GdbMaterial::GdbMaterial() :
	density_(0.0), radLen_(0.0) {
	this->materialName("Non");
}



GdbMaterial::GdbMaterial(	const char* name,
							const double density,
							const double radiationLength) : 
	density_(density), radLen_(radiationLength) {
	this->materialName(name);
}

GdbMaterial::GdbMaterial(const GdbMaterial& material) : 
	density_(material.density()), radLen_(material.radLen()) {
	this->materialName( material.materialName());
}

GdbMaterial::~GdbMaterial(){
}

GdbMaterial& GdbMaterial::operator=(const GdbMaterial& material){
	if(this != &material){
		density_ =	material.density();
		radLen_ =	material.radLen();
		this->materialName(material.materialName());
	}
	return *this;
}



const char* GdbMaterial::materialName() const {
	return materialName_;
}


double GdbMaterial::radLen() const {
	return radLen_;
}

double GdbMaterial::density() const {
	return density_;
}

void GdbMaterial::materialName(const char* name){
	strncpy( materialName_, name, GdbMaterial::maxNameLength);
}

void GdbMaterial::radLen(const double& radLen){
	radLen_ = radLen;
}

void GdbMaterial::density(const double& dens){
	density_ = dens;
}


istream& operator>>(istream& is, GdbMaterial& material){
	string name;
	double  dens(0.0), radLen(0.0);
	if(is.good()){
		is >> name;
		if(is.good()){
			is >> dens;
			if(is.good()){
				is >> radLen;
			} else {
				CsErrLog::Instance()->mes(elInfo, 
					"Material line must have [name] [dens] [radiation length(missing)].");
			}
		} else {
			CsErrLog::Instance()->mes(elInfo, 
			"Material line must have [name] [dens(missing)] [radiation length(missing)].");
		}
		
	} else {
		CsErrLog::Instance()->mes(elInfo, 
		"Material line must have [name(missing)] [dens(missing)] [radiation length(missing)].");		
	}
	material.materialName( name.c_str() );
	material.radLen(radLen);
	material.density(dens);
	return is;
}

ostream& operator<<(ostream& os, const GdbMaterial& material){
	string name = material.materialName();
	os	<< "material" << "\t"	
		<< name	<< "\t"
		<< material.density()		<< "\t"
		<< material.radLen()		;
	return os;
}
