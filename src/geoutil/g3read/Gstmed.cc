#if __GNUG__ >= 2
#  pragma implementation
#endif
/*!
  \file		Gstmed.cc
  \brief	Object for Material entries in Geant 3. (implementation file)
  \author	$Author: benigno $
  \version	$Revisopn$
  \date		$Date: 2003/04/24 13:26:03 $
*/

#include "Gstmed.h"

using namespace std;

//---- Gstmed ---------------------------------------------------------

Gstmed::Gstmed():
	id_(), name_("NoName"), materialId_(0), isvol_(0), ifield_(0), pMaterial_(NULL) {
}

Gstmed::Gstmed(	const unsigned int& id,
				const string& name,
				const unsigned int& materialId,
				const int& isvol,
				const int& ifield) :
	id_(id), name_(name), materialId_(materialId), isvol_(isvol), 
	ifield_(ifield), pMaterial_(NULL) {
}

Gstmed::Gstmed(const Gstmed& gstmed) :
	id_(gstmed.id()), name_(gstmed.name()), materialId_(gstmed.materialId()), 
	isvol_(gstmed.isvol()), ifield_(gstmed.ifield()),
	 pMaterial_(gstmed.pMaterial()) {
}


Gstmed::~Gstmed()
{
}

Gstmed& Gstmed::operator=(const Gstmed& gstmed){
	if(this != &gstmed){
		this->id(gstmed.id());
		this->name(gstmed.name());
		this->materialId(gstmed.materialId());
		this->isvol(gstmed.isvol());
		this->ifield(gstmed.ifield());
		this->pMaterial(gstmed.pMaterial());
	}
	return *this;
}

ostream& operator<<( ostream& os, const Gstmed& med){

	os	<< "ID:\t"	<< med.id()		<< "\t"
		<< "NAME:\t["	<< med.name()	<< "]\t"
		<< "VOL:\t"	<< med.isvol()	<< "\t"
		<< "FIELD:\t"	<< med.ifield()
		<< "Material:\t"	<< med.materialId() << "\t"
		<< "Material*:\t"	<< med.pMaterial();

	return os;
}
