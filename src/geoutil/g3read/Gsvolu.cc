#if __GNUG__ >= 2
#  pragma implementation
#endif
/*!
  \file		Gsvolu.cc
  \brief	volume entries in Geant 3. (implementation file)
  \author	$Author: benigno $
  \version	$Revisopn$
  \date		$Date: 2003/04/24 13:26:03 $
*/

#include "Gsvolu.h"

using namespace std;

//---- Gsvolu ---------------------------------------------------------
Gsvolu::Gsvolu(	
		const string&        name,
		const string&        type,
		const unsigned int& materialId,
		const int&         nparam,
		const vector<double>& param ):
	name_(name), type_(type), materialId_(materialId),
	nparam_(nparam), param_(param), pMedium_(NULL) {
}

Gsvolu::Gsvolu(const Gsvolu& volume):
	name_(volume.name()), type_(volume.type()), 
	materialId_(volume.materialId()),
	nparam_(volume.nparam()), param_(volume.param()),
	pMedium_(volume.pMedium()) {
}

Gsvolu::~Gsvolu()
{
}

ostream& operator<<(ostream& os, Gsvolu& volume ){

	os	<< "NAME:\t["	<< volume.name()	<< "]\t"
		<< "MID:\t"	<< volume.materialId() 
		<< " ("		<< volume.pMedium()	<< ")\t"
		<< "TYPE:\t["	<< volume.type()	<< "]\t"
		<< "N:\t"		<< volume.nparam()	<< "\t";
	
	vector<double> param = volume.param();
	for(int i=0; i< volume.nparam(); i++){
		os	<< param[i]	<< ", ";
	}
	return os;
}
