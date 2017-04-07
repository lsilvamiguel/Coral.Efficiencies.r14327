#if __GNUG__ >= 2
#  pragma implementation
#endif
/*!
  \file		Gsdet.cc
  \brief	Object for GSDET entries in G3 call (implementation file)
  \author	$Author: benigno $
  \version	$Revisopn$
  \date		$Date: 2003/04/24 13:26:02 $
*/

#include "Gsdet.h"

using namespace std; 

//---- Gsdet ---------------------------------------------------------

Gsdet::Gsdet() :
  name_(), mother_(), type_(), 
  hitName_(), hitBit_(), hitFactor_(), hitOffset_(), 
  digitName_(), digitBit_(), user_() {
}

Gsdet::Gsdet( const string& name, const string& mother, const int& type ) : 
  name_(name), mother_(mother), 
  type_(type), pMaterial_(NULL),
  hitName_(), hitBit_(), hitFactor_(), hitOffset_(), 
  digitName_(), digitBit_(), user_() {
}

Gsdet::Gsdet( const Gsdet& gsdet ):
  name_(gsdet.name()), mother_(gsdet.mother()), 
  type_(gsdet.type()), pMaterial_(gsdet.pMaterial()),
  hitName_(gsdet.hitName()), hitBit_(gsdet.hitBit()), 
  hitFactor_(gsdet.hitFactor()), hitOffset_(gsdet.hitOffset()), 
  digitName_(gsdet.digitName()), digitBit_(gsdet.digitBit()), 
  user_(gsdet.user()) {
}

Gsdet::~Gsdet()
{
}

Gsdet& Gsdet::operator=( const Gsdet& gsdet ){
  if( this != &gsdet ){
	this->name(			gsdet.name()		);
	this->mother(		gsdet.mother()		);
	this->type(			gsdet.type()		);
	this->pMaterial(	gsdet.pMaterial()	);
	this->hitName(		gsdet.hitName()		);
	this->hitBit(		gsdet.hitBit()		);
	this->hitFactor(	gsdet.hitFactor()	);
	this->hitOffset(	gsdet.hitOffset()	);
	this->digitName(	gsdet.digitName()	);
	this->digitBit(		gsdet.digitBit()	);
	this->user(			gsdet.user()		);
  }
  return *this;
}
