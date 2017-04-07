#if __GNUG__ >= 2
#  pragma implementation
#endif

/*!
  \file		Gspos.cc
  \brief	GSPOS entry in G3 call object implementation file.
  \author	$Author: tnagel $
  \version	$Revisopn$
  \date		$Date: 2010/04/27 19:13:13 $
*/

#include <algorithm>

#include "Gspos.h"
#include "CsErrLog.h"

using namespace std;
using namespace CLHEP;

//---- Gspos ---------------------------------------------------------
Gspos::Gspos() :
  id_(0), name_(""), unit_(0), mother_(""), center_(), 
  rotation_(), flag_(),  n_(1), motherN_(0),
  pMother_(NULL), pRotation_(NULL), pDetector_(NULL),
  param_(), pVolume_(NULL), pDaughter_(),
  hasDetectorInTree_(false) 
{
  rM_ = HepMatrix(3,3,1);
}

Gspos::Gspos(const string& name, const unsigned int&  unit) :
  id_(0),
  name_(name), unit_(unit), mother_(""), center_(), 
  rotation_(), flag_(),  n_(1), motherN_(0),
  pMother_(NULL), pRotation_(NULL), pDetector_(NULL),
  param_(), pVolume_(NULL), pDaughter_(),
  hasDetectorInTree_(false) {
  rM_ = HepMatrix(3,3,1);
}


Gspos::Gspos(	
  const string&        name,
  const unsigned int&  unit,
  const string&      mother,
  const HepGeom::Vector3D<double>  center,
  const int&         rotation,
  const string&        flag ) :
  id_(0),
  name_(name), unit_(unit), mother_(mother), center_(center), 
  rotation_(rotation), flag_(flag),  n_(1), motherN_(0),
  pMother_(NULL), pRotation_(NULL), pDetector_(NULL),
  param_(), pVolume_(NULL), pDaughter_(),
  hasDetectorInTree_(false) {
  rM_ = HepMatrix(3,3,1);
}

Gspos::Gspos(const Gspos& pos):
  id_(0),
  name_(pos.name()), unit_(pos.unit()),
  mother_(pos.mother()), center_(pos.center()),
  rotation_(pos.rotation()), flag_(pos.flag()), 
  n_(pos.n()), motherN_(pos.motherN()),
  pMother_(pos.pMother()), pRotation_(pos.pRotation()),
  rM_(pos.rM()),
  pDetector_(pos.pDetector()),
  param_(pos.param()), pVolume_(pos.pVolume()), 
  pDaughter_(pos.pDaughter()), 
  hasDetectorInTree_(pos.hasDetectorInTree()) {
}


Gspos::~Gspos(){
}

void Gspos::addToDaughter(Gspos* pPos){
  pDaughter_.push_back( pPos );
}

ostream& operator<<(ostream& stream, Gspos& pos){
  stream
	<<	"ID:"		<<	pos.id()	<<	'\t'
	<<	"NAME:"	<<	pos.name()	<<	"\t"
	<<	"UNIT:"	<<	pos.unit()	<<	"\t" 
	<<	"Number:"	<<	pos.n()		<<	"\n"
	<<	"Mother:"	<< pos.mother()	<<	"\t"
	<<	"Number:" << pos.motherN()	<<	'\t'
	<<	"ADDRESS:"	<< pos.pMother()	<<	'\t'
	<<	"Detector:"	<< pos.pDetector() << "\n\t";
  pos.dumpMother(stream);
  return stream;
}


void Gspos::dumpDaughters(ostream& os, const string tab){

  list<Gspos*>::iterator iPos;

  os	
	<< tab << "+" << this->name() << "\t" 
	<< this->id() << ":" << this->unit() << "\t";

  if( this->pDetector() != NULL ) os << "(Detector)\t";
  os << endl;

  for(	
	iPos = this->refPDaughter().begin();
	iPos != this->refPDaughter().end();
	iPos++){
	(*iPos)->dumpDaughters( os, tab + "|" );
  }
	
}

bool Gspos::checkDetector(){
  this->hasDetectorInTree(false);

  if( this->pDetector() != NULL ) {
	this->hasDetectorInTree(true);
  }

  list<Gspos*>::iterator iPos;
  for(	iPos = this->refPDaughter().begin();
	iPos != this->refPDaughter().end();
	iPos++){
	if( (*iPos)->checkDetector() ) {
	  this->hasDetectorInTree(true);
	}
  }
  return this->hasDetectorInTree();
}

void Gspos::dumpDetectors(ostream& os, const string tab){
  list<Gspos*>::iterator iPos;

  if( this->hasDetectorInTree() ){

	os	<< tab << "+" << this->name() << "\t" 
		<< this->id() << ":" << this->unit() << "\t";
	if( this->pDetector() != NULL ) os << "(Detector)\t";
	os << endl;

	for(	iPos = this->refPDaughter().begin();
		iPos != this->refPDaughter().end();
		iPos++){
		(*iPos)->dumpDetectors( os, tab + "|" );
	}
  }

}

void Gspos::dumpMother( ostream& stream ) const {
  if( this->pMother() !=NULL ){

	stream	<< "<--" 
			<< this->pMother()->id() << "."
			<< this->pMother()->name();

	this->pMother()->dumpMother(stream);

  }	
}

Gspos& Gspos::findPosition(Gspos& pos, bool& result){
  if( ( result = ( (*this) == ( pos ) ) ) ) return *this;

  list<Gspos*> pDaughter = this->pDaughter();
  list<Gspos*>::iterator iPos;

  for(	iPos = pDaughter.begin() ;
	iPos != pDaughter.end(); iPos++){
    if( ( result = ( *(*iPos) == pos) ) ) return *(*iPos);
  }

  for(	iPos = pDaughter.begin() ;
	iPos != pDaughter.end(); iPos++ ){
	Gspos& found = (*iPos)->findPosition( pos , result );
	if( result == true ) return found ;
  }

  result = false;	
  return *this;
}

Gspos& Gspos::findPosition(	const string& name, 
							const unsigned int& unit, 
							bool& result ){	
  result = false;
  if( name == "" ) return *this;
  Gspos target( name, unit );

  return this->findPosition( target, result );
}
