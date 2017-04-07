/*!
  \file		GdbSolidType.cc
  \brief	Geometry DB solid object implementation file
  \author	$Author: miyachi $
  \version	$Revision: 1.3 $
  \date		$Date: 2000/06/21 00:10:52 $
*/
#if __GNUG__ >= 2
#  pragma implementation
#endif

#include "CsErrLog.h"
#include "CsSTD.h"
#include "GdbSolidType.h"
#include "CsCOMGEANT.h"

//---- GdbSolidType ---------------------------------------------------------
GdbSolidType::GdbSolidType() : 
  nParam_(0), 
  parameters_() {
  this->type("NONE");
  this->initParam();
}

GdbSolidType::GdbSolidType(const char* type, const long& nParam) : 
  nParam_(0), 
  parameters_(){

  int n = COMGEANT::nParam( type );
  if( n != int(nParam) ){
	ostrstream out;
	out	<< "The given number of parameters [" << nParam 
		<< "] is not correct. This must be [" << n
		<< "]." ;
	CsErrLog::Instance()->mes( elFatal, out.str());
  }

  this->initParam();	
  this->nParam(nParam);
  this->type(type);
}

GdbSolidType::GdbSolidType(const GdbSolidType& solid) :
  nParam_(solid.nParam()), 
  parameters_(){
  this->type(solid.type());
  this->initParam();	
  for(long i=0; i<solid.nParam(); i++) this->addParam(i, solid.param(i) );
}

GdbSolidType::~GdbSolidType(){
}

void GdbSolidType::initParam(){
  for(long i=0; i<GdbSolidType::maxParameterSize; i++)
	  parameters_[i];
//		this->addParam( i, 0.0);
}

GdbSolidType& GdbSolidType::operator=(const GdbSolidType& solid){
  if(this != &solid){
	  this->type( solid.type());
	  this->nParam( solid.nParam());
	  this->initParam();	
	  for(long i=0; i<solid.nParam(); i++)
			  this->addParam(i, solid.param(i) );
  }
  return *this;
}

void GdbSolidType::addParam(const long& index, const double& param){
  if(this->checkIndex(index)){
	  parameters_[index] = param;
  }
}

double	GdbSolidType::param(const long& index) const {
  if(this->checkIndex(index)){
	  return parameters_[index];
  }
  return -99999999999.9;
}


const char* GdbSolidType::type() const {
  return type_;
}

const char* GdbSolidType::type(const char* type){
  return strncpy( type_, type, GdbSolidType::maxNameLength);
}

void GdbSolidType::nParam(const long nParam){
  if(nParam > GdbSolidType::maxParameterSize) {
#ifdef CsErrLog_h
	  CsErrLog::Instance()->mes(elError,
		  "parameter size must be smaller than maxParameterSize.");
#else
	  cerr	<< "GdbSolidType::nParam\tparameter size must be smaller than maxParameterSize." 
			  << endl;
#endif
	  nParam_ = maxParameterSize;
	  return;
  }
  nParam_ = nParam;
}


long GdbSolidType::nParam() const {
  return nParam_;
}

const double* GdbSolidType::parameters() const {
  return parameters_;
}

bool GdbSolidType::checkIndex(const long& index) const {
  if( index < nParam_ ) {
	  return true;
  }
#ifdef CsErrLog_h
  string mes = "the given index[";
  mes += index;
  mes += "] exceedes the maximun length [";
  mes += nParam_;
  mes += "]";
  CsErrLog::Instance()->mes(elError, mes);
#else
  cerr	<< "the given index[" << index 
		  << "] exceedes the maximum length [" << nParam_ << "]"
		  << endl;
#endif
  return false;
}

ostream& operator<<(ostream& os, const GdbSolidType& solid){
  string type = solid.type();
  long nParam = solid.nParam();

  os	<< "shape\t"	<< type << "\t" << nParam << "\t";
  for(long i = 0; i<nParam; i++) {
	  os << solid.param(i) << "\t";
  }

  return os;
}
