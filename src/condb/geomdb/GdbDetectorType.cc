/*!
  \file		GdbDetectorType.cc
  \brief	Geometry DB detector type object implementation file
  \author	$Author: miyachi $
  \version	$Revision: 1.2 $
  \date		$Date: 2000/06/07 08:36:45 $
*/
#if __GNUG__ >= 2
#  pragma implementation
#endif

#include "CsSTD.h"
#include "GdbDetectorType.h"

GdbDetectorType::GdbDetectorType() :
	detectorType_(0), 
	gateLength_(0.0),
	eff_(1.0),
	bkg_(0.0),
	hasDrift_(false),
	driftVelocity_(0.0),
	driftT0_(0.0),
	doubleHitResolution_(0.0),
	spaceResolution_(0.0),
	timeSlice_(0.0)
	{
	this->name("------");
}

GdbDetectorType::~GdbDetectorType(){
}

GdbDetectorType::GdbDetectorType(	const char* name,
									const int& type,
									const double& gate ) :
	detectorType_(type), 
	gateLength_(gate),
	eff_(1.0),
	bkg_(0.0),
	hasDrift_(false),
	driftVelocity_(0.0),
	driftT0_(0.0),
	doubleHitResolution_(0.0),
	spaceResolution_(0.0),
	timeSlice_(0.0)
	{
	this->name( name );
}

GdbDetectorType::GdbDetectorType(
		const char* name, 
		const int& type, 
		const double& gate,
		const double& eff,
		const double& bkg
	) :
	detectorType_(type),
	gateLength_(gate),
	eff_(eff),
	bkg_(bkg),
	hasDrift_(false),
	driftVelocity_(0.0),
	driftT0_(0.0),
	doubleHitResolution_(0.0),
	spaceResolution_(0.0),
	timeSlice_(0.0)
	{
	this->name( name );
}

GdbDetectorType::GdbDetectorType(const GdbDetectorType& dType) :
	detectorType_(dType.type()), 
	gateLength_(dType.gate()),
	eff_(dType.eff()),
	bkg_(dType.bkg()),
	hasDrift_(dType.hasDrift()),
	driftVelocity_(dType.driftV()),
	driftT0_(dType.T0()),
	doubleHitResolution_(dType.dHitRes()),
	spaceResolution_(dType.sRes()),
	timeSlice_(dType.tSlic())
	{
	this->name(dType.name());
}


GdbDetectorType& GdbDetectorType::operator=(const GdbDetectorType& dType){
	if(this != &dType){
		detectorType_	= dType.type();
		gateLength_		= dType.gate();
		eff_			= dType.eff();
		bkg_			= dType.bkg();
		this->name(		dType.name());
		this->hasDrift(	dType.hasDrift());
		this->driftV(	dType.driftV());
		this->T0(		dType.T0());
		this->dHitRes(	dType.dHitRes());
		this->sRes(		dType.sRes());
		this->tSlic(	dType.tSlic());
	}
	return *this;
}

const char* GdbDetectorType::name(const char* name){
	return  strncpy( detectorName_, name, GdbDetectorType::maxNameLength );		
}


const char* GdbDetectorType::name() const {
	return detectorName_;
}

int GdbDetectorType::type() const {
	return detectorType_;
}

double GdbDetectorType::gate() const{
	return gateLength_;
}

double GdbDetectorType::eff() const{
	return eff_;
}

double GdbDetectorType::bkg() const{
	return bkg_;
}

bool GdbDetectorType::hasDrift() const {return hasDrift_;}

double GdbDetectorType::driftV() const {return driftVelocity_;}
double GdbDetectorType::T0() const {return driftT0_;}

double GdbDetectorType::dHitRes() const {return doubleHitResolution_;}
double GdbDetectorType::sRes() const {return spaceResolution_;}
double GdbDetectorType::tSlic() const {return timeSlice_;}

void GdbDetectorType::hasDrift(const bool& hasDrift) {
	hasDrift_ = hasDrift;
}

void GdbDetectorType::driftV(const double& driftV) { 
	driftVelocity_ = driftV;
} 

void GdbDetectorType::T0(const double& T0) {driftT0_ = T0;} 

void GdbDetectorType::dHitRes(const double& dHitRes) {
	doubleHitResolution_ = dHitRes;
} 

void GdbDetectorType::sRes(const double& sRes) {
	spaceResolution_ = sRes;
} 

void GdbDetectorType::tSlic(const double& tSlic) {timeSlice_ = tSlic;} 

void GdbDetectorType::setDriftInfo(	
	const double& driftV,
	const double& T0,
	const double& dHitRes,
	const double& sRes,
	const double& tSlic) {

	this->hasDrift(true);
	this->driftV(driftV);
	this->T0(T0);
	this->dHitRes(dHitRes);
	this->sRes(sRes);
	this->tSlic(tSlic);

}


ostream& operator<<(ostream& os, const GdbDetectorType& det){
	string name = det.name();
	os	<< "detector"		<< '\t'
		<< name			<< '\t'
		<< det.type()	<< '\t'
		<< det.gate()	<< '\t'
		<< det.eff()	<< '\t'
		<< det.bkg()	;
	return os;
}

