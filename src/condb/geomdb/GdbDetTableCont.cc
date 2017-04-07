/*!
  \file		GdbDetTableCont.cc
  \brief	Geometry DB detector table object implementation file
  \author	$Author: miyachi $
  \version	$Revision: 1.2 $
  \date		$Date: 2000/06/07 08:36:44 $
*/
#if __GNUG__ >= 2
#  pragma implementation
#endif

#include "GdbDetTableCont.h"
#include "CsCOMGEANT.h"

//---- GdbDetTableCont ---------------------------------------------------------
GdbDetTableCont::GdbDetTableCont() : 
	id_(0), unit_(0), type_(0), radLen_(0.0),
	solidtype_(COMGEANT::BOX),
	param_(COMGEANT::nParam(COMGEANT::BOX), 0.0), 
	center_(0.0, 0.0, 0.0),  rotM_(3, 3, 1), 
	wireDistance_(0.0), angle_(0.0), 
	wireNumber_(0), wirePitch_(0.0), 
	efficiency_(0.0), backGround_(0.0), timeGate_(0.0), 
 	hasDrift_(false), 
	driftVelocity_(0.0), driftT0_(0.0), 
	doubleHitResolution_(0.0), 
	spaceResolution_(0.0), 
	timeSlice_(0.0),
	hasDeadspace_(false),
	dId_(0),
	dType_(COMGEANT::BOX),
	dParam_(COMGEANT::nParam(COMGEANT::BOX), 0.0),
	dCenter_(0.0, 0.0, 0.0),
	dRotM_(3,3,1)
	{
	this->name("------");
}
	bool hasDeadspace_;
	vector<double> dParam_;
	HepPoint3D dCenter_;
	HepMatrix  dRotM_;

GdbDetTableCont::GdbDetTableCont(
		const int& id,
		const char* name, 
		const int& unit,
		const int& type,
		const double& radLen,
		const HepPoint3D& size,
		const HepPoint3D& center, 
		const HepMatrix& rotation, 
		const double& wirD,
		const double& angle,
		const int& nWir,
		const double& wirP,
		const double& eff,
		const double& bgk,
		const double& tGate
	 ) :
	id_(id), unit_(unit), type_(type), radLen_(radLen),
	solidtype_(COMGEANT::BOX),
	param_(COMGEANT::nParam(COMGEANT::BOX), 0.0), 
	center_(center),  rotM_(rotation), 
	wireDistance_(wirD), angle_(angle), 
	wireNumber_(nWir), wirePitch_(wirP), 
	efficiency_(eff), backGround_(bgk), timeGate_(tGate), 
 	hasDrift_(false), 
	driftVelocity_(0.0), driftT0_(0.0), 
	doubleHitResolution_(0.0), 
	spaceResolution_(0.0), 
	timeSlice_(0.0),
	hasDeadspace_(false),
	dId_(0),
	dType_(COMGEANT::BOX),
	dParam_(COMGEANT::nParam(COMGEANT::BOX), 0.0),
	dCenter_(0.0, 0.0, 0.0),
	dRotM_(3,3,1)
	{

	this->refParam()[0] = size.x();
	this->refParam()[1] = size.y();
	this->refParam()[2] = size.z();
	
	this->name(name);
}

GdbDetTableCont::GdbDetTableCont(
		const int& id,
		const char* name, 
		const int& unit,
		const int& type,
		const double& radLen,
		const int& solidtype,
		const vector<double>& param,
		const HepPoint3D& center, 
		const HepMatrix& rotation, 
		const double& wirD,
		const double& angle,
		const int& nWir,
		const double& wirP,
		const double& eff,
		const double& bgk,
		const double& tGate
	) :
	id_(id), unit_(unit), type_(type), radLen_(radLen),
	solidtype_(solidtype),
	param_(param), 
	center_(center),  rotM_(rotation), 
	wireDistance_(wirD), angle_(angle), 
	wireNumber_(nWir), wirePitch_(wirP), 
	efficiency_(eff), backGround_(bgk), timeGate_(tGate), 
 	hasDrift_(false), 
	driftVelocity_(0.0),
	driftT0_(0.0), 
	doubleHitResolution_(0.0), 
	spaceResolution_(0.0), 
	timeSlice_(0.0),
	hasDeadspace_(false),
	dId_(0),
	dType_(COMGEANT::BOX),
	dParam_(COMGEANT::nParam(COMGEANT::BOX), 0.0),
	dCenter_(0.0, 0.0, 0.0),
	dRotM_(3,3,1)
	{	
	this->name(name);
}


GdbDetTableCont::GdbDetTableCont(const GdbDetTableCont& tableCont) :
	id_(		tableCont.id()),
	unit_(		tableCont.unit()),
	type_(		tableCont.type()),
	radLen_(	tableCont.radLen()),
	solidtype_(	tableCont.solidtype()),
	param_(		tableCont.param()),
//	size_(		tableCont.size()),
	center_(	tableCont.center()),
	rotM_(		tableCont.rotation()), 
	wireDistance_(	tableCont.wirD()),
	angle_(			tableCont.angle()), 
	wireNumber_(	tableCont.nWir()),
	wirePitch_(		tableCont.wirP()), 
	efficiency_(	tableCont.eff()),
	backGround_(	tableCont.bgk()),
	timeGate_(		tableCont.tGate()), 
 	hasDrift_(	tableCont.hasDrift()), 
	driftVelocity_(	tableCont.driftV()),
	driftT0_(		tableCont.T0()), 
	doubleHitResolution_(	tableCont.dHitRes()), 
	spaceResolution_(		tableCont.sRes()), 
	timeSlice_(				tableCont.tSlic()),
	hasDeadspace_(		tableCont.hasDeadSpace()),
	dId_(				tableCont.dId()),
	dType_(				tableCont.dType()),
	dParam_(			tableCont.dParam()),
	dCenter_(			tableCont.dCenter()),
	dRotM_(				tableCont.dRotM())
	{
	this->name(tableCont.name());
}


GdbDetTableCont::~GdbDetTableCont(){
}

GdbDetTableCont& GdbDetTableCont::operator=(const GdbDetTableCont& tableCont){
	if(this != &tableCont){
		this->id(		tableCont.id());
		this->name(			tableCont.name());
		this->unit(		tableCont.unit());
		this->type(		tableCont.type());
		this->radLen(	tableCont.radLen());
		this->solidtype(	tableCont.solidtype());
		this->param(	tableCont.param());
//		this->size(		tableCont.size());
		this->center(	tableCont.center());
		this->rotation(		tableCont.rotation());
		this->wirD(	tableCont.wirD());
		this->angle(			tableCont.angle());
		this->nWir(	tableCont.nWir());
		this->wirP(		tableCont.wirP()); 
		this->eff(	tableCont.eff());
		this->bgk(	tableCont.bgk());
		this->tGate(		tableCont.tGate());
		this->hasDrift(	tableCont.hasDrift());
		this->driftV(	tableCont.driftV());
		this->T0(		tableCont.T0());
		this->dHitRes(	tableCont.dHitRes());
		this->sRes(		tableCont.sRes());
		this->tSlic(				tableCont.tSlic());
		this->hasDeadSpace( tableCont.hasDeadSpace());
		this->dId(				tableCont.dId()),
		this->dType(				tableCont.dType());
		this->dParam(			tableCont.dParam());
		this->dCenter(			tableCont.dCenter());
		this->dRotM(				tableCont.dRotM());
	}
	return *this;
}

void GdbDetTableCont::name( const char* name ){
	strncpy(name_, name, GdbDetTableCont::MAXCHAR);
}

void GdbDetTableCont::dumpDetTableHeader(ostream& os){
	os	<< setfill( ' ' )
		<< " det"
		<< setw(5)	<< "ID"
		<< setw(6)	<< "NAME"
		<< setw(3)	<< "UNIT"
		<< setw(5)	<< "TYPE"
		<< setw(9)	<< "Rad.Len.";
				
	os	<< setw(4)	<< "N";
	for( unsigned int i=0; i < this->refParam().size(); i++){
		os	<< setw(8)	<< i;
	}

	os	<< setw(11)	<< "Xcm"
		<< setw(11)	<< "Ycm"
		<< setw(11)	<< "Zcm"
		<< setw(7)	<< "RMid"
		<< setw(14)	<< "WirD"
		<< setw(9)	<< "Angle"
		<< setw(6)	<< "N.wir"
		<< setw(9)	<< "W.Pitch"
		<< setw(6)	<< "Effic."
		<< setw(6)	<< "BG"
		<< setw(7)	<< "T-Gate"
		<< setw(11)	<< "Drift-v."
		<< setw(6)	<< "T0"
		<< setw(6)	<< "2hitRes."
		<< setw(7)	<< "Spa.Res."
		<< setw(6)	<< "TimeSlice"
		<< endl;
}
void GdbDetTableCont::setDriftInfo(	
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

void GdbDetTableCont::dumpDetTable(ostream& os){
	string name = this->name();

	int rotMId(1);

	os	<< setfill( ' ' )	<< setiosflags(ios::fixed)
		<< " det"
		<< setw(5)	<< this->id()
		<< setw(6)	<< this->name()
		<< setw(3)	<< this->unit()
		<< setw(5)	<< this->type()
		<< setprecision(2)	
		<< setw(9)	<< this->radLen();

	os	<< setw(4)	<< this->refParam().size();
	for( unsigned int i=0; i < this->refParam().size(); i++){
		os	<< setprecision(3)	<< setw(8)	<< this->refParam()[i];
	}
	
	os	<< setprecision(4)	
		<< setw(11)	<< this->center().x()
		<< setw(11)	<< this->center().y()
		<< setw(11)	<< this->center().z()
		<< setw(7)	<< rotMId
		<< setw(14)	<< this->wirD()
		<< setprecision(3)	
		<< setw(9)	<< this->angle()
		<< setw(6)	<< this->nWir()
		<< setprecision(4)	
		<< setw(9)	<< this->wirP()
		<< setprecision(3)	
		<< setw(6)	<< this->eff()
		<< setw(6)	<< this->bgk()
		<< setprecision(1)	
		<< setw(7)	<< this->tGate();
		
	if(this->hasDrift())
	os	<< setprecision(4)	
		<< setw(11)	<< this->driftV()
		<< setprecision(1)	
		<< setw(6)	<< this->T0()
		<< setw(6)	<< this->dHitRes()
		<< setw(7)	<< this->sRes()
		<< setw(6)	<< this->tSlic();
}		

void GdbDetTableCont::deadZone(	
		const int& type, const vector<double>& param,
		const HepPoint3D& center, const HepMatrix& rot){

		this->hasDeadSpace(true);
		this->dId( this->id() ); // In case there is only one dead space
		this->dType(type);
		this->dParam(param);
		this->dCenter(center);
		this->dRotM(rot);

}

void GdbDetTableCont::dumpDeadSpace(ostream& os){
	if(this->hasDeadSpace()){
		os	<< setfill( ' ' )	<< setiosflags(ios::fixed)
			<< " dead"
			<< setw(6)	<< this->dId()
			<< setw(6)	<< this->name()
			<< setw(3)	<< this->unit()
			<< setw(3)	<< this->dType();
			
		for(unsigned int i = 0; i < this->dParam().size(); i++)
			os	<< setprecision(3)
				<< setw(8)	<< this->dParam()[i];
		
		os	<< setprecision(3)
			<< setw(11)	<< this->dCenter().x()
			<< setw(11)	<< this->dCenter().y()
			<< setw(11)	<< this->dCenter().z();

		int rotDMID(0);		
		os	<< setw(3)	<< rotDMID;

	}
}


ostream& operator<<(ostream& os, GdbDetTableCont& tableCont){
	tableCont.dumpDetTable(os);
	os << endl;
	tableCont.dumpDeadSpace(os);
	os << endl;
	return os;
}

HepPoint3D GdbDetTableCont::size() const {
//	if(this->solidtype() == COMGEANT::TUBE)
//		return HepPoint3D(	this->refParam()[0], 
//						this->refParam()[1],
//						this->refParam()[2]);
//

	vector<double> p = this->param();

	return HepPoint3D( p[0], p[1], p[2]);

}
