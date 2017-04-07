/*!
  \file		GdbWireInfo.h
  \brief	Geometry DB wire information object definition file
  \author	$Author: miyachi $
  \version	$Revision: 1.2 $
  \date		$Date: 2000/06/07 08:36:48 $
*/
#if __GNUG__ >= 2
#  pragma implementation
#endif

#include "CsErrLog.h"
#include "GdbWireInfo.h"

//---- GdbWireInfo ---------------------------------------------------------
GdbWireInfo::GdbWireInfo() :
	id_(0), 
	number_(0), pitch_(0.0), angle_(0.0),
	firstPosition_(0.0),
	flag_(-1) {
}

GdbWireInfo::GdbWireInfo(	const int& id,
							const int& number, 
							const double& pitch,
							const double& angle, 
							const double& firstPosition,
							const int& flag=0) :
	id_(id), number_(number), pitch_(pitch), angle_(angle),
	firstPosition_(firstPosition),
	flag_(flag) {
}

GdbWireInfo::GdbWireInfo(	const int& number, 
							const double& pitch,
							const double& angle, 
							const double& firstPosition,
							const int& flag=0) :
	id_(0), number_(number), pitch_(pitch), angle_(angle),
	firstPosition_(firstPosition),
	flag_(flag) {
}

GdbWireInfo::GdbWireInfo(const GdbWireInfo& wInfo) :
	id_(wInfo.id()),
	number_(wInfo.number()), 
	pitch_(wInfo.pitch()), 
	angle_(wInfo.angle()),
	firstPosition_(wInfo.firstPosition()),
	flag_(wInfo.flag()) {
}
    
GdbWireInfo::~GdbWireInfo(){
}

GdbWireInfo& GdbWireInfo::operator=(const GdbWireInfo& wInfo){
	if(this != &wInfo){
		id_		= wInfo.id();
		number_	= wInfo.number();
		pitch_	= wInfo.pitch();
		angle_	= wInfo.angle();
		firstPosition_ = wInfo.firstPosition();
		flag_	= wInfo.flag();
	}
	return *this;
}


int GdbWireInfo::id() const {
	return id_;
}

int GdbWireInfo::number() const {
	return number_;
}

double GdbWireInfo::pitch() const {
	return pitch_;
}

double GdbWireInfo::angle() const {
	return angle_;
}

double GdbWireInfo::firstPosition() const {
	return firstPosition_;
}

int GdbWireInfo::flag() const {
	return flag_;
}

void GdbWireInfo::dumpDetectorsTable(ostream& os){
	
}


ostream& operator<<(ostream& os, const GdbWireInfo& wire){
	os	<< "wire"					<< "\t"
		<< wire.id()			<< "\t"
		<< wire.number()		<< "\t"
		<< wire.pitch()			<< "\t"
		<< wire.angle()			<< "\t"
		<< wire.firstPosition()	<< "\t"
		<< wire.flag();
	return os;
}

istream& operator>>(istream& is, GdbWireInfo& wire){
	if(is.good()){
		int id; is >> id;
		wire.id(id);
	} else {
		CsErrLog::Instance()->mes(elWarning,
			"wire id number pitch angle fist position [flag]");
		return is;
	}

	if(is.good()){
		int number;
		is >> number;
		wire.number(number);

		if(is.good()){
			double pitch;
			is >> pitch;
			wire.pitch(pitch);
			
			if(is.good()){
				double angle;
				is >> angle;
				wire.angle(angle);

				if(is.good()){
					double fPos;
					is >> fPos;
					wire.firstPosition(fPos);

					if(is.good()){
						int flag;
						is >> flag;
						wire.flag(flag);

					} else {
						CsErrLog::Instance()->mes(elWarning,
							"wire number pitch angle fist position [flag]");
					}
				} else {
					CsErrLog::Instance()->mes(elWarning,
						"wire number pitch angle [fist position] [flag]");
				}
			} else {
				CsErrLog::Instance()->mes(elWarning,
					"wire number pitch [angle] [fist position] [flag]");
			}
		} else {
//			string mes = "wire ";
//			char buf[5];
//			itoa(number, buf, 5);
//			mes.append( buf );
//			mes += "[pitch] [angle] [fist position] [flag]";
			CsErrLog::Instance()->mes(elWarning,
				"wire number [pitch] [angle] [fist position] [flag]");

		}
	} else {
		CsErrLog::Instance()->mes(elWarning,
			"wire [number] [pitch] [angle] [fist position] [flag]");
	}	
	return is;
}

//protected:
void GdbWireInfo::id(const int& id){
	id_ = id;
}


void GdbWireInfo::number(const int& number){
	number_ = number;
}

void GdbWireInfo::pitch(const double& pitch){
	pitch_ = pitch;
}

void GdbWireInfo::angle(const double& angle){
	angle_ = angle;
}

void GdbWireInfo::firstPosition(const double& firstPosition){
	firstPosition_ = firstPosition;
}

void GdbWireInfo::flag(const int& flag){
	flag_ = flag;
}

