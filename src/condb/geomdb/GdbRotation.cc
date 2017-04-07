/*!
  \file		GdbPoint3D.cc
  \brief	Geometry DB rotation matrix object implementation file
  \author	$Author: miyachi $
  \version	$Revision: 1.2 $
  \date		$Date: 2000/06/07 08:36:46 $
*/
#if __GNUG__ >= 2
#  pragma implementation
#endif

#include "GdbRotation.h"

//---- GdbRotation ---------------------------------------------------------
GdbRotation::GdbRotation() :
	xx_(1.0), xy_(0.0), xz_(0.0),
	yx_(0.0), yy_(1.0), yz_(0.0),
	zx_(0.0), zy_(0.0), zz_(1.0) {
}

GdbRotation::GdbRotation(
	const double& xx, const double& xy, const double& xz,
	const double& yx, const double& yy, const double& yz,
	const double& zx, const double& zy, const double& zz
 ) :
	xx_(xx), xy_(xy), xz_(xz),
	yx_(xy), yy_(yy), yz_(yz),
	zx_(xz), zy_(zy), zz_(zz) {
}

GdbRotation::GdbRotation(const GdbRotation& rot) :
	xx_(rot.xx()), xy_(rot.xy()), xz_(rot.xz()),
	yx_(rot.xy()), yy_(rot.yy()), yz_(rot.yz()),
	zx_(rot.xz()), zy_(rot.zy()), zz_(rot.zz()) {
}

GdbRotation::~GdbRotation(){
}

GdbRotation& GdbRotation::operator=(const GdbRotation& m){
	if(this != &m){
		this->xx(m.xx());
		this->xy(m.xy());
		this->xz(m.xz());
		this->yx(m.yx());
		this->yy(m.yy());
		this->yz(m.yz());
		this->zx(m.zx());
		this->zy(m.zy());
		this->zz(m.zz());
	}
	return *this;
}

double GdbRotation::xx() const {return xx_;}
double GdbRotation::xy() const {return xy_;}
double GdbRotation::xz() const {return xz_;}
double GdbRotation::yx() const {return yx_;}
double GdbRotation::yy() const {return yy_;}
double GdbRotation::yz() const {return yz_;}
double GdbRotation::zx() const {return zx_;}
double GdbRotation::zy() const {return zy_;}
double GdbRotation::zz() const {return zz_;}

double GdbRotation::xx(const double& xx) { xx_ = xx; return xx_;}
double GdbRotation::xy(const double& xy) { xy_ = xy; return xy_;}
double GdbRotation::xz(const double& xz) { xz_ = xz; return xz_;}
double GdbRotation::yx(const double& yx) { yx_ = yx; return yx_;}
double GdbRotation::yy(const double& yy) { yy_ = yy; return yy_;}
double GdbRotation::yz(const double& yz) { yz_ = yz; return yz_;}
double GdbRotation::zx(const double& zx) { zx_ = zx; return zx_;}
double GdbRotation::zy(const double& zy) { zy_ = zy; return zy_;}
double GdbRotation::zz(const double& zz) { zz_ = zz; return zz_;}

GdbRotation GdbRotation::operator+(const GdbRotation& rot){
	return GdbRotation(
		this->xx() + rot.xx(), this->xy() + rot.xy(), this->xz() + rot.xz(),
		this->yx() + rot.yx(), this->yy() + rot.yy(), this->yz() + rot.yz(),
		this->zx() + rot.zx(), this->zy() + rot.zy(), this->zz() + rot.zz()
	);
}

GdbRotation GdbRotation::operator-(const GdbRotation& rot){
	return GdbRotation(
		this->xx() - rot.xx(), this->xy() - rot.xy(), this->xz() - rot.xz(),
		this->yx() - rot.yx(), this->yy() - rot.yy(), this->yz() - rot.yz(),
		this->zx() - rot.zx(), this->zy() - rot.zy(), this->zz() - rot.zz()
	);
}

GdbRotation GdbRotation::operator*(const GdbRotation& rot){
	GdbPoint3D x = (*this) * GdbPoint3D( rot.xx(), rot.yx(), rot.zx() );
	GdbPoint3D y = (*this) * GdbPoint3D( rot.xy(), rot.yy(), rot.zy() );
	GdbPoint3D z = (*this) * GdbPoint3D( rot.xz(), rot.yz(), rot.zz() );

	return GdbRotation(
		x.x(),	y.x(),	z.x(),
		x.y(),	y.y(),	z.y(),
		x.z(),	y.z(),	z.z()
	);

}

GdbPoint3D GdbRotation::operator*(const GdbPoint3D& point){
	return GdbPoint3D(
		this->xx()*point.x() + this->xy()*point.y() + this->xz()*point.z(),
		this->yx()*point.x() + this->yy()*point.y() + this->yz()*point.z(),
		this->zx()*point.x() + this->zy()*point.y() + this->zz()*point.z()
	);
}

ostream& operator<<(ostream& os, const GdbRotation& rotM){
	os	<< setprecision(3) ;
	os	<< "rotaion"		<< "\t"	
		<< rotM.xx()	<< "\t"
		<< rotM.xy()	<< "\t"
		<< rotM.xz()	<< "\t"
		<< rotM.yx()	<< "\t"
		<< rotM.yy()	<< "\t"
		<< rotM.yz()	<< "\t"
		<< rotM.zx()	<< "\t"
		<< rotM.zy()	<< "\t"
		<< rotM.zz();
	return os;
}
