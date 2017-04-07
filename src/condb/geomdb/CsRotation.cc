/*!
  \file		CsRotation.cc
  \brief	Rotation matrix object implementation file
  \author	$Author: miyachi $
  \version	$Revision: 1.2 $
  \date		$Date: 2000/06/07 08:36:43 $
*/
#if __GNUG__ >= 2
#  pragma implementation
#endif

#include "CsRotation.h"

//---- CsRotation ---------------------------------------------------------
CsRotation::CsRotation() :
	xx_(1.0), xy_(0.0), xz_(0.0),
	yx_(0.0), yy_(1.0), yz_(0.0),
	zx_(0.0), zy_(0.0), zz_(1.0) {
}

CsRotation::CsRotation(
	const double& xx, const double& xy, const double& xz,
	const double& yx, const double& yy, const double& yz,
	const double& zx, const double& zy, const double& zz
 ) :
	xx_(xx), xy_(xy), xz_(xz),
	yx_(xy), yy_(yy), yz_(yz),
	zx_(xz), zy_(zy), zz_(zz) {
}

CsRotation::CsRotation(const HepMatrix& m) :
	xx_(m(1,1)), xy_(m(1,2)), xz_(m(1,3)),
	yx_(m(2,1)), yy_(m(2,2)), yz_(m(2,3)),
	zx_(m(3,1)), zy_(m(3,2)), zz_(m(3,3)) {
}


CsRotation::CsRotation(const CsRotation& m) :
	xx_(m.xx()), xy_(m.xy()), xz_(m.xz()),
	yx_(m.yx()), yy_(m.yy()), yz_(m.yz()),
	zx_(m.zx()), zy_(m.zy()), zz_(m.zz())	{
}

CsRotation::~CsRotation(){
}

CsRotation& CsRotation::operator=(const CsRotation& m){
	if( this != &m ){
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


CsRotation::operator HepMatrix () const {
	HepMatrix m(3,3);
	m(1,1) = this->xx();m(1,2) = this->xy();m(1,3) = this->xz();
	m(2,1) = this->yx();m(2,2) = this->yy();m(2,3) = this->yz();
	m(3,1) = this->zx();m(3,2) = this->zy();m(3,3) = this->zz();
	return m;
}

ostream& operator<<(ostream& os, const CsRotation& rotM){
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
