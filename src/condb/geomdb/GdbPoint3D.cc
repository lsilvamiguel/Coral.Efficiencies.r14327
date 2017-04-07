/*!
  \file		GdbPoint3D.cc
  \brief	Geometry DB 3-dimensional vector object implementation file
  \author	$Author: miyachi $
  \version	$Revision: 1.2 $
  \date		$Date: 2000/06/07 08:36:46 $
*/
#if __GNUG__ >= 2
#  pragma implementation
#endif

#include "GdbPoint3D.h"

//---- GdbPoint3D ---------------------------------------------------------

GdbPoint3D::GdbPoint3D() :
	x_(0.0), y_(0.0), z_(0.0) {
}

GdbPoint3D::GdbPoint3D( const double& x, const double& y, const double& z ) :
	x_(x), y_(y), z_(z)
	{
}

GdbPoint3D::GdbPoint3D( const GdbPoint3D& point ) : 
	x_(point.x()), y_(point.y()), z_(point.z())
	{
}

GdbPoint3D::~GdbPoint3D(){
}

GdbPoint3D& GdbPoint3D::operator=(const GdbPoint3D& point){
	if( this != & point ){
		this->x( point.x() );
		this->y( point.y() );
		this->z( point.z() );
	}
	return *this;
}

double GdbPoint3D::x() const { return x_;}
double GdbPoint3D::y() const { return y_;}
double GdbPoint3D::z() const { return z_;}

double GdbPoint3D::x(const double& x) { x_ = x ; return x_;}	//!< return x after set x
double GdbPoint3D::y(const double& y) { y_ = y ; return y_;}	//!< return y after set y
double GdbPoint3D::z(const double& z) { z_ = z ; return z_;}	//!< return z after set z

GdbPoint3D GdbPoint3D::operator+(const GdbPoint3D& point){
	return GdbPoint3D(
		this->x() + point.x(),
		this->y() + point.y(),
		this->z() + point.z()
	);
}

GdbPoint3D GdbPoint3D::operator-(const GdbPoint3D& point){
	return GdbPoint3D(
		this->x() - point.x(),
		this->y() - point.y(),
		this->z() - point.z()
	);
}

GdbPoint3D GdbPoint3D::operator*(const double& factor){
	return GdbPoint3D( factor * this->x(), factor * this->y(), factor* this->z());	
}

GdbPoint3D operator*(const double& factor, const GdbPoint3D& point){
	return GdbPoint3D( factor * point.x(), factor * point.y(), factor* point.z());
}

ostream& operator<<( ostream& os, const GdbPoint3D& point ){
	os	<< point.x()	<< '\t'
		<< point.y()	<< '\t'
		<< point.z()	;
	return os;
}
