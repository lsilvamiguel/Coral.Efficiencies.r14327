//-*-Mode: C++;-*-
#ifndef _GdbPoint3D_h_
#define _GdbPoint3D_h_
/*!
  \file		GdbPoint3D.h
  \brief	Geometry DB 3-dimensional vector object definition file
  \author	$Author: miyachi $
  \version	$Revision: 1.2 $
  \date		$Date: 2000/06/07 08:36:46 $
*/

#if __GNUG__ >= 2
#  pragma interface
#endif

#include "GdbSTD.h"
#include "HepODBMS/odbms/HepODBMS.h"

//---- GdbPoint3D -----------------------------------------------------------
/*!
  \class	GdbPoint3D
  \brief	Geometry DB -dimensional vector object. 
			Since CLHEP objects are not persistent capable object or 
			not suitable (warning?), 
			one can store 3-dimentional vector object in DB using this.
*/
class GdbPoint3D {
public:
//! a default constructor
  GdbPoint3D();
/*!
  \fn		GdbPoint3D( const double& x, const double& y, const double& z );
  \param	x is a double	
  \param	y is a double	
  \param	z is a double	
  \brief	construct GdbPoint3D object using the given x, y and z.
*/
  GdbPoint3D( const double& x, const double& y, const double& z );
//! a copy constructor
  GdbPoint3D( const GdbPoint3D& point );

//!< destructor
  ~GdbPoint3D();

//!< assignment operator
  GdbPoint3D& operator=(const GdbPoint3D& point);

  double x() const ;	//!< return x
  double y() const ;	//!< return y
  double z() const ;	//!< return z

  double x(const double& x) ;	//!< return x after set x
  double y(const double& y) ;	//!< return y after set y
  double z(const double& z) ;	//!< return z after set z

//! dump information to ostream
  friend ostream& operator<<( ostream& os, const GdbPoint3D& point );

//! plus operator overloading
  GdbPoint3D operator+(const GdbPoint3D& point);
//! minus operator overloading
  GdbPoint3D operator-(const GdbPoint3D& point);
//! times operator overloading
  GdbPoint3D operator*(const double& factor);
//! scaling operator overloading
  friend GdbPoint3D operator*(const double& factor, const GdbPoint3D& point);

private:
  d_Double x_;	//!< x
  d_Double y_;	//!< y
  d_Double z_;	//!< z
};

#endif
