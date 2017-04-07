#ifndef GdbRotation_h
#define GdbRotation_h
/*!
  \file		GdbPoint3D.h
  \brief	Geometry DB rotation matrix object definition file
  \author	$Author: miyachi $
  \version	$Revision: 1.2 $
  \date		$Date: 2000/06/07 08:36:47 $
*/
#include "GdbSTD.h"
#include "GdbPoint3D.h"
/*!
  \class	GdbRotation
  \brief	Geometry DB rotation matrix object
*/
class GdbRotation {
public:
//!	constructor
  GdbRotation();
/*!
  \fn		GdbRotation(
			  const double& xx, const double& xy, const double& xz,
			  const double& yx, const double& yy, const double& yz,
			  const double& zx, const double& zy, const double& zz
		   );
  \param	xx is a double
  \param	xy is a double
  \param	xz is a double
  \param	yx is a double
  \param	yy is a double
  \param	yxz is a double
  \param	zx is a double
  \param	zy is a double
  \param	zz is a double
  \brief	construct from 9 parameters
*/
  GdbRotation(
	  const double& xx, const double& xy, const double& xz,
	  const double& yx, const double& yy, const double& yz,
	  const double& zx, const double& zy, const double& zz
   );

//! copy constructor
  GdbRotation(const GdbRotation& rot);

//! destructor
  ~GdbRotation();

//! assignment operator
  GdbRotation& operator=(const GdbRotation& rot);

  double xx() const;	//!< return xx
  double xy() const;	//!< return xy
  double xz() const;	//!< return xz
  double yx() const;	//!< return yx
  double yy() const;	//!< return yy
  double yz() const;	//!< return yz
  double zx() const;	//!< return zx
  double zy() const;	//!< return zy
  double zz() const;	//!< return zz

  double xx(const double& xx) ;	//!< set xx
  double xy(const double& xy) ;	//!< set xy
  double xz(const double& xz) ;	//!< set xz
  double yx(const double& yx) ;	//!< set yx
  double yy(const double& yy) ;	//!< set yy
  double yz(const double& yz) ;	//!< set yz
  double zx(const double& zx) ;	//!< set zx
  double zy(const double& zy) ;	//!< set zy
  double zz(const double& zz) ;	//!< set zz

//! dump information to ostream
  friend ostream& operator<<(ostream& os, const GdbRotation& rotM);

//--------------------- operations
  GdbRotation operator+(const GdbRotation& rot);	//!< plus operator
  GdbRotation operator-(const GdbRotation& rot);	//!< minus operator
  GdbRotation operator*(const GdbRotation& rot);	//!< times operator
/*!
  \fn		GdbPoint3D operator*(const GdbPoint3D& point);
  \param	point is a GdbPoint
  \return	a GdbPoint
  \brief	return (*this) * (point)
*/
  GdbPoint3D operator*(const GdbPoint3D& point);


private:
  double xx_;	//!< xx
  double xy_;	//!< xy
  double xz_;	//!< xz
  double yx_;	//!< yx
  double yy_;	//!< yy
  double yz_;	//!< yz
  double zx_;	//!< zx
  double zy_;	//!< zy
  double zz_;	//!< zz
};

#endif \\GdbRotation_h
