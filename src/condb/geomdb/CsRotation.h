#ifndef CsRotation_h
#define CsRotation_h
/*!
  \file		CsRotation.h
  \brief	Rotation matrix object defenition file
  \author	$Author: miyachi $
  \version	$Revision: 1.2 $
  \date		$Date: 2000/06/07 08:36:43 $
*/
#include "GdbSTD.h"
#include <CLHEP/Matrix/Matrix.h>

/*!
	\class	CsRotation
	\brief	Conditional Class for storing rotation matrix
*/
class CsRotation {
public:
//!	default constructor
	CsRotation();
	
/*!
  \fn		CsRotation(
			  const double& xx, const double& xy, const double& xz,
			  const double& yx, const double& yy, const double& yz,
			  const double& zx, const double& zy, const double& zz
			)
  \param	xx is a double
  \param	xy is a double
  \param	xz is a double
  \param	yx is a double
  \param	yy is a double
  \param	yz is a double
  \param	zx is a double
  \param	zy is a double
  \param	zz is a double
  \brief	constructore with nine variables.
*/
	CsRotation(
		const double& xx, const double& xy, const double& xz,
		const double& yx, const double& yy, const double& yz,
		const double& zx, const double& zy, const double& zz
	 );

/*!
  \fn		CsRotation(const HepMatrix& m)
  \param	m is a HepMatrix
  \brief	Copy from the given HepMatrix m
*/
	CsRotation(const HepMatrix& m);

//! copy constructor
	CsRotation(const CsRotation& m);

//!	assignment operator
	CsRotation& operator=(const CsRotation& m);

//! destructor
	~CsRotation();

	HepDouble xx() const {return xx_;}	//!< get xx
	HepDouble xy() const {return xy_;}	//!< get xy
	HepDouble xz() const {return xz_;}	//!< get xz
	HepDouble yx() const {return yx_;}	//!< get yx
	HepDouble yy() const {return yy_;}	//!< get yy
	HepDouble yz() const {return yz_;}	//!< get yz
	HepDouble zx() const {return zx_;}	//!< get zx
	HepDouble zy() const {return zy_;}	//!< get zy
	HepDouble zz() const {return zz_;}	//!< get zz

	HepDouble xx(const double& xx) { xx_ = xx; return xx_;} //!< set and get xx
	HepDouble xy(const double& xy) { xy_ = xy; return xy_;} //!< set and get xy
	HepDouble xz(const double& xz) { xz_ = xz; return xz_;} //!< set and get xz
	HepDouble yx(const double& yx) { yx_ = yx; return yx_;} //!< set and get yx
	HepDouble yy(const double& yy) { yy_ = yy; return yy_;} //!< set and get yy
	HepDouble yz(const double& yz) { yz_ = yz; return yz_;} //!< set and get yz
	HepDouble zx(const double& zx) { zx_ = zx; return zx_;} //!< set and get zx
	HepDouble zy(const double& zy) { zy_ = zy; return zy_;} //!< set and get zy
	HepDouble zz(const double& zz) { zz_ = zz; return zz_;} //!< set and get zz
	
//! dump information to ostream
	friend ostream& operator<<(ostream& os, const CsRotation& rotM);

private:
	HepDouble xx_;	//!< xx
	HepDouble xy_;	//!< xy
	HepDouble xz_;	//!< xz
	HepDouble yx_;	//!< yx
	HepDouble yy_;	//!< yy
	HepDouble yz_;	//!< yz
	HepDouble zx_;	//!< zx
	HepDouble zy_;	//!< zy
	HepDouble zz_;	//!< zz

public:
/*!
  \fn		operator HepMatrix () const
  \brief	cast operator to HepMatrix
*/
	operator HepMatrix () const ;

};

#endif \\CsRotation_h
