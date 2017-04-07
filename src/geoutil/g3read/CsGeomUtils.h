//-*-Mode: C++;-*-
#ifndef _CsGeomUtils_h_
#define _CsGeomUtils_h_
/*!
  \file		CsGeomUtils.h
  \brief	Definition and implementation of several global functions.
  \author	$Author: miyachi $
  \version	$Revisopn$
  \date		$Date: 2000/06/06 10:24:33 $
*/

#if __GNUG__ >= 2
#  pragma interface
#endif

// CLHEP Geometry package
#include <CLHEP/Geometry/Point3D.h>
#include <CLHEP/Geometry/Plane3D.h>
#include <CLHEP/Matrix/Matrix.h>

/*!
  \fn		Hep3Vector operator*( const HepMatrix& M, const Hep3Vector& a)
  \param	M is a HepMatrix
  \param	a is a Hep3Vector
  \return	a Hep3Vector
  \brief	return M*a, where M must be M(3,3) and a be 3-dimentional vector
*/
CLHEP::Hep3Vector operator*( const CLHEP::HepMatrix& M, const CLHEP::Hep3Vector& a){
        CLHEP::Hep3Vector b;
        b.setX(M(1,1)*a.x() + M(1,2)*a.y() + M(1,3)*a.z());
        b.setY(M(2,1)*a.x() + M(2,2)*a.y() + M(2,3)*a.z());
        b.setZ(M(3,1)*a.x() + M(3,2)*a.y() + M(3,3)*a.z());
        return b;
};


#endif
