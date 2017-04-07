//-*-Mode: C++;-*-
#ifndef _CsGeoDbDetector_h_
#define _CsGeoDbDetector_h_
/*!
  \file		CsGeoDbDetector.h
  \brief	Geometry Detector table interface definition.
  \author	$Author: miyachi $
  \version	$Revision: 1.1 $
  \date		$Date: 2000/06/07 08:36:42 $
*/

#if __GNUG__ >= 2
#  pragma interface
#endif

#include "CsSTD.h"
#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Geometry/Point3D.h>

//---- CsGeoDbDetectpr -----------------------------------------------------------
/*!
	\class	CsGeoDbDetectpr [pABC]
	\brief	Interface to access Geometry Information stored in GEO DB.
*/
class CsGeoDbDetector {
public:

  //! default constructor
  CsGeoDbDetector();

  //! destructor
  virtual ~CsGeoDbDetector();

  virtual int id() const = 0; 				//!< get id
  virtual const char* name() const = 0;		//!< get name
  virtual int unit() const = 0;				//!< get unit
  virtual int type() const = 0;				//!< get type
  virtual double radLen() const = 0;		//!< get radiation length
  virtual HepPoint3D center() const = 0;	//!< get center point
  virtual HepPoint3D size() const = 0;		//!< get size
  virtual vector<double> param() const = 0;	//!< get volume parameter
  virtual HepMatrix rotation() const = 0;	//!< get rotation matrix
  virtual bool hasDeadSpace() const = 0;	//!< dose it has dead space?

  virtual double wirD() const = 0;			//!< get wire dimension?
  virtual double angle() const = 0;			//!< get angle
  virtual int nWir() const = 0;				//!< get numner of wire
  virtual double wirP() const = 0;			//!< get wire pitch
  
  virtual double eff() const = 0;			//!< get efficiency
  virtual double bgk() const = 0;			//!< get background
  virtual double tGate() const = 0;			//!< get time gate

  virtual bool hasDrift() const = 0;		//!< dose it have drift?
  virtual double driftV() const = 0;		//!< get drift velocity
  virtual double T0() const = 0;			//!< get t0

  virtual double dHitRes() const = 0;		//!< get double hit resolution
  virtual double sRes() const = 0;			//!< get space slice
  virtual double tSlic() const = 0;			//!< get time slice

  //! dump information to the given ostream
  friend ostream& operator<<(ostream& os, CsGeoDbDetector& detector);

};

#endif
