//-*-Mode: C++;-*-
#ifndef _GdbVolume_h_
#define _GdbVolume_h_
/*!
  \file		GdbVolume.h
  \brief	Geometry DB volume object definition file
  \author	$Author: miyachi $
  \version	$Revision: 1.2 $
  \date		$Date: 2000/06/07 08:36:48 $
*/

#if __GNUG__ >= 2
#  pragma interface
#endif

#include "GdbSolidType.h"
#include "GdbMaterial.h"
#include "GdbWireInfo.h"
#include "GdbRotation.h"
#include "GdbPoint3D.h"

//---- GdbVolume -----------------------------------------------------------
/*!
  \class	GdbVolume
  \brief	Persistent Geometry Volume Objects which is holded by GdbStation.
		  GdbStation is inherited from HepPersObj.
*/
class GdbVolume {
public:
  GdbVolume();	//!< a default constructor

/*!
  \fn		GdbVolume( const GdbPoint3D& center, const GdbRotation& rotation);
  \param	center is a GdbPoint3D
  \param	rotation is a GdbRotation
  \brief	Construct volume with center and rotation matrix.
*/
  GdbVolume( const GdbPoint3D& center, const GdbRotation& rotation);

/*!
  \fn		GdbVolume( const char* volumeName, const GdbSolidType solid, 
			  const GdbMaterial& material, const GdbPoint3D& center);
  \param	volumeName is a pointer to char
  \param	center is a GdbPoint3D
  \param	solid is a GdbSolidType
  \param	material is a GdbMaterial
  \brief	constructor with 4 parameters
*/
  GdbVolume( const char* volumeName, const GdbSolidType solid, 
	const GdbMaterial& material, const GdbPoint3D& center);

/*!
  \fn		GdbVolume( const char* volumeName, const GdbSolidType solid,
			  const GdbMaterial material, const GdbPoint3D& center,
			  const GdbRotation& rotation, const long& nSubVolume,
			  const bool& active, const GdbWireInfo& wire);
  \param	volumeName is a pointer to char
  \param	solid is a GdbSolidType
  \param	material is a GdbMaterial
  \param	center is a GdbPoint3D
  \param	rotation is GdbRotation
  \param	nSubvolume is long
  \param	active is bool
  \param	wire is GdbWIreInfo
  \brief	constructor with 8 parameters
*/
  GdbVolume( const char* volumeName, const GdbSolidType solid,
	const GdbMaterial material, const GdbPoint3D& center,
	const GdbRotation& rotation, const long& nSubVolume,
	const bool& active, const GdbWireInfo& wire);

/*!
  \fn		GdbVolume(const GdbVolume& gVolume);
  \param	gVolume is a GdbVolume
  \brief	Copy constructor
*/
  GdbVolume(const GdbVolume& gVolume);


  virtual ~GdbVolume();	//!< destructor

  GdbVolume& operator=(const GdbVolume& gVolume);	//!< assignment operator
  bool operator<(const GdbVolume& gVolume) const ;	//!< less than operator
  bool operator==(const GdbVolume& gVolume) const ;	//!< equal to operator
  bool operator==(const char* name) const ;			//!< equal to operator 

  const char* volumeName() const ;	//!< get volume name
  GdbSolidType solid() const;		//!< get solid type
  GdbMaterial material() const ;	//!< get material
  GdbPoint3D center() const ;		//!< get center point
  GdbRotation rotMatrix() const ;	//!< ger rotation matrix
  GdbWireInfo wire() const;			//!< get wire information
  bool isActive() const;			//!< get active flag

  void isActive(const bool& active);			//!< set active flag
  void volumeName(const char* volumeName);		//!< set volume name
  void solid(const GdbSolidType& solid);		//!< set solid type
  void material(const GdbMaterial& material);	//!< set material information
  void center(const GdbPoint3D& center);		//!< set center point
  void rotMatrix(const GdbRotation& rotMartix);	//!< set rotation matrix
  void wire(const GdbWireInfo& wire);			//!< set wire information

  long nSubVolume() const;				//!< get number of subvolumes
  void nSubVolume(const long& nSub) ;	//!< set number of subvolumes

//! dump information to ostream
  friend ostream& operator<<(ostream& os, GdbVolume& volume);

//! dump information to ostream
  virtual void dump(ostream& os, const char* pTab = NULL);

/*!
  \struct	ActivePlane
  \brief	active plane information
*/
  struct ActivePlane {
	long 		index;		//!< index 
	GdbPoint3D  center;		//!< center point in HALL
	GdbRotation rotation;	//!< rotation matrix to HALL
  };

private:
  static const int maxNameLength = 8;		//!< max name length
  d_Char volumeName_[maxNameLength];	//!< name
  GdbSolidType solid_;					//!< solid type
  GdbMaterial material_;				//!< material information

  GdbPoint3D center_;					//!< center in MV
  GdbRotation rotMatrix_;				//!< rotation to MV

  d_Long nSubVolume_;					//!< number of subvolumes
  d_Boolean isActive_;					//!< active area flag
  GdbWireInfo wire_;					//!< wire information

};

#endif
