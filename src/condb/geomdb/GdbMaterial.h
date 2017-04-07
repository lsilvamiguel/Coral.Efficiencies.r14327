//-*-Mode: C++;-*-
#ifndef _GdbMaterial_h_
#define _GdbMaterial_h_
/*!
  \file		GdbMaterial.h
  \brief	Geometry DB material description object definition file
  \author	$Author: miyachi $
  \version	$Revision: 1.2 $
  \date		$Date: 2000/06/07 08:36:46 $
*/

#if __GNUG__ >= 2
#  pragma interface
#endif

#include "GdbSTD.h"

//---- GdbMaterial -----------------------------------------------------------
/*!
  \class	GdbMaterial
  \brief	Material description object for Geometry DB.
*/
class GdbMaterial {
public:
  GdbMaterial();	//!< default constructor

/*!
  \fn		GdbMaterial( const char* name, const double density,
			  const double radiationLength);
  \param	name is a pointer to char
  \param	density is a double
  \param	radiationLength is a double
  \brief	constructor with name, density and radiation length of material
*/
  GdbMaterial( const char* name, const double density,
	const double radiationLength);

//! copy constructor
  GdbMaterial(const GdbMaterial& material);

//! destructor
  ~GdbMaterial();

//! assignment operator
  GdbMaterial& operator=(const GdbMaterial& material);


  const char* materialName() const ;	//!< get name
  double radLen() const;				//!< get radiation length
  double density() const;				//!< get density

//! operator >> overload
  friend istream& operator>>(istream& is, GdbMaterial& material);

//! output operator overloading
  friend ostream& operator<<(ostream& os, const GdbMaterial& material);

protected:
  void materialName(const char* name);	//!< set name
  void radLen(const double& radLen);	//!< set radiation length
  void density(const double& dens);		//!< set density

private:
  static const char maxNameLength = 64;	//!< maximum name length
  char materialName_[maxNameLength];	//!< name
  double density_;						//!< densisty
  double radLen_;						//!< radiation length

};

#endif
