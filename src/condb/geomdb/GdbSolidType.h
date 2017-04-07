//-*-Mode: C++;-*-
#ifndef _GdbSolidType_h_
#define _GdbSolidType_h_
/*!
  \file		GdbSolidType.h
  \brief	Geometry DB solid object definition file
  \author	$Author: miyachi $
  \version	$Revision: 1.2 $
  \date		$Date: 2000/06/07 08:36:47 $
*/

#if __GNUG__ >= 2
#  pragma interface
#endif

#include "GdbSTD.h"
#include "HepODBMS/odbms/HepODBMS.h"

//---- GdbSolidType -----------------------------------------------------------
/*!
  \class	GdbSolidType
  \brief	Geometry DB solid object.
			In the object which is consisted by GdbPerStation
			can not use string, and stl container and so on.
			And olso d_String, d_Varray must be avoid to use
			since it is nor allowed to include any d_Varray 
			inside d_Varray. Then in this object, 
			char array and double array sizes are fixed by constant.
*/
class GdbSolidType {
public:

  GdbSolidType();	//!< default constructor
/*!
  \fn		GdbSolidType(const char* type, const long& nParam);
  \param	type is a const char
  \param	nParam is a const long
  \brief	construct GdbSolidType objet with solid type name and 
			number of parameters. 
			(In future relese, this could use COMGEANT naming space
			to determin solid name and number of parameter.) 
*/
  GdbSolidType(const char* type, const long& nParam);

//! copy constructor
  GdbSolidType(const GdbSolidType& solid);

//! destructor
  ~GdbSolidType();
  
//! assignment operator
  GdbSolidType& operator=(const GdbSolidType& solid);

/*!
  \fn		void addParam(const long& index, const double& param);
  \param	index is a long
  \param	param is a double
  \brief	set param as index-th entry
*/
  void addParam(const long& index, const double& param);

  double param(const long& index) const ;	//!< get index-th parameter
  const char* type() const ;				//!< get type name
  long nParam() const ;						//!< get number of parameters
//! get pointer to 1st parameter.
  const double* parameters() const;

//! dump information to ostream
  friend ostream& operator<<(ostream& os, const GdbSolidType& solid);

private:
  static const int maxNameLength = 8;		//!< maximun name length
  d_Char type_[maxNameLength];			//!< type name
  d_Long nParam_;						//!< number of parameters

  static const int maxParameterSize = 28;		//!< max parameter size
  d_Double parameters_[maxParameterSize];	//!< volume parameters

  bool checkIndex(const long& index) const ;	//!< check index

  void initParam();	//!< initialze parameters
  const char* type(const char* type);	//!< set solid type name
  void nParam(const long nParam);	//!< set nuber of paramter

};

#endif
