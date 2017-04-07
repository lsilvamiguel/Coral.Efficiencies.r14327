//-*-Mode: C++;-*-
#ifndef _CsCOMGEANT_h_
#define _CsCOMGEANT_h_
/*!
  \file		CsCOMGEANT.h
  \brief	GEANT 3 type definition namespace, COMGEANT. (definition)
  \author	$Author: miyachi $
  \version	$Revision: 1.2 $
  \date		$Date: 2000/06/07 08:36:42 $
*/

#include "CsSTD.h"

//! Name space for GEANT3 solid type definition
namespace COMGEANT {

/*
  \enum		SOLID
  \brief	Solid type used in GEANT3.
*/
  enum SOLID {
	  BOX	=  1,
	  TRD1	=  2,
	  TRD2	=  3,
	  TRAP	=  4,
	  TUBE	=  5,
	  TUBS	=  6,
	  CONE	=  7,
	  CONS	=  8,
	  SPHE	=  9,
	  PARA	= 10,
	  PGON	= 11,
	  PCON	= 12,
	  ELTU	= 13,
	  HYPE	= 14,
	  GTRA	= 28
  };

//!	return number of parameters needed for the given SOLID
  int nParam(SOLID solidtype);
  
//!	return number of parameters needed for the solid type (string).
  int nParam(const string& type);
  
//!	return SOLID 
  SOLID solidId(const string& type);

//! return solid name
  string solidType( SOLID solidtype );

}

#endif
