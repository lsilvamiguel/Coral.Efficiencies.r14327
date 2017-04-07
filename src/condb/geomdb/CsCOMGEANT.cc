#if __GNUG__ >= 2
#  pragma implementation
#endif
/*!
  \file		CsCOMGEANT.cc
  \brief	GEANT 3 type definition namespace, COMGEANT. (implementation)
  \author	$Author: miyachi $
  \version	$Revision: 1.2 $
  \date		$Date: 2000/06/07 08:36:42 $
*/


#include "CsCOMGEANT.h"

#include <algorithm>
#include <cctype>


int COMGEANT::nParam(SOLID solidtype){
  int rVal;
  switch( solidtype ){
	case BOX:	rVal =  3;	break;
	case TRD1:	rVal =  4;	break;
	case TRD2:	rVal =  5;	break;
	case TRAP:	rVal = 11;	break;
	case TUBE:	rVal =  3;	break;
	case TUBS:	rVal =  5;	break;
	case CONE:	rVal =  5;	break;
	case CONS:	rVal =  7;	break;
	case SPHE:	rVal =  6;	break;
	case PARA:	rVal =  6;	break;
	case PGON:	rVal =  7;	break;
	case PCON:	rVal =  6;	break;
	case ELTU:	rVal =  3;	break;
	case HYPE:	rVal =  4;	break;
	case GTRA:	rVal = 12;	break;
	default:		rVal =  3;	break;				
  }
  return rVal;
}

string COMGEANT::solidType( COMGEANT::SOLID solidtype ){

  switch(solidtype){
	case BOX:	return "BOX";
	case TRD1:	return "TRD1";
	case TRD2:	return "TRD2";
	case TRAP:	return "TRAP";
	case TUBE:	return "TUBE";
	case TUBS:	return "TUBS";
	case CONE:	return "CONE";
	case CONS:	return "CONS";
	case SPHE:	return "SPHE";
	case PARA:	return "PARA";
	case PGON:	return "PGON";
	case PCON:	return "PCON";
	case ELTU:	return "ELTU";
	case HYPE:	return "HYPE";
	case GTRA:	return "GTRA";
	default:		return "BOX";		
  }
}

int COMGEANT::nParam(const string& type){
  return COMGEANT::nParam( COMGEANT::solidId(type) );
}

COMGEANT::SOLID COMGEANT::solidId(const string& type){
  
  string typeUpper = type;
  transform( typeUpper.begin(), typeUpper.end(), typeUpper.begin(), toupper );

  if( typeUpper == "TRD1" )	return COMGEANT::TRD1;
  if( typeUpper == "TRD2" )	return COMGEANT::TRD2;
  if( typeUpper == "TRAP" )	return COMGEANT::TRAP;
  if( typeUpper == "TUBE" )	return COMGEANT::TUBE;
  if( typeUpper == "TUBS" )	return COMGEANT::TUBS;
  if( typeUpper == "CONE" )	return COMGEANT::CONE;
  if( typeUpper == "CONS" )	return COMGEANT::CONS;
  if( typeUpper == "SPHE" )	return COMGEANT::SPHE;
  if( typeUpper == "PARA" )	return COMGEANT::PARA;
  if( typeUpper == "PGON" )	return COMGEANT::PGON;
  if( typeUpper == "PCON" )	return COMGEANT::PCON;
  if( typeUpper == "ELTU" )	return COMGEANT::ELTU;
  if( typeUpper == "HYPE" )	return COMGEANT::HYPE;
  if( typeUpper == "GTRA" )	return COMGEANT::GTRA;
  
  return COMGEANT::BOX;

}
