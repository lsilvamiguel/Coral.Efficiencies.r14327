#if __GNUG__ >= 2
#  pragma implementation
#endif

/*!
   \file    CsTrigHodDbReader.cc
   \brief   Compass Some CDB reader class implementation file.
   \author  Vassili Motchalov
   \version $Revision: 1.1 $
   \date    $Date: 2000/04/05 11:31:29 $
*/

#include "CsTrigHodDbReader.h"

CsTrigHodDbReader::CsTrigHodDbReader( string dbName) : 
CsCondDbReader(dbName) {}


CsTrigHodDbReader::~CsTrigHodDbReader(){
	CsErrLog::Instance()->mes(elDebugging, "TrigHod Reader is going to die.");
}









