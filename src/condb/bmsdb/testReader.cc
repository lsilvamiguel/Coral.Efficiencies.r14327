/*!
   \file    testReader.h
   \brief   Test Program for CsBmsDbReader objects.
   \author  Yoshiyuki Miyachi
   \version $Revision: 1.1 $
   \date    $Date: 2000/06/06 09:48:36 $
*/
#include "CsSTD.h"
#include "CsInit.h"
#include "CsRegistrySing.h"
#include "CsOpt.h"
#include "CsBmsDbReader.h"
#include "CsBMSconstants.h"

int main(int argc, char *argv[]){

//	Current release has serious problem on using CsInit!!
	CsInit* init = CsInit::Instance(argc, argv);
//	CsOpt*  coralOpt = CsOpt::Instance(argc,argv);
	CsRegistrySing *reg_ = CsRegistrySing::Instance();


//	In current schema, handler object has to be create in stack memory with new operator,
//	and it will be deleted when TransactionManager desappear from memory (reg_->callEndMethods()).
//	You should not delete this object by your self.
  CsBmsDbReader* pBmsDbReader = new CsBmsDbReader();

//	Get CsBMSconstants (transient) object at the given time.
  CsTime runTime(1994, 1, 1, 0, 0, 0);
  CsTime beginTime(1999, 1, 1, 0, 0, 0);
  CsTime endTime(2003, 1, 1, 0, 0, 0);
//	CsTime runTime(2000, 4, 1, 0, 0, 0);
//	CsBMSconstants tConsts = pBmsDbReader->getConstants(runTime);
//	print out information
//	tConsts.prntBMSconst();


//	get all constants
//	list<CsBMSconstants> constList = 
//		pBmsDbReader->getAllConstants(beginTime, endTime);
//	for(	list<CsBMSconstants>::iterator it = constList.begin();
//			it != constList.end();
//			it++){
//		(*it).prntBMSconst();
//	}
//
  map< CsInterval, CsBMSconstants> constList = 
	  pBmsDbReader->getAllConstants(beginTime, endTime);

  map< CsInterval, CsBMSconstants>::iterator it;

  for( it = constList.begin(); it != constList.end(); it++){
	  (*it).second.prntBMSconst();
  }

  for( it = constList.begin(); it != constList.end(); it++){
	  cout << (*it).first << endl;
  }



	reg_->callEndMethods();
	
// That is all.
	return 0;
}
