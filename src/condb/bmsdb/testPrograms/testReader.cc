/*!
   \file    testReader.h
   \brief   Test Program for CsBmsDbReader objects.
   \author  Yoshiyuki Miyachi
   \version $Revision: 1.3 $
   \date    $Date: 2000/04/22 07:49:43 $
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
	CsTime runTime(2000, 4, 1, 0, 0, 0);
	CsBMSconstants tConsts = pBmsDbReader->getConstants(runTime);

//	print out information
	tConsts.prntBMSconst();

//	get all constants
	list<CsBMSconstants> constList = pBmsDbReader->getAllConstants(runTime);
	for(	list<CsBMSconstants>::iterator it = constList.begin();
			it != constList.end();
			it++){
		(*it).prntBMSconst();
	}

	reg_->callEndMethods();
	
// That is all.
	return 0;
}
