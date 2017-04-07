/*!
   \file    testUpdater.h
   \brief   Test Program for CsBmsDbUpdater objects.
   \author  Yoshiyuki Miyachi
   \version $Revision: 1.3 $
   \date    $Date: 2000/04/22 07:49:43 $
*/
#include "CsSTD.h"
#include "CsInit.h"
#include "CsRegistrySing.h"
//#include "CsOpt.h"
#include "CsBmsDbUpdater.h"

int main(int argc, char *argv[]){

	CsInit* init = CsInit::Instance(argc, argv);
//	CsOpt*  coralOpt = CsOpt::Instance(argc,argv);
	cout	<< argv[0]
			<< "\t###### Create Register. ######" << endl;
	CsRegistrySing *reg_ = CsRegistrySing::Instance();

//	In current schema, handler object has to be create in stack memory with new operator,
//	and it will be deleted when TransactionManager desappear from memory (reg_->callEndMethods()).
//	You should not delete this object by your self.
	cout	<< argv[0]
			<< "\t###### Create Updater. ######" << endl;
	CsBmsDbUpdater* pBmsUpdater = new CsBmsDbUpdater();

// set interval
//	CsTime beginTime(2000, 1, 1, 0, 0, 0);
//	CsTime   endTime(2010, 1, 1, 0, 0, 0);

// update CDB.	
	cout	<< argv[0]
			<< "\t###### Send Update signal to updater with the interval. ######" << endl;
	pBmsUpdater->updateAll();

// send commit()
	cout	<< argv[0]
			<< "\t###### Send commit signal to updater. ######" << endl;
	pBmsUpdater->commit();

// call end method.
	cout	<< argv[0]
			<< "\t###### Send end method signal to register. ######" << endl;
	reg_->callEndMethods();

// That's all.
	return 0;
}
