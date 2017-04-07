#include "CsSTD.h"
#include "CsInit.h"
#include "CsGeoDbUpdater.h"

int main(int argc, char *argv[]){
	CsInit* init = CsInit::Instance(argc, argv);

	CsRegistrySing *reg_ = CsRegistrySing::Instance();

	CsGeoDbUpdater* pUpdater = new CsGeoDbUpdater();
	
	CsTime beginTime(2000, 1, 1, 0, 0, 0);
	CsTime   endTime(2010, 1, 1, 0, 0, 0);
	
	pUpdater->update( beginTime, endTime);


	pUpdater->commit();


	reg_->callEndMethods();

	exit(0);
}
