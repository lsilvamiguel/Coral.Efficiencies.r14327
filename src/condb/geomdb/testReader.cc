#include "CsSTD.h"
#include "CsInit.h"
#include "CsGeoDbReader.h"

int main(int argc, char *argv[]){
	CsInit* init = CsInit::Instance(argc, argv);

	CsRegistrySing *reg_ = CsRegistrySing::Instance();

	CsGeoDbReader* pReader = new CsGeoDbReader();

	string gversion("v005rel1");
	CsTime runTime(2000, 4, 1, 0, 0, 0);

	list<CsGeoDbDetector*> detectorTable = 
		pReader->getDetectorTable(runTime, gversion);

	list<CsGeoDbDetector*>::iterator iDet;
	for(	iDet =  detectorTable.begin();
		iDet != detectorTable.end();
		iDet++){
		
		cout << *(*iDet) << endl;
		
		
	}

	
	reg_->callEndMethods();

	exit(0);
}
