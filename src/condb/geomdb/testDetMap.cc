#include <CsSTD.h>
#include <CsInit.h>
#include "CsDetectorFile.h"
#include "CsDetectorMap.h"


int main(int argc, char *argv[]){

	CsInit* init = CsInit::Instance(argc, argv);

	CsDetectorMap mapping("mapping.dat");
	
	CsDetectorFile detectorFile( *(init->getDetectorTable()) );
	list<CsDetector> detectorList = detectorFile.getDetector();

	for(list<CsDetector>::iterator 
		idet =  detectorList.begin();
		idet != detectorList.end();
		idet++){

		string detector( mapping.detector((*idet).getName()) );
		cout	<< (*idet).getName() << " <--- "
				<< detector << " in "
				<< mapping.container(detector) << endl;

	}

	return 0;
}
