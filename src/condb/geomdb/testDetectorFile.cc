#include <CsSTD.h>
#include "CsInit.h"
#include "CsDetectorFile.h"
int main(int argc, char *argv[]){

	CsInit* init = CsInit::Instance(argc, argv);

	CsDetectorFile detectorFile( *(init->getDetectorTable()) );
	list<CsDetector> detectorList = detectorFile.getDetector();

	for(list<CsDetector>::iterator 
		idet =  detectorList.begin();
		idet != detectorList.end();
		idet++){
		cout << (*idet).getName() << endl;
	}

	return 0;
}
