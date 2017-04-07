#include "CsSTD.h"
#include "CsInit.h"
#include "CsOpt.h"
#include "GdbDataFile.h"
#include "GdbPerStation.h"

int main(int argc, char *argv[]){
	CsInit* init = CsInit::Instance(argc, argv);
	CsRegistrySing *reg_ = CsRegistrySing::Instance();

	CsGDBupdater* pGdbUpdater = new CsGDBupdater();

	cout << "try open" << endl;

	pGdbUpdater->open();

	string tag = "GeomDB";
	string key = "datafile";
	string datafile;

	if( ! CsOpt::Instance()->getOpt( tag, key, datafile ) ) {
		cerr << " No file is defined."  << datafile << endl;
		exit(1);
	}

	GdbDataFile gDbFile(datafile);

	cout << "-----------------------------------" << endl;
	cout << gDbFile << endl;
	cout << "-----------------------------------" << endl;

//list<HepRef(GdbPerStation)> createPerStationIn(HepCalibDatabase* pDB);
	list<HepRef(GdbPerStation)> 
		listRefStation = 
			gDbFile.createPerStationIn( pGdbUpdater->getDatabase() );

	cout << "Number of station is " << listRefStation.size() << endl;


	for(	list<HepRef(GdbPerStation)>::iterator
		iRefStation  = listRefStation.begin();
		iRefStation != listRefStation.end();
		iRefStation++){
	
		cout << *(*iRefStation) << endl;

	}

	pGdbUpdater->abort();

	reg_->callEndMethods();
	//      delete pGdbUpdater; // Temporaly solution.

	return 0;

	exit(0);
}
