#include "CsSTD.h"
#include "CsInit.h"
#include "CsOpt.h"
#include "GdbDataFile.h"
#include "GdbPerStation.h"

int main(int argc, char *argv[]){
	CsInit* init = CsInit::Instance(argc, argv);

	string tag = "GeomDB";
	string key = "datafile";
	string datafile;

	
	if( ! CsOpt::Instance()->getOpt( tag, key, datafile ) ) {
		cerr << " No file is defined."  << datafile << endl;
		exit(1);
	}

	GdbDataFile gDbFile(datafile);

	cout << gDbFile << endl;

	exit(0);
}
