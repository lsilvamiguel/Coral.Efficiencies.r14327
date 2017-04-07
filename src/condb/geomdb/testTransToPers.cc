#include "CsSTD.h"
#include "CsInit.h"
#include "GdbDBuilder.h"
#include "GdbPerStation.h"

int main(int argc, char *argv[]){
	CsInit* init = CsInit::Instance(argc, argv);

	CsRegistrySing *reg_ = CsRegistrySing::Instance();

	GdbDBuilder builder;
	builder.build();


	ofstream originalTable("original.dat");
	ofstream copiedTable("copied.dat");

	list<GdbDetector>::iterator iDetList;
	for(	iDetList =  builder.detectors().begin();
		iDetList != builder.detectors().end();
		iDetList++){
		
		vector<GdbStation>::iterator iStaList;
		for(	iStaList =  (*iDetList).station().begin();
			iStaList !=  (*iDetList).station().end();
			iStaList++){
			
			originalTable	<< (*iStaList)	<< endl;
					
			cerr << "#------------- COPY START ----------------" << endl;

			GdbPerStation newStation( (*iStaList) );

			cerr << "#------------- COPY DONE ----------------" << endl;

			copiedTable	<< newStation	<< endl;

		}

	}

	reg_->callEndMethods();

	exit(0);
}
