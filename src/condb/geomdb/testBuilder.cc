#include "CsSTD.h"
#include "CsInit.h"
#include "GdbDBuilder.h"

int main(int argc, char *argv[]){
	CsInit* init = CsInit::Instance(argc, argv);

	cout	<< "\n\n\n\n\n\n\n\n\n\n================================================\n"
			<< "====== GDB Detector Builder:(Initialization)" << endl;
	GdbDBuilder builder;


	cout	<< "\n\n\n\n\n\n\n\n\n\n================================================\n"
			<< "===== GDB Detector Builder:(Building)" << endl;
	builder.build();



	cout	<< "\n\n\n\n\n\n\n\n\n\n================================================\n"
			<< "===== GDB Detector Builder:(Building SUccessfully done." << endl;

	cout << builder << endl;

	CsRegistrySing::Instance()->callEndMethods();
	exit(0);
}
