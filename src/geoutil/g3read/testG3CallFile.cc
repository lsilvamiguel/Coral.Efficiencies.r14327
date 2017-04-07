#if __GNUG__ >= 2
#  pragma implementation
#endif

#include "CsInit.h"
#include "CsG3CallFile.h"

int main(int argc, char *argv[]){
	CsInit* init = CsInit::Instance(argc, argv);
	CsRegistrySing *reg_ = CsRegistrySing::Instance();


	std::cout << "Test program for G3call List file reader." << std::endl;


	CsG3CallFile* pG3callList = CsG3CallFile::Instance();
	
	std::cout << (*pG3callList) << std::endl;


	std::cout << "DUMP DETECTOR STRUCUTUR" << std::endl;

	(*pG3callList).HALL().dumpDetectors( std::cout );
	
	std::cout << std::endl;

	std::cout << "END TEST PROGRAM" << std::endl;


	reg_->callEndMethods();

	exit(0);
}
