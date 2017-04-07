#if __GNUG__ >= 2
#  pragma implementation
#endif

#include "CsSTD.h"
#include "CsInit.h"
#include "CsRegistrySing.h"
#include "CsOpt.h"
#include "BmsContainerConf.h"
int main(int argc, char *argv[]){

	CsInit* init = CsInit::Instance(argc, argv);

	BmsContainerConf cofiguration("template.dat");
	
	if( cofiguration.read() ){
		
		cout << cofiguration << endl;
		
	}


	exit(0);

}
