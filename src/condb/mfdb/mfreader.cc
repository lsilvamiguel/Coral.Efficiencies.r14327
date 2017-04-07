#include "CsOpt.h"
#include "CsGeom.h"
#include "CsRegistrySing.h"
#include "CsSTD.h"

main(int argc, char *argv[]){

  CsOpt* coralOpt = CsOpt::Instance(argc,argv);
  CsRegistrySing *reg_ = CsRegistrySing::Instance();
  
  string tag,key,str,detTable_;

  ofstream f("out.dat");
  CsGeom *geom = CsGeom::Instance(new const string("detectors.dat"));
  CsField *field = geom->getCsField();
  float x,y,z,Bx,By,Bz;
  for(int i=0;i<1000;i+=10){
    for(int j=0;j<1000;j+=10){
      for(int k=0;k<10000;k+=100){
    	x=i,y=j,z=k;
	field->getField(x,y,z,Bx,By,Bz);
	f << " x = " << x 
	  << " y = " << y
	  << " z = " << z
	  << " Bx = " << Bx
	  << " By = " << By
	  << " Bz = " << Bz
	  << endl;	
      }
    }
  }
  reg_->callEndMethods();
}













