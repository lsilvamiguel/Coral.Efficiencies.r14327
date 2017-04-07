#include <CsFieldGridSM1.h> 
CsFieldGridSM1::CsFieldGridSM1(){};

CsFieldGridSM1::CsFieldGridSM1(float x,float y,float z,
				 float Bx,float By,float Bz){
  _x  =  x; _y  =  y; _z  = z;
  _Bx = Bx; _By = By; _Bz = Bz;
};
