#include "CsPersFieldSM1.h"

CsPersFieldSM1::CsPersFieldSM1(){
  _noGrid=0;
}

void CsPersFieldSM1::addGrid(CsFieldGridSM1 grid){
  _noGrid++;
  _array[_noGrid-1] = grid;
}

void CsPersFieldSM1::changeGridNumber(int n){
  _array.resize(n);
}

CsFieldGridSM1 CsPersFieldSM1::getGrid(int index){
  return _array[index];
}












