#include "CsPersFieldSol.h"

CsPersFieldSol::CsPersFieldSol(){
  _noGrid=0;
}

void CsPersFieldSol::addGrid(CsFieldGridSol grid){
  _noGrid++;
  _array[_noGrid-1] = grid;
}

void CsPersFieldSol::changeGridNumber(int n){
  _array.resize(n);
}

CsFieldGridSol CsPersFieldSol::getGrid(int index){
  return _array[index];
}











