#include "CsMagFieldSM1.h"

CsMagFieldSM1::CsMagFieldSM1(){
}

void CsMagFieldSM1::addGrid(CsFieldGridSM1 grid){
  _list.push_back(grid);	
}

int CsMagFieldSM1::getNoGrid(){
  return _list.size();
}

