#include "CsMagFieldSol.h"

CsMagFieldSol::CsMagFieldSol(){
}

void CsMagFieldSol::addGrid(CsFieldGridSol grid){
  _list.push_back(grid);	
}

int CsMagFieldSol::getNoGrid(){
  return _list.size();
}

