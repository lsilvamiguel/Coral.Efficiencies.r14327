/*!
   \file    CsRCGraph.cc
   \brief   Compass Graphics Class.
   \version 0.2
   \author  Take-Aki TOEDA
   \author  Alexander Zvyagin
   \date    30 April 1999
*/

#include "CsRCGraph.h"

CsRCGraph* CsRCGraph::instance = 0;

CsRCGraph::CsRCGraph(VoidFuncPtr_t initfuncs[])
{}

CsRCGraph* CsRCGraph::Init(VoidFuncPtr_t initfuncs[])
{
  if(instance == 0){ 
    instance = new CsRCGraph(initfuncs);
  }
  return instance;
}

void CsRCGraph::Close()
{}
