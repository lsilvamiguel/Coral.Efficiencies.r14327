/*!
   \file    CsRCGraph.h
   \brief   Compass Graphics Class.
   \version 0.1
   \author  Take-Aki TOEDA
   \date    21 May 1999
*/

#ifndef CsRCGraph_h
#define CsRCGraph_h

#include "coral_config.h"
#include "CsSTD.h"

#include "Rtypes.h"

/*!
  \class CsRCGraph 
  \brief   Compass Graphics Class.
*/
class CsRCGraph{
  protected:
  CsRCGraph(VoidFuncPtr_t initfuncs[]); 
 public:
  static CsRCGraph *Init(VoidFuncPtr_t initfuncs[]); 
  void Close();

 private:
  static CsRCGraph* instance;
};
#endif
