// $Id:

/*!
   \file    CsRCGRing.h
   \brief   CORAL Event Display Package.
   \version $Revision: 1.2 $
   \author  Take-Aki TOEDA
   \date    $Date: 2010/06/18 10:44:21 $
*/

#ifndef CsRCGRing_h
#define CsRCGRing_h
#include <iostream>

#include "coral_config.h"
#include "TCanvas.h"
#include "TPolyLine3D.h"

#include "CsTypes.h"

class CsRCGRing : public TPolyLine3D
{
 public:
  CsRCGRing( int32 n,float64 *x, int32 color, int32 type);
  virtual ~CsRCGRing();
  int GetType();
  void SetType(int);
  void GetInfo();
  virtual void ExecuteEvent(int32 event, int32 px, int32 py);
 private:
  int32 Type;

  ClassDef(CsRCGRing,1)
};
#endif
