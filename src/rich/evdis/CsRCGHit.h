// $Id:

/*!
   \file    CsRCGHit.h
   \brief   CORAL Event Display Package.
   \version $Revision: 1.2 $
   \author  Take-Aki TOEDA
   \date    $Date: 2010/06/18 10:44:21 $
*/

#ifndef CsRCGHit_h
#define CsRCGHit_h
#include <iostream>

#include "coral_config.h"
#include "TCanvas.h"
#include "TMarker3DBox.h"

#include "CsTypes.h"

/*! \class CsRCGHit
    \brief   CORAL Event Display Package.
    Hit class in graphics 
*/

class CsRCGHit : public TMarker3DBox
{
 public:
  CsRCGHit( int32 type,float64 x, float64 y, float64 z, float64 dx, float64 dy, float64 dz,int32 co);
  virtual ~CsRCGHit();
  int GetType();
  void SetType(int);
  void SetFactorX(float);
  void GetInfo();
  virtual void ExecuteEvent(int32 event, int32 px, int32 py);
 private:
  int32 Type;
  float32 factorX;
  float32 Xcm,Ycm,Zcm;
  float32 DXcm,DYcm,DZcm;

  ClassDef(CsRCGHit,1)
};
#endif
