// $Id:

/*!
   \file    CsRCGCathode.h
   \brief   CORAL Event Display Package.
   \version $Revision: 1.2 $
   \author  Take-Aki TOEDA
   \date    $Date: 2010/06/18 10:44:21 $
*/

#ifndef CsRCGCathode_h
#define CsRCGCathode_h
#include <iostream>

#include "coral_config.h"
#include "TCanvas.h"
#include "TMarker3DBox.h"

#include "CsTypes.h"

/*! \class CsRCGCathode
    \brief   CORAL Event Display Package.
    Hit class in graphics 
*/

class CsRCGCathode : public TMarker3DBox
{
 public:
  CsRCGCathode( int32 type,float64 x, float64 y, float64 z, float64 dx, float64 dy, float64 dz, int32 co);
  virtual ~CsRCGCathode();
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

  ClassDef(CsRCGCathode,1)
};
#endif
