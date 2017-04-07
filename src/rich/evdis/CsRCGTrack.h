// $Id:

/*!
   \file    CsRCGTrack.h
   \brief   CORAL Event Display Package.
   \version $Revision: 1.3 $
   \author  Take-Aki TOEDA
   \date    $Date: 2010/06/18 10:44:21 $
*/

#ifndef CsRCGTrack_h
#define CsRCGTrack_h
#include <iostream>

#include "coral_config.h"
#include "TCanvas.h" //Used kButton1Down
#include "TPolyMarker3D.h"
#include "TList.h"

#include "CsTypes.h"

/*! \class CsRCGTrack
    \brief   CORAL Event Display Package.
    Track class in graphics
*/

class CsRCGTrack : public TPolyMarker3D
{
 public:
  CsRCGTrack();
  CsRCGTrack(int32 n,Float_t *x,Marker_t m,Int_t color, Float_t size);
  virtual ~CsRCGTrack();
  virtual void ExecuteEvent(int32 event, int32 px, int32 py);
  char *GetParticleName(){return Name;}
  int GetCharge(){return Charge;}
  void SetCharge(int);
  void SetName(const char*);
  void GetInfo();
 private:
  Char_t Name[80];
  int32 Charge;
  float factorX;
  int32 Origin;
  int32 Nfirst;

  TList *listTrack;
  ClassDef(CsRCGTrack,1)
};
#endif
