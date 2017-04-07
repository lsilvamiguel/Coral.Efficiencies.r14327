// $Id:

/*!
   \file    CsRCIFtoDisplay.h
   \brief   CORAL Event Display Package.
   \version $Revision: 1.2 $
   \author  Take-Aki TOEDA
   \date    $Date: 2010/06/18 10:44:21 $
*/

#ifndef CsRCIFtoDisplay_h
#define CsRCIFtoDisplay_h

#include "coral_config.h"

#include "TRotMatrix.h"
#include "TMatrix.h"
#ifdef Check
#undef Check
#endif

#include "CsGeant3.h"
#include "CsEvent.h"
#include "CsRichOneDisplay.h"
//#include "CsDecoder.h"
#include "CsDigit.h"
#include "CsRCDisplay.h"
#include "CsRCGHit.h"
#include "CsRCGRing.h"
#include "CsRCGTrack.h"

TRotMatrix *ConvertMatrix(CLHEP::HepMatrix);

/*! \class CsRCIFtoDisplay
    \brief   CORAL Event Display Package.
    Interface class connecting CORAL main part to graphics part
*/

class CsRCIFtoDisplay{
 public:
  CsRCIFtoDisplay();
  ~CsRCIFtoDisplay();
  void Trans(CsRichOneDisplay*);

 private:
  CsRCGHit    *hit;
  CsRCGRing   *ring;
  CsRCGTrack  *track;
  TRotMatrix *RotMat;
  CLHEP::HepMatrix HepRotMat;
};
#endif
