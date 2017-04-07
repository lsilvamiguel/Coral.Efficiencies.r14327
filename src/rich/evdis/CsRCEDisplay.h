// $Id:

/*!
   \file    CsRCEDisplay.h
   \brief   CORAL Event Display Package.
   \version $Revision: 1.2 $
   \author  Take-Aki TOEDA
   \date    $Date: 2010/06/18 10:44:21 $
   
   This is a tempolary class.
*/

#ifndef CsRCEDisplay_h
#define CsRCEDisplay_h

#include "coral_config.h"
#include "TPad.h"

#include "CsTypes.h"

/*! \class CsRCEDisplay
    \brief   CORAL Event Display Package.
    Event display window
    This is tempolary class.
*/

class CsRCEDisplay : public TPad
{
 public:
  CsRCEDisplay();
  virtual ~CsRCEDisplay();
  virtual void ExecuteEvent(int32 event, int32 px, int32 py);
  virtual int32 DistancetoPrimitive(int32 px, int32 py); 
  ClassDef(CsRCEDisplay,0)
};
#endif
