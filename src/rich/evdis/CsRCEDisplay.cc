// $Id:

/*!
   \file    CsRCEDisplay.cc
   \brief   CORAL Event Display Package.
   \version $Revision: 1.2 $
   \author  Take-Aki TOEDA
   \date    $Date: 2010/06/18 10:44:21 $
   
   This is a tempolary class.
*/

#include <iostream>

#include "coral_config.h"
#include "CsRCEDisplay.h"

ClassImp(CsRCEDisplay)

CsRCEDisplay::CsRCEDisplay() : TPad("viewpad", "Event Display",0,0,1,1) {
  SetRightMargin(.01);
  SetLeftMargin(.01);
  SetTopMargin(.01);
  SetBottomMargin(.01);
}

CsRCEDisplay::~CsRCEDisplay()
{}
  
void CsRCEDisplay::ExecuteEvent(int32 event, int32 px, int32 py)
{
   static float32 x,y;
   static float32 xmin,ymin;
   static float32 xrange,yrange;
   static float32 pxold, pyold;
}

int32 CsRCEDisplay::DistancetoPrimitive(int32 px, int32 py)
{

   return 0;
}













