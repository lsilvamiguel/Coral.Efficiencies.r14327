
/*!
  To draw text strings in the event display
*/

#include <iostream>
#include <stdio.h>
#include <math.h>
#include "TDisplay.h"
#include "TOpt.h"
#include "TSetup.h"
#include "TEv.h"
#include "CsEvent.h"
#include "TConstants.h"
#include "higz.h"



void TDisplay::DrawString(std::string s, float x, float y, int icol, float slope)
{
  if(DrOpt[14] < 0) return;

  float xl0=XvfN[0];
  float xh0=XvfN[1];
  float yl0=YvfN[0];
  float yh0=YvfN[1];
  
  float xl=Xvf[0];
  float xh=Xvf[1];
  float yl=Yvf[0];
  float yh=Yvf[1];

  if(x < xl || x > xh || y < yl || y > yh) return;

  // recalclulate x,y to NDC
  float xx,yy;

  float ax=(xh0-xl0)/(xh-xl);
  float ay=(yh0-yl0)/(yh-yl);

  xx = xl0 + ax*(x-xl);
  yy = yl0 + ay*(y-yl);


  ISELNT(0);

  ISTXCI(icol);
  //IGSET("CHHE",0.005);
  //ITX(xx, yy, s.c_str());
  float beta;
  if(slope >= 1.e10) {
    beta = M_PI/2.;
  } else {
    beta  = atan(ay*slope/ax);
  }
  float alpha = beta * 180./ M_PI;
  //float delta = 0.002;
  //xx -= delta * sin(beta);
  //yy += delta * cos(beta);
  IGTEXT(xx,yy,s.c_str(),0.007,alpha,"L");

  ISELNT(10);

}








