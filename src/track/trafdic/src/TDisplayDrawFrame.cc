/*!
  Draw axis
*/

#include "TOpt.h"
#include "TDisplay.h"
#include "higz.h"

void TDisplay::DrawFrame()
{

  ISELNT(0);

  float xl0=XvfN[0];
  float xh0=XvfN[1];
  float yl0=YvfN[0];
  float yh0=YvfN[1];
  
  float xl=Xvf[0];
  float xh=Xvf[1];
  float yl=Yvf[0];
  float yh=Yvf[1];
  
  ISTXCI(1);

  // Draw scale

  IGSET("LAOF",1.5*(xh0-xl0)/100.);
  IGSET("LASI",0.8*(xh0-xl0)/100.);
  IGSET("TMSI",1.0*(xh0-xl0)/100.);
  IGAXIS(xl0,xh0,yl0,yl0,xl,xh,505,"DHS=+" );  // bottom
  IGAXIS(xl0,xh0,yh0,yh0,xl,xh,505,"DHSU-");   // top
  IGAXIS(xl0,xl0,yl0,yh0,yl,yh,505,"DHSU-");   // left
  IGAXIS(xh0,xh0,yl0,yh0,yl,yh,505,"DHS=+" );  // right

}











