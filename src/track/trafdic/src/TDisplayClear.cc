
#include "TOpt.h"
#include "TSetup.h"
#include "TDisplay.h"
#include "higz.h"
#include "TConstants.h"

/*!
  Private methof to clear display  
  and set normalization transformations
  Called by TDisplay::Draw()
*/

void TDisplay::Clear()
{

  //Clear
 
  ICLRWK(0,1);
  ISCLIP(1);
  char c[] = "*"; 
  IGSET(c,0);

  float xl0=0.;
  float xh0=1.;
  float yl0=0.;
  float yh0=TConstants_A4;
   
  AngProj = 0.1*TSetup::Ref().vProj()[Proj];
  
  ISELNT(0);
  
  // srink a bit normalized viewfield 
  float r;
  if(DrOpt[2]==-1) r = 0.04;
  else            r = 0.1;  // if draw to PS file
  xl0+=r; xh0-=r;
  float lv = (xh0-xl0) * TConstants_A4;
  yl0=(TConstants_A4-lv)/2.;
  yh0=yl0+lv;

  // store current normalized viewfield
  XvfN[0]=xl0;
  XvfN[1]=xh0;
  YvfN[0]=yl0;
  YvfN[1]=yh0;

  if(DrOpt[11] == 1) {
    ISFACI(101); // color for background (defined in TDisplay::Init()
    ISFAIS(1); 
    IGBOX(xl0,xh0,yl0,yh0); // fill all view field
    ISFAIS(0); 
  }
  
  //
  //      Normalization transformation (#10)
  //

  float xl=Xvf[0], xh=Xvf[1];
  float yl=Yvf[0], yh=Yvf[1];
  
  ISWN(10,xl,xh,yl,yh);
  ISVP(10,xl0,xh0,yl0,yh0);
  ISELNT(10);

}











