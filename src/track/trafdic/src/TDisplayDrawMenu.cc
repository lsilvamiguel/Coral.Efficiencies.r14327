// $Id: TDisplayDrawMenu.cc,v 1.2 2010/02/03 18:22:23 suhl Exp $

/*!
  Draw menu
  Return user's choise
*/

#include <iostream>
#include <stdio.h>
#include "TEv.h"
#include "TOpt.h"
#include "TSetup.h"
#include "TDisplay.h"
#include "higz.h"

/*
  Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/TDisplayDrawMenu":
   i) Add "Write back" button.
*/

int TDisplay::DrawMenu()

{

  const TEv& ev = TEv::Ref();
  static int iflg=1; // toggle draw/hide menu
  int ich=0,im=0;
  int nbu=0;
  const int len=17;
  const int n = 19;
  float x1,x2,y1,y2;

  char* chu[nbu];

  // Memu items text
  char  chd[n][len];
  char  chv[n][len];
  char  chi[n][len]=
  {
    "     Zoom in    ",
    "    Undo zoom   ",
    "     Zoom out   ",
    "      Move      ",
    "      Locate    ",
    "      Ruler     ",
    "   Track Info   ",
    "   Plane Info   ",
    "   Hit   Info   ",
    "                ", // chi[ 9]
    "                ", // chi[10]
    "  Reset view    ",
    "    Open PS     ", // chi[12]
    "    Options     ",
    "   Next event   ",
    "                ",
    "-               ", // chi[16]
    "   Write Back   ",
    "      Exit      "
  };


  // "context dependent" options

  sprintf(chi[16],"- Proj. %2u/%2zu ", Proj, TSetup::Ref().vProj().size()-1);

  if(DrOpt[2]==-1) strcpy(chi[12],"   Save to PS   ");
  if(DrOpt[2]== 1) strcpy(chi[12],"    Close PS    ");

  if(ev.IsMC() && DrOpt[0] > 0) strcpy(chi[9],    "   MC Hit Info  ");
  if(ev.IsMC() && DrOpt[9] > 0) strcpy(chi[10],   "   MC Trk Info  ");


  for(int i=0; i < n; i++) 
    strcpy(chd[i],"                ");

  ISELNT(0);  // Default Normalization Transformation
  
  if(Xmenu[0] < 0.5){
    x1=Xmenu[0];
    x2=x1+len*Wmenu[0];
  } else {
    x2=Xmenu[0];
    x1=x2-len*Wmenu[0];
  }
  y1 = Xmenu[1];
  y2 = y1 - n*Wmenu[1];

  do{
    if(iflg > 0){
      char s1[] = "Compass", s2[] = "RSTCHWD";
      IGMENU(im,s1,x1,x2,y1,y2,nbu,chu,n,chi,chd,chv,ich,s2);
    }
    else{
      char s1[] = "Compass", s2[] = "TC";
      IGMENU(im,s1,x1,x2,y1,y2,nbu,chu,1,chi,chd,chv,ich,s2);
    }
    if(ich==-1000) iflg*=-1;
  } while(ich == 0);

  IUWK(0,1);
  ISELNT(10); // NT defined by view field size
  return(ich);
}




















