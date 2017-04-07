// $Id: TDisplayDrawMaterialMap.cc 13148 2011-12-28 16:55:25Z kbicker $

/*!
  Draw material map
*/

/*
  Version specific to "lattice" alternative.
  Changes w.r.t. "traffic/TDisplayDrawMaterialMap.cc":
  i) The vertical slice through the map is done on the beam axis.
*/

#include <iostream>
#include <math.h>
#include "TDisplay.h"
#include "TOpt.h"
#include "TSetup.h"
#include "THlx.h"
#include "TConstants.h"
#include "CsGeom.h"
#include "higz.h"

using namespace std;

void TDisplay::DrawMaterialMap()
{

  const TSetup& setup = TSetup::Ref();

  if(DrOpt[12] < 0) return;
  if(DrOpt[12] > 0 && TOpt::ReMode[20] == 0) {
    cout<<"TDisplayDrawMaterialMap ==> Usage of material map is switched off by ReMode[20] option"<<endl;
    return;
  }

  CsMaterialMap *map = CsGeom::Instance()->getCsMaterialMap();

  if (Proj==0 || Proj==1) { // ***** RESTRICT DRAWING to 0 and 90 PROJ. *****
    float par[10];
    par[0]=50;    // number of countour levels
    par[1]=2.107; // line type and color index
    par[2]=Xvf[0];
    par[3]=Xvf[1];
    par[4]=Yvf[0];
    par[5]=Yvf[1];
    par[6]=0.;
    par[7]=0.;
    par[8]=0.;
    par[9]=0.;
    
    int ncx = 300; // number of cells for contour plot alog X
    int ncy = int(ncx*TConstants_A4);

    float dx=fabs(par[2]-par[3])/ncx;
    float dy=fabs(par[4]-par[5])/ncy;
    
    float v[ncy][ncx];
    
    THlx H, Htmp; H(0)=H(1)=H(2)=H(3)=H(4) = 0;
    H(5) = TOpt::dCut[4] ? TOpt::iCut[15]/TOpt::dCut[4] : 0;
    float step(0),rl(0);

    // fill the table
    float xx=par[2]+dx/2;
    for(int ii=0; ii < ncx; ii++){
      Htmp(0) = xx; H.Extrapolate(Htmp); H(0) = xx;
      float yy=par[4]+dy/2;
      for(int jj=0; jj < ncy; jj++){
	if (Proj==0) {  // HORIZONTAL PROJ.: SLICE THROUGH MIDPLANE
	  H(1) = yy;      H(2) = 0;  map->getRadLength(H,true,rl,step);
	}
	if (Proj==1) {  // VERTICAL   PROJ.: SLICE ON THE BEAM AXIS...
	  // ...so taht the outline of the central hole be drawn on the display
	  H(1) = Htmp(1); H(2) = yy; map->getRadLength(H,true,rl,step);
	}
	if(rl > 29800) rl = 30420;        // almost an air :-)
	v[jj][ii]=10*log(30420/rl);       // log scale (something like dB)
	yy+=dy;
      }
      xx += dx; H = Htmp;
    }
    // Draw the table
    //IGTABL(ncx,ncy,(float*)&v[0][0],2,par,"C");

    par[0]=1;
    par[1]=10.; 
    ISPMCI(109);
    IGTABL(ncx,ncy,(const float**)v,2,par,"P");
    ISPMCI(1);

  } else {
  }
  IGSET("BORD",0);
  IUWK(0,1);
 
}
