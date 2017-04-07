// $Id: TDisplayDrawMag.cc 13123 2011-12-20 09:09:57Z kbicker $

/*!
  Draw magnets and mag. field
*/

#include <iostream>
#include "TAlgo.h"
#include "TDisplay.h"
#include "TSetup.h"
#include "TOpt.h"
#include "TConstants.h"
#include "higz.h"

/*
  Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/TDisplayDrawMag":
   i) Bug fix: Drawing SM2 gap as viewed from the side when in the 90 degree
   view.
*/

void TDisplay::DrawMag()
{

  const TSetup& setup = TSetup::Ref();

  float x[5],y[5];
  float x0;
  int im(0);

  if(setup.NMags <  2) return;
  if(setup.NMags == 2) im = 0;
  if(setup.NMags == 3) im = 1;

  if(DrOpt[11] == -1){
    ISPLCI(105);
  }
  if(DrOpt[11] ==  1){
    ISPLCI(105);
    IGSET("BORD",1);
  }

  x0 = setup.MagCenter(im,0);

  //         ***** DRAW MAG FIELD AS CONTOUR LEVELS... *****
  if (DrOpt[8]==1 &&
      TOpt::Graph[3]>0 &&       // ...depending upon Graph[3], which defaults to 3, cf. "TDisplay::DrawModeMenu"
      (Proj==0 || Proj==1)) {   // ...when in the top or side view.
    float par[10];
    par[0]=50;    // Number of contour levels
    par[1]=2.107;
    par[2]=Xvf[0];
    par[3]=Xvf[1];
    par[4]=Yvf[0];
    par[5]=Yvf[1];
    par[6]=0.;
    par[7]=0.;
    par[8]=0.;
    par[9]=0.;
    
    int ncx = 300; // Number of cells for contour plot alog X
    int ncy = int(ncx*TConstants_A4);

    float dx=fabs(par[2]-par[3])/ncx;
    float dy=fabs(par[4]-par[5])/ncy;
    
    float v[ncy][ncx];
    
    double xyz[3]; double fyz[3];
    double f(0),f0(0),f1(0),f2(0);
    
    float xx=par[2]+dx/2;
    for(int ii=0; ii < ncx; ii++){
      float yy=par[4]+dy/2;
      for(int jj=0; jj < ncy; jj++){
	xyz[0]=xx;
	if(Proj == 0) {
	  xyz[1]=yy; xyz[2]=0.; f0=TAlgo::Field(xyz,fyz); f1=fyz[2]; f2=sqrt(fyz[1]*fyz[1] + fyz[2]*fyz[2]);
	}
	if(Proj == 1) {
	  xyz[1]=0.; xyz[2]=yy; f0=TAlgo::Field(xyz,fyz); f1=fyz[1]; f2=sqrt(fyz[1]*fyz[1] + fyz[2]*fyz[2]);
	}
	if (TOpt::Graph[3]==1) f = f1;      // ***** THREE Graph[3] OPTIONS
	if (TOpt::Graph[3]==2) f = f2;
	if (TOpt::Graph[3]==3) f = f0;
	v[jj][ii]=float(fabs(f));
	yy+=dy;
      }
      xx+=dx;
    }
    IGTABL(ncx,ncy,(const float**)v,2,par,"C");

  }
  else {  //         ***** DRAW MAGNETS (OUTLINE of GAP) *****
    //       SM1
    ISLN(1);
    x[0]=x[4]=x[3]=setup.MagCenter(im,0)-TConstants_SM1Siz[0];
    y[0]=y[4]=y[1]=setup.MagCenter(im,1)-TConstants_SM1Siz[1];
    x[1]=x[2]     =setup.MagCenter(im,0)+TConstants_SM1Siz[0];
    y[2]=y[3]     =setup.MagCenter(im,1)+TConstants_SM1Siz[1];
    //y[0]=y[4]=y[1]=setup.MagCenter(im,2)-TConstants_SM1Siz[2];
    //y[2]=y[3]     =setup.MagCenter(im,2)+TConstants_SM1Siz[2];
    if(DrOpt[11] == 1){
      ISFAIS(1);
      ISFACI(104);
      IGRAPH(5,x,y,"F");
      ISFAIS(0);
    } else {
      ISPLCI(105);
      IPLm(5,x,y);
    }
    
    //       SM2
    im++;
    x[0]=x[4]=x[3]=setup.MagCenter(im,0)-TConstants_SM2Siz[0];
    x[1]=x[2]     =setup.MagCenter(im,0)+TConstants_SM2Siz[0];
    if (Proj==1) {
      y[0]=y[4]=y[1]=setup.MagCenter(im,2)-TConstants_SM2Siz[2];
      y[2]=y[3]     =setup.MagCenter(im,2)+TConstants_SM2Siz[2];
    }
    else { 
      y[0]=y[4]=y[1]=setup.MagCenter(im,1)-TConstants_SM2Siz[1];
      y[2]=y[3]     =setup.MagCenter(im,1)+TConstants_SM2Siz[1];
    }
    if(DrOpt[11] == 1){
      ISFAIS(1);
      ISFACI(104);
      IGRAPH(5,x,y,"F");
      ISFAIS(0);
    } else {
      ISPLCI(105);
      IPLm(5,x,y);
    }
  }
  ISPLCI(1);
  IGSET("BORD",0);
  IUWK(0,1);
 
}







