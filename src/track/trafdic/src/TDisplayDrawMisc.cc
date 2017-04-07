// $Id: TDisplayDrawMisc.cc 14069 2015-09-17 20:44:46Z lsilva $

/*!
  Draw RICH, target, muon walls and other "decorations"

  - Get target sub-structure from CsGeom.)
  - Get muWall #2 from TOpt arrays "MuonWall" and "Calo".
  - Derive RICH from "TOpt::RICHPipe"
  - The rest is built-in.
*/

#include <iostream>
#include "higz.h"
#include "CsGeom.h"
#include "TDisplay.h"
#include "TSetup.h"
#include "TOpt.h"


using namespace std;

void TDisplay::DrawMisc()
{

  float x[5],y[5], x0; int im(-1);

  const TSetup &setup = TSetup::Ref();

  if      (DrOpt[11]==-1) ISPLCI(105);
  else if (DrOpt[11]== 1) {
    ISPLCI(105);
    IGSET("BORD",1);
  }

  if (TOpt::Graph[4]==0) {
    //               ***** POLARIZED TARGET SETUP: DRAW TARGET MAGNET *****
    float xt = setup.MagCenter(0,0); 
    x[0]=x[4]=x[3]= 100.0+xt;
    y[0]=y[4]=y[1]=-20.;
    x[1]=x[2]     =-100.0+xt;
    y[2]=y[3]     = 20.;
    if(DrOpt[11] == 1) {
      ISFAIS(1);ISFACI(103);
      IGRAPH(5,x,y,"F");
      ISFAIS(0);
    } else {
      ISPLCI(103);
      IPLm(5,x,y);
    }
  }

  if (DrOpt[12]<0) { //     ********** TARGET SUB-STRUCTURE **********
    const vector<CsTargetCell> cells = CsGeom::Instance()->getTargetCells();
    if (cells.size()) { // ***** IF AVAILABLE: GET SUBSTRUCTURE FROM CsGeom...
      float co = cos(AngProj*3.14159265358979323846/180.);
      float si = sin(AngProj*3.14159265358979323846/180.);
      for (int iCell = 0; iCell<(int)cells.size(); iCell++) {
	const CsTargetCell &cell = cells[iCell];
	float yc = co*cell.x[1]+si*cell.x[2];
	x[0]=x[4]=x[3] = (cell.x[0]+cell.length/2)/10;
	y[0]=y[4]=y[1] = (yc-cell.radius)/10;
	x[1]=x[2]      = (cell.x[0]-cell.length/2)/10;
	y[2]=y[3]      = (yc+cell.radius)/10;
	if (DrOpt[11]==1) {
	  ISFAIS(1); ISFACI(106); IGRAPH(5,x,y,"F"); ISFAIS(0);
	}
	else {
	  ISPLCI(106); IPLm(5,x,y);
	}
      }
    }
    else if (TOpt::Graph[4]==0) {      //  ***** ...ELSE: DEFAULT MUON TARGET...
      float xt = setup.TargetCenter[0]; 
      x[0]=x[4]=x[3]= -5+xt;  // Draw pol. targ. upstream part
      y[0]=y[4]=y[1]=-1.5;
      x[1]=x[2]     = -65+xt;
      y[2]=y[3]     = 1.5;
      if (DrOpt[11]==1) {
	ISFAIS(1);ISFACI(106); IGRAPH(5,x,y,"F"); ISFAIS(0);
      }
      else {
	ISPLCI(106); IPLm(5,x,y);
      }
      x[0]=x[4]=x[3]=65.+xt;   // Draw pol. targ. downstream part
      y[0]=y[4]=y[1]=-1.5;
      x[1]=x[2]     = 5+xt;
      y[2]=y[3]     = 1.5;
      if(DrOpt[11] == 1) {
	ISFAIS(1);ISFACI(106); IGRAPH(5,x,y,"F"); ISFAIS(0);
      }
      else {
	ISPLCI(106); IPLm(5,x,y);
      }
    }
    else {                                //  ***** ...OR: DEFAULT HADRON TARGET
      float xt = setup.TargetCenter[0]; 
      x[0]=x[4]=x[3]=   0.15+xt;
      y[0]=y[4]=y[1]=  -3.0;
      x[1]=x[2]     =  -0.15+xt;
      y[2]=y[3]     =   3.0;
      if(DrOpt[11] == 1) {
	ISFAIS(1);ISFACI(6); IGRAPH(5,x,y,"F"); ISFAIS(0);
      }
      else {
	ISPLCI(106); IPLm(5,x,y);
      }
    }
  }
  
  //       ***** DRAW MUOMN WALL #2 *****
  if (DrOpt[12]<0 && Proj<=1) {  // IF Y or Z PROJ. and NO MATERIAL MAP
    double *data[] =  // Sucessively MF2, HCAL2, ECAL2
      { TOpt::MuonWall+10, TOpt::Calo+10, TOpt::Calo+30 };
    for (int iMat = 0; iMat<3; iMat++) {
      float xMat = *data[iMat], dxMat = *(data[iMat]+3);
      float yMat, dyMat, yHole, dyHole;
      if (Proj==0) {     // HORIZONTAL 
	yMat  = *(data[iMat]+1), dyMat = *(data[iMat]+4);
	yHole = *(data[iMat]+6), dyHole =*(data[iMat]+8);
      }
      else {             // VERTICAL
	yMat  = *(data[iMat]+2), dyMat = *(data[iMat]+5);
	yHole = *(data[iMat]+7), dyHole =*(data[iMat]+9);
      }
      for (int lowUp = -1; lowUp<=1; lowUp += 2) {
	x[0]=x[4]=x[3] = xMat-dxMat;
	y[0]=y[4]=y[1] = yMat+lowUp*dyMat;
	x[1]=x[2]      = x[0]+2*dxMat;
	if (iMat==0 && Proj==0) { // Central hole of muWall MF2 in horizontal...
	  // ...its shape is complex => Use built-in description
	  yHole = 36  ; dyHole = 21; y[2] = yHole+lowUp*dyHole;
	  yHole = 31.5; dyHole = 17; y[3] = yHole+lowUp*dyHole;
	}
	else {
	  y[2]=y[3] = yHole+lowUp*dyHole;
	}
	if (DrOpt[11]==1) {
	  ISFAIS(1); ISFACI(109);
	  IGRAPH(5,x,y,"F");
	  ISFAIS(0);
	}
	else {
	  ISPLCI(109); IPLm(5,x,y);
	}
      }
    }
  }

  //       ***** DRAW HADRON ABSORBER *****
  // - Note that the position and dimension (if not the shape) are taken
  //  from ROOTGeometry
  if (DrOpt[12]<0 && hasAbsorber) {
    if (DrOpt[11]==1) {
      ISFAIS(1); ISFACI(101);
    }
    if      (Proj==0) {    
      if (DrOpt[11] == 1) IGRAPH(13,xA,yA,"F");
      else { ISPLCI(102); IPLm(13,xA,yA); }
    }
    else if (Proj==1) {
      if (DrOpt[11] == 1) IGRAPH(13,xA,zA,"F");
      else { ISPLCI(102); IPLm(13,xA,zA); }
    }
  }

  if (setup.NMags<2) return;
  if (setup.NMags==2) im = 0;
  if (setup.NMags==3) im = 1;
  x0 = setup.MagCenter(im,0);

  //       ***** DRAW RICH *****
  if (DrOpt[11]==1) {
    ISFAIS(1); ISFACI(102);
  }
  x0 += 230.0;
  float dx = 150.; x0 += dx;
  if (TOpt::RICHPipe[0]) {
    x0 = TOpt::RICHPipe[0]; dx = TOpt::RICHPipe[3];
  }
  float xR[5] = { x0+dx,  x0-dx, x0-dx, x0+dx,  x0+dx };
  if      (Proj==0) {
    for (int du = -1; du<=1; du += 2) {
      float yp = TOpt::RICHPipe[1]+du*TOpt::RICHPipe[4];
      float y[5] = {float(du*220.), float(du*160.), yp, yp, float(du*220.) };
      if (DrOpt[11] == 1) IGRAPH(5,xR,y,"F");
      else { ISPLCI(102); IPLm(5,xR,y); }
    }
  }
  else if (Proj==1) {
    for (int du = -1; du<=1; du += 2) {
      float zp = TOpt::RICHPipe[2]+du*TOpt::RICHPipe[4];
      float y[5] = {float(du*220.), float(du*160.), zp, zp, float(du*220.) };
      if (DrOpt[11] == 1) IGRAPH(5,xR,y,"F");
      else { ISPLCI(102); IPLm(5,xR,y); }
    }
  }
  else {
    ISLN(2);
    float y[5] = {-280., -160., 160., 280., -280. };
    ISPLCI(102); IPLm(5,xR,y);
  }

  ISFAIS(0);
  ISLN(1);
  ISPLCI(1);
  IGSET("BORD",0);
  IUWK(0,1);
}








