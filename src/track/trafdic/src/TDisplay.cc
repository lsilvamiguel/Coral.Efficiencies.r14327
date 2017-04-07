#include <iostream>

#include "TROOT.h"
#include "TClass.h"
#include "TObjArray.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoBBox.h"
#include "TGeoTube.h"

#include "CsErrLog.h"
#include "TDisplay.h"
#include "TOpt.h"
#include "TAlgo.h"
#include "TConstants.h"
#include "higz.h"

TDisplay* TDisplay::address = 0; // init static pointer

/*! 
  TDisplay constructor
  calls private member function TDisplay::Init().
  This function put all necessary information into TDisplay data members.

*/

TDisplay::TDisplay():
  rkutraj(false),  // flag to draw Runge-Kutta trajectories
  mmaptraj(false)  // flag to draw calls to material map  
{
  NobjCreated++;
  if(address == 0 ){ //!< if not yet exist 
    address = this;
    if( TOpt::Graph[0] > 0 ){
      // Display initialization
      if(!Init()) assert(false);
    }
  } else {
    std::cout<<"Only one instance of the class TDisplay is allowed"<<std::endl;
    assert(false);
  }

  //         ********** RE-ORDER RECONSTRUCTION ZONES **********
  const TSetup &setup = TSetup::Ref();
  int iZone, nZones = setup.vIplLast().size();
  iZone2Zone = new int[nZones]; iplLasts = new int[nZones];
  for (iZone = 0; iZone<nZones; iZone++) {
    iZone2Zone[iZone] = iZone; iplLasts[iZone] = setup.vIplLast()[iZone];
  }
  for (iZone = 0; iZone<nZones; iZone++) {
    int last = iplLasts[iZone], zone = iZone2Zone[iZone];
    for (int iZp = iZone+1; iZp<nZones; iZp++) if (iplLasts[iZp]<last) {
      iplLasts[iZone] = iplLasts[iZp];     iplLasts[iZp] = last;
      last = iplLasts[iZone];
      iZone2Zone[iZone] = iZone2Zone[iZp]; iZone2Zone[iZp] = zone;
      zone = iZone2Zone[iZone];
    }
  }

  //   ********** RETRIEVE DESCRIPTION OF DRELL-YAN ABSORBER  **********
  hasAbsorber = 0; // 0: no absorber, >0: else
  if (gGeoManager) {
    //  Let's start 1st retrieving the sizes from volumes. We can do this by
    // inspecting the list of TGeoVolume's, where one can access individual
    // volumes independent the tree structure they're embedded in. Next we will
    // retrieve the volume positions from the tree structure via TGeoNode's. But
    // in doing so, we will limit ourselves to a few cases, assuming for the
    // rest, that the voumes are centered in their mother volumes or that they
    // otherwise obey the simplest possible positioning scheme.
    // - "ABSV" is the base volume (BOX) of the absorber
    // - "ALB1|2|3" (TUBE) describe the bore hole filled w/ air
    // - "ABA1" gives the size of the core.
    // - "ABCC" describes the concrete, and correlatively, the length of the
    //  core part which is protruding downstream
    const TObjArray *l = gGeoManager->GetListOfVolumes();
    if (l) {
      TObject *o = l->FindObject("ABSV");
      if(o) {
        TGeoVolume *vol = (TGeoVolume*)o;
        TGeoShape *sh = vol->GetShape();
        if(sh->IsA()->InheritsFrom(TGeoBBox::Class())) {
          hasAbsorber += 1;
          const TGeoBBox *ABSV = (TGeoBBox*)sh;
          double dx = ABSV->GetDX(), dy = ABSV->GetDY(), dz = ABSV->GetDZ();
          // Initializing also bore hole and concrete, expecting them to be
          // // downstream of ABSV
          xA[0]=xA[7]=xA[12]=xA[8]=xA[11] = -dx;
          xA[3]=xA[4] = +dx; 
          yA[0]=yA[1]=yA[12] = +dy; yA[6]=yA[7] = -dy;
          zA[0]=zA[1]=zA[12] = +dz; zA[6]=zA[7] = -dz;
          // // Init bore hole assuming it's justified on upstream edge of "ABSV"
          xA[9]=xA[10] = -dx;
        }
      }
      o = l->FindObject("ABCC");
      if(o) {
        TGeoVolume *vol = (TGeoVolume*)o;
        TGeoShape *sh = vol->GetShape();
        if (sh->IsA()->InheritsFrom(TGeoBBox::Class())) {
          hasAbsorber += 1; //Concrete
          const TGeoBBox *ABCC = (TGeoBBox*)sh;
          xA[1] = ABCC->GetDX(); xA[2]=xA[5]=xA[6] = xA[1];
        }
      }
      o = l->FindObject("ABA1");
      if(o) {
        TGeoVolume *vol = (TGeoVolume*)o;
        TGeoShape *sh = vol->GetShape();
        if (sh->IsA()->InheritsFrom(TGeoBBox::Class())) {
          hasAbsorber += 1;  // Core
          const TGeoBBox *ABA1 = (TGeoBBox*)sh;
          double dy = ABA1->GetDY(), dz = ABA1->GetDZ();
          yA[2] += dy; yA[5] -= dy; yA[3] += dy; yA[4] -= dy;
          zA[2] += dz; zA[5] -= dz; zA[3] += dz; zA[4] -= dz;
        }
      }
      for(int i = 1; i <= 3; ++i) {
        o = l->FindObject(Form("ALB%d",i));
        if(o) {
          TGeoVolume *vol = (TGeoVolume*)o;
          TGeoShape *sh = vol->GetShape();
          if (sh->IsA()->InheritsFrom(TGeoTube::Class())) {
            hasAbsorber += 1;
            const TGeoTube *ALBi = (TGeoTube*)sh;
            double dx = ALBi->GetDz(), R = ALBi->GetRmax();
            xA[9] += 2*dx; xA[10] += 2*dx;
            yA[8]=zA[8]=yA[9]=zA[9] = -R; yA[11]=zA[11]=yA[10]=zA[10] = R;
          }
        }
      }
    }

    // Now the search for nodes "ABSV_1" and "ABCC_1" via "DY14_1"
    l = gGeoManager->GetListOfNodes(); if (l) {
      const TObject *o = l->First(); if (o && !strcmp(o->GetName(),"HALL_1")) {
	TGeoNode *hall = (TGeoNode*)o, *filia, *absv = 0, *dy14 = 0, *abcc = 0;
	for (int d = 0; d<hall->GetNdaughters(); d++) {
	  filia = hall->GetDaughter(d);
	  if (!strcmp(filia->GetName(),"ABSV_1")) { absv = filia; break; }
	}
	if (absv) {
	  hasAbsorber += 1;
	  double local[3] = {0,0,0}, master[3];
	  absv->LocalToMaster(local,master);
	  for (int i = 0; i<13; i++) {
	    xA[i] += master[0]; yA[i] += master[1]; zA[i] += master[2];
	  }
	  for (int d = 0; d<absv->GetNdaughters(); d++) {
	    filia = absv->GetDaughter(d);
	    if (!strcmp(filia->GetName(),"DY14_1")) { dy14 = filia; break; }
	  }
	}
	if (dy14) {
	  // Assuming "DY14_1" is centered on "ABSV_1"
	  for (int d = 0; d<dy14->GetNdaughters(); d++) {
	    filia = dy14->GetDaughter(d);
	    if (!strcmp(filia->GetName(),"ABCC_1")) {
	      hasAbsorber += 1;
	      double local[3] = {0,0,0}, master[3];
	      filia->LocalToMaster(local,master);
	      xA[1] += master[0]; xA[2]=xA[5]=xA[6] = xA[1];
	    }
	  }
	}
      }
      //#define TDisplay_DEBUG_ABSORBER
#ifdef TDisplay_DEBUG_ABSORBER
      printf("%d\nxA =",hasAbsorber);
      for (int i = 0; i<13; i++) printf(" %7.2f",xA[i]);
      printf("\nyA =");
      for (int i = 0; i<13; i++) printf(" %7.2f",yA[i]);
      printf("\nzA =");
      for (int i = 0; i<13; i++) printf(" %7.2f",zA[i]);
      printf("\n");
#endif
    }
    if (hasAbsorber!=8) {
      CsErrLog::mes(elError,"Could not get a complete description of the hadron abosrber\n"
		    " => Will not be able to draw it in the event display!");
      hasAbsorber = 0;
    }
  }
  
}

//! Destructor
TDisplay::~TDisplay()
{
  NobjDestructed++;
  // Close graphics
  if(TOpt::Graph[0] > 0) {
    if(DrOpt[2] == 1) CLOPS(); // close PS file if it was opened
    IGEND();
  }
    address = 0;
}

void TDisplay::point(double xx[], int color)
{
  //
  // Draw dot
  //
  if(TOpt::Graph[0] == 0) return;
  
  float x(0),y(0);
  ISMK(20);ISPMCI(color);
  ISMKSC(0.2);
  x = float(xx[0]);
  float co = cos(AngProj*3.14159265358979323846/180.);
  float si = sin(AngProj*3.14159265358979323846/180.);
  y = co*xx[1] + si*xx[2];
  IPMm(x,y);	
  
}

/*!

  substitutute  HIGZ IPL (draw polyline) functions to 
  get rid of junk lines at big zooms (bug in HIGZ) 

*/
void TDisplay::IPLm(int n, float* xx, float* yy)
{
  if(n < 2) return;

  float x[2],y[2];
  float x1,y1,x2,y2;
  float xl=Xvf[0];
  float xh=Xvf[1];
  float yl=Yvf[0];
  float yh=Yvf[1];
  

  for(int i=0; i < n-1; i++){ // loop over poliline segments

    x1=xx[i];   y1=yy[i];   
    x2=xx[i+1]; y2=yy[i+1]; 
    
    // segment is fully outside the window
    if( (x1 < xl && x2 < xl) || ( x1 > xh && x2 > xh) ||
	(y1 < yl && y2 < yl) || ( y1 > yh && y2 > yh) ){
      continue;
    }

    // segment is fully inside the window
    if( (xl < x1 && x1 < xh) && (yl < y1 && y1 < yh) &&
	(xl < x2 && x2 < xh) && (yl < y2 && y2 < yh) ){ 
      Distort(x1,y1); Distort(x2,y2);
      x[0]=x1; y[0]=y1; x[1]=x2; y[1]=y2;
      IPL(2,x,y); 
      continue;
    }

    // find 4 cross points of the line, defined by segment, with window borders

    float xc[4],yc[4]; // 0 - with xl, 1 - with xh, 2 - with yl, 3 - with yh
    float dx,dy, yp, xp;
    
    dx=x2-x1;
    dy=y2-y1;
    
    // with xl
    xc[0]=xl;
    if(dx == 0.){
      yc[0]=1.E+20;
    } else {
      yp = (y2-y1)/dx;
      yc[0]=y1+yp*(xc[0]-x1);
    }
    // with xh
    xc[1]=xh;
    if(dx == 0.){
      yc[1]=1.E+20;
    } else {
      yp = (y2-y1)/dx;
      yc[1]=y1+yp*(xc[1]-x1);
    }

    // with yl
    yc[2]=yl;
    if(dy == 0.){
      xc[2]=1.E+20;
    } else {
      xp = (x2-x1)/dy;
      xc[2]=x1+xp*(yc[2]-y1);
    }

    // with yh
    yc[3]=yh;
    if(dy == 0.){
      xc[3]=1.E+20;
    } else {
      xp = (x2-x1)/dy;
      xc[3]=x1+xp*(yc[3]-y1);
    }

    int nx=0;
    if(yl < yc[0] && yc[0] < yh) {
      x[nx]=xc[0]; y[nx]=yc[0];
      nx++;
    } 
    if(yl < yc[1] && yc[1] < yh) {
      x[nx]=xc[1]; y[nx]=yc[1];
      nx++;
    } 

    if(xl < xc[2] && xc[2] < xh) {
      x[nx]=xc[2]; y[nx]=yc[2];
      nx++;
    } 

    if(xl < xc[3] && xc[3] < xh) {
      x[nx]=xc[3]; y[nx]=yc[3];
      nx++;
    } 
    if(nx != 2) continue; // segment do not cross window
    
    // ends are fully outside the window
    if( (xl > x1 || x1 > xh) && (xl > x2 || x2 > xh) ||
	(yl > y1 || y1 > yh) && (yl > y2 || y2 > yh) ){
      // ...
      // crosspoint has to be calculated. Not yet done
      IPL(2,x,y); 
      continue;
    }
    
    // other cases (not yet treated) ... to be done
    Distort(x1,y1); Distort(x2,y2);
    x[0]=x1; y[0]=y1; x[1]=x2; y[1]=y2;
    IPL(2,x,y); 
    
  }// end of loop over poliline segments

}


/*!

  substitutute HIGZ IPM (draw polymarker) functions to 
  ignore  polymarkers outside the current view field 
  (to speedup drawing)  

*/
void TDisplay::IPMm(float x, float y)
{
  if(x < Xvf[1] && x > Xvf[0] &&
     y < Yvf[1] && y > Yvf[0]){
    Distort(x,y);
    IPM(1,&x,&y);
  }
}

/*!
  Function to apply
  any kinds of distorstion or coordinate transformation
  to drawing primitives.
  Call is changing x and y
 */
void TDisplay::Distort(float& x, float& y)
{
  int mode = DrOpt[15];
  float x0(x);
  float y0(y);

  float yscale;
  switch(mode)
    {
    case 0:
      // no distorsion
      x=x0; y=y0;
      return;
    case 1:
      yscale = 2000./fabs(x);
      x=x0;
      y=yscale*y0;
      return;
      
    default:
      return;
    }
}

















