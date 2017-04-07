/*!

  Interface to GEANT's Runge-Kutta extrapolation
  
  Input:
  parin(x0,y0,z0,yp,zp,charge/p) - input track parameters
  
  Output:
  parout(x0,y0,z0,yp,zp,charge/p) - output track parameters
  x0 of parout[] have to be preset to x coordinate 
  of the plane you want to extrapolate track to.
  
  return - 0    OK

*/

#include <math.h>
#include <iostream>
#include "TAlgo.h"
#include "TOpt.h"
#include "TConstants.h"

extern "C" void grkuta_(double*, double[], double[], int*);

int TAlgo::Rkutex(const double * const parin, double * const parout)

{
  double sign, step, deltaX, deltaX_prev, vin[7], vout[7];
  int ier=0;


  deltaX = deltaX_prev = fabs(parout[0]-parin[0]); 
  sign=((parout[0]-parin[0]) < 0 ? -1. : +1.);

  if(deltaX <= TConstants_RKuttaMinStep){ // too small extrapolation distance
    for(int i = 0; i < 6; i++) parout[i]=parin[i];
    return(0);
  } else {
    // transform to parameters for grkuta.F
    vin[0]=parin[0]; vin[1]=parin[1]; vin[2]=parin[2];             // x,y,z
    vin[3]=sign/sqrt((parin[3]*parin[3])+(parin[4]*parin[4])+1.);  // COSx 
    vin[4]=parin[3]*vin[3];                                        // COSy
    vin[5]=parin[4]*vin[3];                                        // COSz
    vin[6]=parin[5]*sign;                                          // charge/P

    int iter=0;
    while(deltaX >TConstants_RKuttaMinStep){ // iterations

      step = TOpt::dCut[3] < deltaX ? TOpt::dCut[3] : deltaX ;
      grkuta_(&step, vin, vout, &ier); if(ier != 0) return(ier); // call Runge-Kutta extrapolation subroutine
      deltaX=fabs(parout[0]-vout[0]);      // it will be the next step       
      if(deltaX > deltaX_prev) return(10); // seems it's not getting closer. Stop.
      deltaX_prev=deltaX;

      for(int i = 0; i < 7; i++) vin[i]=vout[i]; // output -> input for next step

    }// end of iterations

    // transform output of grkuta.F back  
    parout[0]=vout[0]; parout[1]=vout[1]; parout[2]=vout[2]; // x,y,z
    parout[3]=vout[4]/vout[3];                               // Yp
    parout[4]=vout[5]/vout[3];                               // Zp
    parout[5]=vout[6]*sign;                                  // charge/P
    return(0);
  }

}








