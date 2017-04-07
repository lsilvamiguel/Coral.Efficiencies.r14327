/*!
  Function to calculate distance from the
  point to the line  segment
  in 2-dimension space. 
  
  \param xpt X coordinate of the point
  \param ypt Y coordinate of the point
  \param xln X coordinates of the line end points (array of 2)
  \param yln Y coordinates of the line end points (array of 2)

  \return smallest  distance to the line segment

*/

#include "TAlgo.h"
#include <math.h>
#include <float.h>

float TAlgo::Point2Line(float xpt, float ypt, float xln[], float yln[])
{
  float a,b,xc,yc,x,y;

  if(fabs(xln[1]-xln[0]) < FLT_EPSILON){
    if(yln[1] < yln[0]){
      y=yln[1]; yln[1]=yln[0]; yln[0]=y;
    }
  } else if(xln[1] < xln[0]){
    x=xln[1]; xln[1]=xln[0]; xln[0]=x;
    y=yln[1]; yln[1]=yln[0]; yln[0]=y;
  }

  if(fabs(xln[1]-xln[0]) < FLT_EPSILON) {
    xc=xln[0];
    yc=ypt;
    if(yc < yln[0]) return((float)sqrt((yln[0]-ypt)*(yln[0]-ypt)+(xln[0]-xpt)*(xln[0]-xpt)));
    if(yc > yln[1]) return((float)sqrt((yln[1]-ypt)*(yln[1]-ypt)+(xln[1]-xpt)*(xln[1]-xpt)));
  } else {
    a=(yln[1]-yln[0])/(xln[1]-xln[0]);
    b=yln[0]-a*xln[0];
    xc=(xpt+a*(ypt-b))/(1.+a*a);
    yc=a*xc+b;
    if(xc < xln[0]) return((float)sqrt((yln[0]-ypt)*(yln[0]-ypt)+(xln[0]-xpt)*(xln[0]-xpt)));
    if(xc > xln[1]) return((float)sqrt((yln[1]-ypt)*(yln[1]-ypt)+(xln[1]-xpt)*(xln[1]-xpt)));
  }
  return((float)sqrt((yc-ypt)*(yc-ypt)+(xc-xpt)*(xc-xpt)));
}
