#ifndef __GeomPlaneScifiG__
#define __GeomPlaneScifiG__

#include "GeomPlane.h"

class GeomPlaneScifiG : public GeomPlane {
  
 public:
  GeomPlaneScifiG(int id,const char* name,int nwir,
		  double x,double y,double z,
		  double dx,double dy,double dz,
		  double ang,double pitch) 
    :   GeomPlane(id,name,nwir,x,y,z,dx,dy,dz,ang,pitch) {}
  
  //  set<CCluster1>& Clusterize(const map<int,CDigit1>& digits);
}; 

#endif

