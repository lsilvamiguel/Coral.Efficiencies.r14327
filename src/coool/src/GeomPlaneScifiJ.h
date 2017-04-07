#ifndef __GeomPlaneScifiJ__
#define __GeomPlaneScifiJ__

#include "GeomPlane.h"

class GeomPlaneScifiJ : public GeomPlane {
  
 public:
  GeomPlaneScifiJ(int id,const char* name,int nwir,
		  double x,double y,double z,
		  double dx,double dy,double dz,
		  double ang,double pitch) 
    :   GeomPlane(id,name,nwir,x,y,z,dx,dy,dz,ang,pitch) {}
  

  std::set<CCluster1>& Clusterize(const std::map<int,CDigit1>& digits);

};

#endif














