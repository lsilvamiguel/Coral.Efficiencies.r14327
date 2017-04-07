#ifndef __GeomPlaneMwpc__
#define __GeomPlaneMwpc__

#include "GeomPlane.h"

class GeomPlaneMwpc : public GeomPlane {
 private:
  static const float fF1_TICK;

 public:
  GeomPlaneMwpc(int id,const char* name,int nwir,
		double x,double y,double z,
		double dx,double dy,double dz,
		double ang,double pitch) 
    :   GeomPlane(id,name,nwir,x,y,z,dx,dy,dz,ang,pitch) {}
  
  //set<CCluster1>& Clusterize(const map<int,CDigit1>& digits);
}; 

#endif
