#ifndef __GeomPlaneSili__
#define __GeomPlaneSili__

#include "GeomPlane.h"
#include "fstream"

class GeomPlaneSili : public GeomPlane {
  
 private:
  float Resolution1;
  float Resolution2;
  float Resolution3;
  bool firstread;

 public:

  GeomPlaneSili(int id,const char* name,int nwir,
		double x,double y,double z,
		double dx,double dy,double dz,
		double ang,double pitch) 
    :   GeomPlane(id,name,nwir,x,y,z,dx,dy,dz,ang,pitch) {firstread=true;}
  
  std::set<CCluster1>& Clusterize(const std::map<int,CDigit1>& digits);
}; 

#endif

