#ifndef __GeomPlaneMumega__
#define __GeomPlaneMumega__

#include "GeomPlane.h"

class GeomPlaneMumega : public GeomPlane {
 private:
  static const float fF1_TICK;
  double Resolution(double wire) const;

 public:
  GeomPlaneMumega(int id,const char* name,int nwir,
		  double x,double y,double z,
		  double dx,double dy,double dz,
		  double ang,double pitch) 
    :   GeomPlane(id,name,nwir,x,y,z,dx,dy,dz,ang,pitch) {}
  
  std::set<CCluster1>& Clusterize(const std::map<int,CDigit1>& digits);

  /// predicate for digit comparison using tot
  struct lessTot {
    bool operator() (const CDigit1& d1, const CDigit1& d2) {
      if(d1.dt.size()<2 || d2.dt.size()<2) return true;
      return d1.dt[1]<d2.dt[1];
    }
  };
  
  //ClassDef(GeomPlaneMumega,0)
}; 

#endif
