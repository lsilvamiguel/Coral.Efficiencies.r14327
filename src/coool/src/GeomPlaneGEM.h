#ifndef __GeomPlaneGEM__
#define __GeomPlaneGEM__

#include "GeomPlane.h"
#include "PlaneGEM.h"

#include "CsGEMPlane.h"

class GeomPlaneGEM : public GeomPlane {
 protected:  
  int    fDefFlag;
  float  fDefPed;
  float  fDefSigma;
  
  CsGEMPlane *fPlane;
  
 public:
  
  GeomPlaneGEM(int id,const char* name,int nwir,
	       double x,double y,double z,
	       double dx,double dy,double dz,
	       double ang,double pitch);

  virtual ~GeomPlaneGEM();
  
  virtual void SetParameters(float ts=3., float tc=5., int s=2);

#if USE_DATABASE==1
  virtual void SetCalibrations(const std::map<unsigned int, PlaneGEM::APVcalib> &c, const CS::Chip::Maps *maps, const std::vector<float> &ct);
#else
  // cannot load calibrations without access to the database, so we need to use
  // default calibrations
  virtual void SetCalibrations(const CS::Chip::Maps *maps);
#endif

  virtual void ResetPlane();

  virtual bool CalcTime(const CDigit1 &digit, double &time);

  virtual std::set<CCluster1>& Clusterize(const std::map<int,CDigit1>& digits);

}; 

#endif

