#ifndef __GeomPlanePGEM__
#define __GeomPlanePGEM__

#include "GeomPlaneGEM.h"

#include "CsPixelGEMPlane.h"

class GeomPlanePGEM : public GeomPlaneGEM {
 private:
  CsPixelGEMPlane *fPlanePix;

  // for the second pixel projections pitch and number of wires may differ
  int    fNwirV;
  double fPitchV;

 public:
  
  GeomPlanePGEM(int id,const char* name,int nwir,
	       double x,double y,double z,
	       double dx,double dy,double dz,
	       double ang,double pitch);

  virtual ~GeomPlanePGEM();

  virtual void AddSubPlane(int id, const char *name, int nwires,
                           double x, double y, double z,
                           double dx, double dy, double dz,
                           double angle, double inpitch);

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

  int    GetNWiresV() const { return fNwirV; }
  double GetPitchV()  const { return fPitchV; }

  void   Pad2Pos(double padX, double padY, double &posX, double &posY);
};

#endif

