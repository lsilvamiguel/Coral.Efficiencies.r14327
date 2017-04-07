#ifndef __GeomPlanePmumega__
#define __GeomPlanePMumega__

#include "PlanePMumega.h"

#include "CsMumegaPlane.h"
#include "CsPixelMumegaPlane.h"
#include "CsRectPixelMumegaPlane.h"

class GeomPlanePMumega : public GeomPlane {
  private :
    CsMumegaPlane *fPlane;
    CsPixelMumegaPlane *fPlanePix;
    CsRectPixelMumegaPlane *fPlanePixMM;

    // for the second pixel projections pitch and number of wirtes may differ
    int fNwirV;
    double fPitchV;

  protected:  
   int    fDefFlag;
   float  fDefPed;
   float  fDefSigma;

  public :

    GeomPlanePMumega(int id, const char* name, Int_t nwires, double x, double y, double z, double dx, double dy, double dz, double ang, double pitch);

  virtual ~GeomPlanePMumega();

  virtual void AddSubPlane(int id, const char *name, Int_t nwires, double x, double y, double z, double dx, double dy, double dz, double angle, double inpitch);
  virtual void SetParameters(float ts = 3., float tc = 5., int s = 2);

#if USE_DATABASE == 1
  virtual void SetCalibrations(const std::map<unsigned int, PlanePMumega::APVcalib> &c, const CS::Chip::Maps *maps, const std::vector<float> &ct);
#else
  // cannot load calibrations without access to the database, so we need to use default calibrations
  virtual void SetCalibrations(const CS::Chip::Maps *maps);
#endif

  virtual void ResetPlane();

  virtual bool CalcTime(const CDigit1 &digit, double &time);

  virtual std::set<CCluster1>& Clusterize(const std::map<int,CDigit1>& digits);

  int GetNWiresV() const { return fNwirV; }
  double GetPitchV() const { return fPitchV; }

  void Pad2Pos(double padX, double padY, double &posX, double &posY);
};

using namespace std;

#endif
