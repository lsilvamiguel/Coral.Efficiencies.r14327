#ifndef __Geometry__
#define __Geometry__

#include <string>
#include <map>

#include "TROOT.h"
#include "TMatrix.h"
#include "GeomPlane.h"

class Geometry : public TObject {

 private:
  static const float fScale;
  const char* fFilename;

  std::map<std::string,GeomPlane*> fAllPlanes;

 public:
  /// constructor 
  Geometry(const char *geofilename);
  ~Geometry() {}

  /// access a given tracker. const qualifers have been removed WARNING !!!!
  GeomPlane* GetPlane(const char* tbname);

  /// converts geant coordinates to coool coordinates
  void Geant2Coool(double &x, double &y, double &z);
  void Geant2Coool(double &x);
  void Geant2Coool(TMatrix &m);

  // returns rotation angle of DRS w/r MRS arouns beam axis
  double DRSAngle(const TMatrix *m); 
  ClassDef(Geometry,0)
};

#endif




