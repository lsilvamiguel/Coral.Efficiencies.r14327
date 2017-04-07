#ifndef __MManGeometry__
#define __MManGeometry__

#include <string>
#include <map>

#include "TROOT.h"
#include "TMatrix.h"
#include "Tracker.h"


class MManGeometry : public TObject {

 private:
  static const float fScale;
  const char* fFilename;

  std::map<std::string,Tracker*> fAllTrackers;

 public:
  /// constructor 
  MManGeometry(const char *geofilename);
  ~MManGeometry() {}

  /// access a given tracker
  const Tracker* GetTracker(const char* tbname) const;

  /// converts geant coordinates to coool coordinates
  void Geant2Coool(double &x, double &y, double &z);
  void Geant2Coool(double &x);
  void Geant2Coool(TMatrix &m);

  // returns rotation angle of DRS w/r MRS arouns beam axis
  double DRSAngle(const TMatrix *m); 
};

#endif




