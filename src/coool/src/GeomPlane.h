#ifndef __GeomPlane__
#define __GeomPlane__

#include <cassert>
#include <map>
#include <set>
#include <string>

//#include "TObject.h"
#include "Containers.h"

class GeomPlane {

 private:

  /// ID
  int fId;
  
  /// name
  std::string fName;
  //fNameC used for SI reading the spacial resolution in GeomPlaneSili
  const char* fNameC;
  
  /// number of wires
  int fNwir;
  
  /// center position
  double fX, fY, fZ;
  
  /// size
  double fDx, fDy, fDz;

  /// angle
  double fAngle, fCos, fSin;

  /// pitch map
  std::map<int, double> fPitchs;

  /// half size
  double fHalfsize;

 protected:

  /// set of clusters (filled in Clusterize)
  std::set<CCluster1> fClusters;

  static const double sqrt12;

 public:
  /// constructor
  GeomPlane(int id,const char* name,int nwir,
	    double x,double y,double z,
	    double dx,double dy,double dz,
	    double ang,double pitch);
  
  /// destructor
  virtual ~GeomPlane();


  /// returns position of the cluster w/r to the center, in cm.
  double Wire2Pos(double wire) const;  

  /// used for variable pitch detectors
  virtual void AddSubPlane(int id,const char* name,int nwir,
		                   double x,double y,double z,
                 		   double dx,double dy,double dz,
                		   double ang,double pitch);
  
  int GetID() const {return fId;}
  //  const char* GetName() const {return fName.c_str();}
  const char* GetName() const {return fNameC;}
  const char* GetNameC() {return fNameC;}
  int GetNWires() const {return fNwir;}
  double GetX() const {return fX;}
  double GetY() const {return fY;}
  double GetZ() const {return fZ;}
  double GetDX() const {return fDx;}
  double GetDY() const {return fDy;}
  double GetDZ() const {return fDz;}
  double GetAngle() const {return fAngle;}
  double GetCos() const {return fCos;}
  double GetSin() const {return fSin;}
  double GetHalfSize() const {return fHalfsize;}
  double GetPitch() const {return (fPitchs.begin())->second;}
  
  /// clustering
  virtual std::set<CCluster1>& Clusterize(const std::map<int,CDigit1>& digits);

  /// \return set of clusters
  const std::set<CCluster1>& GetClusters() const {return fClusters;}
  
  /// \return pitch associated to wire
  double Pitch(double wire) const; 

  /// \return resolution for a cluster at the position wire
  double Resolution(double wire) const;
  
};

#endif






