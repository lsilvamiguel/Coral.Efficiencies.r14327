// Header file RayTrace.h
// for the ray tracing function in RayTrace.cc
// Andreas Mutter, 29.10.2005
// some changes (to avoid memory leaks) applied, 21.11.2005
// andreas.mutter@cern.ch
// added the find lens function
#ifndef RAYTRACE_H
#define RAYTRACE_H



class LightRay
{
 private:
  double sp[3];                       // starting point
  double dir[3];                      // direction vector 
  double wl;                          // wavelength (nm)
  double n;                           // index of refraction of the current medium at wl 
 public:
  LightRay();                         // sets defaults.
  void print();                       // print out all ray parameters
  void norm();                        // normalizes the direction vector
  void SetDir(double d[3]);
  void SetStart(double s[3]);
  void SetWave(double wav);
  void SetInd(double ind);
  double getval(int k);                       // return the kth value
};



class OpticalSurface 
{
 private:
  double dist_to_previous, xshift, yshift;                      // coordinates of the center
  double xtilt, ytilt, ztilt;                                   // tilt in degrees
  double L[3], K[3];                                            // Sellmeier 1 coefficients

 protected:
  double curv;                                                  // curvature
  double asph1, asph2;                                          // aspherical coefficients for r^2 and r^4
  double radius;                                                // radius of the optical element

 public:
  OpticalSurface();                                             // Constructor, sets some defaults
  virtual ~OpticalSurface() {};                                 // default virtual destructor to remove compiler warning
  int debug;

  void SetLoc(double, double, double);                          // sets location 
  void SetTilt(double, double, double);                         // sets tilt
  void SetSell(double, double, double, double, double, double); // set Sellmeier coeffs (for microns!!!)
  void SetSurf(double, double, double);                         // set curvature and asph const
  void SetRadius(double);                                       // set radius of optical element

  int TransformRay(LightRay &photon);                          // transforms ray into local coordinate system
  int TransformRayBack(LightRay &photon);                      // transforms ray back into previous coordinate system

  double RefInd(double wavelength);                             // calculate index of refraction 

  double SurfSag(double x, double y);                           // z coordinate of surface at given x and y

  int CheckVig(LightRay photon);                                // check if photon is vignetted by surface

  int FindNormal(double x, double y, double n[3]);                // calculate normal n at point x, y

  int DoRefraction(LightRay &photon, double normal[3]);         // actually refract the ray according to Snell's law

  virtual int Refract(LightRay &photon) = 0;                    // collection of steps taken at the surface
};

class CoordinateBreak : public OpticalSurface
{
 public:
  int Refract(LightRay &photon);
};

class OpticalPlane : public OpticalSurface 
{
 public:
  int Refract(LightRay &photon);
  virtual int FindIntersection(LightRay &photon);
};
 
class Sphere : public OpticalPlane
{
 public:
  int FindIntersection(LightRay &photon);
};

class Asphere : public OpticalPlane
{
 public:
  int FindIntersection(LightRay &photon);
};

class OptSyst
{
 private:
  OpticalSurface *t[8];
 public:
  OptSyst();
  int RayTrace(double x, double y, double ax, double ay, double e);
  // takes x and y (on the first surface of the field lens) in cm,
  // the x- and y-coordinates of the normalized photon direction,
  // and the photon energy in eV 
  // and returns the PMT pixel in the following numbering scheme:
  //  y
  //  A  12 13 14 15
  //  |   8  9 10 11
  //  |   4  4  6  7
  //  |   0  1  2  3
  //  ---------------> x
  //
  // where x and y are the coordinates in the ZEMAX local frame of the
  // PMT. The PMT pixels were assumed to have a size of 4.5mm x 4.5mm and have no
  // dead zone between them.
  int RayTrace(LightRay &Photon);
  // as above, but with input values as LightRay object, to be used after FindFirstLens or TransPhotonDir
  int FindFirstLens(double xr, double yr, double zr, \
		    double xm, double ym, double zm, \
		    double xd, double yd, double e, \
		    LightRay &Photon);
  // Determines which lens is hit from the coordinates on the mirror,
  // the CsI cathode (in MRS and DRS) and sets the initial values for the 
  // LightRay object to be traced through the telescope.
  // Only needed for "old" ZEBRA files.
  int TransPhotonDir(double xr, double yr, double zr,	     \
		     double xm, double ym, double zm,	     \
		     double xd, double yd, double e,	     \
		     int dn, LightRay &Photon);
  // Transforms photon direction from MRS to the local lens system and sets photon parameters.
  // returns the cathode number (4,6,11,13)
  int CheckEff(LightRay photon,double scale_factor);
  // checks if photon wavelength is within the PMT efficiency region
};



#endif
// end of RayTrace.h
