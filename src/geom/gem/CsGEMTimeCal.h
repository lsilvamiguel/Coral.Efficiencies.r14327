#ifndef CsGEMTimeCal_H
#define CsGEMTimeCal_H

/*
--------------------------------------------------------------------------

 Declaration of classes to store time calibrations for GEM / PixelGEM
 detectors

--------------------------------------------------------------------------

class CsGEMTimeCal :

                 Object holding the parameters for time reconstruction with
                 one amplitude ratio.

--------------------------------------------------------------------------

class CsGEMTimeCals :

                 Object holding two objects of type CsGEMTimeCal, each
                 for one ratio.

--------------------------------------------------------------------------
 v1.0    03/2007         Markus Kr�mer, Thiemo Nagel, Bernhard Ketzer
 v2.0    26.05.2008      Bernhard Ketzer
--------------------------------------------------------------------------

*/

// C++ headers
#include <iostream>
#include <vector>

//-----------------------------------------------------------------------------
// Class definition CsGEMTimeCal
//-----------------------------------------------------------------------------

class CsGEMTimeCal {

 protected:

  bool fvalid;
  bool fnew;

 public:

  virtual       ~CsGEMTimeCal() {};
  bool          IsValid() const { return  fvalid; };
  void          IsValid(const bool _valid) { fvalid=_valid; };
  bool          IsNew() const { return  fnew; };
  void          IsNew(const bool _new) { fnew=_new; };

};

//-----------------------------------------------------------------------------
//Class definition CsGEMTimeCalOld
//-----------------------------------------------------------------------------

class CsGEMTimeCalOld : public CsGEMTimeCal {

 private:

  // Parameters to reconstruct time
  double frmin;
  double frmax;
  double ft0;
  double fsl;
  double ftc;

 public:

  CsGEMTimeCalOld();
  CsGEMTimeCalOld(const double& _a, const double& _b, const double& _c,
                  const double& _d, const double& _e);

  virtual       ~CsGEMTimeCalOld() {};

  double        rmin() const { return frmin; }
  double        rmax() const { return frmax; }
  double        t0() const { return ft0; }
  double        sl() const { return fsl; }
  double        tc() const { return ftc; }
  bool          CalcTime(const float& _a1, const float& _a2, double &_time, double &_timeerr);

};

//-----------------------------------------------------------------------------
//Class definition CsGEMTimeCalNew
//-----------------------------------------------------------------------------

class CsGEMTimeCalNew : public CsGEMTimeCal {

 private:

  // Parameters to reconstruct time
  double fa0;
  double ft0;
  double fr0;
  std::vector<double> fcov;

 public:

  CsGEMTimeCalNew();
  CsGEMTimeCalNew(const double& _a, const double& _t, const double& _r,
                  const double* _cov);
  CsGEMTimeCalNew(const double& _a, const double& _t, const double& _r,
                  const std::vector<double> &_cov);

  virtual       ~CsGEMTimeCalNew() {};

  double                a0() const { return fa0; }
  double                t0() const { return ft0; }
  double                r0() const { return fr0; }
  std::vector<double>   cov() const { return fcov; }
  double                cov(unsigned int i) { if ( i >= fcov.size() ) return 0.; return fcov[i]; }
  bool                  CalcTime(const float& _a1, const float& _a2, const float& _sigma, double &_time, double &_timeerr);

};

//-----------------------------------------------------------------------------
//Class definition CsGEMTimeCals
//-----------------------------------------------------------------------------

class CsGEMTimeCals {

 private:

  std::vector<CsGEMTimeCal*> fCals;
  bool fnew;

 public:

  CsGEMTimeCals();
  CsGEMTimeCals(CsGEMTimeCal *_cal0, CsGEMTimeCal *_cal1);

  virtual ~CsGEMTimeCals() {};

  //  void operator>> (std::istringstream &_in);
  void                Set(CsGEMTimeCal* _cal0, CsGEMTimeCal* _cal1);
  const CsGEMTimeCal* GetCal(unsigned int i) const;
  int                 GetSize()  { return fCals.size(); }
  bool                CalcTime(const std::vector<float>& _amp, const float& _sigma, double &_time, double &_timeerr) const;
  void                SetIsNew(bool _isnew) { fnew = _isnew; }
  bool                IsNew() const { return fnew; }
  bool                IsValid() const;
  void                Print() const;
  void                Clear();

};

std::istream& operator>> (std::istream &_in, CsGEMTimeCals& _cals);

#endif
