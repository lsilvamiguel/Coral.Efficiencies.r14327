#ifndef CsMumegaTimeCal_H
#define CsMumegaTimeCal_H

/*
--------------------------------------------------------------------------

 Declaration of classes to store tima calibrations for Mumega/PixelMumega
 detectors

--------------------------------------------------------------------------

class CsMumegaTimeCal :

                 Object holding the parameters for time reconstruction with
		 one amplitude ratio.

--------------------------------------------------------------------------

class CsMumegaTimeCals :

                 Object holding two objects of type CsMumegaTimeCal, each
		 for one ratio.

--------------------------------------------------------------------------

*/

// C++ headers
#include <vector>
#include <iostream>

//--------------------------------------------------------------------------
// Class definition CsMumegaTimeCal
//--------------------------------------------------------------------------

class CsMumegaTimeCal {

 protected:

  bool fvalid;
  bool fnew;

 public:

  virtual ~CsMumegaTimeCal() {};
  bool    IsValid() const { return fvalid; };
  void    IsValid(const bool _valid) { fvalid = _valid; };
  bool    IsNew() const { return fnew; };
  void    IsNew(const bool _new) { fnew = _new; };

};

//--------------------------------------------------------------------------
// Class definition CsMumegaTimeCalOld
//--------------------------------------------------------------------------

class CsMumegaTimeCalOld : public CsMumegaTimeCal {

 private:

  // Parameters to reconstruct time
  double frmin;
  double frmax;
  double ft0;
  double fsl;
  double ftc;

 public:

  CsMumegaTimeCalOld();
  CsMumegaTimeCalOld(const CsMumegaTimeCalOld &_cal);
  CsMumegaTimeCalOld(const double _a, const double _b, const double _c,
		     const double _d, const double _e);

  virtual ~CsMumegaTimeCalOld() {};

  double rmin() const { return frmin; }
  double rmax() const { return frmax; }
  double t0() const { return ft0; }
  double sl() const { return fsl; }
  double tc() const { return ftc; }
  bool   CalcTime(float _a1, float _a2, double &_time, double &_timeerr);

};

//--------------------------------------------------------------------------
// Class definition CsMumegaTimeCalNew
//--------------------------------------------------------------------------

class CsMumegaTimeCalNew : public CsMumegaTimeCal {

 private:

  double fa0;
  double ft0;
  double fr0;
  std::vector<double> fcov;

 public:

  CsMumegaTimeCalNew();
  CsMumegaTimeCalNew(const CsMumegaTimeCalNew &_cal);
  CsMumegaTimeCalNew(const double _a, const double _t, const double _r,
		     const double* _cov);
  CsMumegaTimeCalNew(const double _a, const double _t, const double _r,
		     const std::vector<double> &_cov);

  virtual ~CsMumegaTimeCalNew() {};

  double              a0() const { return fa0; }
  double              t0() const { return ft0; }
  double              r0() const { return fr0; }
  std::vector<double> cov() const { return fcov; }
  double              cov(unsigned int i) { if ( i >= fcov.size() ) return 0.; return fcov[i]; }
  bool                CalcTime(float _a1, float _a2, float _sigma, double &_time, double &_timeerr);

};

//--------------------------------------------------------------------------
// Class definition CsMumegaTimeCals
//--------------------------------------------------------------------------

class CsMumegaTimeCals {

 private:

  std::vector<CsMumegaTimeCal*> fCals;
  bool fnew;

 public:

  CsMumegaTimeCals();
  CsMumegaTimeCals(CsMumegaTimeCal *_cal0, CsMumegaTimeCal *_cal1);

  virtual ~CsMumegaTimeCals() {};

  void                   Set(CsMumegaTimeCal* _cal0, CsMumegaTimeCal* _cal1);
  const CsMumegaTimeCal* GetCal(unsigned int i) const;
  int GetSize() { return fCals.size(); }
  bool CalcTime(const std::vector<float> _amp, const float _sigma, double &_time, double &_timeerr) const;
  void SetIsNew(bool _isnew) { fnew = _isnew; }
  bool IsNew() { return fnew; }
  bool IsValid() const;
  void Print() const;
  void Clear();

};

std::istream& operator>> (std::istream &_in, CsMumegaTimeCals& _cals);

#endif
