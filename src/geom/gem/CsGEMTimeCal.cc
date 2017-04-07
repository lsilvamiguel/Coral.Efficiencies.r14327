/*
  --------------------------------------------------------------------------

  Implementation of classes to store time calibrations for GEM / PixelGEM
  detectors

  Authors: Markus Kr�mer / Thiemo Nagel / Bernhard Ketzer  March 2007

  --------------------------------------------------------------------------
*/

// Class declarations
#include "CsGEMTimeCal.h"

// C++ headers
#include <cmath>
#include <iostream>

//-----------------------------------------------------------------------------
// Constructors CsGEMTimeCalOld
//-----------------------------------------------------------------------------
CsGEMTimeCalOld::CsGEMTimeCalOld() {
  frmin = frmax = ft0 = fsl = ftc = 0.;
  fvalid = false;
  fnew = false;
};

CsGEMTimeCalOld::CsGEMTimeCalOld(const double& _a, const double& _b, const double& _c,
                                 const double& _d, const double& _e) {
  frmin = _a;
  frmax = _b;
  ft0 = _c;
  fsl = _d;
  ftc = _e;
  fvalid = true;
  fnew = false;
};

//-----------------------------------------------------------------------------
// Time calculation for one ratio with old calibrations
//-----------------------------------------------------------------------------
bool CsGEMTimeCalOld::CalcTime(const float& _a1, const float& _a2, double &_time, double &_timeerr) {

  // Ratio and error on ratio
  double _r=0.;
  if (_a2>0.) {
    _r = (double)(_a1)/(double)(_a2);
  } else
      return false;

  // Limitation to valid regions
  if (_r>0. && ftc/_r>=1. && _r>=frmin && _r<=frmax){
    _time    = ft0 + fsl * log(ftc/_r-1);
    _timeerr = 17.; // 17 ns fixed time error on single ratio

    // Reverse sign and add fixed constant in order to be
    // compatible with new time calculation
    const float TCS_T0 = 40.; // maximum of TCS phase distr.
    _time    = -1.*_time + TCS_T0;
  } else
      return false;

  return true;
}

//-----------------------------------------------------------------------------
// Constructors CsGEMTimeCalNew
//-----------------------------------------------------------------------------
CsGEMTimeCalNew::CsGEMTimeCalNew() {
  fa0 = ft0 = fr0 = 0.;
  fcov.clear();
  fvalid = false;
  fnew = true;
};

CsGEMTimeCalNew::CsGEMTimeCalNew(const double& _a, const double& _t, const double& _r,
                                 const double* _cov) {
  fa0 = _a;
  ft0 = _t;
  fr0 = _r;
  fcov.clear();
  for (int i = 0; i < 9; i++)
    fcov.push_back(_cov[i]);
  fvalid = true;
  fnew = true;
};

CsGEMTimeCalNew::CsGEMTimeCalNew(const double& _a, const double& _t, const double& _r,
                                 const std::vector<double> &_cov) {
  fa0 = _a;
  ft0 = _t;
  fr0 = _r;
  for (std::vector<double>::const_iterator it = _cov.begin();
       it != _cov.end(); ++it) {
    double a = (*it);
    fcov.push_back(a);
  }
  fvalid = true;
  fnew = true;
};

//-----------------------------------------------------------------------------
// Time calculation for one ratio with new calibrations
//-----------------------------------------------------------------------------
bool CsGEMTimeCalNew::CalcTime(const float& _a1, const float& _a2, const float& _sigma, double &_time, double &_timeerr) {

  if (!(_a2>0.))
    return false;

  // Ratio and error on ratio
  double _r  = (double)(_a1)/(double)(_a2);
  double _er = (_sigma + 1./sqrt(12.)) * sqrt(_a1*_a1 + _a2*_a2) / (_a2*_a2);

  // Calculate time and error
  if (!(_r>0.) || fr0/_r<=1.)
    return false;

  _time    = ft0 + fa0 * log(fr0/_r-1);

  double da = log(fr0/_r-1);
  double dt = 1.;
  double dr = fa0/(_r*(fr0/_r-1));
  double dR = (fa0*fr0/(_r*_r*(fr0/_r-1)));
  const int a(0), t(1), r(2);
  _timeerr = sqrt(da*da*fcov[3*a+a]+
                  da*dt*fcov[3*a+t]+
                  da*dr*fcov[3*a+r]+
                  dt*da*fcov[3*t+a]+
                  dt*dt*fcov[3*t+t]+
                  dt*dr*fcov[3*t+r]+
                  dr*da*fcov[3*r+a]+
                  dr*dt*fcov[3*r+t]+
                  dr*dr*fcov[3*r+r]+
                  dR*dR*_er*_er);

  return true;
}

//-----------------------------------------------------------------------------
// Constructor CsGEMTimeCals
//-----------------------------------------------------------------------------
CsGEMTimeCals::CsGEMTimeCals() {
  fnew=false;
};

CsGEMTimeCals::CsGEMTimeCals(CsGEMTimeCal *_cal0, CsGEMTimeCal *_cal1) {
  fnew = true;
  CsGEMTimeCalOld *cal0 = dynamic_cast<CsGEMTimeCalOld*>(_cal0);
  CsGEMTimeCalOld *cal1 = dynamic_cast<CsGEMTimeCalOld*>(_cal1);
  if ( ( cal0 && cal1 ) || ( !cal0 && !cal1 ) ) {
    if (!cal0) fnew = false;
    fCals.push_back(_cal0);
    fCals.push_back(_cal1);
  }
};

//-----------------------------------------------------------------------------
// Set CsGEMTimeCals
//-----------------------------------------------------------------------------
void CsGEMTimeCals::Set(CsGEMTimeCal* _cal0, CsGEMTimeCal* _cal1) {
  Clear();
  fnew = true;
  CsGEMTimeCalOld *cal0 = dynamic_cast<CsGEMTimeCalOld*>(_cal0);
  CsGEMTimeCalOld *cal1 = dynamic_cast<CsGEMTimeCalOld*>(_cal1);
  if ( ( cal0 && cal1 ) || ( !cal0 && !cal1 ) ) {
    if (!cal0) fnew = false;
    fCals.push_back(_cal0);
    fCals.push_back(_cal1);
  }
};

//-----------------------------------------------------------------------------
// Streaming operator for CsGEMTimeCals
//-----------------------------------------------------------------------------
std::istream& operator>> (std::istream &_in, CsGEMTimeCals& _cals) {

  _cals.Clear();

  // Create some variables to perform checks and store data from file
  std::vector<double> _tcal;
  double _datum;
  bool _error = false;
  bool _oldcal = false;

  // Read calibration file
  if (_in.fail()) {
    std::cerr << "CsGEMTimeCals::operator>> (std::istream &_in) : Error reading time calibrations"
              << std::endl;
    _error = true;
  }
  while (_in>>_datum) _tcal.push_back(_datum);
  if (_in.eof()) _in.clear();

  // Perform some checks
  if (_tcal.size()==10) _oldcal = true;
  else if (_tcal.size()==24) _oldcal = false;
  else {
    std::cerr << "CsGEMTimeCals::operator>> (std::istream &_in) : Unexpected size of time calibrations: "
              <<_tcal.size()<< std::endl;
    _error = true;
  }

  // Create time calibration object
  if (!_error) {
    if (_oldcal) {
      CsGEMTimeCalOld *cal0 = new CsGEMTimeCalOld(_tcal[0],_tcal[1],_tcal[2],
                                                  _tcal[3],_tcal[4]);
      CsGEMTimeCalOld *cal1 = new CsGEMTimeCalOld(_tcal[5],_tcal[6],_tcal[7],
                                                  _tcal[8],_tcal[9]);
      _cals.Set(cal0, cal1);
      _cals.SetIsNew(false);
    }
    else {
      std::vector<double> _cov;
      _cov.assign(_tcal.begin()+3,_tcal.begin()+12);
      CsGEMTimeCalNew *cal0 = new CsGEMTimeCalNew(_tcal[0],_tcal[1],_tcal[2],_cov);
      _cov.assign(_tcal.begin()+15,_tcal.begin()+24);
      CsGEMTimeCalNew *cal1 = new CsGEMTimeCalNew(_tcal[12],_tcal[13],_tcal[14],_cov);
      _cals.Set(cal0, cal1);
      _cals.SetIsNew(true);
    }
  }
  return _in;
}

//-----------------------------------------------------------------------------
// Get calibrations for one ratio (i=0,1)
//-----------------------------------------------------------------------------
const CsGEMTimeCal* CsGEMTimeCals::GetCal(unsigned int i) const {
  if (i<fCals.size()) return fCals[i];
  else return NULL;
}

//-----------------------------------------------------------------------------
// Calculate time from 3 amplitudes
//-----------------------------------------------------------------------------
bool CsGEMTimeCals::CalcTime(const std::vector<float>& _amp, const float& _sigma, double &_time, double &_timeerr) const {

  bool tv0=false, tv1=false;
  double t0, t1, te0, te1;

  // Check if calibrations valid
  if (!IsValid()) {
    std::cout<<"CsGEMTimeCals::CalcTime : time calibration not valid!"<<std::endl;
    return false;
  }

  // Old calibrations
  CsGEMTimeCalOld* cal0 = dynamic_cast<CsGEMTimeCalOld*>(fCals[0]);
  CsGEMTimeCalOld* cal1 = dynamic_cast<CsGEMTimeCalOld*>(fCals[1]);

  if (cal0 && cal1) {
    tv0 = cal0->CalcTime(_amp[0], _amp[2], t0, te0);
    tv1 = cal1->CalcTime(_amp[1], _amp[2], t1, te1);
  }

  // New calibrations
  else {
    CsGEMTimeCalNew* cal0 = dynamic_cast<CsGEMTimeCalNew*>(fCals[0]);
    CsGEMTimeCalNew* cal1 = dynamic_cast<CsGEMTimeCalNew*>(fCals[1]);

    if (cal0 && cal1) {
      tv0 = cal0->CalcTime(_amp[0], _amp[2], _sigma, t0, te0);
      tv1 = cal1->CalcTime(_amp[1], _amp[2], _sigma, t1, te1);
    }
  }

  // Calculate weighted mean if times from both ratios valid
  if (tv0 && tv1) {
    _time    = t0/(te0*te0) + t1/(te1*te1);
    _timeerr = 1./(te0*te0) + 1./(te1*te1);
    _time    = _time / _timeerr;
    _timeerr = sqrt(1./_timeerr);
    return true;
  }
  // Only first ratio
  else if (tv0) {
    _time    = t0;
    _timeerr = te0;
    return true;
  }
  // Only second ratio
  else if (tv1) {
    _time    = t1;
    _timeerr = te1;
    return true;
  }

  return false;
}

//-----------------------------------------------------------------------------
// Method to test whether calibrations for both ratios are valid
//-----------------------------------------------------------------------------
bool CsGEMTimeCals::IsValid() const {
  bool test=true;

  unsigned int size=0;
  for (std::vector<CsGEMTimeCal*>::const_iterator it=fCals.begin(); it!=fCals.end(); ++it) {
    test=test&&(*it)->IsValid(); // test will be false as soon as one calibration is not valid
    size++;
  }

  if (size!=2) // valid calibration are always of size 2
    return false;

  return test;
};

//-----------------------------------------------------------------------------
// Print calibrations
//-----------------------------------------------------------------------------
void CsGEMTimeCals::Print() const {
  std::vector<CsGEMTimeCal*>::const_iterator it;
  for (it=fCals.begin();it!=fCals.end();++it) {
    if (CsGEMTimeCalNew* cal = dynamic_cast<CsGEMTimeCalNew*>(*it)) {
      std::cout<<"Set of new calibrations:"<<std::endl;
      std::cout<<"a0 = "<<cal->a0()<<std::endl;
      std::cout<<"t0 = "<<cal->t0()<<std::endl;
      std::cout<<"r0 = "<<cal->r0()<<std::endl;
      std::vector<double> cov = cal->cov();
      for (unsigned int j=0;j<cov.size();j++) {
        std::cout<<cov[j]<<" ";
        if ((j+1)%3==0) std::cout<<std::endl;
      }
    }
    else if (CsGEMTimeCalOld* cal = dynamic_cast<CsGEMTimeCalOld*>(*it)) {
      std::cout<<"Set of old calibrations:"<<std::endl;
      std::cout<<"rmin = "<<cal->rmin()<<std::endl;
      std::cout<<"rmax = "<<cal->rmax()<<std::endl;
      std::cout<<"t0 = "<<cal->t0()<<std::endl;
      std::cout<<"sl = "<<cal->sl()<<std::endl;
      std::cout<<"tc = "<<cal->tc()<<std::endl;
    }
  }
}

//-----------------------------------------------------------------------------
// Clear calibrations
//-----------------------------------------------------------------------------
void CsGEMTimeCals::Clear() {
  // Delete calibration objects
  std::vector<CsGEMTimeCal*>::iterator it = fCals.begin();
  while ( it!=fCals.end() ) {
    delete (*it);
    it = fCals.erase(it);
  }
  fnew = false;
}

