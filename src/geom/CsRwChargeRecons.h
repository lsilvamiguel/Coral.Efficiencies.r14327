#ifndef CsRwChargeRecons_h
#define CsRwChargeRecons_h
#include <iostream>
#include "CsHistograms.h"
#include "CsECAL1.h"

typedef struct  ec1gmm{
       int  zone;
       double  X;
       double  Y;
       double  Zx;
       double  Zy;
       double  E;
       double nhits;
       double  Mom;
       CsHelix hel;
       int    nfirst;
       int    npart;
       int    calobj;
} ec1gmm;

class CsRwChargeRecons {
 public:
  static CsRwChargeRecons* Instance(void);
  CsRwChargeRecons();
  ~CsRwChargeRecons();
  void RwChargeRecons(void);
  
 private:
  
  static CsRwChargeRecons* instance_;       // The Singleton Static pointer

  void  gm_gm_distance(std::vector<ec1gmm> & ec1,  int ng, 
		       double & Dmin, double& Dx, double& Dy, double & x, double& y,
		       double& zx, double & x1, double & y1,double& zy1);
  void  sum_digits(int ibg, double ph_X,double ph_Y, double ph_Z,
				  double Summ_Range,double hits[], int& pl_fired);
  bool gmchk(double Xg, double Yg, double delta);
  bool htchk(int ipl, double Xg, double Yg, int nz, double delta);


  //     coral.options file parameters

  bool RW_MC;
  int RW_histoLevel;
  bool readcal;
  double RW_Summ_Range[3];
  double deltaPR_overlap[3];
  double deltaSpace_overlap[3];
  //
  static const int nplanes=8;

  CsDetector*  RWdet[nplanes];

  double RW_Plane_X[nplanes];
  double RW_Plane_Y[nplanes];
  double RW_Plane_Z[nplanes];
  
  CsECAL1 *calorimeter;
  CsDet *idet1;
  
  //       histograms
  
  CsHist1D  *hi1D[20];
  CsHist2D  *hi2D[10];
};

#endif
