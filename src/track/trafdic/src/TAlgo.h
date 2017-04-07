#ifndef TAlgo_h
#define TAlgo_h
#include <string.h>
#include "CsGeom.h"
#include "TConstants.h"

/*!

  \brief Algorithms

  Utility functions algorithms and wrappers used in tracking

*/

namespace TAlgo {
  
  //! Wrapper to GEANT Runge-Kutta extrapolation
  int Rkutex (const double* const VecIn, double* const VecOut);

  //! GEANT Runge-Kutta, modified for Jacobian calculations. Negligible field gradients are assumed.
  bool RkutaNoGrad (double* SU, double* VO, double& Path);
  
  //! Wrapper to CORAL CsField::getField
  double Field (double p[3], double b[3]);
 
  //! Inverrt 5x5 matrix
  int Inv5(double* a, double* b);

  //! point-to-line distance
  float Point2Line(float xpt, float ypt, float xln[], float yln[]);

  //! Nested loops combinatorics
  bool NLoopComb(int N, int ns[], int i[]);

  //! find track segment in projection
  int FindProj(int igr, int ipass, float X0, int npl, float xpl[], float res[], float sigt[],
	       int fh[], int lh[], int nhit, float Uhit[], float Thit[], 
	       int& Ntk, int hits[], float y0[], float yp[], float tT[], float sT[] );
 
 //! find track segment in space
  int FindSpace(int igr, int ipass, float X0, int npl, int idpl[], int iprj[], float cosa[], float sina[], 
		float res[], float sigt[],  
		int fh[], int lh[], int nhit, float Uhit[], float Thit[],
		int nproj, int Ntk_prj[], float cos_prj[], float sin_prj[], 
		float **Y0_prj, float **Yp_prj, float **tT_prj, float **sT_prj, int& Ntk, int hits[]);

 //! fill internal alignment histograms
  void Alignment(int igr, int npl, int idpl[], int iprj[], float cosa[], float sina[], float res[], 
		 int fh[], int lh[], int nhit, float Uhit[]);

 //! Print object counters
  void PrintObjCounters();

 //! Print iformation about memory alocations of current proccess
  void PrintMem();
  
  //! Print local time
  void PrintTime();
  
};


inline double TAlgo::Field (double p[3], double b[3])
{
  float x,y,z,bx,by,bz;
  x=float(p[1]*10.); y=float(p[2]*10.); z=float(p[0]*10.); 
  CsField* MagField = CsGeom::Instance()->getCsField();
  MagField->getField(x,y,z,bx,by,bz);
  b[0]=bz*10.; b[1]=bx*10.; b[2]=by*10.;
  return(sqrt(bx*bx + by*by + bz*bz));
};

extern "C"
{
  void clops_(void);
  void opps_(const char *s,int);
}

inline void CLOPS(void) {clops_();}
inline void OPPS(const char *s) {opps_(s,strlen(s));}

#endif
