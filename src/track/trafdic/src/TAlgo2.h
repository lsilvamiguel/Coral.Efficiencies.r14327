// $Id: TAlgo2.h,v 1.3 2004/11/17 02:30:29 ybedfer Exp $

#ifndef TAlgo2_h
#define TAlgo2_h

/*!
  \brief Algorithms
  Utility functions, algorithms and wrappers used in tracking 
  (Version specific to "lattice" alternative).
*/

#include "TConstants.h"

class TAlgo2 {

public:
  //! find track segment in projection (when in alternative "lattice")
  int FindProj(int igr, float X0, int npl, float xpl[], float res[], 
	       int fh[], int lh[], int nhit, float Uhit[], 
	       int& Ntk, int hits[], float y0[], float yp[],
	       int mode, float *Y2_trk);

  //! alternative find track segment in projection (when in alternative "lattice")
  int FindProjQ(int igr, float X0, int npl, float xpl[], float res[], 
		int fh[], int lh[], int nhit, float Uhit[], 
		int& Ntk, int hits[], float y0[], float yp[],
		int mode, float *Y2_trk,
		int *idpl, int *href);
  
  //! find track segment in space (when in alternative "lattice")
  int FindSpace(int igr, int mode, float X0,
		int npl, int idpl[], float cosa[], float sina[], float tols[], 
		int fh[], int lh[],
		int nhit, float Uhit[],
		int nproj, int Ntk_prj[], float cos_prj[], float sin_prj[], 
		float Y0_prj[10][TConstants_NTtrack_max], float Yp_prj[10][TConstants_NTtrack_max], float *Y2_prj,
		float *hinfos,
		int &Ntk, int hits[], int hused[]);
 
};
#endif







