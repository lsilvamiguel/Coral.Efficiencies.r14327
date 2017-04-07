// $Id: THlxFindCDA.cc 13148 2011-12-28 16:55:25Z kbicker $

#include "THlx.h"
#include <iostream>

using namespace std;

/*!
  Find closest distance of approach of 2 tracks,
  represented by "this" helix and helix H
  extrapolating them backward.
  At the exit, "this" helix and helix H has
  it's values at the X where they have minimal distance.
*/


extern "C"  // C-like prototypes
{
  void minvar_(float& x, float& y, float& r, float& eps, float& step, int& maxf, 
	       float& a, float& b, float f2minimize_(float& x, int& i) ); // CERNLIB routine
  float f2minimize_(float& x, int& i);
}

namespace HELIX4CDA {
  THlx H1,H2;
}

bool THlx::FindCDA(THlx& H, float X0, float Xmin, float Xmax)
{
  // Defaults: float X0 = 0, float Xmin = -1000, float Xmax = +1000

  THlx& H1 = (*this);
  THlx& H2 = H;

  float x,cda,r,step,eps,Xcross;
  int maxf; 

  // pass 2 helices to "f2minimize_" routine
  HELIX4CDA::H1 = H1;
  HELIX4CDA::H2 = H2;

  // first aproximation of min. position
  x = X0;
  // convergency citeria
  r    = 0.1;
  eps  = 0.01;
  // first step size
  step = 50.;
  // man number of itterations (calls to f2minimize_)
  maxf = 200;

  minvar_(x, cda, r, eps, step, maxf, Xmin, Xmax, f2minimize_); // call CERNLIB minimum search

  if(fabs(Xmin-x)/(1.0+fabs(Xmin)) < 2*eps ||
     fabs(Xmax-x)/(1.0+fabs(Xmax)) < 2*eps) {
    //cout<<"THlx::FindCDA() ==> end up at one of the limits"<<endl;
    return(false); // no local minimum had been found
  }

  // Extrapolate to found min.
  if(! HELIX4CDA::H1.Extrap(x, (*this))) return false;
  if(! HELIX4CDA::H2.Extrap(x, H))       return false;
  
  //cout<<"THlx::FindCDA() ==>  MINIMUM is reached at X = "<< x <<" with CDA = "<<cda<<endl;

  return(true);
}


extern "C" float f2minimize_(float& x, int& i)
{
  static map<float,float> m;
  if(i == 0) m.clear();

  if(m.find(x) == m.end()) {
    THlx H1ext, H2ext;
    float dist;
    HELIX4CDA::H1.Extrap(x, H1ext);
    HELIX4CDA::H2.Extrap(x, H2ext);
    dist = H1ext.Dist(H2ext);
    m[x]=dist;
    //cout<<"f2minimize_ called with \t X = "<<x<<"  \t Dist = "<<dist<<endl;;
    return(dist);
  } else {
    //cout<<"f2minimize_ called with \t X = "<<x<<"  \t Mapped value is used : "<<m[x]<<endl;
    return(m[x]);
  }
}


















