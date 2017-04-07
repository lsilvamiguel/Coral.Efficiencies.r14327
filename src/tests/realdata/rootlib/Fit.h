// $Id: Fit.h,v 1.6 2002/08/17 20:52:04 hpereira Exp $
#ifndef Fit_h
#define Fit_h

//=== _Fit1_ corresponds to Polynomial RTRelation fit, 
//=== otherwise hyperbolic tangent is used.
#define _FIT1_  

#include<TROOT.h>
#include<TObject.h>

class Fit: public TObject { 
  public:
  Fit( ): TObject() { };
  static double t0Fit(   double *x, double *par); //!< Returns t(r), linear
  static double t0FitInv(double *x, double *par); //!< Returns r(t), linear
  static double rtFit(   double *x, double *par); //!< Returns r(t), either polynom or ArcTH (depending on _FIT1_
  static double gridFit( double *x, double *par); //!< Fit RT slice: 2 gaussians, symetric to 0
  static double du2DFit( double *x, double *par); //!< Fit DU (residual) vs anything, linear
  static double du1DFit( double *x, double *par); //!< Fit DU (residual), two gaussians (signal+bg)
  static double du1DFit_GAUS( double *x, double *par); //!< Fit DU (residual), two gaussians (signal+bg)
  static double du1DFit_MWPC( double *x, double *par);   //!< Fit DU (residual) for MWPC like detectors, 2 Erf (leading+trailing)
  static double driftTFit( double *x, double *par); //!< Drift time distribution fit, one Erf on leading edge.
  static void   resFcn( int &npar, double *gin, double &f, double *x, int flag );
  ClassDef(Fit,0)
};


#endif
