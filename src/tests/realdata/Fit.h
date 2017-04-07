// $Id: Fit.h,v 1.1 2002/12/09 08:50:36 hpereira Exp $
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
  static double P1( double *x, double *par); //!< Fit DU (residual) vs anything, linear
  static double GG( double *x, double *par); //!< two gaussians (signal+bg) (6 parameters)
  static double GP0( double *x, double *par );  //!< gaussian + P0, 4 parameters
  static double EE( double *x, double *par );   //!< two Erf
  static double EEG( double *x, double *par );  //!< two Erf(4 pars) + a gaussian (3 pars) 
  static double EEP0( double *x, double *par);  // two Erf + P0 (5 pars) 
  static double driftTFit( double *x, double *par); //!< Drift time distribution fit, one Erf on leading edge.
  static void   resFcn( int &npar, double *gin, double &f, double *x, int flag );
  ClassDef(Fit,0)
};


#endif
