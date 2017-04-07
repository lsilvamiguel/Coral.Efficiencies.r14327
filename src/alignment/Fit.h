// $Id: Fit.h,v 1.4 2006/06/16 15:21:47 conrad Exp $
#ifndef Fit_h
#define Fit_h

/*!
   \file    Fit.h
   \brief   some usefull functions for fits
   \author  Hugo Pereira
   \version $Revision: 1.4 $
   \date    $Date: 2006/06/16 15:21:47 $
*/

#include<TROOT.h>
#include<TObject.h>

/*!
   \class   Fit
   \brief   some usefull functions for fits
*/

class Fit: public TObject { 
  public:
  static double P1( double *x, double *par );   //!< straight line, 2 parameters
  static double GP0( double *x, double *par );  //!< gaussian + P0, 4 parameters
  static double GG( double *x, double *par );   //!< two gaussians, 6 parameters
  static double EE( double *x, double *par );   //!< two Erf
  static double EEP0( double *x, double *par ); //!< two Erf + P0 (5 pars) 
  static double EEG( double *x, double *par );  //!< two Erf(4 pars) + a gaussian (3 pars) 
  ClassDef(Fit,0)
};


#endif
