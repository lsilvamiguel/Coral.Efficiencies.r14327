// $Id: TH2Fit.h,v 1.2 2007/06/01 09:13:52 conrad Exp $

#ifndef TH2Fit_h
#define TH2Fit_h
/*!
   \file    TH2Fit.h
   \brief   to perform minuit fit of a 2D histogram using any x vs y function.
   \author  Hugo Pereira
   \version $Revision: 1.2 $
   \date    $Date: 2007/06/01 09:13:52 $
*/

#include <string>
#include <iostream>

#include <TROOT.h>
#include <TObject.h>
class TH2;
class TF1;
class TMinuit;

/*!
   \class   TH2Fit
   \brief   to perform minuit fit of a 2D histogram using any x vs y function.
*/
class TH2Fit: public TObject {
  public:
  
  /*! \fn TH2Fit( TH2* h, TF1* f, const unsigned int nP)
    creator
    \param f is the fit function.
    \param nP the number of parameters
  */
  TH2Fit( TF1* f, const unsigned int nP );   
    
  /*! \fn int ExecMinuitCommand( unsigned int iP )
    to execute minuit commands
  */
  int ExecMinuitCommand( const char* command, double* pList, int size );

  /*! \fn void Fit( TH2* h, onst double xMin = 1, const double xMax = -1 )
    \brief perform the fit. By default, histogram range is used.
    \param h the histogram to be fitted
    \param xMin min value for fit range
    \param xMax max value for fit range
  */
  bool Fit( TH2* h, const double xMin = 1, const double xMax = -1 );
  

  /*! \fn void FitInverted( TH2* h, onst double xMin = 1, const double xMax = -1 )
    \brief perform the fit. By default, histogram range is used. The fit is inverted: x = f(y) is considered.
    \param h the histogram to be fitted
    \param xMin min value for fit range
    \param xMax max value for fit range
  */
  bool FitInverted( TH2* h, const double xMin = 1, const double xMax = -1 );  
  
  
  private:
  TF1* f_;           //!< function used for the fit
  unsigned int nP_;  //<! number of parameters
  TMinuit *gMinuit_; //!< minuit object
  
  static TH2* h_static_; //!< the histogram to be fitted
  static TF1* f_static_; //!< static copy of f_
  static unsigned int nP_static_;  //<! number of parameters
  static double xMin_;     //<! lower limit for the fit
  static double xMax_;     //<! upper limit for the fit
  static void fcn( int& npar, double *gin, double &res, double *par, int flag );    //!< minuit to be minimised function (\sum (y-f(x))^2 )
  static void fcnInv( int& npar, double *gin, double &res, double *par, int flag ); //!< minuit to be minimised function (\sum (x-f(y))^2 )
  
  
  ClassDef(TH2Fit,0)
};


#endif  
