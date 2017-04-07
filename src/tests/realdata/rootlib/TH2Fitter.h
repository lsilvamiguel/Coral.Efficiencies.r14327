// $Id: TH2Fitter.h,v 1.2 2002/06/08 16:25:02 hpereira Exp $

#ifndef TH2Fitter_h
#define TH2Fitter_h

#include <string>
#include <iostream.h>
#include <TROOT.h>
#include <TObject.h>
#include <TH2.h>
#include <TF1.h>
#include <TMinuit.h>

void fcn( int& npar, double *gin, double &res, double *par, int flag );
void fcnInv( int& npar, double *gin, double &res, double *par, int flag );

class TH2Fitter: public TObject {
  public:
  static TH2Fitter* Instance();   //<! Singleton instanciation  
  
  /*! \fn void SetFunction( TF1* f, const unsigned int nP)
    initialize minuit with number of parameters \param nP.
    \param f is the fit function.
  */
  void SetFunction( TF1* f, const unsigned int nP);
  
  /*! \fn int ExecMinuitCommand( unsigned int iP )
    to execute minuit commands
  */
  int ExecMinuitCommand( const char* command, double* pList, int size );

  /*! \fn void Fit( TH2S* h, TF1* f, 
    const double xMin = 1, const double xMax = -1, 
    const unsigned int nP )
    \brief perform the fit of histogram \param h, using function \param f
    between \param xMin and \param xMax. if \param xMin > \param xMax, 
    histogram range is used.
    \param nP is the number of parameters
    starting parameters are passed through \param f
  */
  bool Fit( TH2S* h, const double xMin = 1, const double xMax = -1 );
  
  /*! \fn void FitInverted( TH2S* h, TF1* f, 
    const double yMin = 1, const double yMax = -1, 
    const unsigned int nP )
    \brief do the same thing as Fit but variables are inverted: x = \param f(y)
    is used to perform the fit between \param yMin and \param yMax.
    if \param yMin > \param yMax, histogram range is used.   
  */
  bool FitInverted( TH2S* h, const double xMin = 1, const double xMax = -1 );  
  
  bool hasF_;
  bool hasNP_;

  TF1* f_;          //<! the function used to the fit
  unsigned int nP_; //<! number of parameters
  double xMin_;     //<! lower limit for the fit
  double xMax_;     //<! upper limit for the fit
  TH2S* h_;         //<! The 2D histogram to be fitted
  static void fcn( int& npar, double *gin, double &res, double *par, int flag );
  static void fcnInv( int& npar, double *gin, double &res, double *par, int flag );
   
  protected:
  TH2Fitter( void ); //<! default Constructor
  
  private:
  static TH2Fitter* instance_;
  TMinuit *gMinuit_; 
  
  ClassDef(TH2Fitter,0)
};


#endif  
