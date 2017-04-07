// $Id: Fit.cc,v 1.11 2002/11/06 14:05:31 hpereira Exp $
#include "Fit.h"
#include "ResManager.h"
#include "DetInfo.h"
#include <math.h>
#include <TROOT.h>
#include <TMath.h>
#include <TMatrix.h>

ClassImp(Fit)

//_______________________________________________________________________________
double Fit::t0Fit(double *x,double *par)
{
  double r = x[0];
  double t = par[0]+fabs(r-par[1])*par[2];
  return t;
}

//_______________________________________________________________________________
double Fit::t0FitInv(double *x,double *par)
{
  double t = x[0] - par[0];
  double r = t*par[1];
  return r;
}

//_______________________________________________________________________________
double Fit::rtFit( double *x, double *par) 
{
  #ifdef _FIT1_   
  // Polynom fit
  double tC = x[0] - par[0];
  double r=tC*par[1]+pow(tC,2)*par[2]+pow(tC,3)*par[3]+pow(tC,4)*par[4];
  #else          
  // Hyperbolic tangent fit
  double cell = 0.5*par[0];
  double tC = x[0] - par[1];
  double pol=tC*par[2]+pow(tC,3)*par[3]+pow(tC,5)*par[4];
  double r = cell * tanh( pol/cell );
  #endif
  
  return r;
}

//_______________________________________________________________________________
double Fit::gridFit( double *x, double *par) 
{
  double r = x[0] - par[2];
  return par[0]*(exp( -0.5*pow( (r-par[1])/par[3], 2 ) ) + exp( -0.5*pow( (r+par[1])/par[3], 2 ) ) );
}

//_______________________________________________________________________________
double Fit::du2DFit( double *x, double *par) 
{
  double u = x[0];
  double du = par[0]+par[1]*u;
  return du;
}

//_______________________________________________________________________________
double Fit::du1DFit( double *x, double *par) 
{
  double xx = x[0];
  return par[0]*exp( -0.5*pow( (xx-par[1])/par[2], 2 ) )+ par[3]*exp( -0.5*pow( (xx-par[4])/par[5], 2 ) );
}

//_______________________________________________________________________________
double Fit::du1DFit_GAUS( double *x, double *par) 
{
  double xx = x[0];
  return par[0]*exp( -0.5*pow( (xx-par[1])/par[2], 2 ) );
}

//_______________________________________________________________________________
double Fit::du1DFit_MWPC( double *x, double *par) 
{
  double xx = x[0];
  double f = par[0]*TMath::Erf((*x-(par[1]-par[2]))/par[3])-par[0]*TMath::Erf((*x-(par[1]+par[2]))/par[3]);
  f+= par[4]*exp( -0.5*pow( (xx-par[5])/par[6], 2 ) );
  return f;
//   return par[0]*(1.0+TMath::Erf((xx-(par[1]-0.5*par[2]))/par[3])*TMath::Erf((par[1]+0.5*par[2]-xx)/par[3]));
}
//_______________________________________________________________________________
double Fit::driftTFit(double *x,double *par)
{
  double t = x[0];
  double N = par[0]+par[1]*TMath::Erf((t-par[2])/par[3]);
  return N;
}
//_____________________________________________________________________________
void Fit::resFcn( int &npar, double *gin, double &f, double *x, int flag )
{
  //=== convert and check pointer to current Manager Object
  ResManager *ResPtr_ = (ResManager*) HManager::ptr_;
  if( !ResPtr_ ) {
    cout << "resFcn - FATAL: Pointer to ResManager is NULL.\n";
    return;
  }
  
  vector< double > res; res.clear();
  for( unsigned int i=0; i < ResPtr_->selDetI_.size(); i++ ) res.push_back( x[i] );
  
  f = 0;
  for( unsigned int i=0; i < ResPtr_->selDetI_.size(); i++ )
  f += pow( ResPtr_->selDetI_[i]->resid_ - sqrt( pow( ResPtr_->_ResCalc( i, res ),2 )+ pow( res[i],2 ) ), 2 ); 
  f *= 100;
  return;
}
