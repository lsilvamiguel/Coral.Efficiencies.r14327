// $Id: Fit.cc,v 1.2 2007/06/01 09:13:52 conrad Exp $
#include "Fit.h"
#include "ResManager.h"
#include "DetectorInfo.h"

#include <math.h>
#include <TROOT.h>
#include <TMath.h>
#include <TMatrix.h>

#include <iostream>

using namespace std;

ClassImp(Fit)


//_______________________________________________________________________________
double Fit::P1( double *x, double *par) 
{
  double u = x[0];
  double du = par[0]+par[1]*u;
  return du;
}


//_______________________________________________________________________________
double Fit::GG( double *x, double *par) 
{
  double xx = x[0];
  return par[0]*exp( -0.5*pow( (xx-par[1])/par[2], 2 ) )+ par[3]*exp( -0.5*pow( (xx-par[4])/par[5], 2 ) );
}

//_______________________________________________________________________________
double Fit::GP0( double *x, double *par) 
{
  double xx = x[0];
  return par[0]*exp( -0.5*pow( (xx-par[1])/par[2], 2 ) )+ par[3];
}

//_______________________________________________________________________________
double Fit::EE( double *x, double *par) 
{
  double xx = x[0];
  double f = par[0]*TMath::Erf((*x-(par[1]-par[2]))/par[3])-par[0]*TMath::Erf((*x-(par[1]+par[2]))/par[3]);
  return f;
}

//_______________________________________________________________________________
double Fit::EEP0( double *x, double *par) 
{
  double xx = x[0];
  double f = par[0]*TMath::Erf((*x-(par[1]-par[2]))/par[3])-par[0]*TMath::Erf((*x-(par[1]+par[2]))/par[3]);
  f+= par[4];
  return f;
}

//_______________________________________________________________________________
double Fit::EEG( double *x, double *par) 
{
  double xx = x[0];
  double f = par[0]*TMath::Erf((*x-(par[1]-par[2]))/par[3])-par[0]*TMath::Erf((*x-(par[1]+par[2]))/par[3]);
  f+= par[4]*exp( -0.5*pow( (xx-par[5])/par[6], 2 ) );
  return f;
}

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
  // Polynom fit
  double tC = x[0] - par[0];
  double r=tC*par[1]+pow(tC,2)*par[2]+pow(tC,3)*par[3]+pow(tC,4)*par[4];
  
  return r;
}

//_______________________________________________________________________________
double Fit::gridFit( double *x, double *par) 
{
  double r = x[0] - par[2];
  return par[0]*(exp( -0.5*pow( (r-par[1])/par[3], 2 ) ) + exp( -0.5*pow( (r+par[1])/par[3], 2 ) ) );
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
  for( unsigned int i=0; i < ResPtr_->dets_.size(); i++ ) res.push_back( x[i] );
  
  f = 0;
  for( unsigned int i=0; i < ResPtr_->dets_.size(); i++ )
  f += pow( ResPtr_->dets_[i]->resid_ - sqrt( pow( ResPtr_->_ResCalc( i, res ),2 )+ pow( res[i],2 ) ), 2 ); 
  f *= 100;
  return;
}
