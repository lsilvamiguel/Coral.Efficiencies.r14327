// $Id: Fit.cc,v 1.4 2008/06/05 12:55:37 rgazda Exp $

/*!
   \file    Fit.cc
   \brief   some usefull functions for fits
   \author  Hugo Pereira
   \version $Revision: 1.4 $
   \date    $Date: 2008/06/05 12:55:37 $
*/

#include "Fit.h"
#include <TROOT.h>
#include <TMath.h>

ClassImp(Fit)


//_______________________________________________________________________________
double Fit::P1( double *x, double *par) 
{
  double u = x[0];
  double du = par[0]+par[1]*u;
  return du;
}

//_______________________________________________________________________________
double Fit::GP0( double *x, double *par) 
{
  double xx = x[0];
  return par[0]*exp( -0.5*pow( (xx-par[1])/par[2], 2 ) )+ par[3];
}

//_______________________________________________________________________________
double Fit::GG( double *x, double *par) 
{
  double xx = x[0];
  return par[0]*exp( -0.5*pow( (xx-par[1])/par[2], 2 ) )+ par[3]*exp( -0.5*pow( (xx-par[4])/par[5], 2 ) );
}

//_______________________________________________________________________________
double Fit::EE( double *x, double *par) 
{
  // double xx = x[0]; // unused variable (jj)
  double f = par[0]*TMath::Erf((*x-(par[1]-par[2]))/par[3])-par[0]*TMath::Erf((*x-(par[1]+par[2]))/par[3]);
  return f;
}

//_______________________________________________________________________________
double Fit::EEP0( double *x, double *par) 
{
  // double xx = x[0]; // unused variable (jj)
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
