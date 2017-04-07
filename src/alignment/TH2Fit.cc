// $Id: TH2Fit.cc,v 1.5 2008/02/20 13:30:13 rgazda Exp $
#include <TMath.h>
#include "TH2Fit.h"
#include <TH2.h>
#include <TF1.h>
#include <TMinuit.h>
using namespace std;
ClassImp( TH2Fit )

//__________________________________
TH2* TH2Fit::h_static_ = 0;
TF1* TH2Fit::f_static_ = 0;
unsigned int TH2Fit::nP_static_ = 0;

// Rem: xMax<xMin means that histogram range is used
double TH2Fit::xMin_ = 1;
double TH2Fit::xMax_ = -1;


//___________________________
TH2Fit::TH2Fit( TF1* f, const unsigned int nP ): 
  TObject(), 
  f_( f ), nP_( nP )
{
  // initialise minuit
  gMinuit_ = new TMinuit( nP_ );
  gMinuit_->mninit(5,6,7);
  
  // set starting parameters 
  double flag = -1;
  int error;
  gMinuit_->mnexcm("SET PRINT", &flag  ,1,error);  // disable all printouts
  for( unsigned int iP=0; iP<nP_; iP++ ) {
    char* pName = new char[15];
    sprintf(pName,"parameter_%i",iP );
    gMinuit_->mnparm(iP,pName, f_->GetParameter( iP ), 1, 0, 0, error);
    if (error) cout << "TH2Fit::Fit - ERROR:Troubles defining parameter" << iP << endl;
    SafeDelete( pName );
  }
  return;
}  
  
//__________________________________________________________________
int TH2Fit::ExecMinuitCommand( const char* command, double* pList, int size )
{  
  int error;
  gMinuit_->mnexcm( command, pList, size, error);
  return error;
}
 
//__________________________________________________________________
bool TH2Fit::Fit( TH2* h, const double xMin, const double xMax ) 
{ 
  
  // copy function and histogram into static vars
  f_static_ = f_;
  nP_static_ = nP_;
  h_static_ = h;
  
  // set limits
  if( xMax < xMin ) {
    cout << "TH2Fit::Fit - INFO: Fit range is histogram range.\n";
    xMin_ =  h_static_->GetXaxis()->GetBinCenter( 1 );
    xMax_ =  h_static_->GetXaxis()->GetBinCenter( h_static_->GetNbinsX() );
  } else {
    xMin_ = xMin;
    xMax_ = xMax;
  }
  
  // set Minuit minimisation function 
  gMinuit_->SetFCN(fcn);
   
  // Do the fit 
  int error;
  double flag0 = 0, flag1 = 1, flag3 = 3;
  gMinuit_->mnexcm("CALL FCN",  &flag1 ,1,error);
  gMinuit_->mnexcm("SET NOGradient", &flag0 ,1,error);
  gMinuit_->mnexcm("MIGRAD",   &flag0 ,0,error);
  gMinuit_->mnexcm("CALL FCN", &flag3 ,1,error);
  
  return true;
}
  
//___________________________________________________________
bool TH2Fit::FitInverted( TH2* h, const double yMin, const double yMax  ) 
{ 

  // copy function and histogram into static vars
  f_static_  = f_;
  nP_static_ = nP_;
  h_static_ = h;
  
  // set limits
  if( yMax < yMin ) {
    cout << "TH2Fit::Fit - INFO: Fit range is histogram range.\n";
    xMin_ =  h_static_->GetXaxis()->GetBinCenter( 1 );
    xMax_ =  h_static_->GetXaxis()->GetBinCenter( h_static_->GetNbinsX() );
  } else {
    xMin_ = yMin;
    xMax_ = yMax;
  }
  
  // set Minuit minimisation function 
  gMinuit_->SetFCN(fcnInv);

  // Do the fit 
  int error;
  double flag0 = 0, flag1 = 1, flag3 = 3;
  gMinuit_->mnexcm("CALL FCN",  &flag1 ,1,error);
  gMinuit_->mnexcm("SET PRINT", &flag0 ,1,error);
  gMinuit_->mnexcm("SET NOGradient", &flag0 ,1,error);
  gMinuit_->mnexcm("MIGRAD",   &flag0 ,0,error);
  gMinuit_->mnexcm("CALL FCN", &flag3 ,1,error);
  return true;
}
  
//_____________________________________________________
void TH2Fit::fcn( int& npar, double *gin, double &res, double *par, int flag )
{
  // Get Histogram and related infos 
  if( !h_static_ ) {
    cout << "TH2Fit::fcn - FATAL: histogram not set" << endl;
    return;
  }
  unsigned int nX = h_static_->GetNbinsX();
  unsigned int nY = h_static_->GetNbinsY();
  
  // Get Function and set parameters 
  if( !f_static_ ) {
    cout << "TH2Fit::fcn - FATAL: function not set" << endl;
    return;
  }
  
  for( unsigned int iP = 0; iP < nP_static_; iP++ )
  f_static_->SetParameter( iP, par[iP] );
    
  // Scan all bins, calculate result 
  double nEnt = 0;
  res = 0;
  for( unsigned int iX = 1; iX <= nX; iX++ )
  for( unsigned int iY = 1; iY <= nY; iY++ ) {
    double x = h_static_->GetXaxis()->GetBinCenter( iX );
    double y = h_static_->GetYaxis()->GetBinCenter( iY );
    double n = h_static_->GetBinContent( iX, iY );
    if( x >= xMin_ && x < xMax_ ) {
      res+= n* pow( y - f_static_->Eval( x ), 2 );
      nEnt+=n;
    }
  }
  res/=nEnt;
  
  return;
}
   
//_____________________________________________________
void TH2Fit::fcnInv( int& npar, double *gin, double &res, double *par, int flag )
{
  // Get Histogram and related infos 
  if( !h_static_ ) {
    cout << "TH2Fit::fcnInv - FATAL: histogram not set" << endl;
    return;
  }
  unsigned int nX = h_static_->GetNbinsX();
  unsigned int nY = h_static_->GetNbinsY();
  
  // Get Function and set parameters 
  if( !f_static_ ) {
    cout << "TH2Fit::fcnInv - FATAL: function not set" << endl;
    return;
  }
  
  for( unsigned int iP = 0; iP < nP_static_; iP++ )
  f_static_->SetParameter( iP, par[iP] );
  
  // Scan all bins, calculate result 
  double nEnt = 0;
  res = 0;
  for( unsigned int iX = 1; iX <= nX; iX++ )
  for( unsigned int iY = 1; iY <= nY; iY++ ) {
    double x = h_static_->GetXaxis()->GetBinCenter( iX );
    double y = h_static_->GetYaxis()->GetBinCenter( iY );
    double n = h_static_->GetBinContent( iX, iY );
    
    // warning xMin_ and xMax here apply to y!
    if( y >= xMin_ && y < xMax_ ) {
      res+= n* pow( x - f_static_->Eval( y, 0, 0 ), 2 );
      nEnt+=n;
    }
  }
  res/=nEnt;
  return;
}
