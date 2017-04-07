// $Id: TH2Fitter.cc,v 1.3 2002/06/25 15:30:23 hpereira Exp $
#include "TH2Fitter.h"

ClassImp( TH2Fitter )

//__________________________________
TH2Fitter* TH2Fitter::instance_ = 0;
TH2Fitter* TH2Fitter::Instance( void ) 
{
  if( instance_ == 0 ) instance_ = new TH2Fitter();
  return instance_;
}

//___________________________
TH2Fitter::TH2Fitter( void ): 
  TObject(), 
  f_( 0 ), hasF_( false ), 
  nP_( 0 ), hasNP_( false ),
  h_( 0 ) , xMin_( 0 ), xMax_( 0 ),
  gMinuit_( 0 ) { }
  
//__________________________________________________________
void TH2Fitter::SetFunction( TF1* f, const unsigned int nP )
{

  f_  = f;  hasF_  = true;
  nP_ = nP; hasNP_ = true;
  if( gMinuit_ != 0 ) delete gMinuit_;

  //=== initialize minuit ===
  gMinuit_ = new TMinuit( nP_ );
  gMinuit_->mninit(5,6,7);

  //=== set starting parameters ===
  for( unsigned int iP=0; iP<nP_; iP++ ) {
    int error;
    char pName[15];
    sprintf(pName,"parameter_%i",iP );
    gMinuit_->mnparm(iP,pName, f_->GetParameter( iP ), 1, 0, 0, error);
    if (error) cout << "TH2Fitter::Fit - ERROR:Troubles defining parameter" << iP << endl;
  }
  return;
}  

//__________________________________________________________________
int TH2Fitter::ExecMinuitCommand( const char* command, double* pList, int size )
{
  if( !(hasF_ && hasNP_ && gMinuit_ != 0 ) ) {
    cout << "TH2Fitter::ExecMinuitCommand - ERROR: TH2Fitter not initialized properly." << endl;
    return 0;
  } 
  
  int error;
  gMinuit_->mnexcm( command, pList, size, error);
  return error;
}
 
//=======================
//=== Fitting methods === 
//=======================

//__________________________________________________________________
bool TH2Fitter::Fit( TH2S* h, const double xMin = 1, const double xMax = -1 ) 
{ 
  //=== check everything is OK
  if( !(hasF_ && hasNP_ && gMinuit_ != 0 ) ) {
    cout << "TH2Fitter::Fit - ERROR: TH2Fitter not initialized properly." << endl;
    return false;
  }
  
  //=== collect histogram and limits
  h_ = h;
  if( xMax < xMin ) {
    cout << "TH2Fitter::Fit - INFO: Fit range is histogram range.\n";
    xMin_ =  h->GetXaxis()->GetBinCenter( 1 );
    xMax_ =  h->GetXaxis()->GetBinCenter( h->GetNbinsX() );
  } else {
    xMin_ = xMin;
    xMax_ = xMax;
  }
  
  //=== set Minuit minimisation function ===
  gMinuit_->SetFCN(fcn);
   
  //=== Do the fit ===
  int error;
  double flag0 = 0, flag1 = 1, flag3 = 3;
  gMinuit_->mnexcm("CALL FCN",  &flag1 ,1,error);
  gMinuit_->mnexcm("SET PRINT", &flag0 ,1,error);
  gMinuit_->mnexcm("SET NOGradient", &flag0 ,1,error);
  gMinuit_->mnexcm("MIGRAD",   &flag0 ,0,error);
  gMinuit_->mnexcm("CALL FCN", &flag3 ,1,error);

  cout << "TH2Fitter::Fit - Successfull. "<< endl;
  
  return true;
}
  
//___________________________________________________________
bool TH2Fitter::FitInverted( TH2S* h, const double yMin = 1, const double yMax = -1  ) 
{ 
  //=== check everything is OK
  if( !(hasF_ && hasNP_ && gMinuit_ != 0 ) ) {
    cout << "TH2Fitter::FitInverteed - ERROR: TH2Fitter not correctly initialised." << endl;
    return false;
  }
  
  //=== collect histogram and limits
  h_ = h;
  if( yMax < yMin ) {
    cout << "TH2Fitter::Fit - INFO: Fit range is histogram range.\n";
    xMin_ =  h->GetXaxis()->GetBinCenter( 1 );
    xMax_ =  h->GetXaxis()->GetBinCenter( h->GetNbinsX() );
  } else {
    xMin_ = yMin;
    xMax_ = yMax;
  }
  
  //=== set Minuit minimisation function ===
  gMinuit->SetFCN(fcnInv);

  //=== Do the fit ===
  int error;
  double flag0 = 0, flag1 = 1, flag3 = 3;
  gMinuit->mnexcm("CALL FCN",  &flag1 ,1,error);
  gMinuit->mnexcm("SET PRINT", &flag0 ,1,error);
  gMinuit->mnexcm("SET NOGradient", &flag0 ,1,error);
  gMinuit->mnexcm("MIGRAD",   &flag0 ,0,error);
  gMinuit->mnexcm("CALL FCN", &flag3 ,1,error);

  cout << "TH2Fitter::FitInverted - Successfull. "<< endl;
  return true;
}
  
//==============================
//=== minimisation functions === 
//==============================
  
//_____________________________________________________
void TH2Fitter::fcn( int& npar, double *gin, double &res, double *par, int flag )
{
  //=== Get Histogram and related infos ===
  TH2S* h = (TH2S*) TH2Fitter::Instance()->h_;
  if( h == 0 ) {
    cout << "fcn - FATAL: No TH2Fitter::Instance->h_ histogram set" << endl;
    return;
  }
  unsigned int nX = h->GetNbinsX();
  unsigned int nY = h->GetNbinsY();
  
  //=== Get Function and set parameters ===
  TF1* f = (TF1*)  TH2Fitter::Instance()->f_;
  if( f == 0 ) {
    cout << "fcn - FATAL: TH2Fitter::Instance->f_ function set" << endl;
    return;
  }
//   
//   if( npar != int(TH2Fitter::Instance()->nP_) ) {
//     cout << "fcn - FATAL: wrong number of parameters" << endl;
//     return;
//   }  
  
  for( unsigned int iP = 0; iP < TH2Fitter::Instance()->nP_; iP++ )
  f->SetParameter( iP, par[iP] );
  
  //=== get Limits ===
  double xMin =   TH2Fitter::Instance()->xMin_;
  double xMax =   TH2Fitter::Instance()->xMax_;
  
  //=== Scan all bins, calculate result ===
  double nEnt = 0;
  res = 0;
  for( unsigned int iX = 1; iX <= nX; iX++ )
  for( unsigned int iY = 1; iY <= nY; iY++ ) {
    double x = h->GetXaxis()->GetBinCenter( iX );
    double y = h->GetYaxis()->GetBinCenter( iY );
    double n = h->GetBinContent( iX, iY );
    if( x >= xMin && x < xMax ) {
      res+= n* pow( y - f->Eval( x ), 2 );
      nEnt+=n;
    }
  }
  res/=nEnt;
  
  return;
}
   
//_____________________________________________________
void TH2Fitter::fcnInv( int& npar, double *gin, double &res, double *par, int flag )
{
  //=== Get Histogram and related infos ===
  TH2S* h = (TH2S*) TH2Fitter::Instance()->h_;
  if( h == 0 ) {
    cout << "fcn - FATAL: No TH2Fitter::Instance->h_ histogram set" << endl;
    return;
  }
  unsigned int nX = h->GetNbinsX();
  unsigned int nY = h->GetNbinsY();
  
  //=== Get Function and set parameters ===
  TF1* f = (TF1*)  TH2Fitter::Instance()->f_;
  if( f == 0 ) {
    cout << "fcn - FATAL: No TH2Fitter::Instance->f_ function set" << endl;
    return;
  }
  
  for( unsigned int iP = 0; iP < TH2Fitter::Instance()->nP_; iP++ )
  f->SetParameter( iP, par[iP] );
  
  //=== get Limits ===
  double yMin =   TH2Fitter::Instance()->xMin_;
  double yMax =   TH2Fitter::Instance()->xMax_;
    
  //=== Scan all bins, calculate result ===
  double nEnt = 0;
  res = 0;
  for( unsigned int iX = 1; iX <= nX; iX++ )
  for( unsigned int iY = 1; iY <= nY; iY++ ) {
    double x = h->GetXaxis()->GetBinCenter( iX );
    double y = h->GetYaxis()->GetBinCenter( iY );
    double n = h->GetBinContent( iX, iY );
    if( y >= yMin && y < yMax ) {
      res+= n* pow( x - f->Eval( y, 0, 0 ), 2 );
      nEnt+=n;
    }
  }
  res/=nEnt;
  return;
}
