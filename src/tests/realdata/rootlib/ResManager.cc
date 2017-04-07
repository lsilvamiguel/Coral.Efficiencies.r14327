// $Id: ResManager.cc,v 1.26 2002/11/06 14:03:40 hpereira Exp $
#include "ResManager.h"
#include "DetInfo.h"
#include "DetFileInfo.h"
#include "Utils.h"
#include "Fit.h"

#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMinuit.h>
#include <iostream.h>

ClassImp( ResManager )

//_______________________________________________________________________________
ResManager::ResManager( void ): HManager() { }
ResManager::ResManager( const char* fileSelection, bool batch = false ):HManager( fileSelection, batch ) {}
ResManager::ResManager( const HManager h ): HManager( h ) {}

//====================================
//=== ShortCuts to DetInfo Methods ===
//====================================
//__________________________________________________________________
void ResManager::FitDU( const char* selection = "*", TCut cut = "T_fnd" )
{
  GetSelection( selection );
  if( tcv_ != 0 )   tcv_->Clear();
  _DivideCanvas( selDetI_.size() );
  for( unsigned int iInf = 0; iInf < selDetI_.size(); iInf++ ) 
  {
    if( tcv_ != 0 ) tcv_->cd( iInf+1 );
    selDetI_[iInf]->FitDU( cut );
    if( tcv_ != 0 ) tcv_->Update();
  }
  
  return;
}

//__________________________________________________________________
void ResManager::FitDU_MWPC( const char* selection = "*", TCut cut = "T_fnd" )
{
  GetSelection( selection );
  if( tcv_ != 0 )   tcv_->Clear();
  _DivideCanvas( selDetI_.size() );
  
  for( unsigned int iInf = 0; iInf < selDetI_.size(); iInf++ ) 
  {
    if( tcv_ != 0 ) tcv_->cd( iInf+1 );
    selDetI_[iInf]->FitDU_MWPC( cut );
    if( tcv_ != 0 ) tcv_->Update();
  }
  
  return;
}

//__________________________________________________________________
void ResManager::FitDUvsU( const char* selection = "*", TCut cut = "T_fnd" )
{
  GetSelection( selection );

  if( tcv_ != 0 ) tcv_->Clear();
  _DivideCanvas( selDetI_.size() );
  
  for( unsigned int iInf = 0; iInf < selDetI_.size(); iInf++ ) 
  {
    if( tcv_ != 0 ) tcv_->cd( iInf+1 );
    selDetI_[iInf]->FitDUvsU( cut );
    if( tcv_ != 0 ) tcv_->Update();
  }
  
  return;
}

//__________________________________________________________________
void ResManager::FitDUvsV( const char* selection = "*", TCut cut = "T_fnd" )
{
  GetSelection( selection );

  if( tcv_ != 0 ) tcv_->Clear();
  _DivideCanvas( selDetI_.size() );
  
  for( unsigned int iInf = 0; iInf < selDetI_.size(); iInf++ ) 
  {
    if( tcv_ != 0 ) tcv_->cd( iInf+1 );
    selDetI_[iInf]->FitDUvsV( cut );
    if( tcv_ != 0 ) tcv_->Update();
  }
  
  return;
}

//__________________________________________________________________
void ResManager::DrawDU_LR( const char* selection = "*", TCut cut = "T_fnd", 
    double OffsetN = 0, 
    double OffsetP = 0 )
{
  GetSelection( selection );

  if( tcv_ != 0 ) tcv_->Clear();
  _DivideCanvas( selDetI_.size() );
  
  for( unsigned int iInf = 0; iInf < selDetI_.size(); iInf++ ) 
  {
    if( tcv_ != 0 ) tcv_->cd( iInf+1 );
    selDetI_[iInf]->DrawDU_LR( cut, OffsetN, OffsetP );
    if( tcv_ != 0 ) tcv_->Update();
  }
  
  return;
}

//__________________________________________________________________
void ResManager::DrawWRSProfiles ( const char* selection = "*", TCut cut = "" )
{
  GetSelection( selection );

  if( tcv_ != 0 )   tcv_->Clear();
  _DivideCanvas( selDetI_.size() );

  for( unsigned int iInf = 0; iInf < selDetI_.size(); iInf++ ) {
    if( tcv_ != 0 ) tcv_->cd( iInf+1 );
    selDetI_[iInf]->DrawWRSProfile( cut );
    if( tcv_ != 0 ) tcv_->Update();
  }
  
  return;
}

#define resFile "res.out"
//_______________________________________________________________________________
bool ResManager::FitResolutions( const char* selection, TCut cut = "T_fnd", const char* detFileName=0 )
{
  //=== get selection matching detFile
  vector< DetInfo* > selTmp = GetSelection( selection );
  selDetI_.clear();
  for( unsigned int i=0; i < selTmp.size(); i++ ) 
  if( selTmp[i]->hasDetFile_ || ( detFileName && selTmp[i]->SetDetFile( detFileName ) ) ) selDetI_.push_back(selTmp[i]);
 
  //=== Set residuals according for selection using cut
  if( tcv_ != 0 )   tcv_->Clear();
  _DivideCanvas( selDetI_.size() );
  for( unsigned int iInf = 0; iInf < selDetI_.size(); iInf++ ) 
  {
    if( tcv_ != 0 ) tcv_->cd( iInf+1 );
    selDetI_[iInf]->OneGaus_ = OneGaus_;
    selDetI_[iInf]->FitDU( cut );
    if( tcv_ != 0 ) tcv_->Update();
  }
  
  //=== dump resolutions to screen
  for( unsigned int i = 0; i < selDetI_.size(); i++ )
  if( selDetI_[i]->resid_ < 0 ) { 
    cout << "ResManager::FitResolutions - ERROR: troubles for \"" << selDetI_[i]->TBName_ << "\" residual.\n";
    return false;
  } else cout << "ResManager::FitResolution - INFO: \"" << selDetI_[i]->TBName_ << "\" -> " << selDetI_[i]->resid_ << "mm.\n";
  
  //=== initialise minuit
  TMinuit *gMinuit = new TMinuit((int)selDetI_.size()); 
  gMinuit->mninit(5,6,7);     // Logical units
  gMinuit->SetFCN(Fit::resFcn);

  cout << "ResManager::FitResolution - INFO: TMinuit initialization performed." << endl;
  
  double dSig = .00001;
  int error;
    
  for( unsigned int i=0; i< selDetI_.size(); i++ ) {
    char pName[16];
    sprintf(pName,"res_%s",selDetI_[i]->TBName_.c_str() );
    gMinuit->mnparm(i,pName,selDetI_[i]->detFI_->res_,dSig, 0,0,error);
    if (error) {
      cout << "ResManager::FitResolution - INFO: Troubles defining parameter " << i << ".\n";
      return false;
    }
  }
  
  double flag0 = 0, flag1 = 1, flag3 = 3;
  gMinuit->mnexcm("CALL FCN",  &flag1 ,1,error);
  gMinuit->mnexcm("SET PRINT", &flag0 ,1,error);
  gMinuit->mnexcm("SET NOGradient", &flag0 ,1,error);
  gMinuit->mnexcm("MIGRAD",   &flag0 ,0,error);
  gMinuit->mnexcm("MINOS",    &flag0 ,0,error);
  gMinuit->mnexcm("CALL FCN", &flag3 ,1,error);
  
  cout << endl << "--------------------------------------------" << endl;
  cout << "AlManager::FitResolution - INFO: TMinuit minimisation performed." << endl;
  cout << "AlManager::FitResolution - INFO: Asuming Fast_Minimisation." << endl;
  
  //=== store Minuit results into vectors
  vector< double > resV;    resV.clear();
  vector< double > resErrV; resErrV.clear();
  for( unsigned int i=0; i<selDetI_.size(); i++ ) {
    double res, resErr;
    gMinuit->GetParameter(i, res, resErr);
    resV.push_back( res );
    resErrV.push_back( resErr );
  }

  //=== dump to screen
  ofstream out( resFile, ios::app );
  out << "===" << fileSelection_ << endl;
  out << "===" << detFileName << endl;
  out << "===" << string( Utils::GetTimeStamp( ) ) << endl;
  for( unsigned int i=0; i<selDetI_.size(); i++ ) {
    cout.form("%8s %10.5f %10.5f %10.5f %10.5f %10.5f\n",
      selDetI_[i]->TBName_.c_str(), 
      resV[i],  resErrV[i], 
      _ResCalc( i, resV ),
      sqrt( pow( _ResCalc( i, resV ), 2 ) + pow( resV[i], 2 ) ),
      selDetI_[i]->resid_  
    );
    out.form("%8s %10.5f %10.5f %10.5f %10.5f %10.5f\n",
      selDetI_[i]->TBName_.c_str(), 
      resV[i],  resErrV[i], 
      _ResCalc( i, resV ),
      sqrt( pow( _ResCalc( i, resV ), 2 ) + pow( resV[i], 2 ) ),
      selDetI_[i]->resid_  
    );
  }
  
  cout.form("%8s %10s %10s %10s %10s %10s\n\n", "TBName", "res", "error", "tr_res", "residF" ,"residI");
  out.form("%8s %10s %10s %10s %10s %10s\n\n", "TBName", "res", "error", "tr_res", "residF" ,"residI");
  out.close();
  delete gMinuit;
  return true;
}
  
//=======================
//=== Private methods ===
//=======================

//_____________________________________________________________________________
double ResManager::_ResCalc( unsigned int it, vector< double > res )
{
  //=== check vector sizes
  if( res.size() != selDetI_.size() ) { 
    cout << "AlManager::_ResCalc - FATAL: resolution and detInfo sizes does not match.\n";
    return -1;
  }
      
  //=== check <it>
  if( it >= selDetI_.size() ) {
    cout << "AlManager::_ResCalc - FATAL: wrong vector index " << it << ".\n";
    return -1;
  }
  
  TMatrix iV_Fast(2,2); 
  for( unsigned int i=0; i < 2; i++ )   
  for( unsigned int j=0; j < 2; j++ )   
  iV_Fast(i,j) = 0;
   
  for( unsigned int i = 0; i < selDetI_.size(); i++ ) {
    if( it == i ) continue;
    TMatrix iR = selDetI_[i]->detFI_->irotM_;
    double r2 = pow( res[i], 2 );
    iV_Fast(0,0) += iR(0,0)*iR(0,0)/r2;
    iV_Fast(0,1) += iR(0,0)*iR(0,1)/r2;
    iV_Fast(1,1) += iR(0,1)*iR(0,1)/r2;
    iV_Fast(1,0) += iR(0,0)*iR(0,1)/r2;
  }
    
  TMatrix V_Fast = TMatrix( TMatrix::kInverted, iV_Fast );
  TMatrix iR = selDetI_[it]->detFI_->irotM_;
  return sqrt( 
      V_Fast(0,0)*iR(0,0)*iR(0,0)
    + V_Fast(1,1)*iR(0,1)*iR(0,1)
    + 2*V_Fast(0,1)*iR(0,0)*iR(0,1)
  );
  
}    

