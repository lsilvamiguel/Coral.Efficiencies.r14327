// $Id: RTManager.cc,v 1.2 2007/06/01 09:13:52 conrad Exp $
 
/*!
   \file    RTManager.cc
   \brief   RT relation Managment Interface Class.
   \author  Hugo Pereira
   \version $Revision: 1.2 $
   \date    $Date: 2007/06/01 09:13:52 $
*/

#include "RTManager.h"
#include "DetectorInfo.h"
#include "DetFileManager.h"
#include "RTInfo.h"
#include "Fit.h"
#include "TH2Fit.h"
#include "Utils.h"

#include <iostream>

#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TText.h>
#include <TGraphErrors.h>

using namespace std;

ClassImp( RTManager )

//_______________________________________________________________________________
void RTManager::DumpRTInfo( const char* selection )
{
  if( !_Select( selection ) ) return;
  
  for( unsigned int i = 0; i < dets_.size(); i++ ) {
    cout << endl << "RTManager::DumpRTInfo \"" << dets_[i]->TBName_ << "\":" << endl;
    cout << "  T0Parameters: " << ( ( dets_[i]->rtInfo_->HasT0Par() ) ? "yes":"no" ) << endl;
    cout << "  RTParameters: " << ( ( dets_[i]->rtInfo_->HasRTPar() ) ? "yes":"no" ) << endl;
    cout << "  RTGrid      : " << ( ( dets_[i]->rtInfo_->HasRTGrid() ) ? "yes":"no" ) << endl;
  }
  return;
}

//__________________________________________________________________
void RTManager::FitRTStraight( const char* selection, const char* cut, double rMin, double rMax, int nEnt )
{
  if( !_Select( selection ) ) return;

  // handle canvas
  cv_->Clear();
  Utils::DivideCanvas( cv_, dets_.size() );
    
  for( unsigned int i = 0; i < dets_.size(); i++ ) {
    cv_->cd( i+1 );
    TH2* h = (TH2*) _Get( "T_rVect-T_duVect:T_tVect", dets_[i], cut, nEnt );   
    cout << "RTManager::FitRTStraight - \"" << dets_[i]->TBName_ << "\" " << ((h)? (int)h->GetEntries():-1) << " entries.\n";
    _FitRTStraight( h, dets_[i]->rtInfo_, rMin, rMax );
    cv_->Update();
  }
}

//__________________________________________________________________
void RTManager::FitRTGrid( const char* selection, const char* cut, unsigned int nBin, int nEnt )
{
  if( !_Select( selection ) ) return;

  // handle canvas
  cv_->Clear();
  Utils::DivideCanvas( cv_, dets_.size() );
    
  for( unsigned int i = 0; i < dets_.size(); i++ ) {
    cv_->cd( i+1 );

    // check rtInfo
    if( dets_[i]->rtInfo_->t0Par_.size() < 2 ){
      cout << "RTManager::FitRTGrid - ERROR: t0 not set for \"" << dets_[i]->TBName_ << "\". use FitRTStraight first.\n";
      continue;
    }

    TH2* h = (TH2*) _Get( "T_rVect-T_duVect:T_tVect", dets_[i], cut, nEnt );   
    cout << "RTManager::FitRTGrid - \"" << dets_[i]->TBName_ << "\" " << ((h)? (int)h->GetEntries():-1) << " entries.\n";
    _FitRTGrid( h, dets_[i]->rtInfo_, nBin );
     cv_->Update();
 }
}

//__________________________________________________________________
void RTManager::FitRTRelation( const char* selection, const char* cut, int nEnt )
{
  if( !_Select( selection ) ) return;

  // handle canvas
  cv_->Clear();
  Utils::DivideCanvas( cv_, dets_.size() );
    
  for( unsigned int i = 0; i < dets_.size(); i++ ) {
    cv_->cd( i+1 );
    
    // check rtInfo
    if( dets_[i]->rtInfo_->t0Par_.size() < 2 ){
      cout << "RTManager::FitRTRelation - ERROR: t0 not set for \"" << dets_[i]->TBName_ << "\". use FitRTStraight first.\n";
      continue;
    }
    
    TH2* h = (TH2*) _Get( "abs(T_rVect-T_duVect):T_tVect", dets_[i], cut, nEnt );   
    cout << "RTManager::FitRTRelation - \"" << dets_[i]->TBName_ << "\" " << ((h)? (int)h->GetEntries():-1) << " entries.\n";
    _FitRTRelation( h, dets_[i]->rtInfo_ );
    cv_->Update();
  }
}

 
    
//__________________________________________________________________
void RTManager::BuildSingleRT( const char* selection, const char* cut, int nEnt )
{
  if( !_Select( selection ) ) return;

  // check single RT relation histogram/RTInfo object
  if( H_single_RT_ )     { SafeDelete(H_single_RT_);     H_single_RT_ = 0; }
  if( H_single_RT_abs_ ) { SafeDelete(H_single_RT_abs_); H_single_RT_abs_ = 0; }
  if( single_rtInfo_ ) SafeDelete( single_rtInfo_ );
  single_rtInfo_ = new RTInfo( selection );

  cv_->Clear();
  
  bool firstCall = true; 
  for( unsigned int i = 0; i < dets_.size(); i++ ) {
  
    if( firstCall || !(H_single_RT_ && H_single_RT_abs_ ) ) {
      
      // algebric values RT relation
      string name = string(selection)+"_rt_single";    
      string title(name+"::"+string( cut ) );
      if( (H_single_RT_ = (TH2*) gDirectory->Get(name.c_str())) ) delete H_single_RT_;
      H_single_RT_ = (TH2*) _Get( "T_rVect-T_duVect:T_tVect", dets_[i], cut, nEnt );   
      if( !H_single_RT_ ) continue;

      H_single_RT_->SetTitle( title.c_str() );
      H_single_RT_->SetName( name.c_str() );
      cout << "RTManager::BuildSingleRT - \"" << dets_[i]->TBName_ << "\" " << ((H_single_RT_)? (int)H_single_RT_->GetEntries():-1) << " entries.\n";
      
      // absolute values RT relation
      name  = string(selection)+"_rt_single_abs";    
      title = name+"::"+string( cut );
      if( (H_single_RT_abs_ = (TH2*) gDirectory->Get(name.c_str())) ) delete H_single_RT_abs_;
      H_single_RT_abs_ = (TH2*) _Get( "abs(T_rVect-T_duVect):T_tVect", dets_[i], cut, nEnt );   
      if( !H_single_RT_abs_ ) continue;

      H_single_RT_abs_->SetTitle( title.c_str() );
      H_single_RT_abs_->SetName( name.c_str() );
    
      firstCall = false;
    
    } else {
     
      // algebric values RT relation
      TH1* h = _GetCloned( H_single_RT_, "T_rVect-T_duVect:T_tVect", dets_[i], cut, nEnt );
      if( h ) H_single_RT_->Add( h ); 
      cout << "RTManager::BuildSingleRT - \"" << dets_[i]->TBName_ << "\" " << ((h)? (int)h->GetEntries():-1) << " entries.\n";
     
      // absolute values RT relation
      h = _GetCloned( H_single_RT_abs_, "abs(T_rVect-T_duVect):T_tVect", dets_[i], cut, nEnt );
      if( h ) H_single_RT_abs_->Add( h ); 
    
    }
    
    // draw algebric values RT relation
    if( !H_single_RT_ ) continue;
    H_single_RT_->Draw("box");
    cv_->Update();
  }
  return;
}

//__________________________________________________________________
void RTManager::WriteRTGridsToDB( const char* selection, const bool useMyT0, const double myT0 )
{
  if( !_Select( selection ) ) return;
  for( unsigned int i = 0; i < dets_.size(); i++ )
  dets_[i]->rtInfo_->WriteRTGridToDB( useMyT0, myT0 );
  return;
}

//__________________________________________________________________
void RTManager::WriteRTParsToDB( const char* selection, const bool useMyT0, const double myT0 )
{
  if( !_Select( selection ) ) return;
  for( unsigned int i = 0; i < dets_.size(); i++ ) 
  dets_[i]->rtInfo_->WriteRTParToDB( useMyT0, myT0 );
  return;
}

//__________________________________________________________________
void RTManager::WriteSingleRTGridsToDB( const char* selection, const bool useMyT0, const double myT0 )
{
  if( !single_rtInfo_ ) {
    cout << "RTManager::WriteSingleRTGridsToDB - no valid rtInfo object.\n";
    return;
  }
  
  string oldTBName = single_rtInfo_->TBName_;
  string oldValid_start = single_rtInfo_->valid_start_;
  string oldValid_stop = single_rtInfo_->valid_stop_;
  
  
  if( !_Select( selection ) ) return;
  for( unsigned int i = 0; i < dets_.size(); i++ ) {
    single_rtInfo_->TBName_ = dets_[i]->TBName_;
    single_rtInfo_->valid_start_ = dets_[i]->rtInfo_->valid_start_;
    single_rtInfo_->valid_stop_  = dets_[i]->rtInfo_->valid_stop_;
    single_rtInfo_->WriteRTGridToDB( useMyT0, myT0 );
  }
  
  single_rtInfo_->TBName_ = oldTBName;
  single_rtInfo_->valid_start_ = oldValid_start;
  single_rtInfo_->valid_stop_  = oldValid_stop;
  
  return;
}

//__________________________________________________________________
void RTManager::WriteSingleRTParsToDB( const char* selection, const bool useMyT0, const double myT0 )
{
  if( !single_rtInfo_ ) {
    cout << "RTManager::WriteSingleRTGridsToDB - no valid rtInfo object.\n";
    return;
  }
  
  string oldTBName = single_rtInfo_->TBName_;
  string oldValid_start = single_rtInfo_->valid_start_;
  string oldValid_stop = single_rtInfo_->valid_stop_;
  
  
  if( !_Select( selection ) ) return;
  for( unsigned int i = 0; i < dets_.size(); i++ ) {
    single_rtInfo_->TBName_ = dets_[i]->TBName_;
    single_rtInfo_->valid_start_ = dets_[i]->rtInfo_->valid_start_;
    single_rtInfo_->valid_stop_  = dets_[i]->rtInfo_->valid_stop_;
    single_rtInfo_->WriteRTParToDB( useMyT0, myT0 );
  }
  
}

//__________________________________________________________________
void RTManager::SetValidity( const char* selection, const char* start, const char* stop )
{
  if( !_Select( selection ) ) return;
  for( unsigned int i = 0; i < dets_.size(); i++ ) 
  dets_[i]->rtInfo_->SetValidity( start, stop ); 
  return;
}

//__________________________________________________________________
void RTManager::_FitRTStraight( TH2* h, RTInfo* rt, double rMin, double rMax )
{
  if( !(h && h->GetEntries() )  ) { cout << "RTManager::_FitRTStraight - Histogram not valid.\n"; return; }
  if( !rt )                       { cout << "RTManager::_FitRTStraight - RTInfo object not valid.\n"; return; } 
   
  h->Draw("box");
  
  // define function used for the fit, define limits 
  if( rMin > rMax ) {
    rMin = h->GetYaxis()->GetBinCenter( 1 );
    rMax = h->GetYaxis()->GetBinCenter( h->GetNbinsY() );
  }
  
  // define starting parameters 
  TF1* f = new TF1("t0Fit", Fit::t0Fit, rMin, rMax, 3);
  double t0Start = h->GetXaxis()->GetBinCenter( 1 );
  double deltaT =  h->GetXaxis()->GetBinCenter( h->GetNbinsX( ) ) - h->GetXaxis()->GetBinCenter( 1 );
  f->SetParameter( 0, t0Start );
  f->SetParameter( 1, 0 );
  f->SetParameter( 2, deltaT/fabs(rMax) );

  TH2Fit fit( f, 3 );
  
  double p0 = 0;
  fit.ExecMinuitCommand("SET PRINT", &p0, 1); 
  fit.FitInverted( h, rMin, rMax );

  double t0 = f->GetParameter( 0 );
  double v = 1.0/f->GetParameter( 2 );
  TF1* fInv = new TF1("t0FitInv", Fit::t0FitInv, t0, h->GetXaxis()->GetBinCenter( h->GetNbinsX() ), 2);
  fInv->SetParameter( 0, t0 - f->GetParameter( 1 ) * f->GetParameter( 2 ) );
  fInv->SetParameter( 1, v );
  fInv->SetLineColor(4);
  fInv->SetLineWidth(1);
  fInv->Draw("same");

  TF1* fInvNeg = new TF1("t0FitInv", Fit::t0FitInv, t0, h->GetXaxis()->GetBinCenter( h->GetNbinsX() ), 2);
  fInvNeg->SetParameter( 0, t0 + f->GetParameter( 1 ) * f->GetParameter( 2 ) );
  fInvNeg->SetParameter( 1, -v );
  fInvNeg->SetLineColor(4);
  fInvNeg->SetLineWidth(1);
  fInvNeg->Draw("same");

  char* text = new char[100];
  sprintf(text,"t0=%9.1f", f->GetParameter( 0 ));
  TText* tt = new TText( );
  tt->SetTextColor( 4 ); 
  tt->SetTextSize(0.06);
  tt->DrawTextNDC( 0.15, 0.15, text );
  SafeDelete( text );

  // save result 
  rt->t0Par_.clear();
  for( unsigned int iP = 0; iP < 3; iP++ ) rt->t0Par_.push_back( f->GetParameter( iP ) );   
  
  return;
}


//__________________________________________________________________
void RTManager::_FitRTGrid( TH2* h, RTInfo* rt, unsigned int nBins )
{
  if( !(h && h->GetEntries() )  ) { cout << "RTManager::_FitRTGrid - Histogram not valid.\n"; return; }
  if( !rt )                       { cout << "RTManager::_FitRTGrid - RTInfo object not valid.\n"; return; } 
     
  // Check t0Pars_
  if( rt->t0Par_.size() < 2 ){
    cout << "RTManager::_FitRTGrid - ERROR: t0 not set for \"" << rt->TBName_ << "\". use FitRTStraight first.\n";
    return;
  }

  h->Draw("box");
  if( !nBins ) nBins = h->GetNbinsX();
  
  cout << "_FitRTGrid - bins " << nBins << endl;
  
  // initialize RTGrid 
  rt->rtGrid_.clear();
  cout << "_FitRTGrid - grid init ok.\n";

  // initialize Fit Function 
  unsigned int nR = h->GetNbinsY();
  TF1* f = new TF1("gridFit", Fit::gridFit,
      h->GetYaxis()->GetBinCenter( 1 ),
      h->GetYaxis()->GetBinCenter( nR ), 4 );
  cout << "_FitRTGrid - TF1 init ok.\n";
   
  // Scan bins 
  unsigned int nTG = 0;
  unsigned int nT = h->GetXaxis()->GetNbins();
  unsigned int n0 = h->GetXaxis()->FindBin( rt->t0Par_[0] ); 
  
  for( unsigned int iB = 0; iB < nBins; iB++ ) {
    #ifdef  START_AT_T0
    unsigned int i1 = (unsigned int) ( (double) iB*(nT-n0)/nBins + n0)+1;
    unsigned int i2 = (unsigned int) ( (double) (iB+1)*(nT-n0)/nBins + n0);
    #else
    unsigned int i1 = (unsigned int) ( (double) iB*nT/nBins)+1;
    unsigned int i2 = (unsigned int) ( (double) (iB+1)*nT/nBins);
    #endif

    TH1* slice = h->ProjectionY("slice", i1, i2 );
    if( ! (slice && slice->GetEntries() ) ) { 
      SafeDelete(slice); 
      continue;
    }
    
    f->SetParameter( 0, slice->GetMaximum() ); 
    f->SetParameter( 1, fabs(slice->GetBinCenter( slice->GetMaximumBin() ) ) ); 
    f->SetParameter( 2, rt->t0Par_[1] );      
    f->SetParameter( 3, 0.3 ); 
      
    slice->Fit(f,"0");   
    double r = fabs( f->GetParameter( 1 ) ) + f->GetParameter( 2 );
    double t = h->GetBinCenter( (i1+i2)/2 );
    double res = fabs( f->GetParameter( 3 )/sqrt( 0.5*slice->GetEntries() ) );
    rt->_AddRTGridPoint( t, r, res );
    
    SafeDelete(slice);
  }
  
  TGraphErrors* tGE = new TGraphErrors();
  tGE->SetMarkerColor(2);
  tGE->SetLineColor(2);
  tGE->SetLineWidth(1);
  tGE->SetMarkerStyle(8);
  tGE->SetMarkerSize(2.9);
  
  for( unsigned int iGr = 0; iGr < rt->rtGrid_.size(); iGr++ ) {
    tGE->SetPoint( iGr, rt->rtGrid_[iGr].t, rt->rtGrid_[iGr].r );
    tGE->SetPointError( iGr, 0, rt->rtGrid_[iGr].res ); 
  }
  
  tGE->Draw("same");
    
  return;
}  

//__________________________________________________________________
void RTManager::_FitRTRelation( TH2 *h, RTInfo *rt )
{
  if( !(h && h->GetEntries() )  ) { cout << "RTManager::_FitRTRelation - Histogram not valid.\n"; return; }
  if( !rt )                       { cout << "RTManager::_FitRTRelation - RTInfo object not valid.\n"; return; } 
    
  // Check t0Pars_
  if( rt->t0Par_.size() < 2 ){
    cout << "RTManager::_FitRTRelation - ERROR: t0 not set for \"" << rt->TBName_ << "\". use FitRTStraight first.\n";
    return;
  }

  h->Draw("box");
   
  unsigned int nP = 5;      
  double p0 = 0, p1 = 1, p2 = 2;
  
  double xMin = h->GetXaxis()->GetBinCenter( 1 );
  double xMax = h->GetXaxis()->GetBinCenter( h->GetNbinsX() );
  double yMax = h->GetYaxis()->GetBinCenter( h->GetNbinsY() );
  TF1* f = new TF1("rtFit", Fit::rtFit , xMin, xMax, nP);
  f->SetLineColor(4);
  f->SetLineWidth(1);
  for( unsigned int iP = 0; iP < nP; iP++ ) f->SetParameter( iP, 0 );
  f->SetParameter( 0, rt->t0Par_[0]-rt->t0Par_[1]*rt->t0Par_[2] ); 
  f->SetParameter( 1, 1.0/rt->t0Par_[2] ); 
  
  TH2Fit fit( f, nP );
  fit.ExecMinuitCommand("FIX", &p1, 1); 
  fit.ExecMinuitCommand("SET PRINT", &p0, 1); 
  fit.Fit( h, xMin, xMax );
  f->Draw("same");
  
  // Update rtInfo
  rt->rtPar_.clear();
  for( unsigned int iP = 0; iP < nP; iP++ ) rt->rtPar_.push_back( f->GetParameter( iP ) );

  return;
}

     
    
