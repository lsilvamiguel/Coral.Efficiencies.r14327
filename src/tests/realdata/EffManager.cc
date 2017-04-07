// $Id: EffManager.cc,v 1.2 2007/06/01 09:13:51 conrad Exp $

/*!
   \file    EffManager.cc
   \brief   Efficiency Managment Interface Class.
   \author  Hugo Pereira
   \version $Revision: 1.2 $
   \date    $Date: 2007/06/01 09:13:51 $
*/

#include "EffManager.h"
#include "DetectorInfo.h"
#include "DetFileManager.h"
#include "Utils.h"

#include <TText.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>

#include <fstream>

#define effFile "eff.out"

using namespace std;

ClassImp( EffManager )

//_______________________________________________________________________________
void EffManager::DrawEff( const char* varSel, const char* selection, TCut cut, const char* opt, int nEnt )
{
  if( !_Select( selection ) ) return;

  char* text = new char[100];

  cout << "EffManager::DrawEff - INFO: Efficiency cut is \"" << cut_ << "\".\n";

  // check if 1D or 2D histos are required
  bool is2D = ( string( varSel ).find(":") != string::npos );
  string stOpt( opt );
  if( is2D ) cout << "EffMAnager::DrawEff - INFO: 2D histos required.\n";
    
  // Update Option in case of 2D Display
  if( stOpt == "" && is2D ) stOpt = string("colz");
  
  // handle canvas
  cv_->Clear();
  Utils::DivideCanvas( cv_, dets_.size() );
  
  // dump header to file
  ofstream out( effFile, ios::app );
  out << "===";
  for( unsigned int i = 0; i < fileSelection_.size(); i++ ) out << " " << fileSelection_[i];
  out << endl;
  out << "===" << string( Utils::GetTimeStamp( ) ) << endl;
  out << "===" << varSel << ":" << cut << endl;
  
    
  for( unsigned int iInf = 0; iInf < dets_.size(); iInf++ ) {
    cv_->cd( iInf+1 );
    
    // get/check reference histogram
    TH1 *hr = _Get( varSel, dets_[iInf], cut, nEnt );
    
    if( !(hr && hr->GetEntries() ) ) {
      cout << "EffManager::DrawEff - \"" << dets_[iInf]->TBName_ << "\". ref histogram empty.\n";
      continue;
    } 
    
    // add efficiency cut to selection cut
    TCut effCut = TCut( cut ) && cut_;

    // get detector histogram 
    TH1 *hd = _GetCloned( hr, varSel, dets_[iInf], (const char*) effCut, nEnt );
    int nE = (int) hd->GetEntries();
        
    // build Efficiency histogram
    string name = string( dets_[iInf]->TBName_ + "_eff");
    string title(name+"("+varSel+")"+"::"+string( cut ) );
    TH1* eff;
    if( (eff = (TH1*) gDirectory->Get(name.c_str())) ) SafeDelete(eff);
    
    eff = ( ( is2D ) ? 
      (TH1*) new TH2D( name.c_str(), title.c_str(), 
        hr->GetNbinsX(), 
        hr->GetXaxis()->GetXmin(),
        hr->GetXaxis()->GetXmax(),
        hr->GetNbinsY(), 
        hr->GetYaxis()->GetXmin(),
        hr->GetYaxis()->GetXmax() ) :
      (TH1*) new TH1D( name.c_str(), title.c_str(), 
        hr->GetNbinsX(), 
        hr->GetXaxis()->GetXmin(),
        hr->GetXaxis()->GetXmax() ) 
    );
    
    double e = ( ( is2D ) ?
      Utils::HDiv( (TH2*) hd, (TH2*) hr, (TH2*) eff ) :
      Utils::HDiv( (TH1*) hd, (TH1*) hr, (TH1*) eff )
    );
    eff->SetEntries(hr->GetEntries() );
    if( ! is2D ) {
      eff->SetMarkerStyle(20);
      eff->SetLineColor( 4 );
      eff->SetMarkerColor(4);
      eff->SetMarkerSize(0.7);
    } else eff->SetMinimum(0.5*e);
    
    eff->Draw(stOpt.c_str() );

    sprintf(text,"%6.1f",100.0*e);
    TText* tt = new TText( eff->GetXaxis()->GetBinCenter( 3 ),0.2, text );
    tt->SetTextColor( 1 ); 
    tt->SetTextSize(0.1);
    tt->Draw();
    if( cv_ ) cv_->Update();
    
    // Dump Entries and Efficiency
    cout << "EffManager::DrawEff - "<< dets_[iInf]->TBName_ << ": " << nE << " entries.\n";
    sprintf(text,"%8s %10i %10.3f\n", dets_[iInf]->TBName_.c_str(), nE, e );
    out<<text;
    
    if( hr ) SafeDelete( hr );
    if( hd ) SafeDelete( hd );
    
  }
  sprintf(text,"%8s %10s %10s\n\n", "", "Entries", "eff" );
  out<<text;
  SafeDelete(text);
  out.close();
} 

//_______________________________________________________________________________
void EffManager::DrawEffUProf( const char* selection, TCut cut, const char* opt, int nEnt )
{
  if( !_Select( selection ) ) return;
  
  // handle canvas
  cv_->Clear();
  Utils::DivideCanvas( cv_, dets_.size() );

  for( unsigned int iInf = 0; iInf < dets_.size(); iInf++ ) {
    cv_->cd( iInf+1 );
    
    // get/check reference histogram
    TH1 *hr = _Get( "T_uLx", dets_[iInf], cut, nEnt );
    if( !(hr && hr->GetEntries() ) ) {
      cout << "EffManager::DrawEff - \"" << dets_[iInf]->TBName_ << "\". ref histogram empty.\n";
      continue;
    }
    
    // add efficiency cut to selection cut
    TCut effCut = TCut( cut ) && cut_;

    // get detector histogram 
    TH1 *hd = _GetCloned( hr, "T_uLx", dets_[iInf], (const char*) effCut, nEnt );
    int nE = (int) hd->GetEntries();
    
    // ref profile        
    hr->SetLineColor(4); 
    hr->Draw();          
    
    // det profile
    hd->SetLineColor(2);          
    hd->Draw("same");
    
    double e = double( hd->GetEntries() )/ hr->GetEntries();
    char text[100];   
    sprintf(text,"eff %6.1f %%",100.0*e);
    TText* tt = new TText( -500, 0.1*hr->GetMaximum(), text );
    tt->SetTextColor( 4 ); 
    tt->SetTextSize(0.06);
    tt->Draw();
    
    if( cv_ ) cv_->Update();
   
    // Dump Entries and Efficiency
    cout << "EffManager::DrawEffUProf - "<< dets_[iInf]->TBName_ << ": " << nE << " entries.\n";
         
  }
}  
