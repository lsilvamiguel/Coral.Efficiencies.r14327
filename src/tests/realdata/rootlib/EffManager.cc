// $Id: EffManager.cc,v 1.29 2002/11/06 14:02:45 hpereira Exp $
#include "EffManager.h"
#include "DetInfo.h"
#include "Utils.h"

#include <fstream.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TText.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>

#define effFile "eff.out"

ClassImp( EffManager )

//_______________________________________________________________________________
EffManager::EffManager( void ):HManager(), cut_("T_fnd") { }
EffManager::EffManager( const char* fileSelection, bool batch = false ):
  HManager( fileSelection, batch ), 
  cut_("T_fnd") {}
EffManager::EffManager( const HManager h ): HManager( h ), cut_("T_fnd") {}

//_______________________________________________________________________________
void EffManager::DrawEff( const char* varSel, const char* selection = "*", TCut cut="1", const char* opt="" )
{
  cout << "EffManager::DrawEff - INFO: Efficiency cut is \"" << cut_ << "\".\n";
  
  //=== check if 1D or 2D histos are required
  bool is2D = ( string( varSel ).find(":") != string::npos );
  string stOpt( opt );
  if( is2D ) cout << "EffMAnager::DrawEff - INFO: 2D histos required.\n";
    
  //=== Update Option in case of 2D Display
  if( stOpt == "" && is2D ) stOpt = string("colz");
  
  //=== Get Selection
  GetSelection( selection );
  vector< DetInfo*> DetOK;
  DetOK.clear();
  for( unsigned int iInf = 0; iInf < selDetI_.size(); iInf++ )
  if( selDetI_[iInf]->T_eff_ ) 
  DetOK.push_back( selDetI_[iInf] );
  
  if( DetOK.empty() ) {
    cout << "EffManager::DrawEff - WARNING: works only with Efficiency tree.\n";
    return;
  }
  
  if( tcv_ )   tcv_->Clear();
  _DivideCanvas( DetOK.size() );
  
  ofstream out( effFile, ios::app );
  out << "===" << fileSelection_ << endl;
  out << "===" << string( Utils::GetTimeStamp( ) ) << endl;
  out << "===" << varSel << ":" << cut << endl;
    
  for( unsigned int iInf = 0; iInf < DetOK.size(); iInf++ ) {
    if( tcv_ ) tcv_->cd( iInf+1 );
    
    TH1 *hr = 0, *hd = 0;
    int nE = 0;
    
    string nameRef = DetOK[iInf]->TBName_+"_ref";
    string nameDet = DetOK[iInf]->TBName_+"_det";
    string varRef = string( varSel ) + ">>" + nameRef;
    if( (hr = (TH1*) gDirectory->Get(nameRef.c_str())) ) delete hr;
    if( (hd = (TH1*) gDirectory->Get(nameDet.c_str())) ) delete hd;
    nE = DetOK[iInf]->T_eff_->Draw(varRef.c_str(), cut,"goff" );
    
    if( is2D ) {
      hr = (TH2S*) gDirectory->Get(nameRef.c_str());
      
      if( !(hr && hr->GetEntries() ) ) {
        cout << "EffManager::DrawEff - \"" << DetOK[iInf]->TBName_ << "\". ref histogram empty.\n";
        continue;
      }
      
      hd = (TH2S*) hr->Clone(nameDet.c_str()); 
      DetOK[iInf]->T_eff_->Project(nameDet.c_str(), varSel,cut_&&cut );
    } else {
      hr = (TH1S*) gDirectory->Get(nameRef.c_str());
      
      if( !(hr && hr->GetEntries() ) ) {
        cout << "EffManager::DrawEff - \"" << DetOK[iInf]->TBName_ << "\". ref histogram empty.\n";
        continue;
      }
      
      hd = (TH1S*) hr->Clone(nameDet.c_str()); 
      DetOK[iInf]->T_eff_->Project(nameDet.c_str(), varSel,cut_&&cut );
    }
    
    string name = string( DetOK[iInf]->TBName_ + "_eff");
    string title(name+"("+varSel+")"+"::"+string( cut ) );
    TH1* eff;
    
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
      Utils::HDiv( (TH2*) hd, (TH2*) hr, (TH2D*) eff ) :
      Utils::HDiv( (TH1*) hd, (TH1*) hr, (TH1D*) eff )
    );

    eff->SetEntries(hr->GetEntries() );
    if( ! is2D ) {
      eff->SetMarkerStyle(20);
      eff->SetLineWidth( 2 );
      eff->SetLineColor( 4 );
      eff->SetMarkerColor(4);
      eff->SetMarkerSize(0.7);
    } else eff->SetMinimum(0.5*e);
    
    eff->Draw(stOpt.c_str() );

    char text[100];   
    sprintf(text,"%6.1f",100.0*e);
    TText* tt = new TText( eff->GetXaxis()->GetBinCenter( 3 ),0.2, text );
    tt->SetTextColor( 1 ); 
    tt->SetTextSize(0.1);
    tt->Draw();
    if( tcv_ ) tcv_->Update();
    
    // Dump Entries and Efficiency
    cout << "EffManager::DrawEff - "<< DetOK[iInf]->TBName_ << ": " << nE << " entries.\n";
    out.form("%8s %10i %10.3f\n", DetOK[iInf]->TBName_.c_str(), nE, e );

    // Cleaning
    if( DetOK[iInf]->T_eff_ ) {
      if( hd ) delete hd;
      if( hr ) delete hr;
    }
  }
  out.form("%8s %10s %10s\n\n", "", "Entries", "eff" );
  out.close();
} 

//_______________________________________________________________________________
void EffManager::DrawEffUProf( const char* selection = "*", TCut cut="T_inActive", const char* opt="" )
{
  cout << "EffManager::DrawEffUProf - INFO: Efficiency cut is \"" << cut_ << "\".\n";
  vector< DetInfo*> DetOK = GetSelection( selection );
  if( DetOK.empty() ) return;
  
  if( tcv_ )   tcv_->Clear();
  _DivideCanvas( DetOK.size() );

  for( unsigned int iInf = 0; iInf < DetOK.size(); iInf++ ) {
    if( tcv_ ) tcv_->cd( iInf+1 );
        
    TH1S *hd = 0, *hr = 0;
    int nE = 0;

    //=== Use Tree
    //=== cleaning (Root HORROR)
    string nameRef = DetOK[iInf]->TBName_+"_ref";
    string nameDet = DetOK[iInf]->TBName_+"_det";
    string varSel( "T_uLx" );
    string varRef = varSel+">>"+nameRef;
    if( (hr = (TH1S*) gDirectory->Get( nameRef.c_str() ) ) ) delete hr;
    if( (hd = (TH1S*) gDirectory->Get( nameDet.c_str() ) ) ) delete hd;
    
    //=== fill histos
    nE = DetOK[iInf]->T_eff_->Draw( varRef.c_str(), cut,"goff" );
    hr = (TH1S*) gDirectory->Get( nameRef.c_str() );
    
    if( !(hr && hr->GetEntries() ) ) {
      cout << "EffManager::DrawEffUProf - \"" << DetOK[iInf]->TBName_ << "\". ref histogram empty.\n";
      continue;
    }
    
    hd = (TH1S*) hr->Clone( nameDet.c_str() ); 
    DetOK[iInf]->T_eff_->Project( nameDet.c_str() , varSel.c_str(), cut_&&cut );
    
            
    hr->SetLineWidth(3); hd->SetLineWidth(3); 
    hr->SetLineColor(4); hd->SetLineColor(2);
    hr->Draw();          hd->Draw("same");
          
    double e = double( hd->GetEntries() )/ hr->GetEntries();
    char text[100];   
    sprintf(text,"eff %6.1f %%",100.0*e);
    TText* tt = new TText( -500, 0.1*hr->GetMaximum(), text );
    tt->SetTextColor( 4 ); 
    tt->SetTextSize(0.06);
    tt->Draw();
    
    if( tcv_ ) tcv_->Update();
   
    // Dump Entries and Efficiency
    cout << "EffManager::DrawAllEffUProf - "<< DetOK[iInf]->TBName_ << ": " << nE << " entries.\n";
         
  }
}  
