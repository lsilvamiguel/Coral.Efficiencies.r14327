// $Id: RTManager.cc,v 1.26 2002/11/06 14:03:09 hpereira Exp $
#include "RTManager.h"
#include "DetInfo.h"
#include "DetFileInfo.h"
#include "RTInfo.h"
#include "CsRTGridPoint.h"
#include "Fit.h"
#include "TH2Fitter.h"

#include <iostream.h>

#include <TFile.h>
#include <TList.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>


ClassImp( RTManager )

//_______________________________________________________________________________
RTManager::RTManager( void ): 
  HManager(), 
  H_single_RT_(0), 
  single_rtInfo_( new RTInfo( "Single" ) ) 
{ }

//_______________________________________________________________________________
RTManager::RTManager( const char* fileSelection, bool batch = false ): 
  HManager( fileSelection, batch ), 
  H_single_RT_(0), 
  single_rtInfo_( new RTInfo( "Single" ) ) { }

//_______________________________________________________________________________
RTManager::RTManager( const HManager h ): 
  HManager( h ), 
  H_single_RT_(0), 
  single_rtInfo_( new RTInfo( "Single" ) ) 
{ }

//===================================
//=== ShortCuts to DetInfo Methods ===
//===================================

//_______________________________________________________________________________
void RTManager::DumpRTInfo( const char* selection = "*" ) 
{
  GetSelection( selection );
  
  for( unsigned int i = 0; i < selDetI_.size(); i++ ) 
  if( selDetI_[i]->HasTBName() ) {  
    cout << endl << " DetInfo \"" << selDetI_[i]->TBName_ << "\":" << endl;
    cout << "  T0Parameters: " << ( ( selDetI_[i]->rtInfo_->HasT0Par() ) ? "yes":"no" ) << endl;
    cout << "  RTParameters: " << ( ( selDetI_[i]->rtInfo_->HasRTPar() ) ? "yes":"no" ) << endl;
    cout << "  RTGrid      : " << ( ( selDetI_[i]->rtInfo_->HasRTGrid() ) ? "yes":"no" ) << endl;
  }
  return;
}

//__________________________________________________________________
void RTManager::DrawT0( const char* selection = "*", TCut cut = "T_fnd" )
{
  GetSelection( selection );

  if( tcv_ != 0 ) tcv_->Clear();
  _DivideCanvas( selDetI_.size() );
  
  TH1S* h;
  for( unsigned int iInf = 0; iInf < selDetI_.size(); iInf++ ) 
  {
    if( tcv_ != 0 ) tcv_->cd( iInf+1 );
    if( !selDetI_[iInf]->T_eff_ ) continue;
    string name = selDetI_[iInf]->TBName_+"_T_tMin_";
    string title = name+":"+string( cut );
    string select = string( "T_tMin>>" )+name;
    if( (h = (TH1S*) gDirectory->Get( name.c_str() )) ) delete h;
    selDetI_[iInf]->T_eff_->Draw( select.c_str(), cut, "goff" );
    h = (TH1S*) gDirectory->Get( name.c_str() );
    
    if( !(h && h->GetEntries() ) ) {
      cout << "RTManager::DrawT0 - \"" << selDetI_[iInf]->TBName_ << "\". histogram empty.\n";
      continue;
    }

    int min = int(h->GetXaxis()->GetXmin());
    int max = int(h->GetXaxis()->GetXmax());
    int bin = max - min;
    h->Rebin( (int)( (h->GetNbinsX()/bin) + 1 ) ); 
    h->Draw();
    if( tcv_ != 0 ) tcv_->Update();
  }
  return;
}

//__________________________________________________________________
void RTManager::FitT0( const char* selection, TCut cut, double tMin, double tMax )
{
  GetSelection( selection );

  if( tcv_ != 0 ) tcv_->Clear();
  _DivideCanvas( selDetI_.size() );
  
  for( unsigned int iInf = 0; iInf < selDetI_.size(); iInf++ ) 
  {
    if( tcv_ != 0 ) tcv_->cd( iInf+1 );
    selDetI_[iInf]->FitT0( cut, tMin, tMax );
    if( tcv_ != 0 ) tcv_->Update();
  }
  return;
}

//__________________________________________________________________
void RTManager::DrawSingleRT( const char* selection = "*", TCut cut = "T_fnd" )
{
  if( H_single_RT_ ) delete H_single_RT_;
  
  GetSelection( selection );

  if( !selDetI_.size() ) return;
  if( tcv_ != 0 ) tcv_->Clear();
  string name(string(selection)+"_T_r_vs_"+( (UseCorrectedVals_) ? "tCor":"tMin") );    
  
  bool firstCall = true; 
  for( unsigned int i = 0; i < selDetI_.size(); i++ ) 
  {
    if( !selDetI_[i]->T_eff_ ) continue;
    if( firstCall || !H_single_RT_ ) {
      firstCall = false;
      string title(name+"::"+string( cut ) );    
      string select = ( (UseCorrectedVals_) ? string("T_rLx:T_tCor>>"):string("T_rLx:T_tMin>>") )+name;
      if( (H_single_RT_ = (TH2S*) gDirectory->Get(name.c_str())) ) delete H_single_RT_;
      selDetI_[i]->T_eff_->Draw(select.c_str(), cut, "goff" );
      H_single_RT_ = (TH2S*) gDirectory->Get(name.c_str());

      if( !H_single_RT_ ) {
        cout << "RTManager::DrawSingleRT - \"" << selDetI_[i]->TBName_ << "\".histogram empty.\n";
        continue;
      }

      H_single_RT_->SetTitle( title.c_str() );
    } else {
      string select = ( (UseCorrectedVals_) ? string("T_rLx:T_tCor>>+"):string("T_rLx:T_tMin>>+") )+name;
      selDetI_[i]->T_eff_->Draw(select.c_str(), cut, "goff" );
    }
    if( !H_single_RT_ ) continue;
    H_single_RT_->Draw("box");
    if( tcv_ != 0 ) tcv_->Update();
  }
  return;
}

//__________________________________________________________________
void RTManager::FitSingleRTStraight( double rMin = 1, double rMax = -1  )
{
  TH2S* h = 0;
  if( !(h=H_single_RT_) ) {
    cout << "RTManager::FitSingleRTStraight - ERROR: call RTManager::DrawSingleRT first.\n";
    return;
  }
   
  h->Draw("box");
  
  //=== define function used for the fit and limits ===
  if( rMin > rMax ) {
    cout << "DetInfo::FitRTStraight - INFO: Fit range is histogram range." << endl;
    rMin = H_single_RT_->GetYaxis()->GetBinCenter( 1 );
    rMax = H_single_RT_->GetYaxis()->GetBinCenter( h->GetNbinsY() );
  }
  
  TF1* f = new TF1("t0Fit", Fit::t0Fit, rMin, rMax, 3);
  
  //=== define starting parameters ===
  double t0Start = h->GetXaxis()->GetBinCenter( 1 );
  double deltaT =  h->GetXaxis()->GetBinCenter( h->GetNbinsX( ) ) - h->GetXaxis()->GetBinCenter( 1 );
  f->SetParameter( 0, t0Start );
  f->SetParameter( 1, 0 );
  f->SetParameter( 2, deltaT/fabs(rMax) );

  TH2Fitter::Instance()->SetFunction( f, 3 );
  TH2Fitter::Instance()->FitInverted( H_single_RT_, rMin, rMax );

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

  char text[100];
  sprintf(text,"t0=%9.1f", f->GetParameter( 0 ));
  TText* tt = new TText( t0,h->GetYaxis()->GetBinCenter( 3 ), text );
  tt->SetTextColor( 2 ); 
  tt->SetTextSize(0.06);
  tt->Draw();


  //=== save result ===
  single_rtInfo_->t0Par_.clear();
  for( unsigned int iP = 0; iP < 3; iP++ ) single_rtInfo_->t0Par_.push_back( f->GetParameter( iP ) );   
  
  return;
}

//__________________________________________________________________
void RTManager::FitSingleRTGrid( unsigned int nBins = 0 )
{
  TH2S* h = 0;
  if( !(h=H_single_RT_) ) {
    cout << "RTManager::FitSingleRTGrid - ERROR: call RTManager::DrawSingleRT first.\n";
    return;
  }
  
  //=== Check t0Pars_
  if( single_rtInfo_->t0Par_.size() < 2 ){
    cout << "RTManager::FitSingleRTGrid - ERROR: use FitSinleRTStraight first." << endl;
    return;
  }

  h->Draw("box");
  if( !nBins ) nBins = h->GetNbinsX();
  
  //=== initialize RTGrid ===
  single_rtInfo_->rtGrid_.clear();
  
  //=== initialize Fit Function ===
  unsigned int nR = h->GetNbinsY();
  TF1* f = new TF1("gridFit", Fit::gridFit,
      h->GetYaxis()->GetBinCenter( 1 ),
      h->GetYaxis()->GetBinCenter( nR ), 4 );
    
  //=== Scan bins ===
  unsigned int nTG = 0;
  unsigned int nT = h->GetXaxis()->GetNbins();
  unsigned int n0 = h->GetXaxis()->FindBin( single_rtInfo_->t0Par_[0] ); 
  
  for( unsigned int iB = 0; iB < nBins; iB++ ) {
    #ifdef  START_AT_T0
    unsigned int i1 = (unsigned int) ( (double) iB*(nT-n0)/nBins + n0)+1;
    unsigned int i2 = (unsigned int) ( (double) (iB+1)*(nT-n0)/nBins + n0);
    #else
    unsigned int i1 = (unsigned int) ( (double) iB*nT/nBins)+1;
    unsigned int i2 = (unsigned int) ( (double) (iB+1)*nT/nBins);
    #endif

    TH1D* slice = h->ProjectionY("slice", i1, i2 );
    if( ! slice->GetEntries() ) { 
      delete slice; 
      continue;
    }
    
    f->SetParameter( 0, slice->GetMaximum() ); 
    f->SetParameter( 1, fabs(slice->GetBinCenter( slice->GetMaximumBin() ) ) ); 
    f->SetParameter( 2, single_rtInfo_->t0Par_[1] );      
    f->SetParameter( 3, 0.3 ); 
      
    slice->Fit(f,"0");   
    double r = fabs( f->GetParameter( 1 ) ) + f->GetParameter( 2 );
    double t = h->GetBinCenter( (i1+i2)/2 );
    double res = fabs( f->GetParameter( 3 )/sqrt( 0.5*slice->GetEntries() ) );
    single_rtInfo_->_AddRTGridPoint( t, r, res );
    
    delete slice;
  }
  
  TGraphErrors* tGE = new TGraphErrors();
  tGE->SetMarkerColor(2);
  tGE->SetLineColor(2);
  tGE->SetLineWidth(1);
  tGE->SetMarkerStyle(8);
  tGE->SetMarkerSize(2.9);
  
  for( unsigned int iGr = 0; iGr < single_rtInfo_->rtGrid_.size(); iGr++ ) {
    tGE->SetPoint( iGr, single_rtInfo_->rtGrid_[iGr].t, single_rtInfo_->rtGrid_[iGr].r );
    tGE->SetPointError( iGr, 0, single_rtInfo_->rtGrid_[iGr].res ); 
  }
  
  tGE->Draw("same");
    
  return;
}  
//__________________________________________________________________
void RTManager::FitRTStraights( const char* selection = "*", TCut cut = "T_fnd", double rMin = 1, double rMax = -1 )
{
  GetSelection( selection );

  if( tcv_ != 0 ) tcv_->Clear();
  _DivideCanvas( selDetI_.size() );
  
  for( unsigned int iInf = 0; iInf < selDetI_.size(); iInf++ ) 
  {
    if( tcv_ != 0 ) tcv_->cd( iInf+1 );
    selDetI_[iInf]->FitRTStraight( cut, rMin, rMax );
    if( tcv_ != 0 ) tcv_->Update();
    
    _WriteT0ToLog( selDetI_[iInf] );
  }
  return;
}

//__________________________________________________________________
void RTManager::FitRTRelations( const char* selection = "*", TCut cut = "T_fnd", const char* detFile = "" )
{
  GetSelection( selection );

  if( tcv_ != 0 ) tcv_->Clear();
  _DivideCanvas( selDetI_.size() );
  
  for( unsigned int iInf = 0; iInf < selDetI_.size(); iInf++ ) 
  {
    if( tcv_ != 0 ) tcv_->cd( iInf+1 );
    if( strlen( detFile ) == 0 || selDetI_[iInf]->SetDetFile( detFile ) ){
      selDetI_[iInf]->FitRTRelation( cut );
      if( tcv_ != 0 ) tcv_->Update();
    }
  }
  return;
}

//__________________________________________________________________
void RTManager::FitRTGrids( const char* selection = "*", TCut cut = "T_fnd", unsigned int nBin = 0  )
{
  GetSelection( selection );

  if( tcv_ != 0 ) tcv_->Clear();
  _DivideCanvas( selDetI_.size() );
  
  for( unsigned int iInf = 0; iInf < selDetI_.size(); iInf++ ) 
  {
    if( tcv_ != 0 ) tcv_->cd( iInf+1 );
    selDetI_[iInf]->FitRTGrid( cut, nBin );
    if( tcv_ != 0 ) tcv_->Update();
  }
  
  return;
}
 
//================================
//=== ShortCuts to I/O Methods ===
//================================

//__________________________________________________________________
void RTManager::WriteRTGridsToDB( const char* selection = "*", 
  const bool useMyT0 = false, 
  const double myT0 = 0 )
{
  GetSelection( selection );

  for( unsigned int iInf = 0; iInf < selDetI_.size(); iInf++ )
  selDetI_[iInf]->rtInfo_->WriteRTGridToDB( useMyT0, myT0 );

}

//__________________________________________________________________
void RTManager::WriteRTParsToDB( const char* selection = "*", 
  const bool useMyT0 = false, 
  const double myT0 = 0 )
{
  GetSelection( selection );

  for( unsigned int iInf = 0; iInf < selDetI_.size(); iInf++ ) 
  selDetI_[iInf]->rtInfo_->WriteRTParToDB( useMyT0, myT0 );

}

//__________________________________________________________________
void RTManager::SetValidity( const char* selection = "*", const char* start = "YYYY-MM-DD-hh:mm:ss", const char* stop = "YYYY-MM-DD-hh:mm:ss" ) 
{
  GetSelection( selection );

  for( unsigned int iInf = 0; iInf < selDetI_.size(); iInf++ ) 
  selDetI_[iInf]->rtInfo_->SetValidity( start, stop ); 
  return;
}
  
//__________________________________________
void RTManager::_WriteT0ToLog( DetInfo* detI )
{
  FILE* out;
  static bool firstCall = true;

  // check detInfo has T0 parameters
  if( !detI->rtInfo_->HasT0Par() ) {
    cout << "RTManager::_WriteT0ToLog - ERROR: DetInfo \"" << detI->TBName_ << "\" has no t0 parameters.\n";
    return;
  }    

  // open file
  if( ( out=fopen( "./t0.log", "a" ) ) == 0 ) {
    cout << "RTManager::_WriteT0ToLog - ERROR: cannot write to file \"./t0.log\"\n";
    return;
  } else if( firstCall ) {
    cout << "RTManager::_WriteT0ToLog - INFO: log file for T0 is \"./t0.log\"\n";
    fprintf( out, "\n");
    firstCall = false;
  }    
  
  fprintf( out, "%s t0Pars ", detI->TBName_.c_str() );
  for( unsigned int ip = 0; ip < detI->rtInfo_->t0Par_.size(); ip++ ) 
  fprintf( out, "%12.3f", detI->rtInfo_->t0Par_[ip] );
  fprintf( out, "\n");
  fclose( out );
}

//=== Check resolution using new RT Relation
//_________________________________________________________
void RTManager::CheckResolutions( const char* selection = "*", TCut cut = "T_fnd", unsigned int nEvt = 1000000000 )
{
  //=== Address for needed branches
  unsigned int T_evt;          // Event number
  double T_duMin;              // track dU/dZ
  double T_tMin;               // track dU/dZ
  double T_rLx;                // distance between helix and closest wire
  
  GetSelection( selection );

  if( tcv_ != 0 ) tcv_->Clear();
  _DivideCanvas( selDetI_.size() );
  
  for( unsigned int iInf = 0; iInf < selDetI_.size(); iInf++ ) 
  {
    DetInfo *det = selDetI_[iInf];
    RTInfo *rt = det->rtInfo_;
    TTree* T_eff_ = det->T_eff_->GetTree();
    
    //=== Check tree
    if( !T_eff_) continue;
    
    //=== Check RT Relation
    if( !( rt->HasRTPar() || rt->HasRTGrid() ) ) continue;
    
    //=== Get Reference Resolution Histogram
    if( tcv_ != 0 ) tcv_->cd( iInf+1 );
    TH1* hRef = 0;
    string name(det->TBName_+"_T_duMin");    
    string title(name+"::"+string( cut ) );
    string select("T_duMin"); select += ">>"+name;
    if( (hRef = (TH1S*) gDirectory->Get(name.c_str())) ) delete hRef;
    T_eff_->Draw(select.c_str(), cut, "goff", nEvt );
    hRef = (TH1S*) gDirectory->Get(name.c_str());
    
    if( !(hRef && hRef->GetEntries() ) ) {
      cout << "RTManager::CheckResolutions - \"" << selDetI_[iInf]->TBName_ << "\". histogram empty.\n";
      continue;
    }

    hRef->SetTitle( title.c_str() );
    hRef->Draw();
  
    //=== Create new histogram
    name = string( det->TBName_+"_T_duMin_new");   
    TH1* hNew = (TH1* ) hRef->Clone( name.c_str() );
    hNew->Reset("ICE");
    hNew->SetTitle( title.c_str() );

    double duMin = hRef->GetXaxis()->GetXmin();
    double duMax = hRef->GetXaxis()->GetXmax(); 
    
    //=== Parse Tree
    T_eff_->SetBranchAddress("T_evt",    &T_evt);
    T_eff_->SetBranchAddress("T_duMin",  &T_duMin);
    T_eff_->SetBranchAddress("T_tMin",   &T_tMin);
    T_eff_->SetBranchAddress("T_rLx",    &T_rLx);
    
    for( unsigned i=0; i<T_eff_->GetEntries() && i < nEvt; i++ ) {

      //=== printout
      if( i % 100 == 0 ) {
        printf("\rEvent %i",i); 
  	    fflush(stdout);
      }

      if( ! T_eff_->Draw("T_evt",cut, "goff", 1, i ) ) continue;
      if( ! T_eff_->GetEntry( i ) ) break;     
      
      bool error;
      double r = rt->GetRfromT( T_tMin, error );
      if( det->hasDetFile_ && r > det->detFI_->wirP_*0.5 ) r = det->detFI_->wirP_*0.5;
      double du = r - T_rLx;
//       cout.form( "Event %6i:T_duMin=%10.4f du=%10.4f [%10.4f-%10.4f]... \n",
//         i,
//         T_duMin,
//         du,
//         duMin, duMax );   
      if( du > duMin && du < duMax ) hNew->Fill( du );
      
    }
    
    hNew->SetLineColor( 2 );
    hNew->Draw("same"); 
     
  }
   
  return;
} 
    
    
