// $Id: DetInfo.cc,v 1.78 2002/11/06 14:02:12 hpereira Exp $
#include "DetInfo.h"
#include "Fit.h"

#include <TText.h>
#include "DetFileInfo.h"
#include "RTInfo.h"
#include "CsRTGridPoint.h"
#include "Utils.h"
#include "TH2Fitter.h"
ClassImp( DetInfo )

//_________________________________________________________
DetInfo::DetInfo( string TBName = "", string detFile = "" ): 
  TObject(),
  TBName_( TBName ),
  OneGaus_( false ),
  UseCorrectedVals_( false ),
  CorrectFromAngle_( false ),
  resid_( -1 ),
  
  //=== RTRelation
  rtInfo_( new RTInfo( TBName ) ),
  
  //=== Tree for everything
  T_eff_( 0 ),

  //=== Geometry
  hasDetFile_( false ),
  detFI_( 0 ),
  detFile_( "" )
{
  if( !TBName.size() ) return;
  if( detFile.size() ) {
    detFile_    =   string(detFile);
    hasDetFile_ =  _ReadDetFile();
  }
}

//_________________________________________________________
DetInfo::DetInfo( int ID, string detFile ):    
  TObject(),
  TBName_( "" ),
  OneGaus_( false ),
  UseCorrectedVals_( false ),
  CorrectFromAngle_( false ),
  resid_( -1 ),

  //=== Tree for everything
  T_eff_( 0 ),

  //=== Geometry
  hasDetFile_( false ),
  detFI_( 0 ),
  detFile_( "" )
{

  //=== Check detFile
  if( !detFile.size() ) { 
    cout << "DetInfo::DetInfo - FATAL: detFile not specified.\n";
    return;
  }
  
  //=== Try loading TBName
  detFile_ = string(detFile);
  if( !_GetTBNameFromDetFile( ID ) ) {
    cout << "DetInfo::DetInfo - FATAL: no match for detID " << ID << " in \"" 
      << detFile << "\".\n";
    return;
  }

  //=== Get DetFileInfos
  hasDetFile_ =  _ReadDetFile();
  
  //=== Create SubObjects
  rtInfo_ = new RTInfo( TBName_ );
  
} 


//_________________________________________________________________________________________
bool DetInfo::SetDetFile( const char* detFile ) 
{
  if( detFile_ == string( detFile ) && hasDetFile_ ) {
    cout << "DetInfo::SetDetFile - Info: \"" << detFile << "\" already loaded." << endl;
    return true;
  }
  
  detFile_ = string( detFile );
  if( ! (hasDetFile_ = _ReadDetFile() ) ) {
    cout << "DetInfo::SetDetFile - Error: bad file <= \"" << detFile << "\"." << endl;
    detFile_ = ("");
    return false;
  }
  
  return true;
}

//_________________________________________________________________________________________
void DetInfo::AddToChain( const char* treeName, const char* fileName )
{
  if( !T_eff_ ) T_eff_ = new TChain( treeName, treeName );
  T_eff_->Add( fileName, -1  );
}

//========================
//=== RTRelation Managment
//______________________________________________________________________________
bool DetInfo::FitT0( TCut cut, double tMin, double tMax )
{
  TH1S *h;

  if( !T_eff_ ) {
    cout << "DetInfo::FitT0 - ERROR: This method works only with tree." << endl;
    return false;
  }

  if( tMin > tMax ) {
    cout << "DetInfo::FitT0 - ERROR: Wrong limits.\n";
    return false;
  }
    
  string select = (UseCorrectedVals_) ? string("T_tCor"):string("T_tMin");
  string name(TBName_+"_"+select);    
  string title(name+"::"+string( cut ) );    
  select += string(">>")+name;
  if( (h = (TH1S*) gDirectory->Get(name.c_str())) ) delete h;
  T_eff_->Draw(select.c_str(), cut, "goff" );
  h = (TH1S*) gDirectory->Get(name.c_str());
  
  if( !( h && h->GetEntries() ) ) {
    cout << "DetInfo::FitT0 - ERROR: empty histogram.\n";
    return false;
  }
  
  h->SetTitle( title.c_str() );
  int min = int (h->GetXaxis()->GetXmin());
  int max = int (h->GetXaxis()->GetXmax());
  int bin = max - min; 
  h->Rebin( (int)( (h->GetNbinsX()/bin) + 1 ) ); 
  
  if( tMin < min ) tMin = min;
  if( tMax > max ) tMax = max;
  
  h->Draw();

  //=== define function used for the fit and limits ===
  TF1* f = new TF1("t0Fit", Fit::driftTFit, tMin, tMax, 4);
  double p0 = double( 0.5*h->GetMaximum() ) ;
  double p1 = double( h->GetMaximum() - p0 );
  double p2 = double( 0.5*(tMax+tMin) );
  double p3 = double( tMax - tMin );
  
  f->SetParameter( 0, p0 );
  
  const unsigned int NTry = 3;
  for( unsigned int t=0; t<NTry; t++ ) {
  
    f->SetParameter( 1, p1 );
    f->SetParameter( 2, p2 );
    f->SetParameter( 3, p3 );
    f->SetLineWidth( 1 );
    f->SetLineColor(4);
 
    // Do the fit
    h->Fit(f,"0","", tMin, tMax );
    
    if( f->GetParameter( 3 )/f->GetParameter(1) >= 0.001 ) break;
    p1*=1.2;
  }
 
  f->Draw("same");
  
  char text[100];
  double t0 = f->GetParameter(2);
  sprintf(text,"%6.1fns",t0);
  TText* tt = new TText( h->GetXaxis()->GetBinCenter( bin/2 ),0.2*h->GetMaximum(), text );
  tt->SetTextColor( 4 ); 
  tt->SetTextSize(0.06);
  tt->Draw();
  
  return true;
}  

//______________________________________________________________________________
bool DetInfo::FitRTStraight( 
  TCut cut="T_fnd", 
  double rMin = 1, double rMax = -1 )
{
  TH2S* h;
  
  if( !T_eff_ )  return false;  
  string name(TBName_+"_T_r_vs_"+( (UseCorrectedVals_) ? "tCor":"tMin") );    
  string title(name+"::"+string( cut ) );    
  string select = ( (UseCorrectedVals_) ? string("T_rLx:T_tCor>>"):string("T_rLx:T_tMin>>") )+name;
  if( (h = (TH2S*) gDirectory->Get(name.c_str())) ) delete h;
  T_eff_->Draw(select.c_str(), cut, "goff" );
  h = (TH2S*) gDirectory->Get(name.c_str());
  
  if( !( h && h->GetEntries() ) ) {
    cout << "DetInfo::FitRTStraight - ERROR: empty histogram.\n";
    return false;
  }
  
  h->SetTitle( title.c_str() );
  h->Draw("");
  
  //=== define function used for the fit and limits ===
  if( rMin > rMax ) {
    cout << "DetInfo::FitRTStraight - INFO: Fit range is histogram range." << endl;
    rMin = h->GetYaxis()->GetBinCenter( 1 );
    rMax = h->GetYaxis()->GetBinCenter( h->GetNbinsY() );
  }
  
  TF1* f = new TF1("t0Fit", Fit::t0Fit, rMin, rMax, 3);
  
  //=== define starting parameters ===
  double t0Start = h->GetXaxis()->GetBinCenter( 1 );
  double deltaT =  h->GetXaxis()->GetBinCenter( h->GetNbinsX( ) ) - h->GetXaxis()->GetBinCenter( 1 );
  f->SetParameter( 0, t0Start );
  f->SetParameter( 1, 0 );
  f->SetParameter( 2, deltaT/fabs(rMax) );

  TH2Fitter::Instance()->SetFunction( f, 3 );
  TH2Fitter::Instance()->FitInverted( h, rMin, rMax );
  
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
  rtInfo_->t0Par_.clear();
  for( unsigned int iP = 0; iP < 3; iP++ ) rtInfo_->t0Par_.push_back( f->GetParameter( iP ) );   
  
  return true;
}

//_______________________________________________________________________________
bool DetInfo::FitRTRelation( TCut cut="T_fnd" )
{
  //=== Check t0Pars_
  if( rtInfo_->t0Par_.size() < 2 ){
    cout << "DetInfo::FitRTRelation - ERROR: use FitRTStraight first." << endl;
    return false;
  }

  #ifndef _FIT1_
  //=== check detFile is present ===
  if( ! hasDetFile_ ) {
    cout << "DetInfo::FitRTRelation - ERROR: detFile needed." << endl ;
    return false;
  }
  #endif

  
  TH2S* h;
  
  if( !T_eff_ ) return false;
  string name(TBName_+"_T_ra_vs_"+( (UseCorrectedVals_) ? "tCor":"tMin") );    
  string title(name+"::"+string( cut ));    
  string select = ( (UseCorrectedVals_) ? string("abs(T_rLx):T_tCor>>"):string("abs(T_rLx):T_tMin>>") )+name;
  if( (h = (TH2S*) gDirectory->Get(name.c_str())) ) delete h;
  T_eff_->Draw(select.c_str(), cut, "goff" );
  h = (TH2S*) gDirectory->Get(name.c_str());
  
  if( !( h && h->GetEntries() ) ) {
    cout << "DetInfo::FitRTRelation - ERROR: empty histogram.\n";
    return false;
  }
 
  h->SetTitle( title.c_str() );
  h->Draw("box");
  
  //=== define function used for the fit and limits ===

  #ifdef _FIT1_
  const unsigned int nP = 5;    

  #else
  const unsigned int nP = 5;      
  #endif

  double xMin = h->GetXaxis()->GetBinCenter( 1 );
  double xMax = h->GetXaxis()->GetBinCenter( h->GetNbinsX() );
  double yMax = h->GetYaxis()->GetBinCenter( h->GetNbinsY() );
  TF1* f = new TF1("rtFit", Fit::rtFit , xMin, xMax, nP);
  f->SetLineColor(4);
  f->SetLineWidth(1);
  
  //=== define starting parameters and initialize TH2Fitter 
  double p1 = 1, p2 = 2;

  for( unsigned int iP = 0; iP < nP; iP++ ) f->SetParameter( iP, 0 );

  #ifdef _FIT1_ 
  // polynomial fit. 
  // First parameter is t0, kept constant
  // Second parameter is drift velocity around 0
  f->SetParameter( 0, rtInfo_->t0Par_[0]-rtInfo_->t0Par_[1]*rtInfo_->t0Par_[2] ); 
  f->SetParameter( 1, 1.0/rtInfo_->t0Par_[2] ); 
  TH2Fitter::Instance()->SetFunction( f, nP );
  TH2Fitter::Instance()->ExecMinuitCommand("FIX", &p1, 1); //Fix parameter 0

  #else
  // hyperbolic tangent fit.
  // First parameter is detector pitch, keps constant 
  // Second parameter is t0, kept constant
  // Third parameter is drift velocity around 0
  f->SetParameter( 0, detFI_->wirP_ ); 
  f->SetParameter( 1, rtInfo_->t0Par_[0]-rtInfo_->t0Par_[1]*rtInfo_->t0Par_[2] ); 
  f->SetParameter( 2, 1.0/rtInfo_->t0Par_[2] ); 
  TH2Fitter::Instance()->SetFunction( f, nP );
  TH2Fitter::Instance()->ExecMinuitCommand("FIX", &p1, 1); //Fix parameter 0
  TH2Fitter::Instance()->ExecMinuitCommand("FIX", &p2, 1); //Fix parameter 1
  #endif

  TH2Fitter::Instance()->Fit( h, xMin, xMax );
  f->Draw("same");
  
  //=== save DetInfo ===
  rtInfo_->rtPar_.clear();
  for( unsigned int iP = 0; iP < nP; iP++ ) rtInfo_->rtPar_.push_back( f->GetParameter( iP ) );
  
  return true;
} 

//_______________________________________________________________________________
bool DetInfo::FitRTGrid( TCut cut="T_fnd", unsigned int nBins = 0 )
{

  //=== Check t0Pars_
  if( rtInfo_->t0Par_.size() < 2 ){
    cout << "DetInfo::FitRTGrid - ERROR: use FitRTStraight first." << endl;
    return false;
  }
  
  TH2S* h;
  
  if( !T_eff_ ) return false;
  string name(TBName_+"_T_r_vs_"+( (UseCorrectedVals_) ? "tCor":"tMin") );    
  string title(name+"::"+string( cut ));    
  string select = ( (UseCorrectedVals_) ? string("T_rLx:T_tCor>>"):string("T_rLx:T_tMin>>") )+name;
  if( (h = (TH2S*) gDirectory->Get(name.c_str())) ) delete h;
  T_eff_->Draw(select.c_str(), cut, "goff" );
  h = (TH2S*) gDirectory->Get(name.c_str());

  if( !( h && h->GetEntries() ) ) {
    cout << "DetInfo::FitRTGrid - ERROR: empty histogram.\n";
    return false;
  }

  h->SetTitle( title.c_str() );
  h->Draw("");  
  if( !nBins ) nBins = h->GetNbinsX();
  
  //=== initialize RTGrid ===
  rtInfo_->rtGrid_.clear();
  
  //=== initialize Fit Function ===
  unsigned int nR = h->GetNbinsY();
  TF1* f = new TF1("gridFit", Fit::gridFit,
      h->GetYaxis()->GetBinCenter( 1 ),
      h->GetYaxis()->GetBinCenter( nR ), 4 );
    
  //=== Scan bins ===
  unsigned int nTG = 0;
  unsigned int nT = h->GetXaxis()->GetNbins();
  unsigned int n0 = h->GetXaxis()->FindBin( rtInfo_->t0Par_[0] ); 

//  #define START_AT_T0
  #ifdef  START_AT_T0
  //=== Add first grid point by hand
  rtInfo_->_AddRTGridPoint( rtInfo_->t0Par_[0], 0, 0 );
  #endif
  
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
    f->SetParameter( 2, rtInfo_->t0Par_[1] );      
    f->SetParameter( 3, 0.3 ); 
      
    slice->Fit(f,"0");   
    double r = fabs( f->GetParameter( 1 ) ) + f->GetParameter( 2 );
    double t = h->GetBinCenter( (i1+i2)/2 );
    double res = fabs( f->GetParameter( 3 )/sqrt( 0.5*slice->GetEntries() ) );
    rtInfo_->_AddRTGridPoint( t, r, res );
    
    delete slice;
  }
  
  TGraphErrors* tGE = new TGraphErrors();
  tGE->SetMarkerColor(2);
  tGE->SetLineColor(2);
  tGE->SetLineWidth(1);
  tGE->SetMarkerStyle(8);
  tGE->SetMarkerSize(2.9);

  for( unsigned int iGr = 0; iGr < rtInfo_->rtGrid_.size(); iGr++ ) {
    tGE->SetPoint( iGr, rtInfo_->rtGrid_[iGr].t, rtInfo_->rtGrid_[iGr].r );
    tGE->SetPointError( iGr, 0, rtInfo_->rtGrid_[iGr].res ); 
  }
  
  tGE->Draw("same");
    
  return true;
}   

//_______________________________________________________________________________
TF1* DetInfo::FitTvsV( TCut cut="T_fnd", double p0 = 0, double p1 = 0 )
{
    
  TH2S* h;

  if( !T_eff_ ) return 0;

  //=== use tree
  string name(TBName_+"_" + ( (UseCorrectedVals_) ?"T_TCor":"T_TMin" ) +"_vs_v");    
  string title(name+"::"+string( cut ) );
  string select( ( (UseCorrectedVals_) ? "T_tCor:T_vLx>>":"T_tMin:T_vLx>>" ) + name);
  if( (h = (TH2S*) gDirectory->Get(name.c_str())) ) delete h;
  T_eff_->Draw(select.c_str(), cut, "goff" );
  h = (TH2S*) gDirectory->Get(name.c_str());
  
  if( !( h && h->GetEntries() ) ) {
    cout << "DetInfo::FitTvsV - ERROR: empty histogram.\n";
    return 0;
  }
   
  h->Draw();
  TProfile* p = (TProfile*) h->ProfileX();
  p->SetErrorOption("s");
  p->SetLineColor( 1 );
  p->SetLineWidth( 1 );
  p->Draw("same");
    
  //=== define function used for the fit and limits ===
  const unsigned int nP = 2;
  double xMin = h->GetXaxis()->GetBinCenter( 1 );
  double xMax = h->GetXaxis()->GetBinCenter( h->GetNbinsX() );
  TF1* f = new TF1("du2DFit", Fit::du2DFit , xMin, xMax, nP);
  f->SetLineColor(2);
  f->SetLineWidth(1);
  
  //=== define starting parameters  
  for( unsigned int iP = 0; iP < nP; iP++ ) f->SetParameter( iP, 0 );
  f->SetParameter( 0, p0 );
  f->SetParameter( 1, p1 );


  TH2Fitter::Instance()->SetFunction( f, nP );
  TH2Fitter::Instance()->Fit( h, xMin, xMax );

   f->Draw("same");
  char text[100];   
  sprintf(text,"t=%g+V*%g",f->GetParameter(0), f->GetParameter(1));
  TText* tt = new TText( h->GetXaxis()->GetBinCenter( 3 ),5, text );
  tt->SetTextColor( 4 ); 
  tt->SetTextSize(0.06);
  tt->Draw();  
  
  return f;
} 

//=======================
//=== Residuals Managment
//_______________________________________________________________________________
bool DetInfo::FitDU( TCut cut="T_fnd" )
{
    
  TH1S* h;
  double du;
  
  if( !T_eff_ ) return false;
  string name(TBName_+"_" + ( (UseCorrectedVals_) ?"T_duCor":"T_duMin") );    
  string title(name+"::"+string( cut ) );
  string select;
  switch( UseCorrectedVals_ + (CorrectFromAngle_<<1) ) {
  case 0: select = string("T_duMin"); break;     // (!UseCorrectedVals_)&&(!CorrectFromAngle_)
  case 1: select = string("T_duCor"); break;     // (UseCorrectedVals_)&&(!CorrectFromAngle_)
  case 2: select = string("T_uLx+cos(atan(T_tuLx))*(T_duMin-T_uLx)"); break;     // (!UseCorrectedVals_)&&(CorrectFromAngle_)
  case 3: select = string("T_uLx+cos(atan(T_tuLx))*(T_duCor-T_uLx)"); break;     // (UseCorrectedVals_)&&(CorrectFromAngle_)
  default:
    cout << "DetInfo::FitDU - ERROR: switches not understood." << endl;
    return false;
  }
  
  cout << "DetInfo::FitDU \"" << TBName_ << "\" - INFO: Selection is \""<< select << "\".";
  select += ">>"+name;
  if( (h = (TH1S*) gDirectory->Get(name.c_str())) ) delete h;
  T_eff_->Draw(select.c_str(), cut, "goff" );
  h = (TH1S*) gDirectory->Get(name.c_str());

  if( !( h && h->GetEntries() ) ) {
    cout << "empty histogram.\n";
    return false;
  }
 
  h->SetTitle( title.c_str() );
  cout << " Entries: " << h->GetEntries() << ".";
  h->Draw();
  
  if( OneGaus_ ) {

    h->Fit( "gaus", "0Q" );     
    h->GetFunction("gaus")->SetLineColor(4);
    h->GetFunction("gaus")->SetLineWidth(1);
    h->GetFunction("gaus")->Draw("same");

    //=== store result for alignment
    resid_ = fabs(h->GetFunction("gaus")->GetParameter(2));
    du =    h->GetFunction("gaus")->GetParameter(1);
    
    char text[100];
    sprintf(text,"%6.1fum",1000.0*resid_);
    TText* tt = new TText( h->GetXaxis()->GetBinCenter( 3 ),0.2*h->GetMaximum(), text );
    tt->SetTextColor( 4 ); 
    tt->SetTextSize(0.06);
    tt->Draw();
          
  } else {

    //=== define function used for the fit and limits ===
    const unsigned int nP = 6;
    double xMin = h->GetXaxis()->GetBinCenter( 1 );
    double xMax = h->GetXaxis()->GetBinCenter( h->GetNbinsX() );
    TF1* f = new TF1("du1DFit", Fit::du1DFit , xMin, xMax, nP);
    f->SetLineColor(4);
    f->SetLineWidth(1);
    
    //=== define starting parameters  
    // first gaussian
    f->SetParameter( 0, h->GetMaximum() );
    f->SetParameter( 1, h->GetBinCenter( h->GetMaximumBin() ) );
    f->SetParameter( 2, 0.5*h->GetRMS( ) );
  
    // second gaussian
    f->SetParameter( 3, 0.2*h->GetMaximum() );
    f->SetParameter( 4, h->GetBinCenter( h->GetMaximumBin() ) );
    f->SetParameter( 5, 2*h->GetRMS( ) );
    
    // do the fit
    h->Fit( f, "0" ); 
    
    f->Draw("same");
  
    // Draw The largest gaussian separately
    TF1* fg = new TF1("gaus", (const char*) "gaus(0)", xMin, xMax );
    fg->SetLineColor(2);
    fg->SetLineWidth(1);
    if( f->GetParameter( 0 ) < f->GetParameter( 3 ) ) {
      fg->SetParameter( 0, f->GetParameter( 0 ) );
      fg->SetParameter( 1, f->GetParameter( 1 ) );
      fg->SetParameter( 2, f->GetParameter( 2 ) );
    } else {
      fg->SetParameter( 0, f->GetParameter( 3 ) );
      fg->SetParameter( 1, f->GetParameter( 4 ) );
      fg->SetParameter( 2, f->GetParameter( 5 ) );
    }
    fg->Draw("same");

    //=== store result for alignment
    resid_ = fabs( ( f->GetParameter(0) > f->GetParameter(3) )? f->GetParameter(2):f->GetParameter(5) );
    du =    ( f->GetParameter(0) > f->GetParameter(3) )? f->GetParameter(1):f->GetParameter(4);
    
    char text[100];
    sprintf(text,"%6.1fum",1000.0*resid_);
    TText* tt = new TText( h->GetXaxis()->GetBinCenter( 3 ),0.2*h->GetMaximum(), text );
    tt->SetTextColor( 4 ); 
    tt->SetTextSize(0.06);
    tt->Draw();
  }
  
  cout << " res=" << resid_ << "mm. ave=" << du << "mm.\n";
  
  if( hasDetFile_ ) 
  for( unsigned int i=0; i < detFList_.size(); i++ ) 
  detFList_[i]->DUtoCenter_CM( -du );
  
  return true;
} 

//_______________________________________________________________________________
bool DetInfo::FitDU_MWPC( TCut cut="T_fnd" )
{
    
  TH1S* h;
  
  if( !T_eff_ ) return false;
  string name(TBName_+"_T_du");    
  string title(name+"::"+string( cut ) );
  string select("T_duMin");
  select += ">>"+name;
  if( (h = (TH1S*) gDirectory->Get(name.c_str())) ) delete h;
  T_eff_->Draw(select.c_str(), cut, "goff" );
  h = (TH1S*) gDirectory->Get(name.c_str());

  if( !( h && h->GetEntries() ) ) {
    cout << "DetInfo::FitDU_MWPC - ERROR: empty histogram.\n";
    return false;
  }

  h->SetTitle( title.c_str() );
  cout << "DetInfo::FitDU_MWPC \"" << TBName_ << "\" - INFO: Entries: " << h->GetEntries() << endl;

  h->Draw();
  
  //=== define function used for the fit and limits ===
  const unsigned int nP = 7;
  double xMin = h->GetXaxis()->GetBinCenter( 1 );
  double xMax = h->GetXaxis()->GetBinCenter( h->GetNbinsX() );
  TF1* f = new TF1("du1DFit_MWPC", Fit::du1DFit_MWPC , xMin, xMax, nP);
  f->SetLineColor(4);
  f->SetLineWidth(1);
  
  //=== define starting parameters  
  f->SetParameter( 0, h->GetMaximum() );
  f->SetParameter( 1, h->GetBinCenter(h->GetMaximumBin() ) );
  f->SetParameter( 2, 0.5*h->GetRMS( ) );
  f->SetParameter( 3, 1.0 );

  f->SetParameter( 4, 0.2*h->GetMaximum() );
  f->SetParameter( 5, h->GetBinCenter(h->GetMaximumBin() ) );
  f->SetParameter( 6, 2*h->GetRMS() );

  h->Fit( f, "0" ); 
  f->Draw("same");

  //=== store result for alignment
  resid_ = f->GetParameter(2);
  double du =    f->GetParameter(1);
  
  if( hasDetFile_ ) {
    for( unsigned int i=0; i < detFList_.size(); i++ ) 
    detFList_[i]->DUtoCenter_CM( -du );
  }
    
    
  return true;
} 

//_______________________________________________________________________________
bool DetInfo::FitDUvsU( TCut cut="T_fnd" )
{
    
  TH2S* h;

  if( !T_eff_ ) return false;
  string name(TBName_+"_" + ( (UseCorrectedVals_) ?"T_duCor":"T_duMin" ) +"_vs_u");    
  string title(name+"::"+string( cut ) );
  string select( ( (UseCorrectedVals_) ? "T_duCor:T_uLx>>":"T_duMin:T_uLx>>" ) + name);
 
  if( (h = (TH2S*) gDirectory->Get( name.c_str() ) ) ) delete h;
  T_eff_->Draw(select.c_str(), cut, "goff" );
  h = (TH2S*) gDirectory->Get(name.c_str());

  if( !( h && h->GetEntries() ) ) {
    cout << "DetInfo::FitDUvsU - ERROR: empty histogram.\n";
    return false;
  }

  cout << "DetInfo::FitDUvsU \"" << TBName_ << "\" - INFO: Entries: " << h->GetEntries();
  cout << ". Select is " << select << endl;
 
  h->Draw();
  
  //=== define function used for the fit and limits ===
  const unsigned int nP = 2;
  double xMin = h->GetXaxis()->GetBinCenter( 1 );
  double xMax = h->GetXaxis()->GetBinCenter( h->GetNbinsX() );
  TF1* f = new TF1("du2DFit", Fit::du2DFit , xMin, xMax, nP);
  f->SetLineColor(2);
  f->SetLineWidth(1);
  
  //=== define starting parameters  
  for( unsigned int iP = 0; iP < nP; iP++ ) f->SetParameter( iP, 0 );

  TH2Fitter::Instance()->SetFunction( f, nP );
  TH2Fitter::Instance()->Fit( h, xMin, xMax );

  f->Draw("same");
  
  return true;
} 

//_______________________________________________________________________________
bool DetInfo::FitDUvsV( TCut cut="T_fnd" )
{
    
  TH2S* h;

  if( !T_eff_ ) return false;
  string name(TBName_+"_" + ( (UseCorrectedVals_) ?"T_duCor":"T_duMin" ) +"_vs_v");    
  string title(name+"::"+string( cut ) );
  string select( ( (UseCorrectedVals_) ? "T_duCor:T_vLx>>":"T_duMin:T_vLx>>" ) + name);
  if( (h = (TH2S*) gDirectory->Get(name.c_str())) ) delete h;
  T_eff_->Draw(select.c_str(), cut, "goff" );
  h = (TH2S*) gDirectory->Get(name.c_str());

  if( !( h && h->GetEntries() ) ) {
    cout << "DetInfo::FitDUvsV - ERROR: empty histogram.\n";
    return false;
  }
  
  h->Draw();
  cout << "DetInfo::FitDUvsV \"" << TBName_ << "\" - INFO: Entries: " << h->GetEntries();
  cout << ". Select is " << select << endl;
  
  //=== define function used for the fit and limits ===
  const unsigned int nP = 2;
  double xMin = h->GetXaxis()->GetBinCenter( 1 );
  double xMax = h->GetXaxis()->GetBinCenter( h->GetNbinsX() );
  TF1* f = new TF1("du2DFit", Fit::du2DFit , xMin, xMax, nP);
  f->SetLineColor(2);
  f->SetLineWidth(1);
  
  //=== define starting parameters  
  for( unsigned int iP = 0; iP < nP; iP++ ) f->SetParameter( iP, 0 );

  TH2Fitter::Instance()->SetFunction( f, nP );
  TH2Fitter::Instance()->Fit( h, xMin, xMax );

  f->Draw("same");
  
  return true;
} 


//_______________________________________________________________________________
bool DetInfo::DrawDU_LR( TCut cut="T_fnd", double OffsetN = 0, double OffsetP = 0 )
{
    
  TH1S *hL, *hR;
  
  if( !T_eff_ )return false;    
  string name(TBName_+"_" + ( (UseCorrectedVals_) ?"T_duCor":"du_TMin" ));    
 
  string title(name+"::"+string( cut ) );
  string select;
  switch( UseCorrectedVals_ + (CorrectFromAngle_<<1) ) {
  case 0: select = string("T_duMin"); break;     // (!UseCorrectedVals_)&&(!CorrectFromAngle_)
  case 1: select = string("T_duCor"); break;     // (UseCorrectedVals_)&&(!CorrectFromAngle_)
  case 2: select = string("T_uLx+cos(atan(T_tuLx))*(T_duMin-T_uLx)"); break;     // (!UseCorrectedVals_)&&(CorrectFromAngle_)
  case 3: select = string("T_uLx+cos(atan(T_tuLx))*(T_duCor-T_uLx)"); break;     // (UseCorrectedVals_)&&(CorrectFromAngle_)
  default:
    cout << "DetInfo::FitDU - \"" << TBName_ << "\" ERROR: switches not understood." << endl;
    return false;
  }
    
  cout << "DetInfo::FitDU - \"" << TBName_ << "\" INFO: Selection is \""<< select << "\". ";
    
  string nameR(name+"R_");     
  char* selectR = new char[256];
  sprintf( selectR, "%s+%f>>%s",select.c_str(), OffsetP, nameR.c_str());
  if( (hR = (TH1S*) gDirectory->Get(nameR.c_str())) ) delete hR;
  T_eff_->Draw(selectR, cut+"T_rLx>0", "goff" );
  hR = (TH1S*) gDirectory->Get(nameR.c_str());

  if( !hR ) {
    cout << "DetInfo::FitDU_LR - ERROR: R histogram not found.\n";
    return false;
  }

  hR->SetTitle( title.c_str() );
  hR->SetLineWidth(1);
  hR->SetLineColor(2);
  delete selectR;
  cout << " Entries (R): " << hR->GetEntries();

  string nameL(name+"L_");    
  char* selectL = new char[256];
  sprintf( selectL, "%s+%f>>%s",select.c_str(), OffsetN, nameL.c_str());
  if( (hL = (TH1S*) gDirectory->Get(nameL.c_str())) ) delete hL;
  T_eff_->Draw(selectL, cut+"T_rLx<0", "goff" );
  hL = (TH1S*) gDirectory->Get(nameL.c_str());

  if( !hL ) {
    cout << "DetInfo::FitDU_LR - ERROR: L histogram not found.\n";
    return false;
  }

  hL->SetTitle( title.c_str() );
  hL->SetLineWidth(1);
  hL->SetLineColor(4);
  delete selectL;
  cout << ". Entries (L): " << hL->GetEntries() << endl;

  //=== Draw histograms
  if( hR->GetMaximum() > hL->GetMaximum() ) {
    hR->Draw();
    hL->Draw("same"); 
  } else {
    hL->Draw();
    hR->Draw("same"); 
  }

  //=== add some text
  char text[100];
  double mean = hR->GetMean();
  sprintf(text,"R>0 %6.1fum",1000*mean);
  TText* tt = new TText( hR->GetXaxis()->GetBinCenter( 3 ),0.4*hR->GetMaximum(), text );
  tt->SetTextColor( 2 ); 
  tt->SetTextSize(0.06);
  tt->Draw();
     
  mean = hL->GetMean();
  sprintf(text,"R<0 %6.1fum",1000*mean);
  tt = new TText( hL->GetXaxis()->GetBinCenter( 3 ),0.2*hL->GetMaximum(), text );
  tt->SetTextColor( 4 ); 
  tt->SetTextSize(0.06);
  tt->Draw();
  
  return true;
}

//_____________________________________________________________________________
bool DetInfo::DrawWRSProfile( TCut cut="" )
{
    
  TH2S* h;

  if( !T_eff_ ) return false;
  string name(TBName_+"_T_u_vs_v_");
  string title(name+"::"+string( cut ) );
  string select("T_vLx:T_uLx");

  T_eff_->Draw(select.c_str(), cut, "goff" );
  h = (TH2S*) gDirectory->Get(name.c_str());

  if( !(h&&h->GetEntries() ) ) {
    cout << "DetInfo::DrawWRSProfile - ERROR: histogram is empty.\n";
    return false;
  }

  h->Draw("color");
  return true;
}

//=======================
//=== Private Methods ===
//=======================

//_______________________________________________________________________________
void DetInfo::_AddDetFileInfo( DetFileInfo* dfi)
{
  if( !dfi ) return;
  detFList_.push_back( dfi );
  if( (!detFI_) || pow( detFI_->xcm_,2 ) + pow( dfi->ycm_,2 ) > pow( detFI_->xcm_,2 ) + pow( dfi->ycm_,2 )
    || ( TBName_.substr(0,2) == "ST" && TBName_[7] == 'b' ) )
  detFI_ = dfi ;
}

//_________________________________________________________________________________________
bool DetInfo::_ReadDetFile( void )
{
  if( ! HasTBName() ) {
    cout << "DetInfo::_GetDetFileInfo - ERROR: TBName not set." << endl;
    return false;
  }  
  
  ifstream in( detFile_.c_str(), ios::in );
  if( !in ) {
    cout << "DetInfo::_GetDetFileInfo - ERROR: cannot read file \"" 
         << detFile_  << "\"." 
         << endl;
    return false;
  }
    
  char line[512];
  unsigned int lineNumber = 0;
  
	while( !in.eof() ) { 
    in.getline( line, sizeof( line ), '\n');
    if( in.eof() || !in.good() ) continue;
    
    if( line[0] == '\0' || line[0] != ' ' ) continue;    
    if( line[1] == ' ' || line[1] == '-' || line[1] == '=' ) continue;
    lineNumber++;      
    
    istrstream s(line);
    string opt;
    s >> opt;

		if ( opt != "det" ) continue;
    
    //=== following lines have been copied from CsGeom.cc ===
    int id;        s >> id;   
    string TBName; s >> TBName; 
    if(  ( TBName == TBName_ ) 
      
      //=== Special case for the straws
      || ( TBName_.substr(0,2) == "ST" && TBName.substr(0,7) == TBName_.substr(0,7) ) 
      
      )
      
    _AddDetFileInfo( new DetFileInfo( detFile_, lineNumber ,line ) ); 
  
  }  

  cout << "DetInfo::ReadDetFile \"" << TBName_ << "\" - INFO: " << detFList_.size() << " matches.\n";
  return (bool) detFList_.size();

}

//_________________________________________________________________________________________
bool DetInfo::_GetTBNameFromDetFile( int ID )
{
   
  ifstream in( detFile_.c_str(), ios::in );
  if( !in ) {
    cout << "DetInfo::_GetDetFileInfo - ERROR: cannot read file \"" 
         << detFile_  << "\"." 
         << endl;
    return false;
  }
    
  char line[512];
  unsigned int lineNumber = 0;
  
	while( !in.eof() ) { 
    in.getline( line, sizeof( line ), '\n');
    if( in.eof() || !in.good() ) continue;
    
    if( line[0] == '\0' || line[0] != ' ' ) continue;    
    if( line[1] == ' ' || line[1] == '-' || line[1] == '=' ) continue;
    lineNumber++;      
    
    istrstream s(line);
    string opt;
    s >> opt;

		if ( opt != "det" ) continue;
    
    //=== following lines have been copied from CsGeom.cc ===
    int id;        s >> id;   
    string TBName; s >> TBName; 
    
    if( id == ID ) { 
      TBName_ = TBName;
      return true;
    }
  }  
  
  return false;
  
}
