// $Id: Macro.cc,v 1.24 2002/11/06 14:04:11 hpereira Exp $
#include "Macro.h"
#include "EffManager.h"
#include "ResManager.h"
#include "RTManager.h"
#include "DetInfo.h"
#include "Utils.h"
#include "DetInfo.h"
#include "CsRTGridPoint.h"

#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>

ClassImp( Macro )

// define output tree and corresponding constants  
const unsigned int nCh_DC = 192;
const unsigned int stp_DC = 8;

//____________________________________________________________________________
// Fit time versus V for all selected rootfiles.
// Write the corresponding slope to a file. 
void Macro::PropTime( const char* fileSelection, 
  const char* cut = "abs(T_vLx)<300", 
  const char* outFile = "PropTime.out" )
{

  ofstream out( outFile, ios::out );
  if( !out ) {
    cout << "Macro::PropTime - FATAL: cannot write to file \"" << outFile << "\".\n";
    return;
  }
     
  //=== Get Selected files
  vector< string > selFiles = _GetFiles( fileSelection ); 
  vector< TCanvas* > Cvs;
  
  for( unsigned int i=0; i<selFiles.size(); i++ ) {
    cout << "Macro::PropTime - INFO: got file \"" << selFiles[i] << "\".\n";
    
    EffManager *e = new EffManager( selFiles[i].c_str() );
    DetInfo *d = e->GetDetInfo();
    if( !d ) continue;
    
    TF1 *f = d->FitTvsV(cut);
    if( !f ) continue;
    
    //=== Parse selected file to get detector TBName
    string TBName = selFiles[i].substr( selFiles[i].rfind(".")+1, selFiles[i].size() );
    
    out << TBName << " ";
    out.form("%9.6f",f->GetParameter(1) );
    out << " ns/mm <==\"" << selFiles[i] << "\".\n";
   
    //=== change canvas name
    char tmp[64];
    sprintf( tmp, "%s_%i", e->GetTCanvas()->GetName(),i );  e->GetTCanvas()->SetName( tmp );
    sprintf( tmp, "%s_%i", e->GetTCanvas()->GetTitle(),i ); e->GetTCanvas()->SetTitle( tmp );
    Cvs.push_back( e->GetTCanvas() );
    delete e;  
  }
  
  out.close();
     
  printf("Press <Return> \n");  
  getchar();
  
  //=== Delete Canvases
  for( unsigned int i = 0; i< Cvs.size(); i++ ) { Cvs[i]->Close(); }
  
  return;
}     

//____________________________________________________________________________
// open all selected files with ResManager
// FitDU for found DetInfo, if any, using <cut>
// Save the result in automaticaly generated ps file, 
// to wich <tag> is added, for comparisons
void Macro::DUFit( const char* fileSelection, 
  const char* cut="", 
  const char* tag="",
  bool oneGaus = false,
  bool UseCorrectedVals = false )
{
  //=== Get Selected files
  vector< string > selFiles = _GetFiles( fileSelection ); 
  vector< TCanvas* > Cvs;

  //=== initialize Macro write status
  Utils::SetWritableStatus( Utils::NO );
    
  // Run EffManager on all files
  for( unsigned int i = 0; i < selFiles.size(); i++ )
  { 
    
    //=== Generate ps name
    string psName = ( UseCorrectedVals ) ?
      
      string( "/afs/cern.ch/user/h/hpereira/eps/coral.rd/duCor.") +
      selFiles[i].substr( selFiles[i].find(".root.")+strlen(".root."), selFiles[i].size()-strlen(".root.")) +
      tag + ".eps" : 
      
      string( "/afs/cern.ch/user/h/hpereira/eps/coral.rd/duMin.") +
      selFiles[i].substr( selFiles[i].find(".root.")+strlen(".root."), selFiles[i].size()-strlen(".root.")) +
      tag + ".eps";
    
      
    cout << "Macro::FitDU - INFO: \"" 
      << selFiles[i] 
      << "\" => \"" 
      << psName 
      << "\".\n";
      
    //=== Set options
    ResManager *r = new ResManager( selFiles[i].c_str() );
    r->OneGaus_ = oneGaus;
    r->UseCorrectedVals_ = UseCorrectedVals;
    r->FitDU( "*", cut );

    //=== Generate postscript file 
    Utils::WStatus wStat = Utils::Writable( psName ); 
    if( wStat == Utils::EXIT ) break; 
    else if( wStat == Utils::YES ) r->GetTCanvas()->SaveAs( psName.c_str() );

    //=== change canvas name
    char tmp[64];
    sprintf( tmp, "%s_%i", r->GetTCanvas()->GetName(),i );  r->GetTCanvas()->SetName( tmp );
    sprintf( tmp, "%s_%i", r->GetTCanvas()->GetTitle(),i ); r->GetTCanvas()->SetTitle( tmp );
    Cvs.push_back( r->GetTCanvas() );
    delete r;
    
  }
     
  printf("Press <Return> \n");  
  getchar();
  
  //=== Delete Canvases
  for( unsigned int i = 0; i< Cvs.size(); i++ ) { Cvs[i]->Close(); }
  
  return;
}

//____________________________________________________________________________
// open all selected files with ResManager
// FitDU for found DetInfo, if any, using <cut>
// Save the result in automaticaly generated ps file, 
// to wich <tag> is added, for comparisons
void Macro::RTGrid( const char* fileSelection, 
  const char* cut="", 
  const char* tag="",
  unsigned int n=0 )
{
  //=== Get Selected files
  vector< string > selFiles = _GetFiles( fileSelection ); 
  vector< TCanvas* > Cvs;

  //=== initialize Macro write status
  Utils::SetWritableStatus( Utils::NO );
    
  // Run EffManager on all files
  for( unsigned int i = 0; i < selFiles.size(); i++ )
  { 
    
    //=== Generate ps name
    string psName =
    string( "/afs/cern.ch/user/h/hpereira/eps/coral.rd/rtGrid.") +
    selFiles[i].substr( selFiles[i].find(".root.")+strlen(".root."), selFiles[i].size()-strlen(".root.")) +
    tag + ".eps";
    
    cout << "Macro::RTGrid - INFO: \"" 
      << selFiles[i] 
      << "\" => \"" 
      << psName 
      << "\".\n";
      
    RTManager *r = new RTManager( selFiles[i].c_str() );
    r->FitRTStraights( "*", cut );
    r->FitRTGrids( "*", cut, n );

    //=== Generate postscript file 
    Utils::WStatus wStat = Utils::Writable( psName ); 
    if( wStat == Utils::EXIT ) break; 
    else if( wStat == Utils::YES ) r->GetTCanvas()->SaveAs( psName.c_str() );

    //=== change canvas name
    char tmp[64];
    sprintf( tmp, "%s_%i", r->GetTCanvas()->GetName(),i );  r->GetTCanvas()->SetName( tmp );
    sprintf( tmp, "%s_%i", r->GetTCanvas()->GetTitle(),i ); r->GetTCanvas()->SetTitle( tmp );
    Cvs.push_back( r->GetTCanvas() );
    delete r;
    
  }
     
  printf("Press <Return> \n");  
  getchar();
  
  //=== Delete Canvases
  for( unsigned int i = 0; i< Cvs.size(); i++ ) { Cvs[i]->Close(); }
  
  return;
}

//______________________________________________________________________________
// open all selected files with EffManager
// DrawEffvsXYs for found detectors, if any, using <cut>
// Save the result in automaticaly generated ps file, 
// to wich <tag> is added, for comparisons
void Macro::EffvsXY( const char* fileSelection, 
  const char* cut="", 
  const char* tag="" )
{
  //=== Get Selected files
  vector< string > selFiles = _GetFiles( fileSelection ); 
  vector< TCanvas* > Cvs;

  //=== initialize Macro write status
  Utils::SetWritableStatus( Utils::NO );
    
  // Run EffManager on all files
  for( unsigned int i = 0; i < selFiles.size(); i++ )
  { 
    
    //=== Generate ps name
    string psName =
    string( "/afs/cern.ch/user/h/hpereira/eps/coral.rd/EffvsXY.") +
    selFiles[i].substr( selFiles[i].find(".root.")+strlen(".root."), selFiles[i].size()-strlen(".root.")) +
    tag + ".eps";
    
    cout << "Macro::DrawEffvsXY - INFO: \"" 
      << selFiles[i] 
      << "\" => \"" 
      << psName 
      << "\".\n";
      
    EffManager *e = new EffManager( selFiles[i].c_str() );
    e->DrawEff( "T_yLx:T_xLx", "*", cut );

    //=== Generate postscript file
    Utils::WStatus wStat = Utils::Writable( psName );
    if( wStat == Utils::EXIT ) break;
    else if( wStat == Utils::YES ) e->GetTCanvas()->SaveAs( psName.c_str() );
   
    //=== change canvas name
    char tmp[64];
    sprintf( tmp, "%s_%i", e->GetTCanvas()->GetName(),i );  e->GetTCanvas()->SetName( tmp );
    sprintf( tmp, "%s_%i", e->GetTCanvas()->GetTitle(),i ); e->GetTCanvas()->SetTitle( tmp );
    Cvs.push_back( e->GetTCanvas() );
    delete e;
    
  }
     
  printf("Press <Return> \n");  
  getchar();
  
  //=== Delete Canvases
  for( unsigned int i = 0; i< Cvs.size(); i++ ) { Cvs[i]->Close(); }
  
  return;
}

//____________________________________________________________________________
// open all selected files with EffManager
// DrawEffvsUVs for found detectors, if any, using <cut>
// Save the result in automaticaly generated ps file, 
// to wich <tag> is added, for comparisons
void Macro::EffvsUV( const char* fileSelection, 
  const char* cut="", 
  const char* tag="" )
{
  //=== Get Selected files
  vector< string > selFiles = _GetFiles( fileSelection ); 
  vector< TCanvas* > Cvs;

  //=== initialize Macro write status
  Utils::SetWritableStatus( Utils::NO );
    
  //=== Run ResManager on all files
  for( unsigned int i = 0; i < selFiles.size(); i++ )
  { 
    
    //=== Generate ps name
    string psName =
    string( "/afs/cern.ch/user/h/hpereira/eps/coral.rd/EffvsUV.") +
    selFiles[i].substr( selFiles[i].find(".root.")+strlen(".root."), selFiles[i].size()-strlen(".root.")) +
    tag + ".eps";
    
    cout << "Macro::DrawEffvsUV - INFO: \"" 
      << selFiles[i] 
      << "\" => \"" 
      << psName 
      << "\".\n";
      
    EffManager *e = new EffManager( selFiles[i].c_str() );    
    e->DrawEff( "T_vLx:T_uLx", "*", cut );

    //=== Generate postscript file
    Utils::WStatus wStat = Utils::Writable( psName );
    if( wStat == Utils::EXIT ) break;
    else if( wStat == Utils::YES ) e->GetTCanvas()->SaveAs( psName.c_str() );
   
    //=== change canvas name
    char tmp[64];
    sprintf( tmp, "%s_%i", e->GetTCanvas()->GetName(),i );  e->GetTCanvas()->SetName( tmp );
    sprintf( tmp, "%s_%i", e->GetTCanvas()->GetTitle(),i ); e->GetTCanvas()->SetTitle( tmp );
    Cvs.push_back( e->GetTCanvas() );
    delete e;
    
  }
     
  printf("Press <Return> \n");  
  getchar();
  
  //=== Delete Canvases
  for( unsigned int i = 0; i< Cvs.size(); i++ ) { Cvs[i]->Close(); }
  
  return;
}

//__________________________________________________________    
void Macro::ProfileToRates( TH2* h, 
    unsigned int nEvt = 1, 
    double tMin=0, double tMax=0 )
{
  if( tMin==0 && tMax==0 ) {
    tMin = h->GetYaxis()->GetXmin();
    tMax = h->GetYaxis()->GetXmax();
  }
  
  cout << "Macro::ProfileToRates - INFO: time limits: [" 
    << tMin << "," 
    << tMax << "]\n";
  
  int nMin = h->GetYaxis()->FindBin( tMin );
  int nMax = h->GetYaxis()->FindBin( tMax );
  
  
  TH1D* h1;
  string name( "rates" );
  if( (h1 = (TH1D*) gDirectory->Get( name.c_str() ) ) ) delete h1;
  h->ProjectionX( name.c_str(), nMin, nMax );
  h1 = (TH1D*) gDirectory->Get( name.c_str() );
  
  if( !( h1 && h1->GetEntries() ) ) { 
    cout << "Macro::ProfileToRates - Histogram is empty.\n";
    return;
  }
  
  h1->SetMarkerStyle(20);
  h1->SetLineWidth( 2 );
  h1->SetLineColor( 4 );
  h1->SetMarkerColor(4);
  h1->SetMarkerSize(0.7);
  h1->Scale(1e6/(nEvt*(tMax-tMin)));
  h1->Draw("p");
}

//=== PRIVATE METHODS ===
//____________________________________________________________________________
vector< string > Macro::_GetFiles( const char* fileSelection )
{
  //=== Get Selected files
  vector< string > selFiles; 
  char command[128];
  
  string tmpFile("/tmp/Macro_FitDU");
  sprintf( command, "ls -1 %s > %s", fileSelection, tmpFile.c_str() );
  system( command );
  
  ifstream in( tmpFile.c_str(), ios::in );
  if( !in ) {
    cout << "Macro::FitDU - FATAL: cannot execute command \"" << command << "\".\n";
    return selFiles;
  } 
  
  char line[128];
  while( !in.eof() ) {
    in.getline( line, sizeof( line ) );
    if( in.eof() || line[0] == '\0' ) continue;
    
    istrstream s(line);
    string name; s >> name;
    if( access( name.c_str(), R_OK ) ){
      cout << "Macro::FitDU - INFO: cannot access file \"" << name << "\".\n";
      continue;
    }
  
    selFiles.push_back( name );
  }

  return selFiles;
}  
