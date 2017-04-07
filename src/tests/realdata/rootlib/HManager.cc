// $Id: HManager.cc,v 1.35 2002/12/09 08:51:06 hpereira Exp $
#include "HManager.h"
#include "DetInfo.h"
#include "Utils.h"
#include "DetFileInfo.h"
#include <functional>
#include <list.h>
#include <strstream.h>

#include <TKey.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TStyle.h>

ClassImp( HManager )
HManager* HManager::ptr_ = 0;

//_______________________________________________________________________________
HManager::HManager( void ):
  TObject(),
  tcv_( 0 ),
  fileSelection_( "" ),
  OneGaus_( true ),
  Batch_( false ),
  UseCorrectedVals_( false ),
  color_( 1 )
{ 
  detInf_.clear();
  ptr_ = this;
}

//_______________________________________________________________________________
HManager::HManager( const char* fileSelection, bool batch = false ):
  TObject(),
  tcv_( 0 ),
  fileSelection_( fileSelection ),
  OneGaus_( true ),
  UseCorrectedVals_( 0 )
{

  //=== some root parameters
  gROOT->SetStyle( "Plain" );

  gStyle->SetPalette( 1 );
  gStyle->SetOptStat( 1 );
  gStyle->SetOptFit( 111 );
  gStyle->SetOptDate( 1 );
  
  tcv_ = new TCanvas( ( (Batch_)?"":"HManager_Cv"),"HManager",  700, 700 );
  if( Batch_ ) tcv_->SetBatch( true );
  detInf_.clear();   
  
  //=== Get Selected files
  vector< string > selFiles = _GetFiles( fileSelection ); 
  for( unsigned int i = 0; i < selFiles.size(); i++ ) {
    TFile *tf = new TFile( selFiles[i].c_str() );
    if( tf == 0 || !tf->IsOpen() ) {
      cout << "HManager::HManager - ERROR: problem with TFile \"" << selFiles[i] << "\".\n";
      continue;
    } 
    
    _LoadTBNames( tf );
  }

  _SortDetInfos();
  cout << "HManager::HManager - INFO: Following TBNames loaded: "<< endl;
  for( unsigned int i = 0; i < detInf_.size(); i++ )  
  cout << "\"" << detInf_[i]->TBName_ 
    << "\" => " 
    << detInf_[i]->T_eff_->GetEntries() 
    << " entries.\n";
  
  ptr_ = this;

}

//_______________________________________________________________________________
void HManager::DumpTBNames( const char* selection = "*" )
{
  GetSelection( selection );
  for( unsigned int iInf = 0; iInf < selDetI_.size(); iInf++ )
  cout << "HManager::DumpTBNames - " << selDetI_[iInf]->TBName_ << endl;
}


//_______________________________________________________________________________
void HManager::Draw( const char* varSel, const char* selection = "*", TCut cut="1", const char* opt="" ) 
{
  
  //=== Get Selection
  GetSelection( selection );
  vector< DetInfo*> DetOK;
  DetOK.clear();
  for( unsigned int iInf = 0; iInf < selDetI_.size(); iInf++ )
  if( selDetI_[iInf]->T_eff_ ) 
  DetOK.push_back( selDetI_[iInf] );
  
  if( DetOK.empty() ) {
    cout << "HManager::Draw - WARNING: works only with Efficiency tree.\n";
    return;
  }
  
  bool same = ( string( opt ).find("same") < strlen( opt ) );
  if( same ) {
    cout << "HManager::Fit INFO: same selected" << endl;
    color_++;
  } else color_ = 1;
  
  if( !same && tcv_ ){
    tcv_->Clear();
    _DivideCanvas( DetOK.size() );
  }
  
  cout << "HManager::Draw - Color is " << color_ << endl;
  for( unsigned int iInf = 0; iInf < DetOK.size(); iInf++ ) {
    if( tcv_ ) tcv_->cd( iInf+1 );
    string name = DetOK[iInf]->TBName_+"_H";
    string title = DetOK[iInf]->TBName_+"_"+varSel+":"+string( cut );
    string varSelH = string(varSel)+">>"+name;
    TH1* h;
    if( (h = (TH1*) gDirectory->Get( name.c_str() ) ) ) delete h;
    int n = (int) DetOK[iInf]->T_eff_->Draw(varSelH.c_str(), cut, "goff" );
    h = (TH1*) gDirectory->Get(name.c_str());
      
    cout << "HManager::Draw - " << DetOK[iInf]->TBName_ << ": " << ((h)? n:-1) << " entries.\n";
    if( !(h && h->GetEntries() ) ) continue;
 
    h->SetTitle( title.c_str() );
    h->SetLineColor( color_ );
    h->SetMarkerColor( color_ );
    h->Draw(opt);
    if( tcv_ ) tcv_->Update();
  }
  
  return;
}  
  
//_______________________________________________________________________________
void HManager::Fit( const char* varSel, const char* selection = "*", TCut cut="1", const char* func="gaus", const char* opt="" ) 
{
  
  //=== Get Selection
  GetSelection( selection );
  vector< DetInfo*> DetOK;
  DetOK.clear();
  for( unsigned int iInf = 0; iInf < selDetI_.size(); iInf++ )
  if( selDetI_[iInf]->T_eff_ ) 
  DetOK.push_back( selDetI_[iInf] );
  
  if( DetOK.empty() ) {
    cout << "HManager::Draw - WARNING: works only with Efficiency tree.\n";
    return;
  }
  
  //=== check if 1D varSel
  if( string( varSel ).find(":") < strlen( varSel ) ) {
    cout << "HManager::Draw - WARNING: works only with 1D histos.\n";
    return;
  }
  
  //=== check if same is selected
  bool same = ( string( opt ).find("same") < strlen( opt ) );
  if( same ) cout << "HManager::Fit INFO: same selected" << endl;
  if( !same && tcv_ ){
    tcv_->Clear();
    _DivideCanvas( DetOK.size() );
  }
  
 
  for( unsigned int iInf = 0; iInf < DetOK.size(); iInf++ ) {
    if( tcv_ ) tcv_->cd( iInf+1 );
    string name = DetOK[iInf]->TBName_+"_H";
    string title = DetOK[iInf]->TBName_+"_"+varSel+":"+string( cut );
    string varSelH = string(varSel)+">>"+name;
    TH1* h;
    if( (h = (TH1*) gDirectory->Get( name.c_str() ) ) ) delete h;
    int n = (int) DetOK[iInf]->T_eff_->Draw(varSelH.c_str(), cut, opt );
    h = (TH1*) gDirectory->Get(name.c_str());
      
    cout << "HManager::Fit - " << DetOK[iInf]->TBName_ << ": " << ((h)? n:-1) << " entries.\n";
    if( !(h && h->GetEntries()  ) ) continue;

    h->SetTitle( title.c_str() );
    h->Fit( func );
    if( tcv_ ) tcv_->Update();
  }
  
  return;
}  


//_______________________________________________________________________________
void HManager::SetDetFile( const char* file ) 
{
  for( unsigned int i = 0; i< detInf_.size(); i++ ) 
  detInf_[i]->SetDetFile( file );
  _SortDetInfos( );
}

//_______________________________________________________________________________
DetInfo* HManager::GetDetInfo( const char* TBName = "" ) 
{
  string TBN(TBName);
  for( unsigned int i = 0; i < detInf_.size(); i++ )
  if( detInf_[i]->TBName_.substr(0,TBN.size()) == TBN ) return detInf_[i];
  return 0;
}
    
//_______________________________________________________________________________
vector<DetInfo*> HManager::GetSelection( const char* selectionList )
{
  
  cout << "HManager::GetSelection -INFO: OneGaus_ is " << OneGaus_ << endl;
  cout << "HManager::GetSelection -INFO: UseCorrectedVals_ is " << UseCorrectedVals_ << endl;

  if( strlen( selectionList ) == 0 ) {
    for( unsigned int iInf = 0; iInf < detInf_.size(); iInf++ ) {
      if( detInf_[iInf] )   detInf_[iInf]->OneGaus_ = OneGaus_;
      if( detInf_[iInf] )   detInf_[iInf]->UseCorrectedVals_ = UseCorrectedVals_;
    }  
    return detInf_;
  }
  
  selDetI_.clear();
  vector< string > selectionVect;
  istrstream in( selectionList );
  do{
    string sel; in >> sel;
    selectionVect.push_back( sel );
  } while( in );
  
  for( unsigned int iSel = 0; iSel < selectionVect.size(); iSel++ ) {
    string selection = selectionVect[iSel];
    if( !selection.size() ) continue;
    for( unsigned int iInf = 0; iInf < detInf_.size(); iInf++ ) {
      if( selection.size() > detInf_[iInf]->TBName_.size() ) continue;
      bool accept = true;
      for( unsigned int i=0; i < selection.size(); i++ ) 
      if( selection[i]!='*' && selection[i] != detInf_[iInf]->TBName_[i] ) {  
        accept = false;
        break;
      }
      if( accept ) {
        if( detInf_[iInf] )   detInf_[iInf]->OneGaus_ = OneGaus_;
        if( detInf_[iInf] )   detInf_[iInf]->UseCorrectedVals_ = UseCorrectedVals_;
        selDetI_.push_back( detInf_[iInf] );
        cout << "HManager::GetSelection - Add \"" << detInf_[iInf]->TBName_ << "\".\n";
      }
    }
  }
  cout << "done.\n";
  return selDetI_;
}

//=== protected functions    
//_______________________________________________________________________________
void HManager::_LoadTBNames( TFile *tf )
{    
  //=== Get Key and check
  TKey* tk =  tf->GetKey("CsEfficiency");
  if( !tk ) {
    cout << "HManager::_LoadTBNames - ERROR: cannot find key \"CsEfficiency\" in \"" 
      << tf->GetName() 
      << "\".\n";
    return;
  }
  
  //=== Get TDirectory and check
  TDirectory *dir = (TDirectory*) tk->ReadObj();
  if( !dir ) {
    cout << "HManager::_LoadTBNames - ERROR: cannot find directory \"CsEfficiency\" in \"" 
      << tf->GetName() 
      << "\".\n";
    return;
  }
  
  TList* list = dir->GetListOfKeys();

  //=== Check all objects in dir ===
  string k = "_eTree_";
  for( TObject* obj = list->First(); obj != 0; obj = list->After( obj ) )  
  if( string( obj->GetName(),8,k.size() ) == k ) {
    string TBName = string( obj->GetName(), 0, 8 ); 
    DetInfo* detI = GetDetInfo( TBName.c_str() );
    if( !detI ) {
      detI = new DetInfo( TBName );
      detInf_.push_back( detI );
      string chainName = string( "CsEfficiency/" ) + obj->GetName();
      detI->AddToChain( chainName.c_str(), tf->GetName() );     
    } else detI->AddToChain( "", tf->GetName() );     

  }  

  return;
}

//____________________________________________________________________________
vector< string > HManager::_GetFiles( const char* fileSelection )
{
  //=== Get Selected files
  vector< string > selFiles; 
  char command[128];
  
  string tmpFile("/tmp/_HManager__GetFiles_");
  tmpFile += string( Utils::GetTimeStamp() );
  sprintf( command, "ls -1 %s > %s", fileSelection, tmpFile.c_str() );
  system( command );
  
  ifstream in( tmpFile.c_str(), ios::in );
  if( !in ) {
    cout << "HManager::_GetFiles - FATAL: cannot execute command \"" << command << "\".\n";
    return selFiles;
  } 
  
  char line[128];
  while( !in.eof() ) {
    in.getline( line, sizeof( line ) );
    if( in.eof() || line[0] == '\0' ) continue;
    
    istrstream s(line);
    string name; s >> name;
    if( access( name.c_str(), R_OK ) ){
      cout << "HManager::_GetFiles- INFO: cannot access file \"" << name << "\".\n";
      continue;
    }
  
    selFiles.push_back( name );
  }

  return selFiles;
}  
  
//_______________________________________________________________________________
void HManager::_DivideCanvas( const unsigned int nPads )
{
  cout << "HManager::_DivideCanvas.\n";
  if( tcv_ == 0 || nPads<=1 ) return;
  unsigned int sqnP = (unsigned int) (sqrt( nPads ));
  unsigned int nC = sqnP;
  unsigned int nL = sqnP;
  
  while( nC*nL < nPads ) 
  if( nC < nL ) nC++;
  else nL++;
   
  tcv_->Divide( nC, nL, 0.005, 0.005, 0 );
  cout << "HManager::_DivideCanvas - Done.\n";
  return;
}

//_______________________________________________
struct HManager::sortDetInfos_ : public binary_function< DetInfo*, DetInfo*, bool > 
{ 

  bool operator() ( DetInfo* di0, DetInfo* di1 ) 
  { if( !( di0->hasDetFile_ && di0->hasDetFile_ ) ) 
    return ( di0->TBName_ < di1->TBName_ ); 
    else return( di0->detFI_->lineNumber_ < di1->detFI_->lineNumber_ );
  } 

};

//_______________________________________________
void HManager::_SortDetInfos( void )
{ 
  //===
  cout << "HManager::_SortDetInfos.\n";
  list<DetInfo*> tmp; tmp.clear();
  for( unsigned int i=0; i<detInf_.size(); i++ ) tmp.push_back( detInf_[i] );
  tmp.sort( HManager::sortDetInfos_() );
  detInf_.clear();
  
  //===
  list<DetInfo*>::iterator Id;
  for( Id = tmp.begin(); Id != tmp.end(); Id++ ) detInf_.push_back( *Id );
  
  return;
}
