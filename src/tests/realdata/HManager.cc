// $Id: HManager.cc,v 1.5 2008/06/27 16:35:20 conrad Exp $
 
/*!
   \file    HManager.cc
   \brief   Calibration Chain Managment Interface Class.
   \author  Hugo Pereira
   \version $Revision: 1.5 $
   \date    $Date: 2008/06/27 16:35:20 $
*/

#include "HManager.h"
#include "DetFileManager.h"
#include "DetectorInfo.h"
#include "Utils.h"
#include "Defs.h"

#include <sstream>
#include <fstream>
#include <list>

#include <TStyle.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TObjString.h>
#include <TH1.h>
#include <TH2.h>
#include <TCut.h>

using namespace std;

//________________________________________
ClassImp(HManager)

// initialise pointer to current HManager object
HManager* HManager::ptr_ = 0;

//________________________________________
HManager::HManager( const char* fileselection, const char* detectorfile ):
  TObject(),
  OneGaus_( true ),
  df_(0),
  detectorFile_( "" ),
  chain_( 0 ),
  nEntries_( 0 ),
  color_( 1 ),
  cv_( 0 )
{
  allDets_.clear();
  allTBNames_.clear();
  
  if( strlen( detectorfile ) )     LoadDetectorFile( detectorfile );
  if( strlen( fileselection ) )    LoadTrees( fileselection ); 
  if( !LoadDetectorsFromString() ) LoadDetectorsFromChain();
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptDate();
  gStyle->SetOptFile();
  cv_ = new TCanvas( "HManager_Cv","HManager",  700, 700 );
  
  ptr_ = this;
}


//________________________________________
DetFileManager* HManager::LoadDetectorFile( const char* name )
{
  if( df_ ) SafeDelete(df_);
  df_ = new DetFileManager( name );
  detectorFile_ = string(name);
  return df_;
}
 
//________________________________________
int HManager::LoadDetectorsFromChain( void )
{
  cout << "HManager::LoadDetectorsFromChain.\n";
  allDets_.clear();

  if( !df_ )    { cout << "HManager::LoadDetectorsFromChain - no detector file loaded.\n"; return 0; }
  if( !chain_ ) { cout << "HManager::LoadDetectorsFromChain - no chain loaded.\n"; return 0; }
  
  // create detectorIds histogram
  string name = "detIds_";
  string var = "T_detVect>>" + name;
  TH1 *h = 0;
  if( (h = (TH1*) gDirectory->Get(name.c_str())) ) SafeDelete(h);
  chain_->Draw( var.c_str(), "1", "goff" );
  h = (TH1*) gDirectory->Get(name.c_str());
  if( !h ) return 0;
  
  vector< DetectorInfo* > dets = df_->GetDetInfo();
  for( unsigned int i=0; i< dets.size(); i++ ) {
    
    int id = dets[i]->id_;
    if( id > h->GetXaxis()->GetXmax() || id < h->GetXaxis()->GetXmin() ) continue;
    int bin = h->FindBin( id );
    int entries = (int)  h->GetBinContent( bin );
    if( entries && allDets_.find( dets[i] ) == allDets_.end() ) {
      allDets_.insert( dets[i] );
      cout << "HManager::LoadDetectors - \"" << dets[i]->TBName_ << "\" " << entries << " entries.\n";
    }
  
  }
  
  return allDets_.size();

}

//________________________________________
int HManager::LoadDetectorsFromString( void )
{
  cout << "HManager::LoadDetectorsFromString.\n";
  allDets_.clear();
  if( !df_ )    { cout << "HManager::LoadDetectorsFromString - no detector file loaded.\n"; return 0; }
  if( !chain_ ) { cout << "HManager::LoadDetectorsFromString - no chain loaded.\n"; return 0; }
  
  vector< DetectorInfo* > dets = df_->GetDetInfo();
  for( unsigned int i=0; i< dets.size(); i++ ) {

    string TBName = dets[i]->TBName_;
    if( allTBNames_.find( TBName ) != allTBNames_.end() && allDets_.find( dets[i] ) == allDets_.end() ) { 
      allDets_.insert( dets[i] );
      cout << "HManager::LoadDetectors - \"" << dets[i]->TBName_ << "\" loaded.\n";
    }
  
  }
  
  return allDets_.size();
}
      
//________________________________________
TChain* HManager::LoadTrees( const char* fileselection )
{
  if( chain_ ) {
    SafeDelete( chain_ );
    chain_ = 0;
  }
  
  AddToChain( fileselection );
  return chain_;
}

//________________________________________
TChain* HManager::AddToChain( const char* fileselection )
{  
    
  vector< string > files = Utils::GetFiles( string(fileselection) );
  
  // try loading tree
  for( unsigned int i=0; i<files.size(); i++ ) {
    string file = files[i];
  
    // First call
    if( !chain_ ) {
      
      // open TFile
      TFile *tf = new TFile( file.c_str() );
      if( tf == 0 || !tf->IsOpen() ) {
        cout << "HManager::AddToChain - ERROR: troubles with TFile \"" << file << "\".\n";
        continue;
      } 
      
      // Get Tree
      const char* treeName = ROOT_PATH "/" TREE_NAME;
      TTree* tree;
	  tf->GetObject(treeName, tree);
      if( !tree ) {
        cout << "HManager::AddToChain - ERROR: Tree \"" << treeName << "\" not found in \"" << file << "\".\n";
        continue;
      }
      
      chain_ = new TChain( treeName, treeName );
      chain_->Add( file.c_str() );

    } else chain_->Add( file.c_str(), -1 ); 
    if( chain_ ) nEntries_ += (unsigned int) chain_->GetEntries();
    cout<<"HManager::AddToChain - file \""<<file.c_str()<<"\" added - "<<nEntries_<<" entrie(s)."<<endl;
  }
  
  // try reading list of recorded TBNames
  for( unsigned int i=0; i<files.size(); i++ ) {
    
    // reopen TFile
    string file = files[i];
    TFile *tf = new TFile( file.c_str() );
    if( tf == 0 || !tf->IsOpen() ) continue;
    
    // try loading TBName_list_
    const char* listName = ROOT_PATH "/" LIST_NAME;
    TObjString* T_list;
	tf->GetObject(listName, T_list);
    if( !T_list ) { 
      cout << "HManager::AddToChain - WARNING: TObjString \"" << listName << "\" not found in \"" << file << "\".\n";
      return chain_;
    }
    
    // dump all TBNames into allTBNames_
    istringstream in( T_list->GetName() );
    while( in ) { 
      string name; 
      in >> name; 
      if( name.size() && allTBNames_.find( name )== allTBNames_.end() ) 
      allTBNames_.insert( name );
    }
  }
    
      
  return chain_;
}

//________________________________________
void HManager::Draw(  const char* var, const char* selection, const char* cut, const char* opt, int nEnt )
{
  if( !_Select( selection ) ) return;
  
  unsigned int nPads = dets_.size();
  if( string( opt ).find("same")==string::npos ) cv_->Clear();
  if( nPads ) Utils::DivideCanvas( cv_, nPads );
  for( unsigned int i=0; i<nPads; i++ ) {
    
    // select main detector
    DetectorInfo *det = dets_[i]->GetMain();
    TH1 *h = _Get( var, det, cut, nEnt );
    if( nPads > 1 ) cv_->cd(i+1);
        
    cout << "HManager::Draw - var \"" 
      << var 
      << "\" det \"" 
      << det->TBName_.c_str() 
      << "\": " 
      << int((h)?h->GetEntries():-1)
      << " entries.\n";

    if( h && h->GetEntries() ) h->Draw(opt); 
    cv_->Update();   
  }
  
}
      
//________________________________________
TH1* HManager::_Get(  const char* var, DetectorInfo *d, const char* cut, int nEnt )
{
  if( !d ) return 0;
  if( !chain_ ) { cout << "HManager::_Get - no tracks loaded.\n"; return 0; }
  if( !nEnt ) nEnt = (int) chain_->GetEntries();

  // Set Cut
  char *buf = new char[64];
  sprintf( buf, "T_detVect==%i",d->id_ );
  TCut detCut = TCut(cut) + TCut(buf);
  SafeDelete(buf);    
  
  // Set histogram title
  string var_H = string(var)+">>"+d->TBName_;
  TH1 *h;
  if( (h = (TH1*) gDirectory->Get(d->TBName_.c_str())) ) SafeDelete(h);
  chain_->Draw( var_H.c_str(), detCut, "goff", nEnt );

  // Retrieve histogram
  h = (TH1*) gDirectory->Get(d->TBName_.c_str());
  return h;
          
}
     
//________________________________________
TH1* HManager::_GetCloned(  const TH1* hIn, const char* var, DetectorInfo *d, const char* cut, int nEnt )
{
  if( !d ) return 0;
  if( !hIn )    { cout << "HManager::_Get - ref histogram is 0.\n"; return 0; }
  if( !chain_ ) { cout << "HManager::_Get - no tracks loaded.\n"; return 0; }
  if( !nEnt ) nEnt = (int) chain_->GetEntries();

  // Set Cut
  char *buf = new char[64];
  sprintf( buf, "T_detVect==%i",d->id_ );
  TCut detCut = TCut(cut) + TCut(buf);
  SafeDelete(buf);    
  
  // create output histogram name/title
  string nameOut  = string( hIn->GetName() )  + "_clone";
  string titleOut = string( hIn->GetTitle() ) + "_clone";
  
  // look for/delete output histogram
  TH1 *hOut;
  if( (hOut = (TH1*) gDirectory->Get(nameOut.c_str())) ) SafeDelete(hOut);
  
  // book output histogram
  hOut = (TH1*) hIn->Clone();
  if( !hOut ) {
    cout << "HManager::_GetCloned() - clone failed.\n";
    return 0;
  }
  
  hOut->SetName( nameOut.c_str() );
  hOut->SetTitle( titleOut.c_str() );

  // project tree in histogram
  chain_->Project( nameOut.c_str(), var, detCut, "", nEnt );
  return hOut;
          
}
                 
//________________________________________
int HManager::_Select( const char* selection, const char* cut, int nEnt )
{
  dets_.clear();
  if( !df_ )    { cout << "HManager::_Select - no detectorfile loaded.\n"; return 0; }
  if( !chain_ ) { cout << "HManager::_Select - no tracks loaded.\n";       return 0; }        
  vector< DetectorInfo* > dets = df_->GetDetSelection( selection );
  cout<<dets.size()<<endl;
  for( unsigned int i=0; i< dets.size(); i++ ) {
      //if( allDets_.find( dets[i] ) != allDets_.end() ) dets_.push_back( dets[i] );
      for( unsigned int j=0; j<dets[i]->sub_.size(); j++) {
	  if( allDets_.find( dets[i]->sub_[j] ) != allDets_.end() ) dets_.push_back( dets[i]->sub_[j] );
      }
  }
  cout << "HManager::_Select - " << dets_.size() << " detectors selected.\n";
  return dets_.size();
}
