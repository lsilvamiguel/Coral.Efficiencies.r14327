// $Id: Tracks.cc,v 1.19 2010/02/03 18:24:05 suhl Exp $
/*!
   \file    Tracks.cc
   \brief   chains and loops over coral alignment tree outputs
   \author  Hugo Pereira
   \version $Revision: 1.19 $
   \date    $Date: 2010/02/03 18:24:05 $
*/

#include <stdio.h>
#include <math.h>
#include "Tracks.h"
#include "DetectorInfo.h"
#include "Defs.h"
#include "Utils.h"

#include <TROOT.h>
#include <TFile.h>
#include <TStyle.h>
#include <TChain.h>
#include <TBranch.h>
#include <TFile.h>
#include <TMath.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TCanvas.h>
#include <TH1.h>

#include <iostream>

using namespace std;

//_______________________________________________
ClassImp(Tracks)
Tracks::Tracks( bool magnets_on  ):
  TObject(),
  T_align_( 0 ),
  magnets_on_(magnets_on),
  ientry_( 0 ),
  nentries_( 0 ),
  cut_("")
{ }

//_______________________________________________
bool Tracks::GetNextEntry( void )
{ 
  // apply cut
  if( cut_.size() ) 
  while( ientry_ < nentries_ && !T_align_->Draw("T_evt",cut_.c_str(), "goff", 1, ientry_ ) ) 
  ientry_++;

  return GetEntry( ientry_++ ); 
}

//_______________________________________________
bool Tracks::GetEntry( unsigned int i ) 
{

  if( T_align_ && i < T_align_->GetEntries() ) {
    T_align_->GetEntry( i ); 
    return true; 
  }
  
  return false;
}

//_______________________________________________
#ifndef __CINT__
bool Tracks::AcceptEntry( vector<DetectorInfo*> dets, unsigned int nPars, double projRange )
{
  unsigned int nP = 0;
  bool diffProj = false;
  double cosT = 0;
  double sinT = 0;
  double sinProjRange=sin( projRange*PI/180.);
  for( int i=0; i< T_nDets; i++ ) 
  for( unsigned int iD=0; iD<dets.size(); iD++ )
  if( dets[iD]->id_ == T_detVect[i] ) { 
    if( !nP )  { 
      cosT = dets[iD]->cosTheta_;
      sinT = dets[iD]->sinTheta_;
    } else if( (!diffProj) && fabs( sinT*dets[iD]->cosTheta_ -  cosT*dets[iD]->sinTheta_ ) < sinProjRange )
    diffProj = true;
    nP++;
    break;
  }
  
  return ( diffProj && nP>nPars ) ;
}
  
//_______________________________________________
void Tracks::DumpEntry( ostream &out, vector<DetectorInfo*> dets )
{ 
  printf( "Tracks::DumpEntry %5i:\n", ientry_ );
  printf( "position  : %10.4f %10.4f %10.4f.\n", T_xLx, T_yLx, T_zLx );
  printf( "angles    : %10.4f %10.4f.\n", T_txLx, T_tyLx );
  printf( "chi2/prob : %10.4f %10.4f.\n", T_chi2/T_ndf, T_prob );
  printf( "nDetectors: %4i.\n", T_nDets );
  for( int i=0; i< T_nDets; i++ )
  for( unsigned int iD=0; iD<dets.size(); iD++ )
  if( dets[iD]->id_ == T_detVect[i] ) {  
    printf( "%3i Det %9s (%4i) - U/DU/Z: %10.4f %10.4f %10.4f.\n",
      i,
      dets[iD]->TBName_.c_str(), 
      T_detVect[i],
      T_uVect[i],
      T_duVect[i],
      dets[iD]->zcm_ );
    break;
  } 
  
  out << endl;
}
#endif      
    
//_______________________________________________
void Tracks::AddCut( const char* cut )
{
  cout << "Tracks::AddCut - adding \"" << cut << "\".\n";
  if( !cut_.size() ) cut_ = string(cut);
  else cut_+= "&&"+string(cut);
}

//_______________________________________________
unsigned int Tracks::AddToChain( const char* fileSelection, bool isBatch )
{
  
  vector< string > files = Utils::GetFiles( string(fileSelection) );
  
  for( unsigned int i=0; i<files.size(); i++ ) {
    // string file = (isBatch) ? Utils::CopyToLocal( files[i] ):files[i]; // jj comment - do not copy anything!
    string file = files[i];

    // First call
    if( !T_align_ ) {
      
      // open TFile
      TFile *tf = new TFile( file.c_str() );
      if( tf == 0 || !tf->IsOpen() ) {
        cout << "Tracks::AddToChain - ERROR: troubles with TFile \"" << file << "\".\n";
        continue;
      } 
      
      // Get "CsAlignment" Key and check
      TKey* tk =  tf->GetKey("CsAlignment");
      if( !tk ) {
        cout << "Tracks::AddToChain - ERROR: cannot find \"CsAlignment\" in \"" 
          << file 
          << "\".\n";
        continue;
      }
      
      // Get TDirectory and check
      TDirectory *dir = (TDirectory*) tk->ReadObj();
      if( !dir ) {
        cout << "Tracks::AddToChain - ERROR: cannot find directory \"CsAlignment\" in \"" 
          << tf->GetName() 
          << "\".\n";
        continue;
      }
        
      TList* list = dir->GetListOfKeys();
      
      // Check all objects in dir ===
      string tag = "T_align_";
      bool found = false;
      for( TObject* obj = list->First(); obj != 0; obj = list->After( obj ) )  { 
        
        // check if object is a tree
        TTree * t = (TTree*) (  (TKey*) obj)->ReadObj();
        if( !t ) continue;
        
        // check tree name
        string name( obj->GetName() );
     
        size_t pos = name.find( tag );
        if( pos == string::npos ) continue;
        found = true;
        treeName_ = string("CsAlignment/")+name;
        break;
      }  
      
      if( !found ) {
        cout << "Tracks::AddToChain - ERROR: no tree found in \"" << file << "\".\n";
        continue;
      }
      
      // book chain and get branches
      T_align_ = new TChain( treeName_.c_str(), treeName_.c_str() );
      T_align_->Add( file.c_str() );
      if( T_align_->GetBranch( "T_evt"     ) ) T_align_->SetBranchAddress("T_evt",      &T_evt);     else cout << "Tracks::AddToChain branch T_evt not found.\n";    
      if( T_align_->GetBranch( "T_trigMsk" ) ) T_align_->SetBranchAddress("T_trigMsk",  &T_trigMsk); else cout << "Tracks::AddToChain branch T_trigMsk not found.\n";
      if( T_align_->GetBranch( "T_zone"    ) ) T_align_->SetBranchAddress("T_zone",     &T_zone);    else cout << "Tracks::AddToChain branch T_zone not found.\n";   
      if( T_align_->GetBranch( "T_cmlt"    ) ) T_align_->SetBranchAddress("T_cmlt",     &T_cmlt);    else cout << "Tracks::AddToChain branch T_cmlt not found.\n";   
      if( T_align_->GetBranch( "T_chi2"    ) ) T_align_->SetBranchAddress("T_chi2",     &T_chi2);    else cout << "Tracks::AddToChain branch T_chi2 not found.\n";   
      if( T_align_->GetBranch( "T_ndf"     ) ) T_align_->SetBranchAddress("T_ndf",      &T_ndf);     else cout << "Tracks::AddToChain branch T_ndf not found.\n";    
      if( T_align_->GetBranch( "T_prob"    ) ) T_align_->SetBranchAddress("T_prob",     &T_prob);    else cout << "Tracks::AddToChain branch T_prob not found.\n";   
      if( T_align_->GetBranch( "T_cop"     ) ) T_align_->SetBranchAddress("T_cop",      &T_cop);     else cout << "Tracks::AddToChain branch T_cop not found.\n";    
      if( T_align_->GetBranch( "T_meanT"   ) ) T_align_->SetBranchAddress("T_meanT",    &T_meanT);   else cout << "Tracks::AddToChain branch T_meanT not found.\n";   
      if( T_align_->GetBranch( "T_nDets"   ) ) T_align_->SetBranchAddress("T_nDets",    &T_nDets);   else cout << "Tracks::AddToChain branch T_nDets not found.\n";  
      
      // detector-wise branches
      if( T_align_->GetBranch( "T_detVect" ) ) T_align_->SetBranchAddress("T_detVect",  T_detVect);  else cout << "Tracks::AddToChain branch T_detVect not found.\n";
      if( T_align_->GetBranch( "T_uVect"   ) ) T_align_->SetBranchAddress("T_uVect",    T_uVect);    else cout << "Tracks::AddToChain branch T_uVect not found.\n";  
      if( T_align_->GetBranch( "T_vVect"   ) ) T_align_->SetBranchAddress("T_vVect",    T_vVect);    else cout << "Tracks::AddToChain branch T_vVect not found.\n";  
      if( T_align_->GetBranch( "T_duVect"  ) ) T_align_->SetBranchAddress("T_duVect",   T_duVect);   else cout << "Tracks::AddToChain branch T_duVect not found.\n"; 
      if( T_align_->GetBranch( "T_dvVect"  ) ) T_align_->SetBranchAddress("T_dvVect",   T_dvVect);   else cout << "Tracks::AddToChain branch T_dvVect not found.\n"; 
      if( T_align_->GetBranch( "T_rVect"   ) ) T_align_->SetBranchAddress("T_rVect",    T_rVect);    else cout << "Tracks::AddToChain branch T_rVect not found.\n";  
      if( T_align_->GetBranch( "T_tVect"   ) ) T_align_->SetBranchAddress("T_tVect",    T_tVect);    else cout << "Tracks::AddToChain branch T_tVect not found.\n";  
      
      // Branches specific to magnets_on_ = false
      if( !magnets_on_ ) {
        if( T_align_->GetBranch( "T_xLx" ))  T_align_->SetBranchAddress("T_xLx",  &T_xLx);  else cout << "Tracks::AddToChain branch T_xLx not found.\n"; 
        if( T_align_->GetBranch( "T_yLx" ))  T_align_->SetBranchAddress("T_yLx",  &T_yLx);  else cout << "Tracks::AddToChain branch T_yLx not found.\n"; 
        if( T_align_->GetBranch( "T_zLx" ))  T_align_->SetBranchAddress("T_zLx",  &T_zLx);  else cout << "Tracks::AddToChain branch T_zLx not found.\n";
        if( T_align_->GetBranch( "T_txLx" )) T_align_->SetBranchAddress("T_txLx", &T_txLx); else cout << "Tracks::AddToChain branch T_txLx not found.\n";
        if( T_align_->GetBranch( "T_tyLx" )) T_align_->SetBranchAddress("T_tyLx", &T_tyLx); else cout << "Tracks::AddToChain branch T_tyLx not found.\n";
      } else {
        if( T_align_->GetBranch( "T_xVect" ))  T_align_->SetBranchAddress("T_xVect",    T_xVect);  else cout << "Tracks::AddToChain branch T_xVect not found.\n"; 
        if( T_align_->GetBranch( "T_yVect" ))  T_align_->SetBranchAddress("T_yVect",    T_yVect);  else cout << "Tracks::AddToChain branch T_yVect not found.\n"; 
        if( T_align_->GetBranch( "T_zVect" ))  T_align_->SetBranchAddress("T_zVect",    T_zVect);  else cout << "Tracks::AddToChain branch T_zVect not found.\n";
        if( T_align_->GetBranch( "T_txVect" )) T_align_->SetBranchAddress("T_txVect",   T_txVect); else cout << "Tracks::AddToChain branch T_txVect not found.\n";
        if( T_align_->GetBranch( "T_tyVect" )) T_align_->SetBranchAddress("T_tyVect",   T_tyVect); else cout << "Tracks::AddToChain branch T_tyVect not found.\n";
        if( T_align_->GetBranch( "T_dzVect" )) T_align_->SetBranchAddress("T_dzVect",   T_dzVect); else cout << "Tracks::AddToChain branch T_dzVect not found.\n";
        if( T_align_->GetBranch( "T_BVect" ))  T_align_->SetBranchAddress("T_BVect",    T_BVect);  else cout << "Tracks::AddToChain branch T_BVect not found.\n";
      }    
        
    } else T_align_->Add( file.c_str(), -1 ); 
    if( T_align_ ) nentries_ = (unsigned int) T_align_->GetEntries();
    printf("Tracks::AddToChain - file \"%s\" added - %i entrie(s).\n", file.c_str(), nentries_ );
  }
  
  return nentries_;
}

