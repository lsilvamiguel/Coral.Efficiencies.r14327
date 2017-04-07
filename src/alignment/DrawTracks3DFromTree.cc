// $Id: DrawTracks3DFromTree.cc,v 1.10 2008/06/05 12:59:21 rgazda Exp $

/*!
   \file    DrawTracks3DFromTree.cc
   \brief   The part of DrawTracks3D which is specific to tracks loaded from alignment output tree
   \author  Hugo Pereira
   \version $Revision: 1.10 $
   \date    $Date: 2008/06/05 12:59:21 $
*/

#include "DrawTracks3D.h"
#include "DetFileManager.h"
#include "DetectorInfo.h"
#include "Tracks.h"
#include "Obj3D.h"
#include <TCanvas.h>
#include <TNode.h>
#include <TGTextEntry.h>

using namespace std;

//________________________________________
void  DrawTracks3D::LoadTracks( const char* fileselection, bool magnets_on )
{
  if( tracks_ ) SafeDelete(tracks_);
  cout << "DrawTracks3D::LoadTracks - magnets: " << ((magnets_on)?"on":"off") << ".\n";
  tracks_ = new Tracks( magnets_on );  
  tracks_->AddToChain( fileselection );
  trackFileSelection_.clear();
  trackFileSelection_.push_back(fileselection);
  return;
}

//________________________________________
void DrawTracks3D::AddToTracks( const char* fileselection )
{
  if( !tracks_ ) {
    cout << "DrawTracks3D::AddToTracks - Tracks object is not defined.\n";
    return;
  }
  tracks_->AddToChain( fileselection );
  
  trackFileSelection_.push_back(fileselection);
  return;
}


//________________________________________
bool DrawTracks3D::DrawTrack( int itr, bool removeOld )
{    
  if( !tracks_ ) { 
    cout << "DrawTracks3D::DrawTrack - ERROR: no tracks loaded.\n"; 
    return false; 
  }
  
  // check if detectors have to be (re)drawn
  if( detSelEntry_ ) {
    string detselection = string( detSelEntry_->GetText() );
    if( detselection != detselection_ ) DrawDetectors();
  } else if( !drawnDetIds_.size() ) DrawDetectors();
      
  if( removeOld ) {
    bool erased = false;
    // remove drawn tracks
    for( unsigned int i=0; i < tracks3D_.size(); i++ ) 
    if( tracks3D_[i] ) { erased = true; SafeDelete( tracks3D_[i] ); }
    tracks3D_.clear();
    
    // remove drawn wires (nodes and shapes only)
    for( unsigned int i=0; i < tWires3D_.size(); i++ ) 
    if( tWires3D_[i] ) { erased = true; tWires3D_[i]->DeleteTObjects(); }
    tWires3D_.clear(); 
    
    // redraw if something was erased
    if( erased ) { mainNode_->Draw(); cv_->Update(); }
    
  }
 
  // Load track parameters
  if( itr >=0 && 
    ( (!tracks_->GetEntry( (unsigned int) itr )) || _GetNMatchingDets() < 2 ) ) {
    if( removeOld ) WriteToLog( "track %i is empty", itr );
    return false; 
  } else {
    bool found = false;
    while( (!found) && tracks_->GetNextEntry() ) found = ( _GetNMatchingDets() >= 2 );
    if( !found ) {
      WriteToLog( "no more tracks" );
      return false;
    }
    itr = int( tracks_->GetIEntry() );
  }
  
  if( _IsDown( bridgedOnly_ ) && !tracks_->T_cop ) { 
    if( removeOld ) WriteToLog( "track %i not bridged", itr );
    return false; 
  }
  
  // create PolyLine3D object to store track points     
  char* buf = new char[64];
  sprintf( buf, "Track_%i", (itr >= 0)? itr:int(tracks_->GetIEntry()) );
  PolyLine3D *pline = new PolyLine3D( string( buf ) );
  SafeDelete( buf );
      
  int nFiredDet = (int) tracks_->T_nDets;  // number of detectors in track
  
  // Store the position at which the first new wire will be booked
  int newWire = tWires3D_.size();  

  // magnets are off    
  if( !tracks_->MagnetsOn() ) { 
    DetectorInfo *first = 0;
    DetectorInfo *last = 0;
    for( int i=0; i<nFiredDet; i++ ) {
      DetectorInfo *det = df_->GetDetInfo( tracks_->T_detVect[i] );
      if( !det || drawnDetIds_.find( det->id_ ) == drawnDetIds_.end() ) continue;
      if( !first ) first = det;
      last = det;
      Obj3D* wire = det->GetWire3D( det->UtoWire( tracks_->T_uVect[i] ) );
      if( wire ) tWires3D_.push_back( wire );
    }
    
    pline->AddPoint( Point( 
      tracks_->T_xLx + tracks_->T_txLx*(first->zcm_ - tracks_->T_zLx ),
      tracks_->T_yLx + tracks_->T_tyLx*(first->zcm_ - tracks_->T_zLx ),
      first->zcm_ ) );
    
    pline->AddPoint( Point( 
      tracks_->T_xLx + tracks_->T_txLx*(last->zcm_ - tracks_->T_zLx ),
      tracks_->T_yLx + tracks_->T_tyLx*(last->zcm_ - tracks_->T_zLx ),
      last->zcm_ ) );
  
  // magnets are on    
  } else {
    for( int i=0; i<nFiredDet; i++ ) {
      DetectorInfo *det = df_->GetDetInfo( tracks_->T_detVect[i] );
      if( !det || drawnDetIds_.find( det->id_ ) == drawnDetIds_.end() ) continue;
      
      pline->AddPoint( Point( 
        tracks_->T_xVect[i] + tracks_->T_txVect[i]*(det->zcm_ - tracks_->T_zVect[i] ),
        tracks_->T_yVect[i] + tracks_->T_tyVect[i]*(det->zcm_ - tracks_->T_zVect[i] ),
        det->zcm_ ) );

      Obj3D* wire = det->GetWire3D( det->UtoWire( tracks_->T_uVect[i] ) );
      if( wire ) tWires3D_.push_back( wire );
    }  
  }        
  
  // Draw Polyline3D  
  pline->MakeNodes( rotNode_, trackColor, trackWidth );
  tracks3D_.push_back( (Obj3D*) pline );
  
  // Draw new wires if required
  if( _IsDown(drawTWires_) )
  for( unsigned int i=newWire; i < tWires3D_.size(); i++ )
  if( tWires3D_[i] && !tWires3D_[i]->Drawn() ) 
  tWires3D_[i]->MakeNodes( rotNode_, tWireColor, wireWidth);
    
  if( removeOld ) {
    WriteToLog( "track %i", itr );
    mainNode_->Draw();
    cv_->Update();
  }
  return true;
}

//________________________________________
void DrawTracks3D::DrawEvent( void )
{  
  
  if( !drawnDetIds_.size() ) DrawDetectors();
        
  // remove drawn tracks
  for( unsigned int i=0; i < tracks3D_.size(); i++ ) 
  if( tracks3D_[i] ) SafeDelete( tracks3D_[i] );
  tracks3D_.clear();
  
  // remove drawn wires (nodes and shapes only)
  for( unsigned int i=0; i < tWires3D_.size(); i++ ) 
  if( tWires3D_[i]  ) tWires3D_[i]->DeleteTObjects(); 
  tWires3D_.clear(); 
 
  // Load first track if needed
  if( eventId_ < 0 ) {
    // bool found = false; // unused variable (jj)
    if( !tracks_->GetEntry( trackId_ ) ) { 
      cout << "DrawTracks3D::DrawDet - no more tracks.\n"; 
      return;
    }
    eventId_ = tracks_->T_evt;
  }
  
  // Load all tracks matching event eventId_
  int nTracksDrawn = 0;  
  while( int(tracks_->T_evt) == eventId_ ) {
    if( _GetNMatchingDets() >= 2 && DrawTrack( trackId_, false ) ) nTracksDrawn++;
    trackId_++;
    if( !tracks_->GetEntry( trackId_ ) ) { 
      cout << "DrawTracks3D::DrawDet - no more tracks.\n"; 
      WriteToLog("no more tracks");
      return;
    }
  }

  if( nTracksDrawn )
  WriteToLog( "evt %i: %i track%s", 
    eventId_, 
    nTracksDrawn,
    (nTracksDrawn > 1)? "s":"" );
  else WriteToLog( "evt %: empty" , eventId_);
  
  eventId_ = tracks_->T_evt;
  mainNode_->Draw();
  cv_->Update();
  return;
}

//_________________________________________ 
void DrawTracks3D::Exit( void ) { exit(0); }

//________________________________________
void  DrawTracks3D::_Init( void )
{
 
  if( drawTWires_ ) drawTWires_->SetState( kButtonDown ); // wires in track
  
  // wires hit in event (disabled when using Tree)
  if( drawHWires_ ) drawHWires_->SetState( kButtonUp ); 
   
  if( bridgedOnly_ ) bridgedOnly_->SetState( kButtonUp );   // tracks with momentum 
  if( drawMagnets_ ) drawMagnets_->SetState( kButtonDown ); // magnets 
  
  // Draw all required detectors
  DrawDetectors();
  return;
}   

//_________________________________________ 
unsigned int DrawTracks3D::_GetNMatchingDets( void )
{
  if( !tracks_ ) return 0;
  unsigned int out = 0;
  for( int i = 0; i < tracks_->T_nDets; i++ ) 
  if( drawnDetIds_.find( tracks_->T_detVect[i] ) != drawnDetIds_.end() ) 
  out ++;
  return out;
}
