// $Id: DrawTracks3D.cc,v 1.20 2008/06/05 13:00:54 rgazda Exp $

/*!
   \file    DrawTracks3D.cc
   \brief   To draw tracks and detectors from alignment output root tree in a 3D display 
   \author  Hugo Pereira
   \version $Revision: 1.20 $
   \date    $Date: 2008/06/05 13:00:54 $
*/

#include "DrawTracks3D.h"

#include "DetFileManager.h"
#include "DetectorInfo.h"
#include "DeadZoneInfo.h"
#include "MagnetInfo.h"
#include "Utils.h"
#include "Obj3D.h"
#include "Point.h"

#include <TNode.h>
//#include <TView.h> //jj - in root 5.16 this class was changed (but it is not used anyway anywhere).
#include <TBRIK.h>
#include <TRootCanvas.h>
#include <TCanvas.h>
#include <TGFrame.h>
#include <TGMenu.h>
#include <TGTextEntry.h>
#include <TGLabel.h>
#include <vector>
using namespace std;
static const int SelProjAll  = 0;
static const int SelProjNone = 1;
static const int SelProjInv  = 2;

//________________________________________
ClassImp(DrawTracks3D)
DrawTracks3D::DrawTracks3D( 
  const char* trackfileselection, 
  const char* detectorfile, 
  bool magnets_on ):
  TObject(),
  df_(0),
  tracks_(0),
  mainNode_( 0 ),
  mainFrame_( 0 ),
  projFrame_( 0 ),
  detectorFile_( "" ),
  eventId_( -1 ),
  trackId_( 0 ),
  
  detZMin_( 0 ),
  detZMax_( 0 ),
  
  // TGObjects
  drawTWires_( 0 ),
  drawHWires_( 0 ),
  drawMagnets_( 0 ),
  bridgedOnly_( 0 ),
  
  log_( 0 ),
  detselection_("FI SI MM GM DC ST DW P MA MB H"),  // list of detectors to be drawn by default  
  detSelEntry_( 0 )

{
  if( detectorfile && strlen( detectorfile ) ) LoadDetectorFile( detectorfile );
  if( trackfileselection && strlen( trackfileselection ) ) LoadTracks( trackfileselection, magnets_on ); 
  
  gROOT->SetStyle("Plain");
  
  // Create canvas
  cv_ = new TCanvas("cv","Compass spectrometer",200,10,1000,500);
  //  view_ = new TView(2); // jj - see comment above
  
  // Create Main shape/node
  TBRIK* mainShape = new TBRIK( "mainS", "mainS", "void", 0.1, 0.1, 0.1 );
  mainNode_  = new TNode( "mainN", "mainN", mainShape );
  mainNode_->SetLineColor( 0 );
  mainNode_->cd();
  
  // Create rotated shape/node
  TBRIK* rotShape = new TBRIK( "rotS", "rotS", "void", 0.1, 0.1, 0.1 );
  rotNode_  = new TNode( "rotN", "rotN", rotShape );
  
  double rot[] = {
     0, 1, 0, 
     0, 0, 1, 
     1, 0, 0 };
  TRotMatrix *rotM = new TRotMatrix( "rotM", "rotM", rot );
  rotNode_->SetMatrix( rotM );
  rotNode_->SetLineColor( 0 );
  rotNode_->cd();

  //! make menu window
  _MakeMainMenu();
  _Init();
  WriteToLog( "Ready.");
}

//________________________________________
DrawTracks3D::~DrawTracks3D( void )
{ 
  if( mainFrame_ ) SafeDelete( mainFrame_ );
  if( projFrame_ ) SafeDelete( projFrame_ );
}

//________________________________________
DetFileManager* DrawTracks3D::LoadDetectorFile( const char* name )
{
  if( df_ ) SafeDelete(df_);
  df_ = new DetFileManager( name );
  detectorFile_ = string(name);
  if( projFrame_ ) SafeDelete( projFrame_ );
  _MakeProjectionMenu();
  return df_;
}
    
//________________________________________
void DrawTracks3D::DrawDetectors( const char* detselection )
{  
  if( !df_ ) { cout << "DrawTracks3D::DrawDetectors - no detectorfile loaded.\n"; return; }
  if( !detselection ) {
    if( detSelEntry_ ) detselection_ = string( detSelEntry_->GetText() );
    detselection = detselection_.c_str();
  }
 
  // get detector selection
  vector< DetectorInfo* > dets = df_->GetDetSelection( detselection );
  if( !dets.size() ) {
    cout << "DrawTracks3D::DrawDetectors - no detectors selected.\n";
    return;
  } 
  
  // remove already drawn detectors
  for( unsigned int i=0; i< detectors3D_.size(); i++ )
  if( detectors3D_[i] ) detectors3D_[i]->DeleteTObjects();
  detectors3D_.clear();
  drawnDetIds_.clear();

  // remove already drawn dead zones
  for( unsigned int i=0; i< deadZones3D_.size(); i++ )
  if( deadZones3D_[i] ) deadZones3D_[i]->DeleteTObjects();
  deadZones3D_.clear();
  

  // remove drawn (track) wires which do not correspond to detector to be drawn
  vector< Obj3D* > tWiresTmp;
  for( unsigned int i=0; i< tWires3D_.size(); i++ ) {
    if( !tWires3D_[i] ) continue;
    bool found = false;
    for( unsigned int j=0; j< dets.size(); j++ )
    if( tWires3D_[i]->TBName_.substr( 0, 8 ) == dets[j]->TBName_ ) { found = true; break; }
    if( !found ) { 
      tWires3D_[i]->DeleteTObjects(); 
    } else tWiresTmp.push_back( tWires3D_[i] );
  }

  // update Wires
  tWires3D_.clear();
  for( unsigned int i=0; i< tWiresTmp.size(); i++ ) tWires3D_.push_back( tWiresTmp[i] );
 
  // remove drawn (hit) wires which do not correspond to detector to be drawn
  vector< Obj3D* > hWiresTmp;
  for( unsigned int i=0; i< hWires3D_.size(); i++ ) {
    if( !hWires3D_[i] ) continue;
    bool found = false;
    for( unsigned int j=0; j< dets.size(); j++ )
    if( hWires3D_[i]->TBName_.substr( 0, 8 ) == dets[j]->TBName_ ) { found = true; break; }
    if( !found ) { 
      hWires3D_[i]->DeleteTObjects(); 
    } else hWiresTmp.push_back( hWires3D_[i] );
  }
 
  // update Wires
  hWires3D_.clear();
  for( unsigned int i=0; i< hWiresTmp.size(); i++ ) hWires3D_.push_back( hWiresTmp[i] );
  
  // update detselection_
  detselection_ = string( detselection );
  
  detZMin_ = 0;
  detZMax_ = 0;
  for( unsigned int i=0; i< dets.size(); i++ ) {
    DetectorInfo *det = dets[i]->GetMain();
    
    int projection = det->projection_;
    if( projection >= 0 && 
      projection < (int) projFlags_.size() && 
      !_IsDown( projFlags_[projection] ) ) continue;
    Obj3D* active = 0;
    
    if( det->zcm_ < detZMin_ || !i ) detZMin_ = det->zcm_;
    if( det->zcm_ > detZMax_ || !i ) detZMax_ = det->zcm_;
    
    // Draw Active area
    if( ( active = det->GetArea3D() ) ) {
      if( !active->Drawn() ) active->MakeNodes( rotNode_, detColor, detWidth );
      if( drawnDetIds_.find( det->id_ ) == drawnDetIds_.end() ) drawnDetIds_.insert( det->id_ );
      detectors3D_.push_back( active );
    }

    
    // Draw DeadZone
    DeadZoneInfo *dead = det->dead_;
    if( dead && (active = dead->GetArea3D() ) ) {
      active->MakeNodes( rotNode_, detColor, detWidth  );
      deadZones3D_.push_back( active );
    }
    
    // redraw wires, if any
    for( unsigned int j=0; j< tWires3D_.size(); j++ ) 
    if( tWires3D_[j] && 
      tWires3D_[j]->TBName_.substr(0,det->TBName_.size() ) == det->TBName_ && 
      !tWires3D_[i]->Drawn() ) 
    tWires3D_[j]->MakeNodes( rotNode_ );
    
    // draw subdetectors   
    for( unsigned int j=0; j<dets[i]->sub_.size(); j++ )
    if( dets[i]->sub_[j] != det && ( active = dets[i]->sub_[j]->GetArea3D() ) ) {
      if( !active->Drawn() ) active->MakeNodes( rotNode_, detColor, detWidth );
      detectors3D_.push_back( active );
    
      // redraw wires, if any
      DetectorInfo *sub = dets[i]->sub_[j];
      for( unsigned int k=0; k< tWires3D_.size(); k++ ) 
      if( tWires3D_[k] && tWires3D_[k]->TBName_.substr(0,sub->TBName_.size() ) == det->TBName_ && !tWires3D_[k]->Drawn()) 
      tWires3D_[k]->MakeNodes( rotNode_ );
    }
    
  }
  
  rotNode_->SetPosition( -0.5*(detZMin_+detZMax_), 0, 0 );
  mainNode_->Draw();
  cv_->Update();
  return;
 
}

//________________________________________
void DrawTracks3D::ToggleMagnets( void )
{
  if( !df_ ) { cout << "DrawTracks3D::ToggleMagnets - no detectorfile loaded.\n"; return; } 
  
  TGCheckButton *btn = (TGCheckButton *) gTQSender;
  vector< MagnetInfo* > mags = df_->GetMagnetInfo();
  for( unsigned int i=0; i< mags.size(); i++ ) {
    Obj3D* active = 0;
    if( !( active = mags[i]->GetArea3D() ) ) continue;
    else if( _IsDown( btn ) && !active->Drawn() ) active->MakeNodes( rotNode_, magnetColor, magnetWidth );
    else if( active->Drawn() )                    active->DeleteTObjects();
  }
  
  mainNode_->Draw();
  cv_->Update();
}

//________________________________________
void DrawTracks3D::ToggleTrackWires( void )
{
  TGCheckButton *btn = (TGCheckButton *) gTQSender;
  
  for( unsigned int i=0; i< tWires3D_.size(); i++ ) 
  if( !tWires3D_[i] ) continue;
  else if( _IsDown( btn) && !tWires3D_[i]->Drawn() ) tWires3D_[i]->MakeNodes( rotNode_, tWireColor, wireWidth );
  else if( tWires3D_[i]->Drawn() )                   tWires3D_[i]->DeleteTObjects();

  mainNode_->Draw();
  cv_->Update();
  
 return;
}

//________________________________________
void DrawTracks3D::ToggleHitWires( void )
{
  TGCheckButton *btn = (TGCheckButton *) gTQSender;
  
  for( unsigned int i=0; i< hWires3D_.size(); i++ ) 
  if( !hWires3D_[i] ) continue;
  else if( _IsDown( btn) && !hWires3D_[i]->Drawn() ) hWires3D_[i]->MakeNodes( rotNode_, hWireColor, wireWidth );
  else if( hWires3D_[i]->Drawn() )                   hWires3D_[i]->DeleteTObjects();

  mainNode_->Draw();
  cv_->Update();
  
 return;
}
 
//________________________________________
void DrawTracks3D::SelectProjections( void )
{ 
  TGButton *btn = (TGButton *) gTQSender;
  int id = btn->WidgetId();
  for( unsigned int i=0; i<projFlags_.size(); i++ ) {
    if( (!_IsDown( projFlags_[i] )) && ( id == SelProjAll || id == SelProjInv ) ) projFlags_[i]->SetState( kButtonDown);
    else if( id == SelProjNone || id == SelProjInv ) projFlags_[i]->SetState( kButtonUp);
  }
     
  return;
}

//________________________________________________________________
void DrawTracks3D::WriteToLog( const char* format, ... )
{
  if( !log_ ) return;
  char* buf = new char[256];
  va_list p;
  va_start(p,format);
  vsprintf(buf, format, p);
	va_end(p);
  log_->SetText( buf );
  SafeDelete( buf );
  return;
}
  
//_________________________________________ 
bool DrawTracks3D::_MakeMainMenu( void )
{
  mainFrame_ = new TGMainFrame( gClient->GetRoot(), 100, 200 );  
  TGVerticalFrame  *menuFrame = new TGVerticalFrame( mainFrame_, 100, 200 );
  
  TGTextButton *button;

  TGLayoutHints *hint  = new TGLayoutHints( kLHintsExpandX, 0,1,0,1);                // justified
  TGLayoutHints *hintl = new TGLayoutHints( kLHintsLeft|kLHintsExpandX, 0,1,0,1);    // aligned left, justified
  TGLayoutHints *hintc = new TGLayoutHints( kLHintsCenterX|kLHintsExpandX, 0,1,0,1); // aligned center, justified
  
  // to draw a set of detectors
  TGLabel *label = new TGLabel( menuFrame, "Select the detectors to be drawn");
  menuFrame->AddFrame( label, hintl );
  
  TGTextEntry *entry = new TGTextEntry( menuFrame, detselection_.c_str() );
  entry->Connect("ReturnPressed()","DrawTracks3D",this,"DrawDetectors()" );
  menuFrame->AddFrame( entry, hint );
  detSelEntry_ = entry;
   
  // To draw a set of detectors
  button = new TGTextButton(menuFrame, "&Draw Detectors" ); 
  button->Connect("Clicked()","DrawTracks3D",this,"DrawDetectors()" );
  menuFrame->AddFrame( button, hint );
    
  button = new TGTextButton(menuFrame, "Next &Track", false );    
  button->Connect("Clicked()","DrawTracks3D",this,"DrawTrack()" );
  menuFrame->AddFrame( button, hint );
  
  // to draw next event
  button = new TGTextButton(menuFrame, "Next &Event");     
  button->Connect("Clicked()","DrawTracks3D",this,"DrawEvent()");
  menuFrame->AddFrame( button, hint );
  
  // to draw/hide magnets
  TGCheckButton *check;
  check = new TGCheckButton(menuFrame, "Draw M&agnets");     
  check->Connect("Clicked()","DrawTracks3D",this,"ToggleMagnets()");
  menuFrame->AddFrame( check, hint );
  drawMagnets_ = check;
   
  // to draw/hide wires on track
  check = new TGCheckButton(menuFrame, "W&ires in Tracks" );     
  check->Connect("Clicked()","DrawTracks3D",this,"ToggleTrackWires()");
  menuFrame->AddFrame( check, hint );
  drawTWires_ = check;
   
  // to draw/hide wires hit in event
  check = new TGCheckButton(menuFrame, "W&ires hit" );     
  check->Connect("Clicked()","DrawTracks3D",this,"ToggleHitWires()");
  menuFrame->AddFrame( check, hint );
  drawHWires_ = check;
  
  // to draw only bridget/all tracks
  check = new TGCheckButton(menuFrame, "&Bridged Tracks Only" );     
  menuFrame->AddFrame( check, hint );
  bridgedOnly_ = check;
  
  // to exit application
  button = new TGTextButton(menuFrame, "E&xit" );    
  button->Connect("Clicked()","DrawTracks3D",this,"Exit()" );
  menuFrame->AddFrame( button, hint );
  
  // to display online comments
  label = new TGLabel( menuFrame, "" );
  menuFrame->AddFrame( label, hintc );
  log_ = label; WriteToLog( "Initialising...");
  
  menuFrame->MapSubwindows();
  menuFrame->Resize(menuFrame->GetDefaultSize());
  menuFrame->MapWindow();    

  mainFrame_->SetWindowName( "menu" );
  mainFrame_->MapSubwindows();
  mainFrame_->Resize(menuFrame->GetDefaultSize());
  mainFrame_->MapWindow();    
  
  return true;
}

//_________________________________________ 
bool DrawTracks3D::_MakeProjectionMenu( void )
{
  if( !df_ ) return false;
  vector< int > projections = df_->GetProjections();
  if( !projections.size() ) return false;
  
  projFrame_ = new TGMainFrame( gClient->GetRoot(), 100, 200 );  
  TGVerticalFrame  *menuFrame = new TGVerticalFrame( projFrame_, 100, 200 );
  
  TGCheckButton *check;
  TGLayoutHints *hint = new TGLayoutHints( kLHintsExpandX, 0,1,0,1);

  for( unsigned int i=0; i< projections.size(); i++ ) {
  
    // to draw/hide detectors corresponding to projection i
    char* buf = new char[64];
    sprintf( buf, "%i degrees",projections[i] );
    check = new TGCheckButton(menuFrame, buf );
    SafeDelete( buf );
    
    check->SetState( kButtonDown );
    check->Connect("Clicked()","DrawTracks3D",this,"DrawDetectors()");
    menuFrame->AddFrame( check, hint );    
    projFlags_.push_back( check );
  }  
  
  TGTextButton *button;
  
  // to select all projections
  button = new TGTextButton(menuFrame, "      &All      ", SelProjAll );    
  button->Connect("Clicked()","DrawTracks3D",this,"SelectProjections()" );
  menuFrame->AddFrame( button, hint );
  
  // to deselect all projections
  button = new TGTextButton(menuFrame, "      &None     ", SelProjNone );    
  button->Connect("Clicked()","DrawTracks3D",this,"SelectProjections()" );
  menuFrame->AddFrame( button, hint );
  
  // to invert selection all projections
  button = new TGTextButton(menuFrame, "     &Invert    ", SelProjInv );    
  button->Connect("Clicked()","DrawTracks3D",this,"SelectProjections()" );

  menuFrame->AddFrame( button, hint );
  menuFrame->MapSubwindows();
  menuFrame->Resize(menuFrame->GetDefaultSize());
  menuFrame->MapWindow();    

  projFrame_->SetWindowName( "Projections" );
  projFrame_->MapSubwindows();
  projFrame_->Resize(menuFrame->GetDefaultSize());
  projFrame_->MapWindow();    
  
  return true;
}
  
