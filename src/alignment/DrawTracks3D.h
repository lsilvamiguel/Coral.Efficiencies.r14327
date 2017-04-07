// $Id: DrawTracks3D.h,v 1.14 2008/02/20 13:45:47 rgazda Exp $
#ifndef DrawTracks3D_h
#define DrawTracks3D_h

/*!
   \file    DrawTracks3D.h
   \brief   To draw tracks and detectors from alignment output root tree in a 3D display 
   \author  Hugo Pereira
   \version $Revision: 1.14 $
   \date    $Date: 2008/02/20 13:45:47 $
*/

#include <iostream>
#include <vector>
#include <string>
#include <cstdarg>
#include <TROOT.h>
#include <TObject.h>
#include <TGButton.h>
#include <set>
class Tracks;
class TNode;
class TCanvas;
class TRootCanvas;
//class TView; //jj
class DetFileManager;
class Obj3D;
class TGMainFrame;
class TGCheckButton;
class TGTextEntry;
class TGLabel;

/*! 
   \class   DrawTracks3D
   \brief   To draw tracks and detectors from alignment output root tree in a 3D display    
*/

class DrawTracks3D: public TObject {
  public:

  /*! \fn  DrawTracks3D( 
    const char* trackfileselection="", 
    const char* detectorfile="", 
    bool magnets_on=false );
    \brief constructor.
    \param trackfileselection input tree(s) obtained from coral. Wildcards are accepted.
    \param detectorfile the name of the det.dat file to be loaded
    \param magnets_on controls the expected structure of the tree. \warning It must match the option used in coral
  */
  DrawTracks3D( 
    const char* trackfileselection= "", 
    const char* detectorfile= "", 
    bool magnets_on=false );

  ~DrawTracks3D(); //!< destructor
  
  /*! \fn DetFileManager* LoadDetectorFile( const char* name )
     \brief Load detector.dat file, returns a pointer to the corresponding DetFileManager object
     \param name the name of the file to be loaded.
  */
  DetFileManager* LoadDetectorFile( const char* name );

  //! Returns pointer to the DetFileManager object, if any, 0 otherwise
  DetFileManager* GetDetectorFile( void ) const { return df_;}
  
  /*! \fn virtual void LoadTracks(  const char* fileselection, bool magnets_on = false )
     \brief load tracks contained in coral alignment tree output, returns pointer to the corresponding Tracks object.
     \param fileselection input tree(s). wildcards are accepted. Trees are chained.
     \param magnets_on controls the expected structure of the tree. \warning It must match the option used in coral
  */
  void LoadTracks( const char* fileselection, bool magnets_on = false );
  
  /*! \fn virtual void AddToTracks( const char* fileselection )
     \brief add tracks contained in coral alignment tree to existing Tracks object, if any. Returns pointer to tracks object.
     \param fileselection input tree(s). wildcards are accepted. Trees are chained.
     \warning added trees must be of same type
  */   
  void AddToTracks( const char* fileselection );
  
  //! Returns pointer to the Tracks object, if any, 0 otherwise.
  Tracks* GetTracks( void ) const { return tracks_; }
   
  /*! \fn void DrawDetectors( const char* detselection = 0 )
     \brief draw detectors in 3D Display using DetFileManager infos.
     \param detselection detector selection, using shorten names, character wise wild cards example "GM01 ST03Y2*b DC01X1__"
     default corresponds to all detectors except BMS and Veto
  */
  void DrawDetectors( const char* detselection = 0 );
  
  /*! \fn Obj3D* DrawTrack( int itr, bool removeOld = true )
     \brief draw a single track read from Tracks object in 3D display. 
     Returns false if track was not drawn for any reason. True otherwise
     \param itr track ID in tree. -1 means next track. Only -1 implemented when running with coral
     \param removeOld if true, delete previously drawn tracks
  */
  bool DrawTrack( int itr = -1, bool removeOld = true );
  void DrawEvent( void ); //!< draw all tracks belonging to same event
  
  void ToggleTrackWires( void ); //!< draw/hide wires in tracks according to check button state
  void ToggleHitWires( void );   //!< draw/hide wires hit in event according to check button state
  void ToggleMagnets( void );    //!< draw/hide magnets according to check button state
  void Exit( void );             //!< exit!
 
  void SelectProjections( void );             //!< change the value for projFlags_ TGCheckButtons
  void WriteToLog( const char* format, ... ); //!< to write online comments to log_ TGLabel
  
  //  TView* GetView( void ) const { return view_; }   //!< returns pointer to TView object (jj
  TCanvas* GetCanvas( void ) const { return cv_; } //!< returns pointer to TCanvas object
 
  private:
  
  //! Returns true if checkbutton is activated
  bool _IsDown( TGCheckButton* btn ) { return ( btn && btn->GetState() == kButtonDown ); }
  void _Init( void );  //!< Initialize all toggle buttons values 
  unsigned int _GetNMatchingDets( void ); //!< returns number of matching detectors in current track
  bool _MakeMainMenu( void );       //!< Makes main menu window returns true if everything was fine
  bool _MakeProjectionMenu( void ); //!< Makes projections menu window returns true if everything was fine
  DetFileManager *df_;     //!< pointer to DetFileManagerObject if any
  Tracks *tracks_;         //!< pointer to Tracks object, if any. (used only when reading alignment tree)
  TCanvas *cv_;            //!< the canvas

  //  TView *view_;        // jj (see $0.cc)
  TNode  *mainNode_;       //!< pointer to 3D Display Main Node
  TNode  *rotNode_;        //!< pointer to 3D Display Rotated Node (z horizontal)
  TGMainFrame *mainFrame_; //!< main frame for menu display
  TGMainFrame *projFrame_; //!< main frame for menu display
  
  std::set< int > drawnDetIds_;             //!< map of drawn detector ids
  std::vector< TGCheckButton* > projFlags_; //!< true/false to show/hide detectors corresponding to a given projection
  std::vector< Obj3D* > detectors3D_;  //!< vector of pointers to the currently drawn detectors
  std::vector< Obj3D* > deadZones3D_;  //!< vector of pointers to the currently drawn dead zones
  std::vector< Obj3D* > tracks3D_;     //!< vector of pointers to the currently drawn tracks
  std::vector< Obj3D* > tWires3D_;     //!< vector of pointers to wires belonging to tracks 
  std::vector< Obj3D* > hWires3D_;     //!< vector of pointers to wires hit in event (meaningfull only when running with coral)
  std::vector< std::string > trackFileSelection_;  //!< selected files from which tracks have been loaded
  std::string detectorFile_; //!< name of the detector.dat file
  int eventId_;         //!< current event id_
  int trackId_;         //!< current track id_

  double detZMin_;   //!< first detector position along the beam
  double detZMax_;   //!< last detector position along the beam

  TGCheckButton* drawTWires_;      //!< to tell if wires belonging to tracks are to be drawn
  TGCheckButton* drawHWires_;      //!< to tell if all hit wires in an event are to be drawn
  TGCheckButton* drawMagnets_;     //!< to tell if magnets are to be drawn
  TGCheckButton* bridgedOnly_;     //!< to tell if only tracks with momentum are to be drawn
  TGLabel* log_;        //!< used to write comments while processing
  std::string detselection_; //!< selected detector string default is "FI SI MM GM DC ST DW P MA MB H"
  TGTextEntry * detSelEntry_; //!< TGObject to change detector selection
  
  static const int tWireColor  = 6; //!< wire objects color (wires in tracks)
  static const int hWireColor  = 4; //!< wire objects color (wires hit in events)
  static const int trackColor  = 2; //!< track objects color
  static const int detColor    = 1; //!< detector objects color
  static const int magnetColor = 2; //!< magnet objects color

  static const int wireWidth   = 1; //!< wire objects line width
  static const int detWidth    = 1; //!< track objects line width
  static const int magnetWidth = 1; //!< detector objects line width
  static const int trackWidth  = 1; //!< magnet objects line width
  
  ClassDef(DrawTracks3D,0)              //!< Root Macro

};


#endif
