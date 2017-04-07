// $Id: DetFileManager.h,v 1.22 2008/03/05 20:36:26 ybedfer Exp $

/*!
   \file    DetFileManager.h
   \brief   read, manipulate, format detectors and dead zones in det.dat file
   \author  Hugo Pereira
   \version $Revision: 1.22 $
   \date    $Date: 2008/03/05 20:36:26 $
*/

#ifndef DetFileManager_h
#define DetFileManager_h
#include <stdio.h>
#include <math.h>
#include <string>
#include <vector>
#include <TROOT.h>
#include <TObject.h>
#include <TMath.h>
#include <functional>

#include "DetectorInfo.h"
#include "DeadZoneInfo.h"
#include "MagnetInfo.h"

class TNode;

/*!
   \class   DetFileManager
   \brief   read, manipulate, format detectors and dead zones in det.dat file
   Fonctionalities, amongst other are: 
   To make X and Y (U and V) detector centers coincide (by translating along the wires)
   To make detector and deadzone centers coincide (with offsets if needed)
   To convert first wire position to the center position
   To update all detectors using the results of aligmment
*/
class DetFileManager: public TObject {

  public:
  /*! \fn DetFileManager( const char* detFileName )
    \brief constructor
    \param detFileName, the detector.dat file name
  */
  DetFileManager( const char* detFileName ); //!< reading from ascii file
  DetFileManager( void );  //!< reading from a root file
  void Init();                  //!< parsing of geometry file. called by ctors

  /*! \fn DetectorInfo* GetDetInfo( const char* selection, bool quiet = false )
    \brief returns pointer to first DetectorInfo object matching detselection 
    \param selection TBName selection string: example: "St03X1*b", "GM01X1__"
  */
  DetectorInfo* GetDetInfo( const char* selection, bool quiet = false );
  
  /*! \fn DetectorInfo* GetDetInfo( int id )
    \brief returns pointer to the DetectorInfo object corresponding to given id
    \param id
    returns 0 if no matching detector is found
  */
  DetectorInfo* GetDetInfo( int id );  
  
  //! Cleaning macro to set all centers together, move the dead zones, sort the detectores ... Dangerous.
  void CleanDetFile( void ); 
  
  /*! \fn bool UpdateFromAlign( const char* file )
    \brief update all DetectorInfo according to alignment output file
    \param file the name of the alignment output file
    returns true if everything is fine.
  */
  bool UpdateFromAlign( const char* file );  
  
  /*! \fn bool SameZ( const char* detFileName, const char* selection="*")
    \brief change zcm_ of all detectors matching selection  according to the values stored in detFileName
    \param selection detector selection
    \param detfilename detector.dat file
    TBName is used for matching. Special care is taken of MM and ST detectors, which have subdetectors
    returns true if everything is fine.
  */
  bool SameZ( const char* detFileName, const char* selection="*");
  
  /*! \fn bool SameSpSli( const char* detFileName, const char* selection="*")
    \brief change spSli_ of all detectors matching selection  according to the values stored in detFileName
    \param selection detector selection
    \param detfilename detector.dat file
    TBName is used for matching. Special care is taken of MM and ST detectors, which have subdetectors
    returns true if everything is fine.
  */
  bool SameSpSli( const char* detFileName, const char* selection="*");
   
  /*! \fn bool SameVel( const char* detFileName, const char* selection="*")
    \brief change vel_ of all detectors matching selection  according to the values stored in detFileName
    \param selection detector selection
    \param detfilename detector.dat file
    TBName is used for matching. Special care is taken of MM and ST detectors, which have subdetectors
    returns true if everything is fine.
  */
  bool SameVel( const char* detFileName, const char* selection="*");
   
  /*! \fn bool SameNames( const char* detFileName, const char* selection="*")
    \brief change name, unit and id of all detectors matching selectionaccording to the values stored in detfilename detector.dat
    \param selection detector selection
    \param detfilename detector.dat file
    TBName is used for matching. Special care is taken of MM and ST detectors, which have subdetectors
    returns true if everything is fine.
  */
  bool SameNames( const char* detFileName, const char* selection="*");

  /*! \fn bool FixWirD( const char* detSelection, double value )
    \brief fix the first wire position of the detectors matching detselection. Move the center position consequently.
    \param detselection detector selection
    \param value the new wire position (mm)
    returns true if everything is fine.
  */
  bool FixWirD( const char* detSelection, double value );
  
  /*! \fn bool MatchCenters( const char* det1, const char* det2 )
    \brief translate detector of TBName det1 and detector of TBName det2 
    along their wire so that their center coincide
    \param det1 TBName of first detector
    \param det2 TBName of second detector
    returns true if everything is fine.
  */
  bool MatchCenters( const char* det1, const char* det2 );
  
  /*! \fn bool MoveCenter( const char* detSelection, double dX, double dY )
    \brief change the center position of detectors matching detSelection
    \param detselection detector selection
    \param dX the offset in along X (horiz, perp to the beam) (mm)
    \param dY the offset in along Y (horiz, perp to the beam) (mm)
    returns true if everything is fine.
  */
  bool MoveCenter( const char* detSelection, double dX, double dY );

  /*! \fn bool MoveDeadZone( const char* detSelection, 
    double xOffset = 0, 
    double yOffset = 0 );
    \brief move the dead zone center of detectors matching detSelection, taking offsets into account (needed for PB )
    \param detselection detector selection
    \param xOffset offset in X (mm)
    \param yOffser offset in Y (mm)
    returns true if everything is fine.
  */
  bool MoveDeadZone( const char* detSelection, 
    double xOffset = 0, 
    double yOffset = 0 );
  
  /*! \fn bool ChangePitch( const char* detSelection, double pitch, bool changeWirD = true ); 
    \brief change the pitch of detectors matching detSelection, recompute the first wire position accordingly
    \param detselection detector selection
    \param pitch the new pitch (mm)
    \param if false, changeWirD do not recompute the first wire position
    returns true if everything is fine.
  */
  bool ChangePitch( const char* detSelection, double pitch, bool changeWirD = true ); 
  
  /*! \fn bool DumpToFile( const char* outFile = 0 )
    \brief rewrites detector.dat file taking all modifications performed into account.
    \param outFile the name of the file to be written. If the file exists, a backup names <outFile>.N is done
    returns true if everything is fine.
  */
  bool DumpToFile( const char* outFile = 0 );
  
  //! dump all main detectors matching detSelection to screen.
  using TObject::Dump; //jj to avoid compilation warning
  void Dump( const char* detSelection );

  //! dump all main+sub detectors matching detSelection to screen.
  void DumpAll( const char* detSelection );

  /*! \fn TBode* Draw3D( const char* detSelection, bool drawSub = true, bool drawWir = false )
    \brief draw detectors matching detSelection in a 3D display, returns main node
    \param detselection detector selection
    \param drawSub if true,   subdetectors are also drawn
    \param drawWir if true, wires are also drawn
  */
  TNode* Draw3D( const char* detSelection = "*", bool drawSub = true, bool drawWir = false );
 
  //! Reorder Detectors and DeadZones matching selection according to zcm and name
  bool Sort( const char* selection );
  
  #ifndef __CINT__  
  //! Returns list of pointer to all _Main_ detectors (DetectorInfo::GetMain())
  std::vector< DetectorInfo* > GetDetInfo( void );    
  
  //! Returns list of pointer to  _ALL_ detectors (Main+subdetectors)
  std::vector< DetectorInfo* > GetAllDetInfo( void );
  
  //! Returns list of pointer to all magnets
  std::vector< MagnetInfo* > GetMagnetInfo( void ) { return MagnetInfo_; }
  
  //! Returns list of detectors matching \param detSelection
  std::vector< DetectorInfo* > GetDetSelection(  std::string detselection ) { return _GetDetSelection(  detselection ); } 

  //! returns list of projections, computed from detector file
  std::vector< int > GetProjections( void ) { return projections_; }
  #endif
  
  std::string detFileName_;   //!< name of the loaded detector file

  char* detFileContents_;     //!< contents of the loaded detector file
  
  //! To sort pointers to DetectorInfo according to zcm and name
  struct sortDet_: public std::binary_function< DetectorInfo*, DetectorInfo*, bool > { 
    bool operator() ( DetectorInfo* di0, DetectorInfo* di1 ) { 
      if( di0->zcm_ - di1->zcm_ < -1e-5 ) return true; 
      else if( fabs(di0->zcm_ - di1->zcm_)<1e-5 && di0->name_ < di1->name_ ) return true;
      else return false;
    } 
  };
   
  //! To sort pointers to DeadZoneInfo according to zcm and name
  struct sortDeadZ_: public std::binary_function< DeadZoneInfo*, DeadZoneInfo*, bool > { 
    bool operator() ( DeadZoneInfo* dzi0, DeadZoneInfo* dzi1 ) 
    { return ( dzi0->zcm_ < dzi1->zcm_ ); } 
  };
  
  struct sortDet_Line_;   //!< to sort pointers to DetectorInfo according to lineNumber
  struct sortDeadZ_Line_; //!< to sort pointers to DeadZoneInfo according to lineNumber
  
  private:
  
  //! Returns list of detectors matching detSelection
  std::vector< DetectorInfo* >  _GetDetSelection(  std::string detselection );

  //! Returns list of detectors matching deadSelection
  std::vector< DeadZoneInfo* > _GetDeadSelection( std::string deadselection );

  std::vector< MagnetInfo* >   MagnetInfo_; //!< vector of pointers to all MagnetInfo objects
  std::vector< DetectorInfo* > DetInfo_;    //!< vector of pointers to all DetectorInfo objects
  std::vector< DeadZoneInfo* > DeadZInfo_;  //!< vector of pointers to all DeadZoneInfo objects

  std::vector< int > projections_;          //!< vector of angles corresponding to different projections
  
  void _SortVects( void );                  //!< To sort the vector of pointers using a list as a transfer
  ClassDef(DetFileManager,1)
  
};

#endif
