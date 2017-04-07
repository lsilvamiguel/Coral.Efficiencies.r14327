// $Id: Align.h,v 1.15 2009/07/21 13:30:19 aaust Exp $
 
/*!
   \file    Align.h
   \brief   Alignment Interface Class.
   \author  Hugo Pereira
   \version $Revision: 1.15 $
   \date    $Date: 2009/07/21 13:30:19 $
*/

#ifndef Align_h
#define Align_h

#include "Defs.h"
#include <TROOT.h>
#include <TObject.h>
#include <iostream>
#include <vector>
#include <string>

class Tracks;
class DetFileManager;
class DetectorInfo;

/*! \class Align 
    \brief Alignment Interface Class.

    This class reads detector table and tracks in root tree format
    obtained from coral. It derives from TObject so that it can be used through root.
    It allows to set all parameters related to alignment as (amongst others)
      - fixing detector parameters
      - tell which parameters has to be fitted
      - add bias to some detectors
      - ...
   
*/
class Align: public TObject {
  public:
  
  /*! \fn  Align( 
    const char* trackfileselection="", 
    const char* detectorfile="", 
    bool magnets_on=false );
    \brief constructor.
    \param trackfileselection input tree(s) obtained from coral. Wildcards are accepted.
    \param detectorfile the name of the det.dat file to be loaded
    \param magnets_on controls the expected structure of the tree. \warning It must match the option used in coral
  */
  Align( 
    const char* trackfileselection="", 
    const char* detectorfile="", 
    bool magnets_on=false );
  
  /*! \fn DetFileManager* LoadDetectorFile( const char* name )
     \brief Load detector.dat file, returns a pointer to the corresponding DetFileManager object
     \param name the name of the file to be loaded.
  */
  DetFileManager* LoadDetectorFile( const char* name );
  
  //! Returns pointer to the DetFileManager object, if any, 0 otherwise
  DetFileManager* GetDetectorFile( void ) { return df_;} 
  
  /*! \fn Tracks* LoadTracks(  const char* fileselection, bool magnets_on = false )
     \brief load tracks contained in coral alignment tree output, returns pointer to the corresponding Tracks object.
     \param fileselection input tree(s). wildcards are accepted. Trees are chained.
     \param magnets_on controls the expected structure of the tree. \warning It must match the option used in coral
  */
  Tracks* LoadTracks(  const char* fileselection, bool magnets_on = false );
  
  /*! \fn Tracks* AddToTracks( const char* fileselection )
     \brief add tracks contained in coral alignment tree to existing Tracks object, if any. Returns pointer to tracks object.
     \param fileselection input tree(s). wildcards are accepted. Trees are chained.
     \warning added trees must be of same type
  */   
  Tracks* AddToTracks( const char* fileselection );
  
  //! Returns pointer to the Tracks object, if any, 0 otherwise.
  Tracks* GetTracks( void ) { return tracks_; }
  
  //! Reload all detectors of DetFileManager object into internal vector
  bool LoadAllDetectors( void );
  
  /*! \fn int UseDetectors( const char* selection )
     \brief selects amongs loaded detectors, detectors to be aligned. returns the number of selected detectors
     \param selection the detectors to select. example: "GM01 ST03Y2*b DC01X1__"
  */
  int UseDetectors( const char* selection );
  
  /*! \fn int ExcludeDetectors( const char* selection )
     \brief exclude detectors from alignment amongst loaded detectors. returns the number of remaining detectors
     \param selection the detectors to select. example: "GM01 ST03Y2*b DC01X1__"
  */
  int ExcludeDetectors( const char* selection );
  
  /*! \fn void ChangeResolution( const char* selection, double res )
     \brief overwrite resolution given in det.dat 
     \param selection loaded detectors whose resolution is to be changed. example: "GM01 ST03Y2*b DC01X1__"
     \param res the new resolution (mm)
  */
  void ChangeResolution( const char* selection, double res );
  
 
  bool AlignU( bool val=true ) { return (alignU_=val); } //!< Controls if alignment of U  (perp to the wires) is to be done
  bool AlignZ( bool val=true ) { return (alignZ_=val); } //!< Controls if alignment of Z  (along the beam) is to be done
  bool AlignT( bool val=true ) { return (alignT_=val); } //!< Controls if alignment of T  (wire orientation) is to be done
  bool AlignP( bool val=true ) { return (alignP_=val); } //!< Controls if alignment of P  (pitch) is to be done
  bool AlignR( bool val=true ) { return (alignR_=val); } //!< Controls if alignment of R0 (prop T0 for drift like detectors) is to be done
  bool AlignL( bool val=true ) { return (alignL_=val); } //!< Controls if alignment of L (Lorentz angle R scaling for drift like detectors) is to be done

  bool AlignOuterST( bool val=false ) { return(alignOuterST_=val); } //!< Controls if alignment of outer straws should be performed (jj)

  bool Iterate( bool val=true ) { return (iterate_=val); } //!< Controls the use of iteration in the matrix inversion
  
  /*! \fn bool SetBias( const char* TBName, const char* biasL )
     \brief add bias to a detector of name \param TBName returns true if Detector was found amongst loaded detectors
     \param biasL bias values in an character string. Order is "<BiasU(mm)> <BiasZ(mm)> <BiasT(deg)> <BiasP(noUnit)>"
  */   
  bool SetBias( const char* TBName, const char* biasL );
   
  /*! \fn bool AddCut( const char* cut )
     \brief Add a selection cut to the Tracks object. Only tree entries matching the cut will be used for alignment
     \param cut the cut to be added. Example "T_chi2/T_ndf<10&&abs(T_duVect)<5". The list of tree variables is found in Tracks 
     object
  */   
  bool AddCut( const char* cut );
  
  //! Sort loaded detectors along z (along the beam)
  bool SortDetectors( void );
  
  /*! \fn void DumpDetectors( const char* selection="*", std::ostream &out = std::cout  )
     \brief dump selected detectors alignment needed parameters to the \param out std::ostream
     \param selection the selection detectors must match to be dumped. example: "GM01 ST03Y2*b DC01X1__"
  */
  void DumpDetectors( const char* selection="*", std::ostream &out = std::cout  );

  /*! \fn bool InitParameters( int nStdDev=1, bool dumpMille=true )
     \brief Init matrix inversion internal parameters. returns true if everything is fine
     \param nStdDev to control the iteration steps.
     \param dumpMille enables millepede output to stdout
  */
  bool InitParameters( int nStdDev=1, bool dumpMille=true );
   
  bool FixU( const char* selection ); //!< Select detectors whose U  (perp to the wire) is not fitted
  bool FixZ( const char* selection ); //!< Select detectors whose Z  (along the beam ) is not fitted
  bool FixT( const char* selection ); //!< Select detectors whose T  (angle) is not fitted
  bool FixP( const char* selection ); //!< Select detectors whose P  (pitch) is not fitted
  bool FixR( const char* selection ); //!< Select detectors whose R0 (prop to T0 for drift like detectors) is not fitted
  bool FixL( const char* selection ); //!< Select detectors whose L  (Lorentz angle R scaling for drift like detectors) is not fitted
  
  /*! \fn  bool Minimize( unsigned int nTracks=0, unsigned int refresh=0 );
     \brief performs the minimisation. Returns true if everything is fine.
     \param nTracks the number of tracks taken from tree to minimize. 0 means all.
     \param refresh to write a line every %refresh events. 0 means nothing written
  */
  bool Minimize( unsigned int nTracks=0, unsigned int refresh=0, const char* reparamdetselection="*" );
  
  /*! \fn bool DumpToFile( const char* filename )
     \brief Dump result of the minimisation to file. Make a backup (<filename>.N) if the file already exist
     returns true if file could be accessed
     \param filename the name of the file
  */ 
  bool DumpToFile( const char* filename, const char* reparamdetselection );
  // bool DumpToFile( string filename );

  bool isBatch_; //!< if true, the input rootfiles are copied to local directory
      
  private:
  //! dump the results of the minimisation to std::ostream \param out
  void _Dump( std::ostream &out = std::cout, const char* reparamdetselection="*"  );
  
  DetFileManager *df_;          //!< pointer to DetFileManager Object. Contains everything related to the detector.dat file
  Tracks *tracks_;              //!< point to Tracks object. Contains chained trees obtained from coral.
  std::string detectorFile_;    //!< name of the loaded detector.dat file
  std::vector<std::string> trackFileSelection_;  //!< vector of all loaded tree filenames
  std::vector< DetectorInfo* > dets_;            //!< vector of selected pointers to DetectorInfo objects used for the alignment minimisation
  std::vector< int > nDetTracks_;                //!< vector of number of tracks passing through the corresponding detector in dets_ (filled only during the minimisation)
  std::string cut_;                              //!< the selection cut applied to the tree entries
  
  bool magnets_on_;  //!< controls the structure of the tree. \warning it must match the corresponding option used by coral
  bool alignU_;      //!< Controls if alignment of U (perp to the wires) is to be done
  bool alignZ_;      //!< Controls if alignment of Z (along the beam) is to be done
  bool alignT_;      //!< Controls if alignment of T (wire orientation) is to be done
  bool alignP_;      //!< Controls if alignment of P (pitch) is to be done
  bool alignR_;      //!< Controls if alignment of R (prop to T0 for drift like detectors ) is to be done
  bool alignL_;      //!< Controls if alignment of L (Lorentz angle R scaling ) is to be done

  bool alignOuterST_; //!< Controls if alignment of outer straws should be performed (jj)

  bool parInit_;     //!< true when InitParameters was called successfuly
  bool iterate_;     //!< controls if millepede makes iterations when inverting the matrix 
  int nStdDev_;      //!< controls the step in matrix inversion iterations
  unsigned int nTracks_; //!< the number of tracks taken from the tree to be used in the minimisation
  float dergb[NGLB];     //!< vector of global derivatives
  float derlc[NPARTRCK]; //!< vector of local derivatives
  float par[NGLB];       //!< vector of parameters
  
  ClassDef(Align,0)
  
};

#endif
