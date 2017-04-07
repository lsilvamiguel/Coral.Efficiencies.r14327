// $Id: HManager.h,v 1.2 2007/06/01 09:13:52 conrad Exp $
 
/*!
   \file    HManager.h
   \brief   Calibration Chain Managment Interface Class.
   \author  Hugo Pereira
   \version $Revision: 1.2 $
   \date    $Date: 2007/06/01 09:13:52 $
*/

#ifndef HManager_h
#define HManager_h

#include <TROOT.h>
#include <TObject.h>
#include <iostream>
#include <vector>
#include <set>
#include <string>

class DetFileManager;
class DetectorInfo;
class TChain;
class TCanvas;
class TH1;
class TH2;
/*! \class HManager 
    \brief Calibration Chain Managment Interface Class.

    This class reads detector table and tracks in root tree format
    obtained from coral. It derives from TObject so that it can be used through root.
*/

class HManager: public TObject {
  public:
  
  /*! \fn  HManager( 
    const char* trackfileselection="", 
    const char* detectorfile="", 
    bool magnets_on=false );
    \brief constructor.
    \param trackfileselection input tree(s) obtained from coral. Wildcards are accepted.
    \param detectorfile the name of the det.dat file to be loaded
  */
  HManager( const char* trackfileselection="", const char* detectorfile=""  );
  
  /*! \fn DetFileManager* LoadDetectorFile( const char* name )
     \brief Load detector.dat file, returns a pointer to the corresponding DetFileManager object
     \param name the name of the file to be loaded.
  */
  DetFileManager* LoadDetectorFile( const char* name );
  
  //! Returns pointer to the DetFileManager object, if any, 0 otherwise
  DetFileManager* GetDetectorFile( void ) { return df_;} 
  
  /*! \fn TChain* LoadTrees(  const char* fileselection )
     \brief load trees contained in coral alignment tree output, returns pointer to the corresponding Tracks object.
     \param fileselection input tree(s). wildcards are accepted. Trees are chained.
  */
  TChain* LoadTrees(  const char* fileselection );
  
  /*! \fn Tracks* AddToChain( const char* fileselection )
     \brief add trees contained in coral alignment tree to existing Tracks object, if any. Returns pointer to tracks object.
     \param fileselection input tree(s). wildcards are accepted. Trees are chained.
     \warning added trees must be of same type
  */   
  TChain* AddToChain( const char* fileselection );
  
  //! Returns pointer to the chain object, if any, 0 otherwise.
  TChain* GetChain( void ) { return chain_; }

  //! To load detectors matching entries in TChain (called in constructor)
  int LoadDetectorsFromChain( void );

  //! To load detectors matching TObjString stored in rootfile
  int LoadDetectorsFromString( void );
  
  /*! \fn void DrawDet( const char* var, const char* detselection, const char* cut="", const char* opt="", int nEnt=0 )
     \brief parse the trees loaded int Tracks object, draw according to the paramaters. Makes a different plot (one plot per pad) for each detector 
     \param var expression of the tree parameters to be drawn
     \param detselection detector selection, using shorten names, character wise wild cards example "GM01 ST03Y2*b DC01X1__"
     \param cut expression of the tree parameters to be used for selection 
     \param opt root options. 
     \param nEnt the number of entries to be parsed. 0 means all
  */
  void Draw(  const char* var, const char* detselection, const char* cut="", const char* opt="", int nEnt=0 );
  bool OneGaus_;              //<! To be passed to all DetInfos

  //! returns pointer to last created HManagerObject
  inline HManager* GetPointer() { return ptr_; }
  static HManager* ptr_;            //!< pointer to last created HManagerObject
  std::vector< DetectorInfo* > dets_;    //!< current vector of selected detector information objects
        
  protected:
  
  //! To check and select detectors matching selection
  int _Select( const char* selection = "*", const char* cut = "1", int nEnt=0 );

  //! To retrieve a 1D histogram associated to a detector from the chain
  TH1* _Get( const char* var, DetectorInfo* d, const char* cut="", int nEnt=0 );

  //! To retrieve a cloned 1D histogram associated to a detector from the chain
  TH1* _GetCloned( const TH1* ref, const char* var, DetectorInfo* d, const char* cut="", int nEnt=0 );

  DetFileManager *df_;            //!< pointer to DetFileManager Object. Contains everything related to the detector.dat file
  std::string detectorFile_;           //!< name of the loaded detector.dat file
  std::vector<std::string> fileSelection_;  //!< vector of all loaded tree filenames
  TChain* chain_;                 //!< coral output chain
  int nEntries_;                  //!< number of entries in chain
  int color_;                     //!< color index for 'same' diagrams 
  TCanvas* cv_;                   //<! Main canvas
 
  std::set< DetectorInfo* > allDets_;   //!< list of detector that have entries in tree
  std::set< std::string > allTBNames_;       //!< list of detector TBNames that have entries in tree
  
  ClassDef(HManager,0)
  
};

#endif
