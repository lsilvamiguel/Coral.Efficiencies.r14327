// $Id: Tracks.h,v 1.16 2009/07/21 13:30:19 aaust Exp $

/*!
  \file    Tracks.h
  \brief   chains and loops over coral alignment tree outputs
  \author  Hugo Pereira
  \version $Revision: 1.16 $
  \date    $Date: 2009/07/21 13:30:19 $
*/

#ifndef Tracks_h
#define Tracks_h
#include <string>
#include <vector>

#include <TROOT.h>
#include <TObject.h>
#include "Defs.h"

class TChain;
class TTree;
class DetectorInfo;

/*! \class Tracks 
    \brief chains and loops over coral alignment tree outputs
    The class derives from root TObject so that it can be used through cint
*/
class Tracks: public TObject {
  public: 
  
  
  /*! \fn Tracks( bool magnets_on = 0 );                  
    Constructor. 
    \param magnets_on monitors the tracks tree structure. 
    \warning magnets_on must match the "main magnet on" entry in the options
    file used to create the tree(s) .
  */
  Tracks( bool magnets_on = false );
  
  //! Set iEntry to 0 to restart loop from first event
  void Init( void ) {ientry_ = 0;}                    
  
  /*! \fn unsigned int AddToChain( const char* fileSelection, bool isBatch=false );  
    \brief adds trees found in fileSelection to current chain if any, or create the chain
    \param fileselection the root files in which the trees are looked for
    \param isBatch if true, the files are first copyed to local directory before being loaded
    returns the total number of entries in the chain
  */
  unsigned int AddToChain( const char* fileSelection, bool isBatch=false );  
  
  //! Get the TChain object, if any, 0 ortherwise
  TChain* GetChain( void ) { return T_align_; }          
  
  //! Add a cut to accept/reject entries while looping \param cut the cut to be used
  void AddCut( const char* cut );
  
  //! returns total number of entries in the chain                        
  unsigned int GetEntries( void ) const { return nentries_; } 
  
  //! returns magnets_on_, which controls the type of the chain
  bool MagnetsOn( void ) const { return magnets_on_; }
  
  /*! \fn bool GetNextEntry( void )                       
     \brief get next entry matching selection cut. Returns false if no entry is found.
  */
  bool GetNextEntry( void ); 
  
  /*! \fn bool GetEntry( unsigned int i )                       
     \brief get entry. Returns false if no entry is found. Returns true if entry exists
     \param i the id of the entry to be retrieved.
  */
  bool GetEntry( unsigned int i );
  
  //! returns index of the current downloaded entry
  inline unsigned int GetIEntry( void ) const { return ientry_; }
    
  #ifndef __CINT__
  /*! \fn bool AcceptEntry( vector<DetectorInfo*> dets, unsigned int nPars, double projRange=2 )
    \brief Check if current entry satisfies ndf criterions i.e nMatchingPlanes-nPars > 0. Check the number of projections (at least two needed)
    \param dets the DetectorInfo list which are checked
    \param nPars number of track parameters
    \param projRange max angles (deg) between 2 detectors to tell they are of same projection
  */   
  bool AcceptEntry( std::vector<DetectorInfo*> dets, unsigned int nPars, double projRange=2 );
  
  //! dump current entries to ostream. only the branches matching detectors in dets are dumped
  void DumpEntry( std::ostream &out, std::vector<DetectorInfo*> dets );  
  #endif
  
  //!  Members linked to tree branches at init.
  unsigned int T_evt;     //!< Event number
  unsigned int T_trigMsk; //!< Event trigger mask
  unsigned int T_zone;    //!< Track zone pattern
  unsigned int T_cmlt;    //!< track multiplicity
  double T_chi2;          //!< track chisquare
  int    T_ndf;           //!< track number of degree of freedom
  double T_prob;          //!< track chisquare propability
  double T_meanT;         //!< track mean time
  double T_cop;           //!< charge over momentum
  int T_nDets;            //!< number of fired detectors
  
  int T_detVect[NPLAN];   //!< list of detectors uniqID contributing to track
  double T_uVect[NPLAN];  //!< u_cluster, coordinate perp to the wire
  double T_vVect[NPLAN];  //!< v_cluster, coordinate along wires (2nd proj for pixel detectors)
  double T_duVect[NPLAN]; //!< u_cluster - u_track 
  double T_dvVect[NPLAN]; //!< v_cluster - v_track     (2D detectors)
  double T_rVect[NPLAN];  //!< u_cluster - u_wire      (drift like detectors only)
  double T_tVect[NPLAN];  //!< t_cluster, cluster time (drift like detectors only)
  
  //!  branches specific to magnets_on_ = false
  double T_xLx;           //!< track X
  double T_yLx;           //!< track Y
  double T_zLx;           //!< track Z
  double T_txLx;          //!< track dX/dZ
  double T_tyLx;          //!< track dY/dZ
  
  //!  branches specific to magnets_on_ = true
  double T_xVect[NPLAN];  //!< track X
  double T_yVect[NPLAN];  //!< track Y
  double T_zVect[NPLAN];  //!< track Z
  double T_txVect[NPLAN]; //!< track dX/dZ
  double T_tyVect[NPLAN]; //!< track dY/dZ
  double T_dzVect[NPLAN]; //!< zoffset between detector and recorded point
  double T_BVect[NPLAN];  //!< local magnetic field along the wires
  
  private:
  
  TChain* T_align_;         //!< Chained trees from coral output
  std::string treeName_;    //!< Name of the tree as read from file
  bool magnets_on_;         //!< Switch between magnet on/off tree structure
  unsigned int ientry_;     //!< internal entry counter
  unsigned int nentries_;   //!< number of entries in tree
  std::string cut_;         //!< the cut used for entry selection
  ClassDef(Tracks,0)
};

#endif
