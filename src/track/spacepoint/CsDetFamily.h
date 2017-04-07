// $Id: CsDetFamily.h,v 1.33 2009/09/14 00:53:51 ybedfer Exp $

/*!
   \file    CsDetFamily.h
   \brief   Detector list to build space points.
   \author  Hugo Pereira
   \version $Revision: 1.33 $
   \date    $Date: 2009/09/14 00:53:51 $
*/

#ifndef CsDetFamily_h
#define CsDetFamily_h

#include "CLHEP/Matrix/Matrix.h"
#include "CsSTD.h"
#include "CsHistograms.h"

class CsDetector;
class CsCluster;
class CsSPMaker;
class CsSpacePoint;
class CsCalSpacePoint;
class CsDriftChamberDetector;
class CsStrawTubesDetector;

class CsDetFamily 
{
 public:
 
  enum spMode { TPOINT, STRAIGHT };  // Modes for fast minimisation
 
  CsDetFamily( int id );
  CsDetFamily( const CsDetFamily& df) { *this = df; }
  CsDetFamily& operator=( const CsDetFamily & );
  virtual ~CsDetFamily() {}
  bool operator==( const CsDetFamily & ) const;

  inline int getID() const { return id_; }                                    //!< detector family index (needed for cuts in option file)
  inline void setName( const std::string &name) { name_ = name; }                  //!< set detector family name              
  inline std::string getName() const { return name_; }                             //!< get detector family name
  virtual std::string getType() const { return "Std"; }                            //!< DUMMY! (used for further devellopment: non drift like detectors)
  
  void addDetector( CsDetector &det, int maxMlt );                            //!< Add one detector to the detFamily, fill corresponding parameters, options, book histograms
  inline std::vector< CsDetector* > getDetectors( void ) const { return d_; }      //!< Detectors in the detFamily
  inline unsigned int detSize( void ) const { return d_.size(); }             //!< Number of detectors in the detFamily
 
  //!< Cuts and parameter settings
  inline void setNClCut( const int cut ) { nClCut_ = cut; }                   //!< Cut on minimum number of clusters/spacepoint
 	inline void setChi2Cut_Fast( const double ccut ) { chi2cut_Fast_ = ccut; }  //!< Cut on fast minimisation chi2
 	inline void setChi2Cut( const double ccut ) { chi2cut_ = ccut; }            //!< Cut on full minimisation chi2
	inline void setZRec( const double z ) { z_ = z; }                           //!< Main position (along the beam) of the detfamily
  
  //!< Tells which mode to use for fast minimisation: STRAIGHT or TPOINT (target pointing)
  inline void setMode( const spMode m ) { 
		mode_ = m; 
    return;
  }
  
  //!< Define Detector family geometry (size and angles)
	inline void setGeometry( 
    const double xMin, const double xMax, 
    const double yMin, const double yMax ) {
    xMin_ = xMin; xMax_ = xMax;
    yMin_ = yMin; yMax_ = yMax;
    geomOK_ = true;
  } 
  
  inline bool   geomOK( void ) const { return geomOK_; }                    //!< Tells if detFamily geometry (max size, max angles) has been set (trough option)
  inline bool   histoBooked( void ) const { return histoBooked_; }          //!< Tells if histograms are booked.
  inline bool   calHistoBooked( void ) const { return calHistoBooked_; }    //!< Tells if calibration histograms (efficiency, resolution) are booked.
  inline bool   LRFilterUsed( void ) const { return useLRFilter_; }         //!< Tells if Left/Right filter is to be used.
  inline bool   removeADet( void ) const { return removeADet_; }            //!< Tells if associate detector is to be removed from cluster combinations. (re-added through cluster LR Association)
  inline bool   useFastMinimisation( void ) const { return useFastMinim_; } //!< Tells if Fast Minimisation (ie. position only, straight tracks or target pointing) is to be done.
  inline bool   useFullMinimisation( void ) const { return useFullMinim_; } //!< Tells if Full Minimisation (ie. position + angles) is to be done.

  inline int    getNClCut( void ) const { return nClCut_; }                 //!< Cut on minimum number of clusters/spacepoint
	inline double getChi2Cut_Fast( void ) const { return chi2cut_Fast_; }     //!< Cut on fast minimisation chi2
	inline double getChi2Cut( void ) const { return chi2cut_; }               //!< Cut on full minimisation chi2
	inline double getZRec( void ) const { return z_; }                        //!< Main position (along the beam) of the detfamily
 	inline spMode getMode( void ) const { return mode_; }                     //!< Mode is straight (1) or target pointing (0) for fast minimisation
  inline std::vector< int > getMaxMlt( void ) const { return maxMlt_; }          //!< Max cluster multiplicity/cluster 
  inline std::vector< CsHist1S* > getResHists( void ) const { return H_du_;}        //!< 1D Residual histograms - when calHistoBooked == true 
  inline std::vector< CsHist2S* > getRes2DHists( void ) const { return H_du_vs_r_;} //!< 2D Residual histograms (du_vs_rAlgebric) - when calHistoBooked == true 
  inline std::vector< int > getADetIndexes( void ) const { return d_A_; }           //!< list of associate detector index( that is d_[d_A_[i]] is associated whith d_[i] )
  
  std::vector< CsSpacePoint* > getSpacePoints( void );              //!< build and return space points, fill monitoring histograms, if booked
  std::vector< CsCalSpacePoint* > getCalibrationSpacePoints( void );   //!< perform internal calibration, return corresponding spacepoints
  void cleanEvent( void );                        //!< clean event; 

  //!< Debuging tools
  void dumpConfig( void );              //!< Write all parameters (detectors, cuts, switches, troubles) to screen
  void dumpGeometry( void );            //!< Write geometry (size, max angles) to screen
  void dumpAssociationTable( void );    //!< Dump detector associations if any.
  void dumpClusters( void );            //!< Write all recorded event clusters infos, needed for spacepoints 
  void dumpClusterSizes( void );        //!< Write number of registered cluster / detector.
  void dumpCmbIndexes( void );          //!< Dump the cluster index combination beeing processed
  void dumpNCmb( void );                //!< Dump the total number of combinations in the event

  //!< Histogram level
  enum histogramLevel { None, Normal, High };
  inline histogramLevel GetHistLevel() const { return hLevel_; }

	
 private:

  //!< parameters
	int id_;							 //!< family id
	std::string name_;					 //!< family name
  double z_;  					 //!< absissa (mm), along the beam, of the plane in which space point coordinates are calculated
	int    nClCut_;				 //!< minimal number of clusters needed in space point
	double chi2cut_Fast_;  //!< cut on fast minimisation chi2 ( position only)
	double chi2cut_; 			 //!< cut on full minimisation chi2 ( position + angles)
	spMode mode_;          //!< mode for minimise_Fast. TPOINT is Target Pointing; Straight is straight tracks;
  bool   useLRFilter_;   //!< should LR filter be used to clean imported cluster lists
  bool   removeADet_;    //!< associate det is removed from minimisation for calibration (RT, resolution, efficiency)
  bool   useFastMinim_;  //!< fast minimisation: position only 
  bool   useFullMinim_;  //!< full minimisation: position+angles
  bool   useClusterAsc_; //!< cluster association means combinations are performed only on half the detectors. 
                         //!< Only clusters with associates are kept.
  unsigned int bigMask_; //!< used to remove ALL associated dets when useClusterAsc_ is true
  
  bool geomOK_;
  double xMin_, xMax_;
  double yMin_, yMax_;
  double tgCenter_;      //!< Target center, taken from CsGeom
  
 	//!< detectors and related infos
	std::vector< CsDetector* >  d_;		            //!< detectors
	std::vector< CsDriftChamberDetector* >  dc_;		//!< drift chamber detectors matching d_
	std::vector< CsStrawTubesDetector* >    st_;		//!< straw tubes detectors matching d_
  std::vector< int > d_A_;           //!< corresponding index of associated detector, if any, -1 otherwise.
  std::vector< CLHEP::HepMatrix > iR_;      //!< inverse rotation matrices associated with d_
  std::vector< double > wirP_;       //!< detector pitches
  std::vector< double > wirD_cor_;   //!< detector first wire position corrected by center offset
  std::vector< int > maxMlt_;        //!< max allowed multiplicity for each detector (read from option file)
  
  //!< imported clusters and related infos
  bool _importClusters( void );         //!< Fills all needed vectors, event/event.
  bool _checkMlt( unsigned int mask );  //!< Check the cluster multiplicity/detector. mask tels which detectors are not checked, as not used for the combinations

  std::vector< CsCluster* > c_;       //!< vector of imported clusters
  std::vector< unsigned int > cd_;    //!< clusters detector index in d_
  std::vector< double > uc_;          //!< clusters position perp to the wires
  std::vector< double > wc_;          //!< clusters position along the beam
  std::vector< double > rc2_;         //!< squared cluster error perp to the wires 

  //!< matrices for fast minimisation ( calculated at _importCluster to spare time )
  void _fillMatrix_Fast( void );
  CLHEP::HepMatrix _getA_Fast( CLHEP::HepMatrix iR, double uc, double wc, double rc2 );
  CLHEP::HepMatrix _getB_Fast( CLHEP::HepMatrix iR, double uc, double wc, double rc2 );
  std::vector< CLHEP::HepMatrix > A_Fast_;   
  std::vector< CLHEP::HepMatrix > B_Fast_;   
 
  //!< matrices for full minimisation
  void _fillMatrix( void );
  CLHEP::HepMatrix _getA( CLHEP::HepMatrix iR, double uc, double wc, double rc2 );
  CLHEP::HepMatrix _getB( CLHEP::HepMatrix iR, double uc, double wc, double rc2 );
  std::vector< CLHEP::HepMatrix > A_;   
  std::vector< CLHEP::HepMatrix > B_;   

  //!< following members are filled only if useClusterAsc_ is true.
  std::vector< int > ic_A_;           //!< associated cluster index, if any. -1 otherwise
  std::vector< CsCluster* > c_A_;      //!< associated cluster, if any, NULL otherwise
  
  //!< for following vectors, XXXX_A_[i] makes sense only when ic_A_[i] != -1
  std::vector< unsigned int > cd_A_;    //!< associated clusters detector index in d_.      
  std::vector< double > uc_A_;          //!< associated cluster position perp to the wires. 
  std::vector< double > wc_A_;          //!< associated clusters position along the beam.   
  std::vector< double > rc2_A_;         //!< squared associated cluster error perp to the wires.           
  std::vector< CLHEP::HepMatrix > A_Fast_A_;   //!< Fast and full minimisation matrices for associated clusters.  
  std::vector< CLHEP::HepMatrix > B_Fast_A_;   
  std::vector< CLHEP::HepMatrix > A_A_;   
  std::vector< CLHEP::HepMatrix > B_A_;   
  
  //!< clusters combinations
  void         _getNextClusterCmb( unsigned int mask ); 
  unsigned int _getCmbSize( unsigned int mask );
  bool hasClusterCmb_;
  int nCmb_;
  std::vector< unsigned int > nc_;   //!< number of detector clusters
  std::vector< unsigned int > fc_;   //!< index of first detector cluster in c_;
  std::vector< unsigned int > ic_;   //!< index of current detector cluster in c_;
  
  std::vector< CsSpacePoint* >    sp_;      //!< list of event spacepoints
  std::vector< CsCalSpacePoint* > sp_cal_;  //!< list of event calibration spacepoints
  
  //!< spacePoint fast minimisation
  bool _minimise_Fast(
    double &x, double &y, 
    double &chi2,
    unsigned int mask );   

  //!< spacePoint full minimisation
  bool _minimise( 
    double &x, double &y, 
    double &tx, double &ty, 
    double &chi2,
    unsigned int mask );       
  
  
  
  //!< monitoring histograms
  histogramLevel hLevel_;     //!< None: no histos; 
                              //!< Normal: monitor histos;
                              //!< High: monitor+internal calibration histos

  bool histoBooked_;
  void _bookHistograms( void ); //!< book monitor Histograms
  CsHist1S* H_spMlt_;           //!< space point multiplicity/event
  CsHist1S* H_spSize_;          //!< space point number of clusters 
  CsHist1S* H_status_;          //!< status for booked space points
    
  CsHist1S* H_chi2_Fast_;       //!< fast minimisation chi2
  CsHist1S* H_chi2_;            //!< full minimisation chi2
  CsHist1S* H_chi2_min_;        //!< full minimisation min chi2

	CsHist2S* H_profile_Fast_;    //!< fast minimisation 2D profile 
	CsHist2S* H_profile_;         //!< full minimisation 2D profile 
	CsHist2S* H_angles_;          //!< fast minimisation 2D angle profiles 
  CsHist2S* H_tpx_;             //!< target pointing correlation in X
  CsHist2S* H_tpy_;             //!< target pointing correlation in y

  //!< calibration histograms
  bool calHistoBooked_;
  void _bookCalHistograms();  //!< Calibration histogram booking. Done at first call to spCalibration()
  
  unsigned int nTBins_;
  
  std::vector< CsHist2S* > H_r_vs_t_;   //!< rt relation histogram vector (r algebric)
  std::vector< CsHist2S* > H_ra_vs_t_;  //!< rt relation histogram vector (r > 0)
  CsHist2S* H_r_vs_t_uniq_;
  CsHist2S* H_ra_vs_t_uniq_;
  
  std::vector< CsHist1S* > H_du_;        //!< residual histogram vector
  std::vector< CsHist2S* > H_du_vs_u_;   //!< correlation histogram vector over the whole chamber
  std::vector< CsHist2S* > H_du_vs_v_;   //!< correlation histogram vector vs distance Along the wire
  std::vector< CsHist2S* > H_du_vs_r_;   //!< correlation histogram vector over one cell

  std::vector< CsHist2S* > H_dt_vs_r_;   //!< correlation histogram vector over one cell
  std::vector< CsHist1S* > H_spMlt_cal_; //!< space point multiplicity/calibration event for each detector
  std::vector< CsHist1S* > H_chi2_cal_;  //!< ref point chi2 
  std::vector< CsHist1D* > H_wir_;       //!< reference profile for efficiency
  std::vector< CsHist1D* > H_wid_;       //!< corresponding found clusters 

};

#endif

 
  
