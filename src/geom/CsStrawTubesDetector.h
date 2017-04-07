// $Id: CsStrawTubesDetector.h,v 1.46 2010/05/27 20:01:55 tnagel Exp $

/*!
   \file    CsStrawTubesDetector.h
   \brief   Compass Straw Tubes like detector Class.
   \author  Benigno Gobbo
   \version $Revision: 1.46 $
   \date    $Date: 2010/05/27 20:01:55 $
*/

#ifndef CsStrawTubesDetector_h
#define CsStrawTubesDetector_h

#include "coral_config.h"
#include "CsHistograms.h"
#include "CsRTRelation.h"

#include <CLHEP/Matrix/Matrix.h>
#include <list>
#include <set>
#include "CsDetector.h"

class CsZone;
class CsCluster;

namespace CS
{
class StrawTubes;
class CsStrawTubesDetector_CalibXray;
class CsStrawTubesDetector_CalibT0;
}

/*! \class CsStrawTubesDetector 
    \brief   Compass Drift Chamber like detector Class.
*/

class CsStrawTubesDetector : public CsDetector {

 public:

  ~CsStrawTubesDetector(void);

  /*! \fn CsStrawTubesDetector( const int row, 
    const int id, const char* name, const int unit,  
    const int type, const double rdLen, const double xsiz,  
    const double ysiz,  const double zsiz, const double xcm,   
    const double ycm, const double zcm, const HepMatrix rotDRS,
    const HepMatrix rotWRS, const double wirD, const double ang, 
    const int nWir, const double wirP, const double eff, const double bkg,
    const double tGate,
    const double vel, const double t0, const double thRes, 
    const double spSli, const double tiSli );
    \brief Constructor for tracking detector types.
    \param row   Detector file raw number
    \param id    Detector identification number
    \param name  Detector name (see Comgeant)
    \param unit  Detector number in station
    \param type  Detector type (see Comgeant)
    \param rdLen Radiation lenght
    \param xsiz  Detector X size (DRS) (mm)
    \param ysiz  Detector Y size (DRS) (mm)
    \param zsiz  Detector Z size (DRS) (mm)
    \param xcm   Detector centre X position (MRS) (mm)
    \param ycm   Detector centre Y position (MRS) (mm)
    \param zcm   Detector centre Z position (MRS) (mm)
    \param rotDRS Rotation matrix of DRS w.r.t MRS
    \param rotWRS Rotation matrix of WRS w.r.t MRS
    \param wirD  1st wire offset (WRS) (mm)
    \param ang   Rotation angle of WRS w.r.t MRS
    \param nWir  Number of wires
    \param wirP  Wires pitch
    \param eff   Detector efficiency
    \param bkg   Detector background
    \param tGate Detector gate
    \param vel   Drift velocity (mm/time_slices)
    \param t0    Time Zero (ns)
    \param thRes Two hits resolution (time slices)
    \param spSli Space resolution (time slices)
    \param tiSli Time Slice length (ns)
  */
  CsStrawTubesDetector( const int    row,
	      const int    id,    const char* name,    const char *TBname,
	      const int    unit,  const int    type,
	      const double rdLen, const double xsiz,  
	      const double ysiz,  const double zsiz,
	      const double xcm,   const double ycm,   
	      const double zcm,   const CLHEP::HepMatrix rotDRS,
	      const CLHEP::HepMatrix rotWRS,
	      const double wirD,  const double ang,   
	      const int    nWir,  const double wirP, 
	      const double eff,   const double bkg,
	      const double tGate,
	      const double vel,   const double t0,
	      const double thRes, const double spSli,
	      const double tiSli );

  bool operator==( const CsStrawTubesDetector& ) const; //!< "equal to" operator
  bool operator<( const CsStrawTubesDetector& ) const;  //!< "less than" operator

//__________________________________________________________________
  /*! \fn double getDistToWire( const double time, bool &error )
    \brief Returns distance to wire corresponding to a given drift time
    using rt relation parameters rtPars_ if any, or 'naive' drift velocity vel_
    \param time the drift time
    \param error is set to true if 
       1- the distance is < 0 
       2- the distance is > pitch/2 
       3- the (local) drift velocity is negative
    \warning this parameter is changed!
  */
  double getDistToWire( const double time, bool &error );

  /*! \fn inline bool hasRTRelation( void ) 
    \brief Returns true if a RT relation grid has been set
    returns false otherwise */
	inline bool hasRTRelation( void ) const { return hasRT_; }

  /*! \fn inline CsRTRelation* getRTRelation( void ) const 
    \brief returns RTrelation or NULL depending on hasRT_
  */
  inline CsRTRelation* getRTRelation( void ) const { return (hasRT_) ? rt_ : NULL; }

  //! book all histograms (they should not be created in the constructor)
  void BookHistograms();

//__________________________________________________________________
  inline bool isMWPC() const { return isMWPC_; }    //!< Returns true if MWPC mode is used
  inline bool hasTDC() const { return false; }      //!< Returns true if TCD informations are available
  inline bool hasDrift() const { return !isMWPC_; } //!< Returns true if Drift informations are available
  inline bool multiHitMCDecoding() const { return multiHitMCDecoding_; } //!< Returns true if multi-hit MC decoding is used

  inline double getVel()   const { return(vel_); }    //!< Returns the drift velocity (mm/time_slices)
  inline double getT0()    const { return(t0_); }     //!< Returns Time Zero (ns)   
  inline double getThRes() const { return(thRes_); }  //!< Returns the two hits resolution (time slices)
  inline double getSpSli() const { return(spSli_); }  //!< Returns the space resolution (time slices)
  inline double getTiSli() const { return(tiSli_); }  //!< Returns the time slice (in ns)          
	
  //! return true if associate detector is present
  inline bool hasAssociateDet() const { return( hasAssociateDet_ ); }
	
  //! set associate detector
  inline void setAssociateDet(CsDetector& det) { hasAssociateDet_ = true; associateDet_ = &det; }
	
  //! get associate detector
  inline CsDetector* getAssociateDet() const { return( associateDet_ ); }

  /*! \fn inline void setAssociationCut( const double cut )
    \brief set association cut
    \param cut The cut to be set. is supposed to be positive
  */	
  inline void setAssociationCut( const double cut ) { associationCut_ = cut; }

  /*! \fn inline double getAssociationCut() const
    \brief return association cut, the maximum distance between 
    associated clusters on shifted planes
  */
  inline double getAssociationCut() const { return( associationCut_ ); }
  
  /*! \fn inline void setLRProbCut( const double cut )
    \brief set Left/Right probability cut
    \param cut the cut to be set. is supposed to be in [0,0.5]
  */	
  inline void setLRProbCut( const double cut ) { LRProbCut_ = cut; }
  
  /*! \fn inline double getLRProbCut() const
    \brief return left/right probability cut. Clusters for which LRProb < LRProbCut are assumed to be 
    fake and can be removed. LRProbCut is supposed to be in [0,0.5]
  */
  inline double getLRProbCut() const { return( LRProbCut_ ); }
  
  /*! \fn inline void setLRMode( const int mode )
    \brief set Left/Right ambiguity raising mode
    \param mode the LRMode to be set: 0 - Target Pointing; 1 - cell overlap; other - nothing
  */
  inline void setLRMode( const int mode ) { LRMode_ = mode; }
  
  /*! \fn inline int getLRMode() const
    \brief return Left/Right ambiguity raising mode:
    0 - Target Pointing; 1 - cell overlap; other - nothing
  */
  inline int getLRMode() const { return( LRMode_ ); }
  
  /*! \fn double getCorrU(const CsCluster &c, const double x, const double y, const double t, bool &error)
    \brief Returns U corrected for propagation given MRS (x,y), and
    track latency tt (fast method).
  */
  double getCorrU(const CsCluster *c,
		  const double x, const double y, const double tt,
		  bool &error);

//   /*! \fn double getCorrT(const CsCluster &c, const double x, const double y, const double t, bool &error)
//     \brief Returns T, drift time corrected for propagation given MRS (x,y), and
//     track latency tt (fast method).
//   */
//   double getCorrT(const CsCluster *c,
// 		  const double x, const double y, const double tt,
// 		  bool &error);

  /*! \fn void updateClusters(double time)
    \brief Updates the U of all clusters w/ time offset. Typically this offset corresponds to the time difference: event time (determined from tracking) minus trigger time (which has been used for the initialisation of U). The function leaves the "time_" attribute of the clusters unchanged.
    \param time  Offset in time (will be subtracted from cluster's time).
  */
  void updateClusters(double time);

  void makeMCDecoding();   //! Decode the MC data for this detector
  void clusterize();       //! Clusterize the digits
  virtual void  DecodeChipDigit (const CS::Chip::Digit &digit); //! Decode raw data
  virtual void readCalibration(time_t timePoint);  //! Read calibration data.

  //LS                                                                                              
  virtual void readMCEffMaps(time_t timePoint);  //! Read Efficiencies for MC.                      
  //double mcEff_;      //!< efficiency                                                             
  std::vector<float> MCEff_, MCEff_err_; //!< MC Efficiency  

  private:
  bool isMWPC_;    //!< if \c true, basic clustering is used. no drift time
  bool alwaysTwo_; //!< if \c true, always to clusters are saved even with 0 drift time
  bool multiHitMCDecoding_; //!< if \c true, multi-hit MC decoding
  //LS                                                                                              
  bool mcEffMapsEnabled_; //!< if \c true, eff. Map enabled    
 
  double vel_;     //!< drift velocity (cm/time slice)
  double t0_;      //!< T0 (ns)
  double thRes_;   //!< two hits resol. (time slices)
  double spSli_;   //!< space slice (time slices)
  double tiSli_;   //!< time slice (ns)

  bool _readRTRelation( void );  //!< to create RT Relation from option file
  bool hasRT_;                   //!< tell RT Relation has been set
  CsRTRelation* rt_;			    	 //!< pointer to RT Relation object
  std::set< int > badChannels_;       //!< bad channels removed from the clustering

  bool decodeCard_; //!<	Decoding set by cards

  struct sortZones_;         //!< used for sorting list of pointers to zones
  bool hasAssociateDet_;     //!<	\c true if any associated detector present (DC shifted planes)  
  CsDetector* associateDet_; //!<	pointer to associate detector, if any
  double associationCut_;    //!< maximum distance along u between associated clusters
  double LRProbCut_;         //!< minimum left/right probability to assume that a cluster is a true cluster
  int LRMode_;               //!< 0 - target pointing method; 1 - cell overlap method; other - nothing

#ifdef CsStraw_OLD_PROPAGATION_CORR
   //=== Propagation time (along the wire) correction === 
  int _getFEOrientation(void);  //!< to set feOrientation_ <= Hard coded, for the moment (DEBUG)
  double propVel_inv_;          //!< Inverse Propagation velocity along the wire [ns/mm] ( ~0.005 ns/mm )
   int feOrientation_;          //!< Propagation time correction is: t_cor = t + v * prop_vel_inv_ * feOrientation_
                                //!<    -1 
                                //!<    +1 
                                //!<     0 means no correction
  double dTdX_, dTdY_;          //!< Propagation corrections
#else
  double propVel_;              //!< Propagation velocity along the wire [mm/ns]
#endif
  int splitMCMin_, splitMCMax_; //!< Domain of split straws in MC

  //_____________
  //!< Histograms
  CsHist1D* H_t_;          //!< time distribution    (full range)
  CsHist1D* H_ch_;         //!< channel distribution ( full range )

  //!< (software) Dead Time: To reject retriggering of discri due to arrival of trailing electons clusters
  double deadTime_;

  public:

  //_________________________
  //!< calibration data class
  class CalibrationData {
  public:
    float t0;                 //!< time offset between TDC time and drift Time
    CsRTRelation* rt;         //!< RT Relation
    std::set< int > badChannels;  //!< badChannelsannels removed for the clustering
    CalibrationData(CsDetector &d) :t0( d.getT0() ) { rt = new CsRTRelation( d ); badChannels.clear(); }
    CalibrationData( const char *s ) {}
    ~CalibrationData( void ) { delete rt; }
        
    friend std::istream& operator>>(std::istream& in,CsStrawTubesDetector::CalibrationData &c);
  };

  private:
    CS::CsStrawTubesDetector_CalibXray  *xray;
    CS::CsStrawTubesDetector_CalibT0    *T0s;
};

#endif // CsStrawTubesDetector_h
