// $Id: CsMWPCDetector.h,v 1.30 2010/01/24 16:10:41 suhl Exp $

/*!
   \file    CsMWPCDetector.h
   \brief   Compass MWPC like detector Class.
   \author  Benigno Gobbo
   \version $Revision: 1.30 $
   \date    $Date: 2010/01/24 16:10:41 $
*/

#ifndef CsMWPCDetector_h
#define CsMWPCDetector_h

#include "coral_config.h"

#include <list>
#include <vector>
#include "CsDetector.h"
#include "CsHistograms.h"

class CsZone;

/*! \class CsMWPCDetector 
    \brief Compass MWPC like detector Class.
*/

class CsMWPCDetector : public CsDetector {

 public:


  /*! \fn CsMWPCDetector( const int row, 
    const int id, const char* name, const int unit,  
    const int type, const double rdLen, const double xsiz,  
    const double ysiz,  const double zsiz, const double xcm,   
    const double ycm, const double zcm, const HepMatrix rotDRS,
    const HepMatrix rotWRS, const double wirD, const double ang, 
    const int nWir, const double wirP, const double eff, const double bkg,
    const double tGate, const double vel, const double spSli, 
    const double tiSli);
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
    \param t0    Time Zero (time slices) 
    \param thRes Two hits resolution (time slices)
    \param spSli Space resolution (time slices)
    \param tiSli Time Slice length (ns)
  */
  CsMWPCDetector( const int    row,
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
              const double tGate, const double vel,
              const double t0,    const double thRes,
              const double spSli, const double tiSli );

  private:
  //CsMWPCDetector( const CsMWPCDetector& ); //!< Copy constructor

  //CsMWPCDetector& operator=( const CsMWPCDetector& ); //!< Assign operator

  public:

  bool operator==( const CsMWPCDetector& ) const; //!< "equal to" operator

  bool operator<( const CsMWPCDetector& ) const; //!< "less than" operator

  //! book all histograms (they should not be created in the constructor)
  void BookHistograms();

  //! Returns true if TCD informations are available
  inline bool hasTDC() const { return( false ); }

  //! Returns true if Drift informations are available
  inline bool hasDrift() const { return( false ); }

  //! Returns true if drift time is simulated in makeMCDecoding
  inline bool hasDriftTimeMC() const { return( driftTimeMC_ ); }

  //! Decode the MC data for this detector
  void makeMCDecoding();

  //! Clusterize the digits
  void clusterize();

  /// Decode raw data
  void          DecodeChipDigit               (const CS::Chip::Digit &digit);
  
  //! Read calibration data
  virtual void readCalibration(time_t timePoint);

  //LS
  virtual void readMCEffMaps(time_t timePoint);  //! Read Efficiencies for MC.
  //double mcEff_;      //!< efficiency
  std::vector<float> MCEff_, MCEff_err_; //!< MC Efficiency

 private:
 
 
  bool            decodeCard_;          //!< Decoding set by cards
  bool            associate_;           //!< Associated clusterisation
  bool            driftTimeMC_;         //!< drift time simulation 

  //LS
  bool mcEffMapsEnabled_; //!< if \c true, eff. Map enabled

  double vel_;                          //!< drift velocity  (mm/time slice)  
  double t0_;                           //!< T0              (time slice)    
  double thRes_;                        //!< two hits resol. (time slices)
  double spSli_;                        //!< space slice (time slices)        
  double tiSli_;                        //!< time slice (ns) 

  float           hittMin_, hittMax_;   //!< Cut on Hit Time (TDC ticks)
  float           cltMin_, cltMax_;     //!< Cut on Cluster Time
  std::vector<float>   calib_data;           //!< Calibration data
  CsHist2D*       hist_;                //!< Time vs channel histogram
  CsHist2D*       histCalib_;           //!< Time vs channel histogram for calib
  CsHist1D*       histCluster_;         //!< Cluster Time

  std::map<std::string, CsHist1D*>   mH1;         //! 1D histogram pointers

  bool            useCalib_;            //!< Use/ Do Not Use calibrations
  std::string          cdbSwitch_;           //!< Calibration source DB
  std::set<int>        badChannels_;         //!< channels to be removed from clustering

  static const float F1bin_;
};

#endif // CsMWPCDetector_h


