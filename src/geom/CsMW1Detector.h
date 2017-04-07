// $Id: CsMW1Detector.h,v 1.7 2010/01/24 16:10:41 suhl Exp $

/*!
   \file    CsMW1Detector.h
   \brief   Compass MW1 detector Class.
   \author  Benigno Gobbo
   \version $Revision: 1.7 $
   \date    $Date: 2010/01/24 16:10:41 $
*/

#ifndef CsMW1Detector_h
#define CsMW1Detector_h

#include "coral_config.h"

#include <CLHEP/Matrix/Matrix.h>
#include <list>
#include <vector>
#include <string>
#include <set>
#include "CsDetector.h"
#include "CsHistograms.h"

class CsZone;

/*! \class CsMW1Detector 
    \brief Compass MWPC like detector Class.
*/

class CsMW1Detector : public CsDetector {

 public:


  /*! \fn CsMW1Detector( const int row, 
    const int id, const char* name, const int unit,  
    const int type, const double rdLen, const double xsiz,  
    const double ysiz,  const double zsiz, const double xcm,   
    const double ycm, const double zcm, const HepMatrix rotDRS,
    const HepMatrix rotWRS, const double wirD, const double ang, 
    const int nWir, const double wirP, const double eff, const double bkg,
    const double tGate );
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
  */
  CsMW1Detector( const int    row,
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
              const double tGate );

  bool operator==( const CsMW1Detector& ) const; //!< "equal to" operator

  bool operator<( const CsMW1Detector& ) const; //!< "less than" operator

  //! book all histograms (they should not be created in the constructor)
  void BookHistograms();

  //! Returns true if TCD informations are available
  inline bool hasTDC() const { return( false ); }

  //! Returns true if Drift informations are available
  inline bool hasDrift() const { return( false ); }

  //! Decode the MC data for this detector
  void makeMCDecoding();

  //! Clusterize the digits
  void clusterize();

  /// Decode raw data
  void          DecodeChipDigit               (const CS::Chip::Digit &digit);
  
  //! Read calibration data
  virtual void readCalibration(time_t timePoint);

  //! redefine CsDetector::Wire2Pos 
  virtual double Wire2Pos(double wire) const;


  //LS
  virtual void readMCEffMaps(time_t timePoint);  //! Read Efficiencies for MC.
  //double mcEff_;      //!< efficiency
  std::vector<float> MCEff_, MCEff_err_; //!< MC Efficiency
                                        
 private:

  //LS
  bool mcEffMapsEnabled_; //!< if \c true, eff. Map enabled                        
  bool            decodeCard_;          //!< Decoding set by cards
  bool            associate_;           //!< Associated clusterisation
  float           hittMin_, hittMax_;   //!< Cut on Hit Time (TDC ticks)
  float           cltMin_, cltMax_;     //!< Cut on Cluster Time
  std::vector<float>   calib_data;           //!< Calibration data
  CsHist2D*       hist_;                //!< Time vs channel histogram
  CsHist2D*       histCalib_;           //!< Time vs channel histogram for calib
  CsHist1D*       histCluster_;         //!< Cluster Time
  bool            useCalib_;            //!< Use/ Do Not Use calibrations
  std::string          cdbSwitch_;           //!< Calibration source DB
  std::set<int>        badChannels_;         //!< channels to be removed from clustering

  static const float F1bin_;
};

#endif // CsMW1Detector_h


