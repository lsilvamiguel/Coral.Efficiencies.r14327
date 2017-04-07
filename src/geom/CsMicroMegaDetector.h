// $Id: CsMicroMegaDetector.h,v 1.48 2010/01/24 16:10:41 suhl Exp $

/*!
   \file    CsMicroMegaDetector.h
   \brief   Compass Micro-Mega like detector Class.
   \author  Benigno Gobbo
   \version $Revision: 1.48 $
   \date    $Date: 2010/01/24 16:10:41 $
*/

// $Log: CsMicroMegaDetector.h,v $
// Revision 1.48  2010/01/24 16:10:41  suhl
// * replacing "unsigned int" or "int" used to store timestamps with "time_t"
//   to fix the strict-aliasing compiler warning
//
// Revision 1.47  2008/01/08 18:09:16  ybedfer
//  ...
//
// Revision 1.46  2007/03/12 10:24:00  ybedfer
//  Time resolution: input via option (can be made to vary from one plane to
// the next), accessor. (Note: In the MC simulation: a built-in of 9 ns is
// still used.)
//
// Revision 1.45  2004/11/20 21:12:59  ybedfer
//  Make use of Colin's "MMLibrary".
//
// Revision 1.44  2003/04/17 09:20:21  benigno
// --> gcc3
//
// Revision 1.43  2003/04/07 13:53:50  cbernet
// *** empty log message ***
//
// Revision 1.42  2003/03/13 15:41:43  cbernet
// added 	hit crates histograms in MM class.
// 	hit ctoff
//
// hit crates is a profile histogram filled only with hits being inside an off time window
// specified with the following line in the option file
//
// MM01V* HitTimeOff [0-1] 300 800
//
// Revision 1.41  2003/02/16 11:45:26  cbernet
// new option :
//
// MM02V* FlagChannels	false/true
//
// when false, channels flagged as bad will be kept. If true or missing, these
// channels are removed at decoding stage
//
// Revision 1.40  2003/02/12 07:59:28  cbernet
// triangle cut in the tot VS t plane can now be defined for MM.
// add this line to the option file :
//
// MM* TriCut  0       0       4       4       400
//
// this cut defines 2 straight lines with slopes 4 and 4.
// 1st 2 numbers are the (t,tot) coordinates of the intersection point
// 400 is tot_max.
//
// this cut is not yet in use.
//
// Revision 1.39  2003/01/23 13:19:46  cbernet
// possibility do define an off trigger time range in the option file :
// MM02* ClusterTime [0-1] -10000 10000
//
// if this cut is passed, new histograms are filled :
// [TBName]_crates		cluster rates
// [TBName]_ctoff         	off trigger cluster time distribution
//
// #ifdef uM_DETECTOR_STUDIES
// is now only used for some cuts. histograms don't depend on defines anymore,
// but on the histogram level
//
// Revision 1.38  2002/12/19 14:27:26  neyret
// Huge clean-ups in calibration code for preparation of the future calib
// MySQLDB
// + some bugs fixed
//
// Revision 1.37  2002/10/18 16:05:57  cbernet
// - created a new 2d histo hit time vs chan for calibration purpose
// (histogramminfg level = high)
//
// - range of histos related to hit time depends on calibration usage.
//
// Revision 1.36  2002/07/01 17:16:46  ybedfer
//  "COMPASS_SETUP == 2001' conditions "uM_PARTIALLY_EQUIPPED".
//
// Revision 1.35  2002/04/22 16:27:55  ybedfer
//  Modifications in the handling of time:
//   - Hit time cut:
//      - from a HitTime range no longer centered on 0,
//      - modifiable by option.
//   - Hit time cut in MC:
//      - according to HitTime range,
//      - on TDC time (instead of MC time).
//   - MC TDC time:
//      - retain TDC closest to trigger (instead of earliest).
//
// Revision 1.34  2002/03/15 11:02:50  ttoeda
// Modify readCalibration and DecodeChipDigit
//
// Revision 1.33  2002/03/08 16:08:24  hpereira
// /src/geom/: almost all detector modified to unify histogram switches. Option format is to be found in /pkopt/hist.opt
// CsMicroMegaDetector.h and .cc: more modifications to allow call to bookHistograms after call to calibrationDB. DeadChannels and partiallyEquiped flags are now set in ReadCalib.
// CsDetector.h/.cc: added method ReadHistLevel(), GetHistLevel() and member hLevel.
// GetHistLevel() Returns either CsDetector::None (ie 0) | CsDetector::Normal (ie 1) | CsDetector::High (ie 2) of type CsDetector::histogramLevel (enum) or unsigned int through cast
//
// Revision 1.32  2001/12/15 16:20:18  ttoeda
// Modified part of reading calibration
//
// Revision 1.31  2001/12/05 21:30:50  ybedfer
//  - Cuts on time (hit, cluster) made non static. Can be modified by
//   options.
//  - Introduced new data members for the cuts necessited by
//   "uM_DETECTOR_STUDIES".
//
// Revision 1.30  2001/11/22 16:32:37  bernet
//
// New function CsDetector::Wire2Pos(int wire) implemented for clustering
// in variable pitch detectors.
//
// this function is now used also in CsTriggerHodoDetector::clusterize()
//
// new histogram added to CsTriggerHodoDetector : cluster position
//
// Revision 1.29  2001/10/26 12:14:25  ybedfer
// *** empty log message ***
//

#ifndef CsMicroMegaDetector_h
#define CsMicroMegaDetector_h

#include "coral_config.h"
#include <CLHEP/Matrix/Matrix.h>
#include <list>
#include <set>
#include "CsDetector.h"
#include "MMLibrary.h"


//#include "CsDecMap.h"
#include "CsHistograms.h"

/*! Leading Edge parity for ChipF1Data.
  1 if leading corresponds to odd data;
  0 if leading corresponds to even data
*/

class CsZone;

/*! \class CsMicroMegaDetector
    \brief Compass Micro-Mega like detector Class.
*/

class CsMicroMegaDetector : public CsDetector {

 public:


  /*! \fn CsMicroMegaDetector( const int row,
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
  CsMicroMegaDetector( const int    row,
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

  ~CsMicroMegaDetector();

  bool operator==( const CsMicroMegaDetector& ) const; //!< "equal to" operator

  bool operator<( const CsMicroMegaDetector& ) const; //!< "less than" operator

  //! Returns true if TCD informations are available
  inline bool hasTDC() const { return( false ); }

  //! Returns the time resolution (improperly called TDC resolution)
  inline double getTDCResol() const { return tRes_; }

  //! Returns true if Drift informations are available
  inline bool hasDrift() const { return( false ); }

  //! Decode the MC data for this detector
  void makeMCDecoding();

  //! Clusterize the digits
  void clusterize();

  /// Decode raw data
  void  DecodeChipDigit (const CS::Chip::Digit &digit);

  /// Read calibration data.
  void  readCalibration(time_t timePoint);

  /// Book histograms
  void BookHistograms();

  //LS
  virtual void readMCEffMaps(time_t timePoint);  //! Read Efficiencies for MC.
  //double mcEff_;      //!< efficiency
  std::vector<float> MCEff_, MCEff_err_; //!< MC Efficiency                                       
 private:

  // will be const next year !
  unsigned int parityLead_;

  /// returns true if data is a leading time
  inline bool IsLeading( const unsigned data ) { return ((data & 1) == (unsigned int) parityLead_) ? true : false; }

  /// returns true if data is a trailing time
  inline bool IsTrailing( const unsigned data ) { return ((data & 1) == (unsigned int) parityLead_) ? false : true; }

  //LS
  bool mcEffMapsEnabled_; //!< if \c true, eff. Map enabled
 
  bool     decodeCard_;       //!< Decoding set by cards
  bool     useCalib_;         //!< Calibrations ON/OFF flag
  bool     associate_;        //!< Associated clusterisation
  int     splitClustersMin_; //!< Cluster splitting min chan
  int     splitClustersMax_; //!< Cluster splitting max chan
  float    clTMin_, clTMax_;  //!< Cut on Cluster Time (ns)
  float    cltoffMin_, cltoffMax_;  //!< Cuts for off-time clusters (ns)
  float    hitTCut_;          //!< Cut on Hit Time (TDC ticks) (derived from previous)
  float    hitTMin_, hitTMax_;//!< Cut on Hit Time (TDC ticks)
  float    hitTOffMin_, hitTOffMax_;//!< Cut on Hit Time for offtime clusters (TDC ticks)
  float    cresolution_;      //!< Cluster resolution ( by default pitch/sqrt(12) )
  double   tRes_;             //!< Time resolution (default = 9 ns, modified by option)

  // Members brought into play upon macro defined
  float    ToTMin_;           //!< Cut on ToT           (ifdef uM_DETECTOR_STUDIES)
  int      clsMax_;           //!< Cut on cluster size  (ifdef uM_DETECTOR_STUDIES)
  int      nclMax_;           //!< Cut on cluster multi (ifdef uM_SINGLE_TRACK)
  int      channel0_,nWtot_;  //!< Channel # and offset (ifdef uM_PARTIALLY_EQUIPPED)


  std::vector<CsHist1D*>  hists1D_;  //!< 1D histograms
  std::vector<CsHist2D*>  hists2D_;  //!< 2D histograms

  /// 1D histograms
  CsHist1D *hch_, *ht_, *hrates_, *htoff_, *cch_, *cs_, *ct_ , *ctot_, *lt_, *tt_, *crates_, *ctoff_;
  
  /// 2D histograms
  CsHist2D *htvsch_, *ctotvst_;

  int      lastChan_;        //!< channel in last digit
  double   *digitData_;      //!< data for digit creation data[0]=chan, data[1]=time
  std::map<int, std::list<CsDigit*>::iterator >  oDigits_;    //!<digits selected for clustering (ordered)

  
  double W2P(double wire);
  double w2pLow, w2pUp, w2pP1, w2pO1, w2pP2, w2pO2, w2pP3, w2pO3; // W2p parameterisation

  std::vector< MM::ChannelCalibExt >   calib_data_ext;
  bool                                 removeBadChannels_;
  
  MM::ECutTriangle *cutTriangle_;

  static const float leadtWght_;
  static const float f1Tick_;
};


//------------------------------------------------------------------------------


#endif // CsMicroMegaDetector_h



