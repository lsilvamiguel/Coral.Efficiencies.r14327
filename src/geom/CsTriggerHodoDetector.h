// $Id: CsTriggerHodoDetector.h,v 1.27 2010/08/13 21:45:44 clamar Exp $

/*!
   \file    CsTriggerHodoDetector.h
   \brief   Compass Trigger Hodoscope like detector Class.
   \author  Benigno Gobbo
   \version $Revision: 1.27 $
   \date    $Date: 2010/08/13 21:45:44 $
*/

#ifndef CsTriggerHodoDetector_h
#define CsTriggerHodoDetector_h

#include "coral_config.h"

#include <CLHEP/Matrix/Matrix.h>
#include <iostream>
#include <vector>

#include "CsDetector.h"
#include "CsHistograms.h"

class CsZone;

/*! \class CsTriggerHodoDetector 
    \brief Compass Trigger Hodoscope like detector Class.
*/

class CsTriggerHodoDetector : public CsDetector {

 public:


  /*! \fn CsTriggerHodoDetector( const int row, 
    const int id, const char* name, const int unit,  
    const int type, const double rdLen, const double xsiz,  
    const double ysiz,  const double zsiz, const double xcm,   
    const double ycm, const double zcm, const HepMatrix rotDRS,
    const HepMatrix rotWRS, const double wirD, const double ang, 
    const int nWir, const double wirP, const double eff, const double bkg,
    const double tGate, const double tRes );
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
    \param tRes  Detector TDC resolution
  */
  CsTriggerHodoDetector( const int    row,
                         const int    id,    const char* name,   const char *TBname, 
                         const int    unit,  const int    type,
                         const double rdLen, const double xsiz,  
                         const double ysiz,  const double zsiz,
                         const double xcm,   const double ycm,   
                         const double zcm,   const CLHEP::HepMatrix rotDRS,
                         const CLHEP::HepMatrix rotWRS,
                         const double wirD,  const double ang,   
                         const int    nWir,  const double wirP, 
                         const double eff,   const double bkg,
                         const double tGate, const double tRes );

  ~CsTriggerHodoDetector();

  virtual void AddSubDetector( const int    row,
                               const int    id,    const char* name,
                               const char *TBname,
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

  private:
  //CsTriggerHodoDetector( const CsTriggerHodoDetector& ); //!< Copy constructor

  //CsTriggerHodoDetector& operator=( const CsTriggerHodoDetector& ); //!< Assign operator
  
  public:

  bool operator==( const CsTriggerHodoDetector& ) const; //!< "equal to" operator

  bool operator<( const CsTriggerHodoDetector& ) const; //!< "less than" operator

  //! Returns true if TCD informations are available
  inline bool hasTDC() const { return( true ); }

  //! Returns true if Drift informations are available
  inline bool hasDrift() const { return( false ); }

  //! Returns the TDC resolution (where applicable)
  inline double    getTDCResol() const { return(tRes_); }

  //! Decode the MC data for this detector
  void makeMCDecoding();

  //! Clusterize the digits
  void clusterize();

  /// Decode raw data
  void          DecodeChipDigit                 (const CS::Chip::Digit &digit);

  /// Search for the raw digits of this detectors in all digits from DDD
  void          DecodeChipDigits                (const CS::Chip::Digits &digits);

  /// Read calibration data.
  virtual void readCalibration(time_t timePoint);

  /// Book histograms
  void BookHistograms();

 private:
  
  double tRes_;                    //!<          TDC resolution
  bool           decodeCard_;      //!<         Decoding set by cards
  bool           ishorizontal_;    //!<         Has horizontal elements
  bool           isup_;            //!<         Is the upper part 
  float  hitTMin_, hitTMax_;       //!< Cut on Hit Time (TDC ticks)

  bool          initDecodingDone;
  bool          calcMeanTime;

  std::vector<CsHist1D*>   hists1D_;        //!<         List of all histos
  CsHist1D       *cabs_;           //!<         cluster abscissa (LWRS)
  CsHist1D       *hittimes_, *hittimesz_;

 public:

  class Calib {
  public:
    int chn;
    float data;
    int flag;
    Calib() : chn(0),data(0),flag(0) {}
    Calib(int chn,float data,int flag) :chn(chn),data(data),flag(flag){}  

    friend std::istream& operator>>(std::istream& in,CsTriggerHodoDetector::Calib &c);
    friend std::istream& operator>>(std::istream& in,std::vector<CsTriggerHodoDetector::Calib> &vc);
  };

 private:
  std::vector<Calib> calib;          //!< Calibration values

  static const float F1bin_;
};

#endif // CsTriggerHodoDetector_h
