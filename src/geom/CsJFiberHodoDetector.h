// $Id: 

/*!
   \file    CsJFiberHodoDetector.h
   \brief   Compass Japanese Fiber Hodoscope like detector Class.
   \author  Benigno Gobbo
   \version $Revision: 1.10 $
   \date    $Date: 2011/02/03 14:48:10 $
*/

#ifndef CsJFiberHodoDetector_h
#define CsJFiberHodoDetector_h

#include "coral_config.h"
#include <list>
#include "CsDetector.h"
#include "CsHistograms.h"

class CsZone;

/*! \class CsJFiberHodoDetector 
    \brief Compass Fiber Hodoscope like detector Class.
*/

class CsJFiberHodoDetector : public CsDetector {

 public:


  /*! \fn CsJFiberHodoDetector( const int row, 
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
  CsJFiberHodoDetector( const int    row,
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
              const double tGate, const double tRes );

  private:
  //CsJFiberHodoDetector( const CsJFiberHodoDetector& ); //!< Copy constructor

  //CsJFiberHodoDetector& operator=( const CsJFiberHodoDetector& ); //!< Assign operator

  public:

  bool operator==( const CsJFiberHodoDetector& ) const; //!< "equal to" operator

  bool operator<( const CsJFiberHodoDetector& ) const; //!< "less than" operator

  //! book all histograms (they should not be created in the constructor)
  void BookHistograms();

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
  virtual void  DecodeChipDigit (const CS::Chip::Digit &digit);
  
  /// Read calibration data.
  virtual void readCalibration(time_t timePoint);

 private:
  
  double    tRes_;                  //!< TDC resolution
  bool      decodeCard_;            //!< Decoding set by cards
  float     hitTMin_, hitTMax_;     //!< Cut on Hit Time (ns)
  double    cluDtMax_;              //!< Max. Hit Time Diff in cluster (sigma)


  CsHist1D *mH1[8];                             //!< 1D histogram pointers
  CsHist2D *mH2[2];                             //!< 2D histogram pointers
  bool use_calib;                               //!< use / don't use calibrations
  bool do_clustering;                           //!< clusterization switch

 public:

  class Calib {
  public:
    int chn;
    float data;
    int flag;
    Calib() : chn(0),data(0),flag(0) {}
    Calib(int chn,float data,int flag) :chn(chn),data(data),flag(flag){}  

    friend std::istream& operator>>(std::istream& in,CsJFiberHodoDetector::Calib &c);
    friend std::istream& operator>>(std::istream& in,std::vector<CsJFiberHodoDetector::Calib> &vc);
  };

 private:
  std::vector<Calib> calib;          //!< Calibration values
};

#endif // CsJFiberHodoDetector_h
