// $Id: CsGEMDetector.h 13196 2012-01-12 17:38:47Z suhl $

/*!
   \file    CsGEMDetector.h
   \brief   Compass GEM like detector Class.
   \author  Benigno Gobbo
   \version $Revision: 13196 $
   \date    $Date: 2012-01-12 18:38:47 +0100 (Thu, 12 Jan 2012) $
*/

#ifndef CsGEMDetector_h
#define CsGEMDetector_h

//----------------------------------------------------------------------------

#include "coral_config.h"

#include <CLHEP/Matrix/Matrix.h>
#include <list>
#include "CsDetector.h"
#include "CsHistograms.h"

class CsGEMPlane;
class CsZone;

/*! \class CsGEMDetector
    \brief Compass GEM like detector Class.
*/

class CsGEMDetector : public CsDetector {

  public:

    enum {ChipsPerPlane=6,ChipChannels=128};

    class APVCal
    {
      public:

        class Channel
        {
          public:
            Channel(int f,double pm,double ps,double cm,double cs) :
                      flag(f),
                      pedestal_mean(pm), pedestal_sigma(ps),
                      calibration_mean(cm), calibration_sigma(cs)
                    {}

            Channel(const char *s)
            {
              if( 5!=sscanf(s,"%d%lf%lf%lf%lf",&flag,&pedestal_mean,&pedestal_sigma,
                                               &calibration_mean,&calibration_sigma) )
                throw CS::Exception("CsGEMDetector::APVCal::Channel::Channel(): bad line \"%s\"",s);
            }

            void Print(std::ostream &o=std::cout,const std::string &prefix="") const
            {
              o<<prefix;
              o << "flag=" << flag
		<< " PED=(" 
		<< std::setw(8) << std::setprecision(3) << pedestal_mean
		<< ","
		<< std::setw(8) << std::setprecision(3) << pedestal_sigma
		<< ") CALIB=("
		<< std::setw(8) << std::setprecision(3) << calibration_mean
		<< ","
		<< std::setw(8) << std::setprecision(3) << calibration_sigma
		<< ")" << std::endl;
            }
	
            int    flag;
            double pedestal_mean;
            double pedestal_sigma;
            double calibration_mean;
            double calibration_sigma;
        };

	int src_id;
	int adc_id;
	int chip_id;

        std::vector<Channel> channels;

        friend std::istream & operator >> (std::istream &,std::vector<APVCal> &c);
    };

 public:


  /*! \fn CsGEMDetector( const int row,
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
    \param spSig space resolution, mm  
    \param eGain effective gain a.u.  
    \param eGSig  gain sigma, a.u.
    \param sWidth signal width, mm
    \param tRes  time resolution, ns    
  */
  CsGEMDetector( const int    row,
	      const int    id,    const char* name,  const char *TBname,
	      const int    unit,  const int    type,
	      const double rdLen, const double xsiz,
	      const double ysiz,  const double zsiz,
	      const double xcm,   const double ycm,
	      const double zcm,   const CLHEP::HepMatrix rotDRS,
	      const CLHEP::HepMatrix rotWRS,
	      const double wirD,  const double ang,
	      const int    nWir,  const double wirP,
	      const double eff,   const double bkg,
              const double tGate, const double spSig, 
              const double eGain, const double eGSig, 
              const double sWidth, const double tRes);

  ~CsGEMDetector();

  private:
  //CsGEMDetector( const CsGEMDetector& ); //!< Copy constructor

  //CsGEMDetector& operator=( const CsGEMDetector& ); //!< Assign operator

  public:

  bool operator==( const CsGEMDetector& ) const; //!< "equal to" operator

  bool operator<( const CsGEMDetector& ) const; //!< "less than" operator


  //! book all histograms (they should not be created in the constructor)
  void BookHistograms();

  //! Returns true if TCD informations are available
  inline bool hasTDC() const { return( false ); }

  //! Returns true if Drift informations are available
  inline bool hasDrift() const { return( false ); }

  //! Decode the MC data for this detector
  void makeMCDecoding();

  //! Amplitude correlation is simulated in MC
  inline bool doAmpCorrelationMC() const {return (ampCorrelationMC_ );}

  //! "Master" detector for amp. correlation in MC
  inline bool isMaster() const {return (isMaster_ );}

  //! Returns "associate" detector plane (e.g. for MC amplitude correlation)
  inline CsDetector* getAssociateDet() const { return( associateDet_ ); }

  //! Sets "associate" detector (e.g. for MC amplitude correlation)
  inline void setAssociateDet(CsDetector &det) { associateDet_ = &det; } 

  //! Returns true if amplitude simulation in MC decoding is used
  inline bool amplitudeMCDecoding() const { return( amplitudeMCDecoding_ ); }

  //! Clusterize the digits
  void clusterize();

  /// Decode raw data (called with all digits)
  virtual void DecodeChipDigits(const CS::Chip::Digits &digits);

  /// Decode raw data
  virtual void DecodeChipDigit(const CS::Chip::Digit &digit);

  //! Read calibration data
  virtual void readCalibration(time_t timePoint);

  //! Map channel calibrations
  bool mapChannelCal();

  //! Map channel calibrations for MC
  bool mapMCChannelCal();

  //! Returns time resolution
  inline double getTRes() const { return(tRes_); }

  //! Correction to the amplitude in order to match counterpart (X vs. Y or V vs. U)
  inline const float *getAmpCorr() const { return fAmpCorr; }

 private:

  bool           decodeCard_;      //!<         Decoding set by cards
  bool          decodeLatch_;      //!< raw data is in latch-all format, do
                                   //   the pedestal correction and
                                   //   zero suppression during the decoding
  bool  amplitudeMCDecoding_;      //!< Amplitude simulation in MC 
  bool     ampCorrelationMC_;      //!< Amplitude correlation in MC
  bool             isMaster_;      //!< True if "master" detector
  CsDetector*  associateDet_;      //!< Pointer to the coupled plane 

  double              spSig_;      //!< space resolution, mm
  double              eGain_;      //!< effective gain a.u.
  double              eGSig_;      //!< gain sigma, a.u.
  double             sWidth_;      //!< signal width, mm
  double               tRes_;      //!< time resolution, ns

  std::vector<APVCal>         apv_chan_cals; //!< pedestals and sigmas for each APV chip

  std::map<std::string, CsHist1D*>   mH1;                     //! 1D histogram pointers
  std::map<std::string, CsHist2D*>   mH2;                     //! 2D histogram pointers
  bool sparse;                                   //! true, if zero supression is ON
  int do_clustering;                            //! clustering procedure switch

  float fThresholdHit_;        //! Threshold for hits in units of sigma
  float fThresholdClu_;        //! Threshold for clusters in units of sigma
  float fLowMult_, fUppMult_;  //! Cuts on multiplicity
  std::vector<float> fAmpRatio13_,fAmpRatio23_; //! Polygon for cut on amplitude ratio
  std::vector<float> fRollOverCor13_1, fRollOverCor23_1; //! Polygon for roll over correction
  std::vector<float> fRollOverCor13_2, fRollOverCor23_2; //! Polygon for roll over correction
  float fMCHitTMin, fMCHitTMax;     //! Cuts on hit time (for Monte Carlo)
  float fAmpCorr[3];  //! Correction to the amplitude in order to match counterpart (X vs. Y or V vs. U)

  CsGEMPlane         *plane;        //! Contains hits and clusters
  int fDebugLevel;                  //! Debug level 
  
  bool apv_chan_cals_mapped;        //! flag to indicate mapping of channel calibrations 

};

namespace CsGEM {
  bool IsInside(float xp, float yp, std::vector<float> x, std::vector<float> y);
} // CsGEM

#endif // CsGEMDetector_h
