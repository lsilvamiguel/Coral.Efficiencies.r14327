// $Id: CsPixelGEMDetector.h 13196 2012-01-12 17:38:47Z suhl $

/*!
   \file    CsPixelGEMDetector.h
   \brief   Compass Pixel GEM class definition.
   \author  Yann.Bedfer@cern.ch
   \version $Revision: 13196 $
   \date    $Date: 2012-01-12 18:38:47 +0100 (Thu, 12 Jan 2012) $
*/

// - CsPixelGEMs (or CsPGs, for short) are the logical entities describing
//  physical pixelGEMs. 3 logical objects per physical one: 1 for the pixelised
//  central core, and 2 for the orthogonal stripped pieces on the outskirts.
// - The CsPG for the pixel core is constructed in 2 steps:
//     I) Constructor proper.
//    II) "CsPG::AddSubDetector".
// - This choice is made by considering how the object could be read from the
//  "detectors.dat", w/o changing much in either detectors.dat's structure or
//  the routine reading it, viz. "CsGeom::readDetTable". It re-cycles the idea
//  already used for constructing CsMicroMegas and other variable pitch
//  detectors:
//   - All info is stored in standard "det" entries (of "detectors.dat").
//   - Extra (w.r.t. simpler, single-pitch, single dim, detectors) info is
//    stored in extra "det" entries, following the 1st one and bearing the
//    same TB name.
//   - Info from these extra lines is processed by a "CsPG::AddSubDetector"
//    paralleling the "CsDetector::AddSubDetector" used for variable picth.
//   - The piece of "CsGeom::readDetTable" software already at play for
//    variable pitch is used to re-direct the info to "CsPG::AddSubDetector".

#ifndef CsPixelGEMDetector_h
#define CsPixelGEMDetector_h

//----------------------------------------------------------------------------

#include "coral_config.h"

#include <CLHEP/Matrix/Matrix.h>
#include <list>
#include "CsDetector.h"
#include "CsHistograms.h"
#include "DaqDataDecoding/Chip.h"
#include "DaqDataDecoding/ChipAPV.h"

class CsGEMPlane;
class CsPixelGEMPlane;
class CsZone;

/*! \class CsPixelGEMDetector
    \brief Compass PixelGEM like detector Class.
*/

class CsPixelGEMDetector : public CsDetector {
  
 public:
  
  int ChipsPerPlane;
  enum {ChipsPerStripPlane=4, ChipsPerPixelPlane=8, ChipChannels=128};
  
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
                throw CS::Exception("CsPixelGEMDetector::Calib::Channel::Channel(): bad line \"%s\"",s);
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

  /*! \fn CsPixelGEMDetector( const int row,
    const int id, const char* name, const int unit,
    const int type, const double rdLen, const double xsiz,
    const double ysiz,  const double zsiz, const double xcm,
    const double ycm, const double zcm, const HepMatrix rotDRS,
    const HepMatrix rotWRS, const double wirD, const double ang,
    const int nWir, const double wirP, const double eff, const double bkg,
    const double tGate );
    \brief Constructor for CsPG: first part, cf. "AddDetector" infra for 2nd part (required for CsPG of the pixel kind).
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
    \param spSig Space resolution, mm  
    \param eGain Effective gain a.u.  
    \param eGSig  Gain sigma, a.u.
    \param sWidth Signal width, mm
    \param tRes  Time resolution, ns    
  */
  CsPixelGEMDetector(const int    row,
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
		     const double sWidth,const double tRes);

  ~CsPixelGEMDetector();

  /*! \fn AddSubDetector( const int row,
    const int id, const char* name, const int unit,
    const int type, const double rdLen, const double xsiz,
    const double ysiz,  const double zsiz, const double xcm,
    const double ycm, const double zcm, const HepMatrix rotDRS,
    const HepMatrix rotWRS, const double wirD, const double ang,
    const int nWir, const double wirP, const double eff, const double bkg,
    const double tGate );
    \brief Complete the construction of PixelGEM 2D-detector.
  */
  void AddSubDetector(const int    row,
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
		      const double tGate);

  //! Returns the first pixel position along V (=WRS+pi/2) (mm)
  double getWirDV() const { return(wirDV_); }
  //! Returns the pixel pitch along V (mm)
  double getWirPV() const { return(wirPV_); } 
  //! Returns number of pixels along V
  int    getNWirV() const { return(nWirV_); }  

 private:

  bool operator==( const CsPixelGEMDetector& ) const; //!< "equal to" operator

  bool operator<( const CsPixelGEMDetector& ) const; //!< "less than" operator

  public:

  //! book all histograms (they should not be created in the constructor)
  void BookHistograms();

  //! Returns true if TCD informations are available
  inline bool hasTDC() const { return( false ); }

  //! Returns true if Drift informations are available
  inline bool hasDrift() const { return( false ); }

  //! Decode the MC data for this detector
  void makeMCDecoding();

  private: // decode the MC data for the two kind
           // of detectors (pixel/strip) and the 
           // two different options in hit
           // creation (simplistic/with amplitude)
  void makeMCPixelsSimple();
  void makeMCPixelsAmplitude();
  void makeMCStripsSimple();
  void makeMCStripsAmplitude();

  public:

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
  void clusterizePixels();
  void clusterizeStrips();

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
  int do_clustering;                            //! clustering procedure switch
  int fClusterConfigMask;                       //! clustering configuration

  float fThresholdHit_;        //! Threshold for hits in units of sigma
  float fThresholdClu_;        //! Threshold for clusters in units of sigma
  float fLowMult_, fUppMult_;  //! Cuts on multiplicity
  std::vector<float> fAmpRatio13_,fAmpRatio23_; //! Polygon for cut on amplitude ratio
  std::vector<float> fRollOverCor13_1, fRollOverCor23_1; //! Polygon for roll over correction
  std::vector<float> fRollOverCor13_2, fRollOverCor23_2; //! Polygon for roll over correction
  float fMCHitTMin, fMCHitTMax;     //! Cuts on hit time (for Monte Carlo)
  float fAmpCorr[3];  //! Correction to the amplitude in order to match counterpart (X vs. Y or V vs. U)
  float fCrossTalkParams[2];   //! Cross talk suppression algorithm parameters
  float fTimeCrossTalkParam;   //! Cross talk in multiplexed analog signal suppression algorithm parameter

  CsGEMPlane              *stripplane;        //! Contains hits and clusters
  CsPixelGEMPlane         *pixelplane;        //! Contains hits and clusters
  int fDebugLevel;                  //! Debug level 
  
  bool apv_chan_cals_mapped;        //! flag to indicate mapping of channel calibrations 
  
#define  wirDU_	  wirD_		//!< Wires distance = diff 1st_wire-center
#define  angU_    ang_		//!< Wires angle
#define  sinAngU_ sinAng_	//!< sin( wires angle)
#define  cosAngU_ cosAng_	//!< cos( wires angle)
#define  nWirU_   nWir_		//!< Number of wires
#define  wirPU_   wirP_		//!< Wire pitch

  double wirDV_;		//!< Wires distance    in 2nd dimension
  double angV_;			//!< Wires angle       in 2nd dimension
  double sinAngV_;		//!< sin( wires angle) in 2nd dimension
  double cosAngV_;		//!< cos( wires angle) in 2nd dimension
  int    nWirV_;		//!< Number of wires   in 2nd dimension
  double wirPV_;		//!< Wire pitch        in 2nd dimension

};

#endif // CsPixelGEMDetector_h
