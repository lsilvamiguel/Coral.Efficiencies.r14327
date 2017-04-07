// $Id: CsRICH1Detector.h,v 1.28 2010/01/24 16:10:41 suhl Exp $

/*!
   \file    CsRICH1Detector.h
   \brief   Prototype for Compass RICH1 detector Class.
   \author  Benigno Gobbo
   \version $Revision: 1.28 $
   \date    $Date: 2010/01/24 16:10:41 $
*/

#ifndef CsRICH1Detector_h
#define CsRICH1Detector_h

#include "coral_config.h"
#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Matrix/Matrix.h>
#include <list>
#include <vector>
#include <string>

#include "CsDet.h"

typedef int iboramultiarray[16][72][72];
typedef float fboramultiarray[16][72][72];

class CsZone;
class CsMCHit;
class CsDigit;
class CsMCDigit;


/*! \class CsRICH1Detector 
    \brief Compass RICH1 detector Class.
*/

class CsRICH1Detector : public CsDet {

 public:


  /*! \fn CsRICH1Detector( const int id );
    \brief Constructor for RICH1 detector.
    \param id    Detector identification number
    \param TBname  Detector technical board name
  */
  CsRICH1Detector( const int id, const std::string &TBname );

  private:
  //CsRICH1Detector( const CsRICH1Detector& ); //!< Copy constructor

  //CsRICH1Detector& operator=( const CsRICH1Detector& ); //!< Assign operator
  
  public:

  bool operator==( const CsRICH1Detector& ) const; //!< "equal to" operator

  /*! \fn inline void clearDigitsList();
    \brief Clear the list of detector digits.
  */
  inline void clearDigitsList() { decodingDone_ = false; myDigits_.clear(); }
 
  /*! \fn void addMCHit( CsMCHit& hit );
    \brief Add a Monte Carlo Hit to this detector.
    \param hit the hit to be added
  */
  void addMCHit( CsMCHit& hit );

  /*! \fn void clearMCHitList();
    \brief Clear the list of detector Monte Carlo hits.
  */
  void clearMCHitList();

  //! Returns list of pointers to the digits associated to this detector
  inline std::list<CsMCHit*> getMyMCHits() { return( myMCHits_ ); } 

  //! Decode the MC data for this detector
  void makeMCDecoding();

  void clusterize() {}

  //! \c true if this detector data was arleady decoded
  bool decoded() const { return( decodingDone_ ); }
  
  //! set this detector as decoded
  void setDecodingDone() { decodingDone_ = true; }

  //! set decoding on this detector
  void setDecode() { decode_ = true; }

  // Returns Cathode number
  inline int getCathode(int address) const { return( address>>20 ); }

  //Returns X-component of the fired pad address
  inline int getPadX(int address) const { return( address & 0x3FF ); }

  //Returns Y-component of the fired pad address
  inline int getPadY(int address) const { return( (address>>10) & 0x3FF ); }

  //Returns pad address
  void getPadADR( CsDigit& digit, int& Cathode, int& ix, int& iy ) const ;

  //encode pad address
  inline  int setPadADR( const int& Cathode, const int& ix, const int& iy ) 
          const { return  (ix + (iy<<10)+ (Cathode<<20)); }

  //! Returns detector centre Z position (MRS) (mm)
  double getZcm() const;

  void setPhotonDet( std::string, double, double, double, CLHEP::HepMatrix );

  void setCathode( int, std::string, std::string, double, double, double, int, int,
		   double, double, double, double, CLHEP::HepMatrix );

  void setMirrNom( std::string, double, double, double, double );

  void setMirrEle( std::string, double, double, double, double, double, double,
                   double, int );
  
  void setMCCFRefInd( double );
  void setMCCFRefIndVS( double );

  virtual void getMCDigits( std::list<CsMCHit*>& );

  inline  int nDetector() { return nDetector_; };

  inline  bool getCalibFlag() const { return calib_flag; };
  inline  float getIndex() const { return index_; };
  inline  float getIndexUV() const { return indexUV_; };
  inline  float getIndexVS() const { return indexVS_; };
  inline  float getThresh() const { return thresh_; };

  void print(); 

  /// Decode raw data
  virtual void       DecodeChipDigits        (const CS::Chip::Digits &digits);
  virtual void       DecodeChipDigit         (const CS::Chip::Digit &digit);

  /// Read calibration data.
  void          readCalibration         (time_t);

  //! True if BORAs multi-arrays were set
  bool borasMultiarrayOK(void) { return _borasOK; }

  //! Get reference to BORAs threshold multi-array
  iboramultiarray &getThresholds(void) { return _threshold; }

  //! Get reference to BORAs pedestals multiarray
  fboramultiarray &getPedestals(void) { return _pedestal; }

  //! Get reference to BORAs sigmas multiarry
  fboramultiarray &getSigmas(void) { return _sigma; }

  //! Get a "geiod"th BORA voltage[i]
  int getBoraVoltage( int geoid, int i ); 

  //! Get a "geiod"th BORA temperature[i]
  int getBoraTemperature( int geoid, int i );  

 /* Andrea ... way to access private array for calibrations */

  bool T0_flag(void) { return _calibT0_flag; } // true if read from DB
  int getPMT_T0( int geoid, int chipid, int chanid );  

 /* Andrea ... end of way to access private array for calibrations */

  //! Set bora calibration data
  bool readBorasTables( int run );

 private:
  
  int            nPhotChmb_;
  int            nDetector_;
  bool           decodingDone_;    //!< \c true id Detector decoded
  bool           decodeCard_;      //!< Decoding set by cards
  bool           decode_;          //!< Decoding set at run time
  std::list<CsMCHit*> myMCHits_;   //!< MC hits associated to this det

  std::vector<float> calib_data;   //!< Calibration data from CDB
  std::vector<int>   calibT0_data;   //!< Calibration data from CDB
  bool          calib_flag;        //!< true if data from CDB
  double index_;
  double indexUV_;
  double indexVS_;
  double thresh_;

  /* Andrea ... private array for calibrations */

  bool              _calibT0_flag;        //!< true if data from CDB
  int  _pmt_t0calarray[160][8][8];

  /* Andrea ... end of private array for calibrations */

  bool              _borasOK;      //!< true if BORAs multi-arrays are set
  iboramultiarray   _threshold;    //!< RICH1 BORAs thresholds
  fboramultiarray   _pedestal;     //!< RICH1 BORAs pedestals
  fboramultiarray   _sigma;        //!< RICH1 BORAs sigmas
  int   _voltage[248][5];          //!< RICH1 BORAs voltages
  int   _temperature[248][5];      //!< RICH1 BORAs temperatures
};

#endif // CsRICH1Detector_h
