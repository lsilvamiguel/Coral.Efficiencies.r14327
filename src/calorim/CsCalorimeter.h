/*!
   \file    CsCalorimeter.h
   \brief   Base class for COMPASS calorimeters
   \version $Revision: 1.64 $
   \author  Vladimir Kolosov
   \author  Alexander Zvyagin
   \author  Denis Murashev
   \date    $Date: 2010/11/29 19:16:06 $
*/

#ifndef CsCalorimeter_h
#define CsCalorimeter_h

#include <cmath>     // for isnormal()
#include <iostream>
#include <string>

#include "coral_config.h"

#include "Reco/Calorimeter.h"
#include "CsDet.h"
#include "CsHist.h"
#include <CsRandom.h>

#include "DaqDataDecoding/Chip.h"
#include "DaqDataDecoding/ChipADC.h"
#include "DaqDataDecoding/ChipSADC.h"
#include "DaqDataDecoding/ChipF1.h"
#include "CsTrigGroupData.h"


class MyPulse;
class CsCalorimeterHist;

class TestHistoSADC;
class TestHistoSADCMore;
class TestHistoSADC_LED;
class DigitizerSADCBase;
class CsDigitizerSADC;
class DevMonitorFEM;
class ManagerShapeTableSADC;

class DigitizerStore
{
  public:
      DigitizerStore ( size_t ncells ) {vsadc.assign(ncells, NULL);}
      void Clear ( void );
  public:
   std::vector < DigitizerSADCBase *> vsadc;
};

class CsCalorimeter : public Reco::Calorimeter, public CsDet
{
  // ============================================================================
  // Types, constants
  // ============================================================================

  public:

    enum CalID {UNKNOWN = 0,ECAL1 = 1, ECAL2 = 2, HCAL1 = 3, HCAL2 =4, ECAL0 = 5 };
    enum ResponseMC {MCExact,MCSmeared};
    enum ROamplitudeFE { UNKNOWN_FEA = 0, FIADC_FEA = 1, SADC_FEA = 2, MSADC_FEA = 3 };

    struct CalorimeterMCData { unsigned int cell_id; double dE;
                               double dT; unsigned int track_id; };

    /*! \brief  Each CsCalorimeter can be divided into number of cells groups
       Sum of amplitudes from the group and timing can be used in trigger logic.
       As by-product timing is used to assign time to reconstructed particle.
    */
    class TrigGroup : public Reco::Calorimeter::SubSet
    {
      // ============================================================================
      // Types, constants
      // ============================================================================
      public:
      enum                GroupInfoType            {CALIB=0,TIME1=1,TIME2=2};

      // =========================================================================
      // Constructors and destructor
      // =========================================================================

      public:

        /// Destructor
        virtual            ~TrigGroup   (void) {}

        ///  Base constructor

                            TrigGroup   (Calorimeter *c, const std::string &name,
                                         const std::vector<size_t> &tg,
                                         int g_l, int g_x, int g_y);

        /// Copy constructor
                            TrigGroup   (const TrigGroup &tg);

      // =========================================================================
      // Operators
      // =========================================================================

      public:

        /// Assignment operator
        TrigGroup          &operator =         (const TrigGroup &tg);

        /// Test for equality
        bool                operator ==        (const TrigGroup &tg);

        /// Print TrigGroup info to output stream
        friend std::ostream     &operator <<        (std::ostream &o,TrigGroup &tg);

      // =========================================================================
      // Methods
      // =========================================================================

      public:

        /// Clear time and amplitude.
        void      Clear            (void) {trig_group_data.Clear();}


      // =========================================================================
      // Data Members
      // =========================================================================

      public:

        /// Info about about trigger group
        CsTrigGroupData          trig_group_data;

    };  // end class TrigGroup


 protected:

  class FitResultSADC
  {
    public:

    double pedestal_;
    double maximum_;
    double ch_max_;
  };

  // ============================================================================
  // Constructors, destructor
  // ============================================================================

  public:

    /// Destructor
    virtual            ~CsCalorimeter           (void) {}

    /*! \brief Construct calorimeter from given COMGEANT file with geometry description.

        \param name             The calorimeter name (EC01P1__,EC01P2__,HC01P1__,HC01P2__)
        \param geom_file        File with geometry description

        For COMGEANT geometry version 6.03 the short GEANT name of the calorimeter will be
        detected automaticaly from geom_file.
    */
                        CsCalorimeter           (const std::string &name,
                                                 const std::string &geom_file);

  private:

    /// You can not use copy constructor (yet)
                        CsCalorimeter           (const CsCalorimeter &c);

  // ============================================================================
  // Operators
  // ============================================================================

  private:

    /// You can not use assignment operator (yet)
    CsCalorimeter      &operator =              (const CsCalorimeter &c);

    /// Print Calorimeter info to output stream
    friend std::ostream     &operator <<             (std::ostream &o,const CsCalorimeter &c);

    /// Get Calorimeter info from input stream
    friend std::istream     &operator >>             (std::istream &in, CsCalorimeter &c);


  // ============================================================================
  // Methods
  // ============================================================================

  public:
    virtual void        Initialize               (void);
    /// Default options initialization.
    virtual void        InitOptions              (void);
    virtual void        ReadOptions              (void);

  public:
    virtual bool       Reconstruction         (void);
    virtual void       EndOfJob               (void) { Reco::Calorimeter::EndOfJob(); }

    bool                AddMCHit                (int detector_number,const void *data);

    ///  Propagate R/O infor from DaqDataDecoding(mapping files) to CsCalorimeters
    virtual void        SetDaqDataDecodingInfoGeneral  ( const std::multimap<CS::Chip::DataID,CS::Chip::Digit*> &daqmap );
    virtual void        SetDaqDataDecodingInfo  ( const CS::Chip::Digit &d );
    void                 SetEventIDInfo  ( void );

    ///  Raw data decoding CsCalorimeters
    virtual void        GetSADCDigitizationOptions            (void);
    ///  Set default values for SADC decoding
    virtual void        SetDefaultSADCDigitizationParameters  (void);
    void                 PrintSADCSettings  (void) const;


    virtual void        DecodeChipDigits            (const CS::Chip::Digits& digits);
    /*! \brief Decoding of a special calibration LED or LASER event.
               This function is not called in CORAL decoding schema so in general we
               can not trust in LED/LASER sadc amplitudes obtained in standard CORAL.
               Actually there ia a way to obtain correct LED/LASER digitization in
               CORAL by replacing DecodeChipSADCDigit(d, false); in CsCalorimeter::DecodeChipDigit
               by something like  DecodeChipSADCDigit(d, isLED); with correct isLED flag.
               But probably speed of decoding will suffer a bit. Need to be measured.
               On top of that DecodeChipDigitsLEDEvent perform decoding of
               LED/LASER system monitoring devices like PINdiod in HCAL1 or 8 FEMs in ECAL1.
    */
    virtual void        DecodeChipDigitsLEDEvent    (const CS::Chip::Digits& digits);
//             bool        ApplyCorrectionsMonitorFEM  ( void );

    virtual void        DecodeChipDigit             (const CS::Chip::Digit& digit);
    virtual void        DecodeChipDigitLEDEvent     (const CS::Chip::Digit& digit);

    virtual void        DecodeChipSADCDigit         (const CS::ChipSADC::Digit& digit, const bool& isLed=false);
    virtual  int        StoreSADCDigit              (const int& icell, const CS::ChipSADC::Digit& digit, const bool& isLed);

    virtual void        DecodeChipADCDigit          (const CS::ChipADC::Digit& digit);

            void        DecodingTestSADC         (void);

    virtual void        TestHistoRawInfo        (void);

 public:
    virtual CsCalorimeter::CalID       GetCalID() const { return UNKNOWN; };


    virtual void        Clear                   (void);
    void                 ClearDigitizersSADC     (void);

    void                 ProcLED                 (void);

    void                CleanCsDigits           (void);

    void                 readCalibration         (time_t timePoint);
    void                 readMCCalibration       (time_t timePoint);

    void                 WriteCalib              (void);

    virtual void         CreateGUI               ();

    /// Draw histograms for one cell
    virtual void       DrawCellHisto           (int icell);

    virtual void        makeMCDecoding          (void);
            void        MakeMCDigitization      (void);

    // overwrite Reco's random number getters
    // use CORAL facilities if Reco us used in CORAL
    virtual double      GetRandomFlat(void) const {
        return CsRandom::flat();
    }
    virtual double      GetRandomGaus(void) const {
        return CsRandom::gauss();
    }

    /// unused dummy function.  Needs to be defined because CsDet declares
    /// an abstract prototype.
    virtual void        clusterize              (void) {};

    /// \return vector of trigger-group data above some energy
    std::vector<CsTrigGroupData> GetTrigGroupData     (double e) const;

    /// \return vector of trigger-group data above some energy at (x,y) point
    std::vector<CsTrigGroupData> GetTrigGroupDataXY   (double e,double x,double y) const;

    /// \return position of the group in trigger_groups vector according to {layer,x,y} numbers
    virtual int        GetTrigGroupNumber       (int layer, int x, int y) const ;

    /*! \brief  CsCalorimeter Test  Function for trigger and timing stuff of the calorimeter */
    void                TrigGroupTest            (void);

    /// Formatting information about NEW trigger groups info for calibration
    int                InputTrGrCalibInfo      (const std::string &s);
    /// Get Formatting information about OLD trigger groups info for for calibration
    int                OutputTrGrCalibInfo     (std::string& s) const;
    /// Formatting information about NEW trigger groups info for Time calibration
    int                InputTrGrTimeCalibInfo  (const std::string &s);
    /// Get Formatting information about OLD trigger groups info for for Time calibration
    int                OutputTrGrTimeCalibInfo (std::string& s) const;

    ///  On-line Time calibration function
    void                CalibrateTrigGroupTime  (void);

    virtual int           InputAnyCalibInfo(const std::string &tag,const std::string &s);
    virtual int           OutputAnyCalibInfo(const std::string &tag, std::string &s,const std::string &comment ="") const;
    virtual int           InputCalibInfo(size_t when, const std::string &s);

    /// Get Formatting information about OLD cells_info for calibration
    virtual int           OutputCalibInfo     ( std::string& s ) const;
    /// Formatting information about NEW cells_info for LEDs
    virtual int           InputLEDInfo      (size_t when, const std::string &s);
    /// Get Formatting information about OLD cells_info for LEDs
    virtual int           OutputLEDInfo     ( std::string& s, const std::string &comment = "" ) const;
    /// Get Formatting information about cells thresholds
    virtual int           InputEcutCellsInfo      (const std::string &s);
    virtual int           OutputEcutCellsInfo     ( std::string& s ) const;
    virtual int           InputTimeCalibInfo      (size_t when, const std::string &s);
    virtual int           InputWidthCalibInfo     (size_t when, const std::string &s);
    virtual int           OutputTimeCalibInfo     ( std::string& s ) const;
    int                    OutputTimeCalibInfo4LED ( std::string& s ) const;
    int                    OutputTimeCalibInfoXY4LED ( std::string& s ) const;

    virtual int           InputSADCInfo  ( const std::string &s);
    virtual int           OutputSADCInfo ( std::string &s ) const;

    virtual int           InputTimeFWHMSADC ( const std::string &s);
    virtual int           InputShapeTableSADC ( const std::string &s);

  public:
    ///  Make CsDigits and store them in CsDet list <*CsDigit> container
    void                MakeCsDigits           ( void );
// This is obsolete function
//    void                SpecialTask_ScaleCalibrationsByDigMaxAd (void);
     void               UpdateFrontEndDependentSettings( void );

    void                SetTCSPhase    ( double tcs_phase ) { tcs_phase_ = tcs_phase; }
    double              GetTCSPhase    ( void ) const {
      if( !std::isnormal(tcs_phase_) )
        throw Reco::Exception( "%s: TCS phase not initialized!", GetName().c_str() );
      return tcs_phase_;
    }
    void                SetTimeInSpill ( double tis ) { time_in_spill_ = tis; }
    double              GetTimeInSpill ( void ) const {
      if( !std::isnormal(time_in_spill_) )
        throw Reco::Exception( "%s: Time in spill not initialized!", GetName().c_str() );
      return time_in_spill_;
    }

  public:
    double GetFFTwidth(const int icell) const;
    ///  Workaround function to distingwish between SADC and MSADC. It is also used in LEDs data processing, so defined as public.
    bool                       IsItMSADCSrcID    (int src_id) const;

    const std::vector<CS::uint16>* GetSADCSample ( int icell ) const;
  protected:

    void                        InitDefaultReadOut    ( void );
    virtual void               InitReadOut    ( void );

    void                        InitArraySADC    ( void );
    DigitizerSADCBase * GetDigitizerSADC2Modify ( size_t icell );
    const DigitizerSADCBase * GetDigitizerSADC ( size_t icell ) const;
    const DigitizerSADCBase * GetLedDigitizerSADC ( size_t icell ) const;

    std::vector<DigitizerSADCBase *> &GetDigitizersSADC( void );
    std::vector<DigitizerSADCBase *> &GetLedDigitizersSADC( void );
    const std::vector<DigitizerSADCBase *> &GetDigitizersSADCs( void ) const;
    const std::vector<DigitizerSADCBase *> &GetLedDigitizersSADCs( void ) const;

    void AddDigitizerSADC( DigitizerSADCBase * dig, size_t icell );
    void AddLedDigitizerSADC( DigitizerSADCBase * dig, size_t icell );

  private:

    void                ReadGeom                (const std::string &det_name,const std::string &geom_file);

    void                ReadShowerProfile       (Reco::CellType& cellType);

  // ============================================================================
  // Attributes, data
  // ============================================================================

 protected:
    // options
    bool skip_decoding_;
    bool skip_led_decoding_;
    bool skip_cs_time_corrections_;
    bool make_tcs_corrections_;
    bool new_sadc_decoding_;
    bool make_sadc_cluster_filtering_;
    bool make_sadc_histo_;
    bool correct_leds_by_fem_signal_;



    /*! \brief First cell number from detectors.dat file for given cell matrix

      The size of this list is equal to size of Calorimeter::matrixes.
    */
    std::vector<size_t>      comgeant_first_cell_n;

    /// \return Number of Trigger Groups in Calorimeter
    size_t              NGroups          (void) const {return trigger_groups.size();}

    /// Trigger groups
    std::vector<TrigGroup>                  trigger_groups;

    /// Trigger groups calibration information
    std::vector<Reco::StatInfo>             group_info[3][2];


    // SADC digits and parameters

    bool sadc_readout_[5000];
    /// Referense to SADC digits
    std::map < int, const std::vector<CS::uint16> *  > sadc_samples_;

    ManagerShapeTableSADC *mantab_sadc_;
    int sadc_format_version_;
    int sadc_decode_version_;
    int sadc_front_;
    int sadc_front_gate_;
    int sadc_ped_min_;
    int sadc_ped_max_;
    int sadc_signal_min_;
    int sadc_signal_max_;
    int sadc_delta_;
    double sadc_clock_;
    double coeff_convert2max_;
    double sadc_max_position_;
    double tcs_phase_;
    double time_in_spill_;
    std::vector < double >  sadc_shape1000_;
    std::vector<double> FFTwidth_;

    double raw_amp_cut_delta_;
    double raw_time_cut_min_;
    double raw_time_cut_max_;
    bool sadc_shape_filter_apply_;
    bool sadc_shape_filter_slc_use_;
    bool sadc_shape_filter_line_fit_use_;

    TestHistoSADC*         test_histo_sadc_;
    TestHistoSADCMore*     test_histo_sadc_more_;
    TestHistoSADC_LED*     test_histo_sadc_led_;

    std::vector < size_t > subsets_list2store_sadc_histo_;

//    std::vector < DigitizerSADCBase * > dig_sadc_;
//     std::map <int,DigitizerStore *> dig_sadc_;
    DigitizerStore* dig_sadc_;

//    std::vector < DigitizerSADCBase * > dig_led_sadc_;
    DigitizerStore* dig_led_sadc_;

    /// The array is initialized after maps reading and serve to identify
    // amplitude read-out front-end attached to the cell
    std::vector < ROamplitudeFE >  read_out_fe_amp_;
    std::vector < int >  msadc_srcids_;

 private:

    /// Level of histograms creation.  Non-positive values mean "DO NOT CREATE ANY HISTOGRAMS".
    int                 hists_book_level;

    /// option to delete CsDigits which must be normaly done im Coral
    bool                clean_cs_digits_;
    bool                add_sadcinfo_to_cs_digits_;
    bool                add_all_sadcinfo_to_cs_digits_;
    bool                add_timeinfo_to_cs_digits_;
    bool                skip_make_cs_digits_;

    /// A test histogram.
    CsHist1F*           hist_MC_hitE;
    CsHist1F*           hist_MC_hitN;
    CsHist1F*           hist_MC_hitT;
    CsHist2F*           hist_MC_hitNcellMCTrackID;

    /// Some more test histograms
    CsCalorimeterHist*  cs_calo_hist;

    ResponseMC          response_mc;

};  // end class CsCalorimeter

////////////////////////////////////////////////////////////////////////////////

#endif // CsCalorimeter_h

