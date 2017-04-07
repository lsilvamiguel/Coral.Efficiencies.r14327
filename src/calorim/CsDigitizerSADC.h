/*!
   \file    CsDigitizerSADC.h
   \brief   Class to decode time and energy from SADC fit
   \version $Revision: 1.21 $
   \author  Vladimir Kolosov
   \author  Alexander Zvyagin
   \author  Denis Murashev
   \date    $Date: 2010/09/29 12:05:22 $
*/

#ifndef CsDigitizerSADC_h
#define CsDigitizerSADC_h


class MyPulse;
#include "TF1.h"

#include <math.h>
#include <vector>
#include <string>

#include "CsTypes.h"
#include "Reco/StatInfo.h"
#include "DigitizerSADCBase.h"

/*!
   \brief   Class to decode time and energy from SADC fit
 */

class CsCalorimeter;

/////////////////////////////////////////////////////////////////////////////////////////////////////////
class CsDigitizerSADC : public DigitizerSADCBase
{
  // ============================================================================
  // Types, constants
  // ============================================================================

  public:
   enum Method { MaxSimple=0, MaxAdvanced=1, SummRange=2, ShapeFit=3, FitPulse=4};
//   enum Type   { SADC=0, MSADC=1};
   enum Shape  { RD_ECAL1=0, RD_ECAL2=1, RD_ECAL1_GAMS=2,
                 RD_ECAL1_MAINZ=3, RD_ECAL1_OLGA=4, LED_ECAL2=5,
                 LED_ECAL1_GAMS=6, LED_ECAL1_MAINZ=7, LED_ECAL1_OLGA=8 };
   enum TimeExtractionMethod {
     Edge, HalfMax
   };

   struct cfd_result_s {

     float time; /*! reconstructed cfd time [clk]*/
     unsigned int ampl : 12; /*! estimated amplitude [ADC channel]*/

   };
   typedef cfd_result_s cfd_result_t;

   struct cfd_option_s {
     unsigned int thr : 12; /*! Amplitude threshold to apply [ADC channel]*/
     unsigned int smin; /*! first sample to analyse */
     unsigned int smax; /*! last sample to analyse */
     unsigned int delay : 12; /*! how long to delay the signal */
     float ampl; /* amplification factor */
   };
   typedef cfd_option_s cfd_option_t;


  private:
//    struct result_s {
//      unsigned char ready; /*! set result ready */
//      double ampl; /*! amplitude in ADC channels */
//      double dampl; /*! amplitude error */
//      double time; /*! signal time */
//      double dtime; /*! error of time */
//      unsigned int max_pos; /*! possition of maximum sample */
//      double base; /*! ped even samples */
//      double base_diff; /*! diffrennce between even and odd samples ped */
//      double dbase; /*! error on base measurement */
//      double chi2; /*! chisquare if fitted */
//      double ndf; /*! if fitted */
//    };
//    typedef result_s result_t;

   const static int                        N_SADC_shapes      = 9;
   const static double                     SADC_shape_factor;
//    const static double                     sadc_clock_        = 12.86;
//    const static double                     TCS_T0             = 40.;

   static std::vector< double >            SADC_shape[N_SADC_shapes];
   static std::vector< double >            D_SADC_shape[N_SADC_shapes];
   static int                              SADC_shape_max_position[N_SADC_shapes];
   static int                              SADC_shape_max_from_start[N_SADC_shapes];
   static double                           SADC_shape_max_from_front[N_SADC_shapes];
   static int                              SADC_shape_sadc_front_gate[N_SADC_shapes];
   static std::vector< double >            SADC_shape_value_ower[N_SADC_shapes];
   static std::vector< double >            SADC_shape_index_ower[N_SADC_shapes];
   static std::vector< double >            SADC_tail[2];


  // =========================================================================
  // Constructors & Destructors
  // =========================================================================

  public:
    virtual            ~CsDigitizerSADC() ;

    CsDigitizerSADC ( CsCalorimeter* parent, Method method, DigitizerSADCBase::Type type, Shape table );

  // =========================================================================
  // Methods
  // =========================================================================

  public:

    //
    // static methods
    //
//     static void PrintSample ( const std::vector<uint16> &sample );
    static void InitStatic ( void );
    static double SADCshape( int is, int table  );


    //
    // const methods
    //
    void  PrintSettings  (void) const;

    double Vshape( int is, double max, double ped ) const;
    double DVshape( int is, double max, double ped ) const;

//     bool   ResultIsReady    ( void ) const { return result_.ready;}
//
//     /// \return amplitude (ADC counts)
//     double GetSignal        ( void ) const { return result_.ampl; }
//
//     /// \return time (ns)
//     double   GetTime        ( void ) const { return result_.time; }
//
//     /// \return pedestal (ADC counts)
//     double   GetPed         ( void ) const { return result_.base; }
//
//     /// \return error of pedestal (ADC counts)
//     double   GetDPed        ( void ) const { return result_.dbase; }
//
//     /// \return position of maximum sample (index to the vector of samples)
//     unsigned GetMaxPosition ( void ) const { return result_.max_pos; }
//
//     /// \return chi2 of fit of signal shape
//     double   GetChi2        ( void ) const;
//
//     const Type &GetType     ( void ) const { return type_; }

    // Need it for development
    std::vector< double > GetCShape ( int size ) const;


    /// \return whether the signal was recognized to be from noise
    bool IsNoise      ( void ) const;
    bool Shape_FilterSLC ( const std::vector<uint16> &smpl ) const; // static at the moment

    /// \return whether the signal was recognized to contain pile-up
    bool IsPileup     ( void ) const;

    /// \return whether the signal was recognized to be from laser/LED
    bool IsLED        ( void ) const;

    /// \return whether saturation could be recovered
    bool SatRecovered ( void ) const;

    //
    // non-const methods
    //
    bool Fit(  const std::vector<uint16> &sample , const unsigned int icell);

    bool  LoadSettings ( const std::vector< double > &settings );

    void GetCcSample ( const std::vector<uint16> &sample );   // non-const !!

//    void  Clear ( void ) { ClearResult(); csample_.clear();}
    void  SetDebugMode ( bool mode ) { debug_ = mode; }

    void  StoreStatInfo  (  const std::vector<uint16> &sample );

    cfd_result_t *CFD(const cfd_option_t &options, const std::vector<unsigned int> &sample);

    void  SetCalibInfo  ( double ped_odd, double ped_even, double dpedd, double cfcv )
      {
				  ped_odd_old_ = ped_odd;
				  ped_even_old_ = ped_even;
				  dped_old_ = dpedd;
				  normref_old_ = cfcv;
      };

    Method GetMethod() {
      return method_;
    };

  protected:
//    bool                FilterByShape ( void );

  private:
   const result_t &FitMaxSimple(  const std::vector<uint16> &sample );
   const result_t &FitMaxAdvanced(  const std::vector<uint16> &sample );
   const result_t &FitSummRange(  const std::vector<uint16> &sample );
   const result_t &FitShape(  const std::vector<uint16> &sample );
   const result_t &FitOverflow(  const std::vector<uint16> &sample );
   const result_t &PulseFit( const std::vector<uint16> &sample, const unsigned int icell, const std::string &DAQDetName );

   double CalcTime(TimeExtractionMethod method, const std::vector<short unsigned int> &sample, const unsigned int max_position);
   double CalcTime_edge(const std::vector<short unsigned int> &sample, const unsigned int max_position);
  public:
   double CalcTime_halfMax(const std::vector<short unsigned int> &sample, const unsigned int max_position);


   std::vector<short unsigned int> Correct_Background_Exp(const std::vector<short unsigned int> &_sample);

   void ClearResult();
  // =========================================================================
  // Data Members
  // =========================================================================

  private:

   Method                           method_;
//   Type                             type_;
   bool                             activated_;

   double sadc_delta_;
   int sadc_ped_min_;
//   int sadc_ped_max_;

   int sadc_signal_min_;
   int sadc_signal_max_;
   int sadc_front_;

   double coeff_convert2max_;  // empiric constant
   Shape                          shape_table_;

//   int                            overflow_amplitude_;

//   bool                           debug_;

//   result_t                        result_;
//   std::vector< double >           csample_;

   int                             front_position_;
   double                          total_summ_;
   double                          signal_slope_;
   double                          fit_range_to_front_;
   double                          fit_range_to_back_;
   double                          fit_cut_slope_;
   double                          fit_cut_offset_;
   double                          fit_cut_debug_;
   double                          fit_base_;
   char*                           fit_ofile_;
   double                          fit_sample_cerr_;
   double                          fit_sample_derr_;
   double                          ac_cut_min_;

//   CsCalorimeter                  *parent_;
   MyPulse                        *mypulse_;
   int                             plane_;

  public:
   bool                            option_store_stat_info_;
   Reco::StatInfo                   stat_ped_odd_;
   Reco::StatInfo                   stat_ped_even_;
   Reco::StatInfo                   stat_dped_;
   Reco::StatInfo                   stat_ped_;
   Reco::StatInfo                   stat_time_;

   double dped_old_;
   double normref_old_;

   bool                            use_shape_filter_slc_;
   bool                            use_shape_filter_line_fit_;
};


#endif
