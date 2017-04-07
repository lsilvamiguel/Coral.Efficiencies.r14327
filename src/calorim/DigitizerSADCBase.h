#ifndef DigitizerSADCBase_h
#define DigitizerSADCBase_h
/*!
   \brief   Class to decode time and energy from SADC fit
 */
#include <vector>
#include <iostream>
#include "CsTypes.h"
#include <typeinfo>
class CsCalorimeter;
class DigitizerSADCBase
{
  public:
   enum Type   { SADC=0, MSADC=1};

  protected:
   const static double                     sadc_clock_;
   const static double                     TCS_T0;
  protected:

   struct result_s {
     unsigned char ready; /*! set result ready */
     double ampl; /*! amplitude in ADC channels */
     double dampl; /*! amplitude error */
     double time; /*! signal time */
     double timeFWHM; /*! FWHM time */
     double timeHFW;    /*! Half Front Width */
     double dtime; /*! error of time */
     unsigned int max_pos; /*! possition of maximum sample */
     double base; /*! ped even samples */
     double base_diff; /*! diffrennce between even and odd samples ped */
     double dbase; /*! error on base measurement */
     double chi2; /*! chisquare if fitted */
     double ndf; /*! if fitted */
     double time1; /*! front time */
     double time2; /*! back front time */
// More development parameters for overfolows test
     int  nover;    /*! amount of regions with overfolows */
     double ampl_overflow;  /*! overfolow amplitude recovered */
     double time_overflow;   /*! overfolow time recovered */
   };
   typedef result_s result_t;
  // =========================================================================
  // Constructors & Destructors
  // =========================================================================

  public:
    virtual            ~DigitizerSADCBase() {}
    DigitizerSADCBase ( CsCalorimeter* parent,  Type type );
    virtual void PrintType();
  public:

    //
    // static methods
    //
    static void PrintSample ( const std::vector<uint16> &sample );

  public:
    std::string GetClassName() const {return typeid(*this).name();}
    bool   ResultIsReady    ( void ) const { return result_.ready;}
    /// Print some parameters
    virtual void  PrintSettings  (void) const;
    /// Fit
    virtual bool Fit(  const std::vector<uint16> &sample , const unsigned int icell) {
      mycell_ = icell;
      return FitBase(sample);
   }
    virtual bool FitBase(  const std::vector<uint16> &sample ) { std::cerr <<" DigitizerSADCBase::Fit(sample) is dummy " << std::endl; return false;}

    /// Clear
    virtual void  Clear ( void ) {
      ClearResult();
      csample_.clear();
      overflow_samples_.clear();
      is_noise_by_shape_anal_=false;
    }
    virtual bool  Check ( void ) const { return true;}

    /// \return false for good physics signal, Bad name. Why only shape? Should be virtual?
    bool   IsNoiseByShape       ( void ) const { return is_noise_by_shape_anal_;}

    /// \return Calorimeter name
    std:: string GetParentName( void ) const;

    /// \return amplitude (ADC counts)
    double GetSignal        ( void ) const { return result_.ampl; }

    /// \return time (ns)
    double   GetTime        ( void ) const { return result_.time; }

    /// TODO implement Time/Ampltude Sigma estimator
    double   GetTimeErr        ( void ) const { return sadc_clock_; }

    /// \return pedestal (ADC counts)
    double   GetPed         ( void ) const { return result_.base; }

    /// \return error of pedestal (ADC counts)
    double   GetDPed        ( void ) const { return result_.dbase; }

    /// \return position of maximum sample (index to the vector of samples)
    unsigned GetMaxPosition ( void ) const { return result_.max_pos; }

    /// \return chi2 of fit of signal shape
    double   GetChi2        ( void ) const { return result_.chi2;}

    const Type &GetType     ( void ) const { return type_; }

    const std::vector< double >  &GetCSSample ( void ) const { return csample_;}

    virtual void GetCcSample ( const std::vector<uint16> &sample );   // non-const !!

    double GetSADCClock ( void ) const { return sadc_clock_;}

    virtual bool OverflowDetected( void ) const { return result_.nover;}
    int GetOverflowDigit( void ) const { return overflow_amplitude_;}
    unsigned int DetectOverflow ( const std::vector<uint16> &sample );

  public:
   double GetTrueCalibPed( void ) const { return (ped_odd_old_truecalib_+ped_even_old_truecalib_)/2.;}
   // Substitution in the code of dynamic pedestals into historic calibrations
   double GetDynamicPed( void ) const { return (ped_odd_old_+ped_even_old_)/2.;}
   void  SetCalibInfo  ( double ped_odd, double ped_even, double dpedd, double cfcv ) {
       ped_odd_old_truecalib_ = ped_odd;
       ped_even_old_truecalib_ = ped_even;
       dped_old_truecalib_ = dpedd;
       normref_old_truecalib_ = cfcv;
   };

  protected:
   const CsCalorimeter *  GetParent( void ) const { return parent_;}
   virtual void ClearResult ( void );

  protected:
   CsCalorimeter                   *parent_;
   Type                                 type_;

   bool                                  debug_;
   result_t                              result_;
   std::vector< double >           csample_;
// Important calibration contant which is set after MSADC/SADC settings It is very important to check correct settings!! TODO Implement CheckBaseSettings()
   int                                     overflow_amplitude_;
   std::vector< int >                      overflow_samples_;
   bool                                  is_noise_by_shape_anal_;

// Calibration data with some waht scrued-up usage
   double ped_odd_old_;
   double ped_even_old_;
   double normref_old_;  // ?? Clarify usage!!
// TRUE Calibration data   from recent development
   double ped_odd_old_truecalib_;
   double ped_even_old_truecalib_;
   double dped_old_truecalib_;
   double normref_old_truecalib_;

// Primary setting for pedestals calculation. Dont see the way to ommit.
   int sadc_ped_max_;
   int mycell_;

};


#endif
