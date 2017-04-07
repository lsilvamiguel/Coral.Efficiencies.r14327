/*!
   \file    CsECAL1.h
   \brief   Class for COMPASS ECAL1 calorimter
   \version $Revision: 1.16 $
   \author  Vladimir Kolosov
   \date    $Date: 2010/05/27 16:17:05 $
*/

#ifndef CsECAL1_h
#define CsECAL1_h

#include <iostream>
#include "coral_config.h"
#include "CsCalorimeter.h"

////////////////////////////////////////////////////////////////////////////////

/*! \brief COMPASS ECAL1 calorimeter class
*/
class CsECAL1 : public CsCalorimeter
{
  // ============================================================================
  // Constructors, destructor
  // ============================================================================

  public:

    /// Destructor
    virtual            ~CsECAL1                 (void) {}

    /// Construct calorimeter with geometry description.
                        CsECAL1                 (const std::string &name,
                                                 const std::string &geom_file);

    /// ECAL1 initialization. Should not be called from constructor.
    void                Initialize              (void);
    /// Default options initialization.
     void               InitOptions            ( void );
    ///  Raw data decoding CsCalorimeters
     void               DecodeChipDigits        (const CS::Chip::Digits &digits);
     void               DecodeChipDigitsLEDEvent(const CS::Chip::Digits &digits);

     void               SetDaqDataDecodingInfo  ( const CS::Chip::Digit &d );

    /// Create GUI
    virtual void         CreateGUI               ();

  private:

    /// You can not use copy constructor (yet)
                        CsECAL1                 (const CsECAL1 &c);

  // ============================================================================
  // Operators
  // ============================================================================

  private:

    /// You can not use assignment operator (yet)
    CsECAL1            &operator =              (const CsECAL1 &c);

  // ============================================================================
  // Methods
  // ============================================================================

  public:
    bool       Reconstruction         (void);
    void       EndOfJob               (void);
    int        InputJOUJOUInfo(const std::string &s);

    ///  Set default values for SADC decoding
    void        SetDefaultSADCDigitizationParameters  (void);
    /// \return -1 in XY is not valid
    int    GetGAMSCellOfColumnRow (int x, int y) const;
    int    GetMaintzCellOfColumnRow (int x, int y) const;
    int    GetOlgaCellOfColumnRow (int x, int y) const;
    ///  \return  cell name
    std::string GetCellName            (int icell) const;
    /// Formatting information about NEW cells_info for calibration
//     int	           InputCalibInfo      (const std::string &s);
    int            InputCalibInfo      (size_t when, const std::string &s);

    /// Get Formatting information about OLD cells_info for calibration
    int            OutputCalibInfo( std::string &s ) const;

    void        ScaleCalib2006  ( size_t when );
    void        ScaleCalib2007  ( size_t when );

    virtual CsCalorimeter::CalID       GetCalID() const { return ECAL1; };

  protected:
    void               InitReadOut    ( void );

  public:

  // ============================================================================
  // Attributes, data
  // ============================================================================
  size_t index_gams_;
  size_t index_maintz_top_;
  size_t index_maintz_bottom_;
  size_t index_olga_saleve_;
  size_t index_olga_jura_;

  double  cell_sizex_gams_;
  double  cell_sizey_gams_;
  double  cell_sizex_maintz_;
  double  cell_sizey_maintz_;
  double  cell_sizex_olga_;
  double  cell_sizey_olga_;

  private:
  /// For backward compatibility  we keep old style SADC timing(2006,2007 data )
  bool new_style_sadc_timing_;

};

////////////////////////////////////////////////////////////////////////////////

#endif // CsECAL1_h
