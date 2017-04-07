#ifndef CsECAL0_h
#define CsECAL0_h

#include <string>
#include <iostream>
#include <cstdio>

#include "CsCalorimeter.h"

////////////////////////////////////////////////////////////////////////////////

/*! \brief COMPASS ECAL0 calorimeter class
*/
class CsECAL0 : public CsCalorimeter
{
  // ============================================================================
  // Constructors, destructor
  // ============================================================================

  public:

    /// Destructor
    virtual            ~CsECAL0                 (void) {}

    /// Construct calorimeter with geometry description.
                        CsECAL0                 (const std::string &name,
                                                 const std::string &geom_file);

  private:

    /// You can not use copy constructor (yet)
                        CsECAL0                 (const CsECAL0 &c);

  // ============================================================================
  // Operators
  // ============================================================================

  private:

    /// You can not use assignment operator (yet)
    CsECAL0            &operator =              (const CsECAL0 &c);

  // ============================================================================
  // Methods
  // ============================================================================

  public:

    /// ECAL0 initialization. Should not be called from constructor.
    void        Initialize              (void);
    /// Default options initialization.
    void               InitOptions            ( void );
    ///  Set default values for SADC decoding
    void        SetDefaultSADCDigitizationParameters  (void);

    void        FillHistoLED4SADC     (void);
    void        EndOfJob              (void);

    /// Fill electron calibrations histograms
    void        CalibrateAtElectronBeam  (void);

    /// Special electron calibrations
    void        CalibrateAtElectronBeamSpecial  (void);

    /// Draw histograms for one cell
    void        DrawCellHisto           (int icell);
    bool        GoodZone4Electrons       ( double x, double y, double e)  const;
    virtual CsCalorimeter::CalID       GetCalID() const { return ECAL0; };

    virtual int StoreSADCDigit(const int& icell, const CS::ChipSADC::Digit& digit, const bool& isLed);

};

////////////////////////////////////////////////////////////////////////////////

#endif // CsECAL0_h
