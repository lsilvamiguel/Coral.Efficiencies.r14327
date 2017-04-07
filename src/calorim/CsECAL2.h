/*!
   \file    CsECAL2.h
   \brief   Class for COMPASS ECAL1 calorimter
   \version $Revision: 1.25 $
   \author  Vladimir Kolosov
   \author  Alexander Zvyagin
   \date    $Date: 2010/07/13 16:18:32 $
*/

#ifndef CsECAL2_h
#define CsECAL2_h

#include <string>
#include <iostream>
#include <cstdio>

#include "CsCalorimeter.h"

////////////////////////////////////////////////////////////////////////////////

/*! \brief COMPASS ECAL1 calorimeter class
*/
class CsECAL2 : public CsCalorimeter
{
  // ============================================================================
  // Constructors, destructor
  // ============================================================================

  public:

    /// Destructor
    virtual            ~CsECAL2                 (void) {}

    /// Construct calorimeter with geometry description.
                        CsECAL2                 (const std::string &name,
                                                 const std::string &geom_file);

  private:

    /// You can not use copy constructor (yet)
                        CsECAL2                 (const CsECAL2 &c);

  // ============================================================================
  // Operators
  // ============================================================================

  private:

    /// You can not use assignment operator (yet)
    CsECAL2            &operator =              (const CsECAL2 &c);

  // ============================================================================
  // Methods
  // ============================================================================

  public:

    /// ECAL2 initialization. Should not be called from constructor.
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
    virtual CsCalorimeter::CalID       GetCalID() const { return ECAL2; };

    virtual int StoreSADCDigit(const int& icell, const CS::ChipSADC::Digit& digit, const bool& isLed);

};

////////////////////////////////////////////////////////////////////////////////

#endif // CsECAL2_h
