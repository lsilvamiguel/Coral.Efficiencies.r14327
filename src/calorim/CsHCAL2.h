/*!
   \file    CsHCAL2.h
   \brief   Class for COMPASS HCAL2 calorimeter
   \version $Revision: 1.20 $
   \author  Vladimir Kolosov
   \author  Alexander Zvyagin
   \date    $Date: 2010/05/27 19:49:03 $
*/

#ifndef CsHCAL2_h
#define CsHCAL2_h

#include <string>

#include "CsCalorimeter.h"

////////////////////////////////////////////////////////////////////////////////

/*! \brief COMPASS HCAL2 calorimeter class
*/
class CsHCAL2 : public CsCalorimeter
{
  // ============================================================================
  // Constructors, destructor
  // ============================================================================

  public:

    /// Destructor
    virtual            ~CsHCAL2                 (void) {}

    /// Construct calorimeter with geometry description.
                        CsHCAL2                 (const std::string &name,
                                                 const std::string &geom_file);

  private:

    /// You can not use copy constructor (yet)
                        CsHCAL2                 (const CsHCAL2 &c);

  // ============================================================================
  // Operators
  // ============================================================================

  private:

    /// You can not use assignment operator (yet)

    CsHCAL2            &operator =              (const CsHCAL2 &c);

  // ============================================================================
  // Methods
  // ============================================================================

  public:

    /// HCAL2 initialization. Should not be called from constructor.
    void                Initialize              (void);
    /// Default options initialization.
     void               InitOptions             (void);
    /// Implementation of raw data decoding in HCAL2
    void        DecodeChipDigit          (const CS::Chip::Digit &digit);
    void        DecodeTriggerGroups      (const CS::Chip::Digits &digits);

    ///  Decode Trigger Group Information
    void           DecodeTriggerGroups  (void);

    virtual CsCalorimeter::CalID       GetCalID() const { return HCAL2; };

  // ============================================================================
  // Attributes, data
  // ============================================================================
  protected:

  ///  Offset for 2x2 trigger groups
  int                    offset_2x2_trig;
  ///  Offset for Layer1 trigger groups
  int                    offset_layer1_trig;
  ///  Offset for Layer2 trigger groups
  int                    offset_layer2_trig;
  ///  Offset for Layer3 trigger groups
  int                    offset_layer3_trig;
  ///  Offset for Layer4 trigger groups
  int                    offset_layer4_trig;
};

////////////////////////////////////////////////////////////////////////////////

#endif // CsHCAL2_h
