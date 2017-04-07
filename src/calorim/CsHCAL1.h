/*!
   \file    CsHCAL1.h
   \brief   Class for COMPASS HCAL1 calorimeter
   \version $Revision: 1.19 $
   \author  Vladimir Kolosov
   \author  Alexander Zvyagin
   \date    $Date: 2010/05/07 15:09:13 $
*/

#ifndef CsHCAL1_h
#define CsHCAL1_h

#include <string>

#include "CsCalorimeter.h"

////////////////////////////////////////////////////////////////////////////////

/*! \brief COMPASS HCAL1 calorimeter class
*/
class CsHCAL1 : public CsCalorimeter
{
  // ============================================================================
  // Constructors, destructor
  // ============================================================================

  public:

    /// Destructor
    virtual            ~CsHCAL1                 (void) {}

    /// Construct calorimeter with geometry description.
                        CsHCAL1                 (const std::string &name,
                                                 const std::string &geom_file);

 private:

    /// You can not use copy constructor (yet)
                        CsHCAL1                 (const CsHCAL1 &c);

  // ============================================================================
  // Operators
  // ============================================================================

 private:

    /// You can not use assignment operator (yet)

    CsHCAL1            &operator =              (const CsHCAL1 &c);

  // ============================================================================
  // Methods
  // ============================================================================

 public:


    /// HCAL1 initialization. Should not be called from constructor.
    void                Initialize              (void);
    /// Default options initialization.
     void               InitOptions             (void);
    /// Implementation of raw data decoding in HCAL1
    void        DecodeChipDigit          (const CS::Chip::Digit &digit);


    ///  Decode Trigger Group Information
    void        DecodeTriggerGroups             (void);

    CsCalorimeter::CalID       GetCalID() const { return HCAL1; };

  // ============================================================================
  // Attributes, data
  // ============================================================================

 private:

  ///  Offset for 2x2 trigger groups
  int                    offset_2x2_trig;
  ///  Offset for Layer1 trigger groups
  int                    offset_layer1_trig;
  ///  Offset for Layer2 trigger groups
  int                    offset_layer2_trig;
};

////////////////////////////////////////////////////////////////////////////////

#endif // CsHCAL1_h
