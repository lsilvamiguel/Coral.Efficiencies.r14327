/*!
   \file    CsTrigGroupData.h
   \brief   Container for data from calorimeter's trigger group
   \version $$
   \author  Vladimir Kolosov
   \date    $$
*/

#ifndef CsTrigGroupData_h
#define CsTrigGroupData_h

#include <iostream>
#include "coral_config.h"

////////////////////////////////////////////////////////////////////////////////

///   Container for trigger-group info for external usage (in trigger logic for example)
class CsTrigGroupData
{
  public:
   ///  Base constructor
  CsTrigGroupData   (int l,int x,int y, double e=0) :
                                       fId_layer(l),fId_x(x),fId_y(y),fAmpl(e)
                    {fSumADC=0;fTime1=0; time1_is_valid=false; fTime2=0; time2_is_valid=false;}

   /// Copy constructor
  CsTrigGroupData   (const CsTrigGroupData &tg) {*this=tg;};
 // =========================================================================
 // Operators
 // =========================================================================

  public:

   /// Assignment operator
  CsTrigGroupData          &operator =         (const CsTrigGroupData &tg);

// =========================================================================
// Methods
// =========================================================================
  public:
   ///  Print info
    void        Print   (void) const { std::cout << " Layer " << fId_layer << " x " << fId_x
                           << " y " << fId_y  << " E " << fAmpl << " SumADC " << fSumADC;
                                              if(time1_is_valid) std::cout << " T1 " << fTime1;
                                              if(time2_is_valid) std::cout << " T2 " << fTime2;
                                                                             std::cout << std::endl;}
    /// Clear time and amplitude.
    void                Clear            (void) {fAmpl=0.; fTime1=0.; fTime2=0.; fSumADC=0.;
                                                time1_is_valid =false; time2_is_valid =false;}

   /// Get layer.
    int                 GetLayer         (void) const {return fId_layer;}

   /// Get x index
    int                 GetX             (void) const {return fId_x;}

   /// Get y index
    int                 GetY             (void) const {return fId_y;}

   /// Get energy.
    double              GetAmpl        (void) const {return fAmpl;}

   /// Set amplitude.
    void                SetAmpl          (double amp) {fAmpl=amp;}

   /// Set Summ of ADC amplitudes.
    void                SetSumADC        (double amp) {fSumADC=amp;}

   /// Add amplitudes to Summ of ADC.
    void                AddADC           (double amp) {fSumADC +=amp;}

   /// \return  summ of ADC amplitudes in trigger group.
    double              GetSumADC        (void) const {return fSumADC;}

   /// Set time.
    void                SetTime          (double tim) {fTime1=tim; fTime2=tim;}

   /// Set time1.
    void                SetTime1         (double tim) {fTime1=tim; time1_is_valid = true;}

   /// Set time2.
    void                SetTime2         (double tim) {fTime2=tim; time2_is_valid = true;}

   /// Get time1.
    double              GetTime1         (void) const {return fTime1;}

   /// Get time2.
    double              GetTime2         (void) const {return fTime2;}

   /// \return validity of Time1
    bool                IsTime1Valid     (void) const {return time1_is_valid;}

   /// \return validity of Time1
    bool                IsTime2Valid     (void) const {return time2_is_valid;}

  // ============================================================================
  // Attributes, data
  // ============================================================================

  private:
   ///  Three integers: (layer,x,y ) identify trigger group in calorimeter
    int        fId_layer;
    int        fId_x;
    int        fId_y;
   /// energy deposit in trigger group
    double     fAmpl;
  /// Summ of ADC amplitudes in the group
    double     fSumADC;
   /// time measurement in trigger group
//     double     fTime;
    double     fTime1;
    double     fTime2;
    bool       time1_is_valid;
    bool       time2_is_valid;
};

////////////////////////////////////////////////////////////////////////////////

#endif // CsTrigGroupData_h

