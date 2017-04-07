#ifndef CompassSoft_Scaler__include
#define CompassSoft_Scaler__include

#include "ChipCol.h"

namespace CS {

////////////////////////////////////////////////////////////////////////////////

/*! \brief COMPASS Scaler

    Some documentation was extracted from "Compass-Note 2001-8"

    \author Alexander Zvyagin
*/
class Scaler: public ChipCol
{
  //============================================================================
  // Types, constants
  //============================================================================

  public:

    class Digit: public Chip::Digit
    {
      public:

                        Digit                   (const DataID &data_id,const DetID &id,uint32 chan,uint32 value);

        void            Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

        const char*     GetNtupleFormat         (void) const {return "chan:val";}
        std::vector<float> GetNtupleData        (void) const;
        
        /// \return Scaler channels
        uint32          GetChannel              (void) const {return channel;}

        /// Amplitude
        int32           GetValue                (void) const {return value;}
 
        void            SetValue                (int v) {value=v;}

      private:

        /// channel number
        uint16          channel;
        
        uint32          value;
    };

    class Map : public Chip::Map
    {
      public:
                        Map                     (const ObjectXML &o);

        void            AddToMaps               (Maps &maps,DaqOption &option) const;

        /// 0-3
        uint16          frame;
    };

  //============================================================================
  // Constructors and destructor
  //============================================================================

  public:
    
                       ~Scaler                  (void) {}
    
    /*! \brief Base constructor.
        \param buf is non-NULL pointer to buffer where first bytes have format of
               Chip::Header.
    */
                        Scaler                  (void const * const buf,bool copy_buf,DaqOption &opt,DaqEvent &ev) :
                                                  ChipCol(buf,copy_buf,opt,ev)
                                                {}

  private:

    /// You can not use copy constructor.
                        Scaler                  (const Scaler &e);
  
  //============================================================================
  // Operators
  //============================================================================

  private:

    /// You can not use assignment operator.
    Scaler             &operator =              (const Scaler &e);

  //============================================================================
  // Methods
  //============================================================================

  public:
  
    /// Print chip information
    void                Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

    /// \return Chip's name. This is "Scaler".
    std::string         GetName                 (void) const {return "Scaler";}

    /// Decode data and \b add new digits to \b digits_list.
    void                Decode                  (const Maps &maps,Digits &digits_list,DaqOption &opt);

  //============================================================================
  // Attributes
  //============================================================================

  private:
};

////////////////////////////////////////////////////////////////////////////////

} // namespace CS

#endif // CompassSoft_Scaler__include
