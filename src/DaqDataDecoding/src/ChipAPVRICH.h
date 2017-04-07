#ifndef CompassSoft_ChipAPVRICH__include
#define CompassSoft_ChipAPVRICH__include

#include <iostream>
#include "ChipAPV.h"

namespace CS {

////////////////////////////////////////////////////////////////////////////////

/*! \brief APV chip for RICH detector

    RICH detector consists from 16 photodetectors. Every photodetector is
    independ from others.

    A photodetector has 72x72 pixels which are organized
    in the following way.
    There are 6x8=48 pixels which form a single pad. In a photon detector
    there are 12 pads along X-axis and 9 pads along Y-axis.

    Every pad has a physical connector. A column of 1x9 pads is read
    by RICHAPV card, which has 9 connectors for 9 pads. So in total
    RICHAPV card has to read 9x48=432 channels. To read full photocathode
    it is needed 12 RICHAPV cards.

    RICHAPV card has 4 APV chips on it. A single APV chip can read 128 channels.
    So in total RICHAPV card may read up to 4x128=512. It is more then enough
    to read 432 pixels. That is why not all 128 channels of APV chip are used.
    The first ten channels and the last ten channels are not readout (4x108=432
    which is exactly equal to the number of pixels we have to read).

    \author Alexander Zvyagin
*/
class ChipAPVRICH: public ChipAPV
{
  //============================================================================
  // Types, constants
  //============================================================================

  public:

    enum {APV_read_channels=108};

    /*! \brief APV chip digit
    
    */
    class Digit: public ChipAPV::Digit
    {
      public:
 
        /// Destructor
        virtual        ~Digit                   (void) {}
 
        /// Base constructor
                        Digit                   (const DataID &data_id,const DetID &id,
                                                 uint16 pd_x, uint16 pd_y);

        /// Print the digit.
        void            Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

        /// \return Digit format for ROOT-ntuple creation.
        const char*     GetNtupleFormat         (void) const;
        
        /// \return vector of values that can be used for filling an Ntuple.
        std::vector<float> GetNtupleData        (void) const;
        
        void            SetPDxy                 (uint16 x,uint16 y) {pd_x=x; pd_y=y;}
        void            SetPixelXY              (int16 x,int16 y) {pixel_x=x; pixel_y=y;}
        
        int16           GetPixelX               (void) const {return pixel_x;}
        int16           GetPixelY               (void) const {return pixel_y;}
        
        uint16          GetPDx                  (void) const {return pd_x;}
        uint16          GetPDy                  (void) const {return pd_y;}
        
        void            SetOldAPVBoard          (bool f) { old_apv_board=f; }
        bool            GetOldAPVBoard          (void) const {return old_apv_board;}
        
        void            SetCardRotation         (char _card_rotation) {card_rotation=_card_rotation;}
        char            GetCardRotation         (void) const {return card_rotation;}
        
      private:
        
        char            card_rotation;
        
        uint16          pd_x;
        uint16          pd_y;
        
        int16           pixel_x;
        int16           pixel_y;

        bool            old_apv_board;
    };

    ////////////////////////////////////////////////////////////////////////////

    /// ChipAPV map
    class Map : public Chip::Map
    {
      public:

        /// Constructor
                        Map                     (const ObjectXML &o);
                        
        void            Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

        void            AddToMaps               (Maps &maps,DaqOption &option) const;

      public:

        char            card_rotation;
        uint16          adc_id;
        uint16          apv_id;
        uint16          pd_x;
        uint16          pd_y;
        bool            old_apv_board;

        ///! Do we always have common mode noise correction in data?
        bool            always_cmc;
    };

    ////////////////////////////////////////////////////////////////////////////

  //============================================================================
  // Constructors and destructor
  //============================================================================

  public:
    
    /// Destructor
                       ~ChipAPVRICH             (void) {}
    
    /*! \brief Base constructor.
        \param buf is non-NULL pointer to buffer where first bytes have format of
               Chip::Header.
    */
                        ChipAPVRICH             (void const * const buf,bool copy_buf,DaqOption &opt,DaqEvent &ev) :
                                                  ChipAPV(buf,copy_buf,opt,ev)
                                                {}

  private:

    /// Copy constructor is disabled.
                        ChipAPVRICH             (const ChipAPVRICH &e);
  
  //============================================================================
  // Operators
  //============================================================================

  private:

    /// Assignment operator is disabled.
    ChipAPVRICH       &operator =               (const ChipAPVRICH &e);

  //============================================================================
  // Static methods
  //============================================================================

  public:
    
    static
    std::pair<int,int>  ch2xy                   (uint8 channel,uint8 card_virtual,bool old_apv_board=false);

  //============================================================================
  // Methods
  //============================================================================

  public:
    
    /// Print chip data
    void                Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

    /// \return Chip's name. This is ChipAPV
    std::string         GetName                 (void) const {return "ChipAPVRICH";}

    /// Decode data and \b add new digits to \b digits_list. Call Scan() method if needed.
    void                Decode                  (const Maps &maps,Digits &digits_list,DaqOption &opt);
    
    /// Scan data, fill \c all_data attribute, detect errors.
    void                Scan                    (DaqOption &opt);

    /// \return all data found by Scan()  method.    
    std::vector<FullData> GetAllData            (void) const {return all_data;}

  //============================================================================
  // Attributes
  //============================================================================

  private:
};

////////////////////////////////////////////////////////////////////////////////

inline std::ostream &operator << (std::ostream &o,const ChipAPVRICH &e)
{
  e.Print(o);
  return o;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace CS

#endif // CompassSoft_ChipAPVRICH__include
