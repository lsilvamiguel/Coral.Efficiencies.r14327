#ifndef CompassSoft_ChipSADC__include
#define CompassSoft_ChipSADC__include

#include "Chip.h"

namespace CS {

////////////////////////////////////////////////////////////////////////////////

/*! \brief COMPASS Sampling ADC chip

    Some documentation was extracted from
    "Multi-channel Sampling Analogue-to-Digital Converter Module i-SADC108032
    for the Mainz Crystal Ball Detector".

    \author Alexander Zvyagin
*/
class ChipSADC : public Chip
{
  //============================================================================
  // Types
  //============================================================================

  public:

    struct HeaderADC
    {
        enum Mode   {UNKNOWN,LATCH,SPARSE,ESPARSE,COMPRESSED};

        uint32      event:12;       /// last 12 bits of the event number
        uint32      size:12;        /// block size in 32-bits word
        uint32      sparse:1;       /// 1-sparse, 0-latch all
        uint32      overflow:1;     /// overflow in multi-hit buffer
        uint32      chip:4;         /// chip ID
        uint32      key:2;          /// Distinguish between latch/sparse/extended sparse modes.
        
        void        Print           (const char *prefix="") const;
        
        Mode        GetMode         (void) const;
    };
    
    struct HeaderADC_compressed
    {
        uint32      event:8;        /// last 8 bits of the event number
        uint32      huff_id:4;      /// reserved
        uint32      size:12;        /// block size in 32-bits word
        uint32      sparse:1;       /// 1-sparse, 0-latch all
        uint32      overflow:1;     /// overflow in multi-hit buffer
        uint32      chip:4;         /// chip ID
        uint32      key:2;          /// Distinguish between latch/sparse/extended sparse modes.
      
        void        Print           (const char *prefix="") const;
    };

    // ADC header of 2008 year versiuon
    struct Data
    {
        uint32      sample0:10;
        uint32      sample1:10;
        uint32      sample2:10;
        uint32      id:2;           /// not used, must be '10'

        void        Print           (const char *prefix="") const;
    };
    
    struct Integral
    {
        uint32      integral:16;    /// integral value
        uint32      not_used_1:7;   /// not used, must be '0'
        uint32      suppression:1;  /// suppression flag
        uint32      not_used_2:2;   /// not used, must be '0'
        uint32      channel:4;      /// channel number
        uint32      not_used_3:2;   /// not used, must be '10'

        void        Print           (const char *prefix="") const;
    };

    struct ChannelInfo
    {
        uint32      sum:16;
        uint32      samples:9;
        uint32      not_used:1;
        uint32      channel:4;
        uint32      id:2;           /// Should be 01

        void        Print           (const char *prefix="") const;
    };

    // ---- Begin of data format 3 (year 2008) description --------------

    struct HeaderADC_v3
    {
        uint32      event:12;       /// last 12 bits of the event number
        uint32      size:12;        /// block size in 32-bits word
        uint32      error:1;        /// error bit
        uint32      adc_id:2;       /// ADC ID
        uint32      hl_port:3;      /// HL port
        uint32      key:2;          /// must be 11
        
        void        Print           (const char *prefix="") const;
    };
    
    struct DataHeader_v3
    {
        uint32      sum:16;         /// Sum (of what?!)
        uint32      samples:9;      /// number of samples
        uint32      zero:1;         /// Must be 0
        uint32      chan_id:4;      ///
        uint32      id:2;           /// Must be 01

        void        Print           (const char *prefix="") const;
    };

    struct Data_v3
    {
        uint32      data1:12;       /// data
        uint32      data2:12;       /// data
        uint32      id:8;           /// Must be 01000000 for the first data and 10000000 for others

        void        Print           (const char *prefix="") const;
    };

    // ---- End of data format 3 (year 2008) description --------------

    class DataID
    {
      public:
        DataID(const Chip::DataID &d) {u.data_id=d;}
        DataID(uint16 srcID,uint16 port,uint16 chip,uint16 chan) {u.s.src_id=srcID; u.s.port=port; u.s.chip=chip; u.s.chan=chan;u.s.none=0;}
        operator Chip::DataID (void) const {return u.data_id;}
        uint16          GetSourceID             (void) const {return u.s.src_id;}
        void            SetSourceID             (uint16 s)   {u.s.src_id=s;}
        uint16          GetPort                 (void) const {return u.s.port;}
        void            SetPort                 (uint16 p)   {u.s.port=p;}
        uint16          GetChip                 (void) const {return u.s.chip;}
        void            SetChip                 (uint16 s)   {u.s.chip=s;}
        uint16          GetChannel              (void) const {return u.s.chan;}
        void            SetChannel              (uint16 s)   {u.s.chan=s;}
        union
        {
            Chip::DataID data_id;
            struct
            {
                Chip::DataID none:24,chan:8,chip:8,port:8,src_id:16;
            } s;
        } u;
        void            Print                   (const std::string &prefix="") const;
    };

    class Digit: public Chip::Digit
    {
      public:

                        Digit                   (const DataID &data_id,const DetID &id)
                                                : Chip::Digit(data_id,id),overflow(false),suppression(false),x(50000),y(50000)
                                                   {}

        void            Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

        const char*     GetNtupleFormat         (void) const {return "x:y";}
        std::vector<float>   GetNtupleData      (void) const;

        const std::vector<uint16>& GetSamples   (void) const {return samples;}
              std::vector<uint16>& GetSamples   (void)       {return samples;}

        const std::vector<uint16>& GetIntegrals (void) const {return integrals;}
              std::vector<uint16>& GetIntegrals (void)       {return integrals;}

        bool            GetOverflow             (void) const {return overflow;}
        void            SetOverflow             (bool s)     {overflow=s;}

        bool            GetSuppression          (void) const {return suppression;}
        void            SetSuppression          (bool s)     {suppression=s;}

        uint16          GetChip                 (void) const {return static_cast<ChipSADC::DataID>(GetDataID()).GetChip();}
        void            SetChip                 (uint16 c)   {static_cast<ChipSADC::DataID>(GetDataID()).SetChip(c);}

        uint16          GetChannel              (void) const {return static_cast<ChipSADC::DataID>(GetDataID()).GetChannel();}
        void            SetChannel              (uint16 c)   {static_cast<ChipSADC::DataID>(GetDataID()).SetChannel(c);}
        
        int16           GetX                    (void) const {return x;}
        void            SetX                    (int16 xx)   {x=xx;}
        
        int16           GetY                    (void) const {return y;}
        void            SetY                    (int16 yy)   {y=yy;}

      private:

        bool            overflow;
        bool            suppression;
        int16           x,y;                    // decoded x,y coordanates
        std::vector<uint16>  samples;
        std::vector<uint16>  integrals;
    };

    class DigitCalib: public Digit
    {
      public:
        struct Box
        {
            Box(int16 _x1,int16 _y1,int16 _x2,int16 _y2) :
                x1(_x1), y1(_y1), x2(_x2), y2(_y2) {}
            void            Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;
            int16 x1,y1,x2,y2;
        };

                        DigitCalib              (const DataID &data_id,const DetID &id) : Digit(data_id,id) {}
        
        const
        std::vector<Box> &   GetBoxes           (void) const {return boxes;}
        std::vector<Box> &   GetBoxes           (void)       {return boxes;}

//        int16           GetX2                   (void) const {return x2;}
//        int16           GetY2                   (void) const {return y2;}

        const std::string& GetRefDet            (void) const {return det_ref;}
        void            SetRefDet               (const std::string &s) {det_ref=s;}

//        AddBox          SetX2Y2                 (int16 xx,int16 yy)   {x2=xx; y2=yy;}

        void            Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

      private:

        std::string     det_ref;
        std::vector<Box>     boxes;
    };

    class Map : public Chip::Map
    {
      public:
                        Map                     (const ObjectXML &o);
                        
        void            Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

        void            AddToMaps               (Maps &maps,DaqOption &option) const;

        bool            IsCalib                 (void) const {return det_ref!="";}

      public:

        uint16          port;
        uint16          chip;
        
        uint16          channel;
        
        int16           x,y;
        
        // Mapping needed for the laster monitoring system
        
        std::string     det_ref;
        std::string     type;

        std::vector<DigitCalib::Box> boxes;
        std::string     data_fmt;
    };


  //============================================================================
  // Constructors and destructor
  //============================================================================

  public:
    
    /// Destructor
                       ~ChipSADC                (void) {Clear();}
    
    /*! \brief Base constructor.
        \sa DaqEvent::DaqEvent(void const * const buf,bool copy_buf)
    */
                        ChipSADC                (void const * const buf,bool copy_buf,DaqOption &opt,DaqEvent &ev);
  
  //============================================================================
  // Methods
  //============================================================================

  private:

    /// \return Chip name.
    std::string         GetName                 (void) const {return "ChipSADC";}

    /// Clear the chip status (do \b not modify the real data buffer).
    void                Clear                   (void);

  public:

    void                Scan                    (DaqOption &opt);
    void                ScanDataVer3            (DaqOption &opt);

    void                Decode                  (const Maps &maps,Digits &digits_list,DaqOption &opt);

  protected:

    void                ScanADCBlockFormat2     (const HeaderADC &adc,const uint32 *start,DaqOption &opt);
    void                ScanADCBlockFormat3     (const HeaderADC &adc,const uint32 *start,DaqOption &opt);

  //============================================================================
  // Data
  //============================================================================
  
  private:
  
    std::list<Digit*>   pre_digits;
};


} // namespace CS

#endif
