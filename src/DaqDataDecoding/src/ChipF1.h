#ifndef CompassSoft_ChipF1__include
#define CompassSoft_ChipF1__include

#include "Chip.h"

namespace CS {

class TriggerTimeConfig;

////////////////////////////////////////////////////////////////////////////////

/*! \brief This is ChipF1 class.

    This is F1 TDC chip.

    Some documentation was extracted from "Compass-Note 2000-8, Compass-Note 2001-8"
    http://hpfr02.physik.uni-freiburg.de/projects/compass/electronics/notes/dataformat-2000-8/format.html
    
    \todo Not all modes are implemented. Please send your request to the maintainer.

    \author Alexander Zvyagin
*/
class ChipF1: public Chip
{
  //============================================================================
  // Types, constants
  //============================================================================

  private:

    enum {ErrorMarker=0xfff};

  public:

    class DataID
    {
      public:
        DataID(const Chip::DataID &d) {u.data_id=d;}
        DataID(uint16 a,uint16 b,uint16 c,uint8 d,uint8 e) {u.s.src_id=a; u.s.mode=b; u.s.geoID_or_port=c;u.s.chip_chan=d; u.s.hit=e;}
        operator Chip::DataID (void) const {return u.data_id;}
        union
        {
            Chip::DataID data_id;
            struct
            {
                Chip::DataID hit:8,chip_chan:8,geoID_or_port:16,mode:16,src_id:16;
            } s;
        } u;
    };

    class Digit: public Chip::Digit
    {
      public:
        virtual        ~Digit                   (void) {}

                        Digit                   (const DataID &data_id,const DetID &id,int32 chan,int16 channel_pos,int32 ampl,double t_unit)
                                                : Chip::Digit(data_id,id), channel(chan), channel_position(channel_pos),
                                                  amplitude(ampl), time_unit(t_unit),time_decoded(-1e6), xy_mode(false), time_reference(0) {}

                        Digit                   (const DataID &data_id,const DetID &id,int16 x,int16 y,int16 channel_pos,int32 ampl,double t_unit)
                                                : Chip::Digit(data_id,id), channel(uint16(x)+(uint16(y)<<16)), channel_position(channel_pos),
                                                  amplitude(ampl), time_unit(t_unit),time_decoded(-1e6), xy_mode(true), time_reference(0) {}

        virtual void    Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

        virtual
        const char*     GetNtupleFormat         (void) const {return "chan:time:pos:unit:tns";}
        virtual
        std::vector<float> GetNtupleData        (void) const;

        /// \return Detector's channel number.
        int32           GetChannel              (void) const {return channel;}

        void            SetChannel              (int32 ch) {channel=ch;}

        /// \return Detector's channel number.
        int32           GetWire                 (void) const {return channel;}

        /// \return Detector's channel X coordinate
        int32           GetX                    (void) const {return channel&0xffff;}

        /// \return Detector's channel Y coordinate
        int32           GetY                    (void) const {return channel>>16;}

        int32           GetChannelPos           (void) const {return channel_position;}

        /// \return time associated with the channel
        int32           GetTime                 (void) const {return amplitude;}

        /// Synonym for GetTime().
        int32           GetAmplitude            (void) const {return amplitude;}
        
        /// Set amplitude (time value)
        void            SetAmplitude            (int32 a) {amplitude=a;}
        
        /// The cost of one time bin
        double          GetTimeUnit             (void) const {return time_unit;}

        /// The cost of one time bin
        void            SetTimeUnit             (double t) {time_unit=t;}
        
        double          GetTimeDecoded          (void) const {return time_decoded;}
        void            SetTimeDecoded          (double t) {time_decoded=t;}
        
        float           GetTimeReference        (void) const {return time_reference;}
        void            SetTimeReference        (float t) {time_reference=t;}
        
      private:

        int32           channel;
        int16           channel_position;
        int32           amplitude;
        double          time_unit;              ///< unit of time measurement
        double          time_decoded;           ///< time in ns with the respect to trigger time
        bool            xy_mode;
        float           time_reference;
    };

    class DigitRICHPMT: public Digit
    {
      public:
                        DigitRICHPMT            (const DataID &data_id,const DetID &id) :
                                                Digit(data_id,id,0,0,0,0),
                                                dcard_x(255),dcard_y(255),type_ul('?'),expert_decoding(false) {}

        void            Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

        const char*     GetNtupleFormat         (void) const {return "x:y:time";}
        std::vector<float> GetNtupleData        (void) const;

        uint16          GetDCardX               (void) const {return dcard_x;}
        uint16          GetDCardY               (void) const {return dcard_y;}

        char            GetType                 (void) const {return type_ul;}

        void            SetDCardX               (uint16 x) {dcard_x=x;}
        void            SetDCardY               (uint16 y) {dcard_y=y;}
        
        void            SetType                 (char ul) {type_ul=ul;}
        
        // The very last stage of the decoding.
        void            Finalize                (void);
        
        bool            IsExpertDecoding        (void) const {return expert_decoding;}
        void            SetExpertDecoding       (bool f) {expert_decoding=f;}

      public:

        
      
      private:
        
        uint16          dcard_x;
        uint16          dcard_y;
        char            type_ul;
        bool            expert_decoding;
    };

    /*! \brief abstract class to represent properties of both debug-header mode and
               debug-data mode.

        One data line is only 32 bits long! But it's structure is rather heavy.
        It uses not only bit-fields, but different bit-field formats for different modes.
        So it is usefull to create a C++ classes structure for this 32-bits word.
    */
    class Data
    {
      public:
        /// empty virtual default destructor to remove compiler warnings
        virtual         ~Data                   (void) {};

        /// operator to assign an integer value
        virtual Data&   operator =              (uint32 d) = 0;
        
        virtual         operator uint32         (void) const = 0;

        /// \return name of the mode
        virtual std::string   GetName           (void) const = 0;
        
        virtual uint16   GetData                (void) const = 0;

        /// PLL attribute (4 bits)
        virtual uint8    GetPLL                 (void) const = 0;

        /// \return port number (4 bits)
        virtual uint8    GetPort                (void) const = 0;

        /// \return channel number (3 bits)
        virtual uint8    GetChannel             (void) const = 0;

        /// \return chip ID (3 bits)
        virtual uint8    GetChip                (void) const = 0;
        
        virtual bool     IsHeader               (void) const = 0;

        /// Print properties.
        virtual void     Print                  (std::ostream &o=std::cout,const std::string &prefix="") const;

        /// Print-operator
        friend std::ostream    &operator <<     (std::ostream &o,const Data &e) {e.Print(o); return o;}
        friend class ChipF1;
    };

    /*! \brief Header debug mode 
    */
    class HeaderDbg: public Data
    {
      public:
      
                        HeaderDbg               (void) {}
                        HeaderDbg               (uint32 d) {*this=d;}
      
        HeaderDbg&      operator =              (uint32 d) {data.all=d; return *this;}

                        operator uint32         (void) const {return  data.all;}

        /// \return it will return name "H-Dbg"
        std::string     GetName                 (void) const {return "H-Dbg";}

        /// \return time (9bits)
        uint16          GetData                 (void) const {return data.s.data;}
        
        /// \return event number (6bits)
        uint8           GetEventNumber          (void) const {return data.s.event;}

        /// \return PLL (4bits)
        uint8           GetPLL                  (void) const {return data.s.pll;}

        /// \return port number (3 bits)
        uint8           GetPort                 (void) const {return data.s.port;}

        /// \return channel number (3 bits)
        uint8           GetChannel              (void) const {return data.s.channel;}

        /// \return chip ID (3 bits)
        uint8           GetChip                 (void) const {return data.s.chip;}
        
        bool            IsTBOBitSet             (void) const {return data.s.tbo;}

        /// \return \b true if it is indead a header
        bool            IsHeader                (void) const {return data.s.zero==0;}

        /// Print properties.
        void            Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

        /// Print-operator
        friend std::ostream    &operator <<     (std::ostream &o,const HeaderDbg &e) {e.Print(o); return o;}

    private:

        union
        {
          struct
          {
            uint32      pll     : 4,            ///< Pll locks for TDC CMC 
                        port    : 4,            ///< Hotlink port number on Catch
                        channel : 3,            ///< channel id
                        chip    : 3,            ///< chip id
                        xor_bit : 1,            ///< xor of the F1 setup
                        data    : 9,            ///< trigger time
                        event   : 6,            ///< event number
                        tbo     : 1,            ///< trigger buffer overflow
                        zero    : 1;            ///< this bit = 0
          } s;
          uint32        all;
        } data;
        friend class ChipF1;
    };

    /*! \brief Data debug mode 
    */
    class DataDbg: public Data
    {
      public:

                        DataDbg                 (void) {}
                        DataDbg                 (uint32 d) {*this=d;}

        /// operator to assign an integer value
        DataDbg&        operator =              (uint32 d) {data.all=d; return *this;}

                        operator uint32         (void) const {return  data.all;}

        /// \return it will return name "D-Dbg"
        std::string     GetName                 (void) const {return "D-Dbg";}

        /// The same as GetData16()
        uint16          GetData                 (void) const {return GetData16bit();}

        /// Get data in 16 bits word
        uint16          GetData16bit            (void) const {return data.s.data;}

        /// Get first 12 bits of data
        uint16          GetData12bit            (void) const {return data.s.data>>4;}

        /// Get last 4 bits of data
        uint8           GetData4bit             (void) const {return data.s.data & 15;}

        /// \return PLL (4bits)
        uint8           GetPLL                  (void) const {return data.s.pll;}

        /// \return port number (3 bits)
        uint8           GetPort                 (void) const {return data.s.port;}

        /// \return channel number (6 bits)
        uint8           GetChipChannel          (void) const {return data.s.chip_chan;}

        /// \return channel number (3 bits)
        uint8           GetChannel              (void) const {return data.s.chip_chan&7;}

        /// \return chip ID (3 bits)
        uint8           GetChip                 (void) const {return data.s.chip_chan>>3;}

        /// \return \b false if it is indead a data.
        bool            IsHeader                (void) const {return data.s.un==0;}

        /// Print properties.
        void            Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

        /// Print-operator
        friend std::ostream    &operator <<     (std::ostream &o,const DataDbg &e) {e.Print(o); return o;}

        union
        {
          struct
          {
            uint32      pll       : 4,          ///< Pll locks for TDC CMC
                        port      : 4,          ///< Hotlink port number on Catch
                        data      :16,          ///< data (hit time)
                        chip_chan : 6,          ///< chip/channel
                        zero      : 1,          ///< this bit = 0
                        un        : 1;          ///< this bit = 1
          } s;
          uint32        all;
        } data;
        friend class ChipF1;
    };

    ////////////////////////////////////////////////////////////////////////////

    class Map : public Chip::Map
    {
      public:
      
                        Map                     (const ObjectXML &o);

      public:

        /// Print the map definition
        void            Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

        void            AddToMaps               (Maps &maps,DaqOption &option) const;

        /*! @brief Read trigger time configuration. */
        void            ReadTTConfig            (TriggerTimeConfig &tt_conf) const;

        int32&          GetX                    (void)       {return wireF;}
        int32&          GetY                    (void)       {return wireL;}
        int32&          GetZ                    (void)       {return wireS;}
        
        int32           GetX                    (void) const {return wireF;}
        int32           GetY                    (void) const {return wireL;}
        int32           GetZ                    (void) const {return wireS;}
        
        int8            GetMode                 (void) const {return mode;}
        int8            GetGeoID                (void) const {return geo_id;}
        int8            GetPort                 (void) const {return port;}
        int8            GetChip                 (void) const {return chip;}
        uint8           GetChipChan             (void) const {return chan_in_chip;}
        uint8           GetHit                  (void) const {return hit;}
        
        int16           GetWireP                (void) const {return wire_position;}
        int16           GetWireX                (void) const {return x;}
        int16           GetWireY                (void) const {return y;}

        /// The cost of one time bin
        double          GetTimeUnit             (void) const {return time_unit;}

        /// The cost of one time bin
        void            SetTimeUnit             (double t) {time_unit=t;}
        
        char            GetType                 (void) const {return type_ul;}

      private:

        /// 'l' for "latch" and 't' for "tdc"
        int8            mode;

        int16           geo_id;

        int16           port;
        
        int16           chip;
        
        uint16          chan_in_chip;
        
        uint16          hit;
        
        int16           x,y;

        char            type_ul;
        
        bool            rich_pmt_expert;

        /// position for splitted wire: 0-normal, 1/-1  up/down
        int16           wire_position;
        
        double          time_unit;
        float           time_reference;         // approximate time measurement point
    };

  //============================================================================
  // Constructors and destructor
  //============================================================================

  public:
    
    /// Destructor
                       ~ChipF1                  (void) {}
    
    /*! \brief Base constructor.
        \sa DaqEvent::DaqEvent(void const * const buf,bool copy_buf)
    */
                        ChipF1                  (void const * const buf,bool copy_buf,DaqOption &opt,DaqEvent &ev);

  private:

     /// Copy constructor
                        ChipF1                  (const ChipF1 &e);
  
  //============================================================================
  // Operators
  //============================================================================

  private:

    /// Assignment operator
    ChipF1             &operator =              (const ChipF1 &e);

  //============================================================================
  // Static methods
  //============================================================================

  public:
  
    /// \return the time difference between the time point and trigger time, taking into account the overolling
    static double       TimeDifference          (int time,double trigger_time,int time_overolling,double time_ref,double cut=30000);

  private:

  //============================================================================
  // Methods
  //============================================================================

  public:
  
    /// Print equipment
    void                Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

    /// Clear the chip status (do \b not modify the real data buffer).
    void                Clear                   (void) {data_all.clear();}

    /// \return Chip name. This is "ChipF1"
    std::string         GetName                 (void) const {return "ChipF1";}

    /// \return \b true if the data are in \c debug mode, otherwise the data are in sparsified mode.
    bool                IsDebugMode             (void) const {return !(GetSLink().GetFormat()&(1<<1));}

    /// \return \b true if the data are in \c latch mode.
    bool                IsLatchMode             (void) const {return GetSLink().GetFormat()&(1<<3);}

    /// \return \b true if the data are in high resolution mode, otherwise the data are in normal resolution.
    bool                IsHighResolutionMode    (void) const {return GetSLink().GetFormat()&(1<<4);}

    /// \return \b true if the data are with leading and trailing edge.
    bool                IsLeadTrail             (void) const {return GetSLink().GetFormat()&(1<<5);}

    /// \return \b true if the data are in HOTL-CMC mode, otherwise the data are in TDC-CMC mode.
    bool                IsHOTLCMC               (void) const {return GetSLink().GetFormat()&(1<<5);}

    /// Decode data and \b add new digits to \b digits_list.
    void                Decode                  (const Maps &maps,Digits &digits_list,DaqOption &opt);

  private:

    /// Decode data and \b add new digits to \b digits_list.
    void                DecodeDebugData         (const Maps &maps,Digits &digits_list,DaqOption &opt);

    /// Decode data and \b add new digits to \b digits_list.
    void                DecodeSparsifiedData    (const Maps &maps,Digits &digits_list,DaqOption &opt);

  protected:

    /*! \brief Scan data, fill \c data_all attribute, detect errors.
    
        The code comes from Wolfgang Kastaun. Some modifications were done by
        Alexander Zvyagin.
        \author Wolfgang Kastaun
        \author Alexander Zvyagin
    */
    void                Scan                    (DaqOption &opt);

  private:

    /*! \brief Scan data, fill \c data_all attribute, detect errors.
    
        The code comes from Wolfgang Kastaun. Some modifications were done by
        Alexander Zvyagin.
        \author Wolfgang Kastaun
        \author Alexander Zvyagin
    */
    void                ScanDebugMode           (DaqOption &opt);

    /*! \brief Scan data, fill \c data_all attribute, detect errors.
    
        \author Alexander Zvyagin
    */
    void                ScanSparsifiedMode      (DaqOption &opt);
    
    void                AddData                 (uint16 time,uint16 mode,uint16 geoID_or_port,uint8 chip_chan);

  //============================================================================
  // Attributes
  //============================================================================

  private:

    std::vector< std::pair<DataID,uint16> >     data_all;         ///< Data for decoding
};

////////////////////////////////////////////////////////////////////////////////

inline std::ostream &operator << (std::ostream &o,const ChipF1 &e)
{
  e.Print(o);
  return o;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace CS

#endif // CompassSoft_ChipF1__include
