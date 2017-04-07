#ifndef CompassSoft_ChipAPV__include
#define CompassSoft_ChipAPV__include

#include <iostream>
#include "Chip.h"
#include "PixelMMDecoding.h"

namespace CS {

////////////////////////////////////////////////////////////////////////////////

/*! \brief COMPASS APV chip

    The format description is available here:
    http://axp01.e18.physik.tu-muenchen.de/~konorov/el_design/sg_adc/sg_adc.html
    
    Some documentation was extracted from the above page.

    \author Alexander Zvyagin
*/
class ChipAPV: public Chip
{
  //============================================================================
  // Types, constants
  //============================================================================

  public:

    enum {APV_READ_CHANNELS=128};

    class DataID
    {
      public:
        DataID(const Chip::DataID &d) {u.data_id=d;}
        DataID(uint16 srcID,uint16 adc_id,uint16 chip_id,uint16 chan) {u.s.src_id=srcID; u.s.adc_id=adc_id; u.s.chip_id=chip_id;u.s.chan=chan;}
        operator Chip::DataID (void) const {return u.data_id;}
        union
        {
            Chip::DataID data_id;
            struct
            {
                Chip::DataID chan:16,chip_id:16,adc_id:16,src_id:16;
            } s;
        } u;
    };

    /*! \brief APV chip digit
    
    */
    class Digit: public Chip::Digit
    {
      public:
 
        /// Destructor
        virtual        ~Digit                   (void) {}

        /// Base constructor
                        Digit                   (const DataID &data_id,const DetID &id,bool sparsifed,
                                                 uint16 chip, uint16 chip_chan,
                                                 uint16 det_chan,
                                                 uint32 addr1,uint32 addr2,uint32 addr3,
                                                 uint32 ampl1,uint32 ampl2,uint32 ampl3,
                                                 uint32 _time_tag,
                                                 int8 chan_position=0,
						 int8 adc_type=0,
                                                 float _amplfactor=1.);

        /// Print the digit.
        virtual
        void            Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

        /// \return Digit format for ROOT-ntuple creation.
        virtual
        const char*     GetNtupleFormat         (void) const;
        
        /// \return vector of values that can be used for filling an Ntuple.
        virtual
        std::vector<float>   GetNtupleData      (void) const;
        
        /// \return Address (already encoded).
        const uint32*   GetAddress              (void) const {return address;}
        
        /// \return Channel number inside the chip (before remapping).
        uint16          GetChipChannel          (void) const {return chip_channel;}

        void            SetChipChannel          (uint16 ch)  {chip_channel=ch;}

        /// \return Detector's channel number.
        uint16          GetChannel              (void) const {return det_channel;}

        /// \return Chip number.    
        uint16          GetChip                 (void) const {return chip;}
        
        void            SetChip                 (uint16 n) {chip=n;}

        /// \return Amplitude
        const uint32*   GetAmplitude            (void) const {return amplitude;}
        
        void            SetAmplitude            (uint32 a1, uint32 a2, uint32 a3) {amplitude[0]=a1; amplitude[1]=a2; amplitude[2]=a3;}

        void            SetAddress              (uint32 a1, uint32 a2, uint32 a3) {address[0]=a1; address[1]=a2; address[2]=a3;}
        
        /// \return Common mode correction amplitudes (Common Noise)
        const uint16*   GetCoNo                 (void) const {return cono;}

        /// Set the common mode correction amplitudes (Common Noise)
        void            SetCoNo                 (uint16 a1,uint16 a2,uint16 a3) {cono[0]=a1;cono[1]=a2;cono[2]=a3;}

        bool            IsSingleFrame           (void) const {return is_single_frame;}
        
        bool            IsSparsifed             (void) const {return is_sparsifed;}
        void            SetSparsifed            (bool b) {is_sparsifed=b;}
        
        /*! \brief Time tage
        
            Time tag is a counter which contains an arrivel time of the trigger to the ADC whith
            a precision of 1 usecond. The time tag counter is reseted by the TCS_RESET signal
            in the beginning of each spill.
        */
        uint32          GetTimeTag              (void) const {return time_tag;}

        void            SetTimeTag              (uint32 t) {time_tag=t;}
        
        void            SetGoodDigit            (uint8 n) {good_digit=n;}
        
        /*! \brief Was the digit good or not.
        
            0 - digit is not good (in the map file it was specified to accept
                all digits, even a bad ones)
            1 - digit is good  (but it can be a bad one in a next event!)
            2 - digit is good, from the mapping file we accept only good digits.
        */
        uint8           GetGoodDigit            (void) const {return good_digit;}
        
        ///! Wire position: -1/+1 for split strips, 0 otherwise (default)
        int8            GetChanPos                (void) const {return chan_position;}

	// ADC type: 1 for RICH type ADC (info only from map file, not in data!), 0 otherwise
	int8            GetADCType               (void) const {return adc_type;}

        // Factor to be applied on amplitudes to compensate 50Ohm caps instead of 100Ohm
        float           GetAmplFactor           (void) const {return amplfactor;}


        void            SetStripConnNb         (int16 n) { strip_conn_number = n;}
        int16           GetStripConnNb         (void) const { return strip_conn_number;}

      private:

        bool            is_single_frame;
        bool            is_sparsifed;
        uint16          chip;
        uint16          chip_channel;
        uint16          det_channel;
        ///! position for split wire: 0-normal, 1/-1  up/down
        int8            chan_position;
	// ADC type: 1 for RICH type ADC, 0 for standard GEM/Si ADC
	int8            adc_type;
        float           amplfactor;  // amplitude multiplication factor (for APV cards with 50Ohm caps)
        uint32          address[3];             ///< Address
        uint32          amplitude[3];           ///< Amplitude
        uint32          time_tag;               ///< Time tag
        uint16          cono[3];                ///< Common mode correction (Common Noise)
        uint8           good_digit;
        int16           strip_conn_number;      ///< Strip connector number (used to identify APV chip in calibrations)
    };

    ////////////////////////////////////////////////////////////////////////////

    /*! \brief APV chip digit for pixel readout

    */
    class DigitPixel: public ChipAPV::Digit
    {
      public:

        /// Destructor
        virtual        ~DigitPixel              (void) {}

        /// Base constructor
                        DigitPixel              (const DataID &data_id,const DetID &id,bool sparsifed,
                                                 uint16 chip, uint16 chip_chan,
                                                 uint16 det_chan,
                                                 uint32 addr1,uint32 addr2,uint32 addr3,
                                                 uint32 ampl1,uint32 ampl2,uint32 ampl3,
                                                 uint32 _time_tag,
                                                 int8 det_orientation=0,
						 int8 adc_type=0);

        /// Print the digit.
        void            Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

        /// \return Digit format for ROOT-ntuple creation.
        const char*     GetNtupleFormat         (void) const;

        /// \return vector of values that can be used for filling an Ntuple.
        std::vector<float>   GetNtupleData      (void) const;

        void            SetDetOrientation       (int x) {det_orientation=x;}
        void            SetMapVersion           (int x) {map_version=x;}
        void            SetPixelXY              (int16 x,int16 y) {pixel_x=x; pixel_y=y;}

        int8            GetDetOrientation       (void) const {return det_orientation;}
        int8            GetMapVersion           (void) const {return map_version;}
        int16           GetPixelX               (void) const {return pixel_x;}
        int16           GetPixelY               (void) const {return pixel_y;}

      private:

        int8            det_orientation;
        int8            map_version;
        int16           pixel_x;
        int16           pixel_y;
    };

    ////////////////////////////////////////////////////////////////////////////

    /*! \brief APV chip digit for pixel readout specific to Pixel MM

    */
    class DigitPixelMM: public ChipAPV::Digit
    {
      public:

        /// Destructor
        virtual        ~DigitPixelMM            (void) {}

        /// Base constructor
                        DigitPixelMM            (const DataID &data_id,const DetID &id,bool sparsifed,
                                                 uint16 chip, uint16 chip_chan,
                                                 uint16 det_chan,
                                                 uint32 addr1,uint32 addr2,uint32 addr3,
                                                 uint32 ampl1,uint32 ampl2,uint32 ampl3,
                                                 uint32 _time_tag,
                                                 int8 det_orientation=0,
						 int8 adc_type=0,
                                                 int8 conn=-1,
                                                 float _amplfactor=1.,
                                                 int _pixelmm_version=1);

        /// Print the digit.
        void            Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

        /// \return Digit format for ROOT-ntuple creation.
        const char*     GetNtupleFormat         (void) const;

        /// \return vector of values that can be used for filling an Ntuple.
        std::vector<float>   GetNtupleData      (void) const;

        void            SetDetOrientation       (int x) {det_orientation=x;}
        void            SetMapVersion           (int x) {map_version=x;}
        void            SetConnNb               (int8 conn) {conn_number=conn;}
        void            SetPixelXYNb            (float x,float y,int16 nb) {pixelmm_x=x; pixelmm_y=y;pixelmm_nb=nb;}
        void            SetLargePixel           (bool largepix) {pixelmm_largepix=largepix;}

        int8            GetDetOrientation       (void) const {return det_orientation;}
        int8            GetMapVersion           (void) const {return map_version;}
        int8            GetConnNb               (void) const {return conn_number;}
        float           GetPixelX               (void) const {return pixelmm_x;}
        float           GetPixelY               (void) const {return pixelmm_y;}
        int16           GetPixelNb              (void) const {return pixelmm_nb;}
        bool            GetLargePixel           (void) const {return pixelmm_largepix;}
        int8            GetPixelMMversion       (void) const {return pixelmm_version;}

      private:

        int8            det_orientation;
        int8            map_version;
        int8            conn_number;  // 1 to 20 max
        bool            pixelmm_largepix;  // true if large pixel
        int16           pixelmm_nb;
        float           pixelmm_x;
        float           pixelmm_y;
        int8		pixelmm_version;
    };

    ////////////////////////////////////////////////////////////////////////////

    class CommonModeCorrection
    {
      public:
                    CommonModeCorrection        (uint32 dd=0) {d.i=dd;}
        uint16      operator[]                  (int i) const {if(i==0)return d.s.CMframe1;
                                                               if(i==1)return d.s.CMframe2;
                                                               if(i==2)return d.s.CMframe3;
                                                               throw Exception("CommonModeCorrection()[]: t=%d",i);}

        union
        {
          struct { uint32 CMframe1:10, CMframe2: 10, CMframe3: 10, zero:1, one:1; } s;
          uint32 i;
        } d;
    };

    ////////////////////////////////////////////////////////////////////////////
  
    /*! \brief ADC header
    */
    class ADCHeader
    {
      public:

        /// Print ADC header
        void            Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

        /// \return data size in 32-bits words, including this header.
        uint32          GetSize                 (void) const {return size;}
        
        /// \return ADC chip ID number.
        uint32          GetID                   (void) const {return id;}
        
        uint8           GetFormat               (void) const {return format_new;}
        
        /// \return format
        bool            IsStandardSparse        (void) const {return GetFormat()&64;}
        
        /// \return \b true if it is single-frame redout, othervise it is tripple-frame readout.
        bool            IsSingleFrame           (void) const {return !(GetFormat()&1);}

        int             GetSparseVer            (void) const {return bool(GetFormat()&2)*(2-IsStandardSparse());}

        /// Pedestal subtraction
        bool            IsPedSubtr              (void) const {return GetFormat()&4;}

        /// Test for common mode correction
        bool            IsComModCor             (void) const {return GetFormat()&8;}
        
        bool            IsReodered              (void) const {return GetFormat()&16;}

        bool            IsDataSkipped           (void) const {return GetFormat()&32;}
        
        /// \return time tag
        uint32          GetTimeTag              (void) const {return time_tag;}

      private:

        /// Block size counts number of 32 bit words including current word of header.
        uint32          size       :16;

        /// ADC ID
        uint32          id         : 8;

        /*! Data format

         Right now only 2 bits of data format are used, the rest is reserved.
         xxxxxxx0 - single frame readout 
         xxxxxxx1 - triple frame readout 
         xxxxxx0x - Latch ALL , no zero suppression 
         xxxxxx1x - Sparse, zero suppressed format
         xxxxx0xx - no pedestals substraction
         xxxxx1xx - pedestals substraction
         xxxx0xxx - no common mode correction aplied
         xxxx1xxx - common mode correction aplied
         xxx0xxxx - reodering off
         xxx1xxxx - reodering on
         xx0xxxxx - all data are readout
         xx1xxxxx - it was too much data, only headers present
         n0xxxxxx - ADC version V1.0
         n1xxxxxx - ADC version V2.0
        */
        uint32          format_old : 8;
        
        /*! \brief Time tag

        Time tag is a counter which counts an arrivel time of the trigger to the ADC with 
        a precision of 1 usecond. The time tag counter is reseted by theTCS_RESET signal 
        in the beginning of each spill. 
        */
        uint32          time_tag   :24;

        uint32          format_new : 8;
    };

    ////////////////////////////////////////////////////////////////////////////

    /*! \brief APV header

        The bits in the APV header are deffined as '1' if the amplitude of the signal bigger or
        equal then 512 and '0' if it is smaller. 
    */
    class APVHeader
    {
      public:
                        APVHeader               (uint32 const *dd) {u.d[0]=dd[0];u.d[1]=dd[1];}
      public:

        void            Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;
        
        /// \return length (in 32 bit words) of the APV data block (including header).
        uint32          GetBlockSize            (void) const {return u.s.data_len;}
        
        /// \return length (in 32 bit words) of the APV header
        uint32          GetSize                 (void) const {return 2;}
        
        /// \return chip number
        uint32          GetChipID               (void) const {return u.s.chip_id;}
        
        CommonModeCorrection const & GetCoNo    (void) const {return cono;}

        union
        {
          struct
          {
            uint32      data_len      :16;      ///< length (in 32 bit words) of the APV data block (including header)
            uint32      event_number  :12;      ///< event number from ADC
            uint32      chip_id       : 4;      ///< chip identifier
            uint32      frame1_good   : 1;      ///< Error bit
            uint32      frame1_gaddr  : 8;      ///< Gray address
            uint32      frame2_good   : 1;      ///< Error bit
            uint32      frame2_gaddr  : 8;      ///< Gray address
            uint32      frame3_good   : 1;      ///< Error bit
            uint32      frame3_gaddr  : 8;      ///< Gray address
            uint32      reserved0     : 2;      ///< 
            uint32      as            : 1;      ///<  APV sync bit
            uint32      gs            : 1;      ///<  global sync bit
            uint32      reserved1     : 1;      ///< 
          } s;

          uint32 d[2];
          
        } u;
        
        CommonModeCorrection    cono;
    };

    ////////////////////////////////////////////////////////////////////////////
    
    /*! \brief Data line structure.
    */
    class Data
    {
      public:

                        Data                    (uint32 data,int sparse) : sparse_ver(sparse) {d.any=data;}

        void            Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

        int             IsSparse                (void) const {return sparse_ver;}
        uint8           GetChannelID            (void) const {return IsSparse() ? d.s.channel_id : d.l.channel_id;}
        uint16          GetFrame3               (void) const {return IsSparse() ? d.s.frame3<<1           : d.l.frame3;  }
        uint16          GetFrame2               (void) const {return IsSparse() ? d.s.frame2<<sparse_ver  : d.l.frame2;  }
        uint16          GetFrame1               (void) const {return IsSparse() ? d.s.frame1<<sparse_ver  : d.l.frame1;  }

        uint16          GetFrameRaw1            (void) const {return IsSparse() ? d.s.frame1 : d.l.frame1;  }
        uint16          GetFrameRaw2            (void) const {return IsSparse() ? d.s.frame2 : d.l.frame2;  }
        uint16          GetFrameRaw3            (void) const {return IsSparse() ? d.s.frame3 : d.l.frame3;  }

        uint16          GetFrameNew1            (void) const;
        uint16          GetFrameNew2            (void) const;
        uint16          GetFrameNew3            (void) const;

      private:
      
        union
        {
          /// LatchAll
          struct { uint32 frame1:10, frame2:10, frame3:10, channel_id:2; } l;
          /// Sparse
          struct { uint32 frame1: 8, frame2: 8, frame3: 9, channel_id:7; } s;
          uint32 any;
        } d;
        int  sparse_ver;        /// can be 0,1,2
    };
    
    ////////////////////////////////////////////////////////////////////////////
    
    /// ChipAPV map
    class Map : public Chip::Map
    {
      public:

        /// Constructor
                        Map                     (const ObjectXML &o);
                        
        void            Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

        virtual void            AddToMaps               (Maps &maps,DaqOption &option) const;

      public:

        uint16          adc_id;
        uint16          chip_id;
        int             det_orientation;
        int             connNb;
        int             map_version;

        ///! Wire position for split wires
	int             wire_position;

        /// Amplitude factor for APV cards with 50Ohm caps on data line
        float           amplfact;

        bool            only_good_digits;
        
        bool            pixel_map;
        bool            pixelmm_map;	// for pixel MM
        bool            pixelmm2s_map;	// for pixel MM, simplified 2nd version
        bool            pixelmm3_map;	// for pixel MM, complete finale version, ELVIA detectors
        bool            rich_adc;

        ///! Do we always have common mode noise correction in data?
        bool            always_cmc;
    };

    ////////////////////////////////////////////////////////////////////////////

    /*! \brief  Data with headers
    */
    class FullData
    {
      public:

                        FullData                (const ADCHeader& _adc,const APVHeader& _apv,const Data& _data,uint16 _channel) :
                                                  adc(_adc), apv(_apv), data(_data), channel(_channel) {}

        const ADCHeader& GetADCHeader           (void) const {return adc;}
        const APVHeader& GetAPVHeader           (void) const {return apv;}
        const Data&      GetData                (void) const {return data;}
        uint16           GetChannel             (void) const {return data.IsSparse() ? data.GetChannelID():channel;}
        void             Print                  (std::ostream &o=std::cout,const std::string &prefix="") const;

      private:

        ADCHeader       adc;
        APVHeader       apv;
        Data            data;
        uint16          channel;
    };
    
  //============================================================================
  // Constructors and destructor
  //============================================================================

  public:
    
    /// Destructor
                       ~ChipAPV                 (void) {}
    
    /*! \brief Base constructor.
        \param buf is non-NULL pointer to buffer where first bytes have format of
               Chip::Header.
    */
                        ChipAPV                 (void const * const buf,bool copy_buf,DaqOption &opt,DaqEvent &ev) :
                                                  Chip(buf,copy_buf,opt,ev)
                                                {}

  private:

    /// Copy constructor is disabled.
                        ChipAPV                 (const ChipAPV &e);
  
  //============================================================================
  // Operators
  //============================================================================

  private:

    /// Assignment operator is disabled.
    ChipAPV            &operator =              (const ChipAPV &e);

  //============================================================================
  // Methods
  //============================================================================

  public:
    
    /// Print chip data
    virtual void        Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

    /// \return Chip's name. This is ChipAPV
    virtual std::string GetName                 (void) const {return "ChipAPV";}

    /// Decode data and \b add new digits to \b digits_list. Call Scan() method if needed.
    virtual void        Decode                  (const Maps &maps,Digits &digits_list,DaqOption &opt);
    
    /// Scan data, fill \c all_data attribute, detect errors.
    virtual void        Scan                    (DaqOption &opt);

    /// \return all data found by Scan()  method.    
    const std::vector<FullData> & GetAllData    (void) const {return all_data;}

  private:

    static bool         InitApvColumn           (void);

  //============================================================================
  // Attributes
  //============================================================================

  protected:
  
    std::vector<FullData>    all_data;
    
    static int          apvColumn[256];
    static bool         init;
};

////////////////////////////////////////////////////////////////////////////////

inline std::ostream &operator << (std::ostream &o,const ChipAPV &e)
{
  e.Print(o);
  return o;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace CS

// Private namespace for PixelGEM detector channel to X,Y mapping (also used in CsPixelGEMDetector.cc)
namespace pixgem {
  const int fXY_size=256;
  typedef struct {int x,y;} fXY_type;
  extern fXY_type fXY[fXY_size];
  bool initconndata();
  std::pair<int,int> ch2xy(int channelnew, int card_virtual);
  std::pair<int,int> detch2xy(int chan);
} // namespace pixgem

#endif // CompassSoft_ChipAPV__include
