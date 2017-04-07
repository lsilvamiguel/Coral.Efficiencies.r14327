#ifndef CompassSoft_ChipGassiplex__include
#define CompassSoft_ChipGassiplex__include

#include <iostream>
#include "ChipCol.h"

namespace CS {

////////////////////////////////////////////////////////////////////////////////

/*! \brief COMPASS Gassiplex chip

    This chip is used in RICH detector readout.

    Some documentation was extracted from COMPASS note 2001-5
    "All you wanted to know about decoding the RICH-1 but were afraid to ask".

    \author Alexander Zvyagin
*/
class ChipGassiplex: public ChipCol
{
  //============================================================================
  // Types, constants
  //============================================================================

  public:
  
    enum            {NTables=2};            // Number of decoiding tables

    /*! \brief Gassiplex chip digit (for RICH detector)
    
    */
    class Digit: public Chip::Digit
    {
      public:
 
        /// Base constructor
                        Digit                   (const DataID &data_id,const DetID &id,
                                                 uint8 _cathode,uint8 _x,uint8 _y,uint16 _ampl,uint16 _table=0)
                                                : Chip::Digit(data_id,id),
                                                  cathode(_cathode),
                                                  ampl(_ampl),
                                                  decoding_table(_table)
                                                { xy[0]=_x; xy[1]=_y; }

        /// Print the digit.
        void            Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

        /// \return Digit format for ROOT-ntuple creation.
        const char*     GetNtupleFormat         (void) const {return "cathode:x:y:ampl";}
        
        /// \return vector of values that can be used for filling an Ntuple.
        std::vector<float>   GetNtupleData      (void) const
                                                {
                                                  std::vector<float> v;
                                                  v.push_back(GetCathode());
                                                  v.push_back(GetX());
                                                  v.push_back(GetY());
                                                  v.push_back(GetAmplitude());
                                                  return v;
                                                }
        
        /// \return Cathde [0-15]
        uint8           GetCathode              (void) const {return cathode;}
        
        /// \return Pixel x-coordinate [0-71]
        uint16          GetX                    (void) const {return xy[0];}

        /// \return Pixel y-coordinate [0-71]
        uint16          GetY                    (void) const {return xy[1];}
        
        /// \return ADC reading (10 bits)
        uint16          GetAmplitude            (void) const {return ampl;}
        
        uint16          GetDecodingTable        (void) const {return decoding_table;}

      private:

        uint8           cathode;                ///< [0-15]
        uint16          xy[2];                  ///< one photocathode has 72x72 pixels
        uint16          ampl;                   ///< 10 bits
        uint16          decoding_table;         ///< decoding table to be used.
    };

    /*! \brief Data line
    */
    class Data
    {
      public:
      
        /// Constructor
                        Data                    (void);

        /// Constructor
                        Data                    (uint8 msb,uint8 chamber,uint8 bora,uint16 chan,uint16 ampl);

        /// Constructor
                        Data                    (uint32 data) {d.all=data;}
                        
                        operator uint32         (void) const {return d.all;}

        /// Print chip data
        void            Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

        /// \return true if the first bit is set to 1 (\sa COMPASS Note-8, page 12).
        bool            IsData                  (void) const {return d.all&(1<<31);}
                        
        /// \return The ADC reading
        uint16          GetAmplitude            (void) const {return d.line.ampl;}
        
        uint16          GetBoraChannel          (void) const {return d.line.chan;}

        /// \return Get channel ID
        uint32          GetChannelID            (void) const {return d.line2.chan;}
        
        uint8           GetBora                 (void) const {return d.line.bora;}
        
        uint8           GetChamber              (void) const {return d.line.chamber;}
        
        uint8           GetMSB                  (void) const {return d.line.msb;}

      private:

        union
        {
          struct
          {
            uint32      ampl    :10;

            /// Channel within BORA
            uint32      chan    : 10;

            /// BORA number 0-23
            uint32      bora    : 5;

            /// Chamber number 0-7
            uint32      chamber : 3;

            /*! The four MSB are for synchronization by the CATHCEs, it
                helps them to recognize a data word from a Header or Trailer words.
                The MSB is 1 while the following three bits are 0.
                Instead, a Header or Trailer word would have bit-31 equal to 0.
                Bits 302-28 are presently set to 0 but they don't necessarily
                need to have that value.
            */
            uint32      msb     : 4;
          } line;

          struct
          {
            uint32      ampl    :10;
            uint32      chan    :18;
            uint32      msb     : 4;
          } line2;

          uint32 all;
        } d;
    };

    class Map : public Chip::Map
    {
      public:
                        Map                     (const ObjectXML &o);
        void            AddToMaps               (Maps &maps,DaqOption &option) const;
        uint16          decoding_table;         ///< decoding table to be used.
    };
    
    class GeoID
    {
      public:
                        GeoID                   (uint8 chamber,uint8 bora) {v=(chamber<<5)+bora;}
                        GeoID                   (uint8 vv) : v(vv) {}

        bool            operator ==             (const GeoID &g) const {return v==g.v;}
        bool            operator !=             (const GeoID &g) const {return !*this==g;}
        bool            operator <              (const GeoID &g) const {return v<g.v;}
                        operator uint8          (void)           const {return v;}

        uint8           GetChamber              (void) const {return v>>5;}
        uint8           GetBora                 (void) const {return v&31;}
        uint8           GetCathode              (void) const {return 15-(GetChamber()*2 + (GetBora()/12?1:0));}
        bool            IsBoraUp                (void) const {return GetBora()/12;}

      private:
        uint8           v;
    };
    
    class PixelCalib
    {
      public:
        enum {ENTRIES=2048};
      public:
                        PixelCalib              (GeoID geo_id,uint16 chan_id,unsigned decoding_table=0);
      public:
        bool            operator <              (const PixelCalib &c) const {return chanID<c.chanID;}
      public:
        void            Print                   (const std::string &preffix="") const;
        uint16          GetChannelID            (void) const {return chanID;}
        uint16          GetCathode              (void) const {return cathode;}
        uint16          GetX                    (void) const {return x;}
        uint16          GetY                    (void) const {return y;}
        void            SetThreshold            (uint16 v) {threshold=v;}
        uint16          GetThreshold            (void) const {return threshold;}
        void            SetSum                  (uint32 v) {sum=v;}
        uint32          GetSum                  (void) const {return sum;}
        void            SetSumSqr               (uint32 v) {sum_sqr=v;}
        uint32          GetSumSqr               (void) const {return sum_sqr;}
        float           GetSigma                (void) const;
      private:
        uint16          chanID;
        uint8           cathode;
        uint8           x,y;
        uint16          threshold;
        uint32          sum;
        uint32          sum_sqr;
    };
    
    class CathodeCalib
    {
      public:
        enum DataLength {SHORT=13,LONG=1021};
        typedef ChipGassiplex::GeoID GeoID;
        /// The buffer must have length at least 12 words +  (1word=32bits)
      public:
                        CathodeCalib            (const uint32 *buf,uint32 &length,unsigned decoding_table=0);
        bool            operator <              (const CathodeCalib &c) const {return geoID<c.geoID;}
      public:
        void            Print                   (const std::string &prefix="") const;
        const std::vector<uint8>& GetVoltage         (void) const {return voltage;}
              std::vector<uint8>& GetVoltage         (void)       {return voltage;}
        const std::vector<uint8>& GetTemperature     (void) const {return temperature;}
              std::vector<uint8>& GetTemperature     (void)       {return temperature;}
        uint32          GetBurstNumber          (void) const {return burst_number;}
        uint32          GetEventNumber          (void) const {return event_number;}
        const GeoID&    GetGeoID                (void) const {return geoID;}
        const std::set<PixelCalib> & GetPixels  (void) const {return pixels;}
              std::set<PixelCalib> & GetPixels  (void)       {return pixels;}
      private:
        GeoID           geoID;              ///< Bora header 8-bits value
        std::vector<uint8>   voltage;
        std::vector<uint8>   temperature;
        std::set<PixelCalib> pixels;
        uint32          burst_number;
        uint32          event_number;
    };
    
  //============================================================================
  // Constructors and destructor
  //============================================================================

  public:
    
    /// Destructor
                       ~ChipGassiplex           (void) {}
    
    /*! \brief Base constructor.
        \param buf is non-NULL pointer to buffer where first bytes have format of
               Chip::Header.
    */
                        ChipGassiplex           (void const * const buf,bool copy_buf,DaqOption &opt,DaqEvent &ev);

    /*! \brief Create fake event.
    */
                        ChipGassiplex           (void);

  private:

    /// Copy constructor is disabled.
                        ChipGassiplex           (const ChipGassiplex &e);
  
  //============================================================================
  // Operators
  //============================================================================

  private:

    /// Assignment operator is disabled.
    ChipGassiplex      &operator =              (const ChipGassiplex &e);

  //============================================================================
  // Methods
  //============================================================================

  public:

    /// Print chip data
    void                Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

    /// \return Chip's name. This is ChipGasspilex
    std::string         GetName                 (void) const {return "ChipGassiplex";}

    /// Decode data and \b add new digits to \b digits_list.
    void                Decode                  (const Maps &maps,Digits &digits_list,DaqOption &opt);

    virtual void        DecodeCalibEvent        (const Maps &maps,DaqOption &opt,...) const;

  private:
  
    static bool         Init                    (void);

  public:
    
    static void         GetXY                   (unsigned decoding_table,const GeoID &geo_id,uint16 chan_id,uint8 &x,uint8 &y);
    
  //============================================================================
  // Data
  //============================================================================
  
  private:
  
    static bool                                         init;
    static std::map<uint32, std::pair<uint8,uint8> >    id_to_xy[NTables];
    static uint16                                       f[NTables][6][8];
};

////////////////////////////////////////////////////////////////////////////////

} // namespace CS

#endif //  CompassSoft_ChipGassiplex__include
