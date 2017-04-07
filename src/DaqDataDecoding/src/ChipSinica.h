#ifndef CompassSoft_ChipSinica__include
#define CompassSoft_ChipSinica__include

#include "Chip.h"

namespace CS {

class TriggerTimeConfig;

////////////////////////////////////////////////////////////////////////////////

/*! \brief This is ChipSinica class.

    This is Academia Scinica chip.

    Some documentation was extracted from "Compass-Note 2000-8, Compass-Note 2001-8"
    http://hpfr02.physik.uni-freiburg.de/projects/compass/electronics/notes/dataformat-2000-8/format.html
    
    \todo Not all modes are implemented. Please send your request to the maintainer.

    \author Yakov Kulinich
*/
class ChipSinica: public Chip
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
        DataID(uint16 a,uint16 b,uint16 c,uint16 d) {u.s.src_id=a; u.s.mode=b; u.s.port=c;u.s.channel=d;}
        operator Chip::DataID (void) const {return u.data_id;}
        union
        {
            Chip::DataID data_id;
            struct
            {
	      Chip::DataID channel:16,port:16,mode:16,src_id:16;  // trigger is fine time of optical trigger to dc5
            } s;
        } u;
    };

    /*
     *  Need this class because we dont just pass time back but also trigger time and temperature
     *  This object gets stored during ScanDbg and then its info added to Digit in Decode
     */
    class Hit
    {
    public:
      Hit(const Hit &d) {hit_time = d.hit_time; trigger_time = d.trigger_time; temperature =d.temperature;}
      Hit(uint16 a,uint16 b,uint16 c) {hit_time=a; trigger_time=b; temperature=c;}
      
      uint16 GetHitTime (void) const {return hit_time;}
      void   SetHitTime (uint16 h) {hit_time = h;};

      uint16 GetTriggerTime (void) const {return trigger_time;}
      void   SetTriggerTime (uint16 t) {trigger_time = t;};      

      double GetTemperature (void) const {return temperature;}
      void   SetTemperature (double t) {temperature = t;};      
 
    private:
      uint16 hit_time, trigger_time;
      double temperature;
    };


    class Digit: public Chip::Digit
    {
      public:
        virtual        ~Digit                   (void) {}

	Digit                   (const DataID &data_id,const DetID &id,int32 wir,int16 wire_pos, double t_unit, int32 overolling)
                                                : Chip::Digit(data_id,id), wire(wir), wire_position(wire_pos),
	                                          hittime(0), triggertime(0),   
	                                          overolling(overolling),  
	                                          temperature(0),
	                                          time_unit(t_unit), time_decoded(-1e6), time_reference(0) {}

        virtual void    Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

        virtual
        const char*     GetNtupleFormat         (void) const {return "wire:wirepos:thit:unit:overolling:tns";}
        virtual
        std::vector<float> GetNtupleData        (void) const;

	/// \set wire number
        void            SetWire                 (int32 w) {wire=w;}

        /// \return Detector's wire number.
        int32           GetWire                 (void) const {return wire;}

        int32           GetWirePos              (void) const {return wire_position;}

        /// Synonym for GetTime().
        int32           GetHitTime              (void) const {return hittime;}
        
        /// Set hit time value
        void            SetHitTime              (int32 h) {hittime=h;}
      
	/// Synonym for GetTime().
        int32           GetTriggerTime          (void) const {return triggertime;}
        /// Set trigger time value
        void            SetTriggerTime          (int32 h) {triggertime=h;}
      
        /// The cost of one time bin
        double          GetTimeUnit             (void) const {return time_unit;}

        /// The cost of one time bin
        void            SetTimeUnit             (double t) {time_unit=t;}
        
	/// Overolling Value
	int32           GetOverolling           (void) const {return overolling;}

	/// 
	void            SetOverolling           (double o) {overolling = o;}

	/// Temperature
	double           GetTemperature         (void) const {return temperature;}

	/// 
	void            SetTemperature          (double t) {temperature = t;}

        double          GetTimeDecoded          (void) const {return time_decoded;}
        void            SetTimeDecoded          (double t) {time_decoded=t;}

	float           GetTimeReference        (void) const {return time_reference;}
        void            SetTimeReference        (float t) {time_reference=t;}
                
      private:

        int32           wire;               
        int16           wire_position;
        int32           hittime;
	int32           triggertime;
	int32           overolling;
	double          temperature;
	double          time_unit;              ///< unit of time measurement
        double          time_decoded;           ///< time
        float           time_reference;
    };

    /*! \brief abstract class to represent properties of both debug-header mode,
         debug-trailer, and debug-data modes.

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
        virtual Data&   operator =              (uint32 d)   = 0;
        
        virtual         operator uint32         (void) const = 0;

        /// \return name of the mode
        virtual std::string   GetName           (void) const = 0;
        
	/// \return data
        virtual uint16   GetData                (void) const = 0;

	/// \return channel number 
        virtual uint8    GetChannel             (void) const = 0;

	/// \is this header?
        virtual bool     IsHeader               (void) const = 0;

	/// \is this trailer?
        virtual bool     IsTrailer              (void) const = 0;

        /// Print properties.
        virtual void     Print                  (std::ostream &o=std::cout,const std::string &prefix="") const = 0;

        /// Print-operator
        friend std::ostream    &operator <<     (std::ostream &o,const Data &e) {e.Print(o); return o;}
        friend class ChipSinica;
    };

    /*! \brief Header debug mode 
    */
    class HeaderDbg: public Data
    {
      public:
      
                        HeaderDbg               (void) {}
                        HeaderDbg               (uint32 d) {*this=d;}
      
	/// ----------- PURE VIRTUALS
        HeaderDbg&      operator =              (uint32 d) {header.all=d; return *this;}

                        operator uint32         (void) const {return  header.all;}

        /// \return it will return name "H-Dbg"
        std::string     GetName                 (void) const {return "H-Dbg";}

        /// \return data number (11 bits)
        uint16          GetData                 (void) const {return header.s.data;}

	/// \return FEM ID (5 bits)
        uint8           GetFEMID                (void) const {return header.s.port;}

	/// \return channel number 
        uint8           GetChannel              (void) const {return -1;}; // dont have channel in header

        /// \return \b true if it is indead a header
        bool            IsHeader                (void) const {return (header.s.zero == 0
								      && header.s.idle == 5);}  // header
        /// \return \b true if it is indead a trailer
        bool            IsTrailer               (void) const {return 0;}  // isnt trailer

	/// ----------- END OF PURE VIRTUAL

	/// \return Port   (5 bits)
        uint8           GetPort                 (void) const {return header.s.port;}
	
	/// \return Temp (4bits)
        uint8           GetTemperature          (void) const {return header.s.temp;}

	/// \return if TBOBit set
        bool            IsTBOBitSet             (void) const {return header.s.tbo;}

	/// \return Idle (4bits)
        uint8           GetIdle                 (void) const {return header.s.idle;}
        
        /// \return event number (8bits)
        uint8           GetEventNumber          (void) const {return header.s.event;}

        /// Print properties.
        void            Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

        /// Print-operator
        friend std::ostream    &operator <<     (std::ostream &o,const HeaderDbg &e) {e.Print(o); return o;}

    private:

        union
        {
          struct
          {
            uint32      temp    : 4,            ///< Temperature 
                        idle    : 4,            ///< always equals 5...
                        port    : 5,            ///< port (femid for dc5)
                        data    : 11,           ///< trigger course time
                        event   : 6,            ///< event number
                        tbo     : 1,            ///< trigger buffer overflow
                        zero    : 1;            ///< this bit = 0
          } s;
          uint32        all;
        } header;
        friend class ChipSinica;
    };
    
    /*! \brief Trailer debug mode 
    */
    class TrailerDbg: public Data
    {
      public:
      
                        TrailerDbg               (void) {}
                        TrailerDbg               (uint32 d) {*this=d;}
      
	/// ----------- PURE VIRTUALS
        TrailerDbg&      operator =              (uint32 d) {trailer.all=d; return *this;}

                        operator uint32          (void) const {return  trailer.all;}

        /// \return it will return name "H-Dbg"
        std::string     GetName                 (void) const {return "T-Dbg";}

        /// \return data number (11 bits)
        uint16          GetData                 (void) const {return trailer.s.data;}

	/// \return channel number 
        uint8           GetChannel              (void) const {return -1;}; // dont have channel in trailer

	/// \return \b true if it is indead a header
        bool            IsHeader                (void) const { return 0;}  // isnt header

        /// \return \b true if it is indead a trailer
        bool            IsTrailer               (void) const {return (trailer.s.zero == 0
								      && trailer.s.idle == 10);}  // trailer
	/// ----------- END OF PURE VIRTUAL

	/// \return Temp (4bits)
        uint8           GetTemperature          (void) const {return trailer.s.temp;}

	/// \return Idle bits (4 bits)
        uint8           GetIdle                 (void) const {return trailer.s.idle;}

        /// Print properties.
        void            Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

        /// Print-operator
        friend std::ostream    &operator <<     (std::ostream &o,const TrailerDbg &e) {e.Print(o); return o;}

    private:

        union
        {
          struct
          {
            uint32      temp    : 4,            ///< Temperature
                        idle    : 4,            ///< always equals 5...
	                data    : 16,           ///< trigger fine time
	                idle2   : 7,            ///< another idle
                        zero    : 1;            ///< this bit = 0
          } s;
          uint32        all;
        } trailer;
        friend class ChipSinica;
    };

    /*! \brief Data debug mode 
    */
    class DataDbg: public Data
    {
      public:

                        DataDbg                 (void) {}
                        DataDbg                 (uint32 d) {*this=d;}

	/// ----------- PURE VIRTUALS
        /// operator to assign an integer value
        DataDbg&        operator =              (uint32 d) {data.all=d; return *this;}

                        operator uint32         (void) const {return  data.all;}

        /// \return it will return name "D-Dbg"
        std::string     GetName                 (void) const {return "D-Dbg";}

        /// \return mode (4bits)
        uint8           GetMode                 (void) const {return data.s.mode;}

        /// \get data - fine time (16 bits)
        uint16          GetData                 (void) const {return data.s.data;}

        /// \return channel number (4 bits)
        uint8           GetChannel              (void) const {return data.s.channel;}

	/// \return \b true if it is indead a header
        bool            IsHeader                (void) const {return data.s.un == 0;} 

	/// \return \b true if it is indead a header
        bool            IsTrailer               (void) const {return data.s.un == 0;} 

	/// ----------- END OF PURE VIRTUAL

        /// Print properties.
        void            Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

	/// \return FEM ID (6 bits)
        uint8           GetIdle                (void) const {return data.s.idle;}

        /// Print-operator
        friend std::ostream    &operator <<     (std::ostream &o,const DataDbg &e) {e.Print(o); return o;}

        union
        {
          struct
          {
            uint32      mode      : 4,          ///< Mode
                        channel   : 4,          ///< FEM Channel id
                        data      :16,          ///< data (hit time)
                        idle      : 6,          ///< dont know
                        zero      : 1,          ///< this bit = 0
                        un        : 1;          ///< this bit = 1
          } s;
          uint32        all;
        } data;
        friend class ChipSinica;
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

        int32&          GetX                    (void)       {return wireF;}
        int32&          GetY                    (void)       {return wireL;}
        int32&          GetZ                    (void)       {return wireS;}
        
        int32           GetX                    (void) const {return wireF;}
        int32           GetY                    (void) const {return wireL;}
        int32           GetZ                    (void) const {return wireS;}
        
        int8            GetMode                 (void) const {return mode;}
	int8            GetGeoID                (void) const {return geo_id;}
	int8            GetPort                 (void) const {return port;}
	int8            GetFEMID                (void) const {return port;}   // femid on dc5 is port
        uint8           GetChannel              (void) const {return channel;}
        uint8           GetHit                  (void) const {return hit;}
        
        int16           GetWireP                (void) const {return wire_position;}
        int16           GetWireX                (void) const {return x;}
        int16           GetWireY                (void) const {return y;}

        /// The cost of one time bin
        double          GetTimeUnit             (void) const {return time_unit;}

        /// The cost of one time bin
        void            SetTimeUnit             (double t) {time_unit=t;}
        
	/// Overolling Value
	int32           GetOverolling           (void) const {return overolling;}

	/// 
	void            SetOverolling           (double o) {overolling = o;}

        char            GetType                 (void) const {return type_ul;}

      private:

        int8            mode;

	int16           geo_id;

        int16           port;
        
        uint16          channel;
        
        uint16          hit;
        
        int16           x,y;

        char            type_ul;

        /// position for splitted wire: 0-normal, 1/-1  up/down
        int16           wire_position;
        
	int32           overolling;
	double          time_unit;
	float           time_reference;         // approximate time measurement point
    };

  //============================================================================
  // Constructors and destructor
  //============================================================================

  public:
    
    /// Destructor
                       ~ChipSinica                  (void) {}
    
    /*! \brief Base constructor.
        \sa DaqEvent::DaqEvent(void const * const buf,bool copy_buf)
    */
                        ChipSinica                  (void const * const buf,bool copy_buf,DaqOption &opt,DaqEvent &ev);

  private:

     /// Copy constructor
                        ChipSinica                  (const ChipSinica &e);
  
  //============================================================================
  // Operators
  //============================================================================

  private:

    /// Assignment operator
    ChipSinica             &operator =              (const ChipSinica &e);

  //============================================================================
  // Static methods
  //============================================================================

  public:
  
    /// \return the time difference between the time point and trigger time, taking into account the overolling
    /// \ChipSinica trigger_time is optical trigger time arrival on front end clock and is thus an int in digits. No averaging/adjusting done so no double like other chips
    /// \if later another way to calculate trigger time is used, change this.
    static double       TimeDifference          (int time,int trigger_time,int time_overolling,double time_ref, int cut=30000);

  private:

  //============================================================================
  // Methods
  //============================================================================

  public:
  
    /// Print equipment
    void                Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

    /// Clear the chip status (do \b not modify the real data buffer).
    void                Clear                   (void) {data_all.clear();}

    /// \return Chip name. This is "ChipSinica"
    std::string         GetName                 (void) const {return "ChipSinica";}

    /// \return \b true if the data are in \c debug mode, otherwise the data are in sparsified mode.
    // had !() before
    bool                IsDebugMode             (void) const {return GetSLink().GetFormat()&(1<<1);}

    /// Decode data and \b add new digits to \b digits_list.
    void                Decode                  (const Maps &maps,Digits &digits_list,DaqOption &opt);

  private:

    /// Decode data and \b add new digits to \b digits_list.
    void                DecodeDebugData         (const Maps &maps,Digits &digits_list,DaqOption &opt);

  protected:

    /*! \brief Scan data, fill \c data_all attribute, detect errors.
    
        The code comes from Wolfgang Kastaun. Some modifications were done by
        Yakov Kulinich.
        \author Wolfgang Kastaun
        \author Yakov Kulinich
    */
    void                Scan                    (DaqOption &opt);

  private:

    /*! \brief Scan data, fill \c data_all attribute, detect errors.
    
        The code comes from Wolfgang Kastaun. Some modifications were done by
        Yakov Kulinich.
        \author Wolfgang Kastaun
        \author Yakov Kulinich
    */
    void                ScanDebugMode           (DaqOption &opt);

    /*! \brief Scan data, fill \c data_all attribute, detect errors.
    
        \author Yakov Kulinich
    */
    
    void                AddData                 (uint16 hittime,uint16 mode,uint16 channel,uint8 port, uint16 triggertime, uint16 temperature);

  //============================================================================
  // Attributes
  //============================================================================

  private:

    std::vector< std::pair<DataID,Hit> >     data_all;         ///< Data for decoding

  public:

    // Maximum data package size from FEM
    static const float FEM_DATA_SIZE_MAX;

    // Maximum allowed temperature before alarm tripped
    static const float FEM_TEMPERATURE_MAX;
};

////////////////////////////////////////////////////////////////////////////////

inline std::ostream &operator << (std::ostream &o,const ChipSinica &e)
{
  e.Print(o);
  return o;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace CS

#endif // CompassSoft_ChipSinica__include
