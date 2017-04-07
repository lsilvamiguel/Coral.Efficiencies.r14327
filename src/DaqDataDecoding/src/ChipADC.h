#ifndef CompassSoft_ChipADC__include
#define CompassSoft_ChipADC__include

#include "ChipCol.h"

namespace CS {

////////////////////////////////////////////////////////////////////////////////

/*! \brief COMPASS ADC chip

    Some documentation was extracted from "Compass-Note 2000-8"
    http://hpfr02.physik.uni-freiburg.de/projects/compass/electronics/notes/dataformat-2000-8/format.html

    \author Alexander Zvyagin
*/
class ChipADC: public ChipCol
{
  //============================================================================
  // Types, constants
  //============================================================================

  public:

    class DataID
    {
      public:
        DataID(const Chip::DataID &d) {u.data_id=d;}
        DataID(uint16 srcID,uint16 geo_id,uint16 chan) {u.s.src_id=srcID; u.s.geo_id=geo_id; u.s.chan=chan;u.s.none=0;}
        operator Chip::DataID (void) const {return u.data_id;}
        union
        {
            Chip::DataID data_id;
            struct
            {
                Chip::DataID none:32,chan:8,geo_id:8,src_id:16;
            } s;
        } u;
    };

    class Digit: public Chip::Digit
    {
      public:

                        Digit                   (const DataID &data_id,const DetID &id,uint16 &chan,int32 ampl)
                                                : Chip::Digit(data_id,id), x(chan), y(0), amplitude(ampl) {}

                        Digit                   (const DataID &data_id,const DetID &id,int32 _x,int32 _y,int32 ampl)
                                                : Chip::Digit(data_id,id), x(_x), y(_y), amplitude(ampl) {}

        void            Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

        const char*     GetNtupleFormat         (void) const {return "x:y:ampl";}
        std::vector<float> GetNtupleData        (void) const {std::vector<float> v; v.push_back(GetX()); v.push_back(GetY()); v.push_back(GetAmplitude()); return v;}

        int32           GetAmplitude            (void) const {return amplitude;}
        void            SetAmplitude            (uint32 a) {amplitude=a;}
        
        int32           GetX                    (void) const {return x;}
        int32           GetY                    (void) const {return y;}
        
      private:

        int32           x,y;
        int32           amplitude;
    };

    class Data
    {
      public:

                        Data                    (const uint32 &d) {u.d=d;}

        int16           GetData                 (void) const {return u.s.ofl ? ((int16)u.s.data-4096) : (int16)u.s.data;}
        
                        operator uint32         (void) const {return u.d;}

        uint8           GetChannel              (void) const {return u.s.channel;}
        
        bool            IsData                  (void) const {return u.s.is_data;}
        
        void            Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;
        
      private:

        union
        {
          struct
          {
            uint32      data       :12,
                        otr        : 1,
                        ofl        : 1,
                        channel    : 6,
                        zero       :10,
                        ped_subst  : 1,
                        is_data    : 1;
          } s;
          uint32 d;
        } u;
    };

    class Map : public Chip::Map
    {
      public:
                        Map                     (const ObjectXML &o);
                        
        void            Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

        void            AddToMaps               (Maps &maps,DaqOption &option) const;

      public:

        /// geographical ID
        uint16          geoID;
        
        uint16          channel;
        
        uint16          x,y;
    };

  //============================================================================
  // Constructors and destructor
  //============================================================================

  public:
    
                       ~ChipADC                 (void) {}
    
    /*! \brief Base constructor.
        \param buf is non-NULL pointer to buffer where first bytes have format of
               Chip::Header.
    */
                        ChipADC                 (void const * const buf,bool copy_buf,DaqOption &opt,DaqEvent &ev) :
                                                  ChipCol(buf,copy_buf,opt,ev)
                                                {}

  private:

    /// You can not use copy constructor.
                        ChipADC                 (const ChipADC &e);
  
  //============================================================================
  // Operators
  //============================================================================

  private:

    /// You can not use assignment operator.
    ChipADC            &operator =              (const ChipADC &e);

  //============================================================================
  // Methods
  //============================================================================

  public:
  
    /// Print chip information
    void                Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

    /// \return Chip's name. This is "ChipADC".
    std::string         GetName                 (void) const {return "ChipADC";}

    /// Decode data and \b add new digits to \b digits_list. Call Scan() method if needed.
    void                Decode                  (const Maps &maps,Digits &digits_list,DaqOption &opt);

  //============================================================================
  // Attributes
  //============================================================================

  private:
};

////////////////////////////////////////////////////////////////////////////////

inline std::ostream &operator << (std::ostream &o,const ChipADC &e)
{
  e.Print(o);
  return o;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace CS

#endif // CompassSoft_ChipADC__include
