#ifndef CompassSoft_ChipHotGeSiCA__include
#define CompassSoft_ChipHotGeSiCA__include

#include "ChipF1.h"

namespace CS {

////////////////////////////////////////////////////////////////////////////////

/*! \brief This is ChipHotGeSiCA class.

    \author Alexander Zvyagin
*/
class ChipHotGeSiCA: public Chip
{
  //============================================================================
  // Types, constants
  //============================================================================

  public:

    class DataID
    {
      public:
        DataID(const Chip::DataID &d) {u.data_id=d;}
        DataID(uint16 srcID,uint8 geo_id,uint8 port,uint8 chan)
            {u.s.src_id=srcID; u.s.geo_id=geo_id; u.s.port=port; u.s.chan=chan; u.s.none=0;}
        operator Chip::DataID (void) const {return u.data_id;}
        operator std::string (void) const;
        union
        {
            Chip::DataID data_id;
            struct
            {
                Chip::DataID none:24,chan:8,port:8,geo_id:8,src_id:16;
            } s;
        } u;
    };

    class Digit: public ChipF1::Digit
    {
      public:
        virtual        ~Digit                   (void) {}

                        Digit                   (const DataID &data_id,const DetID &id,int32 chan,int16 channel_pos,int32 ampl,double t_unit)
                                                : ChipF1::Digit(Chip::DataID(data_id),id,chan,channel_pos,ampl,t_unit) {}

        virtual void    Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;
    };

    class Map : public Chip::Map
    {
      public:
                        Map                     (const ObjectXML &o);

        void            Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

        void            AddToMaps               (Maps &maps,DaqOption &option) const;

        /// The cost of one time bin
        double          GetTimeUnit             (void) const {return time_unit;}

        /// The cost of one time bin
        void            SetTimeUnit             (double t) {time_unit=t;}

        int16           GetWireP                (void) const {return wire_position;}

      public:

        /// geographical ID
        uint16          geoID;
        
        uint16          port;

        /// position for splitted wire: 0-normal
        int16           wire_position;

        double          time_unit;
    };

  //============================================================================
  // Constructors and destructor
  //============================================================================

  public:
    
                       ~ChipHotGeSiCA           (void) {}
    
    /*! \brief Base constructor.
        \param buf is non-NULL pointer to buffer where first bytes have format of
               Chip::Header.
    */
                        ChipHotGeSiCA           (void const * const buf,bool copy_buf,DaqOption &opt,DaqEvent &ev) :
                                                  Chip(buf,copy_buf,opt,ev)
                                                {}

  private:

    /// You can not use copy constructor.
                        ChipHotGeSiCA           (const ChipHotGeSiCA &e);
  
  //============================================================================
  // Operators
  //============================================================================

  private:

    /// You can not use assignment operator.
    ChipHotGeSiCA      &operator =              (const ChipHotGeSiCA &e);

  //============================================================================
  // Methods
  //============================================================================

  public:
  
    /// Print chip information
    void                Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

    /// \return Chip's name. This is "ChipHotGeSiCA".
    std::string         GetName                 (void) const {return "ChipHotGeSiCA";}

    void                Scan                    (DaqOption &opt);

    /// Decode data and \b add new digits to \b digits_list. Call Scan() method if needed.
    void                Decode                  (const Maps &maps,Digits &digits_list,DaqOption &opt);

  //============================================================================
  // Attributes
  //============================================================================

  private:
    
    std::vector< std::pair<DataID,uint16> > all_data;
};


////////////////////////////////////////////////////////////////////////////////

inline std::ostream &operator << (std::ostream &o,const ChipHotGeSiCA &e)
{
  e.Print(o);
  return o;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace CS

#endif // CompassSoft_ChipHotGeSiCA__include
