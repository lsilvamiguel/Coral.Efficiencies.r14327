#ifndef CompassSoft_ChipCol__include
#define CompassSoft_ChipCol__include

#include "Chip.h"

namespace CS {

////////////////////////////////////////////////////////////////////////////////

/*! \brief Chips collection for ChipADC, ChipGassiplex and ChipScalar

    Some documentation was extracted from "Compass-Note 2000-8"
    http://hpfr02.physik.uni-freiburg.de/projects/compass/electronics/notes/dataformat-2000-8/format.html

    \author Alexander Zvyagin
    \todo Scalar data
*/
class ChipCol: public Chip
{
  //============================================================================
  // Types, constants
  //============================================================================

  public:

    enum {ErrorMarker=0x7ff};

    /*! \brief Data header/trailer for ADC,Scaler and RICH chips

    */
    class HeaderTrailer
    {
      public:
      
        /// Base constructor.
                        HeaderTrailer           (uint32 _event_number,uint32 _geoID,bool _is_trailer=true, bool _is_data=false)
                                                {
                                                  d.s.event      = _event_number;
                                                  d.s.geoID      = _geoID;
                                                  d.s.is_trailer = _is_trailer;
                                                  d.s.is_data    = _is_data;
                                                }

        /// Copy constructor.
                        HeaderTrailer           (uint32 v) {d.all=v;}

        /// Assignment operator
        HeaderTrailer   operator =              (uint32 v) {d.all=v; return *this;}

        /// Transfer HeaderTrailer to uint32 value.
                        operator uint32         (void) const {return d.all;}

        /// \return true if it is data, not a trailer/header (\sa COMPASS-Note 8, page 12)
        bool            IsData                  (void) const {return d.s.is_data;}

        /*! \return true if it is a header
        */
        bool            IsHeader                (void) const {return !d.s.is_data && !d.s.is_trailer;}

        /*! \return true if it is a trailer
        */
        bool            IsTrailer               (void) const {return !d.s.is_data && d.s.is_trailer;}

        /*! \return Geographical ID address
    
            \exception Exception if it is data, not a trailer/header
        */
        uint32          GetGeoID                (void) const;

        /*! \return Event number

            \exception Exception if it is data, not a trailer/header
        */
        uint32          GetEventNumber          (void) const;

      private:

        union
        {
          struct
          {
            uint32      event      :20,         ///< Event number
                        geoID      :10,         ///< Geometrical ID
                        is_trailer : 1,         ///< 0 for header and 1 for trailer
                        is_data    : 1;         ///< Zero for header/trailer and 1 for data
          } s;                                  ///< Structure with trailer/header properties
          uint32 all;                           ///< 32-bits word of trailer/header
        } d;                                    ///< Union for 32-bits word and structure.
    };

  //============================================================================
  // Constructors and destructor
  //============================================================================

  public:
    
                       ~ChipCol                 (void) {}

    /*! \brief Base constructor.
        \param buf is non-NULL pointer to buffer where first bytes have format of
               Chip::Header.
    */
                        ChipCol                 (void const * const buf,bool copy_buf,DaqOption &opt,DaqEvent &ev) :
                                                  Chip(buf,copy_buf,opt,ev)
                                                {}

  private:

    /// You can not use copy constructor.
                        ChipCol                 (const ChipCol &e);
  
  //============================================================================
  // Operators
  //============================================================================

  private:

    /// You can not use assignment operator.
    ChipCol            &operator =              (const ChipCol &e);

  //============================================================================
  // Methods
  //============================================================================

  public:
  
    /// Print chip information
    virtual void        Print                   (std::ostream &o=std::cout,const std::string &prefix="") const {Chip::Print(o,prefix);}

    /// \return Chip's name. This is "ChipCol".
    virtual std::string GetName                 (void) const = 0;

    /// Decode data and \b add new digits to \b digits_list. Call Scan() method if needed.
    virtual void        Decode                  (const Maps &maps,Digits &digits_list,DaqOption &opt) = 0;
    
    /*! \brief Add data header/trailer

        \exception Exception if you do not have rights to modify the chip data
        \exception Exception if the trailer/header is not allowed (two headers, two trailers, data after trailer, header after data)
    */
    void                AddHeaderTrailer        (const HeaderTrailer &ht);

    /*! \brief Add data 32-bits word

        \exception Exception if you do not have rights to modify the chip data
        \exception Exception if the trailer was not inserted before.
    */
    void                AddData                 (uint32 data);

    /// Scan data, fill \c all_data attribute, detect errors.
    virtual void        Scan                    (DaqOption &opt);

  private:

    // Check and report about bad event number.
    void                CheckEventNumber          (const HeaderTrailer &ht,DaqOption &opt);

  protected:

    std::vector< std::pair<HeaderTrailer,uint32> > all_data;
    static std::map<int32,bool> geoid_with_bad_events_counter;
};

////////////////////////////////////////////////////////////////////////////////

inline
std::ostream &operator << (std::ostream &o,const ChipCol &e)
{
  e.Print(o);
  return o;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace CS

#endif // CompassSoft_ChipCol__include
