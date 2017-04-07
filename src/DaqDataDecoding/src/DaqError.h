#ifndef CS__DaqError___include
#define CS__DaqError___include

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>

#include "Exception.h"

namespace CS {

class DaqErrorType;

////////////////////////////////////////////////////////////////////////////////

/*! \brief Errors reporting for DAQ

    Rule of the thumb: All objects of DaqError class should be reported
    to an appriate errors list. Never throw an exception with a DaqError object.
    
    \author Alexander Zvyagin
*/
class DaqError
{
  ////////////////////////////////////////////////////////////////////////////
  // Constructors, destructor
  ////////////////////////////////////////////////////////////////////////////

  public:

    /*! The errors severity levels */
    enum SeverityLevel
    {
      OK,                       ///< no warning
      MINOR_PROBLEM,            ///< minor problem
      WARNING,                  ///< warning
      SEVERE_PROBLEM            ///< severe problem warning
    };

    /*! Actions to be taken if a problem occured */
    class Action
    {
      public:
        /// An action type
        enum Type
        {
          NOTHING,                  ///< do nothing, use the data 
          DISCARD_EVENT_ON_GEOID,   ///< discard the rest of the event on this port/GeoID
	  DISCARD_EVENT_ON_PORT,    ///< discard the rest of the event on this port
          DISCARD_EVENT_ON_SRCID,   ///< discard the rest of the event on this SourceID
          DISCARD_EVENT             ///< discard whole event
        };
    };

    /*! Error types */
    enum Type
    {
      UNKNOWN_TYPE,         ///< Error with unknown type
      
      EXCEPTION,            ///< Exception

      TCS_FIFO,             ///< FIFO full bit in TCS errors field was set
      TCS_ECC,              ///< ECC checksum error bit set
      TCS_SYNC,             ///< synchroniation error bit set
      TCS_UNDEF,            ///< other (undefined) bit in this field set      
      SLINK_ERRFLAG,        ///< Error Flag in SLink header was set
      SLINK_ERRNR,          ///< #errorwords field in SLink header was not zero
      SLINK_WRONG_EVENT_NUMBER,

      TDC_SERR1,            ///< Trigger buffer overflow
      TDC_SERR2,            ///< CMC FIFO full
      TDC_SERR3,            ///< transmission error in bit 23..0

      TDC_WRGEOID,          ///< Wrong Geographical ID
      TDC_ERR1,             ///< Data or setup word received instead of header
      TDC_ERR2,             ///< Timeout if no words received
      TDC_ERR3,             ///< Header or trailer with smaller event number received
      TDC_ERR4,             ///< No CMC is connected
      TDC_ERR5,             ///< Timeout if no trailer or too many data words received
      TDC_ERR6,             ///< FIFO full, port is off until end of burst
      TDC_ERR7,             ///< Header or trailer with larger event number received
      TDC_UNDEF8,           ///< Error# 8 :undefined
      TDC_UNDEF9,           ///< Error# 9 :undefined
      TDC_ERR10,            ///< HOTLink transmission error
      TDC_UNDEF11,          ///< Error# 11:undefined
      TDC_ERR12,            ///< Trailer if event data was skipped (too high rate)
      TDC_UNDEF13,          ///< Error# 13:undefined
      TDC_UNDEF14,          ///< Error# 14:undefined
      TDC_UNDEF15,          ///< Error# 15:undefined

      TDC_WRPORT_T,         ///< Trailer Port number not equal to the one in the corresponding header
      TDC_TBO,              ///< Trigger buffer overflow bit in TDC Header word set
      TDC_WRPLL_H,          ///< PLL bits not set on TDC-CMC header
      TDC_ERR3U,            ///< Header or trailer with smaller event number occured,
                            ///< but without an preceeding error word
      TDC_ERR7U,            ///< same with larger event number
      TDC_WRERR,            ///< TDC_ERR3 / 7 error word found , but following event number
                            ///< was not smaller/larger
      TDC_ERR1U,            ///< Data word without header word received, 
                            ///< but without preceeding error word 
      TDC_WRPORT,           ///< Data Port number not equal to the one in the corresponding header
      TDC_WRPLL,            ///< PLL bits not set on TDC-CMC data
      TDC_WRTIME,           ///< Trigger time in two headers differs more than 1.

      TDC_TBO2,             ///< Trigger buffer overflow bit in TDC Header word set
      TDC_TEMPERATURE,      ///< Temperature alert
      TDC_ERR20,            ///< Maximum data package received without header
      TDC_ERR21,            ///< Maximum data package received without trailer
      TDC_ERR22,            ///< Maximum data package received without header or trailer

      ADC_DATA_WO_HEADER,
      ADC_W_GEOID_EW,
      ADC_NO_TRAILER,
      ADC_TRAILER_WO_HEADER,
      ADC_W_GEOID_TRAILER,
      ADC_WR_EV_HEADER,
      ADC_WR_EV_TRAILER,
      
      SADC_WRONG_EVENT,     ///< Wrong event number
      SADC_BAD_H_ADC,       ///< Bad ADC header
      SADC_H_ADC_ERROR,     ///< 'error' bit is set in ADC header
      SADC_H_ADC_BAD_SIZE,  ///< Bad block size in ADC header
      SADC_UNKNOWN_MODE,    ///< Unknown ADC mode
      SADC_BAD_DATA_HEADER, ///< Bad 'data header' word
      SADC_DATA_OUTSIDE,    ///< Data words outside data block
      SADC_MODULE_MISSING,  ///< Module is not connected
      SADC_BAD_DATA,        ///< Bad 'data' word
      SADC_BAD_INTEGRAL,    ///< Bad 'integral' word
      SADC_BAD_INTEGRALS,   ///< Bad 'integral' words are not consistent
      SADC_BAD_CHAN_N,      ///< Channel number is no increasing
      SADC_NO_DATA_LINE,    ///< Data line is not found.
      SADC_UNKNOWN_DATA,    ///< Unknown data line
      SADC_SHORT_SAMPLES,   ///< Too little channel data samples
      SADC_BAD_N_SAMPLES,   ///< Bad number of samples
  
      SLINKM_BAD_SIZE,
      
      APV_H_S_E,
      APV_ADC_OUT,
      APV_APV_OUT,
      APV_H_WR_EV,
      APV_LOCAL_TIMETAG_FAILED,
      APV_LAST_TRIGGER_TICKS_FAILED,
      APV_MISS_ADC,
      
      GASSIPLEX_BAD_MSB,
      GASSIPLEX_BIG_BORA,
      GASSIPLEX_BAD_CHAN,
      GASSIPLEX_BAD_GEOID,
      GASSIPLEX_BAD_CHANID,
      
      ChipHotGeSiCA_WRONG_EVENT1,
      ChipHotGeSiCA_WRONG_EVENT2,
      ChipHotGeSiCA_WRONG_EVENT3,
      ChipHotGeSiCA_WRONG_EVENT4,
      ChipHotGeSiCA_WRONG_PORT_PORT,
      ChipHotGeSiCA_WRONG_CHIP_PORT,
      ChipHotGeSiCA_BAD_EVENT,
      ChipHotGeSiCA_BAD_CHIP,
      ChipHotGeSiCA_UNEXPECTED_DT,
      ChipHotGeSiCA_WRONG_CHIP,
      ChipHotGeSiCA_PORT_TIMEOUT,
      ChipHotGeSiCA_EVENT_PORT_LOST,
      ChipHotGeSiCA_DATA_PORT_LOST,
      ChipHotGeSiCA_SKIP_SPILL, 
      ChipHotGeSiCA_WRONG_EVENT,     
      
      EVENT_CORRUPTED0,
      EVENT_CORRUPTED1,
      EVENT_CORRUPTED2,
      EVENT_SUBEVENT_HEADER_DIFF,
      EVENT_UNKNOWN,
      EVENT_BAD_CHIP_SIZE,
      EVENT_SRCID_TOO_BIG,
      
      EVENT_MISSING_SRCID,
      EVENT_UNKNOWN_SRCID,
      EVENT_WRONG_GEOID_ORDER,
      
      TT_WRONG_HITS_IN_MASTER_TIME,
      TT_BIG_SIGMA,
      WRONG_HITS_IN_TRIGGER_MASK,
      TT_MT_SHIFT_CHANGED,
      
      TCS_PHASE_NO_DIGITS,
      TCS_PHASE_MANY_DIGITS,
      TCS_PHASE_OUT_OF_RANGE,        ///< tcs phase is outside of +/- 200ns

      TIS_NO_DATA,                   ///< no data from FluxScalers for time in spill computation
      
      MAX_TYPE
    };

    /// Types
    enum Arg
    {
        SLINK,
        PORT,
        CHAMBER,
        CHIP,
        BORA,
        CHANNEL,
        CHANNEL_ID,
        SOURCE_ID,
        GEO_ID,
        COUNTER,
        EVENT_N,
        SL_EVENT_N,
        VALUE
    };
    
    typedef double ArgType;
    typedef std::multimap<Arg,ArgType> Args;

  ////////////////////////////////////////////////////////////////////////////
  // Static methods
  ////////////////////////////////////////////////////////////////////////////

  public:

    static const std::string& GetArgName    (Arg t);

  private:

  ////////////////////////////////////////////////////////////////////////////
  // Constructors, destructor
  ////////////////////////////////////////////////////////////////////////////

  public:
  
                       ~DaqError            (void) {delete e;}

                        DaqError            (const DaqError &d) : e(NULL) {*this=d;}
                       
                        DaqError            (const Exception &ee) : type(EXCEPTION) {e=new Exception(ee);}

                        DaqError            (Type t)
                                            : e(NULL), type(t) {}

                        DaqError            (Type t,const Args &_args)
                                            : e(NULL), type(t)
                                            {args=_args;}

                        DaqError            (Type t, Arg a1,ArgType v1)
                                            : e(NULL), type(t)
                                            {args.insert(std::pair<Arg,ArgType>(a1,v1));}

                        DaqError            (Type t, Arg a1,ArgType v1, Arg a2,ArgType v2)
                                            : e(NULL), type(t)
                                            {args.insert(std::pair<Arg,ArgType>(a1,v1));
                                             args.insert(std::pair<Arg,ArgType>(a2,v2));}

                        DaqError            (Type t, Arg a1,ArgType v1, Arg a2,ArgType v2, Arg a3,ArgType v3)
                                            : e(NULL), type(t)
                                            {args.insert(std::pair<Arg,ArgType>(a1,v1));
                                             args.insert(std::pair<Arg,ArgType>(a2,v2));
                                             args.insert(std::pair<Arg,ArgType>(a3,v3));}

                        DaqError            (Type t, Arg a1,ArgType v1, Arg a2,ArgType v2, Arg a3,ArgType v3, Arg a4,ArgType v4)
                                            : e(NULL), type(t)
                                            {args.insert(std::pair<Arg,ArgType>(a1,v1));
                                             args.insert(std::pair<Arg,ArgType>(a2,v2));
                                             args.insert(std::pair<Arg,ArgType>(a3,v3));
                                             args.insert(std::pair<Arg,ArgType>(a4,v4));}

                        DaqError            (Type t, Arg a1,ArgType v1, Arg a2,ArgType v2, Arg a3,ArgType v3,
                                                     Arg a4,ArgType v4, Arg a5,ArgType v5)
                                            : e(NULL), type(t)
                                            {args.insert(std::pair<Arg,ArgType>(a1,v1));
                                             args.insert(std::pair<Arg,ArgType>(a2,v2));
                                             args.insert(std::pair<Arg,ArgType>(a3,v3));
                                             args.insert(std::pair<Arg,ArgType>(a5,v5));}

                        DaqError            (Type t, Arg a1,ArgType v1, Arg a2,ArgType v2, Arg a3,ArgType v3,
                                                     Arg a4,ArgType v4, Arg a5,ArgType v5, Arg a6,ArgType v6)
                                            : e(NULL), type(t)
                                            {args.insert(std::pair<Arg,ArgType>(a1,v1));
                                             args.insert(std::pair<Arg,ArgType>(a2,v2));
                                             args.insert(std::pair<Arg,ArgType>(a3,v3));
                                             args.insert(std::pair<Arg,ArgType>(a5,v5));
                                             args.insert(std::pair<Arg,ArgType>(a6,v6));}

                        DaqError            (Type t, Arg a1,ArgType v1, Arg a2,ArgType v2, Arg a3,ArgType v3,
                                                     Arg a4,ArgType v4, Arg a5,ArgType v5, Arg a6,ArgType v6,
                                                     Arg a7,ArgType v7)
                                            : e(NULL), type(t)
                                            {args.insert(std::pair<Arg,ArgType>(a1,v1));
                                             args.insert(std::pair<Arg,ArgType>(a2,v2));
                                             args.insert(std::pair<Arg,ArgType>(a3,v3));
                                             args.insert(std::pair<Arg,ArgType>(a5,v5));
                                             args.insert(std::pair<Arg,ArgType>(a6,v6));
                                             args.insert(std::pair<Arg,ArgType>(a7,v7));}

  ////////////////////////////////////////////////////////////////////////////
  // Operators
  ////////////////////////////////////////////////////////////////////////////

  public:

    DaqError&           operator =          (const DaqError &d);

                        operator Exception  (void) const;

  ////////////////////////////////////////////////////////////////////////////
  // Methods
  ////////////////////////////////////////////////////////////////////////////

  public:


    const DaqErrorType& GetDaqErrorType     (void) const;
    Type                GetType             (void) const {return type;}
    
    const Args&         GetArgs             (void) const {return args;}
          Args&         GetArgs             (void)       {return args;}
    
    /*! \brief Find and return value of the argument 'arg'
        \arg arg argument name
        \arg number if there are two (or more) arguments of the same type, specify the number of argument you are looking for
        \return an argument value (throw an exception if not found) or \c not_found if it was not set.
    */
    ArgType             GetArg              (const Arg &arg,size_t number=0) const;

    /// Print the error
    void                Print               (std::ostream &o=std::cout,const std::string &prefix="") const;

  ////////////////////////////////////////////////////////////////////////////
  // Attributes, data
  ////////////////////////////////////////////////////////////////////////////

  private:

    Exception*          e;

    /// Error type
    Type                type;

    /// List of arguments
    Args                args;

  ////////////////////////////////////////////////////////////////////////////
  // Static attributes, data
  ////////////////////////////////////////////////////////////////////////////

  private:
    static std::map<Arg, std::string> argument_names;        ///< The explanation of argument names
    static bool         init;                                ///< The class initialization
};

////////////////////////////////////////////////////////////////////////////////

/*! \brief The DAQ error type
*/
class DaqErrorType
{
  ////////////////////////////////////////////////////////////////////////////
  // Types
  ////////////////////////////////////////////////////////////////////////////

  public:

  ////////////////////////////////////////////////////////////////////////////
  // Static methods
  ////////////////////////////////////////////////////////////////////////////

  public:

    static std::map<DaqError::Type,DaqErrorType>& GetDaqErrorTypes(void) {return error_types;}
    static DaqErrorType& GetDaqErrorType(DaqError::Type t);
    static DaqErrorType& GetDaqErrorType(const std::string &name);

  private:

    static bool         Init                    (void);

  ////////////////////////////////////////////////////////////////////////////
  // Constructors, destructor
  ////////////////////////////////////////////////////////////////////////////

  public:

    /// The destructor
                       ~DaqErrorType            (void) {}

    /// Main constructor
                        DaqErrorType            (DaqError::Type t,const std::string &description,
                                                 DaqError::SeverityLevel l,DaqError::Action::Type a);

  ////////////////////////////////////////////////////////////////////////////
  // Methods
  ////////////////////////////////////////////////////////////////////////////

  public:

    const DaqError::Type&          GetType          (void) const {return type;}
    const std::string&             GetName          (void) const {return name;}
    const DaqError::SeverityLevel& GetSeverityLevel (void) const {return level;}
    const DaqError::Action::Type&  GetAction        (void) const {return action;}

    void                           SetSeverityLevel (const DaqError::SeverityLevel& l) {level=l;}
    void                           SetAction        (DaqError::Action::Type a) {action=a;}

  ////////////////////////////////////////////////////////////////////////////
  // Attributes, data
  ////////////////////////////////////////////////////////////////////////////

  public:
  
    DaqError::Type          type;               ///< An error type
    std::string             name;               ///< An error name
    DaqError::SeverityLevel level;              ///< Importance level
    DaqError::Action::Type  action;             ///< Action to be taken if an error happens

  ////////////////////////////////////////////////////////////////////////////
  // Static attrubutes, data
  ////////////////////////////////////////////////////////////////////////////

  private:

    static std::map<DaqError::Type,DaqErrorType> error_types; ///< The container for all errors
    static bool             init;
};

////////////////////////////////////////////////////////////////////////////////

class DaqErrors
{
  public:

    /*! @brief Add a new error 
        @return error index. */
    size_t              Add                     (const DaqError &e);
    
    /*! @brief Remove an error number 'n'. Remove all errors if 'n<0' */
    void                Remove                  (int n=-1);

    /*! @brief Print errors for a given type. */
    void                Print                   (DaqError::SeverityLevel t) const;

    /*! @brief Print errors for a given action. */
    void                Print                   (DaqError::Action::Type t) const;

    void                Print                   (DaqError::Type t) const;

    const std::set<size_t> & Get                (DaqError::Type t) const;
    const std::set<size_t> & Get                (DaqError::SeverityLevel t) const;
    const std::set<size_t> & Get                (DaqError::Action::Type t) const;

  public:

    std::vector<DaqError>                                 errors;
    std::map<DaqError::Type,std::set<size_t> >            errors_type;
    std::map<DaqError::SeverityLevel,std::set<size_t> >   errors_level;
    std::map<DaqError::Action::Type,std::set<size_t> >    errors_action;
};

////////////////////////////////////////////////////////////////////////////////

inline
const DaqErrorType& DaqError::GetDaqErrorType(void) const
{
    return DaqErrorType::GetDaqErrorType(type);
}

////////////////////////////////////////////////////////////////////////////////

} // namespace CS

#endif // CS__DaqError___include
