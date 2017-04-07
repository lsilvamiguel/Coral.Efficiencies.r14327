#ifndef CompassSoft_DaqEvent__include
#define CompassSoft_DaqEvent__include

#include <cstdio>
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include "config.h"
#include "DaqError.h"
#include "Chip.h"
#include "DaqOption.h"
#include "OnlineFilter.h"
#include "TriggerTime.h"

// Old versions of ROOT define this.
#ifdef Check
#undef Check
#endif

namespace CS {

class Event1Run;
typedef void (*DaqEventCallMe)(CS::Chip&);

////////////////////////////////////////////////////////////////////////////////

/*! \brief DaqEvent class was designed for reading and analysis of DAQ
           events in DATE (Daq Test Environment) format in COMPASS experiment.

    DaqEvent provides C++ interface to DAQ event. You can construct event
    either from memory buffer DaqEvent::DaqEvent(void const * const event_buffer)
    or from an input stream DaqEvent::DaqEvent(istream &f), print event to screen
    DaqEvent::Print(), get list of chips in that event with
    DaqEvent::ReadChips() and DaqEvent::GetChips() and of course you
    have access to different fields of a DATE event: event size
    DaqEvent::GetLength(), time DaqEvent::GetTime() and so on.

    \author Alexander Zvyagin
*/
class DaqEvent
{
  //============================================================================
  // Types, constants
  //============================================================================

  public:

    /*! \brief Class for throwing an exception
        DaqError of that type will be thrown if end of data stream will be reached.
    */
    class ExceptionEndOfStream : public Exception
    {
      public:
        ExceptionEndOfStream(void) : Exception("DaqEvent::DaqEvent(): end of stream") {}
    };
    
    enum EventType
    {
      START_OF_RUN         =  1,
      END_OF_RUN           =  2,
      START_OF_RUN_FILES   =  3,
      END_OF_RUN_FILES     =  4,
      START_OF_BURST       =  5,
      END_OF_BURST         =  6,
      PHYSICS_EVENT        =  7,
      CALIBRATION_EVENT    =  8,
      END_OF_LINK          =  9,
      EVENT_FORMAT_ERROR   = 10
    };

    /*! \brief Event header. C++ wrapper for DATE event header structure
    */
    class Header
    {
      public:

        // Minimal amount of bytes that need to be read by
        // DaqEvent::DaqEvent(istream &f) to ensure that the
        // 'GetLength()>=event_header.GetHeaderLength()' assertion is computed
        // correctly.  The range needs to include the 'size' fields of old and
        // new header and the 'version' field of the new header, therefore
        // first 4 header words need to be read at least.
        const static int ReadMinimum = 16;

        enum
        {
          ATTRIBUTE_WORDS       =  2,
          MAX_DETECTOR          = 94,
          SUPER_EVENT_MASK      = 1<<31,
          MASK_LENGTH           = (MAX_DETECTOR/32)+1,
          EVENT_TYPE_MASK       = 0x0000FFFF,
          EVENT_MAGIC_NUMBER    = 0xDA1E5AFE,
          EVENT_SWAPPED         = 0x00020000
        };

        /// Create a fake header.
                        Header                  (uint32 const * const buffer) : buf(buffer) {}
                        Header                  (uint8 const * const buffer) : buf(reinterpret_cast<uint32 const * const>(buffer)) {}

      private:

        struct HeaderDATE_old
        {
          void Print(const std::string &prefix) const;

          int32           size;                   ///< size of event in Bytes
          uint32          magic;                  ///< magic number used for consistency check
          uint32          type;                   ///< event type
          uint32          headLen;                ///< size of header in bytes
          uint32          runNb;                  ///< run number
          uint32          burstNb;                ///< burst number
          uint32          nbInRun;                ///< event number in run
          uint32          nbInBurst;              ///< event number in burst
          uint32          triggerNb;              ///< trigger number for this detector
          uint32          fileSeqNb;              ///< file sequence number for multifiles run */
          uint32          detectorId[MASK_LENGTH];///< detector identification */
          uint32          time;                   ///< time in seconds since 0.00 GMT 1.1.1970
          uint32          usec;                   ///< microseconds
          uint32          errorCode;
          uint32          deadTime;
          uint32          deadTimeusec;
          uint32          typeAttribute[ATTRIBUTE_WORDS]; ///< event type id mask
        };

        struct HeaderDATE_36
        {
          void Print(const std::string &prefix) const;
        
          enum constants {ATTR_SUPER_EVENT = 68};
          uint32          size;
          uint32          magic;
          uint32          headLen;
          uint32          version;
          uint32          type;
          uint32          runNb;                  ///< run number
          uint32          event_id[2];
          uint32          event_trigger_pattern[2];
          uint32          event_detector_pattern;
          uint32          typeAttribute[3];
          uint32          event_ldc_id;
          uint32          event_gdc_id;
          uint32          time;
        };

        const HeaderDATE_old &  GetHeaderOld(void) const {return *reinterpret_cast<const HeaderDATE_old*>(buf);}
        const HeaderDATE_36  &  GetHeader36 (void) const {return *reinterpret_cast<const HeaderDATE_36 *>(buf);}

      public:

        bool            operator ==         (const Header &h) const;
        bool            operator !=         (const Header &h) const {return !(*this==h);}

      public:

        void            Print               (const std::string &prefix="") const;

        uint32          GetVersion          (void) const {return GetHeader36().version;}
        
        uint32          GetMagic            (void) const {return GetVersion()<0xffff ? GetHeaderOld().magic : GetHeader36().magic;}

        uint32          GetLength           (void) const;
        
        uint32          GetHeaderLength     (void) const {return GetVersion()<0xffff ? GetHeaderOld().headLen : GetHeader36().headLen;}
        
        std::pair<time_t,uint32> GetTime     (void) const {return GetVersion()<0xffff ? std::pair<time_t,uint32>(GetHeaderOld().time,GetHeaderOld().usec)
                                                                                      : std::pair<time_t,uint32>(GetHeader36 ().time,0);}
        std::pair<time_t,uint32> GetDeadTime (void) const {return GetVersion()<0xffff ? std::pair<time_t,uint32>(GetHeaderOld().deadTime,GetHeaderOld().deadTimeusec)
                                                                                      : throw "DaqEvent::Header::GetDeadTime(): no information.";}
    
        uint32          GetBurstNumber      (void) const {return GetVersion()<0xffff ? GetHeaderOld().burstNb : (GetHeader36().event_id[1]>>20)&0x00000fff;}
    
        uint32          GetRunNumber        (void) const {return GetVersion()<0xffff ? GetHeaderOld().runNb : GetHeader36().runNb;}
    
        uint32          GetEventNumberInRun (void) const {return GetVersion()<0xffff ? GetHeaderOld().nbInRun : GetHeader36().event_id[0];}

        uint32          GetEventNumberInBurst(void)const {return GetVersion()<0xffff ? GetHeaderOld().nbInBurst : GetHeader36().event_id[1]&0x000fffff;}

        uint32          GetTriggerNumber    (void) const {return GetVersion()<0xffff ? GetHeaderOld().triggerNb : throw "DaqEvent::Header::GetTriggerNumber(): no information.";}

        uint32          GetErrorCode        (void) const {return GetVersion()<0xffff ? GetHeaderOld().errorCode : throw "DaqEvent::Header::GetErrorCode(): no information.";}

        EventType       GetEventType        (void) const {return EventType(GetVersion()<0xffff ? GetHeaderOld().type : GetHeader36().type);}
        
        uint32          GetTrigger          (void) const;

        bool            HaveSubEvents       (void) const;
        
        const uint32 *  GetTypeAttributes   (void) const {return GetVersion()<0xffff ? GetHeaderOld().typeAttribute : GetHeader36().typeAttribute;}

      private:
        
        uint32 const * const buf;
    };
    
  //============================================================================
  // Constructors and destructor
  //============================================================================

  public:

    /// Destructor
                       ~DaqEvent                (void);


  private:
    /// Create a fake event.
                        DaqEvent                (void);

  public:

    /*! \brief Create DaqEvent using given stream of events in the format of DATE library.
        
        If it is not possible to read an event, exception ExceptionEndOfStream will be thrown.
        
        \exception ExceptionEndOfStream
        
        Example:
        \verbatim
          #include <fstream>            // C++ files I/O
          #include "DaqEvent.h"         // Interface to DaqDataDecoding Event
          //...
          ifstream f("file.dat");       // Create an input stream from data file.
          try                           // block to catch an errors/exceptions
          {
            while(true)                 // Cycle on all events in the file
            {
              DaqEvent event(f);        // Create event or throw an exception
                                        // if it is not possible to do this.
              // .... user code is here
            }
          }
          catch( ExceptionEndOfStream ) // Catch the situation when stream is empty
            { cout << "OK. End of stream is detected.\n"; }
          catch( ... )                  // Catch all other situations.
            { cerr << "Error!\n"; }
        \endverbatim
    */
                        DaqEvent                (std::istream &f);

    /*! \brief Create DaqEvent using given stream of events in the format of DATE library.

        If it is not possible to read an event, exception ExceptionEndOfStream will be thrown.

        \exception ExceptionEndOfStreamGetRunNumber

    */
                        DaqEvent                (FILE *f);

    /*! \brief Create DaqEvent using given buffer of DAQ Test Environment (DATE).
        
        \param event_buffer pointer to start of buffer with event information
        
        Event length will be taken from 4-bytes value event_buffer[0]
        
        This constructor will NOT copy any data from the buffer.
        When a request DaqEvent::GetXXX() is issued by a client, appropriate data
        will be extracted directly from the buffer.
        
        Any attempt to access data not in the range [buf,buf+length) will generate
        an exception.
    */

                        DaqEvent                (void const * const event_buffer,bool copy_buf=false);

  private:

    /// Copy constructor.
                        DaqEvent                (const DaqEvent &d);

  //============================================================================
  // Operators
  //============================================================================

  private:
  
    /// Assignment operator
    DaqEvent           &operator =              (const DaqEvent &d);

  //============================================================================
  // Methods
  //============================================================================

  public:

    const
    TriggerTime &       GetTT                   (void) const {return TT;}
    TriggerTime &       GetTT                   (void)       {return TT;}

    /// Print event
    void                Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;
    
    /// Clear the event. Chips list will be erased.
    void                Clear                   (void);

    void                Check                   (void) const;

    void                CheckSLink              (const SLink &slink);
    
    /// \return pointer to event buffer
    const uint32*       GetBuffer               (void) const {return reinterpret_cast<uint32*>(buf);}

    //void                New                     (const uint32 *buf,bool copy_buf=false);
    
    /// \return DaqEvent Header.
    const Header        GetHeader               (void) const {return Header(buf);}

    /// \return Event length in bytes.
    uint32              GetLength               (void) const {return Header(buf).GetLength();}
    
    /// \return Event header length in bytes.
    uint32              GetHeaderLength         (void) const {return Header(buf).GetHeaderLength();}
    
    /*! \return time in seconds since 0.00 GMT 1.1.1970 and (second parameter) microseconds

       Example:
       \verbatim
         DaqEvent *event=...
         time_t  sec = event->GetTime().first;    // seconds since 0.00 GMT 1.1.1970
         uint32 usec = event->GetTime().second;   // microseconds
       \endverbatim
    */
    std::pair<time_t,uint32> GetTime            (void) const {return Header(buf).GetTime();}

    /// \return time
    std::string         GetTimeStr              (void) const;

    /*! \return Get dead time
        \sa GetTime
    */
    std::pair<time_t,uint32> GetDeadTime        (void) const {return Header(buf).GetDeadTime();}
    
    /// \return Run number
    uint32              GetBurstNumber          (void) const {return Header(buf).GetBurstNumber();}

    /// \return Run number
    uint32              GetRunNumber            (void) const {return Header(buf).GetRunNumber();}

    /// \return Event number in run
    uint32              GetEventNumberInRun     (void) const {return Header(buf).GetEventNumberInRun();}

    /// \return Event number in burst
    uint32              GetEventNumberInBurst   (void) const {return Header(buf).GetEventNumberInBurst();}

    /// \return Get trigger number
    uint32              GetTriggerNumber        (void) const {return Header(buf).GetTriggerNumber();}
    
    /// \return Error code
    uint32              GetErrorCode            (void) const {return Header(buf).GetErrorCode();}

    /// \return event type
    EventType           GetType                 (void) const {return Header(buf).GetEventType();}
    
    /*! \brief Read all chips found in the event to the internal list.
        
        New chips will be added to the internal list. Call methot Clear() to
        clear this list.
        
        Example:
        \verbatim
          DaqEvent event(buf);   // 'buf' is a memory pointer to a DATE event
          event.ReadChips(); //
          cout << "It was read "
               << event.GetChips().size()
               << " chips\n";
        \endverbatim
    */
    void                ReadChips               (DaqOption &opt,DaqEventCallMe f=NULL);

    /// The same as previous function, but uses the default options.
    void                ReadChips               (DaqEventCallMe f=NULL);
    
    void                SetTriggerMask          (const DaqOption &options);
    
  public:
    
    void                ChipCreate              (const SLink *sl,bool copy_buf,DaqOption &opt,DaqEventCallMe call_me);

    /*! \return List of chips
        \sa ReadChips
    */
          std::vector<Chip*> &GetChips          (void)       {return chips;}

    /*! \return List of chips
        \sa ReadChips
    */
    const std::vector<Chip*> &GetChips          (void) const {return chips;}

    /*! \return errors related to the event parsing only (do not include chip-related problems)
    
        Errors will be \b added to the list \c el.
    */
    const DaqErrors&    GetDaqErrors            (void) const {return errors;}
          DaqErrors&    GetDaqErrors            (void)       {return errors;}
    
    /*! \brief Get trigger after a trigger mask (see DaqOption) is applied.
        
        Trigger mask is calculated only after a ReadChips() method call!
    */
    uint32              GetTrigger              (void) const {return GetHeader().GetTrigger() & trigger_mask;}
    
    Event1Run *         GetEvent1Run            (void) const {return first_event_of_run;}
    void                SetEvent1Run            (Event1Run *e) {first_event_of_run=e;}
  
    /// \return check that event is good for analysis
    bool                IsGood                  (void) const;
    
    OnlineFilter      & GetOnlineFilter         (void)       {return online_filter;}
    OnlineFilter const& GetOnlineFilter         (void) const {return online_filter;}
    
    const DaqEvent &    GetTopEvent             (void) const {return const_cast<DaqEvent*>(this)->GetTopEvent();}
          DaqEvent &    GetTopEvent             (void)       {return parent==NULL ? *this : parent->GetTopEvent();}

    bool                IsBadDataSource         (uint16 srcID,uint16 geoID=uint16(-1)) const {return srcID_geoID_bad.count((srcID<<16)+geoID)>0;}
    void                AddBadDataSource        (uint16 srcID,uint16 geoID=uint16(-1))       {srcID_geoID_bad.insert((srcID<<16)+geoID);}
    const
    std::set<uint32> &  GetBadDataSources       (void) const {return srcID_geoID_bad;}
    std::set<uint32> &  GetBadDataSources       (void)       {return srcID_geoID_bad;}
    
    void                AddError                (const DaqError &e) {GetTopEvent().errors.Add(e);}

  private:
  
    /// \return Do we have any subevents inside the event?
    bool                HaveSubEvents           (void) const {return GetHeader().HaveSubEvents();}
    
    void                ReadFromStream          (void *stream);

    /// Check that all SrcIDs do exist in the event.
    void                CheckSrcIDs             (DaqOption &opt);

  //============================================================================
  // Attributes, data
  //============================================================================

  private:

    /*! @brief Parent event. It is NULL for a top-level event. */
    DaqEvent           *parent;
  
    /// Do we need to delete the buffer in destructor?
    bool                flag_buffer_delete;

    /// Do we have rights to modify the buffer?
    bool                allow_buf_change;
    
    /// A pointer to the "first event of run" object.
    Event1Run *         first_event_of_run;

    /// Pointer to real data. 
    uint8              *buf;

    /*! \brief Internal list of chips.
        
        This list will be filled by method ReadChips().
    */
    std::vector<Chip*>  chips;
    
    DaqErrors           errors;
    
    /*! @brief List of srcIDs found in the given event
     
        The list is filled by the ReadChips() call.
    */
    std::set<uint16>    srcIDs;

    /*! @brief  List of ((srcID<<16)+geoID) which are known to have problems (we should not add data from those srcIDs+geoIDs). */
    std::set<uint32>    srcID_geoID_bad;

    OnlineFilter        online_filter;
    
    TriggerTime         TT;
    
    uint32              trigger_mask;
};

////////////////////////////////////////////////////////////////////////////////

inline std::ostream &operator << (std::ostream &o,const DaqEvent &e)
{
  e.Print(o);
  return o;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace CS

#endif // CompassSoft_DaqEvent__include
