#ifndef CompassSoft_DaqEventsManager__include
#define CompassSoft_DaqEventsManager__include

#include <cstdio>
#include <iostream>
#include <string>
#include <vector>
#include "config.h"

#include "DaqEvent.h"
#include "DaqOption.h"
#include "Chip.h"
#include "Event1Run.h"

namespace CS {

////////////////////////////////////////////////////////////////////////////////

/*! \brief Manager for raw events reading.

    This manager can read events via
    - DATE library
    - RFIO package ("/castor/dir.../file")
    - direct file

    \code

      try
      {
        CS::DaqEventsManager manager("DaqDataDecoding");

        manager.AddDataSource("file.data");
        manager.Print();

        while(1)
        {
          CS::DaqEvent event = manager.GetEvent();
          event.Print();
        }
      }
      catch( CS::DaqEvent::ExceptionEndOfStream )
      {
        // This is the normal exit from the loop
        cout << "End of data\n";
      }

      // Something is wrong...
      // Print error message and finish this data file.

      catch( ... )
      {
        cerr << "An error!\n";
      }

    \endcode

    \author Alexander Zvyagin
*/
class DaqEventsManager
{
  //============================================================================
  // Types, constants
  //============================================================================

  public:
        
  //============================================================================
  // Constructors and destructor
  //============================================================================

  public:

    /// Destructor
    virtual            ~DaqEventsManager        (void);

    /// Base constructor
                        DaqEventsManager        (const std::string& name="",bool set_signal_handler=true);

  private:

    /// You can not use the copy constructor.
                        DaqEventsManager        (const DaqEventsManager&);

  //============================================================================
  // Operators
  //============================================================================

  private:

    /// You can not use the assignment operator.
    DaqEventsManager&   operator =              (const DaqEventsManager&);

  //============================================================================
  // Methods
  //============================================================================

  public:
  
    /// Print information about manager.
    void                Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;
    
    /// Clear events manager.
    void                Clear                   (void);
    
    /// \return Name of the manager.
    const std::string&  GetName                 (void) const {return name;}
    
    /// Set maximum amount of events to be read.
    void                SetEventsMax            (size_t n) {events_counter_max=n;}
    
    /// \return events amount read by manager
    size_t              GetEventsCounter        (void) const {return events_counter;}

    uint32              GetTriggerMask          (void) const  {return trigger_mask;}
    void                SetTriggerMask          (uint32 mask) {trigger_mask=mask;}

    /*! \brief add data source(s)
        \param src array of character strings with names
        \param n   number of elements in the above array
    */
    void                AddDataSource           (const char *src[],int n);

    /*! \brief add data source
    */
    void                AddDataSource           (const std::string &src) {sources.push_back(src);}
    const std::vector<std::string> &GetDataSources        (void) const {return sources;}
    
    /// \return reference to the event. This method does not read a next event!
    DaqEvent const &    GetEvent                (void) const {if(event==NULL) throw "DaqEventsManager::GetEvent(): no event."; return *event;}
    DaqEvent &          GetEvent                (void)       {if(event==NULL) throw "DaqEventsManager::GetEvent(): no event."; return *event;}
    Event1Run const &   GetEvent1Run            (void) const {return first_event_of_run;}
    Event1Run       &   GetEvent1Run            (void)       {return first_event_of_run;}

	/// \return true if an event is available
	bool                IsEventAvailable        (void) const {return event!=0;}

    /// Decode current event
    bool                DecodeEvent             (void);

    // Read next event
    bool                ReadEvent               (void);

    bool                SetDaqEvent             (uint8 *pExternalBuffer);

    void                OptGoodEventsDeliver    (bool o=true) {good_events_only=o;}
    
    void                SetMapsDir              (const std::string &maps) {maps_dir=maps;}
    const std::string & GetMapsDir              (void) const {return maps_dir;}
    const std::vector<std::string> & GetDetectors         (void) const {return detectors_all;}
    
    Chip::Maps const &  GetMaps                 (void) const {return maps_of_the_run;}
    Chip::Maps       &  GetMaps                 (void)       {return maps_of_the_run;}
    
    DaqOption const &   GetDaqOptions           (void) const {return opts_of_the_run;}
    DaqOption       &   GetDaqOptions           (void)       {return opts_of_the_run;}
    
    Chip::Digits &      GetEventDigits          (void)       {return digits;}
    const Chip::Digits &GetEventDigits          (void) const {return digits;}

    bool                GetTTDecodingSubtraction(void) const {return tt_decoding_and_subtraction;}
    void                SetTTDecodingSubtraction(bool f=true) {tt_decoding_and_subtraction=f;}

  private:

    bool                NextDataSource          (void);
    static void         SignalHandler           (int n);

  //============================================================================
  // Attributes, data
  //============================================================================

  private:
  
    /// Name of the manager
    std::string         name;

    std::string         maps_dir;
    unsigned int        run_number;
    Chip::Maps          maps_of_the_run;
    DaqOption           opts_of_the_run;
    DaqEvent *          event;
    Chip::Digits        digits;
    std::vector<std::string>      detectors_all;
    Event1Run           first_event_of_run;
    
    bool                event_ok;
    bool                good_events_only;
    static bool         flag_end;

    bool                tt_decoding_and_subtraction;
    
    /// Sources
    std::vector<std::string>      sources;
    
    /// Current source number
    size_t              source_number;
    
    /*! \brief Maximum events counter

        \throw ExceptionEndOfStream in the attempt to read events_counter_max+1 event (from all sources)
    */
    size_t              events_counter_max;
    
    /// Events counter.
    size_t              events_counter;
    
    /// Number of events read from current source.
    size_t              events_in_curr_src_counter;
    
    /// Internal buffer.
    void*               buffer;
    
    /// Input stream;
    FILE*               in_stream;

    /// Accepted triggers
    uint32              trigger_mask;

};

////////////////////////////////////////////////////////////////////////////////

} // namespace CS

#endif // CompassSoft_DaqEvent__include
