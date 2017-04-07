#ifndef First1Run___include
#define First1Run___include

#include <map>
#include <iostream>
#include "CatchInfo.h"

namespace CS {

/*! \brief First event of run
    This is an example how to work with the Event1Run
    \code
        Event1Run first_event_of_run; // This is a special event!

        for( all events in a run )
        {
            DaqEvent event = .....;         // read an event somehow
            DaqOption daq_options = ......; // options for the decoding

            if( first_event_of_run.GetRun()!=event.GetRunNumber() )
            {
                if( first_event_of_run.GetRun()!=0 )
                {
                    // The previous run has finished, this is the new one.
                    // We can dump the first event of run to a file.
                    // An example of code:
                    
                    // Print to the standard output
                    first_event_of_run.Print();
                    
                    // Write to a file:
                    char file_name[22];
                    sprintf(file_name,"feor%u.dat",first_event_of_run.GetRun());
                    ofstream f(file_name);
                    if( !f.is_open() )
                        throw "Can not open file for writing first-event-of-run";
                    first_event_of_run.Write(f);
                }

                event.SetEvent1Run(&first_event_of_run);            // Where to store the info
                first_event_of_run.NewRun(event.GetRunNumber());    // Be ready for a next events
            }
            
            // Read the event structure
            event.ReadChips(daq_options);
            
            // If it was a first-event-of-run it is decoded to first_event_of_run
        }
    \endcode
*/
class Event1Run
{
  public:
    /*! Add new catch info
      \arg buf is a pointer to the S-Link followed by header and data
    */
                       ~Event1Run               (void) {Clear();}
                        Event1Run               (unsigned int r=0) : run(r) {}

    void                Print                   (const char *prefix="") const;
    void                Add                     (unsigned int run,const void *buf);
    void                Clear                   (void);
    void                NewRun                  (unsigned int r) {Clear(); run=r;}

    unsigned int        GetRun                  (void) const {return run;}

    const std::map<int,CatchInfo*> &
                        GetCatchesInfo          (void) const {return catch_info;}

    void                Write                   (std::ostream &o) const;
    void                Read                    (std::istream &o);

  private:

    unsigned int        run;
    std::map<int,CatchInfo*> catch_info;
};

} // namespace CS

#endif // First1Run___include
