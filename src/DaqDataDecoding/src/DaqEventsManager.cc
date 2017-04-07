#include <cstdlib>
#include <csignal>
#include <cassert>
#include <unistd.h>
#include "config.h"

#include "DaqError.h"
#include "DaqEventsManager.h"
#include "ChipF1.h"
#include "TriggerTime.h"

#if USE_DATE_LIB
#include "monitoring/monitor.h"
#endif

#if USE_RFIO
#include <shift.h>
#endif

using namespace std;

///////////////////////////////////////////////////////////////////////////////


namespace CS {

bool DaqEventsManager::flag_end=false;

////////////////////////////////////////////////////////////////////////////////

DaqEventsManager::~DaqEventsManager(void)
{
    Clear();
}

////////////////////////////////////////////////////////////////////////////////

DaqEventsManager::DaqEventsManager(const string &s,bool set_signal_handler) :
    name(s),
    event(NULL),
    event_ok(false),
    good_events_only(false),
    buffer(NULL),
    in_stream(NULL),
    trigger_mask(unsigned(-1))
{
    Clear();

    #if USE_DATE_LIB
    int status;
    status = monitorDeclareMp((char*)name.c_str());
    if ( status!=0 )
    {
        string message(monitorDecodeError(status));
        throw Exception("DaqEventsManager::DaqEventsManager(): DATE::monitorDeclareMp(): %s",message.c_str());
    }
    #endif
  
    if( set_signal_handler )
    {
        signal(0,       SignalHandler);
        signal(SIGHUP,  SignalHandler);
        signal(SIGINT,  SignalHandler);
        signal(SIGQUIT, SignalHandler);
        signal(SIGALRM, SignalHandler);
        signal(SIGTERM, SignalHandler);
    }
}

////////////////////////////////////////////////////////////////////////////////

void DaqEventsManager::Clear(void)
{
    sources.clear();
    source_number = size_t(-1);
    SetEventsMax(uint32(-1));  // very big number
    events_counter = 0;
    events_in_curr_src_counter = 0;
    maps_dir="";
    run_number=0;
    digits.Clear();
    detectors_all.clear();
    first_event_of_run.Clear();
    event_ok=false;
    tt_decoding_and_subtraction = true;
    
    if( in_stream!=NULL )
    {
        fclose(in_stream);
        in_stream=NULL;
    }
    free(buffer);
    buffer=NULL;
    
    delete event;
    event=NULL;
}

////////////////////////////////////////////////////////////////////////////////

void DaqEventsManager::SignalHandler(int n)
{
    printf("\n"
           "============================================\n"
           "=== The program has received signal: %3d ===\n"
           "============================================\n\n",n);
    if( flag_end )
    {
        printf("Forcing exit.\n\n");
        abort();
    }
    else
        flag_end = true;
}

////////////////////////////////////////////////////////////////////////////////

void DaqEventsManager::AddDataSource(const char *src[],int n)
{
  for(int i=0; i<n; i++)
    sources.push_back(src[i]);
}

////////////////////////////////////////////////////////////////////////////////

bool DaqEventsManager::DecodeEvent(void)
{
    if( GetEvent().GetBuffer()==NULL )
        return false;

    event_ok=false;
    digits.Clear();

    try
    {
        // make sure everything is initialized properly
        GetEvent().GetTT().InitAndClear( GetDaqOptions().GetTTConfig() );

        GetEvent().ReadChips(opts_of_the_run);

        for( vector<Chip*>::iterator chip=GetEvent().GetChips().begin();
             chip!=GetEvent().GetChips().end(); chip++ )
            (*chip)->Decode(maps_of_the_run,digits,opts_of_the_run);
        if( tt_decoding_and_subtraction )
            GetEvent().GetTT().DecodeAndSubtract(digits,GetDaqOptions(),GetEvent().GetDaqErrors());

        event_ok = GetEvent().IsGood();
        return event_ok;
    }
    catch( const std::exception & e)
    {
        cerr << e.what() << "\n";
    }
    catch( const char * s )
    {
        cerr << s << "\n";
    }
    catch( ... )
    {
        cerr << "Unknown exception.\n";
    }

    printf("Event decoding failed for the trigger mask %d (16 bits are shown).\n",GetEvent().GetTrigger()&0xffff);

    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool DaqEventsManager::ReadEvent(void)
{
    begin:

    if( flag_end || events_counter>=events_counter_max )
        return false;

    if( events_in_curr_src_counter==0 )
        if( !NextDataSource() )
            return false;

    events_in_curr_src_counter++;

    try
    {
        #if USE_DATE_LIB
        // Read event from a DATE stream.
        free(buffer);
        buffer=NULL;
        int status;
        status = monitorGetEventDynamic(&buffer);    // An attempt to get a next event.
        if( status!=0 )
        {
          string message(monitorDecodeError(status));
          Exception("DATE::monitorGetEventDynamic(): %s",message.c_str()).Print(cerr);

          //free(buffer);   // DATE library crashes sometimes. 
          buffer=NULL;
        }

        if( buffer==NULL )
        {
          // There is no more data. Try next data source.
          events_in_curr_src_counter=0;
          goto begin;
        }

        events_counter++;
        delete event;
        event = NULL;
        event = new DaqEvent(buffer);

        #else

        assert(in_stream!=NULL);
        events_counter++;
        delete event;
        event = NULL;
        event = new DaqEvent(in_stream);

        #endif
    }
    catch(...)
    {
        // There are no more data or there is a problem with the current data source.
        // Try next data source.
        
        events_counter--;
        events_in_curr_src_counter=0;
        
        assert(in_stream!=NULL);
        fclose(in_stream);
        in_stream=NULL;
        delete event;
        event = NULL;
        
        goto begin;
    }
    
    
    try
    {
        if( GetEvent().GetBuffer()!=NULL && GetEvent().GetRunNumber()!=run_number )
        {
            run_number = GetEvent().GetRunNumber();
            first_event_of_run.NewRun(run_number);    // Be ready for next events
            maps_of_the_run.Clear();
            opts_of_the_run.Clear();
            detectors_all.clear();
            if( maps_dir.length()>0 )
                if( !Chip::ReadMaps( run_number, maps_dir, maps_of_the_run, opts_of_the_run, detectors_all ) )
                    printf("WARNING!!!! DaqEventsManager::ReadEvent(): mapping files(s) problem!\n");
        }

        GetEvent().SetEvent1Run(&first_event_of_run);  // Where to store the info

        return true;
    }
    catch(...)
    {
        printf("DaqEventsManager::ReadEvent(): Fix the map file problem(s) first!\n");
        return false;
    }
    
    goto begin;
}

////////////////////////////////////////////////////////////////////////////////

bool DaqEventsManager::SetDaqEvent(uint8 *pExternalBuffer) {
	if (event) delete event;
	event = 0;
	event = new DaqEvent(pExternalBuffer);
	
	return ( event != 0 );
}

////////////////////////////////////////////////////////////////////////////////

// \return false if it not possible
bool DaqEventsManager::NextDataSource(void)
{
    events_in_curr_src_counter=0;
    for( source_number++; source_number<sources.size(); source_number++ )
    {
        const char *ss = sources[source_number].c_str();
        cout << "Data stream \"" << ss << "\"\n";
        
        #if USE_DATE_LIB
        int status;
        status = monitorSetDataSource((char*)ss);
        if ( status!=0 )
        {
          cerr << "DaqEventsManager::NextDataSource(): DATE::monitorSetDataSource(): "
               << monitorDecodeError(status) << "\n";
          continue; // Try next data source
        }
        return true;
        #else
        assert(in_stream==NULL);
        in_stream = fopen((char*)ss,(char *)"rb");
        if( in_stream==NULL )
        {
            cerr << "DaqEventsManager::NextDataSource(): Can not open file \"" << ss << "\"\n";
            continue; // Try next data source
        }
        return true;
        #endif
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////

void DaqEventsManager::Print(ostream &o,const string &prefix) const
{
  char s[200];
  snprintf(s,200,"%sDaqEventsManager name .......... %s\n", prefix.c_str(), GetName().c_str());
  o<<s;
  snprintf(s,200,"%sMaximum events to read.......... %zu\n",prefix.c_str(),events_counter_max);
  o<<s;
  snprintf(s,200,"%sTotal events counter ........... %zu\n",prefix.c_str(),events_counter);
  o<<s;
  snprintf(s,200,"%sEvents counter in current src .. %zu\n",prefix.c_str(),events_in_curr_src_counter);
  o<<s;
  snprintf(s,200,"%sTotal number of sources ........ %zu\n",prefix.c_str(),sources.size());
  o<<s;
  snprintf(s,200,"%sCurrent source number is ....... %zu\n",prefix.c_str(),source_number+1);
  o<<s;
  for( size_t i=0; i<sources.size(); i++ )
  {
    snprintf(s,200,"%sSource %3zu:  %s\n",prefix.c_str(),i+1,sources[i].c_str());
    o<<s;
  }
}

////////////////////////////////////////////////////////////////////////////////

} // namespace CS
