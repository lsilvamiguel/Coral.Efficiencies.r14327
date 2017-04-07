#include <cstdio>
#include <exception>
#include <string>
#include <set>
#include <popt.h>
#include <dlfcn.h>
#include <fstream>

#include "config.h"

#include "utils.h"
#include "DaqEventsManager.h"

#if USE_RFIO
#include <shift.h>
#endif

using namespace std;
using namespace CS;

////////////////////////////////////////////////////////////////////////////////

typedef void (*plugin_type)(const DaqEventsManager &);

plugin_type load_plugin(const char *path,const char *name="ddd_call")
{
    if( path==NULL || *path==0 )
        throw "load_plugin(): can not load an empty plugin.";

    printf("Loading plugin %s ...\n",path);
    void *handle=NULL;

    char *error=NULL;

    handle = dlopen (path, RTLD_NOW);
    if (!handle)
        throw dlerror();

    plugin_type ddd_call = (plugin_type) dlsym(handle, name);
    if ((error = dlerror()) != NULL)
        throw error;

    //dlclose(handle);
    printf("Plugin %s has been loaded successfully.\n",path);
    
    return ddd_call;
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc,const char *argv[])
{
    unsigned char error = 0; // success.
    FILE *f_out=NULL;   // output file (if any)

    try
    {
        vector<plugin_type> plugins;
        set<string> dets_digit;
        char *file_out         = strdup("");
        char *maps             = strdup("/afs/cern.ch/compass/detector/maps/");
        char *input_files_list = strdup("");
        int errors_list=0, event_errors=0, chips_stat=0, event_header=0, filter_CAL=0, stat_TT=0;
        const uint32 trigger_mask_all_bits=0xffffffff;
        uint32 trigger_mask=trigger_mask_all_bits;
        DaqEventsManager manager;

        struct poptOption options[] =
        {
            { "events",     '\0',   POPT_ARG_INT,
                                    NULL, 'N',
                                    "Maximum number of events to be analyzed", "NUMBER" },
            { "filter-CAL",'\0',    POPT_ARG_NONE,
                                    &filter_CAL,  0,
                                    "Use filter for calorimter events calibration (random trigger and calibration events).", "" },
            { "out",        '\0',   POPT_ARG_STRING|POPT_ARGFLAG_SHOW_DEFAULT,
                                    &file_out,  0,
                                    "Name of the output raw data file (filtering for random trigger or calibration data).", "PATH" },
            { "input",      '\0',   POPT_ARG_STRING|POPT_ARGFLAG_SHOW_DEFAULT,
                                    &input_files_list,  'l',
                                    "File name with the list of input raw files.", "PATH" },
            { "maps",       '\0',   POPT_ARG_STRING|POPT_ARGFLAG_SHOW_DEFAULT,
                                    &maps,  'm',
                                    "Map file or direcory", "PATH" },
            { "det",        '\0',   POPT_ARG_STRING,
                                    NULL, 'd',
                                    "Print digits from this detector (regular expr.)", "NAME" },
            { "event-errors",'\0',  POPT_ARG_INT|POPT_ARGFLAG_SHOW_DEFAULT,
                                    &event_errors, 0,
                                    "Print the list of encounted problems per event. 0-no print, 1-errors summary, 2-detailed list", "NUMBER" },
            { "errors",     '\0',   POPT_ARG_NONE,
                                    &errors_list, 0,
                                    "Print the list of encounted problems at the end", "" },
            { "stat",       '\0',   POPT_ARG_NONE,
                                    &chips_stat, 0,
                                    "Print chips statistics", "" },
            { "stat-TT",    '\0',   POPT_ARG_NONE,
                                    &stat_TT, 0,
                                    "Print trigger time statistics (for calibration)", "" },
            { "plugin",     '\0',   POPT_ARG_STRING,
                                    NULL, 'p',
                                    "User plugin code (file.so) with function ddd_call(const DaqEvent &)", "PATH" },
            { "trigger",    '\0',   POPT_ARG_INT|POPT_ARGFLAG_SHOW_DEFAULT,
                                    &trigger_mask, 0,
                                    "Trigger mask (skip events if mask does not fit)", "NUMBER" },
            { "sep",        '\0',   POPT_ARG_NONE,
                                    &event_header, 0,
                                    "Separate events on screen", "" },
            POPT_AUTOHELP
            POPT_TABLEEND
        };

        poptContext poptcont=poptGetContext(NULL,argc,argv,options,0);
        poptSetOtherOptionHelp(poptcont,
            "<options...> <data-file> [<data-file> <data-file> ...]\n"
            "  read COMPASS raw data.\n"
            "  Author: Alexander Zvyagin <Alexander.Zvyagin@cern.ch>\n"
        );

        int rc;
        while( (rc=poptGetNextOpt(poptcont))!=-1 )
        {
            switch( rc )
            {
                case 'N':
                    manager.SetEventsMax( atoi(poptGetOptArg(poptcont)) );
                    break;

                case 'd':
                    dets_digit.insert(poptGetOptArg(poptcont));
                    break;

                case 'm':
                    manager.SetMapsDir(poptGetOptArg(poptcont));
                    break;

                case 'p':
                    plugins.push_back(load_plugin(poptGetOptArg(poptcont)));
                    break;

                case 'l':
                {
                    const char *s = poptGetOptArg(poptcont);
                    ifstream f(s);
                    if( !f.is_open() )
                    {
                        printf("%s\n",s);
                        throw "Can not open a list file.";
                    }

                    string fname;
                    while( f>>fname )
                        manager.AddDataSource(fname);
                    break;
                }

                default:
                    fprintf(stderr, "bad argument %s: %s\n",
		            poptBadOption(poptcont, POPT_BADOPTION_NOALIAS),
		            poptStrerror(rc));
                    throw "";
            }
        }

        // Get all data files
        while(true)
        {
            const char *data_file=poptGetArg(poptcont);

            if( data_file==NULL )
                break;
            
            manager.AddDataSource(data_file);
        }

        if( manager.GetDataSources().empty() )
        {
            poptPrintHelp(poptcont,stdout,0);
            return 1;
        }

        manager.Print();

        if( file_out!=NULL && *file_out!=0 )
        {
            f_out = fopen(file_out,"w");
            if( f_out==NULL )
            {
                //printf("---------\"%s\"\n",
                throw "Can not open output file.";
            }
        }

        // This is the main loop over events.
        while( manager.ReadEvent() )
        {
            if( f_out!=NULL )
            {
                int
                    trigger = manager.GetEvent().GetTrigger()&0xffff,
                    type    = manager.GetEvent().GetHeader().GetEventType()&DaqEvent::Header::EVENT_TYPE_MASK;

                bool write_event=true;

                if( filter_CAL )
                    if( !(trigger&(1<<11)) // Random trigger
                        &&
                        type!=DaqEvent::CALIBRATION_EVENT )
                        write_event = false;

                if( write_event )
                {
                    if( 1!=fwrite((void*)manager.GetEvent().GetBuffer(),
                                  manager.GetEvent().GetLength(),1,f_out) )
                    {
                        throw "Can not write to the output file!";
                    }
                }

                continue;
            }

            if( trigger_mask!=trigger_mask_all_bits &&
               !(manager.GetEvent().GetTrigger()&trigger_mask ) )
                continue;

            if( event_header )
            {
                printf("*********************************************\n");
                printf("NEW EVENT. Trigger: ");
                bits_print(cout,manager.GetEvent().GetTrigger());
                printf("\n");
                manager.GetEvent().GetHeader().Print();
                printf("\n");
            }

            if( event_header )
            {
                printf("Decoding started ...\n");
                fflush(stdout);
            }
            bool decoded=manager.DecodeEvent();
            if( event_header )
                printf("Decoding %s\n", decoded?"OK":"failed" );

            for( vector<DaqError>::const_iterator e=manager.GetEvent().GetDaqErrors().errors.begin();
                 e!=manager.GetEvent().GetDaqErrors().errors.end(); e++ )
            {
                try
                {
                    CS::Exception v=*e;
                }
                catch(const std::exception &e)
                {
                    printf("%s\n",e.what());
                }
            }

            //if( !decoded )
            //    continue;
            
            if( dets_digit.size()>0 )
            {
                typedef multimap<CS::DetID,CS::Chip::Digit*>::iterator m_it; // iterator type

                // Loop on all found digits
                for( m_it d_it=manager.GetEventDigits().begin(); d_it!=manager.GetEventDigits().end(); d_it++ )
                    // Loop on all patterns
                    for( set<string>::const_iterator p=dets_digit.begin(); p!=dets_digit.end(); p++ )
                        // Does the detector name match the pattern?
                        if( string_match(d_it->first.GetName(),*p) )
                        {
                            // Yes! Print the digit from the detector.
                            d_it->second->Print();
                        }
            }

            for( size_t i=0; i<plugins.size(); i++ )
                plugins[i](manager);
            
            if( event_errors>0 )
            {
                DaqErrors &e = manager.GetEvent().GetDaqErrors();
                // print event errors
                printf("Event errors: total=%zu  problems(NO,WARN,ERR,CRIT)=(%zu,%zu,%zu,%zu)  actions(NO,GEO,SRC,EVENT)=(%zu,%zu,%zu,%zu)\n",
                       e.errors.size(),
                       e.errors_level[DaqError::OK              ].size(),
                       e.errors_level[DaqError::MINOR_PROBLEM   ].size(),
                       e.errors_level[DaqError::WARNING         ].size(),
                       e.errors_level[DaqError::SEVERE_PROBLEM  ].size(),
                       e.errors_action[DaqError::Action::NOTHING                ].size(),
                       e.errors_action[DaqError::Action::DISCARD_EVENT_ON_GEOID ].size(),
                       e.errors_action[DaqError::Action::DISCARD_EVENT_ON_SRCID ].size(),
                       e.errors_action[DaqError::Action::DISCARD_EVENT          ].size());

                if( event_errors>1 )
                {
                    //printf("Sorry, not implemented yet. Send a mail to Zvyagin.Alexander@gmail.com\n");
                    printf("Errors sorted by LEVEL:\n");
                    printf("==> Error level:  NO ERROR: %zu errors\n",e.errors_level[DaqError::OK              ].size());
                    e.Print(DaqError::OK);
                    printf("==> Error level:  MINOR PROBLEM: %zu errors\n",e.errors_level[DaqError::MINOR_PROBLEM   ].size());
                    e.Print(DaqError::MINOR_PROBLEM);
                    printf("==> Error level:  WARNING: %zu errors\n",e.errors_level[DaqError::WARNING         ].size());
                    e.Print(DaqError::WARNING);
                    printf("==> Error level:  SEVERE_PROBLEM: %zu errors\n",e.errors_level[DaqError::SEVERE_PROBLEM  ].size());
                    e.Print(DaqError::SEVERE_PROBLEM);
                    
                    printf("\nErrors sorted by ACTION:\n");
                    printf("==> Error action:  NO ACTION: %zu errors\n",e.errors_action[DaqError::Action::NOTHING].size());
                    e.Print(DaqError::Action::NOTHING);
                    printf("==> Error action:  DISCARD_EVENT_ON_GEOID: %zu errors\n",e.errors_action[DaqError::Action::DISCARD_EVENT_ON_GEOID].size());
                    e.Print(DaqError::Action::DISCARD_EVENT_ON_GEOID);
                    printf("==> Error action:  DISCARD_EVENT_ON_SRCID: %zu errors\n",e.errors_action[DaqError::Action::DISCARD_EVENT_ON_SRCID].size());
                    e.Print(DaqError::Action::DISCARD_EVENT_ON_SRCID);
                    printf("==> Error action:  DISCARD_EVENT: %zu errors\n",e.errors_action[DaqError::Action::DISCARD_EVENT].size());
                    e.Print(DaqError::Action::DISCARD_EVENT);

                    printf("- - - - - - -    End of event errors report   - - - - - -\n\n");
                }
            }
        }

        if( errors_list )
            CS::Exception::PrintStatistics();

        if( chips_stat )
            CS::Chip::PrintStatistics();

        if( stat_TT )
        {
            manager.GetDaqOptions().GetStatTT().Print();

            int events_max=50;
            printf("Limiting statics to %d events.\n",events_max);
            manager.GetDaqOptions().GetStatTT().SetBufferSize(events_max);
            manager.GetDaqOptions().GetStatTT().Print();
        }

        printf("\n\nEnd of the decoding!\n\n");
    }
    catch(const char *e)
    {
        printf("%s\n",e);
        error = 1;
    }
    catch(const std::exception &e)
    {
        printf("%s\n",e.what());
        error = 1;
    }
    catch(...)
    {
        printf("Unknown exception.\n");
        error = 1;
    }

    if( f_out!=NULL )
        fclose(f_out);

    return error;
}

////////////////////////////////////////////////////////////////////////////////
