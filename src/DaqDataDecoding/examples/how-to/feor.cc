#include <cstdio>
#include <fstream>

#include "DaqEventsManager.h"
#include "Event1Run.h"
#include "ChipGassiplex.h"
#include "utils.h"

using namespace std;
using namespace CS;

///////////////////////////////////////////////////////////////////////////////

int main(int argc,char **argv)
{
    if( argc!=2 )
    {
        printf("Usage: <EXE> data\n");
        return 1;
    }

    try
    {
        DaqEventsManager manager;
        manager.AddDataSource(argv[1]);
        DaqOption daq_options;              // default options
        DaqEvent &event=manager.GetEvent(); // just a simple access to the event

        while( manager.ReadEvent() )
        {
            if( manager.GetEventsCounter()%100==0 )
                cout << "It were read " << manager.GetEventsCounter() << " events.\n";
        
            // To speed up the process we discard not intresting events
            // this check did not work:  event.GetEventNumberInRun()!=1
            if( event.GetEventNumberInBurst() !=1 ||
                event.GetBurstNumber()        !=1 )
                continue;


            event.ReadChips(daq_options);
                        
            if( manager.GetEvent1Run().GetCatchesInfo().size()==0 )
            {
                Exception("Event1Run is empty! Why!? run=%d.\n",event.GetRunNumber()).Print();
                continue;
            }

            manager.GetEvent1Run().Print();
            
            for( map<int,CatchInfo*>::const_iterator it=manager.GetEvent1Run().GetCatchesInfo().begin();
                 it!=manager.GetEvent1Run().GetCatchesInfo().end(); it++ )
            {
                if( it->second->GetSender()==CatchInfo::RICH )
                {
                    it->second->Print();

                    const uint32 *p = (uint32*)it->second->GetData();

                    for( int rest_length=it->second->GetDataSize()/4; rest_length!=0; )
                    {
                        uint32 len=rest_length;

                        try
                        {
                            ChipGassiplex::CathodeCalib cathode_calib(p,len);
                            cathode_calib.Print();
                            
                            for( set<ChipGassiplex::PixelCalib>::const_iterator it=cathode_calib.GetPixels().begin();
                                 it!=cathode_calib.GetPixels().end(); it++ )
                            {
                                // PixelCalib access-functions
                                it->GetChannelID();
                                it->GetCathode();
                                it->GetX();
                                it->GetY();
                                it->GetThreshold();
                                it->GetSigma();
                            }
                            
                        }
                        catch(...)
                        {
                            printf("Problem! Bad words 0-%d:\n",rest_length);
                            for( int i=-10; i<rest_length; i++ )
                            {
                                printf("%4d: ",i);
                                bits_print(cout,*(uint32*)(p+i*4),0,32,"%8x\n");
                            }
                        }
 
                        rest_length-=len;
                        p+=len;
                        
                        assert(rest_length>=0);
                    }
                }
            }
            
            void save_feor(const Event1Run &e);
            save_feor(manager.GetEvent1Run());
        }
    }
    catch( DaqEvent::ExceptionEndOfStream )
    {
        // This is the normal exit from the loop
        cout << "End of data\n";
    }

    // Something is wrong...
    // Print error message and finish this data file.

    catch( const std::exception &e )
    {
        cerr << "exception:\n" << e.what() << "\n";
    }
    catch( const char *s )
    {
        cerr << s << "\n";
    }
    catch( ... )
    {
        cerr << "Oops, unknown exception!\n";
    }
    
    CS::Exception::PrintStatistics();
    CS::Chip::PrintStatistics();

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

void save_feor(const Event1Run &e)
{
    // Print to the standard output
    //e.Print();

    // Write to a file:
    char file_name[22];
    sprintf(file_name,"feor%u.dat",e.GetRun());
    ofstream f(file_name);
    if( !f.is_open() )
        throw "Can not open file for writing first-event-of-run";
    e.Write(f);
}

///////////////////////////////////////////////////////////////////////////////
