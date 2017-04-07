#include "DaqEventsManager.h"
#include "ChipF1.h"
#include "TriggerTime.h"

#include "TFile.h"
#include "TNtuple.h"

using namespace CS;

extern "C"
{
void ddd_call(const CS::DaqEventsManager &m);
}

TFile *f=NULL;
TNtuple *nt=NULL;

void finish(void)
{
    if( f!=NULL )
        f->Write();
}

void ddd_call(const CS::DaqEventsManager &m)
{
    static bool init=true;
    if( init )
    {
        init=false;
        
        f  = TFile::Open("tt.root","RECREATE","",9);
        nt = new TNtuple("tt","trigger time","trig:hits:mt:tt0:tt1:tt2:tt3:tt4:tt5:tt6:tt7:tt8:tt9:tt10:tt11:tt12:tt13:tt14:tt15");
        
        atexit(finish);
    }

    const TriggerTime &tt = m.GetEvent().GetTT();
    

    vector<float> args;
    args.push_back(m.GetEvent().GetTrigger()&0xffff);
    args.push_back(tt.GetShifts().size());
    args.push_back(tt.GetTime(0));

    for( size_t i=0; i<16; i++ )
        args.push_back(0);
    
    for( size_t i=0; i<tt.GetShifts().size(); i++ )
    {
        unsigned n = 3+tt.GetShifts()[i].bit;
        assert( n<args.size() );
        args[n] = tt.GetShifts()[i].tt_mt_diff;
    }

    nt->Fill(&args.front());
}
