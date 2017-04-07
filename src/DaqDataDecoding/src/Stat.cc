#include <cassert>
#include <stdexcept>
#include <cmath>

#include "Stat.h"

namespace CS {

// =============================================================================

void Stat::SetBufferMaxSize(unsigned s)
{
    if( s<2 )
        throw std::runtime_error("Stat::SetBufferMaxSize(): buffer size must be >=2");

    buffer_max_size=s;

    // Number of elemnts to remove from begin and end of the buffer    
    while( buffer.size()>buffer_max_size )
    {
        buffer.erase(buffer.begin());
        buffer.erase(--buffer.end());
    }
}

// =============================================================================

void Stat::Add(float v)
{
    fills_counter++;
    buffer.insert(v);
    
    if( buffer.size()<=buffer_max_size )
        return;

    // Number of elemnts to remove from begin and end of the buffer    
    for( unsigned n = unsigned(buffer_max_size*cut_ratio/2); n>0; n-- )
    {
        buffer.erase(buffer.begin());
        buffer.erase(--buffer.end());
    }
}

// =============================================================================

float Stat::Mean(float events_fraction) const
{
    assert(events_fraction>0 && events_fraction<=1);

    unsigned int
        skip    = (unsigned int)((1-events_fraction)*buffer.size()),
        n_start = skip/2,
        n_end   = buffer.size()-n_start;

    float sum = 0;
    unsigned int i=0, n=0;
    for( std::multiset<float>::const_iterator v=buffer.begin(); v!=buffer.end(); v++ )
    {
        i++;

        if( i<n_start || i>=n_end )
            continue;

        sum += *v;
        n += 1;
    }

    if( n==0 )
        throw std::runtime_error("Stat::Mean(): number of events is 0.");
    
    return sum/n;
}

// =============================================================================

float Stat::Sigma(float events_fraction) const
{
    assert(events_fraction>0 && events_fraction<=1);

    unsigned int
        skip    = (unsigned int)((1-events_fraction)*buffer.size()),
        n_start = skip/2,
        n_end   = buffer.size()-n_start;

    float sum = 0, sum2 = 0;
    unsigned int i=0, n=0;

    for( std::multiset<float>::const_iterator v=buffer.begin(); v!=buffer.end(); v++ )
    {
        i++;
        if( i<n_start || i>=n_end )
            continue;

        sum  += *v;
        sum2 += (*v) * (*v);
        n += 1;
    }

    if( n<2 )
        throw std::runtime_error("Stat::Sigma(): number of events is <2.");

    const float mean = sum/n;
    
    return std::sqrt( sum2/n - mean*mean );
}

// =============================================================================

} // namespace CS



// =============================================================================

#define TEST 0

#if TEST

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TRandom.h"

TROOT root("","");

int main(void)
{
    TFile f("stat.root","RECREATE");
    
    TRandom r;

    CS::Stat stat;
    stat.SetBufferMaxSize(10000);
    unsigned s=stat.Size();
    unsigned fill=0;
    
    for( int i=0; i<100000; i++ )
    {
        float v = r.Gaus(0,1);
    
        stat.Add(v);
        //if( stat.Size()>1 )
        //    printf("i=%d  size=%d   mean=%g   sigma=%g\n",i,stat.Size(),stat.Mean(),stat.Sigma());

        if( s>stat.Size() )
        {
            fill++;
            char buf[22];
            sprintf(buf,"h%d",fill);
            TH1F *h = new TH1F(buf,buf,100,-5,5);

            for( std::set<float>::const_iterator v=stat.GetBuffer().begin(); v!=stat.GetBuffer().end(); v++ )
                h->Fill(*v);
        }

        
        s = stat.Size();
    }

    f.Write();
    
    return 0;
}

#endif
