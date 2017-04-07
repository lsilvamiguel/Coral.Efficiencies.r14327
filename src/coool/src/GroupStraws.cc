#include <cstdio>
#include "TThread.h"
#include "GroupStraws.h"
#include "PlaneStrawTubes.h"
#include "DaqEvent.h"

ClassImp(GroupStraws);

using namespace std;
using namespace CS;

static const unsigned int straw_channels_a_XY[2]={ 96, 64};
static const unsigned int straw_channels_b_XY[2]={222,192};

GroupStraws::GroupStraws(const char* name) :
    Group(name)
{
}

void GroupStraws::Init(void)
{
    Group::Init();
    fRateCounter = 0;
    
    fRates = new TH1F_Ref((fName+"_rates").c_str(),(fName+" rates distribution").c_str(),1128,0,1128,fRateCounter);
    fRates -> GetXaxis() -> SetTitle("channels");
    fRates -> GetXaxis() -> SetTitle("rates per channel (kHz)");
    AddHistogram(fRates);

    fChannels = new TH1F_Ref((fName+"_ch").c_str(),(fName+" channels distribution").c_str(),1128,0,1128,fRateCounter);
    fChannels -> GetXaxis() -> SetTitle("channels");
    fChannels -> GetXaxis() -> SetTitle("Hits counter");
    AddHistogram(fChannels);

    fHits = new TH1F_Ref((fName+"_hits").c_str(),(fName+" hits multiplicity").c_str(),20,0,20,fRateCounter);
    AddHistogram(fHits);

    fTime = new TH1F_Ref((fName+"_t").c_str(),(fName+" drift time").c_str(),400,-14000,-9000,fRateCounter);
    AddHistogram(fTime);

    fTime6mm = new TH1F((fName+"_t_6mm").c_str(),(fName+" drift time for 6mm straws").c_str(),400,-14000,-9000);
    AddHistogram(fTime6mm);
    
    fTime10mm = new TH1F((fName+"_t_10mm").c_str(),(fName+" drift time for 10mm straws").c_str(),400,-14000,-9000);
    AddHistogram(fTime10mm);

#if USE_DATABASE == 1
  setDBpt(fDataBase);
#endif
  OpenReference();
  if (fReferenceDirectory) {
    ((TH1F_Ref*)fHits)->SetReference(fReferenceDirectory);
    ((TH1F_Ref*)fTime)->SetReference(fReferenceDirectory);
    ((TH1F_Ref*)fChannels)->SetReference(fReferenceDirectory);
  }

}

void GroupStraws::EndEvent(const CS::DaqEvent &event)
{
    // fRateCounter++; removed, done by Monitor
    Group::EndEvent();
    
    if(thr_flag)
        TThread::Lock();
    
    unsigned n_hits=0;
    float window=0;
    
    for( vector<const Plane*>::const_iterator pp=fPlanes.begin(); pp!=fPlanes.end(); pp++ )
    {
        const PlaneStrawTubes *p=dynamic_cast<const PlaneStrawTubes *>(*pp);
        if( p==NULL )
        {
            printf("GroupStraws::EndEvent(): Not PlaneStrawTubes found!\n");
            continue;
        }

        float w = ( p->GetOffTrigTVarID().GetMax() - p->GetOffTrigTVarID().GetMin() );
        if( window==0 )
            window = w;
        else
        if( window!=w )
            printf("GroupStraws::EndEvent(): warning: different F1 time window!  %g!=%g\n",window,w);
        
        //printf("%s %d digits!\n",p->GetName(),p->GetDigits().size());
        for( list<CS::Chip::Digit*>::const_iterator it=p->GetDigits().begin(); it!=p->GetDigits().end(); it++ )
        {
            const CS::ChipF1::Digit *d=dynamic_cast<const CS::ChipF1::Digit *>(*it);
            if( d==NULL )
            {
                printf("GroupStraws::EndEvent(): bad digits was found!\n");
                continue;
            }

            const string &s = d->GetDetID().GetName();
            int channel=-1;

            const bool ud = s[6]=='d';
            if      (s[7]=='a') channel = 2*d->GetChannel()+ud;
            else if (s[7]=='b') channel = 2*(straw_channels_a_XY[s[4]=='Y']+d->GetChannel())+ud;
            else if (s[7]=='c') channel = 2*(straw_channels_a_XY[fName[4]=='Y']+straw_channels_b_XY[s[4]=='Y']+d->GetChannel())+ud;

            assert(channel!=-1);
            fChannels->Fill(channel);
            fTime->Fill(d->GetTimeDecoded()/d->GetTimeUnit());
            
            if( s[7]=='b' )
                fTime6mm->Fill(d->GetTimeDecoded()/d->GetTimeUnit());
            else
                fTime10mm->Fill(d->GetTimeDecoded()/d->GetTimeUnit());

            if( d->GetChannelPos()!=0 )
            {
                // Fill physical hole region in another histogram place
                // 80 and 95 are startuing channel numbers of the physical hole
                assert(s[7]=='b');
                int chan=d->GetChannel()-(s[4]=='Y'?80:95);
                if( chan<0 || chan>=32 )
                    printf("Bad channel: for %s %d %d\n",s.c_str(),d->GetChannel(),chan);
                else
                    fChannels->Fill(1000 + (d->GetChannelPos()<0?0:64) + chan*2 + ud);
            }

            if( p->GetTimeVariable().Test(d->GetTimeDecoded()/d->GetTimeUnit()) )
                n_hits ++;
        }
        
        // fill the rates histogram
        const float timewin = window * Plane1V::fF1_TICK * fRateCounter;
        if( timewin > 0 )
            for( int i=0; i<=fChannels->GetXaxis()->GetNbins()+1; i++ )
                fRates->SetBinContent(i,fChannels->GetBinContent(i)/timewin);
    }
    
    fHits->Fill(n_hits);

    if(thr_flag)
        TThread::UnLock();
}
