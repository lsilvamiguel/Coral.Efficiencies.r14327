#include <cstdio>
#include <string>
#include <popt.h>
#include "Chip.h"
#include "ChipF1.h"
#include "ChipAPV.h"
#include "ChipADC.h"
#include "DaqOption.h"
#include "utils.h"

using namespace std;
using namespace CS;

int main(int argc,const char *argv[])
{
    int run=0, chan[]={-1,-1};
    char *det   = strdup("");
    char *schan = strdup("-1,-1");
    char *maps  = strdup("/afs/cern.ch/compass/detector/maps");

    try
    {
        struct poptOption options[] =
        {
            { "maps",       '\0',   POPT_ARG_STRING|POPT_ARGFLAG_SHOW_DEFAULT,
                                    &maps,  0,
                                    "Map file or direcory", "PATH" },
            { "run",        '\0',   POPT_ARG_INT,
                                    &run, 0,
                                    "Run number to be analyuyzed", "RUN" },
            { "det",        '\0',   POPT_ARG_STRING,
                                    &det, 0,
                                    "Detector name(s) (regular expression)", "NAME" },
            { "chan",       '\0',   POPT_ARG_STRING|POPT_ARGFLAG_SHOW_DEFAULT,
                                    &schan, 0,
                                    "Detector channel", "CHANNEL[,POSITION]" },
            POPT_AUTOHELP
            POPT_TABLEEND
        };

        poptContext poptcont=poptGetContext(NULL,argc,argv,options,0);

        poptSetOtherOptionHelp(poptcont,
            "<options...>\n"
            "  Get mapping information.\n"
            "  Author: Alexander Zvyagin <Alexander.Zvyagin@cern.ch>\n"
        );

        if( poptGetNextOpt(poptcont)>=0 )
            throw "Internal problem in options!";

        if( argc<=1 )
        {
            poptPrintHelp(poptcont,stdout,0);
            return 1;
        }


        if( 2!=sscanf(schan,"%d,%d",&chan[0],&chan[1]) )
                        if( 1!=sscanf(schan,"%d",&chan[0]) )
                            throw "Bad channel number.";

        //--------------------------------------------------------------------------


        Chip::Maps maps_of_the_run;
        DaqOption opts_of_the_run;
        vector<string> detectors_all;
        Chip::ReadMaps( run, maps, maps_of_the_run, opts_of_the_run, detectors_all );

        for( Chip::Maps::const_iterator it=maps_of_the_run.begin(); it!=maps_of_the_run.end(); it++ )
        {
            if( !string_match(it->second->GetDetID().GetName(),det) )
                continue;

            const ChipF1::Digit *d_f1=dynamic_cast<const ChipF1::Digit *>(it->second);
            if( NULL!=d_f1 )
            {
                if( (chan[0]==-1 || d_f1->GetX()==chan[0]) &&
                    (chan[1]==-1 || d_f1->GetY()==chan[1]) )
                {
                    d_f1->Print();
                }
            }

            const ChipAPV::Digit *d_apv=dynamic_cast<const ChipAPV::Digit *>(it->second);
            if( NULL!=d_apv )
            {
                if( chan[1]!=-1 )
                    printf("Warning: bad chan=%d,%d for ChipAPV.\n",chan[0],chan[1]);

                if( chan[0]==-1 || d_apv->GetChannel()==chan[0] )
                {
                    d_apv->Print();
                }
            }

            const ChipADC::Digit *d_adc=dynamic_cast<const ChipADC::Digit *>(it->second);
            if( NULL!=d_adc )
            {
                if( (chan[0]==-1 || d_adc->GetX()==chan[0]) &&
                    (chan[1]==-1 || d_adc->GetY()==chan[1]) )
                {
                    d_adc->Print();
                }
            }

        }

        return 0;
    }
    catch(const char *e)
    {
        printf("%s\n",e);
    }
    catch(const std::exception &e)
    {
        printf("%s\n",e.what());
    }
    catch(...)
    {
        printf("Unknown exception.\n");
    }

    return 1;
}

////////////////////////////////////////////////////////////////////////////////
