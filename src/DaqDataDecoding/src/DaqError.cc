#include <cstdio>
#include <sstream>

#include "DaqError.h"

using namespace std;

/*! \file DaqError.cc
    \author Alexander Zvyagin
*/

namespace CS {

map<DaqError::Type,DaqErrorType> DaqErrorType::error_types;
bool DaqErrorType::init=DaqErrorType::Init(); // DaqErrorType class initialisation

map<DaqError::Arg, string> DaqError::argument_names;

////////////////////////////////////////////////////////////////////////////////

bool DaqErrorType::Init(void)
{
    new DaqErrorType(DaqError::UNKNOWN_TYPE,          "DAQ ERROR: Unknown DaqError",
                     DaqError::WARNING,       DaqError::Action::NOTHING                );

    new DaqErrorType(DaqError::EXCEPTION,             "DECODING ERROR: Exception",
                     DaqError::WARNING,       DaqError::Action::NOTHING                );

    //TCS errors
    new DaqErrorType(DaqError::TCS_FIFO,              "TCS ERROR: FIFO full bit set",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_SRCID );
    new DaqErrorType(DaqError::TCS_ECC,               "TCS ERROR: ECC checksum error bit set",
                     DaqError::OK,            DaqError::Action::NOTHING                );
    new DaqErrorType(DaqError::TCS_SYNC,              "TCS ERROR: synchroniation error bit set",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_SRCID  );
    new DaqErrorType(DaqError::TCS_UNDEF,             "TCS ERROR: Undefined error bit set",
                     DaqError::MINOR_PROBLEM, DaqError::Action::NOTHING                );
    new DaqErrorType(DaqError::SLINK_ERRFLAG,         "ERROR: Error flag set",
                     DaqError::MINOR_PROBLEM, DaqError::Action::NOTHING                );
    new DaqErrorType(DaqError::SLINK_ERRNR,           "ERROR: Error counter !=0",
                     DaqError::MINOR_PROBLEM, DaqError::Action::NOTHING                );
    new DaqErrorType(DaqError::SLINK_WRONG_EVENT_NUMBER,"DAQ ERROR: SLink event number differs from the one of the event header.",
                     DaqError::SEVERE_PROBLEM,       DaqError::Action::DISCARD_EVENT  );

    //Catch-special error words
    new DaqErrorType(DaqError::TDC_SERR1,             "FRONTEND ERROR: Trigger buffer overflow",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_GEOID  );
    new DaqErrorType(DaqError::TDC_SERR2,             "HOTLINK/CABLE ERROR: CMC FIFO full",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_GEOID  );
    new DaqErrorType(DaqError::TDC_SERR3,             "HOTLINK/CABLE ERROR: Transm. error in bit 23..0",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_GEOID  );

    //Catch error words	
    new DaqErrorType(DaqError::TDC_WRGEOID,           "FRONTEND ERROR: Wrong Geographical ID",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_GEOID  );
    new DaqErrorType(DaqError::TDC_ERR1,              "FRONTEND ERROR: Data without header received",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_GEOID  );
    new DaqErrorType(DaqError::TDC_ERR2,              "FRONTEND ERROR: Timeout: no words received",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_GEOID  );
    new DaqErrorType(DaqError::TDC_ERR3,              "FRONTEND ERROR: Header or trailer with smaller event number received",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_GEOID  );
    new DaqErrorType(DaqError::TDC_ERR4,              "CATCH ERROR: No CMC connected",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_GEOID  );
    new DaqErrorType(DaqError::TDC_ERR5,              "FRONTEND ERROR: Timeout: no trailer or too many data words",
                     DaqError::MINOR_PROBLEM, DaqError::Action::NOTHING                );
    new DaqErrorType(DaqError::TDC_ERR6,              "HOTLINK/CABLE ERROR: FIFO full, port off till end of burst",
                     DaqError::WARNING,       DaqError::Action::NOTHING                );
    new DaqErrorType(DaqError::TDC_ERR7,              "FRONTEND ERROR: Header or trailer with larger event number received",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_GEOID  );
    new DaqErrorType(DaqError::TDC_UNDEF8,            "CATCH ERROR: undefined Error# 8",
                     DaqError::MINOR_PROBLEM, DaqError::Action::NOTHING                );
    new DaqErrorType(DaqError::TDC_UNDEF9,            "CATCH ERROR: undefined Error# 9",
                     DaqError::MINOR_PROBLEM, DaqError::Action::NOTHING                );
    new DaqErrorType(DaqError::TDC_ERR10,             "HOTLINK/CABLE ERROR: HOTLink transm. error",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_GEOID  );
    new DaqErrorType(DaqError::TDC_UNDEF11,           "CATCH ERROR: undefined Error# 11",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_GEOID  );
    new DaqErrorType(DaqError::TDC_ERR12,             "CATCH ERROR: event data was skipped",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_GEOID  );
    new DaqErrorType(DaqError::TDC_UNDEF13,           "CATCH ERROR: undefined Error# 13",
                     DaqError::MINOR_PROBLEM, DaqError::Action::NOTHING                );
    new DaqErrorType(DaqError::TDC_UNDEF14,           "CATCH ERROR: undefined Error# 14",
                     DaqError::MINOR_PROBLEM, DaqError::Action::NOTHING                );
    new DaqErrorType(DaqError::TDC_UNDEF15,           "CATCH ERROR: undefined Error# 15",
                     DaqError::MINOR_PROBLEM, DaqError::Action::NOTHING                );  

    // ChipF1 format errors  
    new DaqErrorType(DaqError::TDC_WRPORT_T,          "TDC ERROR: trailer Port# != port# in header",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_SRCID  );
    new DaqErrorType(DaqError::TDC_TBO,               "TDC ERROR: Trigger buffer overflow bit set",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_GEOID  );
    new DaqErrorType(DaqError::TDC_WRPLL_H,           "TDC ERROR: PLL bits not set on TDC-CMC header",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_GEOID  );
    new DaqErrorType(DaqError::TDC_ERR3U,             "TDC ERROR: Header/trailer with smaller event#, without error word",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_GEOID  );
    new DaqErrorType(DaqError::TDC_ERR7U,             "TDC ERROR: Header/trailer with larger event#, without error word",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_GEOID  );
    new DaqErrorType(DaqError::TDC_WRERR,             "TDC ERROR: error word 3/7 found , but following event number was equal",
                     DaqError::MINOR_PROBLEM, DaqError::Action::NOTHING                );
    new DaqErrorType(DaqError::TDC_ERR1U,             "TDC ERROR: Data without header, without preceeding error word",
                     DaqError::MINOR_PROBLEM, DaqError::Action::NOTHING                );
    new DaqErrorType(DaqError::TDC_WRPORT,            "TDC ERROR: Data port# != port# in header",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_SRCID  );
    new DaqErrorType(DaqError::TDC_WRPLL,             "TDC ERROR: PLL bits not set on TDC-CMC data",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_GEOID  );
    new DaqErrorType(DaqError::TDC_WRTIME,            "TDC ERROR: trigger time mismatch",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_GEOID  );

 
    // ChipSinica format errors
    new DaqErrorType(DaqError::TDC_TBO2,              "TDC ERROR: Trigger buffer overflow bit set",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_PORT );

    new DaqErrorType(DaqError::TDC_TEMPERATURE,       "TDC ERROR: FEM Temperature High",
                     DaqError::WARNING,       DaqError::Action::NOTHING);

    new DaqErrorType(DaqError::TDC_ERR20,             "TDC ERROR: Maximum data package size without header",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_PORT );

    new DaqErrorType(DaqError::TDC_ERR21,             "TDC ERROR: Maximum data package size without trailer",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_PORT );

    new DaqErrorType(DaqError::TDC_ERR22,             "TDC ERROR: Maximum data package size without header or trailer",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_PORT );
  
  

    //ADC/RICH/Scaler format errors
    new DaqErrorType(DaqError::ADC_DATA_WO_HEADER,    "ADC ERROR:Data without header",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_GEOID  );
    new DaqErrorType(DaqError::ADC_W_GEOID_EW,        "ADC ERROR:extra-word has wrong GeoID",
                     DaqError::MINOR_PROBLEM, DaqError::Action::NOTHING                );
    new DaqErrorType(DaqError::ADC_NO_TRAILER,        "ADC ERROR:missing Trailer",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_GEOID  );
    new DaqErrorType(DaqError::ADC_TRAILER_WO_HEADER, "ADC ERROR:Trailer without header",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_GEOID  );
    new DaqErrorType(DaqError::ADC_W_GEOID_TRAILER,   "ADC ERROR:Trailer with wrong geoID",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_GEOID  );
    new DaqErrorType(DaqError::ADC_WR_EV_HEADER,      "ADC ERROR: Wrong event number in data header",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_GEOID  );
    new DaqErrorType(DaqError::ADC_WR_EV_TRAILER,     "ADC ERROR: Wrong event number in data trailer",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_GEOID  );

    new DaqErrorType(DaqError::SLINKM_BAD_SIZE,       "SLINK MULTIPLEXER ERROR: Bad data block size",
                     DaqError::WARNING,       DaqError::Action::NOTHING                 );

    new DaqErrorType(DaqError::APV_H_S_E,             "APV header sample error",
                     DaqError::WARNING,       DaqError::Action::NOTHING                 );
    new DaqErrorType(DaqError::APV_ADC_OUT,           "APV ERROR: ADC data are outside of the limits",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_SRCID  );
    new DaqErrorType(DaqError::APV_APV_OUT,           "APV ERROR: APV data are outside of the limits",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_SRCID  );
    new DaqErrorType(DaqError::APV_H_WR_EV,           "APV ERROR: header event number mismatches with SLINK ones",
                     DaqError::WARNING,       DaqError::Action::NOTHING  );
    new DaqErrorType(DaqError::APV_LOCAL_TIMETAG_FAILED,"APV ERROR: ADC time tag differs across a single GeSiCA",
                     DaqError::WARNING,       DaqError::Action::NOTHING  );
    new DaqErrorType(DaqError::APV_LAST_TRIGGER_TICKS_FAILED,"APV ERROR: ADC time tag differs across several GeSiCAs",
                     DaqError::WARNING,       DaqError::Action::NOTHING  );
    new DaqErrorType(DaqError::APV_MISS_ADC,          "APV ERROR: ADC is missing",
                     DaqError::WARNING,       DaqError::Action::NOTHING  );

    // ChipGassiplex errors

    new DaqErrorType(DaqError::GASSIPLEX_BAD_MSB,     "GASSIPLEX ERROR: bad MSB",
                     DaqError::WARNING,       DaqError::Action::NOTHING                 );

    new DaqErrorType(DaqError::GASSIPLEX_BIG_BORA,    "GASSIPLEX ERROR: too big BORA number",
                     DaqError::WARNING,       DaqError::Action::NOTHING                 );

    new DaqErrorType(DaqError::GASSIPLEX_BAD_CHAN,    "GASSIPLEX ERROR: bad BORA channel (>575 or %%4=3)",
                     DaqError::WARNING,       DaqError::Action::NOTHING                 );

    new DaqErrorType(DaqError::GASSIPLEX_BAD_GEOID,   "GASSIPLEX ERROR: header_geoID!=data_geoID",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT_ON_SRCID  );

    new DaqErrorType(DaqError::GASSIPLEX_BAD_CHANID,  "GASSIPLEX ERROR: wrong channel ID",
                     DaqError::WARNING,       DaqError::Action::NOTHING                 );

    // ChipHotGeSiCA errors
    new DaqErrorType(DaqError::ChipHotGeSiCA_WRONG_EVENT1,"ChipHotGeSiCA ERROR: event numbers mismatch in Slink and event header",
                     DaqError::WARNING,       DaqError::Action::NOTHING );
    new DaqErrorType(DaqError::ChipHotGeSiCA_WRONG_EVENT2,"ChipHotGeSiCA ERROR: event numbers mismatch in event header/trailer",
                     DaqError::WARNING,       DaqError::Action::NOTHING );
    new DaqErrorType(DaqError::ChipHotGeSiCA_WRONG_EVENT3,"ChipHotGeSiCA ERROR: event numbers mismatch in event header and port header",
                     DaqError::WARNING,       DaqError::Action::NOTHING );
    new DaqErrorType(DaqError::ChipHotGeSiCA_WRONG_EVENT4,"ChipHotGeSiCA ERROR: event numbers mismatch in port header/trailer",
                     DaqError::WARNING,       DaqError::Action::NOTHING );
    new DaqErrorType(DaqError::ChipHotGeSiCA_WRONG_PORT_PORT,"ChipHotGeSiCA ERROR: port numbers mismatch in port header/trailer",
                     DaqError::WARNING,       DaqError::Action::NOTHING );
    new DaqErrorType(DaqError::ChipHotGeSiCA_WRONG_CHIP_PORT,"ChipHotGeSiCA ERROR: port numbers mismatch in chip header/trailer",
                     DaqError::WARNING,       DaqError::Action::NOTHING );
    new DaqErrorType(DaqError::ChipHotGeSiCA_BAD_EVENT,"ChipHotGeSiCA ERROR: Bad event structure",
                     DaqError::WARNING,       DaqError::Action::NOTHING );
    new DaqErrorType(DaqError::ChipHotGeSiCA_BAD_CHIP,"ChipHotGeSiCA ERROR: Bad chip header/trailer structure.",
                     DaqError::WARNING,       DaqError::Action::NOTHING );
    new DaqErrorType(DaqError::ChipHotGeSiCA_UNEXPECTED_DT,"ChipHotGeSiCA ERROR: unexpected data or trailer.",
                     DaqError::WARNING,       DaqError::Action::NOTHING );
    new DaqErrorType(DaqError::ChipHotGeSiCA_WRONG_CHIP,"ChipHotGeSiCA ERROR: wrong chip recieved.",
                     DaqError::WARNING,       DaqError::Action::NOTHING );
    new DaqErrorType(DaqError::ChipHotGeSiCA_WRONG_EVENT,"ChipHotGeSiCA ERROR: wrong event recieved.",
                     DaqError::WARNING,       DaqError::Action::NOTHING );
    new DaqErrorType(DaqError::ChipHotGeSiCA_PORT_TIMEOUT,"ChipHotGeSiCA ERROR: port timeout recieved",
                     DaqError::WARNING,       DaqError::Action::NOTHING );
    new DaqErrorType(DaqError::ChipHotGeSiCA_EVENT_PORT_LOST,"ChipHotGeSiCA ERROR: event lost on port",
                     DaqError::WARNING,       DaqError::Action::NOTHING );
    new DaqErrorType(DaqError::ChipHotGeSiCA_DATA_PORT_LOST,"ChipHotGeSiCA ERROR: truncation of data on port",
                     DaqError::WARNING,       DaqError::Action::NOTHING );
    new DaqErrorType(DaqError::ChipHotGeSiCA_SKIP_SPILL,"ChipHotGeSiCA ERROR: skip spill",
                     DaqError::WARNING,       DaqError::Action::NOTHING );

    // DAQ event errors

    new DaqErrorType(DaqError::EVENT_CORRUPTED0,      "DAQ ERROR: Event is corrupted!",
                     DaqError::SEVERE_PROBLEM,DaqError::Action::DISCARD_EVENT           );

    new DaqErrorType(DaqError::EVENT_CORRUPTED1,      "DAQ ERROR: equipment length is too small or not 4-bytes aligned",
                     DaqError::SEVERE_PROBLEM,DaqError::Action::DISCARD_EVENT           );

    new DaqErrorType(DaqError::EVENT_CORRUPTED2,      "DAQ ERROR: equipment length != SLink event size",
                     DaqError::SEVERE_PROBLEM,DaqError::Action::DISCARD_EVENT           );

    new DaqErrorType(DaqError::EVENT_SUBEVENT_HEADER_DIFF, "DAQ ERROR: event headers in the event and subevent do not match",
                     DaqError::SEVERE_PROBLEM,DaqError::Action::DISCARD_EVENT           );

    new DaqErrorType(DaqError::EVENT_UNKNOWN,         "DAQ ERROR: unknown event type",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT           );

    new DaqErrorType(DaqError::EVENT_BAD_CHIP_SIZE,   "DAQ ERROR: chip data is out of the event",
                     DaqError::WARNING,       DaqError::Action::DISCARD_EVENT           );

    new DaqErrorType(DaqError::EVENT_MISSING_SRCID,   "DAQ ERROR: missing sourceID(s)",
                     DaqError::WARNING,       DaqError::Action::NOTHING                 );

    new DaqErrorType(DaqError::EVENT_UNKNOWN_SRCID,   "DAQ ERROR: unknown sourceID(s)",
                     DaqError::WARNING,       DaqError::Action::NOTHING                 );

    new DaqErrorType(DaqError::EVENT_WRONG_GEOID_ORDER,"DAQ ERROR: wrong geoIDs/ports order",
                     DaqError::WARNING,       DaqError::Action::NOTHING                 );

    new DaqErrorType(DaqError::EVENT_SRCID_TOO_BIG,    "DAQ ERROR: Too big number for a sourceID",
                     DaqError::SEVERE_PROBLEM,       DaqError::Action::NOTHING  );

    new DaqErrorType(DaqError::SADC_WRONG_EVENT,        "SADC ERROR: bad event number in ADC header",
                     DaqError::WARNING,       DaqError::Action::NOTHING                 );

    new DaqErrorType(DaqError::SADC_BAD_H_ADC,      "SADC ERROR: bad ADC header",
                     DaqError::WARNING,       DaqError::Action::NOTHING                 );

    new DaqErrorType(DaqError::SADC_H_ADC_ERROR,    "SADC ERROR: ADC header error bit is set",
                     DaqError::WARNING,       DaqError::Action::NOTHING                 );

    new DaqErrorType(DaqError::SADC_H_ADC_BAD_SIZE, "SADC ERROR: bad block size in ADC header",
                     DaqError::WARNING,       DaqError::Action::NOTHING                 );

    new DaqErrorType(DaqError::SADC_UNKNOWN_MODE, "SADC ERROR: unknown mode in ADC header",
                     DaqError::WARNING,       DaqError::Action::NOTHING                 );

    new DaqErrorType(DaqError::SADC_BAD_DATA,       "SADC ERROR: bad 'data' word",
                     DaqError::WARNING,       DaqError::Action::NOTHING                 );

    new DaqErrorType(DaqError::SADC_BAD_DATA_HEADER,"SADC ERROR: bad 'data header' word",
                     DaqError::WARNING,       DaqError::Action::NOTHING                 );

    new DaqErrorType(DaqError::SADC_DATA_OUTSIDE,"SADC ERROR: Data words outside data block",
                     DaqError::WARNING,       DaqError::Action::NOTHING                 );

    new DaqErrorType(DaqError::SADC_MODULE_MISSING,"SADC ERROR: Module is missing",
                     DaqError::WARNING,       DaqError::Action::NOTHING                 );

    new DaqErrorType(DaqError::SADC_BAD_INTEGRAL,   "SADC ERROR: bad 'integral' word",
                     DaqError::WARNING,       DaqError::Action::NOTHING                 );

    new DaqErrorType(DaqError::SADC_BAD_INTEGRALS,  "SADC ERROR: integral words are not consistent",
                     DaqError::WARNING,       DaqError::Action::NOTHING                 );

    new DaqErrorType(DaqError::SADC_BAD_CHAN_N,     "SADC ERROR: channel numbers are not increasing",
                     DaqError::WARNING,       DaqError::Action::NOTHING                 );

    new DaqErrorType(DaqError::SADC_NO_DATA_LINE,   "SADC ERROR: data line in ADC block is missing",
                     DaqError::WARNING,       DaqError::Action::NOTHING                 );

    new DaqErrorType(DaqError::SADC_UNKNOWN_DATA,   "SADC ERROR: unknown data line",
                     DaqError::WARNING,       DaqError::Action::NOTHING                 );

    new DaqErrorType(DaqError::SADC_SHORT_SAMPLES,  "SADC ERROR: channel data samples are too small",
                     DaqError::WARNING,       DaqError::Action::NOTHING                 );

    new DaqErrorType(DaqError::SADC_BAD_N_SAMPLES,  "SADC ERROR: number of declared an presented samples does not match",
                     DaqError::WARNING,       DaqError::Action::NOTHING                 );

    new DaqErrorType(DaqError::TT_WRONG_HITS_IN_MASTER_TIME, "MASTERTIME ERROR: number of master time hits # should be # for resolution #",
                     DaqError::WARNING,       DaqError::Action::NOTHING  );

    new DaqErrorType(DaqError::TT_BIG_SIGMA,           "MASTERTIME ERROR: Trigger time has too large sigma",
                     DaqError::WARNING,       DaqError::Action::NOTHING  );

    new DaqErrorType(DaqError::WRONG_HITS_IN_TRIGGER_MASK, "MASTERTIME ERROR: Too small or too large number of hits in the trigger mask",
                     DaqError::SEVERE_PROBLEM, DaqError::Action::NOTHING);

    new DaqErrorType(DaqError::TT_MT_SHIFT_CHANGED, "MASTERTIME ERROR: TriggerTime-MasterTime shift changed for the trigger",
                     DaqError::SEVERE_PROBLEM, DaqError::Action::NOTHING);

    new DaqErrorType(DaqError::TCS_PHASE_NO_DIGITS, "There is no 'TCS phase' digit",
                     DaqError::SEVERE_PROBLEM, DaqError::Action::NOTHING);

    new DaqErrorType(DaqError::TCS_PHASE_MANY_DIGITS, "There is more than one 'TCS phase' digit",
                     DaqError::WARNING, DaqError::Action::NOTHING);

    new DaqErrorType(DaqError::TCS_PHASE_OUT_OF_RANGE, "The computed TCS phase lies outside of +/- 200ns",
                     DaqError::SEVERE_PROBLEM, DaqError::Action::NOTHING);

    new DaqErrorType(DaqError::TIS_NO_DATA, "Time in spill cannot be computed because there is no data from any of the FluxScalers",
                     DaqError::SEVERE_PROBLEM, DaqError::Action::NOTHING);

    return true;
}

////////////////////////////////////////////////////////////////////////////////

DaqErrorType::DaqErrorType(DaqError::Type t,const string &description,DaqError::SeverityLevel l,DaqError::Action::Type a)
: type(t), name(description), level(l), action(a)
{
    if( error_types.count(t)!=0 )
        throw Exception("DaqErrorType::DaqErrorType() the id %d is already in use!",int(type));
    error_types.insert(pair<DaqError::Type,DaqErrorType>(t,*this));
}

////////////////////////////////////////////////////////////////////////////////

DaqErrorType& DaqErrorType::GetDaqErrorType(DaqError::Type t)
{
    map<DaqError::Type,DaqErrorType>::iterator it = error_types.find(t);
    if( it==error_types.end() )
        throw Exception("DaqErrorType::GetDaqErrorTypes() can not find id %d!",int(t));
    return it->second;
}

////////////////////////////////////////////////////////////////////////////////

DaqErrorType& DaqErrorType::GetDaqErrorType(const string &s)
{
    for( map<DaqError::Type,DaqErrorType>::iterator it=error_types.begin();
         it!=error_types.end(); it++ )
        if( it->second.GetName()==s )
            return it->second;

    throw Exception("DaqErrorType::GetDaqErrorTypes() can not find error \"%s\"",s.c_str());
}

////////////////////////////////////////////////////////////////////////////////

DaqError& DaqError::operator = (const DaqError &d)
{
    if( this!=&d )
    {
        type    = d.type;
        args    = d.args;

        delete e;
        if( d.e!=NULL )
            e = new Exception(*d.e);
        else
            e = NULL;
    }

    return *this;
}

////////////////////////////////////////////////////////////////////////////////

DaqError::operator Exception(void) const
{
    if( e!=NULL )
        return *e;

    std::ostringstream format;

    vector<ArgType> v;
    format << DaqErrorType::GetDaqErrorType(type).GetName();

    for( Args::const_iterator a=args.begin(); a!=args.end(); a++ )
    {
        format << " " << GetArgName(a->first) << "=%lg";
        v.push_back(a->second);
    }

    switch( v.size() )
    {
        case  0: return Exception(format.str().c_str());
        case  1: return Exception(format.str().c_str(),v[0]);
        case  2: return Exception(format.str().c_str(),v[0],v[1]);
        case  3: return Exception(format.str().c_str(),v[0],v[1],v[2]);
        case  4: return Exception(format.str().c_str(),v[0],v[1],v[2],v[3]);
        case  5: return Exception(format.str().c_str(),v[0],v[1],v[2],v[3],v[4]);
        case  6: return Exception(format.str().c_str(),v[0],v[1],v[2],v[3],v[4],v[5]);
        case  7: return Exception(format.str().c_str(),v[0],v[1],v[2],v[3],v[4],v[5],v[6]);
        case  8: return Exception(format.str().c_str(),v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7]);
        case  9: return Exception(format.str().c_str(),v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8]);
        case 10: return Exception(format.str().c_str(),v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9]);
        case 11: return Exception(format.str().c_str(),v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10]);
        case 12: return Exception(format.str().c_str(),v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11]);
        case 13: return Exception(format.str().c_str(),v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12]);
        case 14: return Exception(format.str().c_str(),v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13]);
        case 15: return Exception(format.str().c_str(),v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14]);
        default: throw Exception("DaqError::operator Exception(): too many (%d) arguments for the message\"%s\"",
                                  v.size(),format.str().c_str());
    }
}

////////////////////////////////////////////////////////////////////////////////

const string& DaqError::GetArgName(DaqError::Arg a)
{
    map<Arg,string> &id=argument_names; // Create a short name.

    if( id.size()==0 )
    {
        id[SLINK            ] = "SLink";
        id[PORT             ] = "port";
        id[CHIP             ] = "chip";
        id[BORA             ] = "BORA";
        id[CHAMBER          ] = "chamber";
        id[CHANNEL          ] = "channel";
        id[CHANNEL_ID       ] = "channelID";
        id[SOURCE_ID        ] = "sourceID";
        id[GEO_ID           ] = "geoID";
        id[COUNTER          ] = "counter";
        id[EVENT_N          ] = "event";
        id[SL_EVENT_N       ] = "SLinkEvent";
        id[VALUE            ] = "value";
    }

    map<Arg,string>::const_iterator it=id.find(a);
    if( it==id.end() )
        throw Exception("DaqError::GetArgName(): unknown type %d",int(a));
    return it->second;
}

////////////////////////////////////////////////////////////////////////////////

DaqError::ArgType DaqError::GetArg(const Arg &arg,size_t number) const
{
    pair<Args::const_iterator,Args::const_iterator> range = args.equal_range(arg);
  
    for( Args::const_iterator a=range.first; a!=range.second; a++ )
        if( number-- == 0 )
            return a->second;
  
    std::ostringstream s;

    Print(s);
    string ss(s.str(),s.tellp());
    if( ss.length()>0 )
      ss[ss.length()-1]=0;  // remove the trailing "\n"
    throw Exception("DaqError::GetArg(): not found: (arg=%s number=%d) args=%s",
                     GetArgName(arg).c_str(),number,ss.c_str());
}

////////////////////////////////////////////////////////////////////////////////

void DaqError::Print(ostream &o,const string &prefix) const
{
    o << prefix << DaqErrorType::GetDaqErrorType(type).GetName();
    for( Args::const_iterator it=args.begin(); it!=args.end(); it++ )
        o << " " << GetArgName(it->first) << "=" << it->second;
    o << "\n";

    if( e!=NULL )
        e->Print(o,prefix);
}

////////////////////////////////////////////////////////////////////////////////

size_t DaqErrors::Add(const DaqError &e)
{
    const size_t n = errors.size(); // The error number.
    errors.push_back(e);

    const DaqErrorType &t = e.GetDaqErrorType();

    errors_type     [t.GetType()            ].insert(n);
    errors_level    [t.GetSeverityLevel()   ].insert(n);
    errors_action   [t.GetAction()          ].insert(n);

    return errors.size();
}

////////////////////////////////////////////////////////////////////////////////

void DaqErrors::Remove(int nn)
{
    if( nn<0 )
    {
        errors.clear();
        errors_type.clear();
        errors_level.clear();
        errors_action.clear();
    }
    else
    {
        size_t n = size_t(nn);

        if( n>=errors.size() )
            throw Exception("DaqError::Remove(): too big index: n=%d  size=%d",n,errors.size());
        
        errors.erase(errors.begin()+n);
        
        for( map<DaqError::Type,set<size_t> >::iterator it=errors_type.begin();
             it!=errors_type.end(); it++ )
            it->second.erase(n);

        for( multimap<DaqError::SeverityLevel,set<size_t> >::iterator it=errors_level.begin();
             it!=errors_level.end(); it++ )
            it->second.erase(n);

        for( multimap<DaqError::Action::Type,set<size_t> >::iterator it=errors_action.begin();
             it!=errors_action.end(); it++ )
            it->second.erase(n);
    }
}

////////////////////////////////////////////////////////////////////////////////

const set<size_t> & DaqErrors::Get(DaqError::SeverityLevel l) const
{
    static const set<size_t> empty;

    map<DaqError::SeverityLevel,set<size_t> >::const_iterator it = errors_level.find(l);
    if( it==errors_level.end() )
        return empty;

    return it->second;
}

////////////////////////////////////////////////////////////////////////////////

void DaqErrors::Print(DaqError::SeverityLevel l) const
{
    for( set<size_t>::const_iterator e=Get(l).begin(); e!=Get(l).end(); e++ )
        errors[*e].Print();
}

////////////////////////////////////////////////////////////////////////////////

const set<size_t> & DaqErrors::Get(DaqError::Action::Type a) const
{
    static const set<size_t> empty;

    map<DaqError::Action::Type,set<size_t> >::const_iterator it = errors_action.find(a);
    if( it==errors_action.end() )
        return empty;

    return it->second;
}

////////////////////////////////////////////////////////////////////////////////

void DaqErrors::Print(DaqError::Action::Type a) const
{
    for( set<size_t>::const_iterator e=Get(a).begin(); e!=Get(a).end(); e++ )
        errors[*e].Print();
}

////////////////////////////////////////////////////////////////////////////////

const set<size_t> & DaqErrors::Get(DaqError::Type t) const
{
    static const set<size_t> empty;

    map<DaqError::Type,set<size_t> >::const_iterator it = errors_type.find(t);
    if( it==errors_type.end() )
        return empty;

    return it->second;
}

////////////////////////////////////////////////////////////////////////////////

void DaqErrors::Print(DaqError::Type t) const
{
    for( set<size_t>::const_iterator e=Get(t).begin(); e!=Get(t).end(); e++ )
        errors[*e].Print();
}

////////////////////////////////////////////////////////////////////////////////

} // namespace CS
