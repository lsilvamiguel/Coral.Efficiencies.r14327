#include <typeinfo>
#include "ChipHotGeSiCA.h"

using namespace std;

namespace CS {

namespace {

class EventHeaderTrailer
{
  public:

                    EventHeaderTrailer      (uint32 d) {u.d=d;}
                    operator uint32         (void) const {return u.d;}
    uint32          GetEvent                (void) const {return u.s.event;}
    uint8           GetGeoID                (void) const {return u.s.geo;}

    union
    {
      struct
      {
        uint32      event      :20,
                    geo        : 4,
                    zeros      : 6,
                    ht         : 1,
                    zero       : 1;
      } s;
      uint32 d;
    } u;
};

class ChipHeaderTrailer
{
  public:

                    ChipHeaderTrailer       (uint32 d) {u.d=d;}
                    operator uint32         (void) const {return u.d;}
    uint32          GetEvent                (void) const {return u.s.event;}
    uint32          GetPort                 (void) const {return u.s.port;}
    bool            IsHeader                (void) const {return u.s.chan==0;}
    bool            IsTrailer               (void) const {return u.s.chan==7;}

    union
    {
      struct
      {
        uint32      chan       : 3,
                    chip       : 3,
                    xor_       : 1,
                    trig_time  : 9,
                    event      : 6,
                    overflow   : 1,
                    null       : 1,
                    port       : 3,
                    zeros      : 4,
                    one        : 1;
      } s;
      uint32 d;
    } u;
};

class PortHeaderTrailer
{
  public:
                    PortHeaderTrailer       (uint32 d) {u.d=d;}
                    operator uint32         (void) const {return u.d;}
    uint32          GetPort                 (void) const {return u.s.port;}
    uint32          GetEvent                (void) const {return u.s.event;}

    union
    {
      struct
      {
        uint32      event      : 6,
                    zeros      :18,
                    port       : 3,
                    mask       : 5;
      } s;
      uint32 d;
    } u;
};

class F1DataLine
{
  public:
                    F1DataLine              (uint32 d) {u.d=d;}
                    operator uint32         (void) const {return u.d;}
    uint32          GetPort                 (void) const {return u.s.port;}
    uint32          GetData                 (void) const {return u.s.data;}
    uint32          GetChannel              (void) const {return u.s.channel;}
    uint32          IsData                  (void) const {return u.s.data_flag;}

    union
    {
      struct
      {
        uint32      data       :16,
                    channel    : 6,
                    _          : 1,
                    data_flag  : 1,
                    port       : 3,
                    zeros      : 4,
                    one        : 1;
      } s;
      uint32 d;
    } u;
};

class Data
{
  public:
    enum ID {None=0,EventHeader=1,EventTrailer=2,Padding=4,PortHeader=8,
             ChipHeader=16,F1Data=32,PortTrailer=64,ChipTrailer=128,Error1=256,Error2=512};
    enum{ Error1Mask=0x7ff };
    
    
    Data(ID _who,uint32 _expected) : id(_who),expected(_expected) {}
    bool Check(ID q) const { return q&expected; }
    
    static ID identify(uint32 data,Chip &chip,DaqOption &opt);
    
    static const string & GetName(ID);
    
    ID id;
    uint32  expected;
  private:
    static map<ID,std::string> names;
};

map<Data::ID,std::string> Data::names;

const string &Data::GetName(ID id)
{
    static bool init=true;
    if(init)
    {
        init=false;

        names[ID(  0)] = "None";
        names[ID(  1)] = "EventHeader";
        names[ID(  2)] = "EventTrailer";
        names[ID(  4)] = "Padding";
        names[ID(  8)] = "PortHeader";
        names[ID( 16)] = "ChipHeader";
        names[ID( 32)] = "F1Data";
        names[ID( 64)] = "PortTrailer";
        names[ID(128)] = "ChipTrailer";
        names[ID(256)] = "Error1";
        names[ID(512)] = "Error2";
    }
    
    map<ID,std::string>::const_iterator it=names.find(id);
    if( it==names.end() )
    {
        static std::string unknown("Unknown");
        return unknown;
    }
    return
        it->second;
}

bool get_bit(unsigned bit,uint32 number)
{ return 1&(number>>bit);

}

// Get the data identification
Data::ID Data::identify(uint32 data,Chip &chip,DaqOption &opt)
{
    uint32 first_bit = get_bit(31,data);
    if( !first_bit )
    {
        if( (data>>20)==Error1Mask )
            return Error1;

        bool is_header = !get_bit(30,data);

        if( is_header )
            return EventHeader;
        else
            return EventTrailer;
    }
    else
    {
        uint32 mask = (data>>27)&15;
        switch( mask )
        {
            case  0:
            {
                // OK, it may be f1-data or chip header/trailer
                F1DataLine f1(data);
                if( f1.IsData() )
                    return F1Data;

                // No, it is not a F1-data!
                ChipHeaderTrailer ht(data);
                if( ht.IsHeader() )
                    return ChipHeader;
                if( ht.IsTrailer() )
                    return ChipTrailer;
                chip.AddError(DaqError(DaqError::ChipHotGeSiCA_BAD_CHIP,
                        DaqError::SOURCE_ID,chip.GetSourceID()),opt);
                return None;
            }
            case  1:    return PortHeader;
            case  2:    return Padding;
            case  3:    return PortTrailer;
            case  4:
            case  5:            
            case  6:             
            case  7:             
            case  8:
            case  9:
            case 10:    return Error2;
            default:    return None;
        }
    }
    
    throw "Data::identify()";
}

}

////////////////////////////////////////////////////////////////////////////////

void ChipHotGeSiCA::Scan(DaqOption &opt)
{
    // Do not scan twice.
    if( is_scaned )
        return;
    else
        is_scaned=true;

    if( GetSLink().IsFirstEventInRun() )
    {
        printf("ChipHotGeSiCA::Scan(): first event of run is skipped!\n");
        return;
    }

    // Do we have any data?
    if( (GetDataEnd()-GetDataStart())==0 )
        return;

    // These variables describe the data flow order. There is more then one data format.
    // We need to look at first words of data to decide which format was used.
    Data
        start           (Data::None,0),
        error1          (Data::None,0),
        error2          (Data::None,0),
        event_header    (Data::None,0),
        event_trailer   (Data::None,0),
        padding         (Data::None,0),
        port_header     (Data::None,0),
        port_trailer    (Data::None,0),
        chip_header     (Data::None,0),
        chip_trailer    (Data::None,0),
        f1_data         (Data::None,0);

    // There are two data formats: full and reduced.
    // Full format is identified by the sequence EventHeader,Padding.
    // Reduced format is identified by the sequence EventHeader,no-Padding.

    // Check that we have at least two words of data.
    if( (GetDataEnd()-GetDataStart())<2 )
    {
        Exception("ChipHotGeSiCA::Scan(): data are too short: %d words",GetDataEnd()-GetDataStart());
        return;
    }

    if( Data::identify(GetDataStart()[0],*this,opt)!=Data::EventHeader )
    {
        Exception("ChipHotGeSiCA::Scan(): First word is not an event header.");
        return;
    }
    
    if( Data::identify(GetDataStart()[1],*this,opt)==Data::Padding )
    {
        // This is a full format

        start               = Data(Data::None,            Data::EventHeader|Data::Error1),
        error1              = Data(Data::Error1,          0xffffffff),
        error2              = Data(Data::Error2,          0xffffffff),
        event_header        = Data(Data::EventHeader,     Data::Padding),
        event_trailer       = Data(Data::EventTrailer,    Data::EventHeader|Data::Error1),
        padding             = Data(Data::Padding,         Data::Padding|Data::PortHeader|Data::EventTrailer|Data::PortTrailer),
        port_header         = Data(Data::PortHeader,      Data::ChipHeader|Data::F1Data|Data::Error2|Data::PortTrailer|Data::Padding),
        port_trailer        = Data(Data::PortTrailer,     Data::PortHeader|Data::Padding),
        chip_header         = Data(Data::ChipHeader,      Data::F1Data|Data::ChipTrailer|Data::Error2),
        chip_trailer        = Data(Data::ChipTrailer,     Data::PortTrailer|Data::Error2|Data::Padding|Data::ChipHeader),
        f1_data             = Data(Data::F1Data,          Data::F1Data|Data::Error2|Data::ChipTrailer|Data::PortTrailer|Data::Padding);
    }
    else
    {
        // Reduced format

        start               = Data(Data::None,            Data::EventHeader|Data::Error1);
        error1              = Data(Data::Error1,          Data::EventHeader|Data::EventTrailer);
        error2              = Data(Data::Error2,          Data::Error1|Data::Error2|Data::F1Data|Data::EventTrailer);
        event_header        = Data(Data::EventHeader,     Data::Error1|Data::Error2|Data::F1Data|Data::EventTrailer);
        event_trailer       = Data(Data::EventTrailer,    Data::Error1|Data::EventHeader);
        f1_data             = Data(Data::F1Data,          Data::Error1|Data::Error2|Data::F1Data|Data::EventTrailer);
    }

    // Final preparations before scanning the data.

    Data
       *previous        = &start,
       *me              = NULL;
  
    const uint32
        *d_event_header = NULL,
        *d_port_header  = NULL,
        *d_chip_header  = NULL;
  
    Data::ID id(Data::None);

    // Loop over data words!

    for( const uint32 *p=GetDataStart(); p<GetDataEnd(); p++ )
    {
        id=Data::identify(*p,*this,opt);

        if(0)
        {
            printf("id=%3d for ",id);
            printf("%3zu  ",p-GetDataStart());
            bits_print(cout,*p,0,32,"%8x\n");
        }

        switch( id )
        {
            case Data::Error1:
                me = &error1;
                //printf("Error1 is ignored...\n");
                break;
 
            case Data::Error2:
            {
                me = &error2;
                uint32 mask = ((*p)>>27)&15;
                uint32 port = F1DataLine(*p).GetPort();
                
                int geo=-1;
                if( d_event_header!=NULL )
                    geo = EventHeaderTrailer(*d_event_header).GetGeoID();

                switch( mask )
                {
                    case  4:
                        AddError(DaqError(DaqError::ChipHotGeSiCA_UNEXPECTED_DT,
                                 DaqError::SOURCE_ID,GetSourceID(),
                                 DaqError::GEO_ID,geo,
                                 DaqError::PORT,port),opt);
                        break;
                    case  5:            
                        AddError(DaqError(DaqError::ChipHotGeSiCA_WRONG_CHIP,
                                 DaqError::SOURCE_ID,GetSourceID(),
                                 DaqError::GEO_ID,geo,
                                 DaqError::PORT,port),opt);
                        break;
                    case  6:             
                        AddError(DaqError(DaqError::ChipHotGeSiCA_WRONG_EVENT,
                                 DaqError::SOURCE_ID,GetSourceID(),
                                 DaqError::GEO_ID,geo,
                                 DaqError::PORT,port),opt);
                        break;
                    case  7:             
                        AddError(DaqError(DaqError::ChipHotGeSiCA_PORT_TIMEOUT,
                                 DaqError::SOURCE_ID,GetSourceID(),
                                 DaqError::GEO_ID,geo,
                                 DaqError::PORT,port),opt);
                        break;
                    case  8:
                        AddError(DaqError(DaqError::ChipHotGeSiCA_EVENT_PORT_LOST,
                                 DaqError::SOURCE_ID,GetSourceID(),
                                 DaqError::GEO_ID,geo,
                                 DaqError::PORT,port),opt);
                        break;
                    case  9:
                        AddError(DaqError(DaqError::ChipHotGeSiCA_DATA_PORT_LOST,
                                 DaqError::SOURCE_ID,GetSourceID(),
                                 DaqError::GEO_ID,geo,
                                 DaqError::PORT,port),opt);
                        break;
                    case 10:
                        AddError(DaqError(DaqError::ChipHotGeSiCA_SKIP_SPILL,
                                 DaqError::SOURCE_ID,GetSourceID(),
                                 DaqError::GEO_ID,geo,
                                 DaqError::PORT,port),opt);
                        break;
                    default: throw "ChipHotGeSiCA::Scan(): Internal problem!";
                }
                break;
            }

            case Data::EventHeader:
                me = &event_header;
                if( d_event_header!=NULL )
                {
                    Exception("ChipHotGeSiCA::Scan(): event header without a previous event trailer.");
                    id=Data::None;  // Stop scan
                    break;
                }
                d_event_header = p;

                if( GetSLink().GetEventNumber()!=EventHeaderTrailer(*d_event_header).GetEvent() )
                {
                    //Exception("ChipHotGeSiCA::Scan(): event numbers mismatches in Slink<->Data");
                    AddError(DaqError(DaqError::ChipHotGeSiCA_WRONG_EVENT1,
                            DaqError::SOURCE_ID,GetSourceID()),opt);
                }
                break;

            case Data::EventTrailer:
                me = &event_trailer;
                if( d_event_header==NULL )
                {
                    Exception("ChipHotGeSiCA::Scan(): event trailer without an event header!");
                    id=Data::None;  // Stop scan
                }
                else
                    if( EventHeaderTrailer(*p).GetEvent()!=EventHeaderTrailer(*d_event_header).GetEvent() )
                    {
                        //Exception("ChipHotGeSiCA::Scan(): event numbers in event header/trailer mismatches!");
                        AddError(DaqError(DaqError::ChipHotGeSiCA_WRONG_EVENT2,
                                 DaqError::SOURCE_ID,GetSourceID()),opt);
                    }
                d_event_header = NULL;
                break;
            
            case Data::Padding:
                me = &padding;
                break;
            
            case Data::PortHeader:
            {
                me = &port_header;
                if( d_port_header!=NULL )
                {
                    Exception("ChipHotGeSiCA::Scan(): port header without a previous port trailer.");
                    id=Data::None;  // Stop scan
                    break;
                }
                d_port_header = p;
                PortHeaderTrailer port_ht(*d_port_header);
                
                if( d_event_header!=NULL )
                    if( (EventHeaderTrailer(*d_event_header).GetEvent()&63)!=port_ht.GetEvent() )
                    {
                        //Exception("ChipHotGeSiCA::Scan(): event numbers mismatches in event_header<->port_header");
                        AddError(DaqError(DaqError::ChipHotGeSiCA_WRONG_EVENT3,
                                 DaqError::SOURCE_ID,GetSourceID(),
                                 DaqError::PORT,port_ht.GetPort()),opt);
                    }
                break;
            }

            case Data::PortTrailer:
            {
                me = &port_trailer;
                if( d_port_header==NULL )
                {
                    Exception("ChipHotGeSiCA::Scan(): port trailer without a port header!");
                    id=Data::None;  // Stop scan
                    break;
                }
                else
                {
                    PortHeaderTrailer port_ht(*p);
                    
                    if( port_ht.GetPort()!=PortHeaderTrailer(*d_port_header).GetPort() )
                    {
                        //Exception("ChipHotGeSiCA::Scan(): port numbers in port header/trailer mismatches!");
                        AddError(DaqError(DaqError::ChipHotGeSiCA_WRONG_PORT_PORT,
                                 DaqError::SOURCE_ID,GetSourceID(),
                                 DaqError::PORT,port_ht.GetPort()),opt);
                    }
                    if( port_ht.GetEvent()!=PortHeaderTrailer(*d_port_header).GetEvent() )
                    {
                        //Exception("ChipHotGeSiCA::Scan(): event numbers in port header/trailer mismatches!");
                        AddError(DaqError(DaqError::ChipHotGeSiCA_WRONG_EVENT4,
                                 DaqError::SOURCE_ID,GetSourceID(),
                                 DaqError::PORT,port_ht.GetPort()),opt);
                    }
                }
                d_port_header = NULL;
                break;
            }
            
            case Data::F1Data:
                me = &f1_data;
                if( d_event_header==NULL )
                {
                    Exception("ChipHotGeSiCA::Scan(): can not add data without an event header!");
                    id=Data::None;  // Stop scan
                }
                else
                {
                    F1DataLine f1(*p);
                    DataID id(GetSourceID(),EventHeaderTrailer(*d_event_header).GetGeoID(),f1.GetPort(),f1.GetChannel());
                    all_data.push_back( pair<DataID,uint16>(id,f1.GetData()) );
                }
                break;
            
            case Data::ChipHeader:
                me = &chip_header;
                if( d_chip_header!=NULL )
                {
                    Exception("ChipHotGeSiCA::Scan(): chip header without a previous chip trailer.");
                    id=Data::None;  // Stop scan
                    break;
                }
                d_chip_header = p;
                break;

            case Data::ChipTrailer:
                me = &chip_trailer;
                if( d_chip_header==NULL )
                {
                    Exception("ChipHotGeSiCA::Scan(): chip trailer without a chip header!");
                    id=Data::None;  // Stop scan
                    break;
                }
                else
                {
                    ChipHeaderTrailer c(*p);
                    if( c.GetPort()!=ChipHeaderTrailer(*d_chip_header).GetPort() )
                    {
                        //Exception("ChipHotGeSiCA::Scan(): port numbers in chip header/trailer mismatches!");
                        AddError(DaqError(DaqError::ChipHotGeSiCA_WRONG_CHIP_PORT,
                                 DaqError::SOURCE_ID,GetSourceID(),
                                 DaqError::PORT,c.GetPort()),opt);
                    }
                }
                break;
            
            default:
                printf("Unknown word id=%d\n",id);
                id = Data::None;
                break;
        }
        
        if( id==Data::None )
        {
            AddError(DaqError(DaqError::ChipHotGeSiCA_BAD_EVENT,
                     DaqError::SOURCE_ID,GetSourceID()),opt);
            break;
        }
        
        if( me!=NULL && !(id==Data::Error1 || id==Data::Error2) )
        {
            if( !previous->Check(me->id) )
            {
                Exception("ChipHotGeSiCA(): srcID=%d Unexpected %s after %s",
                           GetSourceID(),
                           Data::GetName(me->id).c_str(),Data::GetName(previous->id).c_str());
                AddError(DaqError(DaqError::ChipHotGeSiCA_BAD_EVENT,
                         DaqError::SOURCE_ID,GetSourceID()),opt);
                break;
            }
            previous = me;
            //printf("ME=%d expect=%d\n",me->id,me->expected);
        }
    }
    
    if( id!=Data::None && me!=&event_trailer )
    {
        Exception("ChipHotGeSiCA(): srcID=%d  event does not finish with event_trailer. (id=%d)",GetSourceID(),unsigned(id));
    }
}

////////////////////////////////////////////////////////////////////////////////

ChipHotGeSiCA::Map::Map(const ObjectXML &o)
:   Chip::Map(o),
    geoID(0),
    port(0),
    wire_position(0),
    time_unit(-1)
{
    if( version==0 )
        version=1;

    if( GetName()!="ChipHotGeSiCA" )
        throw Exception("ChipHotGeSiCA::Map::Map(): Internal error.");

    if( GetVersion()!=1  )
        throw Exception("ChipHotGeSiCA::Map::Map(): unknown version %d",GetVersion());

    // If "time_unit" has been set in a mapping file, then use it ....
    if( NULL==GetAttribute("time_unit",time_unit) )
    {
        SetAttribute("time_unit","0.12892312");
        GetAttribute("time_unit",time_unit);
    }

    int chanN;
    string name;
    istringstream s(dec_line.c_str());

    s >> name >> source_id >> geoID >> port >> chanF>>chanS>>chanN>>wireF>>wireL>>wireS;

    if( s.fail() )
        throw Exception("ChipHotGeSiCA::Map::Map(): bad format in line: %s",map_line.c_str());

    s >> wire_position;
    if( s.fail() )
        wire_position=0;

    id=DetID(name);

    chanL=chanF+(chanN-1)*chanS;

    Check();
}

////////////////////////////////////////////////////////////////////////////////

void ChipHotGeSiCA::Map::AddToMaps(Maps &maps,DaqOption &options) const
{
    for( size_t n=0; n<GetChanN(); n++ )
    {
        int wire=GetWireF()+n*GetWireS();
        int chan=GetChanF()+n*GetChanS();
        DataID data_id(GetSourceID(),geoID,port,chan);

        Digit *digit = new Digit(data_id,GetDetID(),wire,wire_position,0,time_unit);

        if( maps.end()!=maps.find(data_id) && !IsMultiDigit() )
        {
            Print();
            Exception("ChipF1::Map::AddToMaps():c: map already exists!").Print();
        }

        maps.insert( pair<DataID,Digit*>(data_id,digit) );
    }

    maps.SetWires( GetDetID(), maps.GetWires(GetDetID())+GetChanN() );
}

////////////////////////////////////////////////////////////////////////////////

void ChipHotGeSiCA::Map::Print(ostream &o,const string &prefix) const
{
  Chip::Map::Print(o,prefix);
  o<<prefix;

//  char s[222];
//  sprintf(s,"ChipF1::Map: mode=%c port=%d chip=%d\n",
//             mode,int(port),int(chip));
//  o << s;
}

////////////////////////////////////////////////////////////////////////////////

void ChipHotGeSiCA::Digit::Print(ostream &o,const string &prefix) const
{
    char s[222];

    DataID d = GetDataID();

    sprintf(s,"%s%8s ampl=%5d  wire=%3d  pos=%d time_unit=%gns time_decoded=%g Data=(%s)\n",
               prefix.c_str(),GetDetID().GetName().c_str(),(int)GetAmplitude(),
               GetChannel(),(int)GetChannelPos(),GetTimeUnit(),GetTimeDecoded(),
               std::string(d).c_str());
  o << s;
}

////////////////////////////////////////////////////////////////////////////////

ChipHotGeSiCA::DataID::operator std::string(void) const
{
    char s[222];
    sprintf(s,"src=%d geo=%d port=%d chan=%d",
            unsigned(u.s.src_id), unsigned(u.s.geo_id), unsigned(u.s.port), unsigned(u.s.chan));
    return s;
}

////////////////////////////////////////////////////////////////////////////////

void ChipHotGeSiCA::Decode(const Maps &maps,Digits &digits_list,DaqOption &opt)
{
    if( !is_scaned )
        Scan(opt);
    
    const size_t size=digits_list.size();
    
    for( size_t i=0; i<all_data.size(); i++ )
    {
        const DataID &id = all_data[i].first;
        //cout << std::string(id) << "\n";
        typedef Maps::const_iterator m_it; // Create a short name for map's iterator
        const pair<m_it,m_it> m_range = maps.equal_range(id); // all maps with the given ID
        for( m_it c=m_range.first; c!=m_range.second; c++ )
        {
            const Digit *digit1 = dynamic_cast<Digit*>(c->second);
            if( digit1==NULL )
                throw Exception("ChipHotGeSiCA::Decode(): Internal error");
            Digit *digit2 = new Digit(*digit1);
            digit2->SetAmplitude(all_data[i].second);
            digits_list.insert(pair<DetID,Digit*>(digit2->GetDetID(),digit2));
        }
    }
    
    if( all_data.size()!=(digits_list.size()-size) )
        Exception("ChipHotGeSiCA::Decode():WW: some maps were not found for srcID=%d",GetSourceID());
}

void ChipHotGeSiCA::Print(ostream &o,const string &prefix) const
{
    throw "ChipHotGeSiCA::Print(): no code yet";
}

}
