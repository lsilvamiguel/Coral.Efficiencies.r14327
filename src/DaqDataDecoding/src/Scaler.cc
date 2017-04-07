#include <sstream>

#include "DaqOption.h"
#include "Scaler.h"
#include "utils.h"

using namespace std;

namespace CS {

////////////////////////////////////////////////////////////////////////////////

void Scaler::Decode(const Maps &maps,Digits &digits_list,DaqOption &opt)
{
    if( !is_scaned )
        Scan(opt);

    const size_t size=digits_list.size();

    for( size_t i=0; i<all_data.size(); i++ )
    {
        const size_t frame_size=34;
    
        size_t frame = i/frame_size;
        size_t chan  = i%frame_size;

        uint32 value;
        if( chan<32 )
            value = bits_get_number(all_data[i].second,0,31);
        else if( chan==32 )  // trigger time
            value = bits_get_number(all_data[i  ].second,16,15) +
                   (bits_get_number(all_data[i+1].second,16,15)<<15);
        else // hit_in_scaler
            value = bits_get_number(all_data[i-1].second,0,16) +
                   (bits_get_number(all_data[i  ].second,0,16)<<16);

        const DataID data_id=CreateDataID7(GetSourceID(),(uint8)frame,(uint8)chan,0,0,0,0);

        typedef Maps::const_iterator m_it; // Create a short name for map's iterator
        const pair<m_it,m_it> m_range = maps.equal_range(data_id); // all maps with given source ID
    
        for( m_it c=m_range.first; c!=m_range.second; c++ )
        {
            const Digit *digit1 = dynamic_cast<Digit*>(c->second);
            if( digit1==NULL )
                throw Exception("Scaler::Decode(): Internal error");
            Digit *digit2 = new Digit(*digit1);
            digit2->SetValue(value);
            digits_list.insert(pair<DetID,Digit*>(digit2->GetDetID(),digit2));
        }
    }

    if( size>=digits_list.size() )
        Exception("Scaler::Decode():WW: some maps were not found for srcID=%d",GetSourceID());
}

////////////////////////////////////////////////////////////////////////////////

vector<float> Scaler::Digit::GetNtupleData(void) const
{
    vector<float> v;

    v.push_back(channel);
    v.push_back(value);

    return v;
}

////////////////////////////////////////////////////////////////////////////////

void Scaler::Print(ostream &o,const string &prefix) const
{
    Chip::Print(o,prefix);
}

////////////////////////////////////////////////////////////////////////////////

Scaler::Map::Map(const ObjectXML& o) :
    Chip::Map(o)
{
    if( GetName()!=o.GetName() || GetName()!="Scaler" )
        throw Exception("Scaler::Map::Map(): Internal error.");

    string name;
    istringstream s(dec_line.c_str());

    switch( GetVersion() )
    {
        case 1:
        {
            // Support of the old version
            int det_n;
            s>>det_n;
        }
        case 2:
        {
            s >> name >> source_id >> frame;
            chanF = wireF = 0;
            chanL = wireL = 33;
            chanS = wireS = 1;
            break;
        }
    
        case 3:
        {
            int chanN;
            s >> name >> source_id >> frame >>  chanF>>chanS>>chanN>>wireF>>wireL>>wireS;
            chanL=chanF+(chanN-1)*chanS;
            break;
        }
    
        default:
            throw Exception("Scaler::Map::Map(): unknown version %d",GetVersion());
    }

    if( s.fail() )
        throw Exception("Scaler::Map::Map(): bad format in line: %s",map_line.c_str());

    id=DetID(name);
    Check();
}

////////////////////////////////////////////////////////////////////////////////

void Scaler::Map::AddToMaps(Maps &maps,DaqOption &options) const
{
    size_t size=maps.size();

    // read option "time_in_spill" which contains the list of tbnames which
    // are to be used for time in spill measurement
    string tis;
    if( GetAttribute("time_in_spill", tis) ) {
      istringstream s(tis); string str;
      while (s >> str)
	options.GetTTConfig().tis_tbnames.insert(str);
    }

    for( size_t n=0; n<GetChanN(); n++ )
    {
        const int chan = GetChanF()+n*GetChanS();
        const int wire = GetWireF()+n*GetWireS();
        const DataID data_id=CreateDataID7(GetSourceID(),frame,chan,0,0,0,0);

        if( maps.end()!=maps.find(data_id) && !IsMultiDigit() )
        {
            Print();
            Exception("Scaler::Map::AddToMaps(): map already exists").Print();
        }

        maps.insert( pair<DataID,Digit*>(data_id,new Digit(data_id,GetDetID(),wire,0)) );
    }

    if( (maps.size()-size)!=GetChanN() )
        throw Exception("Scaler::Map::AddToMaps(): internal error in Map version=%d",GetVersion());

    maps.SetWires( GetDetID(), maps.GetWires(GetDetID())+GetChanN() );
}

////////////////////////////////////////////////////////////////////////////////

Scaler::Digit::Digit(const DataID &data_id,const DetID &id,uint32 chan,uint32 _value) :
    Chip::Digit(data_id,id),
    channel(chan),
    value(_value)
{}

////////////////////////////////////////////////////////////////////////////////

void Scaler::Digit::Print(ostream &o,const string &prefix) const
{
    o << prefix;
    char s[222];
    sprintf(s,"%-8s  chan=%2d  value=%d\n",
               GetDetID().GetName().c_str(),GetChannel(),GetValue());
    o << s;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace CS
