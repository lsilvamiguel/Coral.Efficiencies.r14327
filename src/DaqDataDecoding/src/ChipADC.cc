#include <sstream>
#include <unistd.h>
#include <cstdio>

#include "ChipADC.h"

using namespace std;

namespace CS {

////////////////////////////////////////////////////////////////////////////////

void ChipADC::Decode(const Maps &maps,Digits &digits_list,DaqOption &opt)
{
    if( !is_scaned )
        Scan(opt);

    const size_t size=digits_list.size();

    for( size_t i=0; i<all_data.size(); i++ )
    {
        const Data d(all_data[i].second);
        const DataID data_id(GetSourceID(),all_data[i].first.GetGeoID(),d.GetChannel());
        typedef Maps::const_iterator m_it; // Create a short name for map's iterator
        const pair<m_it,m_it> m_range = maps.equal_range(data_id); // all maps with given source ID

        for( m_it c=m_range.first; c!=m_range.second; c++ )
        {
            const Digit *digit1 = dynamic_cast<Digit*>(c->second);
            if( digit1==NULL )
              throw Exception("ChipADC::Decode(): Internal error");
            Digit *digit2 = new Digit(*digit1);
            digit2->SetAmplitude(d.GetData());
            digits_list.insert(pair<DetID,Digit*>(digit2->GetDetID(),digit2));
        }
    }

    if( size>=digits_list.size() )
        Exception("ChipADC::Decode():WW: some maps were not found for srcID=%d",GetSourceID());
}

////////////////////////////////////////////////////////////////////////////////

void ChipADC::Print(ostream &o,const string &prefix) const
{
  Chip::Print(o,prefix);
}

////////////////////////////////////////////////////////////////////////////////

ChipADC::Map::Map(const ObjectXML &o) :
  Chip::Map(o),
  geoID(0),
  channel(0),
  x(0),
  y(0)
{
  if( GetName()!=o.GetName() || GetName()!="ChipADC" )
    throw Exception("ChipADC::Map::Map(): Internal error.");

  istringstream s(dec_line.c_str());

  string name;
  int chanN=-1;

  switch( GetVersion() )
  {
    case 1:
    {
      int det_num;
      s >> det_num >> name >> source_id>>geoID>>chanF>>chanS>>chanN>>wireF>>wireL>>wireS;
      chanL=chanF+(chanN-1)*chanS;
      break;
    }
    
    case 2:
    {
      int det_num;
      s >> det_num >> name >> source_id>>geoID>>channel>>x>>y;
      chanN=1;
      chanF=chanL=channel;
      chanS=1;
      wireF=wireL=x;
      wireS=1;
      break;
    }

    case 3:
    case 4:
    {
      if( GetVersion()==3 )
      {
        int det_n;
        s>>det_n;
      }
      
      if( options=="" )
      {
        s >> name >> source_id >> geoID>>chanF>>chanS>>chanN>>wireF>>wireL>>wireS;
        chanL=chanF+(chanN-1)*chanS;
      }
      else
      if( options=="xy" )
      {
        s >> name >> source_id >> geoID >> channel >> x >> y;
        chanN=1;
        chanF=chanL=channel;
        chanS=1;
        wireF=wireL=x;
        wireS=1;
      }
      else
        throw Exception("ChipADC::Map::Map(): unknown option \"%s\%",options.c_str());

      break;
    }

    default:
      throw Exception("ChipADC::Map::Map(): unknown version %d",GetVersion());
  }

  if( s.fail() )
    throw Exception("ChipADC::Map::Map(): bad format in line: %s",map_line.c_str());

  id=DetID(name);
  Check();
}

////////////////////////////////////////////////////////////////////////////////

void ChipADC::Map::AddToMaps(Maps &maps,DaqOption &options) const
{
  if( GetVersion()<=1 || GetVersion()>=5 )
    throw Exception("ChipADC::Map::AddToMaps(): not supported version %d",GetVersion());

  const DataID data_id(GetSourceID(),geoID,channel);

  if( maps.end()!=maps.find(data_id) && !IsMultiDigit() )
  {
      Print();
      Exception("ChipADC::Map::AddToMaps(): map already exists").Print();
  }

  maps.insert( pair<DataID,Digit*>(data_id,new Digit(data_id,GetDetID(),x,y,0)));
  maps.SetWires( GetDetID(), maps.GetWires(GetDetID())+1 );
}

////////////////////////////////////////////////////////////////////////////////

void ChipADC::Map::Print(ostream &o,const string &prefix) const
{
  Chip::Map::Print(o,prefix);
  o<< prefix << "      ";
  char s[222];
  sprintf(s,"geoID=%3d  channel=%3d   x=%2d  y=%2d\n",geoID,channel,x,y);
  o << s;
}

////////////////////////////////////////////////////////////////////////////////

void ChipADC::Digit::Print(ostream &o,const string &prefix) const
{
  o<<prefix;
  char s[222];
  const DataID &d = reinterpret_cast<const DataID &>(GetDataID());

  sprintf(s,"%-8s xy=(%3d,%3d)  ampl=%5d  Data=(srcID=%d,geoID=%d,chan=%d)\n",
             GetDetID().GetName().c_str(),GetX(),GetY(),GetAmplitude(),
             unsigned(d.u.s.src_id),unsigned(d.u.s.geo_id),unsigned(d.u.s.chan));
  o << s;
}

////////////////////////////////////////////////////////////////////////////////

void ChipADC::Data::Print(ostream &o,const string &prefix) const
{
  o << prefix;
  char s[222];
  sprintf(s,"data=%4d otr=%d ofl=%d channel=%2d zero=%4d ped_subst=%d is_data=%d\n",
             u.s.data, u.s.otr, u.s.ofl, u.s.channel, u.s.zero, u.s.ped_subst, u.s.is_data);
  o << s;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace CS
