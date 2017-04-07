#include <typeinfo>
#include <cstdarg>
#include <cmath>
#include <cassert>
#include <map>
#include <sstream>
#include "ChipGassiplex.h"
#include "DaqEvent.h"
#include "DaqOption.h"
#include "utils.h"

using namespace std;

namespace CS {

template <class T> T sqr(const T &x) {return x*x;}

////////////////////////////////////////////////////////////////////////////////

map<uint32, pair<uint8,uint8> > ChipGassiplex::id_to_xy[ChipGassiplex::NTables];
uint16                          ChipGassiplex::f[ChipGassiplex::NTables][6][8];
bool                            ChipGassiplex::init     = ChipGassiplex::Init();

bool ChipGassiplex::Init(void)
{
  // Version original
  f[0][3][7]=542;  f[0][5][2]=505;  f[0][4][2]=540;
  f[0][4][6]=506;  f[0][3][4]=541;  f[0][0][0]=504;
  f[0][3][6]=470;  f[0][5][3]=433;  f[0][1][2]=468;
  f[0][1][6]=434;  f[0][2][4]=469;  f[0][1][0]=432;
  f[0][2][6]=398;  f[0][0][2]=361;  f[0][5][1]=396;
  f[0][2][5]=362;  f[0][1][4]=397;  f[0][5][0]=360;
  f[0][2][7]=326;  f[0][0][3]=289;  f[0][0][1]=324;
  f[0][3][5]=290;  f[0][4][4]=325;  f[0][4][0]=288;
  f[0][4][7]=254;  f[0][4][3]=217;  f[0][2][2]=252;
  f[0][0][6]=218;  f[0][0][4]=253;  f[0][2][0]=216;
  f[0][5][7]=182;  f[0][1][3]=145;  f[0][3][2]=180;
  f[0][5][6]=146;  f[0][0][5]=181;  f[0][2][1]=144;
  f[0][1][7]=110;  f[0][3][3]= 73;  f[0][1][1]=108;
  f[0][1][5]= 74;  f[0][5][4]=109;  f[0][3][1]= 72;
  f[0][0][7]= 38;  f[0][2][3]=  1;  f[0][4][1]= 36;
  f[0][4][5]=  2;  f[0][5][5]= 37;  f[0][3][0]=  0;

      // Version Novemeber 2001
  f[1][3][7]= 38;  f[1][5][2]=  1;  f[1][4][2]= 36;
  f[1][4][6]=  2;  f[1][3][4]= 37;  f[1][0][0]=  0;
  f[1][3][6]=542;  f[1][5][3]=505;  f[1][1][2]=540;
  f[1][1][6]=506;  f[1][2][4]=541;  f[1][1][0]=504;
  f[1][2][6]=470;  f[1][0][2]=433;  f[1][5][1]=468;
  f[1][2][5]=434;  f[1][1][4]=469;  f[1][5][0]=432;
  f[1][2][7]=398;  f[1][0][3]=361;  f[1][0][1]=396;
  f[1][3][5]=362;  f[1][4][4]=397;  f[1][4][0]=360;
  f[1][4][7]=326;  f[1][4][3]=289;  f[1][2][2]=324;
  f[1][0][6]=290;  f[1][0][4]=325;  f[1][2][0]=288;
  f[1][5][7]=254;  f[1][1][3]=217;  f[1][3][2]=252;
  f[1][5][6]=218;  f[1][0][5]=253;  f[1][2][1]=216;
  f[1][1][7]=182;  f[1][3][3]=145;  f[1][1][1]=180;
  f[1][1][5]=146;  f[1][5][4]=181;  f[1][3][1]=144;
  f[1][0][7]=110;  f[1][2][3]= 73;  f[1][4][1]=108;
  f[1][4][5]= 74;  f[1][5][5]=109;  f[1][3][0]= 72;

  for( unsigned table=0; table<NTables; table++ )
  {
    for( unsigned x=0; x<72; x++ )
      for( unsigned y=0; y<144; y++ )
      {
        unsigned
          //chamber        = x/72 + 4*(y/144),
          BORA_column    = (x/6)%12,
          BORA_line      = y/72,
          BORA,
          connector,
          internal_number;

        if( BORA_line==0 || BORA_line==2 )
        {
          // UP
          BORA            = BORA_column;
          connector       = (y/8)%18;
          internal_number = f[table][x%6][y%8] + connector*4;
        }
        else
        {
          // DOWN
          BORA            = BORA_column+12;
          connector       = 8-(y/8)%9;
          internal_number = f[table][5-x%6][7-y%8] + connector*4;
        }

        //unsigned channel_id = chamber*32768 + BORA*1024 + internal_number;
        int id=BORA*1024 + internal_number;

        if( (id%4)==3 )
          cerr << "ChipGassiplex:init(): bad id: " << id << "\n";

        const map<uint32, pair<uint8,uint8> >::const_iterator it=id_to_xy[table].find(id);
        if( it!=id_to_xy[table].end() )
        {
          fprintf(stderr,"ChipGassiplex:init():  duplication:  id=%d: (%d,%d)  (%d,%d)\n",
                          id, it->second.first,it->second.second, x,y);
        }
        else
          id_to_xy[table][id] = pair<uint8,uint8>(x,y);
      }

    #if 0
    // Checking
    for( size_t i=0; i<576; i++ )
    {
      const map<uint32, pair<uint8,uint8> >::const_iterator it=id_to_xy[table].find(i);
      if( it==id_to_xy[table].end() )
      {
        if( (i%4)!=3 )
          cerr << "ChipGassiplex:init():  unknown channelID: " << i << "\n";
      }
      else
        if( (i%4)==3 )
          cerr << "ChipGassiplex:init():  bad id:  id%%4==3: " << i << "\n";
    }

    if(0)
    {
      for( size_t i=0; i<=253502; i++ )
      {
        const map<uint32, pair<uint8,uint8> >::const_iterator it=id_to_xy[table].find(i);
        int x,y;
        if( it==id_to_xy[table].end() )
          x = y = 512;
        else
        {
          x = it->second.first;
          y = it->second.second;
        }

        //printf(" %d , %d , %d \r\n",i,x,y);
      }
      exit(1);
    }

    #endif
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////

ChipGassiplex::ChipGassiplex(void const * const buf,bool copy_buf,DaqOption &opt,DaqEvent &ev) :
  ChipCol(buf,copy_buf,opt,ev)
{}

////////////////////////////////////////////////////////////////////////////////

void ChipGassiplex::Decode(const Maps &maps,Digits &digits_list,DaqOption &opt)
{
  const DataID data_id=CreateDataID5(GetSourceID(),0,0,0,0);
  typedef Maps::const_iterator m_it; // Create a short name for map's iterator
  const pair<m_it,m_it> m_range = maps.equal_range(data_id); // all maps with given source ID

  // Return if no maps were found.
  if( m_range.first==m_range.second )
  {
    Exception("ChipGassiplex::Decode(): no maps for src=%d",GetSourceID());
    return;
  }
  
  const Digit *digit = dynamic_cast<const Digit*>(m_range.first->second);
  if( digit==NULL )
    throw Exception("ChipGassiplex::Decode(): internal error");
  if( digit->GetDecodingTable()>=NTables )
    throw Exception("ChipGassiplex::Decode(): unknown decoding_table=%d",int(digit->GetDecodingTable()));
  
  if( !is_scaned )
    Scan(opt);

  for( size_t i=0; i<all_data.size(); i++ )
  {
    bool decoded=false;
    const Data d(all_data[i].second);

    const GeoID
      geoID_header (all_data[i].first.GetGeoID()),
      geoID_data   (d.GetChamber(),d.GetBora());

    if( !d.IsData() )
      throw Exception("ChipGassiplex::Decode(): Data word is expected! Internal error.");

    if( d.GetMSB()!=8 )
    {
      AddError(DaqError(DaqError::GASSIPLEX_BAD_MSB,DaqError::SOURCE_ID,GetSourceID(),DaqError::VALUE,d.GetMSB()),opt);
      continue;
    }

    if( d.GetBora()>=24 )
    {
      AddError(DaqError(DaqError::GASSIPLEX_BIG_BORA,
                        DaqError::SOURCE_ID,    GetSourceID(),
                        DaqError::GEO_ID,       all_data[i].first.GetGeoID(),
                        DaqError::BORA,         d.GetBora()),opt);
      continue;
    }

    if( d.GetBoraChannel()>575 || (d.GetBoraChannel()%4)==3 )
    {
      AddError(DaqError(DaqError::GASSIPLEX_BAD_CHAN,
                        DaqError::SOURCE_ID,    GetSourceID(),
                        DaqError::GEO_ID,       all_data[i].first.GetGeoID(),
                        DaqError::BORA,         d.GetBora(),
                        DaqError::CHANNEL,      d.GetBoraChannel()),opt);
      continue;
    }

    // Check that geoID is correct.
    if( geoID_header!=geoID_data )  // Wrong ADC number
    {
      AddError(DaqError(DaqError::GASSIPLEX_BAD_GEOID,
                        DaqError::SOURCE_ID,    GetSourceID(),
                        DaqError::CHAMBER,      geoID_header.GetChamber(),
                        DaqError::CHAMBER,      geoID_data  .GetChamber(),
                        DaqError::BORA,         geoID_header.GetBora(),
                        DaqError::BORA,         geoID_data  .GetBora()),opt);
      continue;
    }

//     const map<uint32, pair<uint8,uint8> >::const_iterator it=id_to_xy[digit->GetDecodingTable()].find(d.GetBoraChannel()+(d.GetBora()<<10));
//     if( it==id_to_xy[digit->GetDecodingTable()].end() )
//     {
//       AddError(DaqError(DaqError::GASSIPLEX_BAD_CHANID,
//                         DaqError::SOURCE_ID,    GetSourceID(),
//                         DaqError::CHAMBER,      d.GetChamber(),
//                         DaqError::BORA,         d.GetBora(),
//                         DaqError::CHANNEL,      d.GetBoraChannel(),
//                         DaqError::CHANNEL_ID,   d.GetChannelID()));
//       continue;
//     }
// 
//     uint16 x=it->second.first, y=it->second.second;
// 
//     if( x>=72 || (geoID_data.IsBoraUp()?(y<72||y>=144):(y>=72)) )
//     {
//       Exception("ChipGassiplex::Decode(): bad xy, srcID=%d chamber=%d bora=%d channel=%d channel_id=%d x=%d y=%d",
//                 (int)GetSourceID(),
//                 (int)d.GetChamber(),
//                 (int)d.GetBora(),
//                 (int)d.GetBoraChannel(),
//                 (int)d.GetChannelID(),
//                 (int)x,y);
//       continue;
//     }
// 
//     if( geoID_data.IsBoraUp() )
//       y-=72;

    if( event->IsBadDataSource(GetSourceID()) )
        continue;

    if( event->IsBadDataSource(GetSourceID(),all_data[i].first.GetGeoID()) )
        continue;
    
    uint8 x,y;
    try
    {
        GetXY(digit->GetDecodingTable(),geoID_data,d.GetBoraChannel(),x,y);
    }
    catch(...)
    {
        // Bad data line, go to the next one.
        continue;
    }

    // The loop on all maps.
    for( m_it c=m_range.first; c!=m_range.second && !decoded; c++ )
    {
      const Digit *digit1 = dynamic_cast<Digit*>(c->second);
      if( digit1==NULL )
        throw Exception("ChipGassiplex::Decode(): Internal error");
      digits_list.insert(pair<DetID,Digit*>(digit1->GetDetID(),
                                            new Digit(digit1->GetDataID(),digit1->GetDetID(),geoID_data.GetCathode(),x,y,d.GetAmplitude())));
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

void ChipGassiplex::DecodeCalibEvent(const Maps &maps,DaqOption &opt,...) const
{
    if(0)
    {
        cout << "ChipGassiplex::DecodeCalibEvent():  " << (GetDataEnd()-buf) << " data words (32bits)\n";
        cout << "SourceID=" << GetSourceID() << "\n";
        for( const uint32 *p=buf; p<GetDataEnd(); p++ )
        {
            printf("%5zu ",p-GetDataStart());
            bits_print(cout,*p,0,32,"%8x\n");
        }
        cout << "=============================================================\n";
    }

    va_list ap;
    va_start(ap,&opt);
    list<CathodeCalib> *cathode_calib = va_arg(ap,list<CathodeCalib>*);
    va_end(ap);

    if( cathode_calib==NULL || (cathode_calib=dynamic_cast<list<CathodeCalib>*>(cathode_calib))==NULL )
        throw Exception("ChipGassiplex::DecodeCalibEvent(): wrong arguments. Arbitrary argument must be list<CathodeCalib>& ");

    const DataID data_id=CreateDataID5(GetSourceID(),0,0,0,0);
    typedef Maps::const_iterator m_it; // Create a short name for map's iterator
    const pair<m_it,m_it> m_range = maps.equal_range(data_id); // all maps with given source ID

    // Return if no maps were found.
    if( m_range.first==m_range.second )
    {
      Exception("ChipGassiplex::DecodeCalibEvent(): no maps for src=%d",GetSourceID());
      return;
    }

    const Digit *digit = dynamic_cast<const Digit*>(m_range.first->second);
    if( digit==NULL )
      throw Exception("ChipGassiplex::DecodeCalibEvent(): internal error");
    if( digit->GetDecodingTable()>=NTables )
      throw Exception("ChipGassiplex::DecodeCalibEvent(): unknown decoding_table=%d",int(digit->GetDecodingTable()));

    
    //  Unknown 20-words offset
    const uint32 *p=GetDataStart()+20;

    while(true)
    {
        int tail_length = GetDataEnd()-p;
        if( tail_length>=CathodeCalib::LONG )
        {
            uint32 length;
//             CathodeCalib(p,length,digit->GetDecodingTable());
//             printf("geoID=%d cathode=%d bora=%d\n",(int)c.GetGeoID(),(int)c.GetGeoID().GetCathode(),(int)c.GetGeoID().GetBora());
//             if( cathode_calib->find(c)!=cathode_calib->end() )
//                 Exception("ChipGassiplex::DecodeCalibEvent(): geoID=%d was found >1 time",int(c.GetGeoID())).Print();
//             size_t n=cathode_calib->size();
            cathode_calib->push_back(CathodeCalib(p,length,digit->GetDecodingTable()));
//             if( cathode_calib->size()!=n+1 )
//                 Exception("ChipGassiplex::DecodeCalibEvent(): failed to add geoID=%d to the calib container",int(c.GetGeoID())).Print();
            p += length;
            continue;
        }

        if( tail_length==0 )
            break;

        Exception("ChipGassiplex::DecodeCalibEvent(): bad data length %u tail=%d",GetLength(),tail_length).Print(cerr);
        break;
    }
}

////////////////////////////////////////////////////////////////////////////////

void ChipGassiplex::Print(ostream &o,const string &prefix) const
{
  Chip::Print(o,prefix);
}

////////////////////////////////////////////////////////////////////////////////

void ChipGassiplex::Digit::Print(ostream &o,const string &prefix) const
{
  o << prefix;
  char s[222];
  sprintf(s,"%-8s cathode=%2d x=%2d y=%2d ampl=%4d\n",
             GetDetID().GetName().c_str(),
             (int)GetCathode(),(int)GetX(),(int)GetY(),(int)GetAmplitude());
  o << s;
}

////////////////////////////////////////////////////////////////////////////////

ChipGassiplex::Data::Data(void)
{
  d.line.ampl    = 0;
  d.line.chan    = 0;
  d.line.bora    = 0;
  d.line.chamber = 0;
  d.line.msb     = 8;
}

////////////////////////////////////////////////////////////////////////////////

ChipGassiplex::Data::Data(uint8 msb,uint8 chamber,uint8 bora,uint16 chan,uint16 ampl)
{
  if( msb>=16 || (msb&8)==0 )
    throw Exception("ChipGassiplex::Data::Data(): bad msb=%u",(unsigned)msb);

  d.line.msb     = msb;
  d.line.chamber = chamber;
  d.line.bora    = bora;
  d.line.chan    = chan;
  d.line.ampl    = ampl;
}

////////////////////////////////////////////////////////////////////////////////

void ChipGassiplex::Data::Print(ostream &o,const string &prefix) const
{
  o << prefix;
  bits_print(o,d.all);
  o << ": ";
  if( !IsData() )
    o << "ChipGassiplex::Data::Print():  this is NOT a data line!\n";
  else
  {
    char s[222];
    sprintf(s,"ampl=%3u BORA_chan=%3u (substr=chan%%4=%d) channelID=%6u  bora=%2u chamber=%1u msb=%d%d%d%db\n",
              (uint32)GetAmplitude(), (uint32)GetBoraChannel(),(uint32)GetChannelID()%4,GetChannelID(),
              (uint32)GetBora(), (uint32)GetChamber(),
              (GetMSB()&8)!=0, (GetMSB()&4)!=0, (GetMSB()&2)!=0, (GetMSB()&1)!=0 );
    o << s;
  }
}

////////////////////////////////////////////////////////////////////////////////

ChipGassiplex::Map::Map(const ObjectXML &o) :
  Chip::Map(o)
{
  if( GetName()!=o.GetName() || GetName()!="ChipGassiplex" )
    throw Exception("ChipGassiplex::Map::Map(): Internal error.");

  if( GetVersion()!=1 && GetVersion()!=2 )
    throw Exception("ChipGassiplex::Map::Map(): unknown version %d",GetVersion());

  // If "decoding_table" has been set in a mapping file, then use it ....
  if( NULL==GetAttribute("decoding_table",decoding_table) )
      decoding_table=0;

  istringstream s(dec_line.c_str());
  string name;

  // Support of the old version
  if( GetVersion()==1 )
    {int det_n; s>>det_n;}

  s >> name >> source_id;

  if( s.fail() )
    throw Exception("ChipGassiplex::Map::Map(): bad format in line: %s",map_line.c_str());

  chanF = wireF = 0;
  chanL = wireL = 72*72*16-1;
  chanS = wireS = 1;

  id=DetID(name);
}

////////////////////////////////////////////////////////////////////////////////

void ChipGassiplex::Map::AddToMaps(Maps &maps,DaqOption &options) const
{
  const DataID data_id=CreateDataID5(GetSourceID(),0,0,0,0);

  if( maps.end()!=maps.find(data_id) )
      Exception("ChipGassiplex::Map::AddToMaps(): map already exists!").Print();

  maps.insert( pair<DataID,Digit*>(data_id,new Digit(data_id,GetDetID(),0,0,0,0,decoding_table)));

  maps.SetWires( GetDetID(), maps.GetWires(GetDetID())+72*72*16 );
}

////////////////////////////////////////////////////////////////////////////////

ChipGassiplex::CathodeCalib::CathodeCalib(const uint32 *buf,uint32 &len,unsigned decoding_table)
: geoID(0,0)
{
    if( decoding_table>=NTables )
        throw Exception("ChipGassiplex::CathodeCalib::CathodeCalib(): too big decoding table number.");


    if( len<SHORT )
        throw Exception("ChipGassiplex::CathodeCalib::CathodeCalib(): too short data length len=%d SHORT=%d!",len,SHORT);

    struct Line1    {uint32 number:20,geoID:8,null:4;};
    const Line1 *line1=(const Line1 *)buf;

    struct LineLast {uint32 number:24,null:7,one:1;};
    const LineLast *burst=(const LineLast *)(buf+10),*event=(const LineLast *)(buf+11);

    struct LineData {uint32 data:8,null:23,one:1;};

    if( line1->null!=0 || burst->null!=0 || burst->one!=1 || event->null!=0 || event->one!=1 )
    {
        bits_print(cout,buf[ 0],0,32," 0 (first): %8x\n");
        bits_print(cout,buf[10],0,32,"10 (event): %8x\n");
        bits_print(cout,buf[11],0,32,"11 (burst): %8x\n");
        throw Exception("ChipGassiplex::CathodeCalib::CathodeCalib(): data structure is bad!");
    }

    //if( line1->number!=((event->number<<12)>>12) )
    //    Exception("ChipGassiplex::CathodeCalib::CathodeCalib():WW: wrong event numbers in header/trailer: %u!=%u",
    //               line1->number,event->number);
    geoID=line1->geoID;
    burst_number = burst->number;
    event_number = event->number;
    
    for( int i=0; i<5; i++ )
    {
        const LineData *data=(const LineData *)(buf+i+1);
        if( data->one!=1 || data->null!=0 )
            Exception("ChipGassiplex::CathodeCalib::CathodeCalib(): bad voltage structure");
        voltage.push_back(data->data);
    }

    for( int i=0; i<4; i++ )
    {
        const LineData *data=(const LineData *)(buf+i+6);
        if( data->one!=1 || data->null!=0 )
            Exception("ChipGassiplex::CathodeCalib::CathodeCalib(): bad temperature structure");
        temperature.push_back(data->data);
    }
    
    // End of the engineering frame
    const uint32 *p=buf+SHORT-1;

    if( ((*p)>>30)==1 )
        len=SHORT;
    else
    {
        // OK, this is the long data format.

        if( len<LONG )
            throw Exception("ChipGassiplex::CathodeCalib::CathodeCalib(): too short data length for the 'long' format. len=%d LONG=%d.",len,LONG);

        // Reading thresholds
        for( int i=0; i<144; i++,p++ )
        {
            const struct Data {uint32 t00:10,t01:10,t10:10,null:1,one:1;} *data=(const Data*)p;
            if( data->null!=0 || data->one!=1 )
                Exception("ChipGassiplex::CathodeCalib::CathodeCalib():WW: bad threshold data");

            const_cast<PixelCalib&>(*pixels.insert(PixelCalib(geoID,(i<<2)+0,decoding_table)).first).SetThreshold(data->t00);
            const_cast<PixelCalib&>(*pixels.insert(PixelCalib(geoID,(i<<2)+1,decoding_table)).first).SetThreshold(data->t01);
            const_cast<PixelCalib&>(*pixels.insert(PixelCalib(geoID,(i<<2)+2,decoding_table)).first).SetThreshold(data->t10);
        }

        // Reading pedestals and noise
        for( int i=0; i<144*3*2; i+=2,p+=2 )
        {
            const struct Data1 {uint32 sum:21,chan:10,one:1;} *data1=(const Data1 *)(p);
            const struct Data2 {uint32 sum_sqr:31,one:1;} *data2=(const Data2 *)(p+1);
            if( data1->one!=1 || data2->one!=1 )
                Exception("ChipGassiplex::CathodeCalib::CathodeCalib():WW: wrong noise/pedestals data");
            try
            {
                set<PixelCalib>::iterator pixel=pixels.find(PixelCalib(0,data1->chan,0));
                assert(pixel!=pixels.end());
                const_cast<PixelCalib*>(&*pixel)->SetSum(data1->sum);
                const_cast<PixelCalib*>(&*pixel)->SetSumSqr(data2->sum_sqr);
            }
            catch(...)
            {}
        }

        len=LONG;
    }
    
    // Check the trailer word
    struct TrailerWord {uint32 event:20,geoID:8,null1:2,one:1,null2:1;};
    const TrailerWord *last_line = (const TrailerWord *)p++;
    if( last_line->null1!=0 || last_line->null2!=0 || last_line->one!=1 )
        Exception("ChipGassiplex::CathodeCalib::CathodeCalib():WW: bad last line");
    if( last_line->geoID!=line1->geoID )
        Exception("ChipGassiplex::CathodeCalib::CathodeCalib():WW: geoID mismatch header!=trailer %u!=%u",
                   (unsigned)line1->geoID,(unsigned)last_line->geoID);
}

////////////////////////////////////////////////////////////////////////////////

ChipGassiplex::PixelCalib::PixelCalib(GeoID geo_id,uint16 chan_id,unsigned decoding_table)
: chanID(chan_id),
  cathode(geo_id.GetCathode()),
  threshold(0),
  sum(0),
  sum_sqr(0)
{
    GetXY(decoding_table,geo_id,chan_id,x,y);
}

////////////////////////////////////////////////////////////////////////////////

float ChipGassiplex::PixelCalib::GetSigma(void) const
{
    return sqrt(sum_sqr/(float)ENTRIES - sqr(sum/(float)ENTRIES));
}

////////////////////////////////////////////////////////////////////////////////

void ChipGassiplex::PixelCalib::Print(const string &prefix) const
{
    printf("%schanID=%5u cathode=%2u x,y=(%2u,%2u) threshold=%4u sum=%7u sum_sqr=%10u sigma=%g\n",
            prefix.c_str(),(unsigned)chanID,(unsigned)cathode,(unsigned)x,(unsigned)y,
            (unsigned)threshold,(unsigned)sum,(unsigned)sum_sqr,GetSigma());
}

////////////////////////////////////////////////////////////////////////////////

void ChipGassiplex::CathodeCalib::Print(const string &prefix) const
{
    printf("%sgeoID=%5d cathode=%2d burst=%u event=%u\n",
           prefix.c_str(),(unsigned)geoID,geoID.GetCathode(),burst_number,event_number);
    for( size_t i=0; i<voltage.size(); i++ )
        printf("%svoltage%zu=%3d\n",prefix.c_str(),i,(unsigned)voltage[i]);
    for( size_t i=0; i<temperature.size(); i++ )
        printf("%stemeprature%zu=%3d\n",prefix.c_str(),i,(unsigned)temperature[i]);
    for( set<PixelCalib>::const_iterator it=pixels.begin(); it!=pixels.end(); it++ )
        it->Print(prefix+"PixelCalib: ");
    printf("\n");
}

////////////////////////////////////////////////////////////////////////////////

void ChipGassiplex::GetXY(unsigned decoding_table,const GeoID &geo_id,uint16 chan_id,uint8 &x,uint8 &y)
{
    if( decoding_table>=NTables )
        throw Exception("ChipGassiplex::GetXY(): too big decoding table number.");

    const map<uint32, pair<uint8,uint8> >::const_iterator it=id_to_xy[decoding_table].find(chan_id+(geo_id.GetBora()<<10));

    if( it==id_to_xy[decoding_table].end() )
        throw Exception("ChipGassiplex::GetXY(): chamber=%u bora=%u chanID=%u",
                         geo_id.GetChamber(),geo_id.GetBora(),chan_id);

    x=it->second.first;
    y=it->second.second;

    if( x>=72 || (geo_id.IsBoraUp()?(y<72||y>=144):(y>=72)) )
        throw Exception("ChipGassiplex::GetXY(): bad xy, chamber=%d bora=%d channel=%d channel_id=%d x=%d y=%d",
                         (int)geo_id.GetChamber(),
                         (int)geo_id.GetBora(),
                         (int)chan_id,
                         x,y);

    if( geo_id.IsBoraUp() )
        y-=72;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace CS
