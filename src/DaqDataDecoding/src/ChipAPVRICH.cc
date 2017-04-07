#include <cassert>
#include "ChipAPVRICH.h"
#include "DaqOption.h"
#include "DaqEvent.h"

using namespace std;

namespace CS {

////////////////////////////////////////////////////////////////////////////////

ChipAPVRICH::Map::Map(const ObjectXML &o) :
  Chip::Map(o),
  card_rotation(false),
  old_apv_board(false),
  always_cmc(false)
{
    if( GetName()!=o.GetName() || GetName()!="ChipAPVRICH" )
        throw Exception("ChipAPVRICH::Map::Map(): Internal error.");

    if( IsOption("old_apv_board") || IsOption("skip10channels") )
        old_apv_board = true;

    string name;

    istringstream s(dec_line.c_str());

    switch( GetVersion() )
    {
        case 1:
            s >> name >> source_id >> adc_id >> apv_id >> pd_x >> pd_y;
            break;
        case 2:
            s >> name >> source_id >> adc_id >> apv_id >> pd_x >> pd_y >> card_rotation;
            break;
        case 3:
            s >> name >> source_id >> adc_id >> apv_id >> pd_x >> card_rotation;
            if( NULL==strchr("JjSs",card_rotation) )
                throw "ChipAPVRICH::Map::Map(): bad card rotation flag!";
            pd_y = 255;
            break;
        default:
            throw Exception("ChipAPVRICH::Map::Map(): bad version number %d",GetVersion());
    }
    
    id=DetID(name);

    if( s.fail() )
        throw Exception("ChipAPVRICH::Map::Map(): bad format in line: %s",map_line.c_str());

  // Do we always have a commmon mode noise correction word in data?
  always_cmc = IsOption("CMC_always");
}

////////////////////////////////////////////////////////////////////////////////

void ChipAPVRICH::Map::AddToMaps(Maps &maps,DaqOption &options) const
{
    // Starting from version 3 one line of mapping describes 4 apv_id
    int number_apv_in_map_line = GetVersion()<3 ? 1 : 4;
    for( int i=0; i<number_apv_in_map_line; i++ )
    {
        const DataID data_id(GetSourceID(),adc_id,apv_id+i,0);

        if( maps.end()!=maps.find(data_id) && !IsMultiDigit() )
        {
            Print();
            Exception("ChipAPVRICH::Map::AddToMaps(): map already exists!").Print();
        }

        Digit *digit = new Digit(data_id,GetDetID(),pd_x,pd_y);
        digit->SetOldAPVBoard(old_apv_board);
        digit->SetCardRotation(card_rotation);
        maps.insert( pair<DataID,Digit*>(data_id,digit) );

        maps.SetWires( GetDetID(), maps.GetWires(GetDetID())+APV_read_channels );
    }

    if( always_cmc )
      options.GetAlwaysCMC4APV().insert(GetSourceID());
}

////////////////////////////////////////////////////////////////////////////////

void ChipAPVRICH::Map::Print(ostream &o,const string &prefix) const
{
    const string prefix2=prefix+"ChipAPVRICH::Map:  ";
    Chip::Map::Print(o,prefix2);
    o<<prefix2;
    char s[222];
    sprintf(s,"adc_id=%d apv_id=%d pd_x=%d pd_y=%d\n",int(adc_id),int(apv_id),int(pd_x),int(pd_y));
    o << s;
}

////////////////////////////////////////////////////////////////////////////////

ChipAPVRICH::Digit::Digit(const DataID &data_id,const DetID &id, uint16 _pd_x, uint16 _pd_y) :
    ChipAPV::Digit(data_id,id, true, 0,0,0,0,0,0,0,0,0,0),
    pd_x(_pd_x),
    pd_y(_pd_y),
    old_apv_board(false)
{}

////////////////////////////////////////////////////////////////////////////////

const char* ChipAPVRICH::Digit::GetNtupleFormat(void) const
{
    static char format[]="x:y:pdx:pdy:chan:chip:chip_chan:ampl1:ampl2:ampl3:addr1:addr2:addr3:time:cono1:cono2:cono3";
    return format;
}

////////////////////////////////////////////////////////////////////////////////

vector<float> ChipAPVRICH::Digit::GetNtupleData(void) const
{
    vector<float> v;

    if( IsSingleFrame() )
        throw Exception("ChipAPV::Digit::GetNtupleFormat(): Single frame output is not supported yet.");
    else
    {
        v.push_back(GetPixelX());
        v.push_back(GetPixelY());
        v.push_back(GetPDx());
        v.push_back(GetPDy());
        
        v.push_back(GetChannel());
        v.push_back(GetChip());
        v.push_back(GetChipChannel());
        v.push_back(GetAmplitude()[0]);
        v.push_back(GetAmplitude()[1]);
        v.push_back(GetAmplitude()[2]);
        v.push_back(GetAddress()[0]);
        v.push_back(GetAddress()[1]);
        v.push_back(GetAddress()[2]);
        v.push_back(GetTimeTag());
        v.push_back(GetCoNo()[0]);
        v.push_back(GetCoNo()[1]);
        v.push_back(GetCoNo()[2]);
    }
    
    return v;
}

////////////////////////////////////////////////////////////////////////////////

namespace {

struct XY_struct {int x,y;} XY[50];

bool initconndata(void)
{
    int ii = 0;

    ii =  0; XY[ii].x = 3; XY[ii].y = 0;
    ii =  2; XY[ii].x = 3; XY[ii].y = 1;
    ii =  4; XY[ii].x = 4; XY[ii].y = 0;
    ii =  6; XY[ii].x = 5; XY[ii].y = 0;
    ii =  8; XY[ii].x = 4; XY[ii].y = 1;
    ii = 10; XY[ii].x = 3; XY[ii].y = 2;
    ii = 12; XY[ii].x = 5; XY[ii].y = 1;
    ii = 14; XY[ii].x = 4; XY[ii].y = 2;
    ii = 16; XY[ii].x = 5; XY[ii].y = 2;
    ii = 18; XY[ii].x = 5; XY[ii].y = 3;
    ii = 20; XY[ii].x = 4; XY[ii].y = 3;
    ii = 22; XY[ii].x = 3; XY[ii].y = 3;
    ii = 24; XY[ii].x = 3; XY[ii].y = 4;
    ii = 26; XY[ii].x = 4; XY[ii].y = 4;
    ii = 28; XY[ii].x = 5; XY[ii].y = 4;
    ii = 30; XY[ii].x = 5; XY[ii].y = 5;
    ii = 32; XY[ii].x = 4; XY[ii].y = 5;
    ii = 34; XY[ii].x = 5; XY[ii].y = 6;
    ii = 36; XY[ii].x = 3; XY[ii].y = 5;
    ii = 38; XY[ii].x = 4; XY[ii].y = 6;
    ii = 40; XY[ii].x = 5; XY[ii].y = 7;
    ii = 42; XY[ii].x = 4; XY[ii].y = 7;
    ii = 44; XY[ii].x = 3; XY[ii].y = 6;
    ii = 46; XY[ii].x = 3; XY[ii].y = 7;

    for( int jj = 1; jj<48; jj+=2 )
    {
        XY[jj].x = 5-XY[jj-1].x;
        XY[jj].y = XY[jj-1].y;
    }
    
    return true;
}

}

////////////////////////////////////////////////////////////////////////////////

/*! \brief A little bit modified code from Damien
*/
pair<int,int> ChipAPVRICH::ch2xy(uint8 channel,uint8 card_virtual,bool old_apv_board)
{
    static bool init=initconndata();

    if( card_virtual>=4 || channel>=128 )
        throw "ChipAPVRICH::ch2xy(): bad input";

    if( !old_apv_board )
    {
        if( channel>117 )
            return pair<int,int>(-1,-1);
        channel = 117 - channel;
    }

    short chincon = -1, connb = -1;

    switch(card_virtual)
    {
        case 0:
            if( channel< 48 ) { connb = 0; chincon = channel +  0; break; }
            if( channel< 96 ) { connb = 1; chincon = channel - 48; break; }
            if( channel<108 ) { connb = 2; chincon = channel - 96; break; }
            break;
        case 1:
            if( channel< 36 ) { connb = 2; chincon = channel + 12; break; }
            if( channel< 84 ) { connb = 3; chincon = channel - 36; break; }
            if( channel<108 ) { connb = 4; chincon = channel - 84; break; }
            break;
        case 2:
            if( channel< 24 ) { connb = 4; chincon = channel + 24; break; }
            if( channel< 72 ) { connb = 5; chincon = channel - 24; break; }
            if( channel<108 ) { connb = 6; chincon = channel - 72; break; }
            break;
        case 3:
            if( channel< 12 ) { connb = 6; chincon = channel + 36; break; }
            if( channel< 60 ) { connb = 7; chincon = channel - 12; break; }
            if( channel<108 ) { connb = 8; chincon = channel - 60; break; }
            break;
        default:
            throw "ch2xy(): bad input, wrong card_virtual";
    }

    if( chincon>=0 && connb >=0 )
        return  pair<int,int>( XY[chincon].x, XY[chincon].y+connb*8 );

    return  pair<int,int>(-1,-1);
}

////////////////////////////////////////////////////////////////////////////////

void ChipAPVRICH::Decode(const Maps &maps,Digits &digits_list,DaqOption &opt)
{
    if( !is_scaned )
        Scan(opt);

    const size_t size=digits_list.size();

    for( vector<FullData>::const_iterator it=all_data.begin(); it!=all_data.end(); it++  )
    {
        const Data &d = it->GetData();

        bool frames_good=true;

        if( !it->GetAPVHeader().u.s.frame1_good ||
            !it->GetAPVHeader().u.s.frame2_good ||
            !it->GetAPVHeader().u.s.frame3_good )
        {
            AddError(DaqError(DaqError::APV_H_S_E,
                              DaqError::SOURCE_ID,GetSourceID(),
                              DaqError::VALUE,(int)it->GetAPVHeader().u.s.frame1_good,
                              DaqError::VALUE,(int)it->GetAPVHeader().u.s.frame2_good,
                              DaqError::VALUE,(int)it->GetAPVHeader().u.s.frame3_good),opt);
            frames_good=false;
        }

        // Check that we do not ignore this source ID
        if( event->IsBadDataSource(GetSourceID()) )
            continue;

        uint32 chan_id = d.IsSparse() ? d.GetChannelID() : it->GetChannel();
        if( !it->GetADCHeader().IsReodered() )
            chan_id =  32 * (chan_id%4) + 8 * (chan_id/4) - 31 * (chan_id/16);
        assert(chan_id<128);
        const DataID data_id(GetSourceID(),it->GetADCHeader().GetID(),it->GetAPVHeader().GetChipID(),0);

        typedef Maps::const_iterator m_it; // Create a short name for map's iterator
        const pair<m_it,m_it> m_range = maps.equal_range(data_id); // all maps with given data ID

        if( it->GetADCHeader().IsSingleFrame() )
            throw Exception("ChipAPVRICH::Decode(): single frame readout is not implemented yet");
        else
        {
            for( m_it c=m_range.first; c!=m_range.second; c++ )
            {
                const Digit *digit1 = dynamic_cast<Digit*>(c->second);
                if( digit1==NULL )
                    throw Exception("ChipAPVRICH::Decode(): Internal error");

                if( !frames_good && digit1->GetGoodDigit()==2 )
                    break;  // do not create any bad digits

                pair<int,int> xy=ch2xy(chan_id,it->GetAPVHeader().GetChipID()%4,digit1->GetOldAPVBoard());

                // Create a new digit!
                Digit *digit2 = new Digit(*digit1);
                digit2->SetAddress(apvColumn[(uint8)it->GetAPVHeader().u.s.frame1_gaddr],
                                   apvColumn[(uint8)it->GetAPVHeader().u.s.frame2_gaddr],
                                   apvColumn[(uint8)it->GetAPVHeader().u.s.frame3_gaddr]);

                digit2->SetChip(it->GetAPVHeader().GetChipID());
                digit2->SetChipChannel(chan_id);
                
                uint32 a0,a1,a2,
                    frame1 = d.GetFrameRaw1(),
                    frame2 = d.GetFrameRaw2(),
                    frame3 = d.GetFrameRaw3();
                    
                if( d.IsSparse() )
                {
                    a0 = frame1*2;
                    a1 = frame2*2;
                    a2 = frame3;
                    a0 = a0<256 ? a0 : (256 + (a0-256)*2);
                    a1 = a1<256 ? a1 : (256 + (a1-256)*2);
                    a2 = a2<256 ? a2 : (256 + (a2-256)*2);
                }
                else
                {
                    a0 = frame1;
                    a1 = frame2;
                    a2 = frame3;
                }
                digit2->SetAmplitude(a0,a1,a2);
                digit2->SetSparsifed(d.IsSparse());
                
                digit2->SetTimeTag(it->GetADCHeader().GetTimeTag());
                digit2->SetCoNo(it->GetAPVHeader().GetCoNo()[0],it->GetAPVHeader().GetCoNo()[1],it->GetAPVHeader().GetCoNo()[2]);
                if( xy==pair<int,int>(-1,-1) )
                    // It seems it is a noise from an unconnected channels!
                    digit2->SetPixelXY(-1,-1);
                else
                {
                    // Normal digit
                    if( digit2->GetCardRotation()=='s' || digit2->GetCardRotation()=='S' )
                        digit2->SetPixelXY(xy.first+digit1->GetPDx()*6, 71-xy.second );
                    else
                        digit2->SetPixelXY(5-xy.first+digit1->GetPDx()*6,xy.second);
                }

                if( !frames_good )
                {
                    assert(digit2->GetGoodDigit()==1);
                    digit2->SetGoodDigit(0);
                }

                digits_list.insert(pair<DetID,Digit*>(digit2->GetDetID(),digit2));
            }
        }
    }

    if( size>=digits_list.size() )
        Exception("ChipAPVRICH::Decode():WW: some maps were not found for srcID=%d",GetSourceID());
}

////////////////////////////////////////////////////////////////////////////////

void ChipAPVRICH::Print(ostream &o,const string &prefix) const
{
    throw "ChipAPVRICH::Print(): No code.";
}

////////////////////////////////////////////////////////////////////////////////

void ChipAPVRICH::Digit::Print(ostream &o,const string &prefix) const
{
    ChipAPV::Digit::Print(o,prefix);
    o << prefix << "pd_x=" << pd_x << "   pd_y=" << pd_y << "\n";
}

////////////////////////////////////////////////////////////////////////////////

void ChipAPVRICH::Scan(DaqOption &opt)
{
    ChipAPV::Scan(opt);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

} // namespace CS
