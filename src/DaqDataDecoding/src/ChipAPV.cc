#include <cassert>
#include <cstdio>
#include <limits>

#include "ChipAPV.h"
#include "ChipF1.h"
#include "TriggerTime.h"
#include "DaqEvent.h"
#include "DaqOption.h"
#include "utils.h"

using namespace std;

namespace CS {

int ChipAPV::apvColumn[256];
bool ChipAPV::init=ChipAPV::InitApvColumn();


////////////////////////////////////////////////////////////////////////////////

void ChipAPV::Scan(DaqOption &opt)
{
  if( is_scaned )
    return;
  else
    is_scaned=true;

  const bool always_cmc = opt.GetAlwaysCMC4APV().count(GetSourceID());

  uint32 lastTrigTick(std::numeric_limits<uint32>::max());

  // save all ADC IDs found in the SLink block to check if all chips at least sent a header
  set<uint32> adcIds;

  for( const uint32* d=GetDataStart(); d<GetDataEnd(); d++ )
  {
    const ADCHeader h_adc = *(ADCHeader*)d;

// char strtmp[200];
// if (GetSourceID()==380) {
// sprintf(strtmp, "%s %d ", GetName().c_str(), GetSourceID());
// h_adc.Print(cerr, strtmp);
// }


    const uint32 *end_adc = d+h_adc.GetSize();

    if( end_adc<=d || end_adc>GetDataEnd() ) {
      AddError(DaqError(DaqError::APV_ADC_OUT,DaqError::SOURCE_ID,GetSourceID()),opt);
      return;
    }

    // this comparison is safe as the time tag is only 24 bits
    if ( opt.GetLocalTriggerTicksSrcIgnore().end()==opt.GetLocalTriggerTicksSrcIgnore().find(GetSourceID()) ) {
        if (lastTrigTick==std::numeric_limits<uint32>::max())
            lastTrigTick = h_adc.GetTimeTag();
        else if (lastTrigTick!=h_adc.GetTimeTag())
                AddError(DaqError(DaqError::APV_LOCAL_TIMETAG_FAILED,DaqError::SOURCE_ID,GetSourceID(),DaqError::PORT,h_adc.GetID()),opt);
    }

    if ( opt.GetLastTriggerTicksSrcIgnore().end()==opt.GetLastTriggerTicksSrcIgnore().find(GetSourceID()) )
        if ( !event->GetTT().SetLastTriggerTicks(h_adc.GetTimeTag()) )
            AddError(DaqError(DaqError::APV_LAST_TRIGGER_TICKS_FAILED,DaqError::SOURCE_ID,GetSourceID(),DaqError::PORT,h_adc.GetID()),opt);

    adcIds.insert(h_adc.GetID());

    for( d+=sizeof(ADCHeader)/sizeof(int32); d<end_adc; d++ ) {
      APVHeader h_apv(d);

      if( h_apv.u.s.event_number!=(0xfff&GetSLink().GetEventNumber()) ) {
          AddError(DaqError(DaqError::APV_H_WR_EV,DaqError::SOURCE_ID,GetSourceID()),opt);
          return;
      }

      const uint32 *end_apv = d+h_apv.GetBlockSize();
      const bool cmc = always_cmc || h_adc.IsComModCor();   // true if cmc is present in data!

      if( cmc )
        h_apv.cono=CommonModeCorrection(*(--end_apv)); // Last word is for common mode correction

      if( end_apv<=d || end_apv>end_adc )
      {
        AddError(DaqError(DaqError::APV_APV_OUT,DaqError::SOURCE_ID,GetSourceID()),opt);
        return;
      }

      uint16 channel=0;
      for( d+=h_apv.GetSize(); d<end_apv; d++ )
        if( !event->IsBadDataSource(GetSourceID()) )
            all_data.push_back(FullData(h_adc,h_apv,Data(*d,h_adc.GetSparseVer()),channel++));

      if( !cmc )
        d--;
    }
    d--;
  }

  map<uint16, set<uint32> >::const_iterator adcsInSrcId;
  if ((adcsInSrcId=opt.GetSrcIDPresentPorts().find(GetSourceID()))!=opt.GetSrcIDPresentPorts().end())
    if (adcIds!=adcsInSrcId->second)
      // find the missing ADCs
      for (set<uint32>::iterator it=adcsInSrcId->second.begin(); it!=adcsInSrcId->second.end(); it++)
        if (adcIds.count(*it)==0)
          AddError(DaqError(DaqError::APV_MISS_ADC,DaqError::SOURCE_ID,GetSourceID(),
                                                   DaqError::PORT,*it),opt);
}

////////////////////////////////////////////////////////////////////////////////

void ChipAPV::Decode(const Maps &maps,Digits &digits_list,DaqOption &opt)
{
    if( !is_scaned )
        Scan(opt);

    const size_t size=digits_list.size();

    for( vector<FullData>::const_iterator it=all_data.begin(); it!=all_data.end(); it++  )
    {
        const Data &d = it->GetData();

// char strtmp[200];
// if (GetSourceID()==380) {
// sprintf(strtmp, "%s %d ", GetName().c_str(), GetSourceID());
// d.Print(cerr, strtmp);
// }

        if( !it->GetAPVHeader().u.s.frame1_good ||
            !it->GetAPVHeader().u.s.frame2_good ||
            !it->GetAPVHeader().u.s.frame3_good )
        {
            AddError(DaqError(DaqError::APV_H_S_E,
                              DaqError::SOURCE_ID,GetSourceID(),
                              DaqError::PORT,it->GetADCHeader().GetID(),
                              DaqError::CHIP,it->GetAPVHeader().GetChipID(),
                              DaqError::VALUE,(int)it->GetAPVHeader().u.s.frame1_good,
                              DaqError::VALUE,(int)it->GetAPVHeader().u.s.frame2_good,
                              DaqError::VALUE,(int)it->GetAPVHeader().u.s.frame3_good),opt);
            continue;
        }

        // Check that we do not ignore this source ID
        if( event->IsBadDataSource(GetSourceID()) )
            continue;

        uint32 chan_id = d.IsSparse() ? d.GetChannelID() : it->GetChannel();
        chan_id =  32 * (chan_id%4) + 8 * (chan_id/4) - 31 * (chan_id/16);
        const DataID data_id(GetSourceID(),it->GetADCHeader().GetID(),it->GetAPVHeader().GetChipID(),chan_id);

        typedef Maps::const_iterator m_it; // Create a short name for map's iterator
        const pair<m_it,m_it> m_range = maps.equal_range(data_id); // all maps with given data ID

        if( it->GetADCHeader().IsSingleFrame() )
            throw Exception("ChipAPV::Decode(): single frame readout is not implemented yet");
        else
        {
            // iterate over all maps matching the data ID
            for( m_it c=m_range.first; c!=m_range.second; c++ )
            {
                const Digit *digit1 = dynamic_cast<Digit*>(c->second);
                if( digit1==NULL )
                    throw Exception("ChipAPV::Decode(): Internal error");
                const DigitPixel *digit1p = dynamic_cast<const DigitPixel*>(digit1);
                const DigitPixelMM *digit1pmm = dynamic_cast<const DigitPixelMM*>(digit1);

                const float af = digit1->GetAmplFactor();
                if ( digit1p==NULL && digit1pmm==NULL ) {
                    Digit *digit2 = new Digit(*digit1);
                    digit2->SetAddress(apvColumn[(uint8)it->GetAPVHeader().u.s.frame1_gaddr],
                    apvColumn[(uint8)it->GetAPVHeader().u.s.frame2_gaddr],
                    apvColumn[(uint8)it->GetAPVHeader().u.s.frame3_gaddr]);

                    // RICH ADC type amplitude decoding (for strip part of PixelGEMs)
                    if ( digit2->GetADCType()==1 )
                    {
                        if( d.IsSparse() )
                            digit2->SetAmplitude((uint32)(af*d.GetFrameNew1()),(uint32)(af*d.GetFrameNew2()),
                                                 (uint32)(af*d.GetFrameNew3()));
                        else
                            digit2->SetAmplitude(d.GetFrameRaw1(),d.GetFrameRaw2(),d.GetFrameRaw3());
                    }

                    // GEM/Si ADC amplitude decoding
                    else
                    {
                        digit2->SetAmplitude(d.GetFrame1(),d.GetFrame2(),d.GetFrame3());
                    }
                    digit2->SetSparsifed(d.IsSparse());

                    digit2->SetTimeTag(it->GetADCHeader().GetTimeTag());
                    digit2->SetCoNo(it->GetAPVHeader().GetCoNo()[0],it->GetAPVHeader().GetCoNo()[1],it->GetAPVHeader().GetCoNo()[2]);

                    digits_list.insert(pair<DetID,Digit*>(digit2->GetDetID(),digit2));
                }
                else if ( digit1p != 0 )
                {
                    // Pixel GEM digit

                    if ( digit1p->GetADCType()!=1 )
                        throw Exception("ChipAPV::Decode(Pixel): only RICH-type ADC used.");

                    // Get pixel coordinates
                    int chan = digit1p->GetChannel();
                    pair<int,int> xy=pixgem::detch2xy(chan);
                    int pixX = xy.first;
                    int pixY = xy.second;

                    // Only create good digits
                    if (pixX!=-1 && pixY!=-1)
                    {
                        // Invert X pads if detector orientation negative (det facing upstream)
                        if( digit1p->GetDetOrientation() < 0 )
                        {
                            pixX = 31 - pixX;
                        }

                        // Create a new digit!
                        DigitPixel *digit2 = new DigitPixel(*digit1p);
                        digit2->SetAddress(apvColumn[(uint8)it->GetAPVHeader().u.s.frame1_gaddr],
                                           apvColumn[(uint8)it->GetAPVHeader().u.s.frame2_gaddr],
                                           apvColumn[(uint8)it->GetAPVHeader().u.s.frame3_gaddr]);

                        digit2->SetChip(it->GetAPVHeader().GetChipID());
                        digit2->SetChipChannel(chan_id);

                        if( d.IsSparse() )
                        {
                            digit2->SetAmplitude(d.GetFrameNew1(),d.GetFrameNew2(),d.GetFrameNew3());
                        }
                        else
                        {
                            digit2->SetAmplitude(d.GetFrameRaw1(),d.GetFrameRaw2(),d.GetFrameRaw3());
                        }
                        digit2->SetSparsifed(d.IsSparse());

                        digit2->SetTimeTag(it->GetADCHeader().GetTimeTag());
                        digit2->SetCoNo(it->GetAPVHeader().GetCoNo()[0],
                                        it->GetAPVHeader().GetCoNo()[1],
                                        it->GetAPVHeader().GetCoNo()[2]);

                        digit2->SetPixelXY(pixX,pixY);

                        digits_list.insert(pair<DetID,Digit*>(digit2->GetDetID(),digit2));
                    }
                }
                else if ( digit1pmm != 0 ) {
                    // Pixel MM digit

                    if ( digit1pmm->GetADCType()!=1 )
                        throw Exception("ChipAPV::Decode(Pixel): only RICH-type ADC used.");

                    // Get pixel coordinates
                    int chan = digit1pmm->GetChannel();
                    int conn = digit1pmm->GetConnNb();
                    int pixnb = pixmm::Instance(digit1pmm->GetPixelMMversion())->GetPixNb(conn,chan);

// if (((ChipAPV::DataID)digit1pmm->GetDataID()).u.s.adc_id==2) digit1pmm->Print(cerr,"digit1pmm: ");
// if (((ChipAPV::DataID)digit1pmm->GetDataID()).u.s.adc_id==2) cerr<<" GetPixelMMversion "<<digit1pmm->GetPixelMMversion()<<" pixnb "<<pixnb<<endl;
                    // Only create good digits
                    if (pixnb>=0) {
                        float pixX, pixY;
                        pixX = pixmm::Instance(digit1pmm->GetPixelMMversion())->GetXPix(pixnb);
                        pixY = pixmm::Instance(digit1pmm->GetPixelMMversion())->GetYPix(pixnb);
                        // Invert X pads if detector orientation negative (det facing upstream)
                        if ( digit1pmm->GetDetOrientation() < 0 ) {
                            pixX = -pixX;
                        }

                        // Create a new digit!
                        DigitPixelMM *digit2 = new DigitPixelMM(*digit1pmm);
                        digit2->SetAddress(apvColumn[(uint8)it->GetAPVHeader().u.s.frame1_gaddr],
                                           apvColumn[(uint8)it->GetAPVHeader().u.s.frame2_gaddr],
                                           apvColumn[(uint8)it->GetAPVHeader().u.s.frame3_gaddr]);

                        digit2->SetChip(it->GetAPVHeader().GetChipID());
                        digit2->SetChipChannel(chan_id);

                        if ( d.IsSparse() ) {
                            digit2->SetAmplitude((uint32)(af*d.GetFrameNew1()),(uint32)(af*d.GetFrameNew2()),
                                                 (uint32)(af*d.GetFrameNew3()));
                        } else {
                            digit2->SetAmplitude(d.GetFrameRaw1(),d.GetFrameRaw2(),d.GetFrameRaw3());
                        }
                        digit2->SetSparsifed(d.IsSparse());

                        digit2->SetTimeTag(it->GetADCHeader().GetTimeTag());
                        digit2->SetCoNo(it->GetAPVHeader().GetCoNo()[0],
                                        it->GetAPVHeader().GetCoNo()[1],
                                        it->GetAPVHeader().GetCoNo()[2]);

                        digit2->SetPixelXYNb(pixX,pixY,pixnb);
                        digit2->SetLargePixel(pixmm::Instance(digit1pmm->GetPixelMMversion())->GetLargePix(pixnb));

                        digits_list.insert(pair<DetID,Digit*>(digit2->GetDetID(),digit2));
// if (digit2 && ((ChipAPV::DataID)digit1pmm->GetDataID()).u.s.adc_id==2) digit2->Print(cerr,"digit2: ");
                    }

                }
            }
        }
    }

    if( all_data.size()>0 && size>=digits_list.size() )
        Exception("ChipAPV::Decode():WW: some maps were not found for srcID=%d",GetSourceID());
}

////////////////////////////////////////////////////////////////////////////////

void ChipAPV::Print(ostream &o,const string &prefix) const
{
  Chip::Print(o,prefix);
}

////////////////////////////////////////////////////////////////////////////////

void ChipAPV::ADCHeader::Print(ostream &o,const string &prefix) const
{
  char s[222];
  sprintf(s,"%sformat=0x%x adc_id=%d len=%d time_tag=%d\n",
             prefix.c_str(),GetFormat(),GetID(),GetSize(),GetTimeTag());
  o << s;
}

////////////////////////////////////////////////////////////////////////////////

void ChipAPV::APVHeader::Print(ostream &o,const string &prefix) const
{
  char s[222];
  sprintf(s,"%sdata_len=%d event_number=%d chip_id=%d\n",
             prefix.c_str(),u.s.data_len,u.s.event_number,u.s.chip_id);
  o << s;
}

////////////////////////////////////////////////////////////////////////////////

void ChipAPV::FullData::Print(ostream &o,const string &prefix) const
{
    adc .Print(o,prefix+"ADC:  ");
    apv .Print(o,prefix+"APV:  ");
    data.Print(o,prefix+"DATA: ");
    o << prefix << "CHAN: " << channel << "\n";
}

////////////////////////////////////////////////////////////////////////////////

void ChipAPV::Data::Print(ostream &o,const string &prefix) const
{
  char s[222];
  sprintf(s,"%sIsSparse=%d channel_id=%d frame1=%d frame2=%d frame3=%d\n",
          prefix.c_str(),IsSparse(),
          GetChannelID(),
          GetFrame1(), GetFrame2(), GetFrame3() );
  o << s;
  sprintf(s,"         new method         frame1=%d frame2=%d frame3=%d\n",
          GetFrameNew1(), GetFrameNew2(), GetFrameNew3() );
  o << s;
}

////////////////////////////////////////////////////////////////////////////////

uint16 ChipAPV::Data::GetFrameNew1(void) const
{
  uint16 a0,
    frame1 = d.s.frame1;

  // 1st/2nd sample 8 bit, least significant bit ignored
  // All samples: two slopes
  a0 = frame1*2;
  a0 = a0<256 ? a0 : (256 + (a0-256)*2);
  return a0;
}
////////////////////////////////////////////////////////////////////////////////

uint16 ChipAPV::Data::GetFrameNew2(void) const
{
  uint16 a1,
    frame2 = d.s.frame2;

  // 1st/2nd sample 8 bit, least significant bit ignored
  // All samples: two slopes
  a1 = frame2*2;
  a1 = a1<256 ? a1 : (256 + (a1-256)*2);
  return a1;
}
////////////////////////////////////////////////////////////////////////////////

uint16 ChipAPV::Data::GetFrameNew3(void) const
{
  uint16 a2,
    frame3 = d.s.frame3;

  // 1st/2nd sample 8 bit, least significant bit ignored
  // All samples: two slopes
  a2 = frame3;
  a2 = a2<256 ? a2 : (256 + (a2-256)*2);
  return a2;
}

////////////////////////////////////////////////////////////////////////////////

ChipAPV::Map::Map(const ObjectXML &o) :
    Chip::Map(o),
    adc_id(0),
    chip_id(0),
    det_orientation(0),
    connNb(-1),
    map_version(0),
    wire_position(0),
    amplfact(1.),
    only_good_digits(true),
    pixel_map(false),
    pixelmm_map(false),
    pixelmm2s_map(false),
    pixelmm3_map(false),
    rich_adc(false),
    always_cmc(false)
{
    map_version = GetVersion();

    // Map for pixel decoding
    pixel_map = IsOption("pixel_gem");
    pixelmm_map = IsOption("pixel_mm");
    pixelmm2s_map = IsOption("pixel_mm_v2s");
    pixelmm3_map = IsOption("pixel_mm_v3");

    // RICH ADC (different amplitude encoding than GEM/Si ADC)
    rich_adc = IsOption("rich_adc");

    // Do we always have a commmon mode noise correction word in data?
    always_cmc = IsOption("CMC_always");

    if( pixel_map )
    {
        if( GetName()!=o.GetName() || GetName()!="ChipAPV" )
            throw Exception("ChipAPV::Map::Map(Pixel): Internal error.");

        string name;
        int chanN;

        istringstream s(dec_line.c_str());

        s>>name>>source_id>>adc_id>>chip_id>>chanF>>chanS>>chanN>>wireF>>wireL>>wireS>>det_orientation;

        id=DetID(name);
        chanL=chanF+(chanN-1)*chanS;

        if( s.fail() )
            throw Exception("ChipAPV::Map::Map(Pixel): bad format in line: %s",map_line.c_str());

        if( chanN<=0 )
          throw Exception("ChipAPV::Map::Map(Pixel): det=\"%s\" chanN must be >0",name.c_str());

        if( GetChanN()!=(unsigned)chanN )
          throw Exception("ChipAPV::Map::Map(Pixel): bad chanN: %s",map_line.c_str());

    }
    else if ( pixelmm_map || pixelmm2s_map || pixelmm3_map )
    {
        if ( GetName()!=o.GetName() || GetName()!="ChipAPV" )
            throw Exception("ChipAPV::Map::Map(Pixel): Internal error.");

        if ( GetVersion()!=3 && GetVersion()!=4 )
          throw Exception("ChipAPV::Map::Map(Strip): unknown version number %d, version=3 or 4 needed with pixel_mm option",GetVersion());

        string name;
        int chanN;

        istringstream s(dec_line.c_str());

        s>>name>>source_id>>adc_id>>chip_id>>chanF>>chanS>>chanN>>connNb>>wireF>>wireL>>wireS>>det_orientation;

        id=DetID(name);
        chanL=chanF+(chanN-1)*chanS;
        if (pixelmm_map) connNb--; // numbering from 1 to 10 in mapping files

        if ( s.fail() )
          throw Exception("ChipAPV::Map::Map(Pixel): bad format in line: %s",map_line.c_str());

        if ( chanN<=0 )
          throw Exception("ChipAPV::Map::Map(Pixel): det=\"%s\" chanN must be >0",name.c_str());

        if ( GetChanN()!=(unsigned)chanN )
          throw Exception("ChipAPV::Map::Map(Pixel): bad chanN: %s",map_line.c_str());

        if( GetVersion()==4 ) {
          s>>amplfact;
          if (s.fail()) amplfact=1.;
        }

    }
    else
    {
        // Read map for strip part

        if( GetName()!=o.GetName() || GetName()!="ChipAPV" )
          throw Exception("ChipAPV::Map::Map(Strip): Internal error.");

        if( GetVersion()!=1 && GetVersion()!=2 && GetVersion()!=3 )
          throw Exception("ChipAPV::Map::Map(Strip): unknown version number %d",GetVersion());

        string name;
        int chanN;

        istringstream s(dec_line.c_str());

        // Support of the old version
        if( GetVersion()==1 )
          { int det_n; s>>det_n;}

        s>>name>>source_id>>adc_id>>chip_id>>chanF>>chanS>>chanN>>wireF>>wireL>>wireS;
        id=DetID(name);
        chanL=chanF+(chanN-1)*chanS;

        if( s.fail() )
          throw Exception("ChipAPV::Map::Map(Strip): bad format in line: %s",map_line.c_str());

        if( chanN<=0 )
          throw Exception("ChipAPV::Map::Map(Strip): det=\"%s\" chanN must be >0",name.c_str());

        if( GetChanN()!=(unsigned)chanN )
          throw Exception("ChipAPV::Map::Map(Strip): bad chanN: %s",map_line.c_str());

        s >> wire_position;
        if( s.fail() )
          wire_position=0;

        if( GetVersion()==3 ) {
          s>>amplfact;
          if (s.fail()) amplfact=1.;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void ChipAPV::Map::AddToMaps(Maps &maps,DaqOption &options) const
{
    const size_t size = maps.size();

    if ( pixel_map )
    {
        for (size_t n=0; n<GetChanN(); n++)
        {
            int chan = GetChanF()+n*GetChanS();
            int wire = GetWireF()+n*GetWireS();
            const DataID data_id(GetSourceID(),adc_id,chip_id,chan);

            if( maps.end()!=maps.find(data_id) && !IsMultiDigit() )
                throw Exception("ChipAPV::Map::AddToMaps(Pixel): map already exists!");

            DigitPixel *digit = new DigitPixel(data_id,GetDetID(),true,chip_id,chan,wire,
                                     0,0,0,  0,0,0, 0,det_orientation,(int8)rich_adc);

            if( IsOption("bad_data_accept") )
            digit->SetGoodDigit(1);

            digit->SetMapVersion(map_version);

            maps.insert( pair<DataID,DigitPixel*>(data_id,digit) );
        }

        if( (maps.size()-size)!=GetChanN() )
            throw Exception("ChipAPV::Map::AddToMaps(Pixel): internal error in Map version=%d",GetVersion());

        maps.SetWires( GetDetID(), maps.GetWires(GetDetID())+GetChanN() );
    }
    else if ( pixelmm_map || pixelmm2s_map || pixelmm3_map )
    {
        register int pixelversion = 1;
        if (pixelmm2s_map) pixelversion = 2;
        if (pixelmm3_map) pixelversion = 3;
        for (size_t n=0; n<GetChanN(); n++) {
            int chan = GetChanF()+n*GetChanS();
            int wire = GetWireF()+n*GetWireS();
            const DataID data_id(GetSourceID(),adc_id,chip_id,chan);

            if( maps.end()!=maps.find(data_id) && !IsMultiDigit() )
                throw Exception("ChipAPV::Map::AddToMaps(PixelMM): map already exists!");

            DigitPixelMM *digit = new DigitPixelMM(data_id,GetDetID(),true,chip_id,chan,wire,
                                                   0,0,0,  0,0,0, 0,det_orientation,(int8)rich_adc,
                                                   connNb, amplfact, pixelversion);

            if( IsOption("bad_data_accept") )
              digit->SetGoodDigit(1);

            digit->SetMapVersion(map_version);

            maps.insert( pair<DataID,DigitPixelMM*>(data_id,digit) );
        }

        if( (maps.size()-size)!=GetChanN() )
            throw Exception("ChipAPV::Map::AddToMaps(Pixel): internal error in Map version=%d",GetVersion());

//         int maxchannel = (connNb+1)*APV_READ_CHANNELS;
//         if (pixelmm2s_map) maxchannel = 768;  // simplified version of pixel scheme
        int maxchannel = pixmm::Instance(pixelversion)->GetNbConn()*APV_READ_CHANNELS;

        if (maxchannel > maps.GetWires(GetDetID())) { maps.SetWires( GetDetID(), maxchannel); }
//        maps.SetWires( GetDetID(), maps.GetWires(GetDetID())+GetChanN() );
    }
    else
    {
        register short stripConn = -1;  // virtual connector number to identify APV chip in pixelMM strips for calibrations
        if ( wire_position == -1 ) {  // short strips, hemisphere -1
          if ( GetWireF() > 255 && GetWireF() < 384 && GetWireL() > 255 && GetWireL() < 384 ) stripConn = 1;
          if ( GetWireF() > 383 && GetWireF() < 512 && GetWireL() > 383 && GetWireL() < 512 ) stripConn = 2;
          if ( GetWireF() > 511 && GetWireF() < 640 && GetWireL() > 511 && GetWireL() < 640 ) stripConn = 3;
        }
        if ( wire_position ==  1 ) {  // short strips, hemisphere +1
          if ( GetWireF() > 255 && GetWireF() < 384 && GetWireL() > 255 && GetWireL() < 384 ) stripConn = 8;
          if ( GetWireF() > 383 && GetWireF() < 512 && GetWireL() > 383 && GetWireL() < 512 ) stripConn = 7;
          if ( GetWireF() > 511 && GetWireF() < 640 && GetWireL() > 511 && GetWireL() < 640 ) stripConn = 6;
        }
        if ( wire_position ==  0 ) {  // long strips
          if ( abs(GetWireS()) == 1 ) {  // new design with non-interleaved long strips
            if ( GetWireF() >=  0 && GetWireF() < 128 && GetWireL() >=  0 && GetWireL() < 128 ) stripConn = 0;
            if ( GetWireF() > 127 && GetWireF() < 256 && GetWireL() > 127 && GetWireL() < 256 ) stripConn = 9;
            if ( GetWireF() > 639 && GetWireF() < 768 && GetWireL() > 639 && GetWireL() < 768 ) stripConn = 4;
            if ( GetWireF() > 767 && GetWireF() < 896 && GetWireL() > 767 && GetWireL() < 896 ) stripConn = 5;
          }
          if ( abs(GetWireS()) == 2 ) {  // old design with interleaved long strips
            if ( GetWireF() == 1 || GetWireL() == 1 ) stripConn = 0;
            if ( GetWireF() == 0 || GetWireL() == 0 ) stripConn = 9;
            if ( GetWireF() == 641 || GetWireL() == 641 ) stripConn = 4;
            if ( GetWireF() == 640 || GetWireL() == 640 ) stripConn = 5;
          }
        }
        for( size_t n=0; n<GetChanN(); n++ ) {
            int chan = GetChanF()+n*GetChanS();
            int wire = GetWireF()+n*GetWireS();
            const DataID data_id(GetSourceID(),adc_id,chip_id,chan);

            if( maps.end()!=maps.find(data_id) && !IsMultiDigit() )
                Exception("ChipAPV::Map::AddToMaps(Strip): map already exists!").Print();

            Digit *digit = new Digit(data_id,GetDetID(),true,chip_id,chan,wire,
                                     0,0,0,  0,0,0, 0,wire_position,(int8)rich_adc, amplfact);

            if( IsOption("bad_data_accept") )
                digit->SetGoodDigit(1);
            digit->SetStripConnNb(stripConn);

            maps.insert( pair<DataID,Digit*>(data_id,digit) );
        }

        if( (maps.size()-size)!=GetChanN() )
            throw Exception("ChipAPV::Map::AddToMaps(Strip): internal error in Map version=%d",GetVersion());

        maps.SetWires( GetDetID(), maps.GetWires(GetDetID())+GetChanN() );
        if( always_cmc )
            options.GetAlwaysCMC4APV().insert(GetSourceID());
    }

    options.AddSrcIDPresentPorts(GetSourceID(), adc_id);
}

////////////////////////////////////////////////////////////////////////////////

void ChipAPV::Map::Print(ostream &o,const string &prefix) const
{
  const string prefix2=prefix+"ChipAPV::Map:  ";
  Chip::Map::Print(o,prefix2);
  o<<prefix2;
  char s[222];
  sprintf(s,"adc_id=%d chip_id=%d det_orientation=%d connNb=%d wire_position=%d\n",
            int(adc_id), int(chip_id), det_orientation, connNb, wire_position);
  o << s;
}

////////////////////////////////////////////////////////////////////////////////

const char* ChipAPV::Digit::GetNtupleFormat(void) const
{
  static char format[]="chan:chip:chip_chan:ampl1:ampl2:ampl3:addr1:addr2:addr3:time:cono1:cono2:cono3:chan_pos";
  if( IsSingleFrame() )
    throw Exception("ChipAPV::Digit::GetNtupleFormat(): Single frame output is not supported yet.");
  else
    return format;
}

////////////////////////////////////////////////////////////////////////////////

ChipAPV::Digit::Digit(const DataID &data_id,const DetID &id,bool sparsifed,
                      uint16 _chip,uint16 _chip_chan,
                      uint16 chan,
                      uint32 addr1,uint32 addr2,uint32 addr3,
                      uint32 ampl1,uint32 ampl2,uint32 ampl3,
                      uint32 _time_tag,int8 _chan_pos,int8 _adc_type, float _amplfactor) :
  Chip::Digit(data_id,id),
  is_single_frame(false),
  is_sparsifed(sparsifed),
  chip(_chip),
  chip_channel(_chip_chan),
  det_channel(chan),
  chan_position(_chan_pos),
  adc_type(_adc_type),
  amplfactor(_amplfactor),
  time_tag(_time_tag),
  good_digit(2),
  strip_conn_number(-1)
{
  address[0]=addr1;
  address[1]=addr2;
  address[2]=addr3;

  amplitude[0]=ampl1;
  amplitude[1]=ampl2;
  amplitude[2]=ampl3;

  SetCoNo(0,0,0);
}

////////////////////////////////////////////////////////////////////////////////

vector<float> ChipAPV::Digit::GetNtupleData(void) const
{
  vector<float> v;

  if( IsSingleFrame() )
    throw Exception("ChipAPV::Digit::GetNtupleFormat(): Single frame output is not supported yet.");
  else
  {
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
    v.push_back(GetChanPos());
  }

  return v;
}

////////////////////////////////////////////////////////////////////////////////

void ChipAPV::Digit::Print(ostream &o,const string &prefix) const
{
  if( IsSingleFrame() )
    throw Exception("ChipAPV::Digit::Print(): single frame readout is not implemented yet");
  else
  {
    char s[300];
    const DataID &d = reinterpret_cast<const DataID &>(GetDataID());

    snprintf(s,300,"%s%-8s  chan=%3d chip=%d chip_chan=%3d [addr,ampl]=[%3u,%3u][%3u,%3u][%3u,%3u]   time_tag=%d  CoNo=(%u,%u,%u) Data=(srcID=%d,adcID=%d,apvID=%d,chan=%d)\n",
              prefix.c_str(),GetDetID().GetName().c_str(),
              GetChannel(), GetChip(), GetChipChannel(),
              GetAddress()[0],GetAmplitude()[0],
              GetAddress()[1],GetAmplitude()[1],
              GetAddress()[2],GetAmplitude()[2],
              GetTimeTag(),GetCoNo()[0],GetCoNo()[1],GetCoNo()[2],
              unsigned(d.u.s.src_id),unsigned(d.u.s.adc_id),unsigned(d.u.s.chip_id),unsigned(d.u.s.chan));
    o << s;
  }
}

////////////////////////////////////////////////////////////////////////////////

ChipAPV::DigitPixel::DigitPixel(const DataID &data_id,const DetID &id, bool sparsified,
                                uint16 chip, uint16 chip_chan,
                                uint16 det_chan,
                                uint32 addr1,uint32 addr2,uint32 addr3,
                                uint32 ampl1,uint32 ampl2,uint32 ampl3,
                                uint32 _time_tag,
                                int8 _det_orientation,
                                int8 _adc_type) :
  ChipAPV::Digit(data_id,id, sparsified,chip,chip_chan,det_chan,addr1,addr2,addr3,ampl1,ampl2,ampl3,_time_tag,0,_adc_type),
  det_orientation(_det_orientation)
{}

////////////////////////////////////////////////////////////////////////////////

const char* ChipAPV::DigitPixel::GetNtupleFormat(void) const
{
    static char format[]="x:y:chan:chip:chip_chan:ampl1:ampl2:ampl3:addr1:addr2:addr3:time:cono1:cono2:cono3";
    return format;
}

////////////////////////////////////////////////////////////////////////////////

vector<float> ChipAPV::DigitPixel::GetNtupleData(void) const
{
    vector<float> v;

    if( IsSingleFrame() )
        throw Exception("ChipAPV::DigitPixel::GetNtupleFormat(): Single frame output is not supported yet.");
    else
    {
        v.push_back(GetPixelX());
        v.push_back(GetPixelY());
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

void ChipAPV::DigitPixel::Print(ostream &o,const string &prefix) const
{
    ChipAPV::Digit::Print(o,prefix);
    o << prefix << " det_orientation=" << int(det_orientation) << "\n";
}

////////////////////////////////////////////////////////////////////////////////

ChipAPV::DigitPixelMM::DigitPixelMM(const DataID &data_id,const DetID &id, bool sparsified,
                                    uint16 chip, uint16 chip_chan,
                                    uint16 det_chan,
                                    uint32 addr1,uint32 addr2,uint32 addr3,
                                    uint32 ampl1,uint32 ampl2,uint32 ampl3,
                                    uint32 _time_tag,
                                    int8 _det_orientation,
                                    int8 _adc_type, int8 _conn_nb, float _amplfactor, int _pixelmm_version) :
  ChipAPV::Digit(data_id,id, sparsified,chip,chip_chan,det_chan,addr1,addr2,addr3,ampl1,ampl2,ampl3,_time_tag,0,_adc_type, _amplfactor),
  det_orientation(_det_orientation),
  conn_number(_conn_nb),
  pixelmm_version(_pixelmm_version)
{
  if (!pixmm::Instance(pixelmm_version)->dataloaded) pixmm::Instance(pixelmm_version)->initconndata();
}

////////////////////////////////////////////////////////////////////////////////

const char* ChipAPV::DigitPixelMM::GetNtupleFormat(void) const
{
    static char format[]="x:y:chan:chip:chip_chan:ampl1:ampl2:ampl3:addr1:addr2:addr3:time:cono1:cono2:cono3:pixnb";
    return format;
}

////////////////////////////////////////////////////////////////////////////////

vector<float> ChipAPV::DigitPixelMM::GetNtupleData(void) const
{
    vector<float> v;

    if( IsSingleFrame() )
        throw Exception("ChipAPV::DigitPixel::GetNtupleFormat(): Single frame output is not supported yet.");
    else
    {
        v.push_back(GetPixelX());
        v.push_back(GetPixelY());
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
        v.push_back(GetPixelNb());
//         v.push_back(GetConnNb());
        v.push_back(pixmm::Instance(pixelmm_version)->GetNConn(GetConnNb()));
    }

    return v;
}

////////////////////////////////////////////////////////////////////////////////

void ChipAPV::DigitPixelMM::Print(ostream &o,const string &prefix) const
{
    ChipAPV::Digit::Print(o,prefix);
    o << prefix << " det_orientation=" << int(det_orientation) << "\n";
    o << prefix << " pixelmm_x=" << pixelmm_x << " pixelmm_y=" << pixelmm_y;
    o << " pixelmm_nb=" << pixelmm_nb<< " pixelmm_largepix=" << pixelmm_largepix;
    o << " pixelmm_version "<< pixelmm_version <<endl;
}

////////////////////////////////////////////////////////////////////////////////

// This is function from Bernhard Ketzer <Bernhard.Ketzer@cern.ch>
// File gemUtils.cxx

//-----------------------------------------------------------------------------
// Initialize array for conversion of APV column code to column address
//-----------------------------------------------------------------------------
bool ChipAPV::InitApvColumn(void)
{
  int A[8], B[8];

  for (int i=0; i<256; i++) apvColumn[i]=-1;

  for (int i=0; i<192; i++) {
    for (int j=0; j<8; j++) A[j] = ( (i+32) >> j ) & 0x1;
    for (int j=0; j<7; j++) B[j] = (A[j+1]?0:1)*A[j] + (A[j]?0:1)*A[j+1];
    B[7] = A[7];

    int this_index = 0;
    for (int j=0; j<8; j++) this_index += ( B[j] << j );

    apvColumn[this_index] = i;

  }

  return true;
}

} // namespace CS


////////////////////////////////////////////////////////////////////////////////
pixgem::fXY_type pixgem::fXY[pixgem::fXY_size];

bool pixgem::initconndata(void) {
  int i=0;

  fXY[i].x=1;  fXY[i].y=0;
  i++; fXY[i].x=5;  fXY[i].y=4;
  i++; fXY[i].x=6;  fXY[i].y=5;
  i++; fXY[i].x=7;  fXY[i].y=6;
  i++; fXY[i].x=2;  fXY[i].y=1;
  i++; fXY[i].x=4;  fXY[i].y=3;
  i++; fXY[i].x=3;  fXY[i].y=2;
  i++; fXY[i].x=7;  fXY[i].y=5;
  i++; fXY[i].x=2;  fXY[i].y=0;
  i++; fXY[i].x=6;  fXY[i].y=4;
  i++; fXY[i].x=8;  fXY[i].y=7;
  i++; fXY[i].x=8;  fXY[i].y=6;
  i++; fXY[i].x=3;  fXY[i].y=1;
  i++; fXY[i].x=5;  fXY[i].y=3;
  i++; fXY[i].x=4;  fXY[i].y=2;
  i++; fXY[i].x=8;  fXY[i].y=5;
  i++; fXY[i].x=9;  fXY[i].y=8;
  i++; fXY[i].x=3;  fXY[i].y=0;
  i++; fXY[i].x=7;  fXY[i].y=4;
  i++; fXY[i].x=9;  fXY[i].y=7;
  i++; fXY[i].x=9;  fXY[i].y=6;
  i++; fXY[i].x=4;  fXY[i].y=1;
  i++; fXY[i].x=6;  fXY[i].y=3;
  i++; fXY[i].x=5;  fXY[i].y=2;
  i++; fXY[i].x=9;  fXY[i].y=5;
  i++; fXY[i].x=4;  fXY[i].y=0;
  i++; fXY[i].x=8;  fXY[i].y=4;
  i++; fXY[i].x=10; fXY[i].y=9;
  i++; fXY[i].x=10; fXY[i].y=8;
  i++; fXY[i].x=5;  fXY[i].y=1;
  i++; fXY[i].x=7;  fXY[i].y=3;
  i++; fXY[i].x=6;  fXY[i].y=2;
  i++; fXY[i].x=10; fXY[i].y=7;
  i++; fXY[i].x=10; fXY[i].y=6;
  i++; fXY[i].x=5;  fXY[i].y=0;
  i++; fXY[i].x=9;  fXY[i].y=4;
  i++; fXY[i].x=10; fXY[i].y=5;
  i++; fXY[i].x=11; fXY[i].y=10;
  i++; fXY[i].x=6;  fXY[i].y=1;
  i++; fXY[i].x=8;  fXY[i].y=3;
  i++; fXY[i].x=7;  fXY[i].y=2;
  i++; fXY[i].x=11; fXY[i].y=9;
  i++; fXY[i].x=6;  fXY[i].y=0;
  i++; fXY[i].x=11; fXY[i].y=8;
  i++; fXY[i].x=11; fXY[i].y=7;
  i++; fXY[i].x=10; fXY[i].y=4;
  i++; fXY[i].x=7;  fXY[i].y=1;
  i++; fXY[i].x=9;  fXY[i].y=3;
  i++; fXY[i].x=8;  fXY[i].y=2;
  i++; fXY[i].x=11; fXY[i].y=6;
  i++; fXY[i].x=11; fXY[i].y=5;
  i++; fXY[i].x=7;  fXY[i].y=0;
  i++; fXY[i].x=12; fXY[i].y=11;
  i++; fXY[i].x=12; fXY[i].y=10;
  i++; fXY[i].x=12; fXY[i].y=9;
  i++; fXY[i].x=8;  fXY[i].y=1;
  i++; fXY[i].x=10; fXY[i].y=3;
  i++; fXY[i].x=9;  fXY[i].y=2;
  i++; fXY[i].x=12; fXY[i].y=8;
  i++; fXY[i].x=8;  fXY[i].y=0;
  i++; fXY[i].x=11; fXY[i].y=4;
  i++; fXY[i].x=12; fXY[i].y=7;
  i++; fXY[i].x=12; fXY[i].y=6;
  i++; fXY[i].x=9;  fXY[i].y=1;
  i++; fXY[i].x=13; fXY[i].y=12;
  i++; fXY[i].x=10; fXY[i].y=2;
  i++; fXY[i].x=12; fXY[i].y=5;
  i++; fXY[i].x=11; fXY[i].y=3;
  i++; fXY[i].x=9;  fXY[i].y=0;
  i++; fXY[i].x=13; fXY[i].y=11;
  i++; fXY[i].x=13; fXY[i].y=10;
  i++; fXY[i].x=13; fXY[i].y=9;
  i++; fXY[i].x=10; fXY[i].y=1;
  i++; fXY[i].x=12; fXY[i].y=4;
  i++; fXY[i].x=11; fXY[i].y=2;
  i++; fXY[i].x=13; fXY[i].y=8;
  i++; fXY[i].x=10; fXY[i].y=0;
  i++; fXY[i].x=13; fXY[i].y=7;
  i++; fXY[i].x=13; fXY[i].y=6;
  i++; fXY[i].x=12; fXY[i].y=3;
  i++; fXY[i].x=11; fXY[i].y=1;
  i++; fXY[i].x=14; fXY[i].y=13;
  i++; fXY[i].x=14; fXY[i].y=12;
  i++; fXY[i].x=13; fXY[i].y=5;
  i++; fXY[i].x=12; fXY[i].y=2;
  i++; fXY[i].x=11; fXY[i].y=0;
  i++; fXY[i].x=14; fXY[i].y=11;
  i++; fXY[i].x=14; fXY[i].y=10;
  i++; fXY[i].x=13; fXY[i].y=4;
  i++; fXY[i].x=12; fXY[i].y=1;
  i++; fXY[i].x=14; fXY[i].y=9;
  i++; fXY[i].x=14; fXY[i].y=8;
  i++; fXY[i].x=13; fXY[i].y=3;
  i++; fXY[i].x=12; fXY[i].y=0;
  i++; fXY[i].x=14; fXY[i].y=7;
  i++; fXY[i].x=14; fXY[i].y=6;
  i++; fXY[i].x=13; fXY[i].y=2;
  i++; fXY[i].x=15; fXY[i].y=14;
  i++; fXY[i].x=14; fXY[i].y=5;
  i++; fXY[i].x=13; fXY[i].y=1;
  i++; fXY[i].x=15; fXY[i].y=13;
  i++; fXY[i].x=14; fXY[i].y=4;
  i++; fXY[i].x=13; fXY[i].y=0;
  i++; fXY[i].x=15; fXY[i].y=12;
  i++; fXY[i].x=14; fXY[i].y=3;
  i++; fXY[i].x=15; fXY[i].y=11;
  i++; fXY[i].x=14; fXY[i].y=2;
  i++; fXY[i].x=15; fXY[i].y=10;
  i++; fXY[i].x=14; fXY[i].y=1;
  i++; fXY[i].x=15; fXY[i].y=9;
  i++; fXY[i].x=14; fXY[i].y=0;
  i++; fXY[i].x=15; fXY[i].y=8;
  i++; fXY[i].x=15; fXY[i].y=7;
  i++; fXY[i].x=15; fXY[i].y=6;
  i++; fXY[i].x=15; fXY[i].y=5;
  i++; fXY[i].x=15; fXY[i].y=4;
  i++; fXY[i].x=15; fXY[i].y=3;
  i++; fXY[i].x=15; fXY[i].y=2;
  i++; fXY[i].x=15; fXY[i].y=1;
  i++; fXY[i].x=15; fXY[i].y=0;
  i++; fXY[i].x=16; fXY[i].y=0;
  i++; fXY[i].x=16; fXY[i].y=1;
  i++; fXY[i].x=16; fXY[i].y=2;
  i++; fXY[i].x=16; fXY[i].y=3;
  i++; fXY[i].x=16; fXY[i].y=4;
  i++; fXY[i].x=16; fXY[i].y=5;
  i++; fXY[i].x=16; fXY[i].y=6;
  i++; fXY[i].x=16; fXY[i].y=7;
  i++; fXY[i].x=16; fXY[i].y=8;
  i++; fXY[i].x=17; fXY[i].y=0;
  i++; fXY[i].x=16; fXY[i].y=9;
  i++; fXY[i].x=17; fXY[i].y=1;
  i++; fXY[i].x=16; fXY[i].y=10;
  i++; fXY[i].x=17; fXY[i].y=2;
  i++; fXY[i].x=16; fXY[i].y=11;
  i++; fXY[i].x=17; fXY[i].y=3;
  i++; fXY[i].x=16; fXY[i].y=12;
  i++; fXY[i].x=18; fXY[i].y=0;
  i++; fXY[i].x=17; fXY[i].y=4;
  i++; fXY[i].x=16; fXY[i].y=13;
  i++; fXY[i].x=18; fXY[i].y=1;
  i++; fXY[i].x=17; fXY[i].y=5;
  i++; fXY[i].x=16; fXY[i].y=14;
  i++; fXY[i].x=18; fXY[i].y=2;
  i++; fXY[i].x=17; fXY[i].y=6;
  i++; fXY[i].x=16; fXY[i].y=15;
  i++; fXY[i].x=19; fXY[i].y=0;
  i++; fXY[i].x=18; fXY[i].y=3;
  i++; fXY[i].x=17; fXY[i].y=7;
  i++; fXY[i].x=17; fXY[i].y=8;
  i++; fXY[i].x=19; fXY[i].y=1;
  i++; fXY[i].x=18; fXY[i].y=4;
  i++; fXY[i].x=17; fXY[i].y=9;
  i++; fXY[i].x=17; fXY[i].y=10;
  i++; fXY[i].x=20; fXY[i].y=0;
  i++; fXY[i].x=19; fXY[i].y=2;
  i++; fXY[i].x=18; fXY[i].y=5;
  i++; fXY[i].x=17; fXY[i].y=11;
  i++; fXY[i].x=17; fXY[i].y=12;
  i++; fXY[i].x=20; fXY[i].y=1;
  i++; fXY[i].x=19; fXY[i].y=3;
  i++; fXY[i].x=18; fXY[i].y=6;
  i++; fXY[i].x=17; fXY[i].y=13;
  i++; fXY[i].x=21; fXY[i].y=0;
  i++; fXY[i].x=17; fXY[i].y=14;
  i++; fXY[i].x=20; fXY[i].y=2;
  i++; fXY[i].x=19; fXY[i].y=4;
  i++; fXY[i].x=21; fXY[i].y=1;
  i++; fXY[i].x=18; fXY[i].y=7;
  i++; fXY[i].x=18; fXY[i].y=8;
  i++; fXY[i].x=18; fXY[i].y=9;
  i++; fXY[i].x=22; fXY[i].y=0;
  i++; fXY[i].x=20; fXY[i].y=3;
  i++; fXY[i].x=19; fXY[i].y=5;
  i++; fXY[i].x=21; fXY[i].y=2;
  i++; fXY[i].x=18; fXY[i].y=10;
  i++; fXY[i].x=22; fXY[i].y=1;
  i++; fXY[i].x=18; fXY[i].y=11;
  i++; fXY[i].x=18; fXY[i].y=12;
  i++; fXY[i].x=20; fXY[i].y=4;
  i++; fXY[i].x=23; fXY[i].y=0;
  i++; fXY[i].x=19; fXY[i].y=6;
  i++; fXY[i].x=22; fXY[i].y=2;
  i++; fXY[i].x=21; fXY[i].y=3;
  i++; fXY[i].x=23; fXY[i].y=1;
  i++; fXY[i].x=18; fXY[i].y=13;
  i++; fXY[i].x=19; fXY[i].y=7;
  i++; fXY[i].x=19; fXY[i].y=8;
  i++; fXY[i].x=24; fXY[i].y=0;
  i++; fXY[i].x=20; fXY[i].y=5;
  i++; fXY[i].x=19; fXY[i].y=9;
  i++; fXY[i].x=23; fXY[i].y=2;
  i++; fXY[i].x=22; fXY[i].y=3;
  i++; fXY[i].x=24; fXY[i].y=1;
  i++; fXY[i].x=21; fXY[i].y=4;
  i++; fXY[i].x=19; fXY[i].y=10;
  i++; fXY[i].x=19; fXY[i].y=11;
  i++; fXY[i].x=25; fXY[i].y=0;
  i++; fXY[i].x=19; fXY[i].y=12;
  i++; fXY[i].x=24; fXY[i].y=2;
  i++; fXY[i].x=23; fXY[i].y=3;
  i++; fXY[i].x=25; fXY[i].y=1;
  i++; fXY[i].x=20; fXY[i].y=6;
  i++; fXY[i].x=20; fXY[i].y=7;
  i++; fXY[i].x=22; fXY[i].y=4;
  i++; fXY[i].x=26; fXY[i].y=0;
  i++; fXY[i].x=21; fXY[i].y=5;
  i++; fXY[i].x=20; fXY[i].y=8;
  i++; fXY[i].x=25; fXY[i].y=2;
  i++; fXY[i].x=24; fXY[i].y=3;
  i++; fXY[i].x=26; fXY[i].y=1;
  i++; fXY[i].x=20; fXY[i].y=9;
  i++; fXY[i].x=20; fXY[i].y=10;
  i++; fXY[i].x=23; fXY[i].y=4;
  i++; fXY[i].x=27; fXY[i].y=0;
  i++; fXY[i].x=20; fXY[i].y=11;
  i++; fXY[i].x=26; fXY[i].y=2;
  i++; fXY[i].x=25; fXY[i].y=3;
  i++; fXY[i].x=27; fXY[i].y=1;
  i++; fXY[i].x=21; fXY[i].y=6;
  i++; fXY[i].x=22; fXY[i].y=5;
  i++; fXY[i].x=24; fXY[i].y=4;
  i++; fXY[i].x=28; fXY[i].y=0;
  i++; fXY[i].x=21; fXY[i].y=7;
  i++; fXY[i].x=21; fXY[i].y=8;
  i++; fXY[i].x=27; fXY[i].y=2;
  i++; fXY[i].x=26; fXY[i].y=3;
  i++; fXY[i].x=28; fXY[i].y=1;
  i++; fXY[i].x=21; fXY[i].y=9;
  i++; fXY[i].x=21; fXY[i].y=10;
  i++; fXY[i].x=25; fXY[i].y=4;
  i++; fXY[i].x=29; fXY[i].y=0;
  i++; fXY[i].x=23; fXY[i].y=5;
  i++; fXY[i].x=28; fXY[i].y=2;
  i++; fXY[i].x=27; fXY[i].y=3;
  i++; fXY[i].x=29; fXY[i].y=1;
  i++; fXY[i].x=22; fXY[i].y=6;
  i++; fXY[i].x=22; fXY[i].y=7;
  i++; fXY[i].x=26; fXY[i].y=4;
  i++; fXY[i].x=30; fXY[i].y=0;
  i++; fXY[i].x=24; fXY[i].y=5;
  i++; fXY[i].x=22; fXY[i].y=8;
  i++; fXY[i].x=29; fXY[i].y=2;
  i++; fXY[i].x=28; fXY[i].y=3;
  i++; fXY[i].x=30; fXY[i].y=1;
  i++; fXY[i].x=22; fXY[i].y=9;
  i++; fXY[i].x=23; fXY[i].y=6;
  i++; fXY[i].x=27; fXY[i].y=4;
  i++; fXY[i].x=31; fXY[i].y=0;
  i++; fXY[i].x=25; fXY[i].y=5;
  i++; fXY[i].x=23; fXY[i].y=7;
  i++; fXY[i].x=23; fXY[i].y=8;
  i++; fXY[i].x=26; fXY[i].y=5;
  i++; fXY[i].x=24; fXY[i].y=6;
  i++; fXY[i].x=24; fXY[i].y=7;
  i++; fXY[i].x=25; fXY[i].y=6;

  assert( i == fXY_size-1 );

  return true;
}

////////////////////////////////////////////////////////////////////////////////
/*! \brief Converts APV channel and position to XY*/
////////////////////////////////////////////////////////////////////////////////
std::pair<int,int> pixgem::ch2xy(int channelnew, int card_virtual) {

  static bool init = false;
  if (!init) initconndata();

  if( card_virtual>=9 || card_virtual<1 || channelnew>=128 || channelnew < 0 ) {
    throw "ChipAPV::ch2xy(): bad input";
  }

  short chincon = -1;
  int sgnX = 0, sgnY =0, offsetX = -1, offsetY = -1;

  short channel = channelnew;

  switch(card_virtual)
    {
      //left bottom
    case 1:
      sgnX = 1;
      sgnY = -1;
      offsetX = 0;
      offsetY = 31;
      chincon = channel + 128;
      if( chincon>=0 )
        return  std::pair<int,int>( sgnX*(fXY[chincon].y)+offsetX, sgnY*(fXY[chincon].x)+offsetY );
      break;

      //left top
    case 2:
      sgnX = 1;
      sgnY = -1;
      offsetX = 0;
      offsetY = 31;
      chincon = channel + 0;
      if( chincon>=0 )
        return  std::pair<int,int>( sgnX*(fXY[chincon].y)+offsetX, sgnY*(fXY[chincon].x)+offsetY );
      break;

      //top left
    case 3:
      sgnX = -1;
      sgnY = -1;
      offsetX = 31;
      offsetY = 31;
      chincon = channel +  128;
      if( chincon>=0 )
        return  std::pair<int,int>( sgnX*(fXY[chincon].x)+offsetX, sgnY*(fXY[chincon].y)+offsetY );
      break;

      //top right
    case 4:
      sgnX = -1;
      sgnY = -1;
      offsetX = 31;
      offsetY = 31;
      chincon = channel + 0;
      if( chincon>=0 )
        return  std::pair<int,int>( sgnX*(fXY[chincon].x)+offsetX, sgnY*(fXY[chincon].y)+offsetY );
      break;

      //right top
    case 5:
      sgnX = -1;
      sgnY = 1;
      offsetX = 31;
      offsetY = 0;
      chincon = channel +  128;
      if( chincon>=0 )
        return  std::pair<int,int>( sgnX*(fXY[chincon].y)+offsetX, sgnY*(fXY[chincon].x)+offsetY );
      break;

      //right bottom
    case 6:
      sgnX = -1;
      sgnY = 1;
      offsetX = 31;
      offsetY = 0;
      chincon = channel +  0;
      if( chincon>=0 )
        return  std::pair<int,int>( sgnX*(fXY[chincon].y)+offsetX, sgnY*(fXY[chincon].x)+offsetY );
      break;

      //bottom right
    case 7:
      sgnX = 1;
      sgnY = 1;
      offsetX = 0;
      offsetY = 0;
      chincon = channel + 128;
      if( chincon>=0 )
        return  std::pair<int,int>( sgnX*(fXY[chincon].x)+offsetX, sgnY*(fXY[chincon].y)+offsetY );
      break;

      //bottom left
    case 8:
      sgnX = 1;
      sgnY = 1;
      offsetX = 0;
      offsetY = 0;
      chincon = channel + 0;
      if( chincon>=0 )
        return  std::pair<int,int>( sgnX*(fXY[chincon].x)+offsetX, sgnY*(fXY[chincon].y)+offsetY );

      break;

    default:
      throw "ChipAPV::ch2xy(): bad input, wrong card_virtual";
      break;
    }

  return  std::pair<int,int>(-1,-1);
}

////////////////////////////////////////////////////////////////////////////////
/*! \brief Converts detector channel to XY*/
////////////////////////////////////////////////////////////////////////////////
std::pair<int,int> pixgem::detch2xy(int chan) {

  // Convert detector channel to format for ch2xy: (ch 0..127, chip position)
  int chip_position = (chan/128) + 1;
  int chdum = 127- (chan%128);
  int pixX=-1, pixY=-1;
  bool mapping_valid = true;

  if (chdum<0 || chdum>127 || chip_position<1 || chip_position>8)
    {
      //      mapping_valid = false;
    }
  else
    {
      std::pair<int,int> xy=ch2xy(chdum,chip_position);
      pixX = xy.first;
      pixY = xy.second;
    }

  return  std::pair<int,int>(pixX,pixY);
}


//------------------------------------------------------------------------------


