/*!
    \file ChipF1.cc
    \author Alexander Zvyagin
*/

#include <cassert>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "ChipF1.h"
#include "DaqError.h"
#include "DaqOption.h"
#include "TriggerTime.h"
#include "DaqEvent.h"
#include "utils.h"

using namespace std;

namespace CS {

////////////////////////////////////////////////////////////////////////////////

ChipF1::ChipF1(void const * const buf,bool copy_buf,DaqOption &opt,DaqEvent &ev)
:   Chip(buf,copy_buf,opt,ev)
{
    Clear();
}

////////////////////////////////////////////////////////////////////////////////

void ChipF1::Decode(const Maps &maps,Digits &digits_list,DaqOption &opt)
{
    // Scan the chip data if it was not done yet.
    if( !is_scaned )
        Scan(opt);

    size_t size=digits_list.size();

    // Loop for all found data lines
    for( vector< pair<DataID,uint16> >::iterator it=data_all.begin();
         it!=data_all.end(); it++ )
    {
        typedef Maps::const_iterator m_it; // Create a short name for map's iterator
        pair<m_it,m_it> m_range = maps.equal_range(it->first); // all maps with given data ID

        if( m_range.first==m_range.second )
        {
            // Map is not found! But may be it is a RICH-PMT?
            DataID d=it->first;
            d.u.s.hit = d.u.s.chip_chan = 255;
            m_range = maps.equal_range(d);
        }

        for( m_it c=m_range.first; c!=m_range.second; c++ )
        {
            const Digit *digit1 = dynamic_cast<Digit*>(c->second);
            if( digit1==NULL )
                throw Exception("ChipF1::Decode(): Internal error");
            
            Digit *digit2=NULL;

            const DigitRICHPMT *digit11 = dynamic_cast<const DigitRICHPMT *>(digit1);
            if( digit11!=NULL )
            {
                DigitRICHPMT *d = new DigitRICHPMT(*digit11);
                d->SetDataID(it->first);
                d->Finalize();
                digit2 = d;
            }
            else
                digit2 = new Digit(*digit1);                     // Clone the digit

            assert( digit11==NULL || NULL!=dynamic_cast<const DigitRICHPMT *>(digit2));

            digit2->SetAmplitude( it->second );                     // Set the measured time
            
            // If high-resolution bit was set, set time_unit to 'high' value,
            // otherwise use settings from the mapping
            //if( GetSLink().GetFormat()&16 )
            //    digit2->SetTimeUnit( ChipF1::GetUnitHigh() );
            digits_list.insert(pair<DetID,Digit*>(digit2->GetDetID(),digit2));
        }

        if( size>=digits_list.size() )
            Exception("ChipF1::Decode():WW: some maps were not found for srcID=%d port_geoID=%d chan=%d",
                       GetSourceID(),int(uint16(it->first>>16)),int(uint8(it->first>>8)));
        
        size=digits_list.size();
    }
}

////////////////////////////////////////////////////////////////////////////////

void ChipF1::Scan(DaqOption &opt)
{
    if( is_scaned )
        return;
    else
        is_scaned=true;

    if (GetSLink().IsFirstEventInRun())
    {
        Exception("ChipF1::Scan():WW: first event of run").Print();
        return;
    }

    if( IsDebugMode() )
        ScanDebugMode(opt);
    else
        ScanSparsifiedMode(opt);
}

////////////////////////////////////////////////////////////////////////////////

// This code comes from Wolfgang Kastaun.
void ChipF1::ScanDebugMode(DaqOption &opt)
{
    uint32 source=GetSLink().GetSourceID();
  
    const bool cmctdc_global=(0==(GetSLink().GetFormat() & 0x40)); //TDC-CMC flag, check also the call: IsDataF1CMC()
    uint32 enrs=GetSLink().GetEventNumber()&0x3f; //last 6 bits of event number
    bool hasAllHeaders = GetSLink().GetEventNumber()%128 == 1; //128th event

    //Debug Mode Checks
    uint32 port=1000;
    uint32 lasterr=0;

    bool expect_trailer=false;
    bool lastwaserror=false;
    bool comptrig=false;
    uint32 trigtime=0;

    for( const uint32 *p=GetDataStart(); p<GetDataEnd(); p++ )
    {
        bool cmctdc = cmctdc_global;

        // If this is an error word, report it and go to the next data line
        ErrorWord errw(*p,ErrorMarker,source,*this,opt);
        if( errw.IsError() )
        {
            lastwaserror=true;
            lasterr=errw.GetErrorCode();
            continue;
        }

        // It is not an error word. So it is one of the header/trailer/data words.

        // Debugmode Data/header/trailer
        HeaderDbg header=*p;    // Can be a trailer also...

        if (header.IsHeader())   // Header or trailer
        {
            // Check do we need to set the 'cmc' flag by force.
            if( opt.IsDataF1CMC(source,header.GetPort()) )
                cmctdc = true;

            if (!lastwaserror)
            { 
	      hasAllHeaders=true;
	      if (comptrig)
		//compass mastertime contains TDC with on purpose different time resolutions
                //and different trigger times, it should not be checked for trigger time 
                if ( !opt.IsMasterTimeSrcID(source) &&
                     !opt.SkipMasterTimeCheckSrcID(source) )       
		    // the maximum trigger time deviation is 2
		    if (abs(TimeDifference(trigtime<<7,header.GetData()<<7,64239,0))>(2<<7)  &&
		        abs(TimeDifference(trigtime<<7,header.GetData()<<7,64600,0))>(1<<7)) 
		      AddError(DaqError(DaqError::TDC_WRTIME,DaqError::SOURCE_ID,source,DaqError::PORT,header.GetPort()),opt);
		comptrig=true;
                trigtime=header.GetData();
            }

            bool b0=true;
            if (lastwaserror) 
            {
                if (1==lasterr) 
                {
                    port=header.GetPort();
                    expect_trailer=false;
                }
                else if (5==lasterr)
                    expect_trailer=true;      
                else if (3==lasterr || 7==lasterr) 
                {
                    if (!expect_trailer)
                        port=header.GetPort();
                }
                else 
                {
                    expect_trailer=true;
                    b0=false;
                }
            }
            else if (!expect_trailer)
                port=header.GetPort();

            if (hasAllHeaders && b0 && expect_trailer && (port!=header.GetPort())) 
                AddError(DaqError(DaqError::TDC_WRPORT_T,
                                  DaqError::SOURCE_ID,source,
                                  DaqError::PORT,header.GetPort(),
                                  DaqError::CHIP,header.GetChip(),
                                  DaqError::CHANNEL,header.GetChannel()),opt);

            expect_trailer=!expect_trailer;

            if( !lastwaserror || 3==lasterr || 7==lasterr )
            {
                if (header.IsTBOBitSet())
                    AddError(DaqError(DaqError::TDC_TBO,
                                      DaqError::SOURCE_ID,source,
                                      DaqError::PORT,    header.GetPort(),
                                      DaqError::CHIP,    header.GetChip(),
                                      DaqError::CHANNEL, header.GetChannel()),opt);

                // For TDC-CMC all PLL bits should be 1 
                if (opt.Use_TDC_WRPLL() && cmctdc && (header.GetPLL()!=0xf))
                    AddError(DaqError(DaqError::TDC_WRPLL_H,
                                      DaqError::SOURCE_ID,source,
                                      DaqError::PORT,    header.GetPort(),
                                      DaqError::CHIP,    header.GetChip(),
                                      DaqError::CHANNEL, header.GetChannel()),opt);
            }

            if (!lastwaserror)
            {
                //check if event numbers are the same
                if (header.GetEventNumber()<enrs)
                    AddError(DaqError(DaqError::TDC_ERR3U,
                                      DaqError::SOURCE_ID,source,
                                      DaqError::PORT,    header.GetPort(),
                                      DaqError::CHIP,    header.GetChip(),
                                      DaqError::CHANNEL, header.GetChannel()),opt);

                if (header.GetEventNumber()>enrs)
                    AddError(DaqError(DaqError::TDC_ERR7U,
                                      DaqError::SOURCE_ID,source,
                                      DaqError::PORT,    header.GetPort(),
                                      DaqError::CHIP,    header.GetChip(),
                                      DaqError::CHANNEL, header.GetChannel()),opt);
            }

            if ( lastwaserror && ( 3==lasterr || 7==lasterr ) ) 
            {
                if (header.GetEventNumber()==enrs)
                    AddError(DaqError(DaqError::TDC_WRERR,
                                      DaqError::SOURCE_ID,source,
                                      DaqError::PORT,    header.GetPort(),
                                      DaqError::CHIP,    header.GetChip(),
                                      DaqError::CHANNEL, header.GetChannel()),opt);
            }
        } // header or trailer end
        else 
        {
            //Debugmode Data
            DataDbg data_word=*p;
            bool good=true;  // Will be true if data line is good.
            
	    if (!expect_trailer && hasAllHeaders) 
            { 
                AddError(DaqError(DaqError::TDC_ERR1U,
                                  DaqError::SOURCE_ID,source,
                                  DaqError::PORT,    header.GetPort(),
                                  DaqError::CHIP,    header.GetChip(),
                                  DaqError::CHANNEL, header.GetChannel()),opt);
                port=data_word.GetPort();
                good=false;
            }
            expect_trailer=true;

            if ((port!=data_word.GetPort()) && hasAllHeaders)
            {
                AddError(DaqError(DaqError::TDC_WRPORT,
                                  DaqError::SOURCE_ID,source,
                                  DaqError::PORT,    data_word.GetPort(),
                                  DaqError::CHIP,    data_word.GetChip(),
                                  DaqError::CHANNEL, data_word.GetChannel()),opt);
                good=false;
            }

            //For TDC-CMC all PLL bits should be 1 :
            if( opt.Use_TDC_WRPLL() )
                if (cmctdc && (data_word.GetPLL()!=0xf))
                {
                    AddError(DaqError(DaqError::TDC_WRPLL,
                                      DaqError::SOURCE_ID,source,
                                      DaqError::PORT,    data_word.GetPort(),
                                      DaqError::CHIP,    data_word.GetChip(),
                                      DaqError::CHANNEL, data_word.GetChannel()),opt);
                    good=false;
                }

            if( good && !event->IsBadDataSource(GetSourceID(),data_word.GetPort()) )
                AddData( data_word.GetData16bit(), 0, data_word.GetPort(), data_word.GetChipChannel() );
        }
        lastwaserror=false;
    } // End of the data loop
}

///////////////////////////////////////////////////////////////////p/////////////

void ChipF1::ScanSparsifiedMode(DaqOption &opt)
{
    // Sparsified mode
    for( const uint32 *p=GetDataStart(); p<GetDataEnd(); p++ )
    {
        if( ErrorWord(*p,ErrorMarker,GetSourceID(),*this,opt).IsError() )
        {
            if( p+1>=GetDataEnd() )
                Exception("ChipF1::ScanSparsifiedMode(): A data/trailer/header is expected after an error!");
            p++;
            // Now (*p) is faulty data or header/trailer
            #warning "ChipF1::Decode(): Skip all faulty data/trailer/header words"
            continue;
        }
      
        const uint16 geoID=(*p)>>22;
        const uint8 chip_chan=((*p)>>16)&63;

        // Check that srcID and geoID are OK.
        if( !event->IsBadDataSource(GetSourceID(),geoID) )
            AddData( uint16(*p), 1, geoID, chip_chan );
    }
}

////////////////////////////////////////////////////////////////////////////////

void ChipF1::AddData(uint16 time,uint16 mode,uint16 geoID_or_port,uint8 chip_chan)
{
    const bool is_latch_mode = IsLatchMode();       // to speed up the code
    uint8 hit_pattern=0;                            // the hit pattern for latch mode
    if( is_latch_mode )
    {
        hit_pattern = time&15;                      // Last 4 bits in latch mode is a hit pattern
        time &= ~15;                                // Set last for bits to zero.
    }

    // Add data for decoding
    if( is_latch_mode )
    {
        for( unsigned int hit=0; hit<4; hit++ )
            if( (1<<hit)&hit_pattern )
            {
                DataID data_id(GetSourceID(),mode,geoID_or_port,chip_chan,hit);
                data_all.push_back( pair<DataID,uint16>(data_id,time) );
            }
    }
    else
    {
        DataID data_id(GetSourceID(),mode,geoID_or_port,chip_chan,0);
        data_all.push_back( pair<DataID,uint16>(data_id,time) );
    }
}

////////////////////////////////////////////////////////////////////////////////

void ChipF1::Print(ostream &o,const string &prefix) const
{
  Chip::Print(o,prefix+"ChipF1 ");
  for( const uint32 *d=GetDataStart(); d<GetDataEnd(); d++ )
  {
    if( IsDebugMode() )
    {
      HeaderDbg h(*d);
      if( h.IsHeader() )
        h.Print(o,prefix);
      else
        DataDbg(*d).Print(o,prefix);
    }
    else
    {
      printf("%s%3zu  %8x\n",prefix.c_str(),d-GetDataStart(),*d);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

void ChipF1::Data::Print(ostream &o,const string &prefix) const
{
  o << prefix << GetName() << "  ";
  char s[222];
  sprintf(s,"ChipChan=%2u=(%2u:%2u) port=%2u PLL=%2u",
            (GetChannel()<<3)+GetChip(),GetChannel(),GetChip(),GetPort(),GetPLL() );
  o<<s;
}

////////////////////////////////////////////////////////////////////////////////

void ChipF1::DataDbg::Print(ostream &o,const string &prefix) const
{
  ChipF1::Data::Print(o,prefix);
  char s[222];
  sprintf(s," data16=%5u    (data12=%4u data4=%2u)\n",
             GetData(),GetData12bit(),GetData4bit());
  o<<s;
}

////////////////////////////////////////////////////////////////////////////////

void ChipF1::HeaderDbg::Print(ostream &o,const string &prefix) const
{
  ChipF1::Data::Print(o,prefix);
  char s[222];
  sprintf(s," time=%3u\n",GetData());
  o<<s;
}

////////////////////////////////////////////////////////////////////////////////

ChipF1::Map::Map(const ObjectXML &o)
:   Chip::Map(o),
    mode('?'),
    geo_id(-1),
    port(-1),
    chip(0),
    chan_in_chip(0),
    hit(0),
    x(0), y(0),
    type_ul(' '),
    rich_pmt_expert(false),
    wire_position(0),
    time_unit(-1),
    time_reference(0)
{
    if( version==0 )
        version=1;

    if( GetName()!="ChipF1" )
        throw Exception("ChipF1::Map::Map(): Internal error.");

    if( GetVersion()<1 || GetVersion()>4  )
        throw Exception("ChipF1::Map::Map(): unknown version %d",GetVersion());

    GetAttribute("time_reference",time_reference);

    // If "time_unit" has been set in a mapping file, then use it ....
    if( NULL==GetAttribute("time_unit",time_unit) )
    {
        SetAttribute("time_unit","0.12892312");
        GetAttribute("time_unit",time_unit);
    }

    if( IsOption("RICH-PMT") || IsOption("RICH-PMT-expert") )
    {
        mode = 'r';

        string name;
        istringstream s(dec_line.c_str());

        s >> name >> source_id >> port >> geo_id >> x >> y >> type_ul;

        if( s.fail() || NULL==strchr("UuLl",type_ul) )
            throw Exception("ChipF1::Map::Map(): bad format (RICH-PMT) in line: %s",map_line.c_str());

        chanF=wireF= 0;
        chanL=wireL=63;
        chanS=wireS= 1;
        
        id=DetID(name);
        
        if( IsOption("RICH-PMT-expert") )
            rich_pmt_expert=true;
    }
    else
    if( IsOption("tdc") || IsOption("latch") )
    {
        string name;
        int chanN;
        istringstream s(dec_line.c_str());

        if( GetVersion()==4 )
        {
            if( IsOption("xy") )
            {
                s >> name >> source_id >> port >> geo_id >> chip >> chanF >> x >> y;
                chanN=1;
                chanS=1;
                wireF=wireL=-1;
                wireS=1;
            }
            else
                s >> name >> source_id >> port >> geo_id >> chanF>>chanS>>chanN>>wireF>>wireL>>wireS;

            if( s.fail() )
                throw Exception("ChipF1::Map::Map(): bad format in line: %s",map_line.c_str());

            s >> wire_position;
            if( s.fail() )
                wire_position=0;

            chanL=chanF+(chanN-1)*chanS;

            mode = IsOption("tdc") ? 't' : 'l';

            id=DetID(name);
        }
        else
        {
            // Support of the old version
            if( GetVersion()==1 )
              { int det_n; s>>det_n;}

            if( GetVersion()==2 && IsOption("latch") )
              s >> name >> source_id >>port>>chip>>chan_in_chip>>chanF>>chanS>>chanN>>wireF>>wireL>>wireS;
            else
            if( GetVersion()==2 && IsOption("xy") )
            {
              s >> name >> source_id >>port>>chip>>chanF>>x>>y;
              chanN=1;
              chanS=1;
              wireF=wireL=-1;
              wireS=1;
            }
            else
              s >> name >> source_id >>port>>chip>>chanF>>chanS>>chanN>>wireF>>wireL>>wireS;

            if( s.fail() )
              throw Exception("ChipF1::Map::Map(): bad format in line: %s",map_line.c_str());

            s >> wire_position;
            if( s.fail() )
              wire_position=0;
        }

        chanL=chanF+(chanN-1)*chanS;

        mode = IsOption("tdc") ? 't' : 'l';

        id=DetID(name);
    }
    else
        throw Exception("ChipF1::Map::Map(): unknown option(s) \"%s\" for line \"%s\"",
                         options.c_str(),map_line.c_str());

    Check();
}

////////////////////////////////////////////////////////////////////////////////

void ChipF1::Map::AddToMaps(Maps &maps,DaqOption &options) const
{
  size_t size=maps.size();
  int factor=1;

  string trigger_mask_data;
  if( NULL!=GetAttribute("trigger_mask_data",trigger_mask_data) )
  {
      size_t n = trigger_mask_data.find('-');
      if( n==string::npos )
        throw Exception("ChipF1::Map::AddToMaps(): Bad attribute \"trigger_mask_data\"=\"%s\"",trigger_mask_data.c_str());
      int
        min = atoi(trigger_mask_data.substr(0,n).c_str()),
        max = atoi(trigger_mask_data.substr(n+1,string::npos).c_str());

      options.GetTTConfig().SetTriggerMaskDataLimit(min,max);
      options.GetTTConfig().trigger_mask_DetID = GetDetID();
      options.GetTTConfig().trigger_mask_srcID = GetSourceID();
      options.GetTTConfig().trigger_mask_port  = port;
      GetAttribute("trigger_time_jitter",options.GetTTConfig().time_jitter);
      GetAttribute("trigger_time_precise",options.GetTTConfig().time_precise);
      GetAttribute("triggers_diff_sigma",options.GetTTConfig().triggers_diff_sigma);
      const TriggerTimeConfig::TTC *ttc = options.GetTTConfig().Find(time_unit);
      if( ttc==NULL )
        throw Exception("ChipF1::Map::AddToMaps(): TT for the trigger_mask time unit %g is not found!",time_unit);
      options.GetTTConfig().trigger_mask_TT_index = ttc->index;
  }
    
  string TCS_phase_DetID;
  if( NULL!=GetAttribute("TCS_phase",TCS_phase_DetID) && TCS_phase_DetID==GetDetID().GetName() )
  {
      // test that the "TCS_phase" detector is only mapped to one channel
      if ( 1!=GetChanN() )
        throw Exception("Too many detector channels mapped to the TCS phase: %d instead of 1",GetChanN());

      options.GetTTConfig().TCS_phase_DetID   = TCS_phase_DetID;
      options.GetTTConfig().TCS_phase_srcID   = GetSourceID();
      options.GetTTConfig().TCS_phase_port    = port;
      options.GetTTConfig().TCS_phase_channel = GetChanF();
  }

  switch( mode )
  {
    case 'r':
    {
        // This is RICH-PMT mode.

        uint16 nn[2]={port,geo_id};
        for( uint16 i=0; i<2; i++ )
        {
            DataID data_id(GetSourceID(),i,nn[i],255,255);
            DigitRICHPMT *digit = new DigitRICHPMT(data_id,GetDetID());
            digit->SetTimeUnit(time_unit);
            digit->SetTimeReference(time_reference);
            digit->SetDCardX(GetWireX());
            digit->SetDCardY(GetWireY());
            digit->SetType(GetType());
            if( rich_pmt_expert )
                digit->SetExpertDecoding(true);

            if( maps.end()!=maps.find(data_id) && !IsMultiDigit() )
            {
                Print();
                Exception("ChipF1::Map::AddToMaps(): map already exists!").Print();
            }

            maps.insert( pair<DataID,Digit*>(data_id,digit) );
        }
        
        break;
    }
    
    case 'l':
    {
      switch( GetVersion() )
      {
        case 1:
        case 4:
        {
          factor=2;

          for( int _chip=0; _chip<8; _chip++ )
            for( int _chan=0; _chan<8; _chan++ )
              for( int _hit=0; _hit<4; _hit++ )
              {
                int wire = _chip*32 + (_chan/4)*16 + _chan%4 + _hit*4;
                try        { wire=CalcWire(wire); }
                catch(...) { continue;            }


                // Create TWO entries in 'maps': one for debug mode (port usage) and another for
                // sparsified mode (for geoID).
                uint16 nn[2]={port,geo_id};
                for( uint16 i=0; i<2; i++ )
                {
                    DataID data_id(GetSourceID(),i,nn[i],(_chip<<3)|_chan,_hit);
                    Digit *digit = new Digit(data_id,GetDetID(),wire,wire_position,0,time_unit);
                    digit->SetTimeReference(time_reference);

                    if( maps.end()!=maps.find(data_id) && !IsMultiDigit() )
                    {
                        Print();
                        Exception("ChipF1::Map::AddToMaps():a: map already exists!").Print();
                    }

                    maps.insert( pair<DataID,Digit*>(data_id,digit) );
                }
              }

          break;
        }

        case 2:
        {
          if( GetChanN()>4 )
            throw Exception("ChipF1::Map::AddToMaps(): too big chanN=%d",GetChanN());

          for( unsigned int n=0; n<GetChanN(); n++ )
          {
            int hit =GetChanF()+n*GetChanS();
            int wire=GetWireF()+n*GetWireS();

            DataID data_id(GetSourceID(),uint16(0),uint16(port),(chip<<3)|chan_in_chip,hit);
            Digit *digit = new Digit(data_id,GetDetID(),wire,wire_position,0,time_unit);
            digit->SetTimeReference(time_reference);

            if( maps.end()!=maps.find(data_id) && !IsMultiDigit() )
            {
                Print();
                Exception("ChipF1::Map::AddToMaps():b: map already exists!").Print();
            }

            maps.insert( pair<DataID,Digit*>(data_id,digit) );
          }
          break;
        }
        
        default:
          throw Exception("ChipF1::Map::AddToMaps(): unknown version %d for latch mode",GetVersion());
      }
      break;
    }
    
    case 't':
    {
        bool cmc_flag = IsOption("cmc");
  
        switch( GetVersion() )
        {
            case 1:
            case 2:
            case 3:
            {
                for( size_t n=0; n<GetChanN(); n++ )
                {
                  int wire=GetWireF()+n*GetWireS();
                  int chan=GetChanF()+n*GetChanS();
                  int _chip = chip==-1 ? (chan>>3) : chip;
                  int _chan = chip==-1 ? (chan&7 ) : chan;
                  DataID data_id(GetSourceID(),uint16(0),uint16(port),(_chip<<3)|_chan,0);

                  Digit *digit;
                  if( IsOption("xy") )
                      digit = new Digit(data_id,GetDetID(),x,y,wire_position,0,time_unit);
                  else
                      digit = new Digit(data_id,GetDetID(),wire,wire_position,0,time_unit);
                  digit->SetTimeReference(time_reference);

                  if( maps.end()!=maps.find(data_id) && !IsMultiDigit() )
                  {
                      Print();
                      Exception("ChipF1::Map::AddToMaps():c: map already exists!").Print();
                  }

                  maps.insert( pair<DataID,Digit*>(data_id,digit) );
                  if( cmc_flag )
                    options.AddDataF1CMC(GetSourceID(),port);
                }
                break;
            }
            
            case 4:
            {
                factor=2;
                for( size_t n=0; n<GetChanN(); n++ )
                {
                    int wire=GetWireF()+n*GetWireS();
                    uint16 chan=GetChanF()+n*GetChanS();
                    if( chan>=64 )
                      throw Exception("ChipF1::Map::AddToMaps(): too big electronic channel=%d",chan);

                    // Create TWO entries in 'maps': one for debug mode (port usage) and another for
                    // sparsified mode (for geoID).
                    uint16 nn[2]={port,geo_id};
                    for( uint16 i=0; i<2; i++ )
                    {
                        DataID data_id(GetSourceID(),i,nn[i],chan,0);

                        Digit *digit;
                        if( IsOption("xy") )
                            digit = new Digit(data_id,GetDetID(),x,y,wire_position,0,time_unit);
                        else
                            digit = new Digit(data_id,GetDetID(),wire,wire_position,0,time_unit);
                        digit->SetTimeReference(time_reference);

                        if( maps.end()!=maps.find(data_id) && !IsMultiDigit() )
                        {
                            Print();
                            Exception("ChipF1::Map::AddToMaps():d: map already exists!").Print();
                        }

                        maps.insert( pair<DataID,Digit*>(data_id,digit) );
                    }
                    if( cmc_flag )
                      options.AddDataF1CMC(GetSourceID(),port);
                }
                break;
            }
            default:
                throw Exception("ChipF1::Map::AddToMaps(): unknown version %d",GetVersion());
        }
        
        ReadTTConfig(options.GetTTConfig());
        
        break;
    }

    default:
        throw Exception("ChipF1::Map::AddToMaps(): unknown mode %c",mode);
    
  }

    if( (mode=='t' || mode=='l') && (maps.size()-size)!=factor*GetChanN() )
        throw Exception("ChipF1::Map::AddToMaps(): internal error in Map mode=%c version=%d",mode,GetVersion());

    maps.SetWires( GetDetID(), maps.GetWires(GetDetID())+GetChanN() );
}

////////////////////////////////////////////////////////////////////////////////

void ChipF1::Map::ReadTTConfig(TriggerTimeConfig &tt_conf) const
{
    // Check for the trigger time decoding options.
    int TT_index=-1, TT_index_recover=-1, TT_overolling=0, TT_channels=0;
    double TT_sigma=0, time_unit=0;
    string sss;

    if( NULL!=GetAttribute("TT_index",TT_index) )
    {
        bool ok = GetAttribute("TT_channels",  TT_channels  )   !=NULL &&
                  //GetAttribute("TT_shift",     TT_shift     )   !=NULL &&
                  //GetAttribute("TT_shift_ref", TT_shift_ref )   !=NULL &&
                  GetAttribute("TT_sigma",     TT_sigma )       !=NULL &&
                  GetAttribute("TT_overolling",TT_overolling)   !=NULL &&
                  GetAttribute("time_unit",    time_unit);
        if(!ok)
            throw Exception("ChipF1::Map::Map(): Bad settings for the trigger time");
        
        GetAttribute("TT_index_recover",TT_index_recover);
            
        TriggerTimeConfig::TTC ttc(id,time_unit,TT_index,TT_index_recover,TT_channels,TT_overolling,TT_sigma);
        ttc.srcID = GetSourceID();
        ttc.port  = port;
        GetAttribute("MT_shift",ttc.MT_shift);
        tt_conf.Add(ttc);
    }
}

////////////////////////////////////////////////////////////////////////////////

void ChipF1::Map::Print(ostream &o,const string &prefix) const
{
  Chip::Map::Print(o,prefix);
  o<<prefix;

  char s[222];
  sprintf(s,"ChipF1::Map: mode=%c port=%d chip=%d\n",
             mode,int(port),int(chip));
  o << s;
}

////////////////////////////////////////////////////////////////////////////////

void ChipF1::Digit::Print(ostream &o,const string &prefix) const
{
    char s[222], s2[44];

    if( xy_mode )
        sprintf(s2,"(x,y)=(%d,%d)",GetX(),GetY());
    else
        sprintf(s2,"wire=%d",GetChannel());

    const DataID &d = reinterpret_cast<const DataID &>(GetDataID());

    sprintf(s,"%s%8s ampl=%5d  %s  pos=%d time_unit=%gns time_decoded=%g Data=(srcID=%d,geoID/port=%d,chip/chan=%d,hit=%d,mode=%d)\n",
               prefix.c_str(),GetDetID().GetName().c_str(),(int)GetAmplitude(),
               s2,(int)GetChannelPos(),GetTimeUnit(),GetTimeDecoded(),
               unsigned(d.u.s.src_id),unsigned(d.u.s.geoID_or_port),unsigned(d.u.s.chip_chan),unsigned(d.u.s.hit),unsigned(d.u.s.mode));
  o << s;
}

////////////////////////////////////////////////////////////////////////////////

void ChipF1::DigitRICHPMT::Print(ostream &o,const string &prefix) const
{
    char s[222];

    const DataID &d = reinterpret_cast<const DataID &>(GetDataID());

    sprintf(s,"%s%8s pixel=(%2d,%2d)  ampl=%5d  Data=(srcID=%d,geoID/port=%d,chip=%d,chan=%d,type=%c)\n",
               prefix.c_str(),GetDetID().GetName().c_str(),
               (int)GetX(),(int)GetY(),(int)GetAmplitude(),
               unsigned(d.u.s.src_id),unsigned(d.u.s.geoID_or_port),
               unsigned(d.u.s.chip_chan)>>3,unsigned(d.u.s.chip_chan)&0x7,GetType());
  o << s;
}

////////////////////////////////////////////////////////////////////////////////

void ChipF1::DigitRICHPMT::Finalize(void)
{
    DataID d = GetDataID();
    
    uint8 chip = d.u.s.chip_chan>>3;
    uint8 chan = d.u.s.chip_chan & 0x7;
    
    assert( chip<=7 && chan<=7 );
    
    // Convert chan to full PMT (two F1 chips)
    if( chip%2 )
        chan += 8;
    
    uint16 x,y;  // pixel coordinates.
    
    // Get (x,y) in the PMT
    x = (15-chan)/4;
    y = 3-chan%4;
    
    if( !expert_decoding )
    {
        // Make the optical mirroring.
        x = 3-x;
        y = 3-y;
    }
    
    // Go to the system of reference of the full D-card
    y += 4*(3-uint8(chip/2));
    
    assert( x<4  && y<16 );
    
    // Apply rotation for A/B cards orientation
    if( GetType()=='u' || GetType()=='U' )
    {
        x =  3-x;
        y = 15-y;
    }

    assert( x<4  && y<16 );
    
    // Now go to the system of refernce of the photocathode
    assert( GetDCardX()<12 && GetDCardY()<3 );
    x += GetDCardX()* 4;    // 48 max
    y += GetDCardY()*16;    // 48 max
    
    assert( x<48 && y<48 );

    if( !expert_decoding )
    {
        // Swapping of Xcoordiante -> to go to coral system of refernce (photonview)
        x = 47 - x;
    }

    SetChannel( (y<<16) + x );

    assert( x==GetX() && y==GetY() );
}

////////////////////////////////////////////////////////////////////////////////

vector<float> ChipF1::DigitRICHPMT::GetNtupleData(void) const
{
  vector<float> v;
  v.push_back(GetX());
  v.push_back(GetY());
  v.push_back(GetTimeDecoded());
  return v;
}

////////////////////////////////////////////////////////////////////////////////

vector<float> ChipF1::Digit::GetNtupleData(void) const
{
  vector<float> v;
  v.push_back(GetChannel());
  v.push_back(GetTime());
  v.push_back(GetChannelPos());
  v.push_back(GetTimeUnit());
  v.push_back(GetTimeDecoded());
  return v;
}

////////////////////////////////////////////////////////////////////////////////

double ChipF1::TimeDifference(int time,double trigger_time,int time_overolling,double time_ref,double cut)
{
  double diff=time-trigger_time-time_ref;

  // Maximum two iterations
  for( int i=0; fabs(diff)>=cut && i<3; i++ )
    if( diff<0 )
      diff += time_overolling;
    else
      diff -= time_overolling;

  return diff+time_ref;
}

////////////////////////////////////////////////////////////////////////////////

}
