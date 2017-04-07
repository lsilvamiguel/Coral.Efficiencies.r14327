/*!
  \file ChipSinica.cc
  \author Alexander Zvyagin
*/

#include <cassert>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "ChipSinica.h"
#include "DaqError.h"
#include "DaqOption.h"
#include "TriggerTime.h"
#include "DaqEvent.h"
#include "utils.h"

using namespace std;

namespace CS {

  // maximum size of data package from FEM: 1header+32data+1trailer words
  const float ChipSinica::FEM_DATA_SIZE_MAX = 34; 
  const float ChipSinica::FEM_TEMPERATURE_MAX=50;

  ////////////////////////////////////////////////////////////////////////////////

  ChipSinica::ChipSinica(void const * const buf,bool copy_buf,DaqOption &opt,DaqEvent &ev)
    :   Chip(buf,copy_buf,opt,ev)
  {
    Clear();
  }

  ////////////////////////////////////////////////////////////////////////////////

  void ChipSinica::Decode(const Maps &maps,Digits &digits_list,DaqOption &opt)
  {
    // Scan the chip data if it was not done yet.
    if( !is_scaned )
      Scan(opt);

    size_t size=digits_list.size();

    int counter = 0;
    // Loop for all found data lines
    for( vector< pair<DataID,Hit> >::iterator it=data_all.begin();
         it!=data_all.end(); it++ )
      {
        typedef Maps::const_iterator m_it; // Create a short name for map's iterator
        pair<m_it,m_it> m_range = maps.equal_range(it->first); // all maps with given data ID

	for( m_it c=m_range.first; c!=m_range.second; c++ )
	  {
	    const Digit *digit1 = dynamic_cast<Digit*>(c->second);
            if( digit1==NULL )
	      throw Exception("ChipSinica::Decode(): Internal error");
            
            Digit *digit2=NULL;
	    digit2 = new Digit(*digit1);                              // Clone the digit

            digit2->SetHitTime(it->second.GetHitTime());              // Set the measured time
            digit2->SetTriggerTime(it->second.GetTriggerTime());      // Set the trigger time
	    digit2->SetTemperature(it->second.GetTemperature());      // Set temperature
	    
	    // printf("digit2 - "); digit2->Print();
            digits_list.insert(pair<DetID,Digit*>(digit2->GetDetID(),digit2));
	  }

        if( size>=digits_list.size() )
	  Exception("ChipSinica::Decode():WW: some maps were not found for srcID=%d port=%d chan=%d",
		    it->first.u.s.src_id, it->first.u.s.port, it->first.u.s.channel ); 
		      
        size=digits_list.size();
      }
  }

  ////////////////////////////////////////////////////////////////////////////////

  void ChipSinica::Scan(DaqOption &opt)
  {
    if( is_scaned )
      return;
    else
      is_scaned=true;

    if (GetSLink().IsFirstEventInRun())
      {
        Exception("ChipSinica::Scan():WW: first event of run").Print();
        return;
      }

    if( IsDebugMode() )
      ScanDebugMode(opt);
    else{
      throw Exception("ChipSinica::Scan(): Wrong Mode");
      exit(1);
    }
  }

  ////////////////////////////////////////////////////////////////////////////////

  // This code comes from Wolfgang Kastaun.
  void ChipSinica::ScanDebugMode(DaqOption &opt)
  {
    uint32 sourceid=GetSLink().GetSourceID();
  
    //Debug Mode Checks
    uint32 numscanned=0;

    uint8 port=-1;    
    uint16 triggertime=-1;  
    
    uint8 headertemperature = -1; uint8 trailertemperature = -1;  // temperature is saved in header and trailer dont ask why
    double FEMtemperature = -1;

    const uint32 *pheader = NULL;
    const uint32 *ptrailer = NULL;

    bool good=false;  // Will be true if data is good
    bool bad=false;   // Will be true if data is good

    for( const uint32 *p=GetDataStart(); p<GetDataEnd(); p++ )
      {
        // Debugmode Data/header/trailer
        HeaderDbg header=*p;    // Can be a trailer also...
	TrailerDbg trailer=*p;
	numscanned++;

	// --- check if header or trailer
	if(header.IsHeader())
	  {
	    if(header.IsTBOBitSet())
	      {
		AddError(DaqError(DaqError::TDC_TBO2,
				  DaqError::SOURCE_ID, sourceid,
				  DaqError::PORT,      header.GetPort()),opt);
	      }
	    port = header.GetPort();
	    headertemperature = header.GetTemperature();
	    if(pheader == NULL) pheader = p;
	    numscanned = 1; // go back down to 1, start to check if we have a full 34 word package header+data+trailer
	  }
	if(trailer.IsTrailer())
	  {
	    triggertime = trailer.GetData();
	    trailertemperature = trailer.GetTemperature();
	    if(ptrailer == NULL) ptrailer = p;
	    if(!pheader){  // Bad Data. give bad flag and add error. Should have header before trailer
	      AddError(DaqError(DaqError::TDC_ERR20,
				DaqError::SOURCE_ID, sourceid,
				DaqError::PORT,      port),opt);
	      bad = true;
	    }
	  }
	// --- header or trailer check end

	// If we have a header and trailer, we loop between them and read data
	if(pheader && ptrailer && numscanned <= FEM_DATA_SIZE_MAX) {
	  // Debugmode Data 
	  FEMtemperature = (headertemperature*16+trailertemperature)*16*0.03125;
	  if(FEMtemperature > FEM_TEMPERATURE_MAX){
	    AddError(DaqError(DaqError::TDC_TEMPERATURE,
			      DaqError::SOURCE_ID, sourceid,
			      DaqError::PORT,      port),opt);
	  }
   
	  // have to loop through again because trailer has trigger info...
	  for( const uint32 *pp=pheader+1; pp<ptrailer; pp++ ){
	    DataDbg data_word=*pp;
	    if( !event->IsBadDataSource(GetSourceID(),port) )
	      AddData( data_word.GetData(), data_word.GetMode(), data_word.GetChannel(), port, triggertime, FEMtemperature);
	  }
	  good = true;
	}
	
	// If we scanned more than allowed size without header or trailer,
	// then something is wrong
	if(numscanned == FEM_DATA_SIZE_MAX && !good){ 
	  if(pheader && !ptrailer)
	    AddError(DaqError(DaqError::TDC_ERR21,
			      DaqError::SOURCE_ID, sourceid,
			      DaqError::PORT,      port),opt);
	  if(!pheader && !ptrailer)
	    AddError(DaqError(DaqError::TDC_ERR22,
			      DaqError::SOURCE_ID, sourceid,
			      DaqError::PORT,      port),opt);
	  bad = true;
	}
	
	// if 
	if(good || bad){
	  pheader = NULL; ptrailer = NULL;
	  port = -1; triggertime = -1;
	  headertemperature = -1; trailertemperature = -1; FEMtemperature = -1;
	  numscanned = 0;
	  good = false; bad = false;
	}

      } // End of the data loop
  }

  ////////////////////////////////////////////////////////////////////////////////

  void ChipSinica::AddData(uint16 hittime,uint16 mode,uint16 channel,uint8 port, uint16 triggertime, uint16 temperature)
  {
    DataID data_id(GetSLink().GetSourceID(),mode,port,channel);  // dont know trigger time yet
    Hit hit(hittime,triggertime, temperature);
    data_all.push_back( pair<DataID,Hit>(data_id,hit) );
  }

  ////////////////////////////////////////////////////////////////////////////////

  void ChipSinica::Print(ostream &o,const string &prefix) const
  {
    Chip::Print(o,prefix+"ChipSinica ");
    for( const uint32 *d=GetDataStart(); d<GetDataEnd(); d++ )
      {
	if( IsDebugMode() )
	  {
	    HeaderDbg  h(*d);
	    TrailerDbg  t(*d);
	    if( h.IsHeader() )
	      h.Print(o,prefix);
	    else if( t.IsTrailer() )
	      t.Print(o,prefix);
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

  void ChipSinica::DataDbg::Print(ostream &o,const string &prefix) const
  {
    o << prefix << GetName() << "  ";
    char s[222];
    sprintf(s,"data16=%5u  FEMchannel=%2u  mode=%2u\n", GetData(),GetChannel(),GetMode() );
    o<<s;
  }

  ////////////////////////////////////////////////////////////////////////////////

  void ChipSinica::HeaderDbg::Print(ostream &o,const string &prefix) const
  {
    // ChipSinica::Data::Print(o,prefix);
    o << prefix << GetName() << "  ";
    char s[222];
    sprintf(s,"Port=%2u  time=%3u  event=%2u  tbo=%1u  temp=%2u\n",GetPort(), GetData(), GetEventNumber(), IsTBOBitSet(), GetTemperature());
    o<<s;
  }

  ////////////////////////////////////////////////////////////////////////////////

  void ChipSinica::TrailerDbg::Print(ostream &o,const string &prefix) const
  {
    // ChipSinica::Data::Print(o,prefix);
    o << prefix << GetName() << "  ";
    char s[222];
    sprintf(s,"TFtime=%3u  temp=%2u\n", GetData(),GetTemperature());
    o<<s;
  }

  ////////////////////////////////////////////////////////////////////////////////

  ChipSinica::Map::Map(const ObjectXML &o)
    :   Chip::Map(o),
	mode('?'),
	geo_id(-1),
	port(0),
	channel(0),
	hit(0),
	x(0), y(0),
	type_ul(' '),
	wire_position(0),
	overolling(0),
	time_unit(-1),
	time_reference(0)
  {
    if( version==0 )
      version=1;

    if( GetName()!="ChipSinica" )
      throw Exception("ChipSinica::Map::Map(): Internal error.");

    if( GetVersion() > 1 )
      throw Exception("ChipSinica::Map::Map(): unknown version %d",GetVersion());

    GetAttribute("time_reference",time_reference);

    // If "time_unit" has been set in a mapping file, then use it ....
    if( NULL==GetAttribute("time_unit",time_unit) )
      {
	SetAttribute("time_unit","1.0714492"); 
        GetAttribute("time_unit",time_unit);
      }
    else
      GetAttribute("time_unit",time_unit);

    // If "overolling" has been set in a mapping file, then use it ....
    if( NULL==GetAttribute("overolling",overolling) )
      {
	SetAttribute("overolling","64239"); 
        GetAttribute("overolling",overolling);
      }
    else
      GetAttribute("overolling",overolling);

    if( IsOption("tdc"))
      {
        string name;
        int chanN;
        istringstream s(dec_line.c_str());

	mode = 't';
 
        if( GetVersion()==1 )
	  {
	    s >> name >> source_id >> port >> geo_id >> chanF >> chanS >> chanN >> wireF >> wireL >> wireS;

            if( s.fail() )
	      throw Exception("ChipSinica::Map::Map(): bad format in line: %s",map_line.c_str());

            s >> wire_position;
            if( s.fail() )
	      wire_position=0;

            chanL=chanF+(chanN-1)*chanS;

            id=DetID(name);
	  }
      }
    else
      throw Exception("ChipSinica::Map::Map(): unknown option(s) \"%s\" for line \"%s\"",
		      options.c_str(),map_line.c_str());

    Check();
  }

  ////////////////////////////////////////////////////////////////////////////////

  void ChipSinica::Map::AddToMaps(Maps &maps,DaqOption &options) const
  {
    size_t size=maps.size();
    int factor=1;
 
    switch( mode )
      {    
      case 't':
	{
	  switch( GetVersion() )
	    {
            case 1:
	      {
                for( size_t n=0; n<GetChanN(); n++ )
		  {
		    int wire=GetWireF()+n*GetWireS();
		    int chan=GetChanF()+n*GetChanS();
		    DataID data_id(GetSourceID(),1,GetPort(),chan);

		    Digit *digit;
		    digit = new Digit(data_id,GetDetID(),wire,wire_position, time_unit, overolling);
		    digit->SetTimeReference(time_reference);

		    maps.insert( pair<DataID,Digit*>(data_id,digit) );
		  }
                break;
	      }
            default:
	      throw Exception("ChipSinica::Map::AddToMaps(): unknown version %d",GetVersion());
	    }
        
	  break;
	}

      default:
        throw Exception("ChipSinica::Map::AddToMaps(): unknown mode %c",mode);
    
      }

    if( (mode=='t') && (maps.size()-size)!=factor*GetChanN() )
      throw Exception("ChipSinica::Map::AddToMaps(): internal error in Map mode=%c version=%d",mode,GetVersion());

    maps.SetWires( GetDetID(), maps.GetWires(GetDetID())+GetChanN() );
  }


  ////////////////////////////////////////////////////////////////////////////////

  void ChipSinica::Map::Print(ostream &o,const string &prefix) const
  {
    Chip::Map::Print(o,prefix);
    o<<prefix;

    char s[222];
    sprintf(s,"ChipSinica::Map: mode=%c  chip=%d\n",
	    mode,int(port));
    o << s;
  }

  ////////////////////////////////////////////////////////////////////////////////

  void ChipSinica::Digit::Print(ostream &o,const string &prefix) const
  {
    char s[222], s2[44];
    
    sprintf(s2,"wire=%d",GetWire());

    const DataID &d = reinterpret_cast<const DataID &>(GetDataID());

    sprintf(s,"%s%8s hittime=%5d  triggertime=%5d  %s  pos=%d  time_unit=%gns overolling=%d time_decoded=%g DataID=(srcID=%d,port=%d,chan=%d,mode=%d) temperature=%f\n",
	    prefix.c_str(),GetDetID().GetName().c_str(),(int)GetHitTime(),(int)GetTriggerTime(),
	    s2,(int)GetWirePos(),GetTimeUnit(),GetOverolling(),GetTimeDecoded(),
	    unsigned(d.u.s.src_id),unsigned(d.u.s.port),unsigned(d.u.s.channel),unsigned(d.u.s.mode),
	    GetTemperature());
    o << s;
  }

  ////////////////////////////////////////////////////////////////////////////////

  vector<float> ChipSinica::Digit::GetNtupleData(void) const
  {
    vector<float> v;
    v.push_back(GetWire());
    v.push_back(GetHitTime());
    v.push_back(GetWirePos());
    v.push_back(GetTimeUnit());
    v.push_back(GetTimeDecoded());
    return v;
  }

  ////////////////////////////////////////////////////////////////////////////////

  double ChipSinica::TimeDifference(int hittime,int trigger_time,int time_overolling, double time_ref, int cut)  // cut = 30000
  {
    int diff=hittime-trigger_time-time_ref;
 
    // Maximum two iterations
    for( int i=0; fabs(diff)>=cut && i<3; i++ ){
      if( diff<0 )
	diff += time_overolling;
      else
	diff -= time_overolling;
    }

    return diff+time_ref;
  }

  ////////////////////////////////////////////////////////////////////////////////

}
