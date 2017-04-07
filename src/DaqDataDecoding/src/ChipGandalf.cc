/*!
    \file ChipGandalf.cc
    \author Boris Iven
*/

#include <cassert>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <utility>

#include "ChipGandalf.h"
#include "DaqError.h"
#include "DaqOption.h"
#include "DaqEvent.h"
#include "utils.h"

namespace CS {

using std::pair;
using std::istringstream;

////////////////////////////////////////////////////////////////////////////////

ChipGandalf::ChipGandalf(void const * const buf,bool copy_buf,DaqOption &opt,DaqEvent &ev)
:   Chip(buf,copy_buf,opt,ev)
{
	//std::cout << "INFO: ChipGandalf::CTOR" << std::endl; //debug
    Clear();
    //Dump(cout, "INFO: "); //debug
}

////////////////////////////////////////////////////////////////////////////////

void ChipGandalf::Decode(const Maps &maps,Digits &digits_list,DaqOption &opt)
{
	//std::cout << "INFO: ChipGandalf::Decode" << std::endl;//debug
    // Scan the chip data if it was not done yet.

    if( !is_scaned )
        Scan(opt);

	
	for( std::list<DigitGADC*>::iterator it=um_digits.begin(); it!=um_digits.end(); ++it )
	{
		
		typedef Maps::const_iterator m_it; // Create a short name for map's iterator
		const pair<m_it,m_it> m_range = maps.equal_range( (*it)->GetDataID() ); // all maps with given Data ID
		
		if( m_range.first==m_range.second ) {
			printf("ChipGandalf::Decode(): No map found for det/chan combination: src %4d, spill/event %4d/%4d \n",GetSLink().GetSourceID(),GetSLink().GetSpillNumber(),GetSLink().GetEventNumber());
			continue;
		}

		for( m_it c = m_range.first; c != m_range.second; ++c )
		{
				const DigitGADC* map_digit = dynamic_cast<DigitGADC*>(c->second);
				if( map_digit )
				{
					//clone the mapped digit
					DigitGADC* final_digit = new DigitGADC(*map_digit);
          
          //copy operation mode and data from unmapped digit
					
					final_digit->setOpMode( (*it)->getOpMode() );
					final_digit->setPulseDet( (*it)->getPulseDet() );
					final_digit->setBaseLErr( (*it)->getBaseLErr() );
					final_digit->setBaseLine( (*it)->getBaseLine() );
					final_digit->setMaxAmplitude( (*it)->getMaxAmplitude() );
					final_digit->setIntegral( (*it)->getIntegral() );
					final_digit->setFrameTime( (*it)->getFrameTime() );
					final_digit->setCoarseTimeMSB( (*it)->getCoarseTimeMSB() );
					final_digit->setCoarseTimeLSB( (*it)->getCoarseTimeLSB() );
					final_digit->setHiResTime( (*it)->getHiResTime() );
					final_digit->setThreshold( (*it)->getThreshold());
					final_digit->setDelay( (*it)->getDelay());
					final_digit->setFracfact( (*it)->getFracfact());

					// TODO: user constructor
					final_digit->SetTimeUnit( map_digit->GetTimeUnit() );
					final_digit->SetDetID(map_digit->GetDetID());
					final_digit->setFrameData( (*it)->getFrameData() );
					
					//add digit to output list
					digits_list.insert(pair<DetID,Digit*>(final_digit->GetDetID(), final_digit));
				}
				else	
					throw Exception("ChipGandalf::Decode(): Internal error");
		    	
		}
	}
  for( std::list<DigitGADC*>::iterator it=um_digits.begin(); it!=um_digits.end(); ++it )    delete(*it);  
  um_digits.clear();
}

///////////////////////////////////////////////////////////////////////////////

void ChipGandalf::Scan(DaqOption &opt)
{
	//std::cout << "INFO: ChipGandalf::Scan" << std::endl;//debug
	
	//!skip first event of run/spill events
	
	if ( isFEOR() || isFEOS() )
	{
		//std::cout << "INFO: ChipGandalf::Scan: Skipped FEOR/FEOS event" << std::endl;//debug
		return;
	}
		
	// check if the format bit is set
	if ( isUndefined() )
	{
		throw Exception("ChipGandalf::Scan: Format not specified !");
	}
	else if ( isTDCReadout() && isGandalfTDC() )
	{
		//TODO
		throw Exception("ChipGandalf::Scan: TDC mode not implemented yet!");
	}
	else if ( (isADCReadout() && (isGandalfADC() || isGandalfScaler())))
	{
		if ( isGandalfADC())
		{
			if (isADCDbgMode() || isADCProcMode())
				scanADC( opt , isADCIlm() );
			else
				throw Exception("ChipGandalf::Scan: Unsupported ADC mode for Gandalf-ADC");
		}
		else //!Gandalf-Scaler
		{
			//TODO
			throw Exception("ChipGandalf::Scan: Gandalf-Scaler mode not implemented yet!");
		}
	}
	else
		throw Exception("ChipGandalf::Scan: Unsupported mode");

}

void ChipGandalf::scanADC( DaqOption &opt, bool interleaved )
{
	
	bool done=false;

	// error message;
	// this string is updated after each step so that it always holds the current possible error
	string msg="unknown";

	const uint32* cdw = GetDataStart();


	// outer loop over datawords; if we have processed data, this loop should be finished after one round
	while( cdw < GetDataEnd() )
	{

		// the digit containing the frame in case we are in dbg mode
		DigitGADC* dbgDigit = NULL;

		// generate digits for each hit
		vector<DigitGADC*> hits;

		const GADCHeadTrail* trail 	= NULL;
		const GADCHeadTrail* head	= NULL;


		//////////////// in dbg mode, first word is header followed by frame data
		if (isADCDbgMode()) {

			head = reinterpret_cast<const GADCHeadTrail*>(cdw);

			msg="could not retrieve channel header";

			if(!head) goto error;


//			//next data word
//			cdw++;
//			msg = "no dataword after header;";
//			if (cdw >= GetDataEnd() ) goto error;

			dbgDigit = new DigitGADC( DataID(GetSourceID(), 0, head->ChID), DetID(""), interleaved? GADC_DEBUG_IL : GADC_DEBUG );
			msg = "could not create new DigitGADC";
			if (!dbgDigit) goto error;

			//next data word
			cdw++;
			std::stringstream tmp;
			tmp << "no dataword after channel header; expecting " << head->WindowLength  << " frames now";
			msg = tmp.str();
			if( cdw >= GetDataEnd() ) goto error;

			//now we expect n(WindowLength) frame words
			const GADCFrameWord* cfw;

			for ( uint32 fw = 0; fw < head->WindowLength; fw++ )
			{
				cfw = reinterpret_cast<const GADCFrameWord*>(cdw);

				msg="could not retrieve frame word";
				if( cfw )
				{
					//TODO: do some integrity checks
					//std::cout << "INFO: ChipGandalf::scanADCDbg: adding frame data #" << 2*fw << " and #" << 2*fw+1 << " to digit..." << std::endl;//debug

					//Add data do digit
					dbgDigit->addFrameData(cfw->DataHigh);
					dbgDigit->addFrameData(cfw->DataLow);
				}
				else goto error;

				//next data word
				cdw++;
				std::stringstream tmp;
				tmp << "missing frame word; expecting " << head->WindowLength-fw << " more frame words";
				msg = tmp.str();
				if( cdw >= GetDataEnd() ) goto error;
			}
		}

		// in debug mode, we could be finished here with the channel data if there are no hits
		trail = reinterpret_cast<const GADCHeadTrail*>(cdw);

		while( (isADCDbgMode() ? trail->Head_Trail != 1 : cdw < GetDataEnd()) )
		{

			// digit where decoded info is stored
			DigitGADC* curDigit=NULL;


			msg="could not retrieve channel/baseline/integral dataword";

			// Channel Baseline Integral
			const GADCProcChIntWord* cpciw = reinterpret_cast<const GADCProcChIntWord*>(cdw);
			if( cpciw )
			{
				//TODO: do some integrity checks

				curDigit = new DigitGADC( DataID(GetSourceID(), 0, cpciw->ChID), DetID(""), interleaved? (isADCDbgMode()?GADC_DEBUG_IL:GADC_PROC_IL) : (isADCDbgMode()?GADC_DEBUG:GADC_PROC) );
				curDigit->setBaseLine(cpciw->Baseline);
				curDigit->setIntegral(cpciw->Integral);
			}
			else goto error;

			msg = "could not create new DigitGADC";
			if (!curDigit)
				goto error;

			//next data word
			cdw++;

			msg = "no dataword after channel/baseline/integral; expecting 2 more datawords";
			// exit if curDigit still NULL or no next word exists
			if (cdw >= GetDataEnd() ) goto error;


			msg = "could not retrieve ct_msb/amplitude dataword";

			//CT_MSB Amplitude
			const GADCProcCTimeAmpWord* cptaw = reinterpret_cast<const GADCProcCTimeAmpWord*>(cdw);
			if( cptaw )
			{
				//TODO: do some integrity checks

				//Add data do digit
				curDigit->setCoarseTimeMSB(cptaw->CTdataMSB);
				curDigit->setMaxAmplitude(cptaw->Amplitude);

			}
			else goto error;

			//next data word
			cdw++;

			msg = "no dataword after ct_msb/amplitude; expecting 1 more dataword";
			if( cdw >= GetDataEnd() ) goto error;

			msg = "could not retrieve ct_lsb/HighRes dataword";

			//CT_LSB HighRes
			const GADCProcCTimeHResWord* cpthw = reinterpret_cast<const GADCProcCTimeHResWord*>(cdw);
			if( cpthw )
			{
				//TODO: do some integrity checks

				//Add data do digit
				curDigit->setCoarseTimeLSB(cpthw->CTdataLSB);
				curDigit->setHiResTime(cpthw->HResTime);
			}
			else goto error;
			

			if (isADCDbgMode())
			{
				//next data word
				cdw++;

				msg = "no dataword after ct_lsb/HighRes; expecting 1 more dataword";
				if( cdw >= GetDataEnd() ) goto error;

				msg = "could not retrieve frameTime dataword";

				//CT_MSB Amplitude
				const GADCDbgFTimeWord* cpftw = reinterpret_cast<const GADCDbgFTimeWord*>(cdw);
				if( cpftw )
				{
					//TODO: do some integrity checks

					//Add data do digit
					curDigit->setFrameTime(cpftw->FrameTime);
					curDigit->setDelay(cpftw->Delay);
					curDigit->setThreshold(cpftw->Threshold);
					curDigit->setFracfact(cpftw->FracFact);
				}
				else goto error;
			}

			// digit is ready to store
			hits.push_back(curDigit);

			//next data word
			cdw++;

			// in debug mode, the last word has to be the trailer
			if (isADCDbgMode()) trail = reinterpret_cast<const GADCHeadTrail*>(cdw);
		}

		if(isADCDbgMode())
		{
			msg = "could not retrieve channel trailer";
			if (trail)
			{

				msg = "channel trailer differs from channel header";
				if ( (*head) != (*trail) )
					goto error;

				if (hits.size()==0) dbgDigit->setPulseDet(0);
				um_digits.push_back(dbgDigit);

				//next data word
				cdw++;

				msg = "reached trailer word, but there are still datawords left";

				// now, all data should be read
				//if( cdw < GetDataEnd() ) goto error;

			}
			else goto error;
		}

		for(vector<DigitGADC*>::iterator hit=hits.begin();hit!=hits.end();++hit)
		{
			if ( (*hit)->getMaxAmplitude()>10 ) um_digits.push_back((*hit));
		  else delete(*hit);
    }
	}

	// if we get until here, everything went through without an error, we can store the digits
	done = true;

	error:
	// we had an error somewhere :(
	if (!done) {
		std::stringstream errormsg;
		errormsg << "ERROR: ChipGandalf::scanADC: " << msg << ";   Source ID: " <<GetSourceID() << "   Spill/Event No: " << GetSLink().GetSpillNumber() << "/"<< GetSLink().GetEventNumber() << std::endl;
		throw Exception(errormsg.str().c_str());
	}
}

////////////////////////////////////////////////////////////////////////////////

void ChipGandalf::Print(ostream &o,const string &prefix) const
{
	std::cout << "ChipGandalf::Print" << std::endl;//debug

}

////////////////////////////////////////////////////////////////////////////////


ChipGandalf::Map::Map(const ObjectXML &o)
:   Chip::Map(o),
	m_mapMode(-1),
	m_port(0),
	m_geoID(0),
	m_x(0),
	m_y(0),
	m_z(0),
	m_r(0),
	m_phi(0),
	time_unit(-1)
{
#ifdef debug
	std::cout << "INFO: ChipGandalf::Map CTOR" << std::endl;//debug
#endif
	int32 chanN;

    if( version==0 )
        version=1;

    if( GetName()!="ChipGandalf" )
        throw Exception("ChipGandalf::Map::Map(): Internal error.");

    if( GetVersion()<1 || GetVersion()>1  )
        throw Exception("ChipGandalf::Map::Map(): unknown version %d",GetVersion());

	string name;
	
    // If "time_unit" has been set in a mapping file, then use it ....
    if( NULL==GetAttribute("time_unit",time_unit) )
    {
        SetAttribute("time_unit","0.12892312");
        GetAttribute("time_unit",time_unit);
    }


    if( IsOption("adc") )
    {
#ifdef debug
    	std::cout << "INFO: ChipGandalf::Map ADC mode" << std::endl;//debug
#endif
    	if( IsOption("xyz") )
    	{
#ifdef debug
    		std::cout << "INFO: ChipGandalf::Map: XYZ-mapping" << std::endl;//debug
#endif
    		m_mapMode = MM_XYZ;
		    istringstream s(dec_line.c_str());

		    s >> name >> source_id >> m_port >> m_geoID >> chanF >> m_x >> m_y >> m_z;
		    
		    chanN = 1;
		    chanL = chanF;
		    
		    wireF = wireL = chanF;
		    wireS = 1;
		    
		    if( s.fail() )
		        throw Exception("ChipGandalf::Map::Map(): bad format (adc) in line: %s",map_line.c_str());
		        
		    //!port will be ignored in adc mode
		    m_port = 0;
		    chanL=chanF+(chanN-1)*chanS;

		    
		    id=DetID(name);
    	
    	}
    	
    	else if( IsOption("rphiz") )
    	{
#ifdef debug
    	    std::cout << "INFO: ChipGandalf::Map: RPHIZ-mapping" << std::endl;//debug
#endif
	   		m_mapMode = MM_RPHIZ;
		    istringstream s(dec_line.c_str());

		    s >> name >> source_id >> m_port >> m_geoID >> chanF >> m_r >> m_phi >> m_z;
		    
		    chanN = 1;
		    chanL = chanF;
		    
		    wireF = wireL = chanF;
		    wireS = 1;
		    
		    if( s.fail() )
		        throw Exception("ChipGandalf::Map::Map(): bad format (adc) in line: %s",map_line.c_str());
		        
		    //!port will be ignored in adc mode
		    m_port = 0;
		    chanL=chanF+(chanN-1)*chanS;

		    
		    id=DetID(name);
    	
    	}
        
        else
        {
#ifdef debug
            std::cout << "INFO: ChipGandalf::Map: normal mapping" << std::endl;//debug
#endif
    	   	m_mapMode = MM_NORMAL;
		    istringstream s(dec_line.c_str());

		    s >> name >> source_id >> m_port >> m_geoID >> chanF >> chanS >> chanN >> wireF >> wireL >> wireS;
		    
		    if( s.fail() )
		        throw Exception("ChipGandalf::Map::Map(): bad format (adc) in line: %s",map_line.c_str());
		        
		    //!port will be ignored in adc mode
		    m_port = 0;
		    chanL=chanF+(chanN-1)*chanS;

		    
		    id=DetID(name);
        }
        
    }
    else if( IsOption("tdc") )
    {
#ifdef debug
       	std::cout << "INFO: ChipGandalf::Map TDC mode" << std::endl;//debug
#endif
    	if( IsOption("xyz") )
    	{
    		std::cout << "INFO: ChipGandalf::Map: XYZ-mapping" << std::endl;//debug
    		m_mapMode = MM_XYZ;
		    istringstream s(dec_line.c_str());

		    s >> name >> source_id >> m_port >> m_geoID >> chanF >> m_x >> m_y >> m_z;
		    
		    chanN = 1;
		    chanL = chanF;
		    
		    wireF = wireL = chanF;
		    wireS = 1;
		    
		    if( s.fail() )
		        throw Exception("ChipGandalf::Map::Map(): bad format (tdc) in line: %s",map_line.c_str());
		    
		    chanL=chanF+(chanN-1)*chanS;
		    id=DetID(name);    	
		}
		else if ( IsOption("rphiz") )
		{
#ifdef debug
    		std::cout << "INFO: ChipGandalf::Map: RPHIZ-mapping" << std::endl;//debug
#endif
			m_mapMode = MM_RPHIZ;
		    istringstream s(dec_line.c_str());

		    s >> name >> source_id >> m_port >> m_geoID >> chanF >> m_r >> m_phi >> m_z;
		    
		    chanN = 1;
		    chanL = chanF;
		    
		    wireF = wireL = chanF;
		    wireS = 1;

		    if( s.fail() )
		        throw Exception("ChipGandalf::Map::Map(): bad format (tdc) in line: %s",map_line.c_str());
		    
		    chanL=chanF+(chanN-1)*chanS;
		    id=DetID(name);
        }
        else
        {
#ifdef debug
       		std::cout << "INFO: ChipGandalf::Map: normal mapping" << std::endl;//debug
#endif
        	m_mapMode = MM_NORMAL;
		    istringstream s(dec_line.c_str());

		    s >> name >> source_id >> m_port >> m_geoID >> chanF >> chanS >> chanN >> wireF >> wireL >> wireS;

		    if( s.fail() )
		        throw Exception("ChipGandalf::Map::Map(): bad format (tdc) in line: %s",map_line.c_str());
		    
		    chanL=chanF+(chanN-1)*chanS;
		    id=DetID(name);        
        }
    }
    else
        throw Exception("ChipGandalf::Map::Map(): unknown option(s) \"%s\" for line \"%s\"",
                         options.c_str(),map_line.c_str());

    Check();//TODO
    
    
}

////////////////////////////////////////////////////////////////////////////////

void ChipGandalf::Map::ReadTTConfig(TriggerTimeConfig &tt_conf) const
{
    // Check for the trigger time decoding options.
    int TT_index=-1;
    double time_unit=0;
    string sss;

    if( NULL!=GetAttribute("TT_index",TT_index) )
    {
        bool ok = GetAttribute("time_unit",    time_unit);
        if(!ok)
            throw Exception("ChipGandalf::Map::Map(): Bad settings for the trigger time");

        TriggerTimeConfig::TTC ttc(id,time_unit,TT_index);
        ttc.srcID = GetSourceID();
        tt_conf.Add(ttc);
    }
}


////////////////////////////////////////////////////////////////////////////////

void ChipGandalf::Map::AddToMaps(Maps &maps,DaqOption &options) const
{
#ifdef debug
	std::cout << "INFO: ChipGandalf::Map::AddToMaps(): start " << std::endl;//debug
#endif

	DataID data_id;
	DigitGADC* digit;

	if( maps.end()!=maps.find(data_id) && !IsMultiDigit() )
	{
	  Print();
	  Exception("ChipADC::Map::AddToMaps(): map already exists").Print();
	}

	if( IsOption("adc") )
	{
		switch( m_mapMode )
		{
			case MM_NORMAL:
				for ( int c = chanF; c <= chanL; c += chanS )
				{
					data_id = DataID(GetSourceID(), m_port, c);
					digit = new DigitGADC( data_id, GetDetID(), GADC_UNKNOWN, CalcWire(c) );
					// TODO: discriminate mapmode with/without time_unit
					digit->SetTimeUnit(time_unit);
					maps.insert( pair<DataID,Digit*>(data_id,  digit) );
#ifdef debug
					std::cout << "INFO: ChipGandalf::Map::AddToMaps(): mapped digit from channel: " << c << " to wire: " << CalcWire(c) << std::endl;//debug
#endif
				}
				break;
			case MM_XYZ:
				data_id = DataID(GetSourceID(), m_port, chanF);
				digit = new DigitGADC( data_id, GetDetID(), GADC_UNKNOWN, chanF, m_x, m_y, m_z );
				maps.insert( pair<DataID,Digit*>(data_id, digit) );
				break;
			case MM_RPHIZ:
				data_id = DataID(GetSourceID(), m_port, chanF);
				digit = new DigitGADC(data_id, GetDetID(), GADC_UNKNOWN, chanF, m_r, m_phi, m_z);
				maps.insert( pair<DataID,Digit*>(data_id, digit) );
				break;
			default:
				throw Exception("ChipGandalf::Map::AddToMaps(): unknown mapping mode %i", m_mapMode ); break;
		}
		
	    ReadTTConfig(options.GetTTConfig());

	}
	else //TDC
	{
		throw Exception( "ChipGandalf::Map::AddToMaps(): TDC mapping not implemented yet!" );
	}
	
	
	maps.SetWires( GetDetID(), maps.GetWires(GetDetID())+1 );
#ifdef debug
	Print();//debug
#endif
}


////////////////////////////////////////////////////////////////////////////////

void ChipGandalf::Map::Print(ostream &o,const string &prefix) const
{
  Chip::Map::Print(o,prefix);
  o<<prefix;

  char s[222];
  sprintf(s,"ChipGandalf::Map: port=%d\n", int(m_port));
  o << s;
}

////////////////////////////////////////////////////////////////////////////////

void ChipGandalf::Digit::Print(ostream &o,const string &prefix) const
{
  o << "ChipGandalf::Digit::Print()";
}


void ChipGandalf::DigitGADC::Print(ostream &o,const string &prefix) const
{
  o << "ChipGandalf::DigitGADC::Print()";
}

/**
 *	If you want to use GandalfDigit in DST, discriminate betweene digit types by vector size
 *	
 *	v.size() == 8: normal digit
 *	v.size() == 9: debug digit; amplitude and integral only meaningfull, if (integral != 0xffee)
 *	v.size() >  9: frame digit; amplitude/integral/times have no meaning. A frame digit is always
 *					followed by a debug digit, which holds the time information.
 *
**/
vector<float> ChipGandalf::DigitGADC::GetNtupleData(void) const
{
  vector<float> v;
  
	// Write processed data for all digits
	v.push_back( getChannel() );
	v.push_back( m_maxAmplitude );
	v.push_back( m_integral );
	v.push_back( GetTimeDecoded() );
	v.push_back( m_coarseTimeMSB );
	v.push_back( m_coarseTimeLSB );
	v.push_back( m_hiResTime );
	v.push_back( GetTimeUnit() );
	
	// check if this is a debug digit
	if (getOpMode()==CS::ChipGandalf::GADC_DEBUG || getOpMode()==CS::ChipGandalf::GADC_DEBUG_IL) {
		
		// write frametime
		v.push_back( getFrameTime() );

		// debug digit with frame info, add number of samples and 
		// the samples to the vector
		if (getNumSamples()>0) {
			v.push_back( getNumSamples() );
			for (int i = 0; i<getNumSamples();i++) {
        			v.push_back( getSample(i) );
	        	}
		}
	}
  return v;
}

////////////////////////////////////////////////////////////////////////////////

}
