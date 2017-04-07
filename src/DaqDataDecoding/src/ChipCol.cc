#include <cstdlib>
#include <cstdio>

#include "DaqEvent.h"
#include "ChipCol.h"

namespace CS {

////////////////////////////////////////////////////////////////////////////////

void ChipCol::Scan(DaqOption &opt)
{
    if( is_scaned )
        return;
    else
        is_scaned=true;

    uint32 source=GetSLink().GetSourceID();
    HeaderTrailer ht(0,0), header(0,0);
    bool hashdr = false;
    enum {DATA,HEADER,TRAILER,ERROR} lastd=TRAILER;
	
    if (GetSLink().IsFirstEventInRun())
        return;
  
    for( const uint32 *p=GetDataStart(); p<GetDataEnd(); p++ )
    {
        if( ErrorWord(*p,ErrorMarker,source,*this,opt).IsError() )
        {
            lastd=ERROR;
            continue;
        }

        ht=*p;

        if (ht.IsData())
        {
            if (lastd==TRAILER) 
                AddError(DaqError(DaqError::ADC_DATA_WO_HEADER,DaqError::SOURCE_ID,source),opt);
            else
                if (hashdr)
                {
                    // Accept data only if srcID and geoID are not ignored
                    if( !event->IsBadDataSource(GetSourceID(),header.GetGeoID()) )
                        all_data.push_back( std::pair<HeaderTrailer,uint32>(header,*p) ); // Add header+data
                }
            lastd=DATA;
        }
        else
        {
            if (ht.IsHeader())
            {
                RegisterGeoID(ht.GetGeoID());

                if ((ht.GetGeoID()==0)&& (lastd==HEADER))
                {
                    if (ht.GetGeoID()!=header.GetGeoID())
                        AddError(DaqError(DaqError::ADC_W_GEOID_EW,
                                          DaqError::SOURCE_ID,source,
                                          DaqError::GEO_ID,header.GetGeoID(),
                                          DaqError::GEO_ID,ht.GetGeoID()),opt);
                }
                else
                {
                    header=*p;
                    hashdr=(lastd!=ERROR);
                    if (lastd==DATA || lastd!=TRAILER) 
                        AddError(DaqError(DaqError::ADC_NO_TRAILER,
                                          DaqError::SOURCE_ID,source, DaqError::GEO_ID,ht.GetGeoID()),opt);
                    lastd=HEADER; 
                    CheckEventNumber(ht,opt);
                }
            }
            else
            {
                if (lastd==TRAILER) 
                    AddError(DaqError(DaqError::ADC_TRAILER_WO_HEADER,
                                      DaqError::SOURCE_ID,source, DaqError::GEO_ID,ht.GetGeoID()),opt);
                lastd=TRAILER; hashdr=false;
                CheckEventNumber(ht,opt);
                if (hashdr && (ht.GetGeoID()!=header.GetGeoID()))
                    AddError(DaqError(DaqError::ADC_W_GEOID_TRAILER,
                                      DaqError::SOURCE_ID,source,
                                      DaqError::GEO_ID,header.GetGeoID(),
                                      DaqError::GEO_ID,ht.GetGeoID()),opt);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void ChipCol::AddData(uint32 data)
{
    if( !CanModify() )
        throw Exception("ChipCol::Add(): you don't have rights to modify the chip data.");
 
    if( (data>>31)!=1 )
        throw Exception("ChipCol::Add(Data): This is not data!");
  
    const HeaderTrailer &ht = *reinterpret_cast<HeaderTrailer*>(buf+GetLength()/4-1);
    if( !(GetDataEnd()>GetDataStart()) && ht.IsTrailer() )
        throw Exception("ChipCol::Add(): I do not expect trailer here");
  
    Chip::Add(&data,sizeof(data)/4);
}

////////////////////////////////////////////////////////////////////////////////

void ChipCol::CheckEventNumber(const HeaderTrailer &ht,DaqOption &opt)
{
    if( GetSLink().GetEventNumber()==ht.GetEventNumber() )
        return;  // good.

    DaqError::Type t = ht.IsTrailer() ? DaqError::ADC_WR_EV_TRAILER : DaqError::ADC_WR_EV_HEADER;
    AddError(DaqError(t,DaqError::SOURCE_ID,GetSourceID(),
                        DaqError::GEO_ID,ht.GetGeoID(),
                        DaqError::SL_EVENT_N,GetSLink().GetEventNumber(),
                        DaqError::EVENT_N,ht.GetEventNumber()),opt);

//   int diff = GetSLink().GetEventNumber() - ht.GetEventNumber();
//   char s_diff[111], sign=diff>0?'+':'-';
//   if( abs(diff)<=3 )
//     sprintf(s_diff,"%+d",diff);
//   else
//   if( abs(diff)<=9 )
//     sprintf(s_diff,"%c4..%c9",sign,sign);
//   else
//   if( abs(diff)<=99 )
//     sprintf(s_diff,"%c10..%c99",sign,sign);
//   else
//     sprintf(s_diff,"%c100..%cinf",sign,sign);
// 
//   AddError(DaqError(Exception("ADC ERROR:%s with wrong event# (src=%d geoID=%d (SLinkEvent-HeaderEvent)=%s)",
//                     ht.IsTrailer()?"Trailer":"Header", GetSourceID(), ht.GetGeoID(), s_diff)));
}

////////////////////////////////////////////////////////////////////////////////

uint32 ChipCol::HeaderTrailer::GetGeoID(void) const
{
    if( !IsData() )
        return d.s.geoID;

    throw Exception("ChipCol::HeaderTrailer::GetGeoID(): this is data.");
}

////////////////////////////////////////////////////////////////////////////////

uint32 ChipCol::HeaderTrailer::GetEventNumber(void) const
{
    if( !IsData() )
        return d.s.event;

    throw Exception("ChipCol::HeaderTrailer::GetEventNumber(): this is header/trailer.");
}

////////////////////////////////////////////////////////////////////////////////

} // namespace CS
