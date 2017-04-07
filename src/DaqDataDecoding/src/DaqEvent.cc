#include "config.h"

#include <algorithm>
#include <iterator>
#include <cassert>

#include "DaqError.h"
#include "DaqEvent.h"
#include "DaqOption.h"
#include "Chip.h"
#include "ChipF1.h"
#include "TriggerTime.h"
#include "DateEquipment.h"
#include "utils.h"

using namespace std;

namespace CS {

////////////////////////////////////////////////////////////////////////////////

void DaqEvent::ReadChips(DaqEventCallMe f)
{
    throw "DaqEvent::ReadChips() do you really need me?! The code does not look good!";
    //DaqOption opt(ObjectXML("DaqOption"));
    //ReadChips(opt,f);
}

////////////////////////////////////////////////////////////////////////////////

string DaqEvent::GetTimeStr(void) const
{
  char s[111];
  time_t t = GetHeader().GetTime().first;
  strftime(s,sizeof(s),"%Y/%m/%d-%T",gmtime(&t));
  sprintf(s+strlen(s)," + %g sec",GetHeader().GetTime().second/1000000.);
  return s;
}

////////////////////////////////////////////////////////////////////////////////

bool DaqEvent::Header::HaveSubEvents(void) const
{
    if( GetVersion()<0xffff )
        return (GetHeaderOld().detectorId[0] & SUPER_EVENT_MASK)
                ||
               (GetHeaderOld().detectorId[2] & SUPER_EVENT_MASK);
    else
    {
        // (event->eventTypeAttribute)[2] & (1<<((ATTR_SUPER_EVENT)&0x1f))
        uint32 n = GetHeader36().typeAttribute[HeaderDATE_36::ATTR_SUPER_EVENT/32];
        return n & (1<<(HeaderDATE_36::ATTR_SUPER_EVENT&0x1f));
    }
}

////////////////////////////////////////////////////////////////////////////////

uint32 DaqEvent::Header::GetTrigger(void) const
{
    if( GetVersion()<0xffff )
        return GetHeaderOld().typeAttribute[1];
    else
    {
        const HeaderDATE_36 &h = GetHeader36();
        return (h.event_trigger_pattern[0] >> 1) + (h.event_trigger_pattern[1] << 31);
    }
}

////////////////////////////////////////////////////////////////////////////////

bool DaqEvent::Header::operator == (const DaqEvent::Header &h) const
{
    return
        GetVersion()            == h.GetVersion()               &&
        GetMagic()              == h.GetMagic()                 &&
        GetBurstNumber()        == h.GetBurstNumber()           &&
        GetRunNumber()          == h.GetRunNumber()             &&
        GetEventNumberInRun()   == h.GetEventNumberInRun()      &&
        GetEventNumberInBurst() == h.GetEventNumberInBurst();
}

////////////////////////////////////////////////////////////////////////////////

uint32 DaqEvent::Header::GetLength(void) const
{
    const uint32 l = GetVersion()<0xffff ? GetHeaderOld().size : GetHeader36().size;
    if( l==0 )
        throw Exception("DaqEvent::GetLength(): can not be zero!");
    return l;
}

////////////////////////////////////////////////////////////////////////////////

void DaqEvent::Header::Print(const std::string &prefix) const
{
    printf("%sversion = %d\n",prefix.c_str(),GetVersion());
    printf("%sevent length = %d bytes, %g words\n",prefix.c_str(),GetLength(),GetLength()/4.);
    printf("%sEvent header length = %d\n",prefix.c_str(),GetHeaderLength());
    printf("%ssub-events = %d\n",prefix.c_str(),HaveSubEvents());
    printf("%srun/burst = %d/%d\n",prefix.c_str(),GetRunNumber(),GetBurstNumber());
    printf("%sevent in: burst/run = %d/%d\n",prefix.c_str(),GetEventNumberInBurst(),GetEventNumberInRun());
    //printf("%strigger number %d\n",prefix.c_str(),GetTriggerNumber());
    printf("%sevent type %d\n",prefix.c_str(),GetEventType());
    //printf("%strigger mask %d\n",prefix.c_str(),GetEventType());
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

void DaqEvent::SetTriggerMask(const DaqOption &options)
{
    for( set<uint16>::const_iterator bit=options.GetTriggerBitsIgnore().begin();
         bit!=options.GetTriggerBitsIgnore().end(); bit++ )
        trigger_mask &= ~(1<<*bit);
}

////////////////////////////////////////////////////////////////////////////////

void DaqEvent::ReadChips(DaqOption &opt,DaqEventCallMe call_me)
{
    SetTriggerMask(opt);

    if( HaveSubEvents() )
    {
        // Loop on all subevents
        for( size_t pos=Header(buf).GetHeaderLength(); pos!=GetLength(); )
        {
            if( pos>=GetLength() )
            {
                Exception("DaqEvent::ReadChips(): Event is damaged. (a) pos>=GetLength():  pos=%d GetLength()=%d",pos,GetLength()).Print();
                AddError(DaqError(DaqError::EVENT_CORRUPTED0));
                return;
            }

            if( GetLength()-pos<sizeof(Header) || pos<=0 )
            {
                Exception("DaqEvent::ReadChips(): Event is damaged.  GetLength()-pos=%d pos=%d.",GetLength()-pos,pos).Print();
                AddError(DaqError(DaqError::EVENT_CORRUPTED0));
                return;
            }

            DaqEvent d(buf+pos);  // New subevent
            d.parent = this;

            if( GetHeader() != d.GetHeader() )
            {
                AddError(DaqError(DaqError::EVENT_SUBEVENT_HEADER_DIFF));
                return;
            }

            d.SetEvent1Run(GetEvent1Run());
            d.online_filter=online_filter;

            if( d.GetLength()>GetLength()-pos )
            {
                Exception("DaqEvent::ReadChips():  Event is damaged. l=%d  GetLength()=%d  pos=%d   GetLength()-pos=%d",
                           d.GetLength(),GetLength(),pos,GetLength()-pos).Print();
                AddError(DaqError(DaqError::EVENT_CORRUPTED0));
                return;
            }

            // Read data from this subevent
            d.ReadChips(opt,call_me);  // call this function again.

            assert( d.errors.errors.empty() );

            // Move chips to the main event
            chips.insert(chips.end(),d.GetChips().begin(),d.GetChips().end());
            for( vector<Chip*>::iterator c=d.GetChips().begin(); c!=d.GetChips().end(); c++ )
                (*c)->SetDaqEvent(*this);

            d.GetChips().clear();
            pos += d.GetLength();   // go to a next subevent

            online_filter = d.online_filter;
        }
    }
    else
    {
        if( (GetHeader().GetEventType()&Header::EVENT_TYPE_MASK)==PHYSICS_EVENT ||
            (GetHeader().GetEventType()&Header::EVENT_TYPE_MASK)==CALIBRATION_EVENT ||
            (GetHeader().GetEventType()&Header::EVENT_TYPE_MASK)==END_OF_BURST )
        {
            for( size_t pos=Header(buf).GetHeaderLength(); pos!=GetLength(); )
            {
                const DateEquipment eq((uint32*)(buf+pos),GetHeader().GetVersion());

                if( pos>=GetLength() || GetLength()-pos<eq.SizeOf() )
                {
                    Exception("DaqEvent::ReadChips(): bad pos:  pos>=GetLength():  pos=%d GetLength()=%d",pos,GetLength()).Print();
                    AddError(DaqError(DaqError::EVENT_CORRUPTED0));
                    return;
                }

                if( GetHeader().GetEventType()&Header::EVENT_SWAPPED )
                    throw Exception("DaqEvent::ReadChips(): event swapping is disabled!");

                if( pos+eq.GetFullLength() > GetLength() )
                {
                    Exception("DaqEvent::ReadChips(): (d) Chip is out of the event: type=%d ChipLen=%d   EventLen=%d  diff=%d",
                                     int(GetHeader().GetEventType()), eq.GetFullLength(),GetLength(),GetLength()-pos-eq.GetFullLength()).Print();
                    AddError(DaqError(DaqError::EVENT_CORRUPTED0));
                    return;
                }

                // ----------------------------------------------------------
                // Check that DATE equipment size and SLink size are correct.
                // ----------------------------------------------------------

                //Full chip data length (header + data) in bytes
                size_t l=eq.GetFullLength();

                if( l<eq.SizeOf() || (l&3) )
                {
                    AddError(DaqError(DaqError::EVENT_CORRUPTED1, DaqError::VALUE,l));
                    return;
                }

                if( l<(eq.SizeOf()+sizeof(SLink)) )
                {
                    pos += l;
                    continue;
                }

                const SLink *slink = reinterpret_cast<const SLink*>(eq.GetBuffer()+eq.SizeOf());
                if( slink->GetEventNumber()!=GetEventNumberInBurst() )
                    AddError(DaqError(DaqError::SLINK_WRONG_EVENT_NUMBER,DaqError::SOURCE_ID,slink->GetSourceID()));

                if( slink->IsSLinkMultiplexer() && opt.FixSLinkMultiplexerSize() )
                    const_cast<SLink*>(slink)->SetEventSize((l-eq.SizeOf())/4);

                if( l!=(slink->GetEventSize()*4+eq.SizeOf()) )
                {
                    AddError(DaqError(DaqError::EVENT_CORRUPTED2,DaqError::SOURCE_ID,slink->GetSourceID(),
                                   DaqError::VALUE,l,     DaqError::VALUE,slink->GetEventSize()*4 ));
                    pos += eq.GetFullLength();
                    continue;
                }

                // ----------------------------------------------------------
                // Create chip!
                // ----------------------------------------------------------

                const size_t chips_amount=chips.size();

                if( slink->IsSLinkMultiplexer() )
                {
                    // First pass to check that the size of the sub-SLinks does
                    // not exceed the size of the multiplexer SLink. This might
                    // happen if the multiplexer does not update the size of
                    // sub-SLink blocks that are truncated (as is the case with
                    // the DC05 Gandalf module in 2015). In this case discard
                    // the complete multiplexer data as there is no way we can
                    // detect which sub-SLink was truncated.
                    size_t pos=2;
                    bool badSizes=false;
                    while( !badSizes && pos<slink->GetEventSize() )
                    {
                        const SLink *sl = (const SLink*) ( ((uint32*)slink)+pos);

                        // Check that the event size of the sub-SLinks is not
                        // equal to zero, because we otherwise run into an
                        // endless loop. Attribute the error to the source ID
                        // of the multiplexer, but set the wrong size 0 as the
                        // error value.
                        if( 0==sl->GetEventSize() )
                        {
                            AddError(DaqError(DaqError::SLINKM_BAD_SIZE, DaqError::SOURCE_ID,slink->GetSourceID(),DaqError::VALUE,sl->GetEventSize()));
                            badSizes=true;
                        }

                        if( pos+sl->GetEventSize() > slink->GetEventSize() )
                        {
                            AddError(DaqError(DaqError::SLINKM_BAD_SIZE, DaqError::SOURCE_ID,slink->GetSourceID(),DaqError::VALUE,sl->GetEventSize()));
                            badSizes=true;
                        }

                        pos+=sl->GetEventSize();
                    }

                    if( !badSizes && pos!=slink->GetEventSize() )
                    {
                        AddError(DaqError(DaqError::SLINKM_BAD_SIZE, DaqError::SOURCE_ID,slink->GetSourceID(),DaqError::VALUE,slink->GetEventSize()));
                        badSizes=true;
                    }

                    pos=2;
                    while( !badSizes && pos<slink->GetEventSize() )
                    {
                        const SLink *sl = (const SLink*) ( ((uint32*)slink)+pos);

                        if( sl->GetEventNumber()!=GetEventNumberInBurst() )
                            AddError(DaqError(DaqError::SLINK_WRONG_EVENT_NUMBER,DaqError::SOURCE_ID,sl->GetSourceID()));

                        ChipCreate(sl,false,opt,call_me);

                        pos+=sl->GetEventSize();
                    }

                    if( !badSizes && chips_amount==chips.size() && slink->GetSourceID()!=opt.GetOnlineFilterSrcID() )
                        AddError(Exception("DaqEvent::ReadChips(): no chips are added!  srcID=%d",slink->GetSourceID()));
                }
                else {
                    ChipCreate(slink,false,opt,call_me);

                    if( chips_amount==chips.size() && slink->GetSourceID()!=opt.GetOnlineFilterSrcID() )
                        AddError(Exception("DaqEvent::ReadChips(): no chips are added!  srcID=%d",slink->GetSourceID()));
                }

                pos += eq.GetFullLength();
            }
        }
        else
        {
            // This is not a physical event. Skip it.
            AddError(DaqError(DaqError::EVENT_UNKNOWN,DaqError::VALUE,GetHeader().GetEventType()&Header::EVENT_TYPE_MASK));
            return;
        }
    }

    if( parent==NULL )
        CheckSrcIDs(opt);
}

////////////////////////////////////////////////////////////////////////////////

void DaqEvent::CheckSrcIDs(DaqOption &opt)
{
    // Check the list of sourceIDs (only for a top-level event!)
    if( parent!=NULL )
        return;

    if( !opt.GetEventSrcIDs().empty() )
    {
        set<uint16> intersection;
        set_intersection(opt.GetEventSrcIDs().begin(), opt.GetEventSrcIDs().end(),
                         srcIDs.begin(), srcIDs.end(),
                         inserter(intersection, intersection.begin()));

        if( intersection.size()!=opt.GetEventSrcIDs().size() || intersection.size()!=srcIDs.size() )
        {
            // There is a difference!
            set<uint16> srcIDs_missing, srcIDs_new;

            set_difference(opt.GetEventSrcIDs().begin(), opt.GetEventSrcIDs().end(),
                           intersection.begin(), intersection.end(),
                           inserter(srcIDs_missing, srcIDs_missing.begin()));

            set_difference(srcIDs.begin(), srcIDs.end(),
                           intersection.begin(), intersection.end(),
                           inserter(srcIDs_new, srcIDs_new.begin()));

            for( set<uint16>::const_iterator it=srcIDs_missing.begin(); it!=srcIDs_missing.end(); it++ )
                AddError(DaqError(DaqError::EVENT_MISSING_SRCID,DaqError::SOURCE_ID,*it));

            for( set<uint16>::const_iterator it=srcIDs_new.begin(); it!=srcIDs_new.end(); it++ )
                AddError(DaqError(DaqError::EVENT_UNKNOWN_SRCID,DaqError::SOURCE_ID,*it));
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void DaqEvent::ChipCreate(const SLink *sl,bool copy_buf,DaqOption &opt,DaqEventCallMe call_me)
{
  try
  {
    // Register the sourceID
    if( sl->GetSourceID()>=Chip::SourceID_MAX )
    {
        AddError(DaqError(DaqError::EVENT_SRCID_TOO_BIG,DaqError::SOURCE_ID,sl->GetSourceID()));
        return;
    }

    GetTopEvent().srcIDs.insert(sl->GetSourceID());

    // Do we need to scan it?
    if( !opt.NeedScanSrcID(sl->GetSourceID()) )
        return;

    if( sl->GetSourceID() == opt.GetOnlineFilterSrcID() )
    {
      OnlineFilter o(*sl);
      online_filter = o;
      return;
    }

    Chip *chip = Chip::Create(sl,copy_buf,opt,GetTopEvent());
    if( chip!=NULL )
    {
        chips.push_back(chip);
        if( call_me!=NULL )
        {
          try
          {
            call_me(*chips.back());
          }
          catch( ... )
          {
            AddError(Exception("DaqEvent::ReadChips(): There was a problem in call_me() for chip=%s src=%d",
                                  chips.back()->GetName().c_str(),chips.back()->GetSourceID()));
          }
        }
    }
  }
  catch( const DaqError &e )
  {
    cerr << "DaqEvent::ReadChips(): Internal problem\n";
    AddError(e);
  }
  catch( const Exception &e )
  {
    e.Print();
    AddError(e);
  }
  catch( const char *e )
  {
    cerr << e << "\n";
    AddError(Exception(e));
  }
  catch( ... )
  {
    AddError(Exception("DaqEvent::ReadChips(): Unknown exception"));
  }
}

////////////////////////////////////////////////////////////////////////////////

bool DaqEvent::IsGood(void) const
{
    return errors.Get(DaqError::SEVERE_PROBLEM).empty() && errors.Get(DaqError::Action::DISCARD_EVENT).empty();
}

////////////////////////////////////////////////////////////////////////////////

void DaqEvent::Clear(void)
{
  for( vector<Chip*>::iterator it=chips.begin(); it!=chips.end(); it++ )
    delete *it;
  chips.clear();
  errors.Remove();
  online_filter.Clear();
  srcID_geoID_bad.clear();
  srcIDs.clear();
}

////////////////////////////////////////////////////////////////////////////////

void DaqEvent::Print(ostream &o,const string &prefix) const
{
  o << prefix << "\n";
  o << prefix << "\n";

  o << prefix
    << "  version="     << GetHeader().GetVersion()
    << "  length="      << GetHeader().GetLength()
    << "  heade_len="   << GetHeader().GetHeaderLength()
    << "  time=("       << GetHeader().GetTime().first << "." << GetHeader().GetTime().second << ")"
    << "  burst="       << GetHeader().GetBurstNumber()
    << "  run="         << GetHeader().GetRunNumber()
    << "  event_run="   << GetHeader().GetEventNumberInRun()
    << "  event_burst=" << GetHeader().GetEventNumberInBurst()
    << "  event_type="  << GetHeader().GetEventType()
    << "  trigger="     << GetHeader().GetTrigger()
    << "  subevents="   << (GetHeader().HaveSubEvents()?"YES":"NO")
    << "\n";

  if( GetHeader().GetVersion()<0xffff )
  {
    o << prefix
      << "  dead_time="         << GetHeader().GetDeadTime().first << "."
                                << GetHeader().GetDeadTime().second << ")"
      << "  trigger_number="    << GetHeader().GetTriggerNumber()
      << "  error_code="        << GetHeader().GetErrorCode()
      << "\n";
  }

  if( chips.size()>0 )
    for( vector<Chip*>::const_iterator it=chips.begin(); it!=chips.end(); it++ )
      (*it)->Print(o,prefix+"chip: ");
  else
  {
    const uint32 *start=GetBuffer()+sizeof(Header)/4;
    for( unsigned i=0; i<(GetLength()-sizeof(Header))/4; i++ )
    {
      o << prefix;
      char s[111];
      sprintf(s,"%3d ",i);
      o<<s;
      bits_print(o,start[i],0,32,"%8x\n");
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

DaqEvent::~DaqEvent(void)
{
    Clear();
    if( flag_buffer_delete )
        free(buf);
}

////////////////////////////////////////////////////////////////////////////////

void DaqEvent::Check(void) const
{
    if( GetHeader().GetMagic()!=Header::EVENT_MAGIC_NUMBER )
        throw Exception("DaqEvent::Check(): Magic number was not found!");

    // Some checks...
    if( GetLength()%4!=0 )
        throw Exception("DaqEvent::Check():  bad event length %d (length%4!=0)",GetLength());

    // Some checks...
    if( GetLength()<GetHeader().GetHeaderLength() )
        throw Exception("DaqEvent::Check(): Bad event length.");
}

////////////////////////////////////////////////////////////////////////////////

DaqEvent::DaqEvent(void const * const event_buffer,bool copy_buf) :
    parent(NULL),
    flag_buffer_delete(false),
    allow_buf_change  (copy_buf),
    first_event_of_run(NULL),
    buf               ((uint8*)event_buffer),      // temporary regardless of copy_buf
    online_filter     (),
    trigger_mask      (0xffffffff)
{
    if( buf==NULL )
        throw Exception("DaqEvent::DaqEvent():  event_buffer=NULL");

    Check();

    if( copy_buf )
    {
        uint32 l = GetLength();
        buf = (uint8*)realloc(NULL,l);
        if( buf==NULL )
            throw Exception("DaqEvent::DaqEvent(): no memory");
        memcpy(buf,event_buffer,l);
        flag_buffer_delete = true;
    }

    Clear();
}

////////////////////////////////////////////////////////////////////////////////

namespace {
// Small class to read from streams.
struct Stream
{
    Stream(FILE *f)    : f1(f),   f2(NULL) {}
    Stream(istream &f) : f1(NULL),f2(&f)   {}

    void Check(void) const {if( (f1==NULL && f2==NULL) || (f1!=NULL && f2!=NULL) )
                              throw "Internal error: bad stream";}
    int Read(void *buf,int bytes);
    static int Read(FILE *f,void *buf,int bytes);
    static int Read(istream &f,void *buf,int bytes);

    FILE    *f1;
    istream *f2;
};

}

////////////////////////////////////////////////////////////////////////////////

DaqEvent::DaqEvent(istream &f) :
  parent(NULL),
  flag_buffer_delete(false),
  allow_buf_change(false),
  first_event_of_run(NULL),
  buf(NULL),
  online_filter(),
  trigger_mask      (0xffffffff)
{
    Stream ff(f);
    ReadFromStream(&ff);
}

////////////////////////////////////////////////////////////////////////////////

DaqEvent::DaqEvent(FILE *f) :
    parent(NULL),
    flag_buffer_delete(true),
    allow_buf_change(false),
    first_event_of_run(NULL),
    buf(NULL),
    online_filter(),
    trigger_mask      (0xffffffff)
{
    Stream ff(f);
    ReadFromStream(&ff);
}

////////////////////////////////////////////////////////////////////////////////

void DaqEvent::ReadFromStream(void *stream)
{
    Stream &f = *reinterpret_cast<Stream*>(stream);

    uint8 buffer[Header::ReadMinimum];
    buf=buffer;
    const Header event_header(buffer);
    int n=0;

    n = f.Read(buffer,Header::ReadMinimum);      // read event header

    if( n==0 )
        throw ExceptionEndOfStream();

    if( n != Header::ReadMinimum )
        throw Exception("DaqEvent::DaqEvent(): can not read event header");

    // Event header was read sucessfully.
    // Now we will read the entire event.

    Check();

    // Event length without header.
    const int event_length_without_header = GetLength() - Header::ReadMinimum;

    buf = (uint8*)malloc(GetLength());
    if( buf==NULL )
        throw Exception("DaqEvent::DaqEvent(): no memory");

    try
    {
        flag_buffer_delete=true;
        memcpy(buf, buffer, Header::ReadMinimum);

        n = f.Read(buf + Header::ReadMinimum, event_length_without_header ); // read rest of the event

        if( n!=int(event_length_without_header) )
            throw Exception("DaqEvent::DaqEvent(): Can not read whole event");

        Clear();
    }
    catch(...)
    {
        free(buf);
        throw;
    }
}

////////////////////////////////////////////////////////////////////////////////

int Stream::Read(void *buf,int bytes)
{
    Check();

    if( f1!=NULL )
        return Read(f1,buf,bytes);

    if( f2!=NULL )
        return Read(*f2,buf,bytes);

    throw "DaqEvent:: internal error in Stream::Read()";
}

////////////////////////////////////////////////////////////////////////////////

int Stream::Read(istream &f,void *buf,int bytes)
{
    f.read((char*)buf,bytes);
    return f.gcount();
}

////////////////////////////////////////////////////////////////////////////////

#if USE_RFIO
#include <shift.h>
#endif

int Stream::Read(FILE *f,void *buf,int bytes)
{
    return fread(buf,1,bytes,f);
}

////////////////////////////////////////////////////////////////////////////////

} // namespace CS
