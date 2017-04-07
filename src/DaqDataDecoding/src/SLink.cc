#include <cstdio>
#include "SLink.h"
#include "Chip.h"
#include "DaqOption.h"

using namespace std;

namespace CS {

////////////////////////////////////////////////////////////////////////////////

SLinkMultiplexer::SLinkMultiplexer(void) :
  event_size     (0),
  source_id      (0),
  event_type     (0),
  error_flag     (0),
  event_number   (0),
  spill_number   (0),
  stat           (0)
{}

////////////////////////////////////////////////////////////////////////////////

bool SLinkMultiplexer::IsSLinkMultiplexer(void) const
{
  return GetSourceID()>=896 && GetSourceID()<=1023;
}

////////////////////////////////////////////////////////////////////////////////

SLink::SLink(void) :
  TCS_error      (8),
  errors_counter (0),
  format         (0)
{}

////////////////////////////////////////////////////////////////////////////////

void SLink::AnalyseErrors(Chip &chip,DaqOption &opt) const
{
  if( GetTCSError()& 1   )
    chip.AddError(DaqError(DaqError::TCS_FIFO,DaqError::SOURCE_ID,GetSourceID()),opt);

  if( !(GetTCSError()& 8)   )
    chip.AddError(DaqError(DaqError::TCS_ECC,DaqError::SOURCE_ID,GetSourceID()),opt);

  if( GetTCSError()&64   )
    chip.AddError(DaqError(DaqError::TCS_SYNC,DaqError::SOURCE_ID,GetSourceID()),opt);

  if( GetTCSError()&0xb6 )
    chip.AddError(DaqError(DaqError::TCS_UNDEF,DaqError::SOURCE_ID,GetSourceID()),opt);

  if( IsError() )
    chip.AddError(DaqError(DaqError::SLINK_ERRFLAG,DaqError::SOURCE_ID,GetSourceID()),opt);

  if( GetErrorsCounter() )
    chip.AddError(DaqError(DaqError::SLINK_ERRNR,
                  DaqError::SOURCE_ID,GetSourceID(),
                  DaqError::COUNTER,GetErrorsCounter()),opt);
}

////////////////////////////////////////////////////////////////////////////////

void SLinkMultiplexer::Print(ostream &o,const string &prefix) const
{
  o << prefix;
  char s[222];
  sprintf(s,"EventSize=%u sourceID=%u EventType=%u Error=%u Event=%u Spill=%u stat=%u\n",
             event_size,source_id,event_type,error_flag,event_number,spill_number,stat);
  o << s;
}

////////////////////////////////////////////////////////////////////////////////

void SLink::Print(ostream &o,const string &prefix) const
{
  SLinkMultiplexer::Print(o,prefix);
  o << prefix;
  char s[222];
  sprintf(s,"  status=%u TCS_error=%u errors_counter=%u format=%u\n",
             status,TCS_error,errors_counter,format);
  o << s;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace CS
