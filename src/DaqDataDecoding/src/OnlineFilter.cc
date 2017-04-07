#include <cstring>

#include "OnlineFilter.h"
#include "SLink.h"
#include "Exception.h"

using namespace std;

namespace CS {

OnlineFilter::OnlineFilter()
  : fVersion(), fDecisions(), fStatistics(), fRawData()
{
  Clear();
}
	
OnlineFilter::OnlineFilter(const SLink &sl)
  : fVersion(), fDecisions(), fStatistics(), fRawData()
{
  const outblock_t   *outblock;
  const outblock_t   *end;

  Clear();

  if( sl.GetEventSize()<4 )
    throw Exception("OnlineFilter::OnlineFilter(): tag does not contain data");

  /* copy raw data */
  fRawData.clear();
  fRawData.reserve(sl.GetEventSize() - 3);
  for(uint32 i = 3; i < sl.GetEventSize(); ++i) {
    fRawData.push_back(((uint32 *)&sl)[i]);
  }

  outblock = reinterpret_cast<const outblock_t *>(&sl + 1);
  end = reinterpret_cast<const outblock_t *>((uint32 *)&sl + sl.GetEventSize());

  /* retrieve general information from the first outblock header */
  DecodeGeneral(outblock);

  /* loop over output blocks and decode the known ones */
  for(; outblock < end; outblock = reinterpret_cast<const outblock_t*>(reinterpret_cast<const char*>(outblock)+outblock->size)) {
    if((char *)outblock + outblock->size > (char *)end) {
      throw Exception("OnlineFilter::OnlineFilter(): corrupt filter tag");
    }
    switch(outblock->type) {
    case ot_decisions:
      DecodeDecisions(outblock);
      break;
    case ot_filter_version:
      DecodeVersion(outblock);
      break;
    case ot_missing_events:
      DecodeMissing(outblock);
      break;
    case ot_eob_filter:
      DecodeStatistics(outblock);
      break;
    default:
      break;
    }
  }

  fTagSeen = true;
}

OnlineFilter &
OnlineFilter::operator=(OnlineFilter &o)
{
  fTagSeen = o.fTagSeen;
  fStatsSeen = o.fStatsSeen;
  fIsMonitoring = o.fIsMonitoring;
  fIsAccepted = o.fIsAccepted;
  fMissing = o.fMissing;
  fDiscarded = o.fDiscarded;
  fTooBig = o.fTooBig;
  fNonPhysics = o.fNonPhysics;
  fNotWritten = o.fNotWritten;
  fRejected = o.fRejected;

  fDecisions = o.fDecisions;
  fStatistics = o.fStatistics;
  fRawData = o.fRawData;

  fStatAccept = o.fStatAccept;
  fStatReject = o.fStatReject;
  fStatDirect = o.fStatDirect;
  fStatTotal = o.fStatTotal;
  fStatTriggers = o.fStatTriggers;
  fStatFilters = o.fStatFilters;

  return *this;
}

void
OnlineFilter::Clear()
{
  fDecisions.clear();
  fRawData.clear();
  fStatistics.clear();

  fTagSeen = false;
  fStatsSeen = false;
  fIsMonitoring = false;
  fIsAccepted = false;
  fMissing = 0;
  fDiscarded = 0;
  fTooBig = 0;
  fNonPhysics = 0;
  fNotWritten = 0;
  fRejected = 0;

  fStatAccept = 0;
  fStatReject = 0;
  fStatDirect = 0;
  fStatTotal = 0;
  fStatTriggers = 0;
  fStatFilters = 0;
}

void
OnlineFilter::DecodeGeneral(const outblock_t *outblock)
{
  fIsMonitoring = outblock->monitoring;
  fIsAccepted = outblock->accepted;
}

void
OnlineFilter::DecodeDecisions(const outblock_t *outblock)
{
  const output_decision_t *decision;
  int                 decisions;

  decision = reinterpret_cast<const output_decision_t *>(outblock + 1);
  decisions = (outblock->size - sizeof(outblock_t)) / sizeof(output_decision_t);

  for(int i = 0; i < decisions; ++i) {
    filter_id_t id=decision[i].filter_id;
    fDecisions[id] =
      FilterDecision(true, decision[i].decided, decision[i].accepted,
		     decision[i].time * 8);
  }
}

void
OnlineFilter::DecodeVersion(const outblock_t *outblock)
{
  const char         *version;

  version = reinterpret_cast<const char *>(outblock + 1);
  if(version[outblock->size - 1 - sizeof(outblock_t)] != 0) {
    throw Exception("OnlineFilter::DecodeVersion(): string not NULL terminated");
  }

  fVersion = version;
}

void
OnlineFilter::DecodeMissing(const outblock_t *outblock)
{
  const output_missing_events_t *missing;

  if(outblock->size != sizeof(outblock_t) + sizeof(output_missing_events_t)) {
    throw Exception("OnlineFilter::DecodeMissing(): block size unknown");
  }

  missing = reinterpret_cast<const output_missing_events_t *>(outblock + 1);
  fMissing = missing->missing;
  fDiscarded = missing->discarded;
  fNotWritten = missing->not_written;
  fNonPhysics = missing->non_physics;
  fTooBig = missing->too_big;
  fRejected = missing->rejected;
}

void
OnlineFilter::DecodeStatistics(const outblock_t *outblock)
{
  const output_eob_filter_t *eob;
  unsigned int        size;

  eob = reinterpret_cast<const output_eob_filter_t *>(outblock + 1);
  size = 4 * eob->max_filter * (6 * (eob->max_trigger + 1) + 2) +
    sizeof(outblock_t) + sizeof(output_eob_filter_t);
  if(outblock->size != size) {
    throw Exception("OnlineFilter::DecodeStatistics: EOB statistics block corrupt");
  }

  fStatAccept = eob->accepted;
  fStatReject = eob->rejected;
  fStatDirect = eob->direct_writeout;
  fStatTotal = fStatAccept + fStatReject + fStatDirect;
  fStatFilters = eob->max_filter;
  fStatTriggers = eob->max_trigger;
  uint32 *data = (uint32 *)(eob + 1);

  for(unsigned int i = 0; i < fStatFilters; ++i) {
    fStatistics[static_cast<filter_id_t>(*data)] =
      FilterStatistics(data + 1, fStatTriggers);
    data += 6 * (fStatTriggers + 1) + 2;
  }

  fStatsSeen = true;
}

OnlineFilter::FilterDecision
OnlineFilter::GetFilterDecision(filter_id_t id) const
{
  map<filter_id_t, FilterDecision>::const_iterator it;

  it = fDecisions.find(id);
  if(it == fDecisions.end()) {
    return FilterDecision();
  } else {
    return it->second;
  }
}

OnlineFilter::FilterStatistics
OnlineFilter::GetFilterStatistics(filter_id_t id) const
{
  map<filter_id_t, FilterStatistics>::const_iterator it;

  it = fStatistics.find(id);
  if(it == fStatistics.end()) {
    return FilterStatistics();
  } else {
    return it->second;
  }
}

void
OnlineFilter::Print(ostream &out) const
{
  out << "OnlineFilter:" << hex << endl;
  out << "  fTagSeen:     " << fTagSeen << endl;
  out << "  fStatsSeen:   " << fStatsSeen << endl;
  out << "  fIsMonitorin: " << fIsMonitoring << endl;
  out << "  fIsAccepted:  " << fIsAccepted << endl;
  out << "  fMissing:     " << fMissing << endl;
  out << "  fDiscarded:   " << fDiscarded << endl;
  out << "  fTooBig:      " << fTooBig << endl;
  out << "  fNonPhysics:  " << fNonPhysics << endl;
  out << "  fNotWritten:  " << fNotWritten << endl;
  out << "  fRejected:    " << fRejected << endl;
  out << "  fVersion:     " << fVersion << '(' << fVersion.size() << ')' << endl;
  out << "  fStatAccept:  " << fStatAccept << endl;
  out << "  fStatReject:  " << fStatReject << endl;
  out << "  fStatDirect:  " << fStatDirect << endl;
  out << "  fStatTotal:   " << fStatTotal << endl;
  out << "  fStatFilters: " << fStatFilters << endl;
  out << "  fStatTrigger: " << fStatTriggers << endl;
}

OnlineFilter::FilterStatistics::FilterStatistics(uint32 *data, int triggers)
{
  if(triggers > MAX_TRIGGERS) {
    triggers = MAX_TRIGGERS;
  }

  memset(this, 0, sizeof(*this));

  not_called = *data++;

  accept_all = *data++;
  reject_all = *data++;
  no_idea_all = *data++;
  total_all = accept_all + reject_all + no_idea_all;

  time_accept_all = *data++;
  time_reject_all = *data++;
  time_no_idea_all = *data++;
  time_total_all = time_accept_all + time_reject_all + time_no_idea_all;

  int i;
  for(i = 0; i < triggers; i++) {
    accept[i] = *data++;
    reject[i] = *data++;
    no_idea[i] = *data++;
    total[i] = accept[i] + reject[i] + no_idea[i];

    time_accept[i] = *data++;
    time_reject[i] = *data++;
    time_no_idea[i] = *data++;
    time_total[i] = time_accept[i] + time_reject[i] + time_no_idea[i];
  }
}

} // namespace CS
