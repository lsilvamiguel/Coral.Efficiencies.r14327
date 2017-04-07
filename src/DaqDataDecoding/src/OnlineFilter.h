#ifndef CompassSoft_OnlineFilter__include
#define CompassSoft_OnlineFilter__include

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include "config.h"

namespace CS {

class SLink;

class OnlineFilter
{
  public:

#include "OnlineFilterTag.h"

  class FilterDecision {

  public:
    FilterDecision() : fHasRun(false), fHasDecided(false) {}
    FilterDecision(bool run, bool dec, bool acc, int time)
      : fHasRun(run), fHasDecided(dec), fHasAccepted(acc), fTimeTaken(time) {}

    bool                HasRun() const {return fHasRun;}
    bool                HasDecided() const {return fHasDecided;}
    bool                HasAccepted() const {return fHasAccepted;}
    int                 TimeTaken() const {return fTimeTaken;}

  private:

    bool                fHasRun;
    bool                fHasDecided;
    bool                fHasAccepted;
    int                 fTimeTaken;

  };

  struct FilterStatistics {

    static const int    MAX_TRIGGERS = 24;

    FilterStatistics() {memset(this, 0, sizeof(*this));}
    FilterStatistics(uint32 *, int);

    /* these methods take a trigger number which defaults to -1, which works
       because the _all variables are declared right before the corresponding
       trigger number arrays. */
    double              GetAcceptRate(int i = -1) const;
    double              GetRejectRate(int i = -1) const;
    double              GetNoIdeaRate(int i = -1) const;
    double              GetTimeAccept(int i = -1) const;
    double              GetTimeReject(int i = -1) const;
    double              GetTimeNoIdea(int i = -1) const;
    double              GetTimeTotal(int i = -1) const;
    double              GetCalledRate() const;

    uint32              not_called;

    uint32              accept_all;
    uint32              accept[MAX_TRIGGERS];
    uint32              reject_all;
    uint32              reject[MAX_TRIGGERS];
    uint32              no_idea_all;
    uint32              no_idea[MAX_TRIGGERS];
    uint32              total_all;
    uint32              total[MAX_TRIGGERS];

    uint32              time_accept_all;
    uint32              time_accept[MAX_TRIGGERS];
    uint32              time_reject_all;
    uint32              time_reject[MAX_TRIGGERS];
    uint32              time_no_idea_all;
    uint32              time_no_idea[MAX_TRIGGERS];
    uint32              time_total_all;
    uint32              time_total[MAX_TRIGGERS];
  };

  public:
                        OnlineFilter            (const SLink &sl);
                        OnlineFilter            (void);
    OnlineFilter       &operator=               (OnlineFilter &);

    void                Clear                   (void);

  private:

                        OnlineFilter            (OnlineFilter &) {}

  public:
  
    bool                FilterTagPresent        (void) const {return fTagSeen;}
    bool                FilterStatisticsPresent (void) const {return fStatsSeen;}
    bool                IsMonitoring            (void) const {return fIsMonitoring;}
    bool                IsAccepted              (void) const {return fIsAccepted;}
    int                 EventsMissing           (void) const {return fMissing;}
    int                 EventsDiscarded         (void) const {return fDiscarded;}
    int                 EventsTooBig            (void) const {return fTooBig;}
    int                 EventsNonPhysics        (void) const {return fNonPhysics;}
    int                 EventsNotWritten        (void) const {return fNotWritten;}
    int                 EventsRejected          (void) const {return fRejected;}
    FilterDecision      GetFilterDecision       (filter_id_t id) const;

    FilterStatistics    GetFilterStatistics     (filter_id_t id) const;
    double              GetAcceptRate           () const;
    double              GetRejectRate           () const;
    double              GetDirectRate           () const;
    uint32              GetStatEvents           () const {return fStatTotal;}

    const std::vector<uint32> &GetRawData       (void) const {return fRawData;}

    void                Print                   (std::ostream &o=std::cout) const;

  private:

    void                DecodeGeneral           (const outblock_t *);
    void                DecodeDecisions         (const outblock_t *);
    void                DecodeVersion           (const outblock_t *);
    void                DecodeMissing           (const outblock_t *);
    void                DecodeStatistics        (const outblock_t *);

  private:

    bool                fTagSeen;
    bool                fStatsSeen;
    bool                fIsMonitoring;
    bool                fIsAccepted;
    int                 fMissing;
    int                 fDiscarded;
    int                 fTooBig;
    int                 fNonPhysics;
    int                 fNotWritten;
    int                 fRejected;
    std::string         fVersion;
    std::map<filter_id_t, FilterDecision> fDecisions;
    std::map<filter_id_t, FilterStatistics> fStatistics;
    uint32              fStatAccept;
    uint32              fStatReject;
    uint32              fStatDirect;
    uint32              fStatTotal;
    uint32              fStatFilters;
    uint32              fStatTriggers;
    std::vector<uint32> fRawData;
};
 
} // namespace CS

inline double CS::OnlineFilter::FilterStatistics::GetAcceptRate(int trigger) const
{
  if(total[trigger]) {
    return (double)accept[trigger] / total[trigger];
  } else {
    return 0;
  }
}

inline double CS::OnlineFilter::FilterStatistics::GetRejectRate(int trigger) const
{
  if(total[trigger]) {
    return (double)reject[trigger] / total[trigger];
  } else {
    return 0;
  }
}

inline double CS::OnlineFilter::FilterStatistics::GetNoIdeaRate(int trigger) const
{
  if(total[trigger]) {
    return (double)no_idea[trigger] / total[trigger];
  } else {
    return 0;
  }
}

inline double CS::OnlineFilter::FilterStatistics::GetTimeAccept(int trigger) const
{
  if(accept[trigger]) {
    return (double)time_accept[trigger] / accept[trigger];
  } else {
    return 0;
  }
}

inline double CS::OnlineFilter::FilterStatistics::GetTimeReject(int trigger) const
{
  if(reject[trigger]) {
    return (double)time_reject[trigger] / reject[trigger];
  } else {
    return 0;
  }
}

inline double CS::OnlineFilter::FilterStatistics::GetTimeNoIdea(int trigger) const
{
  if(no_idea[trigger]) {
    return (double)time_no_idea[trigger] / no_idea[trigger];
  } else {
    return 0;
  }
}

inline double CS::OnlineFilter::FilterStatistics::GetTimeTotal(int trigger) const
{
  if(total[trigger]) {
    return (double)time_total[trigger] / total[trigger];
  } else {
    return 0;
  }
}

inline double CS::OnlineFilter::FilterStatistics::GetCalledRate() const
{
  if(total_all + not_called) {
    return (double)total_all / (total_all + not_called);
  } else {
    return 0;
  }
}

inline double CS::OnlineFilter::GetAcceptRate() const
{
  if(fStatTotal) {
    return (double)fStatAccept / (double)fStatTotal;
  } else {
    return 0;
  }
}

inline double CS::OnlineFilter::GetRejectRate() const
{
  if(fStatTotal) {
    return (double)fStatReject / (double)fStatTotal;
  } else {
    return 0;
  }
}

inline double CS::OnlineFilter::GetDirectRate() const
{
  if(fStatTotal) {
    return (double)fStatDirect / (double)fStatTotal;
  } else {
    return 0;
  }
}

#endif // CompassSoft_OnlineFilter__include
