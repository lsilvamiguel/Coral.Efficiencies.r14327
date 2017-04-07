#ifndef CsDDDStore_h
#define CsDDDStore_h

#include <string>
#include "CsStore.h"

namespace CS {class DaqEvent; class DaqEventsManager;}

class CsDDDStore: public CsStore
{
  public:

    static CsDDDStore*  Instance                (void);

    ~CsDDDStore() {}

    bool                init                    (void);

    bool                scan                    (void);

    bool                next                    (void);

    uint8 *             rawBuffer               (void);

    uint32              getEventInRun           (void) const;

    uint32              getRun                  (void) const;
  
    uint32              getEventInBurst         (void) const;

    uint32              getBurst                (void) const;

    uint32              getTriggerMask          (void) const;
  
    uint32              getErrorCode            (void) const;

    CsTime              getTime                 (void) const;

    int                 getRawBufferLength      (void);

  private:
  
    const CS::DaqEvent&   GetDaqEvent           (void) const;
    CS::DaqEventsManager& GetDaqEventsManager   (void) const;
  
    static CsDDDStore* instance_;
};

#endif // CsDDDStore_h
