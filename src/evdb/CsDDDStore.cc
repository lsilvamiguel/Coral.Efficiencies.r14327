#include "CsDDDStore.h"
#include "CsInit.h"

#include "DaqDataDecoding/DaqEvent.h"
#include "DaqDataDecoding/DaqEventsManager.h"

using namespace std;

CsDDDStore* CsDDDStore::instance_ = NULL;

CsDDDStore* CsDDDStore::Instance()
{
    if( instance_ == NULL )
        instance_ = new CsDDDStore();

    return instance_;
}

bool CsDDDStore::init(void)
{
    CsInit::Instance();
    return true;
}

bool CsDDDStore::scan(void)
{
    CsInit *ci = CsInit::Instance();
    for( list <string*>::const_iterator it=ci->getDateFilesList().begin();
         it!=ci->getDateFilesList().end(); it++ )
        ci->getDaqEventsManager().AddDataSource(**it);

    return true;
}
  
bool CsDDDStore::next(void)
{
    return CsInit::Instance()->getDaqEventsManager().ReadEvent();
}

const CS::DaqEvent &CsDDDStore::GetDaqEvent(void) const
{
    return CsInit::Instance()->getDaqEventsManager().GetEvent();
}

uint8 * CsDDDStore::rawBuffer(void)
{
    return (uint8*)(GetDaqEvent().GetBuffer());
}

int CsDDDStore::getRawBufferLength(void)
{
    printf("CsDDDStore::getRawBufferLength(): no code!\n");
    return 0;
}

uint32 CsDDDStore::getEventInRun(void) const
{
    return GetDaqEvent().GetEventNumberInRun();
}

uint32 CsDDDStore::getRun(void) const
{
    return GetDaqEvent().GetRunNumber();
}
  
uint32 CsDDDStore::getEventInBurst(void) const
{
    return GetDaqEvent().GetEventNumberInBurst();
}

uint32 CsDDDStore::getBurst(void) const
{
    return GetDaqEvent().GetBurstNumber();
}

uint32 CsDDDStore::getTriggerMask(void) const
{
    return GetDaqEvent().GetTrigger();
}
  
uint32 CsDDDStore::getErrorCode(void) const
{
    return GetDaqEvent().GetErrorCode();
}

CsTime CsDDDStore::getTime(void) const
{
    pair<time_t,uint32> t = GetDaqEvent().GetTime();
    return CsTime(t.first,t.second);
}

