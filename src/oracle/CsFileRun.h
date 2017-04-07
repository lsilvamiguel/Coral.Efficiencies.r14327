#ifndef __CS_FILE_RUN_H__
#define __CS_FILE_RUN_H__
#if USE_FileStore

#include "CsDstObject.h"

class CsFileRun : public CsDstRun
{
  std::vector<std::string> fileNames;
  FlatFile                 dstFile;
  unsigned                 curChunkIdx;
  CsDstEvent*              dstEvent;
public:
  CsFileRun(int runnum);
  virtual ~CsFileRun();
  
  virtual bool        Init(std::vector<std::string>& file_names, std::string& chunkName);
  virtual bool        SelectNextChunk(int slot_number, const std::string& chunkName);
  virtual bool        NextEvent();
  virtual uint8*      GetRawEvent() { return NULL; }
  virtual CsDstEvent* GetDstEvent() { return dstEvent; }

  virtual uint32      GetRawEventLength() const { return 0; }
  virtual uint32      GetEventInRun()     const { if(dstEvent) return dstEvent->GetEventNumber();  return 0; }
  virtual uint32      GetEventInBurst()   const { if(dstEvent) return dstEvent->GetEventInBurst(); return 0; }
  virtual uint32      GetBurstNumber()    const { if(dstEvent) return dstEvent->GetBurstNumber();  return 0; }
  virtual uint32      GetTriggerMask()    const { if(dstEvent) return dstEvent->GetTriggerMask();  return 0; }
  virtual uint32      GetErrorCode()      const { if(dstEvent) return dstEvent->GetErrorCode();    return 0; }
  virtual CsTime      GetTime()           const { if(dstEvent) return CsTime(dstEvent->GetTimeInSec(),dstEvent->GetTimeInUSec()); return CsTime(0,0); }
  virtual bool        IsDstExist()        const { return (dstEvent != NULL); }

  virtual bool        UpdateRunLog(const std::string& tbnames,int slot) { return false; }
  virtual bool        CloseDst() { return true; }
  virtual void        Abort() {};
protected:
  virtual bool        OpenFile(std::string& path);
  virtual void        Print(std::ostream& msgout) const {};
};

#endif //#if USE_FileStore
#endif
