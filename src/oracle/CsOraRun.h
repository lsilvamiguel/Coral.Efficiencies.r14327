#ifndef __CS_ORA_RUN_H__
#define __CS_ORA_RUN_H__

#include <string>
#include <vector>
#include <map>

#include <CsTime.h>

#include "CsOraSession.h"
#include "CsDstObject.h"
#include "CsDst_04.h"

struct OraDstEventInfo
{
  uint32               eventNumber;
  uint32               eventSize;
  uint32               eventOffset;
  uint32               triggerMask;
  std::vector<float32> extVar;
};

struct OraDstChunkInfo
{
  uint32                       fileId;
  std::string                  dirName;
  std::string                  fileName;
  uint32                       fileSize;
  uint32                       dstVersion;
  std::vector<std::string>     extVarName;
  std::vector<OraDstEventInfo> eventInfo;
};

struct OraRawEventInfo
{
  uint32 eventNumber;
  uint32 burstNumber;
  uint32 eventInBurst;
  uint32 triggerMask;
  uint32 timeInSec;
  uint32 timeInUsec;
  uint32 errorCode;
  uint32 eventSize;
  uint32 fileOffset;
  bool   eventToRead;
};

struct OraRawChunkInfo
{
  uint32                       fileId;
  std::string                  dirName;
  std::string                  fileName;
  uint32                       fileSize;
  std::vector<OraRawEventInfo> eventInfo;
  OraDstChunkInfo*             dstChunkInfo[3]; // 3 dst slots are available
};

struct OraRunInfo
{
  uint32                       sorFileId;
  std::vector<OraRawChunkInfo> chunkInfo;
};

struct OraEventInfo
{
  OraRawEventInfo* rawEventInfo;
  OraDstEventInfo* dstEventInfo;
};

const int IS_RAW_READING       = 1;
const int IS_DST_READING       = 2;
const int IS_DST_PRODUCTION_04 = 3; // DST version 4 production
const int IS_RUN_LOG_UPDATE    = 4;

class CsOraRun : public CsDstRun
{
  CsOraSession*                           session;

  bool                                    isDstProduction;
  bool                                    isDstReading;
  bool                                    isRawReading;
  bool                                    isRunLogUpdate;

  OraRawChunkInfo*                        curChunkInfo;
  OraRunInfo*                             runInfo;
  std::map<uint32,OraEventInfo>           eventMap;
  std::map<uint32,OraEventInfo>::iterator eventIt;
  bool                                    isNewChunk;

  std::vector<uint8>                      rawBuffer;
  FlatFile                                rawFile;
  FlatFile                                dstFile;

  int                                     slotNumber;
  std::string                             runPeriod;
  int                                     runYear;

 public:
  CsOraRun(int runnum);
  virtual ~CsOraRun();
  
  virtual bool        Init(int todo);
  virtual bool        SelectNextChunk(int slotNumber, const std::string& chunkName = ""); // slot number = [0..2]
  virtual bool        NextEvent();
  virtual uint8*      GetRawEvent();
  virtual CsDstEvent* GetDstEvent();
  virtual bool        CloseDst();
  virtual void        Abort();
  virtual bool        SaveDstEvent(const CsDstEvent*);
  virtual bool        UpdateRunLog(const std::string& runLog, int slot);
  virtual void        Print(std::ostream&) const;

  virtual uint32 GetRawEventLength() const { return eventIt->second.rawEventInfo->eventSize; }
  virtual uint32 GetEventInRun()     const { return eventIt->second.rawEventInfo->eventNumber; }
  virtual uint32 GetEventInBurst()   const { return eventIt->second.rawEventInfo->eventInBurst; }
  virtual uint32 GetBurstNumber()    const { return eventIt->second.rawEventInfo->burstNumber; }
  virtual uint32 GetTriggerMask()    const { return eventIt->second.rawEventInfo->triggerMask; }
  virtual uint32 GetErrorCode()      const { return eventIt->second.rawEventInfo->errorCode; }
  virtual CsTime GetTime()           const 
    { return CsTime(eventIt->second.rawEventInfo->timeInSec,eventIt->second.rawEventInfo->timeInUsec); }
  virtual bool   IsDstExist()        const { return (eventIt->second.dstEventInfo != NULL); }

  // Raw file information subroutines:
  virtual std::string GetRawFileName()    const { return curChunkInfo?curChunkInfo->fileName:""; }
  virtual std::string GetRawFileDirName() const { return curChunkInfo?curChunkInfo->dirName:""; }
  virtual uint32      GetRawFileSize()    const { return curChunkInfo?curChunkInfo->fileSize:0; }
};

#endif
