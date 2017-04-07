#ifndef __CS_FILE_STORE_H__
#define __CS_FILE_STORE_H__
#if USE_FileStore

#include "CsStore.h"
#include "CsFileRun.h"

class CsFileStore : virtual public CsStore, public CsEndOfJob
{
protected:
  std::vector<std::string>  file_names;
public:
  static CsFileStore* Instance();

  virtual bool 		init();        
  virtual bool 		scan();         // Prepare the iteration
  virtual bool 		next();         // Next event

  //=========  Event related methods
  
  virtual uint8* 	rawBuffer()             { return fileRun->GetRawEvent(); }
  virtual int 		getRawBufferLength()    { return fileRun->GetRawEventLength(); }
  virtual uint32 	getEventInRun()   const { return fileRun->GetEventInRun(); }
  virtual uint32 	getRun()          const { return fileRun->GetRunNumber(); }
  virtual uint32 	getEventInBurst() const { return fileRun->GetEventInBurst(); }
  virtual uint32 	getBurst()        const { return fileRun->GetBurstNumber(); }
  virtual uint32 	getTriggerMask()  const { return fileRun->GetTriggerMask(); }
  virtual uint32 	getErrorCode()    const { return fileRun->GetErrorCode(); }
  virtual CsTime 	getTime()         const { return fileRun->GetTime(); }
  virtual const char*   getTBNames()      const { return _allTBNames.c_str(); }
  virtual bool          dst(int slot=0)         { return fileRun->IsDstExist(); }
  //==========
  //! DST readout
  virtual bool downloadDST(CsRecoEvent*, int version = 0);

  virtual bool downloadTBNamesString(std::string &log);
  virtual int  getSelectedSlot() { return dstSlot; }
  virtual bool end(void) { return true; }          //!< Commits transaction and closes DST file.

 protected:
  CsFileStore();
  virtual ~CsFileStore();

  virtual bool ScanDirectory(std::string& dir, int runNum, std::string& chunk);

protected:
   std::string     chunkName;
   int             runNumber;
   int             dstSlot;
   int             dstVersion;

   CsFileRun*      fileRun;
  
   bool            isInit;
   bool            isScanned;

   std::string     _allTBNames;

};
#endif //#if USE_FileStore
#endif

