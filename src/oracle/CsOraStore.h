#ifndef _CS_ORA_STORE_HH_
#define _CS_ORA_STORE_HH_

#include "FlatFile.h"
#include "CsTime.h"
#include "CsStore.h"
#include "CsEndOfJob.h"
#include "CsRecoEvent.h"
#include "CsOraRun.h"
#include "CsDstObject.h"

#include <vector>
#include <math.h>

const int LAST_DST_VERSION = 4;

typedef void (*ExtVarFunc)();

class CsOpt;

class CsOraStore : virtual public CsStore, public CsEndOfJob
{
 public:
  enum  DataType
  {
    DATA_TYPE_ERROR = 0,
    RAW_DATA_TYPE   = 1,
    DST_DATA_TYPE   = 2
  };

  static CsOraStore *Instance();
  
  //------- reader part
  
  virtual bool 		init();        
  virtual bool 		scan();         // Prepare the iteration
  virtual bool 		next();         // Next event

  //=========  Event related methods
  
  virtual uint8* 	rawBuffer()             { return oraRun->GetRawEvent(); }
  virtual int 		getRawBufferLength()    { return oraRun->GetRawEventLength(); }
  virtual uint32 	getEventInRun()   const { return oraRun->GetEventInRun(); }
  virtual uint32 	getRun()          const { return oraRun->GetRunNumber(); }
  virtual uint32 	getEventInBurst() const { return oraRun->GetEventInBurst(); }
  virtual uint32 	getBurst()        const { return oraRun->GetBurstNumber(); }
  virtual uint32 	getTriggerMask()  const { return oraRun->GetTriggerMask(); }
  virtual uint32 	getErrorCode()    const { return oraRun->GetErrorCode(); }
  virtual CsTime 	getTime()         const { return oraRun->GetTime(); }
  virtual const char*   getTBNames()      const { return _allTBNames.c_str(); }
  virtual bool          dst(int slot=0)         { return oraRun->IsDstExist(); }
  //==========
  
  //! DST readout
  bool downloadDST(CsRecoEvent*, int version = 0);

  //! DST production (requires consistent use of RDB tools)
  bool uploadDST(CsRecoEvent*, int fatness, int version);   

  bool downloadTBNamesString(std::string &log);
  bool uploadTBNamesString(std::string &log);

  virtual void          saveAndExit(void);  //!< Transaction commit (i.e. save changes)
  virtual void          abortAndExit(void); //!< Transaction abort (i.e. discard changes)
  virtual bool          end(void);          //!< Commits transaction and closes DST file.

  int  getSelectedSlot() { return dstSlot; }

  static void        setExtVarName(int n, const std::string& name);
  static const char* getExtVarName(int n) { return _extVarNames[n].c_str(); }
  static void        setExtVariable(int n, float32 val)
    { ((n<3 && n>=0)?_extVars[n] = val:0); }
  static float32     getExtVariable(unsigned number)
    { return ((number<3)?_extVars[number]:0); }

  void setPreuploadFunc(ExtVarFunc func) { _preuploadFunc = func; }
  std::string GetCurrentChunkName() const { return (oraRun?oraRun->GetRawFileName():""); }
  std::string GetCurrentRawDir()    const { return (oraRun?oraRun->GetRawFileDirName():""); }

  std::string RawSelectionCriteria(){return rawSelectionCriteria;}

 protected:

            CsOraStore();
   virtual ~CsOraStore();

   bool GetRawSelectionCriteria(CsOpt* myOptions);
   bool GetDstSelectionCriteria(CsOpt* myOptions);
   bool ParseDstCriteriaString();

 protected:
   DataType        dataType;
   std::string     chunkName;
   int             runNumber;
   int             dstSlot;
   int             dstVersion;
   bool            isDstProduction;

   CsOraRun*       oraRun;

   bool            isInit;
   bool            isScanned;

   std::string     _allTBNames;
   std::string     dstSelectionCriteria;
   std::string     rawSelectionCriteria;

   ExtVarFunc                       _preuploadFunc;     
   static std::vector<std::string>  _extVarNames; // Names of the external variables
   static std::vector<float32>      _extVars;     // Values of the external variables
};

#endif
