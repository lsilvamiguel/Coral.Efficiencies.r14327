#ifndef __CS_DST_OBJECT_H__
#define __CS_DST_OBJECT_H__

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <CLHEP/Matrix/Matrix.h>

#include <CsEvent.h>
#include <CsTime.h>
#include "CsStoreMisc.h"
#include "FlatFile.h"

const uint32 ERROR_TYPE           =     0;
const uint32 VECTOR3_TYPE         =    10;
const uint32 MATRIX3_TYPE         =    20;
const uint32 CLUSTER_TYPE_04      =    34;
const uint32 HELIX_TYPE_04        =    44;
const uint32 TRACK_TYPE_04        =    54;
const uint32 VERTEX_TYPE_04       =    64;
const uint32 CALOBJ_TYPE_04       =    74;
const uint32 PARTICLE_TYPE_04     =    84;
const uint32 EVENT_HEADER_TYPE_04 = 10400; // for compatibility with previous version
const uint32 DST_EVENT_TYPE_04    = 10401;
const uint32 DST_CHUNK_TYPE_04    = 10402;

class CsDstObject
{
 public:
  uint32  type;             // unique DST object's type
  int     dstVersion;       // DST version for the object

 public:
  CsDstObject() { type = ERROR_TYPE; dstVersion = 0; }
  virtual ~CsDstObject() { };

  //!                      bool  Load(std::istream&) 
  //  downloads the DST object from the input stream.
  //  returns true if everything is OK, otherwise returns false
  virtual bool Load(std::istream&) = 0;
  virtual bool Load(FlatFile&) = 0;

  //!                 bool  IsEquil(const CssObject*)
  // returns true if both DST objects are equivalented, otherwise false
  virtual bool IsEquil(const CsDstObject*) const = 0;

  //!                      bool  Save(ofstream&)
  // Saves DST object to the output stream.
  // Returns true if no I/O error during writing of the data
  virtual bool Save(std::ofstream&) const;
  virtual bool Save(FlatFile&) const;

  //!                      void Print(std::ostream&)
  // Prints dump of the DST object to output stream
  virtual void Print(std::ostream&) const = 0;

  //!             void GetBuffer(std::vector<uint8>& buffer)
  // Puts DST object's content in "buffer"
  virtual void GetBuffer(std::vector<uint8>& buffer) const = 0;

  //!                    int GetDSTVersion()
  // returns DST version of the object
  int GetDSTVersion() const  { return dstVersion; }
};

class CsDstCluster : public CsDstObject
{
 public:
  CsDstCluster()          { };
  virtual ~CsDstCluster() { };

  virtual void        Init(const CsCluster& cl) = 0;
  virtual CsCluster*  GetCluster() const = 0;
};

class CsDstHelix : public CsDstObject
{
 public:
  CsDstHelix()  {};
  virtual  ~CsDstHelix() {};

  virtual void     Init(const CsHelix& hl) = 0;
  virtual CsHelix* GetHelix() const = 0;
};

class CsDstTrack : public CsDstObject
{
 public:
  CsDstTrack() {};
  virtual ~CsDstTrack() {};

  virtual void     Init(const CsTrack& tr) = 0;
  virtual CsTrack* GetTrack() const = 0;
  virtual int      GetNClusters() const = 0;    // returns number of cluster for the track (if fatness=10)
  virtual void     SetAsBeam() = 0;             // sets track as beam
};

class CsDstVertex : public CsDstObject
{
 public:
  CsDstVertex()  {};
  virtual ~CsDstVertex() {};

  virtual void      Init(const CsVertex& vt) = 0;
  virtual CsVertex* GetVertex() const = 0;
};

class CsDstCalobject : public CsDstObject
{
 public:
  CsDstCalobject() {};
  virtual ~CsDstCalobject() {};

  virtual void  Init(const Reco::CalorimeterParticle& co) = 0;
  virtual Reco::CalorimeterParticle* GetCalobject() const = 0;
};

class CsDstParticle : public CsDstObject
{
 public:
  CsDstParticle() {};
  virtual ~CsDstParticle() {};

  virtual void        Init(const CsParticle& pt) = 0;
  virtual CsParticle* GetParticle() const = 0;
  virtual int16       GetPartType() const = 0;
  virtual std::string GetPartName() const = 0;
};

class CsDst3Vector : public CsDstObject
{
  float64 vec[3];
 public:
  CsDst3Vector();
  virtual ~CsDst3Vector() {};

  virtual void           Init(const float64* v);
  virtual void           Init(const Cs3Vector& v);
  virtual const float64* GetVector() const { return &vec[0]; }
  virtual void           GetBuffer(std::vector<uint8>& buffer) const;

  virtual bool           Load(std::istream& file);
  virtual bool           Load(FlatFile& file);
  virtual bool           IsEquil(const CsDstObject* obj) const;
  virtual void           Print(std::ostream& msgout) const;
};

class CsDst3Matrix : public CsDstObject
{
  float32 cov[9];
 public:
  CsDst3Matrix();
  virtual ~CsDst3Matrix() {};

  virtual void               Init(const float32* m);
  virtual void               Init(const CLHEP::HepMatrix& m);
  virtual CLHEP::HepMatrix*  GetHepMatrix() const;
  virtual void               GetBuffer(std::vector<uint8>& buffer) const;

  virtual bool        Load(std::istream& file);
  virtual bool        Load(FlatFile& file);
  virtual bool        IsEquil(const CsDstObject* obj) const;
  virtual void        Print(std::ostream& msgout) const;
};

class CsDstEventHeader : public CsDstObject
{
 protected:
  uint32  runNumber;
  uint32  eventNumber;
  uint32  eventInBurst;
  uint32  burstNumber;
  uint32  timeInSec;
  uint32  timeInUSec;
  uint32  errorCode;
  uint32  triggerMask;
 public:
  CsDstEventHeader();
  virtual ~CsDstEventHeader() {};

  virtual void    Init(const CsRecoEvent& ev) = 0;
  virtual void    GetHeader(CsRecoEvent& ev) const = 0;

  virtual float32 GetExternalVariable(int n) const { return 0.; }
  virtual void    SetExternalVariable(int n, float32 val) { };

  uint32  GetRunNumber()             const { return runNumber; }
  uint32  GetEventNumber()           const { return eventNumber; }
  uint32  GetEventInBurst()          const { return eventInBurst; }
  uint32  GetBurstNumber()           const { return burstNumber; }
  uint32  GetTimeInSec()             const { return timeInSec; }
  uint32  GetTimeInUSec()            const { return timeInUSec; }
  uint32  GetErrorCode()             const { return errorCode; }
  uint32  GetTriggerMask()           const { return triggerMask; }


  void    SetRunNumber(uint32 val)         { runNumber = val; }
  void    SetEventNumber(uint32 val)       { eventNumber = val; }
  void    SetEventInBurst(uint32 val)      { eventInBurst = val; }
  void    SetBurstNumber(uint32 val)       { burstNumber = val; }
  void    SetTimeInSec(uint32 val)         { timeInSec = val; }
  void    SetTimeInUSec(uint32 val)        { timeInUSec = val; }
  void    SetErrorCode(uint32 val)         { errorCode = val; }
  void    SetTriggerMask(uint32 val)       { triggerMask = val; }
};

class CsDstEvent : public CsDstObject
{
 protected:
  uint32  runNumber;
  uint32  eventNumber;
  uint32  eventInBurst;
  uint32  burstNumber;
  uint32  timeInSec;
  uint32  timeInUSec;
  uint32  errorCode;
  uint32  triggerMask;

  std::vector<float32> extVars;

  bool    isExist;
 public:
  CsDstEvent();
  virtual ~CsDstEvent() { ; }

  virtual void    Init(const CsRecoEvent& ) = 0;
  virtual void    GetEvent(CsRecoEvent&) const = 0;

  virtual uint32  GetRunNumber()             const { return runNumber; }
  virtual uint32  GetEventNumber()           const { return eventNumber; }
  virtual uint32  GetEventInBurst()          const { return eventInBurst; }
  virtual uint32  GetBurstNumber()           const { return burstNumber; }
  virtual uint32  GetTimeInSec()             const { return timeInSec; }
  virtual uint32  GetTimeInUSec()            const { return timeInUSec; }
  virtual uint32  GetErrorCode()             const { return errorCode; }
  virtual uint32  GetTriggerMask()           const { return triggerMask; }


  virtual void    SetRunNumber(uint32 val)         { runNumber = val; }
  virtual void    SetEventNumber(uint32 val)       { eventNumber = val; }
  virtual void    SetEventInBurst(uint32 val)      { eventInBurst = val; }
  virtual void    SetBurstNumber(uint32 val)       { burstNumber = val; }
  virtual void    SetTimeInSec(uint32 val)         { timeInSec = val; }
  virtual void    SetTimeInUSec(uint32 val)        { timeInUSec = val; }
  virtual void    SetErrorCode(uint32 val)         { errorCode = val; }
  virtual void    SetTriggerMask(uint32 val)       { triggerMask = val; }

  virtual float32 GetExternalVariable(int n) const { return 0.; }
  virtual void    SetExternalVariable(int n, float32 val) { };
};

class CsDstRun;

class CsDstChunk : public CsDstObject
{
 protected:
  std::string              chunkName;
  std::string              TBNames;
  std::vector<std::string> extVarNames;
 public:
  CsDstChunk() { };
  virtual ~CsDstChunk() { };

  virtual void        Init(const char* chunkName ) = 0;
  virtual CsDstEvent* GetNextEvent(std::istream& file) const = 0;
  virtual CsDstEvent* GetNextEvent(FlatFile& file) const = 0;
  virtual CsDstEvent* GetNextEvent(std::vector<uint8>& buffer) const = 0;
  virtual const char* GetChunkName() const { return chunkName.c_str(); }
  virtual void        SetExtVarName(int n, const char* name) = 0;
  virtual const char* GetExtVarName(int n) const = 0;
  virtual const char* GetTBNames() const { return TBNames.c_str(); }
  virtual void        SetTBNames(const char* tbNames) { TBNames = tbNames; }
};

class CsDstRun : public CsDstObject
{
 protected:
  int           runNumber;
  time_t        start;
  time_t        stop;
  long          status; // for compatibility with previous version
  std::string   runLog;
  CsDstChunk*   dstChunk;
 public:
  CsDstRun() : runNumber(0), start(0), stop(0), status(0) { dstChunk = NULL; }
  virtual ~CsDstRun() { delete dstChunk; }
  
  virtual bool        SelectNextChunk(int slot_number, const std::string& chunkName) = 0;
  virtual bool        NextEvent()   = 0;
  virtual uint8*      GetRawEvent() = 0;
  virtual CsDstEvent* GetDstEvent() = 0;

  virtual const char* GetRunLog() const { return runLog.c_str(); }
  virtual void        SetRunLog(const char* log) { runLog = log; }
  virtual uint32      GetRawEventLength() const = 0;
  virtual uint32      GetEventInRun()     const = 0;
  virtual uint32      GetRunNumber()      const { return runNumber; }
  virtual uint32      GetEventInBurst()   const = 0;
  virtual uint32      GetBurstNumber()    const = 0;
  virtual uint32      GetTriggerMask()    const = 0;
  virtual uint32      GetErrorCode()      const = 0;
  virtual CsTime      GetTime()           const = 0;
  virtual bool        IsDstExist()        const = 0;

  virtual bool        UpdateRunLog(const std::string& tbnames,int slot) = 0;
  virtual bool        CloseDst() = 0;
  virtual void        Abort() = 0;

  virtual bool        Load(std::istream&);
  virtual bool        Load(FlatFile&);
  virtual bool        IsEquil(const CsDstObject*) const;
  virtual void        GetBuffer(std::vector<uint8>&) const;
};

#endif
