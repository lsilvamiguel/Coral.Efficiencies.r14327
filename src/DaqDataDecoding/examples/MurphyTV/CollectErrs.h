#include "config.h"
#include <stdlib.h>
#include <pthread.h>
#include <set>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <math.h>
using CS::uint32;
using CS::uint16;
using CS::uint8;
using namespace std;

#include "DaqEvent.h"
#include "DaqOption.h"
#include "DaqError.h"
#include "ChipGandalf.h"
using CS::DaqEvent;
using CS::DaqError;
using CS::DaqErrors;
using CS::Exception;
using CS::DaqOption;
using CS::Chip;
using CS::ChipGandalf;

typedef DaqError::Type DaqErrorType; 

/**************************************************************************************************************
This class is used for bookkeeping of the errors. 
**************************************************************************************************************/
class CollectErrs{
  public:
  
  enum oink {MAX_SOURCE=1024,MAX_GEOID=1024,MAX_PORT=16,MAX_VPORT=36,INVALID=2000,R_MAX_TYPE=DaqError::MAX_TYPE};

  private:

  enum iek {
    HISTBUF_SIZE=500,
    HISTBUF_RES=100,
  };
  public:
  class VPortID {    //to store a Portnumber/GeographicID pair
  public:
    VPortID(uint16 port,uint16 geoid){id=geoid | (port<<16);}    
    bool operator == (const VPortID &vp) {return (vp.id==id);}
    uint16 Port() const {return id>>16;}
    uint16 GeoID() const {return id & 0xffff;}      
    uint32 id;
  };
  
  //to store statistics on number of header and data words, total number of errors
  //and number of errors by errortype. 
  class VPortStat {  
  public:
    VPortStat() {NumData=NumHeader=NumAnyErr=0;}
    VPortStat& operator +=(const VPortStat &ps);
    uint32 GetNumErr(DaqError::Type typ) const; 
    map <DaqError::Type,uint32> NumErr;
    uint32 NumData,NumHeader,NumAnyErr;
  };
  typedef map <VPortID,VPortStat> VPortStatMap;

  private:
  class CatchStat { //to store statistics for one SourceID
  private:

  public:
    enum mune {SIZEHIST_SIZE=640, SIZEHIST_RES=20};

    CatchStat();

    //These increment number of header or data words..
    //..in general on this CATCH
    void AddData(uint32 n) {Data+=n;} 
    //..on a given PortNr/GeoID on this CATCH
    void AddHD(uint16 port,uint16 geoid,uint32 header,uint32 data); 
    
    //This adds one Error of type typ, occured at TriggerNr trnr on this CATCH  
    void AddError(uint32 trnr, uint32 spillnb, DaqError::Type typ,uint16 port,uint16 geoid);
    
    //Returns a map with Statistics for the Port/GeoID pairs
    const map <VPortID,VPortStat> &GetVPortStats() {return VPortStats;}

    //Fill the eventsize-histo 
    void FillSizeHist(uint32 size);
 
    //Number of data words   
    uint32 GetData() {return Data;}

    //Number of header words 
    uint32 GetHeader() {return Header;}

    //Total number of errors 
    uint32 GetErrSum() {return ErrSum;}

    //Number of events which contained one ore more errors 
    uint32 GetNrBadEvnts() {return BadEvnts;}

    //Number of errors of type typ 
    uint32 GetErr(DaqError::Type typ) {if (typ<(DaqError::Type)R_MAX_TYPE) return Err[typ]; else return 0;}

    //Number of events which contained one ore more errors of type typ 
    uint32 GetEvntsWithErr(DaqError::Type typ) {if (typ<(DaqError::Type)R_MAX_TYPE) return EvntsWithErr[typ]; else return 0;}

    //TriggerNr when error typ occured the first time 
    uint32 GetFirstErrPos(DaqError::Type typ) {if (typ<(DaqError::Type)R_MAX_TYPE) return FirstErrPos[typ]; else return 0;}

    //Get bin bin of the eventsize histo 
    uint32 GetSizeHist(uint16 bin) {if (bin<SIZEHIST_SIZE) return SizeHist[bin]; else return 0;}

    //Number of spills with events affected by an error on this SrcID
    uint32 GetNbSpillWithErrGen() {return NbSpillWithErrGen;}

    //Last spill with events affected by an error on this SrcID
    uint32 GetLastSpillWithErrGen() {return LastSpillWithErrGen;}

    //Number of spills with events affected by error of type typ 
    uint32 GetNbSpillWithErr(DaqError::Type typ) {if (typ<(DaqError::Type)R_MAX_TYPE) return NbSpillWithErr[typ]; else return 0;}

    //Last spill with events affected by error of type typ 
    uint32 GetLastSpillWithErr(DaqError::Type typ) {if (typ<(DaqError::Type)R_MAX_TYPE) return LastSpillWithErr[typ]; else return 0;}

    //To be called before/after parsing of event
    void BeginEvent();
    void EndEvent();
    
    private:

    VPortStatMap::iterator  FindVPort(uint16 port,uint16 geoid);
    
    uint32 Data;    //Number of data words 
    uint32 Header;  //guess
    uint32 ErrSum;  //Number of errors
    
    bool ErrEx[R_MAX_TYPE];  //which errors were found in this event ?

    uint32 BadEvnts;  //Number of events which had one or more errors

    uint32 SizeHist[SIZEHIST_SIZE];  //eventsize-histo (array instead of ROOThisto to save memory)

    uint32 Err[R_MAX_TYPE];          //Number of errors of each type
    uint32 EvntsWithErr[R_MAX_TYPE]; //Number of events wich contained one ore more error of a given type
    uint32 FirstErrPos[R_MAX_TYPE];  //TriggerNr when a given error type occured the first time

    VPortStatMap VPortStats;  //Statistics for each PortNr/GeoID

    uint32 NbSpillWithErr[R_MAX_TYPE];  //Number of spills with events affected by each error
    uint32 LastSpillWithErr[R_MAX_TYPE];  //Last spill with events affected by each error
    uint32 NbSpillWithErrGen;  //Number of spills with events affected by an error
    uint32 LastSpillWithErrGen;  //Last spill with events affected by an error
  };
          
  public:

  
  CollectErrs(const char *mapfile=0);
  ~CollectErrs();
  
  //hand event to decoding library, read maps if needed 
  void DecodeEvent(DaqEvent *event);    
  
  //count the errors given back by decoding library
  void HandleNewEvent(DaqEvent *event);
  
  //tell decoding library to read its maps and option file.
  void ReadMaps(const char *mapfile=0);


  uint32 GetRunNr() { return RunNr;}              //RunNr of last recorded event
  uint32 GetEventTime() {return EventTime;}      
  uint32 GetEventSpillNb() {return EventSpillNb;}      
  uint32 GetNrRecEvents() {return NrRecEvents;}   //Number of events recorded so far
  uint32 GetParsedSize() {return ParsedSize;}     //Total amount of data parsed so far
  
  
  //unimportant
  uint32 GetHistBufRes() {return HISTBUF_RES;}
  uint32 GetHistBufSize() {return HISTBUF_SIZE;}
  uint32 GetHistPos() {return HistBufPos;}
  uint32 GetErrHist(DaqError::Type typ,uint32 bin) {if (bin<HistBufPos && typ<(DaqError::Type)R_MAX_TYPE) return ErrHist[bin][typ]; else return 0;}
  uint32 GetRecEventsHist(uint32 bin) {if (bin<HistBufPos) return RecEventsHist[bin]; else return 0;}
  uint32 GetEventSizeHist(uint32 bin) {if (bin<HistBufPos) return EventSizeHist[bin]; else return 0;}

  uint32 GetSizeHistSize() {return CatchStat::SIZEHIST_SIZE;}
  uint32 GetSizeHistRes() {return CatchStat::SIZEHIST_RES;}
  uint32 GetSizeHistOfCatch(uint16 source,uint16 bin) {return GetCatchNC(source)->GetSizeHist(bin);}
  
  //Total number of errors
  uint32 GetErrSum();
  
  //Number of errors of type typ
  uint32 GetErr(DaqError::Type typ) {if (typ<(DaqError::Type)R_MAX_TYPE) return Err[typ]; else return 0;}


  uint32 GetLastErrs(DaqError::Type typ) {if (HistBufPos>0 && typ<(DaqError::Type)R_MAX_TYPE) return ErrHist[HistBufPos-1][typ]; else return 0;}
  uint32 GetLastEvents(){if (HistBufPos>0) return RecEventsHist[HistBufPos-1]; else return 0;}

  //Number of events which contained one ore more errors of type typ
  uint32 GetNrEvntsWithErr(DaqError::Type typ) {if (typ<(DaqError::Type)R_MAX_TYPE) return NrEvntsWithErr[typ]; else return 0;}

  //TBNames of detectors connected to SourceID source, mapping files have to be present for this to succeed
  void GetDetsAtCatch(uint16 source,set <string> &detnames);
  
  //return description for error type err
  string GetErrorName(DaqError::Type err);
  
  //Number of errors of type typ on SourceID source
  uint32 GetErrOnCatch(uint16 source,DaqError::Type typ) {return GetCatchNC(source)->GetErr(typ);}

  //Number of events which contained one ore more errors of type typ on SourceID source
  uint32 GetEvntsWithErrOnCatch(uint16 source,DaqError::Type typ) {return GetCatchNC(source)->GetEvntsWithErr(typ);}

  //TriggerNr when error typ occured the first time on SourceID source
  uint32 GetFirstErrPosOnCatch(uint16 source,DaqError::Type typ) {return GetCatchNC(source)->GetFirstErrPos(typ);}

  //Number of errors on SourceID source
  uint32 GetCatchErrSum(uint16 source) {return GetCatchNC(source)->GetErrSum();}

  //Number of events whith one ore more error on SourceID source
  uint32 GetCatchNrBadEvnts(uint16 source) {return GetCatchNC(source)->GetNrBadEvnts();}
  
  //Number of data/header words from SourceID source
  uint32 GetNrDataOnCatch(uint16 source) {return GetCatchNC(source)->GetData();}
  uint32 GetNrHeaderOnCatch(uint16 source) {return GetCatchNC(source)->GetHeader();}
  
  //Number of spills with events affected by errors of type typ on SourceID source, and last spill like that
  uint32 GetNbSpillWithErr(uint16 source,DaqError::Type typ) {return GetCatchNC(source)->GetNbSpillWithErr(typ);}
  uint32 GetLastSpillWithErr(uint16 source,DaqError::Type typ) {return GetCatchNC(source)->GetLastSpillWithErr(typ);}
  
  //Number of spills with events affected by errors on SourceID source, and last spill like that
  uint32 GetNbSpillWithErrGen(uint16 source) {return GetCatchNC(source)->GetNbSpillWithErrGen();}
  uint32 GetLastSpillWithErrGen(uint16 source) {return GetCatchNC(source)->GetLastSpillWithErrGen();}

  //Statistics for all Port/GeoID on SourceID source
  const map <VPortID,VPortStat> &GetVPortStats(uint16 source) {return GetCatchNC(source)->GetVPortStats();}

  uint32 GetRunCount() {return RunCount;}
  bool CatchExists(uint16 source) {if (source<MAX_SOURCE) return (Catches[source]!=0); else return false;}

  //next time DecodeEvent is called all counters are reset to zero
  void ResetAtNextEvent() {Rne=true;}
  //next time a new run starts and DecodeEvent is called all counters are reset to zero
  void ResetAtNewRun(bool doreset=true) {Rnr=doreset;}
  
  //Prints a report to stream out (using only informations about the SourceIDs in srcids,
  //if the set is not empty)
  void PrintReport(ostream &out,set<int> srcids);
  
  //Reset all counters
  void Reset();

  
  //When HandleNewEvent is called by one thread and Information from this class is used by another,
  //these functions have to be called for synchronisation
  void Lock() {pthread_mutex_lock(&sync);}
  void Unlock() {pthread_mutex_unlock(&sync);}

  private:
  
  void Check4NewRun(DaqEvent *event);
  void HandleNewError(uint32 trnr,const DaqError &err);
  CatchStat *GetCatch(uint16 source) {if (!Catches[source]) Catches[source]=new CatchStat(); return Catches[source];}
  CatchStat *GetCatchNC(uint16 source);
 
  void PrintErr(ostream &out,const string &prefix,int nr,int bad_events);
  
  CatchStat defctch;               //Dummy CATCH for internal purposes
  CatchStat *Catches[MAX_SOURCE];  //Statistics for each possible SourceID
  
  uint32 Err[R_MAX_TYPE];                   //number of errors of each type
  uint32 NrEvntsWithErr[R_MAX_TYPE];        //number of events with one or more errors of a given type
  uint32 ErrHist[HISTBUF_SIZE][R_MAX_TYPE]; //number of errors versus trigger number
  uint32 RecEventsHist[HISTBUF_SIZE];       //number of recorded events in each trigger number bin
  uint32 EventSizeHist[HISTBUF_SIZE];       //amount of recorded data   ""
  
  uint32 HistBufPos,NrRecEvents,TriggerNr,RunNr,EventTime,ParsedSize,RunCount,EventSpillNb;

  bool ErrEx[R_MAX_TYPE];  //did this error occur in this event ?
  bool Catch_in_event[MAX_SOURCE]; //was a Catch present in the actual event?
  bool Rne,Rnr; 
  DaqOption opts;     //decoding library options, read in with mapping files (DAQ.xml).
  Chip::Maps maps;    //mapping files, optional, needed for detector name display
  bool mapsread;      
  
  pthread_mutex_t sync;   //needed to sync acces to this class when used multithreaded
         
};

bool operator <(const CollectErrs::VPortID &id1,const CollectErrs::VPortID &id2);

