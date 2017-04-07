// protptype of CDB handler

#ifndef CondDbHandler_h
#define CondDbHandler_h


#include <string>
#include <iostream>
#include <fstream>
#include <iostream>
#include <strstream.h>
#include <map>
#include <list>

#ifndef NO_OBJY
#include <ConditionsDB/ICondDBMgr.h>
#include <ConditionsDB/CondDBObjyDBMgrFactory.h>
#include <ConditionsDB/CondDBObjFactory.h>
#endif

#include "CDB.h"

class CondDbHandler : public CDB{
public:

#ifndef NO_OBJY
  class Interval{
  public:
    Interval() : since_(0,0),till_(0,0){}
    Interval(Time since, Time till) : since_(since),till_(till){}
    Time getSince(void) const{return since_;}
    Time getTill(void) const{return till_;}
  protected:
    Time since_;
    Time till_;
  };

  CondDbHandler(const string &bootFile, const list<string> &detList,Time bTime,Time eTime);
  ~CondDbHandler();
  void open(const string &folderName);
  void open(const list<string> &tbNameList);
  void close();
  void stopTransaction();
  void read(const string &folder, string &data, Time timePoint, const char* keyword=0);

  void close(const string &folder);
  void showFolder(const string &folder) const;
  void showFolderList() const;
  void showFolderSetList() const;
  void dumpCache();
  void dumpCache(const string &outfName);
  void dumpCache(const string &folder, Time timePoint);
  void showCache();

protected:

  bool getData(const string &detector,Interval interval);
  void getSubFolder(vector<string>& folders, const string &folderSet);
  Interval getCachedInterval(const string &foldername);
  bool readObjects(const string & foldername,Interval interval);
  bool checkCachedInterval(const string &foldername,Interval interval);
  void addFolderMap(const string &foldername,Interval interval);
  //vector<string> getFolderListFromMap(string folderser);
  bool clearCache();
  bool clearCache(const string folder);
  string connectString();
  bool isFolder(const string &folder);
  bool isFolderSet(const string &folder);
  void readFromFolder(const string &folder);
  void readFromFolderSet(const string &folder);
  void randmize(vector<string> &folders);
  vector<string> constructFolderName(const list<string> &tbNameList);

  bool addToCache(const string &foldername,ICondDBObject *condObj);
  ICondDBObject *find(const string &foldername,Time t);

  Time _bTime,_eTime;

  multimap<string, ICondDBObject*> _objMap;

  map<string, Interval> _mapFolder;
  vector<string> _folders;
  vector<string> _folderSet;
  long _sizeOfCache;

  //Flags
  bool _openDB;
  bool _transaction;

  //ConditionsDB packages
  ICondDBMgr *_condDBmgr;
  ICondDBFolderMgr *_condFolderMgr;
  ICondDBDataAccess *_condDataAccess;

#else

public:
  CondDbHandler(const string &bootFile, const list<string> &detList,Time bTime,Time eTime);
  ~CondDbHandler(){}

  virtual void read(const string &folder, vector<float> &data, Time timePoint){}
  virtual void read(const string &folder, vector< vector<float> > &data, Time timePoint){}
  virtual void read(const string &folder, string &data, Time timePoint){}

protected:
  Time _bTime,_eTime;
  long _sizeOfCache;

  //Flags
  bool _openDB;
  bool _transaction;

#endif
};
#endif




