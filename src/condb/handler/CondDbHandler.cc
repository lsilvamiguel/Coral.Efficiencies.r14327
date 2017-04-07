#include "CondDbHandler.h"
#include "CondDbHandlerError.h"
#include <unistd.h>
#include <sys/time.h>
#include <iostream>
#ifndef NO_OBJY
#include "oo.h"
#endif

//
//PUBLIC METHODE
//
//

CondDbHandler::CondDbHandler(const string &bootFile, const list<string> &detList, Time bTime, Time eTime) :\
  _bTime(bTime),_eTime(eTime),_openDB(true),_sizeOfCache(0)
{

  if(bTime.first == 0 && bTime.second == 0){
    cout << "CondDbHandler::CondDbHandler Time intervarl is zero!" << endl;
    exit(-1);
  }


#ifndef NO_OBJY

  ooSetRpcTimeout(100); //Default=25s

  #ifdef STANDALONE
  ooNoLock();
  #endif

  _condDBmgr = CondDBObjyDBMgrFactory::createCondDBMgr();
  _condDataAccess = _condDBmgr->getCondDBDataAccess();
  _condFolderMgr = _condDBmgr->getCondDBFolderMgr();
  _condDBmgr->init(bootFile);
  _condDBmgr->startRead();
  _condDBmgr->openDatabase();
  _condFolderMgr->getAllCondDBFolder(_folders);
  _condFolderMgr->getAllCondDBFolderSet(_folderSet);
  _transaction = true;

//    cout << bootFile << endl;
//    for(list<string>::const_iterator it=detList.begin();it!=detList.end();it++){
//      cout << *it << endl;
//    }
//    cout << bTime.first << endl;
//    cout << eTime.first << endl;
  open(detList);
  stopTransaction();

#else
  cout << "You need to recompile with Objectivity/DB" << endl;
#endif
}

#ifndef NO_OBJY

CondDbHandler::~CondDbHandler(){
  if(_openDB)close();
  if(_transaction)stopTransaction();
}

void CondDbHandler::stopTransaction(){
  //cout << "stop transaction" << endl;
  _condDBmgr->abort();
  CondDBObjyDBMgrFactory::destroyCondDBMgr( _condDBmgr );
  _transaction = false;
}

void CondDbHandler::open(const string &folder){

  try{
    if(isFolder(folder)) readFromFolder(folder);
    if(isFolderSet(folder)) readFromFolderSet(folder);
  }
  catch(CondDbHandlerError &e){
    e.Abort();
  }
}

void CondDbHandler::open(const list<string> &tbNameList){
  vector<string> folderList = constructFolderName(tbNameList);
  vector<string>::iterator it;
  try{
    for(it = folderList.begin();it!=folderList.end();it++){
      if(isFolder(*it)) readFromFolder(*it);
      if(isFolderSet(*it)) readFromFolderSet(*it);
    }
  }
  catch(CondDbHandlerError &e){
    e.Abort();
  }

}

void CondDbHandler::close(){
  clearCache();
  _openDB = false;
}

void CondDbHandler::close(const string &folder){
  clearCache(folder);
}


void CondDbHandler::showFolder(const string &folder) const{
  _condFolderMgr->describe(folder);
}

void CondDbHandler::showFolderList() const{
  vector<string>::iterator it = (string*)_folders.begin();
  cout << "Number of Condition Folder : " << _folders.size() << endl;
  cout << "### Begin of Folder List ###" <<endl;
  while(it !=(string*)_folders.end()){
    cout << *it++ << endl;
  }
  cout << "### End of FolderList ###" <<endl;
}

void CondDbHandler::showFolderSetList() const{
  vector<string>::iterator it = (string*)_folderSet.begin();
  cout << "Number of Condition Folder Set: " << _folderSet.size() << endl;
  cout << "### Begin of Folder Set List ###" <<endl;
  while(it !=(string*)_folderSet.end()){
    cout << *it++ << endl;
  }
  cout << "### End of Folder Set List ###" <<endl;
}

void CondDbHandler::dumpCache(){
  cout << "### DUMP CACHE " << endl;
  cout << "Number of Cached Object : " << _objMap.size() << endl;
  map<string,Interval>::iterator it;
  multimap<string,ICondDBObject*>::iterator p;
  CondDBKey since,till;
  string data;
  string description;
  for(it=_mapFolder.begin();it!=_mapFolder.end();it++){
    p = _objMap.find(it->first);
    if(p != _objMap.end()){
      do{
	since = p->second->validSince();
	till  = p->second->validTill();
	p->second->data(data);
	p->second->description(description);
	cout << p->first << " ";
	cout << since <<" ";
	cout << till << " ";
	cout << description << " ";
	cout << data << endl;
	p++;
      } while(p != _objMap.upper_bound(it->first));
    }
  }
}

void CondDbHandler::dumpCache(const string &outfName){
  cout << "### DUMP CACHE " << endl;
  cout << "Number of Cached Object : " << _objMap.size() << endl;
  ofstream f(outfName.c_str());
  map<string,Interval>::iterator it;
  multimap<string,ICondDBObject*>::iterator p;
  CondDBKey since,till;
  string data;
  for(it=_mapFolder.begin();it!=_mapFolder.end();it++){
    p = _objMap.find(it->first);
    if(p != _objMap.end()){
      do{
	since = p->second->validSince();
	till  = p->second->validTill();
	p->second->data(data);
	f << p->first << " ";
	f << since <<" ";
	f << till << " ";
	f << data << endl;
	p++;
      } while(p != _objMap.upper_bound(it->first));
    }
  }
}

void CondDbHandler::dumpCache(const string &folder, Time t){
  cout << "### DUMP CACHE " << endl;
  cout << "Time Point "<< t.first << endl;
  map<string,Interval>::iterator it;
  multimap<string,ICondDBObject*>::iterator p;
  CondDBKey since,till;
  string data;
  //for(it=_mapFolder.begin();it!=_mapFolder.end();it++){
  p = _objMap.find(folder);
  if(p != _objMap.end()){
    do{
      since = p->second->validSince();
      till  = p->second->validTill();
      p->second->data(data);
      p++;
    } while(!((since<t.first &&
	     till>t.first) ||
	      since == t.first) &&
	    p != _objMap.upper_bound(folder) );
    p--;
    cout << p->first << " ";
    cout << since <<" ";
    cout << till << " ";
    cout << data << endl;

  }

}

void CondDbHandler::showCache(){
  cout << "Total cached size : " << _sizeOfCache << endl;
}



void CondDbHandler::read(const string &Folder, string &data, Time tPoint, const char* keyword){
  try{

    string keyw = "";
    if (keyword) keyw = keyword;

    string folder;
    if (Folder.find("/") == Folder.size()) {    // Folder is just the TBname
      folder = "/COMPASS/";
      for(int i=0;i<2;i++) folder = Folder[i];
      folder += "/";
      for(int i=0;i<8;i++) folder += Folder[i];
      folder += keyw;
      folder += "/calibration";
    } else {            // folder is the real /COMPASS/TB/TBname/calibration
      folder = Folder;
    }

    //Remove "__"
    if(folder.find("__") < folder.size()){
      folder.replace(folder.find("__"),2,"");
    }

    tPoint.first -=2*60*60;
    map<string,Interval>::iterator it;
    multimap<string,ICondDBObject*>::iterator p;
    CondDBKey since,till;
    p = _objMap.find(folder);
    //cout << timePoint.secFrEpoch() << endl;
    if(p != _objMap.end()){
    do{
      since = p->second->validSince();
      till  = p->second->validTill();
      if(since<=tPoint.first && till>tPoint.first){
	//	cout << since <<" " << till<< endl;
	break;
      }
      p++;
    } while(p != _objMap.upper_bound(folder));
    if(!(since<=tPoint.first && till>tPoint.first)){
      ostrstream os;
      os<<"CondDBHandler : "<< folder << " : condition does not cached : t="<<tPoint.first<<ends;
      string tmp(os.str());
      throw CondDbHandlerError(tmp);
      //      cout << tmp << endl;
    }

    p->second->data(data);

    }
  }
  catch(CondDbHandlerError &e){
    e.Abort();
  }
}




//
//
//PROTECTED METHODE
//
//

bool CondDbHandler::getData(const string &folder,Interval interval){

//    string fname = detector + ".dat";
//    ifstream f(fname.c_str());
//    string tmp;
//    if(!f.is_open()){
//      tmp = "CondDbHandler:getData : "+ fname + " doesn't exist.";
//      throw tmp;
//    }

//    while(!f.eof()){
//      string item;
//      f>>item;
//      item = "/"+detector+"/"+item;
//      if(!readObjects(item,interval)){
//        ostrstream os;
//        os << "CondDbHandler:getData : "<< item<<" doesn't exist between from "
//  	  <<interval.getSince()<<" to "<<interval.getTill()<<" .";
//        tmp = os.str();
//        throw tmp;
//      }
//    }
  vector<string> folders;
  getSubFolder(folders,folder);
  if(folders.size() == 1) return true;
  else return false;
}

void CondDbHandler::readFromFolder(const string &folder){

  Interval interval(_bTime,_eTime);
  ICondDBObject* oCondObject = 0;
  CondDBKey bValid,eValid,bValidP,eValidP;
  CondDBKey tStep = 1;
  ostrstream os;
  string data;
  //Find first interval
  _condDataAccess->findCondDBObject(oCondObject,folder,interval.getSince().first);
  bValid = oCondObject->validSince();
  eValid = oCondObject->validTill();

  if(data == "") {
    //   throw CondDbHandlerException(folder,CondDBKey);
    os << "CondDbHandler : " << folder;
    os << " : No condition in " << interval.getSince().first;
    string tmp(os.str());
    //cout << tmp << endl;
    //throw CondDbHandlerError(tmp);
    data ="EMPTY";
  }
  oCondObject->data(data);

  addToCache(folder,oCondObject);

  bValidP = bValid;
  eValidP = eValid;
  CondDBKey tPoint = interval.getSince().first;
  //cout << "ev " << eValid << endl;
  //cout << "gt " <<  interval.getTill()<< endl;
  while(eValid<interval.getTill().first){
    //cout << "flag" << endl;
    _condDataAccess->findCondDBObject(oCondObject,folder,tPoint);
    bValid = oCondObject->validSince();
    eValid = oCondObject->validTill();
    oCondObject->data(data);
    if(data == ""){
      os << "CondDbHandler : " << folder;
      os << " : No condition in " << tPoint;
      string tmp(os.str());
      cout << tmp << endl;
    }
    //cout << tPoint << " "<< bValid << " " << eValid << endl;

    if(bValidP!=bValid && eValidP!=eValid){
      addToCache(folder,oCondObject);
    }
    bValidP = bValid;eValidP = eValid;
    //tPoint=tStep+eValidP;
    tPoint=eValidP;
  }

  cout << folder << " T["<<interval.getSince().first <<","<< interval.getTill().first<<") is cached." << endl;

  addFolderMap(folder,interval);

}

void CondDbHandler::readFromFolderSet(const string &folder){
  vector<string> folders;
  getSubFolder(folders,folder);

  if(folders.size() == 0) return;

  randmize(folders);//Randmize reading order
  vector<string>::iterator p;
  p=folders.begin();
  if(p != folders.end()){
    do{
      readFromFolder(*p);
      p++;
    }while(p != folders.end());
  }
}

void CondDbHandler::getSubFolder(vector<string>& folders, const string &folderSet){
  vector<string>::iterator it = (string*)_folders.begin();
  while(it !=(string*)_folders.end()){
    if((*it).find(folderSet,0)<(*it).size()){
      folders.push_back(*it);
    }
    it++;
  }
}

bool CondDbHandler::isFolder(const string &folder){
  vector<string>::iterator it = (string*)_folders.begin();
  while(it !=(string*)_folders.end()){
    if(*it == folder) return true;
    it++;
  }
  return false;
}

bool CondDbHandler::isFolderSet(const string &folder){
  vector<string>::iterator it = (string*)_folderSet.begin();
  while(it !=(string*)_folderSet.end()){
    if(*it == folder) return true;
    it++;
  }
  return false;
}


void CondDbHandler::randmize(vector<string> &folders){
  srand(getpid());
  //cout <<"PID " << getpid() <<endl;
  unsigned int swap1,swap2,roulete;
  vector<string>::iterator p;

//    for(p = folders.begin();p!=folders.end();p++){
//      cout << " randmize : " << *p << endl;
//    }
//    cout << "AFTER" << endl;

  float maxTurn = 1000.;
  roulete = (unsigned int) (maxTurn*rand()/(RAND_MAX+1.0));
  //  cout << "roulete " << roulete << endl;
  for(unsigned int i=0;i<roulete;i++){
    swap1=(unsigned int) ((float)folders.size()*rand()/(RAND_MAX+1.0));
    swap2=(unsigned int) ((float)folders.size()*rand()/(RAND_MAX+1.0));
    string tmp = folders[swap1];
    folders[swap1] = folders[swap2];
    folders[swap2] = tmp;

  }

//      int fSize = folders.size();


//    for(unsigned int i=0;i<diff;i++){
//      string tmp = folders[0];
//      folders.erase(folders.begin());
//      folders.push_back(tmp);
//    }

//    for(p = folders.begin();p!=folders.end();p++){
//      cout << " randmize : " << *p << endl;
//    }

}


CondDbHandler::Interval CondDbHandler::getCachedInterval(const string &foldername){
  map<string,Interval>::iterator it;
  for(it=_mapFolder.begin();it!=_mapFolder.end();it++){
    if(it->first == foldername) return it->second;
  }
  return Interval();
}

bool CondDbHandler::checkCachedInterval(const string &foldername, Interval interval){
  map<string,Interval>::iterator it;
  for(it=_mapFolder.begin();it!=_mapFolder.end();it++){
    //Checking provided interval is in cached intreval
    if(it->first == foldername &&
       (it->second).getSince() <= interval.getSince() &&
       (it->second).getTill() >= interval.getTill()){
      return true;
    }
  }
  return false;
}

void CondDbHandler::addFolderMap(const string & foldername, Interval interval){
  if(_mapFolder.count(foldername) != 0) _mapFolder.erase(foldername);

  _mapFolder.insert(pair<string,Interval>(foldername,interval));
}


bool CondDbHandler::readObjects(const string & foldername, Interval interval){

  for (CondDBKey t=interval.getSince().first; t<interval.getTill().first;t++){
    //Now increment unit is one

    ICondDBObject* oCondObject = 0;

    //    if(! condDataAccess->findCondObject(oCondObject,foldername,t)) return false;
    _condDataAccess->findCondDBObject(oCondObject,foldername,t);

    CondDBKey since = 0;
    CondDBKey till  = 0;

    string objData = "";

    since = oCondObject->validSince();
    till  = oCondObject->validTill();
    oCondObject->data(objData);

    addToCache(foldername,oCondObject);
  }

  addFolderMap(foldername,interval);

  cout << foldername << "["<<interval.getSince().first <<","<< interval.getTill().first<<") is cached." << endl;

  return true;
}

bool CondDbHandler::addToCache(const string &foldername, ICondDBObject *condObj){
  _objMap.insert(pair<string, ICondDBObject*>(foldername,condObj));
  string data;
  condObj->data(data);
  _sizeOfCache+=data.size();
  return true;
}

bool CondDbHandler::clearCache(){
  if(_mapFolder.empty()) return 0;
  multimap<string,ICondDBObject*>::iterator p;
  map<string,Interval>::iterator it;
  for(it=_mapFolder.begin();it!=_mapFolder.end();it++){
    p = _objMap.find(it->first);
    for(;p!=_objMap.upper_bound(it->first);p++){
      CondDBObjFactory::destroyCondDBObject(p->second);
    }
  }
  return 0;
}

bool CondDbHandler::clearCache(const string folder){
  if(_mapFolder.empty()) return 0;
  multimap<string,ICondDBObject*>::iterator p;
  p = _objMap.find(folder);
  if(p!=_objMap.upper_bound(folder)){
    do{
      string data;
      p->second->data(data);
      _sizeOfCache-=data.size();
      CondDBObjFactory::destroyCondDBObject(p->second);
      p++;
    }while(p!=_objMap.upper_bound(folder));
    _mapFolder.erase(folder);
  }else{return -1;}


  return 0;
}

ICondDBObject *CondDbHandler::find(const string &foldername,Time t){
  string tmp;
  ostrstream os;
  //
  //Check folder name is booked
  //
  if(getCachedInterval(foldername).getSince().first > t.first ||
       getCachedInterval(foldername).getTill().first < t.first){
    os<<"CondDbHandler:Provided time is out of interval ";
    tmp = os.str();
    throw CondDbHandlerError(tmp);
    //cout << tmp << endl;
  }

  CondDBKey since,till;
  multimap<string,ICondDBObject*>::iterator it,p;
  p = _objMap.find(foldername);
  for(it = p;it!=_objMap.end();it++){
    since = (it->second)->validSince();
    till  = (it->second)->validTill();
    if(since<=t.first && t.first<=till) return it->second;
  }
  return (ICondDBObject*)0;
}

vector<string> CondDbHandler::constructFolderName(const list<string> &tbNameList){
  //if(tbName.size() == 2){
  //cout << "*** F3" << endl;

  vector<string> folderList;
  list<string>::const_iterator it1;
  vector<string>::const_iterator it2;
  string filter;
  for(it1 = tbNameList.begin();it1!=tbNameList.end();it1++){
    //cout << "it1 CDN " << *it1 << endl;
    filter="/COMPASS/"+(*it1);
    for(it2 = _folders.begin();it2!=_folders.end();it2++){
      //size_t pos = (*it2).find((*it1),0);
      if((*it2).find(filter,0)<(*it2).size()){
	//cout << "it2 CDN " <<*it2 << endl;
	folderList.push_back(*it2);
      }
    }
  }
    //}
  //_folderSet
  return  folderList;
}


//catch CondDbHandlerException
#endif
