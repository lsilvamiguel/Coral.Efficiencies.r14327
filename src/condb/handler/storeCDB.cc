// Basic example of data storage.

#ifndef NO_OBJY  
#include "oo.h"
#include "ConditionsDB/ICondDBMgr.h"
#include "ConditionsDB/CondDBObjyDBMgrFactory.h"
#include "ConditionsDB/CondDBObjFactory.h"
#endif

#include <string>
#include <iostream>
#include <fstream>
#include <strstream>
#include <unistd.h>
#include <sys/param.h>
#include <stdio.h>
#include <time.h>

// forward def
time_t StringToTM(string time);

int main ( int argc, char* argv[] )
{
#ifndef NO_OBJY

  try {

    if(argc<5) exit(-1);
    string bootFile = argv[1];
    string folder  = argv[2];

    string cdbfile;
    cdbfile = argv[3];
    

    string bTime_s = argv[4];
    string eTime_s = argv[5];

    int  bTime = StringToTM(bTime_s);
    int  eTime = StringToTM(eTime_s);

    long time = 0;


#ifdef STANDALONE
    ooNoLock();
#endif

    ICondDBMgr* condDBmgr;
    ICondDBDataAccess* condDataAccess;
    ICondDBFolderMgr*  condFolderMgr; 
    ICondDBObject* condObject;
    condDBmgr = CondDBObjyDBMgrFactory::createCondDBMgr();
    condDBmgr->init(bootFile);

    condDataAccess = condDBmgr->getCondDBDataAccess();
    condFolderMgr  = condDBmgr->getCondDBFolderMgr(); 
    condDBmgr->startUpdate();    
    condDBmgr->openDatabase();

    float element;
    ostrstream os;
    ifstream ifs;
    ifs.open(cdbfile.c_str());
    if(!ifs.is_open()) cout << "File is not opened." << endl;
    os<<ifs.rdbuf()<<ends;
    string data = os.str(); 
    cout << data << endl;
    string description;
    condObject = CondDBObjFactory::createCondDBObject(bTime, eTime,
						    data,description);
    condDataAccess->storeCondDBObject( folder,
					 condObject );
    CondDBObjFactory::destroyCondDBObject(condObject);

//    ostrstream os0;
//    os0<<"/COMPASS/RICH/xxx"<<folderId<<ends;
//    string folder(os0.str());
//    cout << folder << endl;
//    ostrstream os1;
    //os1<<"xxx"<<folderId<<ends;
    //string description(os1.str());
//      for(long long j=0;j<uTime;j+=timeStep){    
//        bTime=j;eTime=j+timeStep;
//        ostrstream os2;
//        for(int k=0;k<10000;k++){
//  	os2 << folderId+j+k+1 << " ";
//        }
//        os2 << ends;
//        string data(os2.str());
	//  	        cout << folder << " " << bTime << " " << eTime << " ";
	//  	        cout << description << " "<<data<<endl;
   // }
    condDBmgr->commit();
    
    return 0;
  }
  catch (CondDBException &e)
    {
      cerr << "*** ConditionsDB exception caught: " << e.getMessage() << "\n"
           << "***   error code: " << e.getErrorCode() << endl;
      
      return 1;  // return failure
    }

  catch(char *pErr){
    cout << "Error:" << pErr << endl; 
  }
  catch(string str){
    cout << "Error:" << str << endl; 
  }
#else
    cout << " You need to compile with Objectivity/DB!" << endl;
#endif

}


time_t StringToTM(string time){
  string year, mon, day, hour, min, sec;

  for(int i=0;i<4;i++) year+=time[i];
  for(int i=5;i<7;i++) mon+=time[i];
  for(int i=8;i<10;i++) day+=time[i];
  for(int i=11;i<13;i++) hour+=time[i];
  for(int i=14;i<16;i++) min+=time[i];
  for(int i=17;i<20;i++) sec+=time[i];

//    cout << year << endl;
//    cout << mon << endl;
//    cout << day << endl;
//    cout << hour << endl;
//    cout << min << endl;
//    cout << sec << endl;

  struct tm t1;
  t1.tm_year = atoi(year.c_str());
  t1.tm_mon  = atoi(mon.c_str());
  t1.tm_mday  = atoi(day.c_str());
  t1.tm_hour = atoi(hour.c_str());
  t1.tm_min = atoi(min.c_str());
  t1.tm_sec = atoi(sec.c_str());
//    cout << "########## Convert ########" << endl;
  
//    cout << t1.tm_year << endl;
//    cout << t1.tm_mon << endl;
//    cout << t1.tm_mday << endl;
//    cout << t1.tm_hour << endl;
//    cout << t1.tm_min << endl;
//    cout << t1.tm_sec << endl;


  t1.tm_year-= 1900;
  t1.tm_mon-=1;

  time_t t1_c = mktime(&t1);

  return t1_c;
}
