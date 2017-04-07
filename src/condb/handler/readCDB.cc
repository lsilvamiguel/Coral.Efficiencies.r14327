#include<iostream>
#include<time.h>
#include<unistd.h>
#include <iomanip.h>
#include "CondDbHandler.h"

#ifndef NO_OBJY
#include "oo.h"
#endif

time_t StringToTM(string time);

int main ( int argc, char* argv[] )
{
#ifndef NO_OBJY
  try {

    ooSetRpcTimeout(1000); //Default=25s
 
    //
    // CDB Handler Test Phase
    //
    if(argc<3) exit(-1);

    string bootPath=argv[1];

    string bTime_s = argv[2];
    string eTime_s = argv[3];

    int  bt = StringToTM(bTime_s);
    int  et = StringToTM(eTime_s);

    CDB::Time bTime(bt,0),eTime(et,0);

    list<string> detList;
    detList.push_back("FI01X1");
    CDB *cdb;

    cdb = new CondDbHandler(bootPath,detList,bTime,eTime);
    
    //Debug
    CondDbHandler *cdbH=dynamic_cast<CondDbHandler*>(cdb);
    if(cdbH!=0){
      cdbH->showFolder("/");
      cdbH->showFolderList();
      cdbH->showFolderSetList();
      cdbH->showCache();
    }

    //Read
    string folder = "/compass/FI01X1__/t0";
    vector<float> data;
    CDB::Time tPoint(bt,0);
    cdb->read(folder,data,tPoint);
    cout << "Data Size " << data.size()<<endl;
    for(unsigned int i=0;i<data.size();i++){
      cout << data[i]<< endl;
    }


    //Termination
    delete cdb;

    return 0; // return success
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

  struct tm t1;
  t1.tm_year = atoi(year.c_str());
  t1.tm_mon  = atoi(mon.c_str());
  t1.tm_mday  = atoi(day.c_str());
  t1.tm_hour = atoi(hour.c_str());
  t1.tm_min = atoi(min.c_str());
  t1.tm_sec = atoi(sec.c_str());

  t1.tm_year-= 1900;
  t1.tm_mon-=1;

  time_t t1_c = mktime(&t1);

  return t1_c;
}
