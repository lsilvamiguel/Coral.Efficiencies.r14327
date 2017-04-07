
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <pthread.h>

#include "CollectErrs.h"

#include "monitor.h"

using CS::uint32;
using CS::DaqError;
using CS::DaqOption;
using CS::ObjectXML;
using CS::Chip;


bool Connect(const char *source){
  static char *table[6]={"All","no","PHY","yes",0};
 
  int a=monitorSetDataSource((char *)source);
  if (0==a) {
    monitorSetNowait();
    a=monitorDeclareMp("errorDumpAll");
    if (a==0) {
	a=monitorDeclareTable(&table[0]);  //We only want physics events
      if (a==0) {
        return true;
      }
      else printf("Date error:%s \n",monitorDecodeError(a));
    } 
    else printf("Date error:%s \n",monitorDecodeError(a));
  } 
  else printf("Date error:%s \n",monitorDecodeError(a));
 
  return false;
}


int main(int argc,char **argv){ 
  char *mfile=0;
  char *ofile=0;
  int nrevents=-1;
  bool qar=false;
  vector <string> datafiles;
  set <int> srcids;

  Chip::SetFillStatistics(false);


  bool se=false; 
  for (int a=1;a<argc;a++) {
    size_t al=strlen(argv[a]);
    char *ra=argv[a];
    char b='s';
    if (argv[a][0]=='-') {
      if (al<2) {se=true;break;}
      if (argv[a][1]=='q') {
        if (al==2) qar=true; else {se=true;break;}
        continue;
      }
      b=argv[a][1];
      if (al<=2) {
        a++;
        if (a>=argc) {se=true;break;}
        ra=argv[a];
      } else ra=argv[a]+2;
    }      
   
    if (b=='s') datafiles.push_back(ra);
    else if (b=='m') mfile=ra;
    else if (b=='o') ofile=ra;
    else if (b=='n') nrevents=atoi(ra);
    else if (b=='i') srcids.insert(atoi(ra));
    else {se=true;break;}

  }
  

  se= se || datafiles.size()==0 || mfile==0;
  if (se) { 
    cout << "usage: errorDumpAll -s <data source> -m <mapping directory> [-o output file]" <<endl
	 << "[-n #events] [-i SrcIDs] [-q]" << endl
         << "-q means exit when run finished (for online)" << endl;
    return 0;
  }

 
 
  uint32 runn=0;
  CollectErrs *trash=new CollectErrs(mfile);
  trash->ResetAtNewRun(false);
  for (vector<string>::iterator fni=datafiles.begin();fni!=datafiles.end();fni++) {
    if (Connect(fni->c_str())) {
      char *ptr;
      int en=0;
      int st;
      uint32 waiting=0;
      bool brk=false;
      while ((en<nrevents || nrevents<0) && !brk) {
      st=monitorGetEventDynamic((void **)&ptr); //get raw event from DATE 
      if (0==st) {
      if (0==ptr) {
    	    usleep(100000); // 100 msec. sleep, prevent too high CPU usage when idle
	    waiting++;
	    if (waiting>108000) brk=true; //quit after 3 hours waiting
	    }
      else { 
        en++;
        DaqEvent *evnt=new DaqEvent(ptr);
        if (runn==0) runn=evnt->GetRunNumber();
        if (runn==evnt->GetRunNumber() || !qar) {
          trash->DecodeEvent(evnt);
          trash->HandleNewEvent(evnt);
        } else brk=true;
        delete evnt;
        free(ptr); 
	}  
      } else {
         cerr << "Date error: " << monitorDecodeError(st);
	 brk=true;
      }
      }
    }
  }
  
  
  if (ofile) {
    ofstream out(ofile);
    if (out) {
      trash->PrintReport(out,srcids);
    }
  }
  else trash->PrintReport(cout,srcids);

  monitorLogout();  
}
