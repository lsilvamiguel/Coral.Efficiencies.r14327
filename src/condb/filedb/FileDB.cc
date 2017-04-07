#include "FileDB.h"
#include <unistd.h>
#include <sys/time.h>

using namespace std;
using namespace Reco;

FileDB::FileDB(string dbPath){
  d = new DataBase(dbPath,DataBase::READ);
}


void FileDB::read(const string &folder, string &data, Time timePoint, const char* keyword){
  tm *t = localtime((time_t*)&timePoint.first);

  string tmp = folder;
  string plane;

  //#define DEBUG_FileDB_read
#ifdef DEBUG_FileDB_read // Debugging the problem reported by Florent in 08/2006...
  // ...viz.: CsDC::readCalib gets empty calib from FileBD
  // Trying to determine what file name does fileDB look into 
  printf("folder: \"%s\" %d %d\n",tmp.c_str(),tmp.find("/"),tmp.size());
  // => Seems that the answer to the "find" query has change over successive
  //   versions of gcc: now (2006) it's = -1
#endif

  bool folderIsTBname;
  folderIsTBname = (int)tmp.find("/")==-1;
  if (folderIsTBname) //       ***** ARG folder IS JUST THE TBname *****
    plane = tmp;
  else {              // ***** ARG folder IS THE CALIB DIR NAME *****
    while (tmp.find("/")<tmp.size()){
      tmp.replace(tmp.find("/"),1," ");
    }
    istringstream is(tmp);
    for (int i = 0; i<3; i++) is>>plane;
  }

  if (keyword) plane += keyword;

  try{
      FILE* calibfd;
      char tmpbuf[8192];
      int nbread = 0;
      calibfd = fopen(d->CalibFilePath(plane, *t).c_str(), "r");
      while (!feof(calibfd)) {
        nbread = fread(tmpbuf, 1, 8192, calibfd);
        if (!nbread) { break;}
        data.append(tmpbuf, nbread);
      }
      fclose(calibfd);
  }
  catch(...) {
   // Dummy catch
//      cout << folder <<" calibrations, valid for ";
//      cout<<t->tm_mday<<"."<<t->tm_mon+1<<"."<<t->tm_year+1900<<" "
//  	<<t->tm_hour<<":"<<t->tm_min<<":"<<t->tm_sec
//  	<<", not found in DB"<<endl;

  }
}






