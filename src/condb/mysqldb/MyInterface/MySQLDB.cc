
#include <cstring>

#include "MySQLDB.h"

using namespace std;

MySQLDB::MySQLDB (const char* server, const char* username, const char* passwd, const char* dbname) :
           CDB(), MySQLDBInterface(server, username, passwd, dbname),
           fEntrytime(""), fSpecETfg(false) {

  fSpecificEntrytime.clear();
}


//----------------------------------------------------------------------------

bool MySQLDB::ConnectDB(void) {
  register bool tmpflag;
  register int retcode = 0;

  tmpflag = connect();
  if (!tmpflag) {
    cerr << "Error in MySQLDB::ConnectDB(): can't connect to DB"<<endl;
    return false;
  }
  if ((retcode = selectDB())) {
    cerr << "Error in MySQLDB::ConnectDB(): can't access DB, error code "<<retcode<<endl;
     return false;
 }
 return true;
}


//----------------------------------------------------------------------------


void MySQLDB::read(const string &folder, string &data, Time timePoint, const char* keyword)
{
  string tmp = folder;
  string plane;
  string str_et = fEntrytime;
  if (tmp.find("/") == tmp.size()) {    // folder is just the TBname
    plane = tmp;
  } else {
    while(tmp.find("/")<tmp.size()){    // folder is the real /COMPASS/TB/TBname/calibration
      tmp.replace(tmp.find("/"),1," ");
    }
    istringstream is(tmp);
    for(int i=0;i<3;i++) is >> plane;
  }

  if (fSpecETfg) {      // there are specific entry times
    typedef std::map<std::string, std::string>::iterator SpecETit;
    SpecETit ispecit = fSpecificEntrytime.find(plane);
    if (ispecit != fSpecificEntrytime.end()) {
      str_et = ispecit->second;
      std::cerr << "MySQLDB::read: using specific entry time "<<str_et<<" for TBName "<<plane<<endl;
    }
  }

  string datestr = toMySQLtime((time_t*)(&timePoint.first));

  string filepath = giveFilepath(plane.c_str(), datestr.c_str(), keyword, str_et.c_str());

  if (filepath == "") return;   // No corresponding entry found in database

  try{
      FILE* calibfd;
      char tmpbuf[8192];
      register int nbread = 0;
      calibfd = fopen(filepath.c_str(), "r");
      if ( !calibfd ) {
        std::cerr<<"MySQLDB::read: can't open calibration file "<<filepath<<endl;
        std::cerr<<"               ... returning empty string"<<endl;
        return;
      }
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




//----------------------------------------------------------------------------

string MySQLDB::decodeTBNameFromFolder(const string& folder) {

  const char *fidx = folder.c_str();
  for (register int ii=0; ii < 3; ii++) {
    fidx = strchr(fidx, '/');
    if (fidx == 0) {
      cerr <<"MySQLDB::decodeTBNameFromFolder: no '/' found in folder string ";
      cerr <<folder<<endl;
      return "";
    }
    fidx++;
    if (*fidx == 0) {
      cerr <<"MySQLDB::decodeTBNameFromFolder: no TBname found in folder string ";
      cerr <<folder<<endl;
      return "";
    }
  }
  return string(fidx, 8);
}

