
/// to have strptime...
#define _XOPEN_SOURCE 1

#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <stdio.h>
#include <exception>
#include <ctype.h>
#include <time.h>

#include "MySQLInterface.h"
#include <mysqld_error.h>

using namespace std;

const string MySQLInterface::fCdbTimeFormat = "yyyy-MM-dd-hh:mm:ss";
const char*  MySQLInterface::fCdbTimeFormatStr = "%Y-%m-%d-%H:%M:%S";
const char*  MySQLInterface::fMySQLTimeFormatStr = "%Y-%m-%d %H:%M:%S";


MySQLInterface::MySQLInterface (const char* server, const char* username,
			        const char* password, const char* dbname)
  : fConnection(0), fConnected(false),
    fResult(0), fResInMem(false), fNRows(0) {

  if ( !mysql_init(&fMysql)) {
    fprintf(stderr, "MySQLInterface::MySQLInterface : %s\n", mysql_error(&fMysql));
    throw mysql_error(&fMysql);
  }
  if(!server || !username || !password || !dbname) {
    cerr<<"cannot create MySQL interface. Bad parameters"<<endl;
    throw "cannot create MySQL interface. Bad parameters";
  }

  fServer = server;
  fUserName = username;
  fPassword = password;
  fDBName = dbname;
  fNumPort = 0;

  cerr<<"creating MySQL interface to server "<<server<<" database "<<dbname<<endl;
}


MySQLInterface::~MySQLInterface() {
  disconnect();
}


//*********************** db interface ******************************

bool MySQLInterface::connect() {

  if (fConnected) {
    cerr << "MySQLInterface::connect: already connected\n";
    return true;
  }

  int trycount = 0;
  while (!fConnected) {
    fConnection = mysql_real_connect(&fMysql, fServer.c_str(),
				     fUserName.c_str(), fPassword.c_str(),
				     "", fNumPort, 0, 0);
    if (fConnection) {
      fConnected = true;
      fRow = 0;
      fResInMem = false;
      cerr<<"connected to DB "<<fDBName<<" @ "<<fServer <<" under "<<fUserName<<endl;
      return true;
    }

    register unsigned int merrno =  mysql_errno(&fMysql);
    if ((merrno != ER_CANT_CREATE_THREAD) && (merrno != ER_CON_COUNT_ERROR)
         && (merrno != ER_TOO_MANY_USER_CONNECTIONS) && (merrno != ER_USER_LIMIT_REACHED)) { break; }
    if (trycount > 10) {
      cerr<<"  MySQLInterface::connect: 10 tries failed, I give up\n";
      break;
    }
    cerr<<"  MySQLInterface::connect: "<<fServer<<" looks overloaded (error: "<<mysql_error(&fMysql)<<")\n";
    cerr<<"    waiting 30 seconds and trying again..."<<endl;
    trycount++;
    sleep(30);
  }

  fConnected = false;
  fprintf(stderr, "MySQLInterface::connect (mysql errno %u): %s\n", mysql_errno(&fMysql), mysql_error(&fMysql));
  return false;
}


bool MySQLInterface::query(const string& request)
{
  if ((!fConnected) && (!connect())) {
    return false;
  }
  if (fResInMem) endQuery();
  fNRows = 0;
  fRow = 0;
  if (mysql_query(fConnection, request.c_str())) {
    fprintf(stderr, "MySQLInterface::query : %s\n", mysql_error(&fMysql));
    return false;
  }
  fResult = mysql_store_result(fConnection);
  if(fResult) {
    fNRows = mysql_num_rows(fResult);
    fRow = mysql_fetch_row(fResult);
    fResInMem = true;
  }
  return true;
}


bool MySQLInterface::getNextRow(void) {
  if (!fConnected) return false;
  if (!fResInMem) return false;
  if (!fResult) return false;
  if (!fRow) return false;
  fRow = mysql_fetch_row(fResult);
  if (fRow) return true;
  return false;
}


void MySQLInterface::endQuery(void) {
  if (!fConnected) return;
  if (!fResInMem) return;
  if (!fResult) return;
  mysql_free_result(fResult);
  fResult = 0;
  fResInMem = false;
}


void MySQLInterface::disconnect(void) {
  if (!fConnected) return;
  endQuery();
  mysql_close(fConnection);
  fConnected = false;
  fRow = 0;
  fResInMem = false;
  cerr<<"disconnected from DB "<<fDBName<<" @ "<<fServer <<" under "<<fUserName<<endl;
}


string MySQLInterface::toMySQLtime(const char* cdbtime) {
  char strtmp[1024];
  struct tm stm;

  stm.tm_sec = stm.tm_min = stm.tm_hour = 0;
  stm.tm_mday = stm.tm_mon = stm.tm_year = 0;
  stm.tm_wday = stm.tm_yday = stm.tm_isdst = 0;
  strptime(cdbtime, fCdbTimeFormatStr, &stm);
  if (strftime(strtmp, 1024, fMySQLTimeFormatStr, &stm)) return strtmp;
    else return "";
}


string MySQLInterface::toMySQLtime(tm* stm) {
  char strtmp[1024];

  if (strftime(strtmp, 1024, fMySQLTimeFormatStr, stm)) return strtmp;
    else return "";
}


string MySQLInterface::toMySQLtime(time_t *timep) {
  tm *t = localtime(timep);
  return toMySQLtime(t);
}


string MySQLInterface::toCDBtime(const char* mysqltime) {
  char strtmp[1024];
  struct tm stm;

  stm.tm_sec = stm.tm_min = stm.tm_hour = 0;
  stm.tm_mday = stm.tm_mon = stm.tm_year = 0;
  stm.tm_wday = stm.tm_yday = stm.tm_isdst = 0;
  strptime(mysqltime, fMySQLTimeFormatStr, &stm);
  if (strftime(strtmp, 1024, fCdbTimeFormatStr, &stm)) return strtmp;
    else return "";
}



