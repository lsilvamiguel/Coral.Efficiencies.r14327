/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/mondb/src/monDB.cc,v $
   $Date: 2002/10/24 15:32:20 $
   $Revision: 1.7 $
   -------------------------------------------------------------------------

   This file is part of Compass monitoring program.

   Authors:
     based on original program of
     Vladimir  Kolosov   ( Kolosov@mx.ihep.su )
     Alexander Zvyagin   ( Alexander.Zviagine@cern.ch, Zvyagin@mx.ihep.su )
     adapted by Colin Bernet & Damien Neyret

   Copyright(C): 1999-2000  V.Kolosov,A.Zvyagin

     This library is free software; you can redistribute it and/or
     modify it under the terms of the GNU Library General Public
     License as published by the Free Software Foundation; either
     version 2 of the License, or (at your option) any later version.

     This library is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
     Library General Public License for more details.

     You should have received a copy of the GNU Library General Public
     License along with this library; if not, write to the Free
     Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <stdio.h>

#include "monDB.h"

// From man-page (man ctime):
//
// struct tm
// {
//         int     tm_sec;         /* seconds */
//         int     tm_min;         /* minutes */
//         int     tm_hour;        /* hours */
//         int     tm_mday;        /* day of the month */
//         int     tm_mon;         /* month */
//         int     tm_year;        /* year */
//         int     tm_wday;        /* day of the week */
//         int     tm_yday;        /* day in the year */
//         int     tm_isdst;       /* daylight saving time */
// };
// The members of the tm structure are:
//
// tm_sec The number of seconds after the minute, normally in the range 0 to 59, but can be up to 61 to allow for leap seconds.
//
// tm_min The number of minutes after the hour, in the range 0 to 59.
//
// tm_hour
//        The number of hours past midnight, in the range 0 to 23.
//
// tm_mday
//        The day of the month, in the range 1 to 31.
//
// tm_mon The number of months since January, in the range 0 to 11.
//
// tm_year
//        The number of years since 1900.
//
// tm_wday
//        The number of days since Sunday, in the range 0 to 6.
//
// tm_yday
//        The number of days since January 1, in the range 0 to 365.
//
// tm_isdst
//        A flag that indicates whether daylight saving time is in effect at the time described.  The value is  positive  if  day-
//        light saving time is in effect, zero if it is not, and negative if the information is not available.

namespace CSMon {

using namespace std;

////////////////////////////////////////////////////////////////////////////////

monDB::monDB (const char* server_, const char* username_, const char* dbname_) :
             fServer(server_), fUserName(username_), fDBName(dbname_),
             fConnection(0), fConnected(false),
             fResult(0), fResInMem(false), fNRows(0) {
  if ( !mysql_init(&fMysql) ) {
    printf("monDB::monDB : %s\n", mysql_error(&fMysql));
    delete this;
  }
}

monDB::~monDB() {
  disconnectDB();
}

////////////////////////////////////////////////////////////////////////////////

bool monDB::connectDB(void)
{
  if (fConnected) {
    cerr << "monDB::connectDB: already connected\n";
    return true;
  }

  fConnection = mysql_real_connect(&fMysql, fServer.c_str(), fUserName.c_str(), "",
				   fDBName.c_str(), 0, 0, 0);
  if (fConnection) {
    fConnected = true;
    fRow = 0;
    fResInMem = false;
    cout<<"connected to DB "<<fDBName<<" @ "<<fServer <<" under "<<fUserName<<endl;
    return true;
  }
  fConnected = false;
  printf("monDB::connectDB : %s\n", mysql_error(&fMysql));
  return false;
}

////////////////////////////////////////////////////////////////////////////////

bool monDB::queryDB(const char* request)
{
  if ((!fConnected) && (!connectDB())) {
    return false;
  }
  if (fResInMem) endQueryDB();
  fNRows = 0;
  fRow = 0;
  if (mysql_query(fConnection, request)) {
    printf("monDB::queryDB : %s\n", mysql_error(&fMysql));
    return false;
  }
  fResult = mysql_store_result(fConnection);
  fNRows = mysql_num_rows(fResult);
  fRow = mysql_fetch_row(fResult);
  fResInMem = true;
  return true;
}

////////////////////////////////////////////////////////////////////////////////

bool monDB::getNextRowDB(void)
{
  if (!fConnected) return false;
  if (!fResInMem) return false;
  if (!fResult) return false;
  if (!fRow) return false;
  fRow = mysql_fetch_row(fResult);
  if (fRow) return true;
  return false;
}

////////////////////////////////////////////////////////////////////////////////

void monDB::endQueryDB(void)
{
  if (!fConnected) return;
  if (!fResInMem) return;
  if (!fResult) return;
  mysql_free_result(fResult);
  fResult = 0;
  fResInMem = false;
}

////////////////////////////////////////////////////////////////////////////////

void monDB::disconnectDB(void)
{
  if (!fConnected) return;
  endQueryDB();
  mysql_close(fConnection);
  fConnected = false;
  fRow = 0;
  fResInMem = false;
  cout<<"disconnected from DB "<<fDBName<<" @ "<<fServer <<" under "<<fUserName<<endl;
}



////////////////////////////////////////////////////////////////////////////////

void monDB::PrintTime(const tm &t,ostream &o)
{
  o << TimeStr(t);
}


////////////////////////////////////////////////////////////////////////////////

const char* monDB::RootRefFile(const char* dettype) {

  if (!dettype) {
    cerr << "monDB::RootRefFile: warning, dettype char* is NULL\n";
    return 0;
  }
  string sqlreq;
  sqlreq = "select runnb,filename from tb_cooolref where tb_cooolref.dettype like '";
  sqlreq += dettype;
  sqlreq += "'";
  if(!queryDB(sqlreq)) {
    cout << " No SQL database answer about reference root file for detector type "<<dettype<<endl;
    return 0;
  }
  if (!getColDB(0)) {
    cout << " No reference root file registered for detector type "<<dettype<<" in SQL database\n";
    cout << "   taking default reference file\n";
    return 0;
  }

  string runnb(getColDB(0));
  string filename = "";
  if (getColDB(1)) filename = getColDB(1);
  endQueryDB();
  string name = "";
  if ( filename != "" ) {
    name += "/afs/cern.ch/compass/detector/monitor/References/";
    name += filename;
  } else {
    name += "/afs/cern.ch/compass/detector/monitor/References/coool_";
    name += runnb;
    name += ".root";
    struct stat statbuf;
    if ( ! stat(name.c_str(), &statbuf) ) {
      // stat succedeed
      if ( S_ISREG(statbuf.st_mode) ) {
        // regular file, we use it instead the file on web server (faster !)
        return name.c_str();
      }
    }
    name = "http://pccoeb03.cern.ch/rootfile.php/runnb=";
    name += runnb;
  }

  return name.c_str();
}


////////////////////////////////////////////////////////////////////////////////

bool operator == (const tm &t1, const tm &t2)
{
  return
    t1.tm_year  == t2.tm_year   &&
    t1.tm_mon   == t2.tm_mon    &&
    t1.tm_mday  == t2.tm_mday   &&
    t1.tm_hour  == t2.tm_hour   &&
    t1.tm_min   == t2.tm_min    &&
    t1.tm_sec   == t2.tm_sec;
}

////////////////////////////////////////////////////////////////////////////////

bool operator > (const tm &t1, const tm &t2)
{
  if( t1.tm_year > t2.tm_year )  return true;
  if( t1.tm_year < t2.tm_year )  return false;

  if( t1.tm_mon  > t2.tm_mon  )  return true;
  if( t1.tm_mon  < t2.tm_mon  )  return false;

  if( t1.tm_mday > t2.tm_mday )  return true;
  if( t1.tm_mday < t2.tm_mday )  return false;

  if( t1.tm_hour > t2.tm_hour )  return true;
  if( t1.tm_hour < t2.tm_hour )  return false;

  if( t1.tm_min  > t2.tm_min  )  return true;
  if( t1.tm_min  < t2.tm_min  )  return false;

  if( t1.tm_sec  > t2.tm_sec  )  return true;
  if( t1.tm_sec  < t2.tm_sec  )  return false;

  return false;
}

////////////////////////////////////////////////////////////////////////////////

string TimeStr(const tm &t)
{
  char s[111];
  sprintf(s,"%4.4d-%2.2d-%2.2d-%2.2d:%2.2d:%2.2d",
             t.tm_year+1900,t.tm_mon+1,t.tm_mday,
             t.tm_hour,t.tm_min,t.tm_sec );
  return string(s);
}

////////////////////////////////////////////////////////////////////////////////

string TimeStrSQL(const tm &t)
{
  char s[111];
  sprintf(s,"%4.4d%2.2d%2.2d%2.2d%2.2d%2.2d",
             t.tm_year+1900,t.tm_mon+1,t.tm_mday,
             t.tm_hour,t.tm_min,t.tm_sec );
  return string(s);
}


////////////////////////////////////////////////////////////////////////////////

} // namespace CSMon

////////////////////////////////////////////////////////////////////////////////









