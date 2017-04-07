/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/DataBase.cc,v $
   $Date: 2010/09/21 10:27:57 $
   $Revision: 1.15 $
   -------------------------------------------------------------------------

   This file is part of cellular calorimeter reconstruction program.

   Authors:
     Vladimir  Kolosov   ( Kolosov@mx.ihep.su )
     Alexander Zvyagin   ( Alexander.Zviagine@cern.ch, Zvyagin@mx.ihep.su )

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

#include "DataBase.h"

using namespace std;

namespace Reco {

////////////////////////////////////////////////////////////////////////////////

DataBase::DataBase(const string &path_, int open_mode_) :
  path          (path_),
  open_mode     (open_mode_),
  creation_flag (false)
{
  if( open_mode&CREATE && !(open_mode&READ) )
    system( string("rm -rf " + path).c_str() );

  if( open_mode&(CREATE|WRITE) )
  {
    mode_t mode=0777;
    int error = mkdir(path.c_str(),mode);
    if( error==0 )
      creation_flag = true;
//#warning TODO  test that directory exists or was created
//     if( error )
//       throw Exception("DataBase::DataBase():  can not create directory \"%s\"   error code %d",path.c_str(),error);
  }
}

////////////////////////////////////////////////////////////////////////////////

string DataBase::TimeStr(const tm &t)
{
  char s[111];
  sprintf(s,"%4.4d-%2.2d-%2.2d-%2.2d:%2.2d:%2.2d",
             t.tm_year+1900,t.tm_mon+1,t.tm_mday,
             t.tm_hour,t.tm_min,t.tm_sec );
  return s;
}

////////////////////////////////////////////////////////////////////////////////

void DataBase::PrintTime(const tm &t,ostream &o)
{
  o << TimeStr(t);
}

////////////////////////////////////////////////////////////////////////////////

void DataBase::ElementBase::Print(ostream &o) const
{
  o << "Element " << GetName();

  if( GetTimeStart()!=0 )
  {
    o << "   start=";
    PrintTime(*GetTimeStart(),o);
  }

  if( GetTimeFinish()!=0 )
  {
    o << "   finish=";
    PrintTime(*GetTimeFinish(),o);
  }
}

////////////////////////////////////////////////////////////////////////////////

void DataBase::ElementBase::SetTimeStart(const tm &t)
{
  if( GetTimeFinish()!=NULL )
    if( t > *GetTimeFinish() )
      throw Exception("DataBase::ElementBase::SetTimeStart(): Attempt to set start time>finish time!\n"
                      "                                       new start=%s  finish=%s",TimeStr(t).c_str(),TimeStr(*GetTimeFinish()).c_str());
  if(time_start==NULL)
    time_start = new tm(t);
  else
    *time_start = t;
}

////////////////////////////////////////////////////////////////////////////////

string DataBase::ElementBase::CreateFileName(void) const
{
  char s[500];

  snprintf(s,500,"%s",name.c_str());

  if( NULL!=time_start )
  {
    const tm &t = *time_start;
    snprintf(s+strlen(s), 500-strlen(s),
             "~~start-%4.4d-%2.2d-%2.2d-%2.2d:%2.2d:%2.2d",
             t.tm_year+1900,t.tm_mon+1,t.tm_mday,
             t.tm_hour,t.tm_min,t.tm_sec);
  }

  if( NULL!=time_finish )
  {
    const tm &t = *time_finish;
    snprintf(s+strlen(s), 500-strlen(s),
             "~~finish-%4.4d-%2.2d-%2.2d-%2.2d:%2.2d:%2.2d",
             t.tm_year+1900,t.tm_mon+1,t.tm_mday,
             t.tm_hour,t.tm_min,t.tm_sec);
  }

  return s;
}

////////////////////////////////////////////////////////////////////////////////

void DataBase::FindAllVersions(const string &name,vector<DataBase::ElementBase> &versions) const
{
  struct dirent **namelist;
  int n = scandir(path.c_str(),&namelist,0,0);
  if (n < 0)
    throw Exception("DataBase::FindAllVersions():  can not scan directory %s",path.c_str());
  else
    while(n--)
    {
      try
      {
        ElementBase e( ElementBase::DecodeFullName(namelist[n]->d_name) );
        if( name==e.GetName() )
          versions.push_back(e);
      }
      catch( std::exception &e )
      {
        cerr << e.what() << "\n";
      }
      catch( const char *s )
      {
        cerr << s << "\n";
      }
      catch( ... )
      {
        cerr << "DataBase::FindAllVersions(): Unknown exception.\n";
      }
    }
}

////////////////////////////////////////////////////////////////////////////////

DataBase::ElementBase DataBase::ElementBase::DecodeFullName(const string &full_name)
{
  ElementBase e("");

  size_t n = full_name.find('~');
  e.name = full_name.substr(0,n);

  if( string::npos!=n )
  {
    n++;
    if( n!=full_name.find('~',n) )
      throw Exception("DataBase::ElementBase::DecodeFullName() error 1:  bad name %s",full_name.c_str());
    n++;

    //   ********** LOOKING FOR START TIME AND FINISH TIME FIELDS **********

    //            ***** I) START TIME *****

    // There has been a series of contradictory changes to the source code over
    // the successive versions in the range [1.8,1.11], in part due to a
    // change in the inner workings of "string::compare". Let's clarify the
    // situation:
    //  - We are expecting file names w/ the following syntax:
    //   "TBddCdcc~~start-<start time>~~finish-<finish time>"
    // (as least this was the initial intent of the author of the code, as
    // one can guess from earliest versions of fileBD, cf.
    // "/afs/cern.ch/compass/scratch/d01/newDB", which dates back from 2002.
    //  - Therefore, for the 1st time field, we check for keyword "start-" w/
    //   "full_name.compare(n,6,"start-")"

    //#define DataBase_LEARN_ABOUT_compare
#ifdef DataBase_LEARN_ABOUT_compare
    printf("\n--------------------%s:\n\n",e.name.c_str());
    string example("TBddCdcc~~start-<start time>~~finish-<finish time>");
    int ok = example.compare(n,6,"start-");
    printf("\"%s\".compare(n,6,\"start-\")=%d\n",example.c_str(),ok);
    ok = full_name.compare(n,6,"start-");
    printf("\"%s\".compare(n,6,\"start-\")=%d\n",full_name.c_str(),ok);
    ok = full_name.compare(n,6,"anything-");
    printf("\"%s\".compare(n,6,\"anything\")=%d\n",full_name.c_str(),ok);
#endif

    //    => "compare" returns 0 if and only if "full_name" fits.
    //  - Other contributors may have some different syntax in mind, which do
    //   NOT fit. Which leads them in turn to require a non zero return value
    //   from "compare". Which is WRONG!
    // Anyway Damien, in version v1.10, has a more elegant wording:
    if (full_name.substr(n,6)=="start-") {
      n += 6; size_t i = full_name.find('~',n);
      e.SetTimeStart(tm());
      DecodeTime(full_name.substr(n,i==string::npos?string::npos:i-n),
		 *e.GetTimeStart());
      if (i!=string::npos) {
        if( string::npos==full_name.find('~',i+1) )
          throw Exception("DataBase::ElementBase::DecodeFullName() error 2:  bad name %s",full_name.c_str());
        n = i+2;
      }
      else
        n = string::npos;
    }

    //            ***** II) FINISH TIME *****

    if (n!=string::npos && full_name.substr(n,7)== "finish-") { // Again reinstating Damien's wording
      n+=7; e.SetTimeFinish(tm());
      DecodeTime( full_name.substr(n,string::npos), *e.GetTimeFinish() );
    }
  }

  return e;
}

////////////////////////////////////////////////////////////////////////////////

void DataBase::ElementBase::DecodeTime(const string &s,tm &t)
{
  char c;
  sscanf(s.c_str(),"%d%c%d%c%d%c%d%c%d%c%d",
          &t.tm_year,   &c,
          &t.tm_mon,    &c,
          &t.tm_mday,   &c,
          &t.tm_hour,   &c,
          &t.tm_min,    &c,
          &t.tm_sec        );
  t.tm_year -= 1900;
  t.tm_mon  -= 1;
}

////////////////////////////////////////////////////////////////////////////////


string DataBase::CalibFilePath (const string &container, const tm &point)
{
  if( !(open_mode&READ) )
    throw Exception("DataBase::CalibFilePath(): data base \"%s\" is not opened for reading",
                     path.c_str());

  vector<ElementBase> vers;
  FindAllVersions(container,vers);
  for( size_t i=0; i<vers.size(); i++ )
  {
    ElementBase &e = vers[i];
    if( (e.GetTimeStart ()==NULL || *e.GetTimeStart ()<=point) &&
        (e.GetTimeFinish()==NULL || *e.GetTimeFinish()>=point) )
    {
      const string file_name = path+'/'+ElementBase(container,e.GetTimeStart(),e.GetTimeFinish()).CreateFileName();
      return file_name;
    }
  }

  throw Exception("DataBase::CalibFilePath():  can not find container \"%s\" in given time point.",
                   container.c_str());
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

} // namespace Reco

////////////////////////////////////////////////////////////////////////////////
