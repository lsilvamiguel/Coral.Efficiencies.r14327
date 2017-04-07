/*
   $Source: RunDataBase.cc,v $ 
   $Date:  $ 
   $Revision:  $ 
   -------------------------------------------------------------------------

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

     Modified by V.Kolosov ( Vladimir.Kolosov@cern.ch, Kolosov@mx.ihep.su )
     from time to run DataBase   Mon Feb 11 2002 
*/

#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>

#include "RunDataBase.h"

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

namespace MN {

using namespace std;
using namespace Reco;

////////////////////////////////////////////////////////////////////////////////

RunDataBase::RunDataBase(const string &path_, int open_mode_) :
  path          (path_),
  open_mode     (open_mode_),
  creation_flag (false),
  structure_id_(PLANE)
{
  if( open_mode&CREATE && !(open_mode&READ) )
    system( string("rm -rf " + path).c_str() );

  if( open_mode&(CREATE|WRITE) )
  {
    mode_t mode=0777;
    int error = mkdir(path.c_str(),mode);
    if( error==0 )
      creation_flag = true;
#warning TODO  test that directory exists or was created
//     if( error )
//       throw Exception("RunDataBase::RunDataBase():  can not create directory \"%s\"   error code %d",path.c_str(),error);
  }
}

////////////////////////////////////////////////////////////////////////////////

string RunDataBase::RunStr(const size_t &t)
{
  char s[111];
  sprintf(s,"%zu",t);
  return s;
}

////////////////////////////////////////////////////////////////////////////////

void RunDataBase::PrintRun(const size_t &t,ostream &o)
{
  o << RunStr(t);
}

/////////////////////////////////////////////////////////////////////////////////

void RunDataBase::CreateStructure(const std::string &container, const std::string &tag, const size_t &start, const size_t &end) const
{
  bool debug = true; 
  if( structure_id_ == PLANE )  
  {
    return;
  }   
  else if( structure_id_ == DETSUBDIRS )  
  {
    string dir = string(path + "/" + tag);
    mode_t mode=0777;
    int error = mkdir(dir.c_str(),mode);
    dir = string(path + "/" + tag + "/" + container);
    error = mkdir(dir.c_str(),mode);
    return;
  }   
  else if( structure_id_ == DETCALIBSUBDIRS )  
  {
    string dir = string(path + "/" + container);
    mode_t mode=0777;
    int error = mkdir(dir.c_str(),mode);
    dir = string(path + "/" + tag + "/" + container+ "/" + tag);
    error = mkdir(dir.c_str(),mode);
    return;
  }
  return;   
}

////////////////////////////////////////////////////////////////////////////////

string RunDataBase::CreateTimeStamp( const size_t &start, const size_t &end) const
{
  if( structure_id_ == PLANE )  
  {
    char s[100];
    if( 0!= start )
    {
      const size_t &t = start;
      sprintf(s,"~~start-Run_%zu",t);
    }

    if( 0!= end )
    {
      const size_t &t = end;
      sprintf(s+strlen(s),"~~finish-Run_%zu",t);
    }
    return string(s);
  }   
  else if( structure_id_ == DETSUBDIRS || structure_id_ == DETCALIBSUBDIRS )  
  {
    char s[100];
    if( 0!= start )
    {
      const size_t &t = start;
      sprintf(s,"%zu-",t);
    }

    if( 0!= end )
    {
      const size_t &t = end;
      sprintf(s+strlen(s),"%zu",t);
    }
    return string(s);
  }   
  else   
  {
    char s[100];
    if( 0!= start )
    {
      const size_t &t = start;
      sprintf(s,"~~start-Run_%zu",t);
    }

    if( 0!= end )
    {
      const size_t &t = end;
      sprintf(s+strlen(s),"~~finish-Run_%zu",t);
    }
    return string(s);
  }   
}

////////////////////////////////////////////////////////////////////////////////

string RunDataBase::CreateFileName(const std::string &container, const std::string &tag, const size_t &start, const size_t &end) const
{
  if( structure_id_ == PLANE )  
  {
    std::string file_name = path + '/' + container + tag + CreateTimeStamp(start,end);
    return file_name;
  }   
  else if( structure_id_ == DETSUBDIRS )  
  {
    std::string file_name = path + '/' + tag + '/' + container + '/' + CreateTimeStamp(start,end);
    return file_name;
  }   
  else if( structure_id_ == DETCALIBSUBDIRS )  
  {
    std::string file_name = path + '/' + container + '/' + tag + '/' + CreateTimeStamp(start,end);
    return file_name;
  }
  else
  {
    std::string file_name = path + '/' + container + CreateTimeStamp(start,end);
    return file_name;
  }
}

////////////////////////////////////////////////////////////////////////////////

string RunDataBase::CreatePathName(const std::string &container, const std::string &tag) const
{
  if( structure_id_ == PLANE )  
  {
    std::string path_name = path;
    return path_name;
  }   
  else if( structure_id_ == DETSUBDIRS )  
  {
    std::string path_name = path + '/' + tag + '/' + container;
    return path_name;
  }   
  else if( structure_id_ == DETCALIBSUBDIRS )  
  {
    std::string path_name = path + '/' + container + '/' + tag;
    return path_name;
  }
  else
  {
    std::string path_name = path;
    return path_name;
  }
}

////////////////////////////////////////////////////////////////////////////////

string RunDataBase::CreateFileNameLast(const std::string &container, const std::string &tag, const size_t &start, const size_t &end) const
{
  if( structure_id_ == PLANE )  
  {
    std::string file_name = path+'/' + container + tag + "~~last";
    return file_name;
  }   
  else if( structure_id_ == DETSUBDIRS )  
  {
    std::string file_name = path + '/' + tag + '/' + container + '/' + "last";
    return file_name;
  }   
  else if( structure_id_ == DETCALIBSUBDIRS )  
  {
    std::string file_name = path + '/' + container + '/' + tag + '/' + "last";
    return file_name;
  }
  else
  {
    std::string file_name = path+'/' + container + tag + "~~last";
    return file_name;
  }
}

///////////////////////////////////////////////////////////////////////////////

bool RunDataBase::FileNameSyntaxOK(const string &file_name, const string &container, const string &tag, std::pair <size_t, size_t> &time_stamps ) const
{
//   bool debug = true;
  bool debug = false;
  if( debug ) cout << " FileNameSyntaxOK ?? for " << file_name << " Mode " << structure_id_ << endl;
  
  bool ok = false;
  size_t &run_start = time_stamps.first;
  size_t &run_finish = time_stamps.second;
  run_start = 0;
  run_finish =0;

  if( structure_id_ == PLANE )  
  {
    if (file_name.substr( 0,(container+tag).size() ) == container+tag )
    {
      if( debug ) cout << " FileNameSyntaxOK ?? file name " << file_name << " match " << container+tag << endl;
    }
    else
    {
//      if( debug ) cout << "FileNameSyntaxOK ?? file name " << file_name << " does not match " << container+tag << endl;
      return false;
    }


    size_t n = file_name.find('~');
    if( debug ) cout <<" (container+tag).size() " << (container+tag).size() << " compare n " << n << endl;
    if( n != (container+tag).size() ) return false;
    if( string::npos!=n )
    {
      n++;
      if( n==file_name.find('~',n) )
      {
        n++;
        if (file_name.substr(n,10)=="start-Run_") {
          n += 10; 
	  size_t i = file_name.find('~',n);
          if( !DecodeRun( file_name.substr(n,i==string::npos?string::npos:i-n), run_start ) ) return false;
          if (i!=string::npos) {
            if( string::npos==file_name.find('~',i+1) )
	    {
              cerr << "DataBase::ElementBase::DecodeFullName() error 2:  bad name " << file_name << endl;
	      return false;
	    }  
            n = i+2;
          }
          else
            n = string::npos;
        }

        if (n!=string::npos && file_name.substr(n,11)== "finish-Run_") 
	{ 
          n+=11;
          if( !DecodeRun( file_name.substr(n,string::npos), run_finish ) ) return false;
        }
      }
      if( run_finish < run_start ) return false;
      ok=true;
    }
    return ok;
  }   
  else if( structure_id_ == DETSUBDIRS )  
  {
    size_t n = file_name.find('-');
    if( n == string::npos ) return false;
    if( debug ) cout << " Check for - n = " << n << " npos " << string::npos << endl;
    size_t ln = strlen(file_name.c_str());
    
    if( debug ) cout << " Decode Run Start " << ln << " 0:" << n <<" "<< file_name.substr(0,n) << endl;
    if( !DecodeRun( file_name.substr(0,n), run_start ) ) return false;
    if( debug ) cout << " Run Start " << run_start << endl;
    n+=1;
    if( ln > n )
    {
      if( debug ) cout << " Decode Run Finish " << ln << " "<< n << ":" << string::npos <<" "<< file_name.substr(n,string::npos) << endl;
      if( !DecodeRun( file_name.substr(n,string::npos), run_finish ) ) return false;
      if( debug ) cout << " Run Finish " << run_finish << endl;
      ok = true;
    }
    else
    {
      run_finish = RUN_MAX;
      ok = true;
    }  
    if( run_finish < run_start ) return false;
    return ok;
  }   
  else if( structure_id_ == DETCALIBSUBDIRS )  
  {
    size_t n = file_name.find('-');
    if( n == string::npos ) return false;
    if( debug ) cout << " Check for - n = " << n << " npos " << string::npos << endl;
    size_t ln = strlen(file_name.c_str());

    if( debug ) cout << " Decode Run Start " << ln << " 0:" <<n<<" "<< file_name.substr(0,n) << endl;
    if( !DecodeRun( file_name.substr(0,n), run_start ) ) return false;
    n+=1;
    if( ln > n )
    {
      if( debug ) cout << " Decode Run Finish " << ln << " "<< n << ":"<< string::npos <<" "<< file_name.substr(n,string::npos) << endl;
      if( !DecodeRun( file_name.substr(n,string::npos), run_finish ) ) return false;
      ok = true;
    }
    else
    {
      run_finish = RUN_MAX;
      ok = true;
    }  
    if( run_finish < run_start ) return false;
    return ok;
  }
  else
  {
    return ok;
  }

}

///////////////////////////////////////////////////////////////////////////////

void RunDataBase::FindAllGoodVersions(const string &name, const string &tag, vector< std::string > &versions, const size_t &point) const
{
//   bool debug = true;
  bool debug = false;
  if( debug ) cout << " FindAllVersions for " << name <<" tag " << tag << endl;
  cout << " FindAllVersions for " << name <<" tag " << tag << endl;
  
  string path_dir = CreatePathName( name, tag);
  struct dirent **namelist;
  if( debug ) cout << " scandir " << path_dir << endl;
  int n = scandir(path_dir.c_str(),&namelist,0,0);
  if (n < 0)
  {
    cerr <<" RunDataBase::FindAllGoodVersions():  can not scan directory " << path_dir << endl;
    return;
  }  
  else
  {
    if( debug ) cout << " scandir " << path_dir << " OK n = " << n << endl;
    while(n--)
    {
      if( debug ) std::cout << " namelist[n]->d_name " << namelist[n]->d_name << std::endl;
      std::pair <size_t, size_t> time_stamps;
      if( ! FileNameSyntaxOK ( namelist[n]->d_name, name, tag, time_stamps ) ) continue;
      if( point  != 0 && ( (time_stamps.first > point) || (time_stamps.second < point ) ) ) continue;
      if( debug )
      {
         cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FileNameSyntaxOK run start " << time_stamps.first <<" point " << point << " run finish " << time_stamps.second  << endl;
      }
      versions.push_back( path_dir + "/" + namelist[n]->d_name );
      if( debug ) std::cout << " Version found!! Name " << namelist[n]->d_name << std::endl;
      std::cout << " Version found!! Name " << namelist[n]->d_name << std::endl;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

bool RunDataBase::DecodeRun(const string &s,size_t &t)
{
  bool debug = false;
  int n = sscanf(s.c_str(),"%zu",&t);
  if( debug ) cout << " DecodeRun n = " << n << endl;
  if( n != 1 ) return false;
  return true;
//   char c;
//   sscanf(s.c_str(),"%d%c%d%c%d%c%d%c%d%c%d",
//           &t.tm_year,   &c,
//           &t.tm_mon,    &c,
//           &t.tm_mday,   &c,
//           &t.tm_hour,   &c,
//           &t.tm_min,    &c,
//           &t.tm_sec        );
//   t.tm_year -= 1900;
//   t.tm_mon  -= 1;
}

////////////////////////////////////////////////////////////////////////////////
// 
// bool operator == (const tm &t1, const tm &t2)
// {
//   return
//     t1.tm_year  == t2.tm_year   &&
//     t1.tm_mon   == t2.tm_mon    &&
//     t1.tm_mday  == t2.tm_mday   &&
//     t1.tm_hour  == t2.tm_hour   &&
//     t1.tm_min   == t2.tm_min    &&
//     t1.tm_sec   == t2.tm_sec;
// }
// 
// ////////////////////////////////////////////////////////////////////////////////
// 
// bool operator > (const tm &t1, const tm &t2)
// {
//   if( t1.tm_year > t2.tm_year )  return true;
//   if( t1.tm_year < t2.tm_year )  return false;
// 
//   if( t1.tm_mon  > t2.tm_mon  )  return true;
//   if( t1.tm_mon  < t2.tm_mon  )  return false;
// 
//   if( t1.tm_mday > t2.tm_mday )  return true;
//   if( t1.tm_mday < t2.tm_mday )  return false;
// 
//   if( t1.tm_hour > t2.tm_hour )  return true;
//   if( t1.tm_hour < t2.tm_hour )  return false;
// 
//   if( t1.tm_min  > t2.tm_min  )  return true;
//   if( t1.tm_min  < t2.tm_min  )  return false;
// 
//   if( t1.tm_sec  > t2.tm_sec  )  return true;
//   if( t1.tm_sec  < t2.tm_sec  )  return false;
// 
//   return false;
// }
// 

} // namespace MN

////////////////////////////////////////////////////////////////////////////////
