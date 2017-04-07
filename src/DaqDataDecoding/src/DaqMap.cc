#include <sstream>
#include <cstdio>
#include <cstdlib>

#include "DaqMap.h"
#include "Exception.h"
#include "ObjectXML.h"

using namespace std;

namespace CS {

////////////////////////////////////////////////////////////////////////////////

DaqMap::DaqMap(const ObjectXML &o) :
  ObjectXML(o),
  version(0)
{
    Attribute *a;
    string s;
    if( NULL!=(a=GetAttribute("runs",s)) )
        SetRunsList(a->second);
    else
        SetAttribute("runs","0-999999999");

    // Read options
    if( NULL!=GetAttribute("option",options) || NULL!=GetAttribute("options",options) );

    GetAttribute("version",version);
}

////////////////////////////////////////////////////////////////////////////////

void DaqMap::Clear(void)
{
    version=0;
    SetName("");
    options="";
    GetAttributes().clear();
    runs_list.clear();
}

////////////////////////////////////////////////////////////////////////////////

void DaqMap::SetRunsList(const string &s)
{
  istringstream str(s.c_str());

  string run;  // Examples:   "1"   "2-5"
  while( str>>run )
  {
    int runF,runL;
    size_t p = run.find('-');
    if( p!=string::npos )
    {
      runF=atoi(run.substr(0  ,p           ).c_str());
      runL=atoi(run.substr(p+1,string::npos).c_str());
    }
    else
      runF=runL=atoi(run.c_str());

    if( runF>runL || runF<0 )
      throw Exception("Chip::Maps::SetRunsList():  RunFirst>RunLast or negative runs:  RunFirst=%d  RunLast=%d",
                       runF,runL);

    if( !runs_list.empty() && runs_list.back().second>=uint32(runF) )
      throw Exception("Chip::Maps::SetRunsList():  Overlapping runs: [%d,%d] [%d,%d]",
                       runs_list.back().first,runs_list.back().second,
                       runF,runL);
    
    runs_list.push_back( pair<uint32,uint32>(uint32(runF),uint32(runL)) );
  }
  
  // Set attributes
  SetAttribute("runs",s);
}

////////////////////////////////////////////////////////////////////////////////

bool DaqMap::TestRunsRange(uint32 run) const
{
  for( vector< pair<uint32,uint32> >::const_iterator it=runs_list.begin(); it!=runs_list.end(); it++ )
      if( run>=it->first && run<=it->second )
        return true;

  return false;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace CS
