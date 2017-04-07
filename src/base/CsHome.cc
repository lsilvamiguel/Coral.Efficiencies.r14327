#include "CsHome.h"
#include "CsErrLog.h"

using namespace std;

static string current_path("/");

////////////////////////////////////////////////////////////////////////////////

CsHome::CsHome(const string &name,const string &title,const string &path)
{
}

////////////////////////////////////////////////////////////////////////////////

string CsHome::GetCurrentPath(const string &package)
{
  if( package=="ROOT" )
  {
  }
  else if( package=="DUMMY" )
    return current_path;
  else 
  {
    CsErrLog::Instance()->mes( elFatal,string("CsHome::GetCurrentPath:  Unknown package name \"")+package+"\"");
  }
  return "";
}

////////////////////////////////////////////////////////////////////////////////

void CsHome::SetCurrentPath(const string &path,const string &package)
{
  if( package=="ROOT" )
  {
  }
  else if( package=="DUMMY" )
    current_path = path;
  else 
  {
    CsErrLog::Instance()->mes( elFatal,string("CsHome::SetCurrentPath:  Unknown package name \"")+package+"\"");
  }
}

////////////////////////////////////////////////////////////////////////////////

CsHome &CsHome::operator = (const CsHome &home)
{
  CsErrLog::Instance()->mes( elFatal,"CsHome::operator=  No body yet.\n");
  
//   if( &home!=this )
//     Rename(home.GetName(),home.GetTitle(),home.GetPath());
  return *this;
}

////////////////////////////////////////////////////////////////////////////////
