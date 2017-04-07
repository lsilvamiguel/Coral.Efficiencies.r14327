/*!
   \file    CsHistograms.cc
   \brief   Compass histograms package.
   \version $Revision: 1.24 $
   \author  Alexander Zvyagin
   \date    $Date: 2010/06/18 10:44:21 $
*/

#include <cassert>
#include <cctype>
#include <string.h>

#include "coral_config.h"

#include "TROOT.h"
#include "TFile.h"

#include "CsHistograms.h"
#include "CsOpt.h"
#include "CsErrLog.h"
#include "CsRegistry.h"
#include "CsHbookProto.h"
#include <string>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

CsHistograms   *CsHistograms::instance                  = NULL;
string          CsHistograms::home_path;
void           *CsHistograms::home                      = NULL;
bool            CsHistograms::is_activated              = false;
string          CsHistograms::top_directory;
CsHistograms::Implementation CsHistograms::implementation = CsHistograms::DUMMY;
unsigned int    CsHistograms::memory_consumed           = 0;
static string   current_path                            = "/";

////////////////////////////////////////////////////////////////////////////////

CsHistograms *CsHistograms::Instance(void)
{
  if( instance==NULL )
  {
    new CsHistograms();
  }
  return instance;
}

////////////////////////////////////////////////////////////////////////////////

void CsHistograms::Init(void)
{
  if( IsActivated() )
    CsErrLog::Instance()->mes(elFatal,"CsHistograms::Init: histograms package was initialised already.");
  else
  {
    is_activated = true;

    if(1)
    {
      // Put CsHistograms to class list for which end() function will be called at the end.
      CsRegistry reg;
      reg.EOJRegistration( Instance() );
    }
    
    // Initialization.
    string name;
    if( !CsOpt::Instance()->getOpt("","histograms package",name) )
      CsErrLog::Instance()->mes( elInfo,"CsHistograms::Init: histograms package is desirable");
    else
    {
      // Check histograms package name.
      if     ( name=="DUMMY" )
      {
        implementation = DUMMY;
      }
      else if( name=="ROOT"  )
      {
        implementation = ROOT;
          if( gROOT==NULL || !gROOT->Initialized() )
          {
            CsErrLog::Instance()->mes( elInfo,"CsHistograms::Init: ROOT package will be used.");
            static TROOT ROOT("","");
          }
      }
      else if( name=="HBOOK" )
      {
        implementation = HBOOK;
        #if USE_HBOOK
          // Nothing here. It is assumed that HBOOK was initialized early.
        #else
          CsErrLog::Instance()->mes( elInfo ,"You can not use option \"histograms package HBOOK\" without HBOOK.");
          CsErrLog::Instance()->mes( elFatal,"CsHistograms::Init: HBOOK package was NOT compiled!");
        #endif
      }
      else
        CsErrLog::Instance()->mes( elFatal, string("CsHistograms::Init: Unknown package name: ")+name );

      // Read histograms home path.
      home_path.erase();
      if( CsOpt::Instance()->getOpt("","histograms home",home_path) && GetImplementation()!=DUMMY )
      {
        CsErrLog::Instance()->mes( elInfo,string("CsHistograms::Init: histograms will be stored in \"")+
                               home_path + "\"" );
        
        // Initialize histograms home.
        switch( CsHistograms::GetImplementation() )
        {
          case CsHistograms::ROOT:
            home = TFile::Open(home_path.c_str(),"RECREATE","histograms",2);
            if( home==NULL || !reinterpret_cast<TFile*>(home)->IsOpen() )
              CsErrLog::Instance()->msg( elFatal, __FILE__,__LINE__,"CsHistograms::Init  ROOT package: can not open \"%s\"",home_path.c_str());
            break;

          #if USE_HBOOK
          case CsHistograms::HBOOK:
          {
            int recl=4096,error;
            hropen(735,"HBOOK",home_path.c_str(),"NQ",recl,error);
            if( error )
              CsErrLog::Instance()->msg( elFatal, __FILE__,__LINE__,"CsHistograms::Init  HBOOK package: can not open \"%s\"",home_path.c_str());
            break;
          }
          #endif

          case CsHistograms::DUMMY:
            break;

          default: 
            CsErrLog::Instance()->mes( elFatal,"CsHistograms::Init:  Interanl error.");
        }
      }
      else
        CsErrLog::Instance()->mes( elInfo,"CsHistograms::Init: histograms will not be saved!!!!!!!!!!");

      // Set top_directory (synonym is "/").
      switch( CsHistograms::GetImplementation() )
      {
        case CsHistograms::ROOT:
          top_directory = gDirectory->GetName();
          break;

        #if USE_HBOOK
        case CsHistograms::HBOOK:
        {
          char dir[]="                                                                               ";
          hcdir(dir,"R");  // Read current directory directory.
          top_directory = dir;
          break;
        }
        #endif

        dummy3:
        case CsHistograms::DUMMY:
          top_directory = "/";
          break;

        default: 
          CsErrLog::Instance()->mes( elFatal,"CsHistograms::Init:  Interanl error.");
      }

      if( top_directory.length()<=0 )
        CsErrLog::Instance()->mes( elFatal,"CsHistograms::Init:  Bad top directory for histograms.");
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

bool CsHistograms::end(void)
{
  if( is_activated )
  {
    switch( CsHistograms::GetImplementation() )
    {
      case CsHistograms::ROOT:
        if( home!=NULL )
        {
          reinterpret_cast<TFile*>(home)->Write("",TObject::kOverwrite);
          reinterpret_cast<TFile*>(home)->Close();
          delete reinterpret_cast<TFile*>(home);
        }
        break;

      #if USE_HBOOK
      case CsHistograms::HBOOK:
      {
//        hcdir("//PAWC"," ");  // Important!
        SetCurrentPath("/");
        int cycle;
        hrout(0,cycle,"T");
        hrend("HBOOK");
        break;
      }
      #endif

      case CsHistograms::DUMMY:
        break;

      default: 
        CsErrLog::Instance()->mes( elFatal,"CsHistograms::end:  Interanl error.");
    }

    home = NULL;
    is_activated = false;
  }
  else
    CsErrLog::Instance()->mes( elWarning,"CsHistograms::end():  This method was called already.");

  return true;
}

////////////////////////////////////////////////////////////////////////////////

string CsHistograms::GetCurrentPath(bool fast)
{
  if(fast)
    return current_path;
  else
  {
    string cur_dir;
    switch( CsHistograms::GetImplementation() )
    {
      case CsHistograms::ROOT:
        cur_dir = gDirectory->GetPath();
        break;

      #if USE_HBOOK
      case CsHistograms::HBOOK:
      {
        char dir[]="                                                                            ";
        hcdir(dir,"R");  // Read current directory directory.
        if( strlen(dir)<3 || dir[0]!='/' || dir[1]!='/' )
          CsErrLog::Instance()->mes( elFatal,"CsHistograms::GetCurrentPath  HBOOK: parse error (test that HBOOK was initialized).");

        cur_dir = dir;
        break;
      }
      #endif

      case CsHistograms::DUMMY:
        cur_dir = current_path;
        break;

      default: 
        CsErrLog::Instance()->mes( elFatal,"CsHistograms::GetCurrentPath  Interanl error.");
    }

    if( cur_dir.compare(0,top_directory.length(),top_directory) ) {
      CsErrLog::Instance()->msg( elFatal,__FILE__,__LINE__,
          "CsHistograms::GetCurrentPath  Can not get current directory:\n"
          "package top directory     = \"%s\"\n"
          "package current directory = \"%s\"\n",
          top_directory.c_str(),cur_dir.c_str());
    }

    cur_dir.replace(0,top_directory.length(),"/");

    // Remove possible double '/' in the path (like in "//path").
    if( cur_dir.length()>=2 )
      if( cur_dir[1]=='/' )
        cur_dir.erase(0,1);

    return cur_dir;
  }
}

////////////////////////////////////////////////////////////////////////////////

namespace { // private namespace of this file

void SetCurrentPathHBOOK(const string &path)
{
//  cout<<"SetCurrentPathHBOOK="<<path<<"\n";

  if( path.empty() )
    return;

  size_t n = (path[0]=='/') ? 1 : 0; // Skip leading '/'

  for( string::size_type i=n; i<path.length(); i++ )
    if( path[i]=='/' || i==path.length()-1 )
    {
      string path_to_go(path,n,i-(path[i]=='/')-n+1);
      n=i+1;

      // Find directory.
      int id=0;
      char type[88]="D", title[88];
      while( hlnext(id,type,title," "), id!=0 )
        if( *type=='D' && path_to_go.length()==strlen(title) )
        {
          // Compare (ignoring case) path_to_go and title.
          bool cmp_result=true;
          for( size_t i=0; i<path_to_go.length(); i++ )
            if( toupper(path_to_go[i])!=toupper(title[i]) )
            {
              cmp_result=false;
              break;
            }

          if( cmp_result )
          {
            hcdir(path_to_go.c_str()            ," ");
            break;
          }
        }

      if( id==0 )    // subdirectory was not found.
        hmdir(path_to_go.c_str(),"S");
    }
}

}

void CsHistograms::SetCurrentPath(const string &path,bool fast)
{
//  cout<<"SetCurrentPath="<<path<<"\n";

  #if 0 && !defined(NDEBUG)
  if( GetCurrentPath(true)!=GetCurrentPath(false) )
  {
    cerr << "GetCurrentPath(true )=\"" << GetCurrentPath(true ) << "\"\n"
         << "GetCurrentPath(false)=\"" << GetCurrentPath(false) << "\"\n";
    exit(1);
  }
  #endif

  if( path==GetCurrentPath(fast) || path.empty() || path=="." )
    return;

  current_path = path;
  
  if( current_path=="/" )
    

  if( current_path.length()>1 && current_path[current_path.length()-1]=='/' )
    current_path.erase(current_path.length()-1); // erase trailing '/'
  
  if( current_path[0]!='/' )
    CsErrLog::Instance()->mes( elFatal,string("CsHistograms::SetCurrentPath  Sorry, the limitation: use absolute paths. dir = ")+current_path );

  if( CsHistograms::GetImplementation()==CsHistograms::DUMMY )
    return;

  size_t n=0;
  if( current_path[0]=='/' )
  {
    n++;
    switch( CsHistograms::GetImplementation() )
    {
      case CsHistograms::ROOT:
      {
        //if( !gDirectory->cd(top_directory.c_str()) )
        TObject *o = gROOT->GetListOfFiles()->FindObject(top_directory.c_str());
        TDirectory *d=NULL;
        if( o!=NULL )
            d = dynamic_cast<TDirectory*>(o);
        if( d==NULL || !d->cd() )
          CsErrLog::Instance()->mes( elFatal,"CsHistograms::SetCurrentPath  ROOT package: can not change directory to \"/\"" );
        break;
      }

      #if USE_HBOOK
      case CsHistograms::HBOOK:
      {
        hcdir(string("//PAWC").c_str()," ");
        hcdir(top_directory.c_str()," ");
        break;
      }
      #endif

      case CsHistograms::DUMMY:
        current_path = "/";
        break;
      
      default:
        CsErrLog::Instance()->mes( elFatal,"CsHistograms::SetCurrentPath  Interanl error.");
    }
  }
  
  if( current_path!="/" )
    switch( CsHistograms::GetImplementation() )
    {
      case CsHistograms::ROOT:
      {
        for( string::size_type i=n; i<current_path.length(); i++ )
          if( current_path[i]=='/' || i==current_path.length()-1 )
            if( n!=i )
            {
              string path_to_go(current_path,n,i-(current_path[i]=='/')-n+1);
              //string path_to_go(current_path,n,i-(current_path[i]=='/'));
              n=i+1;

              TObject *obj = gDirectory->Get(path_to_go.c_str());
              if( obj!=NULL )
              {
                // OK, the object exists.

                // Check that this object is a directory.
                if (!obj->InheritsFrom(TDirectory::Class()))
                  CsErrLog::Instance()->msg( elFatal,__FILE__,__LINE__,
                    "CsHistograms::SetCurrentPath  ROOT package: cur.dir=\"%s\", this is not a directory: \"%s\"",
                    GetCurrentPath(fast).c_str(),path_to_go.c_str() );
                else
                  if( !static_cast<TDirectory*>(obj)->cd() )
                    CsErrLog::Instance()->msg( elFatal,__FILE__,__LINE__,
                      "CsHistograms::SetCurrentPath  ROOT package: can not go from \"%s\", to \"%s\"",
                      GetCurrentPath(fast).c_str(),path_to_go.c_str() );
              }
              else
              {
                if( NULL==gDirectory->mkdir(path_to_go.c_str()) )
                  CsErrLog::Instance()->mes( elFatal,
                                             string("CsHistograms::SetCurrentPath  ROOT package: can not create directory: \"")+path_to_go+"\"" );
                if( !gDirectory->cd(path_to_go.c_str()) )
                  CsErrLog::Instance()->mes( elFatal,
                                             string("CsHistograms::SetCurrentPath  ROOT package: can not change directory to: \"")+path_to_go+"\"" );
              }
            }
      }
      break;

      #if USE_HBOOK
      case CsHistograms::HBOOK:
      {
        hcdir("//PAWC"," ");
        SetCurrentPathHBOOK(current_path);

        hcdir(top_directory.c_str()," ");
        SetCurrentPathHBOOK(current_path);

        break;
      }
      #endif

      default:
        CsErrLog::Instance()->mes( elFatal,"CsHistograms::SetCurrentPath  Interanl error.");
    }
}

////////////////////////////////////////////////////////////////////////////////

CsHistograms::~CsHistograms(void)
{
  if( IsActivated() )
  {
    CsErrLog::Instance()->mes( elWarning,
      "CsHistograms::~CsHistograms:  end() method was not called!  Force call of CsHistograms::end()");
    end();
  }
  instance = NULL;
}

////////////////////////////////////////////////////////////////////////////////

CsHistograms::CsHistograms(void)
{
  if (instance!=NULL)
  {
    CsErrLog::Instance()->mes(elFatal,"CsHistograms::CsHistograms: second instance of histogram package detected.");
  }
  if( !IsActivated() ) {
    Init();
  }

  instance = this;
}

////////////////////////////////////////////////////////////////////////////////
