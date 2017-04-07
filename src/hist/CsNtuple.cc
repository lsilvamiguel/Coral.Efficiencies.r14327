/*!
   \file    CsNtuple.cc
   \brief   Compass Ntuple Class.
   \version $Revision: 1.11 $
   \author  Alexander Zvyagin
   \date    $Date: 2010/06/18 10:44:21 $
*/

#include "coral_config.h"

#include "CsNtuple.h"
#include "CsHistograms.h"
#include "CsErrLog.h"

#include "TNtuple.h"

#if USE_HBOOK
#include "CsHBOOK.h"
#endif

#define Ntuple_dummy reinterpret_cast<CsHomeNamed*>(nt)
#define Ntuple_hbook reinterpret_cast<CsHbookNtuple*>(nt)

using namespace std;

////////////////////////////////////////////////////////////////////////////////

CsNtuple::~CsNtuple(void)
{
  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        break;

    case CsHistograms::HBOOK:
      #if USE_HBOOK
        delete [] Ntuple_hbook->buffer;
        delete Ntuple_hbook;
        break;
      #endif
      goto dummy;

    dummy:
    case CsHistograms::DUMMY:
      delete Ntuple_dummy;
      break;

    default: 
      CsErrLog::Instance()->mes( elFatal,"CsNtuple::~CsNtuple:  Internal error.");
  }
}

////////////////////////////////////////////////////////////////////////////////

CsNtuple::CsNtuple(const string &name,const string &title,const string &format)
{
  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        nt = new TNtuple(name.c_str(),title.c_str(),format.c_str());
        return;

    case CsHistograms::HBOOK:
    {
      #if USE_HBOOK
        nt = new CsHbookNtuple(name,title);
        string fmt(format);
        for( size_t i=0; i<fmt.length(); i++ )
          if( fmt[i]==':' || i+1==fmt.length() )
          {
            if( fmt[i]==':' )
              fmt[i]=',';
            fmt.insert(i+1-(fmt[i]==','),":R");
            i+=2;
            Ntuple_hbook->Nvar ++;
          }
        // cout.form("format OLD: \"%s\n",format.c_str());
        // cout.form("format NEW: \"%s\n",fmt   .c_str());
        
        Ntuple_hbook->buffer = new float32[Ntuple_hbook->Nvar];
        hbnt  (Ntuple_hbook->ID(),Ntuple_hbook->GetHbookTitle().c_str(),"");
        hbname(Ntuple_hbook->ID(),"block",Ntuple_hbook->buffer,fmt.c_str());

        break;
      #endif
      goto dummy;
    }

    dummy:
    case CsHistograms::DUMMY:
      nt = new CsHomeNamed(name,title);
      return;

    default:
      CsErrLog::Instance()->mes( elFatal,"CsNtuple::CsNtuple:  Internal error.");
  }
}

////////////////////////////////////////////////////////////////////////////////

string CsNtuple::GetName(void) const
{
  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
      //        return Ntuple_hbook->GetName();
      return string();  // FIXME

    case CsHistograms::HBOOK:
    {
      #if USE_HBOOK
        return Ntuple_hbook->GetName();
      #endif
      goto dummy;
    }

    dummy:
    case CsHistograms::DUMMY:
      return string();

    default:
      CsErrLog::Instance()->mes( elFatal,"CsNtuple::GetName:  Internal error.");
  }
  
  return string(); // This makes compiler happy.
}

////////////////////////////////////////////////////////////////////////////////

string CsNtuple::GetTitle(void) const
{
  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        return reinterpret_cast<TNtuple*>(nt)->GetTitle();

    case CsHistograms::HBOOK:
    {
      #if USE_HBOOK
        return Ntuple_hbook->GetTitle();
      #endif
      goto dummy;
    }

    dummy:
    case CsHistograms::DUMMY:
      return string();

    default:
      CsErrLog::Instance()->mes( elFatal,"CsNtuple::GetTitle:  Internal error.");
  }
  
  return string(); // This makes compiler happy.
}

////////////////////////////////////////////////////////////////////////////////

void CsNtuple::SetTitle(const string &title)
{
  CsErrLog::Instance()->mes(elFatal,"CsNtuple::SetTitle:: not implemented yet");
}

////////////////////////////////////////////////////////////////////////////////

void CsNtuple::SetName(const string &title)
{
  CsErrLog::Instance()->mes(elFatal,"CsNtuple::SetName:: not implemented yet");
}

////////////////////////////////////////////////////////////////////////////////

CsNtuple &CsNtuple::operator = (const CsNtuple &ntuple)
{
  CsErrLog::Instance()->mes( elFatal,"CsNtuple::operator=  There is no code yet.");
// 
//   switch( CsHistograms::GetImplementation() )
//   {
//     case CsHistograms::ROOT:
//         CsErrLog::Instance()->mes( elFatal,"CsNtuple::operator=  There is no code yet.");
//         break;
// 
//     case CsHistograms::DUMMY:
//       Ntuple_dummy->operator=(*ntuple.h);
//       break;
// 
//     default:
//       CsErrLog::Instance()->mes( elFatal,"CsNtuple::operator=  Internal error.");
//   }
  return *this;
}

////////////////////////////////////////////////////////////////////////////////

int32 CsNtuple::Fill(const float32 *x)
{
  if( x==NULL )
    CsErrLog::Instance()->mes(elFatal,"CsNtuple::Fill: x==NULL");

  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        return reinterpret_cast<TNtuple*>(nt)->Fill(const_cast<float32*>(x));

    case CsHistograms::HBOOK:
      #if USE_HBOOK
        memcpy(Ntuple_hbook->buffer,x,Ntuple_hbook->Nvar*4);
        hfnt(Ntuple_hbook->ID());
        break;
      #endif
      goto dummy;

    dummy:
    case CsHistograms::DUMMY:
      return 0;

    default:
      CsErrLog::Instance()->mes( elFatal,"CsNtuple::Fill:  Internal error.");
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int32 CsNtuple::Fill(float32 x1,float32 x2,float32 x3,float32 x4,
                     float32 x5,float32 x6,float32 x7,float32 x8,float32 x9)
{
  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        return reinterpret_cast<TNtuple*>(nt)->Fill(x1,x2,x3,x4,x5,x6,x7,x8,x9);

    case CsHistograms::HBOOK:
      #if USE_HBOOK
        switch( Ntuple_hbook->Nvar )
        {
          case  9:      Ntuple_hbook->buffer[8] = x9;
          case  8:      Ntuple_hbook->buffer[7] = x8;
          case  7:      Ntuple_hbook->buffer[6] = x7;
          case  6:      Ntuple_hbook->buffer[5] = x6;
          case  5:      Ntuple_hbook->buffer[4] = x5;
          case  4:      Ntuple_hbook->buffer[3] = x4;
          case  3:      Ntuple_hbook->buffer[2] = x3;
          case  2:      Ntuple_hbook->buffer[1] = x2;
          case  1:      Ntuple_hbook->buffer[0] = x1;
          case  0:      break;
          default:      CsErrLog::Instance()->msg( elFatal,__FILE__,__LINE__,
                          "CsNtuple::Fill:  Bad variables amount: %d",Ntuple_hbook->Nvar );
      }
      hfnt(Ntuple_hbook->ID());
      #endif
      goto dummy;

    dummy:
    case CsHistograms::DUMMY:
      return 0;

    default:
      CsErrLog::Instance()->mes( elFatal,"CsNtuple::Fill:  Internal error." );
  }
  return 0;
}


////////////////////////////////////////////////////////////////////////////////

string CsNtuple::GetPath(void) const
{
  CsErrLog::Instance()->mes( elFatal,"CsNtuple::GetPath: Not implemented yet.");
  return string();
}

////////////////////////////////////////////////////////////////////////////////

void CsNtuple::SetPath(const string &path)
{

//   if( CheckPath(path) )
//     CsErrLog::Instance()->mes( elFatal,string("CsHistograms::SetCurrentPath  Your path is incorrect\n")+
//                                        "Path=\"" + path + "\"" );
        CsErrLog::Instance()->mes( elFatal,"CsNtuple::SetHomePath:  Not implemented yet.  (0)");

  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        CsErrLog::Instance()->mes( elFatal,"CsNtuple::SetPath: ROOT: Not implemented yet.");
        break;

    case CsHistograms::HBOOK:
      #if USE_HBOOK
        CsErrLog::Instance()->mes( elFatal,"CsNtuple::SetPath: ROOT:  Not implemented yet.");
        break;
      #endif
      goto dummy;

    dummy:
    case CsHistograms::DUMMY:
      Ntuple_dummy->SetPath(path);
      break;

    default: 
      CsErrLog::Instance()->mes( elFatal,"CsNtuple::GetHomePath:  Internal error.");
  }
}

////////////////////////////////////////////////////////////////////////////////
