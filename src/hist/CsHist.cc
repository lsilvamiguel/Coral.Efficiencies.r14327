/*!
   \file    CsHist.cc
   \brief   Compass histograms classes implementation.
   \version $Revision: 1.23 $
   \author  Alexander Zvyagin
   \date    $Date: 2010/06/18 10:44:21 $
*/


#include <iostream>
#include <cassert>

#include "coral_config.h"
#include "CsErrLog.h"
#include "CsHistograms.h"
#include "CsHist.h"

#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TDirectory.h"

#if USE_HBOOK
#include "CsHbookProto.h"
#include "CsHBOOK.h"
extern int quest_[100];
#endif

using namespace std;


////////////////////////////////////////////////////////////////////////////////
// 
//                         class CsFitParam
// 
////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////

CsFitParam::CsFitParam(float64 v) :
  name         (""),
  value_start  (v),
  value_finish (v),
  error        (1),
  step         (1),
  min          (0),
  max          (0)
{}

////////////////////////////////////////////////////////////////////////////////

CsFitParam::CsFitParam(const string name_,float64 value_start_,float64 step_,float64 min_, float64 max_) :
  name         (name_),
  value_start  (value_start_),
  value_finish (value_start_),
  error        (step_),
  step         (step_),
  min          (min_),
  max          (max_)
{}

////////////////////////////////////////////////////////////////////////////////

void CsFitParam::Print(ostream &o) const
{
  o << "\"" << GetName().c_str() << "\" start=" << (double)GetStart()
    << " ==> finish=" << (double)GetFinish()
    << " error=" << (double)GetError()
    << " step="  <<(double)GetStep()
    << " min=" << (double)GetMax()
    << " max=" <<(double)GetMin() << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

static struct
{
  FitFunction fit_function;
  size_t Nx,Npar;
} fit_data;

static Double_t root_fit(Double_t *x,Double_t *p)
{
  assert( x!=NULL && p!=NULL );
  assert(fit_data.Nx  >0);
  assert(fit_data.Npar>0);
  const vector<Double_t> 
    xx(x,x+fit_data.Nx  ), 
    pp(p,p+fit_data.Npar);
  return fit_data.fit_function(xx,pp);
}

#if USE_HBOOK
//       FUNCTION UFIT(X)
// *    The dimension of DPAR  ||  MUST be 24!
// *                           VV
//       DOUBLE PRECISION DPAR(24),FITFUN 
//       COMMON/HCFITD/DPAR,FITFUN
//       FITFUN = DPAR(1)+DPAR(2)*SQRT(X) +DPAR(3)/X
//       UFIT=FITFUN
//       END
static float32 hbook_fit(const float32 *x)
{
  assert(x!=NULL);
  const vector<Double_t>
    xx(x          ,x          +fit_data.Nx  ),
    pp(hcfitd_.par,hcfitd_.par+fit_data.Npar);
  hcfitd_.fit_result = fit_data.fit_function(xx,pp);
  return hcfitd_.fit_result;
}
#endif

////////////////////////////////////////////////////////////////////////////////

CsHistBase::Dimension::Dimension(size_t bins,float64 min_,float64 max_,const string &label_)
: Nbins(bins), min(min_), max(max_), label(label_)
{
  if( Nbins==0 )
    CsErrLog::Instance()->mes( elFatal,"CsHistBase::Dimension::Dimension:  Nbins==0");

  if( min>=max )
    CsErrLog::Instance()->msg( elFatal,__FILE__,__LINE__,
       "CsHistBase::Dimension::Dimension:  min>=max:  min=%g, max=%g",min,max);
}

////////////////////////////////////////////////////////////////////////////////









////////////////////////////////////////////////////////////////////////////////
// 
//                           class CsHistBase
// 
////////////////////////////////////////////////////////////////////////////////












#define hist_dummy reinterpret_cast<CsHomeNamed*>(h)
#define hist_hbook reinterpret_cast<CsHbookHist*>(h)

////////////////////////////////////////////////////////////////////////////////

CsHistBase::~CsHistBase(void)
{
  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        break;

    case CsHistograms::HBOOK:
      #if USE_HBOOK
        delete hist_hbook;
        break;
      #endif
      goto dummy;

    dummy:
    case CsHistograms::DUMMY:
      delete hist_dummy;
      break;

    default: 
      CsErrLog::Instance()->mes( elFatal,"CsHistBase::~CsHistBase:  Internal error.");
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////

CsHistBase::CsHistBase(const string &name,const string &title) :
h(NULL)
{
  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        break;

    case CsHistograms::HBOOK:
      #if USE_HBOOK
        h = new CsHbookHist(name,title);
        hist_hbook->CsHomeNamed::SetPath(CsHistograms::GetCurrentPath());
        break;
      #endif
      goto dummy;

    dummy:
    case CsHistograms::DUMMY:
      h = new CsHomeNamed(name,title);
      return;

    default: 
      CsErrLog::Instance()->mes( elFatal,"CsHistBase::operator=  Internal error.");
  }
  return;
}

////////////////////////////////////////////////////////////////////////////////

CsHistBase &CsHistBase::operator = (const CsHistBase &hist)
{
  CsErrLog::Instance()->mes( elFatal,"CsHistBase::operator=  not implemented yet.");

//   if( this!=&hist )
//     switch( CsHistograms::GetImplementation() )
//     {
//       case CsHistograms::ROOT:
//           // tempory
//           CsErrLog::Instance()->mes( elFatal,"CsHistBase::operator=  No body yet.");
//           break;
// 
//       case CsHistograms::DUMMY:
//         ((CsHomeNamed*)h)->SetName (((CsHomeNamed*)hist.h)->GetName ());
//         ((CsHomeNamed*)h)->SetTitle(((CsHomeNamed*)hist.h)->GetTitle());
//         break;
// 
//       default: 
//         CsErrLog::Instance()->mes( elFatal,"CsHistBase::operator=  Internal error.");
//     }

  return *this;
}

////////////////////////////////////////////////////////////////////////////////

string CsHistBase::GetName(void) const
{
  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        return reinterpret_cast<TH1*>(h)->GetName();

    case CsHistograms::HBOOK:
      #if USE_HBOOK
        return hist_hbook->GetName();
        break;
      #endif
      goto dummy;

    dummy:
    case CsHistograms::DUMMY:
      return hist_dummy->GetName();
    
    default:
      CsErrLog::Instance()->mes( elFatal,"CsHistBase::GetName  Internal error.");
  }
  return string();  // Just making C++ compiler happy.
}

////////////////////////////////////////////////////////////////////////////////

void CsHistBase::SetName(const string &name)
{
  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        return reinterpret_cast<TH1*>(h)->SetName(name.c_str());
        break;

    case CsHistograms::HBOOK:
      #if USE_HBOOK
        hist_hbook->SetName(name);
        break;
      #endif
      goto dummy;

    dummy:
    case CsHistograms::DUMMY:
      hist_dummy->SetName(name);
      break;
    
    default:
      CsErrLog::Instance()->mes( elFatal,"CsHistBase::GetName  Internal error.");
  }
}

////////////////////////////////////////////////////////////////////////////////

string CsHistBase::GetTitle(void) const
{
  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        return reinterpret_cast<TH1*>(h)->GetTitle();

    case CsHistograms::HBOOK:
    {
      #if USE_HBOOK
        return hist_hbook->GetTitle();
      #endif
      goto dummy;
    }

    dummy:
    case CsHistograms::DUMMY:
      return hist_dummy->GetTitle();
    
    default:
      CsErrLog::Instance()->mes( elFatal,"CsHistBase::GetTitle  Internal error.");
  }
  return string();  // Just making C++ compiler happy.
}

////////////////////////////////////////////////////////////////////////////////

void CsHistBase::SetTitle(const string &title)
{
  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        return reinterpret_cast<TH1*>(h)->SetTitle(title.c_str());
        break;

    case CsHistograms::HBOOK:
      #if USE_HBOOK
        CsErrLog::Instance()->mes(elFatal,"CsHistBase::SetTitle: HBOOK: not implemented yet.");
        break;
      #endif
      goto dummy;

    dummy:
    case CsHistograms::DUMMY:
      hist_dummy->SetTitle(title);
      break;
    
    default:
      CsErrLog::Instance()->mes( elFatal,"CsHistBase::GetName  Internal error.");
  }
}

////////////////////////////////////////////////////////////////////////////////

string CsHistBase::GetPath(void) const
{
  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
      {
        TDirectory const * const dir = reinterpret_cast<TH1*>(h)->GetDirectory();
        if( dir==NULL )
          return "";
        else
        {
          // tempory.
          return dir->GetPath();
//           string p(dir->GetPath());     //  "file:/dir1/dir2"
//           string::size_type n=p.find(':');
//           if( string::npos==
        }
      }

    case CsHistograms::HBOOK:
      #if USE_HBOOK
        return hist_hbook->GetPath();
      #endif
      goto dummy;

    dummy:
    case CsHistograms::DUMMY:
      return hist_dummy->GetTitle();
    
    default:
      CsErrLog::Instance()->mes( elFatal,"CsHistBase::GetPath  Internal error.");
  }
  return string();  // Just making C++ compiler happy.
}

////////////////////////////////////////////////////////////////////////////////

void CsHistBase::SetPath(const string &path)
{
  CsErrLog::Instance()->mes( elFatal,"CsHistBase::SetPath: Not implemented yet.");

//   switch( CsHistograms::GetImplementation() )
//   {
//     case CsHistograms::ROOT:
//         CsErrLog::Instance()->mes( elFatal,"CsHistBase::SetHomePath:  Not implemented yet.");
//         break;
// 
//     case CsHistograms::DUMMY:
//       hist_dummy->SetPath(home_path);
//       break;
// 
//     default: 
//       CsErrLog::Instance()->mes(elFatal,"CsHistBase::SetHomePath  Internal error.");
//   }
}

////////////////////////////////////////////////////////////////////////////////

uint64 CsHistBase::GetEntries(void) const
{
  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        return (uint64) reinterpret_cast<TH1*>(h)->GetEntries();

    case CsHistograms::HBOOK:
      #if USE_HBOOK
        CsHistograms::SetCurrentPath(hist_hbook->GetPath());
        return hnoent(hist_hbook->ID());
      #endif
      goto dummy;

    dummy:
    case CsHistograms::DUMMY:
      return 0;

    default:
      CsErrLog::Instance()->mes( elFatal,"CsHistBase::GetEntries  Internal error.");
  }
  return 0;     // Just making C++ compiler happy.
}

////////////////////////////////////////////////////////////////////////////////

void CsHistBase::SetEntries(uint64 n)
{
  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        reinterpret_cast<TH1*>(h)->SetEntries(n);
        break;

    case CsHistograms::HBOOK:
      #if USE_HBOOK
        CsErrLog::Instance()->mes(elFatal,"CsHistBase::SetEntries: HBOOK: not implemented yet.");
        break;
      #endif
      goto dummy;

    dummy:
    case CsHistograms::DUMMY:
      break;

    default:
      CsErrLog::Instance()->mes( elFatal,"CsHistBase::SetEntries  Internal error.");
  }
}

////////////////////////////////////////////////////////////////////////////////

void CsHistBase::Print(void)
{
  cout << "Histogram \"" << GetName() << "\"\n";
}

////////////////////////////////////////////////////////////////////////////////

CsHistBase::Dimension CsHistBase::GetDim(uint16 n) const
{
  if( n<=0 || n>NDim() )
    CsErrLog::Instance()->msg( elFatal,__FILE__,__LINE__,
      "CsHistBase::Dimension::GetDim():  bad dimension number %d.  Must be 1<=ndim<=%d",n,NDim());

  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
    {
        TAxis const *a=NULL;
        switch(n)
        {
          case  1:  a = reinterpret_cast<TH1*>(h)->GetXaxis(); break;
          case  2:  a = reinterpret_cast<TH1*>(h)->GetYaxis(); break;
          case  3:  a = reinterpret_cast<TH1*>(h)->GetZaxis(); break;
          default:  CsErrLog::Instance()->msg(elFatal,__FILE__,__LINE__,
                      "CsHistBase::GetDim: ROOT: this dimension is not supported: %d",n);
        }
        return Dimension(a->GetNbins(),a->GetXmin(),a->GetXmax(),a->GetName());
    }

    case CsHistograms::HBOOK:
    {
      #if USE_HBOOK
        CsHistograms::SetCurrentPath(hist_hbook->GetPath());
        switch(n)
        {
          case  1:  return Dimension(hist_hbook->nx,hist_hbook->xmin,hist_hbook->xmax);
          case  2:  return Dimension(hist_hbook->ny,hist_hbook->ymin,hist_hbook->ymax);
          default:  CsErrLog::Instance()->msg(elFatal,__FILE__,__LINE__,
                      "CsHistBase::GetDim HBOOK: this dimension is not supported: %d",n);
        }
        break;
      #endif
      goto dummy;
    }

    dummy:
    case CsHistograms::DUMMY:
      return Dimension(1,0,1,"dummy");

    default:
      CsErrLog::Instance()->mes( elFatal,"CsHistBase::Dimension::GetDim()  Internal error.");
  }

  return Dimension(1,0,1,"dummy");
}

////////////////////////////////////////////////////////////////////////////////

float64 CsHistBase::Fit(FitFunction f,vector<float64> &param,const string &options) const
{
  CsErrLog::Instance()->mes(elFatal,"CsHistBase::Fit(): you can not call this method.");
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

float64 CsHistBase::Fit(FitFunction f,vector<CsFitParam> &param,const string &options) const
{
  CsErrLog::Instance()->mes(elFatal,"CsHistBase::Fit(): you can not call this method.");
  return 0;
}

////////////////////////////////////////////////////////////////////////////////













////////////////////////////////////////////////////////////////////////////////
// 
//                           class CsHist1<T>
// 
////////////////////////////////////////////////////////////////////////////////












////////////////////////////////////////////////////////////////////////////////

template<class T>
CsHist1<T>::~CsHist1(void)
{
  if( GetPath()=="" ) // Do not save the histogram.
    switch( CsHistograms::GetImplementation() )
    {
      case CsHistograms::ROOT:
          switch( sizeof(T) )
          {
            case 8: delete reinterpret_cast<TH1D*>(h); break;
            case 4: delete reinterpret_cast<TH1F*>(h); break;
            case 2: delete reinterpret_cast<TH1S*>(h); break;
            default:
              CsErrLog::Instance()->mes( elFatal,"CsHist1::~CsHist1()  Internal error (unknown size)");
          }
          break;

      case CsHistograms::HBOOK:
      {
        #if USE_HBOOK
          CsHistograms::SetCurrentPath(hist_hbook->GetPath());
          hdelet(hist_hbook->ID());
          break;
        #endif
        goto dummy;
      }

      dummy:
      case CsHistograms::DUMMY:
        break;

      default:
        CsErrLog::Instance()->mes( elFatal,"CsHist1::~CsHist1()  Internal error.");
    }
}

////////////////////////////////////////////////////////////////////////////////

template <class T>
CsHist1<T>::CsHist1(const string &name,const string &title,uint32 Xbins,float64 Xmin,float64 Xmax)
: CsHistBase(name,title)
{
  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        switch( sizeof(T) )
        {
          case 8: h = new TH1D(name.c_str(),title.c_str(),Xbins,Xmin,Xmax); break;
          case 4: h = new TH1F(name.c_str(),title.c_str(),Xbins,Xmin,Xmax); break;
          case 2: h = new TH1S(name.c_str(),title.c_str(),Xbins,Xmin,Xmax); break;
          default:
            CsErrLog::Instance()->mes( elFatal,"CsHist1::CsHist1()  Internal error (unknown size)");
        }
        break;

    case CsHistograms::HBOOK:
    {
      #if USE_HBOOK
        hbook1(hist_hbook->ID(),hist_hbook->GetHbookTitle().c_str(),Xbins,Xmin,Xmax);
        hist_hbook->nx   = Xbins;
        hist_hbook->xmin = Xmin;
        hist_hbook->xmax = Xmax;
        break;
      #endif
      goto dummy;
    }

    dummy:
    case CsHistograms::DUMMY:
      break;

    default:
      CsErrLog::Instance()->mes( elFatal,"CsHist1::CsHist1()  Internal error.");
  }
  
  CsHistograms::SetConsumedMemory( CsHistograms::GetConsumedMemory() + sizeof(T)*(Xbins+2) );

  #if debug_histograms
    cout << "CsHist1<T>::CsHist1(): Mem = " << CsHistograms::GetConsumedMemory() << "   " << name << "\n";
  #endif
}

////////////////////////////////////////////////////////////////////////////////

template <class T>
float64 CsHist1<T>::GetMean(void) const
{
  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        return (float64) reinterpret_cast<TH1*>(h)->GetMean();

    case CsHistograms::HBOOK:
      #if USE_HBOOK
        CsHistograms::SetCurrentPath(hist_hbook->GetPath());
        return hstati(hist_hbook->ID(),1,"HIST",0);
      #endif
      goto dummy;

    dummy:
    case CsHistograms::DUMMY:
      return 0;

    default:
      CsErrLog::Instance()->mes( elFatal,"CsHist1::GetMean():  Internal error.");
  }
  return 0;     // Just making C++ compiler happy.
}

////////////////////////////////////////////////////////////////////////////////

template <class T>
float64 CsHist1<T>::GetRMS(void) const
{
  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        return (float64) reinterpret_cast<TH1*>(h)->GetRMS();

    case CsHistograms::HBOOK:
      #if USE_HBOOK
        CsHistograms::SetCurrentPath(hist_hbook->GetPath());
        return hstati(hist_hbook->ID(),2,"HIST",0);
      #endif
      goto dummy;

    dummy:
    case CsHistograms::DUMMY:
      return 0;

    default:
      CsErrLog::Instance()->mes( elFatal,"CsHist1:GetRMS():  Internal error.");
  }
  return 0;     // Just making C++ compiler happy.
}

////////////////////////////////////////////////////////////////////////////////

template <class T>
void CsHist1<T>::Fill(float64 x1,float64 w)
{
  // silently ignore call to Fill() if object isn't instantiated
  if ( !this )
    return;

  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        reinterpret_cast<TH1*>(h)->Fill(x1,w);
        break;

    case CsHistograms::HBOOK:
      #if USE_HBOOK
        CsHistograms::SetCurrentPath(hist_hbook->GetPath());
        hf1(hist_hbook->ID(),x1,w);
        break;
      #endif
      goto dummy;

    dummy:
    case CsHistograms::DUMMY:
      break;

    default:
      CsErrLog::Instance()->mes( elFatal,"CsHist1::Fill()  Internal error.");
  }
}

////////////////////////////////////////////////////////////////////////////////

template <class T>
size_t CsHist1<T>::GetBin(float64 x1) const
{
  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        return reinterpret_cast<TH1*>(h)->FindBin(x1);

    case CsHistograms::HBOOK:
    {
      #if USE_HBOOK
        CsErrLog::Instance()->mes(elFatal,"CsHist1::GetBin: HBOOK: not implemented yet.");
//         int n;
//         hxyij( hist_hbook->ID(), x1,0, &n, NULL);
//         return (size_t)n;
      #endif
      goto dummy;
    }

    dummy:
    case CsHistograms::DUMMY:
      return size_t(-1);

    default:
      CsErrLog::Instance()->mes( elFatal,"CsHist1::GetBin()  Internal error.");
  }
  return size_t(-1);
}

////////////////////////////////////////////////////////////////////////////////

template <class T>
T CsHist1<T>::GetBinContent(size_t bin) const
{
  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        return static_cast<T>( reinterpret_cast<TH1*>(h)->GetBinContent(bin) );

    case CsHistograms::HBOOK:
      #if USE_HBOOK
        CsHistograms::SetCurrentPath(hist_hbook->GetPath());
        return static_cast<T>( hi(hist_hbook->ID(),bin) );
      #endif
      goto dummy;

    dummy:
    case CsHistograms::DUMMY:
      return 0;

    default:
      CsErrLog::Instance()->mes( elFatal,"CsHist1::GetBinContent()  Internal error.");
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

template <class T>
void CsHist1<T>::SetBinContent(size_t bin,T v)
{
  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        switch( sizeof(T) )
        {
          case 8: reinterpret_cast<TH1D*>(h)->SetBinContent(bin,v); break;
          case 4: reinterpret_cast<TH1F*>(h)->SetBinContent(bin,v); break;
          case 2: reinterpret_cast<TH1S*>(h)->SetBinContent(bin,v); break;
          default:
            CsErrLog::Instance()->mes( elFatal,"CsHist1::SetBinContent()  Internal error (unknown size)");
        }
        break;

    case CsHistograms::HBOOK:
      #if USE_HBOOK
        CsErrLog::Instance()->mes(elFatal,"CsHistBase::SetBinContent: HBOOK: not implemented yet.");
      #endif
      goto dummy;

    dummy:
    case CsHistograms::DUMMY:
      break;

    default:
      CsErrLog::Instance()->mes( elFatal,"CsHist1D::SetBinContent  Internal error.");
  }
}

////////////////////////////////////////////////////////////////////////////////

template <class T>
float64 CsHist1<T>::Fit(FitFunction function,vector<CsFitParam> &par,const string &opt,float64 min,float64 max) const
{
  Dimension d=GetDim(1);
  if( min<d.GetMin() || max>d.GetMax() )
    CsErrLog::Instance()->msg(elWarning,__FILE__,__LINE__,"CsHist1::Fit():  histogram range is [%g,%g],  function range is [%g,%g]",
                                d.GetMin(),d.GetMax(), min,max );
  fit_data.fit_function = function;
  fit_data.Nx           = 1;
  fit_data.Npar         = par.size();

  string options(opt);

  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
    {
        options += "RN";  // Range,Nodraw
        TF1 f("fit_function",root_fit,min,max,par.size());

        for( size_t i=0; i<par.size(); i++ )
        {
          f.SetParameter(i,par[i].GetStart());
          f.SetParName  (i,par[i].GetName().c_str());
          f.SetParError (i,par[i].GetError());
          f.SetParLimits(i,par[i].GetMin(),par[i].GetMax());
        }

        reinterpret_cast<TH1*>(h)->Fit("fit_function",options.c_str(),"");

        for( size_t i=0; i<par.size(); i++ )
        {
          par[i].SetFinish(f.GetParameter(i));
          par[i].SetError (f.GetParError(i));
        }

        return f.GetChisquare();
    }

    case CsHistograms::HBOOK:
    {
      #if USE_HBOOK
        options += "RNB";  // range,Nodraw,bondary
        const size_t par_max=24;
        if( par.size()>par_max )
          CsErrLog::Instance()->msg(elFatal,__FILE__,__LINE__,"CsHist1::Fit():  too many (%d) parameters. max=%d",
                                    par.size(),par_max);
        float32 par2[par.size()], step[par.size()], pmin[par.size()], pmax[par.size()], par2_err[par.size()];
        for( size_t i=0; i<par.size(); i++ )
        {
          par2[i] = par[i].GetStart();
          step[i] = par[i].GetStep();
          pmin[i] = par[i].GetMin();
          pmax[i] = par[i].GetMax();
        }
        int ddddd;
        CsHistograms::SetCurrentPath(hist_hbook->GetPath());
        hxyij( hist_hbook->ID(), min, 0, &quest_[10], &ddddd );
        hxyij( hist_hbook->ID(), max, 0, &quest_[11], &ddddd );
        float32 result;
        hfith(hist_hbook->ID(),hbook_fit,options.c_str(),par.size(),par2,step,pmin,pmax,par2_err,result);
        for( size_t i=0; i<par.size(); i++ )
        {
          par[i].SetFinish(par2    [i]);
          par[i].SetError (par2_err[i]);
        }
        return result;
      #endif
      goto dummy;
    }

    dummy:
    case CsHistograms::DUMMY:
      break;

    default:
      CsErrLog::Instance()->mes( elFatal,"CsHist1::Fit()  Internal error.");
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////////////

template <class T>
T CsHist1<T>::GetBinError(size_t bin) const
{
  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        return static_cast<T>( reinterpret_cast<TH1*>(h)->GetBinError(bin) );

    case CsHistograms::HBOOK:
      #if USE_HBOOK
        CsHistograms::SetCurrentPath(hist_hbook->GetPath());
        return static_cast<T>( hie(hist_hbook->ID(),bin) );
      #endif
      goto dummy;

    dummy:
    case CsHistograms::DUMMY:
      return 0;

    default:
      CsErrLog::Instance()->mes( elFatal,"CsHist1::BinContentType()  Internal error.");
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

template <class T>
void CsHist1<T>::SetBinError(size_t bin,CsHist1<T>::BinContentType v)
{
  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        reinterpret_cast<TH1*>(h)->SetBinError(bin,v);
        break;

    case CsHistograms::HBOOK:
      #if USE_HBOOK
        CsErrLog::Instance()->mes(elFatal,"CsHist1::SetBinError(): HBOOK: not implemented yet.");
      #endif
      goto dummy;

    dummy:
    case CsHistograms::DUMMY:
      break;

    default:
      CsErrLog::Instance()->mes( elFatal,"CsHist1::SetBinError()  Internal error.");
  }
}

////////////////////////////////////////////////////////////////////////////////

template <class T>
CsHist1<T> &CsHist1<T>::operator = (const CsHist1 &hist)
{
  CsErrLog::Instance()->mes( elFatal,"CsHist1::operator=  not implemented yet.");

//   if( this!=&hist )
//     switch( CsHistograms::GetImplementation() )
//     {
//       case CsHistograms::ROOT:
//           // tempory
//           CsErrLog::Instance()->mes( elFatal,"CsHistBase::operator=  No body yet.");
//           break;
// 
//       case CsHistograms::DUMMY:
//         ((CsHomeNamed*)h)->SetName (((CsHomeNamed*)hist.h)->GetName ());
//         ((CsHomeNamed*)h)->SetTitle(((CsHomeNamed*)hist.h)->GetTitle());
//         break;
// 
//       default: 
//         CsErrLog::Instance()->mes( elFatal,"CsHistBase::operator=  Internal error.");
//     }

  return *this;
}


////////////////////////////////////////////////////////////////////////////////














////////////////////////////////////////////////////////////////////////////////
// 
//                           class CsHist2<T>
// 
////////////////////////////////////////////////////////////////////////////////












////////////////////////////////////////////////////////////////////////////////

template <class T>
CsHist2<T>::~CsHist2(void)
{
  if( GetPath()=="" ) // Do not save the histogram.
    switch( CsHistograms::GetImplementation() )
    {
      case CsHistograms::ROOT:
          switch( sizeof(T) )
          {
            case 8: delete reinterpret_cast<TH2D*>(h); break;
            case 4: delete reinterpret_cast<TH2F*>(h); break;
            case 2: delete reinterpret_cast<TH2S*>(h); break;
            default:
              CsErrLog::Instance()->mes( elFatal,"CsHist2::~CsHist2  Internal error (unknown size)");
          }
          break;

      case CsHistograms::HBOOK:
      {
        #if USE_HBOOK
          CsHistograms::SetCurrentPath(hist_hbook->GetPath());
          hdelet(hist_hbook->ID());
          break;
        #endif
        goto dummy;
      }

      dummy:
      case CsHistograms::DUMMY:
        break;

      default:
        CsErrLog::Instance()->mes( elFatal,"CsHist2D::~CsHist2D  Internal error.");
    }
}

////////////////////////////////////////////////////////////////////////////////

template <class T>
CsHist2<T>::CsHist2(const string &name,const string &title,uint32 Xbins,float64 Xmin,float64 Xmax,
                                                           uint32 Ybins,float64 Ymin,float64 Ymax)
: CsHistBase(name,title)
{
  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        switch( sizeof(T) )
        {
          case 8: h = new TH2D(name.c_str(),title.c_str(),Xbins,Xmin,Xmax,Ybins,Ymin,Ymax); break;
          case 4: h = new TH2F(name.c_str(),title.c_str(),Xbins,Xmin,Xmax,Ybins,Ymin,Ymax); break;
          case 2: h = new TH2S(name.c_str(),title.c_str(),Xbins,Xmin,Xmax,Ybins,Ymin,Ymax); break;
          default:
            CsErrLog::Instance()->mes( elFatal,"CsHist2::CsHist2  Internal error (unknown size)");
        }
        break;

    case CsHistograms::HBOOK:
    {
      #if USE_HBOOK
        hbook2(hist_hbook->ID(),hist_hbook->GetHbookTitle().c_str(),Xbins,Xmin,Xmax,Ybins,Ymin,Ymax);
        hist_hbook->nx   = Xbins;
        hist_hbook->xmin = Xmin;
        hist_hbook->xmax = Xmax;
        hist_hbook->ny   = Ybins;
        hist_hbook->ymin = Ymin;
        hist_hbook->ymax = Ymax;
        break;
      #endif
      goto dummy;
    }

    dummy:
    case CsHistograms::DUMMY:
      break;

    default:
      CsErrLog::Instance()->mes( elFatal,"CsHist2::CsHist2  Internal error.");
  }

  CsHistograms::SetConsumedMemory( CsHistograms::GetConsumedMemory() + sizeof(T)*(Xbins+2)*(Ybins+2) );

  #if debug_histograms
    cout << "CsHist2<T>::CsHist2(): Mem = " << CsHistograms::GetConsumedMemory() << "   " << name << "\n";
  #endif
}

////////////////////////////////////////////////////////////////////////////////

template <class T>
float64 CsHist2<T>::GetMean(size_t axis) const
{
  if( axis!=1 && axis!=2 )
    CsErrLog::Instance()->msg(elFatal,__FILE__,__LINE__,
          "CsHist2::GetMean():  Bad axis number %d",axis);

  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        return (float64) reinterpret_cast<TH1*>(h)->GetMean(axis);

    case CsHistograms::HBOOK:
      #if USE_HBOOK
        CsErrLog::Instance()->mes(elFatal,"CsHist2::GetMean(): HBOOK: not implemented yet.");
      #endif
      goto dummy;

    dummy:
    case CsHistograms::DUMMY:
      return 0;

    default:
      CsErrLog::Instance()->mes( elFatal,"CsHist2::GetMean():  Internal error.");
  }
  return 0;     // Just making C++ compiler happy.
}

////////////////////////////////////////////////////////////////////////////////

template <class T>
float64 CsHist2<T>::GetRMS(size_t axis) const
{
  if( axis!=1 && axis!=2 )
    CsErrLog::Instance()->msg(elFatal,__FILE__,__LINE__,
          "CsHist2::GetRMS():  Bad axis number %d",axis);

  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        return (float64) reinterpret_cast<TH1*>(h)->GetRMS(axis);

    case CsHistograms::HBOOK:
      #if USE_HBOOK
        CsErrLog::Instance()->mes(elFatal,"CsHist2::GetRMS(): HBOOK: not implemented yet.");
      #endif
      goto dummy;

    dummy:
    case CsHistograms::DUMMY:
      return 0;

    default:
      CsErrLog::Instance()->mes( elFatal,"CsHist2:GetRMS():  Internal error.");
  }
  return 0;     // Just making C++ compiler happy.
}

////////////////////////////////////////////////////////////////////////////////

template <class T>
void CsHist2<T>::Fill(float64 x1,float64 x2, float64 w)
{
  // silently ignore call to Fill() if object isn't instantiated
  if ( !this )
    return;

  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        switch( sizeof(T) )
        {
          case 8: reinterpret_cast<TH2D*>(h)->Fill(x1,x2,w); break;
          case 4: reinterpret_cast<TH2F*>(h)->Fill(x1,x2,w); break;
          case 2: reinterpret_cast<TH2S*>(h)->Fill(x1,x2,w); break;
          default: throw "CsHist2::Fill()  internal error";
        }
      break;

    case CsHistograms::HBOOK:
      #if USE_HBOOK
        CsHistograms::SetCurrentPath(hist_hbook->GetPath());
        hf2(hist_hbook->ID(),x1,x2,w);
        break;
      #endif
      goto dummy;

    dummy:
    case CsHistograms::DUMMY:
      break;

    default:
      CsErrLog::Instance()->mes( elFatal,"CsHist2::Fill  Internal error.");
  }
}

////////////////////////////////////////////////////////////////////////////////

template <class T>
size_t CsHist2<T>::GetBin(float64 x1,float64 x2) const
{
  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        switch( sizeof(T) )
        {
          case 8: return reinterpret_cast<TH2D*>(h)->FindBin(x1,x2);
          case 4: return reinterpret_cast<TH2F*>(h)->FindBin(x1,x2);
          case 2: return reinterpret_cast<TH2S*>(h)->FindBin(x1,x2);
          default: throw "CsHist2::Fit()  internal error";
        }

    {
      #if USE_HBOOK
        CsErrLog::Instance()->mes(elFatal,"CsHist1::GetBin: HBOOK: not implemented yet.");
//         size_t n1,n2;
//         hxyij( hist_hbook->ID(), x1,x2, &n1, &n2 );
//         return n1+n2*();
      #endif
      goto dummy;
    }

    dummy:
    case CsHistograms::DUMMY:
      return size_t(-1);

    default:
      CsErrLog::Instance()->mes( elFatal,"CsHist2::GetBin  Internal error.");
  }
  return size_t(-1);
}

////////////////////////////////////////////////////////////////////////////////

template <class T>
T CsHist2<T>::GetBinContent(size_t bin) const
{
  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        switch( sizeof(T) )
        {
          case 8: return static_cast<T>(reinterpret_cast<TH2D*>(h)->GetBinContent(bin));
          case 4: return static_cast<T>(reinterpret_cast<TH2F*>(h)->GetBinContent(bin));
          case 2: return static_cast<T>(reinterpret_cast<TH2S*>(h)->GetBinContent(bin));
          default: throw "CsHist2::GetBinContent()  internal error";
        }

    case CsHistograms::HBOOK:
    {
      #if USE_HBOOK
        CsHistograms::SetCurrentPath(hist_hbook->GetPath());
        size_t nx = GetDim(1).GetNBins()+2;
        return static_cast<T>(hij(hist_hbook->ID(),bin%nx,bin/nx));
      #endif
      goto dummy;
    }

    dummy:
    case CsHistograms::DUMMY:
      return 0;

    default:
      CsErrLog::Instance()->mes( elFatal,"CsHist2::GetBinContent  Internal error.");
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

template <class T>
void CsHist2<T>::SetBinContent(size_t bin,T v)
{
  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        switch( sizeof(T) )
        {
          case 8: reinterpret_cast<TH2D*>(h)->SetBinContent(bin,v); break;
          case 4: reinterpret_cast<TH2F*>(h)->SetBinContent(bin,v); break;
          case 2: reinterpret_cast<TH2S*>(h)->SetBinContent(bin,v); break;
          default: throw "CsHist2::SetBinContent()  internal error";
        }
        break;

    case CsHistograms::HBOOK:
      #if USE_HBOOK
        CsErrLog::Instance()->mes(elFatal,"CsHist2::SetBinContent(): HBOOK: not implemented yet.");
      #endif
      goto dummy;

    dummy:
    case CsHistograms::DUMMY:
      break;

    default:
      CsErrLog::Instance()->mes( elFatal,"CsHist2::SetBinContent():  Internal error.");
  }
}

////////////////////////////////////////////////////////////////////////////////

template <class T>
T CsHist2<T>::GetBinError(size_t bin) const
{
  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        switch( sizeof(T) )
        {
          case 8: return static_cast<T>( reinterpret_cast<TH2D*>(h)->GetBinError(bin) );
          case 4: return static_cast<T>( reinterpret_cast<TH2F*>(h)->GetBinError(bin) );
          case 2: return static_cast<T>( reinterpret_cast<TH2S*>(h)->GetBinError(bin) );
          default: throw "CsHist2::GetBinError()  internal error";
        }

    case CsHistograms::HBOOK:
    {
      #if USE_HBOOK
        CsHistograms::SetCurrentPath(hist_hbook->GetPath());
        size_t nx = GetDim(1).GetNBins()+2;
        return static_cast<T>( hije(hist_hbook->ID(),bin%nx,bin/nx) );
      #endif
      goto dummy;
    }

    dummy:
    case CsHistograms::DUMMY:
      return 0;

    default:
      CsErrLog::Instance()->mes( elFatal,"CsHist2::GetBinError()  Internal error.");
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

template <class T>
void CsHist2<T>::SetBinError(size_t bin,CsHist2<T>::BinContentType v)
{
  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
        switch( sizeof(T) )
        {
          case 8: reinterpret_cast<TH2D*>(h)->SetBinError(bin,v); break;
          case 4: reinterpret_cast<TH2F*>(h)->SetBinError(bin,v); break;
          case 2: reinterpret_cast<TH2S*>(h)->SetBinError(bin,v); break;
          default: throw "CsHist2::SetBinError()  internal error";
        }
        break;

    case CsHistograms::HBOOK:
      #if USE_HBOOK
        CsErrLog::Instance()->mes(elFatal,"CsHist2::SetBinError: HBOOK: not implemented yet.");
      #endif
      goto dummy;

    dummy:
    case CsHistograms::DUMMY:
      break;

    default:
      CsErrLog::Instance()->mes( elFatal,"CsHist2::SetBinError  Internal error.");
  }
}

////////////////////////////////////////////////////////////////////////////////

template <class T>
float64 CsHist2<T>::Fit(FitFunction function,vector<CsFitParam> &par,const string &opt,
                        float64 Xmin,float64 Xmax,float64 Ymin,float64 Ymax) const
{
  Dimension d1=GetDim(1), d2=GetDim(2);
  if( Xmin<d1.GetMin() || Xmax>d1.GetMax() || Ymin<d2.GetMin() || Ymax>d2.GetMax() )
    CsErrLog::Instance()->msg(elWarning,__FILE__,__LINE__,"CsHist2::Fit():  histogram range is [%g,%g]x[%g,%g],  function range is [%g,%g]x[%g,%g]",
                                d1.GetMin(),d1.GetMax(),d2.GetMin(),d2.GetMax(), Xmin,Xmax,Ymin,Ymax );
  fit_data.fit_function = function;
  fit_data.Nx           = 2;
  fit_data.Npar         = par.size();

  string options(opt);

  switch( CsHistograms::GetImplementation() )
  {
    case CsHistograms::ROOT:
    {
        options += "RN";  // Range,Nodraw
        TF2 f("fit_function",root_fit,Xmin,Xmax,Ymin,Ymax,par.size());

        for( size_t i=0; i<par.size(); i++ )
        {
          f.SetParameter(i,par[i].GetStart());
          f.SetParName  (i,par[i].GetName().c_str());
          f.SetParError (i,par[i].GetError());
          f.SetParLimits(i,par[i].GetMin(),par[i].GetMax());
        }

        switch( sizeof(T) )
        {
          case 8: reinterpret_cast<TH2D*>(h)->Fit("fit_function",options.c_str(),""); break;
          case 4: reinterpret_cast<TH2F*>(h)->Fit("fit_function",options.c_str(),""); break;
          case 2: reinterpret_cast<TH2S*>(h)->Fit("fit_function",options.c_str(),""); break;
          default: throw "CsHist2::Fit()  internal error";
        }

        for( size_t i=0; i<par.size(); i++ )
        {
          par[i].SetFinish(f.GetParameter(i));
          par[i].SetError (f.GetParError(i));
        }

        return f.GetChisquare();
    }

    case CsHistograms::HBOOK:
    {
      #if USE_HBOOK
        options += "RNB";  // range,Nodraw,bondary
        const size_t par_max=24;
        if( par.size()>par_max )
          CsErrLog::Instance()->msg(elFatal,__FILE__,__LINE__,"CsHist2::Fit():  too many (%d) parameters. max=%d",
                                    par.size(),par_max);
        float32 par2[par.size()], step[par.size()], pmin[par.size()], pmax[par.size()], par2_err[par.size()];
        for( size_t i=0; i<par.size(); i++ )
        {
          par2[i] = par[i].GetStart();
          step[i] = par[i].GetStep();
          pmin[i] = par[i].GetMin();
          pmax[i] = par[i].GetMax();
        }
        CsHistograms::SetCurrentPath(hist_hbook->GetPath());
        hxyij( hist_hbook->ID(), Xmin, Ymin, &quest_[10], &quest_[12] );
        hxyij( hist_hbook->ID(), Xmax, Ymax, &quest_[11], &quest_[13] );
        float32 result;
        hfith(hist_hbook->ID(),hbook_fit,options.c_str(),par.size(),par2,step,pmin,pmax,par2_err,result);
        for( size_t i=0; i<par.size(); i++ )
        {
          par[i].SetFinish(par2    [i]);
          par[i].SetError (par2_err[i]);
        }
        return result;
      #endif
      goto dummy;
    }

    dummy:
    case CsHistograms::DUMMY:
      break;

    default:
      CsErrLog::Instance()->mes( elFatal,"CsHist2::Fit():  Internal error.");
  }

  return 0;
}


////////////////////////////////////////////////////////////////////////////////

template <class T>
CsHist2<T> &CsHist2<T>::operator = (const CsHist2 &hist)
{
  CsErrLog::Instance()->mes( elFatal,"CsHist2::operator=  not implemented yet.");

//   if( this!=&hist )
//     switch( CsHistograms::GetImplementation() )
//     {
//       case CsHistograms::ROOT:
//           // tempory
//           CsErrLog::Instance()->mes( elFatal,"CsHistBase::operator=  No body yet.");
//           break;
// 
//       case CsHistograms::DUMMY:
//         ((CsHomeNamed*)h)->SetName (((CsHomeNamed*)hist.h)->GetName ());
//         ((CsHomeNamed*)h)->SetTitle(((CsHomeNamed*)hist.h)->GetTitle());
//         break;
// 
//       default: 
//         CsErrLog::Instance()->mes( elFatal,"CsHistBase::operator=  Internal error.");
//     }

  return *this;
}


////////////////////////////////////////////////////////////////////////////////











// Let instantiate these templates explicitly

template class CsHist1<int16  >;
template class CsHist2<int16  >;
template class CsHist1<float32>;
template class CsHist2<float32>;
template class CsHist1<float64>;
template class CsHist2<float64>;




////////////////////////////////////////////////////////////////////////////////

static void implementation_of_methods(void)
{
  CsHist1D h1D("","",1,1,1);
  CsHist1F h1F("","",1,1,1);
  CsHist1S h1S("","",1,1,1);

  CsHist2D h2D("","",1,1,1,2,2,2);
  CsHist2F h2F("","",1,1,1,2,2,2);
  CsHist2S h2S("","",1,1,1,2,2,2);

  FitFunction *f=NULL;
  vector<CsFitParam> par;

  #define code(h)               \
    h.Fit(*f,par);              \
    h.GetBinContent(0);         \
    h.SetBinContent(0,0);       \
    h.GetBinError(0);           \
    h.SetBinError(0,0);

  #define code1(h)              \
    code(h);                    \
    h.Fill(1);                  \
    h.GetBin(0);                \
    h.GetMean();                \
    h.GetRMS();

  #define code2(h)              \
    code(h);                    \
    h.Fill(1,1);                \
    h.GetBin(0,0);              \
    h.GetMean(1);               \
    h.GetRMS(1);
    

  code1(h1D); code2(h2D); 
  code1(h1F); code2(h2F);
  code1(h1S); code2(h2S);
  
}

////////////////////////////////////////////////////////////////////////////////

#undef hist_root
