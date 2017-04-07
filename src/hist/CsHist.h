/*!
   \file    CsHist.h
   \brief   Interface to histograms.
   \version $Revision: 1.20 $
   \author  Alexander Zvyagin
   \date    $Date: 2010/06/18 10:44:21 $
*/

#ifndef CsHist_h
#define CsHist_h

#include <iostream>
#include <cstdio>
#include "coral_config.h"

#include <vector>
#include <string>

#include "CsTypes.h"
#include "CsHome.h"

////////////////////////////////////////////////////////////////////////////////

/*! \brief Fit function type.

    This is prototype of functions that will be used in a histogram fitting procedure.

    \param x   vector of arguments
    \param par vector of parameters
    
    This is function example:

    \verbatim
      float64 fit_fun(const std::vector<float64> &x,const std::vector<float64> &par)
      {
        // f(x) = a+b*x

        if( x.size()!=1 || par.size()!=2 )
          throw "fit_fun: need 1 argument and 2 parameters");

        return par[0] + par[1]*x[0];
      }
    \endverbatim
*/
typedef float64     (*FitFunction) (const std::vector<float64> &x,const std::vector<float64> &par);

////////////////////////////////////////////////////////////////////////////////

/*! \class CsFitParam
    \author Alexander Zvyagin
    \brief Fit parameter
    
*/
class CsFitParam
{
  // ---------------------------------------------------------------------------
  // Constructors, destructor
  // ---------------------------------------------------------------------------

  public:

    /// Construct parameter from given \c float value.
                        CsFitParam              (float64 v);

    /*! \brief Construct parameter from given properties.
    
        \param name         parameter's name
        \param value_start  initial parameter value
        \param step         initial parameter step size
        \param min          lower bondary
        \param max          upper bondary
    */
    
                        CsFitParam              (const std::string name="",float64 value_start=0,float64 step=1,float64 min=0, float64 max=0);

  // ---------------------------------------------------------------------------
  // Operators
  // ---------------------------------------------------------------------------

  public:

    /// Cast operator - return final parameter value
                        operator float64        (void) {return value_finish;}

    /// Print parameter properties.
    friend std::ostream     &operator <<             (std::ostream &o,const CsFitParam &p) {p.Print(o);return o;}

  // ---------------------------------------------------------------------------
  // Methods
  // ---------------------------------------------------------------------------

  public:

    /// Print parameter properties.
    void                Print                   (std::ostream &o=std::cout) const;

    /// \return Parameter's name
    const std::string       &GetName                 (void) const {return name;}

    /// \return Start value
    float64             GetStart                (void) const {return value_start;}

    /// \return Finish value (after a fitting, for example)
    float64             GetFinish               (void) const {return value_finish;}

    /// \return Parameter error value (after a fitting, for example)
    float64             GetError                (void) const {return error;}

    /// \return Initial step size (for a fitting, for example)
    float64             GetStep                 (void) const {return step;}

    /// \return Lower boundary
    float64             GetMin                  (void) const {return min;}

    /// \return Upper boundary
    float64             GetMax                  (void) const {return max;}
    
    /// \brief Set name
    void                SetName                 (const std::string &v) {name        =v;}

    /// \brief Set initial value
    void                SetStart                (float64 v)       {value_start =v;}

    /// \brief Set final value
    void                SetFinish               (float64 v)       {value_finish=v;}

    /// \brief Set parameter error
    void                SetError                (float64 v)       {error       =v;}

    /// Set step size
    void                SetStep                 (float64 v)       {step        =v;}

    /// Set lower boundary
    void                SetMin                  (float64 v)       {min         =v;}

    /// Set upper boundary
    void                SetMax                  (float64 v)       {max         =v;}

    /// Set both lower and upper boundaries
    void                SetMinMax               (float64 v_min,float64 v_max) {min=v_min; max=v_max;}

  // ---------------------------------------------------------------------------
  // Attributes, data.
  // ---------------------------------------------------------------------------

  private:

    /// Parameter's name
    std::string              name;

    /// Initial value
    float64             value_start;

    /// Final value
    float64             value_finish;

    /// Parameter's error
    float64             error;

    /// Step size
    float64             step;
    
    /// Lower and upper boundaries
    float64             min, max;
};

////////////////////////////////////////////////////////////////////////////////

/*! \class CsHistBase
    \author Alexander Zvyagin
    \brief Base histogram class.
    
    This is abstract class with general histogram properties.
*/
class CsHistBase : public CsHome
{
  // ---------------------------------------------------------------------------
  // Types.
  // ---------------------------------------------------------------------------

  public:

    /*! \class Dimension
        \author Alexander Zvyagin
        \brief A histogram dimension.
    */
    class Dimension
    {
      public:

        /// Destructor.
                       ~Dimension               (void) {}

        /// Base constructor. Two extra bins will be added for underflow and overflow.
                        Dimension               (size_t bins,float64 min_,float64 max_,const std::string &label_="");

        /// \return bins amount. There are two extra bins for underflow and overflow.
        size_t          GetNBins                (void) const { return Nbins; }
        
        /// \return x minimum
        float64         GetMin                  (void) const { return min;   }

        /// \return x maximum
        float64         GetMax                  (void) const { return max;   }

        /*! \brief \return bin number for this value.
           \return 0 for underflow
           \return GetNBins()+1 for overflow
           \return 1...GetNBins() for value \c x inside [min,max).
        */
        size_t          GetBin                  (float64 x) const;

        /// Print dimension information.
        void            Print                   (void) const;

        /// Print dimension information.
        friend std::ostream &operator <<             (std::ostream &o,const Dimension &d);

      private:

        /// Bins amount.
        size_t          Nbins;

        /// Minimum.
        float64         min;

        /// Maximum.
        float64         max;
        
        /// Dimension's label.
        std::string          label;
    };

  // ---------------------------------------------------------------------------
  // Constructors, destructor.
  // ---------------------------------------------------------------------------

  protected:

    /// Destructor.
    virtual            ~CsHistBase              (void);

    /// Empty constructor.
                        CsHistBase              (void) { };

    /// Copy constructor.
                        CsHistBase              (const CsHistBase &h);

    /// Base constructor.
                        CsHistBase              (const std::string &name,const std::string &title);

  // ---------------------------------------------------------------------------
  // Operators
  // ---------------------------------------------------------------------------

  protected:
    
    /// Assignment operator.
    CsHistBase            &operator =              (const CsHistBase &h);

  // ---------------------------------------------------------------------------
  // Methods
  // ---------------------------------------------------------------------------
  
  public:
    
    /// \return Histogram name.
    std::string              GetName                 (void) const;

    /// Change histogram name.
    void                SetName                 (const std::string &name);

    /// \return Histogram title.
    std::string              GetTitle                (void) const;

    /// Change histogram title.
    void                SetTitle                (const std::string &title);

    /// \return Path where histogram will be stored.
    std::string              GetPath                 (void) const;
    
    /*! Set place where a histogram will be stored.

        Example home_path is "dir1/dir2/dir3". If some directories do not exist
        they will be created. The histogram will be moved from current directory
        to home_path.

        If home_path starts from symbol \b / the path will be ralative to
        histograms package home directory, see CsHistograms::GetHomePath().
        Else the path will be relative to current directory, see
        CsHistograms::GetCurrentPath(), CsHistograms::SetCurrentPath().
    */
    void                SetPath                 (const std::string &home_path);

    /// \return Information of histogram dimension \c n (1,2,..).
    Dimension           GetDim                  (uint16 n) const;
    
    /// \return Histogram dimension size.
    virtual uint16      NDim                    (void) const = 0;
    
    /// \return Entries amount.
    uint64              GetEntries              (void) const;
    
    /// Set entries amount.
    void                SetEntries              (uint64 n);

    /*! \return Pointer to raw histogram. This method never returns NULL.
        Use this pointer carefully. It is user responsibility to know what is
        object type on an adrress of the pointer.
        
        If you uses DUMMY histograms package, then you get useless (for you)
        pointer to internal object of dummy class (see CsHist.cc source file).
        If you use "histograms package ROOT", the following rules are applied:
        \arg if your histogram is CsHist1D, the pointer has type TH1D*
        \arg if your histogram is CsHist2D, the pointer has type TH2D*
        \arg if your histogram is CsHist1F, the pointer has type TH1F*
        \arg if your histogram is CsHist2F, the pointer has type TH2F*
        \arg if your histogram is CsHist1S, the pointer has type TH1S*
        \arg if your histogram is CsHist2S, the pointer has type TH2S*
    */
    virtual void       *ThePointer              (void)       {return h;}
    
    /// Constant version of the previous function.
    virtual const void *ThePointer              (void) const {return h;}

    /// Print histogram's information.
    virtual void        Print                   (void);
    
    /*! \brief Histogram fit.
         Available options are
         \arg 'Q' - quiet mode (no print)
         \arg 'V' - verbose mode
         \arg 'L' - use Loglikelihood method (default is chisquare method)
         \arg 'E' - perform better Errors estimation using Minos technique
         \arg 'W' - set all bin errors to 1
    */
    virtual float64     Fit                     (FitFunction f,std::vector<float64> &param,
                                                 const std::string &options="") const;

    /// Almost the same method as the above. \sa CsFitParam
    virtual float64     Fit                     (FitFunction f,std::vector<CsFitParam> &param,
                                                 const std::string &options="") const;

  // ---------------------------------------------------------------------------
  // Attributes
  // ---------------------------------------------------------------------------
    
  protected:
    
    /// The pointer to real histogram object.
    void               *h;
};

////////////////////////////////////////////////////////////////////////////////

/*! \class CsHist1
    \author Alexander Zvyagin
    \brief 1-dimensional histogram.
*/
template <class T>
class CsHist1 : public CsHistBase
{
  // ---------------------------------------------------------------------------
  // Types.
  // ---------------------------------------------------------------------------

  public:
  
    typedef             T                       BinContentType;

  public:

  // ---------------------------------------------------------------------------
  // Constructors, destructor.
  // ---------------------------------------------------------------------------

    /// Destructor.
    virtual            ~CsHist1                 (void);

    /// Copy constructor.
                        CsHist1                 (const CsHist1 &h);

    /// Base constructor.
                        CsHist1                 (const std::string &name,const std::string &title,
                                                 uint32 Xbins,float64 Xmin,float64 Xmax);

  // ---------------------------------------------------------------------------
  // Operators
  // ---------------------------------------------------------------------------

  public:
    
    /// Copy operator.
    CsHist1            &operator =              (const CsHist1 &h);
    
  // ---------------------------------------------------------------------------
  // Methods
  // ---------------------------------------------------------------------------

  public:

    /// \return Histogram dimension. This is number 1.
    uint16              NDim                    (void) const { return 1; }

    /// Clear all bins contents.
    void                Clear                   (void);

    /// Fill histogram.  If the class hasn't been instantiated (i.e. if Fill()
    /// has been called on a NULL pointer), the Fill() request is ignored.
    void                Fill                    (float64 x, float64 weight=1);

    size_t              GetBin                  (float64 x1) const;

    T                   GetBinContent           (size_t bin) const;

    T                   GetBinError             (size_t bin) const;
    
    void                SetBinContent           (size_t bin,T v);

    void                SetBinError             (size_t bin,T v);

    float64             GetMean                 (void) const;

    float64             GetRMS                  (void) const;
    
    /// See CsHistBase::Fit()
    float64             Fit                     (FitFunction f,std::vector<CsFitParam> &param,
                                                 const std::string &options,float64 min,float64 max) const;
    
    /// See CsHistBase::Fit()
    float64             Fit                     (FitFunction f,std::vector<float64> &param,
                                                 const std::string &options,float64 min,float64 max) const
                                                {
                                                  std::vector<CsFitParam> par;

                                                  for( size_t i=0; i<param.size(); i++ )
                                                  {
                                                    char name[11];
                                                    sprintf(name,"par%zu",i+1);
                                                    par.push_back( CsFitParam(name,param[i]) );
                                                  }

                                                  float64 result = Fit(f,par,options,min,max);

                                                  for( size_t i=0; i<param.size(); i++ )
                                                    param[i] = par[i].GetFinish();

                                                  return result;
                                                }

    /// See CsHistBase::Fit()
    float64             Fit                     (FitFunction f,std::vector<float64> &param,
                                                 const std::string &options="") const
                                                {
                                                  const Dimension d=GetDim(1);
                                                  return Fit(f,param,options,d.GetMin(),d.GetMax());
                                                }

    /// See CsHistBase::Fit()
    float64             Fit                     (FitFunction f,std::vector<CsFitParam> &param,
                                                 const std::string &options="") const
                                                {
                                                  const Dimension d=GetDim(1);
                                                  return Fit(f,param,options,d.GetMin(),d.GetMax());
                                                }
};

////////////////////////////////////////////////////////////////////////////////

/*! \class CsHist2
    \author Alexander Zvyagin
    \brief 2-dimensional histogram.
*/
template <class T>
class CsHist2 : public CsHistBase
{
  // ---------------------------------------------------------------------------
  // Types.
  // ---------------------------------------------------------------------------

  public:
  
    typedef             T                       BinContentType;

  // ---------------------------------------------------------------------------
  // Constructors, destructor.
  // ---------------------------------------------------------------------------

  public:

    /// Destructor.
    virtual            ~CsHist2                 (void);

    /// Copy constructor.
                        CsHist2                 (const CsHist2 &h);

    /// Base constructor.
                        CsHist2                 (const std::string &name,const std::string &title,
                                                 uint32 Xbins,float64 Xmin,float64 Xmax,
                                                 uint32 Ybins,float64 Ymin,float64 Ymax);

  // ---------------------------------------------------------------------------
  // Operators
  // ---------------------------------------------------------------------------

  public:
    
    /// Copy operator.
    CsHist2            &operator =              (const CsHist2 &h);
    
  // ---------------------------------------------------------------------------
  // Methods
  // ---------------------------------------------------------------------------

  public:

    /// \return Histogram dimension. This is number 2.
    uint16              NDim                    (void) const { return 2; }

    /// Clear all bins contents.
    void                Clear                   (void);

    /// Fill histogram.  If the class hasn't been instantiated (i.e. if Fill()
    /// has been called on a NULL pointer), the Fill() request is ignored.
    void                Fill                    (float64 x, float64 y, float64 weight=1);

    size_t              GetBin                  (float64 x1,float64 x2) const;

    T                   GetBinContent           (size_t bin) const;

    T                   GetBinError             (size_t bin) const;
    
    void                SetBinContent           (size_t bin,T v);

    void                SetBinError             (size_t bin,T v);

    float64             GetMean                 (size_t axis) const;

    /*! \return Root mean square for the given axis (1,2).
        Does not work for HBOOK.
    */
    float64             GetRMS                  (size_t axis) const;
    
    /// See CsHistBase::Fit()
    float64             Fit                     (FitFunction f,std::vector<CsFitParam> &param,
                                                 const std::string &options,
                                                 float64 Xmin,float64 Xmax,
                                                 float64 Ymin,float64 Ymax) const;
    
    /// See CsHistBase::Fit()
    float64             Fit                     (FitFunction f,std::vector<float64> &param,
                                                 const std::string &options,
                                                 float64 Xmin,float64 Xmax,
                                                 float64 Ymin,float64 Ymax) const
                                                {
                                                  std::vector<CsFitParam> par;

                                                  for( size_t i=0; i<param.size(); i++ )
                                                  {
                                                    char name[11];
                                                    sprintf(name,"par%zu",i+1);
                                                    par.push_back( CsFitParam(name,param[i]) );
                                                  }

                                                  float64 result = Fit(f,par,options,Xmin,Xmax,Ymin,Ymax);

                                                  for( size_t i=0; i<param.size(); i++ )
                                                    param[i] = par[i].GetFinish();

                                                  return result;
                                                }

    /// See CsHistBase::Fit()
    float64             Fit                     (FitFunction f,std::vector<float64> &param,
                                                 const std::string &options="") const
                                                {
                                                  const Dimension d1=GetDim(1), d2=GetDim(2);
                                                  return Fit(f,param,options,d1.GetMin(),d1.GetMax(),d2.GetMin(),d2.GetMax());
                                                }

    /// See CsHistBase::Fit()
    float64             Fit                     (FitFunction f,std::vector<CsFitParam> &param,
                                                 const std::string &options="") const
                                                {
                                                  const Dimension d1=GetDim(1), d2=GetDim(2);
                                                  return Fit(f,param,options,d1.GetMin(),d1.GetMax(),d2.GetMin(),d2.GetMax());
                                                }
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/*!
    \example hist_test.cc
    This is example of how to use histograms package.
*/

typedef CsHist1<int16  > CsHist1S;
typedef CsHist2<int16  > CsHist2S;
typedef CsHist1<float32> CsHist1F;
typedef CsHist2<float32> CsHist2F;
typedef CsHist1<float64> CsHist1D;
typedef CsHist2<float64> CsHist2D;

////////////////////////////////////////////////////////////////////////////////

#endif // CsHist_h
