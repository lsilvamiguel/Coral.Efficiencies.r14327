/*!
   \file    CsNtuple.h
   \brief   Compass Ntuple Class.
   \version $Revision: 1.8 $
   \author  Alexander Zvyagin
   \date    $Date: 2003/04/09 12:27:26 $
*/

#ifndef CsNtuple_h
#define CsNtuple_h

#include "coral_config.h"

# include <string>
#include "CsTypes.h"
#include "CsHome.h"

/*! \class CsNtuple
    \author Alexander Zvyagin
    \brief Ntuple Class.
    
    \example hist_test.cc
    
    In a Ntuple you can store many  (up to 1000) float32 variables per event.
    Only very simple Ntuples are supported. The format is "x1:x2:x3" and so on.
    All variables must have type float32.
*/

class CsNtuple : public CsHome
{
  // ---------------------------------------------------------------------------
  // Constructors, destructor.
  // ---------------------------------------------------------------------------

  public:

    /// Destructor.
    virtual            ~CsNtuple                (void);

  private:

    /// Copy constructor.    
                        CsNtuple                (const CsNtuple &h);

  public:

    /*! \brief Base constructor.

        Example of format: "x:y:z:energy". \b ATTENTION! This format does NOT depend on histograms package
        you are using. It is the same both for HBOOK and ROOT.
    */
                        CsNtuple                (const std::string &name,const std::string &title,const std::string &format);

  // ---------------------------------------------------------------------------
  // Operators
  // ---------------------------------------------------------------------------

  private:

    /// Assignment operator.
    CsNtuple           &operator =              (const CsNtuple &h);


  // ---------------------------------------------------------------------------
  // Methods
  // ---------------------------------------------------------------------------

  public:

    /// Fill ntuple.
    int32               Fill                    (const float32 *x=NULL);

    /// Fill ntuple.
    int32               Fill                    (float32 x1,float32 x2=0,float32 x3=0,float32 x4=0,
                                                 float32 x5=0,float32 x6=0,float32 x7=0,float32 x8=0,float32 x9=0);

    /// \return Ntuple name.
    virtual std::string      GetName                 (void) const;
    
    /// \return Ntuple title.
    virtual std::string      GetTitle                (void) const;
    
    /// \return path to an the ntuple.
    virtual std::string      GetPath                 (void) const;
    
    /// Set ntuple name.
    virtual void        SetName                 (const std::string &name);

    /// Set ntuple title.
    virtual void        SetTitle                (const std::string &title);
    
    virtual void        SetPath                 (const std::string &path);

  // ---------------------------------------------------------------------------
  // Attributes
  // ---------------------------------------------------------------------------

  private:
  
    /// Pointer to real NTuple object.
    void               *nt;

    /// Address of space with data for Fill().
    void               *address;
};

#endif
