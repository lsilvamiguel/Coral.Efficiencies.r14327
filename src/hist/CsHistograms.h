/*!
   \file    CsHistograms.h
   \brief   Compass histograms package.
   \version $Revision: 1.11 $
   \author  Alexander Zvyagin
   \date    $Date: 2003/04/09 12:27:26 $
*/

#ifndef CsHistograms_h
#define CsHistograms_h

#include "coral_config.h"

#include <vector>
#include <string>
#include "CsEndOfJob.h"

/*! \class CsHistograms
   \author Alexander Zvyagin
   \brief Interface to Compass histograms package.

   To use histogram package you have to add two lines in a CORAL options file.
   \arg a) histograms package ROOT
   \arg b) histograms home    file.root
   
   All histograms will be based on ROOT package and will be stored
   in file with name "file.root" 

   Any method of CsHistograms class may be accessed by two ways:
   \arg 1. CsHistograms::Instance()->method_name();
   \arg 2. CsHistograms::method_name();
*/
class CsHistograms : public CsEndOfJob
{
  public:

    /*! \brief Types of histograms package implementations.
    
       In CORAL we will have single interface to histograms
       (see CsHistBase,CsHist1,CsHist2) and different implementations of
       these histograms methods.
       For the time beeing we have two of them:
       \arg DUMMY - empty bodies of all (almost all) methods.
       \arg ROOT  - use ROOT (http://root.cern.ch).
       \arg HBOOK - use HBOOK (http://wwwinfo.cern.ch/asdoc/hbook_html3/hboomain.html).
    */
    enum Implementation {DUMMY=0, ROOT=1, HBOOK=2};

    /* \brief \return Singleton for CsHistograms.
      This static method returns pointer to singleton object of class CsHistograms.
    */
    static CsHistograms*Instance                (void);
    
    /// \return Histogram package implementation.
    static Implementation GetImplementation     (void) { return implementation; }
    
    /*! \return Home path of all histograms.
        This may be either file name or path inside CORAL database.
    */
    static const std::string&GetHomePath             (void) { return home_path; }
    
    static std::string       GetCurrentPath          (bool fast=true);

    /*! \brief Use this function to change current path. If \c path does not exist, it will be created.
        Do nothing if path="",".".
    */
    static void         SetCurrentPath          (const std::string &path,bool fast=true);

    /// \return Was the package activated?
    static bool         IsActivated             (void) {return is_activated;}
    
    /// Get amount of consumed memory by the histogram's package.
    static unsigned     GetConsumedMemory       (void) {return memory_consumed;}

    /// Set amount of consumed memory by the histogram's package.
    static void         SetConsumedMemory       (unsigned int m) {memory_consumed=m;}

    void*               getFile                 (void) {return home;}

  public:

    /// Destructor.
    virtual            ~CsHistograms            (void);

  protected:
  
    /// Protected default constructor.
                        CsHistograms            (void);

    /// Protected copy constructor.
                        CsHistograms            (const CsHistograms &h);

    /// Protected assignment operator.
    CsHistograms       &operator =              (const CsHistograms &h);

  private:

    /// Initialize the package.
    static void         Init                    (void);
    
    /// This method will be called at the very end.
    bool                end                     (void);
    
    /// Pointer to CsHistograms singleton object.
    static CsHistograms*instance;
    
    /// Place where histograms will be stored.
    static std::string       home_path;

    /// Pointer to an object that stores histograms.
    static void        *home;
    
    /// \return Was the package activated?
    static bool         is_activated;

    /// Type of implementation
    static Implementation implementation;
    
    /// Package-dependent top directory.
    static std::string       top_directory;
    
    /// Amount of memory consumed by histogram's package.
    static unsigned int memory_consumed;
    
  friend class CsInit;
};

/*!
    \example coral.options.hist
    This is example of how to use histograms package.
*/

#include "CsHist.h"
#include "CsNtuple.h"

#endif //CsHistograms_h
