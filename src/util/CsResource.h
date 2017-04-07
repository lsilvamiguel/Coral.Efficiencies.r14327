// $Id: CsResource.h,v 1.4 2007/02/05 10:18:53 gobbo Exp $ 

/*!
  \file CsResource.h
  \brief Debugging object to read and store coral voracity
  \author Hugo Pereira
  \version $Revision: 1.4 $
  \date $Date: 2007/02/05 10:18:53 $
*/

#ifndef CsResource_h
#define CsResource_h

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <string>

class CsResource {

  public:
  static CsResource* Instance( void );     //!< Singleton instanciation
  void dumpRCToLog( std::string tag );     //!< Dump Resource To Log File
  void dumpRCDiffToLog( std::string tag ); //!< Dump Resource diff To Log File
  
 
  protected:
  CsResource( void );           //!< The Constructor
  virtual ~CsResource( void );  //!< The Destructor
  
  class CsRCData{  
    public:
    CsRCData( void );
    CsRCData& operator=( const CsRCData& );    
    bool update( void );

    int success_;
    int vmSize_;
    int vmLock_;
    int vmResident_;
    int vmData_;
    int vmStack_;
    int vmExec_;
    int vmLib_;
  };
    
  private:
  static CsResource* instance_;
  std::string logFName_;
  CsRCData* OldRC_;
  CsRCData* NewRC_;
  
  
};

#endif
