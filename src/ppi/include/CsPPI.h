/*!
  \file      CsPPI.h
  \brief     Compass PreProcessor Interface Class
  \author    Markus Kraemer
  \version   $Revision: 1.1 $
  \date      $Date: 2010/02/12 14:50:23 $

  \par       History:
  20091127   Design start of this class
*/

#ifndef __CsPPI_H__
#define __CsPPI_H__

#include <string>
#include <time.h>
#include <vector>
 #include "CsEvent.h"
 #include "DaqDataDecoding/DetID.h"

// using namespace CS;


/************************************************************
 ** CsPPImsg provides message handling interface for CsPPI **
 ************************************************************/
enum error_level {
  DEBUG3 = -3,
  DEBUG2 = -2, 
  DEBUG1 = -1,
  INFO = 0,
  WARNING = 1,
  ERROR = 2,
  FATAL = 3
};

class CsPPImsg {

 protected:

  error_level     fLevel;
  int             fErrCode;
  std::string     fWhat;
  time_t          fWhen;
  std::string     fFile;
  unsigned int    fLine;

 public:

  // Con-/Destruvtor 
  CsPPImsg(const error_level Level, const int Code,
           const std::string& what, 
           const std::string& file,  // use '__FILE__' cpp macro
           const unsigned int line); // use '__LINE__' cpp macro
  ~CsPPImsg();

  // Member access
  error_level             Level() const;
  int                     errcode() const;
  std::string             what() const;
  
};

/*************************************************************
 ** The CsPP class space provides an interface for coral,   **
 ** which provides preprocessing capability. Eache class    **
 ** providing preprocessing capability should inherit from  **
 ** this class.                                             **
 *************************************************************/

class CsPPI {

 protected:

  CS::DetID               fId;
  std::string             fTBName;
  std::vector<CsPPImsg*>  fMSG;
  unsigned int            fNBwarning;
  unsigned int            fNBerror;
  unsigned int            fNBfatal;
  error_level             fMsgThr;

 public:

  // Con-/Destuctors
                                        CsPPI(const CS::DetID &id, 
                                              const std::string &TBname,
                                              const error_level thr = INFO);
                                        CsPPI();
  virtual                               ~CsPPI();

  // Access to member objects
  const CS::DetID&                      GetID() const;
  const char*                           getName() const;     // from DetID
  const std::string&                    GetTBName() const;
  const std::vector<CsPPImsg*>&         GetMSG() const;
  unsigned int                          GetNBfatal() const;
  unsigned int                          GetNBerror() const;
  unsigned int                          GetNBwarning() const;
  const CsPPImsg*                       GetLastError(int code = 0) const;
  const CsPPImsg*                       GetLastFatal() const;
  const CsPPImsg*                       GetLastMSG(int thr = 0) const;
  void                                  ClearMSG(int thr = 0);
  void                                  MSG(CsPPImsg *MSG);
  void                                  MSG(const error_level level, 
                                            const int Code,
                                            const std::string& what, 
                                            const std::string& file,
                                            const unsigned int line);

  // function, which is called during preprocessing 
  // to do the work returns error code
  // to be overloaded by children
  virtual int                           process(); 
  virtual int                           end(); 

};

#endif
