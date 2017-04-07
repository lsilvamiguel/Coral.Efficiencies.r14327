// $Id: CsStore.h,v 1.6 2009/12/22 10:58:04 ybedfer Exp $

#ifndef CsStore_h
#define CsStore_h

/*!
   \file    CsStore.h
   \brief   Interfaces (abstract class) to access event data from 
            Objectivity/DB or DATE stream.
   \author  Massimo Lamanna, Sergei Gerassimov
   \version $Revision: 1.6 $
   \date    26.01.2000

*/

#include <string>
#include "CsTypes.h"
#include "CsTime.h"

class CsRecoEvent;

class CsStore  {

public:
  virtual ~CsStore() {}

  virtual bool init()= 0;         //!< Init call

  virtual bool scan() = 0;         //!< Prepare the iteration
  virtual bool next() = 0;         //!< Next event

  virtual uint8* rawBuffer() = 0;  //!< Makes the raw data buffer available

  virtual uint32 getEventInRun          () const = 0;
  virtual uint32 getRun                 () const = 0;
  virtual uint32 getEventInBurst        () const = 0;
  virtual uint32 getBurst               () const = 0;
  virtual uint32 getTriggerMask         () const = 0;
  virtual uint32 getErrorCode           () const = 0;

  virtual CsTime getTime                () const = 0;

  virtual uint8* dstBuffer()                                       { return NULL; }
  virtual bool   dst(int slot)                                     { return false; }
  virtual bool   uploadDST(CsRecoEvent*, int fatness, int version) { return false; }
  virtual bool   downloadDST(CsRecoEvent*, int version)            { return false; }
  virtual void   saveAndExit()                                     { ; }
  virtual void   abortandExit()                                    { ; }
  virtual bool   end()                                             { return true; }
  virtual int    getRawBufferLength()                              { return 0; }
  virtual bool   uploadTBNamesString(std::string&)                      { return false; }
  virtual bool   downloadTBNamesString(std::string&)                    { return false; }
  virtual int    getSelectedSlot()                                 { return 0; }
};

#endif // CsStore_h
