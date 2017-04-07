// $Id: CsDummyStore.h,v 1.1 2000/03/03 11:05:07 benigno Exp $
#ifndef CsDummyStore_h
#define CsDummyStore_h
/*!
   \file    CsDummyStore.h
   \brief   Dummy interface to switch off Objy and Date dependency...
   \author  Benigno Gobbo 
   \version $Revision: 1.1 $
   \date    $Date: 2000/03/03 11:05:07 $
*/

#include "CsStore.h"

class CsDummyStore : public CsStore  {

 public:

  static CsDummyStore* Instance() { 
    if(i_==0) i_= new CsDummyStore(); 
    return(i_); 
  }

  bool init()              { return(true); } 
  bool scan()              { return(true); } 
  bool next()              { return(true); }
  uint8* rawBuffer()       { return((uint8*)NULL); }
  uint32 getEventInRun()   const { return(0); }
  uint32 getRun()          const { return(0); }
  uint32 getEventInBurst() const { return(0); }
  uint32 getBurst()        const { return(0); }
  uint32 getTriggerMask()  const { return(0); }
  uint32 getErrorCode()    const { return(0); }
  CsTime getTime()         const { CsTime dummy; return(dummy); }

 protected:
  CsDummyStore() {}
  virtual ~CsDummyStore() { i_ = 0; }

 private:
  static CsDummyStore* i_;

};

CsDummyStore* CsDummyStore::i_ = 0;

#endif // CsDummyStore_h
