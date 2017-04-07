// $Id: CsDateStore.h,v 1.14 2010/01/20 23:09:30 tnagel Exp $

#ifndef CsDateStore_h
#define CsDateStore_h

/*!
   \file    CsDateStore.h
   \brief   Engine to access event data from DATE file or DAQ stream
   \author  Sergei.Gerassimov@cern.ch
   \version $Revision: 1.14 $
   \date    $Date: 2010/01/20 23:09:30 $

*/

#include "coral_config.h"

#include "CsTime.h"
#include "CsStore.h"
#include <string>
#include <list>

/*! \class CsDateStore 
    \brief Engine to access event data from DATE file or DAQ stream

    It provides methods to open DATE file or stream
    and retrieve events.
*/

class CsDateStore: public CsStore {

public:

  static CsDateStore* Instance();

  bool init();         //!< Init call

  bool scan();         //!< Prepare the iteration
  bool next();         //!< Next event

  uint8* rawBuffer(); //!< Makes the raw data buffer available

  uint32 getEventInRun()      const ;
  uint32 getRun()             const ;
  uint32 getEventInBurst()    const ;
  uint32 getBurst()           const ;
  uint32 getTriggerMask()     const ;
  uint32 getErrorCode()       const ;
  CsTime getTime()            const ;

 protected:

  CsDateStore();
  virtual ~CsDateStore();

private:

  static CsDateStore* instance_;

  bool init_;
  bool scan_;

  // Added: 20011001, Benigno
  std::list<std::string*> _files;
  std::list<std::string*>::iterator _currentFile;
  // end

  uint8* date_buffer_ptr;           // pointer to the DATE buffer
  // Added 20050109, Benigno
  int   _fileid;
  unsigned int _dateVersionGuess;

  unsigned int _getEvent( void *ptr );  

};

#endif // CsDateStore_h











