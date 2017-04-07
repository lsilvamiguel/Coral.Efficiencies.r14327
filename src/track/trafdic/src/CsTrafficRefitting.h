// $Id: CsTrafficRefitting.h 13148 2011-12-28 16:55:25Z kbicker $

/*!
   \file    CsTrafficRefitfing.h
   \brief Coral interface to TraFFiC's refitting tracks using vertexing info.
   \author  Yann.Bedfer@cern.ch
*/

#ifndef CsTrafficRefitting_h
#define CsTrafficRefitting_h

/*! \class CsTrafficRefitting 
    Coral interface to TraFFiC's refitting tracks using vertexing info.
*/

#include "CsTrack.h"
#include "CsVertex.h"

class CsTrafficRefitting {

 public:

  CsTrafficRefitting();           //!< constructor
  ~CsTrafficRefitting();          //!< destructor

  /*!
    Interface to TraFFiC's tracks refitting.
    \sa TEv::TracksRefit()
  */
  bool doRefitting(CsVertex *pVertex, std::list<CsTrack*> &tracksToBeDeleted);
};

#endif //CsTrafficRefitting_h
