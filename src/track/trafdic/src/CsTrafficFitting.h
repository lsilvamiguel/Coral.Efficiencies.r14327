// $Id: CsTrafficFitting.h 13148 2011-12-28 16:55:25Z kbicker $

/*!
   \file    CsTrafficFitting.h
   \author  Benigno.Gobbo@cern.ch
   \version $Revision: 13148 $
   \date    $Date: 2011-12-28 17:55:25 +0100 (Wed, 28 Dec 2011) $ 

*/

#ifndef CsTrafficFitting_h
#define CsTrafficFitting_h

#include "CsSTD.h"
#include "CsTrkFitting.h"

/*! \class CsTrafficFitting 
  Coral interface to Traffic track fit
*/

class CsTrafficFitting : public CsTrkFitting {

 public:

  CsTrafficFitting();          //!< constructor
  ~CsTrafficFitting();         //!< destructor

  /*! 
    Interface to Traffic track fit
    \sa TEv::TracksFit()
  */
  bool doFitting( std::list<CsTrack*>& tracks, const std::list<CsCluster*>& clusters );
};

#endif //CsTrafficFitting_h
