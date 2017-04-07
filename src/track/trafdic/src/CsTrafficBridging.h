// $Id: CsTrafficBridging.h 13148 2011-12-28 16:55:25Z kbicker $

/*!
   \file    CsTrafficBridging.h
   \author  Benigno.Gobbo@cern.ch
   \version $Revision: 13148 $
   \date    $Date: 2011-12-28 17:55:25 +0100 (Wed, 28 Dec 2011) $ 

*/

#ifndef CsTrafficBridging_h
#define CsTrafficBridging_h

#include "CsSTD.h"
#include "CsTrkBridging.h"

/*! \class CsTrafficBridging 
    Coral interface to Traffic segment bridging
*/

class CsTrafficBridging : public CsTrkBridging {

 public:

  CsTrafficBridging();           //!< constructor
  ~CsTrafficBridging();          //!< destructor

  /*!
    Interface to Traffic segment bridging
    \sa TEv::BridgeSegments()
  */
  bool doBridging( std::list<CsTrack*>& tracks, std::list<CsCluster*>& clusters );
};

#endif //CsTrafficBridging_h
