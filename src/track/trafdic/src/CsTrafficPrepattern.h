// $Id: CsTrafficPrepattern.h 13148 2011-12-28 16:55:25Z kbicker $

/*!
   \file    CsTrafficPrepattern.h
   \author  Benigno.Gobbo@cen.ch
   \version $Revision: 13148 $
   \date    $Date: 2011-12-28 17:55:25 +0100 (Wed, 28 Dec 2011) $ 

*/

#ifndef CsTrafficPrepattern_h
#define CsTrafficPrepattern_h

#include "CsSTD.h"
#include "CsTrkPrepattern.h"

/*! \class CsTrafficPrepattern 
    Coral interface to Traffic Pre-Pattern
*/

class CsTrafficPrepattern : public CsTrkPrepattern {

 public:

  CsTrafficPrepattern();          //!< constructor
  ~CsTrafficPrepattern();         //!< destructor

  /*! 
     Interface to Traffic pre-pattern (track pieces finding)
     \sa TEv::PrePattern()
  */
  bool doPrepattern( const std::list<CsCluster*> &clusters, const std::list<CsZone*> &zones );

  /*!
     Retutns to CORAL found track pieces
  */
  bool getPatterns( std::list<CsTrack*>& tracks );

  /*! 
    Returns to CORAL clusters, not assosiated with found track pieces 
  */
  std::list<CsCluster*> getUnusedClusters(void) const;

};

#endif //CsTrafficPrepattern_h








