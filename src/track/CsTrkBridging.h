// $Id: CsTrkBridging.h,v 1.4 2009/12/22 10:56:48 ybedfer Exp $

/*!
   \file    CsTrkBridging.h
   \brief   Compass Track Bridging Abstract Base Class.
   \author  Benigno Gobbo
   \version $Revision: 1.4 $
   \date    $Date: 2009/12/22 10:56:48 $ 

*/

#ifndef CsTrkBridging_h
#define CsTrkBridging_h

#include "CsTrack.h"


/*! \class CsTrkBridging 
    \brief Compass Track Bridging Abstract Base Class.

    This abstract class is intended to define the mandatory methods for
    all tracks bridging (joining of track segments from different zones) 
    classes. 
*/

class CsTrkBridging {

 public:
  // By making a base class destructor virtual, one ensures that the destructor of any derived class is executed (in addition to, and prior to, that of the base class).
  virtual ~CsTrkBridging() {}

  /*! \fn bool virtual doBridging( list<CsTrack*>& tracks, list<CsCluster*>& clusters )
      \brief This method performs the bridging of track segments belonging
      from two different COMPASS apparatus zones. It returns \c true if
      operations ended correctly, \c false otherwise.
      \warning Track segments are removed from lists if they are used
      for the creation of the new track. New tracks are introduced as
      bridged segments.
      \param tracks List of tracks to be bridged
      \param clusters List of unused digits in the zone
  */
  virtual bool doBridging( std::list<CsTrack*>& tracks, 
			   std::list<CsCluster*>& clusters ) = 0; 

};

#endif //CsTrkBridging_h
