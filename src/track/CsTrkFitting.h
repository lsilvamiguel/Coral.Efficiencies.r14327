// $Id: CsTrkFitting.h,v 1.4 2009/12/22 10:56:48 ybedfer Exp $

/*!
   \file    CsTrkFitting.h
   \brief   Compass Track Parameters Fitting Abstract Base Class.
   \author  Benigno Gobbo
   \version $Revision: 1.4 $
   \date    $Date: 2009/12/22 10:56:48 $ 

*/

#ifndef CsTrkFitting_h
#define CsTrkFitting_h

#include "CsTrack.h"

/*! \class CsTrkFitting 
    \brief Compass Track Parameters Fitting Abstract Base Class.

    This abstract class is intended to define the mandatory methods for
    all track fitting classes. 
*/

class CsTrkFitting {

 public:

  // By making a base class destructor virtual, one ensures that the destructor of any derived class is executed (in addition to, and prior to, that of the base class).
  virtual ~CsTrkFitting() {}

  /*! \fn virtual bool doFitting( list<CsTrack*>& tracks, const list<CsCluster*>& clusters )
      \brief This method performs the track parameters refit of the
      given list of tracks. 
      \warning this method modify the \c CsTrack objects!
      \param tracks The list of tracks to be refitted
      \param clusters The list of unused cluster in the zone
   */ 
  virtual bool doFitting( std::list<CsTrack*>& tracks,
			  const std::list<CsCluster*>& clusters ) = 0;

};

#endif //CsTrkFitting_h
