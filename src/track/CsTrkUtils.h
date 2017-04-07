// $Id: CsTrkUtils.h,v 1.2 2009/12/22 10:56:48 ybedfer Exp $

/*!
   \file    CsTrkUtils.h
   \brief   Compass Tracking Utils Abstract Base Class.
   \author  Y.B.
   \version $Revision: 1.2 $
   \date    $Date: 2009/12/22 10:56:48 $ 

*/

#ifndef CsTrkUtils_h
#define CsTrkUtils_h



/*! \class CsTrkUtils 
    \brief Compass Tracking Utils Abstract Base Class.

    This abstract class was designed to provide for making Dico. But could be
    mad an all purpose simulation abstract base class.
*/

#include "CsCluster.h"
#include "CsZone.h"
#include "CsHelix.h"
#include "CsTrack.h"

class CsTrkUtils {

 public:
  // By making a base class destructor virtual, one ensures that the destructor of any derived class is executed (in addition to, and prior to, that of the base class).
  virtual ~CsTrkUtils() {}

  virtual bool genLattice() =0;

};

#endif //CsTrkUtils_h
