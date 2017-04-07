// $Id: CsSampleBridging.h,v 1.2 1999/12/08 16:16:11 benigno Exp $

/*!
   \file    CsSampleBridging.h
   \brief   Sample of track bridging class
   \author  Benigno Gobbo
   \version $Revision: 1.2 $
   \date    $Date: 1999/12/08 16:16:11 $

*/

#ifndef CsSampleBridging_h
#define CsSampleBridging_h

#include "CsSTD.h"
#include "CsTrkBridging.h"

/*! \class CsSampleBridging 
    \brief Sample of track bridging class

    Use this to develop your classes
*/

class CsSampleBridging : public CsTrkBridging {

 public:

  /*! \fn CsSampleBridging()
    \brief ...
  */
  CsSampleBridging();

  /*! \fn ~CsSampleBridging()
    \brief ...
  */
  virtual ~CsSampleBridging();

  /*! \fn bool doBridging( list<CsTrack>& tracks, list<CsCluster*>& clusters )
    \brief ... 
  */
  bool doBridging( list<CsTrack*>& tracks, list<CsCluster*>& clusters );

 private:



};

#endif //CsSampleBridging_h
