// $Id: CsSampleFitting.h,v 1.2 1999/12/08 16:16:12 benigno Exp $

/*!
   \file    CsSampleFitting.h
   \brief   Sample of track fitting class
   \author  Benigno Gobbo
   \version $Revision: 1.2 $
   \date    $Date: 1999/12/08 16:16:12 $

*/


#ifndef CsSampleFitting_h
#define CsSampleFitting_h

#include "CsSTD.h"
#include "CsTrkFitting.h"

/*! \class CsSampleFitting 
    \brief Sample of track fitting class

    Use this file as starting point for you class development.
*/

class CsSampleFitting : public CsTrkFitting {

 public:

  /*! \fn CsSampleFitting()
    \brief ...
  */
  CsSampleFitting();

  /*! \fn ~CsSampleFitting()
    \brief ...
  */
  virtual ~CsSampleFitting();

  /*! \fn bool doFitting( list<CsTrack>& tracks, const list<CsCluster*>& clusters )
    \brief ... 
  */
  bool doFitting( list<CsTrack*>& tracks, const list<CsCluster*>& clusters );

 private:



};

#endif //CsSampleFitting_h
