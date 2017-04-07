// $Id: CsReconFitting.h,v 1.3 2003/04/17 12:38:32 benigno Exp $

/*!
   \file    CsReconFitting.h
   \brief   Recon track fitting class
   \author  ...
   \version $Revision: 1.3 $
   \date    $Date: 2003/04/17 12:38:32 $

*/

#ifndef CsReconFitting_h
#define CsReconFitting_h

#include "CsSTD.h"
#include "CsTrkFitting.h"

/*! \class CsReconFitting 
    \brief Recon of track fitting class

    Use this file as starting point for you class development.
*/

class CsReconFitting : public CsTrkFitting {

 public:

  /*! \fn CsReconFitting()
    \brief ...
  */
  CsReconFitting();

  /*! \fn ~CsReconFitting()
    \brief ...
  */
  virtual ~CsReconFitting();

  /*! \fn bool doFitting( std::list<CsTrack*>& tracks, std::const list<CsCluster*>& clusters )
    \brief ... 
  */
  bool doFitting( std::list<CsTrack*>& tracks, const std::list<CsCluster*>& clusters );

 private:



};

#endif //CsReconFitting_h
