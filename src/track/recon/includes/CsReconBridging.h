// $Id: CsReconBridging.h,v 1.3 2003/04/17 12:38:32 benigno Exp $

/*!
   \file    CsReconBridging.h
   \brief   Recon track bridging class
   \author  ...
   \version $Revision: 1.3 $
   \date    $Date: 

*/

#ifndef CsReconBridging_h
#define CsReconBridging_h

#include "CsSTD.h"
#include "CsTrkBridging.h"

/*! \class CsReconBridging 
    \brief Recon of track bridging class

    Use this to develop your classes
*/

class CsReconBridging : public CsTrkBridging {

 public:

  /*! \fn CsReconBridging()
    \brief ...
  */
  CsReconBridging();

  /*! \fn ~CsReconBridging()
    \brief ...
  */
  virtual ~CsReconBridging();

  /*! \fn bool doBridging( std::list<CsTrack*>& tracks, std::list<CsCluster*>& clusters )
    \brief ... 
  */
  bool doBridging( std::list<CsTrack*>& tracks, std::list<CsCluster*>& clusters );

 private:



};

#endif //CsReconBridging_h
