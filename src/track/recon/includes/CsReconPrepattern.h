// $Id: CsReconPrepattern.h,v 1.6 2003/04/17 12:38:32 benigno Exp $

/*!
   \file    CsReconPrepattern.h
   \brief   Recon track prepattern class.
   \author  ...
   \version $Revision: 1.6 $
   \date    $Date: 2003/04/17 12:38:32 $

*/

#ifndef CsReconPrepattern_h
#define CsReconPrepattern_h

#include "CsSTD.h"
#include "CsTrkPrepattern.h"

/*! \class CsReconPrepattern 
    \brief Recon of track prepattern class.

    use this class as starting point for your development.
    
*/
class CsReconPrepattern : public CsTrkPrepattern {

 public:
 
  /*! \fn CsReconPrepattern()
    \brief ...
  */
  CsReconPrepattern();

  /*! \fn ~CsReconPrepattern()
    \brief ...
  */
  virtual ~CsReconPrepattern();

  /*! \fn bool doPrepattern( const std::list<CsCluster*> clusters, const std::list<CsZone*> zones )
    \brief ... 
  */
  bool doPrepattern( const std::list<CsCluster*> &clusters, const std::list<CsZone*> &zones );

  /*! \fn bool getPatterns( std::list<CsTrack*>& tracks );
    \brief ... 
  */
  bool getPatterns( std::list<CsTrack*>& tracks );

  /*! \fn std::list<CsCluster*> getUnusedClusters()
    \brief ... 
  */
  std::list<CsCluster*> getUnusedClusters() const;

 private:



};

#endif //CsReconPrepattern_h
