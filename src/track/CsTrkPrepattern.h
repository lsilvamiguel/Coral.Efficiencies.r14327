// $Id: CsTrkPrepattern.h,v 1.5 2009/12/22 10:56:48 ybedfer Exp $

/*!
   \file    CsTrkPrepattern.h
   \brief   Compass Track Prepattern Abstract Base Class.
   \author  Benigno Gobbo
   \version $Revision: 1.5 $
   \date    $Date: 2009/12/22 10:56:48 $ 

*/

#ifndef CsTrkPrepattern_h
#define CsTrkPrepattern_h



/*! \class CsTrkPrepattern 
    \brief Compass Track Prepattern Abstract Base Class.

    This abstract class is intended to define the mandatory methods for
    all Pre pattern recognition classes. 
*/

#include "CsCluster.h"
#include "CsZone.h"
#include "CsHelix.h"
#include "CsTrack.h"

class CsTrkPrepattern {

 public:

  // By making a base class destructor virtual, one ensures that the destructor of any derived class is executed (in addition to, and prior to, that of the base class).
  virtual ~CsTrkPrepattern() {}

  /*! \fn bool doPrepattern( const list<CsCluster*> clusters, const list<CsZone*> zones )
      \brief This method performs the patter recognition on a given
      collection of clusters. Returns \c true if the operation ended correctly;
      returns \c false otherwise.
      \param clusters List of pointers to the clusters to be used in the 
      pattern recognition procedure.
      \param zones List of zones on which perform the procedure
   */
  virtual bool doPrepattern( const std::list<CsCluster*> &clusters, 
			     const std::list<CsZone*>    &zones )  = 0;

  /*! \fn bool getPatterns( list<CsTrack*>& tracks )
      \brief This method returns the list of found track candidates after
      the prepattern procedure. These tracks could be a not completely 
      fitted ones.
      \param The list of track where to add the found patterns
   */
  virtual bool getPatterns( std::list<CsTrack*>& tracks ) = 0;

  /*! \fn list<CsCluster*> getUnusedClusters()
      \brief This method returns the list of clusters not associated to
      any found pattern.
   */
  virtual std::list<CsCluster*> getUnusedClusters() const = 0;

};

#endif //CsTrkPrepattern_h
