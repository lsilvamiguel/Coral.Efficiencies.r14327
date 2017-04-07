// $Id: CsSamplePrepattern.h,v 1.3 2000/10/30 23:27:07 zvyagin Exp $

/*!
   \file    CsSamplePrepattern.h
   \brief   Sample of track prepattern class.
   \author  Benigno Gobbo
   \version $Revision: 1.3 $
   \date    $Date: 2000/10/30 23:27:07 $

*/

#ifndef CsSamplePrepattern_h
#define CsSamplePrepattern_h

#include "CsSTD.h"
#include "CsTrkPrepattern.h"

/*! \class CsSamplePrepattern 
    \brief Sample of track prepattern class.

    use this class as starting point for your development.
    
*/

class CsSamplePrepattern : public CsTrkPrepattern {

 public:

  /*! \fn CsSamplePrepattern()
    \brief ...
  */
  CsSamplePrepattern();

  /*! \fn ~CsSamplePrepattern()
    \brief ...
  */
  virtual ~CsSamplePrepattern();

  /*! \fn bool doPrepattern( const list<CsCluster*> clusters, const list<CsZone*> zones )
    \brief ... 
  */
  bool doPrepattern( const list<CsCluster*> &clusters, const list<CsZone*> &zones );

  /*! \fn list<CsTrack> getPatterns()
    \brief ... 
  */
  bool getPatterns( list<CsTrack*>& tracks );

  /*! \fn list<CsCluster*> getUnusedClusters()
    \brief ... 
  */
  const list<CsCluster*> &getUnusedClusters(void);

 private:



};

#endif //CsSamplePrepattern_h
