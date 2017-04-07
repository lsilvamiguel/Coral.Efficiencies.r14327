/*!
   \file    CsVrtFitting.h
   \brief   Compass Vertex Parameters Fitting Abstract Base Class.
   \author  Alexandre Korzenev
   \version $Revision: 1.8 $
   \date    $Date: 2010/05/24 22:10:56 $ 

*/

#ifndef CsVrtFitting_h
#define CsVrtFitting_h

#include "CsVertex.h"
#include "CsTrack.h"

/*! \class CsVrtFitting 
    \brief Compass Vertex Parameters Fitting Abstract Base Class.

    This abstract class is intended to define the mandatory methods for
    all vertex fitting classes. 
*/

class CsVrtFitting {

 public:
  // By making a base class destructor virtual, one ensures that the destructor of any derived class is executed (in addition to, and prior to, that of the base class).
  virtual ~CsVrtFitting() {}

  /*! \fn virtual bool doFitting(list<CsVertex*> &vrts, map<CsTrack*,bool> &specials,bool reTrackingON, double &T0) = 0
    \brief Performs the vertex/track parameters fit for the given lists of vertex/tracks. 
    \param \e vrts = Reference to list of CsVertex objects to be fitted.
    \param \e specials = Reference to map of scattered-mu ID.
    \param \e reTrackingON: true if expecting to be in the context of a 2nd pass of vertexing, after a re-tracking has taken place.
    \param \e T0: If not a null pointer, on input it hols the T0 used in tracking: if this turns out to deviate too much from eventual best vertex' time, and \e reTrackingON is non true, then doFitting aborts and returns false so that a re-tracking be requested. On output, \e T0 holds the Best Vertex time, which is the recommended reference time to use as T0 in the re-tracking.  
  */ 
  virtual bool doFitting(std::list<CsVertex*> &vrts, std::map<CsTrack*,bool> &specials,
			 bool reTrackingON = false, double *T0 = 0) = 0;

  /*! \fn virtual const CsVertex *getBestPVertex() const
    \brief Returns a point to the Best Primary Vertex, if defined.
   */
  virtual const CsVertex *getBestPVertex() const = 0;
};

#endif //CsVrtFitting_h
