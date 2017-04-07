/*!
   \file    CsVrtPattern.h
   \brief   Compass Vertex Pattern Abstract Base Class.
   \author  Alexandre Korzenev
   \version $Revision: 1.10 $
   \date    $Date: 2010/04/20 00:13:15 $ 

*/

#ifndef CsVrtPattern_h
#define CsVrtPattern_h



/*! \class CsVrtPattern 
    \brief Compass Vertex Pattern Abstract Base Class.

    This abstract class is intended to define the mandatory methods for
    all Pattern recognition classes. 
*/

#include "CsVertex.h"
#include "CsTrack.h"
#include "CsParticle.h"

class CsVrtPattern {

 public:
    // By making a base class destructor virtual, one ensures that the destructor of any derived class is executed (in addition to, and prior to, that of the base class).
  virtual ~CsVrtPattern() {}

  /*! \fn bool doPattern( vector<CsParticle*> &particles ) = 0;
    \brief This method performs the patter recognition on a given
    collection of particles. Returns \c true if the operation ended correctly;
    returns \c false otherwise.
    \param particles = Vector of pointers to the particles to be used in the 
    pattern recognition procedure.
    \param reTrackT0 = if non null on input, allows "doPattern" to request a re-tracking, by returning on output the event time to be used in the re-tracking.
  */
  virtual bool doPattern(std::vector<CsParticle*> &particles, double *reTrackT0 = 0) = 0;

  /*! \fn bool getPatterns( std::list<CsVertex*>& vrts, map<CsTrack*,bool>& specials ) = 0;
    \brief This method returns the list of found vertecies candidates after
    the pattern procedure. 
    \param vrts = The list of vertex candidates.
    \param specials = map of marked tracks (For Primary Vertex Only).
  */
  virtual bool getPatterns(std::list<CsVertex*> &vrts, std::map<CsTrack*,bool> &specials ) = 0;

  /*! \fn virtual const CsVertex *getT0SettingVertex() const
    \brief Returns a pointer to the vertex selected by the PR as the potential Best Primary Vertex and hence assigning the T0 used in a possible refit of all tracks. 
   */
  virtual const CsVertex *getT0SettingVertex() const = 0;
};

#endif //CsVrtPattern_h
