/*! 
  \brief Pair of track segments
  
  Simple class for use in different kinds of track bridging
  procedures.

*/

#include "TTrack.h"

class TTrackPair {
public:

  static unsigned int NobjCreated;
  static unsigned int NobjDestructed;

  TTrackPair() 
  {
    NobjCreated++;
  };

  ~TTrackPair()
  {
    NobjDestructed++;
  };

  TTrackPair(const TTrackPair& tp):
    Chi2(tp.Chi2),
    iL  (tp.iL),
    iR  (tp.iR)
  {
     NobjCreated++;
  };


  double Chi2;                     //!< track candidate Chi2
  std::list<TTrack>::iterator iL;  //!< "upstream"   track segment iterator
  std::list<TTrack>::iterator iR;  //!< "downstream" track segment iterator
  //! "Less" operator
  bool operator < (const TTrackPair& tp) const { return (Chi2 < tp.Chi2);} 
};

