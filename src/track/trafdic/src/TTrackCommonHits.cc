
//
// Method of class TTrack to 
// count number of common hits 
// of "this" track with other track
//

#include "TTrack.h"

/*
// old slow method
int TTrack::CommonHits(const TTrack& t)
{
  int ncommon=0;
  
  for(std::list<int>::iterator i = lHitPat.begin(); i != lHitPat.end(); i++){
    if((*i) < 0) continue;
    for(std::list<int>::const_iterator j = t.lHitPat.begin(); j != t.lHitPat.end(); j++){
      if((*j) < 0) continue;
      if((*i) == (*j)) ncommon++;
    }
  }
  return(ncommon);
}
*/

// new method with use of sIHit sets
// (not yet done)
int TTrack::CommonHits(const TTrack& t)
{
  int ncommon=0;
  
  for(std::list<int>::iterator i = lHitPat.begin(); i != lHitPat.end(); i++){
    if((*i) < 0) continue;
    for(std::list<int>::const_iterator j = t.lHitPat.begin(); j != t.lHitPat.end(); j++){
      if((*j) < 0) continue;
      if((*i) == (*j)) ncommon++;
    }
  }
  return(ncommon);
}
