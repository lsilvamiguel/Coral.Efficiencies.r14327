// $Id: TTrackInsertMissedPlanes.cc 13148 2011-12-28 16:55:25Z kbicker $

#include "TTrack.h"

using namespace std;

/*!
  Insert missed planes references into track's list.

  Corresponding hit references are also inserted, w/ a default value of -1
  -1 ("no hit found") will possibly be replaced by
  -2 ("out of active area") later in TTrack::Refine()
*/

void TTrack::InsertMissedPlanes()
{
  list<int>::iterator ip1 = lPlnRef.begin(), ip2;
  list<int>::iterator ihp = lHitPat.begin(); 
  while(1){
    ip2 = ip1; ip2++;
    if (ip1==lPlnRef.end() || ip2==lPlnRef.end()) break;
    int ipl1 = *ip1; int ipl2 = *ip2;
    if (ipl1+1<ipl2) {              // Insertion is needed
      lPlnRef.insert(++ip1,ipl1+1);   // Insert missing plane number
      lHitPat.insert(++ihp,    -1);   // Insert missing hit ref.
      ip1--; ihp--; // Now points to inserted elements
    }
    else {                         // Nothing to insert
      ip1++; ihp++;
    }
  }

  insPlaneDone = true;
  return;
}
