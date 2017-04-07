#include <iostream>
#include <iomanip>
#include <math.h>
#include "TEv.h"
#include "TTrack.h"
#include "THit.h"
#include "TSetup.h"
#include "TOpt.h"

#ifdef COMPASS_USE_OSPACE_STD
# include <ospace/std/algorithm>
#else
# include <algorithm>  
#endif

using namespace std;

/*!
  Remove hit h from the track
  Remove this track ID from the set of assosiated track ID of the hit.
  Remove plane references, if last or first hit on the track had been removed.
*/

void TTrack::RemoveHit(const THit& h)
{

  const TSetup& setup = TSetup::Ref();
  const TEv& ev = TEv::Ref();

  int ipl = h.IPlane;
  int ih  = h.IHit;
  int iproj = setup.vPlane(ipl).IProj;

  if(find(lHitPat.begin(), lHitPat.end(), ih) ==  lHitPat.end()){
    cout<<"TTrackRemoveHit() ==> Hit # "<<h.IHit<<" do not belongs to the track with ID = "
	<<this->Id<<" ( hit's plane # "<<ipl<<" )"<<endl;
    return;
  }

  list<int>::iterator ipr;
  list<int>::iterator ihp;
  list<int>::iterator idf;

  // check if the hit to remove ("ih") is the only hit of projection "iproj" 
  int nsameprj=0;
  for(ihp == lHitPat.begin(); ihp != lHitPat.end(); ihp++) {
    if((*ihp) < 0) continue; // "hole"
    if(setup.vPlane(ev.vHit(*ihp).IPlane).IProj == iproj ) nsameprj++; 
  }
  if(nsameprj == 0) {
    cout<<"TTrackRemoveHit() ==> Hit # "<<h.IHit<<" has projection # "
	<<iproj<<" which is not in the set of this track's projections. Track ID = "<<this->Id<<endl;
  }
  if(nsameprj == 1) {  // the only hit of this projection
    sProj.erase(iproj); // remove iproj from the  set of used projections
  }

  ipr = lPlnRef.begin(); 
  ihp = lHitPat.begin(); 
  while((*ihp) != ih) {ipr++; ihp++; } // now "ihp" iterator points to "h" hit index
  (*ihp) = -1; // set no hit association on the plane ipl.

  sIHit.erase(ih); // remove hit index from the set
  
  // cut "empty" front tail
  ipr = lPlnRef.begin(); 
  ihp = lHitPat.begin(); 
  while(1) {
    if((*ihp) < 0){
      lPlnRef.erase(ipr);
      lHitPat.erase(ihp);
    } else {
      break;
    }
    ipr++; ihp++;
  }
  // cut "empty" back tail
  ipr = lPlnRef.end(); 
  ihp = lHitPat.end(); 
  while(1) {
    ipr--; ihp--;
    if((*ihp) < 0){
      lPlnRef.erase(ipr);
      lHitPat.erase(ihp);
    } else {
      break;
    }
  }

  const_cast<THit&>(h).eraseTrackID(Id);      // cast away constness and remove hit -> track reference
  NHits--;
}

