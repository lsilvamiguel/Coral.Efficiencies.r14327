#include <iostream>
#include <iomanip>
#include <math.h>
#include "TTrack.h"
#include "THit.h"
#include "TSetup.h"
#include "TOpt.h"

#ifdef COMPASS_USE_OSPACE_STD
# include <ospace/std/algorithm>
#else
# include <algorithm>  
#endif


/*!
  Replace hit h1 of this track by hit h2. Update corresponding crossreferencies
*/

void TTrack::SubstituteHit(const THit& h1, const THit& h2)
{

  const TSetup& setup = TSetup::Ref();
  int ipl1 = h1.IPlane;
  int ih1  = h1.IHit;

  int ipl2 = h2.IPlane;
  int ih2  = h2.IHit;

  if(ipl1 != ipl2){
    std::cout<<"TTrackSubstituteHit() ==> Hit # "<<h2.IHit<<" do not belongs to the same plane as hit # "
	<<h1.IHit<<std::endl;
    return;
  }

  if(find(lHitPat.begin(), lHitPat.end(), ih1) ==  lHitPat.end()){
    std::cout<<"TTrackSubstituteHit() ==> Hit # "<<h1.IHit<<" do not belongs to the track with ID = "
	<<this->Id<<" ( hit's plane # "<<ipl1<<" )"<<std::endl;
    return;
  }

  std::list<int>::iterator ipr = lPlnRef.begin(); 
  std::list<int>::iterator ihp = lHitPat.begin(); 
  while(ihp != lHitPat.end() && (*ihp) != ih1) {ipr++; ihp++;}
  (*ihp) = ih2; // replace hit association on the plane ipl.

  sIHit.erase (ih1); // remove hit index from the set
  sIHit.insert(ih2); // inset new hit index into the set
  const_cast<THit&>(h1).eraseTrackID(Id);      // cast away constness and remove  hit -> track reference
  const_cast<THit&>(h2).addTrackID(Id);        // cast away constness and add new hit -> track reference

}













