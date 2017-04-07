// $Id: TTrackAppend.cc 13148 2011-12-28 16:55:25Z kbicker $

#include "TTrack.h"
#include "TEv.h"

using namespace std;

/*!

  This function makes "this" track longer by track "t". 
  If 2-d argument is omitted (default), this method also 
  establishes (for "this" extended track) corresponding
  THit -> "this" track references and finds corresponding 
  TKine track (in case of MC event).
  \warning Track "t" IS NOT deleted immediately, but just labelled with 
  t.IMark = -1, which means "to be deleted afterwards".

*/

/*
  Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/TTrackAppend":
  i) Guestimates.
  ii) Update "NDFs".
*/

void TTrack::Append(TTrack& t, string mode )
{
  // Check if track t is downstream of "this" track
  if (Hlast(0)>t.Hfirst(0)) {
    cout<<"TTrackAppend == > Track ID "<<this->Id<<"  Hlast (0) = "<<Hlast(0)   <<"  Type = "<<this->Type<<endl
	<<"        Appending track ID "<<t.Id    <<"  Hfirst(0) = "<<t.Hfirst(0)<<"  Type = "<<t.Type    <<endl;
    exit(1);
  }
  
  Type |= t.Type;   // "OR" track's bit patterns of zones
  t.IMark = -1;   // label to delete afterwards
  NHits += t.NHits; NDFs += t.NDFs;

  Hlast = t.Hlast;

  TEv* pEv = TEv::Ptr(); // Pointer to TEv object
  if (pEv==NULL) {
    cout<<"TTrack::Append ==> TEv object must exists!"<<endl; exit(1);
  }
  
  // Update lPlnRef, lHitPat  and list of used projections
  lPlnRef.insert(lPlnRef.end(), t.lPlnRef.begin(), t.lPlnRef.end());
  lHitPat.insert(lHitPat.end(), t.lHitPat.begin(), t.lHitPat.end());
  copy(t.sProj.begin(), t.sProj.end(), inserter(sProj, sProj.begin()));

  if (mode=="") { // Default mode. Must be switched OFF for faster append to temporary track

    Scifi |= t.Scifi; // "OR" track's bit patterns of scifi-only

    //    ********** SET ASSOCIATION NEW THits -> this TTrack **********
    list<int>::iterator i;
    for (i = t.lHitPat.begin(); i!=t.lHitPat.end(); i++){
      if (*i<0) continue; // Missing hit
      //  ***** LOOP OVER HIT REFERENCES ON "t" TRACK *****
      // Cast away constness as the state of the THit object will be changed
      THit &h = const_cast<THit&> (pEv->vHit(*i)); 
      // Add corresponding back reference THit->TTrack
      h.addTrackID(this->Id);
    } // End of loop over hit references
    
    this->FindKine(); // find again corresponding MC track (vKine)

    // Update vector of guestimates
    vGuests.insert(vGuests.end(),t.vGuests.begin(),t.vGuests.end());
    GuestsGroups |= t.GuestsGroups;

    // Reset timing:
    //  - Since previous assignment (if any) does not account for newly
    //   appended segment.
    //  - But, caveat, SigmaTime" will not be reset by any other
    //   modification ("AddHit", "SubHit", etc...) made to the track.
    SigmaTime = -1;
  }

  return;
}
