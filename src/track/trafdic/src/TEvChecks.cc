#include <iostream>
#include "Traffic.h"
#include "TEv.h"
#include "TSetup.h"
#include "TOpt.h"
#ifdef COMPASS_USE_OSPACE_STD
# include <ospace/std/algorithm>
#else
# include <algorithm>  
#endif

using namespace std;

/*!
  Test TEv object information consistency. 
  Print diagnostics and stop the job in the 
  case of found problems.

*/

void TEv::Checks()
{
  
  const TSetup& setup = TSetup::Ref();

  bool err=false;
  { // Check MC tracks <--> Reconstructed tracks association consistency
    for(int ik = 0; ik < int(vecKine.size()); ik++){ // loop over MC tracks
      TKine& t = vecKine[ik];
      if( t.sTrackID().empty() ) continue;
      set<unsigned int, less<unsigned int> >::iterator id;
      for(id = t.sTrackID().begin(); id != t.sTrackID().end(); id++){ //loop over found track IDs
	list<TTrack>::iterator it; int IKine;
	for(it = listTrack.begin(); it != listTrack.end(); it++){ // loop over TTracks
	  TTrack& tr = (*it); IKine = tr.IKine;
	  if(tr.Id == (*id))  goto found;
	} // TTrack loop
	cout<<"TEv::Checks() ==> Invalid reference: Kine track # "
	    <<ik<<" --> rec. track ID = "<< (*id)<<endl;
	err=true; continue;
      found:
	if(IKine != ik){
	  cout<<"TEv::Checks() ==> Invalid reference: Kine track # "
	      <<ik<<" <-- rec. track ID = "<< (*id)<<endl;
	  err=true; continue;
	}
      } //sTrackID loop
    } // Kine track loop
  } // end of check

  { // Reconstructed tracks  <--> MC tracks association consistency 
    list<TTrack>::iterator it; int IKine;
    for(it = listTrack.begin(); it != listTrack.end(); it++){ // loop over TTracks
      TTrack& tr = (*it); IKine = tr.IKine;
      if(IKine < 0) continue;
      if(IKine >= int(vecKine.size())) {
	cout<<"TEv::Checks() ==> for reconstructed track ID "<<tr.Id
	    <<" TKine track index is our of range "<<IKine<<endl;
	err=true; continue;
      }
      TKine& t = vecKine[IKine];
      if (find(t.sTrackID().begin(), t.sTrackID().end(), tr.Id) ==  t.sTrackID().end()){ // Id not found
	cout<<"TEv::Checks() ==> Invalid reference: Reconstructed track ID "
	    <<tr.Id<<" <-- Kine track #  "<< IKine <<endl;
	err=true; continue;
      }
    } // end of loop ovet TTracks
  } // end of check

  if(err){
    this->PrintMC(1);
    this->PrintRecTracks();
    assert(false);
  }

  { // Check Hit <--> Reconstructed tracks association consistency
    for(int ih = 0; ih < int(vecHit.size()); ih++){ // loop over hits
      THit& h = vecHit[ih];
      if( h.sTrackID().empty() ) continue;
      set<unsigned int, less<unsigned int> >::iterator id;
      for(id = h.sTrackID().begin(); id != h.sTrackID().end(); id++){ //loop over found track IDs
	list<TTrack>::iterator it;
	for(it = listTrack.begin(); it != listTrack.end(); it++){ // loop over TTracks
	  if((*it).Id == (*id)) { // track with this ID is found
	    if(find((*it).lHitPat.begin(), (*it).lHitPat.end(), ih) == (*it).lHitPat.end()){ // ih not found
	      cout<<"TEv::Checks() ==> Invalid reference: Hit # "
		  <<ih<<" <-- rec. track ID = "<< (*id)<<endl;
	      assert(false);
	    }
	    goto next_id;
	  }  
	} // end of TTrack loop
	cout<<"TEv::Checks() ==> Invalid reference: Hit # "
	    <<ih<<" --> rec. track ID = "<< (*id)<<endl;
	assert(false);
      next_id:;
      } //sTrackID loop
    } // THits loop
  } // end of check 
  

  return;
}











