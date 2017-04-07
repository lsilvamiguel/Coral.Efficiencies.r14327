// $Id: TKine.cc 13148 2011-12-28 16:55:25Z kbicker $

#include <iostream>
#include "TEv.h"
#include "TKine.h"
#include "TVtxMC.h"

using namespace std;

/*
  (Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/TKine.cc":
  - Secondary track from a halo induced vertex are declared pile-up, up to
  the 4th generation.
 */

bool TKine::isPileup() const
{  

  TEv& ev = TEv::Ref();

  int iv = IVtxOrig;
  if (iv<0) {
    cout<<"TKine::isPileup() ==> Wrong origin vertex index : "<<iv<<endl;
    return false ;
  }
  else if (iv==0)       // Primary track from primary vertex
    return false;
  else {
    const TVtxMC *v; const TKine *t;
    int level = 0; do {
      v = &ev.vVtxMC(iv);// Mother vertex for this track
      if (v->IOrig<0)    // Mother vertex: !primary and w/o mother track => halo
	return true;
      else {
	t = &ev.vKine(v->IOrig); iv = t->IVtxOrig;
	if (iv==0)       // Secondary track from primary vertex
	  return false; 
      }
    } while (++level<=4);
    return false;
  }
}




