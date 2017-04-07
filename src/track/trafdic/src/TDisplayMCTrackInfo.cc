// $Id: TDisplayMCTrackInfo.cc 13148 2011-12-28 16:55:25Z kbicker $

/*!
  Print information about MC track # ik
*/


#include <iostream>
#include <stdio.h>
#include "TDisplay.h"
#include "TDetect.h"
#include "TEv.h"
#include "TSetup.h"
#include "TConstants.h"
#ifdef COMPASS_USE_OSPACE_STD
# include <ospace/std/algorithm>
#else
# include <algorithm>  
#endif

using namespace std;

/*
  Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/TDisplayMCTrackInfo.cc":
  i) Track's time
  ii) Momentum at first measured point
*/

void TDisplay::MCTrackInfo(int ikin)
{  

  TEv& ev = TEv::Ref(); const TSetup& setup = TSetup::Ref();

  // Print argument track info (it usually has been picked from event display).

  string str1, str2; str1 = ev.vKine(ikin).Name();
  
  int iv = ev.vKine(ikin).IVtxOrig;
  const TVtxMC &vtx = ev.vVtxMC(iv);
  if (ev.vKine(ikin).isPileup()) str2 = " pileup "; 
  else {
    int ikin_orig = vtx.IOrig;
    if (ikin_orig<0) {
      if (iv==0)                 str2 = "prim.vtx";
    }
    else                         str2 = ev.vKine(ikin_orig).Name();
  }
  cout<<"MC tr. "<<ikin<<"  "<<str1<<"  from "<<str2<<endl;

  const TKine &k = ev.vKine(ikin);

  //         ********** MOMENTUM at TARGET EXIT/TIME? **********
  // - Momentum at target exit: i.e. after energy is lost in target material
  //   - In fact, one determines the momentum at first measured point. Just to
  //    make things simpler: no need for retrieving target thickness, etc...
  // - Track's time:
  //   - Determine track's time from one of its component hit.
  //   - (Pixel)GEM are excluded, given that their timing can be smeared, when
  //    the "AmpCorrelationMC" feature is activated.
  const vector<THitMC> &vHitMC = ev.vHitMC();
  const vector<int> &vHitMCRef = k.vHitMCRef(); vector<int>::const_iterator i;
  int info; static double time, pFirst; static const TDetect *dFirst;
  for (i = vHitMCRef.begin(), info = 0; i!=vHitMCRef.end(); i++) {
    const THitMC &h = vHitMC[*i]; const TDetect &d = h.DetRef();
    if (!info) {
      const CsMCTrkHit *csh =
	dynamic_cast<CsMCTrkHit*>(const_cast<CsMCHit*>(h.PtrHit()));
      pFirst = csh->getP().mag(); info |= 0x1; dFirst = &d;
    }
    if (d.IType!=26 && d.IType!=28 && d.IType!=29) { // Exclude (pix)GEMs
      info |= 0x2; time = h.DeltaT; break;
    }
  }
  printf("Pmc: %8.3f GeV X0,Y0,Z0: %6.2f,%6.2f,%6.2f cm Yp,Zp: %8.5f,%8.5f",
	 1/k.Pinv(),vtx.V(0),vtx.V(1),vtx.V(2),k.P(1)/k.P(0),k.P(2)/k.P(0));
  if (info) {
    if (info&0x2) printf("  T = %4.1f ns",time);
    else          printf("  T = ?      ");
    printf(" - @ %s(%.1f cm) P = %7.3f GeV\n",
	   dFirst->Name.c_str(),dFirst->X(0),pFirst);
  }
  else printf("\n");
}




