// $Id: TPlane.cc 13148 2011-12-28 16:55:25Z kbicker $

#include "TEv.h"
#include "TPlane.h"

using namespace std;

/*
  Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/TPlane.cc":
  i) "fillPixRefs" for alternative hit references (special pixelGEMs)
*/

// Store reference "this Plane -> Clusters in vecHit" in U coordinate sorted order
void TPlane::addHitRef(int ih)
{
  TEv& ev = TEv::Ref();
  vector<int>::iterator i = vecHitRef.begin(), iend = vecHitRef.end();
  double u = ev.vHit(ih).U; while (i!=iend && ev.vHit(*i).U<u) i++;
  vecHitRef.insert(i,ih);

}

// Store reference "this Plane -> Clusters in vecHit" in X,Y,U,V sorted order. (X,Y,U,V in MRS. U and V defined by 
void TPlane::fillPixRefs(int coord, double cc, double sc)
{

  TEv& ev = TEv::Ref();
  const TDetect &d = TSetup::Ref().vDetect(iPlane);
  double cd = d.Ca, sd = d.Sa, cdc = cd*cc+sd*sc, sdc = sd*cc-cd*sc;
  for (int ih = 0; ih<(int)vecHitRef.size(); ih++) {
    int jh = vecHitRef[ih]; const THit &h = ev.vHit(jh);
    double a = cdc*h.U-sdc*h.V;
    vector<int>::iterator i = vecPixRefs[coord].begin(), iend = vecPixRefs[coord].end();
    while (i!=iend) {
      const THit &hi = ev.vHit(*i); double ai = cdc*hi.U-sdc*hi.V;
      if (ai<a) i++;
      else break;
    }
    vecPixRefs[coord].insert(i,jh);
  }
}

// Store reference "this Plane ->MC hits in vecHitMC" in Y coordinate sorted order
void TPlane::addHitMCRef(int ih)
{
  TEv& ev = TEv::Ref();

  vector<int>::iterator i    = vecHitMCRef.begin(); 
  vector<int>::iterator iend = vecHitMCRef.end  (); 
  while(i != iend && ev.vHitMC(*i).Xyz(1) < ev.vHitMC(ih).Xyz(1)) i++;
  vecHitMCRef.insert(i,ih); 

}


// Reset "Plane -> Hits" references
void TPlane::Reset()
{
  vecHitRef  .erase(vecHitRef  .begin(), vecHitRef  .end());
  for (int i = 0; i<4; i++) {
    vecPixRefs[i].erase(vecPixRefs[i].begin(),vecPixRefs[i].end());
  }
  vecHitMCRef.erase(vecHitMCRef.begin(), vecHitMCRef.end());
}
