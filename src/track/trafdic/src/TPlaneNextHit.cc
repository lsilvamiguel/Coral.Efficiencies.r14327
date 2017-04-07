#include "TEv.h"
#include "TSetup.h"
#include "TPlane.h"
#include "THit.h"
#include "THlx.h"
#include <float.h>


/*!

  For the helix H, propagated  on "this" plane, this function returns (in every consecutive call) 
  following parameters:
  
  \param iHit is hit's index in TEv::vHit
  \param Chi2 estimation of Chi2 increment ( always <= ChiCut )
  \param Dist distance from track to hit 
  
  starting from closest hit.
  This function have to be called in the loop, until it returns "false".
  If the function returns "false", this means: 
 "next hit has a Chi2 > ChiCut"
 
 
 \param flag If set to 1, function returns information for closest hit only (could be called once)

 "Division by 2" method is used to avoid loop over all hits on the plane.
 
 \warning Here it is assumed  that hits on this plane (represented by references vecHitRef)  
 are ordered with increase of measured coordinate (in WRS) 

*/

bool TPlane::NextHit(THlx& H, double ChiCut, int& iHit, double& Chi2, double& Dist, int ifl) const
{
  const TSetup& setup = TSetup::Ref();

  TEv& ev      = TEv::Ref();
  const std::vector<THit>& vHit = ev.vHit();

  // find hit, closest to Helix H  position on the plane

  static bool first = true; 
  static int k,k0,idir;
  static double utrk;
  static int nhits;

  if(first){
    first = false;
    nhits= vHitRef().size();
    if(nhits == 0) goto stop;
    THlx Hrot;
    const double& cosa=setup.vDetect(IDetRef).Ca;
    const double& sina=setup.vDetect(IDetRef).Sa;

    H.Rotate(cosa, sina, Hrot);
    utrk = Hrot(1);

    int kl=0; 
    int kr=nhits-1;
    while((kr-kl) > 1){
      k0=(kl+kr)>>1;
      if(utrk <= vHit[vecHitRef[k0]].U) kr=k0;
      else                              kl=k0;
    }
    if(fabs(utrk-vHit[vecHitRef[kr]].U) < fabs(utrk-vHit[vecHitRef[kl]].U)) k0=kr;
    else                                                                    k0=kl;
    iHit=vecHitRef[k0];
    if(H(0,0) == 0) Chi2 = DBL_MAX; // if helix without cov. matrix
    else            Chi2 = H.HitChi2(vHit[iHit]);
    Dist= utrk - vHit[iHit].U;
    if(Chi2 > ChiCut) {             // closest hit is too far
      goto stop;
    }
    if(ifl != 0) {
      first = true;
      return(true);
    };
    k=k0; idir = -1;          // prepare for next call
    return(true);

  } // end of "first-call" block

  k+=idir; // next hit

  if(idir == -1 && k < 0){
    idir = 1;
    k=k0+1;
    if(k == nhits) goto stop;
    goto getChi2;
  }
  if(idir == 1 && k == nhits) {
    goto stop;
  }

 getChi2:
  iHit = vecHitRef[k]; // ref of current hit 
  Chi2= H.HitChi2(vHit[iHit]);
  Dist= utrk - vHit[iHit].U;

  if(idir == -1 && Chi2 > ChiCut) {
    idir = 1;
    k=k0+1;
    if(k == nhits) goto stop;
    goto getChi2;
  };
  if(idir == 1 && Chi2 > ChiCut) {
    goto stop;
  }

  return(true);

 stop:
  first = true; 
  return(false);

};





