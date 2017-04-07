// $Id: TTrackFindKine.cc 14069 2015-09-17 20:44:46Z lsilva $

#include <iostream>
#include <CsSTD.h>
#include "CsErrLog.h"
#include "CsCluster.h"
#include "TTrack.h"
#include "TEv.h"
#include "TKine.h"
#include "THit.h"
#include "TOpt.h"
#include "TSetup.h"

using namespace std;

/*!
  - Find associated MC-track (in "TEv::vKine"), store its ID in "this->IKine".
  - Associated MC-track is that track that supplies a fraction >=
   "TOpt::dCut[10]" of its hits to the reconstructed track ("TOpt::dCut[10]"
   being typically in the range 80/90%, there is no ambiguity. Otherwise the
   candidate w/ the largest fraction is retained.)
  - Are taken into account in the above fraction: original hits, composite
   hits for which the distance to original hit is reasonably small (cf.
   "TEv::ImportCluster"), and hits that are not original but are part of a
   shower, or belong to muWallA.
  - If no MC track can be associated, "IKine=-1"
  - Max. number of hits from the same MC track ("NHsam") and the corresponding
   track ID ("IKsame") are also stored.
  - If associated MC track is found, set also the back reference of MC-track to
   "this" track.
  - In the case of real data events, this function does nothing.
*/

/*
  Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/TTrackFinKine.cc":
  i) Comment out code for cross-checking the number of hits: where a track
  is processed by the "QNewtonFit" its "NHits" is truncated (to its part
  upstream of the RICH).
  ii) Allow close mirrors to coalesce.
  iii) Handle THit's w/ 2 TTrack references (IKine+IKin2)
  iv) Upon option ("dCut[81]"), special processing of DR hits (note that what we have in view here are DRs w/ lead, as opposed to w/o lead). Firstly, non-original hits are accepted. Secondly, if a shower is detected ("TTrack::HasSHower"), independent evaluation of the DR hits, allowing for non original hits and also mirrors, and based on "dCut[81] rather than dCut[10]".
  v) Non-original hits from MuWallA are accepted, on the ground that I) they will not impact significantly the (backward) track parameters and II) only grant to the track a muID valid to the extent that they are born after in the few last X0s of the absorber, condition which is most probably fulfilled.  
  vi) In-flight pi->mu decay, which failed to pass "dCut[10]" are tentatively rescued.
*/


void  TTrack::FindKine()
{

  TEv& ev = TEv::Ref();
  if(! ev.IsMC() ) return;

  const TSetup& setup = TSetup::Ref();
  float coalesce = TOpt::dCut[27];  // Mirrors whose distance < coalesce*d.Resol

  const vector<THit>&  vHit  = ev.vHit ();
  // cast away constness, as we going to add TKine -> TTrack reference (change the object state)
  vector<TKine>& vKine = const_cast<vector<TKine>&>(ev.vKine());
  
  // Remove old TKine -> this TTrack association 
  // if one had been established before
  if (IKine>=0) vKine[IKine].eraseTrackID(Id);

  int evalDRs = 0, nDRHitsMax, nDRSameMax; if (TOpt::dCut[81]) {
    // Special processing of DR hits:
    // - Note that what we have in view here are DRs w/ lead, as opposed to w/o
    //  lead: only in that case (e.g. 2006/2007 data) should "dCut[81]" be set.
    // - In any case, accept non original hits. This in order to cope w/
    //  hadronic showers, which the PR automatically associates to the hadron
    //  track, and which have a small transverse size and hence entail no
    //  big error on the track parameters. 
    // - If indication an e.m. shower is taking place and yet not firmly
    //  established so that hits form the putative shower are associated to
    //  the track ("TTrack::HasShower==1", cf. "TEv::F2reTrack2RICHWall"), make
    //  an independent evaluation of the DR hits, accepting non original hits
    //  and also mirrors, based on the numerical value of "dCut[81]".
    evalDRs = 1; if (HasShower) evalDRs = 2;
  }
  list<int>::iterator i, j; int nSameMax, ikinMax, nSame2nd, ikin2nd;
  for (i = lHitPat.begin(), nSameMax=ikinMax=nSame2nd=ikin2nd = -1,
	 nDRHitsMax=nDRSameMax = 0; i!=lHitPat.end(); i++) {
    if (*i<0) continue;

    //       *************** LOOP OVER TRACK'S HITS ***************

    const THit &hi = vHit[*i];
    int ikin = hi.IKine, ikin2 = hi.IKin2, niks, ik, nSame, nDRHits, nDRSame;
    if      (ikin==-1) continue;         // Hit non associated
    else if (ikin<-1) continue;          // Mirror
    if (ikin2>=0) niks = 2;       // A 2-track hit
    else          niks = 1;
    for (ik = 0; ik<niks; ik++) {

      // ********** LOOP ON (possibly 2) MC TRACKS ASSOCIATED TO HIT **********

      for (j = lHitPat.begin(), nSame = 0, nDRHits=nDRSame = 0;
	   j!=lHitPat.end(); j++) {
	if (*j<0) continue;

	// ***** LOOP OVER TRACK HITS SEARCHING for SAME MC TRACK *****

	const THit &hj = vHit[*j]; int jkin = hj.IKine;

	int isDR; if (evalDRs) {
	  const TPlane &pj = setup.vPlane(hj.IPlane);
	  const TDetect &dj = setup.vDetect(pj.IDetRef);
	  if (dj.IType==17) { isDR = evalDRs; nDRHits++; }
	  else                isDR = 0;
	}
	else isDR = 0;

	if (jkin==ikin || hj.IKin2==ikin) {
	  if (isDR) nDRSame++;           // DR hits: accept non-original ones...
	  else if (setup.InMuWallA(hj.IPlane)) // ...same for MuWallA hits
	    nSame++;
	  else if (hj.IOrig==hi.IOrig) nSame++;
	}
	else if (isDR==2 &&     // In case of shower, accept also mirrors
		 (jkin==-2-ikin || hj.IKin2==-2-ikin)) nDRSame++;
      }
      if (nSame>nSameMax) {
	nSame2nd = nSameMax; ikin2nd = ikinMax;
	nSameMax = nSame; ikinMax = ikin;
	nDRHitsMax = nDRHits; nDRSameMax = nDRSame;
      }
      else if (nSame>nSame2nd && ikin!=ikinMax) {
	nSame2nd = nSame; ikin2nd = ikin;
      }

      ikin = ikin2;  // Try 2nd track if exists (i.e. niks==2)
    }
  }

  NHsame = nSameMax; IKsame = ikinMax;
  if (evalDRs==1) NHsame += nDRSameMax; 

  if (IKsame>=0 &&
      (evalDRs==2 && double(NHsame)>=TOpt::dCut[10]*(NHits-nDRHitsMax) &&
       double(nDRSameMax)>=TOpt::dCut[81]*nDRHitsMax ||
       double(NHsame)>=TOpt::dCut[10]*NHits)) {
    IKine = IKsame;
    vKine[this->IKine].addTrackID(this->Id); // TKine -> this TTrack association
  }
  else {
    IKine = -1;
    if (ikin2nd>=0 && double(NHsame+nSame2nd)>=TOpt::dCut[10]*NHits) {
      //                ***** IN-FLIGHT PI DECAY *****
      // - In-flight: reco'd track comprises both mother particle and its decay.
      // - In case the mass diff between mother and daughter is small, this is
      //  likely to happen, w/ the track parameters of the chimera closely
      //  approximating the mother's. Such a reco will yield a perfectly valid
      //  particle, that can, e.g., allow to reco in turn its mother resonance
      //  (e.g., a pi decaying into mu combined w/ a K to reco a D0). This, both
      //  in RD and in MC. But, in MC, this reco will only be considered a
      //  genuine one, and not a fake, if the association of the reco'd track w/
      //  the MC mother is granted.
      //  (Note: One may or may not want to discard such a candidate pi, in the
      //  case it accumulates enough %X0 while being a mu and is ID'd as such.
      //  This is irrelevant to the present problem.)
      // - One has to make sure that the parameters of the mother aren't upset
      //  in the process. Therefore, let's restrict ourself to pi decaying into
      //  mu (hence small mass diff).
      // - Basic check would be: hits from 1st and 2nd track come indeed in
      //  succession and don't mix.
      // - Better, which in addition rejects accidental chimeras: check we are
      //  indeed dealing w/ a mother and its daughter.
      // - Even in the pi->mu case, the mother's parameter can be inaccurate.
      //  => Add the condition that 1/P, if !=0, not be too much wrong.
      int pids[] = { vKine[ikinMax].PDGid(), vKine[ikin2nd].PDGid()}, ipi = -1;
      if (abs(pids[0])==211) ipi = 0; if (abs(pids[1])==211) ipi = 1;
      if (ipi>=0 && abs(pids[1-ipi])==13 && pids[0]*pids[1]>0) {
	int ikmu, ikpi, iv; if (ipi==0) { ikpi = ikinMax; ikmu = ikin2nd; }
	else                            { ikpi = ikin2nd; ikmu = ikinMax; }
	if ((iv = vKine[ikmu].IVtxOrig)>=0)
	  if (ev.vVtxMC()[iv].IOrig==ikpi) {
	    double rRD = fabs(Hfirst(5)), rMC = 0; if (rRD) {
	      TKine &pi = vKine[ikpi];
	      const THitMC &h = ev.vHitMC()[pi.vHitMCRef().front()];
	      const CsMCTrkHit *csh =
		dynamic_cast<CsMCTrkHit*>(const_cast<CsMCHit*>(h.PtrHit()));
	      rMC = 1/csh->getP().mag();
	    }
	    if (!rRD || // No reco'd momentum: do associate yet...
		// ...in order to document event display or debugging tools.
		//         Else 2 conditions...
		// ...1st is a cut on (reco'd-MC) in term of sigmas (of reco'd).
		// Turns out to reject some interesting soft pions (from D*)...
		fabs(rRD-rMC)<5*sqrt(Hfirst(5,5)) ||
		// ...Rescued by cut on (reco'd-MC) in term of P_MC at low P_MC.
		rMC>.2 && fabs(rRD-rMC)<.05*rMC) IKine = ikpi;
	  }
      }
    }
  }
  return;
}
