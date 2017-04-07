// $Id: TEvPrePattern2.cc 14094 2015-11-06 15:28:48Z lsilva $

/*!
  \brief Preliminary pattern recognition (when in alternative "lattice")
  Track segments finding in different detector groups
*/

/*
  Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/TEvPrePattern.cc" (in addition to those of
  standard "lattice":
    i) Iterating.
   ii) TWatch.
  iii) Taking into account THit::Status, depending upon iteration and step (
      projection or space).
   iv) (On option) "PR_GEM_QUADRUPLES".
    v) Reevaluation of LR ambiguities after track finding "ReLR"
   vi) Cleaning and pick-up of unused hits either based on "Dico" fitting (case
      of zones 0,1,2) or on straight line fitting (case zones 3,4) track
      segments issued by track finding.
  vii) Refined ideal PR.
 viii) Pattern recognition spanning 2 zones in Z. (options "PR_LOW_P_FORWARD"
      and "PR_2ZONE_ZPROJ").
   ix) PixelGEM/MMs.
*/


#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cassert>
#include "TH2.h"
#include "CsErrLog.h"
#include "CsZone.h"
#include "CsCluster.h"
#include "CsHistograms.h"
#include "CsEvent.h"
#include "CsEventUtils.h"
#include "CsMCUtils.h"
#include "CsDetector.h"
#include "CsGEMDetector.h"
#include "TWatches.h"
#include "Traffic.h"
#include "TOpt.h"
#include "TSetup.h"
#include "TEv.h"
#include "TAlgo2.h"
#include "TConstants.h"
#include "TLattice.h"
// Infra has to be included last: otherwise compiler error in "../../../include/DaqDataDecoding/Chip.h"
#include "TMatrixF.h"
#include "TVectorF.h"

using namespace std;

// Interface
int findMirrorRef(int iref, int iHit, int iMirror, int *href);

void TEv::PrePattern2(CsZone* zone)
{
  const TSetup& setup = TSetup::Ref();

  int igroup;
  for (igroup = 0; igroup<int(setup.vIplFirst().size()); igroup++){
    // loop over all det. groups
    if (setup.Group2Zone(igroup)!=zone) continue; // reconstruction in this group was not requested

    Traffic::Ref().Stopwatch.Start(17+(igroup<3?igroup:3));

    if (TOpt::ReMode[4]!=0) {  // *************** IDEAL PR (MC) ***************

      // Note: The following tries and handles corectly the case of a pi (or K)
      // decaying to a mu, which we would like to reco as a single track (which
      // will be done in non ideal case).

      // iPRpar REQUIREMENT: Select that, expected to be lesser demanding,
      // corresponding to `subsequent' PR iterations
      int nHitsMn = TOpt::iPRpar[igroup*10+2 + 3];

      float incrVSAT = 1.3;// ********** VSAT(=SCIFIS/cGEMS) ENHANCED **********
      float nEffective; // # of Effective (w/ VSAT enhanced) Hits
      // - Enhancement factor used in the real (as opposed to ideal) PR to
      //  compensate for the low redundancy of VSAT. (Note that it is not
      //  certain that the enhancement factor used here correspond exactly to
      //  those used in the real PR. Particularly so for what concerns cGEMs.)
      // - Applied to scifis but also those GEMs (refered to as cGEMs) used in
      //  the VSAT region: pixelGEMs and, if "PR_ENHANCE_CENTRAL_GEMs" defined,
      //  GEMs w/ central zone activated.
      // - PixelMMs are not enhanced, as of 2009/08: waiting for an evaluation
      //  of the impact.
      if      (igroup==0) incrVSAT = nHitsMn/5.;
      else if (igroup==1) incrVSAT = nHitsMn/4.;
      else if (igroup==2) incrVSAT = nHitsMn/4.;
      else                incrVSAT = 0;

      int nMCTs = (int)vecKine.size(), nPats = nMCTs/32+1;
      unsigned int *mctPats; mctPats = new unsigned int[nPats];
      memset((void*)mctPats,0,nPats*sizeof(unsigned int));
      int imct; for (imct = 0; imct<nMCTs; imct++) {
	int ipat = imct/32; unsigned int pat = 1<<imct%32;
	if (mctPats[ipat]&pat) continue; mctPats[ipat] |= pat;

	//       ********** LOOP OVER MC TRACKS **********

	const TKine *ki = &vecKine[imct]; if (!ki->Q()) continue;
	const TKine *kf = 0, *k;  // mu continuation?
	const vector<int> &vtcs = ki->vecVtxMCRef; int pid = ki->PDGid();
	if ((abs(pid)==211/* pi */ || abs(pid)==321/* K */) && vtcs.size()>0) {
	  const vector<int> &vtrks = vecVtxMC[vtcs[0]].vecTrkRef;
	  int nVtrks = (int)vtrks.size();
	  if (nVtrks==1 /* nu's being skipped by COMGeant */) {
	    for (int j = 0; j<nVtrks; j++) {
	      int jmct = vtrks[j];
	      k = &vecKine[jmct]; if (abs(k->PDGid())==13) {
		kf = k; mctPats[jmct/32] |= 1<<jmct%32; break; 
	      }
	    }
	  }
	}
	TTrack t; t.Type = 1<<igroup;
	int ik; for (ik = 0, nEffective = 0, k = ki; ik<2; ik++) {
	  vector<int>::const_iterator ih;
	  for (ih = k->vHitRef().begin(); ih!=k->vHitRef().end(); ih++) {

	    //       ***** LOOP OVER TRACK'S HITS *****

	    if (vecHit[*ih].IKine<0) continue;    // ***** SKIP NON GENUINE HITS
	    if (TOpt::ReMode[4]>1 &&    // ***** SUPERIDEAL PR...
		vecHit[*ih].IKin2>=0) continue;  // ... skip composite hits
	    if (vecHit[*ih].IOrig!=0) continue;  // ***** SKIP NON ORIGINAL HITS
	    int ipl = vecHit[*ih].IPlane;
	    if (ipl<setup.vIplFirst()[igroup] || setup.vIplLast()[igroup]<ipl)
	      continue;                         // ***** HITS from CURRENT GROUP
	    const TPlane  &p = setup.vPlane(ipl);
	    if (p.IFlag==0) continue;                      // ***** TPlane IS ON

	    t.AddHit(vecHit[*ih]);                    // ***** ADD HIT TO TTrack
	    const TDetect &d = setup.vDetect(p.IDetRef);
	    nEffective +=                 // ***** INCREMENT # of EFFECTIVE HITS
		    (d.IType==29 /* PixelGP */) ? 2 : 1;
	    // VSAT(=scifis/cGEMs, but excluding Sis and cMMs) enhanced
	    if (d.IType==22 || d.IType==28) nEffective += incrVSAT;
	    if (d.IType==29)                nEffective += 2*incrVSAT;
	  }
	  
	  if (!kf) break; k = kf;
	}
	if (nEffective<nHitsMn) continue;
	//                                        ***** FIT... *****
	int fitOK = 0;   // ... bound to fail were there any reinteraction
	while (fitOK==0) {
	  if (t.NDFs>5) {                  // #Hits (in fact "TTrack::NDFs">5...
	    if (t.QuickKF(1,0) && t.QuickKF(-1,0)) {
	      t.Hfirst(5,5) = 1.E-20; // ...Fixed momentum fit
	      t.Hfirst(5) = k->Pinv();// ...W/ momentum taken from MC truth bank
	      if (t.FullKF(1) && t.Chi2tot/(t.NDFs-4)<TOpt::dCut[9]) fitOK = 1;
	    }
	  }
	  else if (t.NDFs>4 && t.QuickKF(1,0) &&  // TTrack::NDFs>4: Straigh fit
		   t.Chi2tot/(t.NDFs-4)<TOpt::dCut[9]) fitOK = 1;
	  else break;
	  if (fitOK==0) {                          // Fit failed...
	    list<int>::iterator ipr = t.lPlnRef.end(); ipr--;
	    list<int>::iterator ihp = t.lHitPat.end(); ihp--;
	    const TStation *&s = setup.vPlane(*ipr).Station;
	    while (ihp!=t.lHitPat.begin() && setup.vPlane(*ipr).Station==s) {
	      // ...Erase all ref's in TPlanes of end TStation
	      t.lPlnRef.erase(ipr--);      // erase plane ref.
	      if (*ihp>=0) {
		t.NHits--;
		const TPlane &p = setup.vPlane(*ipr);
		const TDetect &d = setup.vDetect(p.IDetRef);
		t.NDFs     -= (d.IType==29 /* PixelGP */) ? 2
		  : 1;
		nEffective -= (d.IType==29 /* PixelGP */) ? 2
		  : 1;
		// VSAT(=scifis/cGEMs, but excluding Sis) enhanced
		if (d.IType==22 || d.IType==28) nEffective -= incrVSAT;
		if (d.IType==29)                nEffective -= 2*incrVSAT;
		const_cast<THit&>(vHit(*ihp)).eraseTrackID(t.Id);   // Erase hit -> track reference
	      }
	      t.lHitPat.erase(ihp--);      // erase hit ref.
	    }
	    if (nEffective<nHitsMn) break;
	  }
	}
	
	if (fitOK==1 &&           // ***** TTrack Successfully fit...
	    nEffective>=nHitsMn)  // ***** ...AND FULFILLS iPRpar REQUIREMENT...
	  listTrack.push_back(t);             // ...STORE THE Ttrack *****
     
      } // END OF LOOP OVER Kine tracks
      delete mctPats;

      Traffic::Ref().Stopwatch.Stop(17+(igroup<3?igroup:3));
      return;
    }


    //--------------------------------------------------------
  
    TAlgo2 algo;
    // ***** WORKING ARRAYS of FIXED DIMENSION *****
    const int maxpl   = 300, maxhit  = 3000;
    // Max. # of proj. to be used in PR. It is given a fix value, while the
    // total # of projections could be retrieved from "setup.vProj().size()".
    // The rationale being that this is a convenient way to keep this # under
    // control and prevent it from diverging over time. 13 is enough to cope w/
    // all of COMPASS detectors orientations, be they all gathered in a single
    // zone as can happen in case of field off data processing. 13 = XY +/-45
    // (GM,MM,FI) +5/-85 (SI) +/-10 (ST,P*) +15 (MB) +/-20 (DC) +/-30 (DW)
    const int maxproj = 13;
    //                         Index of generic plane in list of active planes
    // (Don't remember why it has to have dimension > maxpl!?)
    int   jpls[1000];
    int   ipls[maxpl];      // Index of active plane in list of all planes
    float xpl [maxpl];               // X of plane
    int   idpl[maxpl];               // Type or index of plane in vDetect
    float res [maxpl];               // resolution of plane
    float tol [maxpl];               // tolerances for track finding of plane
    int    fh [maxpl], lh[maxpl];    // first and after-the-last hit number on plane
    float cosa[maxpl];               // cos(a)
    float sina[maxpl];               // sin(a) of planes
    float Uhit[maxhit];              // Coord of hit
    int   href[maxhit];              // ref to vHit
    int   ifl [2*maxhit];            // Flag used hit
#define PR_GEM_QUADRUPLES
    // GEM_QUADRUPLES (GEM quads) is an attempt at taking advantage of amplitude
    // (amp) correlation in GEMs.
    // - It considers both standard and stripped pieces of pixelGEMs, i.e. "GM"
    //  and "GP..X|Y|U|V", i.e. TDetect::IType==26 and 28. But was designed for
    //  the former and is not well suited for the latter...
    // - It's built around quads of hits mutually consistent, and doubly so: in
    //  amps, X/Y and U/V, and spatial coords X/Y/U/V (cf. TEv::Quadruples).
    //  This yields very reliable combinations, which are then given a bonus. 
    //  The rationale behind this choice is that since we have a high redundancy
    //  of GEMs, the difficult cases are crowded events w/ ambiguities. It
    //  turns out not to fit the VSAT region, where pixelGEMs are sitting and
    //  where redundancy is low.
    // - The bonus (mentioned supra) is twofold:
    //   I) In proj. search: dummy planes are created and filled w/ the relevant
    //    projections of the quads. These dummy planes are two in any case. W/ a
    //    U<V<Y<X ordering of the GEM planes, they are U'(Y,Z) and Y'(U,V). The
    //    bonus comes then from the larger #hits that, the U and X projected,
    //    corresponding tracks can then collect (the discrimination among
    //    candidate tracks being primarily based on #hits). 
    //  II) In space search, the candidate tracks including full quads (4 out of
    //    4 hits) are given a bonus.
    // - The implementation of this algorithm is so that it imposes strict
    //  constraints on the setup of the GEM stations (including the ordering of
    //  the planes w/in the stations): it must be exactly the same for each
    //  station w/in any given zone. If in particular, any GEM plane is turned
    //  off, the whole thing has (cf. TOpt::ReMode[17]) to be disabled.
    // - Other limitations of this present implementation:
    //   - It's in fact enabled only in earlier iterations: the more crowded
    //    ones. (Although, it could in principle be extended).
    //   - 2-zone proj. do not benefit from bonus (I), cf "search_quads" flag.
    //   - Main limitation: the dummy planes would be more useful if they were
    //    Y' and Z' (and if there were 4 of them, btw): would help to make for
    //    the lack of redundancy in VSAT, where Y and Z proj. are the more
    //    likely ones: w/ connection to scifis FI05, in particular.
    // - The impact of the procedure has not been evaluated for a long time. It
    //  is in any case not that big, when measured on MC.
    static float hinfos[maxhit];
    int *ihinfos; ihinfos = (int*)hinfos; if (sizeof(int)!=sizeof(float))
      CsErrLog::mes(elFatal,"Compiler error: sizeof(int) != sizeof(float)");
    static int   hits[(TConstants_NTtrack_max+1)*maxpl*2];  // Found hits

    float X0 = float(setup.iPlane2Detect(setup.vIplFirst()[igroup]).X(0)); // reference plane


    int nPrjs = 0, ndisc = 0, ipr;
    int projind[maxproj]; // Array of indices of projections in this detector group
    int discind[maxproj]; // Array of indices of discarded projections

    int Ntk, nHitsTot, npl, ier;
    int yProj = 0, zProj = 1;  // Indices of Y/Z proj. in "TSetup::vecProj" (=0/1, cf. TSetup::Init)
    int yPrj = -1, zPrj = -1, uPrj = -1; // Indices of Y/Z proj. in "projind"
    int pixProjs[] = {yProj,zProj,-1,-1}; // Projections covered by pixelGEMs
    int pixInds[] = {0,0,0,0}; // Satus of "pixProjs" in "proj(disc)ind"
    double aPixProjs[] = {0,90,45,-45}; const double sqr2 = sqrt(2.);
    double cPixProjs[] = {1,0,sqr2,sqr2}, sPixProjs[] = {0,1,sqr2,-sqr2};
    if (TOpt::iCut[29]&1<<igroup) { // Zones where special pixelGEM projection
      for (int i23 = 2; i23<=3; i23++) {
	double a = TOpt::dCut[80+i23]; aPixProjs[i23] = a; a *= M_PI/180;
	cPixProjs[i23] = cos(a); sPixProjs[i23] = sin(a);
      }
    }
    for (ipr = 0; ipr<(int)setup.vProj().size(); ipr++) {
      if      (abs(floor(setup.vProj()[ipr]/10.+.5)-floor(aPixProjs[2]+.5))<=1)
	pixProjs[2] = ipr;
      else if (abs(floor(setup.vProj()[ipr]/10.+.5)-floor(aPixProjs[3]+.5))<=1)
	pixProjs[3] = ipr;
    }

#ifdef PR_GEM_QUADRUPLES
    // GEM_QUADRUPLES requires XYUV proj's: this must be checked in the code.
    // For the time being, let's simply exclude zones 0x5, which lack at least
    // UV. (Note: revisiting the code, I now find this excluding 0x5 at best
    // useless...)
    int gem = (TOpt::ReMode[17] && igroup && igroup<4) ? 26 : -3;
    int first_gem_prj = -1, first_gem_plane = 10000;
#endif

    //        *************** COUNTING PLANES IN proj... ***************
    // *************** ...AND TOTAL #HITS, SINES AND COSINES ... ***************
#ifdef PR_GEM_QUADRUPLES
    // *************** ... AND WHICH PROJ 1ST gem PLANE IS IN ***************
    // (in view of  reordering array of proj's with 1st gem plane (of a station
    // of 4) belonging to 1st proj in array)
#endif
    float cos_prj[maxproj], sin_prj[maxproj]; // Cos(a), Sin(a) of proj
    for (ipr=nHitsTot = 0; ipr<(int)setup.vProj().size(); ipr++) {

      //         ********** LOOP OVER ALL PROJ'S **********

      int npl = 0, jpr = -1, mhits = 0, nVSATs = 0, nGEMs = 0;
      float cos_mean=0, sin_mean=0, w_mean=0;  // Mean cos/sin in proj
      for (int ipl = setup.vIplFirst()[igroup]; ipl<=setup.vIplLast()[igroup];
	   ipl++){

	//         ***** LOOP OVER PLANES IN GROUP *****

	const TPlane &p = setup.vPlane(ipl);
	if (p.IFlag==0) {                         // PLANE IS OFF
	  if (setup.vDetect(p.IDetRef).IType==gem || // IF PLANE IS CsGEM...
	      setup.vDetect(p.IDetRef).IType==gem+2) // also when it's IType==28
	    CsErrLog::msg(elFatal,__FILE__,__LINE__, // ...EXIT
     "CsGEM \"%s\" turned off while GEM amplitude correlation requested",
			  setup.vDetect(p.IDetRef).Name.c_str());
	  continue;                                  // ELSE: SKIP
	}
	const TDetect &d = setup.vDetect(p.IDetRef);
	if (npl==0) {
	  int yet, i; for (i=yet = 0; i<nPrjs; i++) if (p.IProj==projind[i]) {
	    yet = 1; break;
	  }
	  for (i     = 0; i<ndisc; i++) if (p.IProj==discind[i]) {
	    yet = 1; break;
	  }
	  if (!yet) {                      // Set current projection
	    if (d.IType==29) {
	      // PixelGP is 1st encountered detector in 2nd loop over "ipr"
	      // ...allow it to start the processing a new proj.
	      for (i = 0, jpr = -1; i<4; i++) if (!pixInds[i]) {
		jpr = pixProjs[i]; break;
	      }
	      if (jpr==-1) break;
	    }
	    else jpr = p.IProj;
	  }
	}
	bool currentProj = false; int iPixProj = -1;
	if (d.IType!=29) // If  not a pixelGP,..
	  // ..current detector contributes w/ its own "p.IProj"...
	  currentProj = p.IProj==jpr;
	else { // ...else, instead, by any of 4 proj's specified in "pixProjs"
	  for (int i = 0; i<4; i++) if (pixProjs[i]==jpr) {
	    currentProj = true; iPixProj = i; pixInds[i] = 1; break;
	  }
	}
	if (!currentProj) continue;
	mhits += p.vHitRef().size(); 
	float w = d.Range/d.Resol; // Not strictly correct for pixelGPs...
	if (11<=d.IType && d.IType<=18) w /=2; // Cases of less reliable drifts
	if (iPixProj>=0) {
	  cos_mean += cos(aPixProjs[iPixProj]/180*M_PI)*w;
	  sin_mean += sin(aPixProjs[iPixProj]/180*M_PI)*w; w_mean += w;
	}
	else {
	  cos_mean += d.Ca*w; sin_mean += d.Sa*w; w_mean += w;
	  if      (jpr==zProj)          zPrj = nPrjs;
	  else if (jpr==yProj)          yPrj = nPrjs;
	  else if (fabs(d.Ca-d.Sa)<.02) uPrj = nPrjs;
	}
#ifdef PR_GEM_QUADRUPLES
	if (d.IType==gem || d.IType==gem+2) {
	  nGEMs++;
	  if (ipl<first_gem_plane) {
	    first_gem_plane = ipl; first_gem_prj = nPrjs;
	  }
	}
#endif
	npl++;
	// Allow for VSAT enhanced. (Here only scifis, pixelGEMs and pixelMPMs.
	// One should have also included plain GEMs w/ activated dead zone, but
	// never mind, this only matters in a debugging context, when all planes
	// other than VSATs have been turned off.)
	if (d.IType==22 || d.IType==28 || d.IType==29 || d.IType==32) nVSATs++;
      }
      if (npl<TOpt::iPRpar[igroup*10+3] &&  // NOT ENOUGH PLANES IN THIS PROJ...
	  nVSATs<2 &&                       // ...and not enough VSATs
	  (!nGEMs || !TOpt::ReMode[17])) {
	discind[ndisc++] = jpr; continue;
      }
      if (nPrjs>=maxproj)
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
		      "# of projections in zone %d > %d",igroup,nPrjs,maxproj);
      cos_prj[nPrjs]=cos_mean/w_mean; sin_prj[nPrjs]=sin_mean/w_mean;
      projind[nPrjs++] = jpr; nHitsTot += mhits;

    }  // End of loop over projections


    if (nHitsTot>maxhit) {
      // ********** IF TOO MANY HITS (i.e. > "maxhit") => GIVE UP **********
      CsErrLog::msg(elError,__FILE__,__LINE__,
		    "# of Hits in group %d > %d",igroup,maxhit);
      Traffic::Ref().Stopwatch.Stop(17+(igroup<3?igroup:3));
      return;
    }
    if (vecHit.size()>2*maxhit) {
      CsErrLog::msg(elError,__FILE__,__LINE__,
		    "# of ALL Hits in group %d > %d",igroup,maxhit);
      Traffic::Ref().Stopwatch.Stop(17+(igroup<3?igroup:3));
      return;
    }

#ifdef PR_GEM_QUADRUPLES
    // ***** REORDERING ARRAY OF PROJ'S:
    //   1ST PROJ = 1ST GEM PLANE'S
    //   2ND PROJ = Y PROJ
    if (first_gem_prj>=0) { // If indeed gem's but incorrect ordering
      float temp; int itemp;
      temp = cos_prj[0]; cos_prj[0] = cos_prj[first_gem_prj];
      cos_prj[first_gem_prj] = temp;
      temp = sin_prj[0]; sin_prj[0] = sin_prj[first_gem_prj];
      sin_prj[first_gem_prj] = temp;
      itemp = projind[0]; projind[0] = projind[first_gem_prj];
      projind[first_gem_prj] = itemp;
      int *prjs[] = {&zPrj,&uPrj,&yPrj}; for (int iprj = 0; iprj<3; iprj++) {
	int *prj = prjs[iprj];  
	if      (*prj==0) *prj = first_gem_prj;
	else if (*prj==first_gem_prj) *prj = 0;
      }
      temp = cos_prj[1]; cos_prj[1] = cos_prj[yPrj];  cos_prj[yPrj] = temp;
      temp = sin_prj[1]; sin_prj[1] = sin_prj[yPrj];  sin_prj[yPrj] = temp;
      itemp = projind[1]; projind[1] = projind[yPrj]; projind[yPrj] = itemp;
      for (int iprj = 0; iprj<2; iprj++) {
	int *prj = prjs[iprj]; if (*prj==1) *prj = yPrj;
      }
      yPrj = 1;
    }
#endif

    // ********** PR SPANNING 2 ZONES in Z: VARIOUS CHECKS **********
#ifdef PR_LOW_P_FORWARD
    if (igroup==0 && TOpt::iPRpar[90]) {// Special PR for low P forward tracks..
      if (setup.vIplFirst().size()<2)
	CsErrLog::mes(elFatal,
	  "Low P forward PR enabled (TraF iPRpar[90]), while #zones==1");
      int nPrjsTot = nPrjs+3;           // ...using up 3 extra projections
      if (nPrjsTot>maxproj)
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
	  "# of projections in zone #0 = %d: cannot enable low P forward PR (TraF iPRpar[90])",nPrjs);
      int *prjs[] = {&zPrj,&yPrj,&uPrj}; char ZYU[] = "ZYU";
      for (int jpr = 0; jpr<3; jpr++) if (*prjs[jpr]==-1)
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
	  "No %c proj. while enabling low P forward PR (TraF iPRpar[90])",
		      ZYU[jpr]);
    }
#endif
#ifdef PR_2ZONE_ZPROJ
    if (igroup==1 &&                    // Special 2-zone proj. search in Z...
	(TOpt::iPRpar[100] || TOpt::iPRpar[110])) {
      if (setup.vIplFirst().size()<2)
	CsErrLog::mes(elFatal,
	  "2-zone Z proj. enabled (TraF iPRpar[100|110]), while #zones==1");
      int nPrjsTot = nPrjs +            // ...using up 1 or 2 extra projections
	TOpt::iPRpar[100]*TOpt::iPRpar[110] ? 2 : 1;
      // "PR_2ZONE_ZPROJ", restricted to "igroup==1 && iter==0", and
      // "PR_LOW_P_FORWARD", "igroup==0 && iter==1", are mutually exclusive =>
      // so that the condition supra is safe enough precaution.
      if (nPrjsTot>maxproj)
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
	  "# of projections in zone 1 = %d: cannot enable 2-zone Z proj. (TraF iPRpar[100|110])",nPrjs);
      if (zPrj==-1)
	CsErrLog::mes(elFatal,
	  "No Z proj. while enabling 2-zone Z proj. (TraF iPRpar[100|110])");
    }
#endif

    // ********** ARRAYS TO STORE RESULTS OF TRACK FINDING IN PROJ **********
    static int  Ntk_prj[maxproj];
    static float Y0_prj[maxproj][TConstants_NTtrack_max], Yp_prj[maxproj][TConstants_NTtrack_max];
#ifdef QUADRATIC_PROJ
    static float Y2_prj[TConstants_NTtrack_max];
#else
    float *Y2_prj = 0;
#endif

    int iter, ipl;
    // Iterations:
    // - By default (cf. infra for particular cases for specific zones):
    //   - iter ==0: Narrow routes, large #hits.
    //   - iter ==1: Still narrow routes, lesser #hits
    //   - iter ==2: Enlarged routes, lesser #hits
    // - iter_mn: -1 is meant to allow for 1 extra iter, w/ PR requirements more
    //           demanding than those of iter #0 (cf. "TAlgo2::FindSpace"), in
    //           order to reduce the combinatorics in iter #0.
    // - Special case of zone #2:
    //   - iter ==1: Hodo HI04 allowed in projection search (this is meant
    //              to rescue cases when one of FI07:8XY is inefficient
    //              and therefore tracks FI07[+GM10]+FI08 are eliminated).
    //   - iter ==2: Still narrow routes, with, at the difference w/ iter ==1,
    //    DWs enabled in the projection search AND HI04 disabled.
    //   - iter ==3: Enlarged routes, w/ DWs disabled again in proj. search.
    int nHitsCut = igroup==0 ? 1000 : 2500;
    if (igroup<=2 && TOpt::iCut[33+igroup]) nHitsCut = TOpt::iCut[33+igroup];
    int iter_mn = nHitsTot<nHitsCut ? 0 : -1, iter_mx, iter_enlarged;
    switch (igroup) {
    case 0: 
    case 1: iter_mx = 5; iter_enlarged = 2; break;
    case 2: iter_mx = 4; iter_enlarged = 3; break;
    case 3: iter_mx = 3; iter_enlarged = 2; break;
    case 4: iter_mx = 2; iter_enlarged = 2 /* N.B.: not used */; break;
    default: iter_mx = 2;  iter_enlarged = 2; // Appeasing the compiler...
    }
    static int nfree;
    for (iter = iter_mn; iter<iter_mx; iter++) {

      // *****************************************************************
      // ******************** TRACK FINDING IN PROJ'S ********************
      // *****************************************************************

    Traffic::Ref().Stopwatch.Start(5);

#ifdef PR_GEM_QUADRUPLES
    // 0: Search for quads starting w/ proj P (arbitrarily P is denoted Y)
    // 1: Include counterpart of quads in proj P and P' = Z
    int search_quads = 0; 
    float Zhit[maxhit];   // Z(?) coord of (U,V)(?) component of quad: hit ...
    int   Zref[maxhit];   // ...and ref
    int   nquads[maxpl];  // # of quads per station
#endif
    bool vSAT0x3Enhanced = // In zones 0x3, VSAT enhancement can be delayed...
      // ..."iCut[28]" specifies the PR iterations where it does takes place.
      1<<(iter>=0 ? iter : 0)&TOpt::iCut[28];
    bool vSATYetToBeEnhanced = iter<3 && !vSAT0x3Enhanced && igroup<=1 &&
      // "iCut[16]" = Pattern of zones where VSAT enhanced
      (TOpt::iCut[16]&1<<igroup);
    //         ***** EXCLUDE LOW PROBA DRIFT HITS, in first few iterations *****
    // (N.B.: GEMs spacer hits (status==-1) are independently excluded from the
    // proj. search, cf. infra.)
    int statusMn = iter<2 ? -1 : -2;
    double scifiTMx = 1000; if (igroup==1 && iter>0 && TOpt::dCut[90]) {
      //  ***** EXCLUDE FAR OFFSET SCIFI HITS, in last (iter>0) iterations *****
      //  This is expected to be enabled for 2010-like setups, where high
      // occupancy in pixelGEMs of zone 0x2: in order to disentangle it, the
      // time gate of scifis is typically enlarged so as to match pixelGEM's.
      // The present time cut is applied in later iters, where pixelGEM are
      // discarded form the projection search.
      double extra =// Extra time width compensating for the trigger jitter of
	// the current event being larger than reference trigger jitter.
	CsEvent::Instance()->getExtraTimeWidth();
      scifiTMx = TOpt::dCut[90]+extra;
    }

    int nSearchs = nPrjs;


    // ***** PR SPANNING 2 ZONES in Z: DEFINE EXTRA "prj" *****
#ifdef PR_LOW_P_FORWARD
    if (TOpt::iPRpar[90] && igroup==0 &&// PR for low momentum forward tracks...
	iter==1) { // ...and "iter" when VSAT(=scifis+cGEMs) are enhanced
      // => all VSAT hits are going to be used up => penalising tracks w/
      // low momentum (and hence standing a chance to be reco'd only in
      // later iterations (where routes are wide opened)) traveling through FI03
      // => Have special search in Z proj. extending down to Rich, and hence
      // yielding very reliable projection tracks, for which we will permit
      // "FindS" to enlarge the routes.
      nSearchs += 3;
      cos_prj[nPrjs]   = cos_prj[zPrj]; sin_prj[nPrjs]   = sin_prj[zPrj];
      cos_prj[nPrjs+1] = cos_prj[yPrj]; sin_prj[nPrjs+1] = sin_prj[yPrj];
      cos_prj[nPrjs+2] = cos_prj[uPrj]; sin_prj[nPrjs+2] = sin_prj[uPrj];
    }
#endif
#ifdef PR_2ZONE_ZPROJ                     // Special 2-zone proj. search in Z...
    if ((TOpt::iPRpar[100] || TOpt::iPRpar[110]) && 
	igroup==1 &&
	iter==0) {  // ...and earlier "iter" (and hence no confict with supra)
      nSearchs += 1;
      cos_prj[nPrjs]   = cos_prj[zPrj]; sin_prj[nPrjs]   = sin_prj[zPrj];
      if (TOpt::iPRpar[100]*TOpt::iPRpar[110]) {
	nSearchs += 1;
	cos_prj[nPrjs+1] = cos_prj[zPrj]; sin_prj[nPrjs+1] = sin_prj[zPrj];
      }
    }
#endif
    for (int ipr = 0; ipr<nSearchs; ipr++) {

      // ********************************************************************
      // ************************* LOOP OVER PROJ'S *************************
      // ********************************************************************

#ifdef PR_GEM_QUADRUPLES
      int nqtot = 0, nstations = 0;
#endif
      int statusMx = 0; bool exclIDWs = false; static float uIDWs;
      int nHits = 0, npl = 0, iproj, iplMn, iplMx;
      int iprZ2 = 0;      // Flag set when any of the 2-zone PR is active
      if (ipr<nPrjs) {
	iproj = projind[ipr];
	iplMn = setup.vIplFirst()[igroup]; iplMx = setup.vIplLast()[igroup];
      }
      // ***** PR SPANNING 2 ZONES in Z: PREPARE for PROJ SEARCH *****
#ifdef PR_LOW_P_FORWARD
      else if (iter==1) {                 // Special case of low P forward PR
	int jpr; switch (ipr-nPrjs) {
	case 0:  jpr = zPrj; break;
	case 1:  jpr = yPrj; break;
	case 2:
	default: jpr = uPrj;
	}
	iproj = projind[jpr];
	if (ipr==nPrjs)   // Z proj. extends from target <-> RICH
	  iplMx = setup.vIplLast()[1];
	else iplMx = setup.vIplLast()[igroup];
	iplMn = setup.vIplFirst()[igroup];
	iprZ2 = 1;
      }
#endif
#ifdef PR_2ZONE_ZPROJ
      else {                              // Special 2-zone proj. search in Z
	iproj = projind[zPrj];
	iplMn = setup.vIplFirst()[0]; iplMx = setup.vIplLast()[1];
	statusMx = 1;
	if (TOpt::iPRpar[100] && ipr==nPrjs) iprZ2 = 2; // LAT upstream of RICH
	else                                 iprZ2 = 3; // VSAT
      }
#endif
      for (ipl = iplMn; ipl<=iplMx; ipl++){ // Loop over planes
	const TPlane &p = setup.vPlane(ipl);
	if (p.IFlag==0) continue;                // Plane is off
	const TDetect &d = setup.vDetect(p.IDetRef);

	bool currentProj = false;
	int iPixProj = -1; static double cPix, sPix;
	if (d.IType!=29) // If  not a PixelGEM...
	  // ..current detector contributes w/ its own "p.IProj"...
	  currentProj = p.IProj==iproj;
	else { // ...else, instead, by any of 4 proj's specified in "pixProjs"
	  for (int i03 = 0; i03<=3; i03++) if (pixProjs[i03]==iproj) {
	    currentProj = true; iPixProj = i03;
	    cPix = d.Ca*cPixProjs[i03]+d.Sa*sPixProjs[i03];
	    sPix = d.Sa*cPixProjs[i03]-d.Ca*sPixProjs[i03]; break;
	  }
	}
	if (!currentProj) continue;

	// ***** PR SPANNING 2 ZONES in Z: DETECTORS SELECTION *****
#ifdef PR_LOW_P_FORWARD
	if (iprZ2==1) {                      // Special case of low P forward PR
	  float xRich = 750, xMM01 = 170, xSM1 = 350;
	  if (d.X(0)>xRich) continue;          // Stop before Rich
	  if (d.IType==26) // Since low momenta are looked for: exclude GEMs
	    continue; // (has the advantage to avoid conflict w/ GEM quadruples)
	  if (d.IType==14)                     // Outer Ystraws: out of scope...
	    continue;    // ...for we are after forward (through FI03) tracks...
	  if (d.IType==15 && d.X(0)<xSM1)      // DC00/1: out of scope...
	    continue;    // ...and useless for DC00/1 would provide w/ many hits
	  if (d.IType==22 && d.X(0)>xMM01)     // FI05: out of scope (low P)...
	    continue;    // ...and tracks through FI04 are reco'd standard way
	  if (d.IType==27 && d.X(0)<xMM01)     // MM01 out (forward tracks)
	    continue;
	  //printf("%s\n",d.Name.c_str());
	}
#endif
#ifdef PR_2ZONE_ZPROJ                        // Special 2-zone proj. search in Z
	if (iprZ2==2) {                // LAT upstream of RICH
	  float xRich = 750;
	  if (d.X(0)>xRich) continue;    // Stop before Rich
	  if (d.IType==22) continue;     // Exclude scifis
	  if (d.IType==26 || d.IType==28)// Exclude GEMs => no GEM quadruples
	    continue;
	  if (d.IType==29 || d.IType==32)// Exclude pixelGP/MPs (just as scifis)
	    continue;
	}
	else if (iprZ2==3) {           // VSAT + SAT downstream of RICH
	  float xRich = 750;
	  if (d.IType!=21 && d.IType!=22 && d.IType!=28 && d.IType!=29 &&
	      d.IType!=26) continue;     // Retain 
	  if (d.IType==26 && d.X(0)<xRich) continue;
	}
#endif

	//     ********** LOOP OVER PLANES in PROJ. **********

	// ***** EXCLUDE UNRELIABLE DETECTORS/HITS from PROJ SEARCH *****
	if (igroup==2) {
	  if (d.IType==16) {                    // Zone 0x4: exclude inner DWs..
	    exclIDWs = iter!=2 ||                  // ...if iter != 2...
	      p.vHitRef().size()>d.Nwires/2*3;     // ...or if too many hits
 	    uIDWs = p.IProj==zPrj ? 50 : 90;
	  }
	  else exclIDWs = false;
	  if (d.IType==43 && iter!=1)           // Zone 0x4, all but iter==1...
	    continue;                               // ...Exclude HI04
	}
	else if (igroup==1 && iter>0 && TOpt::iCut[31]) {
	  // Special case of 2010-like setup:      Zone 0x2, last iters...
	  // (Scifi-enhancement is expected to be active, pixelGEMs
	  // w/ their poor timing would then explode the combinatorics.
	  // In order for their useful info to be yet taken into
	  // account, some limited enhancement is activated in
	  // FindSpace (and not FindProjQ) upon iCut[31].)
	  if (d.IType==28 || d.IType==29)           // ...Exclude pixelGEMs
	    continue;
	}
	//else if (d.IType==32 && iter<0)        // Zone #0, in earlier iters...
	// (Let's drop this rejection, which is probably only valid for the
	// earlier MP versions...)
	//continue;                                 // ...Eclude MP (because hampering the bridging SIs w/ GP)
	if (39<=d.IType && d.IType<43) continue; // Exclude hodos(!HI04), MAs
	if (d.IType==17) continue;               // Exclude ST04s (or DRs)
#ifdef PR_SKIP_NON_ENHANCED
	if (vSATYetToBeEnhanced &&     // VSAT(0x3)enhanced due in later iter...
	    (d.IType==22 || d.IType==28 || d.IType==29 || d.IType==32))
	  continue;                    // ...exclude VSAT(=scifis/cGEMs)
#endif

	idpl[npl] = d.IType; xpl[npl] =d.X(0);

#ifdef PR_ENHANCE_CENTRAL_GEMs
	// Enhance GEM w/ central zone activated the same way as scifis are.
	// This is conditioned by "PR_ENHANCE_CENTRAL_GEMs", which is expected
	// to be defined by "make ./src/track/lattice", cf. ../../GNUMakefile...
	if (idpl[npl]==26 && igroup==1 &&    // ...and restricted to zone 0x2...
	    !d.DZSize(0)) {           // ...and requiring DZ is indeed activated
	  // => Tag the plane to inform "FindP" (Use "0x200, to distinguish it
	  //   from the 0x100 used to tags drifts (in a completely unrelated
	  //   part of the code though, viz. internal to "FindS").
	  idpl[npl] |= 0x200;
	  // And make the offset of the plane available in "Y0_prj"
	  Y0_prj[ipr][0] = d.X(1)*d.Ca+d.X(2)*d.Sa;
	}
#endif

	//      ********** ENLARGE ROUTES ABOUT DRIFT DETECTORS **********
	// In order to account for LR wrong assignments, signal propagation,...
	if ((1<<igroup&TOpt::iCut[27]) &&
	    11<=d.IType && d.IType<=15)  // Except ST04, DWs and DRs
	  res[npl] = d.Resol*1.25;
	else
	  res[npl] = d.Resol;

	fh[npl] = nHits;
	const vector<int> *hitRef = iPixProj>=0 ? &p.vPixRefs(iPixProj) :
	  &p.vHitRef();
	vector<int>::const_iterator ihit;
	for (ihit = hitRef->begin(); ihit!=hitRef->end(); ihit++) {

	  //     ***** LOOP OVER HITS (for PROJ. SEARCH) *****

	  THit &h = vecHit[*ihit];
	  if (h.Status>statusMx) continue;       // Exclude used hits
#warning TODO: Relax the condition for recycling low LR probas
	  if (h.Status<statusMn) continue;       // Exclude low proba drift hits
	  if (h.Status==-1)                      // Exclude GEMs spacer hits
	    // (N.B.: This is done in the projection search, where these
	    // would bias the definition of routes, which are opened w/ a width
	    // derived from the detector's resolution and not from the THit
	    // uncertainty. They will be reintroduced in the space search, where
	    // their bias is less harmful given that the routes are already (if
	    // one neglects the effect of track fit in "FindS") defined and
	    // eventually weighted by their enlarged uncertainty in the final
	    // fit.)
	    continue;
	  if (d.IType==22 && fabs(h.Time)>scifiTMx) // Exclude off-time scifis
	    // (This is expected to be enabled for 2010-like setups, cf. supra.)
	    continue;
	  if (exclIDWs && fabs(h.U)<uIDWs) continue; 

	  if (iPixProj>=0) Uhit[nHits] = cPix*h.U-sPix*h.V;
	  else             Uhit[nHits] = h.U;
	  href[nHits++] = *ihit;
	}  // end of loop over hits on the plane
	lh[npl] = nHits;

	if (++npl==maxpl)
	  CsErrLog::msg(elFatal,__FILE__,__LINE__,
	    "# of Planes in Proj (%.3f,%.3f) in Group %d (=%d) > %d",
			cos_prj[ipr],sin_prj[ipr],igroup,npl,maxpl);

#ifdef PR_GEM_QUADRUPLES
	if (iter<=0 && (d.IType==gem || d.IType==gem+2)) {
	  int nqs = 0;
	  const TStation *&s = p.Station;
	  if (search_quads==0 && ipl==s->IPlanes[0] ||
	      search_quads==1 && ipl==s->IPlanes[3]) {
	    fh[npl] = nHits;
	    if      (search_quads==0) {
	      // Search for quads in proj P (probably P==Y)
	      nqs = Quadruples(ipl,&Uhit[nHits],&Zhit[nqtot],
			       &href[nHits],&Zref[nqtot],ifl);
	      int i, j;
	      for (i = nHits, j = nqtot; i<nHits+nqs; i++, j++) {
		href[i] += nHits; Zref[j] += nHits;
	      }
	      nHits += nqs;
	      nquads[nstations] = nqs;
	      idpl[npl] = gem;
	      xpl[npl] = setup.vDetect(p.IDetRef+2).X(0);
	    }
	    else if (search_quads==1) {
	      // Include counterpart of quads in proj P (P?=Z)
	      nqs = nquads[nstations];
	      for (int i = nqtot; i<nqtot+nqs; i++) {
		Uhit[nHits] = Zhit[i];
		href[nHits++] = Zref[i];
	      }
	      idpl[npl] = gem;
	      xpl[npl] = setup.vDetect(p.IDetRef-2).X(0);
	    }
	    res[npl] = d.Resol;
	    lh[npl] = nHits;
	    nqtot += nqs;
	    if (++npl==maxpl) CsErrLog::msg(elFatal,__FILE__,__LINE__,
		 "# of Planes in Proj (%.3f,%.3f) in Group %d (=%d) > %d",
		 cos_prj[ipr],sin_prj[ipr],igroup,npl,maxpl);
	    nstations++;   // Increment nstations provided ITYpe = gem
	  }
	  else if (search_quads==0 && ipl!=s->IPlanes[0] ||
		   search_quads==1 && ipl!=s->IPlanes[3]) {
	    CsErrLog::msg(elFatal,__FILE__,__LINE__,
  "Quads: zone %d pass #%d %s: ipl=%d != s->IPlanes[%d]=%d",igroup,search_quads,
  d.Name.c_str(),ipl,search_quads==0? 0 : 3,s->IPlanes[search_quads==0? 0 : 3]);
	  }
	}
#endif
      }   // End of loop over planes

      if (TOpt::Print[5]) {
	printf("Det. group %1u    Proj %2u   Nplanes %3u   NHits %5u  ",
	       igroup,iproj,npl,nHits); 
	cout<<flush;
      }
    
      // ********** RECONSTRUCTION IN ONE PROJ **********
    
      int Ntk=TConstants_NTtrack_max, ier(0);
      for (int jpl = 0; jpl<npl; jpl++) {           //***** SET TOLERANCES...
	if (igroup==4) {  // Specific case of beam zone
	  // Quick hack to make for the fact that FI15 has .07 of resolution
	  double beam_res = res[jpl]<0.01177 ? 0.01177 : res[jpl];
	  tol[jpl] = beam_res*TOpt::dPRpar[igroup*10+0];
#warning TODO: Cancel enlarged routes in zone#4
	  if (iter>0) // tol[jpl] *= 1.5;
	    tol[jpl] += TOpt::dPRpar[igroup*10+2];
	}
	else {
	  if (iter<0)
	    tol[jpl] = res[jpl]*3;
	  else {
	    if (igroup==0) {
	      // Case of SI mixed w/ detectors of much worse resolution, in
	      // particular pixelGP/MPs.
	      double minRes = TOpt::dCut[84] ? // SI uncertainty correction term
		.010 : .005;
	      if (res[jpl]<minRes)
		tol[jpl] = minRes*TOpt::dPRpar[igroup*10+0];
	      else
		tol[jpl] = res[jpl]*TOpt::dPRpar[igroup*10+0];
	    }
	    else
	      tol[jpl] = res[jpl]*TOpt::dPRpar[igroup*10+0];
	    // => Add extra tolerance to accomodate curvature, independent
	    //   upon detector resolution
	    if      (iprZ2==2) {             // Special 2-zone proj. search in Z
	      tol[jpl] += TOpt::dPRpar[103]; continue;
	    }
	    else if (iprZ2==3) {
	      tol[jpl] += TOpt::dPRpar[113]; continue;
	    }
	    tol[jpl] += TOpt::dPRpar[igroup*10+2]*fabs(cos_prj[ipr]) +
	                TOpt::dPRpar[igroup*10+3]*fabs(sin_prj[ipr]);
	    if (iter>=iter_enlarged || iprZ2==1) 
	      tol[jpl] += TOpt::dPRpar[igroup*10+2]*fabs(cos_prj[ipr]) +
		          TOpt::dPRpar[igroup*10+3]*fabs(sin_prj[ipr]);
	    if (iter>iter_enlarged || iprZ2==1)
	      tol[jpl] += TOpt::dPRpar[igroup*10+2]*fabs(cos_prj[ipr]) +
		          TOpt::dPRpar[igroup*10+3]*fabs(sin_prj[ipr]);
	  }
	}
      }

#ifdef QUADRATIC_PROJ
      float *Y2 = fabs(sin_prj[ipr])<.05 && igroup==1 ? Y2_prj : 0;
#else
      float *Y2 = 0;
#endif
#ifdef PR_GEM_QUADRUPLES
      if (nstations) search_quads++; 
#endif

      // ********** ACTUAL PROJECTION TRACK FINDING **********
	
      int igr = igroup;
#ifdef PR_LOW_P_FORWARD
      if (iprZ2==1) {                   // Special case low P forward PR
	if (ipr==nPrjs) {            // Z proj.
	  igr = 9; 
	  // "PRpar[93]" is taken from options file. Most sensible requirement
	  // being 7, i.e. 2 (resp. 3) out of 3 (FI03,MM02:3) in zone 0x1 and
	  // 6 (resp. 5) out of 4*1.75 = 7
	  // "PRpar[91]", which is serving 2 different purposes (something
	  // that has to be changed at some point), viz. "FindP" and "FindS",
	  // has to be given 2 different values in the present special case
	  TOpt::iPRpar[97] = 2; // => set it (will be reset infra).
	}
	else {                       // Y,U proj.
	  igr = 0;
	  // PRpar are the same as those for a standard PR in group #0:
	  // Usually "iPRpar[93]==3" (understood as 2.5 by "FindP", but this
	  // is irrelevant here since DC have been excluded).
	}
	// In both cases "dPRpar[10k+1]" is taken from options file.
      }
#endif
#ifdef PR_2ZONE_ZPROJ
      if (iprZ2>=2) igr = 8+iprZ2;           // Special 2-zone proj. search in Z
#endif

      ier = algo.FindProjQ(igr, X0, npl, xpl, tol, fh, lh, nHits, Uhit,
			   Ntk, hits, Y0_prj[ipr], Yp_prj[ipr],
			   iter,Y2,idpl,href);

#ifdef PR_2ZONE_ZPROJ
      if (iprZ2>=2) {                     // Special 2-zone proj. search in Z
	// Clean away duplicata between 2-zone and 1-zone Z proj. Keep the
	// 2-zone copies: these are given a bonus in the space search.
	// (The status of this cleaning is not yet clear:
	//  - Exact duplicata are not harmful: they simply consume CPU.
	//  - Approximate ones may be: they can prevent a genuine track,
	//   expected to be that given by the 2-zone search, from being found.
	// So far we exclude only exact duplicata.)
	int it, jt, jBest, nZTrks= Ntk_prj[zPrj], *discTrks = new int[nZTrks];
	for (it = 0, memset(discTrks,0,nZTrks*sizeof(int)); it<Ntk; it++) {
	  float z0 = Y0_prj[nPrjs][it], zp = Yp_prj[nPrjs][it], best = 1e6;
	  for (jt = 0, jBest = -1; jt<nZTrks; jt++) {
	    float z0p = Y0_prj[zPrj][jt], zpp = Yp_prj[zPrj][jt];
	    float dist = (z0-z0p)*(z0-z0p)+1e4*(zp-zpp)*(zp-zpp);
	    if (dist<best) {
	      best = dist; jBest = jt;
	    }
	  }
	  if (jBest>=0 && best<.005) discTrks[jBest] = 1;
	  //#  define PR_DEBUG_DUPLICATA
#  ifdef PR_DEBUG_DUPLICATA
	  if (jBest>=0)
	    printf("%2d %6.2f %8.5f %2d %6.2f %8.5f %f\n",it,z0,zp,
		   jBest,Y0_prj[zPrj][jBest],Yp_prj[zPrj][jBest],best);
	  else printf("%2d %6.2f %8.5f\n",it,z0,zp);
#  endif
	}
	for (it=jt = 0; it<nZTrks; it++) {
	  if (jt!=it) {
	    Y0_prj[zPrj][jt] = Y0_prj[zPrj][it];
	    Yp_prj[zPrj][jt] = Yp_prj[zPrj][it];
	  }
	  if (discTrks[it]==0) jt++;
	}
	Ntk_prj[zPrj] = jt;
	delete[] discTrks;
      }
#endif
      if(ier > 1) {
	Traffic::Ref().Stopwatch.Stop(5);
	Traffic::Ref().Stopwatch.Stop(17+(igroup<3?igroup:3));
	return;
      }
      Ntk_prj[ipr]=Ntk;
     

    } // end of loop over pat. rec.  projections in this group

    Traffic::Ref().Stopwatch.Stop(5);

    if (TOpt::ReMode[5]!=0) {
      Traffic::Ref().Stopwatch.Stop(17+(igroup<3?igroup:3));
      return;   // Space track finding is OFF
    }

    // *****************************************************************
    // ******************** TRACK FINDING IN SPACE ********************
    // *****************************************************************

    Traffic::Ref().Stopwatch.Start(6);

    int nHits = 0; npl=0; Ntk = TConstants_NTtrack_max; 
    int ipl0  = setup.vIplFirst()[igroup]; const TDetect *dFirst;
    int nPats = (setup.vIplLast()[igroup]-setup.vIplFirst()[igroup]+1)/32+1;
    unsigned int driftsPats[nPats];
    for (int ipat = 0; ipat<nPats; ipat++) driftsPats[ipat] = 0;
#ifdef PR_GP_WHICHSIDE_INFO
    // Upon option or if a single-zone job, disable the making use of the which-
    // side info in stripped GPs.
    bool useGPWhichSide = setup.vIplFirst().size()!=1 && !TOpt::ReMode[45];
#endif
    vector<int>::const_iterator ihit;
    for (ipl = ipl0, dFirst = NULL; ipl<=setup.vIplLast()[igroup]; ipl++) {
      //  *************** LOOP on PLANES (for SPACE SEARCH) ***************
      const TPlane &p = setup.vPlane(ipl); int kpl = ipl-ipl0;
      if (p.IFlag==0) { jpls[kpl] = -1; continue; }     // Plane is off
      const TDetect &d = setup.vDetect(p.IDetRef);
      if (!dFirst) dFirst = &d;
#ifdef PR_SKIP_NON_ENHANCED
	if (vSATYetToBeEnhanced &&     // VSAT(0x3)enhanced due in later iter...
	    (d.IType==22 || d.IType==28 || d.IType==29 || d.IType==32)) {
	jpls[kpl] = -1; continue;      // ...exclude VSAT(=scifis/cGEMs)
      }
#endif
      //if (iter<1 && d.IType==17) continue;   // Exclude ST04 upon early iter's
      idpl[npl] = p.IDetRef;
      xpl[npl] = d.X(0); cosa[npl] = d.Ca; sina[npl] = d.Sa; res[npl] = d.Resol;
      if (11<=d.IType && d.IType<=18) driftsPats[kpl/32] |= 1<<kpl%32; 
#warning TODO: Enlarge drift resolution in space search
      else if (d.IType==39) res[npl] *= 1.5; // Enlarge MA res. (hence route)
      else if (d.IType==40) res[npl] *= 2;   // MA02 still larger
      fh[npl] = nHits;
      jpls[kpl] = npl; ipls[npl] = ipl;
      for (ihit = p.vHitRef().begin(); ihit!=p.vHitRef().end(); ihit++) {

	//    ********** LOOP OVER HITS (for SPACE SEARCH) **********

	THit &h = vecHit[*ihit];
	if (h.Status==1) continue;          // Exclude used hits
	if (h.Status<statusMn) continue;    // Exclude low proba drift hits
	Uhit[nHits] = h.U;
#ifdef PR_ST_WHICHSIDE_INFO
	if (d.IType==11 || d.IType==12) // For the time being, no ST04: !tested
	  // The info is an integer (-1, 0, or +1) stored in "double" variable.
	  // Yet we pass it here as a "double", since it's so used in "FindS".
	  hinfos[nHits] = h.ptrCl->getDigitsList().front()->getData()[1];
#endif
	if (d.IType==22) {
	  double time;
	  if (!h.ptrCl->getTime(time)) hinfos[nHits] = 1.e6;
	  else                         hinfos[nHits] = time;
	}
	if (d.IType==26 || d.IType==28) { // ********** GEM or STRIPPED PixelGEM
#if defined PR_GEM_QUADRUPLES && defined PR_USE_GEM_CORR
	  if      (d.IType==gem || d.IType==gem+2) {
	    //                       ***** FILL "hinfos" W/ GEM QUAD CORRELATION
	    // The quad#+1 is stored. 0 meaning not member of any quad.
	    // (If the quadruple search is enabled, all hits for all (pixel)GEMs
	    // in current zone have been processed, cf. TEv::Quadruples.)
	    int quad = ifl[*ihit]; ihinfos[nHits] = quad>=0 ? quad+1 : 0;
	  }
	  else// Reaching over here means GEM quadruples disabled in this zone..
	    ihinfos[nHits] = 0;// ...=> default "hinfos", i.e. hit is in no quad
#else
	  ihinfos[nHits] = 0;
#endif
	}
#ifdef PR_GP_WHICHSIDE_INFO
	if (d.IType==28 || //    ***** STRIPPED pixelGEMs HAVE WHICH-SIDE INFO
	    // - The info is non ambiguous only for those strips that cross the
	    //  pixelised central piece: for all others the cluster might sit
	    //  close to the border between the 2 so-called hemispheres. These
	    //  unambiguous strips correspond to channels [87,168] according to
	    //  Sebastian Uhl, cf. 2009/01 coll. meeting (instead of [88,167]
	    //  I would have expected). Anyway, this residual ambiguity is dealt
	    //  with in "FindS". Therefore, we decide to store the which-side
	    //  info in any case. This in the MS byte of "ihinfos".
	    // - "FindS" can only handle few stripped pixelGEMs configurations.
	    //  => Have to disable it in case of single-zone magnets-OFF job
	    d.IType==31) { //    ***** STRIPPED pixelMMs: SAME THING
	  if (d.IType==31) ihinfos[nHits] = 0; // Case pixelMM: init
	  if (useGPWhichSide) {
	    const vector<double> &info = h.ptrCl->getAllAnalogDataErrors();
	    if (info.size()>=4) {
	      char coord = d.Name[4]; int side;
	      if ((int)info[3]) { // In CsPM, which-info may be absent.
		// CsPG and CsPM classes do not abide by COMPASS convention...
		if (coord=='X' || coord=='V') side = info[3]>0 ?  1 : -1;
		else                          side = info[3]>0 ? -1 :  1;
		side <<=16; ihinfos[nHits] |= side;
	      }
	    }
	  }
	}
#endif
	if (d.IType==29 || d.IType==32) // Pixel GP/MP: info is v-coordinate
	  hinfos[nHits] = h.V;
	href[nHits] = *ihit;
	if(++nHits==maxhit){
	  if (TOpt::Print[0]) 
	    cout<<"TEv::PrePattern ==> (Space).  More then "<<maxhit<<" hits in group "<<igroup<<endl;
	  Traffic::Ref().Stopwatch.Stop(6);
	  Traffic::Ref().Stopwatch.Stop(17+(igroup<3?igroup:3));
	  return;
	}
      }  // end of loop over hits on the plane
      lh[npl] = nHits;

      // ***** SET TOLERANCES *****
      if (igroup==4) {  // Special case of beam zone
	// Set Si's resolution = scifi's. The routes are usually defined by
	// scifi hits (for scifis are situated at both ends of the zone) and
	// therefore are very imprecised compare to real Si resolution.
	double beam_res = res[npl]<0.01177 ? 0.01177 : res[npl];
	tol[npl] = beam_res*TOpt::dPRpar[igroup*10+0];
	if (iter>0) // tol[npl] *= 1.5;
	  tol[npl] += TOpt::dPRpar[igroup*10+2];
      }
      else {
	if (TOpt::ReMode[34] && igroup==3 && d.IType==26) // Special case: GM11
	  tol[npl] = 0.0577* // Set GM11's resolution = MWPC's
	    TOpt::dPRpar[igroup*10+0];
	else  if (igroup==0 && res[npl]<.005)
	  // Case of SI mixed w/ detectors of much worse resolution, in
	  // particular pixelGP/MPs.
	  tol[npl] = .005*TOpt::dPRpar[igroup*10+0];
	else
	  tol[npl] = res[npl]*TOpt::dPRpar[igroup*10+0];
	// => Add extra tolerance to accomodate curvature, independent
	//   upon detector resolution
	tol[npl] += TOpt::dPRpar[igroup*10+2]*fabs(d.Ca) +
	            TOpt::dPRpar[igroup*10+3]*fabs(d.Sa);
	if (iter>=iter_enlarged) // tol[npl] *= 1.5;
	  tol[npl] += TOpt::dPRpar[igroup*10+2]*fabs(d.Ca) +
	              TOpt::dPRpar[igroup*10+3]*fabs(d.Sa);
	if (iter>iter_enlarged) // tol[npl] *= 1.5;
	  tol[npl] += TOpt::dPRpar[igroup*10+2]*fabs(d.Ca) +
	              TOpt::dPRpar[igroup*10+3]*fabs(d.Sa);
      }

      if (++npl==maxpl){
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
		      "Too many (>%d) planes in zone %d",maxpl,igroup);
	Traffic::Ref().Stopwatch.Stop(17+(igroup<3?igroup:3));
	Traffic::Ref().Stopwatch.Stop(6);
	return;
      }
    }   // *************** END LOOP on PLANES for SPACE SEARCH ***************

#ifdef QUADRATIC_PROJ
    float *Y2 = igroup==1 ? Y2_prj : 0;
#else
    float *Y2 = 0;
#endif
#ifdef PR_2ZONE_ZPROJ                        // Special 2-zone proj. search in Z
    // (N.B.: In low P forward PR case, excluded here, "nSearchs==nPrjs+3"
    int mPrjs = (nSearchs==nPrjs+1 || nSearchs==nPrjs+2) ? nSearchs : nPrjs;
#else
    int mPrjs = nPrjs;
#endif
    ier = algo.FindSpace
      (igroup, iter,
       X0, npl, idpl, cosa, sina, tol, fh, lh, nHits, Uhit, 
       mPrjs, Ntk_prj, cos_prj, sin_prj, Y0_prj, Yp_prj, Y2,
       hinfos, Ntk, hits, ifl);
#ifdef PR_LOW_P_FORWARD
    if (ier>0 && nSearchs==nPrjs+3) {        // Special low P forward PR
      int maxplp = maxpl-npl, nplp, jpl, nHitsp, Ntkp;
      float xMM01 = 170;
      for (jpl=nplp=nHitsp = 0; jpl<npl; jpl++) {
	int ipl = ipls[jpl];
	const TPlane &p = setup.vPlane(ipl);
	const TDetect &d = setup.vDetect(p.IDetRef);
	if (d.IType==15) continue;             // DC00/1: out of scope (forward)
	if (d.IType==22 && d.X(0)>xMM01) continue; // FI05: out of scope (low P)
	if (d.IType==27 && d.X(0)<xMM01) continue; // MM01 out (forward tracks)
	//printf("%s\n",d.Name.c_str());
	nHitsp += lh[jpl]-fh[jpl]; idpl[npl+nplp] = idpl[jpl];
	fh[npl+nplp] = fh[jpl]; lh[npl+nplp] = lh[jpl]; tol[nplp] = tol[jpl];
	cosa[npl+nplp] = cosa[jpl]; sina[npl+nplp] = sina[jpl];
	tol[nplp] = tol[jpl] + TOpt::dPRpar[2]*fabs(cosa[jpl]) +
	  TOpt::dPRpar[3]*fabs(sina[jpl]);
	if (++nplp>=maxplp)
	  CsErrLog::msg
	    (elFatal,__FILE__,__LINE__,
	     "Too many (>%d) planes in low P forward PR (iPRpar[90])",maxplp);
      }
      // Call "FindS" with special "igroup=9", and "mode=1" (i.e. no quad fit,
      // which is not consistent w/ the fact that we're targeting low P, but
      // was found to be just as good on MC (2003.01 D*'s) and furthermore
      // able to recover some particular instance of K0's from
      // "/hpss/../03/P1E/DST_mergedDump/evtdump1-30194.raw"
      // => Set "mode>0" PRpar's, that cannot be taken from option file
      // "iPRpar[94]" at least 3, for all of FI03 have probably already been
      // used by the standard "FindS" called just supra. This has to be set for
      // "FindP(igroup==9)" @ a different value: cannot take it from options.
      TOpt::iPRpar[94] = 4;
      TOpt::iPRpar[91] = TOpt::iPRpar[94];  // Needed because of bug in "FindS"
      Ntkp = TConstants_NTtrack_max-Ntk; // Since "Ntk" tracks already built
      ier = algo.FindSpace
	(9, 1, X0, nplp, idpl+npl, cosa+npl, sina+npl, tol, fh+npl, lh+npl, nHitsp, Uhit, 
	 3, Ntk_prj+nPrjs, cos_prj+nPrjs, sin_prj+nPrjs,
	 Y0_prj+nPrjs, Yp_prj+nPrjs, 0,
	 hinfos,Ntkp,hits+ier,ifl);
      if (ier>0) Ntk += Ntkp;
    }
#endif
    Traffic::Ref().Stopwatch.Stop(6);
    //#define PR_DEBUG_CPUTIME
#ifdef PR_DEBUG_CPUTIME
    for (int i = 0; i<nSearchs; i++) printf("%3d",Ntk_prj[i]);
    static double cpuFS = 0;
    printf(" tracks(proj,iter=%d) => %.3f s/evt\n",
	   iter,Traffic::Ref().Stopwatch.SumTime(6)-cpuFS);
    cpuFS = Traffic::Ref().Stopwatch.SumTime(6);
#endif
    if (ier<0)  {
      Traffic::Ref().Stopwatch.Stop(17+(igroup<3?igroup:3)); return;  // Error
    }

    // *************** STORE TRACKS ***************

    int nHitsMn = TOpt::iPRpar[igroup*10 + (iter>0 ? 5 : 2)];
    unsigned int hitsPats[nPats]; // Same # of int's as "driftsPats" supra.
    for (int it = 0, nh = 0; it<Ntk; it++) { // ***** LOOP ON FOUND TRACKS *****
      TTrack t;
      t.Type = 1<<igroup;   // TTrack's bit pattern of zones

      //                                        ***** STORE HITS *****
      int ih0 = nh; float xPrv = -1000;
      int iRichW = igroup==1 ? 1 : 0; int iGM10 = igroup==2 ? 1 : 0;
      bool hasAssociate = false;
      // ***** VSAT DETECTORS: SINGLE OUT TRACKS going through them (almost)only
      // (These tracks are very special: have little redundancy (w/ the possible
      // exception of SI+pixGEM tracks in hadron setup's zone 0x1), are liable
      // to undergo an interaction inside the spectrometer. => They will be
      // flagged (w/ a TTrack attribute called "Scifi" for historical reason).
      // Scifis have in addition excellent time resolution, which will be used
      // in the cleaning phase.)
      int nScifis = 0, nUsedGoodHits = 0, nIsolated = 0;
#ifdef PR_ENHANCE_CENTRAL_GEMs
      int nCGEMs = 0, enhanceCGEMs = 0; static int nScifisCGEMs, nCGEMsCGEMs;
      if      (igroup==1) {
	enhanceCGEMs = 1; nScifisCGEMs = 2; nCGEMsCGEMs = 2;
      }
      else if (igroup==2) {
	enhanceCGEMs = 1; nScifisCGEMs = 0; nCGEMsCGEMs = 4;
      }
#endif
      int nPixGPs = 0, nPixMPs = 0;
      int nSIs = 0;
      float xY0 = X0-1, xZ0 = X0-1; static float xY1, xZ1;
      THit *scifiHits[2]; int ihs[2];
      for (int ipat = 0; ipat<nPats; ipat++) hitsPats[ipat] = 0;
      while (1) {
	int iref = hits[nh++], ref;
	if      (iref==-1) break;               // Next track's hits
	else if (iref<-1)  ref = href[-2-iref];      // Alternative placeholder
	else               ref = href[iref];         // Normal hit
	if (ref>=int(vecHit.size()))            // Check consistency
	  CsErrLog::mes
	    (elFatal,"space ==> something's wrong with hit references");
	if (iref<-1) continue;                  // Skip alternative placeholder
	THit &h = vecHit[ref];
	const TPlane &p = setup.vPlane(h.IPlane);
	const TDetect &d = setup.vDetect(p.IDetRef);
	float x = d.X(0);
	if (hits[nh]==-1 && x-xPrv>250 &&       // Isolated last hit...
	    d.IType<41 && d.IType!=29) { // (excluding hodos and pixelGEMs)
	  // (e.g. case of a ST04 hit separated from the rest by Rich thickness)
	  hits[nh-1] = -2-maxhit;  // ...Flag it as neither hit nor alternative
	  nIsolated++;
	  continue;  	  // => skip it.
	}
	else if (iRichW &&                     // Case of RICHwall ST04/DR/PS...
		 // This is meant to eliminate tracks w/ few trailing hits
		 // in RICH wall detectors: coming after the thickness of RICH
		 // these hits, when few, are not reliable. It eliminates also
		 // tracks w/ few leading RICH wall hits, which must not harm.
		 TOpt::dCut[67]<d.X(0) && d.X(0)<TOpt::dCut[68]) {
	  if (iRichW==1) {
	    if ((int)t.NDFs>nHitsMn &&           // ...if hit not indispensable
		x-xPrv>250) {                    // ...and isolated
	      int nRichWs = 0, nPSs = 0, nGM04s = 0, nShower = 0;
	      int mh = nh-1, jref = hits[mh];
	      while (jref!=-1) {
		if (jref>=0) {
		  int refj = href[jref]; if (refj>=int(vecHit.size()))
		    CsErrLog::mes
		      (elFatal,"space ==> something's wrong with hit references");
		  const THit &hj = vecHit[refj]; int iPj = hj.IPlane;
		  const TPlane &pj = setup.vPlane(iPj);
		  const TDetect &dj = setup.vDetect(pj.IDetRef);
		  if (TOpt::dCut[67]<dj.X(0) && dj.X(0)<TOpt::dCut[68]) {
		    if      (dj.IType==17) {
		      // Make so that one single coord cannot suffice
		      if (TOpt::iCut[26]) {
			double uj = hj.U; int mj; if ((mj = hj.Mirror)!=-1) {
			  uj += vecHit[mj].U; uj /= 2;
			}
			for (int kref = jref+1; kref<nHits; kref++) {
			  const THit &hk = vecHit[href[kref]];
			  if (hk.IPlane!=iPj) break;
			  double uk = hk.U; int mk; if ((mk = hk.Mirror)!=-1) {
			    uk += vecHit[mk].U; uk /= 2;
			  }
			  if (uk>uj+1.1 /* DR's pitch is = 1*/) break;
			  nShower++; uj = uk;
			}
			for (int kref = jref-1; kref>=0; kref--) {
			  const THit &hk = vecHit[href[kref]];
			  if (hk.IPlane!=iPj) break;
			  double uk = hk.U; int mk; if ((mk = hk.Mirror)!=-1) {
			    uk += vecHit[mk].U; uk /= 2;
			  }
			  if (uk<uj-1.1) break;
			  nShower++; uj = uk;
			}
		      }
		      nRichWs++;
		    }
		    else if (dj.IType==2) nPSs++;
		    else if (dj.IType==26) nGM04s++;
		  }
		  else if (dj.IType!=39 && dj.IType!=40)
		    // Break upon encountering any hit outside the RICHwall
		    // complex, except if it's one of MAs
		    break;
		}
		jref = hits[++mh];
	      }
	      if (jref==-1) {// ...and if RICHwall hits are the last of the list
		if (!nGM04s) {
		  if (TOpt::iCut[26]) {                        // Case DR
		    if (nRichWs<5 && nPSs<3 || nShower>nRichWs*3)
		      iRichW = 2;
		  }
		  else if (nRichWs<4 && nPSs<3) iRichW = 2;    // Case ST04
		}
		else if   (nPSs+nGM04s<3)       iRichW = 2;
	      }
	      if (iRichW!=2) iRichW = 0;
	    }
	    else iRichW = 0;
	  }
	  if (iRichW==2) {
	    hits[nh-1] = -2-maxhit; continue;
	  }
	}
	else if (iGM10 && d.IType==26 && d.X(0)>3000) {  // Special case GM10...
	  if (iGM10==1) {
	    if ((int)t.NDFs>nHitsMn) {           // ...if GM10 not indispensable
	      int nGM10s = 1; int mh = nh, jref = hits[mh]; while (jref!=-1) {
		if (jref>=0) {
		  int refj = href[jref]; if (refj>=int(vecHit.size()))
		    CsErrLog::mes
		      (elFatal,"space ==> something's wrong with hit references");
		  THit &hj = vecHit[refj];
		  const TPlane &pj = setup.vPlane(hj.IPlane);
		  const TDetect &dj = setup.vDetect(pj.IDetRef);
		  if (dj.IType==26) nGM10s++;
		  else break;
		}
		jref = hits[++mh];
	      }
	      if (nGM10s<2) iGM10 = 2;           // ...and single => Skip it
              else iGM10 = 0;
	    }
	    else iGM10 = 0;
	  }
	  if (iGM10==2) {
	    hits[nh-1] = -2-maxhit; continue;
	  }
	}
	xPrv = x;
	if (h.Status<=0) {
	  hasAssociate |=
	    h.ptrCl->hasAssociateClusters();  // Quick hack to select Drifts
	  // 1st group: # of scifis
	  // else scifis only track? => Reset t.Scifi
#ifdef PR_ENHANCE_CENTRAL_GEMs
	  if (d.IType==26 && enhanceCGEMs &&         // ***** CENTRAL GEMs:
	      !d.DZSize(0)) {                  // Require DZ indeed activated...
	    double xu = d.X(1)*d.Ca+d.X(2)*d.Sa;
	    if (fabs(h.U-xu)<2.5) // ...and central (in the coord available) hit
	      nCGEMs++;
	  }
#endif
	  if      (d.IType==22) {                    // ***** SCIFIS:...
	    if (nScifis<=1) {// ...When few, remember them, for later evaluation
	      scifiHits[nScifis] = &h; ihs[nScifis] = nh-1;
	    }
	    nScifis++;       // ...Count them
	  }
	  else if (d.IType==28 || d.IType==29) { // ***** GPs: pixels and strips
	    nPixGPs++; if (d.IType==29) nPixGPs++;
	  }
	  else if (d.IType==32) {                // ***** Pixel MPs
	    nPixMPs++;
	  }
	  else if (d.IType==21) nSIs++;          // ***** SIs
	  t.AddHit(h);
	  int kpl = h.IPlane-ipl0; hitsPats[kpl/32] |= 1<<kpl%32;
	  if (p.IProj!=zProj) {
	    if (xY0<X0) xY0 = d.X(0); xY1 = d.X(0);
	  }
	  if (p.IProj!=yProj || d.IType==29) {
	    if (xZ0<X0) xZ0 = d.X(0); xZ1 = d.X(0);
	  }
	}
	else {                               // ***** ``VERY GOOD'' HITS...*****
	  nUsedGoodHits++;
	  //if (d.IType!=22 || igroup!=4) {	
	  hits[nh-1] = -2-maxhit;  // Flag them as neither hit nor alternative
	  //}
          // Special case of scifis in scifi/Si zone (igroup==4, so far)
	  // Owing to the fact that scifi hits, w/ their crude spatial
	  // resolution, can be shared by several genuine tracks: allow the
	  // worst ones, in term of the ordering returned by "TAlgo::FindSpace",
	  // to be rescued if the # of used ``good hits'' correspond to most of
	  // the deficit in total #hits. The tracks rescued will still have to
	  // pass TEv::CleanTrackList...
	  // (This has been added to rescue Evt 151 (1, 152) in
	  //  "dstar.160.full_p-up.2002.04.outpipe.fz.1",
	  // w/ Si re-enabled in projection search, 
	  // w/ iPRpar[40-45] changed from "3 2 6  3 2 6" to "4 3 11  3 2 11".)
	}
      }

      Traffic::Ref().Stopwatch.Start(7);

      int ok = 1, QNstatus = 0; bool straightTrack = false;
      if (xY0<X0 || xZ0<X0 || xY1-xY0<15 || xZ1-xZ0<10) {// ***** UNCONSTRAINED?
	// Reject cases where angle in Y or Z, not constrained: e.g. FI05X+
	// FI06XY(U)V (the status of FI06U being uncertain: X abscissa unknown
	// as of 05/12) but accepting ST03Z1+ST03Z2
	ok = 0; goto pickup;
      }
      if (nScifis>=(int)t.NDFs-1)   // ***** SCIFI-(almost)ONLY TRACKS *****
#warning TODO: Switch from scifi- to VSAT-(almost)only
	t.Scifi = 1<<igroup;  // TTrack's bit pattern of VSAT-only
      else if (nScifis+nPixGPs+nPixMPs+nSIs>=(int)t.NDFs-1)
	t.Scifi = 0x100<<igroup;    // ***** VSAT-(almost)ONLY TRACKS *****
#ifdef PR_ENHANCE_CENTRAL_GEMs   // Enhance Central GEMs:
      if (enhanceCGEMs && nScifis==nScifisCGEMs && nCGEMs==nCGEMsCGEMs)
	// Tracks in Central GEMs (CGEMs) have fewer hits than PR requirements
	// dictate. They only passed the latter thanks to scifis and CGEMs
	// enhancement. And hence are fragile against any subsequent check of
	// these requirements, as in "TEv::CleanTrackList". To protect them... 
	t.Scifi = 1<<igroup; // ...single them out, w/ the "TTrack::Scifi" tag.
#endif
      if (nScifis && igroup!=4 &&   // ***** AVOID STEALING of SCIFI HITS *****
	  (int)t.NDFs-nScifis>nHitsMn &&// If scifi hits not indispensable...
	  (nScifis==1 || nScifis==2 &&  // ...and either single hit or 2 from...
	   // ...2 different stations (checking THit::IPlane's differing by more
	   // than the max. # of planes in a scifi station)
	   scifiHits[1]->IPlane>scifiHits[0]->IPlane+3)) {
	int iS; for (iS = 0; iS<nScifis; iS++) t.SubHit(*scifiHits[iS]);
	t.UseHitTime(); float tT = t.MeanTime, dtT = t.SigmaTime; int nRemoved;
	for (iS=nRemoved = 0; iS<nScifis; iS++)
	  if (dtT>0 && // Track's time IS measured
	      fabs(scifiHits[iS]->Time-tT)<2*dtT) {
	    t.AddHit(*scifiHits[iS]);
	  }
	  else {	    // If hit's time differs from track's: subtract it
	    nRemoved++; hits[ihs[iS]] = -2-maxhit;
	  }
	nScifis -= nRemoved;
      }

      static float QNchi2; // "static" to avoid "warning... uninitialized"
      static int QNmode; static unsigned int QNnh0;
      if (nScifis>=4 &&        // Ensures more than one scifi station
	  //t.NHits-nScifis<6) { // Still 4*FI+7*MM2:3, or 8*MMs as...
	  // ...in mc.2003.01 Evt#42. (N.B.: Best would be to cut on #stations)
	  (igroup!=0 ||
	   // Allow 5*FIs+4*MM03 as in  mc.2003.01 Evt#202
	   // (Best would still be to ...)
#warning TODO: Determine number of TStations right away
	   t.NDFs-nScifis<4 &&
	   // Allow 6*FIs+3MM03 as in D*.2006.03 Evt#145
	   (nScifis<6 || t.NDFs-nScifis<3))) {
	QNmode = -1; QNnh0 = 5;
	if (t.NHits>10) straightTrack = true;
      }
      else if (nPixGPs+nScifis>=4 || (int)t.NDFs==nSIs+nPixGPs+nPixMPs) {
#warning TODO: Tune the taking into account of pixelGMPs
	// Recycling here what's been for case of scifis supra.
	QNmode = -1; QNnh0 = 5;
	if (t.NHits>10) straightTrack = true;
      }
      else if (igroup==0) {
	QNmode = 3; QNnh0 = 6;
      }
      else if (igroup<3) {
	QNmode = 5; QNnh0 = 6;
      }
      else QNmode = -1;

      static int idebug = 0; // 0x2: ReLR, 0x4: Cleaning, 0x8: QN, 
      if (idebug>1) t.DumpHits(0);

      // ***** TOO MANY (VERY GOOD) HITS HAVE BEEN CANCELED *****
      // => GIVE UP...but:
      // ...but avoid throwing away scifis yet-to-be-enhanced (For this, we
      // refrain from declaring ``scifi-track'' hits ``very good'', cf. infra.
      // => The only scifi-hits that are declared ``very-good'' are from mixed
      // scifi/SAT tracks.
      if (nUsedGoodHits) {
#ifndef PR_HADRON_2004
	if (igroup!=4 && (int)t.NDFs<nHitsMn ||
	    // ...but special case of scifi/Si tracks. (usefull?)
	    igroup==4 && (int)t.NDFs<
	    nHitsMn-TOpt::iPRpar[igroup*10+(iter>0?4:1)]-1) {
	  ok = 0; goto pickup;
	}
#else
	// Specific case of hadron setup: "nHitsMn" can be so low that the
	// condition supra fails to reject even the most stupid tracks. Anyway
	// beam intensity is low in hadron data, the probability for 2-beam
	// event is low and the probability for the genuine incident to see its
	// hits stolen by the accidental one negligibly small => no need for
	// the rescue procedure supra.
	if ((int)t.NDFs<nHitsMn) { ok = 0; goto pickup; }
#endif
      }
      // else?.. For scifiSi tracks, what happens next, depends on whether the
      // cleaning step is able to reintroduce these shared good hits (that 
      // have not been AddHit'd but that have not been silenced in the "hits"
      // array). Have to clarify this point.

      // ***** scifiSi TELESCOPE: CANCEL WHEN NOT ENOUGH scifis *****
      // - This in order to be able to loosen the requirement on total #hits,
      //  w/o producing fakes, from 11 to 10 in the case of 2002 setup. Cf. loss
      //  beam reco in "GridKa:valexakh/.../dstar.160.full_p-up.2002.04...fz.3",
      //  evt #9. (Not yet tested on later setups.)
      // - The cut is applied in 2 steps. Here's it's loosened. The full cut
      //  will be checked after the pickup step, where one try to complete the
      //  track by pick-up.
      if (igroup==4 && nScifis<TOpt::iCut[20]-1) { ok = 0; goto pickup; }

      if (idebug>1) {
	if (it==0) {
	  printf("\nFindSpace %d iter %d: %d tracks <\n",igroup,iter,Ntk);
	  for (int iproj = 0; iproj<nPrjs; iproj++)
	    printf("Proj %d: %d\n",iproj,Ntk_prj[iproj]);
	}
	printf("\nNew track %d\n",it);
	t.DumpHits(0);
      }

      if (t.QuickKF(-1,0)) {           // ********** STRAIGHT LINE KF **********
	t.QuickKF(1,0);
	t.IFit = 0x1;  // Set fit type bit pattern to Straight
	float tgx = t.Hfirst(3), tgy = t.Hfirst(4);

	if (TOpt::ReMode[14] && igroup<=3) { // ***** QN Fit AND ReLR *****

	  if (igroup==1) {
	    // Guest value from horizontal angle
	    static double fldint= setup.MagFieldInt[setup.NMags-2];
#  define INIT_QN_1
#  if defined INIT_QN_0
	    t.Haux(5)=(0-t.Hfirst(3))/(1.E-3 * 0.3 * fldint);
	    t.Haux(1)=t.Haux(2)=t.Haux(3)=t.Haux(4)=0;
#  elif defined INIT_QN_1
	    t.Haux(5) = (0-t.Hfirst(3))/(1.E-3 * 0.3 * fldint);
	    t.Haux(0)=t.Haux(1)=t.Haux(2) = 0;
	    t.Haux(3) = t.Hfirst(1)/
	      (t.Hfirst(0)-TSetup::Ref().TargetCenter[0]);
	    t.Haux(4) = t.Hfirst(4);
#  elif defined INIT_QN_2
	    if (igroup==0) {
	      t.Haux(0)=t.Haux(1)=t.Haux(2) = 0;
	      t.Haux(3) = t.Hfirst(3); t.Haux(4) = t.Hfirst(4);
	      t.Haux(5) = (0-t.Hfirst(3))/(1.E-3 * 0.3 * fldint);
	    }
	    else {
	      t.Haux(3) = t.Hfirst(1)/
		(t.Hfirst(0)-TSetup::Ref().TargetCenter[0]);
	      t.Haux(5) = (t.Haux(3)-t.Hfirst(3))/(1.E-3 * 0.3 * fldint);
	      t.Haux(0)=t.Haux(1)=t.Haux(2) = 0;
	      t.Haux(4) = t.Hfirst(4);
	    }
#  endif
	  }
	  float KFchi2tot = t.Chi2tot;
	  if (QNmode>=0) {
	    if (t.QNewtonFit(1,QNmode)) {
	      QNchi2 = t.Chi2tot/(t.NDics-QNnh0);
	      if (QNchi2<5 || QNchi2<5*KFchi2tot/(t.NDFs-5)) {
		QNstatus = 2;
		if (idebug&0x10)
		  printf("%.3f (%.2f,%.2f,%.5f,%.5f,%.3f)",
			 t.Chi2tot/(t.NDFs-4),
			 t.Haux(1),t.Haux(2),t.Haux(3),t.Haux(4),t.Haux(5));
		if (QNchi2>KFchi2tot/(t.NDFs-5) || (idebug&0x10)) {
		  if (!t.QNewtonFit(1,4)) {
		    QNstatus = 0; t.Chi2tot = KFchi2tot;
		  }
		  else QNchi2 = t.Chi2tot/(t.NDics-QNnh0);
		}
		if (idebug&0x10) {
		  if (QNstatus) printf("-> %.3f (%.2f,%.2f,%.5f,%.5f,%.3f)\n",
		    t.Chi2tot/(t.NDFs-4),
		    t.Haux(1),t.Haux(2),t.Haux(3),t.Haux(4),t.Haux(5));
		  else printf("-> !?\n");
		}
	      }
	      else t.Chi2tot = KFchi2tot;
	    }
	  }
#define ReLR
#ifdef ReLR
	  if (hasAssociate) {	    // ********** REEVALUATE LR **********
	    //#  define DEBUG_ReLR  // Watch out: it's RD-uncompatible
	    //#  define NTUPLE_ReLR
#  ifdef NTUPLE_ReLR
	    static float DA_ReLR =  1;
	    static float da_ReLR =  1;
	    static float da_ReLR2=  1; // !?
	    static float hm_ReLR =  0;
	    static float Pf_ReLR =  1;
	    static float Pf_ReLR_new = 1;
	    static float Pf_ReLR_ReA = 1;
	    static ofstream ofile("ReLR.out");
#  else
	    static float DA_ReLR = .0; // Diff Incidence-Target pointing, on a TTrack basis
	    static float da_ReLR = .0; // Diff Incidence-Target pointing, on a TDetect basis
	    static float da_ReLR2= .05;// Incidence diff. used in 2ND CASE.
	    static float hm_ReLR =  2; // Diff hit-mirror in units of resolution
	    static float Pf_ReLR     = .55;// Final (post reevaluation) proba...
	    static float Pf_ReLR_new = .8;// ...if 2nd case: new hit added
	    static float Pf_ReLR_ReA = .55;// ...if ReLR implies ReAssociation
#  endif
	    float tgtptx = t.Hfirst(1)/t.Hfirst(0);
	    float tgtpty = t.Hfirst(2)/t.Hfirst(0);
	    int modified = 0, chi2Status = 1;

	    // ***** COMPARE INCIDENT ANGLE TO TARGET POINTING...
	    if (idebug&0x2) printf("%f %f    %f %f\n",tgtptx,tgx,tgtpty,tgy);
	    if (TOpt::iCut[21]&1<<igroup) { // Target pointing...
	      if (fabs(tgtptx-tgx)<DA_ReLR && fabs(tgtpty-tgy)<DA_ReLR)
		goto clean_track; 	    // ...if no difference => give up
	    }
	    else {                          // Infinity pointing...
	      if (fabs(tgx)<DA_ReLR && fabs(tgy)<DA_ReLR)
		goto clean_track; 	    // ...if no difference => give up
	    }

	    for (int ih = ih0; ih<nh-1; ih++) {

	      // ***** LOOP ON HITS...

	      int iref = hits[ih];
	      if (iref<-1) continue;             // Skip alternative placeholder
	      THit &hi = vecHit[href[iref]];
	      const TPlane &pi = setup.vPlane(hi.IPlane);
	      const TPlane *&associatedPl = pi.Associate;
	      if (associatedPl) {
		CsCluster *ci = hi.ptrCl;
		// ***** in Drift-like detectors
		const TDetect &d = setup.vDetect(pi.IDetRef);

		// ***** Incidence/target pointing for CURRENT det orientation
		float incidence = d.Ca*tgx+d.Sa*tgy;
		float incDiff = (TOpt::iCut[21]&1<<igroup) ?
		  fabs(incidence-d.Ca*tgtptx+d.Sa*tgtpty) : // Target pointing
		  fabs(incidence);                          // Infinity pointing
		if (incDiff<da_ReLR) continue;

		// ***** Init: variables common to the 2 cases infra
		int jh = ih+1, jref = hits[jh];
		while (jref<-1 && jh<nh-1) {  // Skip alternative placeholder...
		  jref = hits[++jh];          // ...and ``very good used hits''
		}
		bool associated = false;
		CsCluster *cj, *ck = NULL, *cl = NULL;
		// Associated clusters (from initial LR)
		list<CsCluster*> associates = ci->getAssociateClusters();
		if (!associates.empty()) {
		  ck = associates.front(); cl = associates.back();
		}
		if (jh<nh-1) {// ***** NEXT TPlane ASSOCIATE of CURRENT 1? *****
		  THit &hj = vecHit[href[jref]]; cj = hj.ptrCl;
		  const TPlane &pj = setup.vPlane(hj.IPlane);
		  if (&pj==associatedPl) {

		    // ***** 1ST CASE: NEXT TPlane IS ASSOCIATE
#  ifdef DEBUG_ReLR
		    float pbi0 = ci->getLRProb(), pbj0 = cj->getLRProb();
		    int action = 0, success = 0;
		    if (idebug&0x2)
		      printf("%d %d  %f %f",hi.IHit,hj.IHit,pbi0,pbj0);
#  endif

		    // ***** REEVALUATE LR
		    CsEventUtils::setLRToMatchAngle(ci,cj,incidence);

		    THit *hmi, *hmj; float dhmi, dhmj;
		    if (hi.Mirror==-1)     { hmi = NULL; dhmi = 0; }
		    else {
		      hmi = &vecHit[hi.Mirror]; // N.B.: Skipping coalesced mir
		      if (hmi->Status==-4) {
			// Restore LR probas that have just been upset. (Note
			// that the setting dhmi = 0, will prevent "hi" to be
			// changed for "hmi" in any case.)
			ci->setLRProb(1); hmi->ptrCl->setLRProb(0); 
			hmi = NULL; dhmi = 0;
		      }
		      else dhmi = fabs(hi.U-hmi->U)/d.Resol;
		    }
		    if (hj.Mirror==-1)     { hmj = NULL; dhmj = 0; }
		    else {
		      hmj = &vecHit[hj.Mirror]; // N.B.: Skipping coalesced mir
		      if (hmj->Status==-4) {
			// Restore LR probas, just as supra.
			cj->setLRProb(1); hmj->ptrCl->setLRProb(0); 
			hmj = NULL; dhmj = 0;
		      }
		      else dhmj = fabs(hj.U-hmj->U)/d.Resol;
		    }
		    float pf; if (iter<=0) pf = cj==ck ? Pf_ReLR : Pf_ReLR_ReA;
		    else pf = .5;
		    if (dhmi>hm_ReLR && ci->getLRProb()<1-pf) {
		      // ***** EXCHANGE CURENT HIT <-> MIRROR
		      if (hmi->Status<statusMn) {
			modified = 1; chi2Status = -1;
			t.SubHit(hi); t.AddHit(*hmi);
			// Mirror is not is the aray of references "href".
			//  => Let it take the place of its image.
			//  => Reassign the "status" attributes so that
			//    "status<statusMn" still means out of "href".
			href[iref] = hmi->IHit; int status = hi.Status;
			hi.status = hmi->Status; hmi->status = status;
		      }
		      else {
			// The mirror just retained has a reference in "href"
			int irefm = findMirrorRef(iref,hi.IHit,hmi->IHit,href);
			if (irefm==-1) {
			  CsErrLog::msg(elError,__FILE__,__LINE__,
                    "ReLR: Track ID=%d Hit %d has Status %d>=%d and yet no ref",
			      t.Id,hmi->IHit,hmi->Status,statusMn);
			}
			else {
			  modified = 1; chi2Status = -1;
			  t.SubHit(hi); t.AddHit(*hmi);
			  ifl[irefm]++; hits[ih] = irefm; ifl[iref]--;
			}
		      }
		    }
		    if (dhmj>hm_ReLR && cj->getLRProb()<1-pf) {
		      // ***** EXCHANGE NEXT HIT <-> MIRROR
		      if (hmj->Status<statusMn) {
			modified = 1; chi2Status = -1;
			t.SubHit(hj); t.AddHit(*hmj);
			// Mirror is not is the aray of references "href"
			href[jref] = hmj->IHit; int status = hj.Status;
			hj.status = hmj->Status; hmj->status = status;
		      }
		      else {
			// The mirror just retained has a reference in "href"
			int jrefm = findMirrorRef(jref,hj.IHit,hmj->IHit,href);
			if (jrefm==-1) {
			  CsErrLog::msg(elError,__FILE__,__LINE__,
                    "ReLR: Track ID=%d Hit %d has Status %d>=%d and yet no ref",
			      t.Id,hmj->IHit,hmj->Status,statusMn);
			}
			else {
			  modified = 1; chi2Status = -1;
			  t.SubHit(hj); t.AddHit(*hmj);
			  ifl[jrefm]++; hits[jh] = jrefm; ifl[jref]--;
			}
		      }
		    }
#  ifdef DEBUG_ReLR
		    // ***** DEBUG: RETRIEVE MC INFO and ASSESS REEVALUATION
		    if (idebug&0x2)
		      printf("  %f %f",ci->getLRProb(),cj->getLRProb());
		    if (hmi &&
			hmi->setTrackID.find(t.Id)!=hmi->setTrackID.end()) {
		      int mhit = hmi->IHit, mkine = hmi->IKine;
		      if (idebug&0x2) printf("  LR1 -> RL1 %d %d",mhit,mkine);
		      action |= 1; if (mkine>=0) success |= 1;
		    }
		    if (hmj &&
			hmj->setTrackID.find(t.Id)!=hmj->setTrackID.end()) {
		      int mhit = hmj->IHit, mkine = hmj->IKine;
		      if (idebug&0x2) printf("  LR2 -> RL2 %d %d",mhit,mkine);
		      action |= 2; if (mkine>=0) success |= 2;
		    }
		    if (cj!=ck && cj!=cl) {
		      if    (idebug&0x2) printf("  Reassociation!\n");
		      action |= 4;
		    }
		    else if (idebug&0x2) printf("\n");
#    ifdef NTUPLE_ReLR
		    if (action) {
		      if (!ofile)
			CsErrLog::mes
			  (elFatal,"Error openning output file \"ReLR.out\"");
		      ofile<<igroup<<" "
			   <<tgtptx<<" "<<tgx<<" "<<tgtpty<<" "<<tgy<<" "
			   <<a<<" "<<d.Ca*tgtptx+d.Sa*tgtpty<<" "
			   <<pbi0<<" "<<pbj0<<" "
			   <<ci->getLRProb()<<" "<<cj->getLRProb()<<" "
			   <<dhmi<<" "<<dhmj<<" "
			   <<action<<" "<<success<<endl;
		    }
#    endif
#  endif
		    ih = jh; associated = true; // Skip next hit
		  }
		}  // End next TPlane is associate of current one
		if (!associated && ck) {

		  // ***** 2ND CASE: ASSOCIATED TPlane NOT in TRACK
		  //  This can result from 2 situations: either we've reached
		  // the end of the track ("jh<nh-1", supra) or next TPlane
		  // in track is farther away ("pj==associatedPl")
		  //  We can still do something: if current hit has associates,
		  // try to incorporate one of them to the track.
		  //  None of them has been incorporated in the 1st place 'cause
		  // the genuine hit @ current TPlane and/or the genuine
		  // associate on next TPlane have been discarded (on view
		  // of their poor LR perfs). This can only occur if the
		  // incidence on the double layer differs significantly from
		  // that used in LR (There's a problem here 'cause the latter
		  // is not so easily known to Traffic: it differs in fact from
		  // one detector to the next. For simplicity a TraF option is
		  // introduced and it's the responsability of the user to
		  // have LR incidences be all equal and match the TraF option).
		  // Therefore we place a tighter cut on diff in incidence than
		  // for the 1ST CASE above. The criterion for modifying the
		  // track here is to find  a pair of associates with high LR
		  // perf, i.e. _local_ compatibility with track incidence. We
		  // demand in addition that the chi2 change be for better.
		  //  If such a pair of associates is found and only then can
		  // the current hit be changed.
		  if (incDiff<da_ReLR2) continue;
#  ifdef DEBUG_ReLR
		  float pbi0 = ci->getLRProb(), pbk0 = ck->getLRProb();
		  int action = 0, success = 0;
		  if (idebug&0x2)
		    printf("%d %d  %f %f",hi.IHit,-1,pbi0,pbk0);
#  endif

		  // ***** REEVALUATE LR
		  CsEventUtils::setLRToMatchAngle(ci,ck,incidence);
		  if (hi.Mirror!=-1) {
		    THit *hmi = &vecHit[hi.Mirror];
		    if (hmi->Status==-4) {
		      // Restore LR probas that have just been upset.
		      ci->setLRProb(1); hmi->ptrCl->setLRProb(0); 
		    }
		  }

		  float pf;    // Cut on Re probability
		  if (ck->getLRProb()<.5) {
		    // ***** ReOrder associates
		    // The fact that we have to reorder associates may
		    // explain why "cl" was not incorporated to the track in
		    // the first place => We can use a loose "pf" cut.
		    cj = cl; pf = Pf_ReLR;
		  }
		  else {
		    cj = ck; pf = Pf_ReLR_new;
		  }

		  // ***** NEW HIT
		  //  No requirement on mirrors interdistance, since the
		  // problem may come from what's happening on current TPlane.
		  if (cj->getLRProb()>pf) {
		    // Determine THit corresponding to "cj"
#  warning TODO: Make sure that no active plane has a disabled associate. And jpls
		    // It can be absent from the array of reference => ...
		    CsCluster *cmj = cj==ck ? cl : ck; // ...Search also mirror
		    int ipa = associatedPl->IPlane-ipl0;
		    int ir, kref, jpa = jpls[ipa], khit; THit *hj;
		    for (ir = fh[jpa], hj = NULL, kref=khit = -1; ir<lh[jpa];
			 ir++) {
		      if (ifl[ir]==0) {
			THit *h = &vecHit[href[ir]];
			if (h->ptrCl==cj) {
			  kref = ir; hj = h; khit = hj->IHit;
			  if (h->Mirror!=-1) {
			    THit *hmj = &vecHit[hj->Mirror];
			    if (hmj->Status==-4) {
			      // Restore LR probas that have just been upset.
			      cj->setLRProb(1); cmj->setLRProb(0); 
			    }
			  }
			  break;
			}
			else if (h->Mirror!=-1) {
			  h = &vecHit[h->Mirror];
			  if (h->ptrCl==cj) {
			    if (h->Status!=-4) {
			      hj = h; khit = hj->IHit;
			      if (h->Status>=statusMn) kref = ir;
			    }
			    else {
			      // Restore LR probas that have just been upset.
			      cj->setLRProb(0); cmj->setLRProb(1);
			    }
			  }
			  break;
			}
		      }
		    }
		    if (khit!=-1) {
		      float oldChi2, chi2; THit *hmi = NULL;
		      if (QNstatus) {
			if (chi2Status==-1) t.QNewtonFit(1,4);
			oldChi2 = t.Chi2tot/(t.NDics-QNnh0);
		      }
		      else {
			if (chi2Status==-1) { t.QuickKF(1,0); t.QuickKF(-1,0); }
			oldChi2 = t.Chi2tot/(t.NDFs-4);
		      }
		      t.AddHit(*hj);   		 // ***** ADD HIT...
		      if (ci->getLRProb()<1-pf) {
			// ***** ReLR EXCHANGE: exchange curent hit <-> mirror
			if (hi.Mirror==-1)
			  CsErrLog::msg(elError,__FILE__,__LINE__,
		           "THit %d,%d,%.3f: associated CsClusters & no Mirror",
					hi.IHit,hi.IPlane,hi.U);
			else {
			  hmi = &vecHit[hi.Mirror];
			  t.SubHit(hi); t.AddHit(*hmi);
			}
		      }
		      if (QNstatus) {
			if (t.QNewtonFit(1,4))
			  chi2 = t.Chi2tot/(t.NDics-QNnh0);
			else
			  chi2 = oldChi2+1;// => Won't pass "chi2<oldChi2" infra
			if (idebug&0x2) t.DumpHits(1);
		      }
		      else {
			t.QuickKF(1,0);        chi2 = t.Chi2tot/(t.NDFs-4);
		      }
		      bool doModify = chi2<oldChi2;
		      if (doModify) {
			if (hmi) {
			  if (hmi->Status<statusMn) {
			    // Mirror is not in the aray of references "href"
			    href[iref] = hmi->IHit; int status = hi.Status;
			    hi.status = hmi->Status; hmi->status = status;
			  }
			  else {
			    // The mirror retained has a reference in "href"
			    int irefm =
			      findMirrorRef(iref,hi.IHit,hmi->IHit,href);
			    if (irefm==-1) {
			      CsErrLog::msg(elError,__FILE__,__LINE__,
                    "ReLR: Track ID=%d Hit %d has Status %d>=%d and yet no ref",
			        t.Id,hmi->IHit,hmi->Status,statusMn);
			      doModify = false;
			    }
			    else {
			      ifl[irefm]++; hits[ih] = irefm; ifl[iref]--;
			    }
			  }
			}
			if (doModify) {
			  modified = 1; chi2Status = 1;
			  if (kref!=-1) ifl[kref]++; 
			  // Set "hitsPats" in order to prevent other hits to
			  // be added to the plane of the new THit
			  hitsPats[ipa/32] |= 1<<ipa%32;
			  QNchi2 = chi2;
#  ifdef DEBUG_ReLR
			  if (idebug&0x2)
			    cout<<ci->getLRProb()<<" "<<ck->getLRProb();
			  if (hmi) {
			    int mhit = hmi->IHit, mkine = hmi->IKine;
			    if (idebug&0x2)
			      printf("  LR1 -> RL1 %d %d",mhit,mkine);
			    action |= 8; if (mkine!=-2) success |= 8;
			  }
			  int kkine = hj->IKine;
			  if (idebug&0x2)
			    printf("  New LR -> RL == %d %d\n",khit,kkine);
			  action |= 16; if (kkine>=0) success |= 32;
#  endif
			}
		      }
		      if (!doModify) {
			t.SubHit(*hj); chi2Status = -1;
			if (hmi) { t.SubHit(*hmi); t.AddHit(hi); }
#  ifdef DEBUG_ReLR
			if (idebug&0x2)
			  cout<<ci->getLRProb()<<" "<<ck->getLRProb()<<endl;
#  endif
		      }
		    }
		  } // End add new hit
#  ifdef DEBUG_ReLR
#    ifdef NTUPLE_ReLR
		  if (action) {
		    float dhmi = 0, dhmk = 0;
		    CsCluster *cmi = ci->getMirrorCluster();
		    if (cmi) dhmi = fabs(ci->getU()-cmi->getU())/10/d.Resol;
		    dhmk = fabs(ck->getU()-cl->getU())/10/d.Resol;
		    if (!ofile)
		      CsErrLog::mes
			(elFatal,"Error openning output file \"ReLR.out\"");
		    ofile<<igroup<<" "
			 <<tgtptx<<" "<<tgx<<" "<<tgtpty<<" "<<tgy<<" "
			 <<a<<" "<<d.Ca*tgtptx+d.Sa*tgtpty<<" "
			 <<pbi0<<" "<<pbk0<<" "
			 <<ci->getLRProb()<<" "<<ck->getLRProb()<<" "
			 <<dhmi<<" "<<dhmk<<" "
			 <<action<<" "<<success<<endl;
		  }
#    endif
#  endif
		}  // End 2nd case
	      }  // End hit from an associated TPlane
	    }  // End loop on hits for ReLR
	    if (modified) {    // ********** ReLR HAS MODIFIED TTrack **********
	      if (!QNstatus) {   // No QN fit
		t.QuickKF(-1,0); // Straight line fit right now...
		t.QuickKF( 1,0);
#  define TRY_QN_AGAIN
#  ifdef TRY_QN_AGAIN
		if (igroup==1) { // ...And try QN fit once again
		  // Guest value from horizontal angle
		  static double fldint= setup.MagFieldInt[setup.NMags-2];
#    if defined INIT_QN_0  // Useful? Isn't it already defined?
		  t.Haux(5)=(0-t.Hfirst(3))/(1.E-3 * 0.3 * fldint);
		  t.Haux(1)=t.Haux(2)=t.Haux(3)=t.Haux(4)=0;
#    elif defined INIT_QN_1
		  t.Haux(5) = (0-t.Hfirst(3))/(1.E-3 * 0.3 * fldint);
		  t.Haux(0)=t.Haux(1)=t.Haux(2) = 0;
		  t.Haux(3) = t.Hfirst(1)/t.Hfirst(0);
		  t.Haux(4) = t.Hfirst(4);
#    elif defined INIT_QN_2
		  t.Haux(3) = t.Hfirst(1)/t.Hfirst(0);
		  t.Haux(5) = (t.Haux(3)-t.Hfirst(3))/(1.E-3 * 0.3 * fldint);
		  t.Haux(0)=t.Haux(1)=t.Haux(2) = 0;
		  t.Haux(4) = t.Hfirst(4);
#    endif
		}
		if (QNmode>=0 && t.QNewtonFit(1,QNmode)) {
		  // If failed, QN does not overwrites anything
		  QNchi2 = t.Chi2tot/(t.NDics-QNnh0); QNstatus = 2;
		}
#  endif
		chi2Status = 1;
	      }
	    }  // End "modified"
	    if (chi2Status!=1 && QNstatus) {
	      bool ret; if (!(ret = t.QNewtonFit(1,4)))
			  ret = t.QNewtonFit(1,QNmode);
	      if (ret) {
		QNchi2 = t.Chi2tot/(t.NDics-QNnh0);
		QNstatus = 2; chi2Status = 1;
	      }
	      else QNstatus = 0;
	    }
	    if (chi2Status==-1) {
	      t.QuickKF(1,0); t.QuickKF(-1,0);
	    }
	  } // End ReLR
#endif
	} // End QN fit and
      clean_track:

	// *************** CLEANING ***************

	if ( TOpt::iCut[36] && nSIs>=TOpt::iCut[36]) {
	  //   ***** CLEANING of SIs: CHECK the CHI2 of each SI TSpacePt *****
	  // - This kind of space point checking is done here for the sole SIs
	  //  at the exclusion of any other detector type, because of:
	  //  - The small distance separating the TPlane's w/in a TStation.
	  //  - The low magnetic field w/in which they are typically embedded (
	  //   hence straight extrapolation from track's 1st point to the SIs).
	  //  - The high resolution of the SIs, combined w/ the loose resolution
	  //   of neighbouring detectors, which leads to wrong hits being easily
	  //   associated, which wrong hits dramatically affects total chi2.
	  // - "iCut[36]" both enables the procedure and sets a threshold on
	  //  the #hits per TStation required for the cleaning to proceed.
	  vector<THit*> badHits; double yp = t.Hfirst(3), zp = t.Hfirst(4);
	  const TStation *sPrv; TSpacePt spt; int ih;
	  int ihLast = nh-2; while (hits[ihLast]<0) ihLast--;
	  for (ih = ih0, sPrv = 0; ih<=ihLast; ih++) {
	    //      ***** LOOP on TRACK'S HITS
	    // - Search for TSpacePt's = sets of hits from a same SI TStation.
	    // - Compute their chi2.
	    // - If bad chi2 (>3), save worst THit into vector of "badHits".
	    int iref = hits[ih]; if (iref<0) continue;   // Bypass alternative
	    THit &h = vecHit[href[iref]]; int ipl = h.IPlane;
	    const TPlane &p = setup.vPlane(ipl);
	    const TStation *s = p.Station;
	    if ((s!=sPrv || ih==ihLast) && sPrv && sPrv->Type==21) {
	      if (ih==ihLast && // Upon last hit in track, care to account...
		  s==sPrv) spt.hs[spt.nHits++] = &h; // ...current hit.
	      int i; double scU, ssU, sU2;
	      for (i = 0, scU=ssU=sU2 = 0; i<spt.nHits; i++) {
		const THit *h = spt.hs[i];
		const TDetect &d = setup.vDetect(h->IPlane);
		double w = h->SigU-TOpt::dCut[85]; // Undo SI worsening.
		if (w<.001) w = .001;// Still, keep a lower bound = 10 microns
		w /= .001; // Rescale
		double c = d.Ca/w, s = d.Sa/w, dx = d.X(0)-sPrv->X0;
		double U = h->U/w-(yp*c+zp*s)*dx;
		spt.sc2 += c*c; spt.ssc += s*c; spt.ss2 += s*s;
		scU += c*U; ssU += s*U; sU2 += U*U;
	      }
	      double det = spt.sc2*spt.ss2-spt.ssc*spt.ssc;
	      spt.y = (scU*spt.ss2-ssU*spt.ssc)/det;
	      spt.z = (ssU*spt.sc2-scU*spt.ssc)/det;
	      spt.chi2 = sqrt((spt.y*spt.y*spt.sc2+spt.z*spt.z*spt.ss2+sU2+
			       2*(spt.y*spt.z*spt.ssc-spt.y*scU-spt.z*ssU))/
			      spt.nHits);
	      spt.chi2 /= .001;
	      if (spt.nHits>=TOpt::iCut[36] && spt.chi2>3) {
		THit *hWorst; double worse;
		for (i = 0, hWorst = 0, worse = 0; i<spt.nHits; i++) {
		  THit *h = spt.hs[i];
		  const TDetect &d = setup.vDetect(h->IPlane);
		  double c = d.Ca, s = d.Sa, dx = d.X(0)-sPrv->X0;
		  double dU = fabs(h->U-spt.y*c-spt.z*s-(yp*c+zp*s)*dx);
		  if (!hWorst || dU>worse) { hWorst = h; worse = dU; }
		}
		badHits.push_back(hWorst);
	      }
	      if (s!=sPrv) spt.reset();
	    }
	    if (s->Type==21) spt.hs[spt.nHits++] = &h;
	    sPrv = s;
	  }
	  int nBads = badHits.size();
	  if ((int)t.NDFs-nBads<TOpt::iPRpar[10*igroup+5]) {
	    ok = 0;                // ***** IF TOO MANY BAD HITS: ERASE TRACK...
	    printf("PR Clean SIs: trk #%d(0x%x,%.2fNDF) removed\n",
		   t.Id,t.Type,t.Chi2tot/(t.NDFs-4));
	    goto pickup;
	  }
	  else {                   // ***** ...ELSE REMOVE BAD HITS...
	    for (int i = 0, ih = ih0; i<(int)badHits.size(); i++) {
	      THit *h = badHits[i]; int ipl = h->IPlane, kpl = ipl - ipl0;
	      int match; for (match = 0; ih<=ihLast; ih++) {
		int iref = hits[ih]; if (iref<0) continue;
		if (&vecHit[href[iref]]==h) { match = 1; break; }
	      }
	      if (!match) CsErrLog::msg(elFatal,__FILE__,__LINE__,
		"Evt #%d Track 3%d Inconsistency: Worst hit %d not in track!\n",
					event,t.Id,h->IHit);
	      t.SubHit(*h); hitsPats[kpl/32] ^= 1<<kpl%32; hits[ih] = -2-maxhit;
 	    }
	    t.QuickKF(1,0); t.QuickKF(-1,0);  // ***** ...AND REFIT
	    if (QNstatus) {
	      if (t.QNewtonFit(1,QNmode)) {
		QNchi2 = t.Chi2tot/(t.NDics-QNnh0); QNstatus = 2;
	      }
	      else QNstatus = 0;
	    }
	  }
	} // End of SI cleaning

	if (TOpt::ReMode[14] && igroup<=2) {  // ********** GROUP<=2 **********

	  if (!QNstatus || (TOpt::ReMode[19]&1<<igroup)) {

	    // ********** KF TTracks (no QN) **********

	    if (t.NDFs<4) {
	      if (t.NDFs+nIsolated<4)
		CsErrLog::msg(elError,__FILE__,__LINE__,
			      "Track %d w/ only %d DFs!",t.Id,t.NDFs);
	      ok = 0;
	    }
#ifdef CUT_KF_UPON_0
	    else if (t.NDFs>4) {
	      float chi2 = t.Chi2tot/(t.NDFs-4);
	      if (iter<=0) {   // ...if 1st iter => cut on chi2
		if (chi2>7.5) ok = 0;
	      }
	      else if (iter==1 /* && igroup==0 */) {
		if (chi2>12.5) ok = 0;
	      }
	      else if (t.Hlast(0)-t.Hfirst(0)>400 &&
		       // Consider only long tracks: otherwise might eliminate
		       // good, low P, tracks (04W28/cdr09002-37277#20974428).
		       // This cannot affect the splitting rescue infra, anyway.
		       chi2>50) ok = 0;
	      if (nScifis>=4 &&     // ***** SCIFI-TRACK: DISCARD BAD chi2 *****
		  // Scifi tracks are high momentum, and therefore not liable
		  // to deviate from straight line.
		  // Loose chi2 cut: so that this lately introduced rescue patch
		  // have no big an impact. Could be made stricter in the
		  // future. Or one could submit scifi-tracks to the straight
		  // track cleaning process, already used in zones 0x18.
		  // One should also introduce an option-driven cut value.
		  // (N.B.: (added lately) Be careful, the cut as is already
		  // throws out mc.2003.01 Evt#42 track through FI03+MM02:3 =>
		  // have a stricter definition of scifi track. Softer than
		  // TTrack::Scifi, for I don't know now what drove me not to
		  // use the latter in the first place...)
		  t.NDFs-nScifis<4 &&
		  chi2>10) ok = 0;
	      //         ********** RESCUE **********
	      if (iter>=1 && ok==0 &&                 // ***** SPLIT TRACK *****
		  (igroup==1 || igroup==2)) {
		float dXMx; static float xPRv; int iSplit, nHs;
		list<int>::iterator ihp;
		for (ihp = t.lHitPat.begin(), dXMx = 300, // I.e. >RICH thickness
		       iSplit = -1, nHs = 0; ihp!=t.lHitPat.end(); ihp++) {
		  const THit &h = vecHit[*ihp];
		  const TPlane &p = setup.vPlane(h.IPlane);
		  const TDetect &d = setup.vDetect(p.IDetRef);
		  if (++nHs>=nHitsMn && d.X(0)-xPrv>dXMx) {
		    dXMx = d.X(0)-xPrv; iSplit = h.IPlane;
		  }
		  xPrv = d.X(0);
		}
		if (iSplit>=0) {
		  t.Shorten(iSplit); if (t.QuickKF(-1,0)) {
		    t.QuickKF(1,0); t.IFit = 0x1; float KFchi2tot = t.Chi2tot;
		    if (KFchi2tot/(t.NDFs-4)<12.5 &&
			KFchi2tot/(t.NDFs-4)<chi2/2) {
		      ok = 1; straightTrack = false;
		      for (int ih = ih0; ih<nh-1; ih++) { // Update "hits" array
			int iref = hits[ih]; if (iref<0) continue;
			if (vecHit[href[iref]].IPlane<iSplit) continue;
			hits[ih] = -2-maxhit;
		      }
		      if (QNmode>-1) {
			if (t.QNewtonFit(1,QNmode)) {
			  QNchi2 = t.Chi2tot/(t.NDics-QNnh0);
			  if (QNchi2<5 || QNchi2<5*KFchi2tot/(t.NDFs-5)) {
			    QNstatus = 2; QNchi2 = t.Chi2tot/(t.NDics-QNnh0);
			  }
			  else QNstatus = 0;
			}
			else QNstatus = 0;
		      }
		    }
		  }
		}
	      }
	    }
#endif
	    if (idebug) {
	      printf("No QN: %.3f ok=%d\n",
		     t.NDFs>4 ? t.Chi2tot/(t.NDFs-4) : 0, ok); 
	      if (ok || (idebug&0x8)) t.DumpHits(0);
	    }
	  }

	  if (QNstatus && (TOpt::ReMode[19]&1<<igroup)==0) {
	    if (QNmode!=0) QNmode = 4;
	    if (QNstatus==-1) {   // => Update "QNchi2" (redundant in principle)
	      t.QNewtonFit(1,QNmode);
	      QNchi2 = t.Chi2tot/(t.NDics-QNnh0); QNstatus = 2;
	    }
	    if (idebug&0x4) { cout<<"QN\n"; t.DumpHits(1); }
	    //#define TEvPP_HISTO_QN
#ifdef TEvPP_HISTO_QN
	    if (TOpt::Hist[2]) {                  // ***** HISTO QN CHI2 *****
	      static bool first = true;
	      static CsHist2D *tQNOK, *tQNKO;
	      if(first){
		first = false;          // What about QN failure?!!
		tQNOK = new CsHist2D("QNOK","`v#^2! OK",100,0,100,2,0,2);
		tQNKO = new CsHist2D("QNKO","`v#^2! KO",100,0,100,2,0,2);
	      }
	      t.FindKine();
	      if (t.NHits==t.NHsame) tQNOK->Fill(QNchi2,igroup);
	      else                   tQNKO->Fill(QNchi2,igroup);
	    }
#endif

	    // Reference to "TLattice". In view of determining whether given
	    // plane has meaningfull guestimates. This is not very satisfying
	    // to have a refernce to "TLattice". Best would have been to have
	    // "TTrack" (which could had guestimates determined by other mean
	    // than "TLattice") return the answer. But for the time being...
	    const TLattice *lat = TLattice::Ptr();

	    int iref, ih;
	    static int   QNnh; // "static" to avoid "warning... uninitialized"
	    if      (QNchi2<30) {

	      // ******************** ALL REASONABLE CHI2's ********************
	      ok = 2;

	      for (ih = ih0, iref = 0; ih<nh-1; ih++) {
		// *************** BUILD HIT PATTERN ***************
		// *************** RAISE ALTERNATIVES AMBIGUITY ***************
		int ialt = iref;  // Remember previous hit: might be alternative
		iref = hits[ih]; if (iref<0) continue;   // Bypass alternative
		if (iref>=nHits) { // Should never happen, cf. "findMirrorRef"
		  CsErrLog::msg(elError,__FILE__,__LINE__,
		    "Alternative: Track ID=%d has reference %d>=nHits=%d",
		    t.Id,iref,nHits); continue;
		}
		THit &h = vecHit[href[iref]];
		int kpl = h.IPlane-ipl0;
		if      (lat->dico_idx[h.IPlane]==-1) continue;
		else if (lat->dico_idx[h.IPlane]==-2) break;
		hitsPats[kpl/32] |= 1<<kpl%32;
		if (ialt>=0 || ialt==-2-maxhit || // ***** ALTERNATIVE EXISTS...
		    // ...and GOOD ENOUGH CHI2 (so that the comparison of...
		    QNchi2>5)    // ...residuals, used hereafter, be meaningful)
		  continue;
		ialt = -2-ialt;
		THit *halt = &vecHit[href[ialt]], *hmalt = NULL;
		if (halt->Status==1) continue;// Hit yet used by very good track
		float guess = t.vGuests[kpl], residu = fabs(guess-h.U),
		  residualt = fabs(guess-halt->U);
		if (halt->Mirror!=-1) {  // ***** ALTERNATIVE HAS MIRROR *****
		  hmalt = &vecHit[halt->Mirror];
		  if (hmalt->Status!=-4) {
		    float residum = fabs(guess-hmalt->U);
		    if (residum<residualt) {
		      THit *htmp = halt; halt = hmalt; hmalt = htmp;
		      residualt = residum;
		      if (halt->Status>=statusMn) {
			// The mirror just retained has a reference in "href"
			ialt =findMirrorRef(ialt,hmalt->IHit,halt->IHit,href);
			if (ialt==-1) {
			  CsErrLog::msg(elError,__FILE__,__LINE__,
             "Alternative: Track ID=%d Hit %d has Status %d>=%d and yet no ref",
					t.Id,halt->IHit,halt->Status,statusMn);
			  continue;
			}
		      }
		    }
		  }
		}
		if (residualt<residu) {
		  t.SubHit(h); t.AddHit(*halt);
		  QNstatus = -1; ifl[iref]--; ifl[ialt]++; // Update flags
		  hits[ih] = ialt;        // Replace hit with alternative
		  hits[ih-1] = -2-maxhit; // Cancel alternative placeholder
		  if (halt->Status<statusMn) {
		    // Retained hit is not is the aray of references "href"
		    href[ialt] = halt->IHit; int status = halt->Status;
		    halt->status = hmalt->Status; hmalt->status = status;
		  }
		}
	      }
	      QNnh = ih;           // Save "nh" in QN range
	      if (QNstatus==-1) {
		t.QNewtonFit(1,QNmode);
		QNstatus = 2;  QNchi2 = t.Chi2tot/(t.NDics-QNnh0);
		if (idebug&0x4) { printf("Alt\n"); t.DumpHits(1); }
		t.IFit &= 0x8;       // Flag straight line fit to be updated
	      }
	    }
	    if (1.5<QNchi2 && QNchi2<30) {

	      // *************** BAD CHI2 ***************

	      int modified = 0; double hBackup[6];
	      for (int i = 1; i<=5; i++) hBackup[i] = t.Haux(i);
	      for (int kter = 0; kter<2; kter++) {
		static TStation *sPrv; static unsigned int sPat;
		for (ih = ih0, iref = 0, sPrv = 0; ih<QNnh; ih++) {
		  int ialt = iref;// Remember previous hit: might be alternative
		  iref = hits[ih]; if (iref<0) continue;   // Ref to alternative
		  if (iref>=nHits) { // Should never happen, cf. "findMirrorRef"
		    CsErrLog::msg(elError,__FILE__,__LINE__,
		      "Cleaning: Track ID=%d has reference %d>=nHits=%d",
		      t.Id,iref,nHits); continue;
		  }
		  THit &h = vecHit[href[iref]];
		  int kpl = h.IPlane-ipl0, jpl = jpls[kpl];
		  unsigned int pat = hitsPats[kpl/32];
		  if ((pat&(1<<kpl%32))==0)  // Erased in a previous "kter"...
		    continue;                // ...Or outside Dico range
		  float residu = (t.vGuests[kpl]-h.U)/res[jpl];
		  // Isolated?
		  const TPlane &p = setup.vPlane(h.IPlane);
		  const TStation *&s = p.Station;
		  if (s!=sPrv) {
		    sPrv = const_cast<TStation*&>(s);
		    int jpat = s->IPat; sPat = pat>>s->JPl;
		    if (jpat<nPats && s->JPl!=32 && s->JPl!=0)
		      sPat |= hitsPats[jpat+1]<<(32-s->JPl);
		    sPat &= s->Pat;
		  }
		  bool isolated =
		    (sPat&s->Pat)==(unsigned int)(1<<kpl%32)>>s->JPl;
		  bool gemUnCorr = false; // GEM correlation: = false <-> OK.
		  if (s->Type==26) {
		    int iGem;   // 1st or 2nd GEM 2-coordinate detector
		    if      ((sPat&0x3)==0x3) iGem = 1;
		    else if ((sPat&0xc)==0xc) iGem = 2;
		    else                      iGem = 0;
		    if (iGem) {
		      static bool unCorr;
		      if (kpl==s->IPlanes[0]-ipl0   && iGem==1 ||
			  kpl==s->IPlanes[0]-ipl0+2 && iGem==2) {
			int ihp; for (ihp = ih+1; ihp<QNnh; ihp++) {
			  int irefp = hits[ihp]; if (irefp<0) continue;
			  if (irefp>=nHits) {
			    // Should never happen, cf. "findMirrorRef"
			    CsErrLog::msg(elError,__FILE__,__LINE__,
		             "Cleaning: Track ID=%d has reference %d>=nHits=%d",
			      t.Id,irefp,nHits); continue;
			  }
			  THit &hp = vecHit[href[irefp]];
			  float ap = hp.ptrCl->getAllAnalogData()[2];
			  float a = h.ptrCl->getAllAnalogData()[2];
			  const TDetect &d = setup.vDetect(p.IDetRef);
			  CsGEMDetector* csG =
			    dynamic_cast<CsGEMDetector*>(d.PtrDet());
			  if (!csG) CsErrLog::msg(elFatal,__FILE__,__LINE__,
						  "ipl = %d = %s not a CsGEM!",
						  kpl+ipl0,d.Name.c_str());
			  const float *ampCorr = csG->getAmpCorr();
			  unCorr = fabs(ap-ampCorr[0]-a*ampCorr[1])>75;
			  break;
			}
			if (ihp==QNnh) {
			  CsErrLog::mes(elFatal,"Inconsistency!");
			}
		      }
		      gemUnCorr = unCorr;
		    }
		  }
		  bool ambiguous = ialt<0 && -2-ialt!=maxhit;
		  bool unraised = h.ptrCl->getLRProb()==.5;
		  if (kter==0) {   // ********** 1ST iter: VERY LARGE RESIDUAL
		    if (fabs(residu)<TOpt::dPRpar[igroup*10+0]*3) continue;
		  }
		  else if         // ********** 2ND iter:...
		    (fabs(residu)<                // ...LARGE RESIDUAL...
		     TOpt::dPRpar[igroup*10+0]
		     /* Allowance for fit imperfection */ * 1.75 &&
		     (fabs(residu)<               // ...OR NOT SO LARGE AND...
		      TOpt::dPRpar[igroup*10+0]*0 ||
		      ifl[iref]<=1 &&                 // ...MULTIHIT
		      !ambiguous &&                   // ...OR AMBIGUOUS
		      !isolated && !unraised &&       // ...OR ISOLATED OR LR=.5
		      !gemUnCorr)) continue;          // ...OR UNCORRELATED GEM
		  if (idebug) printf("Try %d %f %f\n",
				     kpl+ipl0,residu,TOpt::dPRpar[igroup*10+0]);
		  if (!ambiguous && // If non ambiguous need enough "NDics"
		      // to be able to subtract a hit still QN fit
		      t.NDics<=QNnh0+1) continue;
		  const TDetect &d = setup.vDetect(p.IDetRef);
		  float incidence = d.Ca*tgx+d.Sa*tgy;
		  t.SubHit(h);             // ********** SUBTRACT HIT **********
		  vector<float> vguestsp = t.vGuests;
		  if (ambiguous) {              // ***** IF ALSO AMBIGUOUS *****
		    THit &halt = vecHit[href[-2-ialt]]; t.AddHit(halt);

		    // ***** ASSOCIATE? *****
		    // We do not consider the possibilty that the mirror of
		    // the alternative is the genuine hit. But we do look
		    // for the associates of the alternative on the next
		    // plane, and decide upon which to consider on a
		    // re-evaluation of LR.
		    int jh = ih /* init search for kh */, jalt = -1;
		    THit *haj = NULL, *hmaj = NULL;
		    static int jref;  // "static" to avoid compiler's warning
		    const TPlane *pj = p.Associate; static int ijpl, jjpl;
		    if (pj) {
		      ijpl = p.Associate->IPlane-ipl0; jjpl = jpls[ijpl];
		    }
		    CsCluster *ci = h.ptrCl; static CsCluster *cj, *ai, *aj, *al;
		    list<CsCluster*> associates = ci->getAssociateClusters();
		    if (!associates.empty() &&
			(hitsPats[ijpl/32]&(1<<ijpl%32))) {
		      CsCluster *cik = associates.front(),
			*cil = associates.back();
		      for (jh = ih+1; jh<ih+3 && jh<QNnh; jh++) {
			jref = hits[jh]; if (jref<0) continue; // Bypass alter
			if (jref>=nHits) {
			  // Should never happen, cf. "findMirrorRef"
			  CsErrLog::msg(elError,__FILE__,__LINE__,
		            "Cleaning: Track ID=%d has reference %d>=nHits=%d",
		            t.Id,jref,nHits); continue;
			}
			THit &hj = vecHit[href[jref]]; cj = hj.ptrCl;
			if (cj!=cik && cj!=cil) break;
			ai = halt.ptrCl; aj = NULL;
			list<CsCluster*> associates = ai->getAssociateClusters();
			if (associates.empty()) break;
			CsCluster *aik = associates.front(),
			  *ail = associates.back();
			for (int ir = fh[jjpl]; ir<lh[jjpl]; ir++) {
			  haj = &vecHit[href[ir]];
			  if      (haj->ptrCl==aik) aj = aik;
			  else if (haj->ptrCl==ail) aj = ail;
			  if (aj) {
			    CsEventUtils::setLRToMatchAngle(ai,aj,incidence);
			    if (halt.Mirror!=-1) {
			      THit &hmi = vecHit[halt.Mirror];
			      if (hmi.Status==-4) {
				// Restore LR probas that have just been upset.
				ai->setLRProb(1); hmi.ptrCl->setLRProb(0); 
			      }
			    }
			    if (haj->Mirror!=-1) {
			      THit &hmj = vecHit[haj->Mirror];
			      if (hmj.Status==-4) {
				// Restore LR probas that have just been upset.
				aj->setLRProb(1); hmj.ptrCl->setLRProb(0); 
			      }
			    }
			    if (ai->getLRProb()>.5) { // I.e. require ...
			      // alternative to be LR most probable
			      if (aj->getLRProb()<.5) {
				if (haj->Mirror!=-1) {
				  hmaj = &vecHit[haj->Mirror];
				  if (hmaj->Status!=-4) {
				    THit *htmp = haj; haj = hmaj; hmaj = htmp;
				  }
				  else hmaj = NULL;
				}
				else {
				  CsErrLog::msg(elError,__FILE__,__LINE__,
			 "THit %d,%d,%.3f: associated CsClusters & no Mirror",
						haj->IHit,haj->IPlane,haj->U);
				  haj = NULL;
				}
			      }
			      if (haj) {
				if (hmaj && haj->Status>=statusMn) {
				  // Mirror retained has a reference in "href"
				  jalt =
				    findMirrorRef(ir,hmaj->IHit,haj->IHit,href);
				  if (jalt==-1) {
				    CsErrLog::msg(elError,__FILE__,__LINE__,
                "Cleaning: Track ID=%d Hit %d has Status %d>=%d and yet no ref",
				      t.Id,haj->IHit,haj->Status,statusMn);
				    haj = NULL;
				  }
				}
				else  jalt = ir;
			      }
			      if (haj && haj->IHit!=hj.IHit) {
				t.SubHit(hj); t.AddHit(*haj);
			      }
			    }
			    break;
			  }  // End look for associates of alternative hit
			}
			break;   // At least one non alternative encountered
		      }
		    }  // End look for associates of current hit on next plane

		    //                             ***** CHECK ALTERNATIVE *****
		    if (t.QNewtonFit(1,QNmode)) {
		      if (idebug&0x4) t.DumpHits(1);
		      float chi2 = t.Chi2tot/(t.NDics-QNnh0);
		      if (chi2<QNchi2*.8) {
			//                      ***** OK: VALIDATE CHANGES *****
			QNstatus = 2; QNchi2 = chi2; modified = 1;
			for (int i = 1; i<=5; i++) hBackup[i] = t.Haux(i);
			ifl[iref]--; ifl[-2-ialt]++;  // hit<->alternative
			hits[ih] = -2-ialt; hits[ih-1] = -2-maxhit;
			if (jalt>=0) {                // Next plane...
			  ifl[jref]--; ifl[jalt]++;   // ...hit<->alternative
			  hits[jh] = jalt; hits[jh-1] = -2-maxhit;
			  if (haj->Status<statusMn) {
			    // In case not yet imported...
			    if (!hmaj)
			      CsErrLog::msg(elFatal,__FILE__,__LINE__,
			 "THit %d,%d,%.3f has Status %d (<%d), yet no Mirror",
					    haj->IHit,haj->IPlane,haj->U,
					    haj->Status,statusMn);
			    href[jalt] = haj->IHit; int status = haj->Status;
			    haj->status = hmaj->Status; hmaj->status = status;
			  }
			}
			continue;        // ***** IF VALIDATED => continue *****
		      } // End ambiguous case I
		    }

		    // ***** NEXT ANALOGUE PLANE IN SAME(?) STATION *****

		    int kh, kalt = -1, lalt = -1;
		    static int kref, mh, lref; // "static" to avoid compiler's warning
		    THit *hal = NULL, *hmal = NULL;
		    for (kh = jh+1; kh<QNnh; kh++) {
		      kref = hits[kh];
		      int malt = -2-kref;  // Remember previous hit=?alternative
		      if (kref<0) continue;        // Bypass alternative
		      if (kref>=nHits) {
			// Should never happen, cf. "findMirrorRef"
			CsErrLog::msg(elError,__FILE__,__LINE__,
			  "Cleaning: Track ID=%d has reference %d>=nHits=%d",
			  t.Id,kref,nHits); continue;
		      }
		      THit &hk = vecHit[href[kref]];
		      const TPlane &pk = setup.vPlane(hk.IPlane);
		      //if (pk.Station!=s) break;
		      if (pk.IProj!=p.IProj) continue;
		      if (malt<0 || malt==maxhit) break;
		      THit &hak = vecHit[href[malt]];
		      int ikpl = hk.IPlane-ipl0, jkpl = jpls[ikpl];
		      if ((hitsPats[ikpl/32]&(1<<ikpl%32))==0)
			break;  // Must have been erased in a previous kter
		      float residuk   = (vguestsp[ikpl]-hk.U   )/res[jkpl];
		      float residualt = (vguestsp[kpl] -halt.U )/res[jpl];
		      float residuak  = (vguestsp[ikpl]-hak.U)/res[jkpl];
		      if (residuk*residu<0 || residualt*residuak<0) break;
		      t.SubHit(hk); t.AddHit(hak); kalt = malt;

		      // ***** ASSOCIATE? *****
		      const TPlane *pl = pk.Associate; static int ilpl, jlpl;
		      if (pl) {
			ilpl = pk.Associate->IPlane-ipl0; jlpl = jpls[ilpl];
		      }
		      CsCluster *ck = hk.ptrCl, *cl, *ak;
		      list<CsCluster*> associates = ck->getAssociateClusters();
		      if (!associates.empty() &&
			  (hitsPats[ilpl/32]&(1<<ilpl%32))) {
			CsCluster *ckk = associates.front(),
			  *ckl = associates.back();
			for (mh = kh+1; mh<kh+3 && mh<QNnh; mh++) {
			  lref = hits[mh]; if (lref<0) continue; // Bypass alt
			  if (lref>=nHits) {
			    // Should never happen, cf. "findMirrorRef"
			    CsErrLog::msg(elError,__FILE__,__LINE__,
			     "Cleaning: Track ID=%d has reference %d>=nHits=%d",
		              t.Id,lref,nHits); continue;
			  }
			  THit& hl = vecHit[href[lref]]; cl = hl.ptrCl;
			  if (cl!=ckk && cl!=ckl) break;
			  ak = hak.ptrCl; al = NULL;
			  list<CsCluster*> associates = ak->getAssociateClusters();
			  if (associates.empty()) break;
			  CsCluster *akk = associates.front(),
			    *akl = associates.back();
			  for (int ir = fh[jlpl]; ir<lh[jlpl]; ir++) {
			    hal = &vecHit[href[ir]];
			    if      (hal->ptrCl==akk) al = akk;
			    else if (hal->ptrCl==akl) al = akl;
			    if (al) {
			      CsEventUtils::setLRToMatchAngle(ak,al,incidence);
			      if (ak->getLRProb()>.5) {
				if (al->getLRProb()<.5) {
				  if (hal->Mirror!=-1) {
				    hmal = &vecHit[hal->Mirror];
				    if (hmal->Status!=-4) {
				      THit *htmp = hal; hal = hmal; hmal = htmp;
				    }
				    else hmal = NULL;
				  }
				  else {
				    CsErrLog::msg(elError,__FILE__,__LINE__,
			 "THit %d,%d,%.3f: associated CsClusters & no Mirror",
						  hal->IHit,hal->IPlane,hal->U);
				    hal = NULL;
				  }
				  if (hal) {
				    if (hmal && hal->Status>=statusMn) {
				      // Mirror retained w/ reference in "href"
				      lalt =
					findMirrorRef(ir,hmal->IHit,hal->IHit,href);
				      if (lalt==-1) {
					CsErrLog::msg(elError,__FILE__,__LINE__,
                "Cleaning: Track ID=%d Hit %d has Status %d>=%d and yet no ref",
				          t.Id,hal->IHit,hal->Status,statusMn);
					hal = NULL;
				      }
				    }
				    else lalt = ir;
				  }
				  if (hal && hal->IHit!=hl.IHit) {
				    t.SubHit(hl); t.AddHit(*hal);
				  }
				}
			      }
			      break;
			    }
			  }
			  break;   // At least non alternative encountered
			}
		      }
		      break;
		    } // End look for analogue plane
		    //      ***** CHECK ALTERNATIVE on NEXT ANALOGUE PLANE *****
		    float chi2; if (kalt==-1 || t.QNewtonFit(1,QNmode)) {
		      if (idebug&0x4) t.DumpHits(1);
		      chi2 = t.Chi2tot/(t.NDics-QNnh0);
		    }
		    else chi2 = 100;
		    if (chi2<QNchi2*.8) {
		      // ***** OK: VALIDATE CHANGES *****
		      QNstatus = 2; QNchi2 = chi2; modified = 1;
		      for (int i = 1; i<=5; i++) hBackup[i] = t.Haux(i);
		      ifl[iref]--; ifl[-2-ialt]++; // hit<->alternative
		      hits[ih] = -2-ialt; hits[ih-1] = -2-maxhit;
		      if (jalt>=0) {               // Next plane...
			ifl[jref]--; ifl[jalt]++;  // ...hit<-alternative
			hits[jh] = jalt; hits[jh-1] = -2-maxhit;
			if (haj->Status<statusMn) {
			  // In case not yet imported...
			  if (!hmaj)
			    CsErrLog::msg(elFatal,__FILE__,__LINE__,
			 "THit %d,%d,%.3f has Status %d (<%d), yet no Mirror",
					  haj->IHit,haj->IPlane,haj->U,
					  haj->Status,statusMn);
			  href[jalt] = haj->IHit; int status = haj->Status;
			  haj->status = hmaj->Status; hmaj->status = status;
			}
		      }
		      if (kalt>=0) {               // Next analogue plane...
			ifl[kref]--; ifl[kalt]++;  // ...hit<-alternative
			hits[kh] = kalt; hits[kh-1] = -2-maxhit;
			if (lalt>=0) {               // NNanalogue plane...
			  ifl[lref]--; ifl[lalt]++;  // ...hit<-alternative
			  hits[mh] = lalt; hits[mh-1] = -2-maxhit;
			  if (hal->Status<statusMn) {
			    // In case not yet imported...
			    if (!hmal)
			      CsErrLog::msg(elFatal,__FILE__,__LINE__,
			 "THit %d,%d,%.3f has Status %d (<%d), yet no Mirror",
					    hal->IHit,hal->IPlane,hal->U,
					    hal->Status,statusMn);
			    href[lalt] = hal->IHit; int status = hal->Status;
			    hal->status = hmal->Status; hmal->status = status;
			  }
			}
		      }
		      continue;          // ***** IF VALIDATED => continue *****
		    }
		    else {             // ***** !OK: RESTORE INITIAL TRACK *****
		      t.SubHit(halt); QNstatus = 1; 
		      t.vGuests = vguestsp;
		      for (int i = 1; i<=5; i++) t.Haux(i) = hBackup[i];
		      if (jalt>=0) {
			t.SubHit(*haj);
			THit &hj    = vecHit[href[jref]]; t.AddHit(hj);
		      }
		      if (kalt>=0) {
			THit &hak = vecHit[href[kalt]]; t.SubHit(hak);
			THit &hk  = vecHit[href[kref]]; t.AddHit(hk);
			if (lalt>=0) {
			  t.SubHit(*hal);
			  THit &hl = vecHit[href[lref]]; t.AddHit(hl);
			}
		      }
		    }  // End !converging QN fit in ambiguous case II
		  } // End ambiguous case

		  if (t.NDics<=QNnh0+1) { // ***** IF NOT ENOUGH D-HITS... *****
		    t.AddHit(h); t.NDics++;      // ... RESTORE HIT
		    QNstatus = 1; // Not needed if QN fit failed: too bad...
		    continue;                    // ... GIVE UP
		  }

		  if (t.QNewtonFit(1,QNmode)) {

		    // ***** EVALUATE NOW W/ CURRENT HIT SIMPLY SUBTRACTED *****

		    float chi2 = t.Chi2tot, chi2Old = QNchi2*(t.NDics-QNnh0+1);
		    if (idebug&0x4) t.DumpHits(1);
		    if (chi2<(chi2Old-residu*residu)*.8 || chi2<chi2Old*.6) {
		      if ((int)t.NDFs<nHitsMn) {
			// Not enough hits any more => Flag track to be released
			ok = 0; goto next_step;
		      }
		      if (unraised) {
			THit *hm = NULL, *hp;
			CsCluster *mirror = h.ptrCl->getMirrorCluster();
			int jref; float residm;
			if ((jref = iref+1)<nHits) {
			  hp = &vecHit[href[jref]];
			  if (hp->ptrCl==mirror) hm = hp;
			}
			else if ((jref = iref-1)>=0) {
			  hp = &vecHit[href[jref]];
			  if (hp->ptrCl==mirror) hm = hp;
			}
			if (hm &&
			    (residm = fabs((t.vGuests[kpl]-hm->U)/res[jpl]))<
			    TOpt::dPRpar[igroup*10+0]) {
			  t.AddHit(*hm);
			  QNstatus = 1; ifl[jref]++;  // Update flags
			  hits[ih] = iref;         // Replace hit w/ alternative
			  chi2 += residm*residm;       // Approximate chi2
			}
			else {
			  QNstatus = 2; hitsPats[kpl/32] ^= 1<<kpl%32;
			}
		      }
		      else {
			QNstatus = 2; hitsPats[kpl/32] ^= 1<<kpl%32;
		      }
		      ifl[iref]--; QNchi2 = chi2/(t.NDics-QNnh0); modified = 1;
		      if (QNchi2<2)  goto next_step;
		    }
		    else {  // ***** !OK: RESTORE TRACK and HIT *****
		      t.AddHit(h); hitsPats[kpl/32] |= 1<<kpl%32; t.NDics++;
		      QNstatus = 1; t.vGuests = vguestsp;
		      for (int i = 1; i<=5; i++) t.Haux(i) = hBackup[i];
		    }
		  }
		  else {  // ***** QN FAILED: RESTORE HIT *****
		    t.AddHit(h); hitsPats[kpl/32] |= 1<<kpl%32; t.NDics++;
		  }
		}  // End loop on hits
	      } // End of 2 iterations
	    next_step:
	      if (modified && (t.IFit&0x1)) {
		CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "TTrack %d modified by QN cleaning has IFit = %d",t.Id,t.IFit);
		//t.IFit &= 0x8;         // Flag straight line fit to be updated
	      }
              if (ok && QNchi2<20.) {
		if (QNstatus!=2) {         // Update fit guestimates...
		  t.QNewtonFit(1,QNmode);
		  QNchi2 = t.Chi2tot/(t.NDics-QNnh0); QNstatus = 2;
		  if (idebug&0x4) { printf("Clean\n"); t.DumpHits(1); }
		  if (1.75<QNchi2 && QNchi2<5) { 
		    for (ih = ih0, iref = 0; ih<nh-1; ih++) {
		      int ialt = iref;  // Remember previous hit: might be alter
		      iref = hits[ih]; if (iref<0) continue;     // Bypass alter
		      if (iref>=nHits) {
			// Should never happen, cf. "findMirrorRef"
			CsErrLog::msg(elError,__FILE__,__LINE__,
 "Cleaning: Track ID=%d has reference %d>=nHits=%d",t.Id,iref,nHits); continue;
		      }
		      THit &h = vecHit[href[iref]];
		      int kpl = h.IPlane-ipl0;
		      if ((hitsPats[kpl/32]&(1<<kpl%32))==0) continue;
		      if      (lat->dico_idx[h.IPlane]==-1) continue;
		      else if (lat->dico_idx[h.IPlane]==-2) break;
		      if (ialt<0 && ialt!=-2-maxhit) {   // Alternative exists...
			ialt = -2-ialt;
			if (ialt==maxhit) continue;
			THit &halt = vecHit[href[ialt]];
			float residu    = t.vGuests[kpl]-h.U;
			float residualt = t.vGuests[kpl]-halt.U;
			if (fabs(residualt)<fabs(residu)) {
			  t.SubHit(h); t.AddHit(halt);
			  QNstatus = -1; ifl[iref]--; ifl[ialt]++; // Update flags
			  hits[ih] = ialt;      // Replace hit with alternative
			}
		      }
		    }
		    if (QNstatus!=2) {
		      t.QNewtonFit(1,QNmode);
		      QNchi2 = t.Chi2tot/(t.NDics-QNnh0); QNstatus = 2;
		      if (idebug&0x4) { printf("Clean\n"); t.DumpHits(1); }
		    }
		  }
		}
              }
	      else ok = 0;
	      if (modified && (int)t.NDFs==nHitsMn) {
		// ***** ENOUGH HITS? CASE of TTrack's W/ DRIFT HITS *****
		// This amounts to enforcing the ``drifts subdued'' at play in
		// the TAlgo2 methods: there drift hits were given less weight
		// and tracks w/ drift hits making it through these methods
		// must have had more than what's required by iPRpar.
		bool hasDrift = false;
		for (int ipat = 0; ipat<nPats; ipat++)
		  if (hitsPats[ipat]&driftsPats[ipat]) {
		    hasDrift = true; break;
		  }
		if (hasDrift) ok = 0;
	      }
	    }
	    if (ok) {
	      if (!TOpt::dCut[6]) {
		if (QNchi2>=18) ok = 0;
	      }
	      else {
		int oldOK = ok;
		if (igroup==0) {
		  if (QNchi2>=TOpt::dCut[6]) ok = 0;
		}
		else if (igroup==1) {
		  if (t.Hlast(0)<750) {
		    if (QNchi2>=TOpt::dCut[6]) ok = 0;
		  }
		  else if (QNchi2>=15) ok = 0;
		}
		else if (QNchi2>=15) ok = 0;
		if (!ok && QNchi2<20) {
		  int nDrifts, ipat;
		  for (ipat=nDrifts = 0; ipat<nPats; ipat++) {
		    unsigned int pat = hitsPats[ipat]&driftsPats[ipat];
		    for (int ibit = 0; ibit<32; ibit++)
		      if (pat&1<<ibit) nDrifts++;
		  }
		  if (fabs(t.Hlast(3))<.01 && fabs(t.Hlast(4))<.01 &&
		      // Special case of paraxial tracks: could be halo mus, and
		      // then off-time => Do they consist mostly of drift hits.
		      nDrifts>t.NDFs*1/3. && QNchi2<20) ok = oldOK;
		  if (igroup==0 &&
		      // In zone 0x1, consider the case of DC00+DC01 tracks:
		      // these are difficult to track => For the time being, let
		      // us decide to accept them, so that their hits are picked
		      // up by other candidate tracks.
		      nDrifts==(int)t.NDFs && QNchi2<18) ok = oldOK;
		}
	      }
	    }
	    if (idebug) { printf("QN Cleaning done ok = %d!\n",ok); t.DumpHits(1); }
	  }
	}
	if (TOpt::ReMode[14] && (igroup>=3 || straightTrack && ok)) {
	  // *******************************************************************
	  // ********** BEAM, muW and ALL A PRIORI STRAIGHT TRACKS... **********
	  // *******************************************************************
	  float chi2; int jter;
	  for (jter = 0, chi2 = t.Chi2tot/(t.NDFs-4); jter<2; jter++) {
	    if (chi2<TOpt::dCut[2]*.6 &&     // ...BAD ENOUGH CHI2: TRY CLEANING
		(chi2<TOpt::dCut[2]*.4 || (int)t.NDFs>nHitsMn)) break;
	    double x0 = t.Hfirst(0), y0 = t.Hfirst(1), z0 = t.Hfirst(2),
	      yp = t.Hfirst(3), zp = t.Hfirst(4);
	    int worstRefs[4], ih, iworst;
	    float worstResids[4], residCut = TOpt::dPRpar[igroup*10+0]*1.25;
	    for (int iworst = 0; iworst<4; iworst++) {
	      worstResids[iworst] = residCut; worstRefs[iworst] = -1;
	    }
	    for (ih = ih0; ih<nh; ih++) {
	      int iref = hits[ih]; if (iref<0) continue; // Bypass alternatives
	      if (iref>=nHits) { // Should never happen, cf. "findMirrorRef"
		CsErrLog::msg(elError,__FILE__,__LINE__,
 "Cleaning: Track ID=%d has reference to %d>=nHits=%d",t.Id,iref,nHits);
		continue;
	      }
	      THit &h = vecHit[href[iref]];
	      int kpl = h.IPlane-ipl0, jpl = jpls[kpl]; float x = xpl[jpl]-x0;
	      float resid = ((y0+yp*x)*cosa[jpl]+(z0+zp*x)*sina[jpl]-h.U)/h.SigU;
	      //#define PR_DEBUG_STRAIGHT_CLEAN
#ifdef PR_DEBUG_STRAIGHT_CLEAN
	      printf(" %d %.2f",jpl,resid);
#endif
	      for (iworst = 0; iworst<4; iworst++) {
		if (fabs(resid)>worstResids[iworst]) {
		  for (int j = 3; j>=iworst+1; j--) {
		    worstResids[j] = worstResids[j-1];
		    worstRefs[j] = worstRefs[j-1];
		  }
		  worstResids[iworst] = fabs(resid); worstRefs[iworst] = ih;
		  break;
		}
	      }
	    }
#ifdef PR_DEBUG_STRAIGHT_CLEAN
	    printf("\n");
#endif
	    int modified, ref; double bestChi2;
	    for (iworst=modified = 0, ref = -1, bestChi2 = 100; iworst<4;
		 iworst++) {
	      ih = worstRefs[iworst]; if (ih<0) continue;  // ***** LARGE RESIDU
	      int iref = hits[ih]; THit &h = vecHit[href[iref]];
	      t.SubHit(h); t.QuickKF(-1,0);                // => REFIT W/O IT
	      double newChi2 = t.Chi2tot/(t.NDFs-4);
	      if (newChi2<chi2*.75) {                  // ***** IF MPROVEMENT...
		if (newChi2<.5*TOpt::dCut[2]) {
		  modified = 1; ref = ih;      // ...good enough: ratify modif
		  bestChi2 = newChi2; break;
		}
		if (newChi2<bestChi2) {
		  bestChi2 = newChi2; ref = ih;// ...else get max chi2 increment
		}
	      }
	      t.AddHit(h);
	    } // End loop on 4 worst residuals
	    if (ref>=0) { // ***** WORST CHI2 INCREMENT (yet or to be) DISCARDED
	      int iref = hits[ref]; THit &h = vecHit[href[iref]];
	      if (!modified) { t.SubHit(h); t.QuickKF(-1,0); modified = 1; }
	      if ((int)t.NDFs<nHitsMn) {     // If too few hits left over...
		t.UseHitTime(); float dX = dFirst->X(0)-t.Hfirst(0);
		float yFirst = t.Hfirst(1)+dX*t.Hfirst(3);
		float zFirst = t.Hfirst(2)+dX*t.Hfirst(4);
		bool inScifiSi = dFirst->InActive(yFirst,zFirst);
		if ((int)t.NDFs<nHitsMn-1 ||   // ...and much too few...
		    bestChi2>2 ||      // ...or not so good chi2...
		    nScifis>=3 &&      // ...or (if scifis) too bad timing
		    t.DispTime>1.8 ||  // (3sigmas, assuming scifi.Resol=.6)
		    inScifiSi) {       // ...or fully in scifi/Si telescope...
		  ok = 0; break;               // ...Flag track to be released
		}
	      }
	      chi2 = bestChi2; hits[ref] = -2-maxhit; ifl[iref]--;
	      t.QuickKF(1,0);
	      const TPlane &p = setup.vPlane(h.IPlane);
	      const TDetect &d = setup.vDetect(p.IDetRef);
	      if (d.IType==22) {
		nScifis--;
		t.Scifi = nScifis>=(int)t.NDFs-1?1<<igroup:0; // Useful? "nScifis" and "NDFs" being decremented simultaneously, the inequality does not change
	      }
	      if (d.IType==28 || d.IType==29) {
		nPixGPs--; if (d.IType==29) nPixGPs--;
	      }
	      if (d.IType==21) nSIs--;
	    }
	    else t.QuickKF(-1,0);                          // ***** ELSE RESTORE
	    if (!modified) break;
	  } // End loop over 2 cleaning iterations
	  if (ok && TOpt::dCut[2] &&
	      (chi2>TOpt::dCut[2]*2 ||
	       t.Hlast(0)-t.Hfirst(0)<900 && chi2>TOpt::dCut[2]) &&
	      (igroup==4 ||                        // Beam tracks: standalone...
	       straightTrack))                     // ...straight scifi track...
	    ok = 0;                                     // ***** ...CUT BAD CHI2
	  if (idebug) {
	    t.DumpHits(0); if (!ok) printf("Erased...\n");
	  }
	}
      }  // End KF Fit OK
      else ok = 0;               // No KF Fit => Flag track to be released

    pickup:
      Traffic::Ref().Stopwatch.Stop(7);
      if (ok) {
	if (ok==2 &&                                 // ***** QN Fitted...
	    t.Chi2tot/(t.NDics-4)<TOpt::dCut[60]) {  // ...GOOD CHI2...

	  //       *************** PICK-UP HITS ***************

	  // Reference to "TLattice". In view of determining whether given plane
	  // has meaningfull guestimates. This is not very satisfying to have a
	  // reference to "TLattice". Best would have been to have "TTrack" (
	  // which could have guestimates determined by other mean than
	  // "TLattice") return the answer. But for the time being...
	  const TLattice *lat = TLattice::Ptr();    

	  int added = 0;
	  for (int jpl = 0; jpl<npl; jpl++) {
	    int ipl = ipls[jpl];
	    if      (lat->dico_idx[ipl]==-1) continue;
	    else if (lat->dico_idx[ipl]==-2) break;
	    int kpl = ipl-ipl0; if ((hitsPats[kpl/32]&(1<<kpl%32))==0) {
	      // If no hit yet associated with this track on current plane
	      const TPlane &p = setup.vPlane(ipl);
	      const TDetect &d = setup.vDetect(p.IDetRef);
	      if (d.IType==32 && !(p.IFlag&0xc0)) // Isolated MPM...
		continue; //...=> No predicting v-coord => Exclude.
	      if (igroup==1) {
		if (d.X(0)-t.Hlast(0)>350)// No extrap. over full length of RICH
		// (Note: supra rejects also any extrap. over large distance,
		// which may not be a good move, e.g. if high momentum, although
		// high momentum implies easy reco and no need for pickup step.
		  break;
		if (d.IType==39 || d.IType==40) continue;   // Exclude MAs
		if (TOpt::dCut[67]<t.Hlast(0) && t.Hlast(0)<TOpt::dCut[68]) {
		  // Zone #1: Special case of extrapolation over full
		  // length of RICH while track is also defined in ST04
		  // which is not in the dictionary => Exclude PS01/GM04
		  if (750<d.X(0) && d.X(0)<980) continue;
		}
	      }
	      float yy = t.Hfirst(1)+(xpl[jpl]-t.Hfirst(0))*t.Hfirst(3);
	      float zz = t.Hfirst(2)+(xpl[jpl]-t.Hfirst(0))*t.Hfirst(4);
	      if (!setup.vDetect(idpl[jpl]).InActive(yy,zz)) continue;
	      float dUBest = res[jpl]*TOpt::dPRpar[igroup*10+0]* // ...
		// ...times allowance for fit imperfection. Note that the extra
		// route width ("dPRpar[igroup*10+2/3]") is instead disregarded.
		1.75;
	      float dVMax =// For the v coord (of pixel planes): no allowance...
		// (...this, to be on the safe side: this, built-in,
		// reconstruction parameter hasn't been optimised in any way.)
		res[jpl]*TOpt::dPRpar[igroup*10+0];
	      int ir, irBest;
	      for (ir = fh[jpl], irBest = -1; ir<lh[jpl]; ir++) {
		if (ifl[ir]) continue;
		THit &h = vecHit[href[ir]];
		float dU = fabs(t.vGuests[kpl]-h.U);
		if (dU<dUBest) {
		  // Special case of pixel planes: we decide to make it simple:
		  // the v coordinate is not considered in the determination of
		  // best pixel hit: it's only checked that it's OK.
		  if      (p.IFlag&0x50) {  // P- or M-pixel planes
		    // 0x1=GPXY: Y < X < P. 0x4=MPMY|U: X|V < Y|U < MY|U
		    if (fabs(t.vGuests[kpl-2]-h.V)<dVMax) {
		      dUBest = dU; irBest = ir;
		    }
		  }
		  else if (p.IFlag&0xa0) {  // P- or M-pixel planes
		    // 0x2=GPUV: P < U < V. 0x8=MPMX|V: MX|V < X|V < Y|U
		    if (fabs(t.vGuests[kpl-2]-h.V)<dVMax) {
		      dUBest = dU; irBest = ir;
		    }		    
		  }
		  else {
		    dUBest = dU; irBest = ir;
		  }
		}
	      }
	      if (irBest>=0) {
		t.AddHit(vecHit[href[irBest]]); added++; ifl[irBest]++;
	      }
	    }
	  }
	  if (added) {
	    t.QuickKF(-1,0); t.QuickKF(1,0);
	    t.QNewtonFit(1,QNmode);
	    if (idebug) { cout<<"Pick-up\n"; t.DumpHits(0); }
	  }
	}
	else if (ok==1 &&                            // ***** KF Fitted...
		 //                                     ... w/ GOOD chi2
		 // (Note: special (distinct from "dCut[60]" supra) chi2 cut.)
		 t.Chi2tot/(t.NDFs-4)<TOpt::dCut[87] &&
		 //  The pickup procedure is temporarily (until its impact is
		 // evaluated) restricted to zones 0x2 and  0x10.
		 //  Zone 0x10: further restricted to 2nd iteration, when not
		 // enough of scifi hits to pass "iCut[20]".
		 (igroup==4 && iter>0 && nScifis<TOpt::iCut[20] ||
		  //  Zone 0x2: also resticted, so as to focus on tracks in the
		  // VSAT region, lacking FI05, either X or Y. These tracks are
		  // poorly constrained (only constrained over the short
		  // interval [GP02,FI06] while spanning much more) and it's
		  // important to extend to FI05 as much as possible. They are
		  // good candidate for the pickup exercise, since they are
		  // expected to be straight.
		  //  As long as the impact has only been evaluated on 2010
		  // data, condition the whole thing by "iCut[31]". And iter<=0,
		  // else, w/ full scifi enhancement, all scifis are used.
		  //  Only resonably well acertained tracks are considered:
		  // #NDFs>7 (out of 4*GP+5*FI in 2010 setup).
		  igroup==1 && nScifis>3 && TOpt::iCut[31] && iter<=0 && t.NDFs>7 &&
		  (t.Hfirst(0)>750 || // Starting downstream of RICH or...
		   // ...upstream of it, at a scifi, that can then only be FI05
		   // (we want to avoid upsetting tracks w/ hits in 
		   setup.vDetect(t.lPRef().front()).IType==22))) {
	  if (igroup==1) // In the 0x2 case, since many scifi hits, cf. supra,
	    // hence good timing reliability, we'll require time consistency.
	    // (In the 0x10 case, good timing can probably achieved thanks to
	    // to SI, but this requires a good time calibration, which may not
	    // be systematically available: one would have then to condition
	    // time consitency check w/ an option. Not yet done.)
	    t.UseHitTime();
	  double T = t.MeanTime, dT;
	  unsigned int newPats[nPats]; for (int ipat = 0; ipat<nPats; ipat++)
	    newPats[ipat] = 0;
	  double oldChi2 = t.Chi2tot/(t.NDFs-4);
	  int added = 0; for (int jpl = 0; jpl<npl; jpl++) {
	    int ipl = ipls[jpl];
	    int kpl = ipl-ipl0; if (hitsPats[kpl/32]&1<<kpl%32) continue;
	    // No hit associated with this track on current plane
	    const TPlane &p = setup.vPlane(ipl);
	    const TDetect &d = setup.vDetect(p.IDetRef);
	    if (d.IType!=22) continue;
	    float yy = t.Hfirst(1)+(xpl[jpl]-t.Hfirst(0))*t.Hfirst(3);
	    float zz = t.Hfirst(2)+(xpl[jpl]-t.Hfirst(0))*t.Hfirst(4);
	    if (!setup.vDetect(idpl[jpl]).InActive(yy,zz)) continue;
	    float dUBest = res[jpl]*TOpt::dPRpar[igroup*10+0]* //...
		// ...times allowance for fit imperfection. Note that the extra
		// route width ("dPRpar[igroup*10+2/3]") is instead disregarded.
	      1.75;
	    float u = yy*d.Ca+zz*d.Sa, umn = u-2*dUBest, umx = u+2*dUBest;
	    int ir, irBest, mult; for (ir = fh[jpl], irBest = -1, mult = 0;
				       ir<lh[jpl]; ir++) {
	      if (ifl[ir]) continue;
	      THit &h = vecHit[href[ir]];
	      if (igroup==1) {
		dT = // Asumming .5 ns resolution (typical for German scifis).
		  fabs(h.Time-T)/sqrt(.25+h.SigT*h.SigT);
		if (dT>5) continue; // Very bad: skip incrementing "mult".
	      }
	      else dT = 0;
	      if (h.U<umn) continue; if (h.U>umx) break; mult++;
	      if (dT>3)// Mediocre timing: do not retain, but account in "mult".
		continue;
	      float dU = fabs(u-h.U); if (dU<dUBest) {
		dUBest = dU; irBest = ir;
	      }
	    }
	    if (irBest>=0 && mult==1) {
	      t.AddHit(vecHit[href[irBest]]); added++; ifl[irBest]++;
	      nScifis++; newPats[kpl/32] |= 1<<kpl%32;
	    }
	  }
	  if (added) {
	    t.QuickKF(-1,0); double newChi2 = t.Chi2tot/(t.NDFs-4);
	    if (newChi2>2 && newChi2>oldChi2*1.25) {
	      int ipat, kpl; for (ipat=kpl = 0; ipat<nPats; ipat++) {
		for (int j = 0; j<32; j++, kpl++) if (newPats[ipat]&1<<j) {
		  t.SubHit(kpl+ipl0); nScifis--;
		}
	      }
	      t.QuickKF(-1,0);
	    }
	    else {
	      t.QuickKF(1,0);
	      if (idebug) { cout<<"Pick-up\n"; t.DumpHits(0); }
	    }
	  }
	}
	// New that pick-up is performed, re-examine cut on # of scifi hits
	if (igroup==4 && nScifis<TOpt::iCut[20]) ok = 0;
      }
      if (ok && igroup==1 && t.Hfirst(0)<750 && t.Hlast(0)>TOpt::dCut[68]) {
	//     ********** LONG TRACKS in ZONE 0x2 **********
	// There are so many detectors in zone 0x2 that tracks can easily pick
	// up wrong hits. A particular example is that of tracks in the sub-zone
	// SM1<->RICH that are unduly extended beyond the RICH. The case of an
	// extension limited to the RICHwall complex (GM04+PS+DR/ST04) has been
	// treated supra. Let's tackle here the case of tracks extending beyond
	// the RICHwall complex, by requiring a reasonable detector efficiency
	// over the length of the extension downsream of RICH, restricting
	// ourselves to SAT and LAT, and excluding the trickier VSAT (where
	// anyway there is not so much redundancy that it can become a problem).
	float nExpect; int nFound, nVSATs, nMAs, iplF;
	float xL = t.Hlast(0), yL = t.Hlast(1), yp = t.Hlast(3);
	float                  zL = t.Hlast(2), zp = t.Hlast(4);
	int iplL = t.lPlnRef.back(); const TPlane &pL = setup.vPlane(iplL);
	int istL = pL.Station->IStation, ist;
	const vector<TStation> &stations = setup.vStation();
	list<int>::reverse_iterator rp = t.lPlnRef.rbegin();
	for (ist = istL, nExpect = 0, nFound=nVSATs=nMAs = 0, iplF = maxpl;
	     ist>=0; ist--) {
	  //     ***** LOOP on STATIONS in SUB-ZONE SM2... *****
	  // - Starting w/ station of last associated hit.
	  // - Proceeding so, based on stations, allows to extend the efficiency
	  //  evaluation downstream of last associated hit, and yet keep being
	  //  fair since any given particle is expected to stay alive throughout
	  //  the limited distance spanned by a detector station.
	  // - Edge effects (detector efficiency lower close to the central,
	  //  inactivated, hole) is taken into account, on a per station basis.
	  const TStation &s = stations[ist];
	  int spl, giveup, type = s.Type, edge;
	  for (spl = (int)s.IPlanes.size()-1, giveup = 0, edge = -1; spl>=0;
	       spl--) {
	    int ipl = s.IPlanes[spl], kpl = ipl-ipl0, jpl = jpls[kpl];
	    if (jpl<0) continue; float xx = xpl[jpl]; if (xx<750) {
	      giveup = 1; break;  // Limiting ourselves to upstream of RICH
	    }
	    float dx = xx-xL, yy = yL+yp*dx, zz = zL+zp*dx;
	    while (rp!=t.lPlnRef.rend() && *rp>ipl) rp++;
	    int status; if (*rp==ipl) {
	      if (type==39 || type==40) {   // Exclude MAs
		nMAs++; continue;
	      }
	      status = 0x3; if (ipl<iplF) iplF = ipl;
	    } 
	    else {
	      if (type==39 || type==40) continue;
	      status = setup.vDetect(ipl).InActive(yy,zz) ? 0x1 : 0;
	    }
	    if (status) {
	      if (type==22 || 28<=type && type<=29) {
		nVSATs++; if (nVSATs>2) { giveup = 1; break; } continue;
	      }
	      else {
		if (edge<0) { // First plane of station: evaluate edge effect
		  const TDetect &d = setup.vDetect(ipl);
		  yy += yy>d.X(1) ? 2 : -2; zz += zz>d.X(2) ? 2 : -2;
		  edge = setup.vDetect(ipl).InActive(yy,zz) ? 0 : 1;
		  if (!edge) {
		     yy += yy>d.X(1) ? -4 : 4; zz += zz>d.X(2) ? -4 : 4;
		     if (!setup.vDetect(ipl).InActive(yy,zz)) edge = 1;
		  }
		}
		nExpect += edge ? .5 : 1; if (status&0x2) nFound++;
	      }
	    }
	  }  // End loop on planes in station
	  if (giveup) break;
	}  // End loop on stations
	if (nVSATs<=2 && iplF<maxpl) {
	  float eff = nFound/nExpect;
	  if (iter<0 && eff<.75 || eff<.66) {
	    //         ***** TRACK is ``INEFFICIENT''... *****
	    if ((int)t.NDFs-nFound-nMAs>=nHitsMn) {
	      // ...UPSTREAM of RICH OK => STRIP AWAY DOWNSTREAM SEGMENT
	      t.Shorten(iplF);
	      if (QNstatus) {
		if (!t.QNewtonFit(1,QNmode)) QNstatus = 0;
	      }
	      if (!QNstatus) { t.QuickKF(-1,0); t.QuickKF(1,0); }
	    }
	    else ok = 0; // ...ELSE => REMOVE
	  }
	}
      }
      if (ok) {

	//   *************** #PROJs AND #SPACE POINTS ***************

	unsigned int projs, allProjs; int nSpacePts, nProjs, nAllProjs;
	int checkLastYZ = iter<0 &&// Early iteration: it's specific in that few
	  // combinations pass the PR proj. step, and hence those which do pass
	  // the space step haven't had to vie w/ many competitors. They may
	  // then pick up hits belonging in fact to those absent competitors. So
	  // more so as this precisely helps them pass the PR selections.
	  // Particularly difficult is zone 0x2, w/ its 2 sub-zones separated by
	  // the thickness of the RICH. Let's be more demanding there and let's
	  // discard tracks that have no obvious defects, contrary to those
	  // dealt w/ supra, but are unreliable (keeping in mind that those
	  // tracks will get another opportunity in later iterations).
	  // I) Track is deprived of either Y or Z info over a large length (
	  //  here we check the case of tracks extending downstream of the RICH,
	  //  requiring they be constrained in both Y and Z in their downstream
	  //  extension, restricting ourselves to [V]SAT).
	  // II) Track ending in the RICHwall: extensions to RICHwall being
	  //  unreliable when not backed by other hits further downstream.
	  igroup==1 && t.Hlast(0)>TOpt::dCut[68] /* i.e. RICHwall */ ? 1 : 0;
	double lastY = X0, lastZ = X0;
	list<int>::iterator ihp; const TStation *sPrv;
	for (ihp = t.lHitPat.begin(), nSpacePts=nProjs=nAllProjs = 0,
	       projs=allProjs = 0, sPrv = 0; ihp!=t.lHitPat.end(); ihp++) {
	  const THit &h = vecHit[*ihp];
	  const TPlane &p = setup.vPlane(h.IPlane);
	  const TStation *s = p.Station; if (s!=sPrv) {
	    sPrv = s; if (nProjs>=2) { nSpacePts++; nProjs = 0; projs = 0; }
	  }
	  unsigned int proj = 1<<p.IProj;
	  if (p.IFlag&0x30) {   // P-pixel planes: Add two proj....
	    // ...arbitrarily chosen as Y and Z proj.
	    if (!(0x1&projs))     {    projs |= 0x1;     nProjs++; }
	    if (!(0x1&allProjs))  { allProjs |= 0x1;  nAllProjs++; }
	    if (!(0x2&projs))     {    projs |= 0x2;     nProjs++; }
	    if (!(0x2&allProjs))  { allProjs |= 0x2;  nAllProjs++; }
	  }
	  else {
	    if (!(proj&projs))    {    projs |= proj;    nProjs++; }
	    if (!(proj&allProjs)) { allProjs |= proj; nAllProjs++; }
	  }
	  if (checkLastYZ) {
	    const TDetect &d = setup.vDetect(h.IPlane); double x = d.X(0);
	    if (d.Ca>.8) lastY = x; // Considering U/V w/ angle < pi/4 as Y
	    else if (p.IProj==zProj) lastZ = x;
	    if (x>750) {
	      int type = d.IType;
	      if (26<=type && type<=29 /* GM|GP */ || type==22 /* FI */)
		checkLastYZ = 2;
	    }
	  }
	}
	if (nProjs>=2) nSpacePts++;

	if (checkLastYZ) {
	  if (checkLastYZ==2 && (lastY<750 || lastZ<750)) ok = 0;
	}
	if (iter<0 && igroup==1 && // Track ending in the RICHwall
	    TOpt::dCut[67]<t.Hlast(0) && t.Hlast(0)<TOpt::dCut[68]) ok = 0;

	if (ok && (nAllProjs<3 || nSpacePts<3 &&
	     //  ***** DISCARD BAD KF CHI2 w/ FEW PROJ. or SPACE POINTS *****
	     //  The #space points is limited to <= 2 in zones 0x9 in LAT. On
	     // the other hand, bad chi2 w/ only 2 space points most certainly
	     // signals a bad track. Unless it's offset in time, and hence has
	     // bad T0. In LAT, this would correspond to accidentally coincident
	     // far halo tracks, which it might be interesting to reco and ID as
	     // a muon in order to check the muon triggers. Let's try then to
	     // make it easier for zone 0x8 and require few proj. in addition to
	     // few space points before rejecting the track in that case.
	     (igroup!=3 || nAllProjs<3)) &&
	    (t.IFit&0x8)==0 &&
	    (t.Chi2tot/(t.NDFs-3)>5 && iter<2 || // NDFs-3: to make for NDFs=4!
	     //  Upon earlier iterations, let's be stricter: ghost tracks, built
	     // w/ few proj. or space points, can steal hits from the low
	     // momentum tracks that can only be reco'd in later iterations
	     // where PR routes are wider opened. E.g. D*.2006, evt #883. (Note
	     // that one could consider rejecting the few proj./space points
	     // tracks altogether. Only that I don't want to introduce too big a
	     // change at once.)
	     t.Chi2tot/(t.NDFs-3)>3.5 && nAllProjs<3 && iter<1))
	  ok = 0;

	else if (ok) {	  //       ********** VERY GOOD TRACK **********
	  // #hits required is 7. Was 6. This improves overall efficiency.
	  // But, in some cases, produces many more 5-hit scifi tracks,
	  // with shared hits, which can have bridging later failing to
	  // select the good 6-hit segment among the so many alternatives.
	  // A first action to counter this side effect is introduced
	  // in "TEv::BridgeSegments", viz. give a malus to ``uncomplete''
	  // scifi track segments. It is envisioned to also modify the
	  // cleaning algorithm to have a second bridging heat, where tracks
	  // that have been cleaned of their counterparts, can be bridged
	  // a second time...
	  bool veryGoodT = 
	    (QNstatus &&
	     (t.NDics>6 && t.Chi2tot/(t.NDics-5)<1.5 ||
	      t.NDics>9 && t.Chi2tot/(t.NDics-5)<2) ||
	     !QNstatus &&
	     (t.NDFs>6 && t.Chi2tot/(t.NDFs-5)<1.5 ||
	      t.NDFs>9 && t.Chi2tot/(t.NDFs-5)<2)) &&
	    nAllProjs>=3 && nSpacePts>=3;
	  THit *worstHit = 0; if (veryGoodT && t.Scifi==1<<igroup) {
	    // If a scifi-(almost)only track:
	    // - Check time dispersion (<2 sigmas)
	    // - Single out worst timed hit if real badly timed
	    // (Since the processing relies on the detector time resolution, it
	    // cannot be applied to other VSAT-(almost)only tracks.)
	    t.UseHitTime();
	    veryGoodT &= 0<t.DispTime && // DispTime =-1 means it's ill defined
	      t.DispTime < .7; // I.e. 2 sigmas
	    double worstT, tT = t.MeanTime; vector<THit*> hs;
	    for (ihp = t.lHitPat.begin(), worstT = 0; ihp!=t.lHitPat.end();
		 ihp++) {
	      THit &h = vecHit[*ihp]; if (h.SigT<=0) continue; hs.push_back(&h);
	      const TPlane &p = setup.vPlane(h.IPlane);
	      const TDetect &d = setup.vDetect(p.IDetRef);
	      double dt = fabs(h.Time-tT); if (d.IType==22 && dt>worstT) {
		// Restrict search of worst timed hit to scifis hits
		worstT = dt; worstHit = &h;
	      }
	    }
	    double s1, st; int ih; for (ih = 0, s1=st = 0; ih<(int)hs.size();
					ih++) {
	      THit *h = hs[ih]; if (h==worstHit) continue;
	      double w = 1/h->SigT/h->SigT; s1 += w; st += h->Time*w;
	    }
	    if (fabs(worstHit->Time-st/s1)<.7)// Worst not so (<2 sigmas) bad?..
	      worstHit = 0;                   // ...=> Disregard worst hit
	  }

	  if (veryGoodT) {
	    t.IFit |= 0x10;  // Flag as very good track
#warning TODO: Have a TTrack::IFlag for both scifi-only and very-good tracks
	    for (ihp = t.lHitPat.begin(); ihp!=t.lHitPat.end(); ihp++) {
	      //  ********** FLAG HITS of VERY GOOD TRACKS **********
	      THit &h = vecHit[*ihp];
	      const TPlane &p = setup.vPlane(h.IPlane);
	      const TDetect &d = setup.vDetect(p.IDetRef);
#define PR_STRICT_ReUSE_NON_ENHANCED
#ifdef PR_STRICT_ReUSE_NON_ENHANCED
	      if (TOpt::ReMode[48]) {// ***** CASE of SCIFIS/PiXMPS in 2012 DVCS
		if ((d.IType==32 || d.IType==32) &&
		    vSATYetToBeEnhanced && // VSAT enhancement delayed...
		    nScifis+nPixMPs<2) continue; // ...=> Exclude scifis/pixMPs
	      }
#endif
	      if (d.IType==22) {                       // ***** CASE of SCIFIS
		if (&h==worstHit) continue;  // Exclude badly timed hit
#ifdef PR_STRICT_ReUSE_NON_ENHANCED
		if (vSATYetToBeEnhanced && // VSAT enhancement delayed...
		    // ...and few scifis involved: e.g. #22021728 (skip 121)
		    // /hpss/../03/P1E/DST_mergedDump/evtdump1-30194.raw
		    // (Note cGEMs are not considered here (impact not eval'd).)
		    nScifis<2) continue; // ...=> Exclude scifis
#endif
	      }
	      if (d.IType==41 && igroup==2 && iter<2)     // ***** CASE OF HO03 
		// Allow for several tracks to share its hits, given its crude
		// granularity and in order to let the scat'd mu get it.
		// (Beyond iteration iter>1, all high momenta must have already
		// been reco'd: therefore then no longer let HO03 hits through.)
		continue;
	      h.status = 1;                    // ***** FLAG HIT as VERY GOOD...
	      int mirror; if ((mirror = h.Mirror)!=-1) {    // ...and ITS MIRROR
		THit &hm = vecHit[mirror]; if (hm.Status!=-4) hm.status = 1;
	      }
	      // ...=> These flags are to be taken into account by this iter
	    }
	  }  // End block very good track
	}  // End blocks evaluating tracks
      }
      if (ok) {
	listTrack.push_back(t); // ***** STORE THE TRACK
	if (igroup==1 && TOpt::iCut[31] &&
	    nScifis>TOpt::iCut[31] && t.NDFs>7 &&
	    (t.Chi2tot/(t.NDFs-4)<2.0 ||  // Good chi2...
	     // ...w/ some tolerance for long track (1040, in 2010 setup, means
	     // longer than GM06-GM03 or FI06-GM02 but including GM06-GM02 or
	     // FI06-GM01) to make for the effect that MS makes on the total
	     // chi2 by deviating from its straight trajectory a track w/ a
	     // significant #hits at both ends, even when it has, as expected
	     // here, a high momentum.
	     t.Chi2tot/(t.NDFs-4)<2.5 && t.Hlast(0)-t.Hfirst(0)>1040)) {
	  // In 2010 (iCut[31] expected to be set), one wants to withdraw as
	  // many scifi hits as possible before they are fully enhanced (or else
	  // the # of combinations in FindProjQ explodes): this is achieved by
	  // setting here their status =5. Do this only when track has many (out
	  // of 2*FI05+2*FI55+3*FI06=7) scifi hits and their association in
	  // track is thus reliable. Reliability is also ensured by a chi2 cut.
	  // (Note: setting status =5, prevents the hits from being excluded
	  // until the completion of the present iteration: they will pass the
	  // cut "h.Status>statusMx".)
	  list<int>::iterator ihp;
	  for (ihp = t.lHitPat.begin(); ihp!=t.lHitPat.end(); ihp++) {
	    THit &h = vecHit[*ihp]; const TDetect &d = setup.vDetect(h.IPlane);
	    if (d.IType==22 && h.Status==0) h.status = 5;
	  }
	}
	if (TOpt::Hist[6]) {
	  static TH2D *tYP_PP[5], *tZP_PP[5], *tCH_PP[5];
	  static bool book = true;
	  if (book) {             // Book histos
	    book = false;
	    CsHistograms::SetCurrentPath("/Traffic/PrePattern");
	    for (int igr = 0; igr<5; igr++) {
	      char name[] = "tYP_PPi";
	      double xmax = TOpt::dPRpar[1+10*igr]*1.2;
	      if (igr==0) xmax *= 2;
	      sprintf(name,"tYP_PP%d",igr);
	      tYP_PP[igr] = new TH2D(name,name,50,-xmax,xmax,5,0,4);
	      sprintf(name,"tZP_PP%d",igr);
	      tZP_PP[igr] = new TH2D(name,name,50,-xmax,xmax,5,0,4);
	      sprintf(name,"tCH_PP%d",igr);
	      tCH_PP[igr] = new TH2D(name,name,100,0,50.,5,0,4);
	    }
	    CsHistograms::SetCurrentPath("/");
	  }
	  if (igroup<5) {
	    t.FindKine();
	    if (!(t.IFit&0x1)) {
	      t.QuickKF(-1,0);
	      t.QuickKF(1,0);
	      t.IFit |= 0x1;           // Set fit type bit pattern to |Straight
	    }
	    if (t.IKine!=-1) {
	      (tYP_PP[igroup])->Fill(t.Hfirst(3),(double)iter);
	      (tZP_PP[igroup])->Fill(t.Hfirst(4),(double)iter);
	      (tCH_PP[igroup])->Fill(t.Chi2tot/(t.NDFs-2),(double)iter);
	    }
	  }
	}
      }
      else if (iter!=iter_mx-1) {
	// Bad TTrack (and yet iteration to come) => Release all hits
	for (int ih = ih0; ih<nh-1; ih++) {
	  int iref = hits[ih];
	  if (iref>=0) {
	    if (t.Scifi!=1<<igroup)  // No scifis-only, case hit used 2 times...
	      // ...of which one is wrong but may improve in subsequent iter...
	      ifl[iref] = 0;   // ...=> give it a 2nd opportunity
	    else ifl[iref]--;
	  }
	}
      }

    } // End of loop over found tracks

    if (iter!=iter_mx-1) {

      // ******************** TAG ALL USED HITS ********************
      // - In view of excluding them in next iteration.
      // - Except scifi hits if "PR_ReUSE_NON_ENHANCED" & "vSATYetToBeEnhanced".
      // - Exclude also ``alternatives''. There is no a priori reason to
      //  proceed so. But it turned out that it gives better perfs than
      //  excluding only hits proper (one achieves the latter by looping
      //  over all found tracks, and tagging their hits).
      for (int jhit=nfree = 0; jhit<nHits; jhit++) {
	THit &h = vecHit[href[jhit]];
	if (h.Status==5) { // Special case of scifi hits not to be ReUSE'd.
	  h.status = 1; continue;
	}
	const TPlane &p = setup.vPlane(h.IPlane);
	const TDetect &d = setup.vDetect(p.IDetRef);
	if (ifl[jhit]<=0) {   // In fact should be >=0 in any case
	  nfree++; if (d.IType==29) nfree++; continue;
	}
#define PR_ReUSE_NON_ENHANCED
#ifdef PR_ReUSE_NON_ENHANCED
	if (vSATYetToBeEnhanced) {    // Scifis enhanced due in later iter...
	  // (Note cGEMs are not considered here (impact not evaluated).)
	  if (d.IType==22 ||         // ...=> allow scifis to be reused...
	      d.IType==32 && TOpt::ReMode[48]) { // (also pixMPs in 2012 DVCS)
	    // ...Yet, in 2010 (iCut[31] set), one wants to withdraw as many
	    // scifi hits as possible before they are enhanced (otherwise,
	    // the # of combinations in FindProjQ explodes): this is achieved by
	    // setting their status =5, cf. supra
	    nfree++; continue;
	  }
	}
#endif
	if (d.IType==41 && igroup==2 && iter<2) // Case of HO03:...
	  // ...allow for several tracks to share its hits, given its cude
	  // granularity and in order to let the scattered mu pass coral's ID
	  // (Beyond iteration iter>1, all high momenta must have already been
	  // reco'd: therefore then no longer let HO03 hits through.)
	  continue;
	h.status = 1;                              // ***** FLAG USED HIT
	int mirror; if ((mirror = h.Mirror)!=-1) { // ***** AND ITS MIRROR
	  THit &hm = vecHit[mirror]; if (hm.Status!=-4) hm.status = 1;
	}
      }
      if (nfree<TOpt::iPRpar[igroup*10+5] && !vSATYetToBeEnhanced)
	iter = iter_mx;        // Full success => Don't waste time in a 2nd iter
      //else..                    Too many hits left unused => subsquent iter

    }
    } // end of iteration

    if (igroup==1) { // Let's be cautious for the time being...
      //      ********** DISAMBIGUATION in VERTICAL **********
      // Tracks close in vertical (i.e. somewhat less than ST pitch). In the
      // case that stereo angles are quasi horizontal (case of DC,ST,PS), could
      // be that they swapped their vertical projections.
      list<TTrack>::iterator it, jt; unsigned int lastId = listTrack.back().Id;
      it = listTrack.begin(); while (it!=listTrack.end()) {
	TTrack &ti = *it; if (ti.Type!=1<<igroup || (ti.Scifi)) {
	  it++; continue;// Skip previous groups or VSAT-(almost)only tracks
	}
	if (ti.Id>lastId) // From now on, tracks created by dismabiguation...
	  break;          // ...=> exit
	float Li = X0-ti.Hfirst(0), ypi = ti.Hfirst(3), zpi = ti.Hfirst(4);
	float y0i = ti.Hfirst(1)+Li*ypi, z0i = ti.Hfirst(2)+Li*zpi;
	double chi2 = ti.IFit&0x8 ? ti.Chi2aux/(ti.NDics-4) :
	  ti.Chi2tot/(ti.NDFs-4);
	bool doSwap = false; jt = it; jt++; while (jt!=listTrack.end()) {
	  TTrack &tj = *jt;
	  double chj2 = (tj.IFit&0x8) ? tj.Chi2aux/(tj.NDics-4) :
	    tj.Chi2tot/(tj.NDFs-4); if (chi2<3 && chj2<3) {
	      jt++; continue; // Both "ti" and "tj" chi2-ok => next "ti"
	    }
	  float Lj = X0-tj.Hfirst(0), ypj = tj.Hfirst(3), zpj = tj.Hfirst(4);
	  float y0j = tj.Hfirst(1)+Lj*ypj, z0j = tj.Hfirst(2)+Lj*zpj;
	  float dy0 = y0j-y0i, dyp = ypj-ypi, dz0 = z0j-z0i, dzp = zpj-zpi;
	  float hDist = dy0*dy0+1e4*dyp*dyp; bool hClose = hDist<1;
	  float vDist = dz0*dz0+1e4*dzp*dzp; if (vDist<.3) {
	    int quasiYZ, zProjs[2]; bool hOpt = hClose; list<int>::iterator ipr;
	    int ij; TTrack *tij; for (ij = 0, tij = &ti, quasiYZ = 0;
				      ij<2; ij++) {
	      for (ipr = tij->lPlnRef.begin(), zProjs[ij] = 0;
		   ipr!=tij->lPlnRef.end(); ipr++) {
		const TPlane &p = setup.vPlane(*ipr);
		const TDetect &d = setup.vDetect(p.IDetRef);
		if (11<=d.IType && d.IType<=18) { // Drifts: 1/2 per layer
		  if      (p.IProj==zProj) zProjs[ij] += 1;
		  else if (p.IProj!=yProj) quasiYZ |= 1<<ij;
		}
		else {
		  if      (p.IProj==zProj) zProjs[ij] += 2;
		  else if (p.IProj!=yProj) {
		    if (d.IType>2) {
		      // !Y, !Z in non-drift, non-MWPC means 45 deg. stereo...
		      if (hOpt)        // ...close in H, allow for one exception
			hOpt = false;
		      else {           // ...else give-up
			quasiYZ &= 1<<(1-ij); break;
		      }
		    }
		    else quasiYZ |= 1<<ij;
		  }
		}
	      }
	      if (!quasiYZ) break; tij = &tj;
	    }
	    if (!(quasiYZ&0x1)) break; // "ti" unvalid => next one
	    if (quasiYZ && // Quasi YZ => No efficient stereo angles
		// If single Z proj. to swap, the case must have solved supra.
		zProjs[0]>=3 && zProjs[1]>=3) {
	      TTrack ni, nj, *ns[2] = {&ni,&nj}, *nij; list<int>::iterator ihp;
	      for (ij = 0, tij = &ti; ij<2; ij++) {
		for (ipr = tij->lPlnRef.begin(), ihp = tij->lHitPat.begin();
		     ipr!=tij->lPlnRef.end(); ipr++, ihp++) {
		  const TPlane &p = setup.vPlane(*ipr); THit &h = vecHit[*ihp];
		  if (p.IProj==zProj) ns[1-ij]->AddHit(h);
		  else                ns[ij]->AddHit(h);
		}
		tij = &tj;
	      }
	      int ok; double oSChi2, nSChi2;
	      for (ij = 0, tij = &ti, nij = &ni, ok = 1, oSChi2=nSChi2 = 0;
		   ij<2; ij++) {
		double oldChi2, newChi2; if (tij->IFit&0x8) {
		  oldChi2 = tij->Chi2tot/(tij->NDics-4);
		  nij->Haux = tij->Haux; nij->QNewtonFit(1,5);
		  newChi2 = nij->Chi2tot/(nij->NDics-4);
		}
		else {
		  oldChi2 = tij->Chi2tot/(tij->NDFs-4);
		  nij->QuickKF(-1,0);
		  newChi2 = nij->Chi2tot/(nij->NDFs-4);
		}
		if (newChi2>oldChi2 && !hClose || newChi2>1.3*oldChi2) {
		  ok = 0; break;
		}
		oSChi2 += oldChi2; nSChi2 += newChi2; tij = &tj; nij = &nj;
	      }
	      if (ok) {
		if (nSChi2<oSChi2) {
		  ni.Type = 1<<igroup; nj.Type = 1<<igroup;
		  listTrack.push_back(ni); listTrack.push_back(nj);
		  listTrack.erase(jt); doSwap = true; break;
		}
	      }
	    }
	  } // End examining (ti,tj)
	  jt++;
	} // End loop on tj
	if (doSwap) listTrack.erase(it++);
	else it++;
      } // End loop on ti
    } // End vertical disambiguation

    if (igroup==2) { // Let's be cautious for the time being...
      //          ********** MERGING TRACKS **********
      list<TTrack>::iterator it, jt;
      it = listTrack.begin(); while (it!=listTrack.end()) {
	TTrack &ti = *it;
	if (ti.Type!=1<<igroup ||     // Skip previous reco zones or...
	    ti.Chi2tot/(ti.NDFs-4)>3) {// ...bad KF (No QN here since zone==0x4)
	  it++; continue;
	}
	float Li = X0-ti.Hfirst(0), ypi = ti.Hfirst(3), zpi = ti.Hfirst(4);
	float y0i = ti.Hfirst(1)+Li*ypi, z0i = ti.Hfirst(2)+Li*zpi;
	bool merged = false; jt = it; jt++; while (jt!=listTrack.end()) {
	  TTrack &tj = *jt;
	  if (tj.Chi2tot/(tj.NDFs-4)>3) { jt++; continue; }
	  float Lj = X0-tj.Hfirst(0), ypj = tj.Hfirst(3), zpj = tj.Hfirst(4);
	  float y0j = tj.Hfirst(1)+Lj*ypj, z0j = tj.Hfirst(2)+Lj*zpj;
	  float dy0 = y0j-y0i, dyp = ypj-ypi, dz0 = z0j-z0i, dzp = zpj-zpi;
	  float absDist = dy0*dy0+1e4*dyp*dyp+dz0*dz0+1e4*dzp*dzp;
	  if (absDist>.7) { jt++; continue; } // ***** CRUDE CUT(s) ON DISTANCE
	  TTrack *tSave = 0; int ijSave = -1, haveZ = 1; if (absDist>.04) {
	    if (dy0*dy0+1e4*dyp*dyp<.04) {
	      // Try and rescue tracks (e.g. in D*.2006.03, #81, #140):
	      //  - PA-only ones, i.e. short, lacking Z info
	      //  - PA+ST+DW ones (i.e. ending w/ a sequence of drift-like dets)
	      //   w/ only one Z proj. drift plane.	      
	      static double xHO03 = 0; if (xHO03==0) {
		for (int ipl = setup.vIplFirst()[2];
		     ipl<(int)setup.vIplLast()[2]; ipl++) {
		  const TDetect &d = setup.vDetect(ipl);
		  if (d.IType==41 && d.Name.find("HO03")==0) {
		    xHO03 = d.X(0) + /* some margin */10; break;
		  }
		}
		if (xHO03==0) xHO03 = -100000;
	      }
	      int ij, driftT[2], zProjs[2], iplDriftZ[2]; TTrack *tij;
	      for (ij = 0, tij = &ti; ij<2; ij++) {
		list<int>::iterator ipr;
		for (ipr = tij->lPlnRef.begin(), driftT[ij]=zProjs[ij] = 0,
		       iplDriftZ[ij] = -1; ipr!=tij->lPlnRef.end(); ipr++) {
		  const TPlane &p = setup.vPlane(*ipr);
		  const TDetect &d = setup.vDetect(p.IDetRef);
		  if (d.IType>=41) continue; // Exclude hodos (HO03, that is)
		  if (11<=d.IType && d.IType<=18) driftT[ij] = 1;
		  else if (driftT[ij]==1)         driftT[ij] = -1;
		  if (p.IProj==zProj) {
		    zProjs[ij]++;
		    if (11<=d.IType && d.IType<=18) iplDriftZ[ij] = *ipr;
		  }
		}
		tij = &tj;
	      }
	      bool rescuePAT = ti.Hlast(0)<xHO03 && !zProjs[0] ||
		/* */          tj.Hlast(0)<xHO03 && !zProjs[1];
	      int rescueDriftT = -1;
	      if (driftT[0]==1 && zProjs[0]<=1 && zProjs[1]>1) rescueDriftT = 0;
	      if (driftT[1]==1 && zProjs[1]<=1 && zProjs[0]>1) rescueDriftT = 1;
	      if (rescueDriftT>=0 && iplDriftZ[rescueDriftT]>=0) {
		TTrack *tk = rescueDriftT?&tj:&ti;
		if (tk->NHits<8) rescueDriftT = -1;
		else {
		  tSave = new TTrack(*tk); ijSave = rescueDriftT;
		  tk->SubHit(iplDriftZ[rescueDriftT]);
		  if (!tk->QuickKF(1,0) || !tk->QuickKF(-1,0)) {
		    *tk = *tSave; delete tSave; tSave = 0; ijSave = -1;
		    rescueDriftT = -1;
		  }
		  else tk->IFit = 0x1; 
		}
	      }
	      if (!rescuePAT && rescueDriftT<0) { jt++; continue; }
	      haveZ = zProjs[0]*zProjs[1];
	    }
	    else { jt++; continue; }
	  }
	  if (ti.IFit&0x8) { ti.QuickKF(1,0); ti.QuickKF(-1,0); }
	  if (tj.IFit&0x8) { tj.QuickKF(1,0); tj.QuickKF(-1,0); }
	  // Covariant matrix: a fast extrapolation will be used. 2x2 blocks
	  // are tranformed by transfer matrix T via  T*C*tT where T is
	  // T(0,0)=T(1,1)=1, T(0,1)=-L, T(1,0) = 0; L being extrap. length.
	  TMatrixF covi[3] = { TMatrixF(2,2), TMatrixF(2,2), TMatrixF(2,2) };
	  TMatrixF covj[3] = { TMatrixF(2,2), TMatrixF(2,2), TMatrixF(2,2) };
	  int m, k, l; for (m = 0, k = 1; m<2; m++, k++) {
	    covi[m](0,0) = ti.Hfirst(k,k); l = k+2; covi[m](1,1) = ti.Hfirst(l,l);
	    covi[m](0,1)=covi[m](1,0) = ti.Hfirst(k,l);
	    covj[m](0,0) = tj.Hfirst(k,k); l = k+2; covj[m](1,1) = tj.Hfirst(l,l);
	    covj[m](0,1)=covj[m](1,0) = tj.Hfirst(k,l);
	  }
	  covi[2](0,0) = ti.Hfirst(1,2); covi[2](1,1) = ti.Hfirst(3,4);
	  covi[2](0,1) = ti.Hfirst(1,4); covi[2](1,0) = ti.Hfirst(3,2);
	  covj[2](0,0) = tj.Hfirst(1,2); covj[2](1,1) = tj.Hfirst(3,4);
	  covj[2](0,1) = tj.Hfirst(1,4); covj[2](1,0) = tj.Hfirst(3,2);
	  for (m = 0; m<3; m++) {
	    TMatrixF ci(covi[m]);
	    covi[m](0,0) = ci(0,0)+Li*(ci(0,1)+ci(1,0))+Li*Li*ci(1,1);
	    covi[m](0,1) = ci(0,1)+Li*ci(1,1); covi[m](1,0) = ci(1,0)+Li*ci(1,1);
	    TMatrixF cj(covj[m]);
	    covj[m](0,0) = cj(0,0)+Lj*(cj(0,1)+cj(1,0))+Lj*Lj*cj(1,1);
	    covj[m](0,1) = cj(0,1)+Lj*cj(1,1); covj[m](1,0) = cj(1,0)+Lj*cj(1,1);
	  }

	  TMatrixF covij(4,4); for (m = 0; m<2; m++) {
	    for (k = 0; k<2; k++) for (l = 0; l<2; l++) {
	      int i = 2*m+k, j = 2*m+l; covij(i,j) = covi[m](k,l)+covj[m](k,l);
	    }
	  }
	  covij(0,2)=covij(2,0) = covi[2](0,0)+covj[2](0,0);
	  covij(1,3)=covij(3,1) = covi[2](1,1)+covj[2](1,1);
	  covij(3,0)=covij(0,3) = covi[2](0,1)+covj[2](0,1);
	  covij(2,1)=covij(1,2) = covi[2](1,0)+covj[2](1,0);
	  TMatrixF hess(TMatrixF::kInverted,covij);
	  TVectorF d(0,3,dy0,dyp,dz0,dzp,"END");
	  TVectorF tmp = d; tmp *= hess; double dist = tmp*d;
	  //                                     ***** CUT ON TRACKS DISTANCE...
	  if (dist>16 && haveZ ||                    // ...4 sigmas except if...
	      dist>64) {// ...no Z info in any of of ti/j, then: 8 sigmas
	    jt++; continue; if (tSave) delete tSave;
	  }
	  TTrack t;       //      ***** MERGED TRACK *****
	  if (!t.Merge(ti,tj)) { jt++; continue; if (tSave) delete tSave; }
	  if (!t.QuickKF(1,0) || !t.QuickKF(-1,0)) {
	    jt++; continue; if (tSave) delete tSave;
	  }
	  double tli = ti.Hlast(0)-ti.Hfirst(0), tlj = tj.Hlast(0)-tj.Hfirst(0);
	  double tl = t.Hlast(0)-t.Hfirst(0);
	  double chi2ij = (ti.Chi2tot*(1+.25*tl/tli)+tj.Chi2tot*(1+.25*tl/tlj))/
	    (ti.NDFs+tj.NDFs-4);
	  double chi2 = t.Chi2tot/(t.NDFs-4);
	  if (chi2<3.5 || chi2<5 && chi2<chi2ij) {
	    if (dist>16) CsErrLog::msg(elError,__FILE__,__LINE__,
 "Merging tracks #%d,%d (%.1f,%.1fNDF) despite distance = %.1f => %.1fNDF",
			  ti.Id,tj.Id,ti.Chi2tot/(ti.NDFs-4),
			  tj.Chi2tot/(tj.NDFs-4),dist,chi2);
	    t.IFit = 0x1; t.Type = 1<<igroup;// Other attributes (Scifi) skipped
	    listTrack.push_back(t); listTrack.erase(jt); merged = true;
	    break;
	  }
	  else if (tSave) {
	    if (ijSave==0) ti = *tSave;
	    else           tj = *tSave;
	  }
	  jt++;
	}  // End loop on "tj"
	if (merged) { listTrack.erase(it++); continue; }
	it++;
      }  // End loop on "ti"
    }  // End of merging block

    Traffic::Ref().Stopwatch.Stop(17+(igroup<3?igroup:3));
  } // end of loop over det. groups in the zone

#ifdef DEBUG
  static int jdebug = 0;
  if (jdebug) DumpTrackList();
#endif
}

// ***** FIND REFERENCE TO MIRROR GIVEN REFERENCE TO THE ORIGINAL HIT *****
int findMirrorRef(int iref, int iHit, int iMirror, int *href)
{
  // Most of the time mirror is at a +/-1 distance.
  // But this is not guaranteed => Check it
  int dref = iMirror-iHit;
  if (abs(dref)==1) iref += dref;
  else {   // If not, search for it.
    int jref, iref0, iref1;
    if (dref>0) { iref0 = iref+1; iref1 = iref+dref; }
    else        { iref0 = iref+dref; iref1 = iref-1; }
    for (jref = iref0, iref = -1; jref<=iref1; jref++) {
      if (href[jref]==iMirror) { iref = jref; break; }
    }
  }
  return iref;
}
