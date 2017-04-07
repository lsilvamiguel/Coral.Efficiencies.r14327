// $Id: TEvCleanTrackList.cc 14094 2015-11-06 15:28:48Z lsilva $

/*!
   \file    TEvCleanTrackList
   \brief   Clean track list == Retain only one out of several sharing a common hit. The case of SIs is treated independently, cf. "TEv::CleanTrackListSI".
   \author  Y.B
   \version $Revision: 14094 $
   \date    $Date: 2015-11-06 16:28:48 +0100 (Fri, 06 Nov 2015) $ 

*/

#include <cstdio>
#include <iostream>
#include "Coral.h"
#include "CsErrLog.h"
#include "TSetup.h"
#include "TEv.h"
#include "TOpt.h"
#include "TLattice.h"

using namespace std;

void TEv::CleanTrackList(int mode)
{
  int   maxhit = 6000; 
  int   nh = vecHit.size(), jh, jcom, igr, ok, nFailures, nRescues;
  // If too many hits, i.e. more than allocated for, i.e. > "maxhit"
  // => Give up
  if (nh>maxhit) {
    CsErrLog::msg(elError,__FILE__,__LINE__,"%d Hits > maxhit = %d",nh,maxhit);
    return;
  }

  const TSetup& setup = TSetup::Ref();
  int   hit_use[maxhit];
  float hit_worst[maxhit];

  list<TTrack>::iterator it, track_indices[maxhit];

#ifdef DEBUG
  int nDelete = 0;
  //#  define CleanTL_DEBUG
#  ifdef CleanTL_DEBUG
  int hitDebugInfo[maxhit];
  for (jh = 0; jh<nh; jh++) hitDebugInfo[jh] = 0;
#  endif  
#endif
  do { // *************** LOOP UNTIL HIT IN COMMON ARE EXHAUSTED ***************

    // **********  CLEAR USE-COUNTERS

    for (jh = 0; jh<nh; jh++) {
      hit_use[jh] = 0; hit_worst[jh] = -1e6;
    }

    // ********** FOR EACH HIT: AT WHAT TRACK IT IS WORST **********
  
    for (it = listTrack.begin(); it!=listTrack.end(); it++) { // loop on tracks
      TTrack &t = *it;
      float qual;
      if (t.Type==0x10 &&           // Beam tracks: only upon option
	  (TOpt::ReMode[15]&0x10)==0) continue; 
      if (mode==0) {             // ***** "mode==0"... *****
	if (t.IFit&0x8 &&        // ***** QN FITTED TRACKS ONLY, EXCEPT... *****
	    // ...THOSE EXTENDED INTO ZONE 0x4 (for TTrack::Chi2tot refers
	    // only to the fitting the 0x3 part of the track at this point)
	    (t.Type&0x7)!=0x7) {
	  if (t.NDics<=6) {
	    if (t.NDFs>10) continue; // If dico covers very little of track
	    else if (t.NDics<=4) qual = 10;
	    else if (t.NDics==5) qual = 5+t.Chi2tot;
	    else                 qual = t.Chi2tot;
	  }
	  else qual = t.Chi2aux/(t.NDics-5+t.NDics*t.NDics*.03);
	  if (t.IFit&0x10) qual -= 1; // Give a bonus to tracks w/ 3 space points
	}
	else continue;
      }
      else {                     // ***** "mode==1": ALL TRACKS *****
	if (t.Type==0x10 && TOpt::iCut[38]) {
	  qual = (t.Chi2tot+TOpt::iCut[38]/100*t.DispTime*t.DispTime/t.SigmaTime/t.SigmaTime)/
	    (t.NDFs-4);
	}
	else if (t.NDFs<=6) {
	  if      (t.NDFs<=4) qual = 10;
	  else if (t.NDFs==5) qual = 5+t.Chi2tot;
	  else                qual = t.Chi2tot;
	}
	else if (t.Type==0x10 &&   // Special processing of scifi/Si segments...
		 TOpt::iCut[23])   // ...to be enabled for 2003 data, given that
	  // - the low redundancy in 2003 4FIs+8Sis
	  // - coupled to the high multiplicity of off-time Si hits (until good
	  //  Si timing is available),
	  // can have off-time tracks steal hits, particulalrly scifi ones,
	  // from the genuine incident muon: we include the #hits in the
	  // quality function, for the detectors' efficiency is expected to be
	  // lower for off-time tracks. Cf., e.g., evt #45098976 in:
	  // "/hpss/../03/P1E/DST_mergedDump/evtdump1-30194.raw"
	  qual = 2*t.Chi2tot/(t.NDFs-5+t.NHits*t.NHits*.03)-t.NHits;
	else {
	  if (fabs(t.Hfirst(5))<1)
	    qual = t.Chi2tot/(t.NDFs-5+t.NHits*t.NHits*.05);
	  else
	    qual = t.Chi2tot/(t.NDFs-5+t.NHits*t.NHits*.03);
	}
      }
      //float qual = t.Chi2tot/(t.NHits-6+
      //		      (t.NHits*t.NHits)*(t.NHits*t.NHits)*.000004);
      if      (t.Type==0x1 || t.Type==0x2 || t.Type==0x4) qual += 20;
      else if ((t.Type&0x7)==0x7) qual -= 5;

#ifdef CleanTL_DEBUG
      bool kineSet = false;
#endif  
      for (list<int>::iterator ih = t.lHitPat.begin();
	   ih!=t.lHitPat.end(); ih++) {      // ***** LOOP OVER TRACK HITS *****
	if (*ih<0) continue;
	const THit &h = vecHit[*ih];
	if (h.IPlane>=TLattice::IplLast &&  // "IplLast" is outside scope
	    !mode) continue;
	const TPlane  &p = setup.vPlane(h.IPlane);
	const TDetect &d = setup.vDetect(p.IDetRef);
	if (d.IType==41) // Exclude hodos, except HIs: they have poor resolution
	  continue;      // => they are allowed to belong to several tracks.


	hit_use[*ih]++;

	bool isWorse = qual>hit_worst[*ih];
	//#define CleanTL_USE_TIME
#ifdef CleanTL_USE_TIME
	if (hit_use[*ih]>1) {
	  if (d.IType==22) {
	    float hT = h.Time, hsT2 = h.SigT*h.SigT; // Assuming time is always defined for scifis...
	    float tProbs[2]; int i; TTrack *ti; bool haveTime;
	    for (i = 0, ti = &t, haveTime = true; i<2; i++) {
	      if (ti->SigmaTime<0) ti->UseHitTime();
	      if (ti->SigmaTime<0) { haveTime = false; break; }
	      float w2 = ti->SigmaTime*ti->SigmaTime+hsT2;
	      float dt = ti->MeanTime-hT;
	      tProbs[i] = TMath::Prob(dt*dt/w2,1);
	      ti = &(*track_indices[*ih]);
	    }
	    if (haveTime) {
	      static float coeff = 1;
	      isWorse = qual-coeff*tProb[0]>hit_worst[*ih]-coeff*tProb[1];
	    }
	  }
	}
#endif
	//#define CleanTL_USE_MULTIPLE
#ifdef CleanTL_USE_MULTIPLE
	if (hit_use[*ih]>1 && fabs(qual-hit_worst[*ih])<2 && d.IType==22 &&
	    qual<5 /* Excludes in particular single segment tracks */) {
	  //                                            ***** IF AMBIGUITY *****
	  TTrack &tp =  *track_indices[*ih];
	  int nDF  = t .NDFs-5; float chi2  = nDF >0 ? t .Chi2tot/nDF  : 0;
	  int nDFp = tp.NDFs-5; float chi2p = nDFp>0 ? tp.Chi2tot/nDFp : 0;
	  if (nDF>0 && chi2<3 || nDFp>0 && chi2p<3) {
	    const TStation *&st = p.Station;
	    int npl = st->IPlanes.size();
	    int ipli = st->IPlanes[0], iplf = st->IPlanes[npl-1];
	    int nHsSt, nHsStp; list<int>::iterator jh = ih; jh++;
	    for (nHsSt = 1; jh!=t.lHitPat.end(); jh++) {
	      if (*jh<0) continue; THit &hj = vecHit[*jh];
	      if (hj.IPlane>iplf) break; nHsSt++;
	    }
	    jh = ih; jh--; list<int>::iterator jh0 = t.lHitPat.begin(); jh0--;
	    for (         ; jh!=jh0; jh--) {
	      if (*jh<0) continue; THit &hj = vecHit[*jh];
	      if (hj.IPlane<ipli) break; nHsSt++;
	    }
	    for (jh = tp.lHitPat.begin(), nHsStp = 0; jh!=tp.lHitPat.end();
		 jh++) {
	      if (*jh<0) continue; THit &hj = vecHit[*jh];
	      if (hj.IPlane<ipli) continue; if (hj.IPlane>iplf) break;
	      nHsStp++;
	    }
	    if (isWorse) {
	      if (nHsSt >=3 && nHsStp<3) isWorse = false;
	    }
	    else {
	      if (nHsStp>=3 && nHsSt <3) isWorse = true;
	    }
	  }
	}	
#endif
#ifdef CleanTL_DEBUG
	if (hit_use[*ih]>1 && qual<5) {
	  TTrack &tp =  *track_indices[*ih];
	  if (hitDebugInfo[*ih]!=int(t.Id+tp.Id)) {
	    hitDebugInfo[*ih] = t.Id+tp.Id;
	    tp.FindKine(); if (!kineSet) { t.FindKine(); kineSet = true; }
	    if ((isWorse && t.IKine>=0 && t.IKine==h.IKine &&
		 (tp.IKine<0 || tp.IKine!=h.IKine)) ||
		(!isWorse && tp.IKine>=0 && tp.IKine==h.IKine &&
		 (t.IKine<0 || t.IKine!=h.IKine))) {
	      float L = t.Hlast(0)-t.Hfirst(0), Lp = tp.Hlast(0)-tp.Hfirst(0);
	      printf("CleanTL_DEBUG: h %d MC %d\n",h.IHit,h.IKine);
	      printf("t  Id=%4d 0x%x sf 0x%x MC=%3d #hs %3d chi2/N(0x%x) %5.2f L %5.1f Q %5.2f q %5.2f\n",
    t.Id, t.Type, t.Scifi, t.IKine, t.NDFs, t.IFit, t.Chi2tot/(t.NDFs-5),  L,
    qual,          t.Chi2tot/(t.NDFs-5+t.NHits*t.NHits*.03));
	    printf("t' Id=%4d 0x%x sf 0x%x MC=%3d #hs %3d chi2/N(0x%x) %5.2f L %5.1f Q %5.2f q %5.2f\n",
    tp.Id,tp.Type,tp.Scifi,tp.IKine,tp.NDFs,tp.IFit,tp.Chi2tot/(tp.NDFs-5),Lp,
    hit_worst[*ih],tp.Chi2tot/(tp.NDFs-5+tp.NHits*tp.NHits*.03));
	    }
	  }
	}
#endif
	if (isWorse) {
	  hit_worst[*ih] = qual; track_indices[*ih] = it;
	}
      }
    }

    // ********** WORST HIT OF ALL == jcom <-> worst_qual

    float worst_qual;
    for (jh = 0, jcom = -1, worst_qual = -100.; jh<nh; jh++) {
      if ( hit_use[jh]>1 && hit_worst[jh]>worst_qual) {
	worst_qual = hit_worst[jh]; jcom = jh;
      }
    }
  
    if (jcom>=0) {

      // ********** RETRIEVE  THE WORST TRACK CONTAINING IT ...

      it =  track_indices[jcom]; TTrack &t = *it;
      bool goodTrack = t.NDFs>5 || t.Chi2tot/(t.NDFs-5)<5;
#ifdef DEBUG
      nDelete++;
#endif

      // ********** NUMBER OF UNSHARED HITS IN DIFFERENT ZONES => ok

      int nhits[5]; memset((void*)nhits,0,sizeof(nhits));
      int nStations[5]; memset((void*)nStations,0,sizeof(nStations));
      int nPixelGMs = 0, nStationsTot = 0;
      // Projections (excluding worst hit)
      unsigned int proj, projs; // Current, all
      static unsigned int sProjs;
      unsigned int cProjs; int nCProjs;// Contrained proj., i.e. w/ #stations>=2
      list<int>::iterator ih, ip; const TStation *sPrv;
      for (ih = t.lHitPat.begin(), ip = t.lPlnRef.begin(), sPrv = 0,
	     projs=cProjs = 0, nCProjs = 0; ih!=t.lHitPat.end(); ih++, ip++) {
	if (*ih<0) continue;
	THit &h = vecHit[*ih];
	if (hit_use[*ih]<=1 || goodTrack &&  // If bad track, it's enough to...
	    // ...find that it shares several other hits w/o checking that
	    // it is actually the worst track for these hits
	    track_indices[*ih]!=it) {
	  int ipl = *ip; 	    // Which zone is this iplane in?
	  for (igr = 0; igr<(int)setup.vIplFirst().size(); igr++) {
	    if (setup.vIplFirst()[igr]<=ipl &&
		ipl<=setup.vIplLast()[igr]) {
	      nhits[igr]++;
	      const TPlane &p = setup.vPlane(ipl);
	      const TStation *&s = p.Station; if (s!=sPrv) {
		sPrv = s; nStations[igr]++; nStationsTot++;
		if (s->Type==29) {
		  nPixelGMs++; nhits[igr]++;// GP P-pixels yield TW0 hits
		}
		sProjs = 0;
	      }
	      proj = 1<<p.IProj; if (!(proj&sProjs)) {
		sProjs |= proj;
		if (!(proj&projs)) projs |= proj;
		else {
		  if (!(proj&cProjs)) nCProjs++; cProjs |= proj;
		}
	      }
	      break;
	    }
	  } // End loop determining zone of current plane
	}
      }
      if (nCProjs<2) {
	// Single out and erase "unconstrained tracks", an example of which is
	// a XYU space point w/ extra hits along only one proj., which leaves
	// the track unconstrained along one of horiz. or vert. dimensions.
	// This is evaluated w/ "nCProj" = # of constrained projections.
	listTrack.erase(it); continue;
      }
      if (nStationsTot==2 && nPixelGMs==2) {
	// Track w/ 2 pixels from pixel GPs (which corresponds to a grand
	// total of 2 TStations). This kind of track escapes the condition infra
	// on #stations >=2 (because pixel GPs constitute one TStation per se).
	// And that's OK, because pixels give high excellent tracking constrain.
	// (E.g.: a track in zones 0x6, w/ GP03P1 and GP03P2 and no extension
	// downstream is considered reliable enough. A track w/ only GP02P1,P2
	// GP03P1,P2 is borderline, yet accepted.) But a single track w/ only 2
	// pixels is not so much useful: it's probably not a ghost, but the fact
	// that no other hits could be associated to it means that it is
	// probably also off-time. => Let's erase it.
	// (Note: Only at this point of the tracking algorithm is it possible
	// to encounter such a track (in principle). Therefore it's OK to have
	// this rejection taking place here.)
	listTrack.erase(it); continue;
      }
      for (igr=ok=nFailures=nRescues = 0; igr<(int)setup.vIplFirst().size();
	   igr++) {
	int jgr = 1<<igr, kgr = 0x100<<igr;
	if (t.Type&jgr) {
	  if (nStations[igr]<2) {// ***** REQUIRE AT LEAST TWO TStations
	    ok |= jgr; nFailures++;
	  }
	  // ***** DETERMINE WHETHER TTrack HAS ENOUGH HITS PER ZONE *****
	  if (((t.Scifi&jgr) ||            // Scifi-(almost)only
	       (t.Scifi&kgr) && igr!=0) && // VSAT-(almost)only except in zone 0x1 (where the redundancy is somewhat higher)
	      // This last condition infra is probably not useful, given that the VSAT only attibute can in principle not be granted in zone 0x10...
	      igr!=4) {
	    switch (igr) {
	    case 0:
	      if (nhits[igr]<5) {
		ok |= jgr; nFailures++; if (nhits[igr]==4) nRescues++;
	      }
	      break;
	    case 1:
	      if (TOpt::iCut[31]) { // Special VSAT handling in zone 0x2.
		// This is expected to be enabled for setups w/ concomitant GP
		// and FI55/06, and hence some redundancy.
		if (nhits[igr]<8) {
		  ok |= jgr; nFailures++; if (nhits[igr]==7) nRescues++;
		}
	      }
	      else if (TOpt::ReMode[18]&0x8) { // Double bridging requested
		// => One can be more demanding, even if that leads to the
		// loss of the 0x3 or 0x6 track. In case of a 0x6 one, it will
		// get a chance to be reco'd in "DoubleBridge".
		// => Consider it a failure for a FI05*2+GP*4 track to loose
		// one of its FI05 hits.
		// Note that we make the assumption that "nhits[1]==5" means
		// necessarily a FI05*2+GP*4 track.
		if (nhits[igr]<6) {
		  const TDetect &dRemoved = setup.vDetect(vecHit[jcom].IPlane);
		  bool FI05Removed = dRemoved.Name.find("FI05")==0;
		  if (nhits[igr]<5 || FI05Removed) {
		    ok |= jgr; nFailures++;
		    if (nhits[igr]==4 && !FI05Removed) nRescues++;
		  }
		}
	      }
	      else {
		if (nhits[igr]<5) {
		  ok |= jgr; nFailures++; if (nhits[igr]==4) nRescues++;
		}
	      }
	      break;
	    case 2: if (nhits[igr]<4) { ok |= jgr; nFailures++; }
	      break;
	    default:
	      CsErrLog::msg(elFatal,__FILE__,__LINE__,
		"TTrack Id = %d has scifi-only bit pattern = %x!",t.Id,t.Scifi);
	    }
	  }
	  else {
	    if      (nhits[igr]<TOpt::iPRpar[igr*10+2+3]) {
	      ok |= jgr; nFailures++;
	      if (nhits[igr]==TOpt::iPRpar[igr*10+2+3]-1) nRescues++;
	    }
	  }
	}
      }
      // Rescue zone when only 1 hit below par. And spare track if <=1 rescue
      if (nFailures==1 && nRescues==1 &&
	  // Special case where few hits in zone 0x1:
	  // This concerns, e.g., VSAT-(almost)only tracks w/ 4 hits in zone
	  // 0x2 (hence passing rescue condition infra), but which 4 hits
	  // are all GP02, which does not constrain enough the tracking.
	  // => Prevent them from being rescued
	  (!(t.Type&0x2) || nStations[1]>=2)) ok = 0;
      if (ok==0) {

	// ********** IF TRACK W/OUT TOO MANY SHARED HITS
	//            => REBUILD IT W/OUT THE PROBLEMATIC HIT

	THit& h = vecHit[jcom];
	if (TOpt::Print[4])
	  cout << " * TEvCleanTackList: Clean track " << t.GetId() << endl;
	t.SubHit(h); t.Clip();
	if (t.NGroups()==1) {
	  // If the downstream piece of a SM2 yoke track (single group and
	  // "TTrack::IMark>0"): save c/p...
	  double cop = t.IMark>0 ? t.Hlast(5) : 0;
	  t.QuickKF(-1,0); t.QuickKF(1,0);
	  if (t.IMark>0) { t.Hfirst(5)=t.Hlast(5) = cop; } // ...and restore it
	  if (t.IFit&0x8) t.QNewtonFit(1,4);
	}
	else {
	  int iFitSave = t.IFit; t.IFit = 0;
	  if (t.Type!=0x7 &&  // Not much we can do: QN fit will perform badly
	      // on this kind of track, if it ever converges. A possibility
	      // would be to QN fit 1st, SubHit and QN refit. And update chi2 w/
	      // the diff in QN chi2.
	      TOpt::ReMode[14] && h.IPlane<TLattice::IplLast) {
	    // QN Fit: If track outside dico, QN will exit w/o changing
	    // anything to the track's chi2 => The rebuilt track will not
	    // have been updated, but too bad...
	    double chi2Save = t.Chi2tot;
	    bool ret; if (t.IFit&0x8) ret = t.QNewtonFit(1,4);
	    else                      ret = t.QNewtonFit(1,1);
	    if (!ret ||
		// Check still that QN fit is meaningful...
		mode && t.NDics<t.NHits*.8 ||
		// ...and, if dp exists, check QN's p against previous value
		t.Hlast(5,5) && fabs(t.Haux(5)-t.Hlast(5))>5*t.Hlast(5,5)) {
	      t.Chi2tot = chi2Save; t.IFit = 0;
	    }
	  }
	  if (t.IFit!=0x8) t.Chi2tot *= float(t.NDFs-5)/(t.NDFs-4);
	}
	t.UseHitTime();
      }
      if (ok) {

	// ********** IF NOT OR IF REBUILDING UNSUCCESSFUL  => DELETE

	if (TOpt::Print[4])
	  cout << " * TEvCleanTackList: Erase track " << t.GetId() << endl;

	//if (Coral::Instance()->isAMonteCarloEvent()) {
	//  if (t.IKine>=0)                // Erase reference to track in Kine
	//  vecKine[t.IKine].subTrackID(t.Id);
	//
	
	if (t.Type==0x7) {
	  // ***** ...BUT TRY TO RECOVER PART of 0x7 TRACKS *****
	  if      ((ok&0x3)==0x3 || (ok&0x1)) {
	    t.Behead(setup.vIplLast()[0]);
	    t.Type &= 0x1e; t.Scifi &= 0x1e1e;  // Update attributes
	    if (TOpt::ReMode[46]) {
	      t.InsertMissedPlanes(); ok = t.FullKF(-1) ? 0 : 1; t.IFit = 0x21;
	    }
	    else {
	      ok = t.QuickKF(-1,1) ? 0 : 1; t.IFit = 0x3; // Canceling 0x10
	    }
	  }
	  else if (ok==0x4) {
	    t.Shorten(setup.vIplFirst()[2]);
	    ok = 0; t.Type &= 0x13; t.Scifi &= 0x1313;  // Update track's Type
	    // The momentum inherited from the 0x7 track fits its 0x6 piece (
	    // cf. TEv::BridgeSegments2, where upon extending a 0x3 bridging
	    // into the 0x4 zone it's the 0x6 momentum that is retained) and
	    // not its 0x3 piece. => Restore the latter's momentum.
	    t.Hlast(5) = t.Haux(5);
	    if (TOpt::ReMode[46]) {
	      t.InsertMissedPlanes(); ok = t.FullKF(-1) ? 0 : 1; t.IFit = 0x21;
	    }
	    else {
	      ok = t.QuickKF(-1,1) ? 0 : 1; t.IFit = 0x3; // Canceling 0x10
	    }
	    //#define CleanTL_DEBUG_REFIT
#ifdef CleanTL_DEBUG_REFIT
	    double rF = t.Hfirst(5), rL = t.Hlast(5), dr = sqrt(t.Hfirst(5,5));
	    if (fabs(rF-rL)>dr) {
	      double oldChi2 = t.Chi2tot/(t.NDFs-5);
	      t.QuickKF(1,1); t.QuickKF(-1,1);
	      CsErrLog::msg(elError,__FILE__,__LINE__,
  "Track #%d refit after shortening: p(1st,last),chi2 %.3f+-%.3f,%.3f,%.2f -> %.3f+-%.3f,%.3f,%.2f",
     1/rF,dr/rF/rF,1/rL,oldChi2,
     1/t.Hfirst(5),sqrt(t.Hfirst(5,5))/t.Hfirst(5)/t.Hfirst(5),t.Hlast(5),t.Chi2tot/(t.NDFs-5));
	    }
#endif
	    // Preserving 0x10: it applies to track's 0x1-segment and hence it's
	    // still valid. Although: not 100% valid: 0x1-segment may have lost
	    // some of its hits.
	    t.IFit |= 0x6;
	  }
	}
	// ***** ...AND SM2/muWall TRACKS WHEN  *****
	//     - NON SCIFIS-ONLY
	//     - SM1/SM2 COUNTERPART IS SCIFIS-ONLY
	// (This is sensible for scifis-only segments in the 0x2 zone can have
	// very few hits (as few as 4) so that they do not add much to the total
	// chi2 and can hence be bridged w/ almost anything and still yield a
	// reasonable, i.e. passing dCut[9], chi2 when evaluated w/o MS.
	// And this is useful for those track segments in the 0x4 zone can be
	// associated to calorimeter clusters and therefore eliminate the
	// possibility for the associated cluster to correspond to a neutral.)
	else if (t.Type==0x6 && (ok&0x2) &&
		 t.Scifi && (t.Scifi&0x202)==t.Scifi) {
	  t.Behead(setup.vIplLast()[1]);
	  ok = 0; t.Type = 0x4; t.Scifi &= 0x4; // Update track's Type
	  t.QuickKF(1,0); t.QuickKF(-1,0); t.IFit = 0x1; // Canceling 0x10
	}
#define CleanTL_RESCUE_0x3
#ifdef CleanTL_RESCUE_0x3
	else if (t.Type==0x3 && (ok&0x1) &&
		 // Do not care to rescue VSAT-(almost)only 0x2, for these are
		 // not very reliable
		 !(t.Scifi&0x202) &&
		 // Restrict oneself to tracks liable to reach the calorimeter,
		 // for this will be the only usefulness of these tracks, that
		 // are bound to remain momentumless: account for calo clusters.
		 t.Hlast(0)>750) {
	  t.Behead(setup.vIplLast()[0]);
	  if ((int)t.NDFs>TOpt::iPRpar[12]) {
	    // I.e. better than criterion used in 1st iteration of PR.
	    t.QuickKF(1,0); t.QuickKF(-1,0);
	    if (t.Chi2tot/(t.NDFs-4)<TOpt::dCut[9]) {
	      // I.e. a somewhat loose chi2 cut, to compensate for
	      //  multiscattering being neglected in this momentless fit.
	      ok = 0; t.Type = 0x2; t.Scifi &= 0x2; // Update track's Type
	      t.IFit = 0x1; // Canceling 0x10
	    }
	  }
	}
#endif
	if (ok) listTrack.erase(it);
	else t.UseHitTime();
      }
    }
  } while (jcom>=0);

}    

// ***************************************************************************
// ***************************** CleanTrackListSI ****************************
// ***************************************************************************
void TEv::CleanTrackListSI()
{
  // - Dedicated to SI hits.
  // - Called later in the TraFDic flowchart.
  // - Simplified algorithm: We do not proceed as supra by singling out the
  //  worst case first: hits are considered in their TEv ordering.

  //   ************************* INITIALIZATIONS *************************

  const TSetup &setup = TSetup::Ref(); int iplL0 = setup.vIplLast()[0];
  //  Indices of SI spectrometer (as opposed to beam) TStations
  static int iS0SI = -1, iS1SI, iPlLSI; if (iS0SI<0) {
    int iplF0 = setup.vIplFirst()[0];
    const vector<TStation> &stations = setup.vStation();
    int iS; for (iS = 0, iS1SI = -1; iS<(int)stations.size(); iS++) {
      const TStation &station = stations[iS];
      int ipl = station.IPlanes[0]+station.IPlanes.size()-1; // TStation's last
      if (iplL0<ipl) break;    // Consider zone 0x1 only
      if (ipl<iplF0) continue; // Skip upstream of target
      int type = station.Type; if (type==21) {
	if      (iS0SI<0) iS0SI = iS;
	else if (iS1SI<0) iS1SI = iS;
	else {            iS0SI = -1; break; } // No more than 2 SI TStations
	iPlLSI = ipl;
      }
    }
    if (iS0SI<0 || iS1SI<0) CsErrLog::mes(elFatal,
 "# of Si stations downstream of target !=2, while \"TraF ReMode[44]\" booked");
  }

  list<TTrack>::iterator it;
  int ipl; for (ipl = setup.vIplFirst()[0]; ipl<=iplL0; ipl++) {
    const TPlane  &p = setup.vPlane(ipl);
    const TDetect &d = setup.vDetect(p.IDetRef);
    if (d.IType!=21) break;   
    if (d.X(0)<setup.TargetCenter[0]) continue;

    //  *************** LOOP ON SI PLANES DOWNSTREAM of TARGET ***************

    vector<int>::const_iterator ihit, jhit;
    for (ihit = p.vHitRef().begin(); ihit!=p.vHitRef().end(); ihit++) {
      THit &h = vecHit[*ihit];
      if (h.sTrackID().size()<2 ||
	  h.sDigits().size()>=5 && // Allow some sharing if cluster's size large
	  h.sTrackID().size()<3) continue;

      //           ********** LOOP ON ALL SHARED HITS **********

      vector<TTrack*> ts; vector<double> chi2Incrs; int best;
      set<unsigned int>::iterator id;
      for (id = h.sTrackID().begin(), best = -1, it = listTrack.begin();
	   id!=h.sTrackID().end(); id++) {

	//            ***** LOOP ON ASSOCIATED TRACKS *****
	//             => DETERMINE BEST CHI2 INCREMENT

	for (;it!=listTrack.end(); it++) { // Get TTrack corresponding to ID
	  TTrack *t = &(*it); if (t->Id!=*id) continue; if (t->IMark==-1) break;
	  map<int, double>::iterator idx;
	  if ((idx = t->mChi2.find(ipl))==t->mChi2.end()) {
	    if (t->Type==0x1) {
	      chi2Incrs.push_back(-1); ts.push_back(t);
	    }
	    else if (t->Hfirst(5)==0)
	      CsErrLog::msg(elError,__FILE__,__LINE__,
     "Track #%d has no momentum while of type 0x%x",t->Id,t->Type);
	    else if (t->mChi2.empty())
	      CsErrLog::msg(elError,__FILE__,__LINE__,
     "Track #%d(0x%x) momentum but no chi2 map",t->Id,t->Type);
	    else
	      CsErrLog::msg(elError,__FILE__,__LINE__,
     "THit #%d track list and Track #%d chi2 map uncompatible",h.IHit,t->Id);
	  }
	  else {
	    double chi2Incr = (*idx).second;
	    if (best<0 || chi2Incr<chi2Incrs[best]) best = (int)ts.size();
	    chi2Incrs.push_back(chi2Incr); ts.push_back(t);
	  }
	  break;
	} // End loop on tracks in "listTrack"
      } // End loop on associated track
      if (best>=0) {                            // ***** BEST CHI2 INCREMENT...
	double chi2Incr = chi2Incrs[best];
	double chi2 = ts[best]->Chi2tot/(ts[best]->NDFs-5);
	if (chi2>TOpt::dCut[16]) { // ***** ...CORRESPONDS TO BAD TOTAL CHI2?...
	  for (int jt = 0; jt<(int)ts.size(); jt++) if (jt!=best) {
	    TTrack *t = ts[jt]; double chi2Incrp = chi2Incrs[jt];
	    double chi2p = t->Chi2tot/(t->NDFs-5);
	    if (chi2p<TOpt::dCut[16] && chi2p<chi2/2 &&
		chi2Incrp<chi2Incr*2) {
	      best = jt; break; // ***** ...RETAIN GOOD CHI2 ALTERNATIVE INSTEAD
	    }
	  }
	}
	for (int jt = 0; jt<(int)ts.size(); jt++) if (jt!=best) {

	  //   ***** LOOP ON ALL OTHER (non best) ASSOCIATED TRACKS *****

	  TTrack *t = ts[jt]; if (t->IMark<0) continue;
	  double chi2 = t->Chi2tot/(t->NDFs-5);
	  int nH0x1, nPix, nSI0s, nSI1s; list<int>::const_iterator jh;
	  for (jh = t->lHPat().begin(), nH0x1=nPix=nSI0s=nSI1s = 0;
	       jh!=t->lHPat().end(); jh++) {
	    if (*jh<0) continue; if (*jh==*ihit) continue;
	    int jpl = vecHit[*jh].IPlane; if (jpl>iplL0) break;
	    const TPlane &pj = setup.vPlane(jpl); int jS = pj.Station->IStation;
	    nH0x1++; if (pj.IFlag&0x30) /* P-pixel plane */ { nH0x1++; nPix++; }
	    if      (jS==iS0SI) nSI0s++;
	    else if (jS==iS1SI) nSI1s++;
	  }
	  int subPar = 0;
	  if (nPix<2 && nH0x1<8 ||
	      !(t->Scifi&0x100) && nH0x1<TOpt::iPRpar[5])
	    subPar = (t->Type&0x2) ? -2 : -3;
	  else if (nSI0s+nSI1s<3 || nSI0s*nSI1s && nSI0s+nSI1s<=5)
	    subPar = -1;
	  bool ret = true; THit *hbest = 0;
	  if (subPar==-2 && chi2>TOpt::dCut[16])  // ***** SUBPAR + BAD CHI2...
	    t->IMark = -2;                        // ***** ...MARK FOR BEHEADING
	  else if (chi2<TOpt::dCut[16] || subPar) {

	    //   ***** SUBPAR or GOOD CHI2: SEARCH ALTERNATIVE HIT/RESCUE *****

	    static double p0x11, dp0x11;  // ***** Preliminary step:
	    if (t->Type==0x11) {          // ***** Track bridged over target...
	      // ...There's not enough analysing power, even if target in dipole
	      // mode, for any precise determination of the momentum.
	      // => Fix P (and dP) ahead of possible refit...
	      p0x11 = t->Hfirst(5); dp0x11 = t->Hfirst(5,5); // ...and backup 
	      t->Hfirst(5,5) = 1.E-20;
	    }

	    double best; for (jhit = p.vHitRef().begin(), best = 6*d.Resol;
			      jhit!=p.vHitRef().end(); jhit++) {
	      if (jhit==ihit) continue;
	      THit &hp = vecHit[*jhit]; if (!hp.sTrackID().empty()) continue;
	      double diff = fabs(hp.U-h.U); if (diff<best) {
		hbest = &hp; best = diff;
	      }
	    }
	    if (hbest) {
	      const THit *hPrv = t->ReplaceHit(hbest);
	      if ((ret = t->FullKF(1))) {// ***** ALTERNATIVE HIT: REFIT FORWARD
		double chi2p; if ((chi2p = t->Chi2tot/(t->NDFs-5))>chi2*1.5 &&
				  chi2p>TOpt::dCut[16]/2) hbest = 0;
		else {
		  nH0x1++; int iS = p.Station->IStation; subPar = 0;
		  if      (iS==iS0SI) nSI0s++;
		  else if (iS==iS1SI) nSI1s++;
		  if (nPix<2 && nH0x1<8 ||
		      !(t->Scifi&0x100) && nH0x1<TOpt::iPRpar[5])
		    subPar = (t->Type&0x2) ? -2 : -3;
		  else if (nSI0s+nSI1s<3 || nSI0s*nSI1s && nSI0s+nSI1s<=5)
		    subPar = -1;
		}
	      }
	    }
	    if (!hbest && ret &&               // ***** NO or BAD ALTERNATIVE...
		!subPar) {  // (subPar cases are handled w/ infra) 
	      t->CancelHit_Clip(ipl);      // ...=> REMOVE HIT and REFIT FORWARD
	      ret = t->FullKF(1);
	    }
	    if (ret && !subPar) ret = t->FullKF(-1);     // ***** REFIT BACKWARD
	    if (!ret) t->IMark = -3;       // ***** MARK KF FAILURE FOR DELETION
	    else      t->IMark = subPar;   // ***** MARK SUBPAR FOR BEHEADING

	    if (ret && t->Type==0x11) { // Track bridged over target...
	      // ...Momentum may have varied slightly (despite being ''fixed'').
	      // => Resore it, as well as related covariance terms.
	      t->Hfirst(5)=t->Hlast(5) = p0x11;
	      t->Hfirst(5,5)=t->Hlast(5,5) = dp0x11;
	      for (int i = 1; i<5; i++) {
		t->Hfirst(5,i)=t->Hfirst(i,5)=t->Hlast(5,i)=t->Hlast(i,5) = 0;
	      }
	    }
	  }
	} // End loop on worst tracks 
      } // End if best track exists
    } // End loop on hits in plane
    it = listTrack.begin(); while (it!=listTrack.end()) {

      //     **** TRACKS MARKED FOR DELETION or UPDATE *****

      TTrack &t = *it;
      if (t.IMark==-3) listTrack.erase(it++);
      else if (t.IMark<0) {

	//       ***** TRACKS SUB-PAR in 0x1 zone *****

	bool ret = true;
	if (t.IMark==-2) { t.Behead(iplL0); t.Type &= 0x3fe; }
	else {
	  unsigned int projs; int nSpacePts, nProjs; const TStation *sPrv;
	  list<int>::const_iterator jh;
	  for (jh = t.lHPat().begin(), nSpacePts=nProjs = 0, projs = 0,
		 sPrv = 0; jh!=t.lHPat().end(); jh++) {
	    if (*jh<0) continue;
	    int jpl = vecHit[*jh].IPlane; if (jpl>iplL0) break;
	    const TPlane &p = setup.vPlane(vecHit[*jh].IPlane);
	    const TStation *s = p.Station; if (s->Type==21) continue;
	    if (s!=sPrv) {
	      sPrv = s; if (nProjs>=2) { nSpacePts++; nProjs = 0; projs = 0; }
	    }
	    unsigned int proj = 1<<p.IProj;
	    if (!(proj&projs))    { projs |= proj; nProjs++; }
	    if ((p.IFlag&0x30) &&   // P-pixel plane: Add an extra proj...
		// ...arbitrarily chosen as Y or Z proj. (or both...)
		(0x3&projs)!=0x3) { projs |= 0x3;  nProjs++; }
	  }
	  if (nProjs>=2) nSpacePts++;
	  if      (nSpacePts<3 && t.Type==0x1) ret = false;
	  else if (t.Type&0x10) {// Special case of track bridged over target...
	    // ...=> One has to take care not to behead it: would remain nothing
	    // (or very little if beheading up to SI)!
	    // ...We could think of removing only the SIs piece in the case
	    // there are many space points. But there is no method yet to handle
	    // this (viz. removing a central segment and keep head and tail).
	    //  => For the time being, let's simply shorten the track.
	    t.Shorten(setup.vIplFirst()[0]); t.Type &= 0x30;
	  }
	  else if (nSpacePts<2) { t.Behead(iplL0); t.Type &= 0x3e; }
	  else {
	    t.Behead(iPlLSI);
	    CsErrLog::msg(elError,__FILE__,__LINE__,
 "Track #%d (chi2/NDF=%.2f) stripped of its SIs",t.Id,t.Chi2tot/(t.NDFs-5));
	  }
	}
	if (ret) {
	  if     ((t.Type&0x3)==0x3 || (t.Type&0x6)==0x6) {
	    ret = t.FullKF(1); if (ret) ret = t.FullKF(-1);
	  }
	  else if (t.Type&0x10) { // Track bridged over target.
	    // (Possibly now stripped of their 0x1 piece.)
	    // => Fix P (and dP) for the time of the refit.
	    double p0x11 = t.Hfirst(5), dp0x11 = t.Hfirst(5,5);
	    t.Hfirst(5,5) = 1.E-20;
	    ret = t.FullKF(1); if (ret) t.FullKF(-1); if (ret) {
	      t.Hfirst(5)=t.Hlast(5) = p0x11;
	      t.Hfirst(5,5)=t.Hlast(5,5) = dp0x11;
	      for (int i = 1; i<5; i++) {
		t.Hfirst(5,i)=t.Hfirst(i,5)=t.Hlast(5,i)=t.Hlast(i,5) = 0;
	      }
	    }
	  }
	  else {
	    ret = t.QuickKF(1,0); ret = t.QuickKF(-1,0);
#ifdef FRINGE_FIELD_KF
	    if ((t.Type==0x1) &&
		t.QNewtonFit(1,1) &&                 // QN fit successful and...
		(TOpt::dCut[66]<0 || 
		 fabs(1/t.Haux(5))<TOpt::dCut[66])) {// ...low enough momentum
	      t.Hfirst(5) = t.Haux(5);               // => KF full fit
	      t.Hfirst(5,5) = 1.E-4;    // w/ large initial momemtum error
	      ret = t.FullKF(1); if (ret) t.FullKF(-1);
	      //printf("CTSI: Evt #%d, 0x1-track #%d fitted w/ P\n",event,t.Id);
	    }
#endif
	  }
	}
	if (!ret) listTrack.erase(it++);
	else { t.IMark = 0; it++; }
      }
      else it++;
    }
  } // End loop on SI planes
}
