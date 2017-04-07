// $Id: TEvTracksRefit.cc 14069 2015-09-17 20:44:46Z lsilva $

/*!
   \brief Refit tracks using info from vertexing. So far only event time info is
   used. Actual refit taking place depends upon:
    - event's trigger,
    - magnitude of departure from event time = 0,
    - whether event time was already taken into account in TEv::TracksFit2.
   \author  Yann.Bedfer@cern.ch
*/

#include "TH1.h"
#include "TH2.h"
#include "CsErrLog.h"
#include "CsInit.h"
#include "DaqDataDecoding/DaqOption.h"
#include "CsDetector.h"
#include "CsEvent.h"
#include "CsTrack.h"
#include "CsHelix.h"
#include "TEv.h"
#include "TDisplay.h"

using namespace std;

/*
  This method is specific to "lattice" alternative of "traffic".
*/

bool TEv::TracksRefit(CsVertex *pVertex, list<CsTrack*> &tracksToBeDeleted)
{

  static TH2D *hdChi2 = 0, *hdChi2InV;
  static bool first = true; if (first) { // ***** BOOK HISTOS *****
    first = false;
    if (TOpt::Hist[9]) {
      CsHistograms::SetCurrentPath("/Traffic");
      hdChi2    = new TH2D("hdChi2",   "#Delta#chi^{2} upon Refit vs. evt Time",
			100,-5,5,8,0,8);
      hdChi2InV = new TH2D("hdChi2InV","#Delta#chi^{2} vs. evtT (in Vertex)",
			100,-5,5,8,0,8);
      CsHistograms::SetCurrentPath("/");
    }
  }

  //       ***** REQUIRE THAT TRIGGER MATCHES ReMode[37]... *****
  const unsigned int allTrigs = 0xffff; // Although TCS can handle 0x7fffff, it was decided to limit max. #triggers to 16. And take advantage of the fact to make use of the other bits of the trigger TDC module.
  unsigned int evtTrig = ptrEvt()->getAlterMasterTrigger()<0 ? trig_mask :
    1<<ptrEvt()->getAlterMasterTrigger();
  evtTrig &= allTrigs;  // Cut away trailing end bits (i.e. online filter...)
  unsigned int inclMask =// Trigger pattern for which bits we require only an...
    TOpt::ReMode[37]&(~TOpt::ReMode[40]);// ...inclusive match, cf. UpdateDrifts
  if (( evtTrig&TOpt::ReMode[37])!=evtTrig && // Strict match 
      !(evtTrig&inclMask) ||                  // Inclusive match
      // The following can happen in MC: event did not pass the trigger...
      evtTrig==0) // => Exclude the case. No possible trigger offset anyway.
    // Re-fitting is, per "CsAverPattern::doPattern"'s construction, forbidden
    // while re-tracking is on: no need to care about it here.
    return 0;

  //  ********** REQUIRE EVENT TIME TO DIFFER SIGNIFICANTLY FROM 0... **********
  list<CsTrack*> vTracks = pVertex->getTracks();
  double beamTime = vTracks.front()->getMeanTime();
  //#define TrkReF_DEBUG
#ifdef TrkReF_DEBUG
  static int debug = 0; if (debug) {
    list<CsTrack*>::iterator icst;
    for (icst = vTracks.begin(); icst!=vTracks.end(); icst++) {
      CsTrack *cst = *icst; int nHits = cst->getClusters().size();
      printf("%4d  %3d %5.2f  %6.2f  %6.2f\n",cst->getId(),nHits,
	     cst->getChi2()/(nHits-5),
	     1/cst->getHelices()[0].getCop(),cst->getMeanTime());
    }
  }
#endif
  if (fabs(beamTime)<TOpt::dCut[71] &&// Event (=beam's) time ~=0 => No refit...
      (TOpt::ReMode[39]==0 ||      //...unless decoalescence requested AND...
       // ...event time wasn't already considered in the primary fit (presumably
       // because several beam candidates) which prevented the decoalescence
       // from being performed.
       evtTConsidered)) return 0;

  //  ********** HAS EVENT TIME ALREADY BEEN TAKEN INTO ACCOUNT **********
  // This may have been done by "TEv::TracksFit2". If so, it has  been done
  // undiscriminatedly, i.e. (almost) all tracks, whether actually associated
  // to the beam particle from which event time was derived or associated to an
  // accidentally coincident particle which beam track would have been missed
  // (still excepting tracks off-time enough for their T0 to have been derived
  // from their own track's time). This matters only if the beam which was
  // considered is actually associated to an interaction vertex, to the extent
  // that reco'ing accidentals is useful (e.g. because they can be made to
  // account for calorimeter clusters). In all rigor, this kind of
  // not-too-off-time accidental could be reexamined now, and made to benefit
  // from the more rigorous selection performed by the present method. But we
  // decide not to do it.
  double evtT = eventTime ? eventTime : eventTRef;
  if (evtTConsidered &&
      fabs(evtT-beamTime)<TOpt::dCut[71]/2) return 0;
  // (Comment a posteriori: I don't remember why I have introduced the 2nd term
  // in the above condition. In principle, since "eventTime" is set = beam time
  // in TEv::EventTime, and beam time cannot have changed ever since, it is
  // systematically fulfilled.)

  eventTime = beamTime;

  list<CsTrack*> eTracks = ptrEvt()->getTracks();
  list<TTrack>::iterator it;

  const TSetup &setup = TSetup::Ref();
  double XfZ0 = setup.iPlane2Detect(setup.vIplLast() [0]).X(0);
  double XiZ1 = 0, XfZ1 = 0;
  if (setup.vIplFirst().size()>=2) {
    XiZ1 = setup.iPlane2Detect(setup.vIplFirst()[1]).X(0);
    XfZ1 = setup.iPlane2Detect(setup.vIplLast() [1]).X(0);
  }
  //                 ***** MISSING MOMENTUM *****
  // It's later used to check whether a track not associated to pVertex can be
  // considered part of the primary interaction w/o infringing P conservation.
  list<CsTrack*>::iterator icst; double upperPx, sDPx2;
  for (icst = vTracks.begin(), upperPx=sDPx2 = 0; icst!=vTracks.end(); icst++) {
    const CsTrack *cst = *icst;
    map<int, int, less<int> >::const_iterator im =
      mapCsTrId2TrId.find(cst->getId());
    if (im==mapCsTrId2TrId.end()) CsErrLog::msg(elFatal,__FILE__,__LINE__,
  "CsTrack %d has no counterpart in TraFFic's event object",cst->getId());
    unsigned int id = im->second; const TTrack *t;
    for (it = listTrack.begin(), t = 0; it!=listTrack.end(); it++) {
      if ((*it).Id==id) {
	t = &(*it); const THlx &h = t->Hfirst; double cop = h(5);
	if (!cop) CsErrLog::msg(elFatal,__FILE__,__LINE__,
  "CsTrack #%d in pVertex and yet w/ momentumless counterpart #%d in TraFFic",
				cst->getId(),t->Id);
	double p = 1/fabs(cop);
	double ty = h(3), tz = h(4), w2 = 1+ty*ty+tz*tz, w = sqrt(w2);
	double px = p/w, py = ty*px, pz = tz*px;
	// dPx/dtgy = - Py/w^2; dPx/dtgy = -Pz/w^2; dPx/d1/P = -Px.P
	double dpxdty = -py/w2, dpxdtz = -pz/w2, dpxdr = -px*p;
	if (upperPx) upperPx -= px;
	else         upperPx += px;
	sDPx2 += dpxdty*dpxdty*h(3,3)+dpxdty*dpxdty*h(4,4)+dpxdr*dpxdr*h(5,5)+
	  2*(dpxdty*dpxdtz*h(3,4)+dpxdty*dpxdr*h(3,5)+dpxdtz*dpxdr*h(4,5));
	break;
      }
    }
    if (!t) CsErrLog::msg(elFatal,__FILE__,__LINE__,
  "CsTrack %d counterpart in TraFFic (ID=%d) does not exist",cst->getId(),id);
  }
  upperPx += 3*sqrt(sDPx2); // Increment missing Px by 3 sigmas

  //                ********** modified = RETURNED  FLAG **********
  // 0: No modif, 1: At least ONE track refitted (be it really modified or not)
  // (We could have singled out the case where the sole track refitted is an
  // off-time track NOT in vertex. Meaning that, upon return in the vertex
  // package, one would have only to update the list of tracks, without having
  // to refit the primary vertex. But this is too much of a complication for
  // a limited gain in CPU.)
  bool modified = false; it = listTrack.begin(); while (it!=listTrack.end()) {
    TTrack &t = *it;
    if (t.Type&0x10) { it++; continue; }  // Exclude beam tracks

    // N.B.: So far, bypass the case of momentumless tracks. They could also
    //      benefit from the redefinition of the event time, but the code
    //      infra assumes implicitly that momentum is defined (call to
    //      "TTrack::FullKF", and, above all, access to "TTrack::mChi2")
    if (t.Type==0xc) { it++; continue; }  // Exclude muWall momentumless
    if (!t.Hfirst.with_mom()) {           // Exclude ANY momemtumless track
      it++; continue;
    }
    int id = t.Id, csId = -1;

    //     *************** LOOP ON ALL TRACKS w/ MOMENTUM ***************
    
    //    ********** IS TRACK ASSOCIATED TO INTERACTION? **********

    // ***** YES IN ANY CASE IF ASSOCIATED TO THE VERTEX *****
    bool inVertex; for (icst = vTracks.begin(), inVertex = false;
			icst!=vTracks.end(); icst++) {
      const CsTrack *cst = *icst;
      map<int, int, less<int> >::const_iterator im =
	// Note: it's checked supra that "im" does exist.
	mapCsTrId2TrId.find(cst->getId());
      if (im->second==id) { inVertex = true; csId = im->first; break; }
    }

    if (t.NGroups()==1 && !inVertex) {// ***** EXCLUDE SINGLE-ZONE OUT OF VERTEX
      // (For the time being, at least... We intend to include them in a future
      // release. Problem is that we would need a better identification of those
      // single-zone tracks that can benefit from the refit, excluding all
      // tracks already  id'ed as halo muons, but also, possibly, those tracks
      // which have no time measuring hits and coud be better fitted by
      // adjusting their T0.
      it++; continue; }

    //                                                      ***** SM2 YOKE TRACK
    // For the time of the refit, the track's piece downstream of SM2 is
    // stripped away.
    TTrack *yokeTr; const THlx *yokeTrHlast;
    SplitYokeTr(t,yokeTr,yokeTrHlast);// Sets "yokeTr=0" for a standard track

    //        ********** CORRECT DRIFTS FOR EVENT TIME **********
    // (The correction relies on the time info stored in CsCluster's. So that,
    // whatever corrections were made earlier, either affecting only THit, or,
    // if "ReMode[29]&0x4", ported also to "CsCluster::U", the calculation is
    // still correct. Cf. also comments in "TEvTracksFit2.cc".)
    //       ********** and FOR PROPAGATION TIME (+RAY) **********
    // (N.B.: Propagation (+Xray) is applied independent of ReMode[29] for it's
    // not (yet) possible to "getCorrU" for the sole event time correction.)
    bool doCorr = true;
    bool deCoalesce =     // ********** DRIFT DECOALESCENCE **********
      TOpt::ReMode[39];

    // Status: -1: No init yet, 0: Init done,
    //          1: Current track modified
    int tStatus = -1;
    double tT = beamTime;

    // ***** 2 STEPS: I) CORRECTION (and possible deletion)
    //                  + LR ambiguity raise coalesced clusters   => REFIT
    //               II) If bad chi2, try worst hit's alternative => 2ND REFIT
    bool hitsPerZone;

// Commented down. Should skip track, not crash processing:
/*
    if (t.mHef.size()!=t.lHitPat.size() || t.mHeb.size()!=t.lHitPat.size())
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
		    "TTrack %d inconsistent mHef/b(=%d/%d) != lHitPat(=%d)",
		    t.Id,t.mHef.size(),t.mHeb.size(),t.lHitPat.size());
*/
    list<int>::iterator ih; map<int,THlx>::const_iterator im;
    for (ih = t.lHitPat.begin(), hitsPerZone = false; ih!=t.lHitPat.end();
	 ih++) {
      if (*ih<0) continue;
      THit &h = vecHit[*ih];
      int ipl = h.IPlane; const TPlane &p = setup.vPlane(ipl);
      const TDetect &d = setup.vDetect(p.IDetRef);
      CsDetector *det = d.PtrDet(); if (!det->hasDrift()) continue;
      if (tStatus<0) {                   // First detector encountered...
	tStatus = 0;                       // ... => Init track's status
	if (!modified) {
	  // Before first modif takes place, draw step by step graphics.
	  if (TOpt::Graph[0]>0 && TOpt::Graph[6]>0) {
	    cout<<"Refitting is starting"<<endl;
	    TDisplay::Ref().Draw(3);  // draw tracks
	    TDisplay::Ref().Draw(4);  // draw menu
	  }
	  modified = true;
	}
	if (!inVertex) {
	  // ********** NOT IN VERTEX: IS IT ACCIDENTALLY COINCIDENT? **********
	  // - The conditions are slightly relaxed w.r.t. "TEvTracksFit2.cc" (
	  //  no longer any cut on track's chi2), now  that one can tell whether
	  //  track is directly associated to vertex or not. (Of course the
	  //  association can still be indirect, via a decay vertex.)
	  // - One could be more rigorous, by looking for hits in the beam
	  //  telescope, consistent w/ the track: not yet done. In the mean
	  //  time:
	  //   - Determine if off-time. W/ a relaxed cut if it's in addition,
	  //    paraxial.
	  //   - Check for high momenta which could be excluded from primary
	  //    interaction because infringing momentum conservation.
	  double timeCut = fabs(t.Hlast(3))<.05 && fabs(t.Hlast(4))<.05 ?
			 /* paraxial */ 3 : 5;
	  bool accidental =
	    t.DispTime>0 && // This requires 2 time-measuring detectors...
	    // ...eliminating the risk to have track's time set by ghost hit.
	    // (Note: When time undef: SigmaTime<0.)
	    fabs(t.MeanTime-beamTime)/t.SigmaTime>timeCut;
	  if (!accidental) {
	    const THlx &h = t.Hfirst; double cop = h(5);
	    if (cop) {
	      double p = 1/fabs(cop);
	      double ty = h(3), tz = h(4), w2 = 1+ty*ty+tz*tz, w = sqrt(w2);
	      double px = p/w, py = ty*px, pz = tz*px;
	      // dPx/dtgy = - Py/w^2; dPx/dtgy = -Pz/w^2; dPx/d1/P = -Px.P
	      double dpxdty = -py/w2, dpxdtz = -pz/w2, dpxdr = -px*p;
	      double dpx =
		dpxdty*dpxdty*h(3,3)+dpxdty*dpxdty*h(4,4)+dpxdr*dpxdr*h(5,5)+2*
		(dpxdty*dpxdtz*h(3,4)+dpxdty*dpxdr*h(3,5)+dpxdtz*dpxdr*h(4,5));
	      accidental = px>upperPx+3*sqrt(dpx);
	    }
	  }
	  if (accidental) {       // Accidental: set latency "tT"...
	    if (0<t.SigmaTime &&t.SigmaTime<3)
	      tT = t.MeanTime;    // ...= track'time, if good enough resolution,
	    else
	      tT = 0;             // ...= 0, otherwise.
	    if (TOpt::Print[2])
	      printf("TTrack %d considered accidental, refit w/ T0 = %.1f\n",
		     t.Id,tT);
	  }
	}
      }
      // **********  I) CASE COALESCED HIT/MIRROR **********
      // - Try decoalesce (require good chi2)
      // - In any case, apply correction (even if coalesced, because of  X-Ray
      //  and latency (of off-time track, of trigger) corrections)
      // ********** II) STANDARD CASE **********
      //  => Possibly, coalesce hit/mirror that have then gone closer.
      //  => Erase hits that have gone out of time gate ("getCorrU" returns
      //    then an error)
      if (h.Mirror!=-1 && vecHit[h.Mirror].Status==-4) {
	//     ***** COALESCED HIT: TRY RAISE LR  if DISTANCE LARGE ENOUGH *****
	THit &hm = vecHit[h.mirror];
	CsCluster *c = h.PtrClus(), *cm = hm.PtrClus();
	if (TOpt::ReMode[29]&0x4) { // CsCluster's may have been modified =>...
	  if (h.iHit<hm.iHit) det->restoreU(c,cm);  // ...restore them
	  else                det->restoreU(cm,c);
	}
	double u= c->getU()/10, um = cm->getU()/10;
	// Which of next/previous TPlane associated to current one? 
	int jpl = p.Associate ? p.Associate->IPlane : -1;
	if (jpl==ipl+1) im = t.mHeb.find(ipl);
	else            im = t.mHef.find(ipl);
	const THlx &H = (*im).second; double y = H(1), z = H(2);
	double deltaU =
	  sqrt(H(1,1)*d.Ca*d.Ca+2*H(1,2)*d.Ca*d.Sa+H(2,2)*d.Sa*d.Sa)/d.Resol;
	if (t.Chi2tot/(t.NDFs-5)<TOpt::dCut[16] && deltaU<1 && deCoalesce) {
	  if (doCorr) {
	    bool error; double up = det->getCorrU(c,y*10,z*10,tT,error);
	    if (!error) {
	      double ump = det->getCorrU(cm,y*10,z*10,tT,error);
	      if (!error) { u = up/10; um = ump/10; }
	    }
	  }
	  if (fabs(um-u)>2.5*d.Resol) { // Distance > 0.5 * cut
	    double ue = y*d.Ca+z*d.Sa;
	    double du = fabs(u-ue)/d.Resol, dum = fabs(um-ue)/d.Resol;
#ifdef TrkF_DEBUG_DECOALESCE
	    static int jdebug = 0;
	    static int nOKs[3] = {0,0,0}; int ok = 0; char oks[] = "0+-0";
	    if (isMC) {
	      if      (h.IKine ==t.IKine || t.IKine!=-1 && h.IKin2 ==t.IKine)
		ok = 1;
	      else if (hm.IKine==t.IKine || t.IKine!=-1 && hm.IKin2==t.IKine)
		ok = 2;
	    }
#endif
	    if (du<2 && dum>3) {
	      h.u  = u;  h.sigu = sqrt(c->getCov()(1,1))/10.; 
	      h.Rotate();  hm.status = 1; tStatus = 1;
	      if (doCorr && (TOpt::ReMode[29]&0x4)) {
		//#define TrkF_DEBUG_CorrU
#ifdef TrkF_DEBUG_CorrU
		if (d.Name=="ST03X1ub") printf("TRF_CorrU %s %.2f -> %.2f\n",
		      d.Name.c_str(),c->getU(),u*10);
#endif
		c->setU(u*10);   // Update CsCluster: cf. explanations infra
	      }
#ifdef TrkF_DEBUG_DECOALESCE
	      if (isMC) {
		if (ok) nOKs[ok-1]++;
		else    nOKs[2]++;
		if (jdebug)
		  printf("DECOALESCE %c (%d,%d,%d): ID%d %4.1f %s dist,du,dum,deltaU %4.1f %4.1f %4.1f %4.1f\n",
	     	    oks[ok],nOKs[0],nOKs[1],nOKs[2],t.Id,t.Chi2tot/(t.NDFs-5),
	            d.Name.c_str(),(um-u)/d.Resol,du,dum,deltaU);
	      }
#endif
	    }
	    else if (dum<2 && du>3) {
	      hm.u = um; hm.sigu = sqrt(cm->getCov()(1,1))/10.;
	      hm.Rotate(); hm.status = 1; tStatus = 1;
	      t.ReplaceHit(&hm);
	      if (doCorr && (TOpt::ReMode[29]&0x4)) {
#ifdef TrkF_DEBUG_CorrU
		if (d.Name=="ST03X1ub") printf("TRF_CorrU %s %.2f -> %.2f\n",
		  d.Name.c_str(),cm->getU(),um*10);
#endif
		cm->setU(um*10);
	      }
#ifdef TrkF_DEBUG_DECOALESCE
	      if (isMC) {
		if (ok) nOKs[2-ok]++;
		else    nOKs[2]++;
		if (jdebug)
		  printf("DECOALESCE %c (%d,%d,%d): ID%d %4.1f %s dist,du,dum,deltaU %4.1f %4.1f %4.1f %4.1f\n",
		  oks[3-ok],nOKs[0],nOKs[1],nOKs[2],t.Id,t.Chi2tot/(t.NDFs-5),
	          d.Name.c_str(),(um-u)/d.Resol,du,dum,deltaU);
	      }
#endif
	    }
	  }
	}
	if (hm.status!=1 && // If ambiguity has not been raised...
	    doCorr) {
	  bool error; double up = det->getCorrU(c,y*10,z*10,tT,error);
	  if (!error) {
	    double ump = det->getCorrU(cm,y*10,z*10,tT,error); if (!error) {
	      u = (up+ump)/20; if (fabs(h.U-u)>d.Resol/10) {
		// ...check still for X-ray (or the like) correction
		h.u = u; h.sigu = sqrt(c->getCov()(1,1))/10+fabs(ump-up)/20;
		h.Rotate(); tStatus = 1;
		if (TOpt::ReMode[29]&0x4) { // Update CsCluster
#ifdef TrkF_DEBUG_CorrU
		  if (d.Name=="ST03X1ub") printf("TRF_CorrU %s %.2f -> %.2f\n",
		    d.Name.c_str(),c->getU(),h.U*10);
#endif
		  c->setU(h.U*10); c->setSigmaU(h.SigU*10);
		}
	      }
	    }
	  }
	}
      }
      else if (doCorr) {
	im = t.mHef.find(ipl);
	const THlx &H = (*im).second; double y = H(1), z = H(2);
	CsCluster *c = h.PtrClus(); bool error = true;
	double up = det->getCorrU(c,y*10,z*10,tT,error); if (!error) {
	  bool coalesced = false; if (h.Mirror!=-1 && tT>1) {
	    THit &hm = vecHit[h.Mirror];
	    double ump = det->getCorrU(hm.PtrClus(),y*10,z*10,tT,error);
	    if (fabs(ump-up)<30*d.Resol) {   // ***** CORR'ED Us CLOSE ENOUGH...
	      //                                  ...COALESCE HIT + MIRROR *****
	      h.u = (up+ump)/20; h.sigu = fabs(ump-up)/20;
	      h.sigu += sqrt(c->getCov()(1,1))/10; tStatus = 1;
	      hm.u = h.U; hm.sigu = h.SigU; coalesced = true;
	      if (isMC && hm.IKine>=0) {
		hm.Rotate(); t.ReplaceHit(&hm);
		h.status = -4; c = hm.PtrClus();
	      }
	      else {
		h.Rotate(); hm.status = -4;
	      }
	      if (TOpt::ReMode[29]&0x4) { // Update Csluster. N.B.: "c" may have
#ifdef TrkF_DEBUG_CorrU
		if (d.Name=="ST03X1ub") printf("TRF_CorrU %s %.2f -> %.2f\n",
		  d.Name.c_str(),c->getU(),h.U*10);
#endif
		c->setU(h.U*10); c->setSigmaU(h.SigU*10); // been redefined
	      }
	    }
	  }
	  bool doUpdate = !coalesced;
#define TrkF_CorrU_THRESHOLD 20
#ifdef TrkF_CorrU_THRESHOLD // Condition "CorrU" by cut on residual
	  doUpdate &= fabs(c->getU()-up)/10>d.Resol/TrkF_CorrU_THRESHOLD;
#endif
	  if (doUpdate) {
	    h.u = up/10; h.Rotate(); tStatus = 1;     // ***** UPDATE THit *****
	    if (TOpt::ReMode[29]&0x4) {
	      //                        ***** UPON OPTION: RESET CsCluster *****
	      // (For detailed comments on this option, cf. "TEvTrackFit2.cc".)
#ifdef TrkF_DEBUG_CorrU
	      if (d.Name=="ST03X1ub") printf("TRF_CorrU %s %.2f -> %.2f\n",
		d.Name.c_str(),c->getU(),up);
#endif
	      c->setU(up);
	    }
	  }
	}
	else if (fabs(tT)>5) {// "getCorrU" in error while latency is large...
	  // ...try and erase coresponding THit
	  static int nHits[5]; if (!hitsPerZone) { 
	    hitsPerZone = true; list<int>::iterator ihp;
	    for (ihp = t.lHitPat.begin(), memset((void*)nHits,0,5*sizeof(int));
		 ihp!=t.lHitPat.end(); ihp++) {
	      if (*ihp<0) continue;
	      int jpl = vecHit[*ihp].IPlane;
	      for (int igr = 0; igr<(int)setup.vIplFirst().size(); igr++) {
		if (setup.vIplFirst()[igr]<=jpl &&
		    jpl<=setup.vIplLast()[igr]) { nHits[igr]++; break; }
	      }
	    }
	  }
	  int igr; for (igr = 0 ; igr<(int)setup.vIplFirst().size(); igr++) {
	    if (setup.vIplFirst()[igr]<=ipl && ipl<=setup.vIplLast()[igr])
	      break;
	  }
	  if (nHits[igr]>TOpt::iPRpar[10*igr+5]) {
	    t.CancelHit(ipl); nHits[igr]--; tStatus = 1;
	  }
	}
      }  // End correcting propag. (and the like): (de)coalescence, deletion
    }  // End loop on hits
    if (tStatus>0) {      // ********** REFIT TRACK AND EXPORT **********
      CsTrack *cCst;        // ***** CsTrack COUNTERPART of CURRENT TTrack *****
      for (icst = eTracks.begin(), cCst = 0; icst!=eTracks.end(); icst++) {
	CsTrack *cst = *icst;
	if (cst->getId()==csId) { cCst = cst; break; }
	map<int, int, less<int> >::const_iterator im =
	  mapCsTrId2TrId.find(cst->getId());
	if (im==mapCsTrId2TrId.end())
	  CsErrLog::msg(elFatal,__FILE__,__LINE__,
	    "CsTrack %d has no counterpart in TraFFic's event object",
			cst->getId());
	if (im->second==id) { cCst = cst; break; }
      }
      if (cCst==0) CsErrLog::msg(elFatal,__FILE__,__LINE__,
	    "TTrack %d has no counterpart in coral CsEvent object",id);
      double oldChi2 = t.Chi2tot/(t.NDFs-5);          // ***** REFIT TRACK *****
      map<int,double> *mChi2 = 0;
      bool ret = t.FullKF(1); 
      if (ret) {
	if (TOpt::dCut[16] && t.Chi2tot/(t.NDFs-5)>.8*TOpt::dCut[16]) {
	  mChi2 = new map<int,double>;
	  *mChi2 = t.mChi2; // Save chi2 incr. of forward fit
	}
	ret = t.FullKF(-1);
      }
      double chi2 = t.Chi2tot/(t.NDFs-5);
      if (TOpt::ReMode[38] && t.NGroups()>1 && ret &&
	  TOpt::dCut[16] && chi2>.8*TOpt::dCut[16]) {

	//                       ***** CLEANING of TRACKS w/ DRIFTS HITS *****
	// Note: Tracks are subjected to cleaning only if w/ non null "tStatus".
	// Other they haven't been modified since last cleaned in "TracksFit2".
	vector<THit*> badHits, altHits;       // ***** BAD chi2 INCREMENTS *****
	THit *worstHit = 0, *altWorst = 0;   // ***** WORST chi2 INCREMENT *****
	THit *worstHiF = 0, *altWorsF = 0;

	if (fabs(t.Hfirst(5))>.5) {  // ***** LOW P: DEAL 1ST w/ WORST HIT *****
	  double worstChi2Incr, worstChi2IncF = 0;
	  list<int>::const_iterator ih; map<int, double>::iterator idx;
	  for (ih = t.lHPat().begin(), idx = t.mChi2.begin(),
		 worstChi2Incr = .25*t.Chi2tot; ih!=t.lHPat().end(); ih++) {
	    if (*ih<0) continue;
	    double chi2Incr = (*idx).second; if (chi2Incr>worstChi2Incr) {
	      THit &h = vecHit[*ih];
	      if (h.Mirror!=-1 && vecHit[h.Mirror].Status!=-4) {
		worstHit = &h; altWorst = &(vecHit[h.Mirror]);
		worstChi2Incr = chi2Incr;
	      }
	    }
	    idx++;
	  }
	  if (mChi2) {
	    for (ih = t.lHPat().begin(), idx = mChi2->begin(),
		   worstChi2IncF = .25*t.Chi2tot; ih!=t.lHPat().end(); ih++) {
	      if (*ih<0) continue;
	      double chi2Incr = (*idx).second; if (chi2Incr>worstChi2IncF) {
		THit &h = vecHit[*ih];
		if (h.Mirror!=-1 && vecHit[h.Mirror].Status!=-4) {
		  worstHiF = &h; altWorsF = &(vecHit[h.Mirror]);
		  worstChi2IncF = chi2Incr;
		}
	      }
	      idx++;
	    }
	  }
	  if (worstHit) {
	    int ipl = worstHit->IPlane; const TPlane &p = setup.vPlane(ipl);
	    const TDetect &d = setup.vDetect(p.IDetRef); double Xd = d.X(0);
	    CsDetector *det = d.PtrDet();
	    if ((fabs(Xd-XfZ0)<10 || fabs(Xd-XiZ1)<10) &&  // I.e. DC01||DC02
		worstChi2Incr>.25*t.Chi2tot || worstChi2Incr>.33*t.Chi2tot) {
	      //                       ***** WORST IS SIGNIFICANTLY WORST? *****
	      if (doCorr) {
		//           ***** PROPAGATION CORR. of WORST-HIT'S MIRROR *****
		im = t.mHeb.find(ipl);
		const THlx &H = (*im).second; double y = H(1), z = H(2);
		bool error;
		double ump = det->getCorrU(altWorst->PtrClus(),y*10,z*10,tT,error);
		if (!error) altWorst->u = ump/10;
	      }  // End propagation correction for mirror hit
	      altWorst->Rotate(); t.ReplaceHit(altWorst);
	      double newChi2 = 0; if ((ret = t.FullKF(1)) &&// ***** REFIT *****
				      (newChi2 = t.Chi2tot/(t.NDFs-5))<chi2) {
		ret = t.FullKF(-1); newChi2 = t.Chi2tot/(t.NDFs-5);
	      }
	      if (ret && !newChi2) CsErrLog::msg(elFatal,__FILE__,__LINE__,
  "Event #%d track #%d: KF fit successful and yet chi2 not updated",event,t.Id);
	      if (!ret || chi2<newChi2) {     // ***** NOT BETTER: RESTORE *****
		t.ReplaceHit(worstHit);
		ret = t.FullKF(1); if (ret) ret = t.FullKF(-1);
		if (!ret) {
		  CsErrLog::msg(elError,__FILE__,__LINE__,
  "Track #%d: KF fails after original track has been restored => erase.",t.Id);
		  listTrack.erase(it++); continue;
		}
	      }
	      else tStatus |= 0x2;
	      chi2 = t.Chi2tot/(t.NDFs-5);
	    }
	  }  // End candidate "worstHit" found
	  if (!(tStatus&0x2) && worstHiF && worstHiF!=worstHit) {
	    // If the original track was restored after the exchange of the
	    // ``backward'' worst hit, and ``forward'' differs: try it.
	    int ipl = worstHiF->IPlane; const TPlane &p = setup.vPlane(ipl);
	    const TDetect &d = setup.vDetect(p.IDetRef); double Xd = d.X(0);
	    CsDetector *det = d.PtrDet();
	    if ((fabs(Xd-XfZ0)<10 || fabs(Xd-XiZ1)<10) &&  // I.e. DC01||DC02
		worstChi2IncF>.25*t.Chi2tot || worstChi2IncF>.33*t.Chi2tot) {
	      if (doCorr) {
		im = t.mHef.find(ipl);
		const THlx &H = (*im).second; double y = H(1), z = H(2);
		bool error;
		double ump = det->getCorrU(altWorsF->PtrClus(),y*10,z*10,tT,error);
		if (!error) altWorsF->u = ump/10;
	      }  // End propagation correction for mirror hit
	      altWorsF->Rotate(); t.ReplaceHit(altWorsF);
	      double newChi2 = 0; if ((ret = t.FullKF(1)) &&// ***** REFIT *****
				      (newChi2 = t.Chi2tot/(t.NDFs-5))<chi2) {
		ret = t.FullKF(-1); newChi2 = t.Chi2tot/(t.NDFs-5);
	      }
	      if (ret && !newChi2) CsErrLog::msg(elFatal,__FILE__,__LINE__,
  "Event #%d track #%d: KF fit successful and yet chi2 not updated",event,t.Id);
	      if (!ret || chi2<newChi2) {     // ***** NOT BETTER: RESTORE *****
		t.ReplaceHit(worstHiF);
		ret = t.FullKF(1); if (ret) ret = t.FullKF(-1);
		if (!ret) {
	CsErrLog::msg(elError,__FILE__,__LINE__,
  "Track #%d: KF fails after original track has been restored => erase.",t.Id);
		  listTrack.erase(it++); continue;
		}
	      }
	      else tStatus |= 0x8;
	      chi2 = t.Chi2tot/(t.NDFs-5);
	    }
	  }  // End candidate forward "worstHiF" found
	}  // End try and cure worst hit

	if (ret && chi2>.8*TOpt::dCut[16]) {
	  //            ***** chi2 STILL BAD: DEALING NOW  w/ ALL BAD HITS *****
	  double chi2IncrMx = .2*.8*TOpt::dCut[16]*(t.NDFs-5); // 20% chi2 cut
	  list<int>::const_iterator ih; map<int, double>::iterator idx;
	  for (ih = t.lHPat().begin(), idx = t.mChi2.begin();
	       ih!=t.lHPat().end(); ih++) {
	    if (*ih<0) continue;
	    THit &h = vecHit[*ih];
	    if (h.Mirror!=-1 && vecHit[h.Mirror].Status!=-4) {
	      double chi2Incr = (*idx).second;
	      if (chi2Incr>chi2IncrMx) {
		THit &hm = vecHit[h.Mirror];
		int ipl = h.IPlane; const TPlane &p = setup.vPlane(ipl);
		const TDetect &d = setup.vDetect(p.IDetRef);
		CsDetector *det = d.PtrDet();
		im = t.mHef.find(ipl); const THlx &Hf = (*im).second;
		double ye = Hf(1), ze = Hf(2), ue; bool error;
		im = t.mHeb.find(ipl); const THlx &Hb = (*im).second;
		ye += Hb(1); ze += Hb(2); ye /= 2; ze /= 2; ue = ye*d.Ca+ze*d.Sa;
		double ump = doCorr ? 
		  det->getCorrU(hm.PtrClus(),ye*10,ze*10,tT,error)/10 : hm.u;
		double du = h.u-ue, dum = ump-ue; 
		if (fabs(dum)<fabs(du)) {
		  badHits.push_back(&h); altHits.push_back(&hm);
		}
	      }
	    }
	    idx++;
	  }
	  if (badHits.size()) {// ********** ALTERNATIVE HITS LIST... **********
	    for (int i = 0; i<(int)badHits.size(); i++) {
	      THit *h = badHits[i], *hm = altHits[i]; int ipl = h->IPlane;
	      if (doCorr) {
		// ***** PROPAGATION(+XRAY+...) CORR. of HIT'S ALTERNATIVE *****
		const TPlane &p = setup.vPlane(ipl);
		const TDetect &d = setup.vDetect(p.IDetRef);
		CsDetector *det = d.PtrDet();
		im = t.mHef.find(ipl);
		const THlx &H = (*im).second; double y = H(1), z = H(2), ump;
		bool error; ump = det->getCorrU(hm->PtrClus(),y*10,z*10,tT,error);
		if (!error) { hm->u = ump/10; hm->Rotate(); }
	      }
	      t.ReplaceHit(hm);
	    }
	    double newChi2 = 0; if ((ret = t.FullKF(1)) &&  // ***** REFIT *****
				    (newChi2 = t.Chi2tot/(t.NDFs-5))<chi2) {
	      ret = t.FullKF(-1); newChi2 = t.Chi2tot/(t.NDFs-5);
	    }
	    if (ret && !newChi2) CsErrLog::msg(elFatal,__FILE__,__LINE__,
  "Event #%d track #%d: KF fit successful and yet chi2 not updated",event,t.Id);
	    if (!ret || chi2<newChi2) {       // ***** NOT BETTER: RESTORE *****
	      for (int i = 0; i<(int)badHits.size(); i++) {
		THit *h = badHits[i]; t.ReplaceHit(h);
	      }
	      ret = t.FullKF(1); if (ret) ret = t.FullKF(-1);
	      if (!ret) {
	CsErrLog::msg(elError,__FILE__,__LINE__,
  "Track #%d: KF fails after original track has been restored => erase.",t.Id);
		listTrack.erase(it++); continue;
	      }
	    }
	    else tStatus |= 0x4;
	  }  // End bad hits found
	}  // End try and cure all bad hits
	if ((TOpt::ReMode[29]&0x4) && (tStatus&0x6)) {    // ...if requested...
	  //                                     ***** ...UPDATE CsCluster *****
	  if (tStatus&0x4) {
	    for (int i = 0; i<(int)badHits.size(); i++) {
	      THit *hm = altHits[i];
#ifdef TrkF_DEBUG_CorrU
	      CsDet *det = hm->PtrClus()->getDetsList().front();
	      if (det->GetTBName()=="ST03X1ub") printf("TRF_CorrU %s %.2f -> %.2f\n",
                det->GetTBName().c_str(),hm->PtrClus()->getU(),hm->U*10);
#endif
	      hm->PtrClus()->setU(hm->U*10);
	    }
	  }
	  if (tStatus&0x2) {
	    THit *hm = altWorst;
#ifdef TrkF_DEBUG_CorrU
	    CsDet *det = hm->PtrClus()->getDetsList().front();
	    if (det->GetTBName()=="ST03X1ub") printf("TRF_CorrU %s %.2f -> %.2f\n",
                det->GetTBName().c_str(),hm->PtrClus()->getU(),hm->U*10);
#endif
	    hm->PtrClus()->setU(hm->U*10);
	  }
	  if (tStatus&0x8) {
	    THit *hm = altWorsF;
#ifdef TrkF_DEBUG_CorrU
	    CsDet *det = hm->PtrClus()->getDetsList().front();
	    if (det->GetTBName()=="ST03X1ub") printf("TF_CorrU %s %.2f -> %.2f\n",
                det->GetTBName().c_str(),hm->PtrClus()->getU(),hm->U*10);
#endif
	    hm->PtrClus()->setU(hm->U*10);
	  }
	}
      }  // End "if cleaning requested"
      if (mChi2) delete mChi2;
      if (!ret) {                                     // ***** IF FAILS... *****
	CsErrLog::mes(elError,"KF fit failed. Track removed.");
	tracksToBeDeleted.push_back(cCst);      // ...CsTrack IS TO BE DELETED
	if (yokeTr) delete yokeTr;
	listTrack.erase(it++); continue;        // ...ERASE TTrack
      }
      if (hdChi2) {
	double chi2 = t.Chi2tot/(t.NDFs-5);
	hdChi2->Fill(chi2-oldChi2,fabs(beamTime));
	if (inVertex) hdChi2InV->Fill(chi2-oldChi2,fabs(beamTime));
	//#define TrkRf_DEBUG
#ifdef TrkRf_DEBUG
	printf("ID=%3d tT=%.1f(%.1f) %.2f %.2f",
	       t.Id,tT,beamTime,chi2,oldChi2);
	if (chi2>oldChi2+.01) {
	  printf("  KO\n");
	}
	else printf("\n");
#endif
      }
      if ((t.Type&0x1) && t.Hfirst(0)>170) {   // ***** TRY microMegas EXTENSION
	// Idea is that some MM hits may be recovered now that DCs are corrected
	// for T0 offset. And these are important because of their impact on the
	// track's parameters @ vertex. E.g: 07W27/cdr32001-58926, Evt 1055129.
	// For the time being, and for simplicty's sake, let's consider only
	// MM01. And let's build-in most of the info about the setup: the 170 cm
	// supra is mid-way between MM01 and MM02, and we look at tracks reco'd
	// no more upstream than this, which in addition have DC hits, and
	// just corrected one, because we are w/in the "if (tStatus)" block.
	int ih0 = t.lHitPat.front(), ipl0 = -1;
	static const TPlane *p; static const TDetect *d;
	if (ih0>=0) {
	  ipl0 = vecHit[ih0].IPlane; p = &setup.vPlane(ipl0);
	  d = &setup.vDetect(p->IDetRef);
	}
	double oldChi2 = t.Chi2tot/(t.NDFs-5);
	if (ipl0>=0 && p->Station->Type==27 && oldChi2<2) {
	  int doErase, ist; for (ist=doErase = 0;
				 ist<(int)setup.vStation().size(); ist++) {
	    const TStation &s = setup.vStation()[ist]; if (s.Type!=27) continue;
	    //        ***** EXTRAPOLATE to MM01 TStation *****
	    int ipl = s.IPlanes[0]; p = &setup.vPlane(ipl);
	    d = &setup.vDetect(p->IDetRef); double x0 = d->X(0);
	    THlx hlx; hlx(0) = x0; t.Hfirst.Extrapolate(hlx,true);
	    double y0 = hlx(1), z0 = hlx(2), yp = hlx(3), zp = hlx(4);
	    if (!d->InActive(y0,z0)) break;
	    double dy2  = hlx(1,1), dz2  = hlx(2,2), dyz  = hlx(1,2);
	    THit *hBests[8]; int jpl, npl = s.IPlanes.size();
	    for (jpl = 0; jpl<npl; jpl++) hBests[jpl] = 0;
	    int nBests; for (jpl=nBests = 0; jpl<npl; jpl++) {
	      //  ***** LOOP on PLANES in MM01: PICKING-UP HITS *****
	      double ca, sa , u; if (jpl) {
		ipl = s.IPlanes[jpl]; p = &setup.vPlane(ipl);
		if (!p->IFlag) continue; d = &setup.vDetect(p->IDetRef);
		ca = d->Ca; sa = d->Sa;
		double x = d->X(0); u = (y0+yp*(x-x0))*ca + (z0+zp*(x-x0))*sa;
	      }
	      else {
		if (!p->IFlag) continue;
		ca = d->Ca; sa = d->Sa; u = y0*ca+z0*sa;
	      }
	      double du = sqrt(dy2*ca*ca+2*dyz*ca*sa+dz2*sa*sa);
	      double cut = 4*du+3*d->Resol, umn = u-1.5*cut, umx = u+1.5*cut;
	      THit *hbest; int mult; double best;
	      vector<int>::const_iterator ihit;
	      for (ihit = p->vHitRef().begin(), mult = 0, hbest = 0, best = cut;
		   ihit!=p->vHitRef().end(); ihit++) {  // Loop over hits
		THit &h = vecHit[*ihit];
		if (!h.sTrackID().empty())// No need to use "THit::Status" here,
		  // since we have no mirror to worry about
		  continue;
		double U = h.U; if (U<umn) continue; if (U>umx) break;
		mult++; double diff = fabs(U-u);
		if (diff<best) { best = diff; hbest = &h; }
	      }
	      if (hbest && (mult==1 || best<cut/2)) {
		hBests[jpl] = hbest; nBests++;
	      }
	    }
	    if (nBests>=3) {  // ***** REQUIRE AT LEAST 3 HITS *****
	      THlx oldFirst(t.Hfirst), oldLast(t.Hlast);
	      for (jpl = 0; jpl<npl; jpl++) if (hBests[jpl]) {
		t.AddHit(*hBests[jpl]); hBests[jpl]->status = 1;
	      }
	      if (TOpt::ReMode[21]==0) t.InsertMissedPlanes();
	      // One could loosen the variance of 1/p here: don't know if this
	      // would help. On the example quoted supra, it didn't change much.
	      THlx hSave = t.Hfirst;
	      bool ret = t.QuickKF(-1,1); if (ret) ret = t.QuickKF(1,1);
	      if (ret) ret = t.FullKF(1); if (ret) ret = t.FullKF(-1);
	      if (!ret || t.Chi2tot/(t.NDFs-5)>oldChi2*1.5) {
		// ***** REQUIRE NO WORSENING OF CHi2: ELSE UNDO CHANGES *****
		t.Behead(ipl0-1); t.Hfirst = oldFirst; t.Hlast = oldLast;
		t.Chi2tot = oldChi2*(t.NDFs-5);
		for (jpl = 0; jpl<npl; jpl++) if (hBests[jpl])
		  hBests[jpl]->status = 0;
		t.Hfirst = hSave;
		ret = t.FullKF(1); if (ret) ret = t.FullKF(-1);
		if (!ret) {
		  CsErrLog::msg(elError,__FILE__,__LINE__,
  "Track #%d: KF fails after original track has been restored => erase.",t.Id);
		  doErase = 1; break;
		}
	      }
	      else t.UseHitTime();
	    }
	    break; // MM01 encountered
	  } // End loop on TStation's
	  if (doErase) { listTrack.erase(it++); continue; }
	} // End case of track starting @ MM plane
      }

      //#define TrkF_IMPROVE_FIT
#ifdef TrkF_IMPROVE_FIT    // ***** HIGH MOMENTUM in SM1 ONLY: IMPROVE FIT *****
      if ((t.Type&0x7)==0x3 && fabs(t.Hlast(5))<1/60.) {
	t.Hfirst(5,5) = 1.E-4;
	double chi2 = t.Chi2tot/t.NDFs; int iter = 0; while (iter<5) {
	  int ichi2_prev = int(chi2*50);
	  if (!(ret = t.FullKF(1))) break; if (!(ret = t.FullKF(-1))) break; 
	  chi2 = t.Chi2tot/t.NDFs;
	  if (int(chi2*50)== ichi2_prev) break; // converged
	}
	if (!ret) {
	  tracksToBeDeleted.push_back(cCst);      // ...CsTrack IS TO BE DELETED
	  if (yokeTr) delete yokeTr;
	  listTrack.erase(it++); continue;
	}
      }
#endif

      if (yokeTr) {      //    ***** YOKE TRACK *****
	//  Let's append the downstream piece back.
	//  It has not been corrected for the diff. between event and trigger
	// time: too bad, but the purpose of having the downstream piece
	// attached to the track is only that it will ease (compared to a track
	// described by its upstream piece on the one hand and its downstream
	// piece that would have no momentum on the other hand) the task of
	// associating calo clusters to it, so as not to leave orphan clusters
	// that could then be associated to neutrals. Given the coarse
	// granularity of the calos, one can get along w/ the bias due the
	// incorrect T0.
	AddYokeTrTail(t,yokeTr);
	delete yokeTr;
      }

      //                                              ***** UPDATE HELICES *****
      const vector<CsHelix> &helices = cCst->getHelices();
      vector<CsHelix>::const_iterator icsh = helices.begin();
      const_cast<CsHelix&>(*icsh)= t.Hfirst.ExportHelix();
      if (helices.size()-2!=t.lHsm.size())
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
     "TTrack Id=%d has %d smoothed helices, CsTrack counterpart (Id=%d) has %d",
		      t.Id,t.lHsm.size(),cCst->getId(),helices.size()-2);
      list<THlx>::iterator ihsm;
      for (icsh++, ihsm = t.lHsm.begin(); ihsm!=t.lHsm.end(); icsh++, ihsm++) {
	THlx &hsm = *ihsm; double x = hsm(0);	// Loop over "smoothed" helices
	//                               ***** SPECIAL CASE OF SM2 YOKE TRACK...
	// ...Do not update smoothed helices beyond the end point of the
	// upstream piece.
	// (Note: The condition infra (on the distance of the end point w.r.t.
	// smoothing point) is not foolproof! Because the distance between 2
	// successive smoothing points can be arbitrarily small if Monitoring
	// histos are booked (cf. "DstProdMon.cc" in $PHAST/coral) and the
	// upstream piece ends in the upstreammost plane of e.g. a GEM, let's
	// say GMXV. The next smoothing point will then be in GMXU, which Z
	// abscissa is usually set 50 um dowstream in our detectors.dat's (hence
	// the 0.001 margin infra) but need not be so, can be arbitrarily close
	// to GMXV and even, since revision r13067 of "./src/geom/CsGeom.cc",
	// exactly at the same Z as GMXV. It looks like "DstProdMon" tries to
	// have a smoothing point for each and every detector plane, even if
	// they are not more distant apart than 0.000001 (hard to undertstand
	// why: a unique smoothing point per GEM would make no difference => one
	// can only guess it's so for sake of simplifying the source code (and
	// btw it may be that"DstProdMon" does not work properly w/ a setup that
	// takes advantage of the revision r13067 mentioned supra)).
	// => Therefore, one would need to figure another solution to work this
	//   out (other than cutting a priori on the distance of end point
	//   w.r.t. to smoothing point), other solution that may involve changes
	//   also in "DstProdMon"
	//  Fortunately, the case (where the upstream piece of a yoke track ends
	// in a GEM) is very unlikely (as the small # of crashes that can be
	// traced to it testifies).
	// => Therefore, there's no urgency...)
	if (t.IMark>0 && x>(*yokeTrHlast)(0)+0.001) break;
	if (t.GetSmoothed(hsm,x)<0)
	  CsErrLog::mes(elError,
	    "TTrack: Cannot get smooth helix upon refit");
// 	  CsErrLog::msg(elFatal,__FILE__,__LINE__,
// 	    "TTrack Id=%d: Cannot get smooth helix @ x=%.3f upon refit",t.Id,x);
	const_cast<CsHelix&>(*icsh) = hsm.ExportHelix();
      }
      const_cast<CsHelix&>(*icsh) = t.Hlast.ExportHelix();
      if (isMC) t.FindKine();
    }
    it++;
  }  // End loop on TTrack's
  return modified;
}
