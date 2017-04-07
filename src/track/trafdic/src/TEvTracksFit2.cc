// $Id: TEvTracksFit2.cc 14094 2015-11-06 15:28:48Z lsilva $

#include "TH1.h"
#include "TH2.h"
//#define TrkF_DEBUG_DECOALESCE
#ifdef TrkF_DEBUG_DECOALESCE
#  include "CsHistograms.h"
#endif

#include "CsErrLog.h"
#include "CsRandom.h"
#include "CsInit.h"
#include "DaqDataDecoding/DaqOption.h"
#include "CsMCTrkHit.h"
#include "CsDetector.h"
#include "CsEvent.h"
#include "CsEventUtils.h"
#include "CsCluster.h"
#include "TOpt.h"
#include "TSetup.h"
#include "TEv.h"
#include "TLattice.h"

using namespace std;

/*! 
  \brief Fit list of tracks
*/

/*
  Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/TracksFit":
   i) Straight fit of single segment track is moved to lattice specific
     "TEvFitSegments".
  ii) Single-zone track: KF fit only if:
      - QN fit was, successfully, performed.
      - Momentum below cut (this in order to ensure that the
       its determination from fringe field is reliable).
 iii) NDF for the calculation of the reduced chi2 derived from TTrack::NDFs
     instead of TTrack::NHits.
  iv) Track segment downstream of muon wall stripped away if that allows to
     restore good chi2.
   v) Attempt at cleaning away ghost hits for tracks w/ bad chi2.
  vi) Correction for propagation time (as well as trigger jitter and time offset
     of pile-up tracks) in drifts detectors. Together w/ decoalescence of the
     hits. And attempt at eliminating ghost mirrors.
 vii) Foretracking to RICHwall, muWallA and hodos.
viii) Possibility of refit with restricted subset of detectors (Det2Go2Fit).
  ix) Beam:
     - MC, muon setup case, w/ BMS: simulation of BMS by taking the momentum
      from thruth bank and smearing it.
     - In other cases: a momentum taken from options is assigned to the beam
      track (this until, or by default of, BMS reco).
   x) SI disambiguation.
  xi) Tracks extending beyond E+HCAL2+MF2 muFilter (&0x8):
     - Cut on overall chi2, chi2 increment.
     - Cure edge effect w/ alternative KF.
*/

void TEv::TracksFit2()
{

  // ******************** INITIALISATIONS ********************

  const TSetup &setup = TSetup::Ref();

  const unsigned int allTrigs = 0xffff; // Although TCS can handle 0x7fffff, it was decided to limit max. #triggers to 16. And take advantage of the fact to make use of the other bits of the trigger TDC module.
  unsigned int evtTrig = trig_mask&allTrigs; // Cut away trailing bits
  const unsigned int caloTrig = 0x10;
  static unsigned int imloTrigs = 0, imloPat = 0;  // Set infra, if needed
  static unsigned int caloBasedTrigs = caloTrig, highQ2Trig = 0, inclOTrigs = 0;
  static int trigStations[4][2], iplHMn = -1;
  static bool first = true; if (first) {
    first = false;
    if (TOpt::ReMode[33]) {  // ***** FORETRACKING to HODOs (except HGs)
      // Determine the 2 TStation's of the hodo system for each IMLO trigger.
      char IMLO[] = "IMLO"; int imlo; imloTrigs = 0x10f;
      for (imlo = 0; imlo<4; imlo++) {
	trigStations[imlo][0]=trigStations[imlo][1] = -1;
      }
      for (int ist = 0; ist<(int)setup.vStation().size(); ist++) {
	const TStation &s = setup.vStation()[ist];
	int ipl = s.IPlanes[0]; const TDetect &d = setup.vDetect(ipl);
	if (d.Name[0]!='H') continue;
	imlo = -1; for (int i = 0; i<4; i++) if (d.Name[1]==IMLO[i]) imlo = i;
	if (imlo<0) continue;
	int i, i01; for (i = 0, i01 = -1; i<2; i++)
	  if (trigStations[imlo][i]==-1) { i01 = i; break; }
	if (i01==-1) CsErrLog::msg(elFatal,__FILE__,__LINE__,
  "ForeTracking to Hodos requested (TraF ReMode[33]), but # of H%c stations > 2",
				   IMLO[imlo]);
	trigStations[imlo][i01] = ist; if (iplHMn==-1) iplHMn = ipl;
	int isOFF; for (i = 0, isOFF = -1; i<(int)s.IPlanes.size(); i++)
	  if (!setup.vPlane(s.IPlanes[i]).IFlag) isOFF = s.IPlanes[i];
	if (isOFF>=0 && (TOpt::ReMode[33]&0x2))
	  // Note: one may not be able to avoid turning off hodos (e.g. HO03
	  // during 2009 DVCS) => let's not make it a fatal error.
	  CsErrLog::msg(elError,__FILE__,__LINE__,
  "ForeTracking to Hodos requested (TraF ReMode[33]), but hodo plane \"%s\" is OFF",
				    setup.vDetect(isOFF).Name.c_str());
      }
      for (imlo = 0; imlo<4; imlo++) {
	if (trigStations[imlo][0]==-1 || trigStations[imlo][1]==-1)
	    CsErrLog::msg(elError,__FILE__,__LINE__,
  "ForeTracking to Hodos requested (TraF ReMode[33]), but # of H%c stations < 2",
			  IMLO[imlo]);
	else imloPat |= 1<<imlo;
      }
    }
    const CS::DaqOption &opts = CsInit::Instance()->getDaqEventsManager().GetDaqOptions();
    int bit = -1;// Bit 0x200 doesn't necessarily mean high-Q2 trigger: check it
    try { bit = opts.GetTriggerBit("LargeQ2Trigger"); } catch ( ... ) { }
    if (bit>=0) { highQ2Trig = 1<<bit; caloBasedTrigs |= highQ2Trig; }
    bit = -1;    // There may be an Incl.OT, at times,...
    try { bit = opts.GetTriggerBit("InclOuterTrigger"); } catch ( ... ) { }
    if (bit>=0)   inclOTrigs = 1<<bit;
  }
  bool searchMAMu = TOpt::ReMode[32]==2 && // Search for mu in MAs: conditioned
    //  The reliability of "ForeTrack2MA" has not been evaluated yet (as of
    // 05/08). But the only way that this can upset the reco in a definitive way
    // is when a, fake, mu candidate is found in MAs (and tagged as muID'd by
    // "muIDinMW1") while the true mu is not muID'd, if and only if this gets
    // the true mu not to be included in the pVertex. All other cases where
    // fake mu's are reco'd (and later ID'd), whether they lead to 2 mu/mu'
    // vertices or a single, fake, one, are recoverable: would suffice to
    // discard the MA mu's from the analysis.
    //  But there's also a question of CPU: we do not want to waste CPU in MAs
    // when the proba for finding anything there is low => therefore, when
    // ForeTrack2MA is conditional, we consider only 2 cases:
    //  - Pure combination of calo and highQ2 triggers.
    //  - C+Hodo-trigger w/ no reco'd track accounting for the fired hodos. This
    //   can be due to bad reco, but also be that hodo-trig was pulled by a
    //   decay track or a combination of accidentally coincident tracks. Because
    //   of the latter cases, and provided that C-trigger is present (meaning an
    //   interaction did take place), we stand a chance to find a mu in MAs.
    ((evtTrig&caloBasedTrigs)==evtTrig ||
     (evtTrig&imloTrigs));  // (re-evaluated infra)
  int searchHodoMu = TOpt::ReMode[33];

  //    ***** In view of UPDATING DRIFTS THits w/ EVENT TIME ***** 
  //       ***** REQUIRE THAT TRIGGER MATCHES ReMode[36]... *****
  unsigned int inclMask =// Trigger pattern for which bits we require only an...
    TOpt::ReMode[36]&(~TOpt::ReMode[40]);// ...inclusive match, cf. UpdateDrifts
  bool evtTimeRequired = ((evtTrig&TOpt::ReMode[36])==evtTrig ||
			  (evtTrig&inclMask)) &&
    // The following can happen in MC: event did not pass the trigger...
    evtTrig || // => Exclude the case. No possible trigger offset anyway.
    // ...OR THAT ALTERNATIVE MASTER TRIGGER is defined (note that then the
    // event timing is automatically granted, being defined by the time of the
    // alternative trigger, cf. "TEv::SetEventTime")
    ptrEvt()->getAlterMasterTrigger()>=0 ||
    // ...OR 2ND ITERATION of TRACKING GOING ON
    ptrEvt()->reTrackingON();
  double evtT0; if (!evtTimeRequired || (eventTime==0 && eventTRef==0)) {
    // Null "eventTime" means it was not possible to define (or it turned out
    // to be defined strictly identical to zero... but never mind). We will then
    // refrain from decoalescing drift hits.
    evtT0 = 0;         evtTConsidered = false;
  }
  else {
    // If drift CsCluster's have already been updated (by "UpdateDrifts"),
    // "eventTime" has been set = 0 (while the actual value has been backed-up
    // into "eventTRef"). Since "evtT" is used infra in the sole updating drift
    // hits, the, possibly null, "eventTime" is indeed the value to consider.
    evtT0 = eventTime; evtTConsidered = true;
  }


  //        ******************** BEAM TRACKS ********************
  // I.e., in fact, beam telescope
  //  - Either a standard beam track, i.e. restricted to the 0x10 zone of the
  //   beam telescope alone,
  //  - Or a beam track bridged over the target (i.e. 0x11 type). Whether it's
  //   over a target dipole ("MagFlag1[0]>2") or not does not make much
  //   difference: the dipole strength is too low for a precise determination
  //   of a beam particle's momentum (expected to be of the order of 100 GeV).
  BeamsFit2();

  bool ret;
  list<TTrack>::iterator it = listTrack.begin(); while (it!=listTrack.end()) {

    TTrack &t = *it;
    t.UseHitTime();
    if (t.Type==0x10 || t.Type==0x11) {  // Beams: already processed supra.
      it++; continue;
    }

    // ************************* PREPROCESSING(S) *************************
    if (t.NGroups()==0) CsErrLog::mes(elFatal,"TTrack with NGroups == 0!");
    int dofit = 1, fixedP =0;
    if (t.IMark>0 && !(t.Type&0x2)) {
      //                 ***** DOWMSTREAM SEGMENT OF YOKE TRACK *****
      t.Hfirst(5,5) = 1.E-20; fixedP = 1;      
    }
    else if (t.IMark<0) //   ***** TRACK PIECE LEFT OVER FROM BRIDGING *****
      // - One expects to catch here the single-zone segments already merged
      //  w/ (appended to) another track and awaiting the validation of the
      //  merging for their fate to be decided.
      // - "BridgeSegments2" leaves over such track segments in the case of
      //  the bridging over the (downstream) muon absorber (and only in that
      //  case). The reason being that there, the simplistic "TTrack::QuickKF"
      //  method is used which cannot do a decent job in solving the difficult
      //  task of bridging over the thickness of the absorber.
      dofit = 0;
    else if (t.NGroups()==1) {
      // ********** SINGLE SEGMENT, NON-BEAM, NON SM2-YOKE, TRACK **********
      // - Either magnet-OFF data (I), fringe-field track (II) or else
      //  fragments of reconstruction (III).
      // - In all cases, assigning momentum allows to account for multi-
      //  scattering in the fit, and, technically, enables smoothing points,
      //  i.e. best track estimators, that are only computed, at any given
      //  point, by the FullKF.
      //  I) FullKF and later on smoothing can be obtained via option
      //   "ReMode[24]". The momentum can only be the beam line tune, which is
      //   expected to be specified by "iCut[15], dCut[4] and [5]".
      // II) The bending through fringe-field allows to determine momentum,
      //   albeit w/ limited precision. The starting value for the momentum is
      //   retrieved, if availbale, from the dico fit. The whole thing is
      //   assumed to be of little interest in case (I), and therefore ignored
      //   when "ReMode[24]".
      // III) No action explicitly intended. When "ReMode[24]", the fragments
      //   get assigned the momentum of the beam line tune.
      if (TOpt::ReMode[24]==0 || // If smoothing required or...
	  TOpt::ReMode[23]) {    // ...if explictly requested by option...
	t.Hfirst(5) = TOpt::iCut[15]/TOpt::dCut[4];
	t.Hfirst(5,5) = 1.E-20; fixedP = 1;// ...variance ~= 0, to fix momentum
      }
#ifdef FRINGE_FIELD_KF
	// - Require a finite starting value that cannot but have been assigned
	//  by "TTrack::IFit&0x8", cf. "BridgeSegments2".
      else if ((t.IFit&0x8) &&                       // If QN fit and...
	       (TOpt::dCut[66]<0 ||
		fabs(1/t.Haux(5))<TOpt::dCut[66]) && // ...low enough momentum
	       t.InGroup(0)) {                       // ...and in zone 0x1
	t.Hfirst(5) = t.Haux(5);        // => KF full fit
	t.Hfirst(5,5) = 1.E-4;          // ...w/ large initial momemtum error
      }
#endif
      else {
	t.Hfirst(5) = 0.; dofit = 0;   // ...Else nothing...and reset P (to be safe)
      }
    }
    else {  // ********** SOME SPECIAL CASES of MULTI SEGMENT TRACKS **********

      if (t.Type==0xc) {                    // Mu-wall momentum-less
	t.Hfirst(5) = TOpt::iCut[15]/TOpt::dCut[4];
	t.Hfirst(5,5) = 1.E-20; fixedP = 1; // ...variance ~= 0, to fix momentum
      }
      else if (!t.Hfirst.with_mom())
	t.Hfirst = t.Haux;                  // QN only so far
    }

    if (!dofit) { it++; continue; }

    //     ************************************************************
    //     *************************  FITTING *************************
    //     ************************************************************
    // (Note: Beam tracks are fitted elswhere, cf. supra.)

    //   I) Full KF Fit.
    //  II) Cleaning
    // III) Signal propagation and other corrections/cleaning in drifts.
    //  IV) Is ForeTracking for muID required ("searchMAMu/searchHodoMu")?

    //           ***** PREPARATION STEPS *****
    if (TOpt::ReMode[21]==0) t.InsertMissedPlanes();
    if (t.Hfirst(5,5)==0) t.Hfirst(5,5) = 1.E-4;
    TTrack *muFT0x8 = 0;                       // ***** BRIDGED OVER mu-WALL?...
    static int muFKine; if (TOpt::ReMode[1] && t.Type&0x8) {
      // ...Remember partial chi2 (before bridging) and associated 0x8 segment.
      int ia = t.Associate; if (ia<0)            // No associated segment...
	// ...Should never happen, cf. "TEv::BridgeSegments2".
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
  "Inconsistency: Event #%d, Track #%d, extending past muFilter: w/o associate",
		      event,t.Id);
      list<TTrack>::reverse_iterator jt;
      for (jt = listTrack.rbegin(); jt!=listTrack.rend(); jt++) {
	TTrack &tj = *jt; if ((int)tj.Id==ia) { muFT0x8 = &tj; break; }
      }
      if (!muFT0x8)           // TTrack w/ ID = associate's does not exist...
	// ...Should never happen, cf. "TEv::BridgeSegments2".
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
  "Inconsistency: Event #%d, 0x%x-Track #%d's 0x8-associate(#%d) doesn't exist",
		      event,t.Type,t.Id,ia);      
      muFKine = t.IKine;
    }

    //  **** CORRECTING MICROMEGAS CLUSTER POSITION ****
    if (TOpt::ReMode[47]) CorrectMMClusters(t);

    //   ********** FIT BACKWARD/FOREWARD w/ MATERIAL INTO ACCOUNT **********

    double foreChi2 = 0;
    if (!(ret = t.FullKF(1)) ||              // ***** FORWARD  FIT HAS FAILED...
	//                                 ***** ...OR ANOMALOUSLY LARGE CHI2...
	//  This should never happen... But it does.
	// E.g. 07W30/cdr29012-59891#50337789 (because of incompatible 0x3 and
	// 0x6 bridges, which is not yet checked in "TEv::BridgeSegments2") or
	// 07W30/cdr27029-59891#121671929 (bad bridge over muFilter: full chi2
	// is then not yet computed).
	//  N.B.:
	//  - The min. chi2 considered as nomalous is arbitrary. May have to be
	//   tuned?...
	//  - Bad bridges over wall are also considered further down .
	//  In order to be on the safe side, erase such case, except...
	(foreChi2 = t.Chi2tot/(t.NDFs-5))>200) {
      if (muFT0x8) {           // ...if the track is bridged over wall...
	// ...then SHORTHEN TTrack down to last THit upstream of wall
	t.Shorten(setup.vIplFirst()[3]); t.Type &= 0x7;
	CsErrLog::msg(elWarning,__FILE__,__LINE__, // "elError" temporarily...
  "Track #%d of type 0x%x: failing to extend past muFilter: KF fit in error",
		      t.Id,t.Type);
	t.Associate = -1;  // Reference to associated 0x8 track no longer useful
	muFT0x8->IMark = 0;// No longer flag 0x8 piece for erasure
	if (t.Type==0x4 && // I.e. a 0xc track just stripped of its 0x8 tail...
	    t.IMark<=0) {  // ...and not a downstream yoke track segment
	  ret = t.QuickKF(-1,0); if (ret) ret = t.QuickKF(1,0);
	}
	ret = t.FullKF(1);
      }
      if (!ret || (foreChi2 = t.Chi2tot/(t.NDFs-5))>200) {
	CsErrLog::msg(elWarning,__FILE__,__LINE__,
  "Track #%d of type 0x%x: Forward KF failed => Removed.",t.Id,t.Type);
	listTrack.erase(it++); continue;              // ***** ...=> ERASE TRACK
      }
      if (t.Type==0x4 && // I.e. an initially 0xc track stripped of its 0x8 tail
	    t.IMark<=0) {// (yet not a downstream yoke track segment)
	it++; continue;
      }
    }

    map<int,double> mChi2 = t.mChi2; 
    double foreRadLenFr = t.radLenFr;

    ret = t.FullKF(-1); double backChi2;
    if (ret && (backChi2 = t.Chi2tot/(t.NDFs-5))>.8*TOpt::dCut[16] &&
	backChi2>2.5*foreChi2) {
      // Case where fore/backward fits do not agree (while large). Since the
      // case is extraordinary, let's advertise it as a "BasicInfo" => will
      // always be displayed, whatever the configuration of the error logger.
      CsErrLog::msg(elBasicInfo,__FILE__,__LINE__,
  "Event #%d Track #%d(0x%x): KF For/backward: chi2,%%X0 = %.2f,%.1f/%.2f,%.1f",
	event,t.Id,t.Type,foreChi2,foreRadLenFr,backChi2,t.radLenFr);
      if (muFT0x8 && foreRadLenFr>t.radLenFr+2.5)
	// Possible explanation is edge effect: track is travelling on the edge
	// of a piece of absorber, which is seen differently in the fore and
	// backward directions. The case spans a small region of phase space,
	// but this is a region that is difficult to describe in MC. Therefore
	// it's interesting to be as little dependent as possible on the edge
	// effect and try to recover the case.
	ret = t.KFwExtraMS(foreRadLenFr-t.radLenFr);
    }

    //    ********** CLEANING (scifi-tracks/other tracks) TRACKS **********
    // (Do this now that chi2 increments have been set by the BACKWARD KF, which
    // enhances (or so I assume) problems in the upstream, and hence scifi, part
    // of the TTrack's.)

    if (ret && TOpt::ReMode[44]) { // ***** SI DISAMBIGUATION  *****
      // (This is meant for hadron setup, w/ SIs in zone 0x1.)
      // - Correct SI assignment is of utmost importance: otherwise track will
      //  not not pass vertexing.
      // => Let's re-evaluate the SI hits association, w/ "BackTrack2SI".
      // - But let's restrict ourselves to favourable cases: tracks which
      //  already possess SI hits(*) in zone 0x1 and have momentum, hence also
      //  a 0x2 segment.
      // - N.B.: "BackTrack2SI" is bugged: it does not always handle tracks
      //  bridged over the target properly => reject tracks w/ a 0x10 segment.
      // - (*): It would also be useful to re-evaluate the tracks w/o SIs (and
      //  see whether we could assign SIs to them, again because of vertexing
      //  considerations). But this would be better done after all tracks w/
      //  SIs have been checked and corrected (after the current loop on tracks,
      //  thus), using the sole SI hits remaining then free. (Note: nothing
      //  along this line is yet written...)
      if ((t.Type&0x13)==0x3 &&
	  setup.vDetect(t.lPRef().front()).IType==21) { // SIs must be in front
	if (BackTrack2SI(t)) { // Track modified by "BackTrackSI"
	  ret = t.FullKF(1);
	  if (ret) { mChi2 = t.mChi2; ret = t.FullKF(-1); }
	  if (!ret) {    // "BackTrackSI" returning error => Erase track
	    listTrack.erase(it++); continue;
	  }
	}
      }
    }

    bool checkVtxD =    // Check vertex detectors = SIs
      // Tracks are closely packed in the vertex detectors zone => higher
      // probability of mistakes. In order to correct these, better would be to
      // have a dedicated piece of software. Not yet done. In the mean time...
      t.NGroups()>1 && (t.Type&0x1) && TOpt::dCut[76]>0;
    bool scifi0x1Mup =  // Candidate mu' w/ scifis in zone 0x1
      // A ghost hit in the scifis of zone 0x1 can prevent a track from being
      // associated to the pVertex. This is of particular importance if the
      // track turns out to be a candidate mu'. In turn, the latter often ends
      // up in 0x1 scifis, and lends itself easily to chi2-based cleaning
      // because of it high momentum and hence small multiscattering. Example:
      // evt #8 of D*.2006.03    => Let's carefully check the case
      (t.Scifi&0x1) && t.Type==0xf /* i.e. mu ID (in muWall#2) */;

    if (t.NGroups()>1 && t.Type!=0xc && ret &&
	(scifi0x1Mup && t.Chi2tot/(t.NDFs-5)>1.5 || // Case of candidate mu'
	((t.Scifi&0xff) ||               // Case scifi-(almost)only track or...
	 t.Scifi==0x100 && checkVtxD) && // ..or vertex detection
	TOpt::dCut[13] && t.Chi2tot/(t.NDFs-5)>TOpt::dCut[13] ||
	//                                  Standard case
	TOpt::dCut[16] && t.Chi2tot/(t.NDFs-5)>.8*TOpt::dCut[16])) {
      // (Initially, cut was "dCut{16]" (set =5). Lowering it (<4.5) allows to
      // recover a K0 in "hpss/../03/P1E/DST_mergedDump/evtdump1-30194.raw".
      // But be careful, so more so as "dCut[16]" is also used elsewhere...)

      static unsigned int xuSIs = 0, yvSIs = 0, xyuFIs = 0;
      static bool first = true; if (first) {  // ***** INIT *****
	for (int ipl = (int)setup.vIplFirst()[0]; ipl<=(int)setup.vIplLast()[0];
	     ipl++) {
	  const TPlane  &p = setup.vPlane(ipl);
	  const TDetect &d = setup.vDetect(p.IDetRef); if (d.IType==21) {
	    if (d.Ca<.9) xuSIs |= 1<<ipl;
	    else         yvSIs |= 1<<ipl;
	  }
	  if (d.IType==22) xyuFIs |= 1<<ipl;
	}
      }
      int doErase = 0; for (int bf = 0; bf<2; bf++) {
	int oldNDFs = t.NDFs-5; float oldChi2 = t.Chi2tot/oldNDFs;
	map<int,double> &mC = bf ? mChi2 : t.mChi2;
	int worstPl = -1, worstGr = -1; double worstChi2Incr;
	t.WorstHit(mC,worstPl,worstGr,worstChi2Incr);
	// Is the worst hit registered in a zone where track is almost only in
	// VSAT (cf. TTrack::Scifi) region
	bool inScifi = 0x1<<worstGr&t.Scifi;
	//                 ***** WORST IS BADLY WORST? *****
	bool badTrack = scifi0x1Mup &&  // Case: mu' in 0x1 scifis
	  worstChi2Incr>.25*t.Chi2tot && (1<<worstPl&xyuFIs);
	badTrack |= (t.Scifi&0xff) &&   // Case: worst is scifi
	  worstChi2Incr>TOpt::dCut[14] &&
	  inScifi;// Require worst to be taking place in the SCIFI ONLY SEGMENT
	badTrack |= !inScifi &&         // Case: !scifi
	  (worstChi2Incr>.75*t.Chi2tot || // Strict
	   worstChi2Incr>.50*t.Chi2tot && // Looser but require bad chi2
	   t.Chi2tot/(t.NDFs-5)>TOpt::dCut[16]);
	if (badTrack && t.Type==0x6 && t.NDFs<=8 &&
	    setup.vDetect(worstPl).Name.find("FI05")==0) {
	  // Stupid track: bad chi2, few hits and bound to loose one of its FI05
	  // hits if it's to reach a more reasonable chi2, and hence bound to be
	  // left largely unconstrained in one dimension. That kind of tracks,
	  // typically based on GPP, keep creeping up, despite all efforts at
	  // screening them out at an earlier stage.
	  doErase = 1; break;	  // => Let's erase it at last...
	}
	if (badTrack && t.Type==0x3 && setup.InMuWallA(worstPl)) {
	  // If worst is a MA or HG02 hit, let's erase all hits from the muWallA
	  // sub-zone: the track will stand a second chance to get extended down
	  // there in "TEv::ForeTrack2MA".
	  // (In fact, as of 2013/11, any PA02 hits here erased are definitely
	  // lost: "ForeTrack2MA" does not take thme into account.)
	  t.Shorten(setup.FirstMuWallA);
	  t.Clip(true); // Erase trailing plane references
	  ret = t.FullKF(1);                            // ...REFIT
	  if (ret) ret = t.FullKF(-1);
	  if (ret) continue;
	  else { doErase = 1; break; }
	}
	if (badTrack) {
	  int jhit = t.CancelHit_Clip(worstPl);         // => ERASE WORST HIT...
	  t.Hfirst(5,5) = fixedP? 1.E-20: 1.E-4;
	  ret = t.FullKF(1);                            // ...REFIT
	  if (ret) {
	    foreChi2 = t.Chi2tot/(t.NDFs-1); foreRadLenFr = t.radLenFr;
	    mChi2 = t.mChi2;               // Save chi2 incr. of forward fit
	    ret = t.FullKF(-1);
	    if (ret && (backChi2 = t.Chi2tot/(t.NDFs-5))>.8*TOpt::dCut[16] &&
		// Case where fore/backward fits don't agree. Cf. comment supra.
		backChi2>2.5*foreChi2 && muFT0x8 && foreRadLenFr>t.radLenFr+2.5)
	      ret = t.KFwExtraMS(foreRadLenFr-t.radLenFr);
	  }
	  if (ret) {
	    // Refit is successful. Check whether
	    // i) We havent't removed essential positional info: track
	    //   parameters would then still be constrained.
	    bool unconstrained = false; // Are deemed to mean unconstrained...
	    if (t.Hfirst(1,1)>.0225 || t.Hfirst(2,2)>.0225 ||     // ...1.5 mm
		t.Hfirst(3,3)>25.e-6 || t.Hfirst(4,4)>25.e-6) {   // ...5 mrd
	      unconstrained = true;
	    }
	    // ii) We have indeed gained in chi2/NDF.
	    if (unconstrained || t.Chi2tot/(t.NDFs-5)>oldChi2*.9 ||
		t.Scifi && t.Chi2tot/(t.NDFs-5)>oldChi2*.75) {
	      //                      // ...NO ACTUAL or DUBIOUS GAIN => RESTORE
	      t.RestoreHit(worstPl,jhit);
	      ret = t.FullKF(1);
	      if (ret) {
		foreChi2 = t.Chi2tot/(t.NDFs-1); foreRadLenFr = t.radLenFr;
		mChi2 = t.mChi2; // Restore chi2 incr. of forward fit
		ret = t.FullKF(-1);
		if (ret &&
		    (backChi2 = t.Chi2tot/(t.NDFs-5))>.8*TOpt::dCut[16] &&
		    // Case fore/backward fits don't agree. Cf. comment supra.
		    backChi2>2.5*foreChi2 && muFT0x8 &&
		    foreRadLenFr>t.radLenFr+2.5)
		  ret = t.KFwExtraMS(foreRadLenFr-t.radLenFr);
	      }
	      if (!ret) // If track exits "FullKF" w/ an error (hence may be
		    // left in an inconsistent state)...
		break;  // ...break the for/backward loop (track is then lost..)
	    }
	    else {
	      if (checkVtxD && // Checking vertex detectors..
		  (1<<worstPl&(xuSIs|yvSIs))!=0) {
		// ***** ...LOOK FOR BAD CHI2 IN (ALMOST-)SAME ORIENTATION SIs
		unsigned int sis = (1<<worstPl&xuSIs) ? xuSIs : yvSIs;
		int ihit; map<int,double>::iterator idx;
		for (ihit = 0, idx = mC.begin(), worstChi2Incr = 0;
		     ihit<(int)t.NHits; ihit++, idx++) {
		  int ipl = (*idx).first, bpl = 1<<ipl;
		  if (!(bpl&(xuSIs|yvSIs))) break; if (!(bpl&sis)) continue;
		  double chi2Incr = (*idx).second; if (chi2Incr>worstChi2Incr) {
		    worstPl = ipl; worstChi2Incr = chi2Incr;
		  }
		}
		if (worstChi2Incr>TOpt::dCut[76]) {
		  t.CancelHit_Clip(worstPl);  // => ERASE BAD CHI2 HIT...
		  ret = t.FullKF(1);          // ...REFIT
		  if (ret) {
		    mChi2 = t.mChi2; // Save chi2 incr. of forward fit
		    ret = t.FullKF(-1);
		  }
		  if (!ret) // If track exits "FullKF" w/ an error...
		    break;  // ...break the for/backward loop.
		}
	      }
	      break;
	    }  // End successful cleaning
	  }  // End refit OK
	}  // End bad track
      }  // End loop on back/foreward
      if (doErase) { listTrack.erase(it++); continue; }
      if (ret && t.Type==0x3 &&
	  TOpt::dCut[16] && t.Chi2tot/(t.NDFs-5)>.8*TOpt::dCut[16]) {
	//                ***** After 1st attempt at cleaning, STILL BAD CHI2...
	//                              ***** ...try now ISOLATED GM04 REJECTION
	// Isolated hits are in any case supicious. If, in addition, right after
	// RICH and its long base-line, all the more so. Detectors in this zone
	// are already watched after, cf. "TEv::PrePattern2". Except GM04.
	list<int>::reverse_iterator rh = t.lHitPat.rbegin(); if (*rh>=0) {
	  THit &h = vecHit[*rh]; int iplGM04 = h.IPlane;
	  const TPlane &p = setup.vPlane(iplGM04);
	  const TDetect &d = setup.vDetect(p.IDetRef);
	  bool doCleanGM04 = false; if (d.IType==26 && d.Name[3]=='4') {
	    do { rh++; } while (rh!=t.lHitPat.rend() && *rh<0); if (*rh>=0) {
	      const TPlane &pp = setup.vPlane(vecHit[*rh].IPlane);
	      const TDetect &dp = setup.vDetect(pp.IDetRef); if (dp.X(0)<750) {
		// There's indeed an isolated GM04 hit. Let's check that the
		// track can stand (yet another?) loss of one of its hits.
		unsigned int nHitsMn = TOpt::iPRpar[15];
		if      (t.Scifi&0x1)   nHitsMn += 5;
		else if (t.Scifi&0x100) nHitsMn += 7;
		else                    nHitsMn += TOpt::iPRpar[5];
		doCleanGM04 = t.NDFs>=nHitsMn;
	      }
	    }
	  }
	  if (doCleanGM04) {      // Do try cleaning away isolated GM04
	    t.Shorten(iplGM04); ret = t.FullKF(1); if (ret) {
	      mChi2 = t.mChi2; // Save chi2 incr. of forward fit
	      ret = t.FullKF(-1);
	    }
	  }
	}
      }
    }
    if (ret && muFT0x8) {      //    ***** TRACK BRIDGED OVER muFILTER *****
	double sChi20x8 = 0;// Chi2 increment (due to 0x8 segment) per 0x8 hit
	// (Note: the foreward increments are used, for they put more weight
	// on the most downstream hits.)
	int nHits0x8; map<int, double>::const_iterator idx;
	for (idx = mChi2.begin(), nHits0x8 = 0; idx!=mChi2.end(); idx++) {
	  if ((*idx).first<setup.vIplFirst()[3]) continue;
	  sChi20x8 += (*idx).second; nHits0x8++;
	}
	sChi20x8 /= nHits0x8;
	//#define TrkF_DEBUG_muF
#ifdef TrkF_DEBUG_muF
      if (TOpt::Hist[5]) {// Histogramming (or printing) for tuning muFilter cut
	// Since they're expected to seldom be brought into play, we condition
	// them by both a cpp condition and an option.
	static TH2D *h2 = 0; if (!h2)
	  h2 = new TH2D("hmuFChi2I","#muFilter: #chi^{2} increment per 0x8-hit",
			128,0,16,7,-.5,6.5);
	int fake = 6; if (muFKine>=0) {
	  if (muFT0x8->IKine!=muFKine) fake = muFT0x8->IKine>=0 ? 1 : 2;
	  else                         fake = 0;
	  if (vecKine[muFKine].isPileup()) fake += 3; 
	}
	h2->Fill(sChi20x8,fake);
	if (fake==1)
	  printf("TrkF_DEBUG_muF: Track#%d(0x%x,chi2=%.2f) fake%d: Incr %.f/n",
		 t.Id,t.Type,t.Chi2tot/(t.NDFs-5),fake,sChi20x8);
      }
#endif
      if (!ret ||                                // If, fit failed or...
	  TOpt::dCut[17]!=0 &&
	  t.Chi2tot/(t.NDFs-5)>TOpt::dCut[17] || // ...too large chi2/NDF or...
	  TOpt::dCut[86]!=0 &&
	  sChi20x8>TOpt::dCut[86]) {             // ...too large chi2-increment
	//                                          => SHORTHEN TTrack
	t.Shorten(setup.vIplFirst()[3]); t.Type &= 0x7;
	if (ret)
	  CsErrLog::msg(elError,__FILE__,__LINE__, // "elError" temporarily...
 "Evt #%d Trk #%d(0x%x) failing to bridge muFilter: chi2=%.2f + %.2f*#0x8-hits",
			event,t.Id,t.Type,t.Chi2tot/(t.NDFs-5),sChi20x8);
        else
	  CsErrLog::msg(elError,__FILE__,__LINE__, // "elError" temporarily...
 "Evt #%d Trk #%d(0x%x) failing to bridge muFilter: no KF fitting",
			event,t.Id,t.Type);
	if (t.Type==0x4 && // I.e. a 0xc track just stripped of its 0x8 tail...
	    t.IMark<=0) {  // ...and not a downstream yoke track segment
	  ret = t.QuickKF(-1,0); if (ret) ret = t.QuickKF(1,0);
	}
	else {
	  ret = t.FullKF(1);                       // => and REFIT
	  if (ret) {
	    mChi2 = t.mChi2; // Save chi2 incr. of forward fit
	    ret = t.FullKF(-1);
	  }
	}
	muFT0x8->IMark = -2;        // => and NO ERASING 0x8 PIECE (yet, tag it)
      }
      t.Associate = -1; // Reference to associated 0x8 track no longer useful.
    }
    if (!ret) { // Backward fit has failed. Erase track. Next one
      CsErrLog::mes(elWarning,"Backward Kalman fit failed. Track removed.");
      listTrack.erase(it++); continue;
    }
    if (t.Type==0x4 &&// I.e. an initially 0xc track stripped of its 0x8 tail...
	// ...as a result of the block "TRACK BRIDGED OVER muFILTER"....
	t.IMark<=0) { // ...and not a downstream yoke track segment
      it++; continue;
    }

    if ((TOpt::ReMode[29]&0x2) || TOpt::ReMode[38]) {

      //     *************** CORRECTION/CLEANING in DRIFTS ***************

      // TStatus: -1: No drift encountered  => No init (yet),
      //           0: Drift was encountered => When "doCorr": Init done,
      //           1:                       => When "doCorr": Track modified 
      int tStatus = -1; static bool first = true;
      double XfZ0 = setup.iPlane2Detect(setup.vIplLast() [0]).X(0);
      double XiZ1 = 0, XfZ1 = 0;
      if (setup.vIplFirst().size()>=2) {
	XiZ1 = setup.iPlane2Detect(setup.vIplFirst()[1]).X(0);
	XfZ1 = setup.iPlane2Detect(setup.vIplLast() [1]).X(0);
      }

      // ********** PROPAGATION TIME (+XRAY) and LATENCY CORRECTION **********
      bool doCorr = TOpt::ReMode[29]&0x2;
      bool doUpdateCs = (TOpt::ReMode[29]&0x4) &&
	// Do not convey the corrections to the CsCluster counterparts of the
	// THits if a redefinition of the event's timing is required
	// (cf. comments on the "evtTimeRequired" supra) and not granted.
	(!evtTimeRequired || evtTConsidered);
      bool deCoalesce =     // ********** DRIFT DECOALESCENCE **********
	// Do not decoalesce if a redefinition of the event's timing is required
	// (cf. comments on the "evtTimeRequired" supra) and not granted.
	(!evtTimeRequired || evtTConsidered) && TOpt::ReMode[39];
      double tT0 = evtT0;
      double evtT = eventTRef ? eventTRef : eventTime;

      // ***** 2 STEPS: I) CORRECTION (and possible deletion)
      //                  + LR ambiguity raise coalesced clusters   => REFIT
      //               II) If bad chi2, try worst hit's alternative => 2ND REFIT
      bool hitsPerZone;
      if (t.mHef.size()!=t.lHitPat.size() || t.mHeb.size()!=t.lHitPat.size())
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "TTrack %d inconsistent mHef/b(=%d/%d) != lHitPat(=%d)",
	  t.Id,t.mHef.size(),t.mHeb.size(),t.lHitPat.size());
      list<int>::iterator ih; map<int,THlx>::const_iterator im;
      for (ih = t.lHitPat.begin(), hitsPerZone = false; ih!=t.lHitPat.end();
	   ih++) {
	if (*ih<0) continue;
	THit &h = vecHit[*ih];
	int ipl = h.IPlane; const TPlane &p = setup.vPlane(ipl);
	const TDetect &d = setup.vDetect(p.IDetRef);
	CsDetector *det = d.PtrDet(); if (!det->hasDrift()) continue;
	if (tStatus<0) {                     // First detector encountered...
	  tStatus = 0;
	  //     ***** ACCIDENTALLY COINCIDENT? *****
	  if (t.DispTime>0 && // This requires 2 time-measuring detectors...
	      // ...eliminating the risk to have track's time set by ghost hit.
	      fabs(t.MeanTime-evtT)/t.SigmaTime>2 &&        // Off time and...
	      // Note: When time undef: SigmaTime<0
	      t.Chi2tot/(t.NDFs-5)>5) {                     // ...bad chi2
	    if (fabs(t.Hlast(3))<.05 && fabs(t.Hlast(4))<.05 || // Paraxial...
		// (Paraxial is supposed to single out halo muons)
		fabs(t.MeanTime-evtT)/t.SigmaTime>5)        // ...or vastly off
	      // The latter condition very strict: we are not so much interested
	      // in off-time track anyway... (Note: non-paraxial off-time tracks
	      // are only expected for a hadron, as opposed to muon, beam.)
	      tT0 = t.MeanTime-eventTRef;    // ...Set track's latency "tT0"
	  }
	}
	// **********  I) CASE COALESCED HIT/MIRROR **********
	// - Try decoalesce (require good chi2)
	// - In any case, apply correction, even if coalesced, because of X-Ray.
	// ********** II) STANDARD CASE **********
	//  => Possibly, coalesce hit/mirror that have then gone closer.
	//  => Erase hits that have gone out of time gate ("getCorrU" returns
	//    then an error)
	// Note: In the MC case the THit "h" is by construction the genuine one
	//      of the pair, cf. "TEv::ImportCluster". If it is left uncoalesced
	//      or maintained as a result of the decoalescence, it will be
	//      considered a good hit when it comes to evaluate the track's
	//      genuineness. If it is instead replaced by its mirror (which is
	//      by construction the ghost one of the pair), the tracks gets
	//      a bad mark, as it should.
	if (h.Mirror!=-1 && vecHit[h.Mirror].Status==-4) {
	  //   ***** COALESCED HIT: TRY RAISE LR  if DISTANCE LARGE ENOUGH *****
	  THit &hm = vecHit[h.mirror];
	  CsCluster *c = h.PtrClus(), *cm = hm.PtrClus();
	  double u = c->getU()/10, um = cm->getU()/10;
	  // Which of next/previous TPlane associated to current one? 
	  int jpl = p.Associate ? p.Associate->IPlane : -1;
	  if (jpl==ipl+1) im = t.mHeb.find(ipl);
	  else            im = t.mHef.find(ipl);
	  const THlx &H = (*im).second; double y = H(1), z = H(2);
	  double deltaU =
	    sqrt(H(1,1)*d.Ca*d.Ca+2*H(1,2)*d.Ca*d.Sa+H(2,2)*d.Sa*d.Sa)/d.Resol;
	  if (t.Chi2tot/(t.NDFs-5)<TOpt::dCut[16] && deltaU<1 &&
	      deCoalesce) {
	    if (doCorr) {
	      bool error; double up = det->getCorrU(c,y*10,z*10,tT0,error);
	      if (!error) {
		double ump = det->getCorrU(cm,y*10,z*10,tT0,error);
		if (!error) { u = up/10; um = ump/10; }
	      }
	    }
	    if (fabs(um-u)>2.5*d.Resol) { // Distance > 0.5 * cut
	      double ue = y*d.Ca+z*d.Sa;
	      double du = fabs(u-ue)/d.Resol, dum = fabs(um-ue)/d.Resol;
#ifdef TrkF_DEBUG_DECOALESCE
	      static int jdebug = 0;
	      static int nOKs[3] = {0,0,0}; int ok = 0; char oks[] = "0+-0";
	      static TH2D *hDumVsDu[2];
	      static TH1D *hChi2[2], *hUmu[2], *hDeltaU[2];
	      if (isMC) {
		static bool first = true; if (first) {
		  first = false; char hName[] = "hDumVsDu0";
		  CsHistograms::SetCurrentPath("/Traffic/TracksFit");
		  for (int i = 0; i<2; i++) {
		    sprintf(hName,"hDumVsDu%d",i);
		    hDumVsDu[i] = new TH2D(hName,"dum vs. du",12,0,6,12,0,6);
		    sprintf(hName,"hChi2%d",i);
		    hChi2[i]    = new TH1D(hName,"Chi2/NDF",20,0,10);
		    sprintf(hName,"hUmu%d",i);
		    hUmu[i]     = new TH1D(hName,"(u-um)/resol",12,0,6);
		    sprintf(hName,"hDeltaU%d",i);
		    hDeltaU[i]    = new TH1D(hName,"#deltaU/resol",12,0,6);
		  }
		  CsHistograms::SetCurrentPath("/");
		}
		if      (h.IKine ==t.IKine || t.IKine!=-1 && h.IKin2 ==t.IKine)
		  ok = 1;
		else if (hm.IKine==t.IKine || t.IKine!=-1 && hm.IKin2==t.IKine)
		  ok = 2;
	      }
#endif
	      if (du<2 && dum>3) {
		h.u  = u;  h.sigu = sqrt(c->getCov()(1,1))/10.; 
		h.Rotate();  hm.status = 1; tStatus = 1;
		if (doCorr && doUpdateCs) {
		  //#define TrkF_DEBUG_CorrU
#ifdef TrkF_DEBUG_CorrU
		  if (d.Name=="ST03X1ub") printf("TF_CorrU %s %.2f -> %.2f\n",
		      d.Name.c_str(),c->getU(),u*10);
#endif
		  c->setU(u*10);   // Update CsCluster: cf. explanations infra
		}
#ifdef TrkF_DEBUG_DECOALESCE
		if (isMC) {
		  if (ok) {
		    nOKs[ok-1]++; hChi2[ok-1]->Fill(t.Chi2tot/(t.NDFs-5));
		    hDumVsDu[ok-1]->Fill(du,dum); hDeltaU[ok-1]->Fill(deltaU);
		    hUmu[ok-1]->Fill(fabs(um-u)/d.Resol);
		  }
		  else nOKs[2]++;
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
		if (doCorr && doUpdateCs) {
#ifdef TrkF_DEBUG_CorrU
		  if (d.Name=="ST03X1ub") printf("TF_CorrU %s %.2f -> %.2f\n",
		    d.Name.c_str(),cm->getU(),um*10);
#endif
		  cm->setU(um*10);
		}
#ifdef TrkF_DEBUG_DECOALESCE
		if (isMC) {
		  if (ok) {
		    nOKs[2-ok]++; hChi2[2-ok]->Fill(t.Chi2tot/(t.NDFs-5));
		    hDumVsDu[2-ok]->Fill(du,dum); hDeltaU[2-ok]->Fill(deltaU);
		    hUmu[2-ok]->Fill(fabs(um-u)/d.Resol);
		  }
		  else nOKs[2]++;
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
	    bool error; double up = det->getCorrU(c,y*10,z*10,tT0,error);
	    if (!error) {
	      double ump = det->getCorrU(cm,y*10,z*10,tT0,error); if (!error) {
		u = (up+ump)/20; if (fabs(h.U-u)>d.Resol/10) {
		  // ...check still for X-ray (or the like) correction
		  h.u = u; h.sigu = sqrt(c->getCov()(1,1))/10+fabs(ump-up)/20;
		  h.Rotate(); tStatus = 1;
		  if (doUpdateCs) { // Update CsCluster
#ifdef TrkF_DEBUG_CorrU
		    if (d.Name=="ST03X1ub") printf("TF_CorrU %s %.2f -> %.2f\n",
		      d.Name.c_str(),c->getU(),h.U*10);
#endif
		    if (h.iHit>hm.iHit) { c->setU(u+.001); cm->setU(u-.001); }
		    else                { c->setU(u-.001); cm->setU(u+.001); }
		    c->setSigmaU(h.SigU*10); cm->setSigmaU(h.SigU*10);
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
	  double up = det->getCorrU(c,y*10,z*10,tT0,error); if (!error) {
	    bool coalesced = false; if (h.Mirror!=-1 && tT0>1) {
	      THit &hm = vecHit[h.Mirror];
	      double ump = det->getCorrU(hm.PtrClus(),y*10,z*10,tT0,error);
	      if (fabs(ump-up)<30*d.Resol) { // ***** CORR'ED Us CLOSE ENOUGH...
		//                                ...COALESCE HIT + MIRROR *****
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
		if (doUpdateCs) { // Update Csluster. N.B.: "c" may have
#ifdef TrkF_DEBUG_CorrU
		  if (d.Name=="ST03X1ub") printf("TF_CorrU %s %.2f -> %.2f\n",
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
	      h.u = up/10; h.Rotate(); tStatus = 1;   // ***** UPDATE THit *****
	      if (doUpdateCs) {
		//                     ***** UPON OPTION: UPDATE CsCluster *****
		// - This is useful in 2 cases, when one intends:
		//  - To make use of CsClusters outside of TraFDic for, e.g.,
		//   the evaluation of residuals in coral.
		//  - To monitor the results of TraFDic while TraFFiC modularity
		//   ("ReMode[1]") is set !=2, i.e. this->TEv is destructed upon
		//   exiting "CsTrafficFitting", and, together w/ it, the
		//   corrections brought to THit's.
		// - Updating "CsCluster::U" does not modify the time info of
		//  the cluster. And "CsDrift..::getCorrU" is relying mainly on
		//  this quantity to compute the corrected U. Therefore, even if
		//  CsCluster is reset here, one will still be able to correctly
		//  compute further corrections made to U, e.g. when it comes to
		//  account for the event time that will be made available after
		//  the vertexing, while refitting tracks in "TEv::TracksRefit".
		//  There is a little problem yet: "getCorrU" makes also use of
		//  "U". This to some limited extent: only to determine whether
		//  the CsCluster under exam has to be corrected with a >0 or a
		//  <0 sign, by comparing it with its mirror. This is already
		//  carelessly done: if the correction has the cluster move
		//  close to the wire at the cell's center, one should, in all
		//  rigor, coalesce it w/ its mirror, whereas now the correction
		//  can have the cluster and its mirror exchange their ordering.
		//  Which would have next correction be applied w/ an opposite
		//  sign. Nonetheless, we decide to leave this careless handling
		//  as is. The problem is complex anyway and the stake may not
		//  be worth the large effort required to solve it rigorously.
#ifdef TrkF_DEBUG_CorrU
		if (d.Name=="ST03X1ub") printf("TF_CorrU %s %.2f -> %.2f\n",
		  d.Name.c_str(),c->getU(),up);
#endif
		c->setU(up);
	      }
	    }
	  }
	  else if (fabs(tT0)>5) {// "getCorrU" in error while latency is large...
	    // ...try and erase coresponding THit
	    static int nHits[5]; if (!hitsPerZone) { 
	      hitsPerZone = true; list<int>::iterator ihp;
	      for (ihp = t.lHitPat.begin(), memset((void*)nHits,0,sizeof(nHits)); ihp!=t.lHitPat.end(); ihp++) {
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
	      nHits[igr]--; tStatus = 1; t.CancelHit(ipl);
	    }
	  }
	}  // End correcting propag. (and the like): (de)coalescence, deletion
      }  // End loop on hits (skipping hits from non-drift detectors)

      if (tStatus>0) {                    // ***** PROPAGATION CORR: REFIT *****
      //#define TrkF_DEBUG_PROPAG
#ifdef TF_DEBUG_PROPAG
	printf("DEBUG_PROPAG %3d %5.1f  %8.3f %8.3f %5.1f %c\n",
 t.Id,t.Chi2tot/(t.NDFs-5),t.Hlast(3),t.Hlast(4),t.MeanTime,tT0!=0 ? 'c' : ' ');
#endif
	t.Clip();
	ret = t.FullKF(1);
	if (ret) {
	  foreChi2 = t.Chi2tot/(t.NDFs-1); foreRadLenFr = t.radLenFr;
	  mChi2 = t.mChi2; // Save chi2 incr. of forward fit
	  ret = t.FullKF(-1);
	  if (ret && (backChi2 = t.Chi2tot/(t.NDFs-5))>.8*TOpt::dCut[16] &&
	      // Case where fore/backward fits don't agree. Cf. comment supra.
	      backChi2>2.5*foreChi2 && muFT0x8 && foreRadLenFr>t.radLenFr+2.5)
	    ret = t.KFwExtraMS(foreRadLenFr-t.radLenFr);
	}
      }
      double chi2 = t.Chi2tot/(t.NDFs-5);
      if (tStatus>=0 && TOpt::ReMode[38] && t.NGroups()>1 && ret &&
	  TOpt::dCut[16] && chi2>.8*TOpt::dCut[16]) {

	//            ***** CLEANING of TRACKS w/ DRIFTS HITS *****

	vector<THit*> badHits, altHits;             // ***** BAD chi2 INCREMENTS
	THit *worstHit = 0, *altWorst = 0;         // ***** WORST chi2 INCREMENT
	THit *worstHiF = 0, *altWorsF = 0;
	TSpacePt *sptDC = 0;                          // ***** WORST DC TStation

	if (fabs(t.Hfirst(5))>.4) {// ***** LOW P: DEAL 1ST w/ WORST HIT/STATION
	  // - Determine worst (w/ worst chi2 increment) hits (fore/backward).
	  // - Determine worst (in chi2 increment or # of, DC, hits) TStation.
	  const TStation *prvS, *worstS, *leastS;  // Worst TStation
	  static double sWorstS, sIncr;
	  static int nLeastS, pLeastS, nHitsS, nProjsS;
	  static unsigned int projs;
	  double worstChi2Incr, worstChi2IncF = 0; // Worst hit (back/foreward)
	  list<int>::const_iterator ih; map<int, double>::iterator idx;
	  for (ih = t.lHPat().begin(), idx = t.mChi2.begin(), nLeastS = 8,
		 worstS=leastS=prvS = 0, worstChi2Incr = .25*t.Chi2tot;
	       ih!=t.lHPat().end(); ih++) {
	    // First, determine worst hit and worst TStation
	    if (*ih<0) continue; THit &h = vecHit[*ih];
	    double chi2Incr = (*idx).second; idx++;
	    if (chi2Incr>worstChi2Incr) {
	      if (h.Mirror!=-1 && vecHit[h.Mirror].Status!=-4) {
		worstHit = &h; altWorst = &(vecHit[h.Mirror]);
		worstChi2Incr = chi2Incr;
	      }
	    }
	    int ipl = h.IPlane; const TPlane &p = setup.vPlane(ipl);
	    const TStation *s = p.Station;
	    if (s==prvS) {
	      sIncr += chi2Incr; nHitsS++; unsigned int proj = 1<<p.IProj;
	      if (!(proj&projs)) { projs |= proj; nProjsS++; }
	      if (sIncr>.75*t.Chi2tot && (!worstS || sIncr>sWorstS)) {
		sWorstS = sIncr; worstS = s;
	      }
	    }
	    else {
	      if (prvS && prvS->Type==15 && nHitsS<nLeastS) {
		leastS = prvS; nLeastS = nHitsS; pLeastS = nProjsS;
	      }
	      sIncr = chi2Incr; nHitsS = 1; projs = 1<<p.IProj; nProjsS = 1;
	      prvS = s;
	    }
	  }
	  for (ih = t.lHPat().begin(), idx = mChi2.begin(),
		 worstChi2IncF = .25*t.Chi2tot; ih!=t.lHPat().end(); ih++) {
	    // Now determine the worst foreward hit.
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
	  if ((!worstS || worstChi2Incr>.5*t.Chi2tot) && // If worst TStation...
	      // ...not real worse one, cf. definition of "worstS", or is so due
	      // to a single worst hit => try rather TStation w/ fewest hits...
	      leastS && (nLeastS<=3 || nLeastS<=4 && pLeastS<=3) ||
	      pLeastS<=2) // ...Try rather TStation w/ only 2 proj. in any case.
	    worstS = leastS;
	  if (worstS && worstS->Type==15 /* ie. is DC */ && (t.Type&0x3)==0x3) {
	    // ...Let's make another try at the PR through DCs, which is very
	    // difficult for the present case of low momentum, because DCs are
	    // either embedded in high fringe field (DCi, i>=1) or (DC0) it
	    // concerns, mainly, a track segment w/ only 2 space points, viz.
	    // DC00 and yet another drift, hence LR-ambigous, drift detector.
	    // =>  Get track clone, strip it of all its hits in worst TStation
	    //   and try picking-up hits again by back-tracking, using QN fit.
	    TTrack tLAS ; const double xRICH = 750;
	    int nHsDC; THit *hsDC[8]; // <8 planes in DC, cf. "TSetup::Init".
	    int npl = worstS->IPlanes.size(), iplf = worstS->IPlanes.front();
	    int ipl1 = setup.vIplFirst()[0];
	    int nPats = (iplf+npl-1)/32+1; unsigned int hPatsDC[nPats];
	    memset(hPatsDC,0,nPats*sizeof(unsigned int));
	    list<int>::iterator ih; int nDsZ1, nDsZ2;
	    const TLattice *lat = TLattice::Ptr();
	    for (ih = t.lHitPat.begin(), nHsDC=nDsZ1=nDsZ2 = 0;
		 ih!=t.lHitPat.end(); ih++) {
	      if (*ih<0) continue; THit &h = vecHit[*ih]; int ipl = h.IPlane;
	      const TDetect &d = setup.vDetect(ipl); if (d.X(0)>750) break;
	      tLAS.AddHit(h);
	      const TPlane  &p = setup.vPlane(ipl);
	      if (p.Station==worstS) {
		hsDC[nHsDC++] = &h;
		int kpl = ipl-ipl1; hPatsDC[kpl/32] |= 1<<kpl%32;
		// Reset "THit::Status", for we will want to possibly recycle
		//  those hits associated to the current, problematic track.
		h.status = 0;
		int mirror; if ((mirror = h.Mirror)!=-1) {
		  THit &hm = vecHit[mirror]; if (hm.Status!=-4) hm.status = 0;
		}
	      }
	      else if (lat->dico_idx[ipl]>0) {
		if (ipl<=setup.vIplLast()[0]) nDsZ1++;
		else                          nDsZ2++;
	      }
	    }
	    static double chi2LAS; bool ok = nDsZ1>4 && nDsZ2>=4;
	    if (ok) {
	      for (int i = 0; i<6; i++) tLAS.Hfirst(i) = t.Hfirst(i);
	      ok = tLAS.QNewtonFit(1,3);
	    }
	    if (ok) {
	      chi2LAS = tLAS.Chi2tot/(tLAS.NDics-5);
	      // Let's now erase all hits from worst DC TStation and refit, w/
	      // fixed momentum, since there is possibly now too few hits to
	      // contrain its value. (Note that QN fit in fixed momentum mode
	      // computes guestimates for all 0x7 zones, whatever the pattern
	      // assigned to the TTrack::Type attribute, here = 0.)
	      for (int ihDC = 0; ihDC<nHsDC; ihDC++) tLAS.SubHit(*hsDC[ihDC]);
	      ok = tLAS.QNewtonFit(1,2); double chi2LA2;
	      if (ok &&
		  (chi2LA2 = tLAS.Chi2tot/(tLAS.NDics-5))>TOpt::dCut[16])
		ok = false;
	    }
	    if (ok) {
	      list<int>::const_iterator ih; const THlx *HDC;// Get helix @ DC...
	      // ...and covariant matrix (We assume it does not change over the
	      // length of the TStation and use the smooth helix at 1st plane.)
	      for (ih = t.lHPat().begin(), HDC = 0; ih!=t.lHPat().end(); ih++) {
		if (*ih<0) continue; THit &h = vecHit[*ih]; int ipl = h.IPlane;
		const TStation *s = setup.vPlane(ipl).Station; if (s==worstS) {
		  map<int,THlx>::const_iterator im = t.mHeb.find(ipl);
		  if (im!=t.mHeb.end()) HDC = &(*im).second;
		}
	      }
	      if (!HDC) CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "TTrack %d inconsistent mHeb: no TPlane #%d",t.Id,worstS->IPlanes.front());
	      double dy2 = (*HDC)(1,1), dz2 = (*HDC)(2,2), dyz = (*HDC)(1,2);
	      vector<TSpacePt> spts; unsigned int hitsPats[nPats];
	      for (int jpl = 0; jpl<npl; jpl++) {
		int ipl = worstS->IPlanes[jpl];
		const TPlane  &p = setup.vPlane(ipl);
		if (p.IProj!=1) continue; // Consider vertical coordinate planes
		if (!p.IFlag/* OFF */) continue;
		memset(hitsPats,0,nPats*sizeof(unsigned int));
		const TDetect &d = setup.vDetect(ipl);
		double cu = d.Ca, su = d.Sa, du;
		du = 4*sqrt(dy2*cu*cu+2*dyz*cu*su+dz2*su*su+d.Resol*d.Resol);
		double u0 = tLAS.vGuests[ipl-ipl1], umn = u0-du, umx = u0+du;
		for (vector<int>::const_iterator ihit = p.vHitRef().begin();
		     ihit!=p.vHitRef().end(); ihit++) {
		  THit &h = vecHit[*ihit]; if (h.Status) continue;
		  float u = h.U; if (u<umn) continue; if (u>umx) break; 
		  // ***** TRY AND BUILD A SPACE POINT ON ALL COMPATIBLE z-HITS
		  BuildSpacePt(tLAS,&h,hitsPats,dy2,dyz,dz2,spts,HDC);
		}
	      }
	      if (spts.size()!=0) {          // ***** DETERMINE BEST SPACE POINT
		int is, nspts = spts.size(), nHMx, ispt0; static double chi2Mn;
		for (is=nHMx = 0, ispt0 = -1; is<nspts; is++) {
		  TSpacePt &spt = spts[is];
		  if (ispt0<0 ||
		      (spt.nHits>nHMx || spt.nHits==nHMx && spt.chi2<chi2Mn)) {
		    nHMx = spt.nHits; chi2Mn = spt.chi2; ispt0 = is;
		  }
		}	
		const TSpacePt &spt0 = spts[ispt0];
		if (spt0.nHits>=4) {
		  for (int i = 0; i<8; i++) if (spt0.hPat&1<<i) {
		    THit *h = spt0.hs[i]; tLAS.AddHit(*h);
		  }
		  if (tLAS.QNewtonFit(1,2) &&
		      tLAS.Chi2tot/(tLAS.NDics-5)<chi2LAS)
		    sptDC = new TSpacePt(spt0);
		}
	      }
	    }
	    if (ok) {     // If the new PR could indeed be achieved...
	      if (sptDC) {  // ...and it yielded a space point
		for (int i = 0; i<8; i++) if (sptDC->hPat&1<<i) {
		  THit *h = sptDC->hs[i]; int ipl = h->IPlane, kpl = ipl-ipl1;
		  if (hPatsDC[kpl/32]&1<<kpl%32) {
		    t.CancelHit(h->IPlane); hPatsDC[kpl/32] ^= 1<<kpl%32;
		  }
		  t.AnnexHit(*h); h->status = 1;
		  if (doCorr) {   // Propagation correction
		    const TDetect &d = setup.vDetect(ipl);
		    CsDetector *det = d.PtrDet();
		    bool error; double up =
		      det->getCorrU(h->PtrClus(),sptDC->y*10,sptDC->z*10,tT0,error);
		    if (!error) h->u = up/10;
		  }
		}
		for (int ihDC = 0; ihDC<nHsDC; ihDC++) {
		  THit *h = hsDC[ihDC]; int kpl = h->IPlane-ipl1;
		  if (hPatsDC[kpl/32]&1<<kpl%32) t.CancelHit(h->IPlane);
		}
	      }
	      else {        // ...and no z-based space point: subtract z-coords
		for (int ihDC = 0; ihDC<nHsDC; ihDC++) {
		  THit *h = hsDC[ihDC];
		  if (setup.vPlane(h->IPlane).IProj==1) t.CancelHit(h->IPlane);
		  else             // "Status" has just been canceled ...
		    h->status = 1; // ... => restore it
		}
	      }
	      t.Clip(); // Erase planes inserted upstream of, new, first hit
	      t.InsertMissedPlanes();
	      ret = t.FullKF(1); if (ret) ret = t.FullKF(-1);
	      if (ret) {
		double newChi2 = t.Chi2tot/(t.NDFs-5);
		if (!ret || chi2<newChi2) {  // ***** NOT BETTER: RESTORE *****
		  if (sptDC) {
		    for (int i = 0; i<8; i++) if (sptDC->hPat&1<<i) {
		      THit *h = sptDC->hs[i];
		      t.CancelHit(h->IPlane); h->status = 0;
		    }
		    for (int ihDC = 0; ihDC<nHsDC; ihDC++) {
		      THit *h = hsDC[ihDC]; t.AnnexHit(*h); h->status = 1;
		    }
		    t.Clip();
		    delete sptDC; sptDC = 0;
		  }
		  else {
		    for (int ihDC = 0; ihDC<nHsDC; ihDC++) {
		      THit *h = hsDC[ihDC];
		      if (setup.vPlane(h->IPlane).IProj==1) {
			t.AnnexHit(*h); h->status = 1;
		      }
		    }
		  }
		  if (TOpt::ReMode[21]==0) t.InsertMissedPlanes();
		  ret = t.FullKF(1); if (ret) ret = t.FullKF(-1);
		}
		else { chi2 = newChi2; tStatus |= 0x10; }
	      }
	      if (!ret) {
		CsErrLog::msg(elError,__FILE__,__LINE__,
 "Track #%d: KF fails after original track has been restored => erase.",t.Id);
		  listTrack.erase(it++);
		if (sptDC) delete sptDC;
		continue;
	      }
	    } // End new PR could be achieved
	  } // End re-evaluating PR in DCs
	  if (!(tStatus&0x10) && worstHit) {
	    // If the original track was left unchanged after the new PR in DCs,
	    // try and exchange ``backward'' worst hit.
	    int ipl = worstHit->IPlane; const TPlane &p = setup.vPlane(ipl);
	    const TDetect &d = setup.vDetect(p.IDetRef); double Xd = d.X(0);
	    CsDetector *det = d.PtrDet();
	    if ((fabs(Xd-XfZ0)<10 || fabs(Xd-XiZ1)<10) &&  // I.e. DC01||DC02
		worstChi2Incr>.15*t.Chi2tot || worstChi2Incr>.20*t.Chi2tot) {
	      //                       ***** WORST IS SIGNIFICANTLY WORST? *****
	      if (doCorr) {
		//           ***** PROPAGATION CORR. of WORST-HIT'S MIRROR *****
		im = t.mHeb.find(ipl);
		const THlx &H = (*im).second; double y = H(1), z = H(2);
		bool error;
		double ump = det->getCorrU(altWorst->PtrClus(),y*10,z*10,tT0,error);
		if (!error) altWorst->u = ump/10;
	      }  // End propagation correction for mirror hit
	      altWorst->Rotate(); t.ReplaceHit(altWorst);
	      double newChi2 = 0; if ((ret = t.FullKF(1)) &&// ***** REFIT *****
				  (newChi2 = t.Chi2tot/(t.NDFs-5))<chi2) {
		ret = t.FullKF(-1); newChi2 = t.Chi2tot/(t.NDFs-5);
	      }
	      if (!ret || chi2<newChi2) {     // ***** NOT BETTER: RESTORE *****
		t.ReplaceHit(worstHit);
		ret = t.FullKF(1); if (ret) ret = t.FullKF(-1);
		if (!ret) {
		  CsErrLog::msg(elError,__FILE__,__LINE__,
 "Track #%d: KF fails after original track has been restored => erase.",t.Id);
		  listTrack.erase(it++); continue;
		}
		chi2 = t.Chi2tot/(t.NDFs-5);
	      }
	      else { tStatus |= 0x2; chi2 = newChi2; }
	    }
	  }  // End candidate "worstHit" found
	  if (!(tStatus&0x12) && worstHiF && worstHiF!=worstHit) {
	    // If the original track was left unchanged after the exchange of
	    // the ``backward'' worst hit, and ``forward'' differs: try it.
	    int ipl = worstHiF->IPlane; const TPlane &p = setup.vPlane(ipl);
	    const TDetect &d = setup.vDetect(p.IDetRef); double Xd = d.X(0);
	    CsDetector *det = d.PtrDet();
	    if ((fabs(Xd-XfZ0)<10 || fabs(Xd-XiZ1)<10) &&  // I.e. DC01||DC02
		worstChi2IncF>.15*t.Chi2tot || worstChi2IncF>.20*t.Chi2tot) {
	      if (doCorr) {
		im = t.mHef.find(ipl);
		const THlx &H = (*im).second; double y = H(1), z = H(2);
		bool error;
		double ump = det->getCorrU(altWorsF->PtrClus(),y*10,z*10,tT0,error);
		if (!error) altWorsF->u = ump/10;
	      }  // End propagation correction for mirror hit
	      altWorsF->Rotate(); t.ReplaceHit(altWorsF);
	      double newChi2 = 0; if ((ret = t.FullKF(1)) &&// ***** REFIT *****
				  (newChi2 = t.Chi2tot/(t.NDFs-5))<chi2) {
		ret = t.FullKF(-1); newChi2 = t.Chi2tot/(t.NDFs-5);
	      }
	      if (!ret || chi2<newChi2) {     // ***** NOT BETTER: RESTORE *****
		t.ReplaceHit(worstHiF);
		ret = t.FullKF(1); if (ret) ret = t.FullKF(-1);
		if (!ret) {
		  CsErrLog::msg(elError,__FILE__,__LINE__,
 "Track #%d: KF fails after original track has been restored => erase.",t.Id);
		  listTrack.erase(it++); continue;
		}
		chi2 = t.Chi2tot/(t.NDFs-5);
	      }
	      else { tStatus |= 0x8; chi2 = newChi2; }
	    }
	  }  // End candidate forward "worstHiF" found
	}  // End try and cure worst hit/station if low momentum 
	if (ret && chi2>.8*TOpt::dCut[16]) {
	  //            ***** chi2 STILL BAD: DEALING NOW  w/ ALL BAD HITS *****
	  double chi2IncrMx = .2*.8*TOpt::dCut[16]*(t.NDFs-5); // 20% chi2 cut
	  //#define TrkF_DEBUG_MIRRORS
#ifdef TrkF_DEBUG_MIRRORS
	  vector<int> ipls; vector<double> dus, dums, chi2Is;
#endif
	  list<int>::const_iterator ih; map<int, double>::iterator idx;
	  for (ih = t.lHPat().begin(), idx = t.mChi2.begin();
	       ih!=t.lHPat().end(); ih++) {
	    if (*ih<0) continue;
	    THit &h = vecHit[*ih];
	    if (h.Mirror!=-1 && vecHit[h.Mirror].Status!=-4) {
	      double chi2Incr = (*idx).second; bool badIncr = chi2Incr>chi2IncrMx;
#ifdef TrkF_DEBUG_MIRRORS
	      bool doResidual = true;
#else
	      bool doResidual = badIncr;
#endif
	      if (doResidual) {
		THit &hm = vecHit[h.Mirror];
		int ipl = h.IPlane; const TPlane &p = setup.vPlane(ipl);
		const TDetect &d = setup.vDetect(p.IDetRef);
		CsDetector *det = d.PtrDet();
		im = t.mHef.find(ipl); const THlx &Hf = (*im).second;
		double ye = Hf(1), ze = Hf(2), ue; bool error;
		im = t.mHeb.find(ipl); const THlx &Hb = (*im).second;
		ye += Hb(1); ze += Hb(2); ye /= 2; ze /= 2; ue = ye*d.Ca+ze*d.Sa;
		double ump = doCorr ? 
		  det->getCorrU(hm.PtrClus(),ye*10,ze*10,tT0,error)/10 : hm.u;
		double du = h.u-ue, dum = ump-ue; 
#ifdef TrkF_DEBUG_MIRRORS
		dus.push_back(du); dums.push_back(dum); ipls.push_back(ipl);
		chi2Is.push_back(chi2Incr/t.Chi2tot*100);
#endif
		if (badIncr && fabs(dum)<fabs(du)) {
		  badHits.push_back(&h); altHits.push_back(&hm);
		}
	      }
	    }
	    idx++;
	  }
#ifdef TrkF_DEBUG_MIRRORS
	  static int kdebug = 0;
	  if (kdebug) {
	    if (kdebug==1) kdebug = 0;	  
	    for (int i = 0; i<(int)ipls.size(); i++) {
	      const TPlane &p = setup.vPlane(ipls[i]);
	      const TDetect &d = setup.vDetect(p.IDetRef);
	      double du = dus[i], dum = dums[i];
	      if (d.Ca>.707) printf("%s: %6.1f %5.2f %6.1f  %6.1f  %4.1f\n",
	   d.Name.c_str(),d.X(0),du,du*d.Ca/d.Resol,dum*d.Ca/d.Resol,chi2Is[i]);
	      else 
printf("%s:                                  %6.1f %5.2f %6.1f  %6.1f  %4.1f\n",
	   d.Name.c_str(),d.X(0),du,du*d.Sa/d.Resol,dum*d.Sa/d.Resol,chi2Is[i]);
	    }
	  }
#endif
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
		bool error; ump = det->getCorrU(hm->PtrClus(),y*10,z*10,tT0,error);
		if (!error) { hm->u = ump/10; hm->Rotate(); }
	      }
	      t.ReplaceHit(hm);
	    }
	    double newChi2 = 0; if ((ret = t.FullKF(1)) &&  // ***** REFIT *****
				(newChi2 = t.Chi2tot/(t.NDFs-5))<chi2) {
	      ret = t.FullKF(-1); newChi2 = t.Chi2tot/(t.NDFs-5);
	      if (badHits.size()==1 &&                   // ***** SINGLE BAD HIT
		  ret && newChi2<chi2 &&              // ***** CORR'D W/ SUCCESS
		  newChi2>TOpt::dCut[16]) {  // ***** AND CHI2 STILL PREFECTIBLE
		THit *ha = 0; static int ipl;
		const TPlane *p = setup.vPlane(badHits[0]->IPlane).Associate;
		if (p) { 
		  ipl = p->IPlane;             // ***** LOOK AT ASSOCIATED PLANE
		  
		  for (ih = t.lHPat().begin(); ih!=t.lHPat().end(); ih++) {
		    if (*ih<0) continue; THit &h = vecHit[*ih];
		    if (h.IPlane==ipl) { ha = &h; break; }
		    else if (h.IPlane>ipl) break;
		  }
		}
		if (ha) { // Track has hit on associated plane
		  const TDetect &d = setup.vDetect(p->IDetRef);
		  double chi2Incr = (*t.mChi2.find(ipl)).second;
		  if (chi2Incr>chi2IncrMx) {
		    int mirror = ha->Mirror;
		    THit *hp = mirror!=-1 && vecHit[mirror].Status!=-4 ?
		      &vecHit[mirror] : ha; // Try and change it for mirror
		    double dis = altHits[0]->U-hp->U; if (fabs(dis)>d.Pitch/4) {
		      int incr = dis>0?+1:-1, ih =  hp->IHit+incr;
		      THit *hs = (0<=ih && ih<(int)vecHit.size())?&vecHit[ih]:0;
		      if (hs && hs->Status==-4) {
			ih += incr;     // Try and change for next to mmirror
			hs = (0<=ih && ih<(int)vecHit.size())?&vecHit[ih]:0;
		      }
		      if (hs && hs->iPlane==ha->iPlane &&
			  fabs(altHits[0]->U-hs->U)<fabs(dis)) hp = hs;
		      else hp = 0;
		    }
		    if (hp && hp!=ha) { // Attractive alternative found
		      if (doCorr) {
			im = t.mHef.find(ipl); const THlx &H = (*im).second;
			double y = H(1), z = H(2), upp; bool error;
			upp = d.PtrDet()->getCorrU(hp->PtrClus(),y*10,z*10,tT0,error);
			if (!error) { hp->u = upp/10; hp->Rotate(); }
			t.ReplaceHit(hp); chi2 = newChi2;
			if ((ret = t.FullKF(1)) &&
			    (newChi2 = t.Chi2tot/(t.NDFs-5))<chi2) {
			  ret = t.FullKF(-1); newChi2 = t.Chi2tot/(t.NDFs-5);
			}
			if (doUpdateCs) {
			  // Update the CsCluster of the initial bad hit, for
			  // it will be overwritten in the bad hits list
			  THit *hm = altHits[0];
#ifdef TrkF_DEBUG_CorrU
			  if (d.Name=="ST03X1ub") printf("TF_CorrU %s %.2f -> %.2f\n",
			    d.Name.c_str(),hm->PtrClus()->getU(),hm->U*10);
#endif
			  hm->PtrClus()->setU(hm->U*10);
			}
			// Overwrite bad hit list, so that newly cleaned THit
			// - Will be restored, if cleaning fails
			// - Will have its CsCluster counterpart updated.
			badHits[0] = ha; altHits[0] = hp;
		      }
		    }
		  }
		}
	      }
	    }
	    if (!ret || chi2<newChi2) {       // ***** NOT BETTER: RESTORE *****
	      for (int i = 0; i<(int)badHits.size(); i++) {
		THit *h = badHits[i]; t.ReplaceHit(h);
	      }
	      ret = t.FullKF(1); if (ret) ret = t.FullKF(-1);
	    }
	    else tStatus |= 0x4;
	    if (!ret) {
	      CsErrLog::msg(elError,__FILE__,__LINE__,
 "Track #%d: KF fails after original track has been restored => erase.",t.Id);
		listTrack.erase(it++);
		if (sptDC) delete sptDC;
	      continue;
	    }
	  }  // End bad hits found
	}  // End try and cure all bad hits
	if (doUpdateCs && (tStatus&0x1e)) {          // ...if requested...
	  //                                     ***** ...UPDATE CsCluster *****
	  if (tStatus&0x4) {
	    for (int i = 0; i<(int)badHits.size(); i++) {
	      THit *hm = altHits[i];
#ifdef TrkF_DEBUG_CorrU
	      CsDet *det = hm->PtrClus()->getDetsList().front();
	      if (det->GetTBName()=="ST03X1ub") printf("TF_CorrU %s %.2f -> %.2f\n",
                det->GetTBName().c_str(),hm->PtrClus()->getU(),hm->U*10);
#endif
	      hm->PtrClus()->setU(hm->U*10);
	    }
	  }
	  if (tStatus&0x2) {
	    THit *hm = altWorst;
#ifdef TrkF_DEBUG_CorrU
	    CsDet *det = hm->PtrClus()->getDetsList().front();
	    if (det->GetTBName()=="ST03X1ub") printf("TF_CorrU %s %.2f -> %.2f\n",
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
	  if ((tStatus&0x10) && sptDC) {
	    for (int i = 0; i<8; i++) if (sptDC->hPat&1<<i) {
	      THit *h = sptDC->hs[i]; h->PtrClus()->setU(h->U*10);
	    }
	  }
	}
	if (sptDC) delete sptDC;
      }  // End try and improve bad chi2 track
      if (!ret) {                    // ***** FIT FAILED. ERASE TRACK. *****
	CsErrLog::mes(elWarning,"KF fit failed. Track removed.");
	listTrack.erase(it++); continue;
      }
      //#define TrkF_DEBUG_BAD_CHI2
#ifdef TrkF_DEBUG_BAD_CHI2
      {
	static int kdebug = 0;
	if (kdebug && t.Chi2tot/(t.NDFs-5)>TOpt::dCut[16]) {
	  t.QNewtonFit(1,2);
	}
      }
#endif
    }        // End correct for propagation in drifts

    if (TOpt::dCut[17]!=0 && t.NDFs>5 &&  // ********** CUT on CHI2 **********
	t.Chi2tot/(t.NDFs-5)>TOpt::dCut[17]) { 
      //     ********** RESCUE AT LEAST PART OF HIGH CHI2 TRACK **********
      // (That part w/ smallest chi2 contribution.
      //  Note:
      // - Rescued tracks can serve several purpose:
      //   - If 0x1: can be fringe-field fitted, and if not, when it comes to
      //     check the event's exclusivity, can be checked to come from the
      //     p-vertex and if indeed, can signal that some momentum is missing
      //     in the reconstruction.
      //   - Else, even when momentumless, can be associated to calo clusters.
      // - The whole thing is not very satisfactory. In particular, the
      //  sum of chi2 contributions, by which the rescue is conditioned, is
      //  very much biased by the overall fit. (This forced me to introduce an
      //  ad hoc, built-in, coefficient of 1.5 in the quality cut.) Would be
      //  better to rescue several parts, betting that the high chi2 originates
      //  in the bridging of 2 incompatible but otherwise sane tracks.)
      //       ***** INIT RESCUE: Evaluate track quality per zone *****
      int igr, iLast, ngr = setup.vIplFirst().size(); if (ngr>5) ngr = 5;
      int nHits[5], nSpacePts[5], nZs[5]; double sChi2[5];
      const TStation *sPrv; static int nProjs; static unsigned int projs;
      for (igr = 0; igr<5; igr++) {
	sChi2[igr] = 0; nHits[igr]=nSpacePts[igr]=nZs[igr] = 0;
      }
      list<int>::const_iterator ih; map<int, double>::iterator idx;
      for (ih = t.lHPat().begin(), idx = t.mChi2.begin(), iLast = -1, sPrv = 0;
	   ih!=t.lHPat().end(); ih++) {
	if (*ih<0) continue; int ipl = vecHit[*ih].IPlane; if (ipl>iLast) {
	  for (igr = 0; igr<ngr; igr++) if (setup.vIplFirst()[igr]<=ipl &&
					    ipl<=setup.vIplLast()[igr]) {
	    iLast = setup.vIplLast()[igr]; break;
	  }
	}
	nHits[igr]++; sChi2[igr] += (*idx).second;
	const TPlane &p = setup.vPlane(ipl); unsigned int proj = 1<<p.IProj;
	if (p.IFlag&0x30) // P-pixel plane
	  proj |= 0x2;    // ...add Z-proj (whether it's of the XY or UV kind)
	const TStation *&s = p.Station; if (s!=sPrv) {
	  nProjs = 1; projs = proj; if (proj&0x2) nZs[igr]++; sPrv = s;
	}
	else if (!(proj&projs)) {
	  nProjs++;  projs |= proj; if (proj&0x2) nZs[igr]++;
	  if (nProjs==2) nSpacePts[igr]++; // Space pt: require 2 distinct proj.
	}
	idx++;
      }
      int tryRescue = 0; if ((t.Type&0x7)==0x7) {
	for (igr = 0; igr<4 /* Zone 0x10 is ignored */; igr++) {
	  // Out of principle, there's no reason to normalise by (#hits-4). But
	  // we want to give a bonus to track w/ a large #hits.
	  if (nHits[igr]>4)    sChi2[igr] /= nHits[igr]-4;
	  else if (nHits[igr]) sChi2[igr] = TOpt::dCut[17];
	}
	if      (sChi2[0]+sChi2[1]<sChi2[1]+sChi2[2] && // Try and rescue 0x3
		 sChi2[0]+sChi2[1]<TOpt::dCut[17]*1.5 &&
		 nHits[0]>=((t.Scifi&0x1)?5:TOpt::iPRpar[ 5]) &&
		 nHits[1]>=((t.Scifi&0x202)?5:TOpt::iPRpar[15])) {
	  t.Shorten(setup.vIplFirst()[2]); tryRescue = 0x303;
	}
	else if (sChi2[1]+sChi2[2]<sChi2[0]+sChi2[1] && // Try and rescue 06...
		 sChi2[1]+sChi2[2]+sChi2[3]<TOpt::dCut[17]*1.5 && // ...or 0xd?
		 nHits[1]>=((t.Scifi&0x202)?5:TOpt::iPRpar[15]) &&
		 nHits[2]>=((t.Scifi&0x404)?4:TOpt::iPRpar[25]) &&
		 (!(t.Type&0x8) || nHits[3]>=TOpt::iPRpar[35])) {
	  t.Behead(setup.vIplLast()[0]); tryRescue = 0x606;
	}
	if (tryRescue) {                                // Full KF fit
	  ret = t.FullKF(1); if (ret) ret = t.FullKF(-1);
	  if (!ret || t.Chi2tot/(t.NDFs-5)>TOpt::dCut[17]) {
	    if (!ret) CsErrLog::msg(elError,__FILE__,__LINE__,
  "Track #%d: KF fails while rescuing 0x%x piece of 0x7 track => erase.",
		      t.Id,tryRescue);
	    listTrack.erase(it++); continue;
	  }
	  else {
	    t.Type = tryRescue; t.IFit = 0x60; t.Scifi &= tryRescue;
	  }
	}
      }
      // Now try to rescue single zone segments. Use "dCut[16]", which is
      // expected to be more severe, for these tracks which lack the reliability
      // brought by bridging. (This is not very satisfactoy, since "dCut[16]" is
      // rather meant for track cleaning, and can be arbitrarily low. Anyway, I
      // give it a factor 2.) And exclude scifi-(almost)only segments, because
      // their chi2 is not enough constrained and hence not so meaningful.
      else {
	for (igr = 0; igr<5; igr++) {
	  // Out of principle, there's no reason to normalise by (#hits-8). But
	  // we want to give a bonus to track w/ a large #hits. And since
	  // scifi-(almost)only are ignored here, one can ask for #hits>8.
	  if (nHits[igr]>8) sChi2[igr] /= nHits[igr]-8;
	  else if (nHits[igr]) sChi2[igr] = TOpt::dCut[17];
	}
	if ((t.Type&0x3)==0x3) {
	  if      (sChi2[0]<sChi2[1] &&         // Try and rescue 0x1 from 0x3
		   sChi2[0]<TOpt::dCut[17]*1.5 && nHits[0]>=TOpt::iPRpar[ 5]) {
	    t.Shorten(setup.vIplFirst()[1]); tryRescue = 0x1;
	  }                                     // Try and rescue 0x2 from 0x3
	  else if (sChi2[1]<TOpt::dCut[17]*1.5 && nHits[1]>=TOpt::iPRpar[15]) {
	    t.Behead(setup.vIplLast()[0]);   tryRescue = 0x2;
	  }
	  if (t.IMark>0) { // 0X3 piece of yoke track:
	    // The association w/ the dowstream piece has to be cancelled:
	    // - If 0x1 is rescued, because the association is in any case no
	    //  longer reliable. (Also, if there's subsequent refit, the case of
	    //  a 0x5 yoke track is not properly handled by "TracksRefit".)
	    // - If 0x2 is rescued, because the fit of the downstream piece of
	    //  the yoke track, w/ a fixed momentum which would then be =0,
	    //  would yield meaningless results.
	    t.IMark = 0;
	    list<TTrack>::iterator jt = listTrack.begin(); TTrack *tr = 0;
	    while (jt!=listTrack.end()) {
	      if ((*jt).IMark==(int)t.Id) { tr = &(*jt); break; }
	      jt++;
	    }
	    if (tr) tr->IMark = 0;
	  }
	}
	else if ((t.Type&0x6)==0x6) {
	  if      (sChi2[1]<sChi2[2] &&         // Try and rescue 0x2 from 0x6
		   sChi2[1]<TOpt::dCut[17]*1.5 && nHits[1]>=TOpt::iPRpar[15]) {
	    t.Shorten(setup.vIplFirst()[2]); tryRescue = 0x2;
	  }                                     // Try and rescue 0x4 from 0x6
	  else if (sChi2[2]<TOpt::dCut[17]*1.5 && nHits[2]>=TOpt::iPRpar[25]) {
	    t.Behead(setup.vIplLast()[1]);   tryRescue = 0x4;
	    if (t.Type&0x8) // Cancel extension into 0x8, since this extension
	      // can only be granted to tracks w/ momentum. (Note: In fact, the
	      // 0x8 extension can be granted to a momentumless track, cf.
	      // "TEv::BridgeSegments2")
	      t.Shorten(setup.vIplFirst()[3]);
	  }
	}
	else if ((t.Type&0x11)==0x11) { // I.e. track bridged over target
	  if      (sChi2[4]<sChi2[0] &&         // Try and rescue 0x10 from 0x11
		   sChi2[4]<TOpt::dCut[17]*1.5 && nHits[4]>=TOpt::iPRpar[45]) {
	    t.Shorten(setup.vIplFirst()[2]); tryRescue = 0x10;
	  }                                     // Try and rescue 0x1 from 0x11
	  else if (sChi2[0]<TOpt::dCut[17]*1.5 && nHits[0]>=TOpt::iPRpar[ 5]) {
	    t.Behead(setup.vIplLast()[1]);   tryRescue = 0x1;
	  }
	}
	if (tryRescue) {
	  ret = t.QuickKF(-1,0); if (ret) ret = t.QuickKF(1,0);
	  if (!ret || t.Chi2tot/(t.NDFs-4)>2*TOpt::dCut[16]) {
	    if (!ret) CsErrLog::msg(elError,__FILE__,__LINE__,
  "Track #%d: KF fails while rescuing 0x%x piece of 0x%x track => erase.",
				    t.Id,tryRescue,t.Type);
	    listTrack.erase(it++); continue;
	  }
	  else {
	    t.Type = tryRescue; t.IFit = 0x1;
	    if (t.Type==0x1 &&  // Case 0x1: try fit w/ P (from fringe field)...
		TOpt::ReMode[14] && // ...require dico fitting
		nSpacePts[0]>=3 &&  // ...(obviously) >=3 space points
		nZs[0]>=2) { // ...>=2 Z points (as Z planes can be scarce)
	      double oldChi2 = t.Chi2tot;// Save total chi2, as QN overwrites it
	      if (t.QNewtonFit(1,1) && t.Chi2aux/(t.NDics-5)<TOpt::dCut[9] &&
		  // Validate momentum derived by dico fit from bending in
		  // fringe field only if momentum low enough, i.e. < "dCut[66]"
		  fabs(1/t.Haux(5))<TOpt::dCut[66]) {
		t.Hfirst(5) = t.Haux(5);  // => KF full fit
		t.Hfirst(5,5) = 1.E-4;    // ...w/ large initial momemtum error
		ret = t.FullKF(1); if (ret) ret = t.FullKF(-1);
		if (!ret || t.Chi2tot/(t.NDFs-5)>TOpt::dCut[17]) {
		  t.QuickKF(-1,0); t.QuickKF(1,0); t.IFit = 0x1;
		}
		else t.IFit = 0x60;
	      }
	      else { // Track's helices weren't changed. But have to restore...
		t.Chi2tot = oldChi2; t.IFit = 0x1; // ... total chi2 and flag.
	      }
	    }
	  }
	}
      }
      if (!tryRescue) {
	listTrack.erase(it++); continue;
      }
    }
    else t.IFit |= 0x60; // Update IFit flag


    //                                    ***** CONDITIONAL ForeTracking *****
    // - "ForeTrack2MA": Determine whether we have a candidate reco'd muon to
    //  play the role of the scattered muon (before setting out looking for
    //  it). => Require
    //   i) Reco downtream of muFilter (and momentum) => 0xe "Type".
    //  ii) In time.
    // iii) Consistent w/ trigger: we only check it has appropriate hodo hits.
    //     And do not check trigger matrix correlation. Any trigger.
    //  (We restrict the analysis to triggers w/ Calo. For the other ones
    //  the proba to have a MA mu in addition to the mu having fired the
    //  hodos of the trigger is low.)
    // - "ForeTrack2Hs": Select candidates. Requiring 0x6 reco, in time,
    //  1st hodo hit consistent w/ current trigger. Cancel search if a
    //  fully consistent combination is found.
    if ((t.Type&0xe)==0xe && searchMAMu && (evtTrig&imloTrigs) ||
	(t.Type&0x6)==0x6 && (TOpt::ReMode[33]&0x2) &&
	// Not an exclusive combination of C and highQ2
	(evtTrig&caloBasedTrigs)!=evtTrig) {
      t.UseHitTime(); // Usefull? Hasn't this been already executed?
      double evtT = eventTRef ? eventTRef : eventTime;
      if (!evtT || t.SigmaTime>0 &&
	  fabs(t.MeanTime-evtT)<3*t.SigmaTime+1/* 3sigma + 1ns evtT resol. */) {
	list<int>::const_iterator ih, ip; int hodoOKs[4] = {0,0,0,0};
	for (ih = t.lHPat().begin(), ip = t.lPRef().begin();
	     ih!=t.lHPat().end(); ih++, ip++) {
	  if (*ih<0 || *ip<iplHMn) continue;
	  int ist = setup.vPlane()[*ip].Station->IStation;
	  for (int imlo = 0; imlo<4; imlo++) {
	    if (!(1<<imlo&imloPat)) continue;
	    for (int i =0; i<2; i++)
	      if (ist==trigStations[imlo][i]) hodoOKs[imlo] |= 1<<i;
	  }
	}

	for (int imlo = 0; imlo<4; imlo++) {
	  bool isInTrig = (1<<imlo&evtTrig) ||
	    imlo==1 && (evtTrig&0x100) ||    // Inclusive middle
	    imlo==3 && (evtTrig&inclOTrigs); // Inclusive outer (and J/psi?)
	  if (hodoOKs[imlo]==0x3) {
	    // We consider all of IMLO, even if they don't belong to the
	    // current trigger. Again because the idea is to restrict
	    // ourselves to cases w/ a high proba for yet another, MA, mu.
	    searchMAMu = false;
	    if (isInTrig) searchHodoMu = 0;
	  }
	  else if (hodoOKs[imlo]==0x1 && isInTrig) t.Scifi |= 0x10000<<imlo;
	}
      }  // End in-time
    }

    it++;

  }  // End of cleaning current track

  //           ***** ERASE TRACK PIECES LEFT OVER from BRIDGING *****
  // Now that mergings have been evaluated => If validated, the left over piece
  // is flagged by IMark==-1.
  // (Note: One expects to encounter here only those pieces left over from the
  // bridging over the downstream muon absorber. Others must have been dealt w/
  // in "TEv::BridgeSegments2".)
  it = listTrack.begin(); while (it!=listTrack.end()) {
    if ((*it).IMark==-1) { listTrack.erase(it++); }
    else it++;
  }

  if (TOpt::ReMode[44])         // ***** SI DISAMBIGUATION  *****
    // The call to "BackTrack2SI" may have left some tracks sharing a same hit.
    // Following the strategy adopted so far throughout TraFDic, we decide to
    // ``clean'' them. Although it could be argued that SIs are a special case,
    // where tracks heading for the vertex are concentrated in a small region,
    // where sharing a common hit may be frequent and cleaning the tracks may
    // leave some of them w/ too few hits, leading to imprecision in the track
    // parameters or even removal of the track. All this can be dealt with in
    // "CleanTrackListSI" anyway.
    CleanTrackListSI();

  //    *************** EXTEND TTrack's by FORETRACKING... ***************
  // Do this now that TTrack's have been FullKF'd, since extended TTrack's
  // will have to be checked for their chi2 increment w.r.t. their original
  // counterpart.
  ForeTrack2RICHWall();         // ...to RICH WALL 
  //                               ...to MAs
  if (TOpt::ReMode[32]==1 ||          // ...in any case or...
      searchMAMu)                     // ...conditioned: and condition fulfilled
    ForeTrack2MA();
  if (TOpt::ReMode[33] && (!(TOpt::ReMode[33]&0x4) || searchHodoMu))
    ForeTrack2Hs();

  it = listTrack.begin(); while (it!=listTrack.end()) {
    TTrack &t = *it;
    if (!(t.IFit&0x60) ||
	t.IMark>0 && !(t.Type&0x2) /* Skip downstream piece of yoke track*/) {
      it++; continue;
    }

    //             *************** FINALISE TTrack: ***************

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
      if (!ret) { listTrack.erase(it++); continue; }
    }
#endif

    //                                                  ***** "SMOOTHED" HELICES
    for (int i = 0; i<int(sizeof(TOpt::SmoothPos)/sizeof(double)); i++) {
      // Loop over "smoothing" points
      double x = TOpt::SmoothPos[i]; if (x==0.) break;
      THlx Hsm; if (t.GetSmoothed(Hsm,x)<0) continue;
      t.lHsm.push_back(Hsm); // store smoothed
    }

    if (TOpt::ReMode[22]) t.Refine(1);     // ***** HITS MAP: PREPARE FOR, HISTO

    TTrack *tr = 0; if (t.IMark>0) {   // ***** 0X3 PIECE OF YOKE TRACK *****
      // Yoke tracks, i.e. tracks going through the yoke of SM2 have been split
      // into 2 pieces, a 0x2 or 0x3 and a 0x4 one (cf. "TEv::BridgeSegments2"),
      // w/ X-references stored in "TTrack::IMark", in order for the downstream
      // piece not to upset the overall fit. Now that the 0x3 piece has been
      // fitted, let's append the downstream piece. Note that the latter has
      // been fitted w/ a fixed momentum = that of the 0x3 piece obtained from
      // the QuickKF performed in "TEv::BridgeSegments2", i.e. neither updated
      // by 0x3's FullKF, nor corrected for the energy loss in the yoke.
      list<TTrack>::iterator jt = listTrack.begin();
      while (jt!=listTrack.end()) {
	if ((*jt).IMark==(int)(*it).Id) { tr = &(*jt); break; }
	jt++;
      }
      if (!tr) t.IMark = 0;
      else {                   // ***** APPEND THE DOWNSTREAM of SM2 TRACK PIECE
	//  Prepare the hits maps of the downstream piece. They will then be
	// encoded in the hits list of the latter, which will then be merged w/
	// the hits list of the upstream piece.
	tr->Refine(1);
	//  Push upstream piece's last helix into the list of smoothed helices
	// so that the last actual measurement of the momentum is remembered.
	//  Before that, check whether there are any trailing empty placeholder
	// in the hits map. If indeed, shorten the track.
	if (t.lHPat().back()<0) {
	  t.Shorten(setup.vIplFirst()[2]);
	  ret = t.FullKF(1); if (ret) t.FullKF(-1); if (!ret) {
	    CsErrLog::msg(elError,__FILE__,__LINE__,
	      "SM2 Yoke track #%d: Kalman fit failed. Track removed.",t.Id);
	    tr->IMark = 0; listTrack.erase(it++); continue;
	  }
	}
	t.lHsm.push_back(t.Hlast);
	// Save characterstics of the 2 pieces prior to merging.
	int iplTlL = t.lPRef().back(), iplTrF = tr->lPRef().front();
	double lastCop = t.Hlast(5), lastDcop = t.Hlast(5,5);
	double radLenFr = t.radLenFr + tr->radLenFr;
	// Only meaningful fit is that of the 0x3 piece: assign its chi2/NDF to
	// the whole track.
	t.Chi2tot = t.Chi2tot*(t.NDFs+tr->NDFs-5)/(t.NDFs-5);
	t.Append(*tr);
	// Update characteristics of downstream piece.
	t.Hlast(5) = lastCop;
	// Precision on q/P for the downstream piece has been set = infinite (so
	// as to keep momentum fixed to that of the 0x3 piece)...
	t.Hlast(5,5) = lastDcop; // ... update it
	// Insert missed planes in the interval between the end of the 0x3
	// piece and the beginning of the downstream piece, and then update the
	// hits map. The later cannot be done using "TTrack::Refine", given that
	// the arrays of helices have not been defined in the interval.
	if (TOpt::ReMode[21]==0) t.InsertMissedPlanes();
	THlx hs[2]; hs[0] = t.lHsm.back();
	list<int>::iterator ih; list<int>::const_iterator ip; int i12;
	for (ih = t.lHitPat.begin(), ip = t.lPRef().begin(), i12 = 0;
	     ih!=t.lHitPat.end(); ih++, ip++) {
	  int ipl = *ip; if (ipl<=iplTlL) continue; if (iplTrF<=ipl) break;
	  if (*ih!=-1) CsErrLog::msg(elError,__FILE__,__LINE__,
  "Inconsistency when merging ID=%d and ID=%d: hit reference on TPlane %d = %d",
				     t.Id,tr->Id,ipl,*ih);
	  const TPlane &p = setup.vPlane(ipl);
	  const TDetect &d = setup.vDetect(p.IDetRef);
	  THlx &hStrt = hs[i12], &hStop = hs[1-i12]; i12 = 1-i12;
	  hStop(0) = d.X(0); hStrt.Extrapolate(hStop,true);
	  if (!d.InActive(hStop(1),hStop(2))) *ih = -2;
	  else if (TOpt::ReMode[20]<=0 || !setup.InMaterialMap(hStop(0))) {
	    double len    = d.Siz(0)/hStop.DirCos(1); // Detector's thickness
	    double radLen = d.RadLen;                 // Detector's rad. length
	    radLenFr += len/radLen;
	  }
	}
	t.radLenFr = radLenFr + hs[i12].RadLenFr();
	// Correlations between q/P and the rest may have drifted away from 0 
	for (int i = 1; i<5; i++) t.Hlast(i,5) = 0; // ...reset them
	for (int i = 0; i<int(sizeof(TOpt::SmoothPos)/sizeof(double)); i++) {
	  // Add the "smoothing" points located downstream of SM2
	  double x = TOpt::SmoothPos[i]; if (x==0.) break;
	  THlx Hsm; if (tr->GetSmoothed(Hsm,x)<0) continue;
	  Hsm(5) = lastCop; Hsm(5,5) = lastDcop;
	  for (int i = 1; i<5; i++) Hsm(5,i) = 0;
	  t.lHsm.push_back(Hsm); // store smoothed
	}
      }
    }

    t.UseHitTime();                                              // ***** TIMING

    //#define Tracksfit_REFIT_APART    // ***** (IF ENABLED) REFIT in SM1/2 ONLY
#ifdef Tracksfit_REFIT_APART
    if (TOpt::ReMode[41]) // Bridging through SM2 yoke is allowed. Then
      // "Tracksfit_REFIT_APART" becomes a real pain => forbid it.
      CsErrLog::mes(elFatal,"\"Tracksfit_REFIT_APART\" while bridging through SM2 yoke allowed (\"TraF ReMode[41]\")");
    if ((t.Type&0x7)==0x7) { // If track spans both SM1 and SM2...
      //                        ...Refit track independently in the 2 magnets
      TTrack t1 = t;         // I) SM1
      //                        ...SHORTHEN all what's downstream of SM2
      static double pinv;
      t1.Shorten(setup.vIplFirst()[2]); t1.Type &= 0x3;
      if ((ret = t1.QuickKF(1,1))) {
	pinv = (t1.Hfirst(3)-t1.Hlast(3))/(.3E-3*setup.MagFieldInt[1]);
	t1.Hfirst(5) = pinv; t1.Hfirst(5,5) = 1.E-4;
	double chi2 = t1.Chi2tot/t1.NDFs; int iter = 0; while (iter<5) {
	  // ***** REFIT (using same algorithm as in TEv::BridgeSegments2) *****
	  int ichi2_prev = int(chi2*10);
	  t1.Hlast(5) = t1.Hfirst(5); t1.Hlast(5,5) = t1.Hfirst(5,5);
	  if (!(ret = t1.QuickKF(-1,1))) break; 
	  chi2 = t1.Chi2tot/t1.NDFs;
	  if (int(chi2*10)== ichi2_prev) break; // converged
	}
      }
      if (ret) ret = t1.FullKF(1);
      if (ret) {
	if ((ret = t1.FullKF(-1))) { 
	  t.lHsm.push_back(t1.Hfirst); t.lHsm.push_back(t1.Hlast);
	  printf("======= Refitting w/ SM1 only: P %.3f -> %.3f -> %.3f\n",
		 1/t.Hfirst(5),1/pinv,1/t1.Hfirst(5));
	}
      }
      if (!ret)
	CsErrLog::mes(elError,"REFIT_APART SM1 KF refit failed!");
      TTrack t2 = t;         // I) SM2
      //                        ...BEHEAD all what's upstream of SM1
      t2.Behead(setup.vIplLast()[0]); t2.Type &= 0x6;
      if ((ret = t2.QuickKF(-1,1))) {
	pinv = (t2.Hfirst(3)-t2.Hlast(3))/(.3E-3*setup.MagFieldInt[2]);
	t2.Hlast(5) = pinv; t2.Hlast(5,5) = 1.E-4;
	double chi2 = t2.Chi2tot/t2.NDFs; int iter = 0; while (iter<5) {
	  // ***** REFIT (using same algorithm as in TEv::BridgeSegments2) *****
	  int ichi2_prev = int(chi2*10);
	  t2.Hfirst(5) = t2.Hlast(5); t2.Hfirst(5,5) = t2.Hlast(5,5);
	  if (!(ret = t2.QuickKF(1,1))) break; 
	  chi2 = t2.Chi2tot/t2.NDFs;
	  if (int(chi2*10)== ichi2_prev) break; // converged
	}
      }
      if (ret) ret = t2.FullKF(-1);
      if (ret) {
	if ((ret = t2.FullKF(1))) { 
	  t.lHsm.push_back(t2.Hfirst); t.lHsm.push_back(t2.Hlast);
	  printf("======= Refitting w/ SM2 only: P %.3f -> %.3f -> %.3f\n",
		 1/t.Hlast(5),1/pinv,1/t2.Hlast(5));
	}
      }
      if (!ret)
	CsErrLog::mes(elError,"REFIT_APART SM2 KF refit failed!");
    }
#endif
    //#define TrackF_DEBUG_KF
#ifdef TrackF_DEBUG_KF
    static int idebug = 0;
    if (idebug) t.FullKF(1);  // => Get ``backward'' chi2 increments...
#endif
    it++;
  }

  //   ********** CLEAN AWAY DOWNSTREAM PIECES of YOKE TRACKS **********
  it = listTrack.begin(); while (it!=listTrack.end()) {
    if ((*it).IMark==-1) { // Downstream piece got "IMark=-1" while appended
      listTrack.erase(it++); continue;
    }
    else it++;
  }

  if (TOpt::ReMode[26]&0x10) {// ********** BEAM CONTINUED BEYOND SM1 **********
    // Check the reliability of the continuation once again, w/ "FullKF", now
    // that the 0x1X piece of it has been cleaned and FullKF'd.
    double nSigmas = TOpt::iCut[30]>0 /* hadron beam */ ? 3 : sqrt(12.);
    double p0 = TOpt::dCut[4]*TOpt::iCut[15], dp0 = nSigmas*TOpt::dCut[5]*p0*p0;
    for (it = listTrack.begin(); it!=listTrack.end(); it++) {
      TTrack &b = *it;
      if (!(b.Type&0x10)) break; // Tracks of interest, i.e., beams, come first
      int ia = b.Associate; if (ia<0) // Association (cf. "BridgeSegments2")...
	continue;                     // ...required
      TTrack *t; list<TTrack>::iterator jt;
      for (jt = listTrack.begin(), t = 0; jt!=listTrack.end(); jt++) {
	TTrack &tj = *jt; if ((int)tj.Id==ia) { t = &tj; break; }
      }
      if (!t) { // TTrack w/ ID = associate's no longer exists...
	// ...It must have been erased at some point => Cancel the association.
	CsErrLog::msg(elError,__FILE__,__LINE__,
  "Beam track #%d's associate(=#%d) is found to have been erased",b.Id,ia);
	b.Associate = -1; continue;
      }
      if (b.Type!=0x10) {
	CsErrLog::msg(elError,__FILE__,__LINE__,
  "Beam track #%d association(<->#%d) cancelled: it's been bridged over target",
		      b.Id,ia);
	b.Associate = -1;
      }
      else {
	short int iMark = t->IMark; // Backup "TTrack::IMark", cf. infra.
	TTrack t0x1X(b); t0x1X.Append(*t); // Build 0x1X track
	t->IMark = iMark; // Preserve the original 0xX track.
	t0x1X.Hlast = t->Hlast;// 0x1X close to 0xX => need not rescale cov.
	t0x1X.FullKF(-1); double chi20x1X = t0x1X.Chi2tot/(t0x1X.NDFs-5);
	bool ok = chi20x1X<3 || // Tight cut, since we expect a high P track...
	  // ...Yet allow for some reco error in the 2 pieces (likely in 0xX)
	  chi20x1X<4 && t0x1X.Chi2tot<1.5*(b.Chi2tot+t->Chi2tot);
	double p = 1/t0x1X.Hfirst(5), dp = 3*sqrt(t0x1X.Hfirst(5,5))*p*p;
	ok &= p0-dp0<p+dp && p-dp<p0+dp0;
	if (!ok) {
	  CsErrLog::msg(elError,__FILE__,__LINE__,
  "Beam track #%d association(<->#%d) cancelled: chi2/NDF=%.1f, p=%.1f+/-%.1f",
			b.Id,ia,chi20x1X,p,dp/3);
	  b.Associate = -1;
	}
      }
    }
  }


  if (TOpt::Det2Go2Fit[0]!="\0") { // Array of dets to restrict fit to, exists
    if (TOpt::ReMode[41]) // Bridging through SM2 yoke is allowed. Then
      // "RESTRICTED SCOPE ReFIT" becomes a real pain => forbid it.
      CsErrLog::mes(elFatal,"Restricted scope refit (\"TraF Det2Go2Fit[0]\") while bridging through SM2 yoke allowed (\"TraF ReMode[41]\")");

    //            ********** RESTRICTED SCOPE ReFIT **********

    //  Tracks are refit while restricting the list of hits to a subset of dets.
    //  If this subset does not bridge a magnet, momentum is fixed (to its
    // full fit value).
    //  If this subset has too few hits, viz. < 6 => discard track. This
    // limit of 6 is chosen so that it's still possible to refit with
    // a single DC station where a double layer has been turned off for exam.
    // The limit is set so that, for other exercice where several different
    // stations are involved, the refit be more reliable. As an example
    // imagine a refit with MWPCs for a track ending in PS01 => There the
    // refit will rely on the sole 4, closely spaced, coordinates of that det.
    // With the limit supra, such track is discarded.
    //  All hits are kept => NdF of fit is not "Chi2tot/(NDFs - # fit pars)".
    THlx Hextr;
    list<TTrack>::iterator it=listTrack.begin(); 
    while (it != listTrack.end()) {
      if (TOpt::ReMode[3]==0 &&            // ***** If this a JOB w/ BRIDGING...
	  !(*it).Hfirst.with_mom()) { it++; continue; }    // ...REQUIRE MOM
      TTrack t = *it;
      list<int>::iterator ih = t.lHitPat.begin();
      list<int>::iterator ip = t.lPlnRef.begin();
      static double xDet; int nDFs = 0, nHits = 0; t.Type = 0;
      while (ih!=t.lHitPat.end()) {
	if (*ih>=0) {
	  const TPlane &p = setup.vPlane(*ip);
	  // ***** RETAIN ONLY THOSE DET'S FLAGGED "2Go2Fit"... *****
	  // CAVEAT: it DOES NOT update sIHit
	  if (!(p.IFlag&0x2)) {
	    const TDetect &d = setup.vDetect(p.IDetRef);
	    for (int igr = 0; igr<(int)setup.vIplFirst().size(); igr++)
	      if (*ip>=setup.vIplFirst()[igr] && *ip<=setup.vIplLast()[igr])
		// ***** ...AND BUILD CORRESPONDING ZONES PATTERN *****
		t.Type |= 1<<igr;
	    xDet = d.X(0);
	    if (nDFs==0) t.Hfirst(0) = xDet;  // Update start helix
	    nDFs += d.IType==29 ? 2 : 1; nHits++;
	  }
	  else {
	    t.lPlnRef.erase(ip++); t.lHitPat.erase(ih++);
	    continue;
	  }
	}
	else if (*ip>=0) {
	  //  No hit but plane is nonetheless there. Means it has been
	  // introduced w/ "InsertMissedPlanes" together with hit placeholder
	  // (and DC stuff)
	  //  We remove these so that smoothing be done only in the refitted
	  // segment.
	  t.lPlnRef.erase(ip++); t.lHitPat.erase(ih++);
	  continue;
	}
	ip++; ih++;
      }
      t.NHits = nHits; t.NDFs = nDFs;
      if (TOpt::ReMode[21]==0) t.InsertMissedPlanes();
      if   (t.NDFs<6) {         // ***** FEW HITS: ERASE *****
	listTrack.erase(it++); continue;
      }
      else {
	if ((t.Type&0x3)!=0x3 &&
	    (t.Type&0x6)!=0x6)   // ***** NO MAGNET BRIDGING: FIX MOMENTUM *****
	  t.Hfirst(5,5) = 1.E-20;
	// ***** >=1 ZONES: VARIOUS UPDATES... *****
	t.Hlast(0) = xDet;                      // ...End helix
	ip = t.lPlnRef.begin(); ih = t.lHitPat.begin();
	while (1) {                             // ...Cut empty front end
	  if ((*ih)<0) {
	    ip = t.lPlnRef.erase(ip); ih = t.lHitPat.erase(ih); continue;
	  }
	  else break;
	  ip++; ih++;
	}
	ip = t.lPlnRef.end(); ih = t.lHitPat.end();
	while (1) {                             // ...Cut empty tail end
	  ip--; ih--;
	  if (*ih>=0) break;
	}
	ip++; t.lPlnRef.erase(ip,t.lPlnRef.end());
	ih++; t.lHitPat.erase(ih,t.lHitPat.end());
      }

      // ***** REFIT *****

      if (t.Hfirst.with_mom()) {
	if (!t.FullKF(1)) {   // Forward  fit has failed. Erase track. Next one.
	  CsErrLog::mes
	    (elWarning,
	     "Forward Kalman Restricted scope Refit failed. Track removed.");
	  listTrack.erase(it++); continue;
	}
	if (!t.FullKF(-1)) { // Backward  fit has failed. Erase track. Next one.
	  CsErrLog::mes
	    (elWarning,
	     "Backward Kalman Restricted scope Refit failed. Track removed.");
	  listTrack.erase(it++); continue;
	}
      }
      else {
	t.QuickKF(1,0);  // straight line forward fit
	t.QuickKF(-1,0); // straight line backward fit	
      }

      // ***** RECOMPUTE "SMOOTHED" HELICES *****
      t.lHsm.clear();
      for (int i = 0; i<int(sizeof(TOpt::SmoothPos)/sizeof(double)); i++) {
	// Loop over "smoothing" points
	double x = TOpt::SmoothPos[i];
	if (x==0.) break;
	if (t.GetSmoothed(Hextr,x)<0) continue;
	t.lHsm.push_back(Hextr);  // store smoothed
      }

      // ***** COPY RESULT of REFIT to ORIGINAL TRACK *****
      // (Now that smoothing has been performed!)

      TTrack &r = *it;
      //#define DEBUG_RESTRICTED_SCOPE
#ifdef DEBUG_RESTRICTED_SCOPE
      r = t;
#else
      Hextr(0)= r.Hfirst(0); t.Hfirst.Extrapolate(Hextr,true); r.Hfirst = Hextr;
      Hextr(0)= r.Hlast(0);  t.Hlast.Extrapolate(Hextr,true);  r.Hlast = Hextr;
      r.Chi2tot = t.Chi2tot;
      r.lHsm = t.lHsm;
#endif
      it++;
    }
  }  // End Restricted Scope Refit

  if (TOpt::ReMode[49]&0x4) // Drell-Yan option &= 0x4... 
    BackTrackVD();          // ...Hit-to-track association in the vertex detector

  if (!isMC) return;

  // ******************** TTrack <-> MCTrack ASSOCIATION ********************
  // - Do this now that all of tracking is finalised (N.B.: Assuming no
  //  tracking refiniment to take place in a yet to be written 2nd pass of
  //  Traffic)
  // - Extrap. tracks upwards to primary vertex if requested,
  for (it = listTrack.begin(); it!=listTrack.end(); it++) { // Loop over tracks
    (*it).FindKine();
    if (TOpt::ReMode[9]==1 && (*it).InGroup(0)==1){
      if( (*it).Hfirst.with_mom() ) { // track with momentum
	THlx Hextr;
	Hextr(0)=vVtxMC(0).V(0);      // Primary vertex. assumed to come first
	if ((*it).Hfirst.Extrapolate(Hextr)) (*it).Hfirst = Hextr;
      }
    }
  }
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ BeamsFit2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void TEv::BeamsFit2()
{
  // - The fit concerns the beamTelescope piece of the beam track (zone 0x10),
  //  and its possible extension to zone 0x1).
  // - The momentum is taken from option.
  // - It's expected to be either permanent or provisional depending upon:
  //   - BMS options (BMS reco, e.g. carried out by CsBeamRecons", or  BMS
  //    simulation, e.g. carried by TraFFiC, cf. ReMode[25] infra).
  //   - Beam (=BMS+beamTelscope), reco or simulation, success.
  // - Permanent: it's expected to be the case for hadron data where there is no
  //  BMS whatsoever.
  // - Provisional: it's used to take into account multiscattering, while the
  //  final momentum is yet unknown.
  // - Provisional converted into Permanent: if beam (=BMS+beamTelscope) turns
  //  out to fail.
  //   This enables to build a pVertex upon it. 
  //  => The uncertainty, specified in "dCut[5], must correspond to uncertainty
  //    prior to any BMS measurement => hence to the beam energy spread.
  // (NOTA BENE: a value corresponding to the beam momentum is available in
  // "dCut[4]", but this quantity serves several purposes and the decision to
  // actually assign it to beam tracks is taken on view of "dCut[5]" (which
  // specifies momentum spread) => If finite, then do ascribe "dCut[4]", w/
  // charge from "iCut[15] and uncertainty from "dCut[5]"). 
  // - In the MC case, and upon option, this method tries and overwrites the
  //  provisional value w/ one taken from the truth bank and smeared.
  // - In any case, false upon return means the track has to be erased.

  //                 ********** BMS SMEARING...**********
  // ...= Smearing of the momentum of the incident track in MC
  double beamP = TOpt::dCut[4], bMSResolution = beamP?.8/beamP/beamP:.8/160/160;
  // - Date: Mon, 3 Nov 2003 10:30:36 ? (MET)
  //   De: 
  //   Martin von Hodenberg <mvhodenb@axfr01.physik.uni-freiburg.de>
  //   In the reconstruction the error on the momentum is set to 0.8 GeV (0.5%).
  // - The random generation is placed here so that it be performed
  //  systematically and hence does not make the evolution of the random seed
  //  depend upon whether the incident track has been reco'd or not. Can be
  //  useful for debugging purposes, so as to get identical events in various
  //  circumstances.
  // - Likewise, in order for the random generation to proceed identically
  //  whether ''retracking'' is called or not, "BMSSmearing" has been made a TEv
  //  member (set =0 by construction) so that it be saved before TEv is deleted
  //  upon retracking, and reinstated right after w/o any new call to CsRandom.
  if (isMC && !reTracking) BMSSmearing = bMSResolution * CsRandom::gauss();

  list<TTrack>::iterator it = listTrack.begin(); while (it!=listTrack.end()) {

    TTrack &t = *it;
    if (t.Type!=0x10 && t.Type!=0x11) { it++; continue; }
    // Let's revisit timing. This time, calling "UseHitTime" w/ a non null
    // "clean" argument. And let's draw the conclusions from it.
    if (!t.UseHitTime(6))  { listTrack.erase(it++); continue; }
    else if ((int)t.NDFs<TOpt::iPRpar[45] && // If #hits < PR requirement...
	     (int)t.NDFs<TOpt::iPRpar[42] && // (1st PR may be less demanding)
	     // ...and large dispersion: 5sigmas, assuming scifi.Resol=.6
	     t.DispTime>3) {
      listTrack.erase(it++); continue;
    }

    if (TOpt::dCut[5]) { //********** MOMENTUM TAKEN FROM OPTION **********
      t.Hfirst(5)  =t.Hlast(5)   = TOpt::iCut[15]/TOpt::dCut[4];
      t.Hfirst(5,5)=t.Hlast(5,5) = TOpt::dCut[5]*TOpt::dCut[5];
      //  A momentum has been assigned to the beam track, possibly only
      // provisionally, so that multiscattering can be taken into account.
      //  Side benefit: the track can later be subjected to "TTRack::Refine",
      // where pseudo-efficiencies and smoothed pulls are histogrammed.
      if (TOpt::ReMode[21]==0) t.InsertMissedPlanes();
      double save = t.Hfirst(5), dsave = t.Hfirst(5,5);
      t.Hfirst(5,5) = 1.E-20; // Fix momentum for the time of the Kalman fit
      bool ret = t.FullKF(1); if (ret) t.FullKF(-1);
      if (!ret) {
	CsErrLog::msg(elWarning,__FILE__,__LINE__,
		      "Beam track %d: Kalman fit failed => remove.",t.Id);
	listTrack.erase(it++); continue;
      }
      t.IFit |= 0x60; // Flag full KF fit
      t.Hfirst(5)=t.Hlast(5) = save; // Momentum was tightly constrained. Yet it may have varied slightly.
      t.Hfirst(5,5)=t.Hlast(5,5) = dsave; // Restore also momentum spread
      for (int i = 1; i<5; i++) {
	t.Hfirst(5,i)=t.Hfirst(i,5)=t.Hlast(5,i)=t.Hlast(i,5) = 0;
      }
    }

    if (isMC) {    // ********** MC BEAM TRACK BMS SIMULATION **********
      // - It's a posteriori as in the RD case: the present block overwrites
      //  whatever assignment may just been done supra.
      // - It's conditioned by ReMode[25]. Various actions:
      //   - Ideal recos: useful if one wants to factorize as completely as
      //    possible the beam reco efficiency out of coral.
      //   - Momentum retrieve from MC truth bank and smearing.

      //                         ***** IDEAL BMS RECO'S...   *****
      // - The overall beam reco is not ideal: the very beam track (IKine==0)
      //  can indeed be NOT reconstructed.
      // - But BMS is assumed to be good enough that bad scifi tracks (IKine<0)
      //  or even non beam tracks (IKine!=0) are rejected.
      if (TOpt::ReMode[25]&0x4) {
	if (t.IKine) {                          // ***** ...SUPER IDEAL BMS RECO
	  listTrack.erase(it++); continue;  // => Erase all IKine != 0
	}
      }
      else if (TOpt::ReMode[25]&0x2) {
	if (t.IKine<0) {                              // ***** ...IDEAL BMS RECO
	  listTrack.erase(it++); continue;  //  => Erase all IKine < 0
	}
      }
      if ((TOpt::ReMode[25]&0x1) && t.IKine>=0) {
	//                         ***** BMS SIMULATION *****
	// - The momentum is retrieved from the MC truth bank...
	// ...provided associated to MC track is available: require IKine>=0.
	const TKine &k = vecKine[t.IKine];
	double pVrtx = fabs(1/k.Pinv()); // Momentum @ vertex
	// Retrieve the values of the energy at first/last measured points
	// (In fact, either of these may not be available (cf. infra). If so, we
	// take instead the value at the closest available point.) 
	vector<int>::const_iterator ih, ihF, ihL, iend = k.vHitMCRef().end();
	static double minF, minL; double xF = t.Hfirst(0), xL = t.Hlast(0);
	for (ih = k.vHitMCRef().begin(), ihF=ihL = iend; ih!=k.vHitMCRef().end();
	     ih++) {
	  const THitMC &h = vecHitMC[*ih]; // Loop over MC track hits...
	  if (h.IOrig) // ...skipping non original hits, for these do no...
	    continue;  // ...convey any info about original track's energy
	  double x = h.Xyz(0);
	  if (ihF==iend || fabs(x-xF)<minF) { ihF = ih; minF = fabs(x-xF); }
	  if (ihL==iend || fabs(x-xL)<minL) { ihL = ih; minL = fabs(x-xL); }
	}
	static double pFirst, pLast;
	if (ihF!=k.vHitMCRef().end() && ihL!=k.vHitMCRef().end()) {
	  CsMCTrkHit *csh = dynamic_cast<CsMCTrkHit*>
	    (const_cast<CsMCHit*>(vecHitMC[*ihF].PtrHit()));
	  // COMGeant extrapolates beam particle backwards, having it
	  // loosing (instead of gaining) energy in the process...
	  // ... => Restore correct energy loss
	  pFirst = 2*pVrtx-csh->getP().mag();
	  csh = dynamic_cast<CsMCTrkHit*>
	    (const_cast<CsMCHit*>(vecHitMC[*ihL].PtrHit()));
	  pLast = 2*pVrtx-csh->getP().mag();
	}
	else {
	  CsErrLog::mes(elError,"BMS simulation unsuccessful => Momentum values of beam track taken from pVertex w/o any energy loss correction");
	  pFirst=pLast = pVrtx;
	}
	static double bMSRes2 = bMSResolution*bMSResolution;     // ***** SMEARING
	t.Hfirst(5) = (1/pFirst + BMSSmearing)*k.Q();
	t.Hlast(5)  = (1/pLast  + BMSSmearing)*k.Q();
	t.Hfirst(5,5)=t.Hlast(5,5) = bMSRes2;
	//  Backward time travel in COMGeant:
	// - The incident particle is assigned an opposite charge before being
	//  propagated backward in time from interaction point through magnetic
	//  field, cf. $COMGeant/code/src/omgbatch/ompro/ombeam.F".
	// - This in principle. But in practice, it's not always the case.
	// => Let's reverse the sign whenever the charge gotten from the truth
	//   bank disagrees w/ the sign of the beam line "iCut[15]".
	if (t.Hfirst(5)*TOpt::iCut[15]<0) {
	  t.Hlast(5) *= -1; t.Hfirst(5) *= -1;
	}
      } // End BMS simulation
    }  // End MC case
    it++; // next track
  }  // End loop on all tracks
}
