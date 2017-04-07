// $Id: TEvFitSegments.cc 14069 2015-09-17 20:44:46Z lsilva $

#include "CsErrLog.h"
#include "CsInit.h"
#include "DaqDataDecoding/DaqOption.h"
#include "CsDetector.h"
#include "CsEvent.h"
#include "CsCluster.h"
#include "TEv.h"
#include "TSetup.h"

using namespace std;

/*! 
  \brief Fit list of single segment tracks
*/

/*
  Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/TracksFit":
   i) Do not straight fit if yet done.
  ii) Remember QN fit before straight fitting.
 iii) muWall#1 cleaning moved in here, so that bridging can benefit
     from it (Cleaning, at the difference of it traffic's one, includes
     cleaning away tracks, that are relying too much on muWall#1 hits,
     and not only hits).
  iv) Correct for time propagation (and X-ray and trigger jitter and time offset
  of pile-up tracks) in drift detectors.
   v) PixelGEM/MMs.
*/

void TEv::FitSegments()
{
  //    ***** In view of UPDATING DRIFTS THits w/ EVENT TIME ***** 
  //       ***** REQUIRE THAT TRIGGER MATCHES ReMode[36]... *****
  const unsigned int allTrigs = 0xffff; // Although TCS can handle 0x7fffff, it was decided to limit max. #triggers to 16. And take advantage of the fact to make use of the other bits of the trigger TDC module.
  unsigned int evtTrig = trig_mask;
  evtTrig &= allTrigs;  // Cut away trailing end bits (i.e. online filter...)
  unsigned int inclMask =// Trigger pattern for which bits we require only an...
    TOpt::ReMode[36]&(~TOpt::ReMode[40]);// ...inclusive match, cf. UpdateDrifts
  double evtT0 =
    (((evtTrig&TOpt::ReMode[36])==evtTrig || (evtTrig&inclMask)) &&
     // The following can happen in MC: event did not pass the trigger...
     evtTrig || // => Exclude the case. No possible trigger offset anyway.
     //  ***** ...OR THAT ALTERNATIVE MASTER TRIGGER IS DEFINED *****
     ptrEvt()->getAlterMasterTrigger()>=0 ||
     // ...OR 2ND ITERATION of TRACKING GOING ON
     ptrEvt()->reTrackingON()) ? eventTime : 0;


  // Simple fit by straight line of track segments

  const TSetup& setup = TSetup::Ref();

  list<TTrack>::iterator it=listTrack.begin(); 
  while (it!=listTrack.end()) {       // ********** LOOP OVER TRACKS **********
    TTrack& t = (*it);
    t.UseHitTime();

    if (t.NGroups()==0) { CsErrLog::mes(elFatal,"TTrack with NGroups == 0!"); }
    if (t.NGroups()!=1) { it++; continue; }// Not track segment but "long" track
    if (!(t.IFit&0x1)) {                         // No straight fit yet done

      if (!t.QuickKF(1,0)) { // forward  fit has failed. Erase track. Next one.
	CsErrLog::mes(elWarning,
		      "Forward Kalman Straight Fit failed. Track removed.");
	listTrack.erase(it++); 
	continue;
      }
      if (!t.QuickKF(-1,0)) {// backward  fit has failed. Erase track. Next one.
	CsErrLog::mes(elWarning,
		      "Forward Kalman Straight Fit failed. Track removed.");
	listTrack.erase(it++);
	continue;
      }
      t.IFit |= 0x1;
    }

    // ***** ZONE 0X2: REMOVE "RANDOM PICKUP" downstream HCAL1 and muF1 *****

    // Rationale:
    //  - The pick-up hits will automatically grant their host track a large
    //   %X0. => Their association need be very reliable.
    //  - Hits are available in large quantities there, supplied by showers in
    //   the thickness of HCAL and muF. A few of them can easily be picked up
    //   by any track, at random.

    // =>
    //  - Remove hits in PAs and MAs if they are not numerous enough.
    //  - Remove tracks which are essentially built on this kind of hits, to
    //   prevent them from stealing hits from better defined genuine muons.
    //   Cf., e.g., evt #253 in "dstar.2003.01.fz.2".

    // Do this now to prevent this pick-up from spoiling track segments
    // before they can be bridged (over SM1)

    // We will consider only tracks that go through ECAL (not yet), HCAL and
    // muF. Need know the positions of these elements and central aperture.
    double xMF1 = TOpt::MuonWall[0], yMF1[2], zMF1[2];
    double xHC1 = TOpt::Calo[0],     yHC1[2], zHC1[2];
    for (int i = 0; i<2; i++) {
      int j = 2*i-1;
      yMF1[i] = TOpt::MuonWall[6]+j*TOpt::MuonWall[8];
      zMF1[i] = TOpt::MuonWall[7]+j*TOpt::MuonWall[9];
      yHC1[i] = TOpt::Calo[6]    +j*TOpt::Calo[8];
      zHC1[i] = TOpt::Calo[7]    +j*TOpt::Calo[9];
    }
    if (t.Type==0x2 && TOpt::iCut[6] && t.Hlast(0)>xHC1) {

      // Extrapolate to MF1 and HC1(no need for precision given the probably
      // loose positioning of these big blocs) => No error propagation =>
      float ytMF1 = t.Hlast(1)+t.Hlast(3)*(xMF1-t.Hlast(0));
      float ztMF1 = t.Hlast(2)+t.Hlast(4)*(xMF1-t.Hlast(0));
      float ytHC1 = t.Hlast(1)+t.Hlast(3)*(xHC1-t.Hlast(0));
      float ztHC1 = t.Hlast(2)+t.Hlast(4)*(xHC1-t.Hlast(0));
      float xStrt = 0, nWHitsMn = (float)TOpt::iCut[6];
      if      (ytHC1<yHC1[0] || yHC1[1]<ytHC1 ||
	       ztHC1<zHC1[0] || zHC1[1]<ztHC1)
	xStrt = xHC1;
      else if (ytMF1<yMF1[0] || yMF1[1]<ytMF1 ||
	       ztMF1<zMF1[0] || zMF1[1]<ztMF1) {
	xStrt = xMF1; nWHitsMn /= 2;
      }

      if (xStrt>0 && t.Hlast(0)>xStrt) {
	// - Count #hits. Weight = 2.5 on MWPCs in order to retain tracks
	//  through MA1 dead zones. Typically cut on #hits will be 12.
	//  6 MWPCs - 1 (inefficient) * 2.5 > 12 => OK
	//  On track avoiding HC1 but traversing MF1, we would require half
	//  fewer = 6. And 3-MWPC-hit track would fulfill it.
	int nMWAs = 0, nPAs = 0, nHits = 0, ipli = -1; double xplf = 0;
	// - Evaluate quality of track upstream of "xStrt": count # of ``space
	//  points.
	int nSpacePts = 0; TStation *sPrv = 0/*, nStations = 0 */;
	int projs = 0, nProjs = 0, allProjs = 0, nAllProjs = 0;
	list<int>::iterator ipr = t.lPlnRef.begin(), ihp = t.lHitPat.begin(); 
	while (ipr!=t.lPlnRef.end()) {   // Loop over track's plane references
	  if (*ihp<0) { ipr++; ihp++; continue; }
	  int ipl = *ipr; const TDetect &d = setup.iPlane2Detect(ipl);
	  double xD = d.X(0); if (xD>xHC1) {
	    if (xD>xStrt) {
	      if (setup.InMuWallA(ipl)) nMWAs++; // MA or HG02 hit
	      else if (d.IType==1)      nPAs++;  // MWPC hit
	    }
	    if (ipli<0) ipli = ipl; xplf = xD; // Remember first/last plane
	    // (Note that it can have Z < "xStrt". For we have to eliminate all
	    // of MAs (and also MWPCs (why not?)) even if it seems they were
	    // triggered before the track has traversed any significant
	    // thickness of material, so as to avoid any ambiguity (coral users
	    // tend to think that hits in MAs are a muID per se) and because our
	    // evaluation of whether HC1 or MF1 were indeed traversed is here
	    // very crude and the precise evaluation of %X0 could turn out to
	    // oveturn it.)
	    nHits++;  // Count all hits
	  }
	  if (xD<xStrt) {
	    const TPlane &p = setup.vPlane(ipl);
	    const TStation *&s = p.Station; if (s!=sPrv) {
	      sPrv = const_cast<TStation*&>(s); /* nStations++; */
	      if (nProjs>=2) { nSpacePts++; nProjs=projs = 0; }
	    }
	    int jproj = 1<<p.IProj;
	    if ((jproj&projs)==0)        {    projs |= jproj;    nProjs++; }
	    if ((p.IFlag&0x30) && (0x2&projs)==0)
	      /* Pixel plane: Z-proj. */ {    projs |= 0x2;      nProjs++; }
	    if ((jproj&allProjs)==0)     { allProjs |= jproj; nAllProjs++; }
	  }
	  ipr++; ihp++;
	}
	if (nProjs>=2) nSpacePts++;

	// ***** REMOVE TRACKS RELYING TOO MUCH on DOWNSTREAM ABSORBER *****
	if (nSpacePts<2 || nAllProjs<3) {
	  listTrack.erase(it++); continue;
	}
	// MuWallA reconstruction...
	bool isMWA = xplf>xMF1 &&  // ...require hits downstream of MF1
	  // (This in order to tackle the case of tracks that, though not muons,
	  // are able to pick up hits from the PA|MA01s, since these are not that
	  // well shielded against showers and delta rays. Requiring hits
	  // downsteam of the good shielding provided by MF1 will prevent these
	  // tracks from being unduly granted a muID (thanks to the %X0 collected
	  // while traversing HC1 on their way to the PA/MA01s).)
	  nMWAs+nPAs>=TOpt::iCut[6];
	// Special case of PA01+PA02 reconstruction
	int isPA = nPAs>=5;
	if (ipli>=0 && !isMWA && !isPA) {
	  // ***** TOO FEW HITS DOWNSTREAM OF ABSORBER => REMOVE THEM *****
	  if ((int)t.NDFs-nHits<TOpt::iPRpar[15]) {
	    // Very few hits left over. Or, even worse, no hit left at all!
	    listTrack.erase(it++); continue;     // => erase track
	  }
	  t.Shorten(ipli);
	  if (!t.QuickKF(1,0)) {// forward  fit failed. Erase track. Next one.
	    CsErrLog::mes
	      (elWarning,
	       "Forward Kalman Straight Refit failed. Track removed.");
	    listTrack.erase(it++); 
	    continue;
	  }
	  if (!t.QuickKF(-1,0)) {// Backward  fit failed. Erase track. Next.
	    CsErrLog::mes
	      (elWarning,
	       "Backward Kalman Straight Refit failed. Track removed.");
	    listTrack.erase(it++);
	    continue;
	  }
	}
      }  // End check track passing through zone 0x2 absorber
    }  // end of "if iCut[6] ..."

    if (TOpt::ReMode[29]&0x1) {

      // ***** CORRECT PROPAGATION in DRIFTS *****

      // Status: -1: No init yet, 0: Init done,
      //          1: Modified
      int status = -1;
      double X0 = t.Hfirst(0);
      double y0 = t.Hfirst(1)*10, dydx = 10*t.Hfirst(3);
      double z0 = t.Hfirst(2)*10, dzdx = 10*t.Hfirst(4);
      double tT0 = evtT0;
      double evtT = eventTRef ? eventTRef : eventTime;
      //#define FS_DEBUG_PROPAG
#ifdef FS_DEBUG_PROPAG
      double yd, zd;
#endif
      list<int>::iterator ih, ip;
      for (ih = it->lHitPat.begin(), ip = it->lPlnRef.begin();
	   ih!=it->lHitPat.end(); ih++, ip++) {
	if (*ih<0) continue;
	const TPlane &p = setup.vPlane(*ip);
	const TDetect &d = setup.vDetect(p.IDetRef);
	CsDetector *csDet = d.PtrDet(); if (!csDet->hasDrift()) continue;
	if (status<0) {               // First detector encountered...
	  t.UseHitTime(); status = 0; // ... => mean time
	  // It's a halo muon accidentally coincident?
	  if (t.DispTime>0 && // This requires 2 time-measuring detectors...
	      // ...eliminating the risk to have track's time set by ghost hit.
	      fabs(t.MeanTime-evtT)/t.SigmaTime>2 &&       // Off time?
	      // Note: When timing undef: "SigmaTime"<0
	      fabs(dydx)<.05 && fabs(dzdx)<.05 &&     // Paraxial?
	      t.Chi2tot/(t.NDFs-4)>10)
	    tT0 = t.MeanTime-eventTRef;              // => Set latency "tT"
#ifdef FS_DEBUG_PROPAG
	  yd =  y0+dydx*(d.X(0)-X0); zd = z0+dzdx*(d.X(0)-X0);
#endif
	}
	THit &h = vecHit[*ih];
	if (h.Mirror!=-1 && vecHit[h.Mirror].Status==-4) // Coalesced Hit...
	  continue;  // ...keep it coalesced, meaning no propagation effect
	double y = y0+dydx*(d.X(0)-X0), z = z0+dzdx*(d.X(0)-X0);
	bool error; double U = csDet->getCorrU(h.PtrClus(),y,z,tT0,error);
	if (!error) {
	  h.u = U/10; h.Rotate(); status = 1;
	}
	//...else: For the time being, do not delete hit, although...
      }        // End loop on planes
      if (status>0) {
#ifdef FS_DEBUG_PROPAG
	printf("1 %3d %5.1f  %5.1f %5.1f %5.1f %c\n",
	       t.Id,t.Chi2tot/(t.NDFs-4),yd,zd,
	       t.MeanTime,tT0!=0 ? 'c' : ' ');
#endif
	t.QuickKF(1,0);  // straight line forward fit
	t.QuickKF(-1,0); // straight line backward fit
#ifdef FS_DEBUG_PROPAG
	printf("      %5.1f\n",t.Chi2tot/(t.NDFs-4));
#endif
      }
    }        // End correct for propagation in drifts
    it++;
  }

  // In case of MC event:
  // - associate found track with MC track.
  if (!isMC) return;

  for(it = listTrack.begin(); it != listTrack.end(); it++) {  // loop over tracks
#ifdef DEBUG
    TTrack *t = &*it;
#endif
    (*it).FindKine();
  }


  return;
}
