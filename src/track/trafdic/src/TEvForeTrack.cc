// $Id: TEvForeTrack.cc 14094 2015-11-06 15:28:48Z lsilva $

#include "CsErrLog.h"
#include "CsOpt.h"
#include "CsInit.h"
#include "DaqDataDecoding/DaqOption.h"
#include "CsHistograms.h"
#include "CsEventUtils.h"
#include "CsDetector.h"
#include "CsEvent.h"
#include "Traffic.h"
#include "TOpt.h"
#include "TSetup.h"
#include "TEv.h"
// The following headers have to be included last, for they conflict
// w/ "CsChip.h"
#include "TMatrixF.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"

using namespace std;

/*! 
  \brief Various methods for extending TTrack's by ForeTracking, i.e. extrapolate (forward) + pick-up.
  \author  Yann.Bedfer@cern.ch
*/

/*
  - "ForeTrack2RICHWall" = RICHWall in the broadest sense, including PS01,
    "ForeTrack2MA" = MA and HG02,
    "ForeTrack2Hs" = Hodos downstream of SM2.
  - "TTrack::IKine" (and "NHsame") not re-evaluated: it's done in "TracksFit2" (
   and "TracksRefit").
 */


void TEv::ForeTrack2RICHWall()
{
  // ********** EXTEND SM1-Bridged TRACKS TO RICHWall(*) **********
  // (*) RICHWall in the broadest sense of all detectors at exit of RICH,
  // including PS01, or ST04, in addition to DR proper.

  // (PR often fails to extend zone 0x2 track segments into RICH Wall detectors
  // (understood as all LAT detectors in this region), particularly in the
  // bending dimensions (X,U,V), because of the long distance between these and
  // detectors in the SM1<->RICH region and of the thickness of material
  // traversed. This method aims at achieving this extension. It is performed
  // on tracks already (SM1-)bridged and KF'ed. Because:
  //  - We do not want to waste time on unreliable tracks.
  //  - KF, in forward direction, gives a better estimate of the track helix
  //   at the far end of the track than straight fit, and takes MS into
  //   account.)

  // N.B.: This method may erase some tracks.

  if (TOpt::dCut[67]==0) return;
  const TSetup &setup = TSetup::Ref();
  static int ipli, iplf, first = 1, nPlsRW; int ipl, jpl;
  static unsigned int plPat, plYUVPat, plDRPat;
  static float dzYMn, dzZMn, xMn, accYMx, xYMx, accZMx, xZMx;
  static TH2D *hChi2 = 0, *hdChi2 = 0, *hNHits = 0, *hShower = 0;
  // - Determine TPlane <-> proj. association. This is necessitated by 2006
  //  setup with "DR" RICH wall (as opposed to "ST04"), for which we want to
  //  distinguish between 1st, DR01, and 2nd, DR02, station
  static unsigned int *p2Projs;
  static int iDRStation = -1;
  if (first) { 

    // *************** INITIALISATION ***************

    // ***** DETERMINE WHICH PLANES BELONG to RICH WALL *****
    //       => PLANE# RANGE [ipli,iplf] and PLANE# PATTERNs
    //       => MIN and MAX ACEPTANCE
    first = 0;
    for (ipl = setup.vIplFirst()[1], ipli=iplf = -1, nPlsRW = 0,
	   plYUVPat=plPat=plDRPat = 0, dzZMn=dzYMn = 10000, accZMx=accYMx = 0;
	 ipl<=setup.vIplLast()[1]; ipl++) {
      const TPlane &pi = setup.vPlane(ipl); const TStation *&st = pi.Station;
      int spl, nspl = (int)st->IPlanes.size(), sIsWin; ipl += nspl-1;  
      for (spl = 0, sIsWin = -1; spl<nspl; spl++) {
	// Do loop on TStation's and then on planes w/in TStation's, since this
	// what is done infra, in the body of "ForeTrack2RICHWall" code. It
	// will allow to reject cases where RICH wall's range straddles that of
	// one of its TStation's, which would cause the code to derail.
	int kpl = st->IPlanes[spl]; const TPlane &ps = setup.vPlane(kpl);
	const TDetect &d = setup.vDetect(ps.IDetRef);
	if (d.IType==22 || d.IType==26) break;        // Exclude GEMs and scifis
	if (!ps.IFlag) continue; // Plane is OFF
	int pIsWin = TOpt::dCut[67]<d.X(0) && d.X(0)<TOpt::dCut[68] ? 1 : 0;
	if (sIsWin<0) sIsWin = pIsWin;
	if (pIsWin!=sIsWin) {        // TStation cannot straddle RICH wall range
	  const TPlane &psi = setup.vPlane(st->IPlanes[0]);
	  const TDetect &di = setup.vDetect(psi.IDetRef);
	  const TPlane &psf = setup.vPlane(st->IPlanes[nspl-1]);
	  const TDetect &df = setup.vDetect(psf.IDetRef);
	  CsErrLog::msg(elFatal,__FILE__,__LINE__,
   "Station \"%s\"...\"%s\" @ [%.1f,%.1f] cm straddles RICHwall zone (=[%.1f,%.1f] cm) defined by TraF dCut[67-68]",
			di.Name.c_str(),df.Name.c_str(),di.X(0),df.X(0),
			TOpt::dCut[67],TOpt::dCut[68]);
	}
	if (pIsWin!=1) continue;              // Require TPlane to be w/in range
	if (nPlsRW++>=32)
	  CsErrLog::msg(elFatal,__FILE__,__LINE__,
   "Too many(>32) detectors in RICHwall zone [%.1f,%.1f] defined by TraF dCut[67-68]",
		      TOpt::dCut[67],TOpt::dCut[68]);
	if (ipli==-1) ipli = kpl; iplf = kpl; plPat |= 1<<kpl-ipli;
	if (fabs(d.Ca)>.04) plYUVPat |= 1<<kpl-ipli;
	// Planes to later checked for acceptance. The X dimension, more
	// relevant, is considered
	float dzYDim = d.DZSize(0); // Can be null, e.g. for an outer ST04 slice
	if (dzYDim && dzYDim<dzYMn) {                 // Smallest dead zone
	  dzYMn = dzYDim; dzZMn = d.DZSize(1); xMn = d.X(0);
	}
	if (d.Siz(1)>accYMx) {   // Largest Y acceptance
	  accYMx = d.Siz(1);                   xYMx = d.X(0);
	}
	if (d.Siz(2)>accZMx) {   // Z can be distinct: think of ST04 slices
	  accZMx = d.Siz(2);                   xZMx = d.X(0);
	}
	if (d.Name.find("DR")==0) {
	  iDRStation = st->IStation; plDRPat |= 1<<kpl-ipli;
	}
      }
    }
    if (ipli==-1)
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
   "No relevant detector in RICHwall zone [%.1f,%.1f] defined by TraF dCut[67-68]",
		    TOpt::dCut[67],TOpt::dCut[68]);

    //        ***** TPlane <-> PROJ. ASSOCIATION *****
    p2Projs = new unsigned int[iplf-ipli+1];
    memset((void*)p2Projs,0,(iplf-ipli+1)*sizeof(unsigned int));
    for (ipl = ipli; ipl<=iplf; ipl++) {
      const TPlane &p = setup.vPlane(ipl); 
      jpl = ipl-ipli; if (!(1<<jpl&plPat)) continue;
      const TDetect &d = setup.vDetect(p.IDetRef);
      p2Projs[jpl] = 1<<p.IProj;
      if (d.Name[1]=='R' && d.Name[3]=='2') { // Use special proj.'s w/ # >=16
	// (There are certainly fewer than 16 proj. in all of COMPASS...)
	if (d.Name[4]=='X') p2Projs[jpl] |= 0x10000;
	else                p2Projs[jpl] |= 0x20000;
      }
    }

    if (isMC && TOpt::Hist[7]) {
      // ********** HISTOS **********
      CsHistograms::SetCurrentPath("/Traffic/MCmonitor");
      hChi2  = new TH2D("hRWChi2","RichWall tracking perf vs. RW-#chi^{2}",
			10,0,TOpt::dCut[69],25,-1/48.,1+1/48.);
      hdChi2 = new TH2D("hRWdChi2",
			"RichWall tracking perf vs. track's #Delta#chi^{2}",
			10,-TOpt::dCut[17]/5,TOpt::dCut[17]/5,25,-1/48.,1+1/48.);
      hNHits = new TH2D("hRWNHits","RichWall tracking perf vs. #Hits",
			50,-.5,49.5,25,-1/48.,1+1/48.);
      hShower = new TH2D("hRWShower","RichWall showerID",2,-.5,1.5,3,-.5,2.5);
      CsHistograms::SetCurrentPath("/");
    }
  }

  //       *************** PRELIMINARY EVENT STEPS ***************
  // - Reassess "THit::Status":
  // - Determine #hits in all of RichWall => Give up if too many hits.
  list<unsigned int> fishyTrks; int nHitsRW;
  UpdateHitStatus(ipli,iplf,(unsigned long long)plPat,fishyTrks,nHitsRW);
  if (nHitsRW>nPlsRW*TOpt::iCut[25]) {// Don't waste CPU on pathologic events...
    // ...=> Give up if too many hits in the RICHwall zone. Yet tune the cut
    // in proportion to the number of candidate tracks to be checked.
    int nCandidateTrks; list<TTrack>::iterator it; 
    for (it = listTrack.begin(), nCandidateTrks = 0; it!=listTrack.end(); it++)
      if ((*it).Hfirst.with_mom() && // Candidate tracks: require momentum...
	  (*it).Hlast(0)<TOpt::dCut[67]) // ...require end upstream of RICHwall
	nCandidateTrks++;
    if (nHitsRW>nPlsRW*TOpt::iCut[25]+24*nCandidateTrks) {
      CsErrLog::msg(elError,__FILE__,__LINE__,
        "Too many hits (%d) in RichWall for %d candidate tracks: give up",
		    nHitsRW,nCandidateTrks); return;
    }
  }
  //     ***** In view of UPDATING DRIFTS THits w/ EVENT TIME ***** 
  //        ***** REQUIRE THAT TRIGGER MATCHES ReMode[36]... *****
  const unsigned int allTrigs = 0xffff; // Although TCS can handle 0x7fffff, it was decided to limit max. #triggers to 16. And take advantage of the fact to make use of the other bits of the trigger TDC module.
  unsigned int evtTrig = trig_mask&allTrigs; // Cut away trailing bits
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

#ifdef DEBUG
  static int idebug = 0;
#endif

  //         ********** DETECT SHOWERS **********
  // - Very preliminary => Simply count the #hits compatible w/ each given
  //  track. And require this #hits to be > some cut.
  // - If a shower is found to be associated w/ a given track, all extensions
  //  of the track found in the RICHWall zone are bound to be unreliable.
  //  => Therefore we cancel them.
  int *showerHits = new int[iplf-ipli+1], *showerFlags = new int[iplf-ipli+1];

  list<TTrack>::iterator it; int iter, iTrk, nTrkPats = listTrack.size()/32+1;
  unsigned int *trkPats = new unsigned int[nTrkPats];
  for (iter = 0, memset((void*)trkPats,0,nTrkPats*sizeof(unsigned int));
       iter<3; iter++) {
    it = listTrack.begin(); iTrk = 0; while (it!=listTrack.end()) {

    //       ******************** LOOP on TRACKS ********************

    int iPat = iTrk/32; unsigned int bTrk = 1<<iTrk%32; iTrk++;
    if (bTrk&trkPats[iPat]) { it++; continue; }

    TTrack &t = *it;
    if (!iter) {   // ********** 2 ITER'S: FROM BEST -> WORST TRACKS **********
      if ((t.Type&0x7)!=0x3 ||           // ***** SM1(ONLY)-BRIDGED TRACKS *****
	  fabs(t.Pinv())>1/TOpt::dCut[70] ||                // ***** P CUT *****
	  !(t.IFit&0x60)) {                      // ***** FULL KF REQUIRED *****
	it++; trkPats[iPat] |= bTrk; continue;
      }
      if (t.Chi2tot/(t.NDFs-5)<3 || fabs(t.Pinv())>.4) {
	it++; continue;             // ***** 1ST: LOW CHI2, HIGH P TRACKS *****
      }
    }
    else if (iter<2) {
      if (t.Chi2tot/(t.NDFs-5)<3) { // ***** 2ND: LOWER P TRACKS ***** 
	it++; continue;
      }
    }                               // ***** 3RD: REMAINING TRACKS *****
    trkPats[iPat] |= bTrk;

    float x0 = t.Hlast(0);  // Intercepts and slopes taken from last helix
    if (x0>TOpt::dCut[68])         // ***** NO RECO DOWNSTREAM of RichWall *****
      { it++; continue; }
    float y0 = t.Hlast(1), z0 = t.Hlast(2);
    float yp = t.Hlast(3), zp = t.Hlast(4);
    //                                              ***** CHECK ACCEPTANCE *****
    if (fabs(y0+(xYMx-x0)*yp)>accYMx || fabs(z0+(xZMx-x0)*zp)>accZMx ||
	fabs(y0+(xMn -x0)*yp)< dzYMn && fabs(z0+(xMn -x0)*zp)< dzZMn) {
      it++; continue; }
    // ********** INITIAL HITS PATTERN. ARE THERE Z HITS? **********
    // If Z hits, they are made the base for the extension. They are not put
    // into question: the fact that they made it in the PR, means that they are
    // reliable enough: they must be in good enough  agreement with the rest,
    // which indirectly means the corresponding track  did not suffer too large
    // a kick while traversing RICH  thickness.
    static const TPlane *pZ, *pZp, *pYp, *pUp, *pVp;
    static THit *hZ, *hZp;
    list<int>::iterator ih, ip; unsigned int hitPlPat; int nDowns, iplZ, nZs;
    for (ih = t.lHitPat.begin(), ip = t.lPlnRef.begin(), hitPlPat = 0,
	   nDowns=nZs = 0, pZ = NULL, iplZ = -1, hZp = NULL;
	 ih!=t.lHitPat.end(); ih++, ip++) {
      if (*ih<0 || *ip<ipli || iplf<*ip) continue;
      int bpl = 1<<*ip-ipli; if (bpl&plPat) {
	hitPlPat |= bpl;
	if (!(bpl&plYUVPat)) {
	  if (pZ) {
	    pZp = &setup.vPlane(*ip);
	    if (pZp==pZ->Associate) hZp = &vecHit[*ih];
	    else {
	      nZs = 2; break; // Several unrelated Z hits => too complex case
	    }
	  }
	  else {
	    pZ = &setup.vPlane(*ip); hZ = &vecHit[*ih]; iplZ = *ip;
	  }
	}
      }
      else nDowns++;
    }
    //                                      ***** EXCLUDE TRACKS RECO'D... *****
    if (nDowns>0) { it++; continue; }         // ...DOWNSTREAM of RichWall *****
    if (hitPlPat&plYUVPat) { it++; continue; }    // ...in YUV of RichWall *****
    if (nZs<0) { it++; continue; }                 // ...W/ TOO MAY Z HITS *****
#  ifdef DEBUG
    if (idebug) it->DumpHits(idebug);
#  endif
    int iplLast = t.lPlnRef.back(); // TTrack's last TPlane

    // *************** TRACK CANDIDATE for EXTENSION ***************

    vector<TSpacePt> spts; TSpacePt spt, spt4;
    static float y, z, u, c, s, dy2, dz2, dyz, du, dyp2, dzp2, dyzp, dup, umx, umn;
    vector<int>::const_iterator ihZ, ihY, ihU, ihV; int firstPl;
    unsigned int usedPlPatZ;
    for (ipl = iplZ<0?ipli:iplZ, firstPl = 1, usedPlPatZ = 0,
	   memset((void*)showerHits,0,(iplf-ipli+1)*sizeof(int)),
	   memset((void*)showerFlags,0,(iplf-ipli+1)*sizeof(int));
	 ipl<=(iplZ<0?iplf:iplZ); ipl++) {
      int jpl = ipl-ipli, bpl = 1<<jpl;
      if ((bpl&plPat)==0)                   // Doesn't belong to ``RichWall''...
	continue;                             // ...because not a LAT, or is OFF
      if (bpl&plYUVPat) continue;                    // Searching Z planes first

      // ***** LOOP OVER Z PLANES in RichWall *****

      pZ = &setup.vPlane(ipl); usedPlPatZ |= bpl; showerFlags[jpl]++;
      const TDetect &dZ = setup.vDetect(pZ->IDetRef);
      float sZ = dZ.Sa, cZ = dZ.Ca, xZ = dZ.X(0), wZ = dZ.Resol;
      THlx hlx; if (firstPl) {
	// Compute uncertainties upon first detector only. In order to
	// speed things up and since uncertainties are not expected to
	// change much over the length of the RichWall.
	hlx(0) = xZ; t.Hlast.Extrapolate(hlx,true); y = hlx(1); z = hlx(2);
	dy2 = hlx(1,1); dz2 = hlx(2,2); dyz = hlx(1,2); firstPl = 0;
	dyp2 = hlx(3,3); dzp2 = hlx(4,4); dyzp = hlx(3,4);
      }
      else {
	y = y0+yp*(xZ-x0); z = z0+zp*(xZ-x0);
      }
      if (iplZ<0) {
	if (!dZ.InActive(y,z)) continue;                  // Out of active zone
	u = y*cZ+z*sZ; du = sqrt(dy2*cZ*cZ+2*dyz*cZ*sZ+dz2*sZ*sZ);
	umx = u+6*du+4*wZ; umn = u-6*du-4*wZ;
      }
      for (ihZ = pZ->vHitRef().begin(); ihZ!=pZ->vHitRef().end(); ihZ++) {
	if (iplZ<0) {
	  THit *h = &vecHit[*ihZ]; if (h->Status!=0) continue;
	  if (h->U<umn) continue; if (umx<h->U) break;
	  hZ = h; hZp = NULL;
	}

	// ********** FOUND AN EXTENSION in Z **********

	const TStation *&st = pZ->Station;
	if (showerFlags[jpl]==1) {
	  int mirror = hZ->Mirror;
	  showerHits[jpl] += mirror>0 && vecHit[mirror].Status==-4 ? 2 : 1;
	}
	float Z = hZ->U, whZ = hZ->SigU;
	c = cZ/whZ; s = sZ/whZ; u = Z/whZ;
	double sc2Z = c*c, sscZ = s*c, ss2Z = s*s;
	int nhZ = 1; unsigned int hPatZ = 0x1;
	double scuZ  = u*c, ssuZ  = u*s, su2Z  = u*u;
#define ForeTrack_UPDATE_ZP
#ifdef ForeTrack_UPDATE_ZP
	// New value of vertical slope, unless reco'd in RICHWall zone by
	// PrePattern => "zpp" peculiar to current "hZ". "zp" is kept unchanged.
	float zpp = hitPlPat&plPat ? zp : (zp+(Z-z0)/(xZ-x0))/2;
#else
	float zpp = zp;
#endif
	static float xZp, cZp, sZp;
	if (!hZp) {
	  pZp = pZ->Associate; if (pZp && pZp->IFlag) {
	    if (1<<pZp->IPlane-ipli&usedPlPatZ) // Associated plane already
	      // scanned: the combination (hit,associated) has probably
	      continue; // already been checked => skip
	    getAssociatedHit(yp,zpp,cZ,sZ,pZ,hZ,pZp,&hZp);
	  }
	}
	if (hZp) {
	  const TDetect &dZp = setup.vDetect(pZp->IDetRef);
	  cZp = dZp.Ca; sZp = dZp.Sa; xZp = dZp.X(0);
	  // New determination of vertical slope
	  float Zp = hZp->U, whZp = hZp->SigU;
#ifdef ForeTrack_UPDATE_ZP
	  zpp = hitPlPat&plPat ? zp : (zp+(Z-z0)/(xZ-x0)+(Zp-z0)/(xZp-x0))/3;
#endif
	  Zp += (xZ-xZp)*(yp*cZp+zpp*sZp);
	  c = cZp/whZp; s = sZp/whZp; u = Zp/whZp;
	  sc2Z += c*c; sscZ += s*c; ss2Z += s*s; 
	  scuZ += u*c; ssuZ += u*s; su2Z += u*u; nhZ++; hPatZ = 0x3;
	}
	unsigned int projZ = p2Projs[jpl], projY, projU, projV;
	unsigned int sProjPatZ = projZ;

	unsigned int usedPlPatY, usedPlPatU, usedPlPatV;
	vector<int>::const_iterator ipY, ipU, ipV = st->IPlanes.end();
	for (ipY = st->IPlanes.begin(), projY=usedPlPatY = 0;
	     ipY!=st->IPlanes.end(); ipY++) {

	  // ***** LOOP ON PLANES FOR ``Y'' (``Y'' is anything but Z) *****

	  const TPlane *pY = &setup.vPlane(*ipY); int jplY = *ipY-ipli;
	  if (pY->IProj==pZ->IProj || // Instead of "p2Projs[jplY]&sProjPatZ"...
	      // ...in order to exclude also Z=DR01Y and Y=DR02Y
	      projY && p2Projs[jplY]!=projY) continue;
	  if (!pY->IFlag) continue;                      // Plane OFF
	  if (pY->vHitRef().empty()) continue;           // Plane empty
	  const TDetect &dY = setup.vDetect(pY->IDetRef);
	  float cY = dY.Ca, sY = dY.Sa,  xY = dY.X(0), wY = dY.Resol;
	  y = y0+yp*(xY-x0); z = z0+zpp*(xY-x0);
	  if (!dY.InActive(y,z)) continue;               // Out of active zone
	  unsigned int sProjPatZY = sProjPatZ, sProjPatZYU = sProjPatZ;
	  showerFlags[jplY]++;
	  u = y*cY+z*sY; du = sqrt(dy2*cY*cY+2*dyz*cY*sY+dz2*sY*sY);
	  float umxY = u+6*du+4*wY, umnY = u-6*du-4*wY;
#ifdef DEBUG
	  if (idebug) {
	    if (isMC) {
	      printf("t #%d Z %.2f,%d,%d ",t.IKine,hZ->U,hZ->IKine,hZ->IKin2);
	      if (hZp) printf("Z' %.2f,%d,%d ",hZp->U,hZp->IKine,hZp->IKin2);
	    }
	    printf("u %.3f du %.3f => %.3f %.3f\n",u,du,umnY,umxY);
	  }
#endif
	  for (ihY = pY->vHitRef().begin(); ihY!=pY->vHitRef().end(); ihY++) {

	    // ***** LOOP ON ``Y'' HITS *****

	    THit *hY = &vecHit[*ihY], *hYp = NULL; if (hY->Status!=0) continue;
	    float Y = hY->U; if (Y<umnY) continue; if (umxY<Y) break;
	    // Set "projY". Do this only now so that it is NOT required that
	    // the first encountered non Z be efficient
	    projY = p2Projs[jplY]; sProjPatZY |= projY; sProjPatZYU |= projY;
	    usedPlPatY |= 1<<jplY;
	    if (showerFlags[jplY]==1) {
	      int mirror = hY->Mirror;
	      showerHits[jplY] += mirror>0 && vecHit[mirror].Status==-4 ? 2 : 1;
	    }

	    int nhZY = nhZ+1; unsigned int hPatZY = hPatZ | 0x4;
	    Y += (xZ-xY)*(yp*cY+zpp*sY);
	    float whY = hY->SigU; c = cY/whY; s = sY/whY; u = Y/whY;
	    double sc2ZY = sc2Z+c*c, sscZY = sscZ+s*c, ss2ZY = ss2Z+s*s;
	    double scuZY = scuZ+u*c, ssuZY = ssuZ+u*s, su2ZY = su2Z+u*u;
#define ForeTrack_UPDATE_YP
#ifdef ForeTrack_UPDATE_YP
	    // New value of horizontal slope, if "pY" is actually a Y plane
	    float ypp = fabs(sY)<.035 ? (yp+(Y-y0)/(xY-x0))/2 : yp;
#else
	    float ypp = yp;
#endif
	    static float xYp, cYp, sYp;
	    pYp = pY->Associate; if (pYp && pYp->IFlag) {
	      if (1<<pYp->IPlane-ipli&usedPlPatY) // Associated plane already
	      // scanned: the combination (hit,associated) has probably
	      continue; // already been checked => skip
	      getAssociatedHit(ypp,zpp,cY,sY,pY,hY,pYp,&hYp);
	    }
	    if (hYp) {
	      const TDetect &dYp = setup.vDetect(pYp->IDetRef);
	      cYp = dYp.Ca; sYp = dYp.Sa; xYp = dYp.X(0);
	      float Yp = hYp->U, whYp = hYp->SigU;
#ifdef ForeTrack_UPDATE_YP
	      if (fabs(sY)<.035) ypp = (yp+(Y-y0)/(xY-x0)+(Yp-y0)/(xYp-x0))/3;
#endif
	      Yp += (xZ-xYp)*(ypp*cYp+zpp*sYp);
	      c = cYp/whYp; s = sYp/whYp; u = Yp/whYp;
	      sc2ZY += c*c; sscZY += s*c; ss2ZY += s*s; 
	      scuZY += u*c; ssuZY += u*s; su2ZY += u*u; nhZY++; hPatZY |= 0x8;
	    }
	    double detZY = sc2ZY*ss2ZY-sscZY*sscZY;
	    if (fabs(detZY)<1.e-5) continue;
	    float yZY = (scuZY*ss2ZY-ssuZY*sscZY)/detZY;
	    float zZY = (ssuZY*sc2ZY-scuZY*sscZY)/detZY;

	    for (ipU = ipY, projU=usedPlPatU = 0; ipU!=st->IPlanes.end();
		 ipU++) {

	      // ***** LOOP ON PLANES FOR ``U'' *****
	      // Start search from "ipY" since all its predecessors failed
	      // to provide a single hit compatible with track's initial

	      const TPlane *pU = &setup.vPlane(*ipU); int jplU = *ipU-ipli;
	      if ((p2Projs[jplU]&sProjPatZY)==p2Projs[jplU] ||
		  // Forbid DR01X (resp. Y) when already DR02X (resp. Y)
		  projU && p2Projs[jplU]!=projU) continue;
	      projU = p2Projs[jplU];
	      if (!pU->IFlag) continue;                  // Plane OFF
	      ipV = ipU;   // Remember last "ipU"
	      if (pU->vHitRef().empty()) continue;       // Plane empty
	      const TDetect &dU = setup.vDetect(pU->IDetRef);
	      float cU = dU.Ca, sU = dU.Sa, xU = dU.X(0), wU = dU.Resol;
	      y = y0+yp*(xU-x0); z = z0+zp*(xU-x0);// Extrapolate from x0: yp,zp
	      if (!dU.InActive(y,z)) continue;             // Out of active zone
	      showerFlags[jplU]++;
	      u = y*cU+z*sU; du = sqrt(dy2*cU*cU+2*dyz*cU*sU+dz2*sU*sU);
	      dup = sqrt(dyp2*cU*cU+2*dyzp*cU*sU+dzp2*sU*sU)*fabs(xZ-xU);
	      float umxU = u+6*du+8*wU, umnU = u-6*du-8*wU, uZY = yZY*cU+zZY*sU;
#ifdef DEBUG
	      if (idebug) {
		if (isMC) {
		  printf("Y %.2f,%d,%d ",hY->U,hY->IKine,hY->IKin2);
		  if (hYp) printf("Y' %.2f,%d,%d ",hYp->U,hYp->IKine,hYp->IKin2);
		}
		printf("u %f du %f => %f %f   yzuZY %f %f %f\n",
		       u,du,umnU,umxU,yZY,zZY,uZY);
	      }
#endif	      
	      for (ihU = pU->vHitRef().begin(); ihU!=pU->vHitRef().end();
		   ihU++) {

		// ***** LOOP ON ``U'' HITS *****

		THit *hU = &vecHit[*ihU], *hUp = NULL;
		if (hU->Status!=0) continue;
		float U = hU->U; if (U<umnU) continue; if (umxU<U) break;
		sProjPatZYU |= projU;
		usedPlPatU |= 1<<jplU;
		if (showerFlags[jplU]==1) {
		  int mirror = hU->Mirror;
		  showerHits[jplU] += mirror>0 && vecHit[mirror].Status==-4 ? 2 : 1;
		}
		double DUMx = 2*wU, DU = hU->SigU<DUMx ? hU->SigU : DUMx;
		if (U<umnU+4*DU) continue; if (umxU-4*DU<U) break;

		// Evaluate space point
		U += (xZ-xU)*(ypp*cU+zpp*sU);
		if (fabs(uZY-U)>3*DU+3*dup) continue;

		// ********** CANDIDATE ZYU SPACE POINT **********

		float whU = hU->SigU; c = cU/whU; s = sU/whU; u = U/whU;
		spt.sc2 = sc2ZY+c*c, spt.ssc = sscZY+s*c; spt.ss2 = ss2ZY+s*s;
		spt.nHits = nhZY+1; spt.hPat = hPatZY | 0x10;
		spt.scu = scuZY+u*c; spt.ssu = ssuZY+u*s, spt.su2 = su2ZY+u*u;
		static float xUp, cUp, sUp;
		pUp = pU->Associate; if (pUp && pUp->IFlag) {
		  if (1<<pUp->IPlane-ipli&usedPlPatU)// Associated plane already
		    // scanned: the combination (hit,associated) has probably
		    continue; // already been checked => skip
		  getAssociatedHit(ypp,zpp,cU,sU,pU,hU,pUp,&hUp);
		}
		if (hUp) {
		  const TDetect &dUp = setup.vDetect(pUp->IDetRef);
		  cUp = dUp.Ca; sUp = dUp.Sa; xUp = dUp.X(0);
		  float Up = hUp->U, whUp = hUp->SigU;
		  Up += (xZ-xUp)*(ypp*cUp+zpp*sUp);
		  c = cUp/whUp; s = sUp/whUp; u = Up/whUp;
		  spt.sc2 += c*c; spt.ssc += s*c; spt.ss2 += s*s;
		  spt.scu += u*c; spt.ssu += u*s; spt.su2 += u*u;
		  spt.nHits++; spt.hPat |= 0x20;
		}
		else if (st->IPlanes.size()>4 && // If redundancy is possible...
			 spt.nHits<4) continue;  // ...require it ...
		// (...This is the case of ST04, where there are 6 layers (in
		// term of TPlane's there are in fact 3 times as many: we
		// assume that for the future RichWall detectors #TPlane's>4
		// will also always ensure #layers>4).
		double det = spt.sc2*spt.ss2-spt.ssc*spt.ssc;
		y = (spt.scu*spt.ss2-spt.ssu*spt.ssc)/det; spt.y = y;
		z = (spt.ssu*spt.sc2-spt.scu*spt.ssc)/det; spt.z = z;
		spt.chi2 = (y*y*spt.sc2+z*z*spt.ss2+spt.su2+
			    2*(y*z*spt.ssc-y*spt.scu-z*spt.ssu))/(spt.nHits-2);
#ifdef DEBUG
		if (idebug) {
		  if (isMC) {
		    printf("U %.2f,%d,%d ",hU->U,hU->IKine,hU->IKin2);
		    if (hUp)
		      printf("U' %.2f,%d,%d ",hUp->U,hUp->IKine,hUp->IKin2);
		  }
		  printf("3Spt %.3f\n",spt.chi2);
		}
#endif	      
		if (spt.chi2>TOpt::dCut[69]) continue;
		THit *hs[] = { hZ, hZp, hY, hYp, hU, hUp };
		for (int i = 0; i<6; i++) spt.hs[i] = hs[i];
		spts.push_back(spt);
	      }  // End loop on U hits
	    }  // End loop on U planes

	    if (ipV!=st->IPlanes.end()) ipV++;
	    for (projV=usedPlPatV = 0; ipV!=st->IPlanes.end(); ipV++) {

	      // ***** LOOP ON PLANES FOR ``V'' *****
	      // (Starting at last ipU, assuming coordinates are not intermixed
	      // in a givent station)

	      const TPlane *pV = &setup.vPlane(*ipV); int jplV = *ipV-ipli;
	      if ((p2Projs[jplV]&sProjPatZYU)==p2Projs[jplV])
		continue; // Forbid DR01X (resp. Y) when already DR02X (resp. Y)
	      if (projV && p2Projs[jplV]!=projV && st->IStation!=iDRStation)
		CsErrLog::msg(elFatal,__FILE__,__LINE__,
			      "%s: Proj=0x%x != Z,YUVProj=0x%x,0x%x,%d,0x%x",
			      setup.vDetect(pV->IDetRef).Name.c_str(),
			      p2Projs[jplV],projZ,projY,projU,projV);
	      projV = p2Projs[jplV];
	      if (!pV->IFlag) continue;                  // Plane OFF
	      if (pV->vHitRef().empty()) continue;       // Plane empty
	      const TDetect &dV = setup.vDetect(pV->IDetRef);
	      float cV = dV.Ca, sV = dV.Sa, xV = dV.X(0), wV = dV.Resol;
	      y = y0+yp*(xV-x0); z = z0+zp*(xV-x0);// Extrapolate from x0: yp,zp
	      if (!dV.InActive(y,z)) continue;             // Out of active zone
	      showerFlags[jplV]++;
#if defined ForeTrack_2002_2004
#  if defined DEBUG || defined FT_DEBUG_RW
	      float sYV = sY*cV-cY*sV, sZV = sZ*cV-cZ*sV;
	      if (fabs(sZV)<.1 || fabs(sYV)<.1) {  // X-check V!=Z && V!=Y
		CsErrLog::msg(elError,__FILE__,__LINE__,
		  "%s: Proj=0x%x != Z,YProj=0x%x,0x%x: sZV,sYV = %f,%f",
		  dV.Name.c_str(),p2Projs[jplV],projZ,projY,sZV,sYV);
		continue;       // ``Y'' should be anything but Z
	      }
#  endif
#endif
	      u = y*cV+z*sV; du = dy2*cV*cV+2*dyz*cV*sV+dz2*sV*sV;
	      dup = sqrt(dyp2*cV*cV+2*dyzp*cV*sV+dzp2*sV*sV)*fabs(xZ-xV);
	      float umxV = u+6*du+10*wV, umnV = u-6*du-10*wV, uZY = yZY*cV+zZY*sV;
#ifdef DEBUG
	      if (idebug) {
		if (isMC) {
		  printf("Y %.2f,%d,%d ",hY->U,hY->IKine,hY->IKin2);
		  if (hYp) printf("Y' %.2f,%d,%d ",hYp->U,hYp->IKine,hYp->IKin2);
		}
		printf("v %f dv %f => %f %f   yzvZY %f %f %f\n",
		       u,du,umnV,umxV,yZY,zZY,uZY);
	      }
#endif	      
	      for (ihV = pV->vHitRef().begin(); ihV!=pV->vHitRef().end();
		   ihV++) {

		// ***** LOOP ON ``V'' HITS *****

		THit *hV = &vecHit[*ihV], *hVp = NULL;
		if (hV->Status!=0) continue;
		float V = hV->U; if (V<umnV) continue; if (umxV<V) break;
		usedPlPatV |= 1<<jplV;		
		if (showerFlags[jplV]==1) {
		  int mirror = hV->Mirror;
		  showerHits[jplV] += mirror>0 && vecHit[mirror].Status==-4 ? 2 : 1;
		}
		double DVMx = 2*wV, DV = hV->SigU<DVMx ? hV->SigU : DVMx;
		if (V<umnV+4*DV) continue; if (umxV-4*DV<V) break;

		// Evaluate space point
		V += (xZ-xV)*(ypp*cV+zpp*sV);
		if (fabs(uZY-V)>3*DV+3*dup) continue;

		// ********** CANDIDATE ZYV SPACE POINT **********

		static float xVp, cVp, sVp, Vp;
		pVp = pV->Associate; if (pVp && pVp->IFlag) {
		  if (1<<pVp->IPlane-ipli&usedPlPatV)// Associated plane already
		    // scanned: the combination (hit,associated) has probably
		    continue; // already been checked => skip
		  getAssociatedHit(ypp,zpp,cV,sV,pV,hV,pVp,&hVp);
		}
		if (hVp) {
		  const TDetect &dVp = setup.vDetect(pVp->IDetRef);
		  cVp = dVp.Ca; sVp = dVp.Sa; xVp = dVp.X(0);
		  Vp = hVp->U; Vp += (xZ-xVp)*(ypp*cVp+zpp*sVp);
		}
		THit *hs[] = { hZ, hZp, hY, hYp, hV, hVp };

		int ispt, match, nspts = spts.size(), ih;
		for (ispt=match = 0; ispt<nspts; ispt++) {
		  TSpacePt &spt3 = spts[ispt];
		  if (spt3.hPat&0xc0) continue;
		  int sameZY; for (ih = 0, sameZY = 1; ih<4; ih++)
		    if ((hPatZY&1<<ih) && spt3.hs[ih]!=hs[ih]) sameZY = 0;
		  if (!sameZY) continue;
		  // Check whether the current ZYV combination...
		  //  ...is equal to the running ZYU combination => give up
		  for (ih = 4; ih<6; ih++)
		    if (spt3.hs[ih]==hV) { match = 2; break; }
		  if (match==2) break;
		  //  ...has proj. in common w/ running ZYU => skip that ZYU
		  if (spt3.hPat&0x10 &&
		      projV==p2Projs[spt3.hs[4]->IPlane-ipli]) continue;
		  float u = spt3.y*cV+spt3.z*sV;
		  if (fabs(u-V)>3*dV.Resol+3*dup) continue;

		  // ***** CANDIDATE ZYUV SPACE POINT *****

		  float whV = hV->SigU; c = cV/whV; s = sV/whV; u = V/whV;
		  spt4.sc2 = spt3.sc2+c*c, spt4.ssc = spt3.ssc+s*c;
		  spt4.ss2 = spt3.ss2+s*s;
		  spt4.nHits = spt3.nHits+1; spt4.hPat = spt3.hPat | 0x40;
		  spt4.scu = spt3.scu+u*c; spt4.ssu = spt3.ssu+u*s;
		  spt4.su2 = spt3.su2+u*u;
		  if (hVp) {
		    float whVp = hVp->SigU; c = cV/whVp; s = sV/whVp; u = Vp/whVp;
		    spt4.sc2 += c*c; spt4.ssc += s*c; spt4.ss2 += s*s;
		    spt4.scu += u*c; spt4.ssu += u*s; spt4.su2 += u*u;
		    spt4.nHits++; spt4.hPat |= 0x80;
		  }
		  double det = spt4.sc2*spt4.ss2-spt4.ssc*spt4.ssc;
		  y = (spt4.scu*spt4.ss2-spt4.ssu*spt4.ssc)/det; spt4.y = y;
		  z = (spt4.ssu*spt4.sc2-spt4.scu*spt4.ssc)/det; spt4.z = z;
		  spt4.chi2 = (y*y*spt4.sc2+z*z*spt4.ss2+spt4.su2+
			    2*(y*z*spt4.ssc-y*spt4.scu-z*spt4.ssu))/
		    (spt4.nHits-2);
#ifdef DEBUG
		  if (idebug) {
		    if (isMC) {
		      printf("V %.2f,%d,%d ",hV->U,hV->IKine,hV->IKin2);
		      if (hVp)
			printf("V' %.2f,%d,%d ",hVp->U,hVp->IKine,hVp->IKin2);
		    }
		    printf("4Spt %.3f\n",spt4.chi2);
		  }
#endif	      
		  if (spt4.chi2>TOpt::dCut[69]) continue;
		  for (ih = 0; ih<6; ih++) spt4.hs[ih] = spt3.hs[ih];
		  spt4.hs[6] = hV; spt4.hs[7] = hVp; match = 1;
		  spts.push_back(spt4);
		}
		if (match==2) continue;
		if (!match) {
		  float whV = hV->SigU; c = cV/whV; s = sV/whV; u = V/whV;
		  spt.sc2 = sc2ZY+c*c, spt.ssc = sscZY+s*c; spt.ss2 = ss2ZY+s*s;
		  spt.nHits = nhZY+1; spt.hPat = hPatZY | 0x10;
		  spt.scu = scuZY+u*c; spt.ssu = ssuZY+u*s, spt.su2 = su2ZY+u*u;
		  if (hVp) {
		    float whVp = hVp->SigU; c = cVp/whVp; s = sVp/whVp; u = Vp/whVp;
		    spt.sc2 += c*c; spt.ssc += s*c; spt.ss2 += s*s;
		    spt.scu += u*c; spt.ssu += u*s; spt.su2 += u*u;
		    spt.nHits++; spt.hPat |= 0x20;
		  }
		  else if (st->IPlanes.size()>4 &&// If redundancy is possible...
			   spt.nHits<4) continue; // ...require it ...
		  double det = spt.sc2*spt.ss2-spt.ssc*spt.ssc;
		  y = (spt.scu*spt.ss2-spt.ssu*spt.ssc)/det; spt.y = y;
		  z = (spt.ssu*spt.sc2-spt.scu*spt.ssc)/det; spt.z = z;
		  spt.chi2 = (y*y*spt.sc2+z*z*spt.ss2+spt.su2+
			      2*(y*z*spt.ssc-y*spt.scu-z*spt.ssu))/(spt.nHits-2);
#ifdef DEBUG
		  if (idebug) {
		    if (isMC) {
		      printf("V %.2f,%d,%d ",hV->U,hV->IKine,hV->IKin2);
		      if (hVp)
			printf("V' %.2f,%d,%d ",hVp->U,hVp->IKine,hVp->IKin2);
		    }
		    printf("3Spt %.3f\n",spt.chi2);
		  }
#endif	      
		  if (spt.chi2>TOpt::dCut[69]) continue;
		  for (ih = 0; ih<6; ih++) spt.hs[ih] = hs[ih];
		  spts.push_back(spt);
		}
	      }  // End loop on V hits
	    }  // End loop on V planes
	  }  // End loop on Y hits
	}  // End loop on Y planes
	if (iplZ>=0) break;
      }  // End loop on Z hits
    }  // End loop on Z planes

    // ********** RETAIN SPACE POINT w/ BEST PERFS in #HITS and CHI2 **********
    if (spts.size()==0) { it++; continue; }
    int ispt, nspts = spts.size(), nEHMx; double chi2Mn; static int ispt0;
    for (ispt=nEHMx = 0, chi2Mn = TOpt::dCut[69]; ispt<nspts; ispt++) {
      TSpacePt &spt = spts[ispt];
      // Determine type of station
      // ...giving a bonus to non-drift (PS01 or GM04) over drift (DR or ST),
      // to compensate for the fewer coordinates, and, as far as PS01 is
      // concerned, because it's more upstream, w/ therefore less
      // multiscattering and more reliable track
      // extrapolation.
      THit *h; int ih; for (ih = 0, h = 0; ih<8; ih++)
	if (spt.hPat&1<<ih) { h = spt.hs[ih]; break; }
      if (!h) continue; // Should never happen. Yet, to be on the safe side...
      int detType = setup.vDetect(h->IPlane).IType, nEHits;
      if (detType<=2/* PS */ || detType==26/* GM04 */) nEHits = spt.nHits;
      else                                             nEHits = (spt.nHits+1)/2;
      if (nEHits>nEHMx || nEHits==nEHMx && spt.chi2<chi2Mn) {
	nEHMx = nEHits; chi2Mn = spt.chi2; ispt0 = ispt;
      }
    }
    if (nEHMx==0) { /* spts.clear(); */ it++; continue; }
    TSpacePt &spt0 = spts[ispt0];
    bool isDRStation =
      setup.vPlane(spt0.hs[0]->IPlane).Station->IStation==iDRStation;
    if (isDRStation) {    // ********** SPECIAL CASE of DR **********
      // DR have no stereo plane => To ensure track reliabilty we...
      int nYUVs, ih; for (ih = 0, nYUVs = 0; ih<8; ih++) if (spt0.hPat&1<<ih) {
	if (1<<spt0.hs[ih]->IPlane-ipli&plYUVPat) nYUVs++;
      }
      if (spt0.nHits<TOpt::iCut[26] ||     // ***** ...REQUIRE MIN. #HITS AND ...
	  nYUVs<2 || spt0.nHits-nYUVs<2) { // ***** ...AT LEAST 2 HITS per COORD
	/* spts.clear(); */ it++; continue;
      }
    }

    //            ********** SHOWER? **********
    int s1, sHY, sHZ, sH, sj, sj2, sjH;
    for (jpl = 0, s1=sHY=sHZ=sj=sj2=sjH = 0; jpl<=iplf-ipli; jpl++) {
      unsigned int bpl = 1<<jpl; if (!(bpl&plDRPat)) continue;
      s1 += 1; sj += jpl; sj2 += jpl*jpl; 
      int nHs = showerHits[jpl]; sjH += jpl*nHs;
      if (bpl&plYUVPat) sHY += nHs;
      else              sHZ += nHs;
    }
    t.NShowerHits = sHZ+sHY;
    if (sHZ>=12 && sHY>=12 && sHZ+sHY>24) {
      t.HasShower = 1; // 1+1 DR01s + 2+2 DR02s (x 2 for mirrors) in Y & Z, + 1
      float dHdj = float(s1*sjH-sj*(sHZ+sHY))/(s1*sj2-sj*sj);
      if (dHdj>.30) t.HasShower = 2;   // 1+1+1+1+2+2+2+2 would be ~= .381
    }
    if (hShower && t.IKsame>=0) {
      int id = vecKine[t.IKsame].PDGid(), ep = abs(id)==11 ? 1 : 0;
      hShower->Fill(ep,t.HasShower);
    }
    if (t.HasShower==2) { /* spts.clear(); */ it++; continue; }

    //          ********** AMBIGUITY? **********
    bool ambiguity; for (ispt = 0, ambiguity = false; ispt<nspts; ispt++) {
      if (ispt==ispt0) continue; TSpacePt &spt = spts[ispt];
      if ((spt.nHits<spt0.nHits || spt.chi2>spt0.chi2+3) &&
	  // Consider only ``peer'' spt's, i.e. same #hits, w/ not too much
	  // worse chi2.
	  // And still allow for some ``inferior'' ones, w/ smaller (by one
	  // unit, and still >=6) #hits, but better chi2.
	  (spt.nHits<6 || spt.nHits<spt0.nHits-1 || spt.chi2>spt0.chi2-.5))
	continue;
      int ih, nDiffs, nYUVs; for (ih=nDiffs=nYUVs = 0; ih<8; ih++) {
	if (!(spt0.hPat&1<<ih) || !(spt.hPat&1<<ih)) continue;
	THit *h = spt.hs[ih], *h0 = spt0.hs[ih];
	if (1<<h->IPlane-ipli&plYUVPat) nYUVs++;
	if (h==h0 || h->Mirror==h0->IHit) continue;
	const TPlane &p = setup.vPlane(h->IPlane);
	const TDetect &d = setup.vDetect(p.IDetRef);
	if (fabs(h->U-h0->U)>3*d.Resol) nDiffs++;
      }
      ambiguity |= nDiffs>2 && (!isDRStation ||// If DR: require 2 hits/coord
			       nYUVs>=2 && spt.nHits-nYUVs>=2);
    }
    if (ambiguity) { /* spts.clear(); */ it++; continue; }

    // ********** ADD SPACE POINT to TRACK and REFIT **********
    double oldChi2 = t.Chi2tot/(t.NDFs-5);
    bool fillHist = hChi2!=NULL; int nGoodHits = 0;
    if (fillHist) {
      t.FindKine(); fillHist &= t.IKine>=0;
    }
    // It's a halo muon accidentally coincident?
    double evtT = eventTRef ? eventTRef : eventTime;
    double tT0; if (t.DispTime>0 &&// This requires 2 time-measuring detectors..
		    // ...avoiding the risk have track's time set by ghost hit.
		    fabs(t.MeanTime-evtT)/t.SigmaTime>2 &&        // Off time...
		    // Note: When time undef: SigmaTime<0
		    fabs(t.Hlast(3))<.05 && fabs(t.Hlast(4))<.05) // Paraxial...
      tT0 = t.MeanTime-eventTRef;// ...Off-time halo muon: Set track's time "tT"
    else tT0 = evtT0;
    for (int ih = 0; ih<8; ih++) if (spt0.hPat&1<<ih) {
      THit *h = spt0.hs[ih];
      if (fillHist && (h->IKine==t.IKine || h->IKin2==t.IKine)) nGoodHits++;
      if ((1<<h->IPlane-ipli)&hitPlPat) continue;
      if (TOpt::ReMode[29]&0x2) {
	//   ***** TIME PROPAGATION, TRIGGER JITTER, TRACK'S TIME OFFSET *****
	const TPlane &p = setup.vPlane(h->IPlane);
	const TDetect &d = setup.vDetect(p.IDetRef);
	//  N.B.: we need only a crude estimate => non updated
	// values "yp" and "zp" are good enough for that purpose.
	double y = (y0+yp*(d.X(0)-x0))*10, z = (z0+zp*(d.X(0)-x0))*10;
	CsDetector *csDet = d.PtrDet(); if (csDet->hasDrift()) {
	  bool error; double U = csDet->getCorrU(h->PtrClus(),y,z,tT0,error);
	  if (!error) {
	    h->u = U/10; h->Rotate();
	    if (TOpt::ReMode[29]&0x4)
	      // ***** UPON OPTION: RESET CsCluster *****
	      // (This is only useful when one intends to make use of
	      // these clusters outside of TraFDic, e.g.:
	      //  - TraFDic monitoring when TraFFiC modularity != 2
	      //  - Evaluation of residuals in coral, as opposed to TraFDiC)
	      h->PtrClus()->setU(U);
	  }
	}
      }
      t.AnnexHit(*h);
    }

    t.InsertMissedPlanes();
    bool ok = t.FullKF(1); if (ok) {
      double chi2 = t.Chi2tot/(t.NDFs-5);
      ok = chi2<TOpt::dCut[17] && chi2<oldChi2+TOpt::dCut[17]/10;
      if (ok) ok = t.FullKF(-1);
      if (fillHist) {
	if (ok) {
	  hChi2 ->Fill(spt0.chi2,   (double)nGoodHits/spt0.nHits);
	  hdChi2->Fill(chi2-oldChi2,(double)nGoodHits/spt0.nHits);
	  hNHits->Fill(sHY+sHZ,     (double)nGoodHits/spt0.nHits);
	  if (t.NHsame+nGoodHits<TOpt::dCut[10]*t.NHits)
	    CsErrLog::msg(elError,__FILE__,__LINE__,
	      "Event #%d track ID=%d,MC=%d,%d/%d downgraded => %d/%d",event,
		 t.Id,t.IKine,t.NHsame,t.NDFs-spt0.nHits,t.NHsame+nGoodHits,t.NDFs);
	}	  
      }
    }
    if (ok) { // ********** FIT OK => UPDATE Status of THit's **********
      for (int i = 0; i<8; i++) if (spt0.hPat&1<<i) {
	THit *h = spt0.hs[i]; h->status = 1;
	int mirror; if ((mirror = h->Mirror)!=-1) vecHit[mirror].status = 1;
      }
    }
    else {    // ********** FIT FAILS => RESTORE INITIAL TRACK **********
      for (int i = 0; i<8; i++) if (spt0.hPat&1<<i) {
	THit *h = spt0.hs[i];
	if ((1<<h->IPlane-ipli)&hitPlPat) continue;
	if (h->IPlane<iplLast) t.CancelHit(h->IPlane);
      }
      t.Shorten(iplLast+1);
      ok = t.FullKF(1);
      if (ok) {
	ok = t.FullKF(-1);
      }
      if (!ok) {
	CsErrLog::msg
	  (elError,__FILE__,__LINE__,
	   "Full KF fails after initial track (ID=%d,chi2=%f) is restored",
	   t.Id,oldChi2);
	listTrack.erase(it++); continue;
      }
    }
#  ifdef DEBUG
    if (idebug) it->DumpHits(idebug);
#  endif
    /* spts.clear(); */ it++;
  } // End loop on tracks
  } // End of 3 iterations
  delete [] trkPats; delete [] showerHits; delete [] showerFlags;
  //                                                ***** RE-ASSESS FISHY TRACKS
  // (Those are suspicious tracks, w/ associated RICHWall hits, which have been
  // freed for the time of the "ForeTrack2RICHWall" processing. Let's now
  // reassess them:
  //  - If one of their RICHWall has been used, we conclude that they are indeed
  //   ghosts and erase them.
  //  - We determine if a hit (concentrating on RICHWall hits) is used by
  //   checking that its "THit::sTrackID()" is >1.)
  list<unsigned int>::iterator jt;
  for (jt = fishyTrks.begin(); jt!=fishyTrks.end(); jt++) {
    unsigned int id = *jt; for (it = listTrack.begin(); it!=listTrack.end();
				it++) {
      if ((*it).Id!=id) continue;
      TTrack &t = *it; list<int>::iterator ihp, ipr;
      for (ihp = t.lHitPat.begin(), ipr = t.lPlnRef.begin();
	   ihp!=t.lHitPat.end(); ihp++, ipr++) {
	if (*ipr<ipli) continue; if (*ipr>iplf) break; if (*ihp<0) continue;
	THit &h = vecHit[*ihp]; if (h.sTrackID().size()>1) {
	  listTrack.erase(it); break;
	}
      }
      break;
    }
  }
}
void TEv::UpdateHitStatus(int ipli, int iplf, unsigned long long plPat,
			  list<unsigned int> &fishyTrks, int &nHits)
{
  // - Reassess "THit::Status":
  //   - Set it =0 for those THit's free to be used, i.e. not added to TTrack
  //    themselves and not mirror of an added THit.
  //   - Set it =0 also for the hits associated to suspicious tracks, and
  //    remember those (->"lFishyTrks"), so as to re-assess them a posteriori.
  //   "Thit::Status" has been set in "TEv::PrePattern2" and have been upset
  //  in "TEv::FitSegments,CleanTrackList".
  //  (Secondarily, some settings have been assigned by "ImportClusters". But
  //  these concern characteristics that are not at stake here: (as of 05/08)
  //  coalescence ("Status==-4" must have prevented THit from being added to
  //  tracks and must have remained unchanged throughtout PR, and we keep it as
  //  is) and LR probability (the info is used only in PR, may have been
  //  overwritten (in PR, when THit is added) and we feel free to overwrite it
  //  again now).
  //   "TStatus" will have to be taken care of throughout the ForeTrack methods
  //  for at the difference to "THit::sTrackID", it is not (and cannot easily)
  //  be updated by "TTrack::AddHit(and related methods)". For simplicity sake,
  //  we will take care of it even for non-drift THit's)
  // - The method also returns the total number of free hits, be they fully free
  //  or just temporarily freed following the above-described criterion. 
  const TSetup &setup = TSetup::Ref(); const unsigned long long b1 = 0x1;
  int ipl; for (ipl = ipli, nHits = 0; ipl<=iplf; ipl++) {
    const TPlane &p = setup.vPlane(ipl); if (!(b1<<ipl-ipli&plPat)) continue;
    const TDetect &d = setup.vDetect(p.IDetRef);
    vector<int>::const_iterator ihit;
    for (ihit = p.vHitRef().begin(); ihit!=p.vHitRef().end(); ihit++) {
      THit &h = vecHit[*ihit]; const set<unsigned int> &tIDs = h.sTrackID();
      int nIDs = tIDs.size(), status; if (!nIDs) status = 0;
      else if (nIDs==1) {// Associated to a single track => track's reliability?
	// (Note that this, single track, case must be the most common, after
	// multi-track associations have been sorted out in "CleanTrackList".)
	unsigned int id = *tIDs.begin(); TTrack *t; list<TTrack>::iterator it;
	for (it = listTrack.begin(), t = 0; it!=listTrack.end(); it++) {
	  if ((*it).Id==id) { t = &(*it); break; }
	}
	status = 1; if (t) {
	  if (t->Type==0x2) {
	    int nSPts, nProjs;
	    t->Evaluate(setup.vIplFirst()[1],setup.vIplLast()[1],3,nSPts,nProjs);
	    if (nSPts<3 ||
		nSPts<4 && t->Chi2tot/(t->NDFs-4)>3 ||
		nSPts<5 && t->Chi2tot/(t->NDFs-4)>20 ||
		nProjs<3) {
	      status = 0;
	      int yet; list<unsigned int>::iterator jt;
	      for (jt = fishyTrks.begin(), yet = 0; jt!=fishyTrks.end(); jt++) {
		if (*jt==id) { yet = 1; break; }
	      }
	      if (!yet) fishyTrks.push_back(id);
	    }
	  }
	}
	else {
	  CsErrLog::msg(elError,__FILE__,__LINE__,
	    "Inconsistency: Event #%d; Hit #%d associated to non-existent Track",
			event,h.IHit);
	  // W/ such a no longer existing association, one would be tempted to
	  // set "status=0"... but given that this is probably a tricky case,
	  // let's not play the odds.
	}
      }
      else status = 1;
      if (status==0 && d.IType==17) { // Special case of ST04 or DR's...
	int mirror; // ...reassess THit's "Status" together that its mirror's
	if (h.Status==-4)  // Exclude silent component of coalesced cluster
	  status = -4;
	else if ((mirror = h.Mirror)!=-1 &&
		 !vecHit[mirror].sTrackID().empty()) status = 1;
      }
      if      (status==0) { h.status = 0; nHits++; }
      else if (status==1)
	h.status = 1;  // Meaning also GEMs spacer hits info erased
    }
  }
}
int TEv::getAssociatedHit(float yp, float zp,
			  float c, float s, const TPlane *p, THit *h, 
			  const TPlane *pp, THit **hp)
{
  CsCluster *ci = h->ptrCl;
  // Associated clusters (from initial LR)
  list<CsCluster*> associates = ci->getAssociateClusters();
  if (associates.empty()) return 0;
  CsCluster *ck = associates.front();
  double incidence = c*yp+s*zp;	  
  CsEventUtils::setLRToMatchAngle(ci,ck,incidence);
  if ((ci->getLRProb())<.49) { // Allow for cluster/mirror equally good (==50%)
    // Selected Z hit found to be rathor a mirror than a genuine hit
    // (w/ proper handling of the case of coalesced THits. In this case,
    // we have one active THit and one dormant (status==-4), but these still
    // have distinct CsClusters => that of active THit may perform worse in LR.)
    int mirror = h->Mirror; if (mirror==-1) return 0;
    if (vecHit[mirror].Status!=-4) return 0;
  }
  CsCluster *cj = ci->getAssociateClusters().front();
  vector<int>::const_iterator ihit;
  for (ihit = pp->vHitRef().begin(); ihit!=pp->vHitRef().end(); ihit++) {
    // Loop over hits
    THit &hj = vecHit[*ihit];
    if (hj.ptrCl==cj) {
      if      (hj.Status==1) return 0;   // Used => exit
      else if (hj.Status==-4) {          // Silent component of coalesced THit
	*hp = &vecHit[hj.Mirror]; return 1;      // => Take mirror instead
      }
      else if (cj->getLRProb()<.49) {     // Ambiguous and...
	int mirror = hj.Mirror;
	if (mirror!=-1 &&
	    vecHit[mirror].Status!=-4)       // ...not coalesced => exit
	  return 0;
	else { *hp = &hj; return 1; }
      }
      else { *hp = &hj; return 1; }
    }
  }
  return 0;
}

class MWAHitList {
public:
  MWAHitList() { hitPat = 0; nhs = 0; qual = 0; }
  THit *hL[6] /* There may be 2 associated HG02 hits */; int nhs;
  unsigned long long hitPat;
  double qual;  //!< Quality function, based on reduced chi2 of extension, and some bonus
  //! "Less" operator
  bool operator < (MWAHitList &hL) const { return qual>hL.qual; } 
  int mark;
  void Print() {
    printf("%d %6.3f ",nhs,qual);
    for (int i = 0; i<nhs; i++) printf(" %3d %4d",hL[i]->IPlane,hL[i]->IHit);
    printf("  %16Lx\n",hitPat);
  }
};
class Track2MWA {
public:
  list<TTrack>::iterator it;   //!< Reco'd spectrometer track
  int nPA01s;  //!< #PA01 hits in "t"
  unsigned long long hitPlPat;
  TTrack tmwa; //!< Extension into muWallA
  double qual;  //!< Quality function, based on reduced chi2 of extension, and some bonus
  //! "Less" operator
  bool operator < (const Track2MWA &tp) const { return qual<tp.qual; } 
};
void findMWA(vector<THit*> *hSs,
	     const THlx &Hinge, const THlx &H,
	     const TTrack &t,
	     int iZY, int nSigmasF,
	     list<MWAHitList> *mwaHLs);
void TEv::ForeTrack2MA()
{
  // ********** EXTEND TRACKS TO MA DETECTORS (and PA01|2/HG02) **********

  // - PR often fails to extend 0x2 track segments into muWallA detectors.
  //  The present method aims at achieving (or finalising) this extension.
  // - It is performed on tracks already
  //   - (SM1-)bridged,
  //   - KF'ed.
  //   A 3rd condition concerns the RichWall (understood as all LAT detectors in
  //  the sub-zone right after the RICH): the candidate track is required to be
  //  reco'd down to there, if it falls into (the heart of) its acceptance.
  //  (Note still that option ReMode[49] (so-called "DrellYan") allows to bypass
  //  that condition.)
  // - The algorithm is extrapolation and pick-up.
  //   It's in two steps:
  //   I) Extension step.
  //  II) Disambiguation step, selecting best among several extensions sharing
  //     the same picked up hits, based on:
  //     - reduced chi2 AND #hits,
  //     - w/ a bonus given to extensions including hits from HG02, since these
  //       convey some y vs. z correlation info.
  //    (Notes:
  //     - In order to get the best out of HG02, one would like to have their
  //      hit coordinate along y as derived from Jura vs. Saleve TDCs.
  //     - Split MA channels convey also some y vs. z correlation info => could
  //      also be considered for a bonus.)
  // - Hits associated to unreliable tracks (cf. "fishyTrks") are freed for the
  //  time of the processing, as is done in "ForeTrack2RICHWall".

  // N.B.:
  // - Some of the setup info needed by the method is built-in.
  // - This method may erase some tracks.
  // - There three tunable parameters:
  //   - Road width = "dCut[99]",
  //   - Intermediate chi2 cut used to limit the combinatorics = "dCut[65",
  //   - Final chi2 cut, expressed as max. allowed chi2 increment per muWallA
  //    hit = "dCut[94]".
  // - The present algorithm is based on a simplistic extrapolate+pick-up
  //  scheme, not so well adapated to the muWallA context w/ its high material
  //  density and consequent multi-scattering. It's been improved over time by
  //  ad-hoc corrections. But something more innovative would still be wellcome:
  //  e.g. a Kalman-Filter-based hit search.

  const TSetup &setup = TSetup::Ref();
  static int ipli, iplf, first = 1, npls; int ipl, jpl;
  // Bit-Patterns of detector planes of interest.
  static unsigned long long plPat1, plYPat, plPat2, plPatHG, plPatPA1, plPatPA2, plPatAll;
  const unsigned long long b1 = 0x1;
  // Acceptance cuts
  static float dzYMnS, dzYMnJ, dzZMnB, dzZMnT, xMn, accYMx, xYMx, accZMx, xZMx;
  float xMA01 = 0;
  // Acceptance of ''heart of RICHWall''
  static float xPS, accPS[2];
  static float accGM = 15; // Default supplied, since the availibility of GM04 is not checked
  static float xVLAT, dzVLAT[2];
  // PA01|2s
  static const TStation *sPA01, *sPA02;
  static TH2D *hChi2, *hdChi2;
  float chi2IperHit = // Max. allowed chi2 increment for picked-up hits.
    // The setting is preferably loose: in order to provide for:
    // - The imprecision in the description of the distribution of dead
    //  material, that may lead to underestimate the multiscattering affecting
    //  the extrapolation.
    // - The high rate of delta-rays expected in this dense MA environment:
    //   - A hit from delta-rays can outperform genuine hits in the
    //    extrapolate+pick-up algorithm infra, if it happens to sit, while
    //    somewhat offset w.r.t. the particle's true trajectory, in direct
    //    continuation of the original track (reco'd upstream of the dense
    //    material of E/HCal1).
    //   - A loose setting may then accommodate some of those delta-ray hits,
    //    somewhat unduly, but w/o any large impact on the track parameters.
    //    E.g. event #3564 of 2006.03 MC D*'s.
    // Nota Bene: this derivation of the cut value from option "dCut[94]" is
    // also done, independently, in "findMWPAHits".
    TOpt::dCut[61] ? TOpt::dCut[61] : 9;
  int DrellYan =  // Special option for DrellYan
    // It bypasses several of the reconstruction requirements, in particular:
    // - The early discarding of tracks estimated to fall into the DZ of MAs
    //  and HG02 based on straight line extrapolation is bypassed.
    // - The discarding of those that fail to be reco'd in the RICHWall (in the
    //  broad sense, cf. "TEv::ForeTrack2RICHWall"), even though they fall into
    //  the acceptance of its heart (i.e. the highly redundant region of PS01).
    // In addition, when setting >1:
    // - If candidate track turns out to be completely out of the acceptance of
    //  the MA01s, these are not required, and this, even if no PA01 triple.
    // - It also enlarges the road width for the MA02s.
    TOpt::ReMode[49];
  double roadWidth = TOpt::dCut[99] ? TOpt::dCut[99] : 4;
  int nSigmasMA02Y = // Special option to loosen the hit search in MA02Y
    // - The default setting is #sigmas = 3.
    // - This may not be enough because of the structure of the MAs (w/ their
    //  sub-units separated by voids where efficiency is null) and because of
    //  their at times, poor quality of the MAs.
    // - The effect being, presumably, more pronounced...
    //  ... in the vertical, because the track uncertainty is smaller in that
    //    dimension, and then cannot accommodate the intrinsic deficiencies of
    //    the detectors.
    //  ... in station MA02 (as opposed to MA01) because of the uncertainty
    //    accumulated along the way.
    // - The whole thing is on shaky ground.
    //   - Yet some cases show that relaxing the #sigmas to =5 can help, e.g.:
    //     ~/w0/evtDump/2015/evtDump.bad.MAreco-258130.raw
    //   - I have no time to investigate further and want to limit the changes
    //    brought to the source and minimize the amount of work needed to
    //    evaluate its impact.
    // => Therefore, this minimalist option that only loosens the hits search
    //   in MA02Y.
    TOpt::ReMode[50] ? 5 : 3; 

  if (first) { 
    //     **************************************************
    //     ***************   INITIALISATION   ***************
    //     **************************************************

    // ***** DETERMINE WHICH PLANES are MA, PA01|2 or HG02 DETECTORS *****
    //  => PLANE# RANGE [ipli,iplf] and PLANE# PATTERNs
    //  => ACCEPTANCE EXTREMA
    //     ACCEPTANCE at the EXIT of RICH
    //  => REMEMBER THE ACCEPTANCE of PS01
    //   (In order to save CPU, we want to restrict ourselves to tracks w/ good
    //   prospect of success. One idea is to require reco downstream of RICH:
    //   tracks that do not make it there are most probably low momentum ones or
    //   unreliable, and anyway their extrapolation into the MA zone is bound to
    //   be imprecise. On the other hand we would not like to heavily rely on
    //   the RICHWall (in the broad sense, cf. "TEv::ForeTrack2RICHWall") reco.
    //   => Let's disregard only those tracks that fall into the heart of the
    //    RICHWall system, viz. PS01, and this when there is possible overlap w/
    //    VLAT (ST04/DR) and hence reasonable or high redundancy is guaranteed,
    //    and still fail to be reco'd.)
    first = 0; int nYZs012[4] = {0,0,0,0}; float aYMn, aYMx, aZMx;
    for (ipl = setup.vIplFirst()[1], ipli=iplf = -1, aYMn = 3, aYMx=aZMx = 0,
	   npls=plYPat=plPat1=plPat2=plPatHG=plPatPA1=plPatPA2 = 0,
	   xPS=xVLAT = 0, sPA01=sPA02 = 0; ipl<=setup.vIplLast()[1]; ipl++) {
      const TPlane &p = setup.vPlane(ipl); if (!p.IFlag) continue;  // Plane OFF
      const TDetect &d = setup.vDetect(p.IDetRef);
      // Determine the acceptance of the heart of the RICHWall system. This means
      // the overlap of PS01 w/ either the VLAT (ST04 or DR) or GM04, which is
      // not part of the RICHWall system but plays a similar role.
      // In all instances, let's simplify: outside its heart, the RICHWall reco
      // is still quite efficient: therefore if we overestimate its extension, it
      // is not that dramatic.
      bool isVLAT = d.Name.find("DR")==0 || d.Name.find("ST04")==0;
      bool isGM04 = d.Name.find("GM04")==0, isPS01 = d.Name.find("PS01")==0;
      bool isPA01 = d.Name.find("PA01")==0, isPA02 = d.Name.find("PA02")==0;
      if (isVLAT || isGM04 || isPS01 || isPA01 || isPA02) {
	int iYZ = p.IProj;
	if (isPA01 || isPA02) {
	  const TStation **s_p = isPA01 ? &sPA01 : &sPA02; *s_p = p.Station;
	}
	else {
	  if (iYZ<0 || 1<iYZ) // Let's not consider stereo planes (for simplicity)
	    continue;
	  if (isPS01) {
	    accPS[iYZ] = d.Range/2; // Retain U ranges of PS01X/Y
	    if (!xPS) // First PS01X|Y encountered: provide for the other one missing
	      accPS[1-iYZ] = d.Siz(2-iYZ);
	    xPS = d.X(0);
	  }
	  else if (isGM04)// In this region of larger momenta, reco is easier... 
	    accGM = d.Range/2; // ... => let's simplify even further.
	  else { // I.e. VLAT of RICHWall
	    float dzDim = d.DZSize(iYZ);
	    if (!dzDim) // "DZ" can be null, e.g. for an outer ST04 slice
	      continue;
	    if (!xVLAT /* not assigned yet */ || dzDim>dzVLAT[iYZ]) {
	      // Retain the maximum DZ size
	      dzVLAT[iYZ] = dzDim; xVLAT = d.X(0);
	    }
	  }
	  continue;
	} // End case of GM04 and RICHWall's PS01 and VLAT
      }
      bool isHG02 = !d.Name.find("HG02");
      if ((d.IType<39 || 40<d.IType) && !isHG02 && !isPA01 &&!isPA02)
	continue; // Retain MA, PA01|2 and HG02
      if (ipli==-1) ipli = ipl; iplf = ipl;
      unsigned long long bpl = b1<<ipl-ipli;
      if      (isHG02) { plPatHG |= bpl; npls++; continue; }
      else if (isPA01 || isPA02) {
	if (p.IProj==0) plYPat |= bpl;
	if (isPA01) plPatPA1 |= bpl;
	else        plPatPA2 |= bpl;
	npls++; continue;
      }
      else {
	npls++;
	if (p.IProj==0) { // Y proj. is 0th, cf. TSetup::Init
	  plYPat |= bpl;
	  if (d.IType==39) { plPat1 |= bpl; nYZs012[0]++; }
	  else             { plPat2 |= bpl; nYZs012[1]++; }
	}
	else {
	  if (d.IType==39) { plPat1 |= bpl; nYZs012[2]++; }
	  else             { plPat2 |= bpl; nYZs012[3]++; }
	}
      }
      // ***** ACCEPTANCE EXTREMA
      //  The purpose is to discard right away tracks obviously off topic, by
      // checking their straight line extrapolation against the extrema. This,
      // in order to save CPU. Given that we will eventually require(*) hits
      // from MA02s, retaining the extrema (in term of angular acceptance) is
      // very little demanding, and will (most probably) make us loose very few
      // genuine tracks (even though multi-scattering effects are overlooked).
      // (*) Note that the DrellYan option works around this early discarding.
      //  Let's consider only MAs since HG02, in 2 pieces, is tricky to sort out.
      float dzYDim = d.DZSize(0), X = d.X(0), a = dzYDim/X;
      if (xMA01==0) xMA01 = X;
      if (dzYDim && a<aYMn) { // "DZ" can be null, e.g. for an outer ST04 slice
	float Y = d.X(1), Z = d.X(2), dzZDim = d.DZSize(1);
	aYMn = a; dzYMnS = Y-dzYDim; dzYMnJ = Y+dzYDim; xMn = X;
	dzZMnB = Z-dzZDim; dzZMnT = Z+dzZDim;
      }
      float YSize = d.Siz(1), ZSize = d.Siz(2), aY = YSize/X, aZ = ZSize/X;
      if (aY>aYMx) {           // Largest Y acceptance
	aYMx = aY; accYMx = YSize; xYMx = X;
      }
      if (aZ>aZMx) {   // Z can be distinct: think of ST04 slices
	aZMx = aZ; accZMx = ZSize; xZMx = X;
      }
    }
    plPatAll = plPat1|plPat2|plPatHG|plPatPA1|plPatPA2;
    if (!plPat1) CsErrLog::mes(elFatal,"No MA01 detectors in the setup!");
    if (!plPat2) CsErrLog::mes(elFatal,"No MA02 detectors in the setup!");
    char YZ[] = "YZ";
    for (int i = 0; i<4; i++) if (nYZs012[i]>4)
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
		    "Too many (=%d>4) MA0%d%cs",nYZs012[i],i%2+1,i<2?'X':'Y');
    if (npls>24)
      // The code is dimensioned for 3*(PA01+2)+2*4*(MAY+Z)+2*HG02.
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
	"Too many (=%d) detector planes in MA or HG02",npls);
    if (iplf-ipli>=(int)sizeof(unsigned long long)*8)
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
	"Too many (=%d) detector planes in PA01...MA...SM2 zone",iplf-ipli+1);
    //  PS01X|Y and PA01|2X|Y play a definite, albeit minor, role in
    // ForeTrack2MA algorithm. If some of these planes may be DetNameOff'd, at
    // least one plane of each of the pairs should be there. => The source code
    // is hence written assuming it's the case. => Have to require them.
    if (xPS==0)
      CsErrLog::mes(elFatal,"None of PS01X|Y detector planes available!");
    if (!sPA01)
      CsErrLog::mes(elFatal,"None of PA01X|Y detector planes available!");
    if (!sPA02)
      CsErrLog::mes(elFatal,"None of PA02X|Y detector planes available!");

    if (TOpt::Hist[7]) {  // ********** BOOK HISTOS **********
      CsHistograms::SetCurrentPath("/Traffic/RDmonitor");
      const char woHG[] = "MA", wHG[] = "MA/HG02";
      char hT[] = "ForeTrack2MA: #Delta#chi^{2} as a f(#MA/HG02-hits)  ";
      sprintf(hT,"ForeTrack2MA: #chi^{2} as a f(#%s-hits)",
	      plPatHG ? wHG : woHG);
      hChi2  = new TH2D("hMWAChi2",hT,
			10,0,TOpt::dCut[17],25,-.5,24.5);
      sprintf(hT,"ForeTrack2MA: #Delta#chi^{2} as a f(#%s-hits)",
	      plPatHG ? wHG : woHG);
      hdChi2 = new TH2D("hMWAdChi2",hT,
			10,-TOpt::dCut[17]/5,TOpt::dCut[17]/5,25,-.5,24.5);
      CsHistograms::SetCurrentPath("/");
    }
  } // End of initialisation

  //       *******************************************************
  //       *************** PRELIMINARY EVENT STEPS ***************
  //       *******************************************************

  // - Reassess "THit::Status":
  // - Determine #hits in all of RichWall => Give up if too many hits.
  list<unsigned int> fishyTrks; int nHitsMuWallA;
  UpdateHitStatus(ipli,iplf,plPatAll,fishyTrks,nHitsMuWallA);
  if (nHitsMuWallA>npls*40) {   // Just not to waste CPU on pathologic events
      CsErrLog::msg(elError,__FILE__,__LINE__,
		    "Too many hits (%d) in PA012s, MAs and HG02s: give up",nHitsMuWallA);
    return;
  }

#ifdef DEBUG
  static int idebug = 0;
#endif

  list<Track2MWA> lTrack2MWA;
  list<TTrack>::iterator it = listTrack.begin(); while (it!=listTrack.end()) {

    //     **************************************************
    //     ***************   LOOP on TRACKS   ***************
    //     **************************************************

    TTrack &t = *it;
    if ((t.Type&0x3)!=0x3)    // ***** CONSIDER ONLY SM1-BRIDGED TRACKS *****
      // (Allowing for tracks bridged over SM2, having in view those muons that
      // pass through muWallA and SM2 yoke (so-called yoke tracks). These could
      // be assigned a large %X0 in any case, even if they are reco'd w/o any MA
      // hits, and therefore would not have to rely on "ForeTrack2MA" to be
      // later muID'd. But this is not assured (as of 2013/10).)
      { it++; continue; }
    if (!(t.IFit&0x60)) { it++; continue; }            // ***** FULL KF REQUIRED

    float x00, y00, z00, yp0, zp0;            // ***** GET Y/Z INTERCEPTS/SLOPES
    if (t.Type&0x4) {
      THlx hlx; t.GetSmoothed(hlx,setup.vDetect(ipli).X(0));
      x00 = hlx(0);
      y00 = hlx(1);     z00 = hlx(2);     yp0 = hlx(3);     zp0 = hlx(4);
    }
    else { 
      x00 = t.Hlast(0);
      y00 = t.Hlast(1), z00 = t.Hlast(2); yp0 = t.Hlast(3), zp0 = t.Hlast(4);
    }
    //      ***** OFF-TOPIC TRACKS? INTERCEPTS/SLOPES against ACCEPTANCE EXTREMA
    float yMn = y00+(xMn-x00)*yp0, zMn = z00+(xMn-x00)*zp0;
    if (fabs(y00+(xYMx-x00)*yp0)>accYMx || fabs(z00+(xZMx-x00)*zp0)>accZMx ||
	!DrellYan &&
	// We include 1cm tolerance, to make for multi-scattering
	dzYMnS+1<yMn && yMn<dzYMnJ-1 && dzZMnB+1<zMn && zMn<dzZMnT-1) {
      it++; continue; }
    float yMA01 = y00+(xMA01-x00)*yp0, zMA01 = z00+(xMA01-x00)*zp0;
    bool fringy = dzYMnS-4<yMA01 && yMA01<dzYMnJ+4 &&
      dzZMnB-4<zMA01 && zMA01<dzZMnT+4;
    //                            ***** TRACKS W/O ANY GOOD PROSPECT OF SUCCESS?
    if (!DrellYan) {
      double pinv = fabs(t.Pinv());
      if (pinv>.4) { it++; continue; }      // ***** MOMENTUM < 2.5 GeV
      if (pinv>.05 &&                       // ***** MOMENTUM < 20 GeV...
	  // High P at large angle (cf. acceptance check supra) are particularly
	  // interesting in this, muID, context: have a lower probability to be a
	  // pion decay. In addition, they are not that many.
	  //  => let's consider them anyway
	  t.Hlast(0)<750) {
	//  ***** ...AND W/IN HEART of RICHWall COMPLEX? REQUIRE RECO BEYOND RICH
	// (This ''heart'' corresponds to the overlap of PS01 and whatever VLAT
	// is in the setup, and GM04. Cf. explanations supra.)
	float yPS = y00+(xPS-x00)*yp0, zPS = z00+(xPS-x00)*zp0;
	bool winGM = /* approx. */ fabs(yPS)<accGM && fabs(zPS)<accGM;
	bool winPS = fabs(yPS)<accPS[0] && fabs(zPS)<accPS[1];
	bool winVLAT; if (xVLAT) {
	  float yVLAT = y00+(xVLAT-x00)*yp0, zVLAT = z00+(xVLAT-x00)*zp0;
	  winVLAT = dzVLAT[0]<fabs(yVLAT) || dzVLAT[1]<fabs(zVLAT);
	}
	else winVLAT = false;
	if (winGM || winPS && winVLAT) { it++; continue; }
      }
    }
    //   ********** INITIAL HITS PATTERN. # HITS in MuWallA?... **********
    // - "hitPlPat": Initially MA+HG, later to include also PAXs found in the
    //              hit search infra. Caveat: PAU|V not included!
    // - "hitPlPAt": Dedicated to initial PAXs.
    // - PAs are distinct in that they do not fully belong to muWallA: only part
    //  of their sensitive area overlaps the hadron absorbers. => Their hits may
    //  be connected to tracks avoiding completely these. Tracks w/ which we do
    //  want to interfere here. E.g. let's not erase those PA hits initially
    //  associated, at variance w/ want we do w/ initial MA hits, when these are
    //  not fulfilling muID requirements.
    unsigned long long hitPlPat, hitPlPAt;
    THit *hsPA01[3], *hsPA02[3];
    int nMA01s, nMA02s, nHG02s, nPA01s, nPA02s, iplLast = t.lPlnRef.back();
    list<int>::iterator ihp, ipr;
    for (ihp = t.lHitPat.begin(), ipr = t.lPlnRef.begin(),
	   nMA01s=nMA02s=nHG02s=nPA01s=nPA02s = 0, hitPlPat=hitPlPAt = 0;
	 ihp!=t.lHitPat.end(); ihp++, ipr++) {
      if (*ihp<0) continue;
      if (ipli<=*ipr && *ipr<=iplf) {
	unsigned long long bpl = b1<<*ipr-ipli;
	if (bpl&plPatAll) {
	  if (bpl&(plPatPA1|plPatPA2)) {
	    hitPlPAt |= bpl;
	    if (bpl&plPatPA1) hsPA01[nPA01s++] = &vecHit[*ihp];
	    else              hsPA02[nPA02s++] = &vecHit[*ihp];
	  }
	  else {
	    hitPlPat |= bpl;
	    if      (bpl&plPat1) nMA01s++;
	    else if (bpl&plPat2) nMA02s++;
	    else                 nHG02s++;
	  }
	}
      }
    }

    if (nMA01s>=5 && nMA02s>=5 ||
    	nPA01s>=3 && nMA02s>=6) { // Cf. "event/PIDdoMuonID.cc#muIDinMW1"
      // ***** TRACK ALREADY FULFILLING MW1(=muWallA)muID REQUIREMENTS? *****
      // - These requirements, when fulfilled, mean a pretty good job was done
      //  by the PR.
      // - IF NO HG02 IN TSetup OR A HG HIT IS ALSO THERE: WE ARE SATISFIED...
      // - Else: may simply be that HG02 was not efficient, but it can also be
      //  that the PR was not good enough. => then, let's try again.
      if (!plPatHG || nHG02s) { it++; continue; }
    }
    if (hitPlPat) { // ***** ...ELSE REMOVE ALL MA + HG02 HITS AND REFIT *****
      // We remove MA hits so that:
      // - track parameters be no longer biased by them,
      // - MW1 muID be no longer granted.
      // (But let's not do the same w/ PA, for simplicity's sake...)
      for (ipr = t.lPlnRef.begin(); ipr!=t.lPlnRef.end(); ipr++) {
	int ipl = *ipr; if (ipl<ipli) continue; if (ipl>iplf) break;
	unsigned long long bpl = b1<<ipl-ipli;
	if (bpl&hitPlPat) {
	  int ihit = t.CancelHit(ipl); vecHit[ihit].status = 0;
	}
      }
      t.Clip(true); iplLast = t.lPlnRef.back(); hitPlPat = 0;
      bool ok = t.FullKF(1); if (ok) ok = t.FullKF(-1);
      if (!ok) {
	CsErrLog::msg(elError,__FILE__,__LINE__,
	  "Event #%d: Full KF fails when track ID=#%d is stripped of its MAs",
		      event,t.Id);
	listTrack.erase(it++); continue;
      }
    }
#ifdef DEBUG
    if (idebug) it->DumpHits(idebug);
#endif

    // *************** TRACK CANDIDATE for EXTENSION ***************

    THit *hZs01[4], *hYs01[5]/* 4+PA01Y */, *hZs02[6]/* 4+2*HG02 */, *hYs02[5];
    THit *hzs01[4], *hys01[5],              *hzs02[6],               *hys02[5];

    //          ***********************************************
    //          **********     MA01s (and PA01s)     **********
    //          ***********************************************

    double xHC1 =  TOpt::Calo[0]; THlx HHC1, HMA1;
    static float x0, y0, z0, yp, zp, y, z, u, dy2, dz2, dyz, du, umx, umn;
    vector<int>::const_iterator ih; int firstPl, kpl, nExpectsY01, nExpectsZ01;
    vector<THit*> hSZs[6]; int nZs01; // Hit sets in Z01
    for (ipl = ipli, firstPl = 1, nZs01=nExpectsZ01 = 0; ipl<=iplf; ipl++) {
      unsigned long long bpl = b1<<ipl-ipli;
      if ((bpl&plPat1)==0) continue;          // Skip because not a MA01, or OFF
      if (bpl&plYPat) continue;                      // Searching Z planes first

      // ***** LOOP OVER Z PLANES in MA01s *****

      const TPlane &pZ = setup.vPlane(ipl);
      const TDetect &dZ = setup.vDetect(pZ.IDetRef);
      float sZ = dZ.Sa, cZ = dZ.Ca, xZ = dZ.X(0), wZ = dZ.Resol;
      if (firstPl) {   // ***** UPON FIRST ENCOUNTERED ACTIVE PLANE...
	// ...let's compute the track uncertainty. But while doing this, let's
	// restrict the recourse to THlx::Extrapolate in order to save CPU.
	// - Extrapolate once and only once into the MA01 region, taking into
	//  account the, possibly thick due to e.g. E/HCal1, material traversed.
	//  And do this at the first active plane encountered, to provide for
	//  particles entering our region from the beam region, through possibly
	//  interposed frames or support structures.
	// - Then, when the material to be traverse comes down to the thinner
	//  windows of the active areas, do not update the uncertainty.
	if (t.Type&0x4) { 	         // ***** COMPUTE UNCERTAINTIES for MA01
	  t.GetSmoothed(HHC1,xHC1);
	  t.GetSmoothed(HMA1,xZ);
	}
	else {
	  HHC1(0) = xHC1; t.Hlast.Extrapolate(HHC1,true);
	  HMA1(0) = xZ; HHC1.Extrapolate(HMA1,true);
	}
	x0 = HMA1(0); y0 = HMA1(1); z0 = HMA1(2); yp = HMA1(3); zp = HMA1(4);
	y = y0; z = z0;
	dy2 = HMA1(1,1); dz2 = HMA1(2,2); dyz = HMA1(1,2);
      }
      else {
	y = y0+yp*(xZ-x0); z = z0+zp*(xZ-x0);
      }
      if (fringy) {
	// We include a 3 sigma tolerance, to make for multi-scattering
	double ydy = y+(y>0?3:-3)*sqrt(dy2), zdz = z+(z>0?3:-3)*sqrt(dz2);
	if (!dZ.InActive(ydy,zdz)) continue;              // Out of Z active zone
      }
      else if (!dZ.InActive(y,z)) continue;
      firstPl = 0;
      u = y*cZ+z*sZ; du = sqrt(dy2*cZ*cZ+2*dyz*cZ*sZ+dz2*sZ*sZ);
      umx = u+roadWidth*du+4*wZ; umn = u-roadWidth*du-4*wZ;
      for (ih = pZ.vHitRef().begin(); ih!=pZ.vHitRef().end(); ih++) {
	THit *h = &vecHit[*ih]; if (h->Status!=0) continue;
	if (h->U<umn) continue; if (umx<h->U) break;
	double flag = h->ptrCl->getDigitsList().front()->getData()[1];
	if (flag &&                // Non null flag means split channel...
	    flag*y<0) continue;    // ...=> Require it to be on right side
	hSZs[nExpectsZ01].push_back(h);
      }
      if (!hSZs[nExpectsZ01].empty()) nZs01++; nExpectsZ01++;
    }
    if (nZs01>64) {  // Limitation necesitated by "findMWA" method
      it++; continue;
    }
    if (nZs01==0 && nPA01s!=3 &&           // ***** NO Z HITS: GIVE UP UNLESS...
	!((DrellYan&0x2) && nExpectsZ01!=0)) { // ...Special DrellYan &= 0x2 case
      it++; continue;
    }

    vector<THit*> hSYs[6]; int nYs01;  // Hit sets in Y01
    for (ipl = ipli, nYs01=nExpectsY01=kpl = 0; ipl<=iplf; ipl++) {
      unsigned long long bpl = b1<<ipl-ipli;
      if (!(bpl&(plPat1|plPatPA1))) continue; // Skip because not PA/MA01 or OFF
      if (!(bpl&plYPat)) continue;                      // Searching Y planes 

      // ***** LOOP OVER Y PLANES in MA01s *****

      const TPlane &pY = setup.vPlane(ipl);
      const TDetect &dY = setup.vDetect(pY.IDetRef);
      float sY = dY.Sa, cY = dY.Ca, xY = dY.X(0), wY = dY.Resol;
      y = y0+yp*(xY-x0); z = z0+zp*(xY-x0);
      // Uncertainties for MA01 already determined since nZs01!=0
      if (fringy) {
	// Something was found (since nZs01!=0), let's further relax the
	// acceptance tolerance to 4 sigmas.
	double ydy = y+(y>0?3:-3)*sqrt(dy2), zdz = z+(z>0?3:-3)*sqrt(dz2);
	if (!dY.InActive(ydy,zdz)) continue;             // Out of Z active zone
      }
      else if (!dY.InActive(y,z)) continue;
      if (firstPl) {// Exceptional case where no active encountered in Z search and processing has nevertheless been continued
	if (t.Type&0x4) {
	  t.GetSmoothed(HHC1,xHC1);
	  t.GetSmoothed(HMA1,xY);
	}
	else {
	  HHC1(0) = xHC1; t.Hlast.Extrapolate(HHC1,true);
	  HMA1(0) = xY; HHC1.Extrapolate(HMA1,true);
	}
	x0 = HMA1(0); y0 = HMA1(1); z0 = HMA1(2); yp = HMA1(3); zp = HMA1(4);
	y = y0; z = z0;
	dy2 = HMA1(1,1); dz2 = HMA1(2,2); dyz = HMA1(1,2);
      }
      nExpectsY01++;              // Now that it's acounted for in "nExpects"...
      if (bpl&hitPlPAt) continue; // ...skip already associated PA01Y
      u = y*cY+z*sY; du = sqrt(dy2*cY*cY+2*dyz*cY*sY+dz2*sY*sY);
      umx = u+roadWidth*du+4*wY; umn = u-roadWidth*du-4*wY;
      for (ih = pY.vHitRef().begin(); ih!=pY.vHitRef().end(); ih++) {
	THit *h = &vecHit[*ih]; if (h->Status!=0) continue;
	if (h->U<umn) continue; if (umx<h->U) break;
	// for PAs channels are not split, only get channel position for MA01
	double flag = (bpl&plPat1) ? h->ptrCl->getDigitsList().front()->getData()[1] : 0;
	if (flag &&                // Non null flag means split channel...
	    flag*z<0) continue;    // ...=> Require it to be on right side
	hSYs[kpl].push_back(h);
      }
      if (!hSYs[kpl].empty()) nYs01++; kpl++;
    }
    if (nYs01>64) {  // Limitation necesitated by "findMWA" method
      it++; continue;
    }

    //                   ********** REQUIRE MIN. # of, ZY 01-HITS **********
    //                                             ***** REQUIRE MIN. EFFICIENCY
    // I.e. min. found/expected, so as not waste amy more CPU on hopeless cases.
    int nExpects = nExpectsZ01+nExpectsY01, ns01 = nZs01+nYs01;
    if (!(DrellYan&0x2) || nExpects!=0) { // Still exclude special DrellYan null acceptance case
      if (!fringy && ns01+nPA01s<nExpects/2) { it++; continue; }
      //                               ***** REQUIRE ENOUGH FOR THE UPCOMING FIT
      if (nZs01<2 || nYs01<2) { it++; continue; } // Note: is this really mandatory?
    }
    // Also, prepare for the criterion on the total # of 01-hits, including pa01
    // to be used later on.
    int nRequired01;
    if (nPA01s==3 ||  // If enough reliability provided by triple PA01
	(DrellYan&0x2) && nExpects==0) // ...or special Drell-Yan case
      nRequired01 = 0;
    else {
      nRequired01 = TOpt::iCut[6]/2-1;
      if (fringy || nExpects<8) // If on the brim of muWallA acceptance
	// Making for detectors crossed outside active area, in order to ease
	// the reco of tracks on the brim of acceptance. Since we are mostly
	// interested in tracks issued from the target (viz. high-Q2 scattered
	// mu's), this is only done for the 1st group of MAs.
      nRequired01 -= 1;
    }

    // ********** BUILD MWA-TRACKS out of MA01s (and possibly PA01s) **********
    list<Track2MWA> t2mwas; if (ns01) {
      //#define FT_DEBUG_MA 2
#ifdef FT_DEBUG_MA
      printf("\n===== findMWA: Track #%d\n",t.GetId());
#endif
      list<MWAHitList> hLZs01, hLYs01;
      int iZY; vector<THit*> *hSs; list<MWAHitList> *hLs;
      for (iZY = 0, hSs = hSZs, hLs = &hLZs01; iZY<2; iZY++) {
	findMWA(hSs,HHC1,HMA1,t,iZY,3,hLs);
	hSs = hSYs; hLs = &hLYs01;
      }
      list<MWAHitList>::iterator ihLZ, ihLY;
      for (ihLZ = hLZs01.begin(); ihLZ!=hLZs01.end(); ihLZ++) {
	MWAHitList &hLZ = *ihLZ;
	for (ihLY = hLYs01.begin(); ihLY!=hLYs01.end(); ihLY++) {
	  MWAHitList &hLY = *ihLY;
	  Track2MWA t2mwa;
	  TTrack &tmwa = t2mwa.tmwa; tmwa.Hfirst = HHC1;
	  unsigned long long &hPat = t2mwa.hitPlPat; hPat = hitPlPat;
	  MWAHitList *hL; for (iZY = 0, hL = &hLZ; iZY<2; iZY++) {
	    for (int ih = 0; ih<hL->nhs; ih++) {
	      THit *h = hL->hL[ih]; if (!h) continue;
	      unsigned long long bpl = b1<<h->IPlane-ipli; hPat |= bpl;
	      tmwa.AddHit(*h);
	    }
	    hL = &hLY;
	  }
	  if (TOpt::ReMode[21]==0) tmwa.InsertMissedPlanes();
	  bool ok = tmwa.FullKF(1,-1,false); // ********** FIT MA01 (PA01X) HITS
	  if (ok) {
	    // No free parameter in KF (all is constrained by original track)
	    // => NDF = #hits in track extension ("TTRack::NDFs", in fact)
	    // (Note: this is fortunate: OK whatever this #hits, if finite.)
	    ok = tmwa.Chi2tot<chi2IperHit*tmwa.NDFs;
#ifdef FT_DEBUG_MA     
	    static int iZ, iY;
	    if (ihLZ==hLZs01.begin()) iZ = 0; if (ihLY==hLYs01.begin()) iY = 0;
	    printf("Try Z=%d x Y=%d: %d %d... %d... => %.2f\n",
		   iZ++,iY++,tmwa.NDFs,hLZ.hL[0]->iHit,hLY.hL[0]->iHit,
		   tmwa.Chi2tot/tmwa.NDFs);
#endif
	    if (!ok) {
	      int worstPl = -1, worstGr = -1; double worstChi2Incr;
	      tmwa.WorstHit(tmwa.mChi2,worstPl,worstGr,worstChi2Incr);
	      unsigned long long bpl = b1<<worstPl-ipli; if (hitPlPat&bpl) {
		//    ***** BAD CHI2/CHI2 INCREMENT: TRY CLEANING AWAY WORST HIT 
		int worstHit = tmwa.CancelHit_Clip(worstPl);
		hitPlPat &= ~bpl; ns01--; // Let's not update "nY|Zs01": no longer used
		ok = tmwa.FullKF(1,-1,false);
		if (ok) ok = tmwa.Chi2tot<chi2IperHit*tmwa.NDFs;
	      }
	    }
	  }
	  if (ok) {                                     // ***** INCLUDE PA01U|V
	    tmwa.IKine = t.IKine;
	    // PA01 already fully associated?
	    if (nPA01s<3 && (nPA01s<2 || (hitPlPat&plPatPA1)==0)) {
	      int addedPAs = findMWPAHits(tmwa,sPA01,hsPA01,nPA01s);
	      if (addedPAs+nPA01s==0 && (hitPlPat&plPatPA1)) {
		// If no new PA01 to back the PA01X found in the MA/PA search
		// => remove it (but those PA01s associated to "t" are kept).
		for (ihp = tmwa.lHitPat.begin(); ihp!=tmwa.lHitPat.end();
		     ihp++) {
		  if (*ihp<0) continue; THit &h = vecHit[*ihp];
		  unsigned long long bpl = b1<<h.IPlane-ipli;
		  if (!(bpl&plPatPA1)) continue;
		  tmwa.CancelHit_Clip(h.IPlane); hitPlPat &= ~bpl; ns01--;
		  ok = tmwa.FullKF(1,-1,false);
		  if (ok) ok = tmwa.Chi2tot<chi2IperHit*ns01;
		  break;
		}
	      }
	    }
	  } // End include PA01UV
	  if (ok) {
	    //              ********** REQUIRE MIN. # of 01-HITS, INCLUDING PA01
	    if ((int)tmwa.NDFs>=nRequired01) t2mwas.push_back(t2mwa);
	  }
	}
      } // End building MWA extension from 01-hits
      if (t2mwas.empty()) { it++; continue; }
    } // End building MWA-tracks out of 01s
    else { // Else special Drell-Yan case: t2mwas.empty...
      Track2MWA t2mwa;       // ...fill w/ original track 
      t2mwa.tmwa.Hfirst = HHC1; t2mwa.hitPlPat = hitPlPat;
      t2mwas.push_back(t2mwa);
    }

    //          ***********************************************
    //          **********  MA02s (and PA02s,HG02) ************
    //          ***********************************************
    double xMW1 = TOpt::MuonWall[0]; THlx HMW1, HMA2;
    list<Track2MWA>::iterator it2;
    for (it2 = t2mwas.begin(); it2!=t2mwas.end(); it2++) {
      TTrack &tmwa = (*it2).tmwa;
      unsigned long long plPat = plPat2|plPatHG;
      for (int kpl = 0; kpl<5; kpl++) hSZs[kpl].clear(); // Hit sets in Z02
      int nZs02, nExpectsZ02; for (ipl = ipli, nZs02=nExpectsZ02 = 0,
				     firstPl = 1; ipl<=iplf; ipl++) {
	unsigned long long bpl = b1<<ipl-ipli;
	if ((bpl&plPat)==0) continue;      // Skip because not a MA/HG02, or OFF
	if (bpl&plYPat) continue;                    // Searching Z planes first

	// ***** LOOP OVER Z PLANES in MA02s and HG02 *****

	const TPlane &pZ = setup.vPlane(ipl);
	const TDetect &dZ = setup.vDetect(pZ.IDetRef);
	float sZ = dZ.Sa, cZ = dZ.Ca, xZ = dZ.X(0), wZ = dZ.Resol;
	if (firstPl) {           // ***** UPON FIRST ENCOUNTERED ACTIVE PLANE...
	  // (We again try to limit the recourse to THlx::Extrapolate.)
	  x0 = tmwa.Hlast(0);    // ***** ...UPDATE INTERCEPTS/SLOPES 
	  y0 = tmwa.Hlast(1), z0 = tmwa.Hlast(2);
	  yp = tmwa.Hlast(3), zp = tmwa.Hlast(4);
	  //                        ***** ...COMPUTE UNCERTAINTIES for MA02/HG02
	  HMW1(0) = xMW1; tmwa.Hlast.Extrapolate(HMW1,true);
	  HMA2(0) = xZ; HMW1.Extrapolate(HMA2,true);
	  y = HMA2(1); z = HMA2(2);
	  dy2 = HMA2(1,1); dz2 = HMA2(2,2); dyz = HMA2(1,2);
	}
	else {
	  y = y0+yp*(xZ-x0); z = z0+zp*(xZ-x0);
	}
	if (!dZ.InActive(y,z)) continue;                 // Out of Z active zone
	firstPl = 0;
	u = y*cZ+z*sZ; du = sqrt(dy2*cZ*cZ+2*dyz*cZ*sZ+dz2*sZ*sZ);
	if (bpl&plPatHG)       { umx = u+4*du+ 6*wZ; umn = u-4*du- 6*wZ;}
	else if (DrellYan&0x2) { umx = u+3*du+10*wZ; umn = u-3*du-10*wZ;}
	else                   { umx = u+3*du+ 4*wZ; umn = u-3*du- 4*wZ;}
	for (ih = pZ.vHitRef().begin(); ih!=pZ.vHitRef().end(); ih++) {
	  THit *h = &vecHit[*ih]; if (h->Status!=0) continue;
	  if (h->U<umn) continue; if (umx<h->U) break;
	  double flag = (bpl&plPatHG) ? 0 :
	    h->ptrCl->getDigitsList().front()->getData()[1];
	  if (flag &&                 // Non null flag means split channel...
	      flag*y<0) continue;     // ...=> Require it to be on right side
	  hSZs[nExpectsZ02].push_back(h);
	}
	if (!hSZs[nExpectsZ02].empty()) nZs02++; nExpectsZ02++;
      }
      if (!nZs02) continue;
      // ********** FOUND POSSIBLE 02-EXTENSION in Z **********

      for (int kpl = 0; kpl<5; kpl++) hSYs[kpl].clear(); // Hit sets in Y02
      int nYs02, nExpectsY02; for (ipl = ipli, nYs02=nExpectsY02 = 0;
				   ipl<=iplf; ipl++) {
	unsigned long long bpl = b1<<ipl-ipli;
	if (!(bpl&(plPat2|plPatPA2))) continue; // Skip because not PA/MA02 or OFF
	if (!(bpl&plYPat)) continue;                    // Searching Y planes 

	// ***** LOOP OVER Y PLANES in MA02s *****

	const TPlane &pY = setup.vPlane(ipl);
	const TDetect &dY = setup.vDetect(pY.IDetRef);
	float sY = dY.Sa, cY = dY.Ca, xY = dY.X(0), wY = dY.Resol;
	y = y0+yp*(xY-x0); z = z0+zp*(xY-x0);
	if (!dY.InActive(y,z)) continue;                 // Out of Y active zone
	if (bpl&hitPlPAt) continue;             // Skip already associated PA02Y
	u = y*cY+z*sY; du = sqrt(dy2*cY*cY+2*dyz*cY*sY+dz2*sY*sY);
	if (DrellYan&0x2) { umx = u+3*du+10*wY; umn = u-3*du-10*wY;}
	else              { umx = u+4*du+ 4*wY; umn = u-4*du- 4*wY;}
	for (ih = pY.vHitRef().begin(); ih!=pY.vHitRef().end(); ih++) {
	  THit *h = &vecHit[*ih]; if (h->Status!=0) continue;
	  if (h->U<umn) continue; if (umx<h->U) break;
	  // for PAs channels are not split, only get channel position for MA02
	  double flag = (bpl&plPat2) ? h->ptrCl->getDigitsList().front()->getData()[1] : 0;
	  if (flag &&                 // Non null flag means split channel...
	      flag*z<0) continue;     // ...=> Require it to be on right side
	  hSYs[nExpectsY02].push_back(h);
	}
	if (!hSYs[nExpectsY02].empty()) nYs02++; nExpectsY02++;
      }
      if (!nYs02) continue;
      // ********** FOUND POSSIBLE 02-EXTENSION in Y **********

      int ns012 = tmwa.NDFs+nZs02+nYs02;
      //              ********** REQUIRE MIN. OVERALL # of POSSIBLE HITS
      // (Note: we could decide to let through cases slighly (1 hit, e.g) below
      // cut, keeping in mind that one can always disregard them later on.)
      if (nPA01s+nPA02s+ns012+/*possibly 2 PA02U|V*/2<TOpt::iCut[6]) continue;

      // ********** EXTEND MWA-TRACK w/ MA02s (+ possibly PA02s|HG02) **********
      list<MWAHitList> hLZs02, hLYs02;
      int iZY; vector<THit*> *hSs; list<MWAHitList> *hLs;
      for (iZY = 0, hSs = hSZs, hLs = &hLZs02; iZY<2; iZY++) {
	findMWA(hSs,HMW1,HMA2,t,iZY,iZY?3:nSigmasMA02Y,hLs);
	hSs = hSYs; hLs = &hLYs02;
      }
      list<MWAHitList>::iterator ihLZ, ihLY; int nExtensions;
      for (ihLZ = hLZs02.begin(), nExtensions = 0; ihLZ!=hLZs02.end(); ihLZ++) {
	MWAHitList &hLZ = *ihLZ;
	for (ihLY = hLYs02.begin(); ihLY!=hLYs02.end(); ihLY++) {
	  MWAHitList &hLY = *ihLY;
	  Track2MWA t2mwa2;
	  TTrack &tmwa2 = t2mwa2.tmwa; tmwa2 = tmwa;
	  unsigned long long &hPat = t2mwa2.hitPlPat; hPat = (*it2).hitPlPat;
	  MWAHitList *hL; for (iZY = 0, hL = &hLZ; iZY<2; iZY++) {
	    for (int ih = 0; ih<hL->nhs; ih++) {
	      THit *h = hL->hL[ih]; if (!h) continue;
	      unsigned long long bpl = b1<<h->IPlane-ipli; hPat |= bpl;
	      tmwa2.AddHit(*h);
	    }
	    hL = &hLY;
	  }
	  if (TOpt::ReMode[21]==0) tmwa2.InsertMissedPlanes();
	  bool ok =                       // ********** FIT MA02|PA02X|HG02 HITS
	    tmwa2.FullKF(1,-1,false);
	  if (ok) {
	    ok = tmwa2.Chi2tot<chi2IperHit*tmwa2.NDFs;
	    if (!ok && (plPatHG&hPat)) // Relaxed cut for HG02 case.
	      // - HG02 adds to the MWA track reliability, since its two pieces,
	      //  Saleve/Jura, provide some of the vertic/horizontal correlation
	      //  that is otherwise lacking in the muWallA sub-zone.
	      //  (Split MA channels provide also some, btw.)
	      // - This added reliability could be further stengthened by:
	      //   - timing (not yet done)
	      //   - 2D info, would HG02 CsCluster's be assigned a v-coordinate,
	      //    along horizontal axis, derived from the time diff. between
	      //    their Jura and Saleve TDCs (not yet available).
	      ok = tmwa2.Chi2tot<chi2IperHit*ns012*1.5;
	  }
	  if (ok) {                                     // ***** INCLUDE PA02U|V
	    // Cf. the handling of PA01U|V supra
	    if (nPA02s<3 && (nPA02s<2 || (hPat&plPatPA2)==0)) {
	      int addedPAs = findMWPAHits(tmwa2,sPA02,hsPA02,nPA02s);
	      if (addedPAs+nPA02s==0 && (hitPlPat&plPatPA2)) {
		for (ihp = tmwa2.lHitPat.begin(); ihp!=tmwa2.lHitPat.end();
		     ihp++) {
		  if (*ihp<0) continue; THit &h = vecHit[*ihp];
		  unsigned long long bpl = b1<<h.IPlane-ipli;
		  if (!(bpl&plPatPA2)) continue;
		  tmwa2.CancelHit_Clip(h.IPlane); hPat &= ~bpl;
		  ok = tmwa2.FullKF(1,-1,false);
		  if (ok) ok = tmwa.Chi2tot<chi2IperHit*tmwa2.NDFs;
		  break;
		}
	      }
	    }
	  } // End include PA02UV
	  if (ok)                  //  ********** REQUIRE MIN. OVERALL # of HITS
	    // (Note: could decide to let through cases slighly (1 hit, e.g.)
	    // below par: keeping in mind they can anyway be disregard later on.)
	    ok = nPA01s+nPA02s+(int)tmwa2.NDFs>=TOpt::iCut[6];

	  if (ok) {      // ********** EVALUATE EXTENSION'S Q FUNCTION AND STORE
	    t2mwa2.it = it; t2mwa2.nPA01s = nPA01s;
	    int ndfs = tmwa2.NDFs; double qual = tmwa2.Chi2tot/ndfs;
	    double eff = (double)ndfs/nExpects; qual /= eff*eff;
	    if (hitPlPat&plPatHG) // Some bonus for HG02, cf. supra...
	      qual *= .80;
	    // Low momenta tracks experience dramatic MS in the calos/absorbers
	    // and may hence tend to become compatible w/ almost anything.
	    // => Let's cook up something to degrade their Q function (the
	    //  cooking recipe is based on event #22234019 of Catarina's 2012 DY
	    //  file, where 2 two tracks are vying for muWallA muID):
	    qual *= 1+fabs(t.Hfirst(5));  // Penalty for low momenta (useful?)
	    if (t.Hlast(0)<750)// Malus for when no reco'd downstream of RICH...
	      qual *= 1.1;     // ...because then chi2 is less constrained.
	    t2mwa2.qual = qual;
	    lTrack2MWA.push_back(t2mwa2);
	  }
	} // End loop on 02 Y hit lists
      } // End loop on 02 Z hit lists
    } // End loop on Track2MWA stubs found in 01 piece of muWallA
#ifdef FT_DEBUG_MA     
    printf("++++++++++\n");
    list<Track2MWA>::iterator il;
    int itmwa; for (il = lTrack2MWA.begin(), itmwa = 0; il!=lTrack2MWA.end();
		    il++, itmwa++) {
      if (il->it->Id!=t.Id) continue;
      TTrack &tmwa = il->tmwa;
      printf("%2d: %2d",itmwa,tmwa.NDFs);
      list<int>::iterator ihp; for (ihp = tmwa.lHitPat.begin();
				    ihp!=tmwa.lHitPat.end(); ihp++) {
	if (*ihp<0) continue; THit &h = vecHit[*ihp];
	printf(" %2d/%4d",h.IPlane,h.IHit);
      }
      printf(" => %.2f => %.2f\n",tmwa.Chi2tot/tmwa.NDFs,il->qual);
    }
    printf("----------\n");
#endif

    it++;
  } // End loop on listTrack

  lTrack2MWA.sort();  //  *************** SORT BY Q FUNCTION ***************

  //                  *************** LOOP on THE ORDERD LIST ***************
  list<Track2MWA>::iterator il;
  for (il = lTrack2MWA.begin(); il!=lTrack2MWA.end(); il++) {
    Track2MWA &tp = *il; TTrack &tmwa = tp.tmwa;
    list<int>::iterator ihp; int nYs, nZs;
    for (ihp = tmwa.lHitPat.begin(), nYs=nZs = 0; ihp!=tmwa.lHitPat.end();
	 ihp++) {
      if (*ihp<0) continue; THit &h = vecHit[*ihp];
      if (h.Status==0) {
	if (setup.vPlane(h.IPlane).IProj==0) nYs++;
	else                                 nZs++;
      }
      else // Let HG02 hits be shared among several tracks, just as for any
	// hodo. (Note that they contribute to the "nZs" hit count only once)
	if (!(b1<<h.IPlane-ipli&plPatHG)) tmwa.CancelHit(h.IPlane);
    }
    // Check that enough hits remain now that some may have been cancelled
    int nPA01s = tp.nPA01s; if (nPA01s+(int)tmwa.NDFs<TOpt::iCut[6]) continue;
    // Do also the check independently for Y and Z
    if ((nPA01s?1:0)+nYs<TOpt::iCut[6]/3 || nZs<TOpt::iCut[6]/3) continue;
    TTrack &t = *(tp.it);
    bool fillHist = hChi2; static int nGenuineMWAs; static double oldChi2;
    if (fillHist) {
      t.FindKine(); fillHist &= t.IKine>=0;
      nGenuineMWAs = 0; oldChi2 = t.Chi2tot/(t.NDFs-5);
    }
    for (ihp = tmwa.lHitPat.begin(); ihp!=tmwa.lHitPat.end(); ihp++) {
      if (*ihp<0) continue; THit &h = vecHit[*ihp]; h.status = 1; t.AddHit(h);
      if (fillHist && (h.IKine==t.IKine || h.IKin2==t.IKine)) nGenuineMWAs++;
    }
    if (TOpt::ReMode[21]==0) t.InsertMissedPlanes();
    bool ok = t.FullKF(1); if (ok) ok = t.FullKF(-1);
    if (ok) {
      double chi2 = t.Chi2tot/(t.NDFs-5); ok = chi2<TOpt::dCut[17];
      if (fillHist) {
	int nMWAs = tmwa.NDFs;
	hChi2->Fill(chi2,nMWAs); hdChi2->Fill(chi2-oldChi2,nMWAs);
	if (t.NHsame+nGenuineMWAs<TOpt::dCut[10]*t.NHits)
	  CsErrLog::msg(elError,__FILE__,__LINE__,
	    "Event #%d track ID=%d,MC=%d,%d/%d downgraded => %d/%d",event,t.Id,
	       t.IKine,t.NHsame,t.NDFs-nMWAs,t.NHsame+nGenuineMWAs,t.NDFs);
      }
    }
    if (!ok) {           // ********** IF FAILURE => RESTORE INITIAL TRACK
      for (ihp = tmwa.lHitPat.begin(); ihp!=tmwa.lHitPat.end(); ihp++) {
	if (*ihp<0) continue;
	THit &h = vecHit[*ihp]; h.status = 0; t.CancelHit(h.IPlane);
      }
      t.Clip(true /* i.e. downstream */);  // Erase inserted planes
      ok = t.FullKF(1); if (ok) ok = t.FullKF(-1);
      if (!ok) {
	t.IMark = -2; // Flag for deletion
	CsErrLog::msg(elError,__FILE__,__LINE__,
	  "Event #%d: Full KF fails when track ID=#%d is restored",event,t.Id);
      }
    }
#ifdef DEBUG
    else if (idebug) it->DumpHits(idebug);
#endif
  } // End loop on tracks
  //                                                ***** RE-ASSESS FISHY TRACKS
  // (Cf. similar block in "ForeTrack2RICHWall".)
  unsigned long long plPat = plPat2|plPatHG;
  list<unsigned int>::iterator jt;
  for (jt = fishyTrks.begin(); jt!=fishyTrks.end(); jt++) {
    unsigned int id = *jt; for (it = listTrack.begin(); it!=listTrack.end();
				it++) {
      if ((*it).Id!=id) continue;
      TTrack &t = *it; list<int>::iterator ihp, ipr;
      for (ihp = t.lHitPat.begin(), ipr = t.lPlnRef.begin();
	   ihp!=t.lHitPat.end(); ihp++, ipr++) {
	int ipl = *ipr; if (ipl<ipli) continue; if (ipl>iplf) break;
	unsigned long long bpl = b1<<ipl-ipli; if (!(bpl&plPat)) continue;
	int ihit = *ihp; if (ihit<0) continue;
	THit &h = vecHit[ihit]; if (h.sTrackID().size()>1) {
	  listTrack.erase(it); break;
	}
      }
      break;
    }
  }
  //                             ***** ERASE TRACKS NO LONGER PASSING THE KF FIT
  it = listTrack.begin(); while (it!=listTrack.end()) {
    if ((*it).IMark==-2) { listTrack.erase(it++); continue; }
    else it++;
  }
}
double getMWAChi2(THit **hs, int nhs,
		  double xhinge, double hinge, double dhinge);
void histMWAHitList(THit **hs, int nhs, double chi2, double qual, int kine);
void findMWA(vector<THit*> *hSs,               // Hit sets (in muWallA X|Y 1|2)
	     const THlx &Hinge, const THlx &H, // Helix of original track @ hinging point and first muWallA
	     const TTrack &t, // Original track (candidate for MWA extension)
	     int iZY, // TraFDic's Y = horizontal <-> X planes, Z = vertical <-> Yplanes
	     int nSigmasF, // Road w/ for the hit search in first plane
	     list<MWAHitList> *mwaHLs) // stl::List of output hit lists
{
  // Note: Argument track so far only conveys it IKine (one could later also
  // make use of its mean time).

  //               ********** INITIALIZATION **********
  const TSetup &setup = TSetup::Ref();
  // Hinging point (used to constrain Q function)
  double xHinge = Hinge(0), hinge, dHinge2;
  if (iZY) { hinge = Hinge(1); dHinge2 = Hinge(1,1); }
  else     { hinge = Hinge(2); dHinge2 = Hinge(2,2); }
  // Hit pattern:
  // - Is 64 bit long: should be enough in most cases, since we here concentrate
  //  on those sole hits that are compatible w/ original track.
  // - Is used to single out hits already used by better combinations and reject
  //  those combinations that contain too many of them.
  // - HG02 is excluded from this used hits business, for its crude granularity.
  const unsigned long long b1 = 0x1; unsigned long long hitHG02;
  int kpl, lastPl, lastPMA, nPMAs/* #planes but HG02 */, yetHG, ih0s[5], ih0;
  for (kpl=ih0=lastPl=nPMAs=yetHG=lastPMA = 0, hitHG02 = 0; kpl<6; kpl++) {
    ih0s[kpl] = ih0; int size = hSs[kpl].size(); if (!size) continue;
    if (setup.vDetect(hSs[kpl][0]->IPlane).Name.find("HG02")==0) {
      int ih; unsigned long long b; for (ih = 0, b = b1<<ih0; ih<size;
					 ih++, b <<= 1) hitHG02 |= b;
      yetHG = 1;
    }
    ih0 += size; lastPl = kpl; if (!yetHG) { lastPMA = kpl; nPMAs++; }
  }
  //int nhsTot = ih0; unsigned long long bMx = b1<<nhsTot-1; // Useless: cf. infra.
  // The search bases itself on the two pivots planes. These are, by default,
  // 1st two and last two detector planes. Let's single out first the cases
  // where this is not possible.
  if (nPMAs<2) return;
  int nPivotsf, nPivotsl;
  if      (nPMAs>=4) { nPivotsf=nPivotsl = 2; }
  else if (nPMAs==3) { nPivotsf = 2; nPivotsl = 1; }
  else               { nPivotsf=nPivotsl = 1; }
  THlx Hp; // Working helix
  double chi2Cut = TOpt::dCut[65] ? TOpt::dCut[65] : 6;
#if defined FT_DEBUG_MA && FT_DEBUG_MA > 1
  printf("---------- %c\n",iZY?'X':'Y');
  for (kpl = 0; kpl<6; kpl++) {
    int size = hSs[kpl].size(); if (!size) continue;
    printf("%d: %d",kpl,hSs[kpl][0]->IPlane);
    for (int ih = 0; ih<size; ih++) printf(" %d",hSs[kpl][ih]->IHit);
    printf("\n");
  }
  printf("----------\n");
#endif

  int fpl, iPivotf; for (fpl = 0, iPivotf = 0; iPivotf<nPivotsf; fpl++) {
    vector<THit*> &hfs = hSs[fpl]; if (hfs.empty()) continue; iPivotf++;
    const TDetect &df = setup.vDetect(hfs[0]->IPlane); double xf = df.X(0);
    THlx Hf; Hf(0) = xf; H.Extrapolate(Hf,true);
    for (int ihf = 0; ihf<(int)hfs.size(); ihf++) {
      THit *hf = hfs[ihf]; THlx Hhf; Hf.Update(*hf,Hhf);
      double duf2 = // Simplifying out deviation of plane orientation from Y|Z
	iZY ? Hhf(1,1) : Hhf(2,2);
      THit *hL[5]; hL[0] = hf; MWAHitList hLhf; hLhf.mark = 0;
      int lpl, iPivotl; double bestQ;
      for (lpl = lastPMA, iPivotl = 0, bestQ = 0; iPivotl<nPivotsl; lpl--) {
	vector<THit*> &hls = hSs[lpl]; if (hls.empty()) continue; iPivotl++;
	const TDetect &dl = setup.vDetect(hls[0]->IPlane); double xl = dl.X(0);
	THlx Hl; Hl(0) = xl; Hhf.Extrapolate(Hl,true);
	double sl = dl.Sa, cl = dl.Ca, wl = dl.Resol, dUlMx = 2*wl;
	double yl = Hl(1), zl = Hl(2), ul = yl*cl+zl*sl;
	for (int ihl = 0; ihl<(int)hls.size(); ihl++) {
	  THit *hl = hls[ihl]; double U = hl->U;
	  double dU = hl->SigU>dUlMx?dUlMx:hl->SigU, duUf = sqrt(duf2+dU*dU);
	  if (U<ul-nSigmasF*duUf) continue; if (ul+nSigmasF*duUf<U) break;
	  THlx Hhfl; Hl.Update(*hl,Hhfl);
	  double yhl = Hhfl(1), zhl = Hhfl(2), yphl = Hhfl(3), zphl = Hhfl(4);
	  double dul2 = iZY ? Hhfl(1,1) : Hhfl(2,2);
	  hL[1] = hl;
	  unsigned long long hitPat = b1<<ih0s[fpl]+ihf;
	  int kpl, nhs, nHG02s; for (kpl = fpl+1, nhs = 2, nHG02s = 0;
				     kpl<=lastPl; kpl++) {
	    bool isHG02 = kpl>lastPMA;
	    if (kpl>=lpl && !isHG02) // Let's only interpolate bewteen pivots...
	      continue;  // ...and do not extrapolate, but into HG02s
	    vector<THit*> &hks = hSs[kpl]; if (hks.empty()) continue;
	    const TDetect &dk = setup.vDetect(hks[0]->IPlane); double xk = dk.X(0);
	    double sk = dk.Sa, ck = dk.Ca, wk = dk.Resol, dUkMx = 2*wk;
	    double yk = yhl+yphl*(xk-xl), zk = zhl+zphl*(xk-xl), uk = yk*cl+zk*sl;
	    THit *hBest; static double bestR; int ihk; static int bestI;
	    for (ihk = 0, hBest = 0; ihk<(int)hks.size(); ihk++) {
	      THit *hk = hks[ihk]; U = hk->U;
	      dU = hk->SigU>dUkMx?dUkMx:hk->SigU; double duUl = sqrt(dul2+dU*dU);
	      if (U<uk-4*duUl) continue; if (uk+4*duUl<U) break;
	      double r = fabs(U-ul); if (!hBest || r<bestR) {
		hBest = hk; bestR = r; bestI = ihk;
	      }
	    }
	    if (hBest) {
	      hL[nhs++] = hBest; hitPat |= b1<<ih0s[kpl]+bestI;
	      if (isHG02) nHG02s++;
	    }
	  } // End loop on inter-pivot planes
	  double chi2 = getMWAChi2(hL,nhs,xHinge,hinge,dHinge2);
	  double qual = nhs           // Q function based on #hits
	    -(nHG02s>1?1:0)// (two HG02 hits may be associated: count only one)
	    +TMath::Prob(chi2,nhs-1); // W/ chi2 probability bonus of at most 1
	  chi2 /= nhs-1;
	  if (TOpt::Hist[7]) histMWAHitList(hL,nhs,chi2,qual,t.GetIKine());
	  if (chi2<chi2Cut && qual>bestQ) {
	    hLhf.hitPat = hitPat|b1<<ih0s[lpl]+ihl;
	    for (int ih = 0; ih<nhs; ih++) hLhf.hL[ih] = hL[ih]; hLhf.nhs = nhs;
	    hLhf.qual = qual;
	    bestQ = qual;
	  }
	}
      }
      if (bestQ) mwaHLs->push_back(hLhf);
    }
  }

  mwaHLs->sort();     // ********** SORT MWAHitList's... **********

  // ********** ...RETAIN, MUTUALLY (almost) EXCLUSIVE, BEST ONES **********
  list<MWAHitList>::iterator ihL, jhL, ihL0 = mwaHLs->begin();
#ifdef FT_DEBUG_MA
  int i; for (ihL = ihL0, i=0; ihL!=mwaHLs->end(); ihL++, i++) {
    printf("%2d: ",i); ihL->Print();
  }
#endif
  for (ihL = ihL0; ihL!=mwaHLs->end(); ihL++) {
    MWAHitList &hLi = *ihL; unsigned long long hPati = hLi.hitPat&(~hitHG02);
    int match; for (jhL = ihL0, match = 0; jhL!=ihL; jhL++) {
      MWAHitList &hLj = *jhL; if (!hLj.mark) continue;
      unsigned long long overlap = hLj.hitPat&hPati; if (overlap) {
	if (overlap==hPati ||  // "j" is subset of "i" or
	    hLi.nhs<=3) {      //  i is anyway too scare...
	  match = 1; break;    // ... => skip
	}
	int nOverlaps; unsigned long long b;
	//for (b = b1, nOverlaps = 0; b<=bMx; b <<= 1)
	// The above does not work: when "b" reaches 0x8000000000000000, it
	// rolls over...
	int ib; for (ib = 0, b = b1, nOverlaps = 0; ib<64; ib++, b <<= 1)
		  if (b|overlap) nOverlaps++;
	if (nOverlaps>1) { match = 1; break; }
      }
    }
    if (match) continue;
    // Check current hL is not contained entrirely in a combination, that, if
    // it cannot only be but strictly worse, would not be that much worse.
    jhL = ihL; for (jhL++; jhL!=mwaHLs->end(); jhL++) {
      MWAHitList &hLj = *jhL;
      unsigned long long overlap = hLj.hitPat&hLi.hitPat;
      if (overlap==hLi.hitPat && hLj.qual<hLi.qual*1.2) { match = 1; break; }
    }
    if (match) continue;
    hLi.mark = 1; 
  }
  ihL = ihL0; while (ihL!=mwaHLs->end()) {
    MWAHitList &hLi = *ihL; if (!hLi.mark) { mwaHLs->erase(ihL++); continue; }
    ihL++;
  }
}
double getMWAChi2(THit **hs, int nhs,
		  double xHinge, double hinge, double dHinge2)
{
  // Returns total chi2
  const TSetup &setup = TSetup::Ref();

  double w2 = dHinge2;

  double s1 = 1/w2, sx = xHinge/w2, sx2 = xHinge*xHinge/w2;
  double sU = hinge/w2, sU2 = hinge*hinge/w2, sxU = xHinge*hinge/w2;

  for (int ih = 0; ih<nhs; ih++) {
    const THit *h = hs[ih]; const TDetect &d = setup.vDetect()[h->IPlane];
    double xi = d.X(0), dUMx = 2*d.Resol;
    double Ui = h->U, dUi = h->SigU>dUMx ? dUMx : h->SigU, w2 = dUi*dUi;
    s1 += 1/w2; sx += xi/w2; sx2 += xi*xi/w2;
    sU += Ui/w2; sU2 += Ui*Ui/w2; sxU += xi*Ui/w2;
  }
  double D = (sx*sx-s1*sx2), a = (sx*sU-s1*sxU)/D, b = (sx*sxU-sx2*sU)/D;
#if defined FT_DEBUG_MA && FT_DEBUG_MA > 1
  double chi2; int ih; for (ih = -1, chi2 = 0; ih<nhs; ih++) {
    double xi, w2, Ui;
    if (ih>=0) {
      const THit *h = hs[ih]; const TDetect &d = setup.vDetect()[h->IPlane];
      xi = d.X(0); double dUMx = 2*d.Resol;
      Ui = h->U; double dUi = h->SigU>dUMx ? dUMx : h->SigU; w2 = dUi*dUi;
      printf (" %2d/%4d ",h->IPlane,h->IHit);
    }
    else {
      xi = xHinge; Ui = hinge; w2 = dHinge2;
    }
    double ui = a*xi+b, r = Ui-ui, dchi2 = r*r/w2; chi2 += dchi2;
    printf ("%4.1f %5.3f %5.3f",r,dchi2,chi2);
  }
  printf(" => %5.3f/%5.3f(%5.3f) %5.3f\n",
	 chi2,chi2/(nhs-1),(sU2-a*sxU-b*sU)/(nhs-1),nhs+TMath::Prob(chi2,nhs-1));
#endif
  return sU2-a*sxU-b*sU;
}
void histMWAHitList(THit **hs, int nhs, double chi2, double qual, int kine)
{
  static TH2D *hMWAchi2[2][2], *hMWAFoM[2][2]; static int firstMA02 = 0;
  if (!firstMA02) {
    const TSetup &setup = TSetup::Ref();
    for (int ipl = 0; ipl<(int)setup.vDetect().size(); ipl++) {
      if (setup.vDetect(ipl).Name.find("MA02")==0) { firstMA02 = ipl; break; }
    }
    CsHistograms::SetCurrentPath("/Traffic/RDmonitor");
    char hN[] = "hMWA01Ychi2", hT[] = "muWallA: 01Y chi2 - fake/genuine/orig";
    for (int i12 = 1; i12<=2; i12++) for (int iZY = 0; iZY<2; iZY++) {
	sprintf(hN,"hMWA0%d%cchi2",i12,iZY?'X':'Y');
	sprintf(hT,"muWallA: 0%d%c chi2 - fake/genuine/orig",i12,iZY?'X':'Y');
	hMWAchi2[i12-1][iZY] = new TH2D(hN,hT,128,0,8,3,-.5,2.5);
	sprintf(hN,"hMWA0%d%cFoM",i12,iZY?'X':'Y');
	sprintf(hT,"muWallA: 0%d%c FoM - fake/genuine/orig",i12,iZY?'X':'Y');
	hMWAFoM[i12-1][iZY] = new TH2D(hN,hT,128,0,8,3,-.5,2.5);
      }
    CsHistograms::SetCurrentPath("/");
  }
  const TSetup &setup = TSetup::Ref();

  int ih, genuine, orig, i12; static int iZY;
  for (ih = 0, genuine=orig = 1, i12 = 0; ih<nhs; ih++) {
    const THit *h = hs[ih]; if (!h) continue;
    genuine = genuine && (h->IKine==kine || h->IKin2==kine) ? 1 : 0;
    orig = genuine && orig && h->IOrig==0 ? 1 : 0;
    if (!i12) {
      iZY = TSetup::Ref().vPlane(h->IPlane).IProj==0 ? 0 : 1;
      i12 = h->IPlane>=firstMA02 ? 2 : 1;
    }
  }
  if (i12) {
    hMWAchi2[i12-1][iZY]->Fill(chi2,genuine+orig);
    hMWAFoM[i12-1][iZY]->Fill(qual,genuine+orig);
  }
  else CsErrLog::mes(elFatal,"Software code inconsistency");
}
#include "TMatrixD.h"
#include "TVectorD.h"
int TEv::findMWPAHits(TTrack &t,            // Track extension in the MWA sub-zone
		      const TStation *sPA,  // PA TStation
		      THit **hsPA, int nPAs)// PA hits associated to the would-be extended track		      , 
{
  // Search TStation "s" for hits compatible w/ last helix of track "t"
  const TSetup &setup = TSetup::Ref();

  static TH2D *hMWAPAchi2 = 0, *hMWAPAres3, *hMWAPAres2; if (!hMWAPAchi2) {
    CsHistograms::SetCurrentPath("/Traffic/RDmonitor");
    hMWAPAchi2 = new TH2D("hMWAPAchi2","muWallA: PA spacept chi2 - fake/genuine/orig",
			  128,0,16,3,-.5,2.5);
    hMWAPAres3 = new TH2D("hMWAPAres3","muWallA: PA space pt residual^2 - fake/genuine/orig",
			  128,0,16,3,-.5,2.5);
    hMWAPAres2 = new TH2D("hMWAPAres2","muWallA: PA 2-hit residual^2 - fake/genuine/orig",
			  128,0,16,3,-.5,2.5);
    CsHistograms::SetCurrentPath("/");
  }
  int kine = t.GetIKine();

  // Already associated hits?
  int ipli = sPA->IPlanes[0]; vector<THit*> hs[3];
  unsigned int associated; int i, nAssociated;
  for (i=nAssociated = 0, associated = 0; i<nPAs; i++) {
    THit *h = hsPA[i]; int jpl = h->IPlane-ipli;
    hs[jpl].push_back(h); associated |= 1<<jpl; nAssociated++;
  }
  for (list<int>::iterator ihp = t.lHitPat.begin(); ihp!=t.lHitPat.end();
       ihp++) {
    if (*ihp<0) continue; THit &h = vecHit[*ihp]; int jpl = h.IPlane-ipli;
    if (jpl<0) continue; if (2<jpl) break;
    hs[jpl].push_back(&h); associated |= 1<<jpl; nAssociated++;
  }
  if (nAssociated==3) return 0; // All three XUV PA planes associated => exit.

  // ********** SEARCH FOR HITS COMPATIBLE W/ ARG. TRACK (last helix) **********
  THlx &H = t.Hlast;
  double x0 = H(0), y0 = H(1), z0 = H(2), yp = H(3), zp = H(4);
  static double dy2, dz2, dyz;
  int jpl, npls, firstPl; for (jpl=npls = 0, firstPl = 1; jpl<3; jpl++) {
    int ipl = sPA->IPlanes[jpl]; if (1<<jpl&associated) continue;
    const TPlane &p = setup.vPlane(ipl); if (!p.IFlag) continue;  // Plane OFF
    const TDetect &d = setup.vDetect(ipl);
    double xd = d.X(0), y = y0+yp*(xd-x0), z = z0+zp*(xd-x0);
    if (!d.InActive(y,z)) continue;
    if (firstPl) {
      THlx HPA01; HPA01(0) = xd; H.Extrapolate(HPA01,true);
      dy2 = HPA01(1,1); dz2 = HPA01(2,2); dyz = HPA01(1,2); firstPl = 0;
    }
    double sd = d.Sa, cd = d.Ca, wd = d.Resol, dUMx = 2*wd;
    double u = y*cd+z*sd, du2 = dy2*cd*cd+2*dyz*cd*sd+dz2*sd*sd;
    for (vector<int>::const_iterator ih = p.vHitRef().begin();
	 ih!=p.vHitRef().end(); ih++) {
      THit *h = &vecHit[*ih]; if (h->Status) continue;
      double U = h->U, dU = h->SigU>dUMx?dUMx:h->SigU, duU = sqrt(du2+dU*dU);
      if (U<u-3*duU) continue; if (u+3*duU<U) break; hs[jpl].push_back(h);
    }
    if (!hs[jpl].empty()) npls++;
  }
  if (npls==0 ||          // If no new hit candidate for association...
      npls+nAssociated<2) // ...or overall #hits too small
    // (Notes: 2nd condition is fulfilled if 1st one is when "findMWPAHits" is
    // to complement a PA0?X hit found in the MA/PA search of "ForeTrack2MA".
    // A non null npls is required infra, so that dy|z2 be defined.)
    return 0;             // ... => exit.

  //     ********** EVALUATE COMBINATIONS OF COMPATIBLE HITS **********
  // - If any space point (i.e. 3 mutually compatible hits), retain best.
  // - Else retain the 2-hit combination if it's unique.
  int jpl1 = 0; while (hs[jpl1].empty()) jpl1++;
  const TDetect &d1 = setup.vDetect(ipli+jpl1);
  double w2 = d1.Resol*d1.Resol;// Simplifying assumption: same resolution for all
  double c1 = d1.Ca, s1 = d1.Sa, x1 = d1.X(0), c, s, u;
  double y1 = y0+(x1-x0)*yp, z1 = z0+(x1-x0)*zp;
  THit *hsBest[3]; double best; int type;
  int ih1; for (ih1 = 0, best = 9, type = 0;
		ih1<(int)hs[jpl1].size(); ih1++) {
    THit *h1 = hs[jpl1][ih1]; double U1 = h1->U, dU1 = h1->SigU;
    c = c1/dU1; s = s1/dU1; u = U1/dU1;
    double sc21 = c*c, ssc1 = s*c, ss21 = s*s;
    double scu1 = c*u, ssu1 = s*u, su21 = u*u;
    int genuine = kine>0 && (h1->IKine==kine || h1->IKin2==kine) ? 1 : 0;
    int orig = genuine && h1->IOrig==0 ? 1 : 0;

    int jpl2 = jpl1+1; while (hs[jpl2].empty()) jpl2++;
    const TDetect &d2 = setup.vDetect(ipli+jpl2);
    double c2 = d2.Ca, s2 = d2.Sa, incrU2 = (d2.X(0)-x1)*(yp*c2+zp*s2);
    for (int ih2 = 0; ih2<(int)hs[jpl2].size(); ih2++) {
      THit *h2 = hs[jpl2][ih2]; double U2 = h2->U+incrU2, dU2 = h2->SigU;
      c = c2/dU2; s = s2/dU2; u = U2/dU2;
      double sc22 = sc21+c*c, ssc2 = ssc1+s*c, ss22 = ss21+s*s;
      double scu2 = scu1+c*u, ssu2 = ssu1+s*u, su22 = su21+u*u;
#ifdef FT2MA_HISTOS
      genuine = genuine && (h2->IKine==kine || h2->IKin2==kine) ? 1 : 0;
      orig = genuine && orig && h2->IOrig==0 ? 1 : 0;
#endif

      if (jpl2==1) {
	const TDetect &d3 = setup.vDetect(ipli+2);
	double c3 = d3.Ca, s3 = d3.Sa, incrU3 = (d3.X(0)-x1)*(yp*c3+zp*s3);
	for (int ih3 = 0; ih3<(int)hs[2].size(); ih3++) {
	  THit *h3 = hs[2][ih3]; double U3 = h3->U+incrU3, dU3 = h3->SigU;
	  c = c3/dU3; s = s3/dU3; u = U3/dU3;
	  double sc2 = sc22+c*c, ssc = ssc2+s*c, ss2 = ss22+s*s;
	  double scu = scu2+c*u, ssu = ssu2+s*u, su2 = su22+u*u;
	  double ss[4] = {sc2,ssc,ssc,ss2}, su[2] = {scu,ssu};
	  TMatrixD m(2,2,ss), mi(TMatrixD::kInverted,m);
	  TVectorD v(2,su), yz = mi*v; double y = yz[0], z = yz[1]; 
	  double chi2 = y*y*sc2+z*z*ss2+su2+2*(y*z*ssc-y*scu-z*ssu);
#ifdef FT2MA_HISTOS
	  genuine = genuine && (h3->IKine==kine || h3->IKin2==kine) ? 1 : 0;
	  orig = genuine && orig && h3->IOrig==0 ? 1 : 0;
	  hMWAPAchi2->Fill(chi2,genuine+orig);
#endif
	  if (chi2>9) continue;
	  double ry = y-y1, rz = z-z1, r2 = ry*ry+rz*rz;
	  double dr2 = 2*sqrt(ry*ry*(mi(0,0)+dy2)+2*ry*rz*(mi(0,1)+dyz)+rz*rz*(mi(1,1)+dz2));
#ifdef FT2MA_HISTOS
	  hMWAPAres3->Fill(r2/dr2,genuine);
#endif
	  if (r2/dr2<3) {
	    best = r2; hsBest[0] = h1; hsBest[1] = h2; hsBest[2] = h3;
	    type = 0x7;
	  }
	}
      }
      if (type==0) { // No space point yet nor any previous 2-hit combination
	// => Evaluate 2-hit
	double ss[4] = {sc22,ssc2,ssc2,ss22}, su[2] = {scu2,ssu2};
	TMatrixD m(2,2,ss), mi(TMatrixD::kInverted,m);
	TVectorD v(2,su), yz = mi*v; double y = yz[0], z = yz[1]; 
	double ry = y-y1, rz = z-z1, r2 = ry*ry+rz*rz;
	double dr2 = 2*sqrt(ry*ry*(mi(0,0)+dy2)+2*ry*rz*(mi(0,1)+dyz)+rz*rz*(mi(1,1)+dz2));
#ifdef FT2MA_HISTOS
	hMWAPAres2->Fill(r2/dr2,genuine+orig);
#endif
	if (r2/dr2<3) {
	  hsBest[jpl1] = h1; hsBest[jpl2] = h2; type = 1<<jpl1|1<<jpl2;
	}
      }
      else if (type!=0x7) type = -1;
    } // End loop on hits of 2nd plane
  } // End loop on hits of 1st plane
  if (type<=0) return 0; // No combination retained => exit

  //    ********** ADD, NEW, HITS FROM BEST COMBINATION and FIT **********
  int ndfs = t.NDFs; double chi2Tot = t.Chi2tot;
  for (int jpl = 0; jpl<3; jpl++) {
    if (!(1<<jpl&type)) continue; if (1<<jpl&associated)  continue;
    t.AnnexHit(*(hsBest[jpl]));
  }
  int added = t.NDFs-ndfs; if (!added) return 0; // New new hit added => exit
  bool ok = t.FullKF(1,-1,false);
  if (ok) {
    double chi2IperHit = // Max. allowed chi2 increment for picked-up hits.
      // (Nota Bene: this derivation of the cut value from option "dCut[94]" is
      // also done, independently, in tha body of "TEv::ForeTrack2MA".)
      TOpt::dCut[94] ? TOpt::dCut[94] : 9;
    ok = t.Chi2tot-chi2Tot<chi2IperHit*added;
  }
  if (ok) return added;                      // ***** FIT OK: RETURN #HITS ADDED
  //                                                // ***** ELSE: UNDO AND EXIT
  for (int jpl = 0; jpl<3; jpl++) {
    if (!(1<<jpl&type)) continue; if (1<<jpl&associated)  continue;
    t.CancelHit(ipli+jpl);
  }
  t.Clip(false /* i.e. upstream */);
  return 0;
}
void TEv::ForeTrack2Hs()
{
  // ********** EXTEND SAS-TRACKS into ZONE 0x8 **********

  // This is preliminary. A more appropriate "(Fore)Track2Hs" would also:
  // be able to:
  //  - deal w/ cases where the 1st hodo of the pair is missing (or both?),
  //  - select a set of hodo hits consistent w/ the trigger matrix among
  //   several candidates.
  // For the time being, we try to fix the folowing cases:
  //  - O-trigger w/ mu outside the acceptance of the tracking detectors
  //   downstream of the wall.
  //  - I-trigger w/ mu going through the central hole of PB's.
  //  - M- and L-triggers w/ poor detector efficiency in zone 0x8.
  // and we leave unattended the cases where only the 1st hodo is missing.
  const TSetup& setup = TSetup::Ref();

  //         *************** INITIALISATION ***************
  static int trigStations[4][2];
  static TH2D *hEffi[5], *hChi2[5], *hdChi2[5]; static bool doHist = false;
  static unsigned int hodoTrigs = 0, ladderTrig = 0, turnedOffTrigs = 0;
  static unsigned int imloPat;
  static bool first = true; if (first) {
    first = false;
    char IMLO[] = "IMLO0"; int imlo;
    
    //                      ********** TRIGGERS **********
    //               ***** GET THEIR DESCRIPTION from OPTIONs *****
    // - It's typically described in "./src/pkopt/trigger[.<extension>].opt".
    // - Expected are entries like 
    // Correspondance between hodoscope name and trigger mask.
    //			name  mask
    //Trigger mask	I	1
    //Trigger mask	M	2 // M covering both MT and InclMT.
    //etc...
    unsigned int trigWords[4] = {0,0,0,0};
    const map<unsigned int,char> &tcsMasks = CsInit::Instance()->getTCSMasks();
    map<unsigned int,char>::const_iterator im;
    for (im = tcsMasks.begin(); im!=tcsMasks.end(); im++) {
      unsigned int mask = (*im).first; char tag = (*im).second;
      for (imlo = 0; imlo<4; imlo++) {
	if (tag==IMLO[imlo]) {
	  hodoTrigs |= mask; trigWords[imlo] |= mask;
	}
      }
      if (tag=='L') ladderTrig = mask;
    }
    if (!isMC) { //    ***** XCHECK TRIGGERS AGAINST MAPPING *****
      unsigned int hTs = hodoTrigs, lT = ladderTrig; hodoTrigs=ladderTrig = 0;
      const CS::DaqOption &opts = CsInit::Instance()->getDaqEventsManager().GetDaqOptions();
      int bit = -1;// Bits 0x7 don't necessarily equate hodo triggers: check it
      try { bit = opts.GetTriggerBit("InnerTrigger"); } catch ( ... ) { }
      if (bit>=0) hodoTrigs |= 1<<bit;
      bit = -1;
      try { bit = opts.GetTriggerBit("MiddleTrigger"); } catch ( ... ) { }
      if (bit>=0) hodoTrigs |= 1<<bit;
      bit = -1;
      try { bit = opts.GetTriggerBit("MTLAST"); } catch ( ... ) { }
      if (bit>=0) hodoTrigs |= 1<<bit;
      bit = -1; // Ladder trigger: fill also dedicated "ladderTrig"
      try { bit = opts.GetTriggerBit("LadderTrigger"); } catch ( ... ) { }
      if (bit>=0) ladderTrig |= 1<<bit; hodoTrigs |= ladderTrig;
      bit = -1; // Incl. Ladder trigger: fill also dedicated "ladderTrig"
      try { bit = opts.GetTriggerBit("InclLadderTrigger"); } catch ( ... ) { }
      if (bit>=0) ladderTrig |= 1<<bit; hodoTrigs |= ladderTrig;
      bit = -1;
      try { bit = opts.GetTriggerBit("OuterTrigger"); } catch ( ... ) { }
      if (bit>=0) hodoTrigs |= 1<<bit;
      bit = -1;
      try { bit = opts.GetTriggerBit("OTLAST"); } catch ( ... ) { }
      if (bit>=0) hodoTrigs |= 1<<bit;
      bit = -1;
      try { bit = opts.GetTriggerBit("InclMiddleTrigger"); } catch ( ... ) { }
      if (bit>=0) hodoTrigs |= 1<<bit;
      bit = -1;
      try { bit = opts.GetTriggerBit("CalorimeterOuterTrigger"); } catch ( ... ) { }
      if (bit>=0) hodoTrigs |= 1<<bit;
      bit = -1;
      try { bit = opts.GetTriggerBit("JPsiTrigger"); } catch ( ... ) { }
      if (bit>=0) hodoTrigs |= 1<<bit;
      if (hodoTrigs!=hTs) CsErrLog::msg(elFatal,__FILE__,__LINE__,
	"Hodo triggers: mapping (=0x%x) != \"Trigger mask\" options (=0x%x)",
					hodoTrigs,hTs);
      if (ladderTrig!=lT) CsErrLog::msg(elFatal,__FILE__,__LINE__,
	"Ladder trigger: mapping (=0x%x) != \"Trigger mask\" options (=0x%x)",
					ladderTrig,lT);
    }
    //   ***** DETERMINE THE 2 TStation's ASSOCIATED WITH IMLO TRIGGERS *****
    for (imlo = 0, imloPat = 0; imlo<4; imlo++) {
      trigStations[imlo][0]=trigStations[imlo][1] = -1;
      for (im = tcsMasks.begin(); im!=tcsMasks.end(); im++)
	if ((*im).second==IMLO[imlo]) imloPat |= 1<<imlo;
    }
    for (int ist = 0; ist<(int)setup.vStation().size(); ist++) {
      int ipl = setup.vStation()[ist].IPlanes[0];
      const TPlane &p = setup.vPlane(ipl);
      if (ipl<setup.vIplFirst()[2]) continue;
      const TDetect &d = setup.vDetect(p.IDetRef);
      if (d.Name[0]!='H') continue;
      int i; for (i = 0, imlo = -1; i<4; i++)
	       if ((imloPat&1<<i) && d.Name[1]==IMLO[i]) imlo = i;
      if (imlo<0) continue;
      int i01; for (i = 0, i01 = -1; i<2; i++)
	if (trigStations[imlo][i]==-1) { i01 = i; break; }
      if (i01==-1)
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
		      "(TraF ReMode[33]) # of H%c stations > 2",IMLO[imlo]);
      trigStations[imlo][i01] = ist;
    }
    for (imlo = 0, imloPat = 0; imlo<4; imlo++) {
      int i01, turnedOff; for (i01 = 0, turnedOff = 1; i01<2; i01++) {
	int ist = trigStations[imlo][i01];
	if (ist==-1) CsErrLog::msg(elWarning,__FILE__,__LINE__,
	               "(TraF ReMode[33]): # of H%c stations < 2",IMLO[imlo]);
	else {
	  imloPat |= 1<<imlo;
	  const vector<int> &splanes = setup.vStation()[ist].IPlanes;
	  for (int spl = 0; spl<(int)splanes.size(); spl++)
	    if (setup.vPlane(splanes[spl]).IFlag) turnedOff = 0;
	}
      }
      if (turnedOff) {
	turnedOffTrigs |= trigWords[imlo]; hodoTrigs &= ~trigWords[imlo];
      }
    }
    if (TOpt::Hist[7]) {                                   // ***** HISTOS *****
      doHist = true;
      CsHistograms::SetCurrentPath("/Traffic/RDmonitor");
      char hName[] = "hHsdChi2I";
      for (imlo = 0; imlo<5; imlo++) {
	if (!(1<<imlo&imloPat) && imlo<4) continue;
	sprintf(hName,"hHsEffi%c",IMLO[imlo]);
	hEffi[imlo] =
	  new TH2D(hName,"Hodo tracking pseudo-#epsilon vs. #planes",
		   32,-.5,31,20,0,1);
	if (!isMC) continue;
	sprintf(hName,"hHsChi2%c",IMLO[imlo]);
	hChi2[imlo]  =
	  new TH2D(hName,"Hodo tracking perf vs. #chi^{2}",
		   10,0,TOpt::dCut[17],25,-1/48.,1+1/48.);
	sprintf(hName,"hHsdChi2%c",IMLO[imlo]);
	hdChi2[imlo] =
	  new TH2D(hName,"Hodo tracking perf vs. track's #Delta#chi^{2}",
		   10,-TOpt::dCut[17]/5,TOpt::dCut[17]/5,25,-1/48.,1+1/48.);
	CsHistograms::SetCurrentPath("/");
      }
    }
  }

  //       *************** PRELIMINARY STEPS ***************
  // - Reassess "THit::Status": set it =0 for those THit's free to be used, i.e.
  //  not added to TTrack themselves and not mirror of an added THit.
  //   "Thit::Status" has been set in "TEv::PrePattern2" and have been upset
  //  in "TEv::CleanTrackList" (and = 2, for the special case of suspicious
  //  0x8 track segments).
  //   (Secondarily, some of its settings have been assigned by "TSetup::Init".
  //  But these concern characteristics that are not at stake here: (as of
  //  05/08): coalescence ("Status==-4" must have prevented THit from being
  //  added to tracks and must have remained unchanged throughtout PR, and we
  //  keep it as is) and LR probability (the info is used only in PR, may have
  //  been overwritten (in PR, when THit is added) and we feel free to overwrite
  //  it again now).)
  //  ("THits::TStatus" will have to be taken care of throughout this method,
  //  for at the difference to "THit::sTrackID" it is not (and cannot easily) be
  //  updated by "TTrack::AddHit(and related methods)".)
  // => (among others) GEMs spacer info is erased. (This info singled out the
  //   THit's w/in the outer spacer interval, where the bias is relatively mild.
  //   The THit's w/in the much worse inner spacer interval were altogether
  //   disregarded, cf. "ImportClusters".)
  list<TTrack>::iterator it , it0x8;
  for (it = listTrack.begin(), it0x8 = listTrack.end(); it!=listTrack.end();
       it++) {
    // Type 0x8 tracks come last in list. To save CPU, save index of 1st one.
    if ((*it).Type==0x8 || (*it).Type==0xc) { it0x8 = it; break; }
  }
  int ipl, ipli = setup.vIplFirst()[2], iplf = setup.vIplLast()[3];
  for (ipl = ipli; ipl<=iplf; ipl++) {
    const TPlane &p = setup.vPlane(ipl);
    const TDetect &d = setup.vDetect(p.IDetRef);
    vector<int>::const_iterator ihit;
    for (ihit = p.vHitRef().begin(); ihit!=p.vHitRef().end(); ihit++) {
      THit &h = vecHit[*ihit];
      if (h.Status==-4) continue;// Exclude coalesced cluster's silent component
      int nIDs = h.sTrackID().size();
      if      (nIDs==0) h.status = 0;
      else if (nIDs==1) {       // Case one and only associated track...
	unsigned int iID = *(h.sTrackID().begin()); TTrack *t0x8;
	for (it = it0x8, t0x8 = 0; it!=listTrack.end(); it++) {
	  if ((*it).Id==iID) { t0x8 = &(*it); break; }
	}
	if (t0x8 &&
	    (t0x8->Type==0x8 ||  t0x8->Type==0xc)) {  // ...and it's 0x8|c...
	  if (t0x8->Chi2tot/(t0x8->NDFs-4)>TOpt::dCut[16]/2) t0x8->IMark = -2;
	  // ...if in addition it's suspicious (be it that it has been marked
	  // so by "TEv::TracksFit2" or has bad chi2)...
	  if (t0x8->IMark==-2) h.status = 2;// ...release (and yet tag) its hits
	}
	else h.status = 1;
      }
      else h.status = 1;
    }
    if (11<=d.IType && d.IType<=18) {       // Special case STs,DWs...
      // ...reassess status of mirrors
      for (ihit = p.vHitRef().begin(); ihit!=p.vHitRef().end(); ihit++) {
	THit &h = vecHit[*ihit]; int mirror;
	if (h.Status!=-4 && (mirror = h.Mirror)!=-1 &&
	    vecHit[mirror].Status>0) h.status = vecHit[mirror].Status;
      }
    }
  }
  const unsigned int allTrigs = 0xffff; // Although TCS can handle 0x7fffff, it was decided to limit max. #triggers to 16. And take advantage of the fact to make use of the other bits of the trigger TDC module.
  unsigned int evtTrig = trig_mask&allTrigs;     // ...cut away trailing bits
  //   ***** In view of UPDATING DRIFTS THits w/ EVENT TIME ***** 
  //      ***** REQUIRE THAT TRIGGER MATCHES ReMode[36]... *****
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
    evtT0 = eventTime; evtTConsidered = true;
  }
  bool doUpdateCs = (TOpt::ReMode[29]&0x4) &&
    // Do not convey the corrections to the CsCluster counterparts of the THits
    // if a redefinition of the event's timing is required and not granted.
    (!evtTimeRequired || evtTConsidered);
  bool doCorr = TOpt::ReMode[29]&0x2;

  for (it = listTrack.begin(); it!=listTrack.end(); it++) {
    TTrack &t = *it;
    //  We have in mind to fully extrapolate beyond track's end in zones 0xc and
    // will not provide for picking up hits (be they hodo hits) into an already
    // existing track piece. Therefore...
    // ...we require SM2 and reject muWall bridging...
    // ...EXCEPT if no large XX0 despite bridging (w/ a same cut on %X0 as in
    // "$CORAL/src/event/PIDdoMuonID.cc"), in order to still cope w/ cases
    // where the candidate mu is a ladder-triggered one but travels on the edge
    // of the absorber and gets a second chance to be ID'd if it traverses the
    // thickness of material in front of HI05 and gets a hit there, cf.
    // D*.2006.03, #89. Theses cases cover a very small region of phase space
    // but also a region that is difficult to describe in MC: a typical edge
    // effect situation. The second heat w/ the inner provides a way to work
    // around the edge effect problem: it cancels the edge, so to say. It's
    // hence interesting to take it into account. We restrict ourselves to
    // ladder. We do not considered hadron triggers, given that the muID is not
    // that important there that we have to bother w/ small phase space regions.
    int ladderEdge =
      (t.Type&0xe)==0xe && t.RadLenFraction()<30 && (evtTrig&ladderTrig);
    if ((t.Type&0xe)!=0x6 && !ladderEdge) continue;

    //          ***** NO TIME CUT: 
    // - We are mainly interested in the scattered muon. I.e. tracks in coinc.
    //  w/ the primary interaction. Therefore, we could cut on track time w.r.t.
    //  trigger or w.r.t. event time. E.g.:
    //double evtT = eventTRef ? eventTRef : eventTime;
    //double extraTime = ptrEvt()->getExtraTimeWidth();
    //if (t.SigmaTime<0 ||
    //    evtT  && fabs(t.MeanTime-evtT)>3*t.SigmaTime+1 ||
    //    !evtT && fabs(t.MeanTime)>3*t.SigmaTime+extraTime) continue;
    // - Yet, since this "ForeTrack2Hs" turns out no to be CPU intensive, let's
    //  accept all tracks => We'll thus ID more muons among the halo tracks.
    //  (Note that, along the same line, we could also drop the requirement on
    //  the hodo hits infra.)
#ifdef DEBUG
    static int idebug = 0;
    if (idebug) it->DumpHits(idebug);
#endif

    // *************** LOOP on ZONES&=0x6 IN-TIME TRACKS ***************

    float x0 = t.Hlast(0);  // Intercepts and slopes taken from last helix
    float y0 = t.Hlast(1), z0 = t.Hlast(2);
    float yp = t.Hlast(3), zp = t.Hlast(4);
    static float dy2, dyz, dz2;
    list<THit*> hitList; int hodoOK, nHodos, nHodosInList, imlo; THlx hlx; 

    hodoOK = 0; nHodos = 0; nHodosInList = 0;
    if (evtTrig&hodoTrigs) {

      // ********** HODO TRIGGERS: SEARCH for CONSISTENT IMLO HITS **********

      for (imlo = 0; imlo<4; imlo++) {
	if (!(1<<imlo&imloPat)) continue;
	bool isInTrig = (1<<imlo&evtTrig) ||
	  imlo==1 && (evtTrig&0x100) ||  // Inclusive middle
	  imlo==3 && (evtTrig&0x600);    // Inclusive outer (and J/psi?)

	// ***** IMLO CONSISTENT w/ TRIGGER *****
	// Comment out the following instruction: it lowers significantly the #
	// of successfully extended tracks, while not bringing any significant
	// improvements otherwise.
	//if (!isInTrig) continue;

	if ((1<<imlo&0x9) &&        // Inner|Outer-Triggers:
	    (TOpt::ReMode[33]&0x2) && // HI04|HO03 is prerequisite, cf.
	  !(t.Scifi&0x10000<<imlo))   // preselected tracks in "TEv::TracksFit2"
	  continue;
	for (int i01 = 0; i01<2; i01++) {
	  if (i01==0 && (TOpt::ReMode[33]&0x2)
	      && (1<<imlo&0x9))// IO hits already there, cf. preselection
	    continue;
	  const TStation &s = setup.vStation()[trigStations[imlo][i01]];
	  vector<int>::const_iterator ip; bool first;
	  for (ip = s.IPlanes.begin(), first = true; ip!=s.IPlanes.end();
	       ip++) {
	    // ********** LOOP on PLANES in 1ST/2ND HODO TStation **********
	    const TPlane &p = setup.vPlane(*ip);
	    if (!p.IFlag) continue;                      // Plane OFF
	    if (p.vHitRef().empty()) continue;           // Plane empty
	    const TDetect &d = setup.vDetect(p.IDetRef); float x = d.X(0);
	    float y = y0+yp*(x-x0), z = z0+zp*(x-x0);
	    if (!d.InActive(y,z)) continue;       // ***** W/IN ACCEPTANCE *****
#define FT2Hs_CHECK_URANGE
#if defined FT2Hs_CHECK_URANGE && !defined InActive_CHECK_URANGE
	    if (!d.WinRange(y,z)) continue;       // ***** W/IN RANGE *****
#endif
	    if (first) {    // ***** 1st encountered plane: UNCERTAINTIES? *****
	      hlx(0) = x; t.Hlast.Extrapolate(hlx,true);
	      dy2 = hlx(1,1); dz2 = hlx(2,2); dyz = hlx(1,2); first = false;
	    }
	    float sU = d.Sa, cU = d.Ca, wU = d.Resol;// ***** LOOK for HIT *****
	    float u = y*cU+z*sU, du = sqrt(dy2*cU*cU+2*dyz*cU*sU+dz2*sU*sU);
	    float umx = u+3*du+3*wU, umn = u-3*du-3*wU;
	    vector<int>::const_iterator ihit;
	    THit *hBest; static float best, diff;
	    for (ihit = p.vHitRef().begin(), hBest = 0; ihit!=p.vHitRef().end();
		 ihit++) {
	      THit &h = vecHit[*ihit];
	      float U = h.U;
	      if (U<umn) continue; if (U>umx) break; diff = fabs(U-u);
	      if (!hBest || diff<best) { hBest = &h; best = diff; }
	    }
	    if (hBest) {
	      set<unsigned int>::iterator id; int yet;
	      for (id = hBest->sTrackID().begin(), yet = 0;
		   id!=hBest->sTrackID().end(); id++) {
		if (*id==t.Id) { yet = 1; break; }
	      }
	      if (!yet) {
		hitList.push_back(hBest); nHodosInList++;
	      }
	      nHodos++;
	    }
	  }  // End loop on hodo planes in current station
	}  // End loop on the 2 hodo component stations of current IMLO trigger
	if (nHodos>=4 ||            // HM = HM04/5X1 and HM04/5Y1
	    imlo!=1 && nHodos>=2 || // HL: require HL04 AND HL05
	    (1<<imlo&0x9) && nHodos)
	  hodoOK |= 1<<imlo;
      }  // Loop on IMLO triggers
      if (!hodoOK) continue;
    }
    // else: Case no hodo consistency check (e.g. hadron data)...
    // ...=> all tracks are scanned.

    // ********** HODO FOUND: CHECK DETECTORS ON THE WAY, AND BEYOND **********
    // (They must have fired w/ a reasonably good effciency.)
    list<int>::const_iterator ih = t.lHPat().end(); ih--;
    list<int>::const_iterator ip = t.lPRef().end(); ip--;
    static int ipl0; while (ip!=t.lPRef().begin()) {
      if (*ih>=0) { ipl0 = *ip; break; }
      ih--; ip--;
    }
    THlx hlxp(t.Hlast); const TStation *sPrv;
    int nExpected, nFound, nStations, projs, nProjs;
    for (ipl = ipl0+1, nExpected=nFound=nStations=projs=nProjs = 0, sPrv = 0;
	 ipl<=iplf; ipl++) {
      //          ***** LOOP on PLANES from END of TRACK to END of SPECTRO *****
      const TPlane &p = setup.vPlane(ipl); if (p.IFlag==0) continue; // OFF
      const TDetect &d = setup.vDetect(p.IDetRef); float x = d.X(0);
      if (hodoOK!=0 && d.IType>=41 &&
	// This is to exclude hodos, if already scanned. But has to be bypass in
	// the special case of a "ladderEdge" looking at HI05.
	  (!ladderEdge || d.Name.find("HI05"))) continue;
      float y = y0+yp*(x-x0), z = z0+zp*(x-x0);
      if (!d.InActive(y,z)) continue;             // ***** W/IN ACCEPTANCE *****
#if defined FT2Hs_CHECK_URANGE && !defined InActive_CHECK_URANGE
      if (!d.WinRange(y,z)) continue;             // ***** W/IN RANGE *****
#endif
      nExpected++;
      const TStation *&s = p.Station; if (s!=sPrv) {
	//                  ***** NEW TStation: COMPUTE UNCERTAINTIES ANEW *****
	hlx(0) = x; hlxp.Extrapolate(hlx,true);
	dy2 = hlx(1,1); dz2 = hlx(2,2); dyz = hlx(1,2); hlxp = hlx;
      }
      float sU = d.Sa, cU = d.Ca, wU = d.Resol;
      float u = y*cU+z*sU, du = sqrt(dy2*cU*cU+2*dyz*cU*sU+dz2*sU*sU);
      float umx = u+3*du+4*wU, umn = u-3*du-4*wU;
      double incidence = d.Ca*yp+d.Sa*zp; THit *m;
      vector<int>::const_iterator ihit; THit *hBest; static float best, diff;
      for (ihit = p.vHitRef().begin(), hBest=m = 0; ihit!=p.vHitRef().end();
	   ihit++) {                      // ***** LOOK for BEST THit *****
	THit &h = vecHit[*ihit];
	//if (h.Status==1 || h.Status==3 || h.Status==-4) continue;
	if (h.Status==1 || h.Status==3) continue;
	if (&h==m) continue;
	float U = h.U; if (U<umn) continue; if (U>umx) break;
	diff = fabs(U-u); if (hBest && diff>=best) continue;
	int swap = 0, mirror; if ((mirror = h.Mirror)!=-1) {
	  m = &vecHit[mirror];      // ***** IF MIRROR: CONSIDER LR PROBABILITY
	  if (m->Status!=-4 && (!hBest || fabs(m->U-u)<best)) {
	    CsCluster *ci = h.ptrCl, *cj = m->ptrCl;
	    list<CsCluster*> cs = ci->getAssociateClusters(); if (!cs.empty()) {
	      CsEventUtils::setLRToMatchAngle(ci,cs.front(),incidence);
	      if (cj->getLRProb()>ci->getLRProb()*1.1) {
		hBest = m; diff = fabs(m->U-u); swap = 1;
	      }
	    }
	  }
	}
	if (!swap) hBest = &h; best = diff;
      }
      if (hBest) {
	if (d.IType>=41) {     // ***** CASE HODOs: CHECK TIME CONSISTENCY *****
	  // (This condition has been added lately and its impact hasn't been
	  // fully checked. But encountering hodos @ this point is only possible
	  // if "!hodoOK" or "ladderEdge", so that impact is limited to
	  // C-triggers (or ladder or hadron data). In C-triggers case, this
	  // impact is 2-fold:
	  //  i) Add HI05 hit to track, and therefore grant it large %X0.
	  // ii) Allow track to pass "nStations>=2" condition infra. Again this
	  //    concerns (most probably) only HI05, w/ track firing, apart from
	  //    hodo, only GM11 or PBs (close to their dead zone).)
#ifdef ForeT_USE_HODO_TRESOL
	  float hdT = hBest->SigT;
#else
	  // Use built-in time resolution for hodo, for present (as of 05/12)
	  // one is large (1.5ns) and does not provide a strict enough time
	  // cut. 1.5ns is probably over-estimated for HM,HL and even HO(?), or,
	  // for HI does not correspond to the, yet to be computed, time
	  // corrected for light propagation. We can only move on to use the
	  // THit::SigT when the 2 points supra are clarified.
	  float hdT = .75;
#endif
	  float tdT = t.SigmaTime, dT = sqrt(tdT*tdT+hdT*hdT);
	  if (tdT>0 && fabs(t.MeanTime-hBest->Time)>2*dT && // Strict time cut (
	    // 2 sigmas only), since we're here in kind of a rescue procedure...
	      dT<1) // ...and add a condition of overall time resolution, so
	    // that adding this hodo hit (which will get "nStations" infra
	    // incremented) be reliable enough. (We're mostly looking at HI,
	    // so that this latter condition should not be too much demanding.)
	    continue;
	}
	hitList.push_back(hBest); nFound++;
	if (s!=sPrv) nStations++;
	int jproj = 1<<p.IProj;
	if ((jproj&projs)==0) { projs |= jproj; nProjs++; }
      }
      sPrv = s;
    }
    if (!hodoOK &&    // ***** CASE NO CONSISTENCY W/ TRIGGER REQUIRED ... *****
	// (This must correspond to "caloTrig")
	// ... => None of the reliability brought by hodo tight timing
	(nProjs<3 || nStations<2) &&// ...=> TIGHT REQUIREMENTS on FOUND HITS...
	// ...except in the case one of the current hodo-based triggers has
	// been DetNameOff'd: we are then typically doing detector studies,
	// hence accept to trade off some reliability for a more extensive phase
	// space coverage.
	(!(evtTrig&turnedOffTrigs) || nProjs<2))
      continue;
    if (doHist && nExpected) {
      if (hodoOK) {
	for (imlo = 0; imlo<4; imlo++) if (hodoOK&1<<imlo)
	  hEffi[imlo]->Fill(nExpected,(double)nFound/nExpected);
      }
      else hEffi[4]->Fill(nExpected,(double)nFound/nExpected);
    }
    if (hitList.empty()) continue;
    if (2<nExpected &&  // I.e. so few expected hits indicate that we are either
	// close to some boundary of the acceptance, or that the track already
	// came close to the end of zone 0x4 and missed the last few hits
	// because they are simply not there => disregard the efficiency cut.
	nFound<.5*nExpected) {                     // ***** EFFICIENCY CUT *****
      const TDetect &first = setup.vDetect(ipl0);
      const TDetect &last  = setup.vDetect(iplf);
      if (nHodosInList /* and implictly !hitList.empty() */) {
	const TDetect &hodo  = setup.vDetect(hitList.front()->IPlane);
	if (nHodosInList==1)
	  CsErrLog::msg(elError,__FILE__,__LINE__,
 "Cancel adding %s to TTrack #%d: extension %s->%s w/ %d Found/%d Expected",
  hodo.Name.c_str(),t.Id,first.Name.c_str(),last.Name.c_str(),nFound,nExpected);
	else {
	  list<THit*>::iterator ihit = hitList.begin();
	  for (int i = 1; i<nHodosInList; i++) ihit++;
	  const TDetect &hodop  = setup.vDetect((*ihit)->IPlane);
	  CsErrLog::msg(elError,__FILE__,__LINE__,
 "Cancel adding %s..%s to TTrack #%d: extension %s->%s w/ %d Found/%d Expected",
  hodo.Name.c_str(),hodop.Name.c_str(),t.Id,first.Name.c_str(),
			last.Name.c_str(),nFound,nExpected);
	}
      }
      else
	CsErrLog::msg(elError,__FILE__,__LINE__,
 "Cancel extending to TTrack #%d into zone 0x8: %s->%s w/ %d Found/%d Expected",
	  t.Id,first.Name.c_str(),last.Name.c_str(),nFound,nExpected);
      continue;
    }

    bool fillMCHist = doHist && isMC; int nGoodHits = 0;
    if (fillMCHist) {
      t.FindKine(); fillMCHist &= t.IKine>=0;
    }

    // It's a halo muon accidentally coincident?
    double evtT = eventTRef ? eventTRef : eventTime;
    double tT0; if (fabs(t.MeanTime-evtT)/t.SigmaTime>5./3 &&     // Off time...
		    // Note: When time undef: SigmaTime<0
		    fabs(t.Hlast(3))<.05 && fabs(t.Hlast(4))<.05) // Paraxial...
      tT0 = t.MeanTime-eventTRef;// ...Off-time halo muon: Set track's time "tT"
    else tT0 = evtT0;

    list<THit*>::iterator ihit;
    for (ihit = hitList.begin(); ihit!=hitList.end(); ihit++) {
      THit *h = *ihit;
      if (fillMCHist && (h->IKine==t.IKine || h->IKin2==t.IKine)) nGoodHits++;
      if (doCorr) {
	//   ***** TIME PROPAGATION, TRIGGER JITTER, TRACK'S TIME OFFSET *****
	const TPlane &p = setup.vPlane(h->IPlane);
	const TDetect &d = setup.vDetect(p.IDetRef);
	//   N.B.: we need only a crude estimate => non updated
	// values "yp" and "zp" are good enough for that purpose.
	double y = (y0+yp*(d.X(0)-x0))*10, z = (z0+zp*(d.X(0)-x0))*10;
	CsDetector *csDet = d.PtrDet(); if (csDet->hasDrift()) {
	  bool error; double U = csDet->getCorrU(h->PtrClus(),y,z,tT0,error);
	  if (!error) {
	    h->u = U/10; h->Rotate();
	    if (doUpdateCs)                       // ***** RESET CsCluster *****
	      // (This is only useful when one intends to make use of
	      // these clusters outside of TraFDic, e.g.:
	      //  - TraFDic monitoring when TraFFiC modularity != 2
	      //  - Evaluation of residuals in coral, as opposed to TraFDiC)
	      h->PtrClus()->setU(U);
	  }
	}
      }
      t.AnnexHit(*h);                                    // ***** ADD HITS *****
    }
    t.InsertMissedPlanes();
    THlx oldH = t.Hfirst; double oldChi2 = t.Chi2tot/(t.NDFs-5);
    bool ok = t.FullKF(1); if (ok) {                        // ***** REFIT *****
      map<int,double> mChi2 = t.mChi2;
      ok = t.FullKF(-1); double chi2 = t.Chi2tot/(t.NDFs-5); THit *worstHit = 0;
      if (TOpt::ReMode[38] && ok && chi2>oldChi2+TOpt::dCut[17]/10) {
 	//                         ***** CLEANING of TRACKS w/ DRIFTS HITS *****
	double chi2IncrMx = TOpt::dCut[17]/30;// ***** BAD chi2 INCREMENTS *****
	vector<THit*> badHits, altHits; double wIncr;
	map<int,THlx>::const_iterator im;
	list<int>::const_iterator ih; map<int, double>::iterator idx, jdx;
	for (ih = t.lHPat().begin(), idx = t.mChi2.begin(), jdx = mChi2.begin(),
	       wIncr = 0; ih!=t.lHPat().end(); ih++) {
	  if (*ih<0) continue;
	  THit &h = vecHit[*ih]; if (h.IPlane<=ipl0) { idx++; jdx++; continue; }
	  double chi2Incr = (*idx).second, chi2Jncr = (*jdx).second;
	  if      (chi2Incr>wIncr) { wIncr = chi2Incr; worstHit = &h; }
	  else if (chi2Jncr>wIncr) { wIncr = chi2Jncr; worstHit = &h; }
	  if (chi2Incr>chi2IncrMx) {
	    int mirror; if ((mirror = h.Mirror)!=-1 &&
			    vecHit[mirror].Status!=-4) {
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
	      if (fabs(dum)<fabs(du)) {
		badHits.push_back(&h); altHits.push_back(&hm);
	      }
	    }
	  }
	  idx++; jdx++;
	}
	int status = 0; if (badHits.size()) {
	  //     ********** TRY FIRST ALTERNATIVE HITS LIST... **********
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
	  static double newChi2; if ((ok = t.FullKF(1)) &&  // ***** REFIT *****
				     (newChi2 = t.Chi2tot/(t.NDFs-5))<chi2) {
	    ok = t.FullKF(-1); newChi2 = t.Chi2tot/(t.NDFs-5);
	  }
	  if (!ok || chi2<newChi2) {          // ***** NOT BETTER: RESTORE *****
	    for (int i = 0; i<(int)badHits.size(); i++) {
	      THit *h = badHits[i]; t.ReplaceHit(h);
	    }
	    t.Hfirst = oldH;
	    ok = t.FullKF(1); if (ok) ok = t.FullKF(-1);
	  }
	  else status |= 0x4;
	}  // End bad hits found
 	if (doUpdateCs && (status&0x4)) {           // ...if requested...
	  //                                     ***** ...UPDATE CsCluster *****
	  for (int i = 0; i<(int)badHits.size(); i++) {
	    THit *hm = altHits[i]; hm->PtrClus()->setU(hm->U*10);
	  }
	}
	if (ok) ok = chi2<TOpt::dCut[17] && chi2<oldChi2+TOpt::dCut[17]/10;
	if (!ok && worstHit && (nFound-1)>.5*nExpected && wIncr/t.Chi2tot>.4) {
	  //     ***** ...IF HOPELESS STILL, TRY THEN ERASING WORST HIT *****
	  t.CancelHit_Clip(worstHit->IPlane);
	  ok = t.FullKF(1); if (ok) ok = t.FullKF(-1);
	  chi2 = t.Chi2tot/(t.NDFs-5);
	  // We delay the updating "hitList" until it's necessary, i.e.
	  // if the 0x8 extension is validated. Otherwise it's better to keep
	  // the original list, in order to later "TTrack::Shorten" the track.
	}
	else worstHit = 0;
      }  // End try and cure all bad hits
      if (ok) {
	ok = chi2<TOpt::dCut[17] && chi2<oldChi2+TOpt::dCut[17]/10;
	if (worstHit) {
	  for (ihit = hitList.begin(); ihit!=hitList.end(); ihit++) {
	    THit *h = *ihit; if (h==worstHit) { hitList.erase(ihit); break; }
	  }
	}
      }
      if (fillMCHist) {
	if (ok) {
	  if (hodoOK) {
	    for (imlo = 0; imlo<4; imlo++) if (hodoOK&1<<imlo) {
	      hChi2 [imlo]->Fill(oldChi2,(double)nGoodHits/nFound);
	      hdChi2[imlo]->Fill(chi2-oldChi2,(double)nGoodHits/nFound);
	    }
	  }
	  else {
	      hChi2 [4]->Fill(oldChi2,(double)nGoodHits/nFound);
	      hdChi2[4]->Fill(chi2-oldChi2,(double)nGoodHits/nFound);
	  }	    
	  if (t.NHsame+nGoodHits<TOpt::dCut[10]*t.NHits)
	    CsErrLog::msg(elError,__FILE__,__LINE__,
 "Ttrack ID=%d,MC=%d,%d/%d downgraded => %d/%d",t.Id,
	      t.IKine,t.NHsame,nFound+nHodos,t.NHsame+nGoodHits,t.NHits);
	}	  
      }
    }
    if (ok) {// ********** FIT OK: UPDATE THit' Status, TTrack's Type **********
      for (ihit = hitList.begin(); ihit!=hitList.end(); ihit++) {
	THit *h = *ihit;
	int status = // Grant special status=3 to hits previously associated to
	  // 0x8 track segments, so as to tag the latter for deletion (and
	  // accessorily prevent the hits from being associated again).
	  h->Status==2 ? 3 : 1;
	h->status = status;
	int mirror; if ((mirror = h->Mirror)!=-1 && vecHit[mirror].Status!=-4)
	  vecHit[mirror].status = status;
      }
      t.Type |= 0x8;
    }
    else {    // ********** FIT FAILS => RESTORE INITIAL TRACK **********
      //t.Shorten(hitList.front()->IPlane);
      for (ihit = hitList.begin(); ihit!=hitList.end(); ihit++) {
	THit *h = *ihit; t.CancelHit(h->IPlane);
      }
      t.Clip(true); // Erase trailing plane references
      t.Hfirst = oldH; ok = t.FullKF(1); if (ok) ok = t.FullKF(-1);
      if (!ok) {
	CsErrLog::msg(elError,__FILE__,__LINE__,
 "Full KF fails after initial track (ID=%d,chi2=%f) is restored",t.Id,oldChi2);
	listTrack.erase(it++); continue;
      }
    }
#ifdef DEBUG
    if (idebug) {
      t.DumpHits(idebug); if (idebug>1) t.Print(idebug-1);
    }
#endif
  }
  //          ***** CLEAN 0x8-TRACK SEGMENTS *****
  // The foretracking to Hs allows itself to re-use the hits (which it tags w/
  // status=2) of suspicious 0x8 tracks, cf. supra. =>
  //  - Erase those 0x8-tracks w/ re-used hits (they must have been tagged w/
  //   status=3). (And don't try to rescue them by only stripping them of the
  //   re-used hits: they're not that important.)
  //    Update the status of their hits: 3->1, 2->0, 1->1.
  //  - For the rest: update the status of their hits: 3->1, 2->1, 1->1. 
  it = listTrack.begin(); while (it!=listTrack.end()) {
    TTrack &t = *it;
    if (t.IMark!=-2) { it++; continue; }
    //if (t.Type!=0x8) { it++; continue; }
    list<int>::iterator ihp; int doErase;
    for (ihp = t.lHitPat.begin(), doErase = 0; ihp!=t.lHitPat.end(); ihp++) {
      if (*ihp<0) continue; THit &h = vecHit[*ihp];
      if (h.Status==3) { doErase = 1; break; }
      int mirror; if ((mirror = h.Mirror)!=-1 && vecHit[mirror].Status==3) {
	doErase = 1; break;
      }
    }
    for (ihp = t.lHitPat.begin(); ihp!=t.lHitPat.end(); ihp++) {
      if (*ihp<0) continue; THit &h = vecHit[*ihp];
      int status; if (h.Status==2) status = doErase ? 0 : 1;
      else                         status = 1; 
      h.status = status;
      int mirror; if ((mirror = h.Mirror)!=-1 && vecHit[mirror].Status!=-4)
	vecHit[mirror].status = status;
    }
    if (doErase) listTrack.erase(it++);
    else { t.IMark = 0; it++; }
  }
}
