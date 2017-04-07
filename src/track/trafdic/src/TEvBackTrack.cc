// $Id: TEvBackTrack.cc 14069 2015-09-17 20:44:46Z lsilva $

#include "CsInit.h"
#include "DaqDataDecoding/DaqOption.h"
#include "CsEvent.h"
#include "CsDriftChamberDetector.h"
#include "Traffic.h"
#include "TOpt.h"
#include "TSetup.h"
#include "TEv.h"
#include "TLattice.h"
// The following headers have to be included last, for they conflict
// w/ "CsChip.h"
#include "TMatrixF.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TH1.h"
#include "TH2.h"

using namespace std;

/*! 
  \brief Various methods for extending TTrack's by backtracking.
  \author  Yann.Bedfer@cern.ch
*/

/*
  - BackTrackZ2 : Extend 0x4 tracks into zone 0x2.
  - BackTrackSAS: Extend 0x6 tracks into zone 0x1.
  - BackTrackZ1 : Extend 0x2 tracks into zone 0x1.
  - BackTrack2SI:
  - BackTrackVD : Extend &0x3 tracks to the vertex detector
  + Associated methods
  - BuildSpacePt: Build a space point on an argument Y hit, based on a an
                 argument direction (specified by an helix).
  - ...
 */

// Interface
bool Account4RIPipe(TTrack &t);

// ===========================================================================
// ==============================  BackTrackZ2  ==============================
// ===========================================================================
void TEv::BackTrackZ2()
{
  // BackTrackZ2 : Extend 0x4 tracks into zone 0x2.
  // - Backtracking is along a trajectory fitted to the 0x4 track, w/ a
  //  ``vertex'' constraint based on available incident beam tracks.
  // - Then pick up hits from zone 0x2.

  // *************** EXTEND SM2<->muW TRACKS UPSTREAM of SM2 ***************
  // *************** ASSUMING VERTEX @ BEAM INTERCEPT w/ X=0 ***************

  const TSetup &setup = TSetup::Ref();
  // ***** NOTATION:
  // Parameters of zone 0, i.e. target<-> SM1: "_z0"
  // Parameters of zone 1, i.e. SM1<-> SM2:    "_z1"
  int ipli_z0 = setup.vIplFirst()[0];
  int ipli_z1 = setup.vIplFirst()[1], iplf_z1 = setup.vIplLast()[1]; 
  int nPats = (iplf_z1-ipli_z1)/32+1; unsigned int hitsPats[nPats];
  double Xlow = 580;    // Lower limit of X range = FI05 ~= 580

  //    ***** CHECK NOT TOO MANY HITS and... *****
  //       ...AND DETERMINE some PARAMETERS of the PROBLEM
  int iplow, nhits_z1, npls_z1;
  for (iplow = iplf_z1, npls_z1=nhits_z1 = 0; iplow>=ipli_z1;
       iplow--) {
    const TPlane &p = setup.vPlane(iplow); if (!p.IFlag) continue; // Plane OFF
    const TDetect &d = setup.vDetect(p.IDetRef);
    if (11<=d.IType && d.IType<=18) continue;// Skip drift-like (ST04/DR>"Xlow")
    if (39<=d.IType && d.IType<=40) continue;// Skip MAs
    if (d.X(0)<Xlow) break;
    npls_z1++; nhits_z1 += p.vHitRef().size();
  }
  iplow++; // since loop supra went 1 plane too far
  const TLattice *lat = TLattice::Ptr();
  if (!lat || nhits_z1>8*npls_z1) return;

  // ***** BEAM: NO MORE THAN 2 => FIXED Y VERTEX for subsequent fit *****
  TTrack *beams[2]; int nBeams; GetIncidentTracks(2,nBeams,beams);
  if (nBeams== 0 || nBeams>2) return;

  //     ***** FREE (temporarily) HITS from BACHELOR 0x2 TRACKS *****
  // - All unsued hits from zone 0x2 are considered.
  // - To which we would like to add ``badly'' used ones. Along this line, let's
  //  free 0x2 hits associated to unreliable tracks. Are considered unreliable:
  //   - Unbridged => 0x2 tracks.
  //   - Few hits. Or rather not enough constraining set of hits: too few, or
  //    too few effective space points, as in the case of a (FI05X+FI06Y)+
  //    (GP02P1+GP02P2) track, where there are in fact 2.
  //    (But be careful w/ scifi tracks. For the time being at least... Scifi
  //    tracks have low redundancy but also have this extra feature, that adds
  //    to their reliability, that they have a tight time consistency built in.
  //    On the other hand, a 0x2 scifi track, even if it can later be extended
  //    to zone 0x1, is not very useful: scifi tracks are rather expected to be
  //    high momenta spanning all 0x7 zones.)
  //  => Let's flag such tracks w/ TTrack::IMark = -2.
  list<TTrack>::iterator it;
  for (it = listTrack.begin(); it!=listTrack.end(); it++) {
    TTrack &t = *it; if (t.Type!=0x2) continue;
    //  Track could eventually be erased (if any of its hits turns out
    // incorporated to a 0x6 track). Which can pose a problem if it's associated
    // (Note: In principle, as of 10/02, 0x2-tracks can not be associated. Yet,
    // to be on the safe side...)
    //  Association can be flagged by TTrack::Associate or ::IMark. For the
    // latter, require that it be strictly =0, in order to avoid overwriting any
    // other special setting.
    if (t.Associate>=0 || t.IMark!=0) continue;
    if (t.NHits>6) continue; // 6 hits is in any case well enough constraining.
    if ((t.Scifi&0x2) && t.NHits>4) continue; // Be careful w/ scifi tracks...
    list<int>::iterator ih; int nSpacePts = 0; const TStation *sPrv = 0;
    int projs = 0, nProjs = 0, allProjs = 0, nAllProjs = 0;
    for (ih = t.lHitPat.begin(); ih!=t.lHitPat.end(); ih++) {
      if (*ih<0) continue;
      const THit &h = vecHit[*ih]; const TPlane &p = setup.vPlane(h.IPlane);
      const TStation *&s = p.Station; if (s!=sPrv) {
	if (nProjs>=3 || nProjs>=2 && sPrv->Type==22) {
	  nSpacePts++; nProjs=projs = 0;
	}
	sPrv = s;
      }
      int jproj = 1<<p.IProj;
      if ((jproj&projs)==0)         {    projs |= jproj;    nProjs++; }
      if ((p.IFlag&0x30) && (0x2&projs)==0)
	// (Note: This excludes pixel MPs: they don't provide accurate v-info.)
	/* Pixel plane: Z-proj. */  {    projs |= 0x2;      nProjs++; }
      if ((jproj&allProjs)==0)      { allProjs |= jproj; nAllProjs++; }
    }
    if (nProjs>=2 || nProjs>=2 && sPrv && sPrv->Type==22) nSpacePts++;
    if (nAllProjs>=3 && nSpacePts>=3) continue;
    t.IMark = -2; // Flag (temporarily: may be reset in fine) "t" as unreliable.
  }

  double fldint = setup.MagFieldInt[setup.NMags-2]+
    setup.MagFieldInt[setup.NMags-1];

  // When it comes to pick up a scifi hit, we apply a time cut on the time diff
  // of the hit w.r.t. the time of the event.
  // When drift CsCluster's are updated (by "UpdateDrifts"), "eventTime" is set
  // = 0, while the actual value is backed-up into "eventTRef".
  double evtT = eventTRef ? eventTRef : eventTime;

  it = listTrack.begin(); while (it!=listTrack.end()) {
    TTrack &t = *it; if ((t.Type&0x7)!=0x4) { it++; continue; }
    if (t.IMark>0) { it++; continue; } // Bypass 0x4 piece of yoke track

    // ******************** LOOP ON SM2<->muWall TRACKS ********************

    float nSpacePts = 0; int added = 0, nScifis = 0; bool bendingInfo = false;
#ifdef DEBUG
    static int idebug = 0; if (idebug) t.DumpHits(0);
#endif
    double oldChi2 = t.Chi2tot; // Save KF's "Chi2tot" before it's upset by QN
    int iB, iBestB, nBs = nBeams==1 ? 1 : 3; float bestB;
    for (iB = 0, iBestB = -1, bestB = 500; iB<nBs; iB++) {
      // ***** EXTRAPOLATE BEAM to Z beam == Z track *****
      // (do this crudely (i.e. straight) since x=0 already crude assumption)
      TTrack *beam;
      if      (iB<2)      beam = beams[iB];
      else if (iBestB>=0) beam = beams[iBestB];
      else break;
      float y0B = beam->Hfirst(1)-beam->Hfirst(3)*beam->Hfirst(0);
      float z0B = beam->Hfirst(2)-beam->Hfirst(4)*beam->Hfirst(0);
      float ypB = beam->Hfirst(3), zpB = beam->Hfirst(4);
      float xB = (t.Hfirst(2)-t.Hfirst(4)*t.Hfirst(0)-z0B)/(zpB-t.Hfirst(4));
      t.Haux(5) = 	// Guess value from horizontal angle
	(0-t.Hfirst(3))/(1.E-3 * 0.3 * fldint);
      // Let's not fix the vertex Y @ X=xB, for we do not have yet a reliable
      // estimate of the Y angle @ that abscissa (and hence cannot determine,
      // w/in "QNewtonFit", the corresponding value of the Y track parameter
      // @ the origin of the dico. This could be done by iterating on
      // "QNewtonFit", but propably is not worth the effort, given that we
      // we are trying to pick-up hits in a region that is far away from the
      // target (and therefore the vertex) and hence are unsensitive to xB.
      t.Haux(0) = 0; t.Haux(1) = y0B; t.Haux(2) = z0B;
      t.Haux(3) = (t.Hfirst(3)-y0B)/t.Hfirst(0); t.Haux(4)= t.Hfirst(4); 
      if (t.QNewtonFit(1,6) &&
	  // Check resulting chi2 against "upper chi2/NDF cut conditioning QNFit
	  // for BackTracking" ("dCUt[89]"). The cut is not exactly dedicated to
	  // the case, since it's also used in similar (upper cut on QNFit, for
	  // BackTracking) yet distinct (different QN mode) contexts.
	  t.Chi2tot/(t.NDics-(t.NDics>4 ? 4 : 3))<TOpt::dCut[89]) {
	if (iB<2) {
	  if (iBestB==-1 || fabs(xB)<bestB) {
	    iBestB = iB; bestB = fabs(xB);
	  }
	}
      }
      else if (iB==2) iBestB = -1;
    }
    if (iBestB==-1) { it++; continue; }

    vector<TTrack*> tracksToErase;
    int iter; for (iter = 0, memset(hitsPats,0,nPats*sizeof(unsigned int));
		   iter<2; iter++) {

      // *************** 2 ITERATIONS ***************

      if (t.Chi2tot==1e6) break;    // Bad fit
      // Intercept and slopes taken from QNewton output, viz. "Haux"
      // Will be used to determine whether inActive.
      THlx hlx; hlx(0) = Xlow; t.Haux.Extrapolate(hlx,false);
      float x0 = hlx(0);
      float y0 = hlx(1), yp = hlx(3), z0 = hlx(2), zp = hlx(4);
      vector<int>::const_iterator ihit;
      int ipl, joined; for (ipl = iplf_z1, joined = 0; ipl>=iplow; ipl--) {

	// ********** LOOP OVER PLANES in SM1<->SM2 **********

	int idx = lat->dico_idx[ipl]; if (idx<0) continue;
	int jpl = ipl-ipli_z1; if (hitsPats[jpl/32]&(1<<jpl%32)) continue;

	const TPlane &p = setup.vPlane(ipl); if (!p.IFlag) continue;// Plane OFF
	const TDetect &d = setup.vDetect(p.IDetRef);
	if (11<=d.IType && d.IType<=18) continue;// Skip drift (ST04/DR>"iplow")
	if (39<=d.IType && d.IType<=40) continue;// Skip MAs
	float yy = y0+yp*(d.X(0)-x0), zz = z0+zp*(d.X(0)-x0);
	if (!d.InActive(yy,zz)) {
	  if (d.DZSize(0)) continue;
	  else { // Special case of VSAT detectors...
	    // ...These are the most important here (the tracks that stand best
	    // chance to be rescued by "BackTrackZ2" are precisely those at
	    // large momentum and paraxial) and we do not want to miss any
	    // hits there, because of a marginally wrong extrapolation (since
	    // the redundancy is low): let's relax the "InActive" check. We can
	    // do it easily for a VSAT detector, where there is no DZ.
	    double yyp = d.X(1)+(yy-d.X(1))*.9, zzp = d.X(2)+(zz-d.X(2))*.9;
	    if (!d.InActive(yyp,zzp)) continue;
	  }
	}
	bool isScifi = d.IType==22;
	bool isPixel = // Excluding MPM: anyway, as of 2012, none in zone 0x2.
	  p.IFlag&0x30;
	if (iter==0 &&                         // 1st iter: vertical coord's...
	    (fabs(d.Ca)>.1 || d.X(0)<580.)) continue;  // ...downstream of FI05
	float cut = // Facilitate the extrapolation in case of scifi, 2nd iter:
	  // - This is the case when, upon 1st iter, 0x4 track has been extended
	  //  to FI06, and is asked, upon 2nd iter, to extrapolate to FI05X,
	  //  which is far away. Since the time of FI05 hit will be checked, it
	  //  will provide good reliability to the hit,...
	  // - Similarly, extra realiability is obtained in the case of a pixel.
	  (((isScifi || isPixel) && iter) ? 8 : 4.5-iter)*d.Resol;
	int kpl = ipl-ipli_z0;
	float guest = t.vGuests[kpl], umn = guest-1.5*cut, umx = guest+1.5*cut;
	static float guestv; if (isPixel) { // Pixel plane => check v coord
	  // (The pixel det. can be an XY type: ordered Y,X,P => ipl-2, or UV
	  // type: ordered P,U,V => ipl+2)
	  if   (p.IFlag&0x10 /* X/Y type */) guestv = t.vGuests[kpl-2];
	  else                               guestv = t.vGuests[kpl+2];
	}
	THit *hbest; int mult; float best; static TTrack *tbest;
	for (ihit = p.vHitRef().begin(), mult = 0, hbest = 0, best = cut;
	     ihit!=p.vHitRef().end(); ihit++) { // Loop over hits
	  THit &h = vecHit[*ihit]; // ref to Thit vector element
	  TTrack *tp = 0; if (!h.sTrackID().empty()) {
	    // Let's consider unused hits... But also ``badly'' used hits. There
	    // can be many such in zone 0x2, used by very unreliable tracks.
	    if (h.sTrackID().size()==1) {
	      int id = *h.sTrackID().begin(); list<TTrack>::iterator jt;
	      for (jt = listTrack.begin(); jt!=listTrack.end(); jt++) {
		TTrack &tj = *jt; if ((int)tj.Id==id) { tp = &tj; break; }
	      }
	    }
	    if (!tp || tp->Type!=0x2 || tp->IMark!=-2) continue;
	  }
	  if (isScifi &&     // Scifi: not much redundancy in compensation:...
	      fabs(h.Time-evtT)>2) continue;  // ...strict time cut
	  float u = h.U; if (u<umn) continue; if (u>umx) break;
	  if (isPixel && // One uses the "cut" calcuted for pixel's u coord!
	      fabs(h.V-guestv)>cut) continue;
	  mult++; float diff = fabs(u-guest);
	  if (diff<best) { best = diff; hbest = &h; tbest = tp; }
	}
	if (hbest && (mult==1 || best<cut/2)) {
	  if (iter==0) {        // ***** 1ST ITER: LOOK FOR SPACE POINTS *****
	    float Z = hbest->U, cz = d.Ca, sz = d.Sa;
	    const TStation *&s = p.Station;
	    float urange = s->URange, vrange = s->VRange;
	    static THit *ybest, *ubest, *vbest; // To avoid compiler warning
	    float dubest = urange, dvbest = vrange;
	    int nhbest = 0; 

	    vector<int>::const_iterator ipy, ipu0, ipu, ipv0, ipv;
	    vector<int>::const_iterator ihy, ily, ihu, ilu, ihv, ilv;
	    for (ipy = s->IPlanes.begin(); ipy!=s->IPlanes.end(); ipy++) {

	      // ***** LOOP ON PLANES FOR ``Y'' (``Y'' is anything but Z)
	      const TPlane &py = setup.vPlane(*ipy);
	      if (py.IPlane==p.IPlane) continue;
	      if (!py.IFlag) continue;                      // Plane OFF
	      if (py.vHitRef().empty()) continue;           // Plane empty
	      const TDetect &dy = setup.vDetect(py.IDetRef);
	      float cy = dy.Ca, sy = dy.Sa, szy = sz*cy-cz*sy;
	      float Xyz = (dy.X(0)+d.X(0))/2;
	      float Y, U, dU, V, dV;
	      if (fabs(szy)<.1) continue;       // ``Y'' is anything but Z

	      ipu0 = ipy; ipu0++;  // "ipu" iterator for looking for ``U''
	      bool assy = ipu0!=s->IPlanes.end() && py.Associate==&setup.vPlane(*ipu0);
	      // If associate exists and comes next => increment "ipu0"...
	      if (assy) {
		ipu0++;
		// ...and if in addition it's ON => set flag "assy"
		assy &= !py.Associate->IFlag && !
		  py.Associate->vHitRef().empty();
	      }
	      ihy = py.vHitRef().begin(); ily = py.vHitRef().end();

	      bool nexty = true; while (nexty) {

		// ***** LOOP ON ``Y'' HITS *****
		THit &hy = vecHit[*ihy];
		if (!hy.sTrackID().empty()) goto nexty;
		Y = hy.U;

		for (ipu = ipu0; ipu!=s->IPlanes.end(); ipu++) {

		  // ***** LOOP ON PLANES FOR ``U'' *****
		  const TPlane &pu = setup.vPlane(*ipu);
		  if (pu.IPlane==p.IPlane) continue;
		  if (!pu.IFlag) continue;                    // Plane OFF
		  if (pu.vHitRef().empty()) continue;         // Plane empty
		  const TDetect &du = setup.vDetect(pu.IDetRef);
		  float cu = du.Ca, su = du.Sa;
		  float syu = sy*cu-cy*su, szu = sz*cu-cz*su;
		  float dXu = du.X(0)-Xyz;
		  if (fabs(szu)<.1 ||            // ``U'' is anything but Z...
		      fabs(syu)<.1) continue;    // ...and ``Y''
		  float Uyz = (Y*szu-(Z+zp*dXu)*syu)/szy;
		  bool uup = szu/szy>0;

		  ipv0 = ipu; ipv0++;  // "ipv" iterator for looking for ``V''
		  bool assu = ipv0!=s->IPlanes.end() && pu.Associate==&setup.vPlane(*ipv0);
		  // If associate exists and comes next => increment "ipv0"...
		  if (assu) {
		    ipv0++;
		    // ...and if in addition it's ON => set flag "assu"
		    assu &= !pu.Associate->IFlag &&
		      !pu.Associate->vHitRef().empty();
		  }
		  ihu = pu.vHitRef().begin(); ilu = pu.vHitRef().end();

		  bool nextu = true; while (nextu) {

		    // ***** LOOP ON ``U'' HITS *****
		    THit &hu = vecHit[*ihu];
		    if (!hu.sTrackID().empty()) goto nextu;
		    U = hu.U; dU = U-Uyz;

		    if (uup) {
		      if      (dU<-urange) goto nextu;
		      else if (dU> urange) {  // It's over for this TPlane
			ihu = ilu; goto nextpu;   // Try with associate
		      }
		    }
		    else {
		      if      (dU> urange) goto nextu;
		      else if (dU<-urange) {  // It's over for this TPlane
			ihu = ilu; goto nextpu;   // Try with associate
		      }
		    }
		    if (nhbest<=3 && fabs(dU)<dubest) {
		      bool ok = true; if (isScifi) {
			ok &= fabs(hy.Time-evtT)<2;
			ok &= fabs(hu.Time-evtT)<2;
		      }
		      if (ok) {
			nhbest = 3; dubest = fabs(dU);
			ybest = &hy; ubest = &hu;
		      }
		    }

		    for (ipv = ipv0; ipv!=s->IPlanes.end(); ipv++) {

		      // ***** LOOP ON PLANES FOR ``V'' *****
		      const TPlane &pv = setup.vPlane(*ipv);
		      if (pv.IPlane==p.IPlane) continue;
		      if (!pv.IFlag) continue;                // Plane OFF
		      if (pv.vHitRef().empty()) continue;     // Plane empty
		      const TDetect &dv = setup.vDetect(pv.IDetRef);
		      float cv = dv.Ca, sv = dv.Sa;
		      float syv = sy*cv-cy*sv, szv = sz*cv-cz*sv;
		      float dXv = dv.X(0)-Xyz;
		      // no condition: ``V'' comes necessarily next
		      float Vyz = (Y*szv-(Z+zp*dXv)*syv)/szy;
		      bool vup = szv/szy>0;

		      bool assv = pv.Associate && !pv.Associate->IFlag &&
			!pv.Associate->vHitRef().empty();
		      ihv = pv.vHitRef().begin(); ilv = pv.vHitRef().end();

		      bool nextv = true; while (nextv) {

			// ***** LOOP ON ``V'' HITS *****
			THit &hv = vecHit[*ihv];
			if (!hv.sTrackID().empty()) goto nextv;
			V = hv.U; dV = V-Vyz;
			dV -= dU*cu/cv*dXv/dXu;

			if (vup) {
			  if      (dV<-vrange) goto nextv;
			  else if (dV> vrange) {
			    ihv = ilv; goto nextpv;
			  }
			}
			else {
			  if      (dV> vrange) goto nextv;
			  else if (dV<-vrange) {
			    ihv = ilv; goto nextpv;
			  }
			}
			if (nhbest<4 || fabs(dV)<dvbest) {
			  bool ok = true; if (isScifi) {
			    if (fabs(hv.Time-evtT)>2) ok = false;
			  }
			  if (ok) {
			    nhbest = 4;
			    dubest = fabs(dU); dvbest = fabs(dV);
			    ybest = &hy; ubest = &hu; vbest = &hv;
			  }
			}

		      nextv:
			ihv++;
		      nextpv:
			if (ihv==ilv) {
			  if (assv) {       // Continuation on associate plane
			    assv = false;
			    const TPlane *&pva = pv.Associate;
			    ihv = pva->vHitRef().begin();
			    ilv = pva->vHitRef().end();
			    const TDetect &dva = setup.vDetect(pva->IDetRef);
			    cv = dva.Ca; sv = dva.Sa;
			    syv = sy*cv-cy*sv; szv = sz*cv-cz*sv;
			    dXv = dva.X(0)-Xyz;
			    Vyz = (Y*szv-(Z+zp*dXv)*syv)/szy;
			  }
			  else nextv = false;
			}
		      }   // End loop on ``V'' hits
		    }   // End loop on planes for ``V''
		  nextu:
		    ihu++;
		  nextpu:
		    if (ihu==ilu) {
		      if (assu) {       // Continuation on associate plane
			assu = false;
			const TPlane *&pua = pu.Associate;
			ihu = pua->vHitRef().begin();
			ilu = pua->vHitRef().end();
			const TDetect &dua = setup.vDetect(pua->IDetRef);
			cu = dua.Ca; su = dua.Sa;
			syu = sy*cu-cy*su; szu = sz*cu-cz*su;
			dXu = dua.X(0)-Xyz;
			Uyz = (Y*szu-(Z+zp*dXu)*syu)/szy;
		      }
		      else nextu = false;
		    }
		  }   // End loop on ``U'' hits
		  if (nhbest) break;             // Space point found: give up
		}   // End loop on planes for ``U''
	      nexty:
		if (++ihy==ily) {
		  if (assy) {       // Continuation on associate plane
		    assy = false;
		    const TPlane *&pya = py.Associate;
		    ihy = pya->vHitRef().begin(); ily = pya->vHitRef().end();
		    const TDetect &dya = setup.vDetect(pya->IDetRef);
		    cy = dya.Ca; sy = dya.Sa, szy = sz*cy-cz*sy;
		    Xyz = (dya.X(0)+d.X(0));
		  }
		  else nexty = false;
		}
	      }   // End loop on ``Y'' hits
	      if (nhbest) break;               // Space point found: give up
	    }   // End loop on planes for ``Y''
	    if (nhbest) {
	      t.AddHit(*hbest); t.AddHit(*ybest); t.AddHit(*ubest);
	      if (nhbest==4) t.AddHit(*vbest);
	      vector<int>::const_iterator is;
	      for (is = s->IPlanes.begin(); is!=s->IPlanes.end(); is++) {
		int kpl = (*is)-ipli_z1; hitsPats[kpl/32] |= 1<<kpl%32;
	      }
	      joined += nhbest; if (isScifi) nScifis += nhbest;
	      nSpacePts += 1;
	      added += joined;
	      if (tbest) tracksToErase.push_back(tbest);
	    }
	  } // End case iter==0
	  else {              // ***** 2ND ITER: NO SPACE POINT REQUIRED *****
	    t.AddHit(*hbest); joined++; nSpacePts += .334;
	    if (isScifi) {   // Grant 2/3 ``space point'': 2 scifis suffice
	      nScifis++;  nSpacePts += .334;
	    }
	    hitsPats[jpl/32] |= 1<<jpl%32;
	    added += joined;
	    if (tbest) tracksToErase.push_back(tbest);
	  } // End case iter==1
	} // End of Hit found
      } // End loop over planes
      if (joined>=4 ||           // ********** ENOUGH HITS ADDED... **********
	  nScifis &&     // ********** ...or SPECIAL CASE of SCIFIS **********
	  // (Scifi hits are more reliable (because of the tight cut on time).
	  // Therefore one can be less demanding in that case (Hopefully so
	  // since there's almost no redundancy there). What is requisite:
	  //  - A Y proj (which is always the case in the 1st iter) => it 
	  //   ensures that the 0x2 extension is consistent with current &0x4
	  //   track in the Y dimension.
	  //  - 2 other proj's, which are needed to ensure consistency in the
	  //   bending dimension, for, to 1st order, any X proj is compatible
	  //   w/ current track, given the extra DoF brought by momentum.
	  // In all generality, this requisite should only be checked on the
	  // full set of added hits collected during the 2 iterations. This
	  // would allow a FI06YV+FI05XY track. But would imply a rewriting 
	  // of the method => let's keep it as is for the time being.)
	  (iter==0 && nScifis>=3 || iter)) {
	if (nSpacePts<2) {              // Special case only 1 ``space point''
	  // => Determine whether we have >1 hits in the bending dimension.
	  bendingInfo = t.BendingInfo();
	}
	else bendingInfo = 1;
	if (bendingInfo) {
	  if (!t.QNewtonFit(1,7) || t.Chi2tot/(t.NDics-5)>TOpt::dCut[9]) {
	    added = -1; break;
	  }
	}
	else {
	  bool ok = t.QNewtonFit(1,6) && t.Chi2tot/(t.NDics-4)<TOpt::dCut[9];
	  if (ok) {
	    TTrack *beam = beams[iBestB];
	    float y0B = beam->Hfirst(1)-beam->Hfirst(3)*beam->Hfirst(0);
	    float z0B = beam->Hfirst(2)-beam->Hfirst(4)*beam->Hfirst(0);
	    float ypB = beam->Hfirst(3), zpB = beam->Hfirst(4);
	    float xB = (t.Haux(2)-t.Haux(4)*t.Haux(0)-z0B)/(zpB-t.Haux(4));
	    t.Haux(0) = xB; t.Haux(1) = y0B+ypB*xB;
	    ok = t.QNewtonFit(1,6) && t.Chi2tot/(t.NDics-4)<TOpt::dCut[9];
	  }
	  if (!ok) { added = -1; break; }
	}
#ifdef DEBUG
	if (idebug) t.DumpHits(1);
#endif
      }
      else break;                 // ***** ...else BREAK
    } // End of 2 iterations
#warning TODO: Reevaluate LR after extrapolating into target<->SM1
    if (nSpacePts>=2 || nSpacePts>1 && nScifis>=4) {
      t.Hlast(5) = t.Haux(5); t.Hlast(5,5) = 1.E-4;
      if (t.QuickKF(-1,1) && t.Chi2tot/(t.NDFs-5)<TOpt::dCut[9]) {
	t.IFit = 0x2;
	if (t.QuickKF(1,1)) t.IFit |= 0x4;
      }
      if ((t.IFit&0x6)!=0x6) {  // Backtracking...
	t.Behead(iplf_z1);          // ...failed => Discard added THits
	t.QuickKF(-1,0); t.QuickKF(1,0);       // Straight fit
	t.IFit = 0x1;                            // Reset fit type code
	t.Hfirst(5)=t.Hlast(5) = 0;            // Reset momentum
      }
      else {
	t.Type |= 0x2;	            // ...OK     => Update TTrack's Type
	if (!bendingInfo)  // In case not enough bending info: remember that
	  // the result of QN fit, obtained while fixing Y @ X=0 and stored
	  t.IFit |= 0x80;  // in "Haux", can be useful.
	if (nScifis>=added-1) // Scifi-(almost)only, cf. "TEv::PrePattern"
	  t.Scifi |= 0x2;
	for (int i = 0; i<(int)tracksToErase.size(); i++)
	  tracksToErase[i]->IMark = -1;
      }
    }
    else if (added) {
      // ********** IF FAILED, RESET: HITS LIST and CHI2 **********
      t.Behead(iplf_z1); t.Chi2tot = oldChi2;
    }
    it++;
  }
  it = listTrack.begin(); while (it!=listTrack.end()) {
    TTrack &t = *it; if (t.Type==0x2) {
      if (t.IMark==-1) { listTrack.erase(it++); continue; }
      if (t.IMark==-2)// Track was temporarily flagged as unreliable, but...
	t.IMark = 0;  // ...haven't seen any of its hits used => Clear the flag.
    }
    it++;
  }
}

// ===========================================================================
// ==============================  BackTrackSAS ==============================
// ===========================================================================
void TEv::BackTrackSAS()
{
  // ***** EXTEND SAS(0x6) TRACKS, w/ momentum, UPSTREAM OF SM1 (ZONE 0x1) *****

  // - Loop on 0x6 tracks.
  // - In 2 steps, w/ highly reliable (good chi2, high P) tracks first,
  // - Tracks are extrapolated back into zone 0x1, hits being picked up along
  //  the route.
  // - The backtracking + pick-up executed in a loop. At the end of each 
  //  iteration:
  //   - The # of picked ups hits is checked.
  //   - If enough, the track is refitted.
  //   - If chi2 turns out too high, worst hit are cleaned. Then, refit again.
  //   - A new, more precise extrapolation, is derived from the refit.
  // - Dico fit is used:
  //   - To extrapolate back into zone 0x1.
  //   - To evaluate the quality of the extended 0x7 track in the making.
  // - KF fit is used to evaluated the final 0x7 track: "TTrack::QuickKF" or
  //  "FullKF", depending upon option TOpt::ReMode[46].
  // - Special cases:
  //   I) 0x6 track is crossing the RICH pipe (and then badly extrapolates into
  //     zone 0x1).
  //     - As determined from the pipe parameters stored in "TOpt::RICHPipe".
  //     - The special procedure is conditioned by "RICHPipe[5]>0", i.e. the
  //      pipe is made out of heavy stuff.
  //     - The extrapolation is based on that part of the 0x6 track piece that's
  //      upstream of RICH, to which the 0x6 momentum is assigned. => It's not
  //      disturbed by the ms undergone when passing through the pipe.
  //     => Therefore, there's a condition on the # of hits, in fact the # of
  //       space points and # of proj. evaluated by "TTrack::Evaluate", that
  //       this extrapolation be meaningful.

  // Todo: investigate whether to require 0x6 track to start upstream of RICH.

  const TSetup &setup = TSetup::Ref();
  int ipli_z0 = setup.vIplFirst()[0], iplf_z0 = setup.vIplLast()[0];
  int nPats = (iplf_z0-ipli_z0)/32+1; unsigned int hitsPats[nPats];
  double xf_z0 = setup.vDetect(iplf_z0).X(0);
  int ipli_z1 = setup.vIplFirst()[1];
  static int nStrippedGPs = -1; if (nStrippedGPs<0) {
    int ipl; for (ipl = ipli_z0, nStrippedGPs = 0; ipl<ipli_z1; ipl++)
      if (setup.vDetect(ipl).IType==28) nStrippedGPs++;
  }

  // ***** BEAM: TO BE USED in case of TTrack W/O BENDING INFO... *****
  TTrack *beams[2]; int nBeams; GetIncidentTracks(2,nBeams,beams);

  for (int iLoop = 0; iLoop<2; iLoop++) {

    // **********************************************************************
    // ***** 2 LOOPS: 1ST LOOP w/ VERY GOOD TRACKS, i.e. high p, low chi2
    //                2ND LOOP w/ THE REST
    // **********************************************************************

  list<TTrack>::iterator it = listTrack.begin(); while (it!=listTrack.end()) {
    TTrack &t = *it; if ((t.Type&0x7)!=0x6) {
      it++; continue;      // ********** LOOP on ``SAS only'' TTracks **********
    }

    if (fabs(1/t.Hfirst(5))>30 || t.Chi2tot/(t.NDFs-5)<2) {// p>30GeV, chi2/NDF<2
      if (iLoop)  { it++; continue; }
    }
    else {
      if (!iLoop) { it++; continue; }
    }

#ifdef DEBUG
    static int idebug = 0; if (idebug) t.DumpHits(0);
#endif
    double chi2tot0x6 = t.Chi2tot;// Save KF's "Chi2tot" before it's upset by QN
    double oldChi2tot = chi2tot0x6; // Init value of reference chi2
    int iFit0x6 = t.IFit;         // Save also fit type;
    double p0x6 = t.Hlast(5), dp20x6 = t.Hlast(5,5);

    //  ***** RICH pipe: heavy material? is track crossing it? *****
    bool account4RIPipe = Account4RIPipe(t);

    //             ******************** REFIT ********************
    //  "BackTrackSAS" relies on the momentum derived from the SAS reco: in the
    // course of its execution we will need to fit and refit the extended track
    // so as to evaluate the quality of the extension.
    //     => Use, fast, QN fit for this.
    //  Unfortunately QN performs badly on long (spanning zones 0x7) tracks (
    // because does not take into account MS).
    //     => Use special mode =2 of "TTrack:QNewtonFit", where the fit is
    // restricted, if possible, to zone 0x2 while momentum is fixed = SAS's.
    //  Yet, there are cases where the reco in zones 0x6 alone doesn't provide
    // for a reliable value of momentum, refered to as "!bendingInfo" infra. For
    // these: let's not fix momentum and supplement 0x6 info by a contraint on
    // compatibility of the 0x6 track w/ the incident beam candidates.
    //  As long we proceed w/ "QNewtonFit", no need to back up "TTrack::Hfirst"
    // (or "Hlast"), since these are then not modified. "TTrack::Chi2tot" has
    // already been backed up, cf. supra.
    bool bendingInfo; int iBestB = -1;
    if (t.IFit&0x80) bendingInfo = false;  // TTrack W/O BENDING INFO
    // => "guests" are meaningfull => do not refit
    else {
      if (!(bendingInfo = t.BendingInfo())) {
	// Case VSAT-only 0x2 segment w/o FI05X. PR2 produces such . Would
	// possibly be better to discard them before they reach "TEv::BridgeS"
	if ((iBestB = t.Fit2IncidentTrack(nBeams,beams,0))==-1) {
	  it++; continue;
	}
      }
      else {      // Require QN fit, viz. inside dico,etc..
	// Note: TTrack is left undisturbed by the, failed, attempt @ QN
	// fitting. => One need not reset it before moving to the next
	if (account4RIPipe) {
	  // 
	  if (!bendingInfo || !t.QNewtonFit(1,10)) { it++; continue; }
	}
	else if (!t.QNewtonFit(1,2)) { it++; continue; }
      }
    }

    //  ***** EXTRAPOLATE TO UPSTREAM OF SM1 => COVARIANT MATRIX *****
    // As the #tracks candidate for BackTrackSAS is reasonably small, one can
    // spend some CPU on using TTRack's extrapolation method to derive the
    // covariant matrix used later on to open appropriately wide routes.
    THlx hf_z0; hf_z0(0) = xf_z0; t.Hfirst.Extrapolate(hf_z0,true);

    int iter, iterMx = 3, added, projs, nProjs, nScifis, nAllSis;
    int nExpected = 0, nGMPed = 0; // # of expected hits: total or sole GMP hits
    double sCu; double dTMx = 2; // Assuming 1ns scifi resol, taking 2 sigmas
    double tST = t.SigmaTime, dT = sqrt(9*tST*tST+dTMx*dTMx);
    double tMn = t.MeanTime-dT, tMx = tMn+2*dT, tUp, tLow;
    for (iter=added=projs=nProjs=nScifis=nAllSis = 0, sCu = 0,
	   tUp = tMn, tLow = tMx, memset(hitsPats,0,nPats*sizeof(unsigned int));
	 iter<iterMx; iter++) {
      // ***** 3 ITERATIONS w/ NARROWING ROUTES, CONCLUDED by QN FIT *****
      if (t.Chi2tot==1e6) break;    // Bad fit

      // Intercept and slopes taken from QNewton output, viz. "Haux"
      // Will be used to determine whether inActive.
      float x0 = t.Haux(0);
      float y0 = t.Haux(1), z0 = t.Haux(2), yp = t.Haux(3), zp = t.Haux(4);

      // Covariant matrix: a fast extrapolation will be used. 2x2 blocks
      // are tranformed by transfer matrix T via  T*C*tT where T is
      // T(0,0)=T(1,1)=1, T(0,1)=-L, T(1,0) = 0; L being extrapolation length.
      TMatrixF cov[3] = { TMatrixF(2,2), TMatrixF(2,2), TMatrixF(2,2) };
      float xCov = hf_z0(0), L; 
      int m, k, l; for (m = 0, k = 1; m<2; m++, k++) {
	cov[m](0,0) = hf_z0(k,k); l = k+2; cov[m](1,1) = hf_z0(l,l);
	cov[m](0,1)=cov[m](1,0) = hf_z0(k,l);
      }
      cov[2](0,0) = hf_z0(1,2); cov[2](1,1) = hf_z0(3,4);
      cov[2](0,1) = hf_z0(1,4); cov[2](1,0) = hf_z0(3,2);
      float dy2 = cov[0](0,0), dz2 = cov[1](0,0), dyz = cov[2](0,0);

      int prjs = projs; // Set of projections, corresponding to current QN fit.
      int nSis; // Number of Sis in current iteration
      vector<THit*> joinedHits;
      // KF fit may be updated, cf. "FIRST TIME EVER SIs" infra => back it up.
      THlx hBack = t.Hfirst; double cBack = t.Chi2tot;
      int ipl, joined; vector<int>::const_iterator ihit; bool xInfo = sCu>.49;
      for (ipl = iplf_z0, joined=nSis = 0, nExpected=nGMPed = 0; ipl>=ipli_z0;
	   ipl--) {
	// LOOP OVER PLANES in target->SM1
	const TPlane &p = setup.vPlane(ipl); if (p.IFlag==0) continue; // OFF
	const TDetect &d = setup.vDetect(p.IDetRef);
	bool isPPixel = p.IFlag&0x30, isMPixel = p.IFlag&0xc0;
	bool isGMP = isPPixel|isMPixel;
	int jpl = ipl-ipli_z0; if (hitsPats[jpl/32]&(1<<jpl%32)) {
	  nExpected++; if (isPPixel) nExpected++; if (isGMP) nGMPed++;
	  continue;
	}
	if (d.IType==15) // Exclude DCs (No LR, no reassessment of THit::Status)
	  continue;
	if (d.IType==32 && !isMPixel) // => Isolated MPM, cf. "TSetup::Init"...
	  // ...Not covered by present algorithm (because no predicting v-coord)
	  continue; // => Exclude.

	float yy = y0+yp*(d.X(0)-x0), zz = z0+zp*(d.X(0)-x0);
	if (!d.InActive(yy,zz)) continue;      // ***** ACTIVE AREAS ENCOUNTERED
	nExpected++; if (isPPixel) nExpected++;              // ***** COUNT THEM
	if (isGMP) nGMPed++;

	//              ***** EXCLUDE SOME DETECTORS/PLANES depending UPON ITER.
	if (iter==0) {                          // 1st iteration...
	  if (p.IProj!=1 && !isPPixel) continue;// ...coord==Z or pixel detector
	  if (d.IType==28 &&                    // ...bypass stripped GPs if...
	      nStrippedGPs<4) continue; // ...#coord's too small for space point
	}
	bool isScifi = d.IType==22, isSi = d.IType==21;
	if ((L = xCov-d.X(0))>100) {  // Update cov matrix
	  xCov = d.X(0);
	  for (int m = 0; m<3; m++) {
	    TMatrixF c(cov[m]);
	    cov[m](0,0) = c(0,0)-L*(c(0,1)+c(1,0))+L*L*c(1,1);
	    cov[m](0,1) = c(0,1)-L*c(1,1); cov[m](1,0) = c(1,0)-L*c(1,1);
	  }
	  dy2 = cov[0](0,0); dz2 = cov[1](0,0); dyz = cov[2](0,0);
	  //#define DEBUG_EXTRAP_COV
#ifdef DEBUG_EXTRAP_COV
	  // Check fast covariant matrix extrapolation against standard one
	  THlx hlx; hlx(0) = d.X(0); t.Hfirst.Extrapolate(hlx);
	  float dy2p = hlx(1,1), dz2p = hlx(2,2), dyzp = hlx(1,2);
	  printf("\n%f %f %f   %f %f %f\n",dy2,dz2,dyz,dy2p,dz2p,dyzp);
          for (int m = 0; m<3; m++) cov[m].Print();
	  hlx.Print();
#endif
	}
	float cu = d.Ca, su = d.Sa;
	float du = sqrt(dy2*cu*cu+2*dyz*cu*su+dz2*su*su);
	float cut = 5*du+3*d.Resol; if (!xInfo) cut += 4*cu*du;

	int proj = 1<<p.IProj;
	float extra = (proj&prjs) ? .5*d.Ca+.25*fabs(d.Sa) : d.Ca+.5*fabs(d.Sa);
	//float cut2 = nSigmas*d.Resol+nExtras*extra*(d.Resol<.01 ? .01 : d.Resol);
	float guest = t.vGuests[jpl], umn = guest-1.5*cut, umx = guest+1.5*cut;
	static float guestv; int vprojPix; // GMP: consistency along v-axis
	if (isPPixel) {
	  int iplv = (p.IFlag&0x10 /* XY type */) ? ipl-2 : /* UV */ ipl+2;
	  guestv = t.vGuests[iplv-ipli_z0];
	  vprojPix = 1<<setup.vPlane(iplv).IProj; // Tag w/ proj. of v coord
	}
	else if (isMPixel) {
	  // We take the predicted value of the v coord from the neighbouring
	  // MP detector. This is approximation for the Z abscissa, but given
	  // the poor resolution in v of the MPs, it's good enough an approx.
	  int iplv = (p.IFlag&0x40 /* Y|U type */) ? ipl-2 : /* X|V */ ipl+2;
	  guestv = t.vGuests[iplv-ipli_z0];
	  vprojPix = 0;
	}
	else vprojPix = 0;
	THit *hbest; int mult; float best; vector<TSpacePt> spts;
	for (ihit = p.vHitRef().begin(), mult = 0, hbest = 0, best = cut;
	     ihit!=p.vHitRef().end(); ihit++) {  // Loop over hits
	  THit &h = vecHit[*ihit]; // ref to Thit vector element
	  float u = h.U; if (u<umn) continue; if (u>umx) break;
	  if (isScifi) {  // Scifis: check time
	    // (Note: A good scifi candidate must agree both in time and space.)
	    double hT = h.Time; if (hT<tMn || tMx<hT) continue;
	  }
	  if (isPPixel) { // One uses "cut", which was set for pixel's u coord!
	    if (fabs(h.V-guestv)>cut) continue;
	  }
	  else if (isMPixel) { // One neglects the track contribution to the
	    // uncertainty and apply a factor sqrt(12)
	    if (fabs(h.V-guestv)>3.5*h.SigV) continue;
	  }
	  mult++;
	  if ((iter || isPPixel) && (!isSi || nAllSis)) {
	    //           ***** PICK UP INDIVIDUAL HITS *****
	    //  We resort to this simpler solution after the picking up hits has
	    // already been initiated by a space point and hence, we can rely on
	    // a solid backtracking extrapolation into zone 0x1. Pixel planes (
	    // i.e. pixelGP, as of 2012/06) are an exception: don't have enough
	    // coordinates for "BuildSpacePt" to work, in any case. And they
	    // yield automatically a space point because of their very pixel
	    // nature of their hits.
	    //  All hits are considered, even if not free: so as not to bias
	    // the determination of the best candidate.
	    float diff = fabs(u-guest);
	    if (diff<best) { best = diff; hbest = &h; }
	  }
	  else if (h.sTrackID().empty()) {
	    //      ***** BUILD SPACE POINT BASED on VERTICAL HIT *****
	    //  Sace points are expected to be more reliable as a starting
	    // point. Particularly for SIs: since we are to take advantage of
	    // their full resolution, whereas it's artificially worsened
	    // elsewhere in the PR, in order to make for possible misalignment (
	    // misalignment due to thermal expansion in the support frames of
	    // the various COMPASS detectors and not expected to show up w/in
	    // the tightly packed setup of a given single SI station).
	    //  (Note: "BuildSpacePt" requires a hit pattern: useless her, since
	    // the current station of detectors is expected to be virgin
	    // territory. (And it also updates it: not clear why...))
	    unsigned int hPats[nPats];
	    memcpy((void*)hPats,(void*)hitsPats,nPats*sizeof(unsigned int));
	    BuildSpacePt(t,&h,hPats,dy2,dyz,dz2,spts);
	  }
	}
	if ((iter || isPPixel) && (!isSi || nAllSis)) {
	  //   ***** INDIVIDUAL HITS CONTINUATION: DETERMINE BEST HIT *****
	  if (hbest && hbest->sTrackID().empty() && (mult==1 || best<cut/2)) {
	    if (isScifi) {  // Scifis: Update timing
	      double hT = hbest->Time;
	      if (hT<tLow) tLow = hT; if (hT>tUp) tUp = hT; nScifis++;
	    }
	    else if (isSi) nSis++;
	    t.AddHit(*hbest); joined++;	sCu += cu; joinedHits.push_back(hbest);
	    hitsPats[jpl/32] |= 1<<jpl%32;
	    if (!(proj&projs)) nProjs++; projs |= proj;
	    if (vprojPix) {                    // Pixel plane (so far only GPP)
	      joined++;                        // ...increment # of joined DoFs
	      if (!(vprojPix&projs)) nProjs++; // ...update set of proj.'s
	      projs |= vprojPix;               // ...add proj. of v coord
	    }
	  }
	}
	else if (spts.size()!=0) {
	  //  ***** SPACE POINT CONTINUATION: DETERMINE BEST SPACE POINT *****
	  int ispt, nspts = spts.size(), nHMx; double chi2Mn; static int ispt0;
	  for (ispt=nHMx = 0, chi2Mn = TOpt::dCut[69]; ispt<nspts; ispt++) {
	    // Note: "chi2Mn" is set = max. possible chi2 in "BuildSpacePt"
	    // => "ispt0" is always defined at the exit of this loop.
	    TSpacePt &spt = spts[ispt];
	    if (spt.nHits>nHMx || spt.nHits==nHMx && spt.chi2<chi2Mn) {
	      nHMx = spt.nHits; chi2Mn = spt.chi2; ispt0 = ispt;
	    }
	  }
	  const TSpacePt &spt0 = spts[ispt0];
	  if (isSi) {// ***** SPECIAL CASE of SI STATION: Re-EVALUATE CHI2 *****
	    // Uncertainties on hits have been so far artificially worsened (see
	    // comment supra). Re-evaluate the chi2 w/ true unertainties.
	    double y = spt0.y, z = spt0.z, du; int i;
	    for (i = 0, chi2Mn = 0; i<8; i++) if (spt0.hPat&1<<i) {
	      THit *h = spt0.hs[i]; int kpl = h->IPlane;
	      const TDetect &dk = setup.vDetect(kpl);
	      du = (h->U-y*dk.Ca-z*dk.Sa)/(h->SigU-TOpt::dCut[85]);
	      chi2Mn += du*du;
	    }
	    chi2Mn = sqrt(chi2Mn/spt0.nHits);
	  }
	  if (spt0.nHits>=3 &&    // ***** CUT on CHI2 => ADD SPACE POINT *****
	      (chi2Mn<8 && mult==1 || chi2Mn<4)) {
	    for (int i = 0; i<8; i++) if (spt0.hPat&1<<i) {
	      THit *h = spt0.hs[i]; int kpl = h->IPlane, lpl = kpl-ipli_z0;
	      if (1<<lpl%32&hitsPats[lpl/32]) {
		CsErrLog::msg(elError,__FILE__,__LINE__,
 "Adding a 2nd THit on TTrack %d @ TPlane %d",t.Id,h->IPlane); continue;
	      }
	      t.AddHit(*h); joined++; joinedHits.push_back(h);
	      hitsPats[lpl/32] |= 1<<lpl%32;
	      if      (isScifi) nScifis++;
	      else if (isSi) nSis++;
	      const TPlane &pk = setup.vPlane(kpl); int proj = 1<<pk.IProj;
	      if (!(proj&projs)) nProjs++; projs |= proj;
	      const TDetect &dk = setup.vDetect(kpl); sCu += dk.Ca;
	    }
	    if (nSis && nAllSis==0) // If it's first time SIs are associated...
	      break;   // ...break before trying to extend the track any further
	  } // Valid space point exists
	} // End Best hit or best Space point
      } // End loop over planes
      if (nSis && nSis<2) { // The region of Sis facing the interacion vertex
	// can be crowded and it's hence easy to pick up fakes...
	//             ***** => DO NOT ALLOW SINGLE Si *****
	for (vector<THit*>::iterator ih = joinedHits.begin();
	     ih!=joinedHits.end(); ih++) {
	  int ipl = (*ih)->IPlane; if (setup.vDetect(ipl).IType==21) {
	    t.SubHit(ipl); joinedHits.erase(ih); joined--; break;
	  }
	}
      }
      //     ********** EVALUATE HITS ADDED IN ZONE 0x1: REFIT **********
      if (joined) {
	if (!bendingInfo) {  // ***** TTrack INITIALLY W/O BENDING INFO *****
	  if (!(bendingInfo = t.BendingInfo())) {     // IF STILL TRUE...
	    if (iter) { added = -1; break; }            // ...IF 2ND ITER: EXIT
	    int nBs; TTrack ** bs; if (iBestB>=0) { 
	      nBs = 1;      bs = &beams[iBestB];
	    }
	    else {
	      nBs = nBeams; bs = beams;
	    }
	    if (t.Fit2IncidentTrack(nBs,bs,1)==-1) {
	      added = -1; break;                        // ...ELSE: NEED BEAM
	      // (Can happen despite "IFit&080": "t" has changed)
	    }
	  }
	  else { t.IFit &= 0xff7f; }
	}
	if (bendingInfo) {
	  int QNmode; // In order to evaluate the quality of the extension into
	  // LAS just obtained supra, we need, QN, refit.
	  //  In case SIs are involved, the combination of a long track (even if
	  // restricted to zones 0x3) and high resolution SI, puts QN to the
	  // limit. => Let's then use special mode =9 of "TTrack::QNewtonFit".
	  // (Note that this implies freeing the momentum, we will later on need
	  // validate the  extension by a full, KF, fit.
	  if      (account4RIPipe) QNmode = 10;
	  else if (nAllSis+nSis>2) QNmode =  9;
	  else                     QNmode =  2; 
	  if (t.QNewtonFit(1,QNmode)) {   // ********** REFIT **********
	    //  Check resulting chi2 against "upper chi2/NDF cut conditioning
	    // QNFit for BackTracking" ("dCUt[89]").
	    //  Because the BackTracking hit search is little demanding and
	    // therefore not very reliable, all the reliability comes from the
	    // goodness of chi2. => Let's be strict. And also, let's try to
	    // clean the track extension so that chi2 be even lower than "upper
	    // cut". This cleaning should in principle be brought into play by
	    // checking chi2 against a lower chi2 cut. In order to keep tne # of
	    // cuts low, let's recycle "dCUt[89]", times some factor <1.
	    double chi2 = t.Chi2tot/(t.NDics-4), chi2Ref = TOpt::dCut[89];
	    double cleaningThr = // For very good track, upon 1st iteration,...
	      // ...when ambiguity can be raised in next iter, be more demanding
	      // Cf. W40/22001-52613 evt #1050430
	      chi2Ref*(iLoop==0 && iter==0 ? .75 : .9);
	    bool unDo = chi2>20;                // ***** REQUIRE REASONABLE CHI2
	    if (!unDo && chi2>cleaningThr) {
	      //            ***** MEDIUM CHI2: CLEANING *****
	      // E.g. W40/22001-52613 evt #1049860. (Note: best would be to
	      // subtract hits w/ large chi2 increment.)
	      list<int>::iterator ih; float worsts[2]; THit *hWorsts[2];
	      for (ih = t.lHitPat.begin(), hWorsts[0]=hWorsts[1] = 0,
		     worsts[0]=worsts[1] = 4.5; ih!=t.lHitPat.end(); ih++) {
		if (*ih<0) continue; THit &h = vecHit[*ih];
		int ipl = h.IPlane; if (ipl>iplf_z0) break;
		const TPlane  &p = setup.vPlane(ipl); int proj = 1<<p.IProj;
		const TDetect &d = setup.vDetect(p.IDetRef);
		float residu = fabs(h.U-t.vGuests[ipl-ipli_z0])/h.SigU;
		if (p.IFlag&0xf0) { // PixelGP/MP
		  int iplv = (p.IFlag&0x50 /* X/Y type */) ? ipl-2 : ipl+2;
		  float res2 = fabs(h.V-t.vGuests[iplv-ipli_z0])/h.SigV;
		  if (res2>residu) residu = res2;
		}
		if (residu>worsts[0]) {
		  worsts[1] = worsts[0]; hWorsts[1] = hWorsts[0];
		  worsts[0] = residu;    hWorsts[0] = &h;
		}
		else if (residu>worsts[1]) {
		  worsts[1] = residu;    hWorsts[1] = &h;
		}
	      }
	      int iw, nDFsBack = t.NDFs; float chi2s[2]; 
	      int iwStrt = worsts[0]>worsts[1]*2 ? 0 : 1;
	      for (iw = iwStrt, chi2s[0]=chi2s[1] = 0; iw>=0; iw--) {
		THit *h = hWorsts[iw]; if (!h) continue;
		t.SubHit(*h);
		if (t.QNewtonFit(1,QNmode)) chi2s[iw] = t.Chi2tot/(t.NDics-4);
		else                        chi2s[iw] = -1;
		if (iw) t.AddHit(*h);
	      }
	      int erased = 0;
	      if (chi2s[1]>0 && (chi2s[1]<chi2s[0] || chi2s[0]<0)) {
		t.SubHit(*hWorsts[1]); t.AddHit(*hWorsts[0]);
		chi2s[1]=chi2s[0] = -1; hWorsts[0] = hWorsts[1];
		if (t.QNewtonFit(1,QNmode)) chi2s[0] = t.Chi2tot/(t.NDics-4);
		else                        chi2s[0] = -1;
		erased = 1;
	      }
	      vector<THit*>::iterator jh = joinedHits.begin();
	      while (jh!=joinedHits.end()) {
		if (*jh==hWorsts[erased]) jh = joinedHits.erase(jh);
		else jh++;
	      }
	      if (chi2s[0]<0 || chi2Ref<chi2s[0]) unDo = true;
	      else if (0<chi2s[0] && (chi2s[0]<chi2*.9 || chi2>chi2Ref)) {
		// Validate cleaning if nice gain in chi2 or original chi2 is
		// high (higher than exact upper chi2/NDF cut).
		list<int>::iterator ip;
		for (ih = t.lHitPat.begin(), ip = t.lPlnRef.begin(),
		       nProjs=projs = 0; ih!=t.lHitPat.end(); ih++, ip++) {
		  if (*ip>iplf_z0) break; if (*ih<0) continue;
		  const TPlane  &p = setup.vPlane(*ip); int proj = 1<<p.IProj;
		  if (!(proj&projs)) { nProjs++; projs |= proj; }
		  if (p.IFlag&0x30) { // Pixel plane (P-pixel, that is)
		    if (!(0x2&projs)) { nProjs++; projs |= 0x2; }
		  }
		}
		joined -= nDFsBack-t.NDFs;
		int ipl = hWorsts[0]->IPlane, jpl = ipl-ipli_z0;
		hitsPats[jpl/32] &= ~(1<<jpl%32);
		const TDetect &d = setup.iPlane2Detect(ipl); sCu -= d.Ca;
		if      (d.IType==22) nScifis--;// Note: no update of "tUp/tLow"
		else if (d.IType==21) nSis--;
		if (joined==0) break;
	      }
	    }
	    if (!unDo && nSis && nAllSis==0) {// ***** FIRST TIME EVER SIs *****
	      // When SIs, QN fit is not constrained by fixed momentum taken
	      // from SAS reco. => Have to check the LAS vs. SAS compatibility.
	      //  =>                                                ***** KF FIT
	      int fastFit; if (TOpt::ReMode[46]) {
		t.InsertMissedPlanes(); fastFit = 0; unDo = !t.FullKF(-1);
	      }
	      else {                    fastFit = 1; unDo = !t.QuickKF(-1,1); }
	      if (!unDo) {
		double newChi2 = t.Chi2tot/(t.NDFs-5);
		t.IFit = 0;// Remember that track KF-refitted after modification
		unDo = newChi2>chi2Ref*1.5;                // ***** CUT BAD CHI2
	      }
	      if (fastFit && !unDo && (t.Chi2tot-oldChi2tot)/joined>chi2Ref) {
		//       ***** BAD CHI2 INCREMENT of EXTENSION: CHECK w/ FULL KF
		oldChi2tot = t.Chi2tot;    // => New reference chi2 
		unDo = !t.FullKF(-1) || t.Chi2tot/(t.NDFs-5)>chi2Ref*0.75;
	      }
	      else oldChi2tot = t.Chi2tot; // => New reference chi2 
	      if (!unDo) {
		if (iter!=iterMx)          // Update cov matrix for nex iter
		  t.Hfirst.Extrapolate(hf_z0,true);
		int ipli_t = t.lPlnRef.front(), ipl;
		THlx H1 = t.Hfirst, H2, &Hi = H1, &Hj = H2;
		for (ipl = ipli_t-1; ipl>=ipli_z0; ipl--) {
		  int jpl = ipl-ipli_z0; if (hitsPats[jpl/32]&(1<<jpl%32))
					   continue;
		  const TDetect &d = setup.vDetect(ipl);
		  Hj(0) = d.X(0); Hi.Extrapolate(Hj,false);
		  t.vGuests[ipl-ipli_z0] = Hj(1)*d.Ca+Hj(2)*d.Sa;
		  THlx &tmp = Hi; Hi = Hj; Hj = tmp;
		}
		H1 = t.Hfirst; Hi = H1; Hj = H2;
		for (ipl = ipli_t+1; ipl<=iplf_z0; ipl++) {
		  int jpl = ipl-ipli_z0; if (hitsPats[jpl/32]&(1<<jpl%32))
					   continue;
		  const TDetect &d = setup.vDetect(ipl); if (d.IType!=21) break;
		  Hj(0) = d.X(0); Hi.Extrapolate(Hj,false);
		  t.vGuests[ipl-ipli_z0] = Hj(1)*d.Ca+Hj(2)*d.Sa;
		  THlx &tmp = Hi; Hi = Hj; Hj = tmp;
		}
	      }
	    }
	    if (unDo) {// ***** BAD CHI2: UNDO LATEST EXTENSION and BREAK  *****
	      if (joined) {
		for (int i = 0; i<(int)joinedHits.size(); i++) {
		  THit *h = joinedHits[i]; int ipl = h->IPlane; t.SubHit(ipl);
		  const TDetect &d = setup.vDetect(ipl);
		  if (d.IType==22) nScifis--;
		}
		t.Clip(); t.Hfirst = hBack; t.Chi2tot = cBack;
	      }
	      break;
	    }
	    nAllSis += nSis;  // Update total #SIs
	  }
	  else { added = -1; break; }
	}
	added += joined; // Update total # of added hits after latest iter
#ifdef DEBUG
	if (idebug) t.DumpHits(1);
#endif
      }
      else break;
    }  // End of 3 iterations
    //       ********** REQUIRE minimum # of HITS and PROJS... **********
    // Several considerations to take into account:
    // - Basic requirement: less demanding than std PR, because I) it need be
    //  so: we want to rescue those tracks that fail to pass PR), II) it can be
    //  so: the compatibility between picked up hits and SAS track adds to the
    //  reliability.
    int nHitsMn = // PR is typically 8 => 7. If so, we allow here for one
      // failing MM in a secondary from a V0 decay downstream of 1st MM.
      TOpt::iPRpar[5]-1;
    // - Next requirement: Provide for reco'ing VSAT tracks, where little
    //  redundancy: 6FIs (2002:7) | 8SIs+1(+1)GP (2008) | 8SIs+2(1+1)G/MP. We
    //  based ourself on "nExpected" and "nGMPed", while pure scifi case is
    // dealt w/ infra.
    if (nExpected<=10) { // Case of 2SIs+GP: allowing, e.g, for one failing SI
      // and GP, each removing 2 hits due to correlation, or downstream decay.
      // And let's be more demanding in case of a RICH-pipe crossing.
      nHitsMn = account4RIPipe ? 7 : 6;
    }
    if (nExpected<=12 &&   // Case of 2SIs+GP+MP: distinguish it from 12*MMs...
	nGMPed>=3)         // ...by requiring in addition some associated SIs.
      nHitsMn = 7;
    if ((added>=nHitsMn ||
	 ((nScifis>=4 && added>=4 ||             // ...Special case of scifis
	   nScifis>=3 && added>=5) &&             // E.g mc.2003.01 Evt #95
	  tUp-tLow<2*dTMx)) && 
	nProjs>=3) {
      if (t.Hlast(5,5)==0) t.Hlast(5,5) = 1.E-4; t.IFit = 0;
      int fastFit; bool ret; if (TOpt::ReMode[46]) {
	fastFit = 0;                               // FullKF systematic
	t.InsertMissedPlanes(); ret = t.FullKF(-1);
      }
      else {
	fastFit = 1;                               // QuickKF
	ret = t.QuickKF(-1,1);
      }
      if (ret &&         // ********** REFIT and CHECK CHI2 **********
	  t.Chi2tot/(t.NDFs-5)<((t.Scifi&0xff) ? TOpt::dCut[16]*1.5 : TOpt::dCut[17])) {
	if (fastFit) ret = t.FullKF(1);
	else         ret = t.QuickKF(1,1);
	if (ret) {
	  if (nScifis>=added-1) // Scifi-(almost)only, cf. "TEv::PrePattern"
	    t.Scifi |= 0x1;
	  for (list<int>::iterator ih = t.lHitPat.begin(); ih!=t.lHitPat.end();
	       ih++) {
	    if (*ih<0) continue; THit &h = vecHit[*ih];
	    if (h.IPlane>iplf_z0) break;
	    h.status = 1;                            // ***** FLAG ZONE 0x1 HITS
	  }
	  t.IFit |= 0x6; added = 0; t.Type |= 0x1;
	}
      }
    }
    if (added) {// ***** MODIF FAILED => RESET HITS LIST, HELICES, CHI2 *****
      t.Behead(iplf_z0);
      if (t.IFit==0) { // TTrack has been KF-refitted while in modified state...
	t.Hlast(5) = p0x6; t.Hlast(5,5) = dp20x6;   // ***** ...=> RESET HELICES
	bool ret; if (TOpt::ReMode[46]) {
	  ret = t.FullKF(-1);    if (ret) ret = t.FullKF(1);
	}
	else {
	  ret = t.QuickKF(-1,1); if (ret) ret = t.QuickKF(1,1);
	}
	if (ret) t.IFit = 0x6;
	else t.IMark = -1;
      }
      else {
	t.Chi2tot = chi2tot0x6; t.IFit = iFit0x6;  // ***** ReINSTATE CHI2, IFIT
      }
    }
    it++;
  } // End loop on tracks
  } // End of 2 loops: 1st w/ very good tracks, 2nd w/ the rest
  list<TTrack>::iterator it = listTrack.begin(); while (it!=listTrack.end()) {
    if ((*it).IMark==-1) { listTrack.erase(it++); continue; }
    it++;
  }
}

// ===========================================================================
// ==============================  BackTrackZ1  ==============================
// ===========================================================================
void TEv::BackTrackZ1(int selection)
{
  // ********** EXTEND SM1<->RICH TRACKS (momentum-less) **********
  // ********** ASSUMING VERTEX @ X WHERE Ztrack==Zbeam **********

  // special drift

  const TSetup &setup = TSetup::Ref();
  int ipli = setup.vIplFirst()[0], iplf = setup.vIplLast()[0]; 
  // ***** NO MORE THAN 64 TPlane's in SPAN! *****
  if (iplf-ipli>63) return;
  int nPats = (iplf-ipli)/32+1;
  unsigned int plPats[nPats], hitsPats[nPats];

  //      ******************** PRELIMINARY STEPS ********************

  // - Determine #hits per TStation => Exclude TStation w/ too many hits.
  // - Reassess "THit::Status": set it =0 for those THit's free to be used, i.e.
  //  not added to TTrack themselves and not mirror of an added THit.
  //   "THit::Status" has been set in "TEv::PrePattern2" and may have been upset
  //  in "TEv::CleanTrackList".
  //   (Secondarily, some of its settings have been assigned by "ImportClusters".
  //  But these concern characteristics that are not at stake here: (as of
  //  05/08): coalescence ("Status==-4" must have prevented THit from being
  //  added to tracks and must have remained unchanged throughtout PR, and we
  //  keep it as is) and LR probability (the info is used only in PR, may have
  //  been overwritten (in PR, when THit is added) and we feel free to overwrite
  //  it again now).
  //  ("THits::Status" will have to be taken care of throughout this method,
  //  for at the difference to "THit::sTrackID" it is not (and cannot easily) be
  //  updated by "TTrack::AddHit(and related methods)". For simplicity sake,
  //  we will take care of it even for non-drift THit's)
  int ipl = ipli, kpl, nSs = 0, nPlsS, nHitsS;
  for (int iPat = 0; iPat<nPats; iPat++) plPats[iPat] = 0xffffffff;
  while (ipl<=iplf) {
    const TPlane &p = setup.vPlane(ipl); const TStation *&s = p.Station;
    for (kpl = 0, nPlsS=nHitsS = 0; kpl<(int)s->IPlanes.size(); kpl++) {
      ipl = s->IPlanes[kpl]; const TPlane &pk = setup.vPlane(ipl);
      int jpl = ipl-ipli;
      if (!pk.IFlag) { plPats[jpl/32] &= ~(1<<jpl%32); continue; }// Plane's OFF
      vector<int>::const_iterator ihit; int nHitsP;
      for (ihit = pk.vHitRef().begin(), nHitsP = 0; ihit!=pk.vHitRef().end();
	   ihit++) {
	THit &h = vecHit[*ihit]; int status = h.sTrackID().empty() ? 0 : 1;
	if (status==0 && s->Type==15) { // Special case of DC...
	  int mirror; // ...reassess THit's "Status" together with its mirror's
	  if (h.Status==-4)  // Exclude silent component of coalesced cluster
	    status = -4;
	  else if ((mirror = h.Mirror)!=-1 &&
		   !vecHit[mirror].sTrackID().empty()) status = 1;
	}
	if      (status==0) { h.status = 0; nHitsP++; }
	else if (status==1) 
	  h.status = 1;  // Meaning also GEMs spacer hits info erased
      }
      if (nHitsP>(s->Type==15 ? 24 /* Larger #hits for DCs?? */ :24))
	plPats[jpl/32] &= ~(1<<jpl%32);// ***** DISCARD PLANE W/ TOO LARGE #HITS
      else {
	nHitsS += nHitsP; nPlsS++;
      }
    }
    if (nHitsS>(s->Type==15 ? 20 /* Larger #hits for DCs */ :16)*nPlsS) {
      //                 ***** DISCARD ALL PLANES OF TStation W/ TOO LARGE #HITS
      for (kpl = 0; kpl<(int)s->IPlanes.size(); kpl++) {
	int jpl = s->IPlanes[kpl]-ipli; plPats[jpl/32] &= ~(1<<jpl%32);
      }
    }
    else nSs++;
    ipl++;
  }
  const TLattice *lat = TLattice::Ptr();     
  if (!lat || nSs<2) return;

  //            *************** REQUIRE  BEAM ***************
  // (No more than 2 beams. Most appropriate (cf. infra) will be used to fix
  // TTrack's Y @ ``vertex'', and QN fit w/ this constraint in order to
  // extrapolate over SM1 despite lack of momentum and predict positions in
  // zone 0x1.)
  TTrack *beams[2]; int nBeams; GetIncidentTracks(2,nBeams,beams);
  if (nBeams== 0 || nBeams>2) return;

  static bool first = true; static TH2D *hChi2 = 0; static int iplp_RICH = -1;
  if (first) {
    first = false;
    if (isMC && TOpt::Hist[7]) {                          // ***** HISTOS *****
      CsHistograms::SetCurrentPath("/Traffic/MCmonitor");
      hChi2  = new TH2D("h0x1Chi2","0x1 BackTracking perf vs. RW-#chi^{2}",
			10,0,TOpt::dCut[69],25,-1/48.,1+1/48.);
      CsHistograms::SetCurrentPath("/");
    }
    double xRICH = 750;
    for (int ipl = setup.vIplFirst()[1]; ipl<setup.vIplLast()[1]; ipl++) {
      iplp_RICH = ipl; if (setup.vDetect(ipl).X(0)>xRICH) break;
    }
    if (iplp_RICH<0) CsErrLog::mes(elFatal,
      "Inconsistency: No detector plane downstream of RICH in zone 0x2");
  }

  // ***** INITIALISATION: In view of UPDATING DRIFTS THits w/ EVENT TIME ***** 
  //        ***** REQUIRE THAT TRIGGER MATCHES ReMode[36]... *****
  //                                                           ***** GET TRIGGER
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
  // When it comes to pick up a scifi hit, we apply a time cut on the time diff
  // of the hit w.r.t. the time of the event.
  // When drift CsCluster's are updated (by "UpdateDrifts"), "eventTime" is set
  // = 0, while the actual value is backed-up into "eventTRef".
  double evtT = eventTRef ? eventTRef : eventTime;

  list<TTrack>::iterator it = listTrack.begin(); while (it!=listTrack.end()) {
    TTrack &t = *it;
    if ((t.Type&0x7)!=0x2 || t.Hfirst(0)>900)  { it++; continue; }
    if (!(t.Scifi&0x202) &&
	(t.Chi2tot/(t.NDFs-4)>TOpt::dCut[17] ||
	 t.Hfirst(0)<500 &&  t.Hlast(0)<750 &&
	 // Stricter chi2 cut for track in SM1-RICH zone and fully outside field
	 t.Chi2tot/(t.NDFs-4)>TOpt::dCut[88])) { it++; continue; }
    if (selection && !(t.IFit&0x10))           { it++; continue; }

    // ******************** LOOP OVER SM1<->SM2 tracks... ********************
    // ******************** ...UPSTREAM of RICH           ********************

    // Upstream of RICH, 'cause otherwise, extrapolation will be uncertain. The
    // cut x0<xRICH will be supplemented by the QN fit with code 8.
    int added = 0; int nScifis = 0, nSIs = 0;
    //#define BTZ1_DEBUG 2
#if defined DEBUG || defined BTZ1_DEBUG
#  ifdef BTZ1_DEBUG
    static int idebug = 1;
#  else
    static int idebug = 0;
#  endif
    if (idebug) t.DumpHits(0);
#endif
    double oldChi2 = t.Chi2tot; // Save KF's "Chi2tot" before it's upset by QN

    // ********** FIT TO BEST INCIDENT TRACK **********
    // (Chi2, from the fit w/ fix Y @ vertex, is not selective => Use "xB"
    // in its stead => Cannot use "TTrack::Fit2IncidentTrack".)
    int iB, iBestB, nBs = nBeams==1 ? 1 : 3; float cutB = 250, bestB;
    double fldint = setup.MagFieldInt[setup.NMags-2];
    for (iB = 0, iBestB = -1, bestB = cutB; iB<nBs; iB++) {
      // ***** EXTRAPOLATE BEAM to Z beam == Z track *****
      // (do this crudely (i.e. straight) since x=0 already crude assumption)
      TTrack *beam;
      if      (iB<2)      beam = beams[iB];
      else if (iBestB>=0) beam = beams[iBestB];
      else break;
      float y0B = beam->Hfirst(1)-beam->Hfirst(3)*beam->Hfirst(0);
      float z0B = beam->Hfirst(2)-beam->Hfirst(4)*beam->Hfirst(0);
      float ypB = beam->Hfirst(3), zpB = beam->Hfirst(4);
      float xB = (t.Hfirst(2)-t.Hfirst(4)*t.Hfirst(0)-z0B)/(zpB-t.Hfirst(4));
      t.Haux(5) = 	// Guess value from horizontal angle
	(0-t.Hfirst(3))/(1.E-3 * 0.3 * fldint);
      // Let's fix the vertex Y @ X=xB, Since we do not have yet a reliable
      // estimate of the Y angle @ that abscissa (and hence cannot determine,
      // w/in "QNewtonFit", the corresponding value of the Y track parameter
      // @ the origin of the dico => do it in 2 steps.
      t.Haux(0) = 0; t.Haux(1) = y0B; t.Haux(2) = z0B;
      t.Haux(3) = (t.Hfirst(3)-y0B)/t.Hfirst(0); t.Haux(4)= t.Hfirst(4); 
      if (t.QNewtonFit(1,8)) { // First fit sets the Y' track parameter
	xB = (t.Haux(2)-t.Haux(4)*t.Haux(0)-z0B)/(zpB-t.Haux(4));
	t.Haux(0) = 0; t.Haux(1) = y0B+ypB*xB;
	bool ret = t.QNewtonFit(1,8); double chi2 = t.Chi2tot/(t.NDics-4);
	if (ret &&
	    // Check resulting chi2 against "upper chi2/NDF cut conditioning
	    // QNFit for BackTracking" ("dCUt[89]"). The cut is not exactly
	    // dedicated to the case, since it's also used in similar (upper cut
	    // on QNFit, for BackTracking) yet distinct (different QN mode)
	    // contexts.
	    (chi2<TOpt::dCut[89] ||
	     // Allow some more margin for low momentum, hence high MS, tracks
	     chi2<TOpt::dCut[89]*1.2 && fabs(t.Haux(5))>.5)) {
	  if (iB<2) {
	    if (iBestB==-1 || fabs(xB)<bestB) {
	      iBestB = iB; bestB = fabs(xB);
	    }
	  }
	}
	else if (iB==2) iBestB = -1;
      }
      else if (iB==2) iBestB = -1;
    }
    if (iBestB==-1) { it++; t.Chi2tot = oldChi2; continue; }

    // Segment upstream of RICH: used infra to determine the uncertainty of the
    // extrapolation into zone 0x1.
    TTrack tp = t; tp.Shorten(iplp_RICH); tp.QuickKF(-1,0);

    int iter; for (iter = 0, memset(hitsPats,0,nPats*sizeof(unsigned int));
		   iter<2; iter++) {
      //            ********************************************
      //            *************** 2 ITERATIONS ***************
      //            ********************************************

      if (t.Chi2tot==1e6) break;    // Bad fit

      // Intercept and slopes taken from QNewton output, viz. "Haux"
      // Will be used to determine whether inActive.
      float x0 = t.Haux(0);
      float y0 = t.Haux(1), z0 = t.Haux(2), yp = t.Haux(3), zp = t.Haux(4);

      // Covariant matrix: a fast extrapolation will be used. 2x2 blocks
      // are tranformed by transfer matrix T via  T*C*tT where T is
      // T(0,0)=T(1,1)=1, T(0,1)=-L, T(1,0) = 0; L being extrapolation length.
      TMatrixF cov[3] = { TMatrixF(2,2), TMatrixF(2,2), TMatrixF(2,2) };
      float xCov = t.Hfirst(0), L; 
      int m, k, l; for (m = 0, k = 1; m<2; m++, k++) {
	cov[m](0,0) = tp.Hfirst(k,k); l = k+2; cov[m](1,1) = tp.Hfirst(l,l);
	cov[m](0,1)=cov[m](1,0) = tp.Hfirst(k,l);
      }
      cov[2](0,0) = tp.Hfirst(1,2); cov[2](1,1) = tp.Hfirst(3,4);
      cov[2](0,1) = tp.Hfirst(1,4); cov[2](1,0) = tp.Hfirst(3,2);
      float dy2 = cov[0](0,0), dz2 = cov[1](0,0), dyz = cov[2](0,0);

      unsigned int hPats[nPats]; int ipl, joined;
      for (ipl = iplf, joined = 0, 
	     memcpy((void*)hPats,(void*)hitsPats,nPats*sizeof(unsigned int));
	   ipl>=ipli; ipl--) {

	// ********** LOOP OVER PLANES in SM1->target **********

	int jpl = ipl-ipli; unsigned int bpl = 1<<jpl%32;
	if (!(plPats[jpl/32]&bpl)) continue; if (hitsPats[jpl/32]&bpl) continue;
	const TPlane &p = setup.vPlane(ipl);
	const TDetect &d = setup.vDetect(p.IDetRef);
	float cu = d.Ca, su = d.Sa;
	if (iter==0) { if (cu>.05) continue; } // 1st iter: coord==Z
	else { if (d.IType==15) continue; }   // 2nd iter => No LR => No DC
	float yy = y0+yp*(d.X(0)-x0), zz = z0+zp*(d.X(0)-x0);

	if (!d.InActive(yy,zz)) continue;          // ***** IN ACTIVE AREA *****
	bool isScifi = d.IType==22, isSI = d.IType==21;
	float du = 0;// ***** UNCERTAINTY of GUESTIMATE for CURRENT PLANE  *****
	if (iter && added>2) { // Case plane belongs to already added TStation?
	  const TStation *&s = p.Station; int kpl, nPlsS;
	  for (kpl=nPlsS = 0; kpl<(int)s->IPlanes.size(); kpl++) {
	    int lpl = s->IPlanes[kpl], mpl = lpl-ipli;
	    if (hitsPats[mpl/32]&1<<mpl%32) {
	      const TPlane &pl = setup.vPlane(ipl);
	      const TDetect &dl = setup.vDetect(pl.IDetRef);
	      float cl = dl.Ca, sl = dl.Sa, dul = dl.Resol/(cl*cu+sl*su);
	      du += 1/dul/dul; nPlsS++;
	    }
	  }
	  if (nPlsS>2 && nPlsS>s->IPlanes.size()*.6) // We restrict ourself...
	    // ...to cases where the internal consistency w/in the TStation is
	    // already ascertained (for we do not want to preclude that one at
	    // least of the already associated hits is wrong) => therefore we
	    // require at least 3 coordinates and at least 5 hits in the case of
	    // an 8 plane DC.
	    du = 1/sqrt(du);
	  else
	    du = 0;
	}
	if (!du) {             // Case plane belongs to still excluded TStation
	  if ((L = xCov-d.X(0))>100) {  // Update cov matrix
	    xCov = d.X(0);
	    for (int m = 0; m<3; m++) {
	      TMatrixF c(cov[m]);
	      cov[m](0,0) = c(0,0)-L*(c(0,1)+c(1,0))+L*L*c(1,1);
	      cov[m](0,1) = c(0,1)-L*c(1,1); cov[m](1,0) = c(1,0)-L*c(1,1);
	    }
	    dy2 = cov[0](0,0); dz2 = cov[1](0,0); dyz = cov[2](0,0);
	    //#define DEBUG_EXTRAP_COV
#ifdef DEBUG_EXTRAP_COV
	    // Check fast covariant matrix extrapolation against standard one
	    THlx hlx; hlx(0) = d.X(0); t.Hfirst.Extrapolate(hlx);
	    float dy2p = hlx(1,1), dz2p = hlx(2,2), dyzp = hlx(1,2);
	    printf("\n%f %f %f   %f %f %f\n",dy2,dz2,dyz,dy2p,dz2p,dyzp);
	    for (int m = 0; m<3; m++) cov[m].Print();
	    hlx.Print();
#endif
	  }
	  du = sqrt(dy2*cu*cu+2*dyz*cu*su+dz2*su*su);
	}
	//	float cut = (4.5-iter)*d.Resol, best, diff; 
	float cut = 4*du+3*d.Resol, best, diff; THit *hbest; int mult;
	float guest = t.vGuests[jpl], umn = guest-cut, umx = guest+cut;
	vector<TSpacePt> spts; static int hZStatus;
	vector<int>::const_iterator ihit;
	for (ihit = p.vHitRef().begin(), best = cut, hbest = 0, mult = 0;
	     ihit!=p.vHitRef().end(); ihit++) {
	  THit &h = vecHit[*ihit];
	  if (iter==0) {  // Special case of 1st iter (in vertical, Z, proj)...
	    // ...in order to provide for the reco of e^+e^- pairs, wich don't
	    // separate much if at all in the Z dimension, let's allow space
	    // points to be built on an already used (presumbaly by the other
	    // member of the e-pair) hit. The status of the hit will be set =0
	    // before the call to "BuildSpacePt", and then reset =1. And the
	    // hit will be removed from the space points. Don't do this if the
	    // station has only three coordinates: case of FIs.
	    if (d.IType==22 && h.Status ||
		h.Status==-4) // Exclude silent component of coalesced cluster
	      continue;
	  }
	  else if (h.Status) continue;

	  //                                   ***** LOOP OVER UNUSED HITS *****

	  float u = h.U; if (u<umn) continue; if (u>umx) break; 
	  if (isScifi &&  // Scifi: not much redundancy in compensation:...
	      fabs(h.Time-evtT)>2) continue;   // ...strict time cut
	  if (iter &&         // ***** 2ND ITER: ACCEPT INDIVIDUAL HITS... *****
	      (!isSI ||    // ...Special case of SIs: individual hits only if...
	       // ...there is enough already added hits in zone 0x1 for the
	       // extrapolation to be reliable This would typically require 2
	       // TStations, or, given that 1st iter can only add space points,
	       // 2 such. Yet, one has to provide for high momentum tracks,
	       // which would typically encounter few stations.
	       //  For the time being (09/05) and since no much time to dedicate
	       // to the problem (be it only for testing), I try to modify as
	       // little as possible existing algotithm. The aim being to reco
	       // #7342720 in 2008/W33, evtDump.R.vs.Q-69853.raw, which picks up
	       // ghost SIs. => Let's require that there be already some SIs
	       nSIs)) {
	    diff = fabs(u-guest); if (diff<best) { hbest = &h; best = diff; }
	  }
	  else {     	           // ***** 1ST ITER: REQUIRE SPACE POINTS *****
	    // (w/ incidence derived from "Haux")
	    mult++;
	    hZStatus = h.Status; h.status = 0;
	    BuildSpacePt(t,&h,hPats,dy2,dyz,dz2,spts);
	    if (hZStatus) {// If hit on which space point is built a used one...
	      // ...update #hits and hPat and set = 0 the relevant hit pointer.
	      // (Unfortunatly, now way to change the chi2 at this point.)
	      int nspts = spts.size(); for (int ispt = 0; ispt<nspts; ispt++) {
		TSpacePt &spt = spts[ispt];
		for (int i = 0; i<8; i++) {
		  if (spt.hs[i]!=&h) continue;
		  spt.hs[i] = 0; spt.nHits--; spt.hPat ^= 1<<i; break;
		}
	      }
	      h.status = hZStatus;
	    }
#ifdef BTZ1_DEBUG
#  if BTZ1_DEBUG > 1
	    if (idebug) {
	      printf("BTZ1: ipl %d => %d:",ipl,spts.size());
	      for (int i = 0; i<(int)spts.size(); i++)
		printf(" %d,%d,0x%x",i,spts[i].nHits,spts[i].hPat);
	      printf("\n");
	    }
#  endif
#endif
	  }
	}
	if (iter) {
	  if (hbest) {
	    t.AddHit(*hbest); joined++; hitsPats[jpl/32] |= 1<<jpl%32;
	    if (isScifi) nScifis++; if (isSI) nSIs++;
	  }
	}
	else if (spts.size()!=0) {
	  //                              ***** DETERMINE BEST SPACE POINT *****
	  int ispt, nspts = spts.size(), nHMx; double chi2Mn; static int ispt0;
	  for (ispt=nHMx = 0, chi2Mn = TOpt::dCut[69]; ispt<nspts; ispt++) {
	    // Note: "chi2Mn" is set = max. possible chi2 in "BuildSpacePt"
	    // => "ispt0" is always defined at the exit of this loop.
	    TSpacePt &spt = spts[ispt];
	    if (spt.nHits>nHMx || spt.nHits==nHMx && spt.chi2<chi2Mn) {
	      nHMx = spt.nHits; chi2Mn = spt.chi2; ispt0 = ispt;
	    }
	  }
	  const TSpacePt &spt0 = spts[ispt0];
	  if (spt0.nHits>=3 &&    // ***** GOOD ENOUGH BEST SPACE POINT... *****
	      (chi2Mn<8 && mult==1 || chi2Mn<4)) {
	    //                                 ***** ...=> ADD SPACE POINT *****
	    bool fillHist = hChi2!=NULL; int nGoodHits = 0;
	    if (fillHist) { t.FindKine(); fillHist &= t.IKine>=0; }
	    for (int i = 0; i<8; i++) if (spt0.hPat&1<<i) {
	      THit *h = spt0.hs[i]; int kpl = h->IPlane-ipli;
	      if (1<<kpl%32&hitsPats[kpl/32]) {
		CsErrLog::msg(elError,__FILE__,__LINE__,
			      "Adding a 2nd THit on TTrack %d @ TPlane %d",
			      t.Id,h->IPlane); continue;
	      }
	      if (fillHist && (h->IKine==t.IKine || h->IKin2==t.IKine))
		nGoodHits++;
	      if (TOpt::ReMode[29]&0x1) {// ***** PROPAGATION CORRECTION...*****
		// ...+X-ray+latency required at the PR level.
		const TPlane &p = setup.vPlane(h->IPlane);
		const TDetect &d = setup.vDetect(p.IDetRef);
		CsDetector *csDet = d.PtrDet(); if (csDet->hasDrift()) {
		  // Space point. N.B.: we need only a crude estimate => non
		  // updated values "yp","zp" are good enough for that purpose.
		  double y = y0+yp*(d.X(0)-x0), z = z0+zp*(d.X(0)-x0);
		  bool error; double U = csDet->getCorrU(h->PtrClus(),y*10,z*10,evtT0,error);
		  if (!error) {
		    h->u = U/10; h->Rotate();
		    if (TOpt::ReMode[29]&0x4)  // UPON OPTION: RESET CsCluster
		      // (This is only useful when one intends to make use of
		      // these clusters outside of TraFDic, e.g.:
		      // - TraFDic monitoring when TraFFiC modularity != 2
		      // - Evaluating residuals in coral, as opposed to TraFDiC)
		      h->PtrClus()->setU(U);
		  }
		}
	      }
	      t.AddHit(*h); joined++; hitsPats[kpl/32] |= 1<<kpl%32;
	      if (isScifi) nScifis++; if (isSI) nSIs++;
	    }
	    if (fillHist) hChi2->Fill(spt0.chi2,(double)nGoodHits/spt0.nHits);
	  } // Valid space point exists
	} // Space point exists
      } // End loop over planes
      added += joined;
      if (joined>=3 ||       // ***** ENOUGH HITS ADDED in target<->SM1... *****
	  iter==0 && nScifis || // ...or A SCIFI SPACE POINT @ iter ==0
	  nScifis>4) {          // ...or TOTAL # of SCIFI HITS (they're in time)
	if (!t.QNewtonFit(1,7) ||              // ...=> REFIT
	    t.Chi2tot/(t.NDics-5)>TOpt::dCut[9]) {
	  added = -1; break;
	}
#if defined DEBUG || defined BTZ1_DEBUG
	if (idebug) t.DumpHits(1);
#endif
      }
      else break;           // ***** ...else BREAK *****
    } // End of 2 iterations
    if (added>=7 || added>=5 && nScifis>=5) {
      if (t.QNewtonFit(1,7)) {  // ********** REFIT **********
	//  Check resulting chi2 against "upper chi2/NDF cut conditioning QNFit
	// for BackTracking" ("dCUt[89]").
	//  Because the BackTracking hit search is little demanding and
	// therefore not very reliable, all the reliability comes from the
	// goodness of chi2. => Let's be strict. Yet, since BackTrackZ1 is
	// expected to deal w/ low momenta, let's not be as strict as in
	// "BackTrackSAS".
	double chi2 = t.Chi2tot/(t.NDics-5), chi2Ref = TOpt::dCut[89];
	if (chi2>20) added = -1;       // ***** BAD CHI2: EXIT *****
	else if (chi2>chi2Ref) {    // ***** MEDIUM CHI2: CLEANING *****
	  // E.g. 07W27/32001-58926 evt #1048947.
	  list<int>::iterator ih; float worsts[2]; THit *hWorsts[2];
	  for (ih = t.lHitPat.begin(), hWorsts[0]=hWorsts[1] = 0,
		 worsts[0]=worsts[1] = 4.5; ih!=t.lHitPat.end(); ih++) {
	    if (*ih<0) continue; THit &h = vecHit[*ih];
	    int ipl = h.IPlane; if (ipl>iplf) break;
	    const TPlane  &p = setup.vPlane(ipl); int proj = 1<<p.IProj;
	    const TDetect &d = setup.vDetect(p.IDetRef);
	    float residu = fabs(h.U-t.vGuests[ipl-ipli])/d.Resol;
	    if (p.IFlag&0xf0) { // PixelGP/MP
	      int iplv = (p.IFlag&0x50 /* X/Y type */) ? ipl-2 : ipl+2;
	      float res2 = fabs(h.V-t.vGuests[iplv-ipli])/h.SigV;
	      if (res2>residu) residu = res2;
	    }
	    if (residu>worsts[0]) {
	      worsts[1] = worsts[0]; hWorsts[1] = hWorsts[0];
	      worsts[0] = residu;    hWorsts[0] = &h;
	    }
	    else if (residu>worsts[1]) {
	      worsts[1] = residu;    hWorsts[1] = &h;
	    }
	  }
	  int iw, nDFsBack = t.NDFs; float chi2s[2]; 
	  int iwStrt = worsts[0]>worsts[1]*2 ? 0 : 1;
	  for (iw = iwStrt, chi2s[0]=chi2s[1] = 0; iw>=0; iw--) {
	    THit *h = hWorsts[iw]; if (!h) continue;
	    t.SubHit(*h);
	    if (t.QNewtonFit(1,7)) chi2s[iw] = t.Chi2tot/(t.NDics-5);
	    else                   chi2s[iw] = -1;
	    if (iw) t.AddHit(*h);
	  }
	  if (chi2s[1]>0 && (chi2s[1]<chi2s[0] || chi2s[0]<0)) {
	    t.SubHit(*hWorsts[1]); t.AddHit(*hWorsts[0]);
	    chi2s[1]=chi2s[0] = -1; hWorsts[0] = hWorsts[1];
	    if (t.QNewtonFit(1,7)) chi2s[0] = t.Chi2tot/(t.NDics-5);
	    else                   chi2s[0] = -1;
	  }
	  if (chi2s[0]<0 || chi2Ref*1.5<chi2s[0]) added = -1;
	  else if (chi2s[0]>0) {
	    added--;
	    int ipl = hWorsts[0]->IPlane, jpl = ipl-ipli;
	    hitsPats[jpl/32] &= ~(1<<jpl%32);
	    const TDetect &d = setup.iPlane2Detect(ipl);
	    if (d.IType==22) nScifis--;
	    // What about THit.status
	  }
	}
      }
      else added = -1;
    }
    if (added>=7 || added>=5 && nScifis>=5) {
      //                                    ***** Quick KF FIT: CHECK CHI2 *****
      t.Hlast(5) = t.Haux(5); t.Hlast(5,5) = 1.E-4; t.IFit = 0;
      if (t.QuickKF(-1,1) &&
	  t.Chi2tot/(t.NDFs-5)<TOpt::dCut[9]) {
	t.IFit |= 0x2;
	if (t.QuickKF(1,1)) t.IFit |= 0x4;
      }
      if ((t.IFit&0x6)==0x6) { // ***** CHECK ``TOPOLOGICAL'' REQUIREMENTS *****
	int nSpacePts = 0; TStation *sPrv = 0;
	int projs = 0, nProjs = 0, allProjs = 0, nAllProjs = 0;
	for (list<int>::iterator ih = t.lHitPat.begin(); ih!=t.lHitPat.end();
	     ih++) {
	  const THit &h = vecHit[*ih];
	  const TPlane &p = setup.vPlane(h.IPlane); if (h.IPlane>iplf) break;
	  const TStation *&s = p.Station; if (s!=sPrv) {
	    sPrv = const_cast<TStation*&>(s); /* nStations++; */
	    if (nProjs>=2) { nSpacePts++; nProjs=projs = 0; }
	  }
	  int jproj = 1<<p.IProj;
	  if ((jproj&projs)==0)    { projs |= jproj;    nProjs++; }
	  if ((jproj&allProjs)==0) { allProjs |= jproj; nAllProjs++; }
	}
	if (nProjs>=2) nSpacePts++;
	if (nAllProjs<3) t.IFit = 0;      // # of proj.: strict requirements
	else if (nSpacePts<2) t.IFit = 0; // # of space points: at least 2...   
	else if (nSpacePts<3) {             
	  if (nScifis<5) {    // ...if not scifi: is it due to acceptance?
	    //  The criterion adopted here is to require a "SpacePt" for
	    // each TStation on TTrack's path, w/ TStation being declared on
	    // path if 2 at least of its TPlane's are on path.
	    //  However we provide for decays occuring downstream of the 1st
	    // TStation, allowing 2 TStation's to fire out of 3 if these 2 are
	    // the last 2 and then only if good efficiency: 7 hits
	    float x0 = t.Hfirst(0);
	    float y0 = t.Hfirst(1), z0 = t.Hfirst(2);
	    float yp = t.Hfirst(3), zp = t.Hfirst(4);
	    nSs = 0; int nPls = 0;
	    unsigned int last2Pats[nPats];
	    memset(last2Pats,0,nPats*sizeof(unsigned int));
	    while (ipl<=iplf) {
	      const TPlane &p = setup.vPlane(ipl);
	      const TStation *&s = p.Station;
	      for (kpl=nPlsS = 0; kpl<(int)s->IPlanes.size(); kpl++) {
		ipl = s->IPlanes[kpl]; const TPlane &pk = setup.vPlane(ipl);
		int jpl = ipl-ipli; if (!(plPats[jpl/32]&1<<jpl%32)) continue;
		const TDetect &d = setup.vDetect(p.IDetRef);
		float yy = y0+yp*(d.X(0)-x0), zz = z0+zp*(d.X(0)-x0);
		if (d.InActive(yy,zz)) {
		  nPls++; nPlsS++; if (nSs) last2Pats[jpl/32] |= 1<<jpl%32;
		}
	      }
	      if (nPlsS>=3) nSs++; ipl++;
	    }
	    if (nSs>2) {
	      int i; bool last2only; for (i = 0, last2only = true; i<nPats; i++)
		last2only &= (last2Pats[i]&hitsPats[i])==hitsPats[i];
	      if (!last2only || added<.9*nPls) t.IFit = 0;
	    }
	  }
	}
      }
      if ((t.IFit&0x6)!=0x6) {    // Backtracking...
	t.Behead(iplf);              // ...failed => Discard any added THits
	t.Hfirst(5)=t.Hlast(5) = 0;           // Reset momentum
	t.QuickKF(-1,0); t.QuickKF(1,0);      // Straight fit
	t.IFit = 0x1;                         // Reset fit type code
      }
      else {	                     // ...OK     => Update Type
	t.Type |= 0x1; if (nScifis>=added-1) t.Scifi |= 0x1;
	list<int>::const_iterator ih, ip;
	for (ih = t.lHPat().begin(), ip = t.lPRef().begin();
	     ih!=t.lHPat().end(); ih++, ip++) {
	  if (*ip>iplf) break;
	  if (*ih>=0) {
	    THit &h = vecHit[*ih]; h.status = 1;
	    int mirror; if ((mirror = h.Mirror)!=-1) vecHit[mirror].status = 1;
	  }
	}
      }
    }
    else if (added) {
      t.Behead(iplf);                // ...failed => discard any added THits
      t.Chi2tot = oldChi2;                       //  and reset chi2
    }
    else t.Chi2tot = oldChi2;
    it++;
  }
}
// ===========================================================================
// ==============================  BackTrackVD  ==============================
// ===========================================================================
void TEv::BackTrackVD(){

  // "BackTrackVD": Hit-to-track association in the vertex detector (VD).
  //  It is activated by the Drell-Yan option "TraF ReMode[49]&=0x4", see
  // "TEv::TrackFit2".
  //  Track candidates in the Drell-Yan acceptance ("TOpt::DY_InAcceptance==1"),
  // or all tracks ("TOpt::DY_InAcceptance==0"), are extrapolated backward to
  // the first plane (most downstream) with hits of the VD. For each
  // extrapolated track in the active area of the VD, all hits fulfilling
  // |Time(pivot_hit) - Time(track)| < "TOpt::DY_VD_time" are considered as
  // candidates for the hit-to-track association. By default the time cut is set
  // to 5 ns.
  //  The track parameters are updated for each of the considered hits and,
  // thereafter, the track tries to pick-up the closest hit to the track in the
  // remaining VD projections. Here the search window is limited to 1 sigma (
  // sigma = detector resolution + error of the extrapolation).
  //  The combination of hits is considered to be a good one if and only if the
  // updated track intersects at least one beam upstream of the VD. The
  // intersection is verified in the y coordinate. Moreover, the updated track
  // is extrapolated to the z coordinate of the intersection. Hits are kept if
  // the x coordinate of the beam is compatible with the extrapolated x(track) 
  // within 3 sigmas of the track parameter.
  //  Finally, an extra time cut is applied to the "pivot_hit". The time
  // difference between this hit and the beam must not be greater than 3 ns. All
  // the surviving combinations of VD hits are compared. 
  //  Priority is given to the 3-hits scenario (we consider only combinations of
  // 2 or 3 VD hits). Those combinations originating a chi2/ndf increment
  // greater than ("TOpt::DY_VD_Chi2" *  #hits) are discarded. By default this
  // option is set to a value of 0.3 .The choosen combination of VD hits is the
  // one producing the lowest chi2/ndf. 

  const TSetup &setup = TSetup::Ref();
  static int ipli, iplf; int ipl;
  
  for (ipl = setup.vIplLast()[4]+1, ipli=iplf = -1; ipl<=setup.vIplFirst()[0]-1; ipl++) {
    
    const TPlane &p = setup.vPlane(ipl); if (!p.IFlag) continue;  // Plane OFF
    const TDetect &d = setup.vDetect(p.IDetRef);
    if (d.IType != 22) continue;        // Retain the vertex detector

    if (ipli==-1) 
      ipli = ipl; 
    
    iplf = ipl;
  }


  int nBeams; 
  TTrack *beams[100]; GetIncidentTracks(100,nBeams,beams); // Get all the beams in the event

  if(nBeams == 0) 
    return;

  double y, z, u, dy2, dz2, dyz, du, umx, umn;
  double sZ, cZ, xZ, wZ, Range, deltaU, chi2;
  int nMA02, nMB, nPB, hit_cand, flag, flag1, flag2;

  double ref;
  bool ok = 0;
  bool badfit = 0;

  string strx1 = "FI35X1";
  string strx2 = "FI35X2";
  string strv1 = "FI35V1";
  string strv2 = "FI35V2";
  string stru1 = "FI35U1";
  string stru2 = "FI35U2";

  list<int>::iterator ihp, ipr;
  list<TTrack>::iterator it = listTrack.begin(); while (it!=listTrack.end()) {

    // ********** LOOP on TRACKS **********
    TTrack &t = *it;

    badfit = 0;

    if (!(t.Type&0x3)) { it++; continue; }  // Retains tracks which are bridged over the SM1
    if (!(t.IFit&0x60)) { it++; continue; } // Retains fitted tracks
    if(fabs(1./t.Hfirst(5)) < 2.5){it++; continue;}

    ref =  t.Chi2tot/(t.NDFs-5);

    if(TOpt::DY_InAcceptance){ //Retains tracks in the DY acceptance

      for (ihp = t.lHitPat.begin(), ipr = t.lPlnRef.begin(), nMA02=0, nMB=0, nPB=0; 
	   ihp!=t.lHitPat.end(); ihp++, ipr++) {
	
	if(*ihp<0) continue;
	const TPlane  &p = setup.vPlane(*ipr);
	const TDetect &d = setup.vDetect(p.IDetRef);
	
	if(d.IType == 40)
	  nMA02++;
	if(d.IType == 18 && d.X(0) > 4000)
	  nMB++;
	if(d.IType == 1 && d.X(0) > 4000)
	  nPB++;
      }
      
      if(!(nMA02 > 4 || nPB > 4 || nMB > 6)){it++; continue;}
    }

    // *************** TRACK CANDIDATE for EXTENSION ***************

    vector<int>::const_iterator ih; 

    double BestChi2_3 = ref+(3.0*TOpt::DY_VD_Chi2);  // Maximum increment in chi2 after adding 3 hits
    double BestChi2_2 = ref+(2.0*TOpt::DY_VD_Chi2);  // Maximum increment in chi2 after adding 2 hits

    THit *hBest_0 = 0;
    THit *hBest_1 = 0;
    THit *hBest_2 = 0;

#define NHITS_MX 150
    THit *hit_array[NHITS_MX];
    hit_cand = 0; flag = 0;

    int NHits1 = 0;
    int NHits2 = 0;

    for (ipl = iplf; ipl>=iplf-3; ipl--) { // Loop over planes in the vertex detector

      if(badfit) break;

      if(ipl == iplf-2){
	hit_cand = 0; 
	flag = 0; 
      }

      const TPlane &pZ = setup.vPlane(ipl);
      const TDetect &dZ = setup.vDetect(pZ.IDetRef);

      if(dZ.Name.find(strx1) == 0 || dZ.Name.find(strv1) == 0) 
	NHits1 = pZ.vHitRef().size();

      else if(dZ.Name.find(strx2) == 0 || dZ.Name.find(strv2) == 0) 
	NHits2 = pZ.vHitRef().size();

      if((NHits1 + NHits2) >= NHITS_MX) continue;

      THlx hlx;
      hlx(0) = dZ.X(0); 
      t.Hfirst.Extrapolate(hlx, true);
      y = hlx(1); z = hlx(2);
      dy2 = hlx(1,1); dz2 = hlx(2,2); 
      
      if(!(dZ.InActive(y,z))){
	// Check active area within 1 sigma of y and z
	if (!(dZ.InActive(y+sqrt(dy2), z+sqrt(dz2))) && !(dZ.InActive(y+sqrt(dy2), z-sqrt(dz2))) && !(dZ.InActive(y-sqrt(dy2), z+sqrt(dz2))) && !(dZ.InActive(y-sqrt(dy2), z-sqrt(dz2)))) {
	  if(dZ.Name.find(strx1) == 0 || dZ.Name.find(strv1) == 0 || ((dZ.Name.find(strx2) == 0 && flag == 0) ||  (dZ.Name.find(strv2) == 0 && flag == 0)))
	    continue;
	}
      }
          
      for (ih = pZ.vHitRef().begin(); ih!=pZ.vHitRef().end(); ih++) { // Loop over hits in the downstream plane
	
	THit *h = &vecHit[*ih]; 
	if(fabs(h->Time - t.MeanTime) > TOpt::DY_VD_time) continue; // Removes hits with |Time(pivot_hit) - Time(track)|  > 5 ns (by default)
	
	hit_array[hit_cand] = h; // All hits surviving the time cut are kept as candidates 
	hit_cand++;
	flag = ipl;
	
      }// End loop over hits in the most downstream plane
      
      if(dZ.Name.find(strx1) == 0 || dZ.Name.find(strv1) == 0) continue; // Search for the remaining hits in the X or V projection (there are 2 planes for each projection)
      if(flag == 0) continue; // No pivot_hit was found
    
      THit *hit_array1[NHITS_MX], *hit_array2[NHITS_MX];

      NHits1 = 0;
      NHits2 = 0;
      	
      for(int i = 0; i < hit_cand; i++){ // Loop over all pivot_hit candidates
	
	t.AddHit(*hit_array[i]);
	ok = t.FullKF(-1);
	
	if (!ok) {
	 
	  t.Behead(iplf);
	  ok = t.FullKF(-1);
	  
	  if(!ok){
	    listTrack.erase(it++); 
	    badfit = 1;
	    break;
	  }
	  else
	    continue;
	}
      
	flag1 = 0; flag2 = 0;

	for (int ipl1 = ipl-1; ipl1>=ipli; ipl1--) { // Loop over the remaining planes in the vertex detector
	  
	  const TPlane &pZ = setup.vPlane(ipl1);
	  const TDetect &dZ = setup.vDetect(pZ.IDetRef);
	  
	  if(dZ.Name.find(stru1) == 0){
	    NHits1 = 0;
	    NHits2 = 0;
	  }
	  
	  if(dZ.Name.find(strv1) == 0 || dZ.Name.find(stru1) == 0)
	    NHits1 = pZ.vHitRef().size();

	  else if(dZ.Name.find(strv2) == 0 || dZ.Name.find(stru2) == 0)
	    NHits2 = pZ.vHitRef().size();

	  if((NHits1 + NHits2) >= NHITS_MX) continue;

	  sZ = dZ.Sa, cZ = dZ.Ca, xZ = dZ.X(0), wZ = dZ.Resol;
	
	  THlx hlx1;
	  hlx1(0) = xZ; 
	  t.Hfirst.Extrapolate(hlx1, true);
	  y = hlx1(1); z = hlx1(2);
	  dy2 = hlx1(1,1); dz2 = hlx1(2,2); dyz = hlx1(1,2);
	  
	  u = y*cZ+z*sZ; du = sqrt(dy2*cZ*cZ+2*dyz*cZ*sZ+dz2*sZ*sZ);
	  umx = u+1*du+1*wZ; umn = u-1*du-1*wZ; 
	  
	  if(!(dZ.InActive(y,z))) continue;
	  
	  if(ipl1 == ipl-1 || ipl1 == ipl-3)
	    Range = umx-umn; // Maximum window to search for the closest hit in the U and V projections
	  
	  for (ih = pZ.vHitRef().begin(); ih!=pZ.vHitRef().end(); ih++) { // Loop over hits
	    
	    THit *h = &vecHit[*ih]; 
	    if(fabs(h->Time - hit_array[i]->Time) > 3) continue; // Removes hits with |Time(hit) - Time(pivot_hit)| > 3 ns

	    if(h->U<umn) continue; 
	    if(umx<h->U) break;

	    deltaU = fabs(h->U-u); 
	    
	    if(deltaU < Range && ipl1 > ipli+1){ // Selects the closest hit to the track in the second plane
	      
	      Range = deltaU;
	      hit_array1[i] = h;
	      flag1 = 1;
	    }
	    else if(deltaU < Range && ipl1 <= ipli+1){ // Selects the closest hit to the track in the most upstream plane
	      Range = deltaU;
	      hit_array2[i] = h;
	      flag2 = 1;
	    }
	  }// End loop over hits
	}// End 2nd loop over planes
      
	if(flag1 && flag2){ // Updates the track parameters in case of 3 hits in the vertex detector
	  
	  t.AddHit(*hit_array1[i]);
	  t.AddHit(*hit_array2[i]);
	  
	  ok = t.FullKF(-1);
	
	  if (!ok) {
	    t.Behead(iplf);
	    ok = t.FullKF(-1);
	  	    
	    if(!ok){
	      listTrack.erase(it++); 
	      badfit = 1;
	      break;
	    }
	    else
	      continue;
	  }
        
	  double y0B, z0B, ypB, zpB, xT, yT, zT;
	  int yintersection = 0;
	  
	  for (int it = 0; it < nBeams; it++){
	  
	    if((fabs(0 - beams[it]->MeanTime) - fabs(0 - hit_array2[i]->Time)) > 0 || (fabs(0 - beams[it]->MeanTime) - fabs(0 - hit_array2[i]->Time)) < -3) continue; // The time of the most upstream hit must not differ from the time of the beam by more than 3 ns
	    
	    y0B = beams[it]->Hfirst(1)-beams[it]->Hfirst(3)*beams[it]->Hfirst(0);
	    z0B = beams[it]->Hfirst(2)-beams[it]->Hfirst(4)*beams[it]->Hfirst(0);
	    ypB = beams[it]->Hfirst(3);
	    zpB = beams[it]->Hfirst(4);
	    
	    xT = (t.Hfirst(2)-t.Hfirst(4)*t.Hfirst(0)-z0B)/(zpB-t.Hfirst(4)); // z coordinate (in PHAST notation) of the intersection (using the y coordinate) between the beam and the muon track
	    yT = y0B+ypB*xT; // x coordinate (in PHAST notation) of the intersection
	    zT = z0B+zpB*xT; // y coordinate (in PHAST notation) of the intersection
	    
	    if(xT > -305 && xT < -100 && sqrt((yT*yT)+(zT*zT)) < 2.0){ // Allowed domain for the intersection 
	      THlx helix;
	      helix(0) = xT; 
	      t.Hfirst.Extrapolate(helix, true);
	      
	      if((yT < (helix(1) + 3*sqrt(helix(1,1)))) && (yT > (helix(1) - 3*sqrt(helix(1,1))))){ // Require a 3sigma compatibility between the x coordinate of the beam and track
		yintersection = 1;
		break;
	      }
	    }
	  }
		  
	  if(yintersection != 0){

	    chi2 = t.Chi2tot/(t.NDFs-5);
	    
	    if(chi2 < BestChi2_3){ // Retains the 3-hit combination with the lowest chi2/ndf
	      BestChi2_3 = chi2;
	      hBest_0 = hit_array[i];
	      hBest_1 = hit_array1[i];
	      hBest_2 = hit_array2[i];
	    }
	  }
	}
	else if((flag1 && !flag2) || (!flag1 && flag2)){ // Updates the track parameters in case of 2 hits in the vertex detector
	  
	  if(flag1)
	    t.AddHit(*hit_array1[i]);
	  else 
	    t.AddHit(*hit_array2[i]);

	  ok = t.FullKF(-1);
	
	  if (!ok) {
	    t.Behead(iplf);
	    ok = t.FullKF(-1);
	    
	    if(!ok){
	      listTrack.erase(it++); 
	      badfit = 1;
	      break;
	    }
	    else
	      continue;
	  }
      
	  double y0B, z0B, ypB, zpB, xT, yT, zT;
	  int yintersection = 0;
	  
	  for (int it = 0; it < nBeams; it++){

	    if((fabs(0 - beams[it]->MeanTime) - fabs(0 - hit_array[i]->Time)) > 0 || (fabs(0 - beams[it]->MeanTime) - fabs(0 - hit_array[i]->Time)) < -3) continue;  // The time of the most upstream hit must not differ from the time of the beam by more than 3 ns

	    y0B = beams[it]->Hfirst(1)-beams[it]->Hfirst(3)*beams[it]->Hfirst(0);
	    z0B = beams[it]->Hfirst(2)-beams[it]->Hfirst(4)*beams[it]->Hfirst(0);
	    ypB = beams[it]->Hfirst(3);
	    zpB = beams[it]->Hfirst(4);
	    
	    xT = (t.Hfirst(2)-t.Hfirst(4)*t.Hfirst(0)-z0B)/(zpB-t.Hfirst(4)); // z coordinate (in PHAST notation) of the intersection (using the y coordinate) between the beam and the muon track
	    yT = y0B+ypB*xT; // x coordinate (in PHAST notation) of the intersection
	    zT = z0B+zpB*xT; // y coordinate (in PHAST notation) of the intersection
	    
	    if(xT > -305 && xT < -100 && sqrt((yT*yT)+(zT*zT)) < 2.0){ // Target domain for the intersection 

	      THlx helix;
	      helix(0) = xT; 
	      t.Hfirst.Extrapolate(helix, true);
	      
	      if((yT < (helix(1) + 3*sqrt(helix(1,1)))) && (yT > (helix(1) - 3*sqrt(helix(1,1))))) { // Require a 3sigma compatibility between the x coordinate of the beam and track
		yintersection = 1;
		break;
	      }
	    }
	  }	
	
	  if(yintersection != 0){
	  
	    chi2 = t.Chi2tot/(t.NDFs-5);

	    if(chi2 < BestChi2_2){ // Retains the 2-hit combination with the lowest chi2/ndf
	      
	      if(ipl == iplf-3)
		hBest_1 = 0;

	      BestChi2_2 = chi2;
	      if(flag1){
		hBest_0 = hit_array[i];
		hBest_1 = hit_array1[i];
	      }
	      else{
		hBest_0 = hit_array[i];
		hBest_2 = hit_array2[i];
	      }
	    }
	  }
	}
    	  
	t.Behead(iplf);
	ok = t.FullKF(-1);
	
	if(!ok){
	  listTrack.erase(it++); 
	  badfit = 1;
	  break;
	}
	
      }// End loop over all pivot_hits 
      
      if(badfit) break;

      if(hBest_0 && hBest_1 && hBest_2) // Do not search for a 2-hit combination if we already have a good 3-hit combination
	break;

    }//End main loop over planes
  
    if(badfit) continue;
  
    if(hBest_0 && hBest_1 && hBest_2){ // Updates the track parameters with the best combination of 3 hits    
      t.AddHit(*hBest_0);
      t.AddHit(*hBest_1);
      t.AddHit(*hBest_2);
      t.InsertMissedPlanes();
    }
    else if((hBest_0 && hBest_1) || (hBest_0 && hBest_2)){ // Updates the track parameters with the best combination of 2 hits (in the 2 most downstream planes)
	
      if(hBest_0 && hBest_1){
	t.AddHit(*hBest_0);
	t.AddHit(*hBest_1);
	t.InsertMissedPlanes();
      }
      else{ // Updates the track parameters with the best combination of 2 hits (upstream + downstream planes)
	t.AddHit(*hBest_0);
	t.AddHit(*hBest_2);
	t.InsertMissedPlanes();
      }
    }

    ok = t.FullKF(1);
       
    if(ok) {
      ok = t.FullKF(-1);
      
      if (!ok) {
	
	t.Behead(iplf);
	ok = t.FullKF(1);
	if(ok) ok = t.FullKF(-1);

	if(!ok){
	  listTrack.erase(it++); 
	  badfit = 1;
	  continue;
	}
      }
    }
    else{
    
      t.Behead(iplf);
      ok = t.FullKF(1);
      if(ok) ok = t.FullKF(-1);
      
      if(!ok){
	listTrack.erase(it++); 
	badfit = 1;
	continue;
      }
    }
 
    if((hBest_0 && hBest_1 && hBest_2) || (hBest_0 && !hBest_1 && hBest_2) || (hBest_0 && hBest_1 && !hBest_2)) 
      t.UseHitTime(); // Recalculates the time of the track after a successful addition of VD hits
    
    it++;    
  }// Endl loop over tracks
} 

// ===========================================================================
// =========================    GetIncidentTracks    =========================
// ===========================================================================
void TEv::GetIncidentTracks(int nBeamsMx, int &nBeams, TTrack **beams)
{
  // Get array of incident tracks, of max. dimension "nBeamsMx"
  // Partial (enforced only when >1 candidates) rejection of tracks w/
  // continuation into the spectrometer

  list<TTrack>::iterator ib = listTrack.begin(); nBeams = 0;
  while (ib!=listTrack.end()) {
    TTrack &b = *ib;
    if (b.Type==0x10 &&
	!(b.IFit&0x100) /* Flag tagging beams w/ continuation into spectro*/ ) {
      if  (nBeams<nBeamsMx) beams[nBeams] = &b; nBeams++;
    }
    ib++;
  }
  if (nBeams>1) {
    // If >1 beams, have they a continuation into the spectrometer?
    // In which case, they did not give rise to an interaction and cannot be
    // of any help in the reco of any other spectrometer track
    //   => Tag them and remove them from the array

    static int prvEvt = -1; int evt = ptrEvt()->getEventNumberInRun();
    if (evt==prvEvt)  // Tagging already done, and accounted for supra => return
      return;
    prvEvt = evt;

    int sign = TOpt::iCut[15]; double energyCut = .9*TOpt::dCut[4];// y cut = .9
    const TSetup& setup = TSetup::Ref();
    bool dipoleField; if (setup.MagScale[0]==0) // Scale == 0 => No field.
      dipoleField = false;
    else dipoleField = setup.MagFlag1[0]>=3;	// Dipole: OD(4) and SMC(3)
    double dYCut = dipoleField ? .25 : .1, dZCut = .1, dYZdXCut = .0010;
    if (!sign || !energyCut)
      CsErrLog::mes(elFatal,"beam's q/P (\"TraF iCut[15]/dCut[4]\") unspecified");
    nBeams = 0; ib = listTrack.begin(); while (ib!=listTrack.end()) {
      TTrack &beam = *ib;
      if (beam.Type!=0x10) { ib++; continue; }
      const THlx &bH = beam.H('d');
      bool continuation = false; list<TTrack>::iterator it = listTrack.begin();
      while (it!=listTrack.end()) {
	TTrack &t = *it; double cop = t.Hfirst(5);
	if (!(t.Type&0x7) ||// Would be better if continuation is reco'd in both
	    // SM1, for precise extrapolation to target, and SM2, for a precise
	    // momentum. But it's too demanding and we are only here trying to
	    // select candidate beams to serve as a support for back-tracking.
	    sign/cop<energyCut-sqrt(t.Hfirst(5,5))/cop/cop)
	  { it++; continue; }
	if (fabs(t.Hfirst(4)-bH(4))>dYZdXCut) {          // dZ/dX Fast rejection
	  it++; continue;
	}
	THlx tH; tH(0) = bH(0); t.Hfirst.Extrapolate(tH,true);
	if (fabs(tH(1)-bH(1))>dYCut || fabs(tH(2)-bH(2))>dZCut ||  // Y, Z, ...
	    fabs(tH(3)-bH(3))>dYZdXCut) {             // ...dY/dX Fast rejection
	  it++; continue;
	}
	TMatrixDSym m(4); TVectorD v(4), w(4);
	int i, j, k, l; for (i=0, k = 1; i<4; i++, k++) {
	  v(i)=w(i) = tH(k)-bH(k); m(i,i) = tH(k,k)+bH(k,k);
	  for (j = i+1, l = k+1; j<4; j++, l++) {
	    m(i,j)=m(j,i) = tH(k,l)+bH(k,l);
	  }
	}
	m.Invert(); w *= m; double deviation = w*v;
	if ((continuation = deviation<17)) break;  // 4(+eps)-sigma cut
	it++;
      }
      if (continuation)            // Beam HAS continuation =>...
	beam.IFit |= 0x100;              // ...Flag it
      else {                             // ...Disregard it
	if  (nBeams<nBeamsMx) beams[nBeams] = &beam; nBeams++;
      }
      ib++;
    }
  }
}

// ===========================================================================
// =========================     FitIncidentTrack    =========================
// ===========================================================================
int TTrack::Fit2IncidentTrack(int nBeams, TTrack **beams,
			      bool mode)  // 0: No "Haux" on input, 1: "Haux" already exists
{
  // QN Fit w/ Y vertex fixed @ intersection w/ incident track in Z dimension

  const TSetup &setup = TSetup::Ref();
  double fldint = setup.MagFieldInt[setup.NMags-2];
  int iBestB = -1; if (0<nBeams && nBeams<=2) {
    int iB, nBs = nBeams==1 ? 2 : 3;
    float cutB = TOpt::dCut[17], bestB, chi2;
    for (iB = 0, bestB = cutB; iB<nBs; iB++) {
      // ***** EXTRAPOLATE BEAM to Z beam == Z track *****
      // Do it twice: 1st needed to get a reasonably good estimate of
      // "Haux" in order for QN to set correctly its internal fixed value
      // of Y from the constraint set on Y @ X such that Zt==Zb.
      // So that:
      //  - Newly added, i.e. added after input "Haux" was last determined,
      //   are taken into account.
      //  - "Haux" can be undetermined on input ("mode==0").
      TTrack *beam;
      if      (iB<nBs-1)  beam = beams[iB];
      else if (iBestB>=0) beam = beams[iBestB];
      else break;
      float y0B = beam->Hfirst(1)-beam->Hfirst(3)*beam->Hfirst(0);
      float z0B = beam->Hfirst(2)-beam->Hfirst(4)*beam->Hfirst(0);
      float ypB = beam->Hfirst(3), zpB = beam->Hfirst(4);
      float xB; if (mode==0 && iB<nBs-1) {
	// Init "Haux" from "Hfirst"
	xB = (Hfirst(2)-Hfirst(4)*Hfirst(0)-z0B)/(zpB-Hfirst(4));
	Haux(5) = 	// Guess value from horizontal angle
	  (0-Hfirst(3))/(1.E-3 * 0.3 * fldint);
	Haux(0) = xB; Haux(1) = y0B+ypB*xB; Haux(2) = z0B+zpB*xB;
	Haux(3) = (Hfirst(3)-y0B)/Hfirst(0); Haux(4)= Hfirst(4);
      } 
      else {
	xB = (Haux(2)-Haux(4)*Haux(0)-z0B)/(zpB-Haux(4));
	Haux(0) = xB; Haux(1) = y0B+ypB*xB;
      }
      if (QNewtonFit(1,6) && (chi2 = Chi2tot/(NDics-4))<cutB) {
	if (iB<nBeams && chi2<bestB) {
	  iBestB = iB; bestB = chi2;
	}
      }
      else if (iB>=nBeams) iBestB = -1;
    }
  }
  return iBestB;
}

// ===========================================================================
// =========================       BuildSpacePt      =========================
// ===========================================================================
void TEv::BuildSpacePt(TTrack &t, THit *hZ, unsigned int *usedPlPatZ,
		       float dy2, float dyz, float dz2,
		       vector<TSpacePt> &spts,
		       const THlx *Haux,  // To be used instead of t.Haux.
		       bool free)         // Built on free hits only
{
  // ********** Found an extension in Z: Build Space point upon it **********

  float y, z, x0, y0, z0, yp, zp;
  if (Haux) {
    x0 = (*Haux)(0);
    y0 = (*Haux)(1); z0 = (*Haux)(2); yp = (*Haux)(3); zp = (*Haux)(4);
  }
  else {
    x0 = t.Haux(0);
    y0 = t.Haux(1); z0 = t.Haux(2); yp = t.Haux(3); zp = t.Haux(4);
  }
  static const TPlane *pZ, *pZp, *pYp, *pUp, *pVp;
  TSpacePt spt, spt4;
  vector<int>::const_iterator ihZ, ihY, ihU, ihV;

  const TSetup &setup = TSetup::Ref();
  int iplZ = hZ->IPlane; pZ = &setup.vPlane(iplZ);
  int ipli = setup.vIplFirst()[0];
  const TDetect &dZ = setup.vDetect(pZ->IDetRef);
  float sZ = dZ.Sa, cZ = dZ.Ca, xZ = dZ.X(0), wZ = dZ.Resol;
  bool isScifi = dZ.IType==22;
  // When it comes to pick up a scifi hit, we apply a time cut on the time diff
  // of the hit w.r.t. the time of the event.
  // When drift CsCluster's are updated (by "UpdateDrifts"), "eventTime" is set
  // = 0, while the actual value is backed-up into "eventTRef".
  double evtT = eventTRef ? eventTRef : eventTime;

  float Z = hZ->U; THit *hZp = 0; const TStation *&st = pZ->Station;
  float c = cZ/wZ, s = sZ/wZ, u = Z/wZ;
  double sc2Z = c*c, sscZ = s*c, ss2Z = s*s;
  int nhZ = 1; unsigned int hPatZ = 0x1;
  double scuZ  = u*c, ssuZ  = u*s, su2Z  = u*u;
  static float xZp, cZp, sZp;
  pZp = pZ->Associate; if (pZp && pZp->IFlag) {
    int jplZp = pZp->IPlane-ipli; if (1<<jplZp%32&usedPlPatZ[jplZp/32])
      // Associated plane already scanned: the combination (hit,associated) has
      return; // probably already been checked => skip
    getAssociatedHit(yp,zp,cZ,sZ,pZ,hZ,pZp,&hZp);
    if (hZp) {
      const TDetect &dZp = setup.vDetect(pZp->IDetRef);
      cZp = dZp.Ca; sZp = dZp.Sa; xZp = dZp.X(0);
      float Zp = hZp->U; Zp += (xZ-xZp)*(yp*cZp+zp*sZp);
      c = cZp/wZ; s = sZp/wZ; u = Zp/wZ;
      sc2Z += c*c; sscZ += s*c; ss2Z += s*s; 
      scuZ += u*c; ssuZ += u*s; su2Z += u*u; nhZ++; hPatZ = 0x3;
      usedPlPatZ[jplZp/32] |= 1<<jplZp%32;
    }
  }
  int projZ = pZ->IProj, projY, projU, projV;
  unsigned int sProjPatZ = 1<<projZ;

  unsigned int usedPlPatY[2], usedPlPatU[2], usedPlPatV[2];
  vector<int>::const_iterator ipY, ipU, ipV = st->IPlanes.end();
  for (ipY = st->IPlanes.begin(), projY = -1,
	 memset((void*)usedPlPatY,0,2*sizeof(unsigned int));
       ipY!=st->IPlanes.end(); ipY++) {

    // ***** LOOP ON PLANES FOR ``Y'' (``Y'' is anything but Z) *****

    const TPlane *pY = &setup.vPlane(*ipY); int jplY = *ipY-ipli;
    if ((1<<pY->IProj&sProjPatZ) ||
	projY!=-1 && pY->IProj!=projY) continue;
    if (!pY->IFlag) continue;                      // Plane OFF
    if (pY->vHitRef().empty()) continue;           // Plane empty
    const TDetect &dY = setup.vDetect(pY->IDetRef);
    float cY = dY.Ca, sY = dY.Sa,  xY = dY.X(0), wY = dY.Resol;
    y = y0+yp*(xY-x0); z = z0+zp*(xY-x0);
    if (!dY.InActive(y,z)) continue;                  // Out of active zone
    unsigned int sProjPatZY = sProjPatZ, sProjPatZYU = sProjPatZ;
    // Determine route to look for Y hits: guess value from "Haux" , width w/
    // contributions from:
    //  - Detector's resolution,
    //  - Uncertainty of measured 0x2-track extrap'd (approx.) to current X,
    //  - Uncertainty due bending in SM1 set = 2% of bending effect @ current X.
    //   There's no real reason backing this choice of 2%, only that we start
    //   from uncertainty on 1/p being ~1% in COMPASS and allow for twice more,
    //   to make for the very peculiar way (fixed Y vertex w/ not so well
    //   defined vertex) bending is determined in the present case. In any case
    //   this turns out to solve the case of evt #10 of "dstar.2003.01"
    u = t.vGuests[jplY]; float du = 25*(dy2*cY*cY+2*dyz*cY*sY+dz2*sY*sY);
    du += 25*wY*wY; du = sqrt(du);// Quadratic sum. But not below, to simplify...
    du += fabs(t.Hfirst(3)-t.Haux(3))*.02*(t.Hfirst(0)-xY)*cY;
    float umxY = u+du, umnY = u-du;
    //#define BTSPT_DEBUG
#if defined DEBUG || defined BTSPT_DEBUG
#  ifdef BTSPT_DEBUG
    static int idebug = 1;
#  else
    static int idebug = 0;
#  endif
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

      THit *hY = &vecHit[*ihY], *hYp = NULL;
      if (free && hY->Status) continue;
      float Y = hY->U; if (Y<umnY) continue; if (umxY<Y) break;
      if (isScifi &&  // Scifi: not much redundancy in compensation:...
	  fabs(hY->Time-evtT)>2) continue;   // ...strict time cut
      // Set "projY". Do this only now so that it is NOT required that
      // the first encountered non Z be efficient
      projY = pY->IProj; sProjPatZY |= 1<<projY; sProjPatZYU |= 1<<projY;
      usedPlPatY[jplY/32] |= 1<<jplY/32;
      int nhZY = nhZ+1; unsigned int hPatZY = hPatZ | 0x4;
      Y += (xZ-xY)*(yp*cY+zp*sY);
      c = cY/wY; s = sY/wY; u = Y/wY;
      double sc2ZY = sc2Z+c*c, sscZY = sscZ+s*c, ss2ZY = ss2Z+s*s;
      double scuZY = scuZ+u*c, ssuZY = ssuZ+u*s, su2ZY = su2Z+u*u;
      static float xYp, cYp, sYp;
      pYp = pY->Associate; if (pYp && pYp->IFlag) {
	int jplYp = pYp->IPlane-ipli; if (1<<jplYp%32&usedPlPatY[jplYp/32])
	  // Associated plane already scanned: the combination (hit,associated)
	  continue; // has probably already been checked => skip
	getAssociatedHit(yp,zp,cY,sY,pY,hY,pYp,&hYp);
	if (hYp) {
	  const TDetect &dYp = setup.vDetect(pYp->IDetRef);
	  cYp = dYp.Ca; sYp = dYp.Sa; xYp = dYp.X(0);
	  float Yp = hYp->U; Yp += (xZ-xYp)*(yp*cYp+zp*sYp);
	  c = cYp/wY; s = sYp/wY; u = Yp/wY;
	  sc2ZY += c*c; sscZY += s*c; ss2ZY += s*s; 
	  scuZY += u*c; ssuZY += u*s; su2ZY += u*u; nhZY++; hPatZY |= 0x8;
	}
      }
      double detZY = sc2ZY*ss2ZY-sscZY*sscZY;
      if (fabs(detZY)<1.e-5) continue;
      float yZY = (scuZY*ss2ZY-ssuZY*sscZY)/detZY;
      float zZY = (ssuZY*sc2ZY-scuZY*sscZY)/detZY;

      for (ipU = ipY, projU = -1,
	 memset((void*)usedPlPatU,0,2*sizeof(unsigned int));
	   ipU!=st->IPlanes.end(); ipU++) {

	// ***** LOOP ON PLANES FOR ``U'' *****
	// Start search from "ipY" since all its predecessors failed
	// to provide a single hit compatible with track's initial

	const TPlane *pU = &setup.vPlane(*ipU); int jplU = *ipU-ipli;
	if ((1<<pU->IProj&sProjPatZY) ||
	    projU!=-1 && pU->IProj!=projU) continue;
	if (!pU->IFlag) continue;                  // Plane OFF
	ipV = ipU;   // Remember last "ipU"
	if (pU->vHitRef().empty()) continue;       // Plane empty
	const TDetect &dU = setup.vDetect(pU->IDetRef);
	float cU = dU.Ca, sU = dU.Sa, xU = dU.X(0), wU = dU.Resol;
	float Y = hY->U, Z = hZ->U, xYU = xY, xZU = xZ;
	if (hYp) { Y += hYp->U; Y /= 2; xYU += xYp; xYU /= 2; }
	if (hZp) { Z += hZp->U; Z /= 2; xZU += xZp; xZU /= 2; }
	xYU -= xU; xZU -= xU; double det = cY*sZ-cZ*sY;
	y = (Y*sZ-Z*sY-yp*(xYU*cY*sZ-xZU*cZ*sY)-zp*sY*sZ*(xYU-xZU))/det;
	z = (Z*cY-Y*cZ-zp*(xZU*sZ*cY-xYU*sY*cZ)-yp*cZ*cY*(xZU-xYU))/det;
	if (!dU.InActive(y,z)) continue;              // Out of active zone
	// Determine route to look for U hits:
	//  - Guess value is now based on already found Z's and Y's, and
	//   slopes "yp" and "zp"
	//  - Width should have a contribution from the uncertainty on the
	//   slopes. Instead,in order to simplify, take 4 times "wU"...
	u = y*cU+z*sU; du = wU*3*(1+fabs(t.Haux(5)));
	float umxU = u+du, umnU = u-du, uZY = yZY*cU+zZY*sU;
#if defined DEBUG || defined BTSPT_DEBUG
	if (idebug) {
	  if (isMC) {
	    printf("Y %.2f,%d,%d ",hY->U,hY->IKine,hY->IKin2);
	    if (hYp) printf("Y' %.2f,%d,%d ",hYp->U,hYp->IKine,hYp->IKin2);
	  }
	  printf("u %f du %f => %f %f   yzuZY %f %f %f\n",
		 u,du,umnU,umxU,yZY,zZY,uZY);
	}
#endif	      
	for (ihU = pU->vHitRef().begin(); ihU!=pU->vHitRef().end(); ihU++) {

	  // ***** LOOP ON ``U'' HITS *****

	  THit *hU = &vecHit[*ihU], *hUp = NULL;
	  if (free && hU->Status) continue;
	  float U = hU->U; if (U<umnU) continue; if (umxU<U) break;
	  if (isScifi &&  // Scifi: not much redundancy in compensation:...
	      fabs(hU->Time-evtT)>2) continue;   // ...strict time cut
	  projU = pU->IProj; sProjPatZYU |= 1<<projU;
	  usedPlPatU[jplU/32] |= 1<<jplU%32;		
	  // Evaluate space point
	  U += (xZ-xU)*(yp*cU+zp*sU);
	  if (fabs(uZY-U)>du) continue; // It's duplicating above cut in fact

	  // ********** CANDIDATE ZYU SPACE POINT **********

	  c = cU/wU; s = sU/wU; u = U/wU;
	  spt.sc2 = sc2ZY+c*c, spt.ssc = sscZY+s*c; spt.ss2 = ss2ZY+s*s;
	  spt.nHits = nhZY+1; spt.hPat = hPatZY | 0x10;
	  spt.scu = scuZY+u*c; spt.ssu = ssuZY+u*s, spt.su2 = su2ZY+u*u;
	  static float xUp, cUp, sUp;
	  pUp = pU->Associate; if (pUp && pUp->IFlag) {
	    int jplUp = pUp->IPlane-ipli; if (1<<jplUp%32&usedPlPatU[jplUp/32])
	      // Associated plane already scanned: the combination (hit,
	      continue; // associated) has probably already been checked => skip
	    getAssociatedHit(yp,zp,cU,sU,pU,hU,pUp,&hUp);
	  }
	  if (hUp) {
	    const TDetect &dUp = setup.vDetect(pUp->IDetRef);
	    cUp = dUp.Ca; sUp = dUp.Sa; xUp = dUp.X(0);
	    float Up = hUp->U; Up += (xZ-xUp)*(yp*cUp+zp*sUp);
	    c = cUp/wU; s = sUp/wU; u = Up/wU;
	    spt.sc2 += c*c; spt.ssc += s*c; spt.ss2 += s*s;
	    spt.scu += u*c; spt.ssu += u*s; spt.su2 += u*u;
	    spt.nHits++; spt.hPat |= 0x20;
	  }
	  else if (st->IPlanes.size()>4 && // If redundancy is possible...
		   spt.nHits<4) continue;  // ...require it ...
	  double det = spt.sc2*spt.ss2-spt.ssc*spt.ssc;
	  y = (spt.scu*spt.ss2-spt.ssu*spt.ssc)/det; spt.y = y;
	  z = (spt.ssu*spt.sc2-spt.scu*spt.ssc)/det; spt.z = z;
	  spt.chi2 = (y*y*spt.sc2+z*z*spt.ss2+spt.su2+
		      2*(y*z*spt.ssc-y*spt.scu-z*spt.ssu))/(spt.nHits-2);
#if defined DEBUG || defined BTSPT_DEBUG
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

      for (projV = -1, memset((void*)usedPlPatV,0,2*sizeof(unsigned int));
	   ipV!=st->IPlanes.end(); ipV++) {

	// ***** LOOP ON PLANES FOR ``V'' *****
	// (Starting at last ipU, assuming coordinates are not intermixed
	// in a givent station)

	const TPlane *pV = &setup.vPlane(*ipV); int jplV = *ipV-ipli;
	if (1<<pV->IProj&sProjPatZYU) continue;
	if (projV!=-1 && pV->IProj!=projV) {
	  const TDetect &dV = setup.vDetect(pV->IDetRef);
	  string stDets = dV.Name.substr(0,4); stDets += string(" = ");
	  vector<int>::const_iterator jp;
	  for (jp = st->IPlanes.begin(); jp!=st->IPlanes.end(); jp++) {
	    const TPlane *pj = &setup.vPlane(*jp);
	    const TDetect &dj = setup.vDetect(pj->IDetRef);
	    if (jp!=st->IPlanes.begin()) stDets += string(",");
	    stDets += dj.Name;
	  }
	  CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "%s: Proj=%d != Z,Y,U,VProj=%d,%d,%d,%d\n"
 " => Check your \"detectors.dat\" for badly misoriented detector planes in station %s",
	    setup.vDetect(pV->IDetRef).Name.c_str(),pV->IProj,
	    projZ,projY,projU,projV,stDets.c_str());
	}
	if (!pV->IFlag) continue;                  // Plane OFF
	if (pV->vHitRef().empty()) continue;       // Plane empty
	const TDetect &dV = setup.vDetect(pV->IDetRef);
	float cV = dV.Ca, sV = dV.Sa, xV = dV.X(0), wV = dV.Resol;
	float Y = hY->U, Z = hZ->U, xYV = xY, xZV = xZ;
	if (hYp) { Y += hYp->U; Y /= 2; xYV += xYp; xYV /= 2; }
	if (hZp) { Z += hZp->U; Z /= 2; xZV += xZp; xZV /= 2; }
	xYV -= xV; xZV -= xV; double det = cY*sZ-cZ*sY;
	y = (Y*sZ-Z*sY-yp*(xYV*cY*sZ-xZV*cZ*sY)-zp*sY*sZ*(xYV-xZV))/det;
	z = (Z*cY-Y*cZ-zp*(xZV*sZ*cY-xYV*sY*cZ)-yp*cZ*cY*(xZV-xYV))/det;
	if (!dV.InActive(y,z)) continue;              // Out of active zone
#if defined DEBUG || defined BTBSP_DEBUG || defined ForeTrack_DEBUG
	float sYV = sY*cV-cY*sV, sZV = sZ*cV-cZ*sV;
	if (fabs(sZV)<.08 || fabs(sYV)<.08) {  // X-check V!=Z && V!=Y
	  CsErrLog::msg(elError,__FILE__,__LINE__,
			"%s: Proj=%d != Z,YProj=%d,%d: sZV,sYV = %f,%f",
			dV.Name.c_str(),pV->IProj,projZ,projY,sZV,sYV);
	  continue;       // ``Y'' should be anything but Z
	}
#endif
	u = y*cV+z*sV; du = wV*3*(1+fabs(t.Haux(5)));
	float umxV = u+du, umnV = u-du, uZY = yZY*cV+zZY*sV;
#if defined DEBUG || defined BTSPT_DEBUG
	if (idebug) {
	  if (isMC) {
	    printf("Y %.2f,%d,%d ",hY->U,hY->IKine,hY->IKin2);
	    if (hYp) printf("Y' %.2f,%d,%d ",hYp->U,hYp->IKine,hYp->IKin2);
	  }
	  printf("v %f dv %f => %f %f   yzvZY %f %f %f\n",
		 u,du,umnV,umxV,yZY,zZY,uZY);
	}
#endif	      
	for (ihV = pV->vHitRef().begin(); ihV!=pV->vHitRef().end(); ihV++) {

	  // ***** LOOP ON ``V'' HITS *****

	  THit *hV = &vecHit[*ihV], *hVp = NULL;
	  if (free && hV->Status) continue;
	  float V = hV->U; if (V<umnV) continue; if (umxV<V) break;
	  projV = pV->IProj;
	  usedPlPatV[jplV/32] |= 1<<jplV%32;		
	  // Evaluate space point
	  V += (xZ-xV)*(yp*cV+zp*sV);
	  if (fabs(uZY-V)>du) continue; // It's duplicating above cut in fact

	  // ********** CANDIDATE ZYV SPACE POINT **********

	  static float xVp, cVp, sVp, Vp;
	  pVp = pV->Associate; if (pVp && pVp->IFlag) {
	    int jplVp = pVp->IPlane-ipli; if (1<<jplVp%32&usedPlPatV[jplVp/32])
	      // Associated plane already scanned: the combination (hit,
	      continue; // associated) has probably already been checked => skip
	    getAssociatedHit(yp,zp,cV,sV,pV,hV,pVp,&hVp);
	    if (hVp) {
	      const TDetect &dVp = setup.vDetect(pVp->IDetRef);
	      cVp = dVp.Ca; sVp = dVp.Sa; xVp = dVp.X(0);
	      Vp = hVp->U; Vp += (xZ-xVp)*(yp*cVp+zp*sVp);
	    }
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
		projV==setup.vPlane(spt3.hs[4]->IPlane).IProj) continue;
	    float u = spt3.y*cV+spt3.z*sV;
	    if (fabs(u-V)>3*dV.Resol) continue;

	    // ***** CANDIDATE ZYUV SPACE POINT *****

	    c = cV/wV; s = sV/wV; u = V/wV;
	    spt4.sc2 = spt3.sc2+c*c, spt4.ssc = spt3.ssc+s*c;
	    spt4.ss2 = spt3.ss2+s*s;
	    spt4.nHits = spt3.nHits+1; spt4.hPat = spt3.hPat | 0x40;
	    spt4.scu = spt3.scu+u*c; spt4.ssu = spt3.ssu+u*s;
	    spt4.su2 = spt3.su2+u*u;
	    if (hVp) {
	      c = cV/wV; s = sV/wV; u = V/wV;
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
#if defined DEBUG || defined BTSPT_DEBUG
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
	    c = cV/wV; s = sV/wV; u = V/wV;
	    spt.sc2 = sc2ZY+c*c, spt.ssc = sscZY+s*c; spt.ss2 = ss2ZY+s*s;
	    spt.nHits = nhZY+1; spt.hPat = hPatZY | 0x10;
	    spt.scu = scuZY+u*c; spt.ssu = ssuZY+u*s, spt.su2 = su2ZY+u*u;
	    if (hVp) {
	      c = cVp/wV; s = sVp/wV; u = Vp/wV;
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
}

// ===========================================================================
// =========================       BackTrack2SI      =========================
// ===========================================================================
bool TEv::BackTrack2SI(TTrack &t)
{
  // *************** REASSESS SI HITS (in ZONE 0x1) ***************
  // - The routine is expecting argument track to indeed own SI hits (it's
  //  checked infra).
  // - Also assumed: SIs come as first in the list of detectors of zone 0x1.
  // - There's a track selection: (tentatively) only proceed w/ those tracks
  //  that are liable to have ghost SIs. The selection scheme is based on the
  //  average chi2 contribution from SIs compared to average from other
  //  detectors. It's (as of 2008/01) expressed awkwardly, given that it largely
  //  depends upon the setting of "TOpt::dCut[84]", which resets uncertainties
  //  in SI hits to make for misalignment (or upon misalignment, if "dCut[84]"
  //  happens to be disabled). It besides accepts all cases where there are only
  //  few (hence unreliable) SI hits.
  // - 2nd selection: retain only tracks where we have all what's needed to
  //  do something. The reassessment algorith is relying on a QN fit of the
  //  track's segment in zone 0x1 stripped of its SIs. 
  // - All SI hits are considered, *including* those already used (in some other
  //  other track). The desambiguation will eventually be performed by
  //  "TEv::CleanTrackListSI".
  // - Is non satisfactory the handling of cases where 2 tracks are close to one
  //  another, typically the case of an e+e- pair. There, the 2 tracks should be
  //  processed simultaneously, and the decision (as to whether a given hit
  //  belongs to a particular track) taken on view of the sum of the chi2 of the
  //  2 tracks.
  // - Returns true if the track has been modified.

  //          *************** INITIALIZATIONS... ***************
  // ...I) On a per job basis
  const TSetup &setup = TSetup::Ref();
  int iplF = setup.vIplFirst()[0], iplL = setup.vIplLast()[0], ipl;
  int nPats = (iplL-iplF)/32+1;	unsigned int hPats[nPats];
  //  Indices of SI spectrometer (as opposed to beam) TStations: 2 TStations are
  // expected (although "BackTrack2SI" probably works w/ more, or less, than 2).
  //  Indices of the upstream-most tracking TStation other than SI (we will
  // require that track has hits in any of these, before proceeding w/ the SIs).
  const vector<TStation> &stations = setup.vStation();
  static int iS0SI = -1, iS1SI, iS0MM, iPlLSI, iPl0MM;
  static unsigned int iS0TBs; // Indices of the 1st TStation of each type.
  if (iS0SI<0) {
    int iS; unsigned int types;
    for (iS = 0, iS1SI=iS0MM = -1, iS0TBs=types = 0; iS<(int)stations.size();
	 iS++) {
      const TStation &station = stations[iS];
      int ipl = station.IPlanes[0]+station.IPlanes.size()-1; // TStation's last
      if (iplL<ipl) break; if (ipl<iplF) /* Skip upstream of target */ continue;
      int type = station.Type; if (type==21) {
	if      (iS0SI<0) iS0SI = iS;
	else if (iS1SI<0) iS1SI = iS;
	else {            iS0SI = -1; break; } // No more than 2 SI TStations
	iPlLSI = ipl;
      }
      else if (type==15 /* DC */ || type==22 /* FI */ ||
	       26<=type && type<=32 /* G/MM/P */) {
	unsigned int b = 1<<type-15; if (!(b&types)) {
	  iS0TBs |= 1<<iS; types |= b;
	  if (type==27) {  // TStation::Type 27 covers both MMs and MPs
	    iS0MM = iS; iPl0MM = ipl;
	  }
	}
      }
    }
    if (iS0SI<0 || iS1SI<0) CsErrLog::mes(elFatal,
 "# of Si stations downstream of target !=2, while \"TraF ReMode[44]\" booked");
  }
  // ...II) On a per track basis
  int paraxial; // Paraxial means track travelling in the VSAT throughout...
    // ...It's flagged by requiring VSAT hits ("TTrack::Scifi") and high enough
    // momentum that it's not bend away from (V)SAT by SM1. It implies a lower
    // redundancy => Few constraints on trajectory once track's stripped of SIs.
    // Therefore, in this particular, we extend QN fit down to RICH.
  if (t.Scifi&0x101 && fabs(t.Hfirst(5))<.15)
    paraxial = fabs(t.Hfirst(5))<.02 ? 2 : 1;
  else
    paraxial = 0;
  double xRICH = 750;

#define BTSI_DEBUG
#ifdef BTSI_DEBUG
  static int idebug = 0;
  if (idebug) t.DumpHits(0);
#endif

  //           *************** TRACK SELECTION ***************
  //   I) CHECK that SI HITS are CONSTRAINED (#hits, #projs/#stations)
  //  II) CHECK chi2 OF SI HITS
  //     - Their contribution to the overall track chi2
  //     - (if numerous enough) Their ``internal'' chi2.
  //     - (if numerous enough) Their ``space point'' chi2's.
  // III) ENOUGH HITS that the TRACK STRIPPED of its SIs CAN BE FITTED
  double sChi2SI; int nSIs, nSIStations, nH0x1, nPix, nHs, lastSI;
  const TDetect *worstD; double worstChi2Incr;  const TStation *sPrv;
  double sptChi2s[2] = {0,0};
  static double sc2, ssc, ss2, scU, ssU, sU2, x0; static int s1;
  double yp = t.Hfirst(3), zp = t.Hfirst(4);
  list<int>::iterator ih; map<int, double>::iterator idx;
  for (ih = t.lHitPat.begin(), idx = t.mChi2.begin(), sChi2SI = 0,
	 nSIs=nSIStations=nH0x1=nPix=nHs = 0, sPrv = 0, lastSI = -1,
	 worstChi2Incr = t.Chi2tot*.25, worstD = 0; idx!=t.mChi2.end(); ih++) {
    if (*ih<0) continue; THit &h = vecHit[*ih];
    ipl = (*idx).first; double chi2Incr = (*idx).second; idx++;
    if (ipl!=h.IPlane) CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "Event #%d Track #%d(0x%x): Inconsisten hit and chi2 maps",event,t.Id,t.Type);
    const TDetect &d = setup.vDetect(ipl);
    if (chi2Incr>worstChi2Incr) { worstChi2Incr = chi2Incr; worstD = &d; }
    if  (paraxial) { if (d.X(0)>xRICH) break; }
    else           { if (ipl>iplL) break; }
    if (d.IType==21) {
      sChi2SI += chi2Incr; nSIs++; lastSI = ipl;
      const TStation *st = setup.vPlane(ipl).Station; if (st!=sPrv) {
	if (sPrv) {
	  if (s1>2) {
	    double det = sc2*ss2-ssc*ssc;
	    double Y = (scU*ss2-ssU*ssc)/det, Z = (ssU*sc2-scU*ssc)/det;
	    double chi2 = sqrt((Y*Y*sc2+Z*Z*ss2+sU2+2*(Y*Z*ssc-Y*scU-Z*ssU))/
			       s1);
	    sptChi2s[sPrv->IStation==iS0SI ? 0 : 1] = chi2/.001;
	  }
	  else sptChi2s[sPrv->IStation==iS0SI ? 0 : 1] = -1;
	}
	sc2=ssc=ss2=scU=ssU=sU2= 0; s1 = 0; x0 = d.X(0);
	sPrv = st; nSIStations++;
      }
      double w = h.SigU-TOpt::dCut[85]; // Undo SI artificial worsening.
      if (w<.001) w = .001; // Still, keep a lower bound = 10 microns
      w /= .001; // Rescale
      double c = d.Ca/w, s = d.Sa/w, dx = d.X(0)-x0, U = h.U/w-(yp*c+zp*s)*dx;
      s1++;
      sc2 += c*c; ssc += s*c; ss2 += s*s; scU += c*U; ssU += s*U; sU2 += U*U;
    }
    else {
      if (ipl<=iplL) {
	nH0x1++;
	if (setup.vPlane(ipl).IFlag&0x30) { // P-pixels: provide also v-axis.
	  nH0x1++; nPix++;
	}
      }
      nHs++;
    }
  }
  if (sPrv && sPrv->Type==21) {
    if (s1>2) {
      double det = sc2*ss2-ssc*ssc;
      double Y = (scU*ss2-ssU*ssc)/det, Z = (ssU*sc2-scU*ssc)/det;
      double chi2 = sqrt((Y*Y*sc2+Z*Z*ss2+sU2+2*(Y*Z*ssc-Y*scU-Z*ssU))/s1);
      sptChi2s[sPrv->IStation==iS0SI ? 0 : 1] = chi2/.001;
    }
    else sptChi2s[sPrv->IStation==iS0SI ? 0 : 1] = -1;
  }
  if (nSIs==0) {
    // Argument tracks are expected to own SI hits, cf. call to "BackTrack2SI"
    // in "TEv::TracksFit2".
    CsErrLog::msg(elError,__FILE__,__LINE__,
		  "Track #%d starts w/ a SI TPlane has no SI hit",t.Id);
    return false;
  }
  bool unConstrained = nSIs<4 || nSIStations>1 &&
    (sptChi2s[0]<0 || sptChi2s[1]<0) /* I.e. #hits in any station <3 */;

  if (paraxial && // Paraxial case: there are few hits apart from SIs. If one...
      // ...of those few turns out to be the worst hit, given that it's upstream
      // of SIs, it hence, in a backward KF fit, cannot owe them its bad chi2
      // increment, but can induce bad chi2 in them on the contrary. It's the
      // very hit that one should try to clean away instead of checking the
      // SIs. Such a cleaning the worst hit is attempted later on in the track
      // reconstruction (cf. "TracksFit2"). No new re-evaluation of the SIs,
      // that would take place thereafter, is foreseen, as of 2010/02. Yet the
      // present one, will probably go awry, because of this very worst hit...
      worstD && worstD->IType!=21 &&
      !unConstrained) return false;  // ...=> Let's give up. 

  double sChi2SAT = t.Chi2tot-sChi2SI; int nSATs = t.NDFs-nSIs;
  bool reAssess =                      // Do re-assess SIs if...
    sChi2SI/nSIs>2.5*sChi2SAT/nSATs || // ...significant contribution to chi2
    // (Note: The presence of random pickup from the SIs will upset the whole
    // fit. So that we have to resort to this kind of very loose chi2 cut. But
    // even this is not completely satisfactory given, that its impact will
    // depend upon the setting of "TOpt::dCut[85]", which resets uncertainties
    // in SI hits to make for misalignment (or upon misalignment, if
    // "dCut[84]" happens to be disabled).)
    sptChi2s[0]>1.8 || sptChi2s[1]>1.8 || // ...or if large ``space pt'' chi2's
    // (Note: This second criterion no longer depends upon  the setting of
    // "TOpt::dCut[85]". It was introduced lately. May utlimately superseed the
    // first one, instead of merely supplementing it.)
    unConstrained;  // ...or SI unconstrained (meaning possible random pickup).
  bool checkIntChi2 =  // ...else let's further examine ``internal'' SI chi2
    sChi2SI/nSIs>3 &&  // (if contribution to chi2 sizeable and...
    nSIs>=6;           // ...#SIs large enough that a SI track can be defined)
  if (!reAssess && !checkIntChi2) return false;

  TTrack *t0x1 = 0; static double chi2SAT; int TBS0 = 0;
  if (nH0x1>8 ||
      paraxial && nHs>=7-paraxial /* e.g. MP00XY+GP01XY+FI05XY */) {
    t0x1 = new TTrack();
    for (ih = t.lHitPat.begin(); ih!=t.lHitPat.end(); ih++) {
      if (*ih<0) continue; THit &h = vecHit[*ih]; int ipl = h.IPlane;
      const TDetect &d = setup.vDetect(ipl);
      if (d.IType!=21 && !reAssess) {
	t0x1->QuickKF(1,0); double intChi2 = t.Chi2tot/(t.NHits-4);
	if (TOpt::dCut[85]) {
	  double f = (d.Resol+TOpt::dCut[85])/d.Resol; intChi2 *= f*f;
	}
	if (intChi2<3) return false;
	else           reAssess = true;
      }
      if (paraxial) { if (d.X(0)>xRICH) break; }
      else          { if (ipl>setup.vIplLast()[0]) break; }
      t0x1->AddHit(h);
      // Is there any hit in the upstream-most, non-SI tracking detectors?
      int iS = setup.vPlane(ipl).Station->IStation; if (1<<iS&iS0TBs) TBS0 = 1; 
    }
    for (int i = 0; i<6; i++) t0x1->Hfirst(i) = t.Hfirst(i);
    bool t0x1OK = t0x1->QNewtonFit(1,2); if (t0x1OK) {
      double chi20x1 = t0x1->Chi2tot/(t0x1->NDics-4); t0x1->Behead(lastSI);
      if ((t0x1OK = t0x1->QNewtonFit(1,2))) {
	chi2SAT = t0x1->Chi2tot/(t0x1->NDics-4);
	t0x1OK = chi2SAT<3 && chi2SAT<chi20x1*1.2 ||
	  paraxial && chi2SAT<6;
      }
    }
    if (!t0x1OK) { delete t0x1; t0x1 = 0; }
  }
  else if (t.Chi2tot/(t.NDFs-5)<TOpt::dCut[16]) return false;

  int nActives = -1; if (t0x1) {
    //    ********** BACKTRACK, PICKUP SIs and FIT **********

#ifdef DEBUG
    if (idebug) t0x1->DumpHits(1);
#endif
    // Intercept and slopes taken from QNewton output, viz. "Haux"
    // Will be used to determine whether inActive.
    float x0 = t0x1->Haux(0), y0 = t0x1->Haux(1), z0 = t0x1->Haux(2);
    float yp = t0x1->Haux(3), zp = t0x1->Haux(4);

    // Covariant matrix: a straight extrapolation will be used.
    bool getKFCov = fabs(t.Hfirst(5))>.06; if (getKFCov) {
      // At low momentum, get the cov matrix from a full KF fit, so that the
      // impact of MS is taken into account.
      t0x1->QuickKF(1,0);
      t0x1->Hlast(5) = t.Hlast(5);
      t0x1->Hlast(5,5) = 1.E-20; // Fix momentum for the time of the Kalman fit
      getKFCov = t0x1->FullKF(-1);
    }
    if (!getKFCov) {
      t0x1->QuickKF(-1,0);
      t0x1->Hfirst(5) = t.Hfirst(5);
      t0x1->Hfirst(5) = t.Hfirst(5); t0x1->Hfirst(5,5) = t.Hfirst(5,5);
    }
    bool initErr = true; static float dy2, dyz, dz2;

    int added, iter, iterMx = 1; unsigned int sPat, actives;
    for (iter = 0, added = 0, actives=sPat = 0; iter<=iterMx; iter++) {
      // ***** 2 ITERATIONS: so as to consider the (2) SI TSations in 2
      // different ways, viz.: i) building a space point on a picked up Y hit,
      // ii) picking up individual hits. Method (i) is first executed, on all
      // TStations. If no space point found for a given TSation (which mishap is
      // remembered in "sPat"), method (ii) is then tried. 
      for (int station = iS0MM; station>=iS0SI; station--) {
	if (station==iS0MM) {// Upstream-most MM...
	    // ...only if no other upstream-most TStation yet and we will use
	    // only one iteration (w/ method based on individual hits).
	  if (TBS0 || iter) continue;
	}
	else if (station!=iS0SI && station!=iS1SI) continue;
	if (sPat&1<<station) continue; // Already 
	const TStation &s = stations[station];
	int iplFS = s.IPlanes.front(), iplLS = s.IPlanes.back();
	bool doSpacePt = iter==0 && station!=iS0MM;
	THit *altHits[8]; // MM|P stations, including pixels, may have up to 8 planes
	int joined, free, kpl;
	for (ipl = iplLS, free=joined = 0, kpl = -1,
	       memset((void*)altHits,0,sizeof(altHits)); ipl>=iplFS; ipl--) {
	  const TPlane  &p = setup.vPlane(ipl); if (!p.IFlag/* OFF */) continue;
	  const TDetect &d = setup.vDetect(ipl);

	  //            ********** LOOP OVER SI (or MM|P) TPlanes **********

	  int jpl = ipl-iplF; kpl++;
	  float cu = d.Ca, su = d.Sa;
	  float yy = y0+yp*(d.X(0)-x0), zz = z0+zp*(d.X(0)-x0);

	  if (!d.InActive(yy,zz)) continue;              // ***** IN ACTIVE AREA
	  if (d.IType==21) actives |= 1<<jpl;
	  if (doSpacePt && p.IProj!=1) continue;
	  //       ***** UNCERTAINTY of GUESTIMATE for CURRENT PLANE  *****
	  if (initErr) {  // Update cov matrix
	    THlx hSI; hSI(0) = d.X(0); t0x1->Hfirst.Extrapolate(hSI);
	    dy2 = hSI(1,1); dz2 = hSI(2,2); dyz = hSI(1,2); initErr = false;
	  }
	  float du = sqrt(dy2*cu*cu+2*dyz*cu*su+dz2*su*su);
	  float cut;                     // ***** CUT on RESIDUAL HIT-GUESTIMATE
	  // The cutting scheme (cut value and its dependence upon context, viz.
	  // whether in X/U or not, multiplicity) is an empirical recipe, only
	  // adjusted so far on a few events from 08W33 "evtDump.LR-69853.raw".
	  if (station==iS0SI) {
	    cut = 4*du+3*d.Resol; if (cu>.2 /* I.e. X||U planes */) cut += 2*du;
	  }
	  else cut = 3*du+3*d.Resol;
	  float guest = t0x1->vGuests[jpl], umn = guest-cut, umx = guest+cut;
	  THit *hbest; int mult; float best; vector<TSpacePt> spts;
	  vector<int>::const_iterator ihit;
	  for (ihit = p.vHitRef().begin(), mult = 0, hbest = 0, best = cut;
	       ihit!=p.vHitRef().end(); ihit++) {
	    THit &h = vecHit[*ihit];//  ***** LOOP ON ALL HITS (also used ones)
	    float u = h.U; if (u<umn) continue; if (u>umx) break; 
	    mult++; float diff = fabs(u-guest);
	    if (diff<cut && doSpacePt) {
	      //      ***** BUILD SPACE POINT BASED on VERTICAL HIT *****
	      memset((void*)hPats,0,nPats*sizeof(unsigned int));
	      BuildSpacePt(*t0x1,&h,hPats,dy2,dyz,dz2,spts,0,false);
	    }
	    //         ***** else INDIVIDUAL HITS: DETERMINE BEST HIT *****
	    else if (diff<best) {
	      if (hbest && best<cut/3) altHits[kpl] = hbest;
	      best = diff; hbest = &h;
	    }
	    else if (diff<cut/3 &&
		     (!altHits[kpl] || diff<fabs(altHits[kpl]->U-guest)))
	      altHits[kpl] = &h;
	  }
	  if (doSpacePt && spts.size()!=0) {// ********** SPACE POINT **********
	    //                ***** Re-EVALUATE CHI2 *****
	    // (For uncertainties on SI hits have been artificially worsened.)
	    int ispt, nspts = spts.size(); double x0 = d.X(0);
	    for (ispt = 0; ispt<nspts; ispt++) {
	      TSpacePt &spt = spts[ispt];
	      int i; double sc2, ssc, ss2, scU, ssU, sU2, s1;
	      for (i = 0, sc2=ssc=ss2=scU=ssU=sU2=s1 = 0; i<8; i += 2)
		if (spt.hPat&1<<i) {
		  THit *h = spt.hs[i]; int kpl = h->IPlane;
		  const TDetect &dk = setup.vDetect(kpl);
		  double w = h->SigU-TOpt::dCut[85]; // Undo SI worsening.
		  if (w<.001) w = .001;// Still, keep a lower bound = 10 microns
		  w /= .001; // Rescale
		  double c = dk.Ca/w, s = dk.Sa/w, dxk = dk.X(0)-x0;
		  double U = h->U/w-(yp*c+zp*s)*dxk; s1 += 1; 
		  sc2 += c*c; ssc += s*c; ss2 += s*s;
		  scU += c*U; ssU += s*U; sU2 += U*U;
		}
	      double det = sc2*ss2-ssc*ssc;
	      double Y = (scU*ss2-ssU*ssc)/det, Z = (ssU*sc2-scU*ssc)/det;
	      double chi2 = sqrt((Y*Y*sc2+Z*Z*ss2+sU2+2*(Y*Z*ssc-Y*scU-Z*ssU))/
				 s1);
	      spt.chi2 = chi2/.001;
	    }
	    //         ***** DETERMINE BEST SPACE POINT *****
	    const TSpacePt *spt0; int nHMx; double chi2Mn;
	    for (ispt=nHMx = 0, chi2Mn = 5, spt0 = 0; ispt<nspts; ispt++) {
	      const TSpacePt &spt = spts[ispt];
	      if (spt.chi2<5 && (spt.nHits>nHMx ||
				 spt.nHits==nHMx && spt.chi2<chi2Mn)) {
		nHMx = spt.nHits; chi2Mn = spt.chi2; spt0 = &spt;
	      }
	    }
	    if (spt0) {  // ***** => ADD SPACE POINT *****
	      for (int i = 0; i<8; i++) if (spt0->hPat&1<<i) {
		THit *h = spt0->hs[i]; joined++;
		if (h->sTrackID().empty()) free++;
		else if (h->sTrackID().size()==1 &&
			 *h->sTrackID().begin()==t.Id) free++;
		t0x1->AddHit(*h);
	      }
	    } // Valid space point exists
	  }
	  else if (hbest && (mult==1 || best<cut/2)) {
	    // ***** INDIVIDUAL HITS:ACCOUNT for MULTIPLICTY (of hits w/in cuts)
	    joined++;
	    if (hbest->sTrackID().empty()) free++;
	    else if (hbest->sTrackID().size()==1 &&
		     *hbest->sTrackID().begin()==t.Id) free++;
	    t0x1->AddHit(*hbest);
	  } // End Best hit or best Space point
	} // End loop over planes
	int jpl; for (jpl = 0, nActives = 0; jpl<iPlLSI-iplF; jpl++)
		   if (1<<jpl&actives) nActives++;
	//                                         ***** EVALUATE PICKED UP HITS
	if (joined>=3 &&                // Require one space point...
	    (free>=2 ||                 // ...and 2 free hits unless...
	     // (Note: 2 free hits because that corresponds to the case where
	     // current track is close to another track in one dimension (hence
	     // ambiguity for 2 coord's, (X,U) or (Y,V), given that coord's w/in
	     // these pairs are only separated by 5 degrees) and isolated in
	     // the other dimension (hence no ambiguity in complementary pair.)
	     station==iS0SI && // ...station==0 (=> anyway more crowded) and...
	     (added || // ...SIs added (=> backtracking more reliable) except...
	      nActives!=8))) { // ...if not w/in acceptance of both stations.
	  if (t0x1->QNewtonFit(1,2)) {  // ***** REFIT w/ FIXED P
	    added += joined; 
	    double bestChi2 = t0x1->Chi2tot/(t0x1->NDics-4); int best, kpl;
	    int ambig; for (kpl = 0, best = -2, ambig = -1; kpl<8;
			    kpl++) {              // ***** REFIT ALTERNATIVES...
	      // ...and determine:
	      //  - Best alternative on view of chi2. (Alternatives are examined
	      //   one at a time, and not in combination, for simplicity sake.)
	      //  - Whether there is ambiguity, i.e. best close to next best. If
	      //   a single ambiguity, remember it.
	      THit *halt = altHits[kpl]; if (halt) {
		if (best==-2) best = -1;
		const THit *hPrv = t0x1->ReplaceHit(halt); if (hPrv) {
		  t0x1->QNewtonFit(1,2);
		  double chi2 = t0x1->Chi2tot/(t0x1->NDics-4);
		  if (chi2<bestChi2) {
		    if (chi2>bestChi2*0.95) ambig = ambig==-1 ? kpl : -2;
		    bestChi2 = chi2; best = kpl;
		  }
		  else if (bestChi2>chi2*0.95) ambig = ambig==-1 ? kpl : -2;
		  t0x1->ReplaceHit(hPrv);
		}
	      }
	    }
	    if (ambig>=0 && station==iS0SI) {        // ***** CASE TStation 0...
	      //                         ***** ...RAISE AMBIGUITY by FULL KF FIT
	      TTrack tp(t); tp.Behead(lastSI);
	      for (ih = t0x1->lHitPat.begin(); ih!=t0x1->lHitPat.end(); ih++) {
		if (*ih<0) continue; THit &h = vecHit[*ih];
		if (setup.vDetect(h.IPlane).IType!=21) break;
		tp.AddHit(h);
	      }
	      if (best>=0 && best!=ambig) tp.ReplaceHit(altHits[best]);
	      tp.FullKF(-1); double chi2 = tp.Chi2tot/(tp.NDFs-5);
	      tp.ReplaceHit(altHits[ambig]);
	      tp.FullKF(-1); if (tp.Chi2tot/(tp.NDFs-5)>chi2) {
		ambig = -2; if (best==ambig) best = -1;
	      }
	    }
	    if (bestChi2>6 && (!paraxial ||// ***** CUT ON CHI2 & CHI2 INCREMENT
			       bestChi2>3) &&
		// (Note: the cut on the chi2 increment is tricky, for the
		// increment due to SIs could be large in case of misalignment,
		// particularly so if not alleviated by "TOpt::dCut[85]".)
		bestChi2*(t0x1->NDics-4)>
		chi2SAT*(nHs-4)+(3+(paraxial==1?2:0))*added) {
	      // Erase from "t0x1" track those hits just picked up, so as to
	      // either clean the ground for next or preserve previous iter.
	      t0x1->SubHits(iplFS,iplLS); added -= joined; continue;
	    }
	    else if (best!=-2) {// Were ambiguities considered, hence fitted...
	      if (best>=0)  t0x1->ReplaceHit(altHits[best]);
	      if (ambig>=0) t0x1->ReplaceHit(altHits[ambig]);
	      if (!t0x1->QNewtonFit(1,2)) {           // ...=> QN refit
		// Should never occur. Yet, to be on the safe side...
		t0x1->SubHits(iplFS,iplLS); added -= joined; continue;
	      }
	    }
	    sPat |= 1<<station;
#ifdef DEBUG
	    if (idebug) t0x1->DumpHits(1);
#endif
	  }
	  else          // Picked-up hits yield track that cannot be fitted
	    joined = -1;// => Flag the case ("joined=-1"), that hits be erased.
	  if (joined && station==iS0MM)
	    printf("Evt #%d: %d MM hits to TTrack #%d\n",event,joined,t.Id);
	}
	else if (joined)// Too few (or questionable) hits picked-up: unreliable.
	  joined = -1;  // => Flag the case ("joined=-1"), that hits be erased.
#ifdef DEBUG
	if (idebug) t0x1->DumpHits(0);
#endif
	if (joined<=0) {   // ***** HIT PICK-UP UNSUCCESSFUL IN CURRENT TStation
	  if (station==iS0MM) break;
	  if (added==0 && iter==iterMx) {
	    if    (station==iS0SI)// No hit picked up in "station==1" either...
	      break;                  //  ...=> Give up.
	    else if (nActives>1)  // Unsuccessful while w/in active of >1 SIs...
	      break; // ...=> Give up. (So as not to waste CPU in "station==0".)
	  }
	  if (!joined) continue;
	  // Erase from auxiliary "t0x1" track those hits just picked up, so as
	  // to either clean the ground for next or preserve previous iteration.
	  t0x1->SubHits(iplFS,iplLS);
#ifdef DEBUG
	  if (idebug) t0x1->DumpHits(0);
#endif
	}
      } // End loop on 2 (or 3, if MM) stations
    } // End of 2 iterations
    if (added) {                       // ***** SUCCESS! => ADD SI HITS TO TRACK
      // Do it by scanning in parallel current track "t" and auxiliary "tx0x1".
      int modified;// Set "modified" if any diff => current track needs a refit.
      list<int>::iterator jh, jp;
      for (ih = t0x1->lHitPat.begin(), modified = 0, jh = t.lHitPat.begin(),
	     jp = t.lPlnRef.begin(); ih!=t0x1->lHitPat.end(); ih++) {
	if (*ih<0) continue; THit &hi = vecHit[*ih]; int ipl = hi.IPlane; 
	if (iPl0MM<ipl) break;
	int jpl = -1; while (jp!=t.lPlnRef.end() && (jpl = *jp)<ipl) {
	  if (*jh>=0) {
	    const THit &hj = vecHit[*jh];
	    const_cast<THit&>(hj).eraseTrackID(t.Id);   // Cast away constness and cancel hit -> track reference
	    if (!modified) modified = 1; t.NHits--; t.NDFs--; *jh = -1;
	  }
	  if (ih==t0x1->lHitPat.begin()) {
	    // Remove all references upstream of the 1st SI hit retained.
	    jp = t.lPlnRef.erase(jp);      // Erase plane ref.
	    jh = t.lHitPat.erase(jh);      // Erase hit ref.
	  }
	  else { jp++; jh++; }
	}
	if (jpl<0) break; // Should never happen.
	if (jpl==ipl) {
	  if (*jh!=*ih) {
	    if (!modified) modified = 1;
	    if (*jh>=0) {
	      const THit &hj = vecHit[*jh];
	      const_cast<THit&>(hj).eraseTrackID(t.Id); // Cast away constness and cancel hit -> track reference
	    }
	    else { t.NHits++; t.NDFs++; }
	    *jh = *ih;
	    const_cast<THit&>(hi).addTrackID(t.Id);     // Cast away constness and add    hit -> track reference
	  }
	}
	else {
	  jp = t.lPlnRef.insert(jp,ipl);      // Insert new plane ref.
	  jh = t.lHitPat.insert(jh,*ih);      // Insert new hit ref.
	  modified = 2; t.NHits++; t.NDFs++;
	  const_cast<THit&>(hi).addTrackID(t.Id);       // Cast away constness and add    hit -> track reference
	}
	jp++; jh++;
      }
      while (jp!=t.lPlnRef.end()) {
	// Now turn to those SI hits of current track "t" not in "t0x1".
	// (Note: given that "added!=0", there is at least on hit in "t0x1",
	// which we must have then encountered supra. => No bothering here about
	// plane reference that would be left stranded upstream of 1st hit.)
	if (*jh>=0) {
	  const THit &hj = vecHit[*jh]; int jpl = hj.iPlane;
	  if (iPlLSI<jpl) break;
	  const_cast<THit&>(hj).eraseTrackID(t.Id);   // Cast away constness and cancel hit -> track reference
	  if (!modified) modified = 1; t.NHits--; t.NDFs--; *jh = -1;
	}
	jp++; jh++;
      }
      delete t0x1;
      if (modified) {
	if (modified==2) t.InsertMissedPlanes(); return true;
      }
      else return false;
    }
    else { delete t0x1; t0x1 = 0; }
  }

  //                  ***** NOTHING WE COULD DO... *****
  if (t0x1) 
    CsErrLog::msg(elFatal,__FILE__,__LINE__,
		  "Evt #%d, track #%d: inconsistency",event,t.Id);
  if (nSIs) {
    if ((nH0x1>6 || nPix>=2) &&// ...Yet if enough hits for a decent 0x1 segment
	(sChi2SI/nSIs>3*sChi2SAT/nSATs ||           // and SIs real bad...
	 nActives>=0 && (2*nSIs<nActives || nSIs<3) ||  // ...or SIs real few
	 nActives==-1 && nSIs<=3)) {  // ...or SI few and we couldn't check them
      t.Behead(lastSI);                      // ***** STRIP TRACK of its SI HITS
      CsErrLog::msg(elError,__FILE__,__LINE__,
 "Track #%d (chi2/NDF=%.2f, chi2 per SI=%.2f) stripped of its SIs",
		    t.Id,t.Chi2tot/(t.NDFs-5),sChi2SI/nSIs);
      return true;
    }
    return false;
  }
  return false;
}
// ===========================================================================
// =========================      Account4RIPipe     =========================
// ===========================================================================
bool Account4RIPipe(TTrack &t) {
  if (TOpt::RICHPipe[5]) return false;

  //  If the RICH pipe is of heavy material (the case up to 2011), the
  // tracks undergoes high ms when passing through. => A fit w/o ms taken into
  // account, such as "QNewtonFit" or "QuickKF" cannot work.
  //  A way out is to restrict the fittting to that part of the track that's
  // upstream of RICH in the building stage of the "BackTrackSAS" method
  // (which is possible if the track has 2 space points at least upstream of
  // RICH) and finalize w/ "FullKF".
  //  A better solution would be to rely on "FullKF" throughout, but this
  // would require a complete rewriting of the whole method.
  double *p = TOpt::RICHPipe, l = p[3], r = p[4], r2 = r*r; int io[2];
  const THlx &Hf = t.H('u');
  float x0 = Hf(0), y0 = Hf(1), z0 = Hf(2), yp = Hf(3), zp = Hf(4);
  for (int ud = 0; ud<2; ud++) {
    double dx = p[0]+(2*ud-1)*l-x0, y = y0+dx*yp-p[1], z = z0+dx*zp-p[2];
    io[ud] = y*y+z*z<r2 ? -1 : 1;
  }
  if (io[0]*io[1]>0) return false;
  else {
    const TSetup &setup = TSetup::Ref();
    int ipli_z1 = setup.vIplFirst()[1];
    static int ipla_RICH = -1; if (ipla_RICH<0) {
      double xRICH = 750;
      for (int ipl = ipli_z1; ipl<setup.vIplLast()[1]; ipl++) {
	if (setup.vDetect(ipl).X(0)>xRICH) break; ipla_RICH = ipl;
      }
      if (ipla_RICH<0) CsErrLog::mes(elFatal,
        "Inconsitency: No detector plane upstream of RICH in zone 0x2");
    }
    int nSpacePts, nProjs;
    t.Evaluate(ipli_z1,ipla_RICH,0x1,nSpacePts,nProjs);
    if (nSpacePts<2 || nProjs<3) return false; 
  }
  return true;
}
