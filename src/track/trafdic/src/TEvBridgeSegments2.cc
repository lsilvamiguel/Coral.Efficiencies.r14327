// $Id: TEvBridgeSegments2.cc 14069 2015-09-17 20:44:46Z lsilva $

#include <iostream>
#include <string.h>
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "CsMCTrkHit.h"
#include "CsOpt.h"
#include "CsErrLog.h"
#include "TOpt.h"
#include "TSetup.h"
#include "TEv.h"
#include "Traffic.h"
#include "TWatches.h"
#include "TVtxMC.h"
#include "TLattice.h"
#include "TDisplay.h"
//#define BridgeS_DEBUG_RECYLEFIT
#ifdef BridgeS_DEBUG_RECYLEFIT
#  include "TMatrixDSym.h"
#endif

using namespace std;

/*! 
  \brief Traffic bridging methods
*/
/*!

  Bridge track segments through the magnets and muFilter 

*/

// - "TEv::BridgeSegments2": Bridging...
//   ... over the magnets SM1 and SM2,
//   ... (upon option) over the target,
//   ... over muFilterB.
// - TEv::BridgeOverTarget: dedicated to bridging over the target, called by
//  "BridgeSegments2"
// - TEv::DoubleBridge: simultaneous bridging over SM1 and SM2, called by
//  "BridgeSegments2".
// - TEv::UpdateKF: dedicated to the FullKF fit of track pairs, upon option
//  "ReMode[46]==0x1".

/*
  Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/TEvBridgeSegments.cc":
   i) Selection of pairs of segments to be tried for bridging made dependent
     upon, reconstructed, momentum, and including a, loose, cut on the
     left-right y coord.
  ii) TWatches to measure time spent fitting.)
 iii) Tag tracks passing the bridging selection and successfully fitted, so
     that they be not refitted a second time. QN fit, which is uncomplete, is
     given a specific tag, which is to be used to force refit.
  iv) Clean tracks list (so as to discard fakes produced by an enhanced
     "PrePattern").
   v) Backtrack w/ pickup for SM2-bridged, SM1->SM2 and SM2-muFilter TTrack's
  vi) Cleaning of THit's contributing worst to chi2 for scifi only track
     segments.
 vii) PixelGEM/MMs.
viii) etc..etc..
 */

// Interface
bool IsYokeTrk(const TTrack &tl, const TTrack &tr);
bool upDownCompatible(int imag, float pinv, float y0, float z0,
		      double ddip, double dz, double dy);
//#define NTUPLE_BRIDGE_CUTS
#ifdef NTUPLE_BRIDGE_CUTS
void fillNtupleBrigeCuts(TTrack &tl, TTrack &tr, float pimc,
			 int imag, float pinv, float y0, float z0,
			 double ddip, double dz, double dy, double dazi);
#endif

void TEv::BridgeSegments2()
{
  if (TOpt::ReMode[3]>0) return; // Bridging is OFF

  const TSetup& setup = TSetup::Ref();

  int print = TOpt::Print[4]; 

  static int hist = 0;   // Histogramming
  // Monitoring of the convergence of the momentum fit: # of iters and, in MC
  // case, reco'd-gen'd at various stages of the fit.
  static TH2D *hBSiter, *hBSerr;
  // Bridging over muWallA
  static TH2D *muW_dy, *muW_dz, *muW_dazi, *muW_ddip;

  // VSAT region => little redundancy => tracks there may turn out to be below
  // par, i.e. have a pattern of hits that does not constrain enough the chi2.
  // "vSATPar" records what #hits may be considered as signaling such cases (
  // which will then be examined further and possibly granted a malus).
  static int vSATPar[3], vSATType[3];
  static bool first = true; if (first) {
    first = false;

    //    ******************** INITIALISATIONS ********************

   if (TOpt::Hist[4]) {                                  // ***** HISTO BOOKING
      hist = 1;
      if (!gDirectory->cd("/Traffic")) {
	TDirectory *traffic = gDirectory->mkdir("/Traffic","TraFFiC");
	traffic->cd();
      }
    }
    if (hist) {
      hBSiter = new TH2D("hBSiter","BridgeS: #iters for imag#timesgood/bad#timesOK/!OK",
			 16,-.5,15.5,12,-.5,11.5);
      hBSerr = new TH2D("hBSerr","BridgeS: 1/P-1/PMC for SM1(</>40 GeV)/2#timesseed/dico/KF)",
			1000,-.1,.1,9,-.5,8.5);
    }
    if (hist) {
      muW_dy   = new TH2D("muW_dy",  "muW dy",  100,-15,15,5,0,.25);
      muW_dz   = new TH2D("muW_dz",  "muW dz",  100,-15,15,5,0,.25);
      muW_dazi = new TH2D("muW_dazi","muW dazi",100,-15,15,5,0,.25);
      muW_ddip = new TH2D("muW_ddip","muW ddip",100,-15,15,5,0,.25);
    }
    //              **********  VSAT: WHAT #HITS may be BELOW PAR (in zones 0x3)
    // - This is setup dependent => infra a few known setups are considered. 
    // - Note: Considering the sole #hits may not suffice (to tell when a
    //  pattern of hits is below par): here we only define the #hits (or #DFs)
    //  below which below par patterns may occur. 
    int iZ; for (int iZ = 0; iZ<3; iZ++) { vSATPar[iZ] = 0; vSATType[iZ] = 0; }
    int nZs0x7 = setup.vIplFirst().size()>2 ? 3 : setup.vIplFirst().size();
    int nFIs, nGPs, ipl;
    for (iZ = 0; iZ<nZs0x7; iZ++) {
      for (ipl = setup.vIplFirst()[iZ], nFIs=nGPs = 0;
	   ipl<=setup.vIplLast()[iZ]; ipl++) {
	const TPlane  &p = setup.vPlane(ipl);
	const TDetect &d = setup.vDetect(p.IDetRef);	    
	if      (d.IType==22) nFIs++;                // Scifi
	else if (d.IType==28 || d.IType==29) nGPs++; // PixelGEM
      }
      if      (iZ==0) {
	if (nFIs==6) // mu setup in [2002,2007]: Any pair of XYU space points
	  // may give good chi2. Define par as being at least that, viz. 6 hits
	  vSATPar[0] = 6;
	else // I haven't thought about it yet => Let par also be #hits = 6
	  vSATPar[0] = 6;
      }
      else if (iZ==1) {
	// #hits in FI05<2 is definitely below par, because no constraint in
	// either Y or Z dimension or, in the [2006,2007] small lever arm. This
	// can happen when total #hits is less than # of available detectors.
	// => Set "vSATPar[1]" accordingly (i.e. so as to check the condition on
	// #hits in FI05). Apart from this, total #hits < 5 does not allow for
	// any Y/Z correlation => also below par but somewhat less dramatic.
	if      (nFIs==7)            // mu setup in [2006,2007] (w/ FI55).
	  vSATPar[1] = 7;
	else if (nFIs==5)            // mu setup in [2002,2004]
	  vSATPar[1] = 5;
	else if (nFIs==2 && nGPs) {  // hadron setup in [2008,2009]
	  vSATPar[1] = 6;  // Par requires 2FIs+2pixels => #DFs>=6...
	  // PixGEMs provide for an automatic Y/Z correlation. => Remember it.
	  vSATType[1] = 1;
	}
      }
      else {
	if      (nFIs==2 && nGPs) {  // hadron setup in [2008,2009]
	  vSATPar[2] = 6; vSATType[2] = 1;
	}
      }
    }
  }

  list<TTrack>::iterator il,ir;

#define BridgeS_ScifiCLEANING
#ifdef BridgeS_ScifiCLEANING
  // Determine # of track segments in zones 0x1 and 0x2
  int nTrks_0x1 = 0, nTrks_0x2 = 0, nTrks_0x10 = 0;
  for (il = listTrack.begin(); il!=listTrack.end(); il++) {
    if      (il->Type==0x1)  nTrks_0x1++;
    else if (il->Type==0x2)  nTrks_0x2++;
    else if (il->Type==0x10) nTrks_0x10++;
  }
#endif

  TTrackPair2 tp;
  list<TTrackPair2> lTrackPair1, lTrackPair2, lTrackPair3;
  list<TTrackPair2>::iterator itp1,itp2;

  int targetField = 0, targetBridging = TOpt::ReMode[26];
#define BridgeS_OVER_TARGET
#ifdef BridgeS_OVER_TARGET
  int imag0 = (targetBridging&0x1b) ? 0 : 1;
  list<TTrackPair2> lTrackPair0; list<TTrackPair2>::iterator itp0;
#else
  const int imag0 = 1;
#endif
  for (int imag = 2; imag>=imag0; imag--) {

    // ********** LOOP OVER MAGNETS TRACKS ARE BRIDGED THROUGH **********

    TTrack tr , tl, tx; unsigned int unusedId = tx.Id;

    static double x0; double fldint(0); list<TTrackPair2> *lTP;
    if      (imag==1) {
      x0 = setup.MagCenter(1,0); fldint = setup.MagFieldInt[1];
      lTP = &lTrackPair1;
    }
    else if (imag==2) {
      x0 = setup.MagCenter(2,0); fldint = setup.MagFieldInt[2]; 
      lTP = &lTrackPair2;
    }
#ifdef BridgeS_OVER_TARGET
    else if (imag==0) { //   ***** SPECIAL CASE of BRIDGING OVER TARGET *****
      //  I.e. looking for continuations of beam telescope tracks into zone 0x1
      // by bridging over target.
      // - Special cases:
      //   - 0x08 = Catching downstream interaction.
      //    The bridging is used to solve cases where the 0x7-reco of the
      //    scattered particle displays an uncompatibility between LAS/SM1 and
      //    SAS/SM2.
      //   - 0x10 = Continuation beyond SM1.
      //    No actual bridging is performed. The association is instead saved
      //    into an attribute of the, 0x10, track ("TTrack::Associate", later to
      //    be ported to corresponding CsTrack). Having in view vertexing, and,
      //    w/in vertexing, the determination of the best vertex candidate: in
      //    that context, vertices built on beam tracks w/o associate will take
      //    precedence.
      //   - In both cases, we restrict ourself to tracks VSAT-(almost)only in
      //    zone 0x1, because a beam is likely to end up in VSAT while in 0x1.
      // - In any case, one has to try to preserve low Q^2 events, which can
      //  easily be confused w/ tracks not interacting in the target.
      if (targetBridging&0x8) {
	// 0x08 = Catching downstream interaction: we further restrict ourself
	// to events where there are few candidates, be it on the beam telescope
	// or the 0x1-zone side, in order to save CPU and also because the rate
	// of fake continuations is then kept low.
	int nVSATs0x1, nTs0x10; // # of VSAT-(almost)only 0x1 and of 0x10 tracks
	for (il = listTrack.begin(), nVSATs0x1=nTs0x10 = 0;
	     il!=listTrack.end(); il++) {
	  if (((*il).Type&0x1) && ((*il).Scifi&0x101)) nVSATs0x1++;
	  if ((*il).Type&0x10) nTs0x10++;
	}
	if (nTs0x10>4 ||         // Note: I put arbitrary cuts here...
	    nVSATs0x1>30)        // ...the 0x1 cut allows for combinatorics
	  targetBridging &= 0x13;
	if (!targetBridging) continue;
      }
      x0 = setup.MagCenter(0,0);
      // Target: nothing(0), solenoid(1) or dipole(2)?
      if (setup.MagScale[0]==0) // Scale == 0 => No field.
	targetField = 0;
      else
	// Solenoid case subdivides in 2 cases: OD(1) and SMC(2).
	// Dipole                             : OD(4) and SMC(3)
	targetField = setup.MagFlag1[0]>=3 ? 2 : 1;
      if (targetField<2) fldint = 0;
      else               fldint = setup.MagFieldInt[0]; 
      lTP = &lTrackPair0;
    }
#endif
    else lTP = 0; // In order to comply w/ compiler...

    float pimc(0), pinv(0), pinv_prev; double chi2, chi2_prev; int iter;

    int prvRId, backupRId; static double prvCop;
    TTrack prvRT; double Hnoise[15];
    for (ir = listTrack.begin(), prvRId=backupRId = -1; ir!=listTrack.end();
	 ir++) {

      // ***** LOOP OVER RIGHT TRACK PIECES *****

      switch (imag) {
      case 0: if (!(ir->Type&0x1))     continue;
#ifdef BridgeS_OVER_TARGET
	if ((targetBridging&0x18) && !(*ir).Scifi) // Over target, optional...
	  continue;              // ...consider only VSAT-(almost)only tracks
	if (targetBridging==0x10) {    // 0x10: Case continuation beyond SM1...
	  // ...Consider only SM1-bridged (in previous iteraction) tracks.
	  int match; for (itp1 = lTrackPair1.begin(), match = 0;
	       itp1!=lTrackPair1.end(); itp1++) {
	    TTrackPair2 &tp1 = *itp1; TTrack &t = *tp1.iL;
	    if (t.Id!=(*ir).Id) continue; if (tp1.Chi2>3) break;
	    match = 1; break; // ...and w/ pretty good chi2
	  }
	  if (!match) continue;
	}
#endif
	break;
      case 1: if ((ir->Type&0x3)!=0x2) continue; break;
      case 2: if ((ir->Type&0x6)!=0x4) continue; break;
      case 3: if (ir->Type!=0x8)       continue; break;
      }
      tr = *ir; // copy right track

      list<TTrackPair2>::iterator minITP; static double minChi2; 
      for (il = listTrack.begin(), minITP = lTP->end(); il!=listTrack.end();
	   il++) {

	// ***** LOOP OVER LEFT TRACK PIECES *****

	if (il==ir) continue; // the same
	switch (imag) {
	case 0: if (il->Type!=0x10)      continue; break;
	case 1: if (il->Type!=0x1)       continue; break;
	case 2: if ((il->Type&0x6)!=0x2) continue; break;
	case 3: if ((il->Type&0xc)!=0x4) continue; break;
	}
	tl = *il; // copy left track

	bool good_pair = false; // True for correct pair of well reco'd tracks
	if (tl.IKine>=0 && tr.IKine>=0) {
	  good_pair = tl.IKine==tr.IKine; int iv;
	  if (!good_pair && (iv = vecKine[tr.IKine].IVtxOrig)>=0 &&
	      vecVtxMC[iv].IOrig==tl.IKine) {
	    int idl = vecKine[tl.IKine].PDGid(), idr = vecKine[tr.IKine].PDGid();
	    if ((abs(idl)==211/* pi */ || abs(idl)==321/* K */) &&
		abs(idr)==13/* mu */ && idl*idr>0) good_pair = true;
	  }
	}
	if (good_pair) {
	  double pivtx = vecKine[tl.IKine].Pinv();
	  const TKine &k = vKine(tl.IKine);
	  // Momentum at target exit: i.e. after energy is lost in target In
	  // fact, one gets the momentum at the 1st entry the MC hit list. Just
	  // to make things simpler.
	  const THitMC &hmc = vHitMC()[k.vHitMCRef().front()];
	  const CsMCTrkHit *csh =
	    dynamic_cast<CsMCTrkHit*>(const_cast<CsMCHit*>(hmc.PtrHit()));
	  pimc = (pivtx>0?1:-1)/csh->getP().mag();
	}
	int jmag = imag==1 && fabs(pimc)>.025 ? 0 : imag; 
	if (TOpt::ReMode[6] && !good_pair) continue; // Ideal bridging

	// Rough track par. comparisson
	double dazi= tr.Hfirst.azi() - tl.Hlast.azi();
	double ddip= tr.Hfirst.dip() - tl.Hlast.dip();
	double dz = tr.Hfirst(2) + (x0-tr.Hfirst(0))*tan(tr.Hfirst.dip())
	  -        (tl.Hlast (2) + (x0-tl.Hlast (0))*tan(tl.Hlast.dip()));
	double dy = tr.Hfirst(1) + (x0-tr.Hfirst(0))*tr.Hfirst(3)
	  -        (tl.Hlast (1) + (x0-tl.Hlast (0))*tl.Hlast (3));
	double dt = tl.SigmaTime<0 || tr.SigmaTime<0 ? 0 : 
	  (tr.MeanTime-tl.MeanTime)/
	  sqrt(tr.SigmaTime*tr.SigmaTime+tl.SigmaTime*tl.SigmaTime);
	
	if (print>0) {
	  cout <<"\n BridgeSegments2 ==> "<<tl.Id<<"("<<tl.IKine<<")   "<<tr.Id<<"("<<tr.IKine<<")";
	  if(good_pair) cout<<"  Pmc = "<<1./pimc;
	  cout<<"  mag = "<<imag<<" FldInt = "<<fldint<<" x0 = "<<x0<<endl;
	  cout<<"Delta Azi (mrad) "<<dazi*1000.<<"   Delta Dip (mrad) = "<<ddip*1000.
	      <<"    Delta Z = "<<dz<<endl;
	}

	// y and z @ ``bridging'' center: y0, z0
	float y0 = tr.Hfirst(1) + (x0-tr.Hfirst(0))*tr.Hfirst(3);
	float z0 = tr.Hfirst(2) + (x0-tr.Hfirst(0))*tan(tr.Hfirst.dip());
	//Estimate momentum
	if (fldint) {
	  // Radius of curvature, rho, in a uniform field of length L, w/
	  // deviation theta, incidence alpha et exit angle beta
	  // rho * tg(theta/2) * (cos(alpha)+cos(beta)) = L
	  double thly = atan(tl.Hlast(3)), thry = atan(tr.Hfirst(3));
	  pinv = 2*cos((thly+thry)/2)*sin((thly-thry)/2)/.3E-3/fldint;
	  pinv *= cos(tl.Hlast.dip());
	  if (hist && good_pair) hBSerr->Fill(pinv-pimc,3*jmag);
	}
	else if (imag==0) {
	  if (targetField==1)
	    pinv = TOpt::iCut[15]/TOpt::dCut[4]; // Case solenoid
	  else                                        // Case no field...
	    // ...might be useful if "lowRedundancy", and hence MS,...
	    pinv = TOpt::iCut[15]/TOpt::dCut[4]; // ...is ever enabled for that case
	}

#ifdef DEBUG
	static int idebug = 0;
	if (idebug>1) { tl.DumpHits(0); tr.DumpHits(0); }
#endif

	//   ********************************************************
	//   ************************* CUTS *************************
	//   ********************************************************
#ifdef NTUPLE_BRIDGE_CUTS
	if (good_pair)
	  fillNtupleBrigeCuts(tl,tr,pimc,imag,pinv,y0,z0,ddip,dz,dy,dazi);
#endif
	if (TOpt::ReMode[6]==0) { // If not ideal bridging
	  if (TOpt::dCut[15]!=0 &&            // ***** TIME DIFF CUT... *****
	      (tr.Scifi&0xff) && ((tl.Scifi&0xff) // ...If scifi-only segments
				  || imag==0)) {  // ...Or briding over target
	    // Let's also require reasonable time dispersion, in order to
	    // avoid cases where wrongly associated hit biases track's time. 
	    if (tl.DispTime<4*tl.SigmaTime &&
		tr.DispTime<4*tr.SigmaTime &&
		fabs(dt)>TOpt::dCut[15]) continue; 
	  }
	  if (!upDownCompatible(imag,pinv,y0,z0,ddip,dz,dy)) continue;
	}

	if (print>0) cout<<"P estim. =  "<<1./pinv<<endl;

	if (TOpt::ReMode[6]==0)	  // ********** CUT on MOMENTUM **********
	  if (fabs(1./pinv)<0.1 || 1000<fabs(1/pinv)) continue;

	//       ********** BUILD MERGED TRACK **********
	if (TOpt::ReMode[46] &&         // FullKF fit requested...
	    !(TOpt::ReMode[46]&0x2) &&  // ...but not systematically...
	    TOpt::ReMode[47] && // ...and recycling fit from previous "(tl,tr)"
	    (int)tr.Id==prvRId && backupRId!=prvRId) {
	  // If "tr" of current "(tl,tr)" pair same as previous, let's prepare
	  // for recycling the FullKF fit of the "tr" segment.
	  prvRT = tr;
	  prvRT.Id = unusedId;// So that THit::setTrackID's don't change when the prvRT object is deleted.
	  backupRId = prvRId; // So as to remember current "tr" is backed up.
	  prvRT.Hlast(5) = prvCop; prvRT.Hlast(5,5) = 1.E-20;
	  prvRT.QuickKF(-1,1);
	  const THlx &H = prvRT.Hfirst, &H_MS = tx.mHub[tr.lPlnRef.front()];
	  int i,j,k; for (i = 1, k = 0; i<=5; i++) for (j = 1; j<=i; j++, k++) {
	      Hnoise[k] = H_MS(i,j)-H(i,j);
	    }
	}
	tx = tl; tx.Append(tr,"tmp"); tx.Chi2tot = 1e6;
#ifdef DEBUG
	if (idebug) tx.DumpHits(idebug);
#endif
	if (imag) {
	  if   (good_pair) Traffic::Ref().Stopwatch.Start(imag+8);
	  else             Traffic::Ref().Stopwatch.Start(imag+10);
	}

	//       ********** INITIAL VALUES **********
	TTrackPair2 *minTP = minITP!=lTP->end() ? &*minITP : 0;
	iter = 0; // Iteration counter (count KF fits only)
	chi2 = 1.E10; tx.Hfirst(5) = pinv;
	if (imag==0 && targetField!=2) tx.Hfirst(5,5) = 1.E-20;
	else                           tx.Hfirst(5,5) = 1.E-4;
	int fastFit = 1;                         // Default use QuickKF
	tp.WorstPl = -1;
#ifdef BridgeS_ScifiCLEANING
	// ***** LOW REDUNDANCY (mainly scifi) TTrack's *****
	// (Restrict to case of zones 0x1 segments, because, otherwise, too
	// few hits in play to do anything valuable. This means one accepts
	// that case when, e.g., one of the 4 or 5 hits of a 0x2 scifi
	// segment being fake has the reco of the track fail (still one
	// has the possibility to remove that fake hit by cleaning the full
	// TTrack in later stage, provided that the reco makes it there).
	bool lowRedundancy = TOpt::dCut[11] &&
	  (tl.NHits==5 && imag==1 && nTrks_0x1 *nTrks_0x2<=50 ||
	   tr.NHits==5 && imag==0 && targetField &&
	   nTrks_0x10*nTrks_0x1<=50);
#else
	bool lowRedundancy = false;
#endif
	if (TOpt::ReMode[6]>1) { // ***** SUPER IDEAL: INIT VALUE(S) <- MC *****
	  // - Based on ReMode[6], several actions are foreseen:
	  //   =2: Take momentum P from associated TKine object (which cannot
	  //      but exist in case ReMode[6]>0, cf. supra), i.e. the MC truth
	  //      @ the vertex of origin.
	  //   =3: Take also the rest of the helix parameters: this is
	  //      ineffective, since the other track parameters are usually
	  //      pretty accurately computed from the straight line fit of "tl",
	  //      and anyway are not used in what follows: "Hlast" of "tr" is
	  //      used instead.
	  //   =4: Take P from the MC helix at the last measured detector.
	  tx.Hfirst(5) = pimc;
	  if (TOpt::ReMode[6]>2) {   // Hyper ideal: all of helix is init'd
	    TKine &k = vecKine[tr.IKine];
	    int ok = 0; if (TOpt::ReMode[6]>3) { // ...init'd from CsMCTrkHit
	      const vector<int> &ihs = k.vHitMCRef();
	      int idetf = setup.vDetect(tl.lPlnRef.back()).IDet;
	      for (int ih = (int)ihs.size()-1; ih>=0; ih--) {
		const THitMC &h = vecHitMC[ihs[ih]];
		if (h.DetRef().IDet==idetf && h.IOrig==0) {
		  CsMCHit *ch = const_cast<CsMCHit*>(h.PtrHit());
		  CsMCTrkHit *th = dynamic_cast<CsMCTrkHit*>(ch);
		  if (th) {
		    const CLHEP::Hep3Vector P = th->getP(); double p = P.mag();
		    if (p) {
		      ok = 1;
		      tx.Hfirst(5) = k.Q()/p;// Later to be assigned to "Hlast".
		      //tx.Hfirst(3) = P[0]/P[2]; tx.Hfirst(4) = P[1]/P[2];
		      //tx.Hfirst(1) = th->getX()/10; tx.Hfirst(2) = th->getY()/10;
		      goto success;
		    }
		  }
		  break;
		}
	      }
	    }
	    if (!ok) {
	      float px = (float)k.P(0), py = (float)k.P(1), pz = (float)k.P(2);
	      tx.Hfirst(3) = py/px; tx.Hfirst(4) = pz/px;
	      const TVtxMC &vtx = vVtxMC(k.IVtxOrig); tx.Hfirst(0) = vtx.V(0);
	      tx.Hfirst(1) = vtx.V(1); tx.Hfirst(2) = vtx.V(2);
	    }
	  }
	}
	if (imag==0 && targetField<2) {
	  // ***** BRIDGING OVER TARGET: CASE NO BENDING MAGNET *****
	  if (!tx.QuickKF(-1,0)) continue;
	  // We are after high momenta => no need to make for MS by normalising
	  // chi2 by more than (#hits-4), as is done in standard case.
	  if ((chi2 = tx.Chi2tot/(tx.NDFs-4))<7) {
	    tx.Hlast(5) = TOpt::iCut[15]/TOpt::dCut[4]; 
	    for (int iter = -1; iter<=1; iter += 2) {
	      if (iter==-1) tx.Hlast  *= 1.E4; // Scale before fit, w/ some...
	      else          tx.Hfirst *= 1.E2; // ...unchecked cooking recipe
	      tx.Hlast(5,5) = 1.E-20; // Fix momentum
	      if (!tx.QuickKF(iter,1)) break;
	      if ((chi2 = tx.Chi2tot/(tx.NDFs-4))<5.5) {
		// A built-in cut for the time being...
		tp.IFit = 0x2; goto success;
	      }
	    }
	  }
	  goto next_pair;
	}

	// Option "ReMode[46]" for bypassing fast fit (in order not to miss
	// tracks undergoing high MS in RICH pipe): 2 versions:
	//  I) If &= 0x2, systematically switches to "FullKF".
	// II) Otherwise, switches to "FullKF" only if track goes through RICH
	//   pipe (and while bridging over SM1). A "FullKF" is then performed a
	//   posteriori on all successfully bridged pairs.
	if (TOpt::ReMode[46]) {
	  if (TOpt::ReMode[46]&0x2)              // FullKF systematic
	    fastFit = 0;
	  else if (imag==1 && tr.CrossRPipe())   // FullKF if crossing RICH pipe
	    fastFit = 0;
	  if (!fastFit) tx.InsertMissedPlanes();
	}

	tp.IFit = 0;
	if ((TOpt::ReMode[14]&0x8) &&                // Dico fit optional? 
	    imag==1 &&        // Dico fit only for SM1
	    !lowRedundancy) { // Exclude low redundancy TTrack's

	  //          ***** DICO FIT ATTEMPT*****

	  if (tx.QNewtonFit(1,3)) { // I.e. not QN misfit (outside dico,etc..)
	    if (fastFit) {
	      if (tx.Chi2tot==1e6) goto next_pair;
	      chi2 = tx.Chi2tot/(tx.NDics-6); 
	      tp.IFit = 0x8; pinv = fabs(tx.Haux(5));
	      if (hist && good_pair) hBSerr->Fill(pinv-pimc,3*jmag+1);
#ifdef DEBUG
	      if (idebug) tx.DumpHits(1);
#endif
	      goto success;
	    }
	    else if (tx.Chi2tot<1.e5 /* I.e. !=1e6 */) {
	      // - Let's seed FullKF w/ the dico momentum. In case the track
	      // undergoes high MS (in, e.g., the heavy RICH pipe) QN fit won't
	      // converge and we fall back to the initial guess for the seed.
	      // - Since we're about to pursue w/ the iterative KF momentum fit
	      // let's not update "chi2": QN's chi2 cannot be compared w/ KF's.
	      double chi2QN = tx.Chi2tot/(tx.NDics-6); 		
	      if (chi2QN<TOpt::dCut[7]) tx.Hfirst(5) = tx.Haux(5);
	    }
	  }
	}

	//        ********** ITERATIVE KF MOMENTUM FIT **********

	while (iter < 15) {                      // ***** ITER LOOP
	  pinv_prev = tx.Hfirst(5); chi2_prev = chi2;
	  tx.Hlast(5,5) = tx.Hfirst(5,5);

	  bool ret; if (fastFit) {               // ***** FIT CANDIDATE *****
	    tx.Hlast(5) = tx.Hfirst(5);
	    ret = tx.QuickKF(-1,1);
	  }
	  else {
	    double r = tx.Hfirst(5);
	    if (iter) // Energy loss correction
	      // P @ last helix = P @ first helix (found in previous iter) - ELoss
	      tx.Hlast(5) = r/(1-tx.ELoss()*fabs(r));
	    else
	      tx.Hlast(5) = r;
	    ret = tx.FullKF(-1); 
	  }
	  if (!ret) goto next_pair;

	  chi2 = tx.Chi2tot/tx.NDFs;
#ifdef BridgeS_ScifiCLEANING
	  if (lowRedundancy &&           // ***** LOW REDUNDANCY: CLEANING *****
	      iter==0 &&                               // UPON 1ST iter
	      chi2>TOpt::dCut[11] &&                   // NOT A VERY GOOD FIT...
	      // ...But still reasonable: one does not want this piece
	      // of code to be resorted to too often (because it consumes
	      // CPU (although I don't know how much)) and its reliability
	      // hasen't been fully evaluated).
	      chi2<TOpt::dCut[7]/3) {
	    // => There might well be A GHOST HIT PICKED UP?: Search
	    //   for it by looking for large chi2 increment. Which implies
	    //   full KF, for only full KF prevents MS to blur the picture
	    //   (and, in any case, only it does provide chi2 increments)
	    if (!fastFit || tx.FullKF(-1)) {
	      if (fastFit && // Meaning FullKF has just been performed on "tx".
		  // And this, w/o inserted missed planes which disqualifies...
		  (int)tr.Id==prvRId)// "tx" as a source of already fitted "tr".
		prvRId = -1;
	      int worstPl = -1, worstGr = -1; double worstChi2Incr;
	      tx.WorstHit(tx.mChi2,worstPl,worstGr,worstChi2Incr);
	      if (worstChi2Incr>TOpt::dCut[12] &&      // WORST IS BADLY WORST?
		  // Is this taking place in the SCIFI ONLY TRACK SEGMENT?
		  worstGr==0) {
		//                                      => ERASE IT...
		// (Don't update THit->TTrack reference since it hasn't been set
		// by "TTrack::Append".)
		int jhit = tx.CancelHit(worstPl,false);
		tx.Hfirst(5) = pinv;                            // ...INIT AGAIN
		if (imag==0 && targetField!=2) tx.Hfirst(5,5) = 1.E-20;
		else                           tx.Hfirst(5,5) = 1.E-4;
		if (fastFit==0) ret = tx.FullKF(-1);            // ...REFIT
		else            ret = tx.QuickKF(-1,1);
		if (!ret) goto next_pair;
		chi2 = tx.Chi2tot/tx.NDFs;
		if (chi2>TOpt::dCut[11]) {                   // ...UNSUCCESSFUL?
		  tx.RestoreHit(worstPl,jhit,false);         // => RESTORE...
		  //tp.WorstPl = -1; // Useless: already set =1, upon iter==0
		  tx.Hfirst(5) = pinv;
		  if (imag==0 && targetField!=2) tx.Hfirst(5,5) = 1.E-20;
		  else                           tx.Hfirst(5,5) = 1.E-4;
		  if (fastFit==0) ret = tx.FullKF(-1);       // ...REFIT
		  else            ret = tx.QuickKF(-1,1);
		  if (!ret) goto next_pair;
		  chi2 = tx.Chi2tot/tx.NDFs;
		}
	      }
	      else tp.WorstPl = worstPl;
	    }
	  }  // End of scifi cleaning
#endif

	  pinv = fabs(tx.Hfirst(5));

	  if (TOpt::ReMode[6]==0) {
#warning TODO: Lower cut on momentum
	    if (pinv>10/3.) {                 // ***** LOW CUT on MOMENTUM *****
	      if (print>1) cout<<"Too small momentum: "<<1./tx.Hfirst(5)<<endl;
	      goto next_pair;
	    }
	  }
	  double chi2Cut;                  // ***** CHI2 CUT... *****
	  if (fastFit) {
	    if      (pinv<.2) chi2Cut = TOpt::dCut[7];
	    else if (pinv<.8) chi2Cut = TOpt::dCut[7]*1.5;   // Looser for low P
	    else if (iter>2)  chi2Cut = TOpt::dCut[7]*1.5;
	    //                       @ very low P, looser for earlier iterations
	    else if (iter>1)  chi2Cut = TOpt::dCut[7]*2;
	    else if (iter)    chi2Cut = TOpt::dCut[7]*3;
	    else              chi2Cut = TOpt::dCut[7]*20;
	  }
	  else {
	    if (pinv<.2 || iter) chi2Cut = TOpt::dCut[97];
	    else {                                        // Upon 1st iter...
	      if (pinv<.8) chi2Cut = TOpt::dCut[97]*1.5;  // ...looser for low P
	      else         chi2Cut = TOpt::dCut[97]*2;
	    }
	  }
	  if (chi2>chi2Cut) {
	    if (print>1)
	      printf("Too bad chi2/NDFs %.2f, P = %.2f\n",chi2,1/pinv);
	    goto next_pair;
	  }

	  if (print >0)
	    cout<<"P = "<<1./tx.Hfirst(5)<<"\t  Err = "<<sqrt(tx.Hfirst(5,5))
		<<"\t  Pinv-PinvMC = "<<tx.Hfirst(5)-pimc<<endl;

	  //                                        ***** CONVERGENCE TEST *****
	  int ichi2      = int(chi2     *10); // Keep only 1 decimal digit
	  int ichi2_prev = int(chi2_prev*10);
	  if (ichi2==ichi2_prev &&
	      chi2<TOpt::dCut[7]*1.5 /* Skip early, loose, iter's */) {
	    if (print>1) cout<<"---> Converged!\n"; break;
	  }
	  
	  iter++;
	}; // <---------- end of iteration loop 

	if (iter==15) {
	  printf("TEv::BridgeSegments2 ==> Max #iterations(=%d) reached:"
		 "mag #%d, %s pair, 1/P: %.4f->%.4f, chi2: %.3f->%.3f\n",
		 iter,imag,good_pair?"good":"wrong",pinv_prev,tx.Hfirst(5),
		 chi2_prev,chi2);
	  goto next_pair;
	}

      success:                                             // ***** CONVERGENCE!
	tp.iL = il; tp.iR = ir;
	if (fastFit && TOpt::ReMode[46] && imag!=0) {
	  //      ***** CASE FULL KF FIT REQUIRED and NOT DONE YET *****
	  //  Since this FullKF consumes CPU, particularly so when in
	  // conjunction w/ ROOTGeometry, let's try to limit its use.
	  // - ReMode[47]&0x1: Do not automatically refit full tracks, but, in a
	  //  sequence of "tl,tr" pairs w/ same "tr", refit only the "tl" piece.
	  // - ReMode[47]&0x2: Do not systematically refit after cleaning.
	  // - In any case, let's skip this block while BridgeOverTarget: Full
	  //  KF is useful to put on an equal footing tracks going through
	  //  vastely different thicknesses of material, while tracks bound to
	  //  undergo BoT, and *interesting*, all follow approx. the same path.
	  if (!UpdateKF(imag,targetField,tx,tp,
			prvRId,prvRT,prvCop,Hnoise,
			minTP,minChi2,
			chi2))
	    goto next_pair;
	}
        else {// else, in all other cases, assign also "TTrackPair2::Hfirst" 
	  if (tp.IFit==0x8) {// QN fit: fitted helix is in "TTrack::Haux"
	    tp.Hfirst = tx.Haux; tp.NDics = tx.NDics;
	  }
	  else tp.Hfirst = tx.Hfirst;
	  chi2 = tx.Chi2tot/(tx.NDFs-5);
	}
	if (!TOpt::ReMode[6] /* !ideal */) {   // ***** FINAL CHI2 CUT (dCut[9])
	  if (chi2>TOpt::dCut[9]*
	      (pinv>.2 ? (pinv>1. ? 3 : 2) : 1)) // Looser for low momenta
	    goto next_pair;
	}
	if (!tp.IFit /* i.e. !QN fitted nor FullKF upon fastFit supra */)
	  tp.IFit = fastFit ? 0x2 : 0x20;

	//       ****************************************
	//       ***** STORE TRACK PAIR CANDIDATES *****
	//       ****************************************

	// Since the effect of multiscattering is at this point not yet taken
	// into account, chi2/NDF is disadvantageous to long tracks
	//          ***** => GIVE THESE A BONUS *****
#warning TODO: Bonus based on #DFs instead of #hits
	tp.Chi2 = !fastFit ? chi2 :
	  tx.Chi2tot/(tx.NDFs-6+tx.NHits*tx.NHits*.01);
	tp.Chi2tot = tx.Chi2tot;
	{ //   ***** VSAT TRACKS BELOW PAR: GIVE THEM A MALUS *****
	  // - VSAT region has little redundancy. Chi2 may be a bad estimator
	  //  of track quality. Determine if this is the case, give then malus.
	  // - So far as scifis are involved, would be good to consider timing,
	  //  although in principle this is already done in "TAlgo::FindSpace".
	  int nDFsl = tl.NDFs, nDFsr = tr.NDFs;
	  if (tp.WorstPl>=0) {
	    // P-pixel planes contribute for 2 units of DF.
	    int cor = (setup.vPlane(tp.WorstPl).IFlag&0x30) ? 2 : 1;
	    if (tp.WorstPl<setup.vIplFirst()[imag]) nDFsl -= cor;
	    else                                    nDFsr -= cor;
	  }
	  int belowParT; if (imag==1 && nDFsl<vSATPar[0]) {
	    belowParT = 1; if (nDFsl<vSATPar[0]-1) belowParT++;
	  }
	  else if (imag==2 && nDFsr<vSATPar[2]) belowParT = 1;
	  else /* Case of zone 0x2: init w/ 0 and update infra */ belowParT = 0;
	  if (imag==1 && nDFsr<vSATPar[1] || imag==2 && nDFsl<vSATPar[1]) {
	    // #hits in FI05<2 is definitely below par, because no constraint in
	    // either Y or Z dimension or, in the [2006,2007] small lever arm.
	    TTrack *t = imag==1 ? &tr : &tl;
	    list<int>::iterator ipr = t->lPlnRef.begin(); ipr++;
	    const TPlane &p = setup.vPlane(*ipr);
	    const TDetect &d = setup.vDetect(p.IDetRef);
	    if (d.X(0)>750 /* i.e. RICH abscissa */) belowParT++;
	    if (vSATType[1])// Hadron 2008/09 = "vSATPar[1]")...
	      // ...6 slots available and required. (So that a missing FI05
	      // yields "belowParT=2", one missing pixel hit yields "1", even if
	      // compensated by one non VSAT extra hit.)
	      belowParT++;
	    else // Else require only 5. Meaning 5<=#hits<7 in [2006,2007], w/ 2
	      // FI05 hits, yields "belowParT=0".
	      if (imag==1 && nDFsr<5 || imag==2 && nDFsl<5) belowParT++;
	  }
	  if (belowParT) {                // ***** GIVE MALUS to BELOW PAR TRACK
	    // Awkward wise for the time being, so as to keep w/ earlier version
	    // in the muon case.
	    int incr = vSATType[1] ? 1 : 2;
	    if (belowParT>1) {
	      incr += 1; if (belowParT>2) incr += 1;
	    }
	    tp.Chi2tot += incr*tx.NDFs*1.; tp.Chi2 += incr;
	  }
	}

	lTP->push_back(tp);
	if (minTP==&tp) minITP = lTP->rbegin().base();

      next_pair:
	if (hist) {
	  if (tp.IFit) {
	    if (tp.IFit==0x8) hBSiter->Fill(0.,    4*imag+(good_pair?0:2));
	    else              hBSiter->Fill(iter+1,4*imag+(good_pair?0:2));
	    if (good_pair && (tp.IFit&0x22)) hBSerr->Fill(pinv-pimc,3*jmag+2);
	  }
	  else
	    hBSiter->Fill(iter+1,4*imag+(good_pair?1:3));
	}
	if (imag) {
	  if   (good_pair) Traffic::Ref().Stopwatch.Stop(imag+8);
	  else             Traffic::Ref().Stopwatch.Stop(imag+10);
	}
      }  // End of loop over right track pieces
    }  // End of loop over left  track pieces
    tr.Id = ++TTrack::TrackCounter;// Make IDs of temporary tracks unique
    tl.Id = ++TTrack::TrackCounter;//   before destructor is called
    tx.Id = ++TTrack::TrackCounter;// (to prevent removal of it's IDs from sTrackIDs of TKine)

    lTP->sort(); //  *************** SORT MERGED PAIRS BY Chi2 ***************

  } // End of loop on "imag"



  //      ******************** BUILD GLOBAL TRACKS ********************

  //     ***** CHI2 CUTS *****
  // (Note that the default values rather correspond to the "QuickKF" case.)
  double goodChi2 =        TOpt::dCut[91] ? TOpt::dCut[91] : 3;
  double badChi2 =         TOpt::dCut[93] ? TOpt::dCut[93] : 20;
  double badVSATChi2 =     TOpt::dCut[92] ? TOpt::dCut[92] : 6;
  double chi20x7IncrCut =  TOpt::dCut[95] ? TOpt::dCut[95] : 3;
  double chi20x7TightCut = TOpt::dCut[96] ? TOpt::dCut[96] : 2;

  double xRICH = 750;
  list<TTrack>::iterator it; 
  //#define BridgeS_DEBUG_STEAL
  //  Patch to debug the stealing of the 0x1 segment of a high momentum track by
  // a lower momentum one that takes advantage of the high uncertainty
  // multiscattering imparts to its trajectory to fit any available 0x1 segment.
  //  This stealing process is particularly worrisome in the Primakoff case,
  // where e+/e- tracks can easily steal from the scattered pion or muon.
#ifdef BridgeS_DEBUG_STEAL
  TTrackPair2 *tpSteal = 0; int ikSteal = 0;
#endif

  if (TOpt::dCut[94]) {
    for (itp2 = lTrackPair2.begin(); itp2!=lTrackPair2.end(); itp2++) {
      TTrackPair2 &tp2 = *itp2; TTrack &tl = *tp2.iL, &tr = *tp2.iR;

      //   *************** FIRST LOOP OVER SM2 TRACK PAIRS ***************
      // - This first loop embeds a loop on SM1 pairs in search for prepending
      //  a 0x1 piece to 0x6 bridged pairs.
      // - It is meant to give reliable 0x7 tracks the edge over less reliable
      //  0x3 ones, in the competition to secure a 0x1 piece.
      // - This is particularly useful when the "FullKF" is used to sort out
      //  candidate track pairs: then a low momentum particle, tracked in 0x2,
      //  can easily steal a 0x1 piece from a higher momentum one, because of
      //  the versatility given to it by the uncertainty due to multiscattering.
      // - Reliability is evaluated on chi2, momentum, abscissa of starting
      //  point xFirst.
      //   - The cut on chi2 is recylcing the "dCut[91]" option used to assess
      //    0x6 track.
      //   - The, lower, cut on momentum makes so that the extraplation of the
      //    0x6 through the SM1 field is precise: dCut[94].
      //   - Cut on xFirst: same rationale as for momentum, while also rejecting
      //    V0 decays, far downstream ones, that is; xFirst<xRICH.
      // - If an extension is found that fails to pass all cuts (on P and T
      //  consistency and chi2), but does so only marginally, give up the whole
      //  bridging, including the 0x6 step, so that the 0x7 gets a second chance
      //  to be built, w/ this time less strict cuts, but also while competing
      //  w/ everybody, including 0x3-only tracks.

      //                                      ********** STANDARD PAIR SELECTION
      if (tl.Type!=0x2 || tl.IMark!=0)  // Already merged or in a better pair
	continue;
      if (tr.IMark!=0) continue;        // Already merged or in a better pair
      if ((tl.Scifi&0x202) && (tr.Scifi&0x404)) {// Hint at precise timing?...
	//                                  *****  CHECK TIMING COMPATIBILITY...
	// "tl/tr", being "Scifi" flagged, do possess time info: no need to
	// check their "SigmaTime" is >0 and hence indeed meaningful.
	double sl = tl.SigmaTime, sr = tr.SigmaTime, dT = sqrt(sl*sl+sr*sr);
	if (dT<.5 &&// Check that we have indeed good timing (cut corresponds...
	    // ...to 2 scifis at least in both "tl" and "tr").
	    fabs(tr.MeanTime-tl.MeanTime)>3*dT)// ...TIGHT TIME CONSISTENCY CUT
	  continue;// (Note: In 2nd loop, there's a rescue (of the 0x2) attempt.
      }
      //                                                ********** TAG BEST PAIR
      // (Note: there is no TTrack::Id ==0, cf. "PrePattern2".)
      tl.IMark = tr.Id; tr.IMark = tl.Id;
      //                                          ********** RELIABILITY CUTS...
      if (IsYokeTrk(tl,tr)) continue;               // ...Bypass yoke tracks
      if (tl.Hfirst(0)>xRICH) continue;             // ...xFirst cut
      double chi20x6 = tp2.Chi2, p0x6 = fabs(1/tp2.Hfirst(5));
      if (chi20x6>goodChi2) continue;               // ...Chi2 cut
      if (p0x6<TOpt::dCut[94]) continue;            // ...Momentum cut
      double cop2 = tp2.Hfirst(5), d2cop2 = tp2.Hfirst(5,5);

      int extend0x1; // =0: No extension, =2: Extension, =1: Ambiguity
      for (itp1 = lTrackPair1.begin(), extend0x1 = 0; itp1!=lTrackPair1.end();
	   itp1++) {
	TTrackPair2 &tp1 = *itp1; TTrack &tl1 = *tp1.iL, &tr1 = *tp1.iR;

	// ********** LOOK FOR A 0x1 (BEFORE SM1) TO PREPEND **********

	if (tp1.iR!=tp2.iL) continue;   // Must have common 0x2-track
	if (tl1.Type!=0x1) continue;    // Already prepended
	//           ***** FOUND A CONTINUATION *****
	//                           ***** DO NOT SPOIL GOOD CHI2 0x6 w/ BAD 0x3
	// - Although the current FIRST loop over SM2 pairs is precisely meant
	//  to favour 0x7 combinations, let's still be cautious, since the
	//  extension of 0x6 into 0x1 can later be achieved by "BackTrackSAS".
	// - As to decide between 0x3 and 0x6: currently considered 0x6 have
	//  already been carefully selected, let's validate it whether
	//  incompatible, very good, 0x3 combinations are available or not.
	double chi20x3 = tp1.Chi2; if (chi20x3>2*badVSATChi2) break;
	//                                              ***** DISAGREEMENT IN P?
	double cop1 = tp1.Hfirst(5), d2cop1 = tp1.IFit!=0x8 ? tp1.Hfirst(5,5) :
	  4e-4; // QN fit has no dcop: let's assume 2%
	double cop21 = fabs(cop2-cop1), dcop = sqrt(d2cop1+d2cop2);
	int pMismatch, tMismatch;
	if (cop21<3*dcop)        // A 3-sigma, somewhat strict, cut.
	  // In the current FIRST loop over SM2 pairs, where 0x6 tracks are
	  // exempted from competing w/ 0x2 ones, we must be strict.
	  pMismatch = 0;
	else if (cop21<6*dcop || // A 6-sigma, intermediate, cut, or...
		 // ...some more margin to make for possible tracking errors, w/
		 // chi2 as a telltale, that might bias the 0x3 momentum.
		 cop21<9*dcop && chi20x3>badVSATChi2/2)
	  pMismatch = 1;
	else pMismatch = 2;
	if (pMismatch>=2) break;
	tl.UseHitTime(); double dT0x6 = tl.SigmaTime>0, dT0x1 = tl1.SigmaTime;
	tMismatch = 0; if (dT0x6>0 && dT0x1>0) {
	  double qT = fabs(tl.MeanTime-tl1.MeanTime)/
	    sqrt(dT0x6*dT0x6+dT0x1*dT0x1);
	  if      (qT>5) tMismatch = 2;
	  else if (qT>3) tMismatch = 1;
	}
	if (tMismatch>=2) break;

	//                                        ********** ATTEMPT 0x7 MERGING
	TTrack t0x7(tl1); t0x7.Append(tl); t0x7.Append(tr); 
	extend0x1 = 1;
	if (tp1.WorstPl>=0) t0x7.SubHit(tp1.WorstPl);
	if (tp2.WorstPl>=0 && tp2.WorstPl!=tp1.WorstPl) t0x7.SubHit(tp2.WorstPl);
	t0x7.Hlast(5) = cop2; t0x7.Hlast(5,5) = d2cop2;
	//                                                         ***** FIT 0x7
	int chi2OK; if (TOpt::ReMode[46]) {
	  t0x7.InsertMissedPlanes(); chi2OK = t0x7.FullKF(-1) ? 2 : 0;
	}
	else                         chi2OK = t0x7.QuickKF(-1,1) ? 2 : 0;
	double chi20x7; if (chi2OK) {                          // ***** CHI2 CUT
	  chi20x7 = t0x7.Chi2tot/(t0x7.NDFs-5);
	  if      (chi20x7>goodChi2+2*chi20x7TightCut) chi2OK = 0;
	  else if (chi20x7>goodChi2+  chi20x7TightCut) chi2OK = 1;
	}
	if (chi2OK>=2) {                             // ***** CHI2 INCREMENT CUT
	  // We have here Reasonably low total chi2,..
	  // Let's check the chi2 increment when going from 0x3 or 0x6 to 0x7.
	  double chi2_34 = (chi20x3+tr.Chi2tot/(tr.NDFs-4))/2;
	  double chi2_16 = (tl1.Chi2tot/(tl1.NDFs-4)+chi20x6)/2;
	  if (chi20x7>chi2_34+chi20x7TightCut &&
	      chi20x7>chi2_16+chi20x7TightCut) chi2OK = 1;
	  // Also check the quality of 0x3: a bad tracking there may be what
	  // lure the 0x6 into matching it.
	  if (chi20x3>badVSATChi2) chi2OK = 1;
	}
	if (chi2OK>=2 && !pMismatch && !tMismatch) {      // ***** PUSH_BACK 0x7
	  // Note Adding, removing elements in "listTrack", because it's a list,
	  // does not invalidate any of the TTrackPair2 iterator members.
	  t0x7.IFit = TOpt::ReMode[46] ? 0x20 : 0x2;
	  t0x7.Hlast(5)   = t0x7.Hfirst(5); // Store momentum @ last point Helix
	  t0x7.Hlast(5,5) = t0x7.Hfirst(5,5);
	  t0x7.Haux = tp1.Hfirst;           // Save 0x3 momentum
	  listTrack.push_back(t0x7);
	  // Have to assign the new instatiated 0x7 track the "TTrack::Id" of
	  // its parent 0x1 track, in order to keep track of its taking part to
	  // "tp0", i.e. BoT, pairs.
	  TTrack &t = listTrack.back(); t.Id = tl1.Id;
	  tl1.IMark=tl.IMark=tr.IMark = -1;   // ***** FLAG 0x1,2,4 FOR DELETION
#ifdef BridgeS_DEBUG_STEAL
	  if (tl1.IKine==1 || tl1.IKsame==1 && tl1.NHsame>.75*tl1.NHits &&
	      tpSteal==0) {
	    if (tl.IKine!=1) { tpSteal = &tp1; ikSteal = tl.IKine; }
	  }
#endif
#ifdef BridgeS_OVER_TARGET
	  if (targetBridging&0x11) {
	    int nDFs0x1 = tl1.NDFs; double chi20x1 = tl1.Chi2tot/(nDFs0x1-4);
	    BridgeOverTarget(lTrackPair0,t,nDFs0x1,chi20x1);
	  }
#endif
	  extend0x1 = 2;
	}
	else {
	  tl.IMark=tr.IMark = 0;// ***** NO EXTENSION | AMBIGUITY: RESET "IMark"
	  if (!chi2OK) extend0x1 = 0;
	}
	break;
      }
      if (extend0x1==0) {   // ********** NO 0X1 EXTENSION => ACTUAL 0x6 MERGING
	tl.IFit = tp2.IFit;
	tl.Chi2tot = tp2.Chi2tot; // (Note: this is probably useless...)
	tl.Hfirst = tp2.Hfirst;
	tl.Append(tr);  // Append 0x4 to 0x2
	if (tp2.WorstPl>=0) tl.SubHit(tp2.WorstPl);
	//                                       ***** ASSIGN P/dP TO LAST HELIX
	tl.Hlast(5) = cop2; tl.Hlast(5,5) = d2cop2;
	tl.IMark = 0; tr.IMark = -1; // "IMark" flags: tl/tr to be kept/deleted
      }
    }
    for (it = listTrack.begin(); it!=listTrack.end(); it++) {
      //                                       ********** RESET >0 "IMark" FLAGS
      short int &iMark = (*it).IMark; if (iMark>0) iMark = 0;
    }
  } // End option for 1st loop on SM2 pairs

  for (itp1 = lTrackPair1.begin(); itp1!=lTrackPair1.end(); itp1++) {
    TTrackPair2 &tp1 = *itp1; TTrack &tl = *tp1.iL, &tr = *tp1.iR;

    //     *************** LOOP OVER SM1 TRACK PAIRS ***************

#ifdef BridgeS_DEBUG_STEAL
    if (tl.IKine==1 || tl.IKsame==1 && tl.NHsame>.75*tl.NHits ||
	tpSteal && tl.Id==(*(tpSteal->iL)).Id) {
      if (tr.IKine!=1 && tpSteal==0) { tpSteal = &tp1; ikSteal = tl.IKine; }
      if (tr.IKine==1 && tpSteal) {
	TTrack &tSl = *(tpSteal->iL), &tSr = *(tpSteal->iR);
	double chi2S = tpSteal->Chi2, pS = 1/tpSteal->Hfirst(5);
	printf("BridgeS_DEBUG_STEAL:"
	       "#%3d(%5.1fGeV,mc=%2d+%2d) %5.2fNDF < "
	       "#%3d(%5.1fGeV,mc=%2d+%2d) %5.2fNDF: %.2f/%.2f\n",
	       tSl.Id,pS,             tSl.IKsame,tSr.IKine,chi2S,
	       tl.Id, 1/tp1.Hfirst(5),ikSteal,   tr.IKine,tp1.Chi2,
	       chi2S/tp1.Chi2,fabs(tpSteal->Hfirst(5)/tp1.Hfirst(5)));
      }
    }
#endif

    //                                           ********** BASIC PAIR SELECTION
    if (tl.Type!=0x1 ||          // Already extended
	tl.IMark==-1) continue;  // Already prepended supra, or appended in BoT
    if (tr.Type!=0x2 || tr.IMark==-1) continue; // Already bridged/appended
    int nDFs0x1 = tl.NDFs; // ***** SAVE SOME OF THE ORIGINAL TRACK'S ATTRIBUTES
    if (tp.WorstPl>=0 && tp.WorstPl<setup.vIplFirst()[0])
      nDFs0x1 -= (setup.vPlane(tp.WorstPl).IFlag&0x30) ? 2 : 1;
    double xL0x1 = tl.Hlast(0), chi20x1 = tl.Chi2tot/(nDFs0x1-4);
    tl.IFit = tp1.IFit;                                         // ***** MERGING
    tl.Chi2tot =// Needed by "CleanTrackList" and infra when considering muWallA
      tp1.Chi2tot;
    if (tp1.IFit==0x8) {   // QN fit
      tl.Haux = tp1.Hfirst; tl.NDics = tp1.NDics; tl.Chi2aux = tp1.Chi2tot;
    }
    else tl.Hfirst = tp1.Hfirst; // KF fit
    tl.Append(tr);               // Append 0x2 to 0x1
#ifdef MULTIPLE_SCIFI_BRIDGING
    bool unreliableScifiTrack = false;
    if  ((tl.Scifi==0x3 || tl.Scifi==0x201) && tl.NDFs<11) {
      tr.IMark = -2; unreliableScifiTrack = true;
    }
#endif
    int worstHit1 = -1; if (tp1.WorstPl>=0) worstHit1 = tl.SubHit(tp1.WorstPl);
    double cop1 = tp1.Hfirst(5), d2cop1 = tp1.IFit!=0x8 ? tp1.Hfirst(5,5) :
      1e-4; // QN fit has no dcop: let's assume 1%
    //                                           ***** ASSIGN P/dP TO LAST HELIX
    tl.Hlast(5) = cop1; tl.Hlast(5,5) = d2cop1;

    // ********** LOOK IF THERE IS CONTINUATION AFTER SM2 **********

    for (itp2 = lTrackPair2.begin(); itp2!=lTrackPair2.end(); itp2++) {
      TTrackPair2 &tp2 = *itp2; TTrack &tl2 = *tp2.iL, &tr2 = *tp2.iR;
      //                                              ***** LOOP SM2 TRACK PAIRS
      if (tp2.iL!=tp1.iR) continue;  // Must have common 0x2-track
      if (tr2.IMark==-1) continue;   // Already appended
      if (tr2.IMark>0) continue;     // Already 0x4 piece of yoke track
      //           ***** FOUND A CONTINUATION *****
      //                                                ***** DISAGREEMENT IN P?
      double cop2 = tp2.Hfirst(5), d2cop2 = tp2.Hfirst(5,5);
      double cop21 = fabs(cop2-cop1), dcop = sqrt(d2cop1+d2cop2);
      int pMismatch; if (cop21>3*dcop) {
	pMismatch = 1;
	if (cop21>5*dcop) // A loose cut, in order to make for possible...
	  pMismatch = 2;  // ...tracking error, that can be fixed later on.
      }
      else pMismatch = 0;
      if ((tr.Scifi&0x202) &&  // VSAT-(almost)only? chi2 not very meaningful...
	  pMismatch>1 ||
	  tr.Hlast(0)<xRICH && // That the reconstruction ends prematurely
	  // upstream of RICH is very well possible: the tracking of particles
	  // through the large thickness of the RICH and beyond being tricky,
	  // cf. "TEv::PrePattern2". Remains that such cases are not the best
	  // candidates for a 0x4 extension. => Let's disregard them when they
	  // in addition show a 0x6 vs. 0x3 mismatch in momentum.
	  pMismatch>1) {
	// ...=> search SM2 track pairs for a better bridging of zone 0x4 track
	bool better; list<TTrackPair2>::iterator jtp2;
	for (jtp2 = lTrackPair2.begin(), better = false;
	     jtp2!=lTrackPair2.end(); jtp2++) {
	  if (jtp2->iR!=tp2.iR) continue;
	  if (jtp2->iL==tp2.iL) { break; }
	  if (jtp2->iL->Type!=0x2 ||       // Already extended
	      jtp2->iL->IMark==-1 ||       // Already appended
	      jtp2->iL->IMark>0)           // Already 0x2 piece of yoke track
	    continue;
	  better = true; break;
	}
	if (better) continue;   // If better bridging found => continue
      }

      //                                ***** TRACK GOING THROUGH SM2 YOKE *****
      // The distribution of material and field in SM2 yoke is not properly
      // described (as of 2008/05).  Extending 0x3 track through the yoke will
      // spoil it (spoil the determination of its momentum) => Cancel such an
      // extension.
      // (Note:
      //  - We do not reject yoke extension earlier, in the pairing step, for
      //   the 0x6-pairing is perfectly valid. Allowing it to take place, even
      //   though we do not retain it in the final end, helps to solve
      //   adequately the pairing of other tracks.)
      //  - Therefore let's reject yoke extensions now. The case, including
      //   0x6-only yoke tracks, will be re-evaluated in any case in
      //   "TracksFit2".)
      bool isYokeTrk = IsYokeTrk(tl2,tr2); if (isYokeTrk) {
	if (!TOpt::ReMode[41] ||// ReMode[41]: no bridging SM2 yoke requested...
	    // ...or 0x3 piece ending upstream of RICH (meaning that tracking
	    // hasn't been able to follow its trajectory from RICH down to SM2:
	    // making it doubtful the two (0x3 and 0x4) piece are related)
	    tr.Hlast(0)<750)
	  continue;
      }
      // else: Now that the case of yoke tracks has been dealt with one could...
      // ... discard "itp2" continuations w/ incompatible momentum, e.g.
      // if (dcop && cop21>10*dcop) continue;
      // but this would validate "itp1", whereas, infra, one considers the
      // possibility to validate "itp2" instead.

      //                       ***** DO NOT SPOIL VERY GOOD 0x6 w/ BAD 0x3 *****
      // - Opt for favouring 0x6 at the expense of 0x7, since the extension
      //  of 0x6 into 0x1 can later be achieved by "BackTrackSAS".
      // - As to decide between 0x3 and 0x6, the case is less straightforward:
      //   In disfavour of 0x3: bad chi2, w/ high P requiring better chi2. When
      //  deciding for a particular numeric cut, one has to remember that a 0x3
      //  track falling into the acceptance of SM2 (sine qua no possible 0x6)
      //  has P typically >20 GeV, and hence must have relatively low chi2, even
      //  when MS is not accounted (which is the case here w/ "QuickKF").
      //   VSAT-(almost)only is a special case, w/ few hits => even lower chi2.
      double chi20x3 = tl.Chi2tot/(tl.NDFs-5);
      int bad0x3; if (chi20x3>badChi2 ||
	  // >40 GeV, i.e. the medium to high range of tracks making it to SM2.
	  fabs(cop1)<.025 && chi20x3>badChi2/2) {
	bad0x3 = 1;
	if (chi20x3>badChi2*1.5 || fabs(cop1)<.025 && chi20x3>badChi2)
	  bad0x3 = 2;
      }
      else if (chi20x3-(chi20x1+tr.Chi2tot/(tr.NDFs-4))/1>3 && (tl.Scifi&0x1))
	// To rescue a 0x6 track in #8475622 of 2010/evtDump.FI.GP.2-85821.raw.
	// To solve this kind of case (where 0x3 is not that bad), would be
	// better to remember the 0x2<->0x4 association, and decide later. E.g.,
	// for the case supra, until CleanTrackList has stripped the 0x3 track
	// of its 0x1 ghost head.
	bad0x3 = 1;
      else bad0x3 = 0;
      bool badVSAT0x3 = (tl.Scifi&0x101) && (tr.Scifi&0x202) &&
	chi20x3>badVSATChi2;
      double chi20x6 = tp2.Chi2tot/(tr.NDFs+tr2.NDFs-5);
      double dT0x2 = tr.SigmaTime, dT0x4 = tr2.SigmaTime;
      bool tMismatch = false; int good0x6 = 0; if (dT0x2>0 && dT0x4>0) {
	if (fabs(tr.MeanTime-tr2.MeanTime)>3*sqrt(dT0x2*dT0x2+dT0x4*dT0x4)) {
	  good0x6 = 0; tMismatch = true;
	}
	else if (chi20x6<goodChi2) good0x6 = 1;
      }
      if (good0x6 && ((bad0x3==2 || badVSAT0x3) ||
		      bad0x3 && pMismatch)) {
	tr.IMark = 0;   // In any case, free 0x2 track
	if (chi20x3>TOpt::dCut[17]*1.5)
	  tl.IMark = -1;      // If 0x1 track was bad, delete 0x3
	else {                // If not, restore it
	  if (worstHit1>=0) { tl.AddHit(vecHit[worstHit1]); worstHit1 = -1; }
	  tl.Shorten(setup.vIplFirst()[1]); tl.Type = 0x1; tl.Scifi &= 0x101;
	  tl.Hfirst(5)=tl.Hlast(5) = 0;
	  tl.QuickKF(-1,0); tl.QuickKF(1,0); tl.IFit = 0x1;
	  tl.FindKine();
	}
	break;
      }
      tl.Append(tr2);                                        // ***** APPEND 0x4
      int worstHit2 = -1; if (tp2.WorstPl>=0 && tp2.WorstPl!=tp1.WorstPl)
			    worstHit2 = tl.SubHit(tp2.WorstPl);
      //                                                           ***** FIT 0x7
      tl.Haux = tp1.Hfirst;                      // Save 0x3 momentum
      tl.Hlast(5) = cop2; tl.Hlast(5,5) = d2cop2; // Initial P: from 0x6 fit
      int ok; if (TOpt::ReMode[46]) {
	tl.InsertMissedPlanes(); if (tl.FullKF(-1)) { ok = 1; tl.IFit = 0x20; }
	else                                        { ok = 0; tl.IFit = 0; }
      }
      else if (tl.QuickKF(-1,1))                    { ok = 1; tl.IFit = 0x2; }
      else                                          { ok = 0; tl.IFit = 0; }

      // Mismatch between SM1 and SM2 fits:
      //   i) Decay that occurs in zone 0x2, leading a good 0x3 track to be
      //     spoiled by its 0x6 extension.
      //  ii) High multiscattering, e.g. through HCAL and/or muAbsober for mu.
      // iii) Interaction downstream of vertex detectors (FI03/4 and the like),
      //     while 0x1 track hits from these. => The small lever arm of the LAS
      //     makes it possible that the scattering angle is absorbed by the fit
      //     in the bending due to SM1 and the multiscattering, yielding a
      //     fairly good 0x3 chi2.
      //  iv) Purely accidental.
      // - Identified:
      //   - I) For scifi-(almost)only tracks, by bad overall chi2 or bad chi2
      //       increment when going from 0x3+0x4 to 0x7
      //   - II) Then by good chi2 for 0x3/0x6 pieces, and much worse for 0x7.
      // - Still chi2 is normalised by (#hits-5) + quad terms in #hits so as to
      //  allow for much multiscattering if the 0x7 track is long.
      // - Several rescue procedures:
      //   I) Downtream interaction: we look a bridging (of the 0x1 segment) w/
      //     the beam telescope. If found, halt the building a 0x7 track.
      //     Leave untouched the 0x2 and 0x4 pieces (to be bridged later on in
      //     the "BridgeSegments2" flow-chart. Possibly build an extended 0x11
      //     beam track.
      //  II) Strip away the 0x6 piece: various requirements being successively
      //     tried
      double chi2_34 = (chi20x3+tr2.Chi2tot/(tr2.NDFs-4))/2;
      double chi2_16 = (chi20x1+chi20x6)/2;
      double chi20x7 = tl.Chi2tot/(tl.NDFs-5);
      double yL = tr.Hlast(1), zL = tr.Hlast(2), r2Last = yL*yL+zL*zL;
      bool forwardTrack = (tl.Scifi&0x101) &&// VSAT-(almost)only in zone 0x1...
	(tl.Scifi&0x202 ||         // ...and VSAT-(almost)only also in 0x2 or...
	 // ...at least reasonably paraxial, which precludes large bend in SM1
	 // and hence a low momentum track spanning zones 0x3
	 r2Last<64); // somewhat arbitray: tuned on 08W33/69853#7345348
	 //  => Cannot but have high momentum => expected to have low chi2,
	 //    even w/o multiscattering taken into account.
      bool highIncrement = // Ambiguous forward track case: reasonably low total
	// chi2, but large chi2 increment when going from 0x3 or 0x6 to 0x7...
	ok && forwardTrack && chi20x7<=TOpt::dCut[17] &&
	(chi20x7>chi2_34+2*chi20x7IncrCut && chi20x7>chi2_16+chi20x7IncrCut ||
	chi20x7>chi2_34+15);
      bool msProne = !pMismatch && !tMismatch && fabs(cop2)>.04 &&
	TOpt::dCut[17]<chi20x7 && chi20x7<TOpt::dCut[17]*1.5 &&
	// Large #hits, and hence of detectors traversed, possibly high MS (
	// energy loss being anyway neglected) impacting on the trajectory in
	// zone 0x4. The cut on #hits>26 means traversing more than all GEMs in
	// zone 0x2 (=24), i.e. traveling in the overlap of GEMs and MWPCs.
	tr.NHits>26 && tr2.NHits>25;
      if (ok &&
	  (highIncrement || msProne) && // ***** ...PERFORM FULL KF on 0x7 TRACK
	  !TOpt::ReMode[46]) {          //       ...when QUICK KF was enabled     
	if (tl.FullKF(-1)) { // Backward KF, since SM2 momentum more reliable
	  double chi2 = tl.Chi2tot/(tl.NDFs-5);
	  CsErrLog::msg(elError,__FILE__,__LINE__,
 "Evt #%d Track #%d = (0x3,%.2f GeV,%.2fNDF) + (0x6,%.2f GeV,%.2fNDF) "
 "Full KF -> %.2fNDF",event,tl.Id,1/cop1,chi20x3,1/cop2,chi20x6,chi20x7);
	  if (chi2<TOpt::dCut[16]*1.5) {  // If good chi2, i.e. below lower
	    // cut for track cleaning, or no more than 50% higher (hoping that
	    // future track cleaning will be able to fix that)...
	    //                                  ***** ... VALIDATE THE 0X7 TRACK
	    if (highIncrement) highIncrement = false;
	    // For the time being, do not apply the following
	    //t->Hlast(5) = t->Hfirst(5); // Account for improved P estimate
	    //t->IFit = 0x20; // Remember FullKF has been executed
	    //t->Chi2tot = // Full KF may give an undue advantage to current...
	    // ...tracks when it comes to CleanTrackList => make a trade-off.
	    // (t->Chi2tot+chi20x7*(t->NDFs-5))/2
	    if (msProne) chi20x7 = TOpt::dCut[17]-.01;
	    tl.IFit = 0x20;
	  }
	  else { ok = 0; tl.IFit = 0; }
	}
	else { ok = 0; tl.IFit = 0; }
	if (!ok) CsErrLog::msg(elError,__FILE__,__LINE__,
 "Evt #%d Track #%d = (0x3,%.2f GeV,%.2fNDF) + (0x6,%.2f GeV,%.2fNDF) "
 "Full KF fit fails",event,tl.Id,1/cop1,chi20x3,1/cop2,chi20x6);
      }
      if (ok && forwardTrack && (chi20x7>TOpt::dCut[17] || highIncrement)) {
#ifdef BridgeS_OVER_TARGET
	if (chi20x6<TOpt::dCut[16] && // If 0x6 also fine...
	    //                       ...I) Could it be a downstream interaction?
	    // (Restricting ourself to scifi-(almost)only and upon option.)
	    (targetBridging&0x8)) {
	  int nPars = targetField<2 ? 4 : 5;
	  for (itp0 = lTrackPair0.begin(); itp0!=lTrackPair0.end(); itp0++) {
	    // Loop over target track pairs
	    TTrackPair2 &tp0 = *itp0; TTrack &t0x10 = *tp0.iL;
	    if (tp0.iR!=tp1.iL) continue;        // Must have common 0x1-track
	    if (t0x10.Type&0x1) continue;        // Already extended
	    double chi20x10 = t0x10.Chi2tot/(t0x10.NDFs-4);
	    double chi20x11 = tp0.Chi2tot/(t0x10.NDFs+nDFs0x1-nPars);
	    if (chi20x11<(chi20x1+chi20x10)*1.5 ||// Small chi2 increment...
		// ...Or very small absolute chi2 (in particular to make for the
		// fact that chi2 may be understimated in zone 0x1 because of
		// option "dCut[84]"). The actual cut is arbitrary. Would need
		// a re-assessment, in particular in the muon setup case.
		TOpt::dCut[84] && chi20x11<.5) {
	      // Good enough 0x11 bridging => get vertex' X from Z intersect
	      THlx hl = t0x10.Hlast, hr = tr.Hfirst;
	      double x0l = hl(0), z0l = hl(2), zpl = hl(4);
	      double x0r = hr(0), z0r = hr(2), zpr = hr(4);
	      double xV = (z0r-zpr*x0r-z0l+zpl*x0l)/(zpl-zpr);
	      if (hl.FindCDA(hr,xV,tl.Hfirst(0)+10,x0r)) xV = hl(0);
	      else                                       xV = tl.Hfirst(0)+10;
	      if (xV>tl.Hfirst(0)) {
		ok = -1; tr.IMark=tr2.IMark = 0; // Release 0x6 segments
		if (xV>xL0x1) { // Vertex downstream end of 0x1 segment?..
		  // ...=> Do actually bridge over target
		  tl.Shorten(setup.vIplFirst()[1]);
		  tl.Type = 0x1; tl.Scifi &= 0x101; 
		  t0x10.Append(tl); 
		  t0x10.Hfirst = tp0.Hfirst;
		  bool ret = targetField?t0x10.QuickKF(1,1):t0x10.QuickKF(1,0);
		  if (ret)
		    ret =   targetField?t0x10.QuickKF(-1,1):t0x10.QuickKF(-1,0);
		  if (!ret) {
		    t0x10.Shorten(setup.vIplFirst()[0]); t0x10.Type = 0x10;
		    t0x10.Scifi = 0; tl.IMark = 0;
		  }
		  else {
		    CsErrLog::msg(elWarning,__FILE__,__LINE__,
 "Downstream interaction: @ Z~=%.2f, between #%d(0x11) and #%d(0x7)",
				xV,t0x10.Id,tr.Id);
		    t0x10.IFit = tp0.IFit;
		  }
		}
		else {
		  tl.IMark = -1; // Flag 0x1 track for deletion
		  CsErrLog::msg(elWarning,__FILE__,__LINE__,
 "Downstream interaction: @ Z~=%.2f, between #%d(0x10) and #%d(0x7)",
				xV,t0x10.Id,tr.Id);
		}
	      }
	    }
	  } // End loop over target track pairs
	} // End checking for downstream interaction
	if (ok>0) ok = 0;
#else
	ok = 0;
#endif
      }
#define BridgeS_CUT_MISMATCH
#ifdef BridgeS_CUT_MISMATCH
      if (ok>0 && (chi20x3<TOpt::dCut[16]/2 && chi20x7>TOpt::dCut[9]*1.5 ||
		 // Case of extremely bad 0x7 chi2. E.g. evt #4 of
		 // "/castor/~guskov/zebra_august_07_pi16/zebradat.0.fz"
		 chi20x3<TOpt::dCut[16]*2 && chi20x7>TOpt::dCut[9]*2.5)) {
	THlx &hrL = tr.Hlast; float xL = hrL(0);
	if (hrL(0)<750) ok = 0;            // ...IIa) 0x2-track ends before RICH
	else {
	  float reducedChi2 = TOpt::ReMode[46] ? tl.Chi2tot/(tl.NDFs-5) :
	    tl.Chi2tot/(tl.NDFs-5+tl.NHits*tl.NHits*.01);
	  if (reducedChi2>TOpt::dCut[9]*1.5) {       // Consider reduced chi2...
	    float xMAL = TOpt::MuonWall[0]-xL;// ...Check position @ muWallA
	    // against central hole in absorber (the latter = dead zone of MAs ~= 70x40)
	    if (fabs(hrL(1)+xMAL*hrL(3))>TOpt::MuonWall[8] ||
		fabs(hrL(2)+xMAL*hrL(4))>TOpt::MuonWall[9]) {
	      if (reducedChi2>TOpt::dCut[9]*2) // ...IIb) If outside: higher cut
		ok = 0;
	    }
	    else ok = 0;                       // ...IIc) If w/in: standard cut
	  }
	}
      }
#endif
      if (ok>0 && !isYokeTrk) {
#ifdef MULTIPLE_SCIFI_BRIDGING
	if (unreliableScifiTrack) tr2.IMark = -2;
#endif
	tl.Hlast(5)   = tl.Hfirst(5);   // Store momentum @ last point Helix
	tl.Hlast(5,5) = tl.Hfirst(5,5);
      }
      else if (ok>=0) {
	int used0x1 = 0; if (ok==0 &&
			     chi20x6<TOpt::dCut[16]/2 && (tl.Scifi&0x1)) {
	  // At this step (0x7 fit fails) while 0x6 is OK, one can question the
	  // validity of the 0x1 segment if it's scifi-(almost)only: such track
	  // being by construction unreliable given the low vSAT redundancy.
	  // => Let's check whether its hits aren't used by other bridged tracks
	  // (which are bound to be more reliable since they've shown up earlier
	  // in the chi2-ordered list of track pairs).
	  // (Note: this block has been added lately, to address a specific
	  // reco inefficiency and has not been tuned. Its basic idea (checking
	  // for used hits) may be worth generalising.)
	  list<int>::iterator ih; int nUsed, nH0x1; list<TTrack>::iterator it;
	  set<unsigned int, less<unsigned int> >::iterator iID;
	  int iplMx = setup.vIplLast()[0];
	  for (ih = tl.lHitPat.begin(), nUsed=nH0x1 = 0; ih!=tl.lHitPat.end();
	       ih++) {
	    if (*ih<0) continue;
	    THit &h = vecHit[*ih]; if (h.IPlane>iplMx) break; nH0x1++;
	    int used; for (iID = h.sTrackID().begin(), used = 0;
			   iID!=h.sTrackID().end(); iID++) {
	      for (it = listTrack.begin(); it!=listTrack.end(); it++) {
		TTrack &ti = *it; if (ti.Id!=*iID) continue;
		if (ti.Type&0x2) used = 1; break;
	      }
	      if (used) { nUsed++; break; }
	    }
	  }
	  if (nUsed>nH0x1-4) used0x1 = nUsed==nH0x1 ? 2 : 1;
	}
	if (ok==0 &&
	    // ***** 0x6 REASONABLY GOOD while 0x3 BAD ENOUGH...
	    (chi20x6<TOpt::dCut[16]*2 && TOpt::dCut[17]<chi20x3 ||
	     // ...OR 06 VERY GOOD while 0x3 QUESTIONABLE...
	     chi20x6<TOpt::dCut[16]/2 && 
		     (used0x1 && TOpt::dCut[16]/2<chi20x3 || used0x1==2))) {
	  tr.IMark=tr2.IMark = 0;           // ...Release 0x6 segments
	  tl.Shorten(setup.vIplFirst()[1]); // ...Rescue 0x1
	  if (worstHit1>=0 &&            // If 0x3's worst hit was subtracted...
	      tp1.WorstPl<setup.vIplFirst()[1]) { // ...and belong to zone 0x1
	    tl.AddHit(vecHit[worstHit1]); worstHit1 = -1; // => Restore it
	  }
	  tl.Type = 0x1; tl.Scifi &= 0x101; tl.QuickKF(1,0); tl.QuickKF(-1,0);
	  tl.IFit = 0x1;
	  tl.FindKine();
	}
	else {  // ***** ELSE RETAIN THE 0x3 PIECE
	  tl.Shorten(setup.vIplFirst()[2]); tl.Type &= 0x13; tl.Scifi &= 0x1313;
	  if (worstHit2>=0 &&            // If 0x6's worst hit was subtracted...
	      tp2.WorstPl<setup.vIplFirst()[2]) { // ...and belong to zone 0x2
	    tl.AddHit(vecHit[worstHit2]); worstHit2 = -1; // => Restore it
	  }
	  tl.Chi2tot = tp1.Chi2tot;
	  tl.Hlast = tr.Hlast; // "tr" has kept the 0x2 end of the 0x3 piece
	  tl.Hlast(5) = cop1; tl.Hlast(5,5) = d2cop1; // Complete it w/ P/dP.
	  tl.Hfirst = tp1.Hfirst; tl.IFit = tp1.IFit;
	  tl.FindKine();
	  if (ok && isYokeTrk) { // Yoke track (0x4 piece thereof)
	    // Remember association 0x3 <-> 0x4
	    tr2.IMark = tl.Id; tl.IMark = tr2.Id;
	    // Assign to it the momentum of 0x3 piece, w/ infinite precision, so
	    // as to insulate against any further modification.
	    tr2.Hfirst(5)=  tr2.Hlast(5)   = tl.Hlast(5);
	    tr2.Hfirst(5,5)=tr2.Hlast(5,5) = 1.E-20;;
	  }
	  else tr2.IMark = 0;
	}
      }
      break;
    }
#ifdef BridgeS_OVER_TARGET
    if (targetBridging&0x11)
      BridgeOverTarget(lTrackPair0,tl,nDFs0x1,chi20x1);
#endif    
  }  // End loop over SM1 track pairs

  for (itp2 = lTrackPair2.begin(); itp2!=lTrackPair2.end(); itp2++) {

    //    *************** LOOP OVER SM2 TRACK PAIRS ***************

    TTrackPair2 &tp2 = *itp2; TTrack &tl = *tp2.iL, &tr = *tp2.iR;
    if (tl.Type!=0x2 ||           // Already extended
	tl.IMark==-1) continue;   // Already appended

    if (tr.Type!=0x4)    // The condition below NEVER occurs.
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
		    "Evt #%d Trk #%d member of SM2 pairs has type 0x%x",
		    event,tr.Id,tr.Type);
    if (tr.IMark==-1 ||          // Already appended
	tr.IMark>0) continue;    // Already 0x4 piece of yoke track

    if ((tl.Scifi&0x202) && (tr.Scifi&0x404)) { // If hint at precise timing...
      //       *****  CHECK TIMING COMPATIBILITY... *****
      // - Fake association is more likely to happen with (0x2+0x4) than supra
      //  where we deal w/ 0x3+0x4.
      double sl = tl.SigmaTime, sr = tr.SigmaTime, dT =
	sqrt(sl*sl+sr*sr); // "tl/tr", being "Scifi" flagged do possess time
      
      if (dT<.5 && // Check that we have indeed good timing (cut corresponds...
	  // ...to 2 scifis at least in both "tl" and "tr").
	  fabs(tr.MeanTime-tl.MeanTime)>6*dT) { // Loose cut on compatibility
	//     ***** ...IF TURNS OUT to be BAD: LOOK for ALTERNATIVE *****
	bool better; double bestChi2 = tp2.Chi2;
	list<TTrackPair2>::iterator jtp2 = itp2;
	for (jtp2++, better = false; jtp2!=lTrackPair2.end(); jtp2++) {
	  TTrackPair2 &tp2p = *jtp2; TTrack &trp = *tp2p.iR;
	  if (tp2p.iL!=tp2.iL) continue;   // Same 0x2 track segment
	  if (trp.IMark==-1) continue;     // Already appended
	  sr = trp.SigmaTime;        // ...=> Let's require good timing again
	  if (!(trp.Scifi&0x404) || (dT = sqrt(sl*sl+sr*sr))>.5) continue;
	  if (fabs(trp.MeanTime-tl.MeanTime)<4*dT &&// "tl/trp" time compatible
	      // ...and yet tracks pair's chi2 close to best (Note: would
	      // probably be better to cut on chi2 increment (current-best)
	      // based on #scifis.)
	      (*jtp2).Chi2<bestChi2*1.5) {
	    better = true; break;
	  }
	}
	if (better) continue;   // If better bridging found => continue
      }
    }

    //               ***** SM2 YOKE TRACKs *****
    // - Cf. definition supra.
    // - The yoke tracks reco'd in SAS alone (as opposed to those also reco'd in
    //  LAS) have no unbiased momentum. (The sign can even be wrong if the track
    //  travels into the region of the yoke where, in reality, field lines close
    //  the magnetic circuit, while, in field maps, they don't.)
    // - Yet, they are genuinely bridged, or almost genuinely if one considers
    //  that the details of the field map matters. Therefore let's still keep
    //  those that will not interfere much w/ the rest of the reco, e.g.
    //  vertexing (but RICH reco will be wrong...), i.e. halo tracks.
    // - For the time being, let's keep them as is. Will think of a mean to
    //  tag them later on...
    if (IsYokeTrk(tl,tr)) {
      bool haloTrack = fabs(tl.Hfirst(3))<.03 && fabs(tl.Hfirst(4))<.03;
      if (!haloTrack) continue;
    }

    tl.IFit = tp2.IFit;
    tl.Chi2tot = // Would be needed by "CleanTrackList"(0) if bridging...
      tp2.Chi2tot; // ...resulted from QN fit. Improbable, but let's dot it...
    tl.Hfirst = tp2.Hfirst;
    tl.Append(tr);  // Append 0x4 to 0x2
    if (tp2.WorstPl>=0) tl.SubHit(tp2.WorstPl);
    tl.Hlast(5)  = tp2.Hfirst(5);// Store momentum in last point Helix
    tl.Hlast(5,5)= tp2.Hfirst(5,5);
  }  // End loop over SM2 track pairs

#ifdef BridgeS_OVER_TARGET
  if (targetBridging&0x2) {

    //      ********** BRIDGE BEAM w/ ZONE 0x1 **********

    // This is meant for hadron data, where interaction often takes place
    // downstream of main target. We limit the bridging to those segments
    // that have not beam bridged over SM1 (0x1 tracks), the rest having
    // been dealt with supra (upon option).
    // (But could be written so as to accomodate various other cases.)
    for (itp0 = lTrackPair0.begin(); itp0!=lTrackPair0.end(); itp0++) {
      TTrackPair2 &tp0 = *itp0; TTrack &tl = *tp0.iL, &tr = *tp0.iR;

      // ***** LOOP OVER TARGET TRACK PAIRS
      if (tl.Type!=0x10) continue; // Already extended
      if (tr.Type!=0x1) continue;  // Already appended
      float chi2 = itp0->Chi2tot/(tl.NDFs+tr.NDFs-4);
      if      (chi2>7) continue; // Strict chi2 cut 
      else if (chi2>3) {         // Try cleaning away FI03...
	// ...for it might be that it is the origin of the interaction and track
	// has picked up a hit from one of the outgoing track.
	static int iplFI03i = -1; static float xFI03ante, xFI03post;
	if (iplFI03i==-1) {
	  iplFI03i = -2; if (!(targetBridging&0x4)) continue;
	  for (int ipl = setup.vIplFirst()[0]; ipl<setup.vIplLast()[0]; ipl++) {
	    const TPlane  &p = setup.vPlane(ipl);
	    const TDetect &d = setup.vDetect(p.IDetRef);	    
	    if (d.Name.find("FI03")==0) {
	      if (iplFI03i<0) { iplFI03i = ipl; xFI03ante = d.X(0)-.01; }
	      xFI03post = d.X(0)+.01;
	    }
	  }
	}
	if (iplFI03i>=0) { // Used a flag summarising all what's requisite
	  if (tr.Hlast(0)<xFI03ante ||         // Must go up to FI03...
	      xFI03post<tr.Hlast(0)) continue; // ...but not beyond
	  TTrack trp(tr); trp.Shorten(iplFI03i);
	  tl.Append(trp); tl.QuickKF(-1,0);
	  if (tl.Chi2tot/(tl.NDFs-4)>3) {                  // Cleaning failed...
	    tl.Shorten(setup.vIplFirst()[0]);                   // ...restore...
	    tl.Type = 0x10; tl.QuickKF(1,0); tl.QuickKF(-1,0);  // ...and update
	    tl.IFit = 0x1;
	    continue;
	  }
	}
	else continue;
      }
      else {
	tl.Hfirst = itp0->Hfirst;
	tl.Append(tr);  // Append track type 0x1 to track type 0x10
	if (itp0->WorstPl>=0) tl.SubHit(itp0->WorstPl);
      }
      tl.QuickKF(1,0); tl.IFit = 0x1;
      if (!TOpt::dCut[4])
	CsErrLog::mes(elFatal,
    "Bridging over target of 0x11 tracks requested while dCut[4]==0");
      tl.Hfirst(5)=tl.Hlast(5) = TOpt::iCut[15]/TOpt::dCut[4];
    }
  }
#endif

  // ********** CLEAN TRACK LIST

  // Erase tracks pieces and fit single-segment target->SM1 tracks
  it = listTrack.begin(); while (it!=listTrack.end()) {// ***** Track loop *****
    if (it->IMark<0)
      listTrack.erase(it++);     // Erase appended track piece
    else if (it->Type==0x1 && it->NDFs<=5 ||  // Single-segment
	     it->Type==0x2 && it->NDFs<=4)    // 4:5-(scifi)hit tracks:
      // - they are unreliable and one cannot do anything about that
      // - they cannot be reasonably QN fitted
      // - they have chi2~=0 => one cannot perform any meaningful selection
      // - they cannot help in the cleaning process for they convey absolutely
      //  no PR info.
      listTrack.erase(it++);     // => erase 
    // else if #hits==6: the corresponding tracks are not so useful either but
    // we keep them for they may in the cleaning process to come.
    else if (it->Type==0x1 &&                  // Single-segment target->SM1
	     TOpt::ReMode[14] && !(it->IFit&0x8)) {
#warning TODO: Which init value fr QN refit?
      it->QNewtonFit(1,1);                        // ...Update QN Fit
      it++;
    }
    else
      it++;
  }   // end of track loop

#ifdef DEBUG
  if (TOpt::Graph[0]>0 && TOpt::Graph[6]>1) {
    cout<<"Before cleaning 0"<<endl;
    TDisplay::Ref().Draw(3); TDisplay::Ref().Draw(4);
  }
  static int idebug = 0;
  if (idebug) DumpTrackList(); if (idebug>1) DumpHits();
#endif
  if (TOpt::ReMode[15] && TOpt::ReMode[14])
    // When several tracks share a common hit within QN reach: retain only best
    CleanTrackList(0);

  //                                            ********** FITTING... **********
  Traffic::Ref().Stopwatch.Start(13); // Time used by fitting w/in bridging step
  it = listTrack.begin(); while (it!=listTrack.end()) {
    TTrack &t = *it;
    if (((t.Type&0x3)==0x3 || (t.Type&0x6)==0x6) && // 2(or more)-ZONE TRACKS...
	!(t.IFit&0x22)) {                           // ...not yet KF fitted
      if (t.IFit==0x8) {
	t.Hlast(5) = t.Haux(5); t.Hlast(5,5) = 1.E-4;
      }
      bool ret; if (TOpt::ReMode[46]) {
	t.InsertMissedPlanes(); ret = t.FullKF(-1);    t.IFit = 0x20;
      }
      else {                    ret = t.QuickKF(-1,1); t.IFit = 0x2; }
      if (ret) {          // ***** KF BACKWARD *****
	if ((t.IFit&0x8) && t.Chi2tot/(t.NDFs-5)>TOpt::dCut[16]) {
	  // If QN fit => may be useful to iterate...
	  int iter = 0; float chi2 = 1.3*t.Chi2aux/(t.NDics-5);
	  while (ret && iter<2 && t.Chi2tot/(t.NDFs-5)>chi2) {
	    t.Hlast(5) = t.Hfirst(5); t.Hlast(5,5) = t.Hfirst(5,5);
	    iter++;
	    ret = TOpt::ReMode[46] ? t.FullKF(-1) : t.QuickKF(-1,1);
	  }
	}
      }
      if (!ret) { listTrack.erase(it++); continue; }
    }
    it++;
  }
  Traffic::Ref().Stopwatch.Stop(13);

#ifdef DEBUG
  if (TOpt::Graph[0]>0 && TOpt::Graph[6]>1) {
    cout<<"After cleaning 0"<<endl;
    TDisplay::Ref().Draw(3); TDisplay::Ref().Draw(4);
  }
  if (idebug) DumpTrackList(); if (idebug>1) DumpHits();
#endif

  if (TOpt::ReMode[15])
    // When several tracks share a common hit: retain only best
    CleanTrackList(1);

#ifdef DEBUG
  if (TOpt::Graph[0]>0 && TOpt::Graph[6]>1) {
    cout<<"After cleaning 1"<<endl;
    TDisplay::Ref().Draw(3); TDisplay::Ref().Draw(4);
  }
  if (idebug) DumpTrackList(); if (idebug>1) DumpHits();
#endif

  //              ********** TRACK'S TIME **********

  it = listTrack.begin(); while (it!=listTrack.end()) {
    TTrack &t = *it;
    int clean; // ***** HIT TIME CLEANING *****
    //  The idea is to prevent badly timed hits from upsetting track's time.
    // This ahead of backtracking and vertexing. A positive by-product being
    // getting rid of these badly timed hits, which are also probably ghosts.
    //  W/in any single zone, scifi timing consistency is already pretty much
    // ensured by the PR algorithm (cutting on time spread among scififs).
    //  => Let's consider only bridged tracks
    if ((t.Type&0x3)==0x3 || (t.Type&0x6)==0x6) {
      //  We need precisely timed hits => means scifi and hodo hits. And we need
      // also precisely reliably timed track. These requirements are enforced in
      // "TTrack::UseHitTime" (by taking into account both hit's and track's
      // time uncertainties). Yet, in order to be on the safe side and avoid
      // e.g. throwing away a hodo hit that could be helpful in the scattered mu
      // ID, a two-gear approach is considered for the cleaning process. First
      // gear is for common tracks: cleaning is strictly conditionned by a tight
      // cut on the time distance (in unit of sigma) between the faulty hit and
      // the track. 2nd gear is for very well timed tracks, which we single out
      // on view of "TTRack::Scifi"
      clean = t.Scifi&0xff ? 6 : 15;
    }
    else clean = 0;
    if (!t.UseHitTime(clean)) { // (Note: Have to foresee a better management of this timing, w/ flags telling track's version and timing's version.)
      listTrack.erase(it++); continue;
    }
    it++;
  }

  if (TOpt::ReMode[18]) {

    // ******************** BACKTRACKING ********************

    it = listTrack.begin(); while (it!=listTrack.end()) {
      TTrack &t = *it;      // ***** ERASE BAD SCIFI-(ALMOST)ONLY TRACKS *****
      if (t.Scifi==0x3 && t.Type==0x3 &&           // I) in ZONES 0x3...
	  // These tracks are:
	  //  - Cannot but have high momentum => expected to have low chi2, even
	  //   w/o multiscattering taken into account.
	  //  - Not very reliable, because having very few hits.
	  //  - Potentially harmful, because possibly stealing hits from other
	  //   tracks, in particuliar from that of the scattered muon (when the
	  //   latter failed to be reco'd because of scifi inefficiency (e.g.
	  //   in zone 0x2) and can only be recovered by backtracking infra...
	  t.Chi2tot/(t.NDFs-5)>TOpt::dCut[17]) {       // ...when HIGH CHI2
	listTrack.erase(it++); continue;
      }
      else if (t.Scifi==0x2 && t.Type==0x2 &&       // II) in ZONE 0x2...
	       t.NHits<=6) {                           // ...when VERY FEW HITS 
	// (A pure scifi-only 0x2-track can have at most 5 hits in 2002:2004.
	// And 7 hits in 2006.)
	listTrack.erase(it++); continue;
      }
      it++;
    }

    if (TOpt::ReMode[18]&0x4)
      BackTrackZ2();  // ********** BACKTRACK 0x4-TRACKS INTO ZONE0x2 **********

    if (TOpt::ReMode[18]&0x8)
      DoubleBridge();

    if (TOpt::ReMode[18]&0x3) {  // ********** BACKTRACK INTO ZONE0x1 **********
      // ***** REMOVE ``SCIFI 0x1 TRACKS'' *****
      // i.e. ZONE1_only  segments w/ MANY SciFi hits, in order to free them
      // ***** REMOVE ``BAD 0x1 TRACKS'' *****
      int nBTZ1s = 0; // ***** DETERMINE #TRACKS CANDIDATES TO BackTrackZ1 *****

       // #tracks candidate to "BackTrackZ1"
      it = listTrack.begin(); while (it!=listTrack.end()) {
	TTrack &t = *it; if (t.Type==0x1) {            // ZONE_1 tracks
#define BackTrack_CHI2MAX 3
	  bool goodTrk =  (t.IFit&0x8) &&
	    t.Chi2aux/(t.NDics-5)<BackTrack_CHI2MAX +
	    2*fabs(t.Haux(5)) + // Relax chi2 cut for low momenta (cf. d*.2006.02 evt #87)
	    (t.NDics>12?1:0) && // Relax chi2 cut if #hits > 12, i.e. > #MMs: it could the be a track w/ SIs => the case is re-evaluated infra.
	    (t.NDics>=8 && t.Hlast(0)<230 || // At this point, allow for SI+GP tracks to pass through w/ a reasonable efficiency requirement: the requirement on #hits is re-evaluated infra.
	     t.NDics>=10);
	  int nScifis = 0, nPGs = 0, nSIs = 0, nStations = 0, nSpacePts = 0;
	  float nAbscissae = 0;
	  TStation *sPrv = 0; double xPrv = -1000; int projPrv = 0;
	  int projs = 0, nProjs = 0, allProjs = 0, nAllProjs = 0;
	  if (goodTrk) {
	    list<int>::iterator ih = t.lHitPat.begin();
	    while (ih!=t.lHitPat.end()) {
	      if (*ih>=0) {  // loop over the hits
		const THit &h = vecHit[*ih];
		const TPlane  &p = setup.vPlane(h.IPlane);
		const TStation *&s = p.Station; if (s!=sPrv) {
		  sPrv = const_cast<TStation*&>(s); nStations++;
		  if (nProjs>=2) nSpacePts++; nProjs=projs = 0;
		}
		int proj = 1<<p.IProj;
		if ((proj&projs)==0)    { projs |= proj;    nProjs++; }
		if ((proj&allProjs)==0) { allProjs |= proj; nAllProjs++; }
		const TDetect &d = setup.vDetect(p.IDetRef);
		if (d.IType==22) nScifis++; if (d.IType==21) nSIs++;
		if (d.IType==29) {
		  nPGs += 2;
		  if ((0x2&projs)==0)    { projs |= 0x2;      nProjs++; }
		  if ((0x2&allProjs)==0) { allProjs |= 0x2;   nAllProjs++; }
		}
		if (d.IType==28) nPGs++;
		nAbscissae += proj==projPrv && d.X(0)-xPrv<1 ? .51 : 1;
		xPrv = d.X(0); projPrv = proj;
	      }
	      ih++;
	    }
	    if (nProjs>=2) nSpacePts++;
	    goodTrk &= nStations>=3 && nAllProjs>=3 &&
	      // And, in order to exclude cases where a 3-(or 4-)ple plus a
	      // segment in only one projection, like, e.g., DC01XYUV+MM01:3Y
	      nSpacePts>=3;
	    if (nAbscissae<10 && nPGs<2)
	      //  i) #hits replaced by #abscissae, considering a DC double-layer
	      //    corresponds to less than 2 abscissae.
	      // ii) Enforce the #hits>=10 condition (that was evaded in case
	      //   zLast<< supra) if it turns out track is not a SI(assumed)+PG
	      //   track (anyway it must have been already been enforced
	      //   indirectly by the requirements on # of space points...)
	      goodTrk = false;
	    if (t.Chi2aux/(t.NDics-5)>BackTrack_CHI2MAX+2*fabs(t.Haux(5)) &&
		// High chi2 track: could only make it through here because it
		// has large #hits, and hence could be a track w/ SIs. SIs...
		//  ...Are difficult to fit (misalignment),
		//  ...Are re-evaluated later on (in "TrackFit2"),
		//  ...Invest the track w/ high sensitivity to vertexing, which
		//   will then help answer such question as: exclusive event? or
		//   track from an accidentally coincident interaction.
		nSIs<5) // => Let's be tolerant w/ SI tracks.
	      goodTrk = false;
	  }
#define BackTrack_SCIFIS_CHI2MAX 3
#define BackTrack_SCIFIS_NMAX 2
	  if (!goodTrk ||
	      nScifis>=BackTrack_SCIFIS_NMAX && t.NHits-nScifis<10 &&
	      t.Chi2aux/(t.NDics-5)>BackTrack_SCIFIS_CHI2MAX)
	    listTrack.erase(it++);
	  else it++;
	}
	else {
	  if (t.Type==0x2 && t.Hfirst(0)<900 &&  // Candidates to "BackTrackZ1"
	      (t.Scifi==0x2 ||
	       t.Chi2tot/(t.NDFs-4)<TOpt::dCut[9])) {
	    // (In view of skipping it, if #tracks turns out to be too high:
	    // no momentum required => this # can bee very large indeed
	    nBTZ1s++;
	  }
	  it++;
	}
      }
#ifdef DEBUG
      if (idebug>1) DumpHits(); 
#endif
      //      ***** Re-ASSESS THits' STATUS IN ZONE 0x1 *****

      //  Set "THit::Status=0" for those THits free to be used, i.e. not
      // associated to a TTrack, directly or as the mirror of an associated
      // counterpart.
      //  "Status" collects aspects of THits. Primarily, whether they are used
      // or not. This is set in "PrePattern2" and has not been updated (at the
      // difference to "THit::sTrackID") after several tracks have been erased
      // since Secondarily, some other aspects are assigned by "ImportClusters".
      // But they are no longer relevant, now that "PrePattern2" was completed,
      // and we feel free to overwrite them. Except coalescence related
      // "Status==-4", which must have prevented the THit from being associated
      // and must have remained unchanged throughtout PR, and we want to keep it
      // as is. (Note: Another still relevant aspect may be the GEM spacer flag.
      // But cannot be taken into account here, since no special "Status" value
      // singles it out.)
      int ipli = setup.vIplFirst()[0], iplf = setup.vIplLast()[0];
      for (int ipl = ipli; ipl<=iplf; ipl++) {
	const TPlane &p = setup.vPlane(ipl); if (!p.IFlag/* Is OFF */) continue;
	vector<int>::const_iterator ihit;
	for (ihit = p.vHitRef().begin(); ihit!=p.vHitRef().end(); ihit++) {
	  THit &h = vecHit[*ihit]; int status = h.sTrackID().empty() ? 0 : 1;
	  int mirror;	    // Case of DC: take mirror's status into account...
	  if (h.Status==-4) // ...Exclude silent component of coalesced cluster
	    status = -4;
	  else if ((mirror = h.Mirror)!=-1 &&
		   !vecHit[mirror].sTrackID().empty()) status = 1;
	  h.status = status;
	}
      }

      if (TOpt::ReMode[18]&0x1)
	BackTrackSAS();        // ***** BACKTRACK SAS, w/ momentum, TRACKS *****
      if (TOpt::ReMode[18]&0x2) {
#define BackTrack_0x2TRACKS_NMAX 10
	if (nBTZ1s<=BackTrack_0x2TRACKS_NMAX)
	  BackTrackZ1(0);      // ***** BACKTRACK 0x2-TRACKS to zone 0x1 *****
	else {
	  CsErrLog::msg(elError,__FILE__,__LINE__,
	    "# of 0x2 Tracks (=%d) > %d => Restrict BackTracking to zone 0x1",
			nBTZ1s,BackTrack_0x2TRACKS_NMAX);
	  BackTrackZ1(1); // BackTrack only ``very good'' tracks to zone 0x1 
	}
      }
    }
  }

  // ******************** MAKE BACKWARD/FORWARD GLOBAL FIT ********************

  Traffic::Ref().Stopwatch.Start(13);
  it = listTrack.begin();
  while (it!=listTrack.end()) {
    // All track must have been KF backward fitted by now. Except cases when
    // track has been updated (e.g. cleaning) after fit was performed
    TTrack &t = *it;
    if ((t.Type==0x3 || t.Type==0x6 || t.Type==0x7) &&   // Bridged track...
	(!(t.IFit&0x22) ||              // ...not yet KF fitted
	 t.Chi2tot/(t.NDFs-6)>TOpt::dCut[9]/10)) {       // ...or bad chi2...
      if (t.IFit&0x8) t.Hlast(5) = t.Haux(5);
      t.Hlast(5,5) = 1.E-4;
      bool ret; if (TOpt::ReMode[46]) {
	t.InsertMissedPlanes(); ret = t.FullKF(-1);    t.IFit = 0x20;
      }
      else {                    ret = t.QuickKF(-1,1); t.IFit = 0x2; }
      if (!ret) {// Best would be to recover the pieces (i.e. 1-zone segments)
	// before discarding the track. But, if one agrees that the overall
	// tracks are, at the present stage, fairly reliable (an assumption that
	// may need a more scrutiny), that "QuickKF" fails is very improbable.
	// => Therefore, we don't loose much by erasing the track altogether.  
	listTrack.erase(it++); continue;
      }
    }
    if (t.Type==0x3 || t.Type==0x6 || t.Type==0x7) {    // Bridged track...
      if (t.Chi2tot/(t.NDFs-6)>TOpt::dCut[9]/10 ||      // ...w/ bad chi2...
	  t.Hfirst(5)<5) {                              // ...or low momentum...
	if (it->IFit&0x8) {
	  t.Hfirst(5) = t.Haux(5); t.Hlast(5,5) = 1.E-4;
	}
	bool ret; if (TOpt::ReMode[46]) {
	  ret = t.FullKF(1);    t.IFit |= 0x40;
	}
	else {
	  ret = t.QuickKF(1,1); t.IFit |= 0x4;
	}
	if (!ret) {
	  listTrack.erase(it++); continue;
	}
      }
    }
    //#define DISPLAY_SINGLE_SEGMENTS
#ifdef DISPLAY_SINGLE_SEGMENTS
    if (t.IFit&0x08) it->Hfirst = it->Haux;
#else
    if ((t.Type==0x3 || t.Type==0x6 || t.Type==0x7) &&
	!it->Hfirst.with_mom()) it->Hfirst = it->Haux;
#endif
    it++;
  }   // end of track loop

  if (setup.vIplFirst().size()>2) { // Implies == 5 or 6, cf. "TSetup::Init".
    // ******************** BRIDGING OVER MUON WALL#2 ********************
    //  The bridging over the muFilter (in order to connect 0x4 and 0x8 track
    // segments) is based on hinging points, at which 0x4 and 0x8 extrapolations
    // are required to meet. (Note: The definition of the muFilter here meant is
    // an extended one, including H/ECAL2. The simple hinging point connection
    // is followed, in "TracksFit2", by a more rigorous fitting, using the
    // material maps.)
    //  The hinging points are derived from the material maps. There are 4,
    // corresponding to 4 different regions defined in the yz plane.
    //  The inner region is delimited by the central hole of HCAL2
    double xHC2  = TOpt::Calo[10];
    double y1Low = TOpt::Calo[16]-TOpt::Calo[18], y1Up = y1Low+2*TOpt::Calo[18];
    double z1Low = TOpt::Calo[17]-TOpt::Calo[19], z1Up = z1Low+2*TOpt::Calo[19];
    //  Then next region is that delimited by the central hole in MF2
    double xMF2 =  TOpt::MuonWall[10];
    double y2Low = TOpt::MuonWall[16]-TOpt::MuonWall[18], y2Up = y2Low+2*TOpt::MuonWall[18];
    double z2Low = TOpt::MuonWall[17]-TOpt::MuonWall[19], z2Up = z2Low+2*TOpt::MuonWall[19];
    //  The next 2 regions are delimited by the outer reaches of ECAL2
    double xEC2  = TOpt::Calo[30];
    double y3Low = TOpt::Calo[31]-TOpt::Calo[34], y3Up = y3Low+2*TOpt::Calo[34];
    double z3Low = TOpt::Calo[31]-TOpt::Calo[35], z3Up = z3Low+2*TOpt::Calo[35];
    TTrack tr, tl, tx; list<TTrack>::iterator il, ir;
    for (il = listTrack.begin(); il!=listTrack.end(); il++) {

      //        ***** LOOP OVER LEFT TRACK PIECES *****
      if (!((*il).Type&0x4)) continue;
      tl = *il; // Copy left track
      // Extrapolate to H/ECAL2 to determine in which region the bridging is
      // supposed to take place. And assign corresponding hinging point to "x0".
      double dx = xHC2-tl.Hlast(0);
      double yl =  tl.Hlast(1)+dx*tl.Hlast(3), zl =  tl.Hlast(2)+dx*tl.Hlast(4);
      dx = xMF2-tl.Hlast(0);
      double yl2 = tl.Hlast(1)+dx*tl.Hlast(3), zl2 = tl.Hlast(2)+dx*tl.Hlast(4);
      static double x0; if (y1Low<yl && yl<y1Up && z1Low<zl && zl<z1Up)
	x0 = setup.MuFilterHinges[0];
      else if (y2Low<yl2 && yl2<y2Up && z2Low<zl2 && zl2<z2Up)
	x0 = setup.MuFilterHinges[1];
      else if (xEC2) {
	dx = xEC2-tl.Hlast(0);
	yl = tl.Hlast(1)+dx*tl.Hlast(3); zl = tl.Hlast(2)+dx*tl.Hlast(4);
	if (y3Low<yl && yl<y3Up && z3Low<zl && zl<z3Up)
	  x0 = setup.MuFilterHinges[2];
	else
	  x0 = setup.MuFilterHinges[3];
      }
      else
	x0 = setup.MuFilterHinges[2];

      for (ir = listTrack.begin(); ir!=listTrack.end(); ir++) {

	//        ***** LOOP OVER RIGHT TRACK PIECES *****
	if ((*ir).Type!=0x8)   continue;
	tr = *ir; // Copy right track

	if (TOpt::ReMode[6]>0) { // Ideal bridging
	  // - Require well reco'd tracks, w/ same associated MC.
	  // - Allow still for a pi or K to decay to a muon.
	  bool good_pair = false;
	  if (tl.IKine>=0 && tr.IKine>=0) {
	    good_pair = tl.IKine==tr.IKine; int iv;
	    if (!good_pair && (iv = vecKine[tr.IKine].IVtxOrig)>=0 &&
		vecVtxMC[iv].IOrig==tl.IKine) {
	      int idl = vecKine[tl.IKine].PDGid(),
		idr = vecKine[tr.IKine].PDGid();
	      if ((abs(idl)==211/* pi */ || abs(idl)==321/* K */) &&
		  abs(idr)==13/* mu */ && idl*idr>0) good_pair = true;
	    }
	  }
	  if (!good_pair) continue;
	}

	// Track par. comparisson
	double dazi= tr.Hfirst.azi() - tl.Hlast.azi();
	double ddip= tr.Hfirst.dip() - tl.Hlast.dip();
	double dz  = (tr.Hfirst(2) + (x0-tr.Hfirst(0))*tan(tr.Hfirst.dip())) 
	          -  (tl.Hlast (2) + (x0-tl.Hlast (0))*tan(tl.Hlast.dip()));
	double dy  = tr.Hfirst(1) + (x0-tr.Hfirst(0))*tr.Hfirst(3) 
	          - (tl.Hlast (1) + (x0-tl.Hlast (0))*tl.Hlast (3));
	static float sigmady = 0.3, sigmadz = 1,
	  sigmadazi = 500*1.5, sigmaddip = 500*2.5;


	float chi2 = sqrt(ddip*ddip/sigmaddip/sigmaddip+
			  dazi*dazi/sigmadazi/sigmadazi+
			  dz*dz/sigmadz/sigmadz+
			  dy*dy/sigmady/sigmady)/4;
	if (hist) {
	  muW_dy->Fill(dy,tl.Pinv());
	  muW_dz->Fill(dz,tl.Pinv());
	  muW_dazi->Fill(dazi*500,tl.Pinv());
	  muW_ddip->Fill(ddip*500,tl.Pinv());
	}
	if (chi2>6) continue;
	if (tl.Type==0x4 &&      // Left side track is momentumless and...
	    tl.Hfirst(0)>3000)   // ...and has 1st point @ DW04 or beyond
	  // => It's most probably a short segment that the reco failed to
	  //   recognise as a part of a longer track.
	  // => Try to have the longer track win over it in any case.
	  chi2 += 3;         //  => Give it a malus.
	tp.Chi2 = chi2; tp.iL = il; tp.iR = ir; tp.Hfirst = tl.Hfirst;
	lTrackPair3.push_back(tp);

      }  // End of loop over right track pieces
    }  // End of loop over left  track pieces
    tr.Id = ++TTrack::TrackCounter;// make IDs of temporary tracks unique
    tl.Id = ++TTrack::TrackCounter;// before destructors is called
    tx.Id = ++TTrack::TrackCounter;// (to prevent removal of it's IDs from sTrackIDs of TKine)

    //      ********** SORT and BUILD ALL-OUT GLOBAL TRACKS **********

    lTrackPair3.sort();

    list<TTrackPair2>::iterator itp3;
    for (itp3 = lTrackPair3.begin(); itp3!=lTrackPair3.end(); itp3++) {

      //        ***** LOOP OVER MUON WALL#2 TRACK PAIRS

      TTrackPair2 &tp = *itp3; TTrack &tl = *(tp.iL), &tr = *(tp.iR);
      if (tl.IMark==-1 ||        // NEVER occurs
	  (tl.Type&0x8))         // Already extended
	continue;
      if (tr.Type!=0x8 ||        // NEVER occurs
	  tr.IMark==-1)          // Already appended
	continue;
      tl.Haux = tl.Hlast;        // Save Helix before wall
      tl.Append(tr);             // Append 0x10 segment to track
      tl.Hlast(5)   = tl.Haux(5);
      tl.Hlast(5,5) = tl.Haux(5,5);
      tl.Associate = tr.Id;
    }
    // Let's not erase tracks pieces left over behind wall now. Keep them
    // instead, in order to:
    //  - have them at hand when it comes to (in "TrackFit2") evaluate the chi2
    //   increment corresponding to the extension downstream of the wall,
    //  - erase them only when that increment turns out to be reasonably small.
  }

  Traffic::Ref().Stopwatch.Stop(13);

}
// ~~~~~~~~~~~~~~~~~~~~~~~~ BridgeOverTarget ~~~~~~~~~~~~~~~~~~~~~~~~~
void TEv::BridgeOverTarget(list<TTrackPair2> &lTrackPair0, TTrack &t, int nDFs0x1, double chi20x1)
{
  // ********** CONTINUATION BEFORE TARGET FOR 0x3 or 0X7 TRACKS **********

  if (TOpt::iCut[15]*TOpt::dCut[4]*TOpt::dCut[5]==0)// Beam q, P, dP unspecifed
    return; // => Cannot evaluate momentum compatibility => Give up
  double cop1 = t.Hfirst(5), d2cop1 = t.IFit!=0x8 ? t.Hfirst(5,5) : 1e-4;
  double copb = TOpt::iCut[15]/TOpt::dCut[4], d2copb = TOpt::dCut[5]*TOpt::dCut[5];
  double cop1b = fabs(cop1-copb), dcop = sqrt(d2cop1+d2copb);
  //#define BridgeS_DEBUG_BoT
#ifdef BridgeS_DEBUG_BoT
  printf("BoT: B %9.2f+/-%8.2f, Trk(0x%x) %9.2f+/-%8.2f => %8.1f %s %8.1f\n",
	 1/copb,sqrt(d2copb)/copb/copb,t.Type,1/cop1,sqrt(d2cop1)/cop1/cop1,
	 fabs(1/cop1-1/copb),cop1b>5*dcop?">5*":"<5*",sqrt(d2cop1/cop1/cop1/cop1/cop1+d2copb/copb/copb/copb/copb));
  printf("BoT: B %9.6f+/-%8.6f, Tkk(0x%x) %9.6f+/-%8.6f => %8.6f %s %8.6f\n",
	 copb,sqrt(d2copb),            t.Type,cop1,sqrt(d2cop1),
	 cop1b,              cop1b>5*dcop?">5*":"<5*",dcop);
  static int bdebug = 1;
  if (bdebug)
#endif
  if (cop1b>5*dcop) return;

  list<TTrackPair2>::iterator itp0;
  for (itp0 = lTrackPair0.begin(); itp0!=lTrackPair0.end(); itp0++) {
    //                    ***** LOOP OVER TARGET TRACK PAIRS *****
    TTrackPair2 &tp0 = *itp0; TTrack &tl0 = *tp0.iL, &tr0 = *tp0.iR;
    if (tr0.Id!=t.Id) continue;         // Must have common track type 0x1
    if (tl0.Type&0x1) continue;         // Already extended
    if (TOpt::ReMode[26]&0x10) {  // 0x10: Beam continuation beyond SM1...
      // ...The idea is to flag cases where the spectrometer track is a
      // mere continuation of a, non-interacting, incident track.
      //  In doing so, one need take care to preserve low Q^2 events,
      // which can easily be confused w/ the cases supra. A condition on
      // Q^2 being impossible at the present stage of the reco, one can
      // only require good chi2 and reasonably large 0x3/7 momentum. Good
      // chi2 is requested from the 0x11 fit, both in absolute value
      // and in increment. For the increment, we decide to check it on the
      // total, as opposed to reduced, chi2, w/o any real evaluation of
      // the pros and cons of the 2 solutions.
      //  The reliability of the continuation hypothesis will be checked
      // again, in "TracksFit2" and later in the vertex package. Here we
      // do it only to minimise CPU consumption.
      if (tl0.Associate>=0) continue;     // Already has an associate
      double chi20x11 = tp0.Chi2tot/(tl0.NDFs+nDFs0x1-4);
      double sChi2 = tl0.Chi2tot+chi20x1*(nDFs0x1-4);
      bool ok = chi20x11<4 && tp0.Chi2tot<5*sChi2 ||
	// Allow some of the 0x11 chi2 to arise from either 0x1 or 0x10.
	chi20x11<5.5 && tp0.Chi2tot<2*sChi2;
      if (ok) {// Require p+/-5sigma. 5 sigmas: to account for the fact...
	// ...that p is retrieved from the crude "QuickKF".
	double p = 1/t.Hlast(5);
	double dr = (t.IFit&0x8) ? p/30: sqrt(t.Hlast(5,5)), dp = 5*dr*p*p;
	double nSigmas = TOpt::iCut[30]>0 /* hadron beam */ ? 3 : sqrt(12.);
	double p0 = TOpt::dCut[4]*TOpt::iCut[15],
	  dp0 = nSigmas*TOpt::dCut[5]*p0*p0;
	if (p0-dp0<p+dp && p-dp<p0+dp0) tl0.Associate = t.Id;
      }
    }
    else {   
      tl0.Append(t);                 // Prepend track type 0x10
      if (tp0.WorstPl>=0 && // Have also to check that this very "WorstPl" has
	  // not already been removed (because specified in TTrackPair2 "tp1").
	  // Since "tp1" is not accessible here, let's instead check that
	  // "WorstPl" belong to the 0x10 zone.
	  tp0.WorstPl<TSetup::Ref().vIplFirst()[0]) tl0.SubHit(tp0.WorstPl);
      tl0.Hlast(5)   = t.Hlast(5);   // Store mom @ last point
      tl0.Hlast(5,5) = t.Hlast(5,5);
      tl0.IFit = 0; // No backward fit yet both based on 0x3 (or 0x7) momentum and including the 0x10 segment
    }
    break;
  }
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ IsYokeTrk ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool IsYokeTrk(const TTrack &tl, const TTrack &tr)
{
  //                 ***** SM2 YOKE TRACKs *****
  const TSetup& setup = TSetup::Ref();
  // xSM2I/O = Abscissa of entrance and exit.
  // y/zSM2Low/up = Dimensions of SM2 gap.
  double xSM2I =   setup.MagCenter(2,0)-200, xSM2O =  setup.MagCenter(2,0)+200;
  double ySM2Low = setup.MagCenter(2,1)-100, ySM2Up = setup.MagCenter(2,1)+100;
  double zSM2Low = setup.MagCenter(2,2)- 50, zSM2Up = setup.MagCenter(2,2)+ 50;
  const THlx &hFirst = tr.H('d'), &hLast = tl.H('u');
  double ySM2I = hLast(1)+(xSM2I-hLast(0))*hLast(3);
  if (ySM2I<ySM2Low || ySM2Up<ySM2I) return true;
  else {
    double ySM2O = hFirst(1)+(xSM2O-hFirst(0))*hFirst(3);
    if (ySM2O<ySM2Low || ySM2Up<ySM2O) return true;
    else {
      double zSM2I = hLast(2) +(xSM2I-hLast(0))*hLast(4);
      if (zSM2I<zSM2Low || zSM2Up<zSM2I) return true;
      else {
	double zSM2O = hFirst(2)+(xSM2O-hFirst(0))*hFirst(4);
	if (zSM2O<zSM2Low || zSM2Up<zSM2O) return true;
      }
    }
  }
  return false;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DoubleBridge ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#include "CsPixelGEMDetector.h"
/*! 
  \brief Bridge a 0x1 and a 0x4 segments over SM1 and SM2 and pick up hits in zone 0x2. 
  \author  Yann.Bedfer@cern.ch

  - The routine is meant to supplement the bad redundancy in the VSAT region of
   zone 0x2. Hence, the 0x1 segment is required to be an (almost)VSAT-only one.
  - There's some provision for the 0x4 segment to have been bridged over the
   mu-wall (yielding a 0xc segment).
  - Since the routine is (expected to be) CPU intensive and intended to only
   correct for what remains an priori small deficiency of the standard tracking
   algorithm, many precautions are taken to limit its use to only useful cases:
   - Severe selection of input 0x1 and 0x4 segments.
   - Evaluation (w/ full evaluation meaning full KF fit and pick-up in addition
   to quick KF fit) of the sole best candidate pairs.

  - NOTA BENE: The source code relies on a number of assumption concerning
   the experimental setup. In particular, the ordering of the coordinate planes
   in a pixelGEM station is assumed to be U<V<Y<X.
*/

void TEv::DoubleBridge()
{
  static TH2D *hDBdz, *hDBddip;
  static int hist = 0; if (!hist) {  //            ***** BOOK HISTOS *****
    if (TOpt::Hist[12]) {
      hist = 1;
      if (!gDirectory->cd("/Traffic")) {
	TDirectory *traffic = gDirectory->mkdir("/Traffic","TraFFiC");
	traffic->cd();
      }
      hDBdz   = new TH2D("hDBdz",  "DbleBridge: dz - 0/1 = good/bad pairs",
		       200,-2,2,2,-.5,1.5);
      hDBddip = new TH2D("hDBddip","DbleBridge: ddip - 0/1 = good/bad pairs",
		       200,-2,2,2,-.5,1.5);
    }
    else hist = -1;
  }

  //    ***** CHECK NOT TOO MANY 0x1 AND 0X4|c TRACK SEGMENTS *****
  // - The double bridging is based on track extrapolation, hence CPU intensive.
  //  We want to avoid wasting CPU on shower events, which easily produce
  //  orphan track segments.
  list<TTrack>::iterator it, jt; int nTs0x1, nTs0x4;
  for (it = listTrack.begin(), nTs0x1=nTs0x4 = 0; it!=listTrack.end(); it++) {
    const TTrack &t = *it;
    if (t.Type==0x1) nTs0x1++;
    if ((t.Type&0x7)==0x4 /* Allow for 0xc */) nTs0x4++;
  }
  // Concerning the cut on the max. #segments: What we have in mind in
  // particular are Primakoff events (where the 0x2 zone, w/ its low redundancy,
  // would be missing because of inefficiency). Therefore we start with the very
  // simple event topology: 1 0x1 + 1 0x4. To which we have to add some margin
  // to make for possible accidentally coincident beam particles (let's say 2,
  // at most) XOR background reaction (e.g. pi0 conversion) => Let's cut @ 9.
  if (nTs0x1*nTs0x4>9 || nTs0x1*nTs0x4==0) return;

  //            ***** INITIALIZATION *****
  const TSetup &setup = TSetup::Ref();
  double xSM1 = setup.MagCenter(1,0), xSM2 = setup.MagCenter(2,0);
  double xSM12 = (xSM1+xSM2)/2;
  double fldint = setup.MagFieldInt[1]+setup.MagFieldInt[2];
  // Indices of Y/Z and U/V = +/-45deg. proj. This is used infra in the
  // book-keeping of proj. , when dealing w/ G|MP which are expected to be at
  // either of 0, 90, +/45 deg., cf. "TSetup::Init".
  static unsigned int yzProj = 3 /* Cf. "TSetup::Init" */, uvProj = 0;
  // Amplitude correlation of GPXYUVs
  static map<int,int> ids;
  static vector<float> aYX0s, aYX1s, aYX2s, aUV0s, aUV1s, aUV2s;

  static bool first = true; if (first) {
    first = false;
    int iproj, iProjU, iProjY = 1;  // Indices of Y/Z and U/V = +/-45deg.
    for (int iproj=iProjU = 0; iproj<(int)setup.vProj().size(); iproj++) {
      unsigned int proj = 1<<iproj; double alpha = setup.vProj()[iproj]/10.;
      if      (fabs(alpha-45)<2) { uvProj |= proj; iProjU = iproj; }
      else if (fabs(alpha+45)<2)   uvProj |= proj;
    }
    // Amplitude correlation of GPXYUVs:
    // - Store the coefficients of the polynomial parameterisations aU(aV) and
    //  aY(aX) in static vectors.
    // - Map the vector index on the pixelGEM station index.
    // - The pixelGEM coordinate planes are expected to be ordered U<V<Y<X.
    //  Standard stations have all of UVXY. UV or XY can be missing (case of UV
    //  of GP01 in 2008/9 setup). All other combinations are not provided for.
    const vector<TStation> &stations = setup.vStation();
    for (int is = 0; is<(int)stations.size(); is++) {
      const TStation &s = stations[is]; if (s.Type!=28) continue;
      int id = ids.size(); ids[is] = id;
      int jpl = 0;
      const TPlane &pU = setup.vPlane(s.IPlanes[0]);
      if (pU.IProj==iProjU) {
	const TDetect &dU = setup.vDetect(pU.IDetRef);
	CsPixelGEMDetector *pg = dynamic_cast<CsPixelGEMDetector*>(dU.PtrDet());
	if (!pg)
	  CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "TStation #%d(Type=28)'s TPlane, \"%s\", cannot be cast into a CsPixelGEM",
			is,dU.Name.c_str());
	const float *ampCorr = pg->getAmpCorr();
	aUV0s.push_back(ampCorr[0]); aUV1s.push_back(ampCorr[1]);
	aUV2s.push_back(ampCorr[2]);
	jpl = 2;
      }
      else { // Case: UV is missing (as is e.g. the case, in 2008/9 setup, w/...
	// ...GP01 (which, btw, is not concerning us here, since in zone 0x1))
	aUV0s.push_back(0); aUV1s.push_back(0); aUV2s.push_back(0);
      }
      const TPlane &pY = setup.vPlane(s.IPlanes[jpl]);
      if (pY.IProj==iProjY) {
	const TDetect &dY = setup.vDetect(pY.IDetRef);
	CsPixelGEMDetector *pg = dynamic_cast<CsPixelGEMDetector*>(dY.PtrDet());
	if (!pg)
	  CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "TStation #%d(Type=%d)'s TPlane, \"%s\", cannot be cast into a CsPixelGEM",
			is,s.Type,dY.Name.c_str());
	const float *ampCorr = pg->getAmpCorr();
	aYX0s.push_back(ampCorr[0]); aYX1s.push_back(ampCorr[1]);
	aYX2s.push_back(ampCorr[2]);
      }
      else {
	aYX0s.push_back(0); aYX1s.push_back(0); aYX2s.push_back(0);
      }
    }
  }

  TTrack tx;
  TTrackPair2 tp;
  list<TTrackPair2> lTrackPair; list<TTrackPair2>::iterator itp;

  for (it = listTrack.begin(); it!=listTrack.end(); it++) {
    TTrack &ti = *it;
    if ((ti.Type&0x7)!=0x4) continue;            // ***** LOOP ON 0x4|c SEGMENTS
    //                                                ***** 0X4|c SEGMENT FILTER
    // - Require the segment to fulfill PR cut on #hits, meaning it has not lost
    //  some it hits in TEv::CleanTrackList, and hence has not had to compete w/
    //  some other neighbouring track. So as to minimize the risk of a fake
    //  bridging w/ a 0x1 segment that would belong to this other track.
    int nHitsMn = TOpt::iPRpar[25]; double chi2Mx = TOpt::dCut[87];
    if (ti.Type==0xc) {
      nHitsMn += TOpt::iPRpar[25];
      // chi2/NDF: the case of 0xc track is problematic, because such a track
      // cannot but have had its momentum assigned peremptorily = beam momentum.
      // Which leads to underestimating the MC it undergoes in the mu-wall.
      chi2Mx *= 1.2;      // => Give it some margin
    }
    if ((int)ti.NDFs<nHitsMn || ti.Chi2tot/(ti.NDFs-4)>chi2Mx) continue;
    // - Track segment could eventually be erased (if indeed double bridged).
    //  Which can pose a problem when it's associated
    //   Association can be flagged by TTrack::Associate or ::IMark. For the
    //  latter, require that it be strictly =0, in order to avoid overwriting
    //  any other special setting.
    if (ti.Associate>=0 || ti.IMark!=0) continue;
    // - Require a min. # of space points and projections.
    list<int>::iterator ih; const TStation *sPrv;
    int nPls, mPls; double data[4];
    int nSpacePts, nProjs, nAllProjs; unsigned int projs, allProjs;
    for (ih = ti.lHitPat.begin(), nSpacePts=nProjs=nAllProjs=nPls=mPls = 0,
	   projs=allProjs = 0, sPrv = 0; ih!=ti.lHitPat.end(); ih++) {
      if (*ih<0) continue;
      const THit &h = vecHit[*ih]; const TPlane &p = setup.vPlane(h.IPlane);
      const TStation *&s = p.Station; if (s!=sPrv) {
	if (sPrv) {
	  if      (sPrv->Type==29) {             // Pixels (GP..P) station...
	    // ...One SpacePt per pixel (note: one pixel plane <-> one station)
	    if (nProjs>=2 /* It's either 0 or 2 */) {
	      nSpacePts++; nProjs = 0; projs = 0;
	    }
	  }
	  else if (sPrv->Type==28) {             // Strips (GP..XYUV) station...
	    // Let's try to build 1 SpacePt per detector, by relying on the
	    // amplitude correlation. This in order to be able qualify tracks
	    // such as GP03X+Y+U+V + FI08X+Y in the 2009 setup.
	    int id = ids[sPrv->IStation];
	    if (nPls==2) {
	      double aU = data[0], aV = data[1];
	      double deltaUV = aUV0s[id]+aUV1s[id]*aV+aUV2s[id]*aV*aV-aU;
	      if (fabs(deltaUV)<TOpt::ReMode[17]) {
		nSpacePts++; nProjs = 0; projs = 0;
	      }
	    }
	    if (mPls==2) {
	      double aY = data[2], aX = data[3];
	      double deltaYX = aYX0s[id]+aYX1s[id]*aX+aYX2s[id]*aX*aX-aY;
	      if (fabs(deltaYX)<TOpt::ReMode[17]) {
		nSpacePts++; nProjs = 0; projs = 0;
	      }
	    }
	    nPls=mPls = 0;
	  }
	  else if (sPrv->Type==22) {             // Scifis station...
	    // Let's try to build 2 space points per station, by relying on the
	    // time correlation. This in order to be able qualify tracks such as
	    // GP03X+Y+U+V + FI08X+Y in the 2009 setup.
	    if (nPls==2 && fabs(data[0]-data[1])<2) {
	      nSpacePts++; nProjs = 0; projs = 0;
	    }
	    nPls = 0;
	  }
	  if (nProjs>=3) { nSpacePts++; nProjs = 0; projs = 0; }
	}
	sPrv = s;
      }
      unsigned int proj = 1<<p.IProj;
      if      (s->Type==22) data[nPls++] = h.Time;
      else if (s->Type==28) {
	if (proj&uvProj) data[nPls++] = h.ptrCl->getAllAnalogData()[2];
	else             data[2+mPls++] = h.ptrCl->getAllAnalogData()[2];
      }
      if ((proj&projs)==0)        {    projs |= proj;    nProjs++; }
      if ((proj&allProjs)==0)     { allProjs |= proj; nAllProjs++; }
      if (p.IFlag&0x30) { // P-pixels: provide also v-axis.
	if   (p.IFlag&0x10) proj ^= yzProj;
	else                proj ^= uvProj;
	if ((proj&projs)==0)      {    projs |= proj;    nProjs++; }
	if ((proj&allProjs)==0)   { allProjs |= proj; nAllProjs++; }
      }
    }
    if (sPrv) { // Account for last TStation
      if      (sPrv->Type==29) {              // Pixels (GP..P) station
	if (nProjs>=2) { nSpacePts++; nProjs = 0; }
      }
      else if (sPrv->Type==28) {              // Strips (GP..XYUV) station
	int id = ids[sPrv->IStation];
	if (nPls==2) {
	  double aU = data[0], aV = data[1];
	  double deltaUV = aUV0s[id]+aUV1s[id]*aV+aUV2s[id]*aV*aV-aU;
	  if (fabs(deltaUV)<TOpt::ReMode[17]) { nSpacePts++; nProjs = 0; }
	}
	if (mPls==2) {
	  double aY = data[2], aX = data[3];
	  double deltaYX = aYX0s[id]+aYX1s[id]*aX+aYX2s[id]*aX*aX-aY;
	  if (fabs(deltaYX)<TOpt::ReMode[17]) { nSpacePts++; nProjs = 0; }
	}
      }
      else if (sPrv->Type==22) {              // Scifis station
	if (nPls==2 && fabs(data[0]-data[1])<2) { nSpacePts++; nProjs = 0; }
      }
      if (nProjs>=3) nSpacePts++;
    }
    if (nAllProjs<3 || nSpacePts<3) continue;
    for (jt = listTrack.begin(); jt!=listTrack.end(); jt++) {
      TTrack &tj = *jt;
      // MC case
      if (tj.Type!=0x1) continue;                  // ***** LOOP ON 0x1 SEGMENTS
      //                                                ***** 0x1 SEGMENT FILTER
      bool bonnePaire = ti.IKine>=0 && tj.IKine>=0 && tj.IKine==ti.IKine;
      // - Require (almost)VSAT-only. Otherwise, the particle avoids the very
      //  beam region in zone 0x2, which is where the redundancy is low and
      //  prevents the 0x2 segment of the track from being reconstructed.
      if (!(tj.Scifi&0x101)) continue;
      // - Do not require the segment to fulfill PR cut on #hits, since we're
      //  in the VSAT region. Yet, to be on the safe side, ask for some minimal
      //  #hits.
      //   Although chi2 is most of the time not much constraining for VSAT
      //  track, still cut on chi2/NDF,so as to also account for setups w/ SIs.
      if ((int)ti.NDFs<5 /* Out of a 6 in e.g. the [2002,2007] setup */ ||
	  tj.Chi2tot/(ti.NDFs-4)>TOpt::dCut[87])
      // - 0x1 segments cannot be associated. Yet, to on the safe side...
      if (tj.Associate>=0 || tj.IMark!=0) continue;
      // - Require a, loose, min. # of space points and projections.
      //  (Cf. supra for detailed comments.)
      for (ih = tj.lHitPat.begin(), nSpacePts=nProjs=nAllProjs=nPls=mPls = 0,
	     projs=allProjs = 0, sPrv = 0; ih!=tj.lHitPat.end(); ih++) {
	if (*ih<0) continue;
	const THit &h = vecHit[*ih]; const TPlane &p = setup.vPlane(h.IPlane);
	const TStation *&s = p.Station; if (s!=sPrv) {
	  if (sPrv) {
	    if      (sPrv->Type==29) {                  // Pixel GP
	      if (nProjs>=2) { nSpacePts++; nProjs = 0; projs = 0; }
	    }
	    else if (sPrv->Type==28) {                  // Strip GP
	      int id = ids[sPrv->IStation];
	      if (nPls==2) {
		double aU = data[0], aV = data[1];
		double deltaUV = aUV0s[id]+aUV1s[id]*aV+aUV2s[id]*aV*aV-aU;
		if (fabs(deltaUV)<TOpt::ReMode[17]) {
		  nSpacePts++; nProjs = 0; projs = 0;
		}
	      }
	      if (mPls==2) {
		double aY = data[2], aX = data[3];
		double deltaYX = aYX0s[id]+aYX1s[id]*aX+aYX2s[id]*aX*aX-aY;
		if (fabs(deltaYX)<TOpt::ReMode[17]) {
		  nSpacePts++; nProjs = 0; projs = 0;
		}
	      }
	      nPls=mPls = 0;
	    }
	    else if (sPrv->Type==22) {                  // Scifis station
	      if (nPls==2 && fabs(data[0]-data[1])<2) {
		nSpacePts++; nProjs = 0; projs = 0;
	      }
	      nPls = 0;
	    }
	    if (nProjs>=3) { nSpacePts++; nProjs = 0; projs = 0; }
	  }
	  sPrv = s;
	}
	unsigned int proj = 1<<p.IProj;
	if      (s->Type==22) data[nPls++] = h.Time;
	else if (s->Type==28) {
	  if (proj&uvProj) data[nPls++] = h.ptrCl->getAllAnalogData()[2];
	  else             data[2+mPls++] = h.ptrCl->getAllAnalogData()[2];
	}
	if ((proj&projs)==0)        {    projs |= proj;    nProjs++; }
	if ((proj&allProjs)==0)     { allProjs |= proj; nAllProjs++; }
	if (p.IFlag&0x30) { // P-pixels: provide also v-axis.
	  if   (p.IFlag&0x10) proj ^= yzProj;
	  else                proj ^= uvProj;
	  if ((proj&projs)==0)      {    projs |= proj;    nProjs++; }
	  if ((proj&allProjs)==0)   { allProjs |= proj; nAllProjs++; }
	}
      }
      if (sPrv) { // Account for last TStation
	if      (sPrv->Type==29) {                    // Pixel GP
	  if (nProjs>=2) { nSpacePts++; nProjs = 0; }
	}
	else if (sPrv->Type==28) {                    // Strip GP
	  int id = ids[sPrv->IStation];
	  if (nPls==2) {
	    double aU = data[0], aV = data[1];
	    double deltaUV = aUV0s[id]+aUV1s[id]*aV+aUV2s[id]*aV*aV-aU;
	    if (fabs(deltaUV)<TOpt::ReMode[17]) { nSpacePts++; nProjs = 0; }
	  }
	  if (mPls==2) {
	    double aY = data[2], aX = data[3];
	    double deltaYX = aYX0s[id]+aYX1s[id]*aX+aYX2s[id]*aX*aX-aY;
	    if (fabs(deltaYX)<TOpt::ReMode[17]) { nSpacePts++; nProjs = 0; }
	  }
	}
	else if (sPrv->Type==22) {                    // Scifis station
	  if (nPls==2 && fabs(data[0]-data[1])<2) { nSpacePts++; nProjs = 0; }
	}
	if (nProjs>=3) nSpacePts++;
      }
      if (nAllProjs<3 || nSpacePts<3) continue;
      // MC case
      int goodPair = ti.IKine>=0 && tj.IKine>=0 && tj.IKine==ti.IKine ? 1 : 0;
      if (TOpt::ReMode[6]>0 && !goodPair) continue; // "Ideal" bridging...
      if (TOpt::ReMode[6]==0) {                     // ...else
	//                          ***** CUT ON COMPATIBILITY of TRACK SEGMENTS
	double dazi= tj.Hfirst.azi() - ti.Hlast.azi();
	double ddip= tj.Hfirst.dip() - ti.Hlast.dip();
	double dz = tj.Hfirst(2) + (xSM12-tj.Hfirst(0))*tan(tj.Hfirst.dip())
	  -        (ti.Hlast (2) + (xSM12-ti.Hlast (0))*tan(ti.Hlast.dip()));
	static int idebug = 0;
	if (idebug && isMC) {
	  printf("Tracks #%04d+%04d %s: dz,ddip,dazi = %6.2f,%6.2f,%6.2f\n",
		 ti.Id,tj.Id,goodPair?"OK":"KO",dz,1000*ddip,1000*dazi);
	}
	if (hist>0) {
	  hDBdz->Fill(dz,goodPair); hDBddip->Fill(1000*ddip,goodPair);
	}
	if (fabs(dz)>.5 || fabs(ddip*1000)>.5) continue;
	if (TOpt::dCut[15]!=0 &&            // ***** TIME DIFF CUT... *****
	    (tj.Scifi&0x1) && (ti.Scifi&0x4)) {// ...Scifi-(almost)only segments
	  double dt = (ti.MeanTime-tj.MeanTime)/
	    sqrt(ti.SigmaTime*ti.SigmaTime+tj.SigmaTime*tj.SigmaTime);
	  // Let's make the time cut systematic, contrary to what's done in the
	  // standard "BridgeSegments2" method.
	  if (fabs(dt)>TOpt::dCut[15]) continue; 
	}
      }
      double thiy = atan(ti.Hlast(3)), thjy = atan(tj.Hfirst(3));
      double pinv = 2*cos((thiy+thjy)/2)*sin((thjy-thiy)/2)/.3E-3/fldint;
      pinv *= cos(ti.Hlast.dip());
      //       ********** BUILD MERGED TRACK **********
      tx = tj; tx.Append(ti,"tmp"); ti.IMark = 0; tx.Chi2tot = 1e6;
      //       ********** INITIAL VALUES **********
      double chi2 = 1.E10; tx.Hfirst(5) = pinv; tx.Hfirst(5,5) = 1.E-6;
      int iter = 0, converged = 0; while (iter < 15) { // ***** ITER LOOP *****

	double pinv_prev = tx.Hfirst(5); double chi2_prev = chi2;
	tx.Hlast(5) = tx.Hfirst(5); tx.Hlast(5,5) = tx.Hfirst(5,5);

	if (!tx.QuickKF(-1,1)) break;    // ***** FIT CANDIDATE *****
	chi2 = tx.Chi2tot/tx.NDFs;
	if (chi2>(iter==0 ? 2 : 1)*TOpt::dCut[7]) break;
	//                                        ***** CONVERGENCE TEST *****
	int ichi2      = int(chi2     *10); // Keep only 1 decimal digit
	int ichi2_prev = int(chi2_prev*10);
	if (ichi2==ichi2_prev) { converged = 1; break; }
	  
	iter++;
      }
      if (!converged) continue;
      //       ********** CUT ( FINAL CHI2) **********
      if (TOpt::ReMode[6]==0) {
	pinv = fabs(tx.Hfirst(5));
	if (chi2>TOpt::dCut[9]*(pinv>.2 ? 2 : 1)) // Looser for low momenta
	  continue;
      }
      //       ***** STORE TRACK PAIR CANDIDATES *****
      tp.Chi2 = tx.Chi2tot/(tx.NDFs-6+tx.NHits*tx.NHits*.01); 
      tp.iL = jt; tp.iR = it; tp.Chi2tot = tx.Chi2tot;
      tp.Hfirst = tx.Hfirst;
      tp.IFit = 0; tp.WorstPl = -1;
      tp.NDics = 0; // Not used when IFit!=0x8, but let's assign dummy value to it, to avoid compiler warning
      lTrackPair.push_back(tp);
    }
  }
  tx.Id = ++TTrack::TrackCounter;// Make IDs of temporary track unique before destructor is called (to prevent removal of it's IDs from sTrackIDs of TKine)

  lTrackPair.sort(); //   ********** SORT MERGED PAIRS BY Chi2 **********

  //      ******************** BUILD GLOBAL TRACKS ********************
  for (itp = lTrackPair.begin(); itp!=lTrackPair.end(); itp++) {
    TTrackPair2 &tp = *itp; TTrack &t = *tp.iL, &tr = *tp.iR;
    if (t.IMark<0) continue;    // 0x1 segment already examined
    if (tr.IMark<0) continue;   // 0x4 segment already examined
    TTrack tx(t); tx.Append(tr); tx.Hfirst = tp.Hfirst;
    // We intend to examine only the best pair any given segment belongs to...
    t.IMark=tr.IMark = -2; // ...let's flag both segments as already examined.
    tx.InsertMissedPlanes();
    bool ret = tx.FullKF(1); if (ret) ret = tx.FullKF(-1);
    if (!ret) continue;
    THlx H; double chi2;
    int iplF = setup.vIplFirst()[1], iplL = setup.vIplLast()[1];
    double dTMx = 2; // Timing consistency: assuming 1ns scifi resol, taking 2 sigmas
    tx.UseHitTime();double tST = tx.SigmaTime, dT = sqrt(9*tST*tST+dTMx*dTMx);
    double tMn = t.MeanTime-dT, tMx = tMn+2*dT;
    int expected, added, nPGs; vector<THit*> hitList; 
    list<int>::iterator ip; map<int,THlx>::const_iterator imf, imb;
    for (ip = tx.lPlnRef.begin(), imf = tx.mHef.begin(), imb = tx.mHeb.begin(),
	   expected=added=nPGs = 0; ip!=tx.lPlnRef.end(); ip++, imf++, imb++) {
      int ipl = *ip;
      if (ipl<iplF) continue; if (ipl>iplL) break;
      const TPlane &p = setup.vPlane(ipl); if (p.IFlag==0) continue; // OFF
      const TDetect &d = setup.vDetect(ipl);
      int type = d.IType; // Type of detector...
      bool isScifi = type==22;
      // ...We're interested in (V)SAT: let's skip the rest
      if (!isScifi && type<26 || 32<type /* !G|MM|P */) continue;
      const THlx &Hleft = (*imf).second;
      if (!d.InActive(Hleft(1),Hleft(2))) continue;
      expected++; if (type==28) nPGs++; if (type==29) nPGs += 2;
      bool isPixel = type==29 || type==32; if (isPixel) expected++; 
      const THlx &Hright = (*imb).second;
      Hleft.Update(Hright,H,chi2);
      double y = H(1), z = H(2), dy2 = H(1,1), dyz = H(1,2), dz2 = H(2,2);
      double ca = d.Ca, sa = d.Sa, u = y*ca+z*sa, v = z*ca-y*sa;
      double du = sqrt(dy2*ca*ca+2*dyz*ca*sa+dz2*sa*sa);
      double dv = sqrt(dy2*sa*sa-2*dyz*ca*sa+dz2*ca*ca);
      double cut = 5*du+3*d.Resol, cvt =  5*du+3*d.Resol;
      double umn = u-1.5*cut, umx = u+1.5*cut;
      vector<int>::const_iterator ihit; THit *hBest; double best; int mult;
      for (ihit = p.vHitRef().begin(), best = cut, hBest = 0, mult = 0;
	   ihit!=p.vHitRef().end(); ihit++) {
	THit &h = vecHit[*ihit];
	// All hits are considered, even if not free: so as not to bias the
	// determination of the best candidate.
	double U = h.U; if (U<umn) continue; if (U>umx) break; 
	if (isPixel && fabs(h.V-v)>cvt) continue;
	if (isScifi) {  // Scifis: check time
	  double hT = h.Time; if (hT<tMn || tMx<hT) continue;
	}
	mult++;
	double diff = fabs(U-u); if (diff<best) { best = diff; hBest = &h; }
      }
      if (hBest && hBest->sTrackID().empty() && (mult==1 || best<cut/2)) {
	hitList.push_back(hBest); added++; if (isPixel) added++;
      }
    }
    if (added<4 || added<.75*expected && // And, in order to rescue FI5XY+GP02
	// tracks of the [2008,2009] setup, where only one GP detector is
	// failing, leading to 2, correlated, planes missing and to an overall
	// efficiency of only added=.666*expected
	(expected>6 || nPGs<4)) continue;
    for (int i = 0; i<(int)hitList.size(); i++) tx.AnnexHit(*hitList[i]);
    if (!tx.FullKF(1) || t.Chi2tot/(t.NDFs-5)>TOpt::dCut[17]) continue;
    tx.Type |= 0x2;       // Account for newly added zone 0x2
    listTrack.push_back(tx); t.IMark=tr.IMark = -1;
  }
  //           ***** ERASE TRACK PIECES SUUCESSFULLY BRIDGED  *****
  it = listTrack.begin(); while (it!=listTrack.end()) {
    if ((*it).IMark==-1) { listTrack.erase(it++); }
    else {
      if ((*it).IMark==-2) (*it).IMark = 0;
      it++;
    }
  }
  
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ upDownCompatible ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// - Check compatibility of up- and down-stream segments, based on diffs
//  in dip angle and positions in z and also y.

bool upDownCompatible(int imag, float pinv, float y0, float z0,
		      double ddip, double dz, double dy)
{
#ifdef BridgeS_OVER_TARGET
  if (imag==0) {
    if (fabs(ddip*1000.)>15.) return false;
    if (fabs(dz)>6) return false;
    if (fabs(dy)>7.5) return false;
  }
  else
#endif
  {
    //#define STRICT_BRIDGE_CUTS
#ifdef STRICT_BRIDGE_CUTS
    if (imag==2) {
      if (fabs(pinv)<.075) {
	if (fabs(ddip*1000.)>15.) return false;
	if (fabs(dz)>6) return false;
	if (fabs(dy)>7.5) return false;
      }
      else {
	if (fabs(ddip*1000.)>40) return false;
	if (fabs(dz)>20) return false;
	if (pinv>0) { if (dy<-15 || dy>30) return false; }
	else        { if (dy<-30 || dy>15) return false; }
      }
    }
    else {
      if (fabs(pinv)<.20) {
	if (fabs(y0)<30 && fabs (z0)<30) {
	  if (fabs(ddip*1000.)>7) return false;
	  if (fabs(dz)>.6) return false;
	  if (pinv>0) { if (dy<-2.5 || dy>4  ) return false; }
	  else        { if (dy<-4   || dy>2.5) return false; }
	}
	else {
	  if (fabs(ddip*1000.)>35) return false;
	  if (fabs(dz)>2.5) return false;
	  if  (pinv>0) { if (dy<-2.5 || dy>4  ) return false; }
	  else         { if (dy<-4   || dy>2.5) return false; }
	}
      }
      else {
	if (fabs(y0)<40 && fabs (z0)<40) {
	  if (fabs(ddip*1000.)>55) return false;
	  if (fabs(dz)>2.5) return false;
	  if (pinv>0) { if (dy<-4   || dy>7.0) return false; }
	  else        { if (dy<-7.0 || dy>4  ) return false; }
	}
	else {
	  if (fabs(pinv)<1.1) {
	    if (fabs(ddip*1000.)>180) return false;
	  }
	  else
	    if (fabs(ddip*1000.)>360) return false;
	  if (fabs(dz)>7) return false;
	  if (fabs(pinv)<2) {
	    if (pinv>0) { if (dy<-6   || dy>8.5) return false; }
	    else        { if (dy<-8.5 || dy>6  ) return false; }
	  }
	  else if (fabs(pinv)<2.5) {
	    if (pinv>0) { if (dy<-4  || dy>12 ) return false; }
	    else        { if (dy<-12 || dy>4  ) return false; }
	  }
	  else {
	    if (pinv>0) { if (dy<0   || dy>20 ) return false; }
	    else        { if (dy<-20 || dy>0  ) return false; }
	  }
	}
      }
    }  // End selection for imag==1
#else
    if (imag==2) {
      if (fabs(pinv)<.075) {
	if (fabs(ddip*1000.)>15.) return false;
	if (fabs(dz)>6) return false;
	if (fabs(dy)>7.5) return false;
      }
      else {
	if (fabs(ddip*1000.)>40) return false;
	if (fabs(dz)>20) return false;
	if (pinv>0) { if (dy<-15 || dy>30) return false; }
	else        { if (dy<-30 || dy>15) return false; }
      }
    }
    else {
      //#  define SMALL_SM1_GAP
#  ifndef SMALL_SM1_GAP
      if (fabs(pinv)<.20) {
	if (fabs(y0)<30 && fabs (z0)<30) {
	  if (fabs(ddip*1000.)>14) return false;
	  if (fabs(dz)>.6) return false;
	  if (pinv>0) { if (dy<-2.5 || dy>4  ) return false; }
	  else        { if (dy<-4   || dy>2.5) return false; }
	}
	else {
	  if (fabs(ddip*1000.)>70) return false;
	  if (fabs(dz)>2.5) return false;
	  if  (pinv>0) { if (dy<-2.5 || dy>4  ) return false; }
	  else         { if (dy<-4   || dy>2.5) return false; }
	}
      }
      else {
	if (fabs(y0)<40 && fabs (z0)<40) {
	  if (fabs(ddip*1000.)>110) return false;
	  if (fabs(dz)>2.5) return false;
	  if (fabs(pinv)<1.8) {
	    if (pinv>0) { if (dy<-4   || dy>7.0) return false; }
	    else        { if (dy<-7.0 || dy>4  ) return false; }
	  }
	  else {
	    if (pinv>0) { if (dy<-4   || dy>8.0) return false; }
	    else        { if (dy<-8.0 || dy>4  ) return false; }
	  }
	}
	else {
	  if (fabs(pinv)<1.1) {
	    if (fabs(ddip*1000.)>360) return false;
	  }
	  else
	    if (fabs(ddip*1000.)>720) return false;
	  if (fabs(dz)>7) return false;
	  if (fabs(pinv)<2) {
	    if (pinv>0) { if (dy<-6   || dy>8.5) return false; }
	    else        { if (dy<-8.5 || dy>6  ) return false; }
	  }
	  else if (fabs(pinv)<2.5) {
	    if (pinv>0) { if (dy<-4  || dy>12 ) return false; }
	    else        { if (dy<-12 || dy>4  ) return false; }
	  }
	  else {
	    if (pinv>0) { if (dy<0   || dy>20 ) return false; }
	    else        { if (dy<-20 || dy>0  ) return false; }
	  }
	}
      }
# else
      if (fabs(pinv)<.133) {
	if (fabs(y0)<30 && fabs (z0)<30) {
	  if (fabs(ddip*1000.)>14) return false;
	  if (fabs(dz)>.6) return false;
	  if (pinv>0) { if (dy<-2.5 || dy>4  ) return false; }
	  else        { if (dy<-4   || dy>2.5) return false; }
	}
	else {
	  if (fabs(ddip*1000.)>70) return false;
	  if (fabs(dz)>2.5) return false;
	  if  (pinv>0) { if (dy<-2.5 || dy>4  ) return false; }
	  else         { if (dy<-4   || dy>2.5) return false; }
	}
      }
      else {
	if (fabs(y0)<40 && fabs (z0)<40) {
	  if (fabs(ddip*1000.)>110) return false;
	  if (fabs(dz)>2.5) return false;
	  if (pinv>0) { if (dy<-4   || dy>7.0) return false; }
	  else        { if (dy<-7.0 || dy>4  ) return false; }
	}
	else {
	  if (fabs(pinv)<.733) {
	    if (fabs(ddip*1000.)>360) return false;
	  }
	  else
	    if (fabs(ddip*1000.)>720) return false;
	  if (fabs(dz)>7) return false;
	  if (fabs(pinv)<1.33) {
	    if (pinv>0) { if (dy<-6   || dy>8.5) return false; }
	    else        { if (dy<-8.5 || dy>6  ) return false; }
	  }
	  else if (fabs(pinv)<1.67) {
	    if (pinv>0) { if (dy<-4  || dy>12 ) return false; }
	    else        { if (dy<-12 || dy>4  ) return false; }
	  }
	  else {
	    if (pinv>0) { if (dy<0   || dy>20 ) return false; }
	    else        { if (dy<-20 || dy>0  ) return false; }
	  }
	}
      }
#  endif
    } // End SM2 else SM1
#endif
  } // End BoT else SM2/1
  return true;
}

#ifdef NTUPLE_BRIDGE_CUTS
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ fillNtupleBrigeCuts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Ntuple for tuning the compatibility cuts

void fillNtupleBrigeCuts(TTrack &tl, TTrack &tr, float pimc,
			 int imag, float pinv, float y0, float z0,
			 double ddip, double dz, double dy, double dazi)
{
  const TSetup &setup = TSetup::Ref(); TEv &ev = TEv::Ref();
  const TKine &kine = ev.vKine(tl.GetIKine());
  const vector<THit> &hits = ev.vHit();

  int nh_shower, nhits, nhsame;
  vector<int>::const_iterator ih;
  for (ih = kine.vHitRef().begin(), nh_shower=nhits = 0;
       ih!=kine.vHitRef().end(); ih++) {
    const THit &h = hits[*ih]; if (!setup.vPlane(h.IPlane).IFlag) continue;
    if (h.IOrig==0) nhits++;
    else            nh_shower++;
  }
  nhsame = tl.GetNHsame()+tr.GetNHsame();
  int pileup = kine.isPileup() ? 1 : 0;
  static ofstream ofile; static bool open = true;
  if (open) {
    open = false;
    string home_path;
    if (CsOpt::Instance()->getOpt("","histograms home",home_path)) {
      int i = home_path.find(".root");
      if (i<0) CsErrLog::mes(elFatal,
		 "No \"histograms home\" specified in options file");
      else {
	home_path.replace(i,5,".bridge.out");
	ofile.open(home_path.c_str());
	if (ofile)
	  CsErrLog::mes(elInfo,string("Bridging Ntuple output in \"")+
			home_path+"\"");
	else
	  CsErrLog::mes(elFatal,string("Error opening Bridging Ntuple output \"")+
			home_path+"\"");
      }
    }
  }
  ofile << imag       << " " << pimc       << " " << pinv << " "
	<< nhits      << " " << nh_shower  << " " << nhsame << " "
	<< y0         << " " << z0         << " "
	<< tl.H('d')(0) << " " << tr.H('u')(0) << " "
	<< tl.H('u')(0) << " " << tr.H('d')(0) << " "
	<< ddip*1000 << " " << dz << " "
	<< dazi*1000 << " " << dy << " " << pileup << endl;
}
#endif

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ UpdateKF ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// - ReMode[47]&0x1 = Recycle Fit, i.e. do not automatically refit full tracks,
//  but, in a sequence of "tl,tr" pairs w/ same "tr", refit only the "tl" piece.
// - ReMode[47]&0x2 = No systematic refit after cleaning, i.e. when cleaning the
//  "(tl,tr)" pair (of its worst hit), do not refit systematically, but only the
//  pair of minimum chi2 for the current "tr".

bool TEv::UpdateKF(int imag, int targetField,
		   TTrack &tx, TTrackPair2 &tp,
		   int &prvRId, TTrack &prvRT, double &prvCop, double *Hnoise,
		   TTrackPair2 *minTP, double &minChi2,
		   double &chi2)
{
  const TSetup &setup = TSetup::Ref();
  TTrack &tl = *tp.iL, &tr = *tp.iR;

  bool recycleFit = // Do recyle "tr" fit...
    TOpt::ReMode[47] &&  // ...if requested
    (int)tr.Id==prvRId;  // ...and current "tr" is same as previous

  static map<int,double> prvMChi2;

  tx.InsertMissedPlanes();                     // ***** PREPARE for FullKF REFIT
  tx.Hlast(5) = tp.IFit==0x8 ? tx.Haux(5) : tx.Hfirst(5); 
  tx.Hlast(5,5) = 1.E-6;

  //    ***** TRY AND RECYCLE PART of the FULL FIT of PREVIOUS (l,r) COMBINATION
#ifdef BridgeS_DEBUG_RECYLEFIT
  static int doPrint = 0; TTrack txp = tx; double coP = tx.Hlast(5);
#endif
  if (recycleFit) {
    double cop = tx.Hlast(5), r = cop/prvCop;
    if (0.6<r && r<1.4) { // Extra condition: new (tl,tr) has similar P
      prvRT.Hlast(5) = cop; prvRT.Hlast(5,5) = 1.E-20;
      prvRT.QuickKF(-1,1);
      THlx &H = prvRT.Hfirst; int i,j,k;
      for (i = 1, k = 0; i<=5; i++) for (j = 1; j<=i; j++, k++) {
	  H(i,j) += Hnoise[k]*r*r;
	}
      int ipl = tr.lPlnRef.front();
#ifdef BridgeS_DEBUG_RECYLEFIT
      // Let's do the recycling fit, and compare w/ the complete fit, cf. infra.
      TTrack txp = tx; txp.mChi2 = prvMChi2; txp.mHub[ipl] = H;
      txp.FullKF(-1,ipl-1); recycleFit = false;
      if (doPrint) { H.Print(); txp.Hfirst.Print(); }
#else
      tx.mChi2 = prvMChi2; tx.mHub[ipl] = H;
      if (!tx.FullKF(-1,ipl-1)) return false;
#endif
    }
    else recycleFit = false;
  }

  if (!recycleFit && !tx.FullKF(-1))  // ***** NO RECYCLING: COMPLETE FullKF FIT
    return false;

  chi2 = tx.Chi2tot/(tx.NDFs-5); float oldChi2 = chi2;

#ifdef BridgeS_DEBUG_RECYLEFIT
  if (doPrint) {
    tx.Hfirst.Print();
    double diff, di[5]; int i,j; TMatrixDSym cov(5);
    printf("%.3f %.3f %6.2f%%\n diff:",1/prvCop,1/coP,coP/prvCop*100);
    for (i = 1; i<=5; i++) {
      di[i-1] = txp.Hfirst(i)-tx.Hfirst(i); double d = tx.Hfirst(i,i);
      printf("  %d %6.2f",i,di[i-1]/sqrt(d));
      cov(i-1,i-1) = d;
      for (j = 1; j<i; j++) {
	cov(i-1,j-1)=cov(j-1,i-1) = tx.Hfirst(i,j);
      }
    }
    TMatrixDSym c = TMatrixDSym(TMatrixDSym::kInverted,cov);
    for (i = 0, diff = 0; i<5; i++) {
      diff += di[i]*di[i]*c(i,i);
      for (j = 0; j<i; j++) diff += 2*di[i]*di[j]*c(i,j);
    }
    printf(" => %6.2f\n",sqrt(diff));
    printf("chi2: %6.2f  %6.2f\n",txp.Chi2tot/(txp.NDFs-5),chi2);
  }
#endif

  if (imag==0 && targetField<2) {// BridgeOverTarget, w/o transverse field
    // In FullKF, w/ ROOTG (or dE/dX MatMap), P evolves even if P is ''fixed'' (
    // because P can only be fixed by init'ing d(1/P)<<1, which insulates P
    // against update by the KF fit, but not against dE/dx). => Let's reset 1/P,
    // and d(1/P).
    // (Note that, at this stage, the actual value assigned to 1/P, and d(1/P),
    // does not matter so much, since it's not the final assignment. But later
    // on, in "TracksFit2", the question will rise again: should dCut[4] be
    // assigned to "Hfirst" and integrated dE subtracted from "Hlast"? what
    // value for d(1/P)?)
    tx.Hfirst(5) = TOpt::iCut[15]/TOpt::dCut[4];
    tx.Hfirst(5,5) = 1.E-20;
  }
  tp.Hfirst = tx.Hfirst;

  if ((int)tr.Id!=prvRId) {  // New "tr": store its map of chi2 increments
    prvCop = tx.Hlast(5); prvRId = tr.Id; prvMChi2 = tx.mChi2;
    int iplMx = setup.vIplFirst()[imag];
    map<int,double>::iterator idx = prvMChi2.begin();
    while (idx!=prvMChi2.end() && (*idx).first<=iplMx) idx++;
    prvMChi2.erase(prvMChi2.begin(),--idx);
  }

  int worstPl = -1; if (tp.WorstPl==-1) {                   // ***** CLEANING...
    // - The idea is to clean away ghost hits degrading the whole chi2 and
    //  penalizing an otherwise genuine track when it comes to compete w/ other
    //  combinations making use of the same either up- or down-stream segment.
    // - A single hit (the worst one) may only be cleaned away.
    // - This cleaning is particularly welcome in the FullKF context, where such
    //  isolated ghost hit defects cannot be easily 'diluted' by noise from
    //  multi-scattering, whereas otherwise more profound but more largely
    //  spread defects can.
    // - The cleaning is conditioned by a cut on the worst chi2 increment.
    // - The cleaned track is refitted via a partial FullKF, minimizing CPU.
    //  Upon option, this is only done if the "(tl,tr)" pair is the optimum of
    //  the current "tr".
    // (Nota bene: cleaning is hence not done when the FullKF is already used
    // supra, in the iterative adjustment of momentum.)
    int worstGr = -1; double worstChi2Incr;
    tx.WorstHit(tx.mChi2,worstPl,worstGr,worstChi2Incr);
    TTrack &worstT = worstGr==imag ? tr : tl; bool enoughHits;
    if (worstT.NDFs>6)      enoughHits = true;
    else if (worstT.NDFs<4) enoughHits = false;
    else { // We may have, after cleaning, a not enough contraining hit pattern:
      // (e.g. case of a MP01X+Y+FI04X+U 0x1 track in 2012 DVCS). Which would
      // result in an unduly low chi2/NDF for the track fit, which would in turn
      // give it an edge over its competitors.
      // - Unconstrained cases are handled in "BridgeSegments2", by giving a
      //  malus to track below par, cf. "belowParT". But this strategy, w/ steps
      //  forward (cleaning) and steps backward (malus), is difficult to tune.
      // - Let's try here something more clear-cut and refrain from cleaning
      //  when it's bound to yield a badly unconstrained track.
      // => To that end, let's require at least two coordinates measured twice,
      //  - of which one is vertical (i.e. Z in TraFFiC convention),
      //  - understanding "coordinate" in a broader sense: either Y, Z, U or V.
      //   and "twice" in a stricter sense: counting the 2 associated planes of
      //   a drift-like TStation only ounce.
      // Note: Loose hits patterns get another opportunity in "BackTrackSAS",
      // where the momentum is already constrained by the 0x6 track and the
      // requirements on the 0x1 piece can hence be loosened.
      static unsigned int proj; // "proj" is always defined when used in the for blcok infra => let's avoid compiler warning w/ a "static" declaration
      unsigned int projPrv, allProjs, match; const TStation *sPrv;
      list<int>::iterator ih;
      for (ih = worstT.lHitPat.begin(), projPrv=allProjs=match = 0, sPrv = 0;
	   ih!=worstT.lHitPat.end(); ih++) {
	if (*ih<0) continue;
	const TPlane &p = setup.vPlane(vecHit[*ih].IPlane);
	if (p.IFlag&0x30) {
	  match = 0x3; break; // A pixel ensures that conditions are met
	}
	else {
	  unsigned int proj = 1<<p.IProj; if (!(proj&0x3)) {
	    int angle = setup.vProj()[p.IProj];
	    if      (abs(angle)<225)     proj = 0x1;
	    else if (abs(angle-450)<225) proj = 0x4;
	    else if (abs(angle-900)<225) proj = 0x2;
	    else if (abs(angle+450)<225) proj = 0x8;
	    else                         proj = 0;
	  }
	}
	const TStation *&s = p.Station; if (s==sPrv && proj==projPrv) continue;
	sPrv = s; projPrv = proj;
	if (proj&allProjs) { match |= proj==0x2 ? 0x2 : 0x1; allProjs ^= proj; }
	allProjs |= proj;
      }
      enoughHits = match==0x3;
    }
    static double f = 12.25; // This presently built-in parameter could be made an option in a future version...
    if (worstChi2Incr>f && enoughHits) {
      int jhit = tx.CancelHit_Clip(worstPl,false); bool unconstrained;
      if (!(TOpt::ReMode[47]&0x2)) {// Systematic refit after cleaning
	if (tx.FullKF(-1,worstPl)) {                // ***** ...PARTIAL RE-REFIT
	  chi2 = tx.Chi2tot/(tx.NDFs-5);
	  if (imag==0 && targetField<2) { // BoT, w/o transverse field
	    tx.Hfirst(5) = TOpt::iCut[15]/TOpt::dCut[4];
	    tx.Hfirst(5,5) = 1.E-20;
	  }
	  unconstrained = // Are deemed to mean unconstrained...
	    tx.Hfirst(1,1)>0.0225 || tx.Hfirst(2,2)>0.0225 ||// ...1.5mm
	    tx.Hfirst(3,3)>25.e-6 || tx.Hfirst(4,4)>25.e-6;  // ...5 mrd
	}
	else {
	  chi2 = oldChi2; unconstrained = false;
	}
      }
      else {
	chi2 = (tx.Chi2tot-worstChi2Incr)/(tx.NDFs-5); unconstrained = false;
      }
      if (chi2>oldChi2*.9 || // If cleaning not significant or...
	  // ...leave the track unconstrained
	  unconstrained) {                          // ***** RESTORE INITIAL FIT
	tx.RestoreHit(worstPl,jhit,false); // Restore track (which apart from restoring TTrack::NHits/NDFs is probably useless...)
	chi2 = oldChi2; tx.Chi2tot = chi2*(tx.NDFs-5); // Restore chi2
	worstPl = -1;
      }
      else {
	bool transient = // If track already very good or cleaning poorly
	  // significant, let's not strip it definitively of its worst hit but
	  // let's not restore the old chi2 so as not to give an unfair
	  // advantage to its competitors.
	  oldChi2<TOpt::dCut[16] || chi2>oldChi2*.8;
	if (transient) {
	  tx.RestoreHit(worstPl,jhit,false);
	  tx.Chi2tot = oldChi2*(tx.NDFs-5); // Restore chi2 of track  ("tx") but not that of track pair ("tp")  
	}
	else {
	  tp.Hfirst = tx.Hfirst; // Update "tp.Hfirst" w/ cleaned track helix.
	  tp.WorstPl = worstPl;  // Store "worstPl" in "tp": that it can later be definitely erased.
	}
      }
      //else "chi2" and "tp.Hfirst" are left unchanged. 
    }
  }
  bool updateMinTP = !minTP || chi2<minChi2;
  if (updateMinTP &&                  // If current "(tl,tr)" is new minimum...
      worstPl>=0 && tp.WorstPl>=0 &&  // ...and is a cleaned track...
      (TOpt::ReMode[47]&0x2)) {       // ...yet not refitted => do refit
    if (tx.FullKF(-1,worstPl)) {
      chi2 = tx.Chi2tot/(tx.NDFs-5);
      if (chi2>TOpt::dCut[97]) return false;
      if (imag==0 && targetField<2) { // BoT, w/o transverse field
	tx.Hfirst(5) = TOpt::iCut[15]/TOpt::dCut[4];
	tx.Hfirst(5,5) = 1.E-20;
      }
      bool unconstrained =
	tx.Hfirst(1,1)>0.0225 || tx.Hfirst(2,2)>0.0225 ||// ...1.5mm
	tx.Hfirst(3,3)>25.e-6 || tx.Hfirst(4,4)>25.e-6;  // ...5 mrd
      if (chi2>minChi2 ||    // If refit turns out to yield no new minimum...
	  unconstrained)     // ...or an unconstrained track (cf. supra)...
	updateMinTP = false; // ...back off
    }
    else updateMinTP = false;
  }
  if (updateMinTP) {
    minTP = &tp; minChi2 = chi2;
  }
  tp.IFit = 0x20;
  return true;
}
