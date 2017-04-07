// $Id: TTrackUseTime.cc 14094 2015-11-06 15:28:48Z lsilva $

#include "CsCluster.h"
#include "TOpt.h"
#include "TSetup.h"
#include "TEv.h"
#include "TTrack.h"

using namespace std;

/*!
  \file   TTrackUseTime.cc
  \brief  Set track mean time and uncertainty (if track has time measurements). Calculate also time dispersion.
  \param \e clean = if !=0, "clean" ouliers, defined as scifi or hodo hits outside of ( mean +/- \e clean * uncertainty ). "Clean" means reset uncertainty to very large value and recompute track's mean time and uncertainty accordingly. And, if in addition track's chi2 is bad, remove worst outlier.
*/

/*
  Notes:
  i) The "clean" option has to be enabled not too early: only at a stage of
    reconstruction when hits can no longer be shared by several tracks, i.e.
    after "TEv::CleanTrackList" has been executed (for the obvious reason that
    the operation affects the "cleaned" hits).
 ii) Dispersion takes as weights for scifis and hodos, the detector resolution
    instead of hit uncertainty. This in order to work around the possible
    resetting of the latter by the clean option and so that "DispTime" remains
    always a good estimator of the time consistency of the track, even after
    cleaning.
iii) Returns false when cleaning goes wrong and track has to be erased.
*/
bool TTrack::UseHitTime(int clean)
{
  const TSetup& setup = TSetup::Ref();
  const TEv &ev = TEv::Ref();

  list<int>::const_iterator ihp; int n1; double s1, st, z1, zt, zt2;
  for (ihp = lHitPat.begin(), n1 = 0, s1=st=z1=zt=zt2 = 0;
       ihp!=lHitPat.end(); ihp++) {
    // ***** LOOP OVER HIT REFERENCES *****
    if (*ihp<0) continue;
    const THit &h = ev.vHit(*ihp);
    if (h.SigT>0) { // Hit does have time
      double w = 1/h.SigT/h.SigT;
      s1 += w; st += h.Time*w;
      n1++;
      const TDetect &d = setup.vDetect(h.IPlane);
      if (d.IType==22 /* scifi */|| d.IType>=41 /* hodos */)
	w = 1/d.TResol/d.TResol;
      z1 += w; zt += h.Time*w; zt2 += h.Time*h.Time*w;
    }
  }

  if (n1) {
    MeanTime = st/s1; SigmaTime = sqrt(1/s1);
    if (n1>1) DispTime = sqrt(zt2*z1-zt*zt)/z1;
    else      DispTime = -1;
  }
  else {
    // Time may have been set at an earlier stage of the reco when
    // TTrack had some time info, which it subsequently lost...
    SigmaTime=DispTime = -1;  //  ...=> Reset "SigmaTime", "DispTime"
  }

  if (0<SigmaTime && SigmaTime<.5 && // Precisely enough timed
      clean) {

    //                    ********** CLEANING **********
    // (Note: In principle scifi hit timing is checked in the PreP. But this
    // only applies to scifi-only or almost-only segments. When segments are
    // concatenated via bridging we may still face inconsistent hit times.
    // E.g. 07W27/cdr32001-58926, event #1055129.)

    int nHits[5], nScifis[5]; // We include a fifth zone: to be on the safe side
    const vector<int> &iplF = setup.vIplFirst(), &iplL = setup.vIplLast();
    int nZones = iplF.size(), zone, iplWorst; double worstDev; bool modified;
    for (ihp = lHitPat.begin(), s1=st=z1=zt=zt2 = 0,
	   memset((void*)nHits,0,sizeof(nHits)),
	   memset((void*)nScifis,0,sizeof(nScifis)),
	   iplWorst = -1, worstDev = 0; ihp!=lHitPat.end(); ihp++) {
      // ***** LOOP OVER HIT REFERENCES *****
      if (*ihp<0) continue;
      const THit &h = ev.vHit(*ihp);
      const TPlane &p = setup.vPlane(h.IPlane); int ipl = h.IPlane;
      const TDetect &d = setup.vDetect(p.IDetRef);
      int iz; for (iz=zone = 0; iz<nZones; iz++) {
	if (iplF[iz]<=ipl && ipl<=iplL[iz]) { zone = iz; break; }
      }
      nHits[zone]++;
      if (h.SigT<=0) continue;
      double time = h.Time, sigt = h.SigT;
      if (d.IType==22 /* scifi */|| d.IType>=41 /* hodos */) {
	if (d.IType==22) nScifis[zone]++;
	double dev = fabs(time-MeanTime)/sqrt(sigt*sigt+SigmaTime*SigmaTime);
	if (dev>clean) {   // ***** FOUND BADLY TIMED HIT... *****
	  if (dev>worstDev) {
	    iplWorst = ipl; worstDev = dev;
	  }
	  THit &hp = const_cast<THit&>(ev.vHit(*ihp));
	  sigt = 100; hp.setSigT(sigt);
	  hp.PtrClus()->setTime(hp.Time,-2);
	}
	double w = 1/sigt/sigt; s1 += w; st += time*w;
	w = d.TResol;           z1 += w; zt += time*w; zt2 += time*time*w;
      }
      else {
	double w = 1/sigt/sigt;
	s1 += w; st += time*w;
	z1 += w; zt += time*w; zt2 += time*time*w;
      }
    }
    if (iplWorst>=0) {            // ***** DISCARDING WORST HIT? *****
      MeanTime = st/s1; SigmaTime = sqrt(1/s1);
      bool badChi2Track = NDFs<5 || Chi2tot/(NDFs-5)>TOpt::dCut[16];
      if (badChi2Track) {
	// Processing depends upon whether bad or good track:
	// if the latter => simply redefine the sigma of the badly timed hit (so
	// that it does not weight any longer on track's mean time: neither for
	// the present nor for later timings of the track. 
	const TDetect &d = setup.vDetect(iplWorst);
	int iz; for (iz=zone = 0; iz<nZones; iz++) {
	  if (iplF[iz]<=iplWorst && iplWorst<=iplL[iz]) { zone = iz; break; }
	}
	int nHs = nHits[zone]-1;
	if (d.IType==22 && (1<<zone&TOpt::iCut[16])) { // Enhance scifis
	  // (Note that this is done crudely, and not in accordance w/ what's
	  // in "Talgo2::FindSpace".)
	  if (nScifis[zone]>1) nHs += 2*(nScifis[zone]-1);
	}
	if (nHs>=TOpt::iPRpar[10*zone+5]) {
	  // Enough hits otherwise => we can safely discard the hit
	  int ihit = CancelHit_Clip(iplWorst,true);
	  nScifis[zone]--; if (!nScifis[zone]) Scifi &= ~(1<<zone);
	  if (n1>1) {	                          // Update dispersion
	    double w = d.TResol; const THit &h = ev.vHit(ihit);
	    z1 -= w; zt -= h.Time*w; zt2 -= h.Time*h.Time*w;
	    DispTime = sqrt(zt2*z1-zt*zt)/z1;
	  }
	  else DispTime = -1;
	  iplWorst = -2; // Track is to be refit
	}
      }
    }
    if (iplWorst==-2) {       // ***** WORST HIT DISCARDED: REFIT
      bool ret = true;
      IFit &= ~0x8; // Do not update QN fit
      if (IFit&0x2) {
	ret = QuickKF(-1,1); if (ret && (IFit&0x4))  ret = QuickKF(1,1);
      }
      if (IFit&0x20) {
	ret = FullKF(-1);    if (ret && (IFit&0x40)) ret = FullKF(1);
      }
      return ret;
    }
  }

  return true;
}
