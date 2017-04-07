// $Id$

// Various TEv method for handling yoke tracks.

#include "CsErrLog.h"
#include "TEv.h"
#include "TSetup.h"

using namespace std;

void TEv::SplitYokeTr(TTrack &t, TTrack *&yokeTr, const THlx *&yokeTrHlast)
{
  if (t.IMark<=0) { yokeTr = 0; yokeTrHlast = 0; return; } 
  const TSetup &setup = TSetup::Ref();
  yokeTr = new TTrack(t); yokeTr->Behead(setup.vIplLast()[1]);
  const TPlane &pFirst = setup.vPlane(vHit(yokeTr->lHPat().front()).IPlane);
  yokeTr->Hfirst(0) = setup.vDetect(pFirst.IDetRef).X(0);
  t.Shorten(setup.vIplFirst()[2]);
  const TPlane &pLast = setup.vPlane(vHit(t.lHPat().back()).IPlane);
  double zLast = setup.vDetect(pLast.IDetRef).X(0);
  // The upstream piece's last helix was pushed into the list of smoothed
  // helices (cf. "TEv::TracksFit2): retrieve it.
  list<THlx>::iterator ihsm;
  for (ihsm = t.lHsm.begin(), yokeTrHlast = 0; ihsm!=t.lHsm.end(); ihsm++) {
    THlx &hsm = *ihsm; if (fabs(hsm(0)-zLast)<.01) {
      t.Hlast = hsm; yokeTrHlast = &hsm;
    }
  }
  if (!yokeTrHlast)
    CsErrLog::msg(elFatal,__FILE__,__LINE__,
		  "SM2 yoke track %d has no helix @ the end point of its upstream piece",
		  t.Id);
}
void TEv::AddYokeTrTail(TTrack &t, TTrack *yokeTr)
{
  // Reunite the upstream (arg. "&t") and downstream piece "*yokeTr".
  if (t.IMark<=0)
    CsErrLog::msg(elFatal,__FILE__,__LINE__,
		  "Track %d has is not a yoke track (IMark=%d)",t.Id,t.IMark);
  const TSetup &setup = TSetup::Ref();

  // Save characteristics of the 2 pieces prior to merging.
  THlx hs[2]; hs[0] = t.Hlast;
  int iplTlL = t.lPRef().back(), iplTrF = yokeTr->lPRef().front();
  double lastCop = t.Hlast(5), lastDcop = t.Hlast(5,5);
  double radLenFr = t.radLenFr;
  // Only meaningful fit is that of the 0x3 piece: assign its chi2/NDF to
  // the whole track.
  t.Chi2tot = t.Chi2tot*(t.NDFs+yokeTr->NDFs-5)/(t.NDFs-5);
  t.Append(*yokeTr);
  t.Hlast(5) = lastCop; t.Hlast(5,5) = lastDcop;
  // Insert missed planes in the interval between the end of the 0x3
  // piece and the beginning of the downstream piece, and then update the
  // hits map. The later cannot be done using "TTrack::Refine", given that
  // the arrays of helices have not been defined in the interval.
  if (TOpt::ReMode[21]==0) t.InsertMissedPlanes();
  list<int>::iterator ih; list<int>::const_iterator ip; int i12;
  for (ih = t.lHitPat.begin(), ip = t.lPRef().begin(), i12 = 0;
       ih!=t.lHitPat.end(); ih++, ip++) {
    int ipl = *ip; if (ipl<=iplTlL) continue;
    if (ipl<iplTrF && *ih!=-1) CsErrLog::msg(elError,__FILE__,__LINE__,
 "Inconsistency when extending SM2 yoke ID=%d: hit reference on TPlane %d = %d",
					     t.Id,ipl,*ih);
    const TPlane &p = setup.vPlane(ipl);
    const TDetect &d = setup.vDetect(p.IDetRef);
    THlx &hStrt = hs[i12], &hStop = hs[1-i12]; i12 = 1-i12;
    hStop(0) = d.X(0); hStrt.Extrapolate(hStop,true);
    if (!d.InActive(hStop(1),hStop(2))) {
      if (ipl<iplTrF) *ih = -2;
    }
    else if (TOpt::ReMode[20]<=0 || !setup.InMaterialMap(hStop(0))) {
      double len    = d.Siz(0)/hStop.DirCos(1); // Detector's thickness
      double radLen = d.RadLen;                 // Detector's rad. length
      radLenFr += len/radLen;
    }
  }
  t.radLenFr = radLenFr + hs[i12].RadLenFr();
}
