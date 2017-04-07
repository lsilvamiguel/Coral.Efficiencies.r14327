// $Id: TTrackAddHit.cc 14069 2015-09-17 20:44:46Z lsilva $

#include <iostream>
#include <iomanip>
#include <cmath>
#include "CsErrLog.h"
#include "TOpt.h"
#include "TSetup.h"
#include "TEv.h"
#include "THit.h"
#include "TTrack.h"

#ifdef COMPASS_USE_OSPACE_STD
# include <ospace/std/algorithm>
#else
# include <algorithm>  
#endif

using namespace std;

/*!
  Bunch of methods to add and remove THit's.
  - Add(remove) hit and plane numbers in lPlnRef and lHitPat of the track.
  - Add(remove) this track ID to the set of associated track ID of the hit.

  - Basic methods: "AddHit", "SubHit".
  - Extra methods: "Merge", "AnnexHit", "ReplaceHit", "SubHit", "CancelHit",
    "CancelHit_Clip", "RestoreHit", "Shorten", "Behead", "Clip", "SingleDiff".
  - CAVEAT:
     - All these methods do not update "sIHit"
     - Do not update "sProj" (i.e. TTrack's set of proj.). In fact, "sProj" is
      completely disregarded by TraFDic, for it could only be maintained, cf.
      case of a method removing THit's, at the expense of much CPU => Better
      recalculate the # of proj. each time it's needed.
  - Pixel planes (flagged w/ TPlane::IFlag&0x30): update "NDFs".
  - Method "BendingInfo".
*/


void TTrack::AddHit(const THit &h)
{
  int ipl = h.IPlane, ih  = h.IHit;

  list<int>::iterator ipr = lPlnRef.begin();
  list<int>::iterator ihp = lHitPat.begin();
  while (ipr!=lPlnRef.end() && *ipr<ipl) { ipr++; ihp++; }
  if (ipr!=lPlnRef.end() && *ipr==ipl) {
    if (*ihp>=0) {
      if (TOpt::ReMode[4]!=0) return; // Skip warning message if ideal PR
      int ikine = h.IKine<-1 ? -2-h.IKine : h.IKine;
      CsErrLog::msg(elError,__FILE__,__LINE__,
		    "Another hit (#%d, TKine #%d) for track ID=%d on plane #%d",
		    ih,ikine,Id,ipl); return;
    }
    else *ihp = ih;
  }
  else {
    lPlnRef.insert(ipr,ipl);      // Insert new plane ref.
    lHitPat.insert(ihp,ih);       // Insert new hit ref.
  }

  NHits++;
  NDFs += (TSetup::Ref().vPlane(ipl).IFlag&0x30) ? 2 : 1;

  const_cast<THit&>(h).addTrackID(Id); // Cast away constness and add hit->track reference
}

bool TTrack::Merge(const TTrack &ti, const TTrack &tj)
{
  //  This track is being built out of the merging the 2 tracks supplied as
  // arguments.
  //  It is checked that the 2 latter tracks have no TPlane in common. There ia
  // an exceptions though: a HO03 is allowed, provided it has a same THit on
  // the 2 tracks.
  //  The check is done before the actual building the new track is undertaken,
  // in a distinct, initial loop over the plane references. Thus preventing the
  // new track from being referenced in the "setTrackID" of the THit's.
  // (Note: Good example to check the algorithm can be dstar.2006.03#140.)
  TEv &ev = TEv::Ref();
  list<int>::const_iterator ipr = ti.lPlnRef.begin(), jpr = tj.lPlnRef.begin();
  list<int>::const_iterator ihp = ti.lHitPat.begin(), jhp = tj.lHitPat.begin(); 
  list<int>::const_iterator iend = ti.lPlnRef.end(), jend = tj.lPlnRef.end();
  while (ipr!=ti.lPlnRef.end() || jpr!=tj.lPlnRef.end()) {
    //    ********** DOUBLE LOOP ON ti AND tj TPlane/THit's **********
    if      (ipr!=iend && jpr!=jend && *ipr==*jpr) {
                               // ***** REQUIRE NO TPlane COMMON TO ti and tj...
      if (TSetup::Ref().vDetect(*ipr).IType!=41) // ...except if it's hodo plane
	return false;
      ipr++; jpr++;
    }
    if      (ipr!=iend && (jpr==jend || *ipr<*jpr)) ipr++;
    else if (jpr!=jend && (ipr==iend || *ipr>*jpr)) jpr++;
  }
  ipr = ti.lPlnRef.begin(); jpr = tj.lPlnRef.begin();
  int nCommons = 0;
  while (ipr!=ti.lPlnRef.end() || jpr!=tj.lPlnRef.end()) {
    //    ********** DOUBLE LOOP ON ti AND tj TPlane/THit's **********
    // (Cast away constness and add hit->track reference)
    if      (ipr!=iend && jpr!=jend && *ipr==*jpr) {
      lPlnRef.push_back(*ipr); lHitPat.push_back(*ihp);
      const THit &h = ev.vHit(*ihp); const_cast<THit&>(h).addTrackID(Id);
      ipr++; ihp++; jpr++; jhp++; nCommons++;
    }
    else if (ipr!=iend && (jpr==jend || *ipr<*jpr)) {
      lPlnRef.push_back(*ipr); lHitPat.push_back(*ihp);
      const THit &h = ev.vHit(*ihp); const_cast<THit&>(h).addTrackID(Id);
      ipr++; ihp++;
    }
    else if (jpr!=jend && (ipr==iend || *ipr>*jpr)) {
      lPlnRef.push_back(*jpr); lHitPat.push_back(*jhp);
      const THit &h = ev.vHit(*jhp); const_cast<THit&>(h).addTrackID(Id);
      jpr++; jhp++;
    }
  }
  NHits = ti.NHits+tj.NHits-nCommons;
  NDFs  = ti.NDFs +tj.NDFs -nCommons; // Common hits cannot be of a pixel type
  return true;
}

void TTrack::AnnexHit(const THit &h)
{
  // - The difference w.r.t. to "AddHit" is that here the plane to which arg
  //  THit belongs may already have been inserted in "lPlnRef". "AnnexHit" is to
  //  be used in the latest steps of the track reco, after "InsertMissedPlanes"
  //  has been performed.
  // - NOTA BENE: For revisons > 12334, and since TTrack::FullKF can be used
  //  in such early step as Bridging, "AnnexHit" is identical to "AddHit".
  int ipl = h.IPlane, ih  = h.IHit;

  list<int>::iterator ipr = lPlnRef.begin();
  list<int>::iterator ihp = lHitPat.begin();
  while (ipr!=lPlnRef.end() && *ipr<ipl) { ipr++; ihp++; }
  if (ipr!=lPlnRef.end() && *ipr==ipl) {
    if (*ihp>=0) {
      if (TOpt::ReMode[4]!=0) return; // Skip warning message if ideal PR
      int ikine = h.IKine<-1 ? -2-h.IKine : h.IKine;
      CsErrLog::msg(elError,__FILE__,__LINE__,
		    "Another hit (#%d, TKine #%d) for track ID=%d on plane #%d",
		    ih,ikine,Id,ipl); return;
    }
    else *ihp = ih;
  }
  else {
    lPlnRef.insert(ipr,ipl);      // Insert new plane ref.
    lHitPat.insert(ihp,ih);       // Insert new hit ref.
  }

  NHits++;
  NDFs += (TSetup::Ref().vPlane(ipl).IFlag&0x30) ? 2 : 1;

  const_cast<THit&>(h).addTrackID(Id); // Cast away constness and add hit->track reference
}

const THit *TTrack::ReplaceHit(const THit *h)
{
  int ipl = h->IPlane, ih  = h->IHit;

  list<int>::iterator ipr = lPlnRef.begin(); 
  list<int>::iterator ihp = lHitPat.begin(); 
  while (ipr!=lPlnRef.end() && *ipr<ipl) {ipr++; ihp++; }
  if (ipr!=lPlnRef.end() && *ipr==ipl && *ihp>=0) {
    TEv &ev = TEv::Ref(); const THit *hp = &(ev.vHit(*ihp));
    const_cast<THit*>(hp)->eraseTrackID(Id); // Cast away constness and cancel hit->track reference
    *ihp = ih;
    const_cast<THit*>(h)->addTrackID(Id);    // Cast away constness and add    hit->track reference
    return hp;
  }
  else {
      CsErrLog::msg(elError,__FILE__,__LINE__,
	"Attempting to replace an inexisting hit for track ID=%d on plane #%d",
		    Id,ipl); return 0;
  }
}

void TTrack::SubHit(const THit &h)
{
  int ipl = h.IPlane;

  list<int>::iterator ipr = lPlnRef.begin(); 
  list<int>::iterator ihp = lHitPat.begin(); 
  while (ipr != lPlnRef.end() && *ipr<ipl) { ipr++; ihp++; }
  if (ipr==lPlnRef.end() || *ipr!=ipl || *ihp<0) {
    CsErrLog::msg(elError,__FILE__,__LINE__,
		  "Erasing non associated hit (#%d, plane #%d) in track ID=%d",
		  h.IHit,ipl,Id); return;
  }
  lPlnRef.erase(ipr);      // Erase plane ref.
  lHitPat.erase(ihp);      // Erase hit ref.

  NHits--;
  NDFs -= (TSetup::Ref().vPlane(ipl).IFlag&0x30) ? 2 : 1;

  const_cast<THit&>(h).eraseTrackID(Id); // Cast away constness and erase hit->track reference
}

int TTrack::SubHit(int ipl)
{
  list<int>::iterator ipr = lPlnRef.begin(); 
  list<int>::iterator ihp = lHitPat.begin(); 
  while (ipr!=lPlnRef.end() && *ipr<ipl) { ipr++; ihp++; }
  if (ipr==lPlnRef.end() || *ipr!=ipl || *ihp<0) {
    CsErrLog::msg(elError,__FILE__,__LINE__,
		  "Erasing non associated hit from plane #%d in track ID=%d",
		  ipl,Id); return -1;
  }

  // first erase the hit -> track reference
  TEv &ev = TEv::Ref(); const THit &h = ev.vHit(*ihp);
  const_cast<THit&>(h).eraseTrackID(Id); // Cast away constness and erase hit->track reference
  int ihit = *ihp;

  // then erase all references of the plane/hit to the track
  lPlnRef.erase(ipr);      // Erase plane ref.
  lHitPat.erase(ihp);      // Erase hit ref.

  NHits--;
  NDFs -= (TSetup::Ref().vPlane(ipl).IFlag&0x30) ? 2 : 1;
  return ihit;
}

void TTrack::SubHits(int ipli, int iplf)
{
  TEv &ev = TEv::Ref(); const TSetup &setup = TSetup::Ref();
  list<int>::iterator ipr = lPlnRef.begin();
  list<int>::iterator ihp = lHitPat.begin();
  while (ipr!=lPlnRef.end() && *ipr<ipli) { ipr++; ihp++; }
  while (ipli<=*ipr && *ipr<=iplf) {
    // Erase all ref's comprised between "ipli" and "iplf" (included)
    lPlnRef.erase(ipr++);      // Erase plane ref.
    if (*ihp>=0) {
      NHits--; 
      const THit &h = ev.vHit(*ihp);
      const_cast<THit&>(h).eraseTrackID(Id); // Cast away constness and erase hit->track reference
      NDFs -= (setup.vPlane(h.IPlane).IFlag&0x30) ? 2 : 1;
    }
    lHitPat.erase(ihp++);      // Erase hit ref.
  }
}

int TTrack::CancelHit(int ipl, bool updateTID)
{
  // - The difference w.r.t. "SubHit" is that the plane ref. corresponding to
  //  the erased hit is kept.
  // - Also, upon option, the hit->track reference may be left not updated:
  //  useful when dealing w/ temporary tracks.
  // - It is meant to be used at a later stage of the track reconstruction,
  //  where references to all planes w/in track's range have been inserted
  //  in preparation for a FullKF fit.
  // - The method does not modify the meaning of iterators defined on either
  //  the plane or the hit list, which allows it to be used w/in a loop
  //  based on such iterators.
  list<int>::iterator ipr = lPlnRef.begin(); 
  list<int>::iterator ihp = lHitPat.begin(); 
  while (ipr != lPlnRef.end() && *ipr<ipl) {ipr++; ihp++; }
  if (ipr==lPlnRef.end() || *ipr!=ipl || *ihp<0) {
    return -1;
  }
  TEv &ev = TEv::Ref(); const THit &h = ev.vHit(*ihp);
  if (updateTID)
    const_cast<THit&>(h).eraseTrackID(Id); // Cast away constness and erase hit->track reference
  int ihit = *ihp; *ihp = -1;

  NHits--;
  NDFs -= (TSetup::Ref().vPlane(ipl).IFlag&0x30) ? 2 : 1;
  return ihit;
}

int TTrack::CancelHit_Clip(int ipl, bool updateTID)
{
  // - Erase hit on plane "ipl". Return its index.
  // - Upon option, update the hit->track reference.
  // - If hit is first(last), erase all plane references up(down) to next hit.
  int ihit; if (ipl==lPlnRef.back()) {
    list<int>::reverse_iterator ipr = lPlnRef.rbegin(), ihp = lHitPat.rbegin(); 
    while (ipr!=lPlnRef.rend() && *ipr>ipl) { ipr++; ihp++; }
    if (ipr==lPlnRef.rend() || *ipr!=ipl || *ihp<0) return -1;
    ihit = *ihp; *ihp = -1;
    while (ipr!=lPlnRef.rend()) {
      // Erase all ref's up to first (==first encountered) THit (excluded)
      if (*ihp>=0) break;
      else {
	lPlnRef.erase(--ipr.base()); lHitPat.erase(--ihp.base());// Erase hit and plane ref.
      }
    }
  }
  else {
    list<int>::iterator ipr = lPlnRef.begin(), ihp = lHitPat.begin(); 
    while (ipr!=lPlnRef.end() && *ipr<ipl) { ipr++; ihp++; }
    if (ipr==lPlnRef.end() || *ipr!=ipl || *ihp<0) return -1;
    ihit = *ihp; *ihp = -1;
    if (ipl==lPlnRef.front()) {
      while (ipr!=lPlnRef.end()) {
	// Erase all ref's up to first (==first encountered) THit (excluded)
	if (*ihp>=0) break;
	else {
	  ipr = lPlnRef.erase(ipr); ihp = lHitPat.erase(ihp);// Erase hit and plane ref.
	}
      }
    }
  }

  const THit &h = TEv::Ref().vHit(ihit);
  if (updateTID)
    const_cast<THit&>(h).eraseTrackID(Id); // Cast away constness and erase hit->track reference
  NHits--;
  NDFs -= (TSetup::Ref().vPlane(ipl).IFlag&0x30) ? 2 : 1;

  return ihit; // Return index of erased hit
}

void TTrack::RestoreHit(int ipl, int ihit, bool updateTID)
{
  // - Insert hit on plane "ipl".
  // - Upon option, refrain from updating the hit->track reference.
  // - If hit is first(last), insert plane references up(down) to next hit.
  //  => TO be used in conjonction w/ "CancelHit_Clip".
  if      (ipl>lPlnRef.back()) {
    int back = lPlnRef.back();
    list<int>::reverse_iterator ipr = lPlnRef.rbegin(), ihp = lHitPat.rbegin(); 
    lPlnRef.insert(ipr.base(),ipl); lHitPat.insert(ihp.base(),ihit); // Insert new hit and plane ref.
    ipr++; ihp++;
    int jpl = ipl-1; while (jpl>back) {
      lPlnRef.insert(ipr.base(),jpl--); lHitPat.insert(ihp.base(),-1);
      ipr++; ihp++;
    }
  }
  else if (ipl<lPlnRef.front()) {
    int front = lPlnRef.front();
    list<int>::iterator ipr = lPlnRef.begin(), ihp = lHitPat.begin(); 
    lPlnRef.insert(ipr,ipl); lHitPat.insert(ihp,ihit); // Insert new hit and plane ref.
    int jpl = ipl+1; while (jpl<front) {
      lPlnRef.insert(ipr,jpl++); lHitPat.insert(ihp,-1);
    }
  }
  else {
    list<int>::iterator ipr = lPlnRef.begin(), ihp = lHitPat.begin(); 
    while (ipr!=lPlnRef.end() && *ipr<ipl) { ipr++; ihp++; }
    if (ipr==lPlnRef.end() || *ipr!=ipl) {
      CsErrLog::msg(elError,__FILE__,__LINE__,
		    "Trying to restore hit on missing plane #%d in track ID=%d",
		    ipl,Id); return;
    }
    *ihp = ihit;
  }

  const THit &h = TEv::Ref().vHit(ihit);
  if (updateTID)
    const_cast<THit&>(h).addTrackID(Id); // Cast away constness and add hit->track reference
  NHits++;
  NDFs += (TSetup::Ref().vPlane(ipl).IFlag&0x30) ? 2 : 1;
  return;
}

void TTrack::Shorten(int ipl)
{
  TEv &ev = TEv::Ref(); const TSetup &setup = TSetup::Ref();
  list<int>::iterator ipr = lPlnRef.end(); ipr--;
  list<int>::iterator ihp = lHitPat.end(); ihp--;
  while (*ipr>=ipl) {
    // Erase all ref's down to TPlane ipl (included)
    lPlnRef.erase(ipr--);      // Erase plane ref.
    if (*ihp>=0) {
      NHits--; 
      const THit &h = ev.vHit(*ihp);
      const_cast<THit&>(h).eraseTrackID(Id); // Cast away constness and erase hit->track reference
      NDFs -= (setup.vPlane(h.IPlane).IFlag&0x30) ? 2 : 1;
    }
    lHitPat.erase(ihp--);      // Erase hit ref.
  }
  while (ipr!=lPlnRef.begin()) {
    // Erase all ref's down to last (==first encountered) THit (excluded)
    if (*ihp>=0) break;
    else {
      lPlnRef.erase(ipr--);      // Erase plane ref.
      lHitPat.erase(ihp--);      // Erase hit ref.
    }
  }
  int jpl = *ipr; const vector<int> &iplFirst = setup.vIplFirst();
  for (int zone = 0; zone<int(iplFirst.size()); zone++) {
    if (jpl>=iplFirst[zone]) continue;
    Type &= ~(1<<zone);
  }
}

void TTrack::Behead(int ipl)
{
  TEv &ev = TEv::Ref(); const TSetup &setup = TSetup::Ref();
  list<int>::iterator ipr = lPlnRef.begin();
  list<int>::iterator ihp = lHitPat.begin();
  while (*ipr<=ipl) {
    // Erase all ref's up to TPlane ipl (included)
    lPlnRef.erase(ipr++);      // erase plane ref.
    if (*ihp>=0) {
      NHits--;
      const THit &h = ev.vHit(*ihp);
      const_cast<THit&>(h).eraseTrackID(Id); // Cast away constness and erase hit->track reference
      NDFs -= (setup.vPlane(h.IPlane).IFlag&0x30) ? 2 : 1;
    }
    lHitPat.erase(ihp++);      // erase hit ref.
  }
  while (ipr!=lPlnRef.end()) {
    // Erase all ref's up to first (==first encountered) THit (excluded)
    if (*ihp>=0) break;
    else {
      lPlnRef.erase(ipr++);      // erase plane ref.
      lHitPat.erase(ihp++);      // erase hit ref.
    }
  }
  int jpl = *ipr; const vector<int> &iplLast = setup.vIplLast();
  for (int zone = 0; zone<int(iplLast.size()); zone++) {
    if (jpl<=iplLast[zone]) continue;
    Type &= ~(1<<zone);
  }
}

void TTrack::Clip(bool downstream)
{
  // Erase planes inserted up(resp. down)stream of first (resp. last) hit
  if (downstream) {
    list<int>::reverse_iterator ipr = lPlnRef.rbegin(), ihp = lHitPat.rbegin();
    while (ipr!=lPlnRef.rend()) {
      // Erase all ref's down to last (==first encountered) THit (excluded)
      if (*ihp>=0) break;
      else {
	lPlnRef.erase(--ipr.base()); lHitPat.erase(--ihp.base());// Erase hit and plane ref.
      }
    }
  }
  else {
    list<int>::iterator ipr = lPlnRef.begin();
    list<int>::iterator ihp = lHitPat.begin();
    while (ipr!=lPlnRef.end()) {
      // Erase all ref's up to first (==first encountered) THit (excluded)
      if (*ihp>=0) break;
      else {
	lPlnRef.erase(ipr++);      // erase plane ref.
	lHitPat.erase(ihp++);      // erase hit ref.
      }
    }
  }
}

bool TTrack::BendingInfo()
{
  // Determine whether there are enough hits in the bending dimension,
  // meaning all of non Z projections, upstream of SM2 to provide precise
  // enough info to determine bending.
  const TSetup &setup = TSetup::Ref();
  list<int>::iterator ipr = lPlnRef.begin();
  list<int>::iterator ihp = lHitPat.begin();
  const TStation *sPrv = 0; int projs =0, nBendingProjs = 0;
  while (ipr!=lPlnRef.end()) {
    if (*ihp>=0) {
      int ipl = *ipr; if (ipl<=setup.vIplLast()[1]) {
	const TPlane &p = setup.vPlane(ipl);
	const TStation *&s = p.Station; if (s!=sPrv) {
	  if ((projs&0x2)!=projs)  // Z proj. is 0x2, cf. "TSetup::Init"
	    nBendingProjs++;
	  sPrv = s; projs = 0;
	}
	projs |= 1<<p.IProj;
      }
    }
    ipr++; ihp++;
  }
  if ((projs&0x2)!=projs) nBendingProjs++;
  return nBendingProjs>1;
}

int TTrack::SingleDiff(const TTrack *t, int &typeDiff) const
{
  // - Check that the diff w.r.t. argument TTrack concerns a single plane.
  // - If indeed:
  //   - Return that plane's index.
  //   - Return also the type of diff: either missing, different or extra.
  // - Else:
  //   - Return -1.
  //   - "typeDiff" meaningless, or even undefined.
  list<int>::const_iterator ipr1 = t->lPlnRef.begin(), ipr2 = lPlnRef.begin();
  list<int>::const_iterator ihp1 = t->lHitPat.begin(), ihp2 = lHitPat.begin();
  int ipl = -1;
  while (ipr1!=t->lPlnRef.end()) {
    while (*ipr2<*ipr1 && ipr2!=lPlnRef.end()) {
      if (*ihp2>=0) {
	if (ipl>=0) return -1; ipl = *ipr2; typeDiff = 2;
      }
      ipr2++; ihp2++;
    }
    if      (*ipr2!=*ipr1 || *ihp1>=0 && *ihp2<0) {
      if (ipl>=0) return -1; ipl = *ipr1; typeDiff = 0;
    }
    else if (*ihp1<0 && *ihp2>=0) {
      if (ipl>=0) return -1; ipl = *ipr1; typeDiff = 2;
    }
    else if (*ihp2!=*ihp1) {
      if (ipl>=0) return -1; ipl = *ipr1; typeDiff = 1;
    }
    ipr1++; ihp1++;
    if (ipr2!=lPlnRef.end()) { ipr2++; ihp2++; }
  }
  while (ipr2!=lPlnRef.end()) {
    if (*ihp2>=0) {
      if (ipl>=0) return -1; ipl = *ipr2; typeDiff = 2;
    }
    ipr2++; ihp2++;
  }
  return ipl;
}
