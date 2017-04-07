// $Id: TTrack.cc 14069 2015-09-17 20:44:46Z lsilva $

#include "TTrack.h"
#include "TEv.h"
#include "THit.h"

using namespace std;

/*
  Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/TTrack.cc":
  i) Member object "IFit" of class "Ttrack".
  ii) Member objects "Haux" (auxilliary helix), "vGuests" (guestimates) and
  "GuestsGroups".
  iii) Method "CrossRPipe".
*/

// Constructor
TTrack::TTrack():
Id(++TTrack::TrackCounter),
Type   (0),
Scifi  (0),
IMark  (0),
IFit   (0),
fFitDone(false),
bFitDone(false),
insPlaneDone(false),
IKine (-1),
NHits  (0),
NDFs   (0),
NDics  (0),
NHsame (0),
IKsame(-1),
Chi2tot(9999),
MeanTime(9999),
SigmaTime(-1),
DispTime(-1),
radLenFr(-1),
ptrTrk(NULL),
NShowerHits (0),
HasShower (0),
Associate (-1),

Hfirst     (), 
Hlast      (), 
Haux       (),

sProj      (), 
lPlnRef    (), 
lHitPat    (),
vGuests    (),
GuestsGroups(0),
Chi2aux(0)
{
}


/*!

  Copy track "t" to "this" track and
  establishes corresponding THit -> "this" track references
  and TKine -> "this" track reference (in case of MC event)

*/
TTrack::TTrack(const TTrack& t):
Id(++TTrack::TrackCounter),
Type     (t.Type),
Scifi    (t.Scifi),
IMark    (t.IMark),
IFit     (t.IFit),
fFitDone (t.fFitDone),
bFitDone (t.bFitDone),
insPlaneDone(t.insPlaneDone),
IKine    (t.IKine),
NHits    (t.NHits),
NDFs     (t.NDFs),
NDics    (t.NDics),
NHsame   (t.NHsame),
IKsame   (t.IKsame),
Chi2tot  (t.Chi2tot),
MeanTime (t.MeanTime),
SigmaTime(t.SigmaTime),
DispTime (t.DispTime),
radLenFr (t.radLenFr),
ptrTrk   (t.ptrTrk),
NShowerHits (t.NShowerHits),
HasShower (t.HasShower),
Associate (t.Associate),

Hfirst   (t.Hfirst),
Hlast    (t.Hlast),
Haux     (t.Haux),

sProj    (t.sProj),
lPlnRef  (t.lPlnRef),
lHitPat  (t.lHitPat),
vGuests  (t.vGuests),
GuestsGroups(t.GuestsGroups),
Chi2aux(t.Chi2aux)
{


  TEv* pEv = TEv::Ptr(); // pointer to TEv object
  if( pEv == NULL) {
    cout<<"TTrack::TTrack(const TTrack&) ==> TEv object must exists!"<<endl; 
    assert(false);
  }
  
  // Set association THit -> TTrack
  
  list<int>::iterator i;
  for(i = lHitPat.begin(); i != lHitPat.end(); i++){     // loop over hit references on the track
    if((*i) < 0) continue; // missing hit
    // cast away constness as the state of the THit object will be changed
    THit&    h = const_cast<THit&> (pEv->vHit(*i)); 
    // add corresponding back reference THit->TTrack
    h.addTrackID(this->Id);
  } // end of loop over hit references on the track
  
  
  // Set corresponding TKine -> this TTrack association 
  
  if (this->IKine >= 0) { // if there is association with MC track
    // cast away constness as the state of the TKine object will be changed
    TKine& k = const_cast<TKine&> (pEv->vKine(this->IKine)); 
    k.addTrackID(this->Id);  // add reference to the track
  }
  
  
}


/*!

 \warning This "light-weight" assignment operator do the fast copy of track contents 
 (including "unique"track ID :-)
 and DO NOT establishes corresponding THit -> "this" track references
 and TKine -> "this" track reference. To be used ONLY for saving 
 of track information into temporary track objects 

*/

TTrack& TTrack::operator = (const TTrack& t)
{
  Id          = t.Id;
  Type        = t.Type;
  Scifi       = t.Scifi;
  IMark       = t.IMark;
  IFit        = t.IFit;
  fFitDone    = t.fFitDone;
  bFitDone    = t.bFitDone;
  insPlaneDone = t.insPlaneDone;
  IKine       = t.IKine;
  NHits       = t.NHits;
  NDFs        = t.NDFs;
  NDics       = t.NDics;
  NHsame      = t.NHsame;
  IKsame      = t.IKsame;
  Chi2tot     = t.Chi2tot;
  MeanTime    = t.MeanTime;
  SigmaTime   = t.SigmaTime;
  DispTime    = t.DispTime;
  radLenFr    = t.radLenFr;
  NShowerHits = t.NShowerHits;
  HasShower   = t.HasShower;
  Associate   = t.Associate;

  Hfirst      = t.Hfirst;
  Hlast       = t.Hlast;

  sProj       = t.sProj;
  lPlnRef     = t.lPlnRef;
  lHitPat     = t.lHitPat;

  return (*this);

}

/*!

  Delete "this" track and remove
  corresponding THit -> "this" track references
  and TKine -> "this" track reference (in case of MC event)

*/

TTrack::~TTrack()
{
  TEv* pEv = TEv::Ptr(); // pointer to TEv object
  if( pEv != NULL) {

    // Remove assosiation THit -> TTrack
    
    list<int>::iterator i;
    for(i = lHitPat.begin(); i != lHitPat.end(); i++){     // loop over hit references on the track
      if((*i) < 0) continue; // missing hit
      // cast away constness as the state of the THit object will be changed
      THit&    h = const_cast<THit&> (pEv->vHit(*i)); 
      // erase corresponding back reference THit->TTrack
      h.eraseTrackID(this->Id);
    } // end of loop over hit referencies on the track
    
    // Remove corresponding TKine -> this TTrack association 
    if (this->IKine >= 0) { // if there is association with MC track
      // cast away constness as the state of the TKine object will be changed
      TKine& k = const_cast<TKine&> (pEv->vKine(this->IKine)); 
      k.eraseTrackID(this->Id);  // Erase reference to the track
    }
    
  } 
  
}

// ===========================================================================
// ==============================    Evaluate   ==============================
// ===========================================================================
#include "CsGEMDetector.h"
#include "CsPixelGEMDetector.h"
void TTrack::Evaluate(int ipli, int iplf, unsigned int mode, 
		      int &nSpacePts, int &nAllProjs)
{
  // - # of spacePts and the # of proj.
  // - In track's segment comprised w/in "[ipli,iplf]".
  // - One space point per set of TStation's providing at least 3 proj.
  // - Upon option, one space point per:
  //   - (mode&0x1) Scifi TStation of 2 proj., provided time consistency.
  //   - (mode&0x2|4) GM|P XY or UV TStation, provided amplitude correlation.
  int fi = (mode&0x1) ? 22 : 0;
  int gm = (mode&0x2) ? 26 : 0;
  int gp = (mode&0x4) ? 28 : 0;

  const TSetup &setup = TSetup::Ref(); const TEv &ev = TEv::Ref();

  //            ***** INITIALIZATION *****
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
    // Amplitude correlation of GM|PXY|UVs
    const vector<TStation> &stations = setup.vStation();
    for (int is = 0; is<(int)stations.size(); is++) {
      const TStation &s = stations[is]; if (s.Type!=26 && s.Type!=28) continue;
      int id = ids.size(); ids[is] = id;
      int jpl = 0;
      const TPlane &pU = setup.vPlane(s.IPlanes[0]);
      if (pU.IProj==iProjU) {
	const TDetect &dU = setup.vDetect(pU.IDetRef);
	const float *ampCorr; if (s.Type==26) {
	  CsGEMDetector *g = dynamic_cast<CsGEMDetector*>(dU.PtrDet());
	  if (!g) CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "TStation #%d(Type=26)'s TPlane, \"%s\", cannot be cast into a CsGEM",
			is,dU.Name.c_str());
	  ampCorr = g->getAmpCorr();
	}
	else {
	  CsPixelGEMDetector *g = dynamic_cast<CsPixelGEMDetector*>(dU.PtrDet());
	  if (!g) CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "TStation #%d(Type=28)'s TPlane, \"%s\", cannot be cast into a CsPixelGEM",
			is,dU.Name.c_str());
	  ampCorr = g->getAmpCorr();
	}
	aUV0s.push_back(ampCorr[0]); aUV1s.push_back(ampCorr[1]);
	aUV2s.push_back(ampCorr[2]);
	jpl = 2;
      }
      else { // Can be the case of a lone CsPixelGEM, e.g. GP01XY in 2008,9.
	aUV0s.push_back(0); aUV1s.push_back(0); aUV2s.push_back(0);
      }
      const TPlane &pY = setup.vPlane(s.IPlanes[jpl]);
      if (pY.IProj==iProjY) {
	const TDetect &dY = setup.vDetect(pY.IDetRef);
	const float *ampCorr; if (s.Type==26) {
	  CsGEMDetector *g = dynamic_cast<CsGEMDetector*>(dY.PtrDet());
	  if (!g) CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "TStation #%d(Type=26)'s TPlane, \"%s\", cannot be cast into a CsGEM",
			is,dY.Name.c_str());
	  ampCorr = g->getAmpCorr();
	}
	else {
	  CsPixelGEMDetector *g = dynamic_cast<CsPixelGEMDetector*>(dY.PtrDet());
	  if (!g) CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "TStation #%d(Type=28)'s TPlane, \"%s\", cannot be cast into a CsPixelGEM",
			is,dY.Name.c_str());
	  ampCorr = g->getAmpCorr();
	}
	aYX0s.push_back(ampCorr[0]); aYX1s.push_back(ampCorr[1]);
	aYX2s.push_back(ampCorr[2]);
      }
      else {
	aYX0s.push_back(0); aYX1s.push_back(0); aYX2s.push_back(0);
      }
    }
  }

  list<int>::iterator ih; const TStation *sPrv;
  int nPls, mPls; double data[4]; int nProjs; unsigned int projs, allProjs;
  for (ih = lHitPat.begin(), nSpacePts=nProjs=nAllProjs=nPls=mPls = 0,
	 projs=allProjs = 0, sPrv = 0; ih!=lHitPat.end(); ih++) {
    if (*ih<0) continue;
    const THit &h = ev.vHit(*ih); int ipl = h.IPlane;
    if (ipl<ipli || iplf<ipl) continue;
    const TPlane &p = setup.vPlane(ipl);
    const TStation *&s = p.Station; if (s!=sPrv) {
      if (sPrv) {
	if      (sPrv->Type==29) {             // Pixels (GP..P) station...
	  // ...One SpacePt per pixel (note: one pixel plane <-> one station)
	  if (nProjs>=2 /* It's either 0 or 2 */) {
	    nSpacePts++; nProjs = 0; projs = 0;
	  }
	}
	else if (sPrv->Type==gm ||             // GM station...
		 sPrv->Type==gp) {             // Strips (GP..XYUV) station...
	  // Let's try to build 1 SpacePt per detector, by relying on the
	  // amplitude correlation.
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
	    double deltaYX = aYX0s[id]+aYX1s[id]*aX+aUV2s[id]*aX*aX-aY;
	    if (fabs(deltaYX)<TOpt::ReMode[17]) {
	      nSpacePts++; nProjs = 0; projs = 0;
	    }
	  }
	  nPls=mPls = 0;
	}
	else if (sPrv->Type==fi) {             // Scifis station...
	  // Check the time correlation.
	  if (nPls>=2 && fabs(data[0]-data[1])<2) {
	    nSpacePts++; nProjs = 0; projs = 0;
	  }
	  nPls = 0;
	}
	if (nProjs>=3) { nSpacePts++; nProjs = 0; projs = 0; }
      }
      sPrv = s;
    }
    unsigned int proj = 1<<p.IProj;
    if      (s->Type==fi) data[nPls++] = h.Time;
    else if (s->Type==gm || s->Type==gp) {
      if (proj&uvProj) data[nPls++] = h.PtrClus()->getAllAnalogData()[2];
      else             data[2+mPls++] = h.PtrClus()->getAllAnalogData()[2];
    }
    if ((proj&projs)==0)        {    projs |= proj;    nProjs++; }
    if ((proj&allProjs)==0)     { allProjs |= proj; nAllProjs++; }
    if (p.IFlag&0x30) { // P-pixel planes: 2nd projection.
      if   (p.IFlag&0x10) proj ^= yzProj;
      else                proj ^= uvProj;
      if ((proj&projs)==0)      {    projs |= proj;    nProjs++; }
      if ((proj&allProjs)==0)   { allProjs |= proj; nAllProjs++; }
    }
  } // End loop on hits
  if (sPrv) { // Account for last TStation
    if      (sPrv->Type==29) {                 // Pixels (GP..P) station
      if (nProjs>=2) { nSpacePts++; nProjs = 0; }
    }
    else if (sPrv->Type==gm ||                 // GM station...
	     sPrv->Type==gp) {                 // Strips (GP..XYUV) station...
      int id = ids[sPrv->IStation];
      if (nPls==2) {
	double aU = data[0], aV = data[1];
	double deltaUV = aUV0s[id]+aUV1s[id]*aV+aUV2s[id]*aV*aV-aU;
	if (fabs(deltaUV)<TOpt::ReMode[17]) { nSpacePts++; nProjs = 0; }
      }
      if (mPls==2) {
	double aY = data[2], aX = data[3];
	double deltaYX = aYX0s[id]+aYX1s[id]*aX+aUV2s[id]*aX*aX-aY;
	if (fabs(deltaYX)<TOpt::ReMode[17]) { nSpacePts++; nProjs = 0; }
      }
    }
    else if (sPrv->Type==fi) {                 // Scifis station
      if (nPls==2 && fabs(data[0]-data[1])<2) { nSpacePts++; nProjs = 0; }
    }
    if (nProjs>=3) nSpacePts++;
  }
}

bool TTrack::CrossRPipe() {
  // Return whether track crosses RICH pipe, as determined by straight extrapolation from its first helix.
  if (TOpt::RICHPipe[5]!=4) // Cf. "RCH1HOLE" volume in "$COMGeant/data/geom/geom_general.ffr"
    // 4 corresponds to the original, steel, pipe...
    // 6 corresponds to the light Al pipe installed in 2012: let's disregard it
    // else no info available or type code is unforeseen: let's also disregard the case
    return false;

  double *p = TOpt::RICHPipe, l = p[3], r = p[4], r2 = r*r; int io[2];
  const THlx &Hf = H('u');
  float x0 = Hf(0), y0 = Hf(1), z0 = Hf(2), yp = Hf(3), zp = Hf(4);
  for (int ud = 0; ud<2; ud++) {
    double dx = p[0]+(2*ud-1)*l-x0, y = y0+dx*yp-p[1], z = z0+dx*zp-p[2];
    io[ud] = y*y+z*z<r2 ? -1 : 1;
  }
  if (io[0]*io[1]>0) return false;
  else return true;
}

void TTrack::WorstHit(map<int,double> &mChi2,
		      int &worstPl, int &worstGr, double &worstChi2Incr)
{
  // Return plane#, zone# and value of worst chi2 increment 
  const TSetup &setup = TSetup::Ref();

  int jgr, kgr;
  static int igrs[6] = {0,1,2,3,4,5}, ngrs;
  static bool first = true; if (first) {  // ***** INIT *****
    // Reorder the spectrometer zones in increasing plane index order
    first = false; ngrs = (int)setup.vIplFirst().size();
    if (ngrs>6) CsErrLog::mes(elFatal,"More than 5 zones!");
    for (jgr = 0; jgr<ngrs; jgr++) for (kgr = jgr+1; kgr<ngrs; kgr++) {
      int lgr = igrs[jgr], mgr = igrs[kgr];
      if (setup.vIplFirst()[mgr]<setup.vIplFirst()[lgr]) {
	igrs[jgr] = mgr; igrs[kgr] = lgr;
      }
    }
  }
  int ihit; map<int,double>::iterator idx;
  for (ihit = 0, jgr = 0, worstChi2Incr = 0, idx = mChi2.begin();
       ihit<(int)NHits; ihit++, idx++) {
    //      ***** WORST HIT? in BACK/FOREWARD KF FIT *****
    double chi2Incr = (*idx).second; int ipl = (*idx).first; 
    while (ipl>setup.vIplLast()[igrs[jgr]] && jgr<ngrs) jgr++;
    if (jgr==6) CsErrLog::msg(elError,__FILE__,__LINE__,
       "Track #%d worst plane has unphysical index = %d!",Id,ipl);
    if (chi2Incr>worstChi2Incr) {
      worstPl = ipl; worstGr = igrs[jgr]; worstChi2Incr = chi2Incr;
    }
  }
}
