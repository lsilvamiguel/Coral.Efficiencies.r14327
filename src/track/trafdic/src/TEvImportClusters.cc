// $Id: TEvImportClusters.cc 14069 2015-09-17 20:44:46Z lsilva $

#include "TH1.h"
#include "CsErrLog.h"
#include "CsHistograms.h"
#include "CsEvent.h"
#include "CsEventUtils.h"
#include "CsMCDigit.h"
#include "CsMCHit.h"
#include "CsMCTrkHit.h"
#include "CsMCTrack.h"
#include "CsMCUtils.h"
#include "TEv.h"
#include "TDigit.h"
#include "TSetup.h"
#include "TOpt.h"
#ifdef COMPASS_USE_OSPACE_STD
# include <ospace/std/set>
#else
# include <set>
#endif

using namespace std;
using namespace CLHEP;

/*! 
  \brief Import CsClusters into Traffic Event TEv.
*/

/*
  Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/TEvImportClusters":
    i) Set the ".Status" and ".Mirror" member data of THit objects.
   ii) Move in LR Filter: THit's having CsCluster counterpart w/ low LR proba
      are tagged w/ ".Status<-1".
  iii) LR pairs with close components are coalesced into a single hit w/ average
      position and enlarged uncertainty (mirror THit still created, w/ same
      caracteristics but ".Status=-4").
   iv) Maintain MC-track/hit association when mirror, by setting "iKine" member
      data = "-2-jKine", where "jKine" is "iKine" of genuine hit.
    v) 2-MC-track hits: "iKine" and "iKin2".  (Note that this design is not
      foolproof; it does not distinguish between the case where there two
      particles contributing, the 2nd of which being ill-defined and the case of
      a single particle: it yields in both case "IKine>=0, IKin2==-1".)
   vi) Special handling of GEM ``spacer hits''.
  vii) Various SI specifics: time cut, enlargement of time and position uncertainties.
 viii) Pixel hits. 
   ix) Two-step import, w/ beam telescope imported alone in first step.
*/

void TEv::ImportClusters(const list<CsCluster*> &lClust,
			 int beamTelescope) // =1: Beam only, =-1: All but beam, =0: Everything
{
  const TSetup& setup = TSetup::Ref();

  bool hist = TOpt::Hist[2] && isMC;

  static TH1D *h[5];
  if (hist) {
    static bool first = true;
    if (first) {
      first=false;
      CsHistograms::SetCurrentPath("/Traffic/import_clusters");
      h[0] = new TH1D("h1","Total # of CsClusters",    100, 0, 5000);
      h[1] = new TH1D("h3","# of !associated clusters",100, 0, 100);
      h[2] = new TH1D("h5","MCTrack mixing in cluster (sMCtrack size)",20,0,20);
      h[3] = new TH1D("h6","Iorig flag mixing in cluster (sIorig size)",
			  20, 0, 20);
      h[4] = new TH1D("h7","# of MCHits in cluster (sMChit size)",     20,0,20);
    }
    CsHistograms::SetCurrentPath("/");
    if (beamTelescope<=0) // Avoid double counting: fill histo only in 2nd step
      h[0]->Fill(lClust.size()+0.5);
  }

  //              ***** BEAM TELESCOPE *****
  // Upon physics data (i.e. 5 zones are defined, of which 5th one is beam
  // telescope), the clusters are imported in 2 steps:
  //  - First step: the beam telescope is imported alone, and processed. It can
  //   can then be taken into account in the definition of an event time that
  //   may differ significantly from the trigger time.
  //  - Second step, after CsClusters from drift detectors have been possibly
  //   updated w/ the newly defined event time: all other zones are imported.
  // The implementation of this two-step scheme assumes that the beam telescope
  // comes first in the list of detectors.
  static int nClsImported_m1;   // # of imported clusters !associated
  int nZones = (int)setup.vIplFirst().size(), iplMn, iplMx; if (beamTelescope) {
    if (nZones==6) {
      if (!TOpt::ReMode[49])
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
	  "Two-step import requested while a 6th reconstruction zone defined"
	  " and yet Drell-Yan \"ReMode[49]\"not specified");
    }
    else if (nZones!=5) 
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
	"Two-step import requested while only %d reconstruction zones defined",
		    nZones);
    if      (beamTelescope==1) {
      iplMn = setup.vIplFirst()[4]; iplMx = setup.vIplLast()[4];
      nClsImported_m1 = 0;
    }
    else {
      if (beamTelescope!=-1)
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
	  "Inconsistency: argument \"beamTelescope\" = %d",beamTelescope);
      if (nZones==6)  // Special Drell-Yann case
	// The 6th zone is the upstream-most, beam telescope aside.
	iplMn = setup.vIplFirst()[5];
      else
	iplMn = setup.vIplFirst()[0];
      iplMx = setup.vIplLast()[3];
    }
  }
  else {
    iplMn = 0; int zone; for (zone = 0,iplMx = 0; zone<nZones; zone++) {
      int iplLast = setup.vIplLast()[zone]; if (iplLast>iplMx) iplMx = iplLast;
    }
    nClsImported_m1 = 0;
  }

  //   ***** SPECIFICS of SOME DETEECTOR TYPES: GEMs and Sis
  bool gemSpacers = setup.GEMSpacers();
  double siTSigCut = TOpt::dCut[77], siTCut = TOpt::dCut[78];
  double siTSigScale = TOpt::dCut[79] ? TOpt::dCut[79] : 1;
  int siCut = (siTSigCut?0x1:0) | (siTCut?0x2:0);
  double extraTime = // If current event's trigger jitter is larger than reference's
    ptrEvt()->getExtraTimeWidth();
  siTCut += extraTime; // In case "siTCut" is enabled, it will use this enlarged cut
  double siTT0     = // Wrong master trigger was selected by Decoding library
    ptrEvt()->getTriggerTimeCorr() +
    TOpt::dCut[80];  // T0 offset (Mainly used to catch off-time tracks in alignment procedure. Could also make for T0 mis-calibration.)
  double siSigCorr = // Position uncertainties correction term: for when ...
    TOpt::dCut[85];  // ...badly aligned, e.g.

  CsCluster *mirror, *prvCluster; int prvPlane;
  list<CsCluster*>::const_iterator ic;
  for (ic = lClust.begin(), prvCluster = 0, prvPlane = -1; ic!=lClust.end();
       ic++) { 

    // *****************************************************************
    // ******************** MAIN LOOP ON CsClusters ********************
    // *****************************************************************
    CsCluster* clus = *ic;

    mirror = clus->getMirrorCluster();  // ********** MIRRORS **********

    // Check if the cluster already was imported
    map<CsCluster*, int, less<CsCluster*> >::const_iterator im;
    im = mapCsCl2Hit.find(clus);
    if(im != mapCsCl2Hit.end()) continue; // already there, skip it

    // ******************** CsCluster -> TPlane -> TDetect ********************

    list<CsDetector*> lDet = clus->getDetsList();
    if (lDet.empty()) {
      CsErrLog::mes(elError,"CsCluster has no reference to CsDetectors"); break;
    }
    if (lDet.size()>1) {
      CsErrLog::mes(elError,"CsCluster belonging to more then 1 detector");
      continue;
    }
    CsDetector *csDet = lDet.front();
    int iPlane = setup.Id2iPlane(csDet->GetID());
    if (iPlane<iplMn || iplMx<iPlane) continue;
    const TDetect &det = setup.iPlane2Detect(iPlane);

    double u = clus->getU()/10;
    HepMatrix m = clus->getCov(); double sigu = sqrt(m(1,1))/10.; 

    int status = 0; if (det.IType==26 && gemSpacers) {
      //          *************** GEM SPACERS ***************
      static const float *gemBounds; if (det.IPlane!=prvPlane)
	gemBounds = setup.GEMSpacerBounds(det.IDet);
      int i, winBounds; for (i=winBounds = 0 ; i<3; i++) {
	int k = 2*i, l = k+1;
	if (u<gemBounds[k]) break;
	if (u<gemBounds[l]) {
	  winBounds++; k += 6; l += 6; 
	  if (gemBounds[k]<u && u<gemBounds[l]) winBounds++;
	}
      }
      if (winBounds==2) continue; // Cluster w/in inner spacer interval -> Skip 
      if (winBounds==1) {         // Cluster w/in outer spacer interval -> ...
	sigu += TOpt::dCut[74]; status = -1;  // ...Enlarge uncertainty and Flag
	if (TOpt::ReMode[29]&0x4) clus->setSigmaU(sigu*10);
      }
    }

    //         *************** SI ***************
    // It's placed at this point of the flowchart, so that time cut can reject
    // current cluster as early as possible.
    if (det.IType==21) {
      if (iPlane<setup.vIplLast()[0])                          // ***** POSITION
	sigu += siSigCorr; // (Only those of zone 0x1)
      double time, sigt;                                           // ***** TIME
      if (!clus->getTime(time) || // These cases should never occur...
	  !clus->getTimeError(sigt)) continue;
      if (time>999)             // ***** TIMELESS HIT (they are assigned t=1000)
	// Complete uncertainty => keep but...
	status = -2; // ...to be used only in later iterations of "PrePattern2"
      else {
	if ((siCut&0x1) && sigt>0 &&
	    fabs(time-siTT0)>siTSigCut*sigt+extraTime)  // ***** OUT OF TIME HIT
	  continue;	        // => Can be safely discarded
	else if ((siCut&0x2) && fabs(time-siTT0)>siTCut)
	  // ***** LARGE TIME & LARGE TIME UNCERTAINTY
	  // In the so-called "cind" clusterisation mode: there are many hits w/
	  // time in the 10 ns range and beyond but still compatible w/ 0 due to
	  // their large uncertainty. Let's discard them, assuming uncertainty
	  // is overestimated. This assumption is backed by comparison w/ the
	  // "ratio" clusterisation mode (where many of these hits get rejected)
	  // and by the time distribution of hits associated to tracks when the
	  // absolute time cut is relaxed.
	  continue;
      }
    }

    //       ******************** MC ONLY BLOCK ********************

    list<CsDigit*> lDig = clus->getDigitsList(); list<CsDigit*>::iterator id;
    if (lDig.empty()) {
      CsErrLog::mes(elError,"CsCluster has empty list<CsDigit*>"); continue;
    }

    int ikine = -1, ikin2 = -1, iorig = -1;
    // "jkine" = "ikine" OR "-2-ikine" if ikine <-1
    static int jkine;  // "static" to avoid "warning... uninitialized"
    
    set<CsMCHit*,less<CsMCHit*> > sMChit;// MCHits associated to current cluster
    set<CsMCHit*, less<CsMCHit*> >::iterator iMChit;

    if (isMC) {

      // *************** HITS: FIND POINTER->CsMCTrack & ORIGIN ***************

      // Case 1 ORIGINAL component:  Fill "THit.iKine", iOrig=0
      // Case 1 or 2 ORIGINAL components mixed w/ NON-ORIGINAL hits:
      //   - Check whether CsCluster position info is consistent (4 sigmas)
      //    with original trajectory.
      //   - If indeed, fill "THit.iKine", AND, possibly, "THit.iKin2", iOrig=0
      // Case 1 or more NON-ORIGINAL from a unique MC track: Fill "THit.iKine",
      //  iOrig = this->Origin if unique, else = 1;
      // Other cases: No TKine, iOrig=-1
      // (The rationale backing this choice:
      //  - Try to be as fair as possible when it comes to evaluate
      //   tracking efficiency
      //  => Assign TKine (i.e. declare composite THit as good hit) as often as
      //    possible: Composite are not of the same quality as singletons
      //    but they still convey meaningful info on track position and it is
      //    often too demanding to require tracking to distinguish them from the
      //    latter.
      // N.B.: The problem is reevaluated later on, for mirror hits.)

      float c = det.Ca, s = det.Sa;

      set<CsMCTrack*,less<CsMCTrack*> > sMCtrack;
      set<int,less<int> > sIorig; 
      // (These are sets because we want to count only once like entries, for
      // purpose of statistics (otherwise only part of the info is later used).)
      
      // Plain array for assigning TKine and Origin according to rules supra.
      const int nRefsMx = 3; CsMCTrack *kineRefs[nRefsMx];
      int nKineRefs, iHOrigs[nRefsMx], nHs[nRefsMx]; float uHs[nRefsMx];
      for (id = lDig.begin(), nKineRefs = 0; id!=lDig.end(); id++) {

	// ***** LOOP OVER CsDigits IN THE CsCluster *****

	list<CsMCHit*> lHit; list<CsMCHit*>::iterator ih;
	CsMCDigit* mcd = dynamic_cast<CsMCDigit*>(*id);
	if(mcd == NULL){
	  CsErrLog::mes(elError,"Cannot cast CsDigit* to CsMCHit* (MC)"); break;
	}
	else lHit = mcd->getHits();
	if(lHit.empty()){
	  CsErrLog::mes(elError,"CsMCDigit has empty list<CsMCHit*>"); break;
	}
	for (ih = lHit.begin(); ih!=lHit.end(); ih++) { 

	  // ***** LOOP OVER MCHits OF THE DIGIT *****

	  if (*ih==NULL) {
	    CsErrLog::mes(elError,"CsMCDigit w/ NULL pointer to CsMCHit");break;
	  }
	  if ((*ih)->getMCTrack()==NULL) {
	    CsErrLog::mes(elError,"CsMCHit w/ NULL pointer to CsMCTrack");break;
	  }
	  sMChit.insert(*ih);                      // ***** FILL sets
	  CsMCTrack *pMCTr = (*ih)->getMCTrack(); sMCtrack.insert(pMCTr);

	  int iHOrig = (*ih)->getOrigin();  sIorig.  insert(iHOrig);
	  int iRef, match; for (iRef = 0, match = -1; iRef<nKineRefs; iRef++)
	    if (iRef<nRefsMx && pMCTr==kineRefs[iRef]) match = iRef;
	  if (match==-1 ||   // New entry...
	      iHOrig==0) {   // ...or original contribution to an existing entry
	    if (iHOrig==0) {
	      float uH;
	      CsMCTrkHit *tHit = dynamic_cast<CsMCTrkHit*>(*ih);
	      if (tHit==0) {
		CsErrLog::mes(elError,"CsMCHit w/o CsMCTrkHit"); continue;
	      }
	      else {
		uH = tHit->getX()*c+tHit->getY()*s;
		if   (nKineRefs<nRefsMx && (match==-1 || iHOrigs[match]!=0)) {
		  if (match==-1) {
		    uHs[nKineRefs] = uH; nHs[nKineRefs] = 1;
		    iHOrigs[nKineRefs] = 0;
		  }
		  else {
		    uHs[match] = uH;     nHs[match] = 1;
		    iHOrigs[match] = 0;
		  }
		}
		else {
		  uHs[match] += uH;      nHs[match]++;
		}
	      }
	    }
	    if (match==-1) {
	      if (nKineRefs<nRefsMx) {
		iHOrigs[nKineRefs] = iHOrig; kineRefs[nKineRefs++] = pMCTr;
	      }
	      else nKineRefs++; 
	    }
	  }
	}// end of loop over MC hits of the digit 
      }  // End of loop over MCDigits in the cluster

      // ***** ORIGIN of the CsCluster? / HOW MANY TKine REFERENCES? *****

      int iRef, nOriginals;
      int nRefsUpLimit = // We base our answer on the sole TKine references...
	// ...for which info was collected (into arrays[nRefsMx]). The rest of
	// the info is lost, leading to possibly underestimating the quality of
	// the THit, and hence that of the reco perfs later on.
	nKineRefs<=nRefsMx?nKineRefs:nRefsMx;
      for (iRef=nOriginals = 0, iorig = -2; iRef<nRefsUpLimit; iRef++) {
	if (iHOrigs[iRef]==0) {
	  float uC = uHs[iRef]/nHs[iRef]/10;
	  if (fabs(u-uC)>4*sigu) iHOrigs[iRef] = -1;
	  else nOriginals++;
	  iorig = 0;
	}
	if      (iorig==-2)                        iorig = iHOrigs[iRef];
	else if (iHOrigs[iRef]!=iorig &&
		 iorig!=0)  // Original hit takes precedence
	  iorig = -1;
      }

      if (hist) {
	h[2]->Fill(sMCtrack.size()+0.5); h[3]->Fill(sIorig.size()+0.5);
	h[4]->Fill(sMChit.size()+0.5);
      }

      if (sMCtrack.size()==0)
	CsErrLog::mes(elError,"No CsCluster association to CsMCTrack");
      
      // ********** FIND TKine TRACK NUMBER **********
      // Several such, in fact 1 or 2: 2 tracks being allowed only if no mirror

      if (nKineRefs==1 || nOriginals==1 ||    // ***** ONLY 1 TKine REF
	  nOriginals==2 && !mirror) {         // ***** or 2 IF NO MIRROR
	for (iRef = 0; iRef<nRefsUpLimit; iRef++) {
	  if (nKineRefs==1 || iHOrigs[iRef]==iorig) {
	    CsMCTrack *pMCTr = kineRefs[iRef];
	    for (int it = 0; it < int(vecKine.size()); it++) {
	      if (pMCTr==vecKine[it].PtrTrk()) { 
		if   (iRef==0)
		  ikine = it;
		else // A 2nd TKine ref. is added to the THit. Nota bene that...
		  // ...for this to happen in the case of an original THit (i.e
		  // assigned THit::IOrig=0, and what case only really matters
		  // in the later evaluation of reco perfs where non original
		  // THit's are considered as non genuine), this 2nd TKine (as
		  // well as the 1st one) has to be consistent w/ the cluster's
		  // centroid (condition "|u-uC|>4*sigu" required for condition
		  // "iHOrigs[iRef]==iorig==0" to be fulfilled). Hence is nicely
		  // settled the case of 2 neighbouring GEM (e.g.) clusters
		  // which had parly overlapped each other, have been split
		  // apart by the GEM clusterisation method but still retain,
		  // small, contributions from each other.
		  ikin2 = it;
		break;
	      }
	    }
	  }
	}
      }
      if (ikine >= 0 ) {        // for debuging
	if( fabs(1./vecKine[ikine].Pinv()) < TOpt::dCut[18] ||
	    fabs(1./vecKine[ikine].Pinv()) > TOpt::dCut[19] ) continue;
      }
      jkine = ikine;  // Save "ikine" in "jkine"

      if (iorig==-2) iorig = -1;  // ********** FINALISE Iorig FLAG **********

    } // End of MC block


    // ******************** CREATE THit object ********************

    THit h;
    h.iPlane = iPlane;
    int ihit = vecHit.size(); h.iHit = ihit;       // Ref. to itself 
    // ********** SET MIRROR REFERENCE **********
    const TPlane &p = setup.vPlane(iPlane);
    if (mirror) {
      double up = mirror->getU()/10;
      if (prvCluster && prvCluster==mirror) {
	// This only works when all clusters are imported in one go (and not
	// track by track, via "ImportTracks" and then the rest)... except if
	// there are more than one (cluster,pair) per cell.
	vecHit.back().mirror = ihit; h.mirror = ihit-1;
      }
      else if (prvPlane==iPlane && up<=u) {
	// This corresponds to the case when mirror has already been imported in
	// a previous pass when "ImportClusters" was called by "ImportTracks".
	//  => Scan vHitRef in search for mirror THit
	// (It's slow but should never be brought into play if one cares to
	// use Traffic's secret trap ReMode[1])
	vector<int>::const_iterator ih; bool match;
	for (ih = p.vHitRef().begin(), match = false; ih!=p.vHitRef().end();
	     ih++) {
	  THit &hp = vecHit[*ih]; if (hp.ptrCl==mirror) {
	    match = true; h.mirror = *ih; hp.mirror = ihit;
	  }
	}
	if (!match)
	  CsErrLog::msg
	    (elError,__FILE__,__LINE__,
	     "TPlane %d THit %d,%.3f has mirror %.3f and yet no THit mirror",
	     iPlane,ihit,u,up);
      }
      // ***** COALESCE (drift) THit WITH ITS MIRROR WHEN CLOSE DISTANCE *****
      const float xSAS = 1300;
      float coalesce = det.X(0)<xSAS ?
	TOpt::dCut[26] /*LAS*/ : TOpt::dCut[27] /*SAS*/;
      if (p.IFlag && // Exclude OFF plane (for it's anyway useless and
	  // OFF probably means one wants to study the response of the detector,
	  // which is more faithfully maintained if no coalescence.
	  fabs(up-u)<coalesce*det.Resol) {
	sigu += fabs(up-u)/2; u = (up+u)/2;
	// The 1st encountered member of the coalesced pair is checked. If
	// genuine, gets a proba of 1, while its counterpart gets 0, and later
	// gets a TStatus of -4 (meaning not to be used in PR). If fake, proba's
	// are inverse, current member gets a TStatus of 4 and its counterpart
	// is later left unchanged.
	if (h.mirror==-1) {
	  if (!isMC || isMC && clus->isGenuine()) {
	    clus->setLRProb(1); mirror->setLRProb(0);
	  }
	  else {                       // If !genuine...
	    clus->setLRProb(0); mirror->setLRProb(1);
	    h.status = -4;             // ...=> Assign special value to status
	  }
	}
	else {                       // Counterpart already encountered 
	  if (vecHit[h.Mirror].Status!=-4)           // If conterpart genuine...
	    h.status = -4;             // ...=> Assign special value to status
	}
      }
    }
    //else h.mirror = -1; // Set by the constructor

    h.u = u; h.sigu = sigu;
    h.v = clus->getV()/10.; h.sigv = sqrt(m(2,2))/10.; 
    h.ptrCl = clus;

    // ********** SET Status **********
    if (status)  // Take into account GEMs spacer flag and SI timeless flag
      h.status = status;
    if (h.Status==0) {          // May already be set, cf. supra and coalescence
      //      if (clus->getLRProb()<.5)
      if (clus->getLRProb()<.45) { // ***** GIVE MALUS to LOW LR PROBAS *****
	//  I) LR proba = 1 for non drift detectors
	// II) Outer detectors are less crowded => looser cut on outer STs and
	//    outer parts of DCs, DWs and MDs.  DRs (or ST04) are special case
	//    crowded as they can be by photon conversions.
	//III) Low proba: "status=-2", i.e. used only in later iterations of PR
	// IV) Very low proba: "status=-3", i.e. not used in PR
	int typ = det.IType; float prob = clus->getLRProb(); bool inner;
	if      (typ==11 || typ==12) inner = true;  // Inner ST 
	else if (typ==15 || typ==16 || typ==18) {   // DC, DW, MD
	  float wire = (h.U-det.Uorig)/det.Pitch+.5;
	  inner = det.Nwires*.25<wire && wire<det.Nwires*.75;
	}
	else inner = true;
	if (prob<.35 || inner) h.status = -2;
	if (TOpt::ReMode[28] && (prob<.25 || inner && prob<.35)) h.status = -3;
      }
    }
    if (isMC) {// ******************** YET ANOTHER MC BLOCK ********************

      // If this cluster is mirror hit in drift detector:
      // - Remove association with TKine. But remember it by setting "ikine"
      //  = "-2-ikine", in order to document the debug info.
      // Note concerning coealesced clusters: In that case the dormant component
      // (i.e. that w/ "Status=-4") is chosen to be always the fake mirror, and
      // it receives here a <0 "TKine" attribute. If the cluster is de-coalesced
      // at some point, the decision to promote one or the other mirror is taken
      // (cf. TEv::TracksFit2 or TracksRefit) on view the U coordinate of the
      // associated CsClusters. Therefore, the THit eventually retained will get
      // the proper fake or genuine qualification. If the cluster remains
      // coalesced, the THit retained will be declared genuine, as should be.
      if (csDet->hasDrift() && !clus->isGenuine()) {
	if (ikine!=-1) ikine = -2-ikine;
	// - Do not do the same for ikin2: to simplify
      }
      h.iKine  = ikine; h.iKin2  = ikin2; h.iOrig  = iorig;
    }

    // ********** CsDigit info -> Current THit's set<TDigits> **********
    for (id = lDig.begin(); id!=lDig.end(); id++) {
      // ***** LOOP OVER CsDigits IN THE CsCluster
      int iwire    = (*id)->getAddress();
      int ndata    = (*id)->getDataSize(); // Associated data size
      double* data = (*id)->getData();     // Pointer to associated data array
      vector<float> info;
      if(ndata <= 0 && data == NULL){
	/*if (!isMC){
	  cout<<"TEv::ImportClusters ==> No associated information for CsDigit"<<endl
	      <<"Det. plane # "<<iPlane<< " THit # "<<h.IHit<<" Digit from wire # "<<iwire<<endl;
	      }*/
      } else {
	for(int j = 0; j < ndata; j++) info.push_back(float(data[j]));
      }
      TDigit td(iwire,info);  // create TDigit object
      h.setDigits.insert(td); // store it into THit object "h"
    }

    if (det.TResol<0)              // No time measurement. Cf. "TSetupInit".
      h.sigt = -1;
    else {                         // Detector  does measure time
      if (clus->getTime(h.time)) {   // CsCluster does have time measurement
	bool hasSigT = clus->getTimeError(h.sigt); assert(hasSigT);
	//if (h.sigt<0) h.sigt = det.TResol; // This is enforced by infra...
	if (h.sigt<det.TResol)
	  h.sigt =  det.TResol;        // Lower limit = det. time resolution
#define ImpClusters_SIs_SUBDUED
#ifdef ImpClusters_SIs_SUBDUED
	if (det.IType==21) // Sis'in 2007 (at least):
	  // i) Uncertainties are overoptimistic. ii) Calibrations are careless.
	  // => Reduce their weight
	  h.sigt *= siTSigScale;
#endif
      }
      else {                         // CsCluster does'nt have time measurement
	h.sigt = -1;
	CsErrLog::msg(elError,__FILE__,__LINE__,
		      "CsCluster from \"%s\" should have, but has no time",
		      det.Name.c_str());
      }
    }

    if (det.IType!=44) // Exclude cylindrical hodos, for their description...
      // ...is not yet understood by TraFDic, and their processing by
      // "THit::Rotate" triggers the following kind of error message:
      // "THit::Rotate ==> non-realistic Sig V = 2.4e+17  from HH02R"
      // Means that those detectors have to be turned off in tracking.
      h.Rotate();

    vecHit.push_back(h);   // *************** STORE the THit ***************

    // ********** REFERENCE TPlane -> THit in coord U sorted order **********
    TPlane &pl = const_cast<TPlane&>(setup.vPlane(iPlane));
    pl.addHitRef(h.iHit);

    if (isMC) {// ******************** YET ANOTHER MC BLOCK ********************

      // *************** TKine->THit ASSOCIATION **********

      if (jkine>=0 && ikine!=-1) {
	vecKine[jkine].vecHitRef.push_back(h.iHit);   // Reference TKine->THit
	if (ikin2>=0) // 2nd track associated to hit
	  vecKine[ikin2].vecHitRef.push_back(h.iHit);// ...2-track hit
      }
      //                                       ***** THitMC->THit association
      if (sMChit.empty())
	CsErrLog::mes(elError,"CsCluster has no reference to CsMCHit");
      for (iMChit = sMChit.begin(); iMChit!=sMChit.end(); iMChit++){
	// LOOP over "MOTHER" MCHits of THIS CLUSTER
	map<CsMCHit*, int, less<CsMCHit*> >::const_iterator im;
	im=mapCsMCHit2HitMC.find(*iMChit);
	if(im == mapCsMCHit2HitMC.end()){
	  CsErrLog::mes(elError,"CsMCHit not imported by TEv::GetMCInfo");
	  continue;
	}
	THitMC& hmc = vecHitMC[(*im).second];
	hmc.setHitRef.insert(h.iHit);
      } 

      if      (ikine==-1) nClsImported_m1++;
    } // End of  yet another MC block

    mapCsCl2Hit[h.ptrCl]=h.iHit;    // Fill CsCluster* -> THit index map

    prvCluster = h.ptrCl; prvPlane = iPlane; // Remember cluster/plane in next loop
    
  }   // ******************** END OF LOOP OVER CsClusters ********************

  //    ***** FILL ARRAYS of ALTERNATIVE HIT REFERENCES in PIXEL PLANES *****
  // (Note: The pixel planes herein considered are the P-pixel ones, as opposed
  // to M-pixels of MPs. As of 2012/07, they comprise only pixel pieces of GPs.)
  for (int zone = 0; zone<(int)setup.vIplFirst().size(); zone++) {
    double sqr2 = sqrt(2.), cs[4] = {1,0,sqr2,sqr2}, ss[4] = {0,1,sqr2,-sqr2};
    if (1<<zone&TOpt::iCut[29]) { // Zones where special pixel projection
      double a = TOpt::dCut[82]/180*M_PI; cs[2] = cos(a); ss[2] = sin(a);
      a        = TOpt::dCut[83]/180*M_PI; cs[3] = cos(a); ss[3] = sin(a);
    }
    for (int ipl = setup.vIplFirst()[zone]; ipl<=setup.vIplLast()[zone];
	 ipl++) {
      TPlane &p = const_cast<TPlane&>(setup.vPlane(ipl));
      if (!(p.IFlag&0x30)) continue;
      for (int coord = 0; coord<4; coord++)     // Fill XYUV-wise references
	p.fillPixRefs(coord,cs[coord],ss[coord]);
    }
  }

  if (hist && beamTelescope<=0) h[1]->Fill(nClsImported_m1+0.5);
}
