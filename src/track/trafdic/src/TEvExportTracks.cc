// $Id: TEvExportTracks.cc 13148 2011-12-28 16:55:25Z kbicker $

/*!
  Export tracks reconstructed in TrafDic to coral
*/

#include "CsErrLog.h"
#include "CsEvent.h"
#include "CsTrack.h"
#include "CsHelix.h"
#include "CsCluster.h"
#include "CsZone.h"
#include "CsMCTrack.h"
#include "TOpt.h"
#include "TSetup.h"
#include "TEv.h"

using namespace std;

/*
  Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/TEvExportTracks.cc":
  i) Set CsTrack's shower attributes from corresponding TTrack's
*/

void TEv::ExportTracks(list<CsTrack*>& lCsTrk)
{
  if (TOpt::ReMode[30]>0) return; // Export nothing

  const TSetup& setup = TSetup::Ref();
  CsEvent* csev = ptrEvt();
  assert(csev != NULL);
  map<CsDetector*,int> det2bit = csev->getDetectorMap();
  if(det2bit.empty()){
    cout<<"TEv::ExportTracks ==> CsEvent::getDetectorMap() returns empty map !"<<endl;
    assert(false);
  }
  map<CsDetector*,int>::const_iterator id2b;

  list<TTrack>::iterator it;
  for (it = listTrack.begin(); it!=listTrack.end(); it++) {

    // ************************* LOOP OVER TTracks *************************
    vector<CsHelix> helices; list<CsCluster*> clusters; list<CsZone*> zones;
    TTrack& t = *it;

    //                 ********** PREPARE HELICES **********
    if (!t.Hfirst.empty())
      helices.push_back(t.Hfirst.ExportHelix());     // First helix
    list<THlx>::iterator ihsm;
    for (ihsm = t.lHsm.begin(); ihsm!=t.lHsm.end(); ihsm++)
      helices.push_back((*ihsm).ExportHelix());      // Smoothed helices, if any
    if (!t.Hlast .empty())
      helices.push_back(t.Hlast. ExportHelix());     // Last helix

    //                ********** PREPARE CLUSTERS **********
    list<int>::iterator ih;
    for (ih = t.lHitPat.begin(); ih!=t.lHitPat.end(); ih++) {
      //                                      ***** LOOP OVER TTrack's HITS LIST
      if (*ih<0) continue;
      CsCluster* pcl = vHit(*ih).PtrClus();
      if (pcl==NULL) {
	cout<<"TEv::ExportTracks ==> Track ID "<<t.Id<<" Hit # "<<(*ih)<<" has NULL pointer to CsCluster"<<endl;
	assert(false);
      }
      clusters.push_back(pcl);
    }

    // ******************** INSTANTIATE CsTrack ******************************
    // Prepare zones
    for (int igr = 0; igr<int(setup.vIplFirst().size()); igr++) {
      if (t.InGroup(igr)!=0) zones.push_back(setup.Group2Zone(igr));
    }
    CsTrack *pCsTrk = new CsTrack(helices,clusters,zones,t.Chi2tot);

    //               ********** EXTRA INFO **********
    pCsTrk->setNDFs(t.NDFs);
    if (t.SigmaTime>0) pCsTrk->setMeanTime(t.MeanTime, t.SigmaTime);
    pCsTrk->setXX0(t.radLenFr);
    pCsTrk->setShower(t.NShowerHits,t.HasShower);

    //          ********** FILL HIT BIT PATTERNS **********
    unsigned int expect[CSTRACK_MAPSIZE], hited[CSTRACK_MAPSIZE];
    for (int i = 0; i<CSTRACK_MAPSIZE; i++) { expect[i]=hited[i] = 0; }
    list<int>::const_iterator ipl = t.lPlnRef.begin();
    list<int>::const_iterator ihp = t.lHitPat.begin(); 
    while (ipl!=t.lPlnRef.end()) { // Loop over track's plane references
      const TDetect& d = setup.iPlane2Detect(*ipl);
      id2b = det2bit.find(d.PtrDet());
      if (id2b==det2bit.end()) {
	cout<<"TEv::ExportTracks ==> CsDetect* of "<<d.Name<<" not in  CsEvent's map<CsDetectors*,int>"<<endl;
	assert(false);
      }
      int word = int(((*id2b).second)/32);
      int bit =      ((*id2b).second)%32 ;
      if (word>=CSTRACK_MAPSIZE) {
	cout<<"TEv::ExportTracks ==> int array["<<CSTRACK_MAPSIZE<<
	  "] is not enough to store hit pattern's bit # "
	    << (*id2b).second <<"  for TPlane # "<<(*ipl)<<endl; 
	assert(false);
      }
      if((*ihp) >= -1 ) expect[word] = expect[word] | 1<<bit;
      if((*ihp) >=  0 ) hited [word] = hited [word] | 1<<bit;
      ipl++; ihp++;
    } // End of loop over track's plane references

    // store hit pattern bitmasks

    //if(t.Pinv() != 0) {cout<<1./t.Pinv()<<" Gev  ";for(int i=0; i<10; i++) {cout<<expect[i]<<" ";} cout<<endl;}
    pCsTrk->setExpectedDetsBitmap(expect);
    //for(int i=0; i<10; i++) cout<<hited[i] <<" "; cout<<endl;
    pCsTrk->setFiredDetsBitmap   (hited );
    //cout<<endl;

    // store the pointer
    int id = pCsTrk->getId();
    map<int, int, less<int> >::const_iterator im;
    im = mapCsTrId2TrId.find(id);
    if(im != mapCsTrId2TrId.end()){ // already there. Something is wrong.
      cout<<"TEv::ExportTracks ==> CsTrack ID = "<<(*im).first
	  <<"  already had been maped to TTrack ID = "<<(*im).second<<endl;
      assert(false);
    }
    mapCsTrId2TrId [id] = t.Id; // map CsTrack ID to TTrack ID
    t.ptrTrk = pCsTrk;          // TTrack -> CsTrack back reference

    if (t.IKine>=0) {          // ***** EXPORT MC <-> RECO ASSOCIATION *****
      CsMCTrack *pCsMCTrk =
	const_cast<CsMCTrack*>(vecKine[t.IKine].PtrTrk());
      pCsMCTrk->addAssociatedTrackID(t.PtrTrk()->getId());
      pCsTrk->setAssociatedMCTrack(pCsMCTrk);
    }

    lCsTrk.push_back(pCsTrk);   // store new track pointer
  } // End of loop over TTracks

  //            ***** ASSOCIATED RECONSTRUCTED TRACK *****
  //  TTrack::Association (so far) concerns tracks (presumed) non-interacting
  // and continued beyond SM1, i.e. 0x1X (X=3,7,f), and their 0xX and 0x10
  // pieces: 0x10.Associate = 0xX.
  // (In the vertex package (where the 0x10 piece may have been assigned a
  // momentum), it will be checked whether it's not a Q^2=0 event or a non
  // interacting beam. If the latter, it will be discarded/dwongraded!?.)
  list<CsTrack*>::iterator icst, jcst;
  for (icst = lCsTrk.begin(); icst!=lCsTrk.end(); icst++) {
    CsTrack *cst = *icst; int csId = cst->getId();
    map<int,int,less<int> >::const_iterator im = mapCsTrId2TrId.find(csId);
    if (im==mapCsTrId2TrId.end()) // Some CsTracks may not orginate from TraFFiC
      continue;
    int id = im->second; const TTrack *t;
    for (it = listTrack.begin(), t = 0; it!=listTrack.end(); it++) {    
      TTrack &ti = *it; if ((int)ti.Id==id) { t = &ti; break; }
    }
    if (!t) CsErrLog::msg(elFatal,__FILE__,__LINE__,
  "CsTrack #%d's TTrack counterpart(=#%d) does not exist",csId,id);
    if (t->Type!=0x10) // Association (so far) only for beams...
      continue;       // ...nevertheless, to be on the safe side.
    int ia = t->Associate, csIa; if (ia<0) continue;
    for (im = mapCsTrId2TrId.begin(), csIa = -1; im!=mapCsTrId2TrId.end(); im++)
      if (im->second==ia) { csIa = im->first; break; }
    if (csIa<0) // No CsTrack w/ ID corresponding to associate. While the
      // associate's integrity has been checked in "TracksFit2" = Fatal.
      CsErrLog::msg(elFatal,__FILE__,__LINE__,
  "TTrack #%d's associate(=#%d) has no CsTrack counterpart",id,ia);
    cst->setAssociatedTrack(csIa);
  }
}
