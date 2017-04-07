// $Id: TEvImportTracks.cc 13130 2011-12-20 14:36:05Z kbicker $

#include "TEv.h"
#include "CsTrack.h"

#ifdef COMPASS_USE_OSPACE_STD
# include <ospace/std/algorithm>
#else
# include <algorithm>  
#endif

using namespace std;

void TEv::ImportTracks(list<CsTrack*>& lCsTrk, string mode)
{
  if( mode != "cp" && mode != "mv" ){
    cout<<"TEv::ImportTracks: Wrong mode : "<<mode<<endl;
    assert(false);
  }

  const TSetup& setup = TSetup::Ref();

  list<CsTrack*>::iterator it = lCsTrk.begin();
  while(it != lCsTrk.end()){ // loop over CsTracks
    CsTrack* pTrk = (*it); // pointer to this CsTrack
    // Check if this track  already was imported
    map<int, int, less<int> >::const_iterator im;
    im = mapCsTrId2TrId.find(pTrk->getId());
    if(im != mapCsTrId2TrId.end()) { // already exists
      if(mode == "mv"){
	delete pTrk;        // delete corresponding CsTrack
	lCsTrk.erase(it++); // delete entry to the list
      } else {
	it++;
      }
      continue; 
    }

    // import clusters, assosiated with this track
    list<CsCluster*> lCsCl = pTrk->getClusters();
    ImportClusters(lCsCl,0);

    TTrack t; // create new track
    t.Type = 0;

    // fill bitmask of det. groups (Type) of the track 
    list<CsZone*> lZone;
    lZone = pTrk->getZones();
    // as N groups could be >= N zones
    for(int igroup = 0; igroup < int(setup.vIplFirst().size()); igroup++){ //loop over all det. groups
      CsZone* pZone = setup.Group2Zone(igroup);
      if( find(lZone.begin(), lZone.end(), pZone) != lZone.end() ) t.Type |= 1<<igroup;
    } 
    // store hits
    list<CsCluster*>::iterator ic;
    map<CsCluster*, int, less<CsCluster*> >::const_iterator im1;
    for(ic = lCsCl.begin(); ic != lCsCl.end(); ic++){ // loop over track's clusters
      im1 = mapCsCl2Hit.find(*ic);
      if(im1 == mapCsCl2Hit.end()) {
	cout<<"TEv::ImportTracks ==> Something is wrong: the cluster not yet imported to TEv object"<<endl;
	assert(false);
      }
      int ih = (*im1).second;
      t.AddHit(vHit(ih));
    }
    // store track parameters
    vector<CsHelix> vH = pTrk->getHelices();
    if(vH.size() == 1) t.Hfirst.ImportHelix(vH[0]);
    if(vH.size() >  1) {
      t.Hfirst.ImportHelix(vH.front());
      t.Hlast. ImportHelix(vH.back ());
    }

    //Chi2, time (if any), etc.
    t.Chi2tot  = pTrk->getChi2();
    double mt;
    if(pTrk->hasMeanTime()) t.MeanTime = pTrk->getMeanTime();
    else                    t.MeanTime = 8888.0;

    t.radLenFr = pTrk->getXX0();

    if(mode == "mv"){
      delete pTrk;        // delete corresponding CsTrack
      lCsTrk.erase(it++); // delete entry to the list
    } else {
      t.ptrTrk = pTrk;   // TTrack -> CsTrack reference (valid only in "copy" mode) 
      it++;
    }

    // do TTrack <--> TKine assosiation (if MC event)
    t.FindKine();
    // calculate track mean time, error, dispersion 
    t.UseHitTime();

    listTrack.push_back(t); // store new TTrack


  } // end of loop over CsTracks

}




