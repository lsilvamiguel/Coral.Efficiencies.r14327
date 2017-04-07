/*
------------------------------------------------------------------------

 Implementation fo classes to store/cluster events from COMPASS Mumega
 detectors

 Author: Bernhard Ketzer    02/06/2009     v0

------------------------------------------------------------------------
*/

#include <cassert>
#include <cstdlib>   // for abs()

#include "CsMumegaPlane.h"

#include <TApplication.h>

#include <fstream> // for debug 


ClassImp(CsMumegaPlane)


CsMumegaPlane::CsMumegaPlane() {

  // Create a CsMumegaPlane object
  fChannels.clear();
  fHits.clear();
  fClusters.clear();
  fIsClusterized = kFALSE;
  fIsHitSorted = kFALSE;
  fIsHitAmpSorted = kFALSE;
  fIsClusterSorted = kFALSE;
  fKey = 0;
  fDisplayCanvas = 0;
  fDisplayHist = 0;
  fName = "";

}


CsMumegaPlane::CsMumegaPlane(std::string _name) {

  // Create a CsMumegaPlane object
  fChannels.clear();
  fHits.clear();
  fClusters.clear();
  fIsClusterized = kFALSE ;
  fIsHitSorted = kFALSE;
  fIsHitAmpSorted = kFALSE;
  fIsClusterSorted = kFALSE;
  fKey = 0;
  fDisplayCanvas = 0;
  fDisplayHist = 0;
  fName = _name;  

}


CsMumegaPlane::~CsMumegaPlane() {
  Clear();
  ClearChannels();
  fTimeCals.Clear();
  if (fDisplayCanvas != NULL) fDisplayCanvas->Close();
  if (fDisplayHist != NULL) fDisplayHist->Delete();

}


void CsMumegaPlane::AddChan(Int_t _channr, Int_t _hem, Int_t _flag, Float_t _ped, Float_t _sigma, Int_t _apvchip, Int_t _apvchan, int _stripconn) {
  // Channel Id and channel calibration objects
  CsMumegaChanId id = CsMumegaChanId(_channr, _hem);
  id.SetAPVInfo(_apvchip, _apvchan);
  id.SetSNConnNb(_stripconn);
//   id.SetSConnNb(_channr, _hem);
  CsMumegaChanCal cal = CsMumegaChanCal(_flag, _ped, _sigma);

  // Return if channel with same Id already exists
  if (GetChan(id) != NULL) {
    std::cout<<"CsMumegaPlane::AddChan() : Channel with same Id already exists ! Det.chan = "<<_channr<<", hemisphere = "<<_hem<<std::endl;
    return;
  }

  // Create new channel
  CsMumegaChan* chan = new CsMumegaChan(id, cal);

  // Loop over all existing channels to update list of active neighbours
  std::map<CsMumegaChanId, CsMumegaChan*>::iterator itchan;
  for (itchan = fChannels.begin(); itchan != fChannels.end(); itchan++) {

    bool neighbour = false;

    // Direct neighbours
    if (abs(itchan->first.GetDetectorChannel() - _channr) <= 1) neighbour = true; //original
    //   if (abs(itchan->first.GetDetectorChannel() - _channr) == 1 && (itchan->first.GetHemisphere() == _hem ) ) neighbour = true;  //u adjacent and same hem 
    //   if (abs(itchan->first.GetDetectorChannel() - _channr) == 0 ) neighbour = true; //same wire on different hem 

       // for PixelMumega -RECTANGULAR PIXELS- (names start with "MP") do not connect the two hemispheres of strips [384, 511] (128 strips), they are divided by the pixels
    if ( GetName().substr(0,2) == "MP" ) { // strip plane of PixelMumega detector
      if ( ( (_channr >= 384 && _channr <= 511) // either to be added or the existing channel have to be in the range
	     || (itchan->first.GetDetectorChannel() >= 384 && itchan->first.GetDetectorChannel() <= 511) )
	   && itchan->first.GetHemisphere() != _hem ) { // the hemispheres of added and existing differ
	// then strips have to be in different hemispheres, otherwise mapping is broken
	neighbour = false;
      }
      }


    // Other criteria for neighbours may be added here, e.g. if separated by flagged channels

    // Add to list if active neighbour
    if (neighbour) {
      if (itchan->second->GetCal()->GetFlag() != 0) chan->AddActiveNeighbour(itchan->second);
      if (_flag != 0) itchan->second->AddActiveNeighbour(chan);
    }
    
  }
  // Add new channel to map
  AddChan(chan);
}


void CsMumegaPlane::AddHit(Int_t _detchan, Int_t _hem, Float_t _amp0, Float_t _amp1, Float_t _amp2) {
  std::vector<float> amp;
  amp.push_back(_amp0);
  amp.push_back(_amp1);
  amp.push_back(_amp2);
  AddHit(_detchan, _hem, amp);
}


void CsMumegaPlane::AddHit(Int_t _detchan, Int_t _hem, std::vector<Float_t> _amp) {
  // Get pointer to corresponding CsMumegaChan
  CsMumegaChanId id = CsMumegaChanId(_detchan, _hem);
  const CsMumegaChan* chan = GetChan(id);

  // Channel not found
  if (chan == NULL) {
    std::cout<<"CsMumegaPlane::AddHit() : chan "<<_detchan<<", hem "<<_hem<<" does not exist in "<<fName<<" ! Ignoring it"<<std::endl;
    return;
  }

  // Get channel calibrations
  const CsMumegaChanCal *cal = chan->GetCal();

  // Check flag
  int flag = cal->GetFlag();
  if (flag == 0) return;

  // Check amplitude
  float sigma = cal->GetPedSigma();
  float thr = fPar.GetThrHit();
  int sample = fPar.GetSample();
  if ( !(_amp[sample] > thr*sigma) ) return; 
  
  // Create new hit
  CsMumegaHit *hit = new CsMumegaHit(chan, _amp);
  /*  // Calculate hit time
  if (!fTimeCals.IsValid()) {
    std::cout<<"csMumegaPlane::AddHit() : time calibration not valid for plane "<<fName<<std::endl;
  } else {
    double time, timeerr;
    if (fTimeCals.CalcTime(_amp, sigma, time, timeerr))
      hit->SetTime(time, timeerr);
      }*/
  fHits.push_back(hit);
}


struct CsMumegaPlane::fCompareHits :
  public std::binary_function<CsMumegaHit*, CsMumegaHit*, bool> {
  bool operator() (CsMumegaHit *h1, CsMumegaHit *h2) {
    if (!h1) return true;
    if (!h2) return false;
    return (h1->GetChan()->GetId() < h2->GetChan()->GetId());
  }
};


void CsMumegaPlane::SortHits() {
  // Sort list of hits by channel Id
  fHits.sort( fCompareHits() );

  // Set flag indicating that event has been sorted
  fIsHitSorted = kTRUE;
  fIsHitAmpSorted = kFALSE;
}


struct CsMumegaPlane::fCompareHitAmps0 :
  public std::binary_function<CsMumegaHit*, CsMumegaHit*, bool> {
  bool operator() (CsMumegaHit *h1, CsMumegaHit *h2) {
    if (!h1) return false;
    if (!h2) return true;
    return ((h1->GetAmp())[0] > (h2->GetAmp())[0]);
  }
};


struct CsMumegaPlane::fCompareHitAmps1 :
  public std::binary_function<CsMumegaHit*, CsMumegaHit*, bool> {
  bool operator() (CsMumegaHit *h1, CsMumegaHit *h2) {
    if (!h1) return  false;
    if (!h2) return true;
    return ((h1->GetAmp())[1] > (h2->GetAmp())[1]);
  }
};


struct CsMumegaPlane::fCompareHitAmps2 :
  public std::binary_function<CsMumegaHit*, CsMumegaHit*, bool> {
  bool operator() (CsMumegaHit *h1, CsMumegaHit *h2) {
    if (!h1) return false;
    if (!h2) return true;
    return ((h1->GetAmp())[2] > (h2->GetAmp())[2]);
  }
};


void CsMumegaPlane::SortHitAmps() {
  // Sort list of hits by channel Id
  int sample = fPar.GetSample();
  if (sample == 0) fHits.sort( fCompareHitAmps0() );
  else if (sample == 1) fHits.sort( fCompareHitAmps1() );
  else fHits.sort( fCompareHitAmps2() );

  // Set flag indicating that event has been sorted
  fIsHitAmpSorted = kTRUE;
  fIsHitSorted = kFALSE;
}


struct CsMumegaPlane::fCompareClusters :
  public std::binary_function<CsMumegaCluster*, CsMumegaCluster*, bool> {
  bool operator() (CsMumegaCluster *c1, CsMumegaCluster *c2) {
    if (!c1) return true;
    if (!c2) return false;
    return (c1->GetPosition() < c2->GetPosition());
  }
};


void CsMumegaPlane::SortClusters() {
  // Sort list of clusters by cluster center
  fClusters.sort( fCompareClusters() );
  // Set flag indicating that event has been sorted
  fIsClusterSorted = kTRUE;
}


void CsMumegaPlane::Clear(Option_t *option) {
  // Delete all hits in list
  ClearHits();

  // Delete all clusters in list
  ClearClusters();
}


int CsMumegaPlane::ClearChannels(Option_t *option) {
  // return if channels are cleared
  if (!fChannels.size()) return 0;

  // Delete all objects in map
  std::map<CsMumegaChanId, CsMumegaChan*>::iterator it;
  for (it = fChannels.begin(); it != fChannels.end(); it = fChannels.begin()) {
    CsMumegaChan* ptr = it->second;
    fChannels.erase(it);
    delete ptr;
    if (!fChannels.size()) return 0;
  }

  // should never be reached
  fChannels.clear();

  return 1;
}


void CsMumegaPlane::ClearHits(Option_t *option) {
  // Delete all objects in list
  std::list<CsMumegaHit*>::iterator it;
  for (it = fHits.begin(); it != fHits.end(); it++) {
    delete *it;
  }
  fHits.clear();
  fIsHitSorted = kFALSE;
  fIsHitAmpSorted = kFALSE;
}


void CsMumegaPlane::ClearClusters(Option_t *option) {
  // Delete all objects in list
  std::list<CsMumegaCluster*>::iterator it;
  for (it = fClusters.begin(); it != fClusters.end(); it++) {
    delete *it;
  }
  fClusters.clear();
  fIsClusterized = kFALSE;
  fIsClusterSorted = kFALSE;
}


void CsMumegaPlane::Reset(Option_t *option) {
  // Static function to reset all static objects for this event
}

void CsMumegaPlane::SetHeader(Int_t i, Int_t run, Int_t date)
{
   fHeader.Set(i, run, date);
}

void CsMumegaPlane::Clusterize() {
  // Clear list of clusters
  ClearClusters();

  // Clustering method
  int method = fPar.GetClusMeth();
  switch(method) {

    // Simple clustering
  case 0:
    {
      SimpleClustering();
      break;
    }
    // Full clustering
  case 1:
    {
      FullClustering();
      break;
    }
    // Full clustering, always fit
  case 2:
    {
      FullClustering();
      break;
    }

  default:
    FullClustering();
  }

  // Calculate cluster properties
  float thrhit = fPar.GetThrHit(); // threshold on single hit amplitude
  int sample = fPar.GetSample(); // sample used for clustering
  int share = fPar.GetShareHits(); // hit sharing ?
  std::vector<Float_t> time_cal_offset = fPar.GetTimeCalOffset();
  std::vector<Float_t> time_cal_slope = fPar.GetTimeCalSlope();
  std::vector<Float_t> time_cal_time_offset = fPar.GetTimeCalTimeOffset();
  std::list<CsMumegaCluster*>::iterator itclus;
  itclus = fClusters.begin();
  while (itclus != fClusters.end()) {
    if ((*itclus)->IsOldClus()) {delete(*itclus); itclus = fClusters.erase(itclus); continue;}
    // Cluster center and amplitude
    if (method == 2 || (*itclus)->ContainsFlaggedChannels()) {
      // (*itclus)->Fit(thrhit);
      std::cout<<"CsMumegaPlane::Clusterize() : Fitting of clusters with flagged channels not yet implemented !"<<std::endl;
      std::cout<<"                       Using CoG method instead..."<<std::endl;
      (*itclus)->CalcAmps(sample, thrhit, share);
      (*itclus)->CalcCoG(sample, thrhit, share);
    } else {
      (*itclus)->CalcAmps(sample, thrhit, share);
      (*itclus)->CalcCoG(sample, thrhit, share);
    }

    //    std::cout<<GetName()<<" offset "<<time_cal_offset[0]<<" slope "<<time_cal_slope[0]<<" time offset "<<time_cal_time_offset[0]<<std::endl;

	// Cluster time
    if (GetPar()->IsTimeCal()) (*itclus)->CalcTime(time_cal_offset, time_cal_slope, time_cal_time_offset);
    else (*itclus)->SetTime(0., 1.e9);
    
    ++itclus;
    
  }

    // Select only clusters above threshold
    SelectClusters();

    // Set flag indicating that event has been clusterized
    fIsClusterized = kTRUE;

    return;
}



void CsMumegaPlane::FullClustering() {

  // Sort hits by amplitude
  if (!fIsHitAmpSorted) SortHitAmps();

  // Mark cross talk hits
  if (GetPar()->IsTxtCal()) FindXTalk();

  // Iterators
  std::list<CsMumegaHit*>::iterator ithit;
  std::list<CsMumegaCluster*>::iterator itclus;

  // Group adjacent hits to cluster candidates
  for (ithit = fHits.begin(); ithit != fHits.end(); ithit++) {
    //    std::cout<<"hit "<<(*ithit)->getChan()->GetId()->GetDetectorChannel()<<std::endl;

    bool assigned = false;

    // Loop over existing clusters
    for (itclus = fClusters.begin(); itclus != fClusters.end(); itclus++) {
      // Check if hit belongs to existing cluster, update its properties
      if ((*itclus)->AddHit(*ithit)) {
	assigned = true;
	(*ithit)->AddCluster(*itclus);
	(*ithit)->IncNrClusters();
	//if ((*ithit)->GetNrClusters()>1) (*itclus)->IncrNSharedHits();
      }
    }

    // Create new cluster if hit was not assigned to any existing cluster
    if (!assigned) {
      CsMumegaCluster* cl = new CsMumegaCluster();
      cl->AddHit(*ithit);
      fClusters.push_back(cl);
      (*ithit)->AddCluster(cl);
      (*ithit)->IncNrClusters();
      // if ((*ithit)->GetNrClusters()>1) cl->IncrNSharedHits();
    }
  }


  for (ithit = fHits.begin(); ithit != fHits.end(); ithit++) {      // Set numbers of shared hits in each cluster (done here because not needed elsewhere yet)
    if ((*ithit)->GetNrClusters()>1){
      std::list<CsMumegaCluster*> cllist = (*ithit)->GetClusters();
      std::list<CsMumegaCluster*>::iterator ithitcl;
      for (ithitcl = cllist.begin(); ithitcl != cllist.end(); ithitcl++){
	(*ithitcl)->IncrNSharedHits();
      }
    }
  }

  // Cluster manipulations
  if (fPar.GetMergingAmpratioThreshold() < 1) MergeClusters(); // merge small clusters
  // SplitClusters(); // split large clusters

  return;
}


void CsMumegaPlane::SimpleClustering() {
  // Iterators
  std::list<CsMumegaHit*>::iterator ithit;

  // Loop over hits in plane
  for (ithit = fHits.begin(); ithit != fHits.end(); ithit++) {
    CsMumegaCluster* cl = new CsMumegaCluster();
    cl->AddHit(*ithit);
    fClusters.push_back(cl);
    (*ithit)->IncNrClusters();
     }

  return;
}

//-----------------------------------------------------------------------------
// Function which returns true if point xp,yp lies inside the
// polygon defined by the points in vectors x and y, false otherwise
// NOTE that the polygon must be a closed polygon (1st and last point
// must be identical) (taken from ROOT)
// Used in CsMumegaPlane::SelectClusters() for cuts on amplitude ratios
//-----------------------------------------------------------------------------
bool CsMumegaPlane::IsInside(float xp, float yp, std::vector<float> x, std::vector<float> y) {
	double xint;
	int i;
	int inter = 0;
	int np = x.size();
	for (i=0;i<np-1;i++) {
		if (y[i] == y[i+1]) continue;
		if (yp < y[i] && yp < y[i+1]) continue;
		if (y[i] < yp && y[i+1] < yp) continue;
		xint = x[i] + (yp-y[i])*(x[i+1]-x[i])/(y[i+1]-y[i]);
		if (xp < xint) inter++;
	}
	if (inter%2) return true;
	return false;
}

void CsMumegaPlane::SelectClusters() {

  float thrclus = fPar.GetThrClus(); // threshold on cluster amplitude
  int sample = fPar.GetSample(); // sample used for clustering
  double tcsphase = fHeader.GetTCSPhase();  //TCS Phase

  // Loop over clusters
  std::list<CsMumegaCluster*>::iterator itclus = fClusters.begin();
  while ( itclus!=fClusters.end() ) {
    bool remove = false; // remove current cluster

    double a0 = (*itclus)->GetAmp()[0];
    double a1 = (*itclus)->GetAmp()[1];
    double a2 = (*itclus)->GetAmp()[2];
    double amp = (*itclus)->GetAmp()[sample];

    // Remove cluster if below threshold
    double noise = (*itclus)->GetNoise();
       if ( !(amp > thrclus*noise) ) 
      {remove = true;
	}
      
   
    // Cut on amplitude ratios
    if (fPar.DoAmpCuts()){

      if ( !(IsInside(a1/a2, a0/a2, fPar.GetCoorda1a2(), fPar.GetCoorda0a2())) ) {remove = true;
}

      if ( !(IsInside(tcsphase, a1/a2, fPar.GetCoordtcs(), fPar.GetCoorda1a2tcs())) ) {remove = true;
}
    }

    if (fPar.IsTimeCal()){ //time cuts possible only if there is a valid time calibration
      
      if ( !((*itclus)->GetTime() - tcsphase > fPar.GetTimeCuts()[0] ) ) {remove = true;
}
      if ( !((*itclus)->GetTime() - tcsphase < fPar.GetTimeCuts()[1] ) ) {remove = true;
}
    }

    if (remove) {
      delete (*itclus);
      itclus = fClusters.erase(itclus);
    } else
      itclus++;
  }
  return;
}



void CsMumegaPlane::MergeClusters() { 

  //  std::list<CsMumegaHit*>::iterator ithit;    

  Float_t ampratiothr = fPar.GetMergingAmpratioThreshold(); 
  int sample = fPar.GetSample();

  std::list<CsMumegaCluster*>::iterator iclus;    // iterator for list of clusters of the event
  std::list<CsMumegaHit*>::iterator ihit;         // iterator for list of hits of the clusters 
  std::list<CsMumegaCluster*>::iterator iclhit;   // iterator for list clusters containing shared hits
  std::list<CsMumegaHit*>::iterator ihitclhit;    // iterator for list of hits belonging to clusters with shared hits
  std::list<CsMumegaHit*>::iterator idglist;      // iterator for list of hits of the merged cluster
  
  for (iclus = fClusters.begin(); iclus != fClusters.end(); iclus++) {    // loop on clusters of the event
    if ((*iclus)->GetNSharedHits() < 1) continue;                    // cluster must contain shared hits
    if ((*iclus)->IsOldClus() ||(*iclus)->IsNewClus()) continue;     // cluster must not be an "old" one (already merged into another cluster) nor a "new" one (created from merged clusters) (ie : clusters can only be merged once)
  
    std::list<CsMumegaHit*> dgclus = (*iclus)->GetHits();     // get hits of clusters with shared hits
    CsMumegaHit *dgmax = dgclus.front();   // hit with maximum amplitude
   
    for ( ihit = dgclus.begin(); ihit != dgclus.end(); ihit++){
   
      if ((*ihit)->GetNrClusters() > 1) {    //  hit must be shared between several clusters
	double ampratioshared = (*ihit)->GetAmp()[sample] / dgmax->GetAmp()[sample];    // compute amplitude ratio between the shared hits and the hit with maximum amplitude
	
	if (ampratioshared > ampratiothr){    //  amplitude ratio must be over threshold
	
	  std::list<CsMumegaCluster*> cldgcl = (*ihit)->GetClusters();     // get clusters the shared hits belongs to
	  std::list<CsMumegaHit*> dglist;                                  // list to store the hits of the new cluster
	  Bool_t flags = false;                                            // bool to store the "flags" value of the new cluster 
	  Int_t NSharedHits = 0;                                           // int to store the number of shared hits of the new cluster
	 
	  for ( iclhit = cldgcl.begin(); iclhit != cldgcl.end(); iclhit++){   // loop on the clusters the shared hits belongs to
	    flags |= (*iclhit)->ContainsFlaggedChannels();                 // "compute" flags value of the "new" cluster
 
	    std::list<CsMumegaHit*> cldgcldg = (*iclhit)->GetHits();   // get hits of these clusters
	   
	    for (ihitclhit = cldgcldg.begin(); ihitclhit != cldgcldg.end();  ihitclhit++){  // loop on them
	    
	      if (dglist.size()==0) {dglist.push_back((*ihitclhit));  // store the hits in the list 
	      }
	      else {
		Bool_t include = true;
		for (idglist = dglist.begin(); idglist != dglist.end(); idglist++){
		  if( (*idglist) == (*ihitclhit) ) include = false;
		}
		if (include){
		  dglist.push_back((*ihitclhit));
		}
	      }
	      
	    }    // end of loop on the hits of the clusters the shared hits belongs to

	    (*iclhit)->OldClus();                         // tag cluster as "old" cluster 
	    (*iclhit)->DeacrNSharedHits();                // decrease its number of shared it
	    NSharedHits += (*iclhit)->GetNSharedHits();   // compute the number of shared hits of the "new" cluster
	    (*ihit)->DeacNrClusters(); 
	    
	  }    // end of loop on the clusters the shared hits belongs to
	  
	  (*ihit)->GetClusters().clear(); // clear list of clusters of the shared hit

	  CsMumegaCluster *cl = new CsMumegaCluster(dglist, NSharedHits, flags);  // create a new cluster with the hits of the old ones 
	  fClusters.push_back(cl);                                 // add it to the plane
	  cl->SortHitAmps();                                       // sort its hits by amplitude
	  cl->NewClus();                                           // tag it as "new" cluster
	  dglist.clear();                                          // clear temporary list of hits
	  (*ihit)->AddCluster(cl);                                 // add the new cluster to the former shared hit
	  (*ihit)->IncNrClusters();                                // and increase its nb of clusters
	  if ((*ihit)->GetNrClusters()>1) cl->IncrNSharedHits();   // in case the hit is still shared, increase the nb of shared hits of the cluster (might be useful later for pixels, after a few adaptations)

	}    // end if amplitude ratio is above threshold

      }    // end if hit is a shared hit

    }    // end loop on clusters with shared hits

  }    // end loop on clusters of the event
  
}


void CsMumegaPlane::GetKey(Int_t event, Int_t key, Int_t keysym, TObject* _selected) {
  // Only for keyboard events
  if (event != kKeyPress) return;
  // printf("symbol = %c (%x)\n", key, keysym);
  fKey = keysym;
}


void CsMumegaPlane::Display(Int_t sample) {
  // Create histogram if it does not exist
  if (fDisplayHist == NULL)
    fDisplayHist = new TH1F(fName.c_str(), fName.c_str(), 768, -0.5, 767.5);
  else
    fDisplayHist->Reset();

  // Fill hits
  std::list<CsMumegaHit*>::iterator ithit;
  for (ithit = fHits.begin(); ithit != fHits.end(); ithit++) {
    int strip = (*ithit)->GetChan()->GetId()->GetDetectorChannel();
    float amp = (*ithit)->GetAmp()[sample];
    fDisplayHist->Fill(strip, amp);
  }

  // Cluster positions
  int ncl = fClusters.size();
  int icl;
  float xcl[ncl], ycl[ncl];
  std::list<CsMumegaCluster*>::iterator itclu;
  for (itclu = fClusters.begin(), icl = 0; itclu != fClusters.end(); itclu++, icl++) {
    xcl[icl] = (*itclu)->GetPosition();

    // Strips belonging to cluster
    std::list<CsMumegaHit*> clhits = (*itclu)->GetHits();
    ycl[icl] = 0.;
    float amp = 0.;
    for (ithit = clhits.begin(); ithit != clhits.end(); ithit++) {
      amp = (*ithit)->GetAmp()[sample];
      if (amp > ycl[icl]) ycl[icl] = amp;
    } // end of loop over strips of clusters
  } // end of loop over clusters

  // Indicate cluster positions
  TPolyMarker *pm = new TPolyMarker(ncl, xcl, ycl);
  pm->SetMarkerStyle(23);
  pm->SetMarkerColor(kRed);
  pm->SetMarkerSize(1.3);

  // A valid non-batch gApplication object must exist, the default created with "new TCanvas" is a batch one
  if (!gApplication)
    new TApplication("monitor", 0, 0);

  // Open canvas if not yet open (!connect signal to SetKey method)
  if (fDisplayCanvas == NULL) {
    fDisplayCanvas = new TCanvas(fName.c_str(), fName.c_str(), 500, 500);
    fDisplayCanvas->Connect("ProcessedEvent(Int_t, Int_t, Int_t, TObject*)", "CsMumegaPlane", this, "GetKey(Int_t, Int_t, Int_t, TObject*)");
  }

  // Plot
  fDisplayCanvas->cd();
  fDisplayHist->Draw();
  pm->Draw();
  fDisplayCanvas->Modified();
  fDisplayCanvas->Update();

  // Check for user input in canvas
  gROOT->SetInterrupt(kFALSE);
  bool endLoop = false;
  if (fKey == 0x73) fKey = 0; // stop at this plane if "s" had been pressed
  do {
    // Process ROOT events
    gSystem->ProcessEvents();

    // check keystrokes in fDisplayCanvas
    switch (fKey) {
    case 0x6e:   // "n" : next event (plot all planes until this plane comes again)
      fKey = 0;
      endLoop = true;
      break;
      case 0x63:   // "c" : continue without stopping
	endLoop = true;
	break;
    case 0x73:   // "s" : step to next plane (set fKey to zero for all planes)
      fKey = 0;
      endLoop = true;
      break;
    }
  } while (!gROOT->IsInterrupted()&&!endLoop);

  // Clean up
  pm->Delete();
}


void CsMumegaPlane::PrintChannels() {
  std::cout<<"============================================================"<<std::endl;
  std::cout<<fName<<" "<<fKey<<std::endl;

  // Channels
  std::map<CsMumegaChanId, CsMumegaChan*>::iterator itchan;
  std::cout<<"  List of channels"<<std::endl;
  std::cout<<"    Ch  Hem  Flag Ped  Sig Neighbours"<<std::endl;
  for (itchan = fChannels.begin(); itchan != fChannels.end(); itchan++) {
    std::list<CsMumegaChan*>::iterator itnb;
    std::list<CsMumegaChan*> nb = itchan->second->GetActiveNeighbours();
    std::cout<< "    "
	     <<itchan->second->GetId()->GetDetectorChannel()<<" "
	     <<itchan->second->GetId()->GetPosition()<<" "
	     <<itchan->second->GetCal()->GetFlag()<<" "
	     <<itchan->second->GetCal()->GetPedMean()<<" "
	     <<itchan->second->GetCal()->GetPedSigma()<<" ";
    for (itnb = nb.begin(); itnb != nb.end(); itnb++) {
      std::cout<<(*itnb)->GetId()->GetDetectorChannel()<<" "
	       <<(*itnb)->GetId()->GetPosition()<<" ";
    }
    std::cout<<std::endl;
  }
}

void CsMumegaPlane::PrintClusters(Int_t sample) {
  std::cout<<"============================================================"<<std::endl;
  std::cout<<fName<<" "<<fKey<<std::endl;

  // Cluster positions
  std::list<CsMumegaCluster*>::iterator itclu;
  std::cout<<std::endl;
  std::cout<<"  List of clusters"<<std::endl;
  std::cout<<"    Pos Hem Amp0 Amp1 Amp2 Hits"<<std::endl;
  for (itclu = fClusters.begin(); itclu != fClusters.end(); itclu++) {
    std::cout<<"   "<<(*itclu)->GetPosition()
	     <<" "<<(*itclu)->GetHemisphere()
	     <<" "<<(*itclu)->GetAmp()[0]
	     <<" "<<(*itclu)->GetAmp()[1]
	     <<" "<<(*itclu)->GetAmp()[2]<<" ";

    // Strips belonging to cluster
    std::list<CsMumegaHit*> clhits = (*itclu)->GetHits();
    std::list<CsMumegaHit*>::iterator ithit;
    float amp = 0.;
    for (ithit = clhits.begin(); ithit != clhits.end(); ithit++) {
      amp = (*ithit)->GetAmp()[sample];
      std::cout<<(*ithit)->GetChan()->GetId()->GetDetectorChannel()
	       <<" "<<amp<<" ";
    } // end of loop over strips of cluster
    std::cout<<std::endl;
  } // end of loop over clusters
  std::cout<<std::endl;
}


void CsMumegaPlane::FindXTalk() { 
  const std::vector<float> &ctalktcuti = GetPar()->GetTimeCrossTalkHitRatioi();
  const std::vector<float> &ctalktcutip1 = GetPar()->GetTimeCrossTalkHitRatioip1();

  int sample = GetPar()->GetSample();

  // loop over hits
  std::list<CsMumegaHit*>::iterator ithit = fHits.begin();
  while (ithit != fHits.end() ) {
    // loop over hits a second time
    std::list<CsMumegaHit*>::iterator ithit2 = ithit; ithit2++;
    while ( ithit2 != fHits.end() ) {
      //      std::cout<<"chiphit1 "<<(*ithit)->GetChan()->GetId()->GetAPVChip()<<" "<<"chiphit2 "<<(*ithit2)->GetChan()->GetId()->GetAPVChip()<<std::endl;
      // Check for neighbouring channels
      /*      if ( abs((*ithit)->GetChan()->GetId()->GetDetectorChannel() - (*ithit2)->GetChan()->GetId()->GetDetectorChannel()) <= 1 &&
	   !(fabs((*ithit)->GetChan()->GetId()->GetPosX() - (*ithit2)->GetChan()->GetId()->GetDetectorChannel()) <= 1
	     && fabs((*ithit)->GetChan()->GetId()->GetPosY() - (*ithit2)->GetChan()->GetId()->GetPosY()) <= 1) ) {
	if ( (*ithit2)->GetAmp()[sample] < ctalkr * (*ithit2)->GetAmp()[sample] ) {
	  //cout << "Found crosstalk !";
	  (*ithit2)->SetXTalk(1.);
	  //cout << "marked" << endl;
	}
	if ( (*ithit)->GetAmp()[sample] < ctalkr * (*ithit2)->GetAmp()[sample] ) {
	  //cout << "Found crosstalk !";
	  (*ithit)->SetXTalk(1.);
	  //cout << "marked" << endl;
	}
	}*/

      // check for cross talk in time
      // both hits have to be on the same APV chip
      // make sure, we are not running with invalid pedestals or MC data

      const int &conn1 = (*ithit)->GetChan()->GetId()->GetConnNb();
      const int &conn2 = (*ithit2)->GetChan()->GetId()->GetConnNb();

      if ( conn1 != -1 && conn1 == conn2 ) {
	// now let us calculate the position of the hit in the multiplexed APV analog signal
	unsigned int posmux1 = (*ithit)->GetChan()->GetId()->GetAPVChannel();
	unsigned int posmux2 = (*ithit2)->GetChan()->GetId()->GetAPVChannel();
	
	// the multiplexer has a multiplicity of 6
	for (unsigned int i = 0; i < 6; i++) {
	  posmux1 = 32*(posmux1%4) + 8*(posmux1/4) - 31*(posmux1/16);
	  posmux2 = 32*(posmux2%4) + 8*(posmux2/4) - 31*(posmux2/16);
	}
	// neighbouring in the multiplexed signal
	if ( abs((int)posmux1 - (int)posmux2) == 1 ) {

	  if ( posmux1 > posmux2 ) {
	    std::vector<Float_t> vampi = (*ithit2)->GetAmp();
	    std::vector<Float_t> vampip1 = (*ithit)->GetAmp();
	    float sigmai = (*ithit2)->GetChan()->GetCal()->GetPedSigma();
	    float sigmaip1 = (*ithit)->GetChan()->GetCal()->GetPedSigma();
	    float thr = fPar.GetThrHit();
	    
	    Bool_t erase1 = false;
	    Bool_t erase2 = false;
	    
	    if (vampi[sample] -  ctalktcuti[conn1] * vampip1[sample] > thr*sigmai){
	      for (int i = 0; i<vampi.size(); i++){
		if (vampi[i] >= ctalktcuti[conn1] * vampip1[i]){
		  vampi[i] = vampi[i] - ctalktcuti[conn1] * vampip1[i];}
		else vampi[i] = 0.;
	      }
	      (*ithit2)->SetAmp(vampi);
	      (*ithit2)->SetMultiplexAmpRatio((vampi[sample]) / (vampip1[sample]));
	    }
	    else {
	      erase1 = true;
	    }

	    if (vampip1[sample] -  ctalktcutip1[conn1] * vampi[sample] > thr*sigmaip1){
	      for (int i = 0; i<vampi.size(); i++){
		if (vampip1[i] >= ctalktcutip1[conn1] * vampi[i]){
		  vampip1[i] = vampip1[i] - ctalktcutip1[conn1] * vampi[i];}
		else vampip1[i] = 0.;
	      }
	      (*ithit)->SetAmp(vampip1);
	    }
	    else {
	      erase2 = true;
	    }

	    if (erase1 || erase2){
	      if (erase1) {delete(*ithit2);ithit2 = fHits.erase(ithit2);}
	      if (erase2) {delete(*ithit);ithit = fHits.erase(ithit);}
	      continue;
	    }

	  } else if ( posmux2 > posmux1 ) {
	    std::vector<Float_t> vampi = (*ithit)->GetAmp();
	    std::vector<Float_t> vampip1 = (*ithit2)->GetAmp();
	    float sigmai = (*ithit)->GetChan()->GetCal()->GetPedSigma();
	    float sigmaip1 = (*ithit2)->GetChan()->GetCal()->GetPedSigma();
	    float thr = fPar.GetThrHit();
	    
	    Bool_t erase1 = false;
	    Bool_t erase2 = false;
	        
	    if (vampi[sample] - ctalktcuti[conn1] * vampip1[sample] > thr*sigmai){
	      for (int i = 0; i<vampi.size(); i++){
		if (vampi[i] >= ctalktcuti[conn1] * vampip1[i]){
		  vampi[i] = vampi[i] -  ctalktcuti[conn1] * vampip1[i];
		}
		else vampi[i] = 0.;
	      }
	      (*ithit)->SetAmp(vampi);
	      (*ithit)->SetMultiplexAmpRatio((vampi[sample]) / (vampip1[sample]));
	    }

	    else {
	      erase1 = true;
	    }

	    if (vampip1[sample] - ctalktcutip1[conn1] * vampi[sample] > thr*sigmaip1){
	      for (int i = 0; i<vampi.size(); i++){
		if (vampip1[i] >= ctalktcutip1[conn1] * vampi[i]){
		  vampip1[i] = vampip1[i] -  ctalktcutip1[conn1] * vampi[i];
		}
		else vampip1[i] = 0.;
	      }
	      (*ithit2)->SetAmp(vampip1);
	    }
	    else {
	      erase2 = true; 
	    }
	    
	    if (erase1 || erase2){
	      if (erase2) {delete(*ithit2);ithit2 = fHits.erase(ithit2);}
	      if (erase1) {delete(*ithit);ithit = fHits.erase(ithit);}
	      continue;
	    }
	  }
	}
      } // end of same APV chips
      ithit2++;
    } // end of loop over second hit
    // iterate to next hit
    ithit++;
  } // end loop over first hit
}
