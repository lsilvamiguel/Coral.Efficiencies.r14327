/*
-------------------------------------------------------------------------

 Implementation of classes to store/analyze events from COMPASS rectangular PixelMumega
 detectors

 Author : Bernhard Ketzer    05/06/2009
 Modified by : Damien Neyret  1/02/2011
-------------------------------------------------------------------------
*/

#include <cstdlib>   // for abs()

#include "CsRectPixelMumegaPlane.h"

#include <TApplication.h>


ClassImp(CsRectPixelMumegaPlane)

  CsRectPixelMumegaPlane::CsRectPixelMumegaPlane() {
  // Create a CsRectPixelMumegaPlane object
  fChannels.clear();
  fHits.clear();
  fClusters.clear();
  fIsClusterized = kFALSE;
  fIsHitAmpSorted = kFALSE;
  fIsClusterSorted = kFALSE;
  fKey = 0;
  fDisplayCanvas = 0;
  fDisplayHist = 0;
  fName = "";

}


CsRectPixelMumegaPlane::CsRectPixelMumegaPlane(std::string _name) {
  fChannels.clear();
  fHits.clear();
  fClusters.clear();
  fIsClusterized = kFALSE;
  fIsHitAmpSorted = kFALSE;
  fIsClusterSorted = kFALSE;
  fKey = 0;
  fDisplayCanvas = 0;
  fDisplayHist = 0;
  fName = _name;
}


CsRectPixelMumegaPlane::~CsRectPixelMumegaPlane() {
  Clear();
  ClearChannels();
  //fTimeCals.Clear();
  if (fDisplayCanvas != NULL) fDisplayCanvas->Close();
  if (fDisplayHist != NULL) fDisplayHist->Delete();
}


void CsRectPixelMumegaPlane::AddChan(Int_t _chandet, Int_t _pixnb, Float_t _xpos, Float_t _ypos, Int_t _flag, Float_t _ped, Float_t _sigma, Int_t _apvchip, Int_t _apvchan, Int_t _connnb) {

  // Channel Id and channel calibration objects
  CsMumegaChanId id = CsMumegaChanId(_chandet, _pixnb, _xpos, _ypos);
  id.SetAPVInfo(_apvchip, _apvchan);
  id.SetRPConnNb(_connnb);
  CsMumegaChanCal cal = CsMumegaChanCal(_flag, _ped, _sigma);

  // Return if channel with same Id already exists
  if (GetChan(id) != NULL) {
    std::cout<<"CsRectPixelMumegaPlane::AddChan() : Channel with same Id already exists ! Det.chan = "
	     <<_chandet<<", X = "<<_xpos<<", Y = "<<_ypos<<std::endl;
    return;
  }

// std::cout<<"AddChan: _chandet "<<_chandet<<" _pixnb "<<_pixnb<<" _xpos "<<_xpos<<" _ypos "<<_ypos<<" _apvchip "<<_apvchip<<" _apvchan "<<_apvchan<<std::endl;
// std::cout<<"AddChan: id GetDetectorChannel "<<id.GetDetectorChannel()<<" GetPosition "<<id.GetPosition()
//   <<" GetPosX "<<id.GetPosX()<<" GetPosY "<<id.GetPosY()<<" GetType "<<id.GetType()<<std::endl;
  // Create new channel
  CsMumegaChan* chan = new CsMumegaChan(id, cal);

  //std::cout<<"_xpos"<<_xpos<<"_ypos"<<_ypos<<std::endl;

    // get pitchx & pitchy of the channel
  float _pitchx = chan->GetId()->GetPitchX();
  float _pitchy = chan->GetId()->GetPitchY();

  // update list of active neighbours
  // loop over all existing channels
  std::map<CsMumegaChanId, CsMumegaChan*>::iterator itchan;
  for (itchan = fChannels.begin(); itchan != fChannels.end(); itchan++) {
    bool neighbour = false;

    // get x and y pos of (possible) neighbour channel
    float n_xpos = itchan->second->GetId()->GetPosX();
    float n_ypos = itchan->second->GetId()->GetPosY();
    float n_pitchx = itchan->second->GetId()->GetPitchX();
    float n_pitchy = itchan->second->GetId()->GetPitchY();

    // Direct neighbours ?
    //  if ( fabs(n_xpos - _xpos) < 1.1*n_pitchx && fabs(n_ypos - _ypos) < 1.1*n_pitchy){

    if ( fabs(n_xpos - _xpos) < 0.55*(_pitchx + n_pitchx) && fabs(n_ypos - _ypos) < 0.55*(_pitchy + n_pitchy) /*&& ((fabs(n_xpos - _xpos) == 0) || (fabs(n_ypos - _ypos) <= 0.5*fabs(_pitchy - n_pitchy)))*/){ 
      //last condition added on 02/23/2012 to avoid "corner neighbouring"
      neighbour = true;
      // std::cout<<"pixnb "<<itchan->second->GetId()->GetPosition()<<"n_xpos "<<n_xpos<<"n_ypos "<<n_ypos<<"pitchx "<<n_pitchx<<"pitchy"<<n_pitchy<<std::endl;
    }
//     else if ( abs(n_xpos - _xpos)== 1 && abs(n_ypos - _ypos) == 1 && (GetPar()->GetClusConfigMask()&2) )
//       neighbour = true;

    // Other criteria for neighbours may be added here, e.g. if separated by flagged channels

    // Add to list if active neighbour
    if (neighbour) {
      if (itchan->second->GetCal()->GetFlag() != 0) chan->AddActiveNeighbour(itchan->second); //commenté le 21/12/11
      if (_flag != 0) itchan->second->AddActiveNeighbour(chan);
    }
    
 
  }
  //if (GetName() == "MP01MX__") std::cout<<GetName()<<" -- PIXELPLANE -- channb : "<<id.GetDetectorChannel()<<" ped : "<<cal.GetPedMean()<<" sigma : "<<cal.GetPedSigma()<<std::endl;
  // Add new channel to map
  AddChan(chan);
}


void CsRectPixelMumegaPlane::AddHit(Int_t _detchan, Int_t _pixnb, Float_t _px, Float_t _py, Float_t _amp0, Float_t _amp1, Float_t _amp2) {
  std::vector<float> amp;
  amp.push_back(_amp0);
  amp.push_back(_amp1);
  amp.push_back(_amp2);
  AddHit(_detchan, _pixnb, _px, _py, amp);
}


void CsRectPixelMumegaPlane::AddHit(Int_t _detchan, Int_t _pixnb, Float_t _px, Float_t _py, std::vector<Float_t> _amp) {
  // Get pointer to corresponding CsMumegaChan
  CsMumegaChanId id = CsMumegaChanId(_detchan, _pixnb, _px, _py);
  const CsMumegaChan* chan = GetChan(id);

  // Channel not found
  if (chan == NULL) {
    std::cout<<"CsRectPixelMumegaPlane::AddHit() : detchan "<<_detchan<<", pixnb "<<_pixnb<<", Xpos "<<_px
	     <<", Ypos "<<_py<<" does not exist in "<<fName<<" ! Ignoring it"
	     <<std::endl;
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
  if ( !(_amp[sample] >thr*sigma) ) return;

  CsMumegaHit *hit = new CsMumegaHit(chan, _amp);

  /*  // Calculate hit time
  if ( !fTimeCals.IsValid() ) {
    std::cout<<"CsRectPixelMumegaPlane::AddHit : time calibration not valid for plane "<<fName<<std::endl;
  } else {
    double time, etime;
    if ( fTimeCals.CalcTime(_amp, sigma, time, etime) )
      hit->SetTime(time, etime);
  }
  */
  // Create new hit
  fHits.push_back(hit);
}


struct CsRectPixelMumegaPlane::fCompareHitAmps0 : public std::binary_function<CsMumegaHit*, CsMumegaHit*, bool> {
  bool operator() (CsMumegaHit *a1, CsMumegaHit *a2) {
    if (!a1) return false;
    if (!a2) return true;
    return (a1->GetAmp()[0] > a2->GetAmp()[0]);
  }
};


struct CsRectPixelMumegaPlane::fCompareHitAmps1 : public std::binary_function<CsMumegaHit*, CsMumegaHit*, bool> {
  bool operator() (CsMumegaHit *a1, CsMumegaHit *a2) {
    if (!a1) return false;
    if (!a2) return true;
    return (a1->GetAmp()[1] > a2->GetAmp()[2]);
  }
};


struct CsRectPixelMumegaPlane::fCompareHitAmps : public std::binary_function<CsMumegaHit*, CsMumegaHit*, bool> {
  bool operator() (CsMumegaHit *a1, CsMumegaHit *a2) {
    if (!a1) return false;
    if (!a2) return true;
    return (a1->GetAmp()[2] > a2->GetAmp()[2]);
  }
};


void CsRectPixelMumegaPlane::SortHitAmps() {
  if (fIsHitAmpSorted) return;

  // Sort list of hits by amplitude
  if ( fPar.GetSample() == 0 ) fHits.sort(fCompareHitAmps0());
  else if ( fPar.GetSample() == 1 ) fHits.sort(fCompareHitAmps1());
  else if ( fPar.GetSample() == 2 ) fHits.sort(fCompareHitAmps());
  else std::cout << "Not a valid sample selection : " << fPar.GetSample()
		 << std::endl;

#if 0
  cout << "Sortedd : ";
  for ( std::list<CsMumegaHit*>::iterator i = fHits.begin(); i != fHits.end(); i++ ) {
    cout << (*i)->Amp3() << " ";
  }
  cout << endl;
#endif

  // Set flag indicating that event has been sorted
  fIsHitAmpSorted = kTRUE;
}


void CsRectPixelMumegaPlane::Clear(Option_t *option) {
  // Delete all hits in list
  ClearHits();

  // Delete all clusters in list
  ClearClusters();
}


void CsRectPixelMumegaPlane::ClearChannels(Option_t *option) {
  // Delete all objects in map
  std::map<CsMumegaChanId, CsMumegaChan*>::iterator it;
  for ( it = fChannels.begin(); it != fChannels.end(); it++ ) {
    delete ( it->second );
  }
  fChannels.clear();
}


void CsRectPixelMumegaPlane:: ClearClusters(Option_t *option) {
  // Delete all objects in list
  std::list<CsRectPixelMumegaCluster*>::iterator itclu;
  for ( itclu = fClusters.begin(); itclu != fClusters.end(); itclu++ ) {
    delete *itclu;
  }

  fClusters.clear();
  fIsClusterized = kFALSE;
  fIsClusterSorted = kFALSE;
}

void CsRectPixelMumegaPlane::ClearHits(Option_t *option)
{
    // Delete all objects in list
    std::list<CsMumegaHit*>::iterator it;
    for (it=fHits.begin(); it!=fHits.end(); it++) {
	delete *it;
    }
    fHits.clear();
    fIsHitAmpSorted = kFALSE;
}


void CsRectPixelMumegaPlane::SetHeader(Int_t i, Int_t run, Int_t date) {
  fHeader.Set(i, run, date);
}


void CsRectPixelMumegaPlane::Clusterize() {

  // Clear list of clusters
  ClearClusters();

  // Clustering method
  int method = fPar.GetClusMeth();
  switch (method) {

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
  int share = ( fPar.GetClusConfigMask() >> 6) & 3; // hit sharing ?
  std::vector<Float_t> time_cal_offset = fPar.GetTimeCalOffset();
  std::vector<Float_t> time_cal_slope = fPar.GetTimeCalSlope();
  std::vector<Float_t> time_cal_time_offset = fPar.GetTimeCalTimeOffset();
  const std::vector<float> &corrections = fPar.GetPosCorrs();
  std::list<CsRectPixelMumegaCluster*>::iterator itclus;

  itclus = fClusters.begin();
  while (itclus != fClusters.end()) {
    if ((*itclus)->IsOldClus()) {delete(*itclus); itclus = fClusters.erase(itclus); continue;}
    // Cluster center and amplitude

    // Position correction
    //if (GetPar()->IsPosCorrCal()) (*itclus)->CorrectPos(corrections);

    (*itclus)->CalcAmps(sample, thrhit, share);
    (*itclus)->CalcCoG(sample, thrhit, share);
    
    //   std::cout<<GetName()<<" offset "<<time_cal_offset[9]<<" slope "<<time_cal_slope[9]<<" time offset "<<time_cal_time_offset[9]<<std::endl;

	// Cluster time
    //if(GetName()=="MP01MX__")std::cout<<GetName()<<std::endl;
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


void CsRectPixelMumegaPlane::FullClustering() {
  // Sort hits by amplitude (descending)
  if ( !fIsHitAmpSorted ) SortHitAmps();

  // mark cross talk hits
  if (GetPar()->IsTxtCal()) FindXTalk();

  // Iterators
  std::list<CsMumegaHit*>::iterator ithit;
  std::list<CsRectPixelMumegaCluster*>::iterator itclus;

  // Group adjacent hits to cluster candidates
  for ( ithit = fHits.begin(); ithit != fHits.end(); ithit++ ) {
    //    std::cout << "hit " << (*ithit)->GetChan()->GetId()->GetDetectorChannel() << std::endl;

    bool assigned = false;

    // Loop over existing clusters
    for (itclus = fClusters.begin(); itclus != fClusters.end(); itclus++ ) {
      // Check if hit belongs to existing cluster, update its properties
      if ((*itclus)->AddHit(*ithit, GetPar()->GetClusConfigMask(), GetPar()->GetThrHit(), GetPar()->GetSample())) {
	assigned = true;
	(*ithit)->AddCluster(*itclus);
	(*ithit)->IncNrClusters();
      }
    }

    // Create new cluster if hit was not assigned to any existing cluster
    if ( !assigned ) {
      CsRectPixelMumegaCluster* cl = new CsRectPixelMumegaCluster();
      cl->AddHit(*ithit);
      fClusters.push_back(cl);
      (*ithit)->AddCluster(cl);
      (*ithit)->IncNrClusters();
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
  // if (fPar.GetMergingAmpratioThreshold() < 1) MergeClusters(); // merge small clusters
  // SplitClusters(); // split large clusters

  return;
}


void CsRectPixelMumegaPlane::SimpleClustering() {
  // Iterators
  std::list<CsMumegaHit*>::iterator ithit;

  // Loop over hits in plane (assumed to be sorted by strip number)
  for ( ithit = fHits.begin(); ithit != fHits.end(); ithit++ ) {
    // Add cluster
    CsRectPixelMumegaCluster *cl = new CsRectPixelMumegaCluster();
    cl->AddHit( *ithit );
    fClusters.push_back(cl);
    (*ithit)->IncNrClusters();
  } // end of loop over hits in plane
}


//-----------------------------------------------------------------------------
// Function which returns true if point xp,yp lies inside the
// polygon defined by the points in vectors x and y, false otherwise
// NOTE that the polygon must be a closed polygon (1st and last point
// must be identical) (taken from ROOT)
// Used in CsRectPixelMumegaPlane::SelectClusters() for cuts on amplitude ratios
//-----------------------------------------------------------------------------
bool CsRectPixelMumegaPlane::IsInside(float xp, float yp, std::vector<float> x, std::vector<float> y) {
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


void CsRectPixelMumegaPlane::SelectClusters() {
  float thrclus = fPar.GetThrClus(); // threshold on cluster amplitude
  int sample = fPar.GetSample(); // sample used for clustering
  double tcsphase = fHeader.GetTCSPhase();  //TCS Phase

  // Loop over clusters
  std::list<CsRectPixelMumegaCluster*>::iterator itclus = fClusters.begin();
  while ( itclus != fClusters.end() ) {
    bool remove = false; // remove current cluster

    double a0 = (*itclus)->GetAmp()[0];
    double a1 = (*itclus)->GetAmp()[1];
    double a2 = (*itclus)->GetAmp()[2];
    double amp = (*itclus)->GetAmp()[sample];
    

    // Remove cluster if below threshold
    double noise = (*itclus)->GetNoise();
        if ( !(amp > thrclus*noise ) )
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
      //std::cout<<GetName()<<" cutinf : "<<fPar.GetTimeCuts()[0]<<" time : "<<(*itclus)->GetTime() - tcsphase<<" cutsup : "<<fPar.GetTimeCuts()[1]<<std::endl;
  
    }
  
    if ( remove ) {
      delete (*itclus);
      itclus = fClusters.erase(itclus);
    } else
      itclus++;
  }
  return;
}


void CsRectPixelMumegaPlane::MergeClusters() { 

  
  //  std::list<CsMumegaHit*>::iterator ithit;    



  Float_t ampratiothr = fPar.GetMergingAmpratioThreshold(); 
  int sample = fPar.GetSample();

  std::list<CsRectPixelMumegaCluster*>::iterator iclus;    // iterator for list of clusters of the event
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

	  CsRectPixelMumegaCluster *cl = new CsRectPixelMumegaCluster(dglist, NSharedHits, flags, fPar);  // create a new cluster with the hits of the old ones 
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



void CsRectPixelMumegaPlane::GetKey(Int_t event, Int_t key, Int_t keysym, TObject* _selected) {

  // Only for keyboard events
  if ( event != kKeyPress ) return;
  fKey = keysym;
}


void CsRectPixelMumegaPlane::Display(Int_t sample, bool wait) {
  gStyle->SetPalette(1);
  gStyle->SetPadRightMargin(0.14);
  gStyle->SetOptStat(kFALSE);

  // Create histogram if it does not exist
  if ( fDisplayHist == NULL )
    fDisplayHist = new TH2F(fName.c_str(), fName.c_str(), 128, -25.6, 25.6, 40, -25, 25);
  else
    fDisplayHist->Reset();

  // Fill hits
  std::list<CsMumegaHit*>::iterator ithit;
  for ( ithit = fHits.begin(); ithit != fHits.end(); ithit++) {
    register float px = (*ithit)->GetChan()->GetId()->GetPosX();
    register float py = (*ithit)->GetChan()->GetId()->GetPosY();
    register float amp = (*ithit)->GetAmp()[sample];
    register int type = (*ithit)->GetChan()->GetId()->GetType();
    if (type == 2) {  // small rectangular pixel = 2*subpixel
      fDisplayHist->Fill(px, py+0.625, amp);
      fDisplayHist->Fill(px, py-0.625, amp);
    }
    if (type == 3) {  // large rectangular pixel=  5*subpixel
      fDisplayHist->Fill(px, py, amp);
      fDisplayHist->Fill(px, py-1.25, amp);
      fDisplayHist->Fill(px, py-2.5, amp);
      fDisplayHist->Fill(px, py+1.25, amp);
      fDisplayHist->Fill(px, py+2.5, amp);
    }
  }

  // Cluster positions
  int ncl = fClusters.size();
  int icl;
  float xcl[ncl], ycl[ncl];
  std::list<CsRectPixelMumegaCluster*>::iterator itclu;
  for ( itclu = fClusters.begin(), icl = 0; itclu != fClusters.end(); itclu++, icl++ ) {
    xcl[icl] = (*itclu)->GetPositionX();
    ycl[icl] = (*itclu)->GetPositionY();
  } // end of loop over clusters

  // Indicate cluster positions
  TPolyMarker *pm = new TPolyMarker(ncl, xcl, ycl);
  pm->SetMarkerStyle(8);
  pm->SetMarkerColor(kBlack);
  pm->SetMarkerSize(1.);


  // a valid non-batch gApplication object must exist, the default created with "new TCanvas" is a batch one
  if ( !gApplication )
    new TApplication("monitor", 0, 0);

  // Open canvas if not yet open ( !connect signal to SetKey method )
  if ( fDisplayCanvas == NULL ) {
    fDisplayCanvas = new TCanvas(fName.c_str(), fName.c_str(), 500, 500);

    if (wait)
      fDisplayCanvas->Connect("ProcessedEvent(Int_t, Int_t, Int_t, TObject*)", "CsRectPixelMumegaPlane", this, "GetKey(Int_t, Int_t, Int_t, TObject*)");
  }

  // Plot
  fDisplayCanvas->cd();
  fDisplayHist->SetOption("COLZ");
  fDisplayHist->Draw();
  pm->Draw();
  fDisplayCanvas->Modified();
  fDisplayCanvas->Update();


  if (wait) {

    gROOT->SetInterrupt(kFALSE);
    bool endLoop = false;
    if ( fKey == 0x73 ) fKey = 0; // stop at this plane if "s" had been pressed
    do {
      // Process ROOT events
      gSystem->ProcessEvents();

      // Check keystrokes in fDisplayCanvas
      switch (fKey) {
      case 0x6e: // "n" : next event ( plot all planes until this plane comes again )
	fKey = 0;
	endLoop = true;
	break;
      case 0x63: // "c" : continue without stopping
	endLoop = true;
	break;
      case 0x73: // "s" : step to next plane ( set fKey to zero for all planes )
	fKey = 0;
	endLoop = true;
	break;
      }
    } while ( !gROOT->IsInterrupted() && !endLoop );
  } // end if waits

  // Clean up
  pm->Delete();
}


void CsRectPixelMumegaPlane::PrintChannels(Int_t sample) {
  std::cout << "=============================================================" << std::endl;
  std::cout << fName << " " << fKey << std::endl;
  std::cout << "  List of channels" << std::endl;
  std::cout << "  detchan  posX posY flag pedmean pedsigma" << std::endl;

  // Channels
  std::map<CsMumegaChanId, CsMumegaChan*>::iterator itchan;
  for ( itchan = fChannels.begin(); itchan != fChannels.end(); itchan++ ) {
    std::list<CsMumegaChan*>::iterator itnb;
    std::list<CsMumegaChan*> nb = itchan->second->GetActiveNeighbours();
    std::cout << "   "
	      << itchan->second->GetId()->GetDetectorChannel() << " "
	      << itchan->second->GetId()->GetPosX() << " "
	      << itchan->second->GetId()->GetPosY() << " "
	      << itchan->second->GetCal()->GetFlag() << " "
	      << itchan->second->GetCal()->GetPedMean() << " "
	      << itchan->second->GetCal()->GetPedSigma() << " ";
    for ( itnb = nb.begin(); itnb != nb.end(); itnb++ ) {
      std::cout << (*itnb)->GetId()->GetDetectorChannel() << " "
		<< (*itnb)->GetId()->GetHemisphere() << " ";
    }
    std::cout << std::endl;
  }
}


void CsRectPixelMumegaPlane::PrintChannelsSimple(Int_t sample) {
  std::cout << "=============================================================" << std::endl;
  std::cout << fName << " " << fKey << std::endl;
  std::cout << "  List of channels" << std::endl;
  std::cout << "  detchan position posX posY type APVchip APVchannel" << std::endl;

  // Channels
  std::map<CsMumegaChanId, CsMumegaChan*>::iterator itchan;
  for ( itchan = fChannels.begin(); itchan != fChannels.end(); itchan++ ) {
    std::list<CsMumegaChan*>::iterator itnb;
    std::list<CsMumegaChan*> nb = itchan->second->GetActiveNeighbours();
    std::cout << "   "
	      << itchan->second->GetId()->GetDetectorChannel() << " "
	      << itchan->second->GetId()->GetPosition() << " "
	      << itchan->second->GetId()->GetPosX() << " "
	      << itchan->second->GetId()->GetPosY() << " "
	      << itchan->second->GetId()->GetType() << " "
	      << itchan->second->GetId()->GetAPVChip() << " "
	      << itchan->second->GetId()->GetAPVChannel() << " ";
    std::cout << std::endl;
  }
}


void CsRectPixelMumegaPlane::PrintClusters(Int_t sample) {
  std::cout << "=============================================================" << std::endl;
  std::cout << fName << " " << fKey << std::endl;

  // Cluster positions
  std::list<CsRectPixelMumegaCluster*>::iterator itclu;
  std::cout << std::endl;
  std::cout << "  List of clusters" << std::endl;
  std::cout << "  X/U Y/V Amp0 Amp1 Amp2 Time" << std::endl;
  for ( itclu = fClusters.begin(); itclu != fClusters.end(); itclu++ ) {
    std::cout << "- " << (*itclu)->GetPositionX()
	      << " " << (*itclu)->GetPositionY()
	      << " " << (*itclu)->GetAmp()[0]
	      << " " << (*itclu)->GetAmp()[1]
	      << " " << (*itclu)->GetAmp()[2]
	      << " " << (*itclu)->GetTime() << std::endl << "    ";

    // Pads belonging to cluster
    std::list<CsMumegaHit*> clhits = (*itclu)->GetHits();
    std::list<CsMumegaHit*>::iterator ithit;
    float amp = 0.;
    for ( ithit = clhits.begin(); ithit != clhits.end(); ithit++ ) {
      amp = (*ithit)->GetAmp()[sample];
      std::cout <<" ch "<< (*ithit)->GetChan()->GetId()->GetDetectorChannel()
		<< " amp " << amp << ",";
    } // end of loop over strips of cluster
    std::cout << std::endl;
  } // end of loop over clusters
  std::cout << std::endl;
}


void CsRectPixelMumegaPlane::FindXTalk() { 
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
