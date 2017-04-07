/*
--------------------------------------------------------------------------

 Implementation of classes to store clusters from Mumegas
 detectors

 Author : Bernhard Ketzer    03/06/2009   v0

--------------------------------------------------------------------------
*/

#include <cassert>
#include <cstdlib>   // for abs()

// Class declarations
#include "CsMumegaCluster.h"
#include "CsMumegaHit.h"

CsMumegaCluster::CsMumegaCluster() {
  fPos = 0.;
  fPosErr = 0.;
  fNoise = 0.;
  fXTalk = 0.;
  fHasTime = false;
  fHemisphere = 0.;
  fContainsFlaggedChannels = kFALSE;
  fNewClus = false;
  fOldClus = false;
  fNSharedHits = 0;
}


CsMumegaCluster::CsMumegaCluster (const std::list<CsMumegaHit*> &_hits, Int_t _NSharedHits, Bool_t _flags, int _sample, float _thrhit, int _share/*, const CsMumegaPlanePar &_par*/) {
  // Save reference to hits in cluster
  fHits = _hits;

  // does the cluster contain flagged channels
  fContainsFlaggedChannels = _flags;

  fNSharedHits = _NSharedHits;

  //std::vector<Float_t> _time_cal = _par.GetTimeCal();

  // calculate cluster amplitude from hits
  CalcAmps(_sample, _thrhit, _share);

  // calculate cluster position from hits
  CalcCoG(_sample, _thrhit);

  // calculate cluster time
  //CalcTime(_time_cal);
  fNewClus = false;
  fOldClus = false;
}


bool CsMumegaCluster::AddHit(CsMumegaHit* _hit, Int_t _clusconfigmask, Float_t _thr, Int_t _sample) {
  std::list<CsMumegaChan*>::iterator itchan;
  std::list<CsMumegaHit*>::iterator ithit;
  bool insert = false;

  // Add if list of hits in cluster is empty
  if (fHits.size() == 0) insert = true;

  // Check if new hit should be added to existing list of hits
  else {

    // Get list of active neighbours for new hit
    std::list<CsMumegaChan*> neighbours = _hit->GetChan()->GetActiveNeighbours();
    for (itchan = neighbours.begin(); itchan != neighbours.end(); itchan++) {

      // Check if any hit in cluster is a neighbour of the new hit
      for (ithit = fHits.begin(); ithit != fHits.end(); ithit++) {
	if ((*itchan) == (*ithit)->GetChan()) {
	  insert = true;
	  //	  std::cout<<"CsMumegaCluster::AddHit() : hit "
	  //		   <<_hit->GetChan()->GetId()->GetDetectorChannel()
	  //		   <<"neighbour of "
	  //		   <<(*ithit)->GetChan()->GetId()->GetDetectorChannel()
	  //		   <<std::endl;

	  // Check if there is a gap between neighbours
	  if (abs((*itchan)->GetId()->GetDetectorChannel() -
		  (*ithit)->GetChan()->GetId()->GetDetectorChannel()) > 1)
	    fContainsFlaggedChannels = true;
	}
      }
    }
  }


  // Add hit to cluster
  if (insert) {
    fHits.push_back(_hit);
  }

  return insert;
}


void CsMumegaCluster::CalcAmps(int _sample, float _thr, int _share) {
  std::vector<float> sum(3,0.);
  double sumt = 0.;
  double wsumxt = 0.;
  double ssumsq = 0.;
  double ssum;
  double dcenter = 0.;

  // calculate position by center-of-gravity method
  // Loop over hits
  std::list<CsMumegaHit*>::iterator ithit;
  for (ithit = fHits.begin(); ithit != fHits.end(); ithit++) {

    // Detector channel, hemisphere, sigma of hit
    float sigma = (*ithit)->GetChan()->GetCal()->GetPedSigma();

    // get the fraction of the hit amplitude that should be assigned to this cluster
    float frac = 1.;
    if (_share == 2)
      frac = 1. / (*ithit)->GetNrClusters();

    // Sums
    for (int i = 0; i < 3; i++) {
      sum[i] += (*ithit)->GetAmp()[i]*frac; // amplitude sum for each of the samples
    }
    sumt += frac*((*ithit)->GetAmp()[_sample] - _thr*sigma);
    wsumxt += frac*((*ithit)->GetAmp()[_sample] - _thr*sigma) * (*ithit)->GetXTalk();
    ssumsq += sigma*sigma;

  }

  // Cluster noise
  ssum = sqrt(ssumsq/(double)fHits.size());

  // Cluster cross talk
  wsumxt /= sumt;

  fAmp = sum;
  fNoise = ssum;
  fXTalk = wsumxt;

  return;
}


void CsMumegaCluster::CalcCoG(int _sample, float _thr, int _share) {

  std::vector<float> sum(3,0.);
  float sigma;
  int detchan;
  double sumt = 0.;
  double wsumt = 0.;
  double dcenter = 0.;
  double hem;
  double whemt = 0.;
  double center = 0.;

  // Loop over hits
  std::list<CsMumegaHit*>::iterator ithit;
  for (ithit = fHits.begin(); ithit != fHits.end(); ithit++) {

    // Detector channel, hemisphere, sigma of hit
    detchan = (*ithit)->GetChan()->GetId()->GetDetectorChannel();
    hem = (*ithit)->GetChan()->GetId()->GetHemisphere();
    sigma = (*ithit)->GetChan()->GetCal()->GetPedSigma();

    // get the fractionof the hit amplitude that should be assigned to this cluster
    float frac = 1.;
    if (_share == 2)
      frac = 1. / (*ithit)->GetNrClusters();

    sumt += frac*((*ithit)->GetAmp()[_sample] - _thr*sigma);
    wsumt += frac*((*ithit)->GetAmp()[_sample] - _thr*sigma) * detchan;
    whemt += frac*((*ithit)->GetAmp()[_sample] - _thr*sigma)* hem;
  }

  if (sumt <= 0.) {
    fPos = 0.;
    fPosErr = 9999999.;
    fHemisphere = 0.;
    return;
  }

  // Center of gravity
  center = wsumt / sumt;

  // Error on center of gravity (assuem 50 um experimental resolution)
  //  dcenter = 0.015625;
  dcenter = 0.0833; // 1/sqrt(12) ^ 2
  //  double sumtsq = sum * sum;
  //  for (int j = 0; j <= nhit; j++) {
  //    dcenter += pow((x[j]*sumt - wsumt), 2)/pow(sumtsq, 4) *
  //      pow(ey[j], 2);
  //  }
  dcenter = sqrt(dcenter);

  // Weighed hemisphere
  whemt = whemt / sumt;

  // Update cluster properties
  fPos = center;
  fPosErr = dcenter;
  fHemisphere = whemt;

  return;
}


void CsMumegaCluster::CalcTime(const std::vector<Float_t> &_timecaloffset, const std::vector<Float_t> &_timecalslope, const std::vector<Float_t> &_timecaltimeoffset ) {

  fHasTime = false;
  fTime = 0.;
  fTimeErr = 1.e9;

  //double thit, ethit;

  /*
  // Loop over hits
  std::list <CsMumegaHit*>::iterator ithit;
  for (ithit = fHits.begin(); ithit != fHits.end(); ithit++) {
    if ( (*ithit)->GetTime(thit, ethit) ) {
      fHasTime = true;
      fTime += thit / (ethit*ethit);
      fTimeErr += 1. / (ethit*ethit);
    }
  }
  */

  CsMumegaHit *mainhit = fHits.front();
  double a1 = mainhit->GetAmp()[1];
  double a2 = mainhit->GetAmp()[2];
  int chipno = mainhit->GetChan()->GetId()->GetConnNb();

  fTime = ((a1/a2) - _timecaloffset[chipno]) / (_timecalslope[chipno]) - _timecaltimeoffset[chipno];
  fTimeErr = 20.;

  fHasTime = true;

  //std::cout<<"chipno "<<chipno<<"offset : "<<_timecaloffset[chipno]<<" - slope : "<<_timecalslope[chipno]<<" - time offset : "<< _timecaltimeoffset[chipno]<<std::endl;

  return;
}


bool CsMumegaCluster::GetTime(double &_time, double &_etime) const {
  _time = fTime;
  _etime = fTimeErr;

  return fHasTime;
}


struct CsMumegaCluster::fCompareHitAmps0 :
  public std::binary_function<CsMumegaHit*, CsMumegaHit*, bool> {
  bool operator() (CsMumegaHit *h1, CsMumegaHit *h2) {
    if (!h1) return false;
    if (!h2) return true;
    return ((h1->GetAmp())[0] > (h2->GetAmp())[0]);
  }
};


struct CsMumegaCluster::fCompareHitAmps1 :
  public std::binary_function<CsMumegaHit*, CsMumegaHit*, bool> {
  bool operator() (CsMumegaHit *h1, CsMumegaHit *h2) {
    if (!h1) return  false;
    if (!h2) return true;
    return ((h1->GetAmp())[1] > (h2->GetAmp())[1]);
  }
};


struct CsMumegaCluster::fCompareHitAmps2 :
  public std::binary_function<CsMumegaHit*, CsMumegaHit*, bool> {
  bool operator() (CsMumegaHit *h1, CsMumegaHit *h2) {
    if (!h1) return false;
    if (!h2) return true;
    return ((h1->GetAmp())[2] > (h2->GetAmp())[2]);
  }
};

void CsMumegaCluster::SortHitAmps(int sample) {
  // Sort list of hits by channel Id
  if (sample == 0) fHits.sort( fCompareHitAmps0() );
  else if (sample == 1) fHits.sort( fCompareHitAmps1() );
  else fHits.sort( fCompareHitAmps2() );

  // Set flag indicating that event has been sorted
  fIsHitAmpSorted = kTRUE;
}



Bool_t CsMumegaCluster::operator< (const CsMumegaCluster &cluster) const {
  return (fPos < cluster.GetPosition());
}
