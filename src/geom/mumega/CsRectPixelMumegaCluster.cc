
// include C++ headers
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <cstdlib>   // for abs()

#include "CsRectPixelMumegaCluster.h"
#include "CsRectPixelMumegaPlane.h"


CsRectPixelMumegaCluster::CsRectPixelMumegaCluster() {
  fSizeX = 0;
  fSizeY = 0;
  fPosX = 0.;
  fPosY = 0.;
  fPosXErr = 0.;
  fPosYErr = 0.;
  fEtaX = -1.;
  fEtaY = -1.;
  fNoise = 0.;
  fXTalk = 0.;
  fHasTime = false;
  fContainsFlaggedChannels = kFALSE;
  fNewClus = false;
  fOldClus = false;
  fNSharedHits = 0;
}


CsRectPixelMumegaCluster::CsRectPixelMumegaCluster(const std::list<CsMumegaHit*> &_hits, Int_t _NSharedHits, Bool_t _flags, const CsRectPixelMumegaPlanePar &_par) {
  
  fHasTime = false;

  // Save reference to hits in cluster
  fHits = _hits;

  // does the cluster contain flagged channels ?
  fContainsFlaggedChannels = _flags;

  fNSharedHits = _NSharedHits;

  // get clusterization parameters
  Int_t _sample = _par.GetSample();
  Float_t _thrhit = _par.GetThrHit();
  Int_t _share = (_par.GetClusConfigMask() >> 6) & 3;
  //std::vector<Float_t> _time_cal = _par.GetTimeCal();

  // Calculate cluster amplitude from hits
  CalcAmps(_sample, _thrhit, _share);

  // Calculate cluster position from hits
  CalcCoG(_sample, _thrhit, _share);

  // Calculate cluster time
  //CalcTime(_time_cal);
  fNewClus = false;
  fOldClus = false;

}


void CsRectPixelMumegaCluster::CalcCoG(int _sample, float _thr, int _share) {
  //std::cout<<"_share = "<<_share<<std::endl;
  std::vector<float> sum(3,0.);
  double sumt = 0.;
  double wsumx = 0.;
  double wsumy = 0.;
  double dcenter = 0.;
  float centerx = 0.;
  float centery = 0.;

  // Calculate position by center-of-gravity method
  // Loop over hits
  std::list<CsMumegaHit*>::iterator ithit;
  for ( ithit = fHits.begin(); ithit != fHits.end(); ithit++ ) {

    // Detector channel, hemisphere, sigma of hit
    float posX = (*ithit)->GetChan()->GetId()->GetPosX();
    float posY = (*ithit)->GetChan()->GetId()->GetPosY();
    float sigma = (*ithit)->GetChan()->GetCal()->GetPedSigma();
    //std::cout<<"Nrclus  hit to add ter"<<(*ithit)->GetNrClusters()<<std::endl;
    // Get the fraction of the hit amplitude that should be assigned to this cluster
    float frac = 1.;
    if (_share == 2)
      frac = 1. / (*ithit)->GetNrClusters();

    sumt += frac * ((*ithit)->GetAmp()[_sample] - _thr*sigma);
    wsumx += frac * ((*ithit)->GetAmp()[_sample] - _thr*sigma) * posX;
    wsumy += frac * ((*ithit)->GetAmp()[_sample] - _thr*sigma) * posY;
  }

  if (sumt <= 0.) {
    fPosX = 0.;
    fPosY = 0.;
    fPosXErr = 9999999.;
    fPosYErr = 9999999.;
    return;
  }

  // Center of gravity
  centerx = wsumx / sumt;
  centery = wsumy / sumt;

  // Error on center of gravity
  dcenter = 0.0833; // s / sqrt(12) ^ 2
  dcenter = sqrt(dcenter);

  // Update cluster properties
  fPosX = centerx;
  fPosY = centery;
  fPosXErr = dcenter;
  
  if (fHits.front()->GetChan()->GetId()->GetType()==2)  fPosYErr = 2*dcenter;
  else if (fHits.front()->GetChan()->GetId()->GetType()==3)  fPosYErr = 5*dcenter;

// std::cout<<"CalcCoG: _sample "<<_sample<<" centerx "<<centerx<<" wsumx "<<wsumx<<" centery "<<centery<<" wsumy "<<wsumy
//     <<" sumt "<<sumt<<" fPosX "<<fPosX<<" fPosY "<<fPosY<<std::endl;

  return;
}

void CsRectPixelMumegaCluster::CorrectPos(const std::vector<Float_t> &poscorrs) { //correction of cluster position from CoG

  //  const std::vector<Float_t> poscorrs = GetPar()->GetPosCorr();

  //test
  /*for (int i =0; i<poscorrs.size(); i++){
    std::cout<<"pos_corr"<<poscorrs[i]<<std::endl;
    }*/

  Float_t maxX = -25.6;
  Float_t minX = 25.6;
  Float_t maxY = -25.0;
  Float_t minY = 25.0;

  std::list<CsMumegaHit*>::iterator ihit;
  for (ihit = fHits.begin(); ihit != fHits.end(); ihit++){

    if ((*ihit)->GetChan()->GetId()->GetPosX() > maxX) maxX = (*ihit)->GetChan()->GetId()->GetPosX();
    if ((*ihit)->GetChan()->GetId()->GetPosX() < minX) minX = (*ihit)->GetChan()->GetId()->GetPosX();
    if ((*ihit)->GetChan()->GetId()->GetPosY() > maxY) maxY = (*ihit)->GetChan()->GetId()->GetPosY();
    if ((*ihit)->GetChan()->GetId()->GetPosY() < minY) minY = (*ihit)->GetChan()->GetId()->GetPosY();

  }
  // std::cout<<"u "<<u<<" minU "<<minU<<std::endl;
  //Eta variable
  fEtaX = fPosX - minX;
  fEtaY = fPosY - minY;
    
  // switch (size_u) 
  Float_t  param_x_0 = 0, param_x_1 = 0, param_x_2 = 0, param_x_3 = 0;
  Float_t  param_y_0 = 0, param_y_1 = 0, param_y_2 = 0, param_y_3 = 0;  
  
   
  switch (fSizeX) { // X correction
  case 2 : //if cluster is 2 pixels wide
    param_x_0 = poscorrs[0]; //order 0 coefficient
    param_x_1 = poscorrs[1];  //order 1 coefficient
    param_x_2 = poscorrs[2];  //order 2 coefficient
    param_x_3 = poscorrs[3]; //order 3 coefficient
    break;
  case 3 : //if cluster is 3 pixels wide
    param_x_0 = poscorrs[4];  //order 0 coefficient
    param_x_1 = poscorrs[5];  //order 1 coefficient
    param_x_2 = poscorrs[6];  //order 2 coefficient
    param_x_3 = poscorrs[7]; //order 3 coefficient
    break;
  case 4 : //if cluster is 4 pixels wide
    param_x_0 = poscorrs[8];  //order 0 coefficient
    param_x_1 = poscorrs[9];  //order 1 coefficient
    param_x_2 = poscorrs[10];  //order 2 coefficient
    param_x_3 = poscorrs[11]; //order 3 coefficient
    break;
  case 5 : //if cluster is 5 pixels wide
    param_x_0 = poscorrs[12];  //order 0 coefficient
    param_x_1 = poscorrs[13];  //order 1 coefficient
    param_x_2 = poscorrs[14];  //order 2 coefficient
    param_x_3 = poscorrs[15]; //order 3 coefficient
    break;
  default : //other cases (eg: 1 pixel)
    param_x_0 = 0.; 
    param_x_1 = 0.; 
    param_x_2 = 0.; 
    param_x_3 = 0.;
  }
  
  //std::cout<<"size_u "<<fSizeX<<"params pos corr "<<param_x_0<<" "<<param_x_1<<" "<<param_x_2<<" "<<param_x_3<<std::endl;

  /*
  switch (fSizeX) {
  case 2 :
    param_x_0 = 0.00446557; param_x_1 = 0.0570062; param_x_2 = -0.686702; param_x_3 = 1.10376;
    break;
  case 3 :
    param_x_0 = 0.0423281; param_x_1 = -0.286459; param_x_2 = 0.674426; param_x_3 = -0.570413;
    break;
  case 4 :
    param_x_0 = -0.00196475; param_x_1 = 0.147756; param_x_2 = -0.356912; param_x_3 = 0.201105;
    break;
  case 5 :
    param_x_0 = -0.214278; param_x_1 = 0.884343; param_x_2 = -1.13435; param_x_3 = 0.46416;
    break;
  default :
    param_x_0 = 0.; param_x_1 = 0.; param_x_2 = 0.; param_x_3 = 0.;
  }
  */

  /*
  switch (fSizeY) {
  case 2 :
    param_y_0 = 0; param_y_1 = 0; param_y_2 = 0; param_y_3 = 0;
    break;
  case 3 :
    param_y_0 = 0; param_y_1 = 0; param_y_2 = 0; param_y_3 = 0;
    break;
  case 4 :
    param_y_0 = 0; param_y_1 = 0; param_y_2 = 0; param_y_3 = 0;
    break;
  case 5 :
    param_y_0 = 0; param_y_1 = 0; param_y_2 = 0; param_y_3 = 0;
    break;
  default :
    param_y_0 = 0.; param_y_1 = 0.; param_y_2 = 0.; param_y_3 = 0.;
    }*/

  /*  if (meanpixsize==1) {
    param_y_0 = 0.117; 
    param_y_1 = -0.0868;
  }
  else if (meanpixsize==0)  {
    param_y_0 = 0.308; 
    param_y_1 = -0.0935;
    }*/

  //  param_y_0 = 0.308 - (0.308 - 0.117) * meanpixsize;
  // param_y_1 = -0.0935 - (-0.0935 + 0.0868) * meanpixsize;

  // apply correction
  fPosX = 10 * (0.1*fPosX + param_x_0 + param_x_1*fEtaX + param_x_2*fEtaX*fEtaX + param_x_3*fEtaX*fEtaX*fEtaX);
  //fPosY = 10 * (0.1*fPosY + param_y_0 + param_y_1*fEtaY + param_y_2*fEtaY*fEtaY + param_y_3*fEtaY*fEtaY*fEtaY);

}

bool CsRectPixelMumegaCluster::AddHit(CsMumegaHit* _hit, Int_t _clusconfmask, Float_t _thr, Int_t _sample) {
  // std::cout<<"Nrclus  hit to add "<<_hit->GetNrClusters()<<std::endl;

  std::list<CsMumegaChan*>::iterator itchan;
  std::list<CsMumegaHit*>::iterator ithit;
  bool insert = false; // should the current hit be added to the cluster ?
  bool changeflagged = false; // contains the expanded cluster flagged channels
  bool xplus = true; // does the cluster's x/y size increase ?
  bool yplus = true;

  // read the clusterization configuration mask
  bool clpossizesplit = (_clusconfmask >> 4) & 1; // cut on the distance of hit to current CoG
  bool clpossizemethod = (_clusconfmask >> 5) & 1; // cut radial (true) or else (false)
  Int_t share = (_clusconfmask >> 6) & 3; // hit sharing between clusters

  // Add if list of hits in cluster empty
  if ( fHits.size() == 0 )
    insert = true;
  else { // Check if new hit should be added to existing list of hits
    // Get list of active neighbours for new hit
    std::list<CsMumegaChan*> neighbours = _hit->GetChan()->GetActiveNeighbours();
    for ( itchan = neighbours.begin(); itchan != neighbours.end(); itchan++ ) {
      // Check if any hit in cluster is a neighbour of the new hit
      for ( ithit = fHits.begin(); ithit != fHits.end(); ithit++ ) {
	Float_t distxn = fabs( (*itchan)->GetId()->GetPosX() - (*ithit)->GetChan()->GetId()->GetPosX() );
	Float_t distyn = fabs( (*itchan)->GetId()->GetPosY() - (*ithit)->GetChan()->GetId()->GetPosY() );
	Float_t pitchxn = ((*itchan)->GetId()->GetPitchX() + (*ithit)->GetChan()->GetId()->GetPitchX() ) / 2.; //04/01/12 : changement + -> -
	Float_t pitchyn = ((*itchan)->GetId()->GetPitchY() + (*ithit)->GetChan()->GetId()->GetPitchY() ) / 2.;
	
	Float_t distx = fabs(_hit->GetChan()->GetId()->GetPosX() - (*ithit)->GetChan()->GetId()->GetPosX() );
	Float_t disty = fabs(_hit->GetChan()->GetId()->GetPosY() - (*ithit)->GetChan()->GetId()->GetPosY() );
	Float_t mpitchy = fabs(_hit->GetChan()->GetId()->GetPitchY() - (*ithit)->GetChan()->GetId()->GetPitchY() ) / 2.;

	if ((*itchan) == (*ithit)->GetChan()) {
	  insert = true;
	  
	  // std::cout<<"X neighb hit"<<(*ithit)->GetChan()->GetId()->GetPosX()<<"Y neighb hit"<<(*ithit)->GetChan()->GetId()->GetPosY()<<std::endl;
	  
	  // Check if there is a gap between neighbours
	  if ( distxn > pitchxn || distyn > pitchyn )
	    changeflagged = true;
	}
	if (!(distx > 0.)) xplus = false; // does the cluster's x size increase ?
	if (!(disty > mpitchy)) yplus = false; // does the cluster's y size increase ?
      } 
    }


    // hit and cluster are neighbours, but does the hit fulfill all conditions
    if (insert) {
      // calculate current cluster center and distance to current hit
      // std::cout<<"Nrclus  hit to add bis"<<_hit->GetNrClusters()<<std::endl;
      CalcAmps(_sample, _thr, share);
      CalcCoG(_sample, _thr, share);
      //   std::cout<<"fPosX "<<fPosX<<std::endl;
      //std::cout<<"fPosY "<<fPosY<<std::endl;
      float distx = fabs( _hit->GetChan()->GetId()->GetPosX() - fPosX );
      float disty = fabs( _hit->GetChan()->GetId()->GetPosY() - fPosY );
      float pitchx = _hit->GetChan()->GetId()->GetPitchX();
      float pitchy = _hit->GetChan()->GetId()->GetPitchY();
      float dists = distx*distx + disty*disty;
      float pitchs = pitchx*pitchx + pitchy*pitchy;

      // cluster position size splitting
      if (clpossizesplit) {
	if ( clpossizemethod && dists > pitchs ) insert = false;
	if ( clpossizemethod && distx > 0.81 ) insert = false;
	if ( !clpossizemethod && (distx > pitchx || disty > pitchy) ) insert = false;	
      }

      /*  //no Y clustering
	  if (_hit->GetChan()->GetId()->GetPosY() != fPosY) insert = false;*/

	   //limited Y clustering
      //if (disty > (_hit->GetChan()->GetId()->GetPitchY())) insert = false;
      
    }
  }

  // Add hit to cluster
  if ( insert) {

    fHits.push_back(_hit);
    if (changeflagged)
      fContainsFlaggedChannels = true;

    if(xplus) fSizeX++;
    if(yplus) fSizeY++;
   
    //   std::cout<<"end cluster"<<std::endl;
  }

  return insert;
}


Bool_t CsRectPixelMumegaCluster::operator< (const CsRectPixelMumegaCluster &cluster) const {
  if (fPosX == cluster.GetPositionX())
    return (fPosY < cluster.GetPositionY());

  return (fPosX < cluster.GetPositionX());
}


