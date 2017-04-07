/*
-------------------------------------------------------------------------

 Implementation of classes to store/analyze events from COMPASS PixelMumega
 detectors

 Author : Bernhard Ketzer    05/06/2009
-------------------------------------------------------------------------
*/

#include <cstdlib>   // for abs()

#include "CsPixelMumegaPlane.h"

#include <TApplication.h>

ClassImp(CsPixelMumegaPlane)

  CsPixelMumegaPlane::CsPixelMumegaPlane() {
  // Create a CsPixelMumegaPlane object
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


CsPixelMumegaPlane::CsPixelMumegaPlane(std::string _name) {
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


CsPixelMumegaPlane::~CsPixelMumegaPlane() {
  Clear();
  ClearChannels();
  fTimeCals.Clear();
  if (fDisplayCanvas != NULL) fDisplayCanvas->Close();
  if (fDisplayHist != NULL) fDisplayHist->Delete();
}


void CsPixelMumegaPlane::AddChan(Int_t _channr, Int_t _xpos, Int_t _ypos, Int_t _flag, Float_t _ped, Float_t _sigma, Int_t _apvchip, Int_t _apvchan) {
  // Channel Id and channel calibration objects
  int pos = _xpos + (_ypos<<10);
  CsMumegaChanId id = CsMumegaChanId(_channr, pos);
  id.SetAPVInfo(_apvchip, _apvchan);
  CsMumegaChanCal cal = CsMumegaChanCal(_flag, _ped, _sigma);

  // Return if channel with same Id already exists
  if (GetChan(id) != NULL) {
    std::cout<<"CsPixelMumegaPlane::AddChan() : Channel with same Id already exists ! Det.chan = "
	     <<_channr<<", X = "<<_xpos<<", Y = "<<_ypos<<std::endl;
    return;
  }

  // Create new channel
  CsMumegaChan* chan = new CsMumegaChan(id, cal);

  // update list of active neighbours
  // loop over all existing channels
  std::map<CsMumegaChanId, CsMumegaChan*>::iterator itchan;
  for (itchan = fChannels.begin(); itchan != fChannels.end(); itchan++) {
    bool neighbour = false;

    // get x and y pos of (possible) neighbour channel
    Int_t n_xpos = (Int_t) itchan->second->GetId()->GetPosX();
    Int_t n_ypos = (Int_t) itchan->second->GetId()->GetPosY();

    // Direct neighbours ?
    if ( (abs(n_xpos - _xpos) + abs(n_ypos - _ypos)) <= 1)
      neighbour = true;
    else if ( abs(n_xpos - _xpos)== 1 && abs(n_ypos - _ypos) == 1 && (GetPar()->GetClusConfigMask()&2) )
      neighbour = true;

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


void CsPixelMumegaPlane::AddHit(Int_t _detchan, Int_t _px, Int_t _py, Float_t _amp0, Float_t _amp1, Float_t _amp2) {
  std::vector<float> amp;
  amp.push_back(_amp0);
  amp.push_back(_amp1);
  amp.push_back(_amp2);
  AddHit(_detchan, _px, _py, amp);
}


void CsPixelMumegaPlane::AddHit(Int_t _detchan, Int_t _px, Int_t _py, std::vector<Float_t> _amp) {
  // Get pointer to corresponding CsMumegaChan
  CsMumegaChanId id = CsMumegaChanId(_detchan, _px, _py);
  const CsMumegaChan* chan = GetChan(id);

  // Channel not found
  if (chan == NULL) {
    std::cout<<"CsPixelMumegaPlane::AddHit() : chan "<<_detchan<<", Xpos "<<_px
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
  /*
  // Calculate hit time
  if ( !fTimeCals.IsValid() ) {
    std::cout<<"CsPixelMumegaPlane::AddHit : time calibration not valid for plane "<<fName<<std::endl;
  } else {
    double time, etime;
    if ( fTimeCals.CalcTime(_amp, sigma, time, etime) )
      hit->SetTime(time, etime);
  }
  */
  // Create new hit
  fHits.push_back(hit);
}


struct CsPixelMumegaPlane::fCompareHitAmps0 : public std::binary_function<CsMumegaHit*, CsMumegaHit*, bool> {
  bool operator() (CsMumegaHit *a1, CsMumegaHit *a2) {
    if (!a1) return false;
    if (!a2) return true;
    return (a1->GetAmp()[0] > a2->GetAmp()[0]);
  }
};


struct CsPixelMumegaPlane::fCompareHitAmps1 : public std::binary_function<CsMumegaHit*, CsMumegaHit*, bool> {
  bool operator() (CsMumegaHit *a1, CsMumegaHit *a2) {
    if (!a1) return false;
    if (!a2) return true;
    return (a1->GetAmp()[1] > a2->GetAmp()[2]);
  }
};


struct CsPixelMumegaPlane::fCompareHitAmps : public std::binary_function<CsMumegaHit*, CsMumegaHit*, bool> {
  bool operator() (CsMumegaHit *a1, CsMumegaHit *a2) {
    if (!a1) return false;
    if (!a2) return true;
    return (a1->GetAmp()[2] > a2->GetAmp()[2]);
  }
};


void CsPixelMumegaPlane::SortHitAmps() {
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


void CsPixelMumegaPlane::Clear(Option_t *option) {
  // Delete all hits in list
  ClearHits();

  // Delete all clusters in list
  ClearClusters();
}


void CsPixelMumegaPlane::ClearChannels(Option_t *option) {
  // Delete all objects in map
  std::map<CsMumegaChanId, CsMumegaChan*>::iterator it;
  for ( it = fChannels.begin(); it != fChannels.end(); it++ ) {
    delete ( it->second );
  }
  fChannels.clear();
}


void CsPixelMumegaPlane:: ClearClusters(Option_t *option) {
  // Delete all objects in list
  std::list<CsPixelMumegaCluster*>::iterator itclu;
  for ( itclu = fClusters.begin(); itclu != fClusters.end(); itclu++ ) {
    delete *itclu;
  }

  fClusters.clear();
  fIsClusterized = kFALSE;
  fIsClusterSorted = kFALSE;
}

void CsPixelMumegaPlane::ClearHits(Option_t *option)
{
    // Delete all objects in list
    std::list<CsMumegaHit*>::iterator it;
    for (it=fHits.begin(); it!=fHits.end(); it++) {
	delete *it;
    }
    fHits.clear();
    fIsHitAmpSorted = kFALSE;
}


void CsPixelMumegaPlane::SetHeader(Int_t i, Int_t run, Int_t date) {
  fHeader.Set(i, run, date);
}


void CsPixelMumegaPlane::Clusterize() {

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
  std::vector<Float_t> time_cal = fPar.GetTimeCal();
  const std::vector<float> &corrections = fPar.GetPosCorrs();
  std::list<CsPixelMumegaCluster*>::iterator itclus;
  for ( itclus = fClusters.begin(); itclus != fClusters.end(); itclus++ ) {
    (*itclus)->CalcAmps(sample, thrhit, share);
    (*itclus)->CalcCoG(sample, thrhit, share);

    // Position correction
    (*itclus)->CorrectPos(corrections);

    // Cluster time
    //(*itclus)->CalcTime(time_cal);
  }

  // Select only clusters above threshold
  SelectClusters();

  // Set flag indicating that event has been clusterized
  fIsClusterized = kTRUE;

  return;
}


void CsPixelMumegaPlane::FullClustering() {
  // Sort hits by amplitude (descending)
  if ( !fIsHitAmpSorted ) SortHitAmps();

  // mark cross talk hits
  FindXTalk();

  // Iterators
  std::list<CsMumegaHit*>:: iterator ithit;
  std::list<CsPixelMumegaCluster*>::iterator itclus;

  // Group adjacent hits to cluster candidates
  for ( ithit = fHits.begin(); ithit != fHits.end(); ithit++ ) {
    //    std::cout << "hit " << (*ithit)->GetChan()->GetId()->GetDetectorChannel() << std::endl;

    bool assigned = false;

    // Loop over existing clusters
    for (itclus = fClusters.begin(); itclus != fClusters.end(); itclus++ ) {
      // Check if hit belongs to existing cluster, update its properties
      if ((*itclus)->AddHit(*ithit, GetPar()->GetClusConfigMask(), GetPar()->GetThrHit(), GetPar()->GetSample())) {
	assigned = true;
	(*ithit)->IncNrClusters();
      }
    }

    // Create new cluster if hit was not assigned to any existing cluster
    if ( !assigned ) {
      CsPixelMumegaCluster* cl = new CsPixelMumegaCluster();
      cl->AddHit(*ithit);
      fClusters.push_back(cl);
      (*ithit)->IncNrClusters();
    }
  }

  return;
}


void CsPixelMumegaPlane::SimpleClustering() {
  // Iterators
  std::list<CsMumegaHit*>::iterator ithit;

  // Loop over hits in plane (assumed to be sorted by strip number)
  for ( ithit = fHits.begin(); ithit != fHits.end(); ithit++ ) {
    // Add cluster
    CsPixelMumegaCluster *cl = new CsPixelMumegaCluster();
    cl->AddHit( *ithit );
    fClusters.push_back(cl);
    (*ithit)->IncNrClusters();
  } // end of loop over hits in plane
}


void CsPixelMumegaPlane::SelectClusters() {
  float thrclus = fPar.GetThrClus(); // threshold on cluster amplitude
  int sample = fPar.GetSample(); // sample used for clustering

  // Loop over clusters
  std::list<CsPixelMumegaCluster*>::iterator itclus = fClusters.begin();
  while ( itclus != fClusters.end() ) {
    bool remove = false; // remove current cluster

    // Remove cluster if below threshold
    double amp = (*itclus)->GetAmp()[sample];
    double noise = (*itclus)->GetNoise();
    if ( !(amp > thrclus*noise ) )
      remove = true;

    if ( remove ) {
      delete (*itclus);
      itclus = fClusters.erase(itclus);
    } else
      itclus++;
  }
  return;
}


void CsPixelMumegaPlane::GetKey(Int_t event, Int_t key, Int_t keysym, TObject* _selected) {

  // Only for keyboard events
  if ( event != kKeyPress ) return;
  fKey = keysym;
}


void CsPixelMumegaPlane::Display(Int_t sample, bool wait) {
  gStyle->SetPalette(1);
  gStyle->SetPadRightMargin(0.14);
  gStyle->SetOptStat(kFALSE);

  // Create histogram if it does not exist
  if ( fDisplayHist == NULL )
    fDisplayHist = new TH2F(fName.c_str(), fName.c_str(), 32, -0.5, 31.5, 32, -0.5, 31.5);
  else
    fDisplayHist->Reset();

  // Fill hits
  std::list<CsMumegaHit*>::iterator ithit;
  for ( ithit = fHits.begin(); ithit != fHits.end(); ithit++) {
    float px = (*ithit)->GetChan()->GetId()->GetPosX();
    float py = (*ithit)->GetChan()->GetId()->GetPosY();
    float amp = (*ithit)->GetAmp()[sample];
    fDisplayHist->Fill(px, py, amp);
  }

  // Cluster positions
  int ncl = fClusters.size();
  int icl;
  float xcl[ncl], ycl[ncl];
  std::list<CsPixelMumegaCluster*>::iterator itclu;
  for ( itclu = fClusters.begin(), icl = 0; itclu != fClusters.end(); itclu++, icl++ ) {
    xcl[icl] = (*itclu)->GetPositionX();
    ycl[icl] = (*itclu)->GetPositionY();
  } // end of loop over clusters

  // Indicate cluster positions
  TPolyMarker *pm = new TPolyMarker(ncl, xcl, ycl);
  pm->SetMarkerStyle(8);
  pm->SetMarkerColor(kBlack);
  pm->SetMarkerSize(1.3);


  // a valid non-batch gApplication object must exist, the default created with "new TCanvas" is a batch one
  if ( !gApplication )
    new TApplication("monitor", 0, 0);

  // Open canvas if not yet open ( !connect signal to SetKey method )
  if ( fDisplayCanvas == NULL ) {
    fDisplayCanvas = new TCanvas(fName.c_str(), fName.c_str(), 500, 500);

    if (wait)
      fDisplayCanvas->Connect("ProcessedEvent(Int_t, Int_t, Int_t, TObject*)", "CsPixelMumegaPlane", this, "GetKey(Int_t, Int_t, Int_t, TObject*)");
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


void CsPixelMumegaPlane::PrintChannels(Int_t sample) {
  std::cout << "=============================================================" << std::endl;
  std::cout << fName << " " << fKey << std::endl;

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


void CsPixelMumegaPlane::PrintClusters(Int_t sample) {
  std::cout << "=============================================================" << std::endl;
  std::cout << fName << " " << fKey << std::endl;

  // Cluster positions
  std::list<CsPixelMumegaCluster*>::iterator itclu;
  std::cout << std::endl;
  std::cout << "  List of clusters" << std::endl;
  std::cout << "    X/U Y/V Amp0 Amp1 Amp2 Hits" << std::endl;
  for ( itclu = fClusters.begin(); itclu != fClusters.end(); itclu++ ) {
    std::cout << "   " << (*itclu)->GetPositionX()
	      << " " << (*itclu)->GetPositionY()
	      << " " << (*itclu)->GetAmp()[0]
	      << " " << (*itclu)->GetAmp()[1]
	      << " " << (*itclu)->GetAmp()[2] << " ";

    // Pads belonging to cluster
    std::list<CsMumegaHit*> clhits = (*itclu)->GetHits();
    std::list<CsMumegaHit*>::iterator ithit;
    float amp = 0.;
    for ( ithit = clhits.begin(); ithit != clhits.end(); ithit++ ) {
      amp = (*ithit)->GetAmp()[sample];
      std::cout << (*ithit)->GetChan()->GetId()->GetDetectorChannel()
		<< " " << amp << " ";
    } // end of loop over strips of cluster
    std::cout << std::endl;
  } // end of loop over clusters
  std::cout << std::endl;
}


void CsPixelMumegaPlane::FindXTalk() {
  float ctalkr = GetPar()->GetCrossTalkHitRatio();
  float ctalkt = GetPar()->GetTimeCrossTalkHitRatio();
  int sample = GetPar()->GetSample();

  // loop over hits
  std::list<CsMumegaHit*>::iterator ithit = fHits.begin();
  while (ithit != fHits.end() ) {
    // loop over hits a second time
    std::list<CsMumegaHit*>::iterator ithit2 = ithit; ithit2++;
    while ( ithit2 != fHits.end() ) {
      // Check for neighbouring channels
      if ( abs((*ithit)->GetChan()->GetId()->GetDetectorChannel() - (*ithit2)->GetChan()->GetId()->GetDetectorChannel()) <= 1 &&
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
      }

      // check for cross talk in time
      // both hits have to be on the same APV chip
      // make sure, we are not running with invalid pedestals or MC data
      if ( (*ithit)->GetChan()->GetId()->GetAPVChip() != -1 &&
	   (*ithit)->GetChan()->GetId()->GetAPVChip() == (*ithit2)->GetChan()->GetId()->GetAPVChip() ) {
	// now let us calculate the position of the hit in the multiplexed APV analog signal
	unsigned int posmux1 = (*ithit)->GetChan()->GetId()->GetAPVChannel();
	unsigned int posmux2 = (*ithit2)->GetChan()->GetId()->GetAPVChannel();

	// the multiplexer has a multiplicity of ?
	for (unsigned int i = 0; i < 6; i++) {
	  posmux1 = 32*(posmux1%4) + 8*(posmux1/4) - 31*(posmux1/16);
	  posmux2 = 32*(posmux2%4) + 8*(posmux2/4) - 31*(posmux2/16);
	}

	// neighbouring in the multiplexed signal
	if ( abs((int)posmux1 - (int)posmux2) == 1 ) {
	  if ( posmux1 > posmux2 ) {
	    if ( (*ithit)->GetAmp()[sample] < ctalkt * (*ithit2)->GetAmp()[sample] )
	      (*ithit)->SetXTalk(1.);
	  } else if ( posmux2 > posmux1 ) {
	    if ( (*ithit2)->GetAmp()[sample] < ctalkt * (*ithit)->GetAmp()[sample] )
	      (*ithit2)->SetXTalk(1.);
	  }
	}
      } // end of same APV chips

      ithit2++;
    } // end of loop over second hit

    // iterate to next hit
    ithit++;
  } // end loop over first hit
}
