/*
--------------------------------------------------------------------------

 Implementation of classes to store/cluster events from COMPASS GEM
 detectors

 Author: Bernhard Ketzer   27/11/2002    v0
                           20/05/2008    v1

--------------------------------------------------------------------------
*/

// Class declarations
#include "CsGEMPlane.h"

// C++ headers
#include <cstdlib>
#include <iostream>

// ROOT headers
#include <TApplication.h>
#include <TPolyMarker.h>
#include <TROOT.h>
#include <TSystem.h>

ClassImp(CsGEMPlane)

//-----------------------------------------------------------------------------
// CsGEMPlane constructor
//-----------------------------------------------------------------------------
CsGEMPlane::CsGEMPlane()
{
    // Create a CsGEMPlane object.
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

//-----------------------------------------------------------------------------
// CsGEMPlane constructor with name
//-----------------------------------------------------------------------------
CsGEMPlane::CsGEMPlane(std::string _name)
{
    // Create a CsGEMPlane object.
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
    fName = _name;
}

//-----------------------------------------------------------------------------
// CsGEMPlane destructor
//-----------------------------------------------------------------------------
CsGEMPlane::~CsGEMPlane()
{
  Clear();
  ClearChannels();
  fTimeCals.Clear();
  if (fDisplayCanvas!=NULL) fDisplayCanvas->Close();
  if (fDisplayHist!=NULL) fDisplayHist->Delete();
}

//-----------------------------------------------------------------------------
// AddChan member function
//-----------------------------------------------------------------------------
void CsGEMPlane::AddChan(int _channr, int _hem, int _flag, float _ped, float _sigma, int _apvchip, int _apvchan)
{
    // Channel Id and channel calibration objects
    CsGEMChanId id = CsGEMChanId(_channr, _hem);
    id.SetAPVInfo(_apvchip, _apvchan);
    CsGEMChanCal cal = CsGEMChanCal(_flag, _ped, _sigma);

    // Return if channel with same Id already exists
    if (GetChan(id)!=NULL) {
        std::cout<<"CsGEMPlane::AddChan() : Channel with same Id already exists! Det.chan="<<_channr<<", hemisphere="<<_hem<<std::endl;
        return;
    }

    // Create new channel
    CsGEMChan* chan = new CsGEMChan(id, cal);

    // Loop over all existing channels to update list of active neighbours
    std::map<CsGEMChanId,CsGEMChan*>::iterator itchan;
    for (itchan=fChannels.begin(); itchan!=fChannels.end(); ++itchan) {

        bool neighbour = false;

        // Direct neighbours
        if (abs(itchan->first.GetDetectorChannel()-_channr)<=1) neighbour = true;

        // for PixelGEM (names start with "GP") do not connect the two hemispheres of
        // strips [87,168] (82 strips), they are divided by the pixels
        if ( GetName().substr(0,2)=="GP" ) { // strip plane of PixelGEM detector
            if ( ( (_channr>=87 && _channr<=168) // either the to be added or the existing channel have to be in the range
                   || (itchan->first.GetDetectorChannel()>=87 && itchan->first.GetDetectorChannel()<=168) )
                 && itchan->first.GetHemisphere()!=_hem ) { // the hemispheres of added and existing differ
                // then strips have to be in different hemispheres, otherwise mapping is broken
                neighbour = false;
            }
        }

        // Other criteria for neighbours may be added here, e.g. if separated by flagged channels

        // Add to list if active neighbour
        if (neighbour) {
            if (itchan->second->GetCal()->GetFlag()!=0) chan->AddActiveNeighbour(itchan->second);
            if (_flag!=0) itchan->second->AddActiveNeighbour(chan);
        }
    }

    // Add new channel to map
    AddChan(chan);
}

//-----------------------------------------------------------------------------
// GetChan member function
//-----------------------------------------------------------------------------
const CsGEMChan* CsGEMPlane::GetChan(const CsGEMChanId& _id) const {
    std::map<CsGEMChanId,CsGEMChan*>::const_iterator it = fChannels.find(_id);

    if ( it == fChannels.end() )
        return NULL;

    return it->second;
}

//-----------------------------------------------------------------------------
// AddHit member function
//-----------------------------------------------------------------------------
void CsGEMPlane::AddHit(int _detchan, int _hem, float _amp0, float _amp1, float _amp2)
{
  std::vector<float> amp;
  amp.push_back(_amp0);
  amp.push_back(_amp1);
  amp.push_back(_amp2);
  AddHit(_detchan, _hem, amp);
}

//-----------------------------------------------------------------------------
// AddHit member function
//-----------------------------------------------------------------------------
void CsGEMPlane::AddHit(int _detchan, int _hem, const std::vector<float>& _amp)
{
    // Get pointer to corresponding CsGEMChan
    CsGEMChanId id = CsGEMChanId(_detchan, _hem);
    const CsGEMChan* chan = GetChan(id);

    // Channel not found
    if (chan==NULL) {
        std::cout<<"CsGEMPlane::AddHit() : chan "<<_detchan<<", hem "<<_hem<<" does not exist in "<<fName<<"! Ignoring it"<<std::endl;
        return;
    }

    // Get channel calibrations
    const CsGEMChanCal *cal = chan->GetCal();

    // Check flag
    int flag = cal->GetFlag();
    if (flag ==0) return;

    // Check amplitude
    float sigma = cal->GetPedSigma();
    float thr = fPar.GetThrHit();
    int sample = fPar.GetSample();
    if ( !(_amp[sample] > thr*sigma) ) return;

    // Create new hit
    CsGEMHit *hit = new CsGEMHit(chan, _amp);
    // Calculate hit time
    if (!fTimeCals.IsValid()) {
        std::cout<<"CsGEMPlane::AddHit() : time calibration not valid for plane "<<fName<<std::endl;
    } else {
        double time, timeerr;
        if (fTimeCals.CalcTime(_amp, sigma, time, timeerr))
            hit->SetTime(time, timeerr);
    }
    fHits.push_back(hit);
}

//-----------------------------------------------------------------------------
// Comparison function to sort list of hits by channel Id
//-----------------------------------------------------------------------------
struct CsGEMPlane::fCompareHits :
    public std::binary_function<CsGEMHit*, CsGEMHit*, bool>
{
    bool operator() (CsGEMHit *h1, CsGEMHit *h2)
        {
            if (!h1) return true;
            if (!h2) return false;
            return (h1->GetChan()->GetId() < h2->GetChan()->GetId());
        }
};

//-----------------------------------------------------------------------------
// SortHits member function
//-----------------------------------------------------------------------------
void CsGEMPlane::SortHits()
{
    // Sort list of hits by channel Id
    fHits.sort( fCompareHits() );

    // Set flag indicating that event has been sorted
    fIsHitSorted = kTRUE;
    fIsHitAmpSorted = kFALSE;
}

//-----------------------------------------------------------------------------
// Comparison function to sort list of hits by ascending amplitudes
//-----------------------------------------------------------------------------
struct CsGEMPlane::fCompareHitAmps0 :
    public std::binary_function<CsGEMHit*, CsGEMHit*, bool>
{
    bool operator() (CsGEMHit *h1, CsGEMHit *h2)
        {
            if (!h1) return false;
            if (!h2) return true;
            return ((h1->GetAmp())[0] > (h2->GetAmp())[0]);
        }
};

//-----------------------------------------------------------------------------
// Comparison function to sort list of hits by ascending amplitudes
//-----------------------------------------------------------------------------
struct CsGEMPlane::fCompareHitAmps1 :
    public std::binary_function<CsGEMHit*, CsGEMHit*, bool>
{
    bool operator() (CsGEMHit *h1, CsGEMHit *h2)
        {
            if (!h1) return false;
            if (!h2) return true;
            return ((h1->GetAmp())[1] > (h2->GetAmp())[1]);
        }
};

//-----------------------------------------------------------------------------
// Comparison function to sort list of hits by ascending amplitudes
//-----------------------------------------------------------------------------
struct CsGEMPlane::fCompareHitAmps2 :
    public std::binary_function<CsGEMHit*, CsGEMHit*, bool>
{
    bool operator() (CsGEMHit *h1, CsGEMHit *h2)
        {
            if (!h1) return false;
            if (!h2) return true;
            return ((h1->GetAmp())[2] > (h2->GetAmp())[2]);
        }
};

//-----------------------------------------------------------------------------
// SortHits member function
//-----------------------------------------------------------------------------
void CsGEMPlane::SortHitAmps()
{
    // Sort list of hits by channel Id
    int sample = fPar.GetSample();
    if (sample==0) fHits.sort( fCompareHitAmps0() );
    else if (sample==1) fHits.sort( fCompareHitAmps1() );
    else  fHits.sort( fCompareHitAmps2() );

    // Set flag indicating that event has been sorted
    fIsHitAmpSorted = kTRUE;
    fIsHitSorted = kFALSE;
}

//-----------------------------------------------------------------------------
// Comparison function to sort list of pointers
//-----------------------------------------------------------------------------
struct CsGEMPlane::fCompareClusters :
  public std::binary_function<CsGEMCluster*, CsGEMCluster*, bool>
{
  bool operator() (CsGEMCluster *c1, CsGEMCluster *c2)
  {
    if (!c1) return true;
    if (!c2) return false;
    return (c1->GetPosition() < c2->GetPosition());
  }
};

//-----------------------------------------------------------------------------
// SortClusters member function
//-----------------------------------------------------------------------------
void CsGEMPlane::SortClusters()
{
  // Sort list of clusters by cluster center
  fClusters.sort( fCompareClusters() );

  // Set flag indicating that event has been sorted
  fIsClusterSorted = kTRUE;
}

//-----------------------------------------------------------------------------
// Clear lists of hits and clusters
//-----------------------------------------------------------------------------
void CsGEMPlane::Clear()
{
    // Delete all hits in list
    ClearHits();

    // Delete all clusters in list
    ClearClusters();
}

//-----------------------------------------------------------------------------
// Clear map of channels
//-----------------------------------------------------------------------------
void CsGEMPlane::ClearChannels()
{
    // Delete all objects in map
    std::map<CsGEMChanId,CsGEMChan*>::iterator it;
    for (it=fChannels.begin(); it!=fChannels.end(); ++it)
        delete (it->second);
    fChannels.clear();
}

//-----------------------------------------------------------------------------
// Clear list of hits
//-----------------------------------------------------------------------------
void CsGEMPlane::ClearHits()
{
    // Delete all objects in list
    std::list<CsGEMHit*>::iterator it=fHits.begin();
    while (it!=fHits.end()) {
        delete *it;
        it=fHits.erase(it);
    }

    fIsHitSorted = kFALSE;
    fIsHitAmpSorted = kFALSE;
}

//-----------------------------------------------------------------------------
// Clear list of clusters
//-----------------------------------------------------------------------------
void CsGEMPlane::ClearClusters()
{
  // Delete all objects in list
  std::list<CsGEMCluster*>::iterator it=fClusters.begin();
  while (it!=fClusters.end()) {
    delete *it;
    it=fClusters.erase(it);
  }

  fIsClusterized = kFALSE;
  fIsClusterSorted = kFALSE;
}

//-----------------------------------------------------------------------------
// Set header
//-----------------------------------------------------------------------------
void CsGEMPlane::SetHeader(unsigned int date, unsigned int run, unsigned int evtNrRun, unsigned int evtNrSpill) {
   fHeader.Set(date, run, evtNrRun, evtNrSpill);
}

//-----------------------------------------------------------------------------
// Clusterize member function: switches between different methods depending
// on method set in fPar.
//-----------------------------------------------------------------------------
void CsGEMPlane::Clusterize()
{

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
    std::list<CsGEMCluster*>::iterator itclus;
    for (itclus=fClusters.begin(); itclus!=fClusters.end(); ++itclus) {
        // Cluster center and amplitude
        if (method==2 || (*itclus)->ContainsFlaggedChannels()) {
            // (*itclus)->Fit(thrhit);
            std::cout<<"CsGEMPlane::Clusterize() : Fitting of clusters with flagged channels not yet implemented!"<<std::endl;
            std::cout<<"                           Using CoG method instead..."<<std::endl;
            (*itclus)->CalcAmps();
            (*itclus)->CalcCoG();
        } else {
            (*itclus)->CalcAmps();
            (*itclus)->CalcCoG();
        }

        // Cluster time
        (*itclus)->CalcTime();
    }

    // Select only clusters above threshold
    SelectClusters();

    // Set flag indicating that event has been clusterized
    fIsClusterized = kTRUE;

    return;
}

//-----------------------------------------------------------------------------
// Full clustering method
//-----------------------------------------------------------------------------
void CsGEMPlane::FullClustering()
{

  // Sort hits by amplitude
  if (!fIsHitAmpSorted) SortHitAmps();

  // mark cross talk hits
  FindXTalk();

  // Iterators
  std::list<CsGEMHit*>::iterator ithit;
  std::list<CsGEMCluster*>::iterator itclus;

  // Group adjacent hits to cluster candidates
  for (ithit=fHits.begin(); ithit!=fHits.end(); ++ithit) {
    //    std::cout<<"hit "<<(*ithit)->GetChan()->GetId()->GetDetectorChannel()<<std::endl;

    bool assigned = false;

    // Loop over existing clusters
    for (itclus=fClusters.begin(); itclus!=fClusters.end(); ++itclus) {

      // Check if hit belongs to existing cluster, update its properties
      if ((*itclus)->AddHit(*ithit)) {
        assigned = true;
        (*ithit)->IncNrClusters();
      }
    }

    // Create new cluster if hit was not assigned to any existing cluster
    if (!assigned) {
      CsGEMCluster* cl = new CsGEMCluster(&fPar);
      cl->AddHit(*ithit);
      fClusters.push_back(cl);
      (*ithit)->IncNrClusters();
    }
  }

  // Cluster manipulations
  // MergeClusters(); // merge small clusters
  // SplitClusters(); // split large clusters

  return;
}

//-----------------------------------------------------------------------------
// Simple clustering: 1 hit = 1 cluster
//-----------------------------------------------------------------------------
void CsGEMPlane::SimpleClustering() {

    // Iterators
    std::list<CsGEMHit*>::iterator ithit;

    // Loop over hits in plane
    for (ithit=fHits.begin(); ithit!=fHits.end(); ++ithit) {

        CsGEMCluster* cl = new CsGEMCluster(&fPar);
        cl->AddHit(*ithit);
        fClusters.push_back(cl);
        (*ithit)->IncNrClusters();
    }

    return;
}

//-----------------------------------------------------------------------------
// Select good clusters and discard bad ones
//-----------------------------------------------------------------------------
void CsGEMPlane::SelectClusters() {
    float thrclus = fPar.GetThrClus(); // threshold on cluster amplitude
    int sample = fPar.GetSample(); // sample used for clustering

    // Loop over clusters
    std::list<CsGEMCluster*>::iterator itclus = fClusters.begin();
    while ( itclus!=fClusters.end() ) {
        //Remove cluster if below threshold
        double amp = (*itclus)->GetAmp()[sample];
        double noise = (*itclus)->GetNoise();

        if ( !(amp>thrclus*noise) ) {
            delete (*itclus);
            itclus = fClusters.erase(itclus);
        } else
            ++itclus;
    }
    return;
}

//-----------------------------------------------------------------------------
// Method to get key and keysymbol of TCanvas event
//-----------------------------------------------------------------------------
void CsGEMPlane::GetKey(int event, int /*key*/, int keysym, TObject* /*_selected*/)
{
  // Only for keyboard events
  if(event != kKeyPress) return;
  //  printf("symbol = %c (%x)\n", key, keysym);
  fKey = keysym;
}

//-----------------------------------------------------------------------------
// CsGEMPlane::Display member function to plot a single event
//-----------------------------------------------------------------------------
void CsGEMPlane::Display(bool wait)
{

  // Create histogram if it doesn't exist
  if (fDisplayHist==NULL)
    fDisplayHist = new TH1F(fName.c_str(), fName.c_str(), 768, -0.5, 767.5);
  else
    fDisplayHist->Reset();

  // Fill hits
  std::list<CsGEMHit*>::iterator ithit;
  for (ithit=fHits.begin(); ithit!=fHits.end(); ++ithit) {
    int strip = (*ithit)->GetChan()->GetId()->GetDetectorChannel();
    float amp = (*ithit)->GetAmp()[fPar.GetSample()];
    fDisplayHist->Fill(strip,amp);
  }

  // Cluster positions
  int ncl = fClusters.size();
  int icl;
  float xcl[ncl], ycl[ncl];
  std::list<CsGEMCluster*>::iterator itclu;
  for (itclu=fClusters.begin(),icl=0; itclu!=fClusters.end(); ++itclu, ++icl) {
    xcl[icl] = (*itclu)->GetPosition();

    // Strips belonging to cluster
    std::list<CsGEMHit*> clhits = (*itclu)->GetHits();
    ycl[icl] = 0.;
    float amp = 0.;
    for (ithit=clhits.begin(); ithit!=clhits.end(); ++ithit) {
        amp = (*ithit)->GetAmp()[fPar.GetSample()];
        if (amp>ycl[icl]) ycl[icl] = amp;
    } // end of loop over strips of cluster
  } // end of loop over clusters

  // Indicate cluster positions
  TPolyMarker *pm = new TPolyMarker(ncl, xcl, ycl);
  pm->SetMarkerStyle(23);
  pm->SetMarkerColor(kRed);
  pm->SetMarkerSize(1.3);

  // a valid non-batch gApplication object must exist, the default created with "new TCanvas" is a batch one
  if ((!gApplication) || (gApplication && gApplication->TestBit(TApplication::kDefaultApplication)))
    new TApplication("monitor", 0, 0);

  // Open canvas if not yet open (!connect signal to SetKey method)
  if (fDisplayCanvas==NULL) {
    fDisplayCanvas = new TCanvas(fName.c_str(), fName.c_str(), 500, 500);

    if (wait)
      fDisplayCanvas->Connect("ProcessedEvent(int,int,int,TObject*)","CsGEMPlane",this,"GetKey(int,int,int,TObject*)");
  }

  // Plot
  fDisplayCanvas->cd();
  fDisplayHist->Draw();
  pm->Draw();
  fDisplayCanvas->Modified();
  fDisplayCanvas->Update();

  // Check for user input in canvas
  if (wait) {
    gROOT->SetInterrupt(kFALSE);
    bool endLoop = false;
    if (fKey==0x73) fKey=0; // stop at this plane if "s" had been pressed
    do {
      // Process ROOT events
      gSystem->ProcessEvents();

      // Check keystrokes in fDisplayCanvas
      switch (fKey) {
        case 0x6e:   // "n": next event (plot all planes until this plane comes again)
          fKey = 0;
          endLoop = true;
          break;
        case 0x63:   // "c": continue without stopping
          endLoop = true;
          break;
        case 0x73:   // "s": step to next plane (set fKey to zero for all planes)
          fKey = 0;
          endLoop = true;
          break;
      }
    } while (!gROOT->IsInterrupted() && !endLoop);
  } // end if waits

  // Clean up
  pm->Delete();
}

//-----------------------------------------------------------------------------
// PrintChannels member function to print all channels of a CsGEMPlane
//-----------------------------------------------------------------------------
void CsGEMPlane::PrintChannels()
{
  std::cout << "=============================================================" << std::endl;
  std::cout << fName << " " << fKey << std::endl;

  // Channels
  std::map<CsGEMChanId,CsGEMChan*>::iterator itchan;
  std::cout<<"  List of channels"<<std::endl;
  std::cout<<"    Ch  Hem  Flag Ped  Sig  Neighbours"<<std::endl;
  for (itchan=fChannels.begin(); itchan!=fChannels.end(); ++itchan) {
    std::list<CsGEMChan*>::iterator itnb;
    std::list<CsGEMChan*> nb = itchan->second->GetActiveNeighbours();
    std::cout<<"    "
             <<itchan->second->GetId()->GetDetectorChannel()<<" "
             <<itchan->second->GetId()->GetPosition()<<" "
             <<itchan->second->GetCal()->GetFlag()<<" "
             <<itchan->second->GetCal()->GetPedMean()<<" "
             <<itchan->second->GetCal()->GetPedSigma()<<" ";
    for (itnb=nb.begin(); itnb!=nb.end(); ++itnb) {
      std::cout<<(*itnb)->GetId()->GetDetectorChannel()<<" "
               <<(*itnb)->GetId()->GetPosition()<<" ";
      }
    std::cout<<std::endl;
  }
}

//-----------------------------------------------------------------------------
// PrintClusters member function to print clusters of plane object
//-----------------------------------------------------------------------------
void CsGEMPlane::PrintClusters()
{
  std::cout << "=============================================================" << std::endl;
  std::cout << fName << " " << fKey << std::endl;

  // Cluster positions
  std::list<CsGEMCluster*>::iterator itclu;
  std::cout<<std::endl;
  std::cout<<"  List of clusters"<<std::endl;
  std::cout<<"    Pos Hem Amp0 Amp1 Amp2 Hits"<<std::endl;
  for (itclu=fClusters.begin(); itclu!=fClusters.end(); ++itclu) {
    std::cout << "    " << (*itclu)->GetPosition()
              << " " << (*itclu)->GetHemisphere()
              << " " << (*itclu)->GetAmp()[0]
              << " " << (*itclu)->GetAmp()[1]
              << " " << (*itclu)->GetAmp()[2] << " ";

    // Strips belonging to cluster
    std::list<CsGEMHit*> clhits = (*itclu)->GetHits();
    std::list<CsGEMHit*>::iterator ithit;
    float amp = 0.;
    for (ithit=clhits.begin(); ithit!=clhits.end(); ++ithit) {
      amp = (*ithit)->GetAmp()[fPar.GetSample()];
      std::cout << (*ithit)->GetChan()->GetId()->GetDetectorChannel()
                << " " << amp << " ";
    } // end of loop over strips of cluster
    std::cout<<std::endl;
  } // end of loop over clusters
  std::cout<<std::endl;
}

void CsGEMPlane::FindXTalk() {
    float ctalkt = GetPar()->GetTimeCrossTalkHitRatio();
    int sample = GetPar()->GetSample();

    // std::cout << "Begin crosstalk identification." << std::endl;
    // loop over hits
    std::list<CsGEMHit*>::iterator ithit = fHits.begin();
    while ( ithit != fHits.end() ) {
        // loop over hits a second time
        std::list<CsGEMHit*>::iterator ithit2 = ithit; ++ithit2;
        while ( ithit2 != fHits.end() ) {
            // check for cross talk in time
            // both hits have to be on the same APV chip
            // make sure, we are not running with invalid pedestals or MC data
            // and we do not need this for standard GEMs
            if ( (*ithit)->GetChan()->GetId()->GetAPVChip()!=-1 &&
                 (*ithit)->GetChan()->GetId()->GetAPVChip()==(*ithit2)->GetChan()->GetId()->GetAPVChip() ) {
                // now lets calculate the position of the hit in the multiplexed APV analog signal
                unsigned int posmux1 = (*ithit)->GetChan()->GetId()->GetAPVChannel();
                unsigned int posmux2 = (*ithit2)->GetChan()->GetId()->GetAPVChannel();

                // the multiplexer has a periodicity of 7
                for (unsigned int i=0; i<6; i++) {
                    posmux1 = 32*(posmux1%4) + 8*(posmux1/4) - 31*(posmux1/16);
                    posmux2 = 32*(posmux2%4) + 8*(posmux2/4) - 31*(posmux2/16);
                }

                // neighbouring in multiplexed signal
                if ( abs((int)posmux1-(int)posmux2)==1 ) {
                    if (posmux1>posmux2) {
                        if ( (*ithit)->GetAmp()[sample] < ctalkt * (*ithit2)->GetAmp()[sample] )
                            (*ithit)->SetXTalk(1.);
                    } else if (posmux2>posmux1) {
                        if ( (*ithit2)->GetAmp()[sample] < ctalkt * (*ithit)->GetAmp()[sample] )
                            (*ithit2)->SetXTalk(1.);
                    }
                }
            } // end of same APV chips

            ++ithit2;
        } // end of loop over second hit

        // iterate to next hit
        ++ithit;
    } // end loop over first hit;
}
