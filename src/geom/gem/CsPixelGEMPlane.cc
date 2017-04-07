/*
--------------------------------------------------------------------------

 Implementation of classes to store/analyze events from COMPASS PixelGEM/Sil
 detectors

 Author: Markus Kr�mer       Juli 2006
         Bernhard Ketzer     June 2008

--------------------------------------------------------------------------
*/

// Class declarations
#include "CsPixelGEMPlane.h"

// C++ headers
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

// ROOT headers
#include <TApplication.h>
#include <TPolyMarker.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

ClassImp(CsPixelGEMPlane)

//-----------------------------------------------------------------------------
// CsPixelGEMPlane constructor
//-----------------------------------------------------------------------------
CsPixelGEMPlane::CsPixelGEMPlane()
{
  // Create a CsPixelGEMPlane object.
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


//-----------------------------------------------------------------------------
// CsPixelGEMPlane constructor
//-----------------------------------------------------------------------------
CsPixelGEMPlane::CsPixelGEMPlane(std::string _name)
{
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

//-----------------------------------------------------------------------------
// CsPixelGEMPlane destructor
//-----------------------------------------------------------------------------
CsPixelGEMPlane::~CsPixelGEMPlane()
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
void CsPixelGEMPlane::AddChan(int _channr, int _xpos, int _ypos, int _flag, float _ped, float _sigma, int _apvchip, int _apvchan)
{
    // Channel Id and channel calibration objects
    CsGEMChanId id = CsGEMChanId(_channr, _xpos, _ypos);
    id.SetAPVInfo(_apvchip, _apvchan);
    CsGEMChanCal cal = CsGEMChanCal(_flag, _ped, _sigma);

    // Return if channel with same Id already exists
    if (GetChan(id)!=NULL) {
      std::cout<<"CsPixelGEMPlane::AddChan() : Channel with same Id already exists! Det.chan="
               <<_channr<<", X="<<_xpos<<", Y="<<_ypos<<std::endl;
      return;
    }

    // Create new channel
    CsGEMChan* chan = new CsGEMChan(id, cal);

    // update list of active neighbours
    // loop over all existing channels
    std::map<CsGEMChanId,CsGEMChan*>::iterator itchan;
    for (itchan=fChannels.begin(); itchan!=fChannels.end(); ++itchan) {
        bool neighbour = false;

        // get x and y pos of (possible) neighbour channel
        int n_xpos = itchan->second->GetId()->GetPosX();
        int n_ypos = itchan->second->GetId()->GetPosY();

        // Direct neighbours?
        if ( (abs(n_xpos - _xpos)+abs(n_ypos - _ypos)) <= 1 )
            neighbour = true;
        else if ( abs(n_xpos - _xpos)==1 && abs(n_ypos - _ypos)==1 && (GetPar()->GetClusConfigMask()&2) )
            neighbour = true;

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
const CsGEMChan* CsPixelGEMPlane::GetChan(const CsGEMChanId& _id) const {
    std::map<CsGEMChanId,CsGEMChan*>::const_iterator it = fChannels.find(_id);

    if ( it == fChannels.end() )
        return NULL;

    return it->second;
}

//-----------------------------------------------------------------------------
// AddHit member function
//-----------------------------------------------------------------------------
void CsPixelGEMPlane::AddHit(int _detchan, int _px, int _py, float _amp0, float _amp1, float _amp2)
{
  std::vector<float> amp;
  amp.push_back(_amp0);
  amp.push_back(_amp1);
  amp.push_back(_amp2);
  AddHit(_detchan, _px, _py, amp);
}

//-----------------------------------------------------------------------------
// AddHit member functions
//-----------------------------------------------------------------------------
void CsPixelGEMPlane::AddHit(int _detchan, int _px, int _py, const std::vector<float>& _amp)
{
    // Get pointer to corresponding CsGEMChan
    CsGEMChanId id = CsGEMChanId(_detchan, _px, _py);
    const CsGEMChan* chan = GetChan(id);

    // Channel not found
    if (chan==NULL) {
        std::cout<<"CsPixelGEMPlane::AddHit() : chan "<<_detchan<<", Xpos "<<_px<<", Ypos "<<_py<<" does not exist in "<<fName<<"! Ignoring it"<<std::endl;
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

    CsGEMHit *hit = new CsGEMHit(chan, _amp);

    // Calculate hit time
    if (!fTimeCals.IsValid()) {
        std::cout<<"CsPixelGEMPlane::AddHit() : time calibration not valid for plane "<<fName<<std::endl;
    } else {
        double time, etime;
        if ( fTimeCals.CalcTime(_amp, sigma, time, etime) )
            hit->SetTime(time, etime);
    }
    // Create new hit
    fHits.push_back(hit);
}

//-----------------------------------------------------------------------------
// Comparison function to sort list of pointers
//-----------------------------------------------------------------------------
struct CsPixelGEMPlane::fCompareHitAmps0 : public std::binary_function<CsGEMHit*, CsGEMHit*, bool> {
    bool operator() (CsGEMHit *a1, CsGEMHit *a2) {
        if (!a1) return false;
        if (!a2) return true;
        return (a1->GetAmp()[0] > a2->GetAmp()[0]);
    }
};
struct CsPixelGEMPlane::fCompareHitAmps1 : public std::binary_function<CsGEMHit*, CsGEMHit*, bool> {
    bool operator() (CsGEMHit *a1, CsGEMHit *a2) {
        if (!a1) return false;
        if (!a2) return true;
        return (a1->GetAmp()[1] > a2->GetAmp()[1]);
    }
};
struct CsPixelGEMPlane::fCompareHitAmps : public std::binary_function<CsGEMHit*, CsGEMHit*, bool> {
    bool operator() (CsGEMHit *a1, CsGEMHit *a2) {
        if (!a1) return false;
        if (!a2) return true;
        return (a1->GetAmp()[2] > a2->GetAmp()[2]);
    }
};


//-----------------------------------------------------------------------------
// Sort hits by fAmp3
//-----------------------------------------------------------------------------
void CsPixelGEMPlane::SortHitAmps()
{
  if (fIsHitAmpSorted) return;

  // Sort list of hits by amplitude
  if      ( fPar.GetSample()==0 ) fHits.sort(fCompareHitAmps0());
  else if ( fPar.GetSample()==1 ) fHits.sort(fCompareHitAmps1());
  else if ( fPar.GetSample()==2 ) fHits.sort(fCompareHitAmps());
  else                            std::cout << "Not a valid sample selection: " << fPar.GetSample() << std::endl;

#if 0
  std::cout << "Sorted: ";
  for (std::list<CsGEMHit*>::iterator i=fHits.begin(); i!=fHits.end(); i++) {
    std::cout << (*i)->Amp3() << " ";
  }
  std::cout << std::endl;
#endif

  // Set flag indicating that event has been sorted
  fIsHitAmpSorted = kTRUE;
}

//-----------------------------------------------------------------------------
// Clear lists of hits and clusters
//-----------------------------------------------------------------------------
void CsPixelGEMPlane::Clear()
{
    // Delete all hits in list
    ClearHits();

    // Delete all clusters in list
    ClearClusters();
}

//-----------------------------------------------------------------------------
// Clear map of channels
//-----------------------------------------------------------------------------
void CsPixelGEMPlane::ClearChannels()
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
void CsPixelGEMPlane::ClearHits()
{
    // Delete all objects in list
    std::list<CsGEMHit*>::iterator it=fHits.begin();
    while (it!=fHits.end()) {
        delete *it;
        it=fHits.erase(it);
    }

    fIsHitAmpSorted = kFALSE;
}

//-----------------------------------------------------------------------------
// Clear list of clusters
//-----------------------------------------------------------------------------
void CsPixelGEMPlane::ClearClusters()
{
  // Delete all objects in list
  std::list<CsPixelGEMCluster*>::iterator it=fClusters.begin();
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
void CsPixelGEMPlane::SetHeader(unsigned int date, unsigned int run, unsigned int evtNrRun, unsigned int evtNrSpill) {
   fHeader.Set(date, run, evtNrRun, evtNrSpill);
}

//-----------------------------------------------------------------------------
// Clusterize member function: switches between different methods depending
// on method set in fPar.
//-----------------------------------------------------------------------------
void CsPixelGEMPlane::Clusterize()
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
    std::list<CsPixelGEMCluster*>::iterator itclus;
    for (itclus=fClusters.begin(); itclus!=fClusters.end(); ++itclus) {
        (*itclus)->CalcAmps();
        (*itclus)->CalcCoG();

        // Position correction
        (*itclus)->CorrectPos();

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
// FullClustering member function, shared pad clustering.
// The hits of the plane are first sorted
// by amplitudes (descending). Than beginning by the first hit
// it is controled, if there is a cluster next to it. If there is
// the hit is assumed to be part of the cluster, else a new cluster
// is created and the hit is asumed to be part of the new one.
//-----------------------------------------------------------------------------
void CsPixelGEMPlane::FullClustering()
{
    // Sort hits by amplitude (descending)
    if (!fIsHitAmpSorted) SortHitAmps();

    // mark cross talk hits
    FindXTalk();

    // Iterators
    std::list<CsGEMHit*>::iterator ithit;
    std::list<CsPixelGEMCluster*>::iterator itclus;

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
            CsPixelGEMCluster* cl = new CsPixelGEMCluster(&fPar);
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
// SimpleClustering member function, primitive clustering: 1 hit = 1 cluster
//-----------------------------------------------------------------------------
void CsPixelGEMPlane::SimpleClustering()
{
  // Iterators
  std::list<CsGEMHit*>::iterator ithit;

  // Loop over hits in plane (assumed to be sorted by strip number)
  for (ithit=fHits.begin(); ithit!=fHits.end(); ++ithit) {
    // Add cluster
    CsPixelGEMCluster *cl = new CsPixelGEMCluster(&fPar);
    cl->AddHit( *ithit );
    fClusters.push_back(cl);
    (*ithit)->IncNrClusters();
  } // end of loop over hits in plane
}

//-----------------------------------------------------------------------------
// Select good clusters and discard bad ones
//-----------------------------------------------------------------------------
void CsPixelGEMPlane::SelectClusters() {
    float thrclus = fPar.GetThrClus(); // threshold on cluster amplitude
    int sample = fPar.GetSample(); // sample used for clustering

    // Loop over clusters
    std::list<CsPixelGEMCluster*>::iterator itclus = fClusters.begin();
    while ( itclus!=fClusters.end() ) {
        bool remove = false; // remove current cluster

        //Remove cluster if below threshold
        double amp = (*itclus)->GetAmp()[sample];
        double noise = (*itclus)->GetNoise();
        if ( !(amp>thrclus*noise) )
            remove=true;

        if ( remove ) {
            delete (*itclus);
            itclus = fClusters.erase(itclus);
        } else
            ++itclus;
    }
    return;
}

//-----------------------------------------------------------------------------
// method to gett key and print keysymbol of TCanvas event
//-----------------------------------------------------------------------------
void CsPixelGEMPlane::GetKey(int event, int /*key*/, int keysym, TObject* /*_selected*/)
{
  // Only for keyboard events
  if(event != kKeyPress) return;
  //  printf("symbol = %c (%x)\n", key, keysym);
  fKey = keysym;
}


//-----------------------------------------------------------------------------
// CsPixelGEMPlane::Display member function to plot a single event
//-----------------------------------------------------------------------------
void CsPixelGEMPlane::Display(bool wait)
{
  gStyle->SetPalette(1);
  gStyle->SetPadRightMargin(0.14);
  gStyle->SetOptStat(kFALSE);

  // Create histogram if it doesn't exits
  if (fDisplayHist==NULL)
    fDisplayHist = new TH2F(fName.c_str(), fName.c_str(), 32, -0.5, 31.5, 32, -0.5, 31.5);
  else
    fDisplayHist->Reset();

  // Fill hits
  std::list<CsGEMHit*>::iterator ithit;
  for (ithit=fHits.begin(); ithit!=fHits.end(); ++ithit) {
    float px  = (*ithit)->GetChan()->GetId()->GetPosX();
    float py  = (*ithit)->GetChan()->GetId()->GetPosY();
    float amp = (*ithit)->GetAmp()[fPar.GetSample()];
    fDisplayHist->Fill(px, py, amp);
  }

  // Cluster positions
  int ncl = fClusters.size();
  int icl;
  float xcl[ncl], ycl[ncl];
  std::list<CsPixelGEMCluster*>::iterator itclu;
  for (itclu=fClusters.begin(),icl=0; itclu!=fClusters.end(); ++itclu, ++icl) {
    xcl[icl] = (*itclu)->GetPositionX();
    ycl[icl] = (*itclu)->GetPositionY();
  } // end of loop over clusters

  // Indicate cluster positions
  TPolyMarker *pm = new TPolyMarker(ncl, xcl, ycl);
  pm->SetMarkerStyle(8);
  pm->SetMarkerColor(kBlack);
  pm->SetMarkerSize(1.3);

  // a valid non-batch gApplication object must exist, the default created with "new TCanvas" is a batch one
  if ((!gApplication) || (gApplication && gApplication->TestBit(TApplication::kDefaultApplication)))
    new TApplication("monitor", 0, 0);

  // Open canvas if not yet open (!connect signal to SetKey method)
  if (fDisplayCanvas==NULL) {
    fDisplayCanvas = new TCanvas(fName.c_str(), fName.c_str(), 500, 500);

    if (wait)
      fDisplayCanvas->Connect("ProcessedEvent(int,int,int,TObject*)","CsPixelGEMPlane",this,"GetKey(int,int,int,TObject*)");
  }

  // Plot
  fDisplayCanvas->cd();
  fDisplayHist->SetOption("COLZ");
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
// PrintChannels member function to print all channels of a CsPixelGEMPlane
//-----------------------------------------------------------------------------
void CsPixelGEMPlane::PrintChannels()
{
  std::cout << "=============================================================" << std::endl;
  std::cout << fName << " " << fKey << std::endl;

  // Channels
  std::map<CsGEMChanId,CsGEMChan*>::iterator itchan;
  std::cout<<"  List of channels"<<std::endl;
  std::cout<<"    Ch  X  Y  Flag Ped  Sig  Neighbours"<<std::endl;
  for (itchan=fChannels.begin(); itchan!=fChannels.end(); ++itchan) {
    std::list<CsGEMChan*>::iterator itnb;
    std::list<CsGEMChan*> nb = itchan->second->GetActiveNeighbours();
    std::cout<<"    "
             <<itchan->second->GetId()->GetDetectorChannel()<<" "
             <<itchan->second->GetId()->GetPosX()<<" "
             <<itchan->second->GetId()->GetPosY()<<" "
             <<itchan->second->GetCal()->GetFlag()<<" "
             <<itchan->second->GetCal()->GetPedMean()<<" "
             <<itchan->second->GetCal()->GetPedSigma()<<" ";
    for (itnb=nb.begin(); itnb!=nb.end(); ++itnb) {
      std::cout<<(*itnb)->GetId()->GetDetectorChannel()<<" "
               <<(*itnb)->GetId()->GetHemisphere()<<" ";
      }
    std::cout<<std::endl;
  }
}

//-----------------------------------------------------------------------------
// PrintClusters member function to print clusters of plane object
//-----------------------------------------------------------------------------
void CsPixelGEMPlane::PrintClusters()
{
  std::cout << "=============================================================" << std::endl;
  std::cout << fName << " " << fKey << std::endl;

  // Cluster positions
  std::list<CsPixelGEMCluster*>::iterator itclu;
  std::cout<<std::endl;
  std::cout<<"  List of clusters"<<std::endl;
  std::cout<<"    X/U Y/V Amp0 Amp1 Amp2 Hits"<<std::endl;
  for (itclu=fClusters.begin(); itclu!=fClusters.end(); ++itclu) {
    std::cout << "    " << (*itclu)->GetPositionX()
              << " " << (*itclu)->GetPositionY()
              << " " << (*itclu)->GetAmp()[0]
              << " " << (*itclu)->GetAmp()[1]
              << " " << (*itclu)->GetAmp()[2] << " ";

    // Pads belonging to cluster
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

//---------------------------------------------------------------------------------
// Function to identify crosstalk and set flag for crosstalk hits
//---------------------------------------------------------------------------------
void CsPixelGEMPlane::FindXTalk() {
    float ctalkr = GetPar()->GetCrossTalkHitRatio();
    float ctalkt = GetPar()->GetTimeCrossTalkHitRatio();
    int sample = GetPar()->GetSample();

    // std::cout << "Begin crosstalk identification." << std::endl;
    // loop over hits
    std::list<CsGEMHit*>::iterator ithit = fHits.begin();
    while ( ithit != fHits.end() ) {
        // loop over hits a second time
        std::list<CsGEMHit*>::iterator ithit2 = ithit; ++ithit2;
        while ( ithit2 != fHits.end() ) {
            // check for neighbouring channels
            if ( abs((*ithit)->GetChan()->GetId()->GetDetectorChannel() - (*ithit2)->GetChan()->GetId()->GetDetectorChannel()) <= 1 &&
                 !(fabs((*ithit)->GetChan()->GetId()->GetPosX() - (*ithit2)->GetChan()->GetId()->GetPosX()) <= 1
                   && fabs((*ithit)->GetChan()->GetId()->GetPosY() - (*ithit2)->GetChan()->GetId()->GetPosY()) <= 1) ) {
                if ( (*ithit2)->GetAmp()[sample] < ctalkr * (*ithit)->GetAmp()[sample] ) {
                    // std::cout << "Found crosstalk!";
                    (*ithit2)->SetXTalk(1.);
                    // std::cout << "marked" << std::endl;
                }
                if ( (*ithit)->GetAmp()[sample] < ctalkr * (*ithit2)->GetAmp()[sample] ) {
                    // std::cout << "Found crosstalk!";
                    (*ithit)->SetXTalk(1.);
                    // std::cout << "marked" << std::endl;
                }
            }

            // check for cross talk in time
            // both hits have to be on the same APV chip
            // make sure, we are not running with invalid pedestals or MC data
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

