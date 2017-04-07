#ifndef CsPixelMumegaPlane_H
#define CsPixelMumegaPlane_H

/*
----------------------------------------------------------------------------

 Declaration of classe to store/analyze events from COMPASS PixelMumega
 detectors

 Author : Bernhard Ketzer        05/06/2009
----------------------------------------------------------------------------

class CsPixelMumegaPlane :

                 Object representing a plane of hit strips in one event.
		 The data member fhits is a list of pointers to CsMumegaHit
		 objects.
		 The data member fClusters is a list of pointers to
		 CsPixelMumegaCluster objects, each containing a list
		 of references to the CsMumegaHit objects it is made of.

----------------------------------------------------------------------------

class CsPixelMumegaPlaneHeader :

                 Header object
*/

// ROOT headers
#include "TMath.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPolyMarker.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TLine.h"
#include "RQ_OBJECT.h"
#include "KeySymbols.h"
#include "TH1.h"
#include "TH2.h"

// C++ headers
#include <list>
#include <set>
#include <map>
#include <iomanip>
#include <cassert>

// Mumega headers
#include "CsMumegaTimeCal.h"
#include "CsMumegaChan.h"
#include "CsMumegaHit.h"
#include "CsPixelMumegaCluster.h"

class CsPixelMumegaPlaneHeader {

 private:

  Int_t fEvtNum; // Event number
  Int_t fRun; // Run number
  Int_t fDate; // Date

 public:

  CsPixelMumegaPlaneHeader() : fEvtNum(0), fRun(0), fDate(0) { };

  virtual ~CsPixelMumegaPlaneHeader() { };

  void Set(Int_t _i, Int_t _r, Int_t _d) { fEvtNum = _i; fRun = _r; fDate = _d; };
  Int_t GetEvtNum() const { return fEvtNum; };
  Int_t GetRun() const { return fRun; };
  Int_t GetDate() const { return fDate; };
  void Print() {
    std::cout << std::endl << std::left << std::setw(20) << "EvtNum"
	      << std::setw(20) << "Run" << std::setw(20) << "Date"
	      << std::endl << std::setw(20) << fEvtNum
	      << std::setw(20) << fRun << std::setw(20) << fDate << std::endl;
  }

};


class CsPixelMumegaPlanePar {

 private:

  Int_t fClusMeth; // Clustering method : 0_primitive, 1_full, 2_always-fit
  Int_t fClusConfigMask; // Clustering configuration
  Float_t fThrClus; // Threshold on cluster amplitude
  Float_t fThrHit; // Threshold on single hit amplitude
  Float_t fCrossTalkR; // Maximal ratio of one hit to be marked cross talk
  Float_t fTimeCrossTalkR; // Maximal ratio of one hit to be marked cross talk in time
  Int_t fSample; // Sample to use for clustering

  std::vector<Float_t> fTimeCal;

  std::vector<Float_t> fPosCorrs;

 public:

  CsPixelMumegaPlanePar() : fClusMeth(1), fClusConfigMask(179), fThrClus(5.), fThrHit(3.), fCrossTalkR(0.3), fTimeCrossTalkR(0.), fSample(2) { }

  virtual ~CsPixelMumegaPlanePar() { }

  void Set(Int_t _meth = 1, Int_t _mask = 179, Float_t _thrclus = 5., Float_t _thrhit = 3., Float_t _ctalkr = 0.3, Float_t _ctalkt = 0., Int_t _sample = 2) {
    fClusMeth = _meth;
    fClusConfigMask = _mask;
    fThrClus = _thrclus;
    fThrHit = _thrhit;
    fCrossTalkR = _ctalkr;
    fTimeCrossTalkR = _ctalkt;
    fSample = _sample;
  }
  void SetClusMeth(Int_t _meth = 1) { fClusMeth = _meth; }
  void SetClusConfigMask(Int_t _mask = 179) { fClusConfigMask = _mask; }
  void SetThrClus(Float_t _thrclus = 5.) { fThrClus = _thrclus; }
  void SetThrHit(Float_t _thrhit = 3.) { fThrHit = _thrhit; }
  void SetCrossTalkHitRatio(Float_t _ctalkr = 0.3) { fCrossTalkR = _ctalkr; }
  void SetTimeCrossTalkHitRatio(Float_t _ctalkt = 0.) { fTimeCrossTalkR = _ctalkt; }
  void SetSample(Int_t _sample = 2) { fSample = _sample; }
  bool SetPosCorrs(const std::vector<Float_t> &_poscorrs) {
    fPosCorrs.clear();
    if (_poscorrs.size() == 16) {
      fPosCorrs = _poscorrs;
      return true;
    } else
      return false;
  }

  Int_t GetClusMeth() const { return fClusMeth; }
  Int_t GetClusConfigMask() const { return fClusConfigMask; }
  Float_t GetThrClus() const { return fThrClus; }
  Float_t GetThrHit() const { return fThrHit; }
  Float_t GetCrossTalkHitRatio() const { return fCrossTalkR; }
  Float_t GetTimeCrossTalkHitRatio() const { return fTimeCrossTalkR; }
  Int_t GetSample() const { return fSample; }

  const std::vector<Float_t>& GetPosCorrs() const { return fPosCorrs; }
  std::vector<Float_t> GetTimeCal() const {return fTimeCal; }
};


class CsPixelMumegaPlane {

 protected:

  CsPixelMumegaPlaneHeader fHeader; // plane header
  CsPixelMumegaPlanePar fPar; // plane parameters
  CsMumegaTimeCals fTimeCals; // time calibrations
  std::map<CsMumegaChanId, CsMumegaChan*> fChannels; // map of all channels including calibrations and active neighbours
  std::list<CsMumegaHit*> fHits; // list with pointers to hits
  std::list<CsPixelMumegaCluster*> fClusters; // list with pointers to clusters
  Bool_t fIsHitAmpSorted;
  Bool_t fIsClusterSorted;
  Bool_t fIsClusterized;
  std::string fName; // Name of the plane
  struct fCompareHitAmps0; // Needed to sort hits
  struct fCompareHitAmps1; // Needed to sort hits
  struct fCompareHitAmps; // Needed to sort hits
  TCanvas* fDisplayCanvas; //! ROOT canvas to display single event
  TH2F* fDisplayHist; //! ROOT histogram to display single event
  Int_t fKey; //! keystroke in fDisplayCanvas

 public:

  CsPixelMumegaPlane();
  CsPixelMumegaPlane(std::string name);

  virtual ~CsPixelMumegaPlane();

  virtual void Clear(Option_t *option = ""); // Clear all entries of plane
  void ClearChannels(Option_t* _option = "");
  void ClearHits(Option_t* _option = "");
  void ClearClusters(Option_t* _option = ""); // Clear all clusters
  void SetHeader(Int_t i, Int_t run, Int_t date);
  void SetPar(Int_t _meth = 1, Int_t _mask = 179, Float_t _thrclus = 5., Float_t _thrhit = 3., Float_t _ctalkr = 0.3, Int_t _sa = 2) { fPar.Set(_meth, _mask, _thrclus, _thrhit, _ctalkr, _sa); }
  void AddChan(Int_t _channr, Int_t _xpos, Int_t _ypos, Int_t _flag, Float_t ped, Float_t _sigma, Int_t _apvchip = -1, Int_t _apvchan = -1);
  void AddChan(CsMumegaChan* _chan) { fChannels[*(_chan->GetId())] = _chan; }
  void SetTimeCals(const CsMumegaTimeCals& _timecals) { fTimeCals = _timecals; }
  void SetTimeCals(const CsMumegaTimeCals* _timecals) { fTimeCals = *_timecals; }
  void AddHit(Int_t _detchan, Int_t _px, Int_t _py, std::vector<Float_t> _amp);
  void AddHit(Int_t _detchan, Int_t _px, Int_t _py, Float_t _amp0, Float_t _amp1, Float_t _amp2);
  void Clusterize();
  void SimpleClustering();
  void FullClustering();
  void SelectClusters();
  // Add a cluster to plane, defaut is no flagged pad
  void SortHitAmps(); // sort hits by amplitude
  Bool_t IsHitAmpSorted() const { return fIsHitAmpSorted; }
  Bool_t IsClusterSorted() const { return fIsClusterSorted; }
  Bool_t IsClusterized() const { return fIsClusterized; }
  Int_t GetNchannel() const { return fChannels.size(); }
  Int_t GetNhit() const { return fHits.size(); };
  Int_t GetNcluster() const { return fClusters.size(); };
  std::string GetName()const { return fName; }
  CsPixelMumegaPlaneHeader* GetHeader() { return &fHeader; }
  CsPixelMumegaPlanePar* GetPar() { return &fPar; }
  const CsMumegaTimeCals* GetTimeCals() const { return &fTimeCals; }
  const std::map<CsMumegaChanId, CsMumegaChan*> &GetChannels() const { return fChannels; }
  const std::list<CsMumegaHit*> &GetHits() const { return fHits; };
  const std::list<CsPixelMumegaCluster*> &GetClusters() const { return fClusters; };
  const CsMumegaChan* GetChan(CsMumegaChanId &_id) const {
    return (fChannels.find(_id) != fChannels.end() ?
	    fChannels.find(_id)->second : NULL);
  }
  void Display(Int_t sample = 2, bool wait = true);
  void GetKey(Int_t event, Int_t key, Int_t keysystem, TObject* _selected);
  void PrintChannels(Int_t _sample = 2);
  void PrintClusters(Int_t _sample = 2);
  virtual void FindXTalk();

  ClassDef(CsPixelMumegaPlane, 2)
    };

#endif
