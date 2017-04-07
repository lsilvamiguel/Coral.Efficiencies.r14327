#ifndef CsRectPixelMumegaPlane_H
#define CsRectPixelMumegaPlane_H

/*
----------------------------------------------------------------------------

 Declaration of classe to store/analyze events from COMPASS rectangular PixelMumega
 detectors

 Author : Bernhard Ketzer        05/06/2009
 Modifed by : Damien Neyret      1/02/2011
----------------------------------------------------------------------------

class CsRectPixelMumegaPlane :

                 Object representing a plane of hit strips in one event.
		 The data member fhits is a list of pointers to CsMumegaHit
		 objects.
		 The data member fClusters is a list of pointers to
		 CsRectPixelMumegaCluster objects, each containing a list
		 of references to the CsMumegaHit objects it is made of.

----------------------------------------------------------------------------

class CsRectPixelMumegaPlaneHeader :

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
#include "CsRectPixelMumegaCluster.h"

class CsRectPixelMumegaPlaneHeader {

 private:

  Int_t fEvtNum; // Event number
  Int_t fRun; // Run number
  Int_t fDate; // Date
  Float_t fTCSPhase; // TCSPhase time

 public:

  CsRectPixelMumegaPlaneHeader() : fEvtNum(0), fRun(0), fDate(0), fTCSPhase(-100) { };

  virtual ~CsRectPixelMumegaPlaneHeader() { };

  void Set(Int_t _i, Int_t _r, Int_t _d, Float_t _t = -100) { fEvtNum = _i; fRun = _r; fDate = _d; fTCSPhase = _t; };
  Int_t GetEvtNum() const { return fEvtNum; };
  Int_t GetRun() const { return fRun; };
  Int_t GetDate() const { return fDate; };
  Float_t GetTCSPhase() const { return fTCSPhase; };
  void Print() {
    std::cout << std::endl << std::left << std::setw(20) << "EvtNum"
	      << std::setw(20) << "Run" << std::setw(20) << "Date"
	      << std::endl << std::setw(20) << fEvtNum
	      << std::setw(20) << fRun << std::setw(20) << fDate << std::endl;
  }

};


class CsRectPixelMumegaPlanePar {

 private:

  Int_t fClusMeth; // Clustering method : 0_primitive, 1_full, 2_always-fit
  Int_t fClusConfigMask; // Clustering configuration
  Float_t fThrClus; // Threshold on cluster amplitude
  Float_t fThrHit; // Threshold on single hit amplitude
  Float_t fCrossTalkR; // Maximal ratio of one hit to be marked cross talk (not implemented yet)
  std::vector<Float_t> fTimeCrossTalkRi; // value of min ai/aip1 ratio, used to correct amp of hits suspected of XTalk in time
  std::vector<Float_t> fTimeCrossTalkRip1; // value of min aip1/ai ratio, used to correct amp of hits suspected of XTalk in time
  Int_t fSample; // Sample to use for clustering
  std::vector<Float_t> fTimeCalOffset; // Parameters for cluster time computation
  std::vector<Float_t> fTimeCalSlope; // Parameters for cluster time computation
  std::vector<Float_t> fTimeCalTimeOffset; // Parameters for cluster time computation

  std::vector<Float_t> fTimeCuts; // Parameters for cluster time cuts
  std::vector<Float_t> fPosCorrs; // Parameters for correction of CoG position
  std::vector<Float_t> fCoorda1a2, fCoorda0a2, fCoordtcs, fCoorda1a2tcs; // Parameters for Cluster Amplitude Ratio Cuts

  Bool_t TimeCal, TxtCal, PosCorrCal; 
  Float_t fMergingAmpratioThreshold; //ratio between hits w/ highest and lowest amplitude in a cluster with (at least) a shared hit, above which the clusters are merged

  Bool_t fDoAmpCuts;

 public:
  
  CsRectPixelMumegaPlanePar() : fClusMeth(1), fClusConfigMask(179), fThrClus(5.), fThrHit(3.), fSample(2),   TimeCal(false),  TxtCal(false), PosCorrCal(false), fMergingAmpratioThreshold(0.35), fDoAmpCuts(false) { }

  virtual ~CsRectPixelMumegaPlanePar() { }

  void Set(Int_t _meth = 1, Int_t _mask = 179, Float_t _thrclus = 5., Float_t _thrhit = 3., Float_t _ctalkr = 0.3, Int_t _sample = 2, Float_t _MergingAmpratioThreshold = 0.35) {
    fClusMeth = _meth;
    fClusConfigMask = _mask;
    fThrClus = _thrclus;
    fThrHit = _thrhit;
    fCrossTalkR = _ctalkr;
    fSample = _sample;
    fMergingAmpratioThreshold = _MergingAmpratioThreshold;
  }
  void SetClusMeth(Int_t _meth = 1) { fClusMeth = _meth; }
  void SetClusConfigMask(Int_t _mask = 179) { fClusConfigMask = _mask; }
  void SetThrClus(Float_t _thrclus = 5.) { fThrClus = _thrclus; }
  void SetThrHit(Float_t _thrhit = 3.) { fThrHit = _thrhit; }
  void SetCrossTalkHitRatio(Float_t _ctalkr = 0.3) { fCrossTalkR = _ctalkr; }
  void SetTimeCrossTalkHitRatioi(const std::vector<Float_t> &_ctalkt) {fTimeCrossTalkRi = _ctalkt; }
  void SetTimeCrossTalkHitRatioip1(const std::vector<Float_t> &_ctalkt) {fTimeCrossTalkRip1 = _ctalkt; } 
  void SetMergingAmpratioThreshold( Float_t _ampratio = 0.35){fMergingAmpratioThreshold = _ampratio;}
  void SetSample(Int_t _sample = 2) { fSample = _sample; }
  bool SetPosCorrs(const std::vector<Float_t> &_poscorrs) {
    fPosCorrs.clear();
    if (_poscorrs.size() == 16) {
      fPosCorrs = _poscorrs;
      return true;
    } else
      return false;
  }
  void SetTimeCal(const std::vector<Float_t> &_time_cal_offset, const std::vector<Float_t> &_time_cal_slope, const std::vector<Float_t> &_time_cal_time_offset) {fTimeCalOffset = _time_cal_offset, fTimeCalSlope = _time_cal_slope, fTimeCalTimeOffset = _time_cal_time_offset;};
  void SetTimeCuts(const std::vector<Float_t> &_time_cuts ) {fTimeCuts = _time_cuts;};
  void SetAmpRatioCutsParams(const std::vector<Float_t> &_coord_a1a2, const std::vector<Float_t> &_coord_a0a2,const std::vector<Float_t> &_coord_tcs, const std::vector<Float_t> &_coord_a1a2tcs ){fCoorda1a2 = _coord_a1a2; fCoorda0a2 = _coord_a0a2, fCoordtcs = _coord_tcs; fCoorda1a2tcs = _coord_a1a2tcs; };
  void TimeCalDefault() {TimeCal = false;}
  void TxtCalDefault() {TxtCal = false;}
  void PosCorrCalDefault() {PosCorrCal = false;}
  void DisableAmpCuts(){fDoAmpCuts = false;}  
  
  void EnableTimeCal() {TimeCal = true;}
  void EnableTxtCal() {TxtCal = true;}
  void EnablePosCorrCal() {PosCorrCal = true;}
  void EnableAmpCuts(){fDoAmpCuts = true;}

  Int_t GetClusMeth() const { return fClusMeth; }
  Int_t GetClusConfigMask() const { return fClusConfigMask; }
  Float_t GetThrClus() const { return fThrClus; }
  Float_t GetThrHit() const { return fThrHit; }
  Float_t GetCrossTalkHitRatio() const { return fCrossTalkR; }
  std::vector<Float_t> GetTimeCrossTalkHitRatioi() const { return fTimeCrossTalkRi; }
  std::vector<Float_t> GetTimeCrossTalkHitRatioip1() const { return fTimeCrossTalkRip1; }
  Int_t GetSample() const { return fSample; }

  const std::vector<Float_t>& GetPosCorrs() const { return fPosCorrs; }

  const std::vector<Float_t>& GetTimeCalOffset() const {return fTimeCalOffset; }
  const std::vector<Float_t>& GetTimeCalSlope() const {return fTimeCalSlope; }
  const std::vector<Float_t>& GetTimeCalTimeOffset() const {return fTimeCalTimeOffset; }

  const std::vector<Float_t>& GetTimeCuts() const {return fTimeCuts; }
  Float_t GetMergingAmpratioThreshold() {return fMergingAmpratioThreshold; }
  const std::vector<Float_t>& GetCoorda1a2() const {return fCoorda1a2;}
  const std::vector<Float_t>& GetCoorda0a2() const {return fCoorda0a2;}
  const std::vector<Float_t>& GetCoordtcs() const {return fCoordtcs;}
  const std::vector<Float_t>& GetCoorda1a2tcs() const {return fCoorda1a2tcs;}
  
  Bool_t IsTimeCal(){return TimeCal;}
  Bool_t IsTxtCal(){return TxtCal;}
  Bool_t IsPosCorrCal(){return PosCorrCal;}

  Bool_t DoAmpCuts(){return fDoAmpCuts;}
};


class CsRectPixelMumegaPlane {

 protected:

  CsRectPixelMumegaPlaneHeader fHeader; // plane header
  CsRectPixelMumegaPlanePar fPar; // plane parameters
  CsMumegaTimeCals fTimeCals;
  std::map<CsMumegaChanId, CsMumegaChan*> fChannels; // map of all channels including calibrations and active neighbours
  std::list<CsMumegaHit*> fHits; // list with pointers to hits
  std::list<CsRectPixelMumegaCluster*> fClusters; // list with pointers to clusters
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

  CsRectPixelMumegaPlane();
  CsRectPixelMumegaPlane(std::string name);

  virtual ~CsRectPixelMumegaPlane();

  virtual void Clear(Option_t *option = ""); // Clear all entries of plane
  void ClearChannels(Option_t* _option = "");
  void ClearHits(Option_t* _option = "");
  void ClearClusters(Option_t* _option = ""); // Clear all clusters
  void SetHeader(Int_t i, Int_t run, Int_t date);
  void SetPar(Int_t _meth = 1, Int_t _mask = 179, Float_t _thrclus = 5., Float_t _thrhit = 3., Float_t _ctalkr = 0.3, Int_t _sa = 2) { fPar.Set(_meth, _mask, _thrclus, _thrhit, _ctalkr, _sa); }
  void AddChan(Int_t _chandet, Int_t _pixnb, Float_t _xpos, Float_t _ypos, Int_t _flag, Float_t ped, Float_t _sigma, Int_t _apvchip = -1, Int_t _apvchan = -1, Int_t _connnb = -1);
  void AddChan(CsMumegaChan* _chan) { fChannels[*(_chan->GetId())] = _chan; }
  void SetTimeCals(const CsMumegaTimeCals& _timecals) { fTimeCals = _timecals; }
  void SetTimeCals(const CsMumegaTimeCals* _timecals) { fTimeCals = *_timecals; }
  void AddHit(Int_t _detchan, Int_t _pixnb, Float_t _px, Float_t _py, std::vector<Float_t> _amp);
  void AddHit(Int_t _detchan, Int_t _pixnb, Float_t _px, Float_t _py, Float_t _amp0, Float_t _amp1, Float_t _amp2);
  void Clusterize();
  void SimpleClustering();
  void FullClustering();
  bool IsInside(float xp, float yp, std::vector<float> x, std::vector<float> y);
  void SelectClusters();
  void MergeClusters();

  // Add a cluster to plane, defaut is no flagged pad
  void SortHitAmps(); // sort hits by amplitude
  Bool_t IsHitAmpSorted() const { return fIsHitAmpSorted; }
  Bool_t IsClusterSorted() const { return fIsClusterSorted; }
  Bool_t IsClusterized() const { return fIsClusterized; }
  Int_t GetNchannel() const { return fChannels.size(); }
  Int_t GetNhit() const { return fHits.size(); };
  Int_t GetNcluster() const { return fClusters.size(); };
  std::string GetName()const { return fName; }
  CsRectPixelMumegaPlaneHeader* GetHeader() { return &fHeader; }
  CsRectPixelMumegaPlanePar* GetPar() { return &fPar; }
  const CsMumegaTimeCals* GetTimeCals() const { return &fTimeCals; }
  const std::map<CsMumegaChanId, CsMumegaChan*> &GetChannels() const { return fChannels; }
  const std::list<CsMumegaHit*> &GetHits() const { return fHits; };
  const std::list<CsRectPixelMumegaCluster*> &GetClusters() const { return fClusters; };
  const CsMumegaChan* GetChan(CsMumegaChanId &_id) const {
    return (fChannels.find(_id) != fChannels.end() ?
	    fChannels.find(_id)->second : NULL);
  }
  void Display(Int_t sample = 2, bool wait = true);
  void GetKey(Int_t event, Int_t key, Int_t keysystem, TObject* _selected);
  void PrintChannels(Int_t _sample = 2);
  void PrintChannelsSimple(Int_t _sample = 2);
  void PrintClusters(Int_t _sample = 2);
  virtual void FindXTalk();

  ClassDef(CsRectPixelMumegaPlane, 2)
};

#endif
