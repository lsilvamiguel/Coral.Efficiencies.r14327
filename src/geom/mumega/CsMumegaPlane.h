#ifndef CsMumegaPlane_H
#define CsMumegaPlane_H

/*
-------------------------------------------------------------------------

 Declaration of classes to store/analyse events from COMPASS PMumegas
 detectors

-------------------------------------------------------------------------

class MumegaPlane :

                 Object representing a plane of hits in one event
		 The data member fHits is a list of pointers to CsMumegaHit
		 objects
		 The data member fClusters is a list of pointers to
		 CsMumegaCluster objects, each comtaining a list of
		 references to the csMumegaHit objects it is made of.

-------------------------------------------------------------------------

class CsMumegaPlaneHeader :

                 Header object containing event number, run number, date
		 (not used for the time being)

-------------------------------------------------------------------------

CsMumegaPlanePar :

                 Object containing thresholds for hits and clusters,
		 clustering parameters and switches

-------------------------------------------------------------------------
v1.0     29/05/2009     by     Bernhard Ketzer
-------------------------------------------------------------------------
*/

// ROOT headers
#include "TMath.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPolyMarker.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TLine.h"
#include "RQ_OBJECT.h"
#include "KeySymbols.h"

// C++ headers
#include <iostream>
#include <cmath>
#include <list>
#include <vector>
#include <functional>
#include <map>
#include <string>

// Mumega headers
#include "CsMumegaTimeCal.h"
#include "CsMumegaChan.h"
#include "CsMumegaHit.h"
#include "CsMumegaCluster.h"

//-------------------------------------------------------------------------
// CsMumegaPlaneHeader class declaration
//-------------------------------------------------------------------------

class CsMumegaPlaneHeader {

 private:

  Int_t fEvtNum; // Event number
  Int_t fRun; // Run Number;
  Int_t fDate; // Date
  Float_t fTCSPhase; // TCSPhase time

 public:

  CsMumegaPlaneHeader() : fEvtNum(0), fRun(0), fDate(0), fTCSPhase(-100)  { };

  virtual ~CsMumegaPlaneHeader() { };

  void Set(Int_t _i = 0, Int_t _r = 0, Int_t _d = 0, Float_t _t = -100) { fEvtNum = _i; fRun = _r; fDate = _d; fTCSPhase = _t; }

  Int_t GetEvtNum() const { return fEvtNum; }
  Int_t GetRun() const { return fRun; }
  Int_t GetDate() const { return fDate; }
  Float_t GetTCSPhase() const { return fTCSPhase; }

};

//-------------------------------------------------------------------------
// CsMumegaPlanePar class declaration
//-------------------------------------------------------------------------

class CsMumegaPlanePar {

 private:

  Int_t fClusMeth; // Clustering method : 0-primitive, 1-full, 2-always_fit
  Int_t fClusConfigMask; // Clustering configuration
  Int_t fShareHits; // 0-no_sharing, 1-to_all, 2-fraction
  Float_t fThrClus; // Threshold on cluster amplitude
  Float_t fThrHit; // Threshold on single hit amplitude 
  std::vector<Float_t> fTimeCrossTalkRi; // Maximal ratio of one hit to be marked cross talk in time, for each chip
  std::vector<Float_t> fTimeCrossTalkRip1; // Peak of ai/aip1 ratio, used to correct amp of hits suspected of XTalk in time
  Int_t fSample; // sample to use for clustering
  std::vector<Float_t> fTimeCalOffset; // Parameters for cluster time computation
  std::vector<Float_t> fTimeCalSlope; // Parameters for cluster time computation
  std::vector<Float_t> fTimeCalTimeOffset; // Parameters for cluster time computation
  std::vector<Float_t> fTimeCuts;
  std::vector<Float_t> fCoorda1a2, fCoorda0a2, fCoordtcs, fCoorda1a2tcs; // Parameters for Cluster Amplitude Ratio Cuts


  Bool_t TimeCal, TxtCal;
  Float_t fMergingAmpratioThreshold; //ratio between hits w/ highest and lowest amplitude in a cluster with (at least) a shared hit, above which the clusters are merged
  
  Bool_t fDoAmpCuts;

 public:

  CsMumegaPlanePar() : fClusMeth(1), fClusConfigMask(51), fShareHits(2), fThrClus(5.), fThrHit(3.), fSample(2), TimeCal(false), TxtCal(false), fMergingAmpratioThreshold(0.35), fDoAmpCuts(false)   { }

  virtual ~CsMumegaPlanePar() { }


  void Set(Int_t _meth = 1, Int_t _mask = 51, Int_t _share = 2, Float_t _thrclus = 5., Float_t _thrhit = 3., Int_t _sample = 2, Float_t _MergingAmpratioThreshold = 0.35) {
    fClusMeth = _meth; fClusConfigMask = _mask; fShareHits = _share;
    fThrClus = _thrclus; fThrHit = _thrhit; fSample = _sample;
    fMergingAmpratioThreshold = _MergingAmpratioThreshold;
  }

  void SetClusMeth(Int_t _meth = 1) { fClusMeth = _meth; }
  void SetClusConfigMask(Int_t _mask = 51) { fClusConfigMask = _mask; }
  void SetShareHits(Int_t _share = 2) { fShareHits = _share; }
  void SetThrClus(Float_t _thrclus = 5.) { fThrClus = _thrclus; }
  void SetThrHit(Float_t _thrhit = 3.) { fThrHit = _thrhit; }
  void SetTimeCrossTalkHitRatioi(const std::vector<Float_t> &_ctalkt) {fTimeCrossTalkRi = _ctalkt; }
  void SetTimeCrossTalkHitRatioip1(const std::vector<Float_t> &_ctalkt) {fTimeCrossTalkRip1 = _ctalkt; } 
  void SetMergingAmpratioThreshold( Float_t _ampratio = 0.35){fMergingAmpratioThreshold = _ampratio;}
  void SetSample(Int_t _sample = 2) { fSample = _sample; }
  void SetTimeCal(const std::vector<Float_t> &_time_cal_offset, const std::vector<Float_t> &_time_cal_slope, const std::vector<Float_t> &_time_cal_time_offset) {fTimeCalOffset = _time_cal_offset, fTimeCalSlope = _time_cal_slope, fTimeCalTimeOffset = _time_cal_time_offset;};
  void SetTimeCuts(const std::vector<Float_t> &_time_cuts ) {fTimeCuts = _time_cuts;};
  void SetAmpRatioCutsParams(const std::vector<Float_t> &_coord_a1a2, const std::vector<Float_t> &_coord_a0a2,const std::vector<Float_t> &_coord_tcs, const std::vector<Float_t> &_coord_a1a2tcs ){fCoorda1a2 = _coord_a1a2; fCoorda0a2 = _coord_a0a2, fCoordtcs = _coord_tcs; fCoorda1a2tcs = _coord_a1a2tcs; };
  void TimeCalDefault() {TimeCal = false;}
  void TxtCalDefault() {TxtCal = false;}
  void DisableAmpCuts(){fDoAmpCuts = false;}
  
  void EnableTimeCal() {TimeCal = true;}
  void EnableTxtCal() {TxtCal = true;}
  void EnableAmpCuts(){fDoAmpCuts = true;}
  
  Int_t GetClusMeth() const { return fClusMeth; }
  Int_t GetClusConfigMask() const { return fClusConfigMask; }
  Int_t GetShareHits() const { return fShareHits; }
  Float_t GetThrClus() const { return fThrClus; }
  Float_t GetThrHit() const { return fThrHit; }
  std::vector<Float_t> GetTimeCrossTalkHitRatioi() const { return fTimeCrossTalkRi; }
  std::vector<Float_t> GetTimeCrossTalkHitRatioip1() const { return fTimeCrossTalkRip1; }
  Int_t GetSample() const { return fSample; }
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
  Bool_t DoAmpCuts(){return fDoAmpCuts;}
};

//-------------------------------------------------------------------------
// CsMumegaPlane class declaration
//-------------------------------------------------------------------------

class CsMumegaPlane {

 private:

  CsMumegaPlaneHeader fHeader; // plane header
  CsMumegaPlanePar fPar; // plane parameters
  CsMumegaTimeCals fTimeCals; // time calibrations
  std::map<CsMumegaChanId, CsMumegaChan*> fChannels; // map of all channels including calibrations and active neighbours
  std::list<CsMumegaHit*> fHits; // list with pointers to hits
  std::list<CsMumegaCluster*> fClusters; // list with pointers to clusters
  Bool_t fIsHitSorted; // list of hits sorted by channel Id
  Bool_t fIsHitAmpSorted; // list of hits sorted by amplitude
  Bool_t fIsClusterSorted; // cluster sorted by channel Id
  Bool_t fIsClusterized;
  std::string fName; // name of the plane
  struct fCompareHits; //! to compare list of pointers
  struct fCompareHitAmps0; //! to compare list of pointers
  struct fCompareHitAmps1; //! to compare list of pointers
  struct fCompareHitAmps2; //! to compare list of pointers
  struct fCompareClusters; //! to compare list of pointers
  TCanvas* fDisplayCanvas; //! ROOT canvas to display single event
  TH1F* fDisplayHist; //! ROOT histogram to display single event
  Int_t fKey; //! keystroke in fDisplayCanvas


 public:

  CsMumegaPlane();
  CsMumegaPlane(std::string name);

  virtual ~CsMumegaPlane();
  
  void         Clear(Option_t* _option = "");
  int          ClearChannels(Option_t* _option = "");
  void         ClearHits(Option_t* _option = "");
  void         ClearClusters(Option_t* _option = "");
  static void  Reset(Option_t* _option = "");
  void         SetHeader(Int_t i, Int_t run, Int_t date);
  void         SetPar(Int_t _meth = 1, Int_t _mask = 51, Int_t _share = 2, Float_t _thrclus = 5., Float_t _thrhit = 3., Int_t _sa = 2) { fPar.Set(_meth, _mask, _share, _thrclus, _thrhit, _sa); }
  void         AddChan(Int_t _detchan, Int_t _hem, Int_t _flag, Float_t _ped, Float_t _sigma, Int_t _apvchip = -1, Int_t _apvchan = -1, int _stripconn = -1 );
  void         AddChan(CsMumegaChan* _chan) { fChannels[*(_chan->GetId())] = _chan; }
  void         SetTimeCals(const CsMumegaTimeCals& _timecals) { fTimeCals = _timecals; }
  void         SetTimeCals(const CsMumegaTimeCals* _timecals) { fTimeCals = *_timecals; }
  void         AddHit(Int_t _detchan, Int_t _hem, Float_t _amp0, Float_t _amp1, Float_t amp2);
  void         AddHit(Int_t _detchan, Int_t _hem, std::vector<Float_t> _amp);
  void         Clusterize();
  void         SimpleClustering();
  void         FullClustering();
  void         FindXTalk();
  void         MergeClusters();
  bool         IsInside(float xp, float yp, std::vector<float> x, std::vector<float> y);
  void         SelectClusters();
  void         SortHits(); // sort hits by channel Id
  void         SortHitAmps(); // sort hits by amplitude
  void         SortClusters();
  Bool_t       IsHitSorted() const { return fIsHitSorted; }
  Bool_t       IsClusterSorted() const { return fIsClusterSorted; }
  Bool_t       IsClusterized() const { return fIsClusterized; }
  Int_t        GetNchannel() const { return fChannels.size(); }
  Int_t        GetNhit() const { return fHits.size(); }
  Int_t        GetNcluster() const { return fClusters.size(); }
  std::string  GetName() { return fName; }
  CsMumegaPlaneHeader* GetHeader() { return &fHeader; }
  CsMumegaPlanePar* GetPar() { return &fPar; }
  const CsMumegaTimeCals* GetTimeCals() const { return &fTimeCals; }
  const std::map<CsMumegaChanId, CsMumegaChan*> &GetChannels() const { return fChannels; }
  const std::list<CsMumegaHit*> &GetHits() const { return fHits; }
  const std::list<CsMumegaCluster*> &GetClusters() const { return fClusters; }
  const CsMumegaChan* GetChan(CsMumegaChanId _id) const { return (fChannels.find(_id) != fChannels.end() ? fChannels.find(_id)->second : NULL); }
  void         Display(Int_t _sample = 2);
  void         GetKey(Int_t _event, Int_t _key, Int_t _keysym, TObject* _selected);
  void         PrintChannels();
  void         PrintClusters(Int_t _sample = 2);

  ClassDef(CsMumegaPlane, 2)
};

#endif
