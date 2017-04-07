#ifndef CsMumegaHit_H
#define CsMumegaHit_H

/*
-------------------------------------------------------------------------

 Declaration of classes to store hits from COMPASS Mumega
 detectors

-------------------------------------------------------------------------

class CsMumegaHit :

                 Object corresponding to the raw data of
		 a single channel as given by the APV chip in
		 multi (3-sample) mode. Its data members are
		 CsMumegaChane* fChan
		 vector<float> fAmp

-------------------------------------------------------------------------
v1.0       04/06/2009      by     Bernhard Ketzer
-------------------------------------------------------------------------
*/

// ROOT headers
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"

// C++ headers
#include <iostream>
#include <cmath>
#include <list>
#include <vector>
#include <functional>
#include <string>

// Mumega headers
#include "CsMumegaChan.h"
#include "CsMumegaCluster.h"

class CsMumegaHit {

 private:

  const CsMumegaChan* fChan; // detector channel
  std::vector<Float_t> fAmp; // amplitudes of three samples
  bool fHasTime; // is a valid time set
  Double_t fTime; // time of hit, calculated from three samples
  Double_t fTimeErr; // error on time of hit
  Double_t fXTalk; // cross talk probability
  Int_t fNrClu; // number of clusters in which the hit is included
  Double_t fMultiplexAmpRatio; // ratio ai / ai+1 in multiplexed signal
  std::list<CsMumegaCluster*> fClusters; //clusters the hit belongs to

 public:

  CsMumegaHit() { }
  CsMumegaHit(const CsMumegaChan* _chan, std::vector<Float_t> _amp);

  virtual ~CsMumegaHit() { }

  const CsMumegaChan* GetChan() const { return fChan; }
  std::vector<Float_t> GetAmp() const { return fAmp; }
  void SetAmp(std::vector<Float_t> _amp) { fAmp = _amp; }
  void SetTime(Double_t _time, Double_t _etime);
  bool GetTime(Double_t &_time, Double_t &_etime) const;
  void SetXTalk(Double_t _xtalk) { fXTalk = _xtalk; }
  Double_t GetXTalk() const { return fXTalk; }
  void IncNrClusters() { fNrClu++; }
  void DeacNrClusters() { fNrClu--; }
  Int_t GetNrClusters() const { return fNrClu; }
  void SetMultiplexAmpRatio(Float_t _ampratio = -1) { fMultiplexAmpRatio = _ampratio; }
  Float_t GetMultiplexAmpRatio() const { return fMultiplexAmpRatio; }
  std::list<CsMumegaCluster*> GetClusters(){ return fClusters;}
  void AddCluster(CsMumegaCluster* _clus) { fClusters.push_back(_clus);}
};

#endif
