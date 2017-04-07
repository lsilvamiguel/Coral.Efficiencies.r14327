#ifndef CsMumegaCluster_H
#define CsMumegaCluster_H

/*
---------------------------------------------------------------------------

 Declaration of classes to store/analyze clusters from COMPASS Mumegas
 detectors

---------------------------------------------------------------------------

class CsMumegaCluster :

                 Object corresponding to a cluster of channels

---------------------------------------------------------------------------
v1.0      03/06/2009    by    Bernhard Ketzer
---------------------------------------------------------------------------
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
//#include "CsMumegaHit.h"

//---------------------------------------------------------------------------
// CsMumegaCluster class declaration
//---------------------------------------------------------------------------

class CsMumegaHit;


class CsMumegaCluster {

 protected:

  std::vector<Float_t> fAmp; // amplitudes of three samples
  Float_t fNoise; // cluster noise
  Float_t fXTalk; // cross talk probability

  Bool_t fHasTime; // is a valid time set
  Double_t fTime; // time of cluster
  Double_t fTimeErr; // error of time

  Bool_t fContainsFlaggedChannels; // 1 if cluster extends over flagged channels
  std::list<CsMumegaHit*> fHits; // array with hit objects belonging to cluster  Bool_t fIsHitSorted; // list of hits sorted by channel Id
  Bool_t fIsHitAmpSorted; // list of hits sorted by amplitude
  struct fCompareHits; //! to compare list of pointers
  struct fCompareHitAmps0; //! to compare list of pointers
  struct fCompareHitAmps1; //! to compare list of pointers
  struct fCompareHitAmps2; //! to compare list of pointers

  Bool_t fNewClus; // 1 if cluster was created from old clusters that has been merged
  Bool_t fOldClus; // 1 if cluster has been merged into a new cluster

  Int_t fNSharedHits; //number of hits shared by several clusters

 private:

  // strip specific properties
  Float_t fPos; // position of cluster in plane
  Float_t fPosErr; // error of cluster position in plane
  Float_t fHemisphere; // average hemisphere


 public:

  CsMumegaCluster();
  CsMumegaCluster(const std::list<CsMumegaHit*> &_hits, Int_t _NSharedHits, Bool_t _flags,
		  int _sample = 2, float _thrhit = 3., int _share = 2);

  virtual ~CsMumegaCluster() { }

  virtual std::vector<Float_t> GetAmp() const { return fAmp; }
  virtual Float_t GetNoise() const { return fNoise; }
  virtual Float_t GetXTalk() const { return fXTalk; }
  virtual int GetSize() const { return fHits.size(); }
  virtual double GetTime() const { return fTime; }
  virtual bool GetTime(double &_time, double &_etime) const;
  virtual Bool_t ContainsFlaggedChannels() const { return fContainsFlaggedChannels; }
  virtual void CalcAmps(int _sample = 2, float _thr = 3., int _share = 2);
  virtual void CalcTime(const std::vector<Float_t> &_timecaloffset, const std::vector<Float_t> &_timecalslope, const std::vector<Float_t> &_timecaltimeoffset);
  virtual const std::list<CsMumegaHit*> &GetHits() { return fHits; }
  virtual void SetTime(Float_t _time, Float_t _time_err) {fTime = _time, fTimeErr = _time_err;}

  // strip specific methods
  virtual void CalcCoG(int _sample = 2, float _thr = 3., int _share = 2);
  virtual bool AddHit(CsMumegaHit* _hit, Int_t _clusconfigmask = 51, Float_t _thr = 3., Int_t _sample = 2);
  Float_t GetPosition() const { return fPos; }
  Float_t GetPositionErr() const { return fPosErr; }
  Float_t GetHemisphere() const { return fHemisphere; }

  void NewClus(){fNewClus = true;}; 
  void OldClus(){fOldClus = true;}; 

  Bool_t IsNewClus(){return fNewClus;};
  Bool_t IsOldClus(){return fOldClus;};

  void IncrNSharedHits(){fNSharedHits++;};
  void DeacrNSharedHits(){fNSharedHits--;};

  Int_t GetNSharedHits(){return fNSharedHits;};

  void SortHitAmps(int _sample = 2); // sort hits by amplitude

  Bool_t operator< (const CsMumegaCluster &cluster) const;

};

#endif
