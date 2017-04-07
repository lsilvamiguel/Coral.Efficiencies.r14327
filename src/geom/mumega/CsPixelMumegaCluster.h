/*
---------------------------------------------------------------------------

 Declaration of classes to store/analyze events from COMPASS PixelMumega
 detectors

 Author : Bernhard Ketzer       10/06/2009

---------------------------------------------------------------------------

class CsPixelMumegaCluster :

                 Object corresponding to a cluster of pads.

---------------------------------------------------------------------------
*/

#ifndef CsPixelMumegaCluster__H
#define CsPixelMumegaCluster__H

#include "CsMumegaCluster.h"
#include "CsMumegaHit.h"

class CsPixelMumegaPlanePar;

class CsPixelMumegaCluster : public CsMumegaCluster {

 private:

  Float_t fPosX; // Position of cluster in plane
  Float_t fPosY;
  Float_t fPosXErr; // Error of cluster position in plane
  Float_t fPosYErr;

  int Categorize();


 public:

  CsPixelMumegaCluster();
  CsPixelMumegaCluster(const std::list<CsMumegaHit*> &_hits, Bool_t _flags, const CsPixelMumegaPlanePar &_par);

  virtual ~CsPixelMumegaCluster() { }

  Float_t GetPositionX() const { return fPosX; }
  Float_t GetPositionY() const { return fPosY; }
  Float_t GetPositionXErr() const { return fPosXErr; }
  Float_t GetPositionYErr() const { return fPosYErr; }
  virtual bool AddHit(CsMumegaHit* _hit, Int_t _clusconfigmask = 179, Float_t _thr = 3., Int_t _sample = 2);
  virtual void CalcCoG(int _sample, float _thr, int _share);
  virtual void CorrectPos(const std::vector<float> &_corrs);

  Bool_t operator< (const CsPixelMumegaCluster &cluster) const;
};

#endif
