/*
---------------------------------------------------------------------------

 Declaration of classes to store/analyze events from COMPASS rectangular PixelMumega
 detectors

 Author : Bernhard Ketzer       10/06/2009
 Modifed by : Damien Neyret     1/02/2011

---------------------------------------------------------------------------

class CsRectPixelMumegaCluster :

                 Object corresponding to a cluster of rectangular pixels

---------------------------------------------------------------------------
*/

#ifndef CsRectPixelMumegaCluster__H
#define CsRectPixelMumegaCluster__H

#include "CsMumegaCluster.h"
#include "CsMumegaHit.h"

class CsRectPixelMumegaPlanePar;

class CsRectPixelMumegaCluster : public CsMumegaCluster {

 private:

  Float_t fPosX; // Position of cluster in plane
  Float_t fPosY;
  Float_t fPosXErr; // Error of cluster position in plane
  Float_t fPosYErr;
  Float_t fEtaX; // Eta variables for cluster position correction
  Float_t fEtaY;
  Int_t fSizeX;
  Int_t fSizeY;

 public:

  CsRectPixelMumegaCluster();
  CsRectPixelMumegaCluster(const std::list<CsMumegaHit*> &_hits, Int_t _NSharedHits, Bool_t _flags, const CsRectPixelMumegaPlanePar &_par);
  

  virtual ~CsRectPixelMumegaCluster() { }

  Float_t GetPositionX() const { return fPosX; }
  Float_t GetPositionY() const { return fPosY; }
  Float_t GetPositionXErr() const { return fPosXErr; }
  Float_t GetPositionYErr() const { return fPosYErr; }
  virtual bool AddHit(CsMumegaHit* _hit, Int_t _clusconfigmask = 179, Float_t _thr = 3., Int_t _sample = 2);

  virtual void CalcCoG(int _sample, float _thr, int _share);
  void CorrectPos(const std::vector<Float_t> &poscorrs);
  Bool_t operator< (const CsRectPixelMumegaCluster &cluster) const;
};

#endif
