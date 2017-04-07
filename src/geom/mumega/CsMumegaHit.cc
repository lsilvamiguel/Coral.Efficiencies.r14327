/*
--------------------------------------------------------------------------

 Implementation of classes to store hits from Mumega
 detectors

 Author : Bernhard Ketzer   05/06/2009    v0

--------------------------------------------------------------------------
*/

// Class declarations
#include "CsMumegaHit.h"

CsMumegaHit::CsMumegaHit(const CsMumegaChan* _chan, std::vector<Float_t> _amp) {
  fChan = _chan;
  fAmp = _amp;
  fHasTime = false;
  fXTalk = 0.;
  fNrClu = 0;
}


bool CsMumegaHit::GetTime(Double_t &_time, Double_t &_etime) const {
  _time = fTime;
  _etime = fTimeErr;

  return fHasTime;
}


void CsMumegaHit::SetTime(Double_t _time, Double_t _etime) {
  fHasTime = true;
  fTime = _time;
  fTimeErr = _etime;
}
