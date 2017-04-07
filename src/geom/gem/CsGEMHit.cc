/*
--------------------------------------------------------------------------

 Implementation of classes to store hits from GEM
 detectors

 Author: Bernhard Ketzer   27/11/2002   v0
                           20/05/2008   v1

--------------------------------------------------------------------------
*/

// Class declarations
#include "CsGEMHit.h"

//-----------------------------------------------------------------------------
// CsGEMHit constructor
//-----------------------------------------------------------------------------
CsGEMHit::CsGEMHit(const CsGEMChan* _chan, const std::vector<float>& _amp) :
    fChan(_chan), fAmp(_amp), fHasTime(false), fTime(0.), fTimeErr(1e12), fXTalk(0.), fNrClu(0)
{
}

bool CsGEMHit::GetTime(double &_time, double &_etime) const
{
    _time  = fTime;
    _etime = fTimeErr;

    return fHasTime;
}

void CsGEMHit::SetTime(const double& _time, const double& _etime)
{
    fHasTime = true;
    fTime    = _time;
    fTimeErr = _etime;
}
