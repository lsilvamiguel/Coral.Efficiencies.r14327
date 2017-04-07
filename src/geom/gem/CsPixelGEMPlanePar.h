#ifndef CsPixelGEMPlanePar_H
#define CsPixelGEMPlanePar_H

// GEM headers
#include "CsGEMPlanePar.h"

// C++ headers
#include <vector>

//-----------------------------------------------------------------------------
// CsGEMPlanePar class declaration
//-----------------------------------------------------------------------------
class CsPixelGEMPlanePar : public CsGEMPlanePar {

  private:

    float               fCrossTalkR;        // Maximal ratio of one hit to be marked cross talk

    std::vector<float>  fPosCorrs;

  public:

    CsPixelGEMPlanePar() : CsGEMPlanePar(), fCrossTalkR(0.3) {
      SetClusConfigMask(179);
    }

    virtual ~CsPixelGEMPlanePar() { }

    virtual void    SetClusConfigMask(int _mask=179)             { CsGEMPlanePar::SetClusConfigMask(_mask); }
    virtual void    SetCrossTalkHitRatio(float _ctalkr=0.3)      { fCrossTalkR=_ctalkr; }
    virtual bool    SetPosCorrs(const std::vector<float> &_poscorrs) {
      fPosCorrs.clear();
      if (_poscorrs.size()==48) {
        fPosCorrs=_poscorrs;
        return true;
      } else
        return false;
    }

    virtual float                     GetCrossTalkHitRatio() const { return fCrossTalkR; }
    virtual const std::vector<float>& GetPosCorrs()          const { return fPosCorrs; }

};

#endif

