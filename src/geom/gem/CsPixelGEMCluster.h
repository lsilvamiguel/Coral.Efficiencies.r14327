/*
--------------------------------------------------------------------------

 Declaration of classes to store/analyze events from COMPASS PixelGEM/Sil
 detectors

 Author: Markus Kr�mer / Thiemo Nagel     March 2007

--------------------------------------------------------------------------

class CsPixelGEMCluster :

        Object corresponding to a cluster of pads.

--------------------------------------------------------------------------
*/

#ifndef CsPixelGEMCluster_H
#define CsPixelGEMCluster_H

#include "CsGEMCluster.h"
#include "CsPixelGEMPlanePar.h"

class CsPixelGEMPlanePar;

//-----------------------------------------------------------------------------
// CsPixelGEMCluster class declaration
//-----------------------------------------------------------------------------
class CsPixelGEMCluster : public CsGEMCluster {

  private:

    float              fPosX;                    // position of cluster in plane
    float              fPosY;
    float              fPosXErr;                 // error of cluster position in plane
    float              fPosYErr;

    int                Categorize();

  public:
                  CsPixelGEMCluster() : // default constructor to make ROOT happy
                    CsGEMCluster(), fPosX(-1.), fPosY(-1.), fPosXErr(1e12), fPosYErr(1e12) { }
                  CsPixelGEMCluster(const CsPixelGEMPlanePar* _par) :
                    CsGEMCluster(_par), fPosX(-1.), fPosY(-1.), fPosXErr(1e12), fPosYErr(1e12) { }
                  CsPixelGEMCluster(const std::list<CsGEMHit*>& _hits, bool _flags,
                                    const CsPixelGEMPlanePar* _par);

    virtual      ~CsPixelGEMCluster() { }

    float         GetPositionX()                         const { return fPosX; }
    float         GetPositionY()                         const { return fPosY; }
    float         GetPositionXErr()                      const { return fPosXErr; }
    float         GetPositionYErr()                      const { return fPosYErr; }
    virtual bool  AddHit(CsGEMHit* _hit);
    virtual void  CalcCoG();
    virtual void  CorrectPos();

    bool operator< (const CsPixelGEMCluster &cluster)                const;

};

#endif // CsPixelGEMCluster_H

