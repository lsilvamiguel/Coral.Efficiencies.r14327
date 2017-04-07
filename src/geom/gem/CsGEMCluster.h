#ifndef CsGEMCluster_H
#define CsGEMCluster_H

/*
--------------------------------------------------------------------------

 Declaration of classes to store/analyze clusters from COMPASS GEM
 detectors

--------------------------------------------------------------------------

class CsGEMCluster :

                 Object corresponding to a cluster of channels.

--------------------------------------------------------------------------
 v1.0     19/11/2002    by    Bernhard Ketzer

 v2.0     26/05/2008    by    Bernhard Ketzer
--------------------------------------------------------------------------
*/

// GEM headers
#include "CsGEMHit.h"
#include "CsGEMPlanePar.h"

// C++ headers
#include <list>
#include <vector>

//-----------------------------------------------------------------------------
// CsGEMCluster class declaration
//-----------------------------------------------------------------------------
class CsGEMCluster {

  protected:
    std::vector<float>   fAmp;                     // amplitudes of three samples
    float                fNoise;                   // cluster noise
    float                fXTalk;                   // cross talk probablility

    bool                 fHasTime;                 // is a valid time set
    double               fTime;                    // time of cluster
    double               fTimeErr;                 // error of time

    bool                 fContainsFlaggedChannels; // 1 if cluster extends over flagged channels
    std::list<CsGEMHit*> fHits;                    // array with hit objects belonging to cluster

    const CsGEMPlanePar* fClusterParams;           // parameters for clusterization

  private:
    // strip specific properties
    float                fPos;                     // position of cluster in plane
    float                fPosErr;                  // error of cluster position in plane
    float                fHemisphere;              // average hemisphere

  public:
                                CsGEMCluster() : // default constructor to make ROOT happy
                                    fNoise(0.), fXTalk(0.), fHasTime(false), fTime(0.), fTimeErr(1e12),
                                    fContainsFlaggedChannels(false), fClusterParams(0), fPos(-1.), fPosErr(1e12),
                                    fHemisphere(0) { }
                                CsGEMCluster(const CsGEMPlanePar* _par) :
                                    fNoise(0.), fXTalk(0.), fHasTime(false), fTime(0.), fTimeErr(1e12),
                                    fContainsFlaggedChannels(false), fClusterParams(_par), fPos(-1.), fPosErr(1e12),
                                    fHemisphere(0) { }

                                CsGEMCluster(const std::list<CsGEMHit*>& _hits, bool _flags,
                                             const CsGEMPlanePar* _par);

    virtual                    ~CsGEMCluster() { }

    virtual const std::vector<float>&   GetAmp()                               const { return fAmp; }
    virtual float                       GetNoise()                             const { return fNoise; }
    virtual float                       GetXTalk()                             const { return fXTalk; }
    virtual int                         GetSize()                              const { return fHits.size(); }
    virtual bool                        GetTime(double &_time, double &_etime) const;
    virtual bool                        ContainsFlaggedChannels()              const { return fContainsFlaggedChannels; }
    virtual void                        CalcAmps();
    virtual void                        CalcTime();
    virtual const std::list<CsGEMHit*>& GetHits()                              const { return fHits; }

    // strip specific methods
    virtual void                        CalcCoG();
    virtual bool                        AddHit(CsGEMHit* _hit);
    float                               GetPosition()                          const { return fPos; }
    float                               GetPositionErr()                       const { return fPosErr; }
    float                               GetHemisphere()                        const { return fHemisphere; }

    bool operator< (const CsGEMCluster &cluster)                     const;

};

#endif

