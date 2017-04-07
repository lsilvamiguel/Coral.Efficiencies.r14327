#ifndef CsGEMPlanePar_H
#define CsGEMPlanePar_H

/*
--------------------------------------------------------------------------

 Declaration of classes to store/analyze events from COMPASS GEM
 detectors

--------------------------------------------------------------------------

class CsGEMPlanePar

                 Object containing thresholds for hits and clusters,
                 clustering parameters and switches

--------------------------------------------------------------------------
v1.0     19/11/2002    by    Bernhard Ketzer

v2.0     27/05/2008    by    Bernhard Ketzer
--------------------------------------------------------------------------
*/

#include <cmath>

//-----------------------------------------------------------------------------
// CsGEMPlanePar class declaration
//-----------------------------------------------------------------------------
class CsGEMPlanePar {

  private:

    int                 fClusMeth;          // Clustering method: 0: primitive, 1: full, 2: always fit
    int                 fClusConfigMask;    // Clustering configuration
    int                 fShareHits;         // 0 : no sharing; 1 : to all; 2 : fraction
    float               fThrClus;           // Threshold on cluster amplitude
    float               fThrHit;            // Threshold on single hit amplitude
    float               fTimeCrossTalkR;    // Maximal ratio of one hit to be marked cross talk in time
    int                 fSample;            // Sample to use for clustering
    std::vector<float>  fClusSizeRes;       // Size dependent spatial resolution of clusters

  public:

    CsGEMPlanePar() : fClusMeth(1), fClusConfigMask(51),fShareHits(2),fThrClus(5.),fThrHit(3.),fTimeCrossTalkR(.0),fSample(2) {
      fClusSizeRes.push_back(1./sqrt(12.));
    }

    virtual ~CsGEMPlanePar() { }

    virtual void    SetClusMeth(int _meth=1)                                { fClusMeth         = _meth;        }
    virtual void    SetClusConfigMask(int _mask=51)                         { fClusConfigMask   = _mask;        }
    virtual void    SetShareHits(int _share=2)                              { fShareHits        = _share;       }
    virtual void    SetThrClus(float _thrclus=5.)                           { fThrClus          = _thrclus;     }
    virtual void    SetThrHit(float _thrhit=3.)                             { fThrHit           = _thrhit;      }
    virtual void    SetTimeCrossTalkHitRatio(float _ctalkt=.0)              { fTimeCrossTalkR   = _ctalkt;      }
    virtual void    SetSample(int _sample=2)                                { fSample           = _sample;      }
    virtual void    SetClusSizeRes(const std::vector<float>& _clusSizeRes)  { fClusSizeRes      = _clusSizeRes; }

    virtual int                         GetClusMeth()               const { return fClusMeth;       }
    virtual int                         GetClusConfigMask()         const { return fClusConfigMask; }
    virtual int                         GetShareHits()              const { return fShareHits;      }
    virtual float                       GetThrClus()                const { return fThrClus;        }
    virtual float                       GetThrHit()                 const { return fThrHit;         }
    virtual float                       GetTimeCrossTalkHitRatio()  const { return fTimeCrossTalkR; }
    virtual int                         GetSample()                 const { return fSample;         }
    virtual const std::vector<float>&   GetClusSizeRes()            const { return fClusSizeRes;    }
};

#endif

