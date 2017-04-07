#ifndef CsGEMHit_H
#define CsGEMHit_H

/*
--------------------------------------------------------------------------

 Declaration of classes to store hits from COMPASS GEM
 detectors

--------------------------------------------------------------------------

class CsGEMHit :

                 Object corresponding to the raw data of
                 a single channel as given by the APV chip in
                 multi (3-sample) mode. Its data members are
                 CsGEMChan* fChan
                 vector<float> fAmp

--------------------------------------------------------------------------
v1.0     19/11/2002    by    Bernhard Ketzer

v2.0     26/05/2008    by    Bernhard Ketzer
--------------------------------------------------------------------------
*/

// GEM headers
#include "CsGEMChan.h"

// C++ headers
#include <vector>

//-----------------------------------------------------------------------------
// CsGEMHit class declaration
//-----------------------------------------------------------------------------
class CsGEMHit {

 private:

    const CsGEMChan*     fChan;      //! detector channel
    std::vector<float>   fAmp;       //  amplitudes of three samples
    bool                 fHasTime;   //  is a valid time set
    double               fTime;      //  time of hit, calculated from three samples
    double               fTimeErr;   //  error on time of hit
    double               fXTalk;     //  cross talk propability
    int                  fNrClu;     //  number of clusters this hit contributes to

 public:

                                CsGEMHit() : fChan(0), fHasTime(false), fTime(0.), fTimeErr(1e12), fXTalk(0.), fNrClu(0) { }
                                CsGEMHit(const CsGEMChan* _chan, const std::vector<float>& _amp);

    virtual                    ~CsGEMHit() { }

    const CsGEMChan*            GetChan()                                           const { return fChan; }
    const std::vector<float>&   GetAmp()                                            const { return fAmp; }
    void                        SetAmp(const std::vector<float>& _amp)                    { fAmp=_amp; }
    void                        SetTime(const double& _time, const double& _etime);
    bool                        GetTime(double &_time, double &_etime)              const;
    void                        SetXTalk(const double& _xtalk)                            { fXTalk = _xtalk; }
    double                      GetXTalk()                                          const { return fXTalk; }
    void                        IncNrClusters()                                           { fNrClu++; }
    int                         GetNrClusters()                                     const { return fNrClu; }

};

#endif

