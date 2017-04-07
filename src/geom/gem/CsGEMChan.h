#ifndef CsGEMChan_H
#define CsGEMChan_H

/*
--------------------------------------------------------------------------

 Declaration of class to store channel calibrations and neighbouring
 active channels for GEM detectors

--------------------------------------------------------------------------

class CsGEMChanId :

                 Object holding detector channel number and position
                 - the detector channel is given by the mapping
                 - position = hemisphere for strip-like detectors
                 - position = x + (y<<10) for pixel detectors

--------------------------------------------------------------------------
class CsGEMChanCal :

                 Object holding flag, pedestal, sigma for one channel

--------------------------------------------------------------------------

class CsGEMChan :

                 Object holding CsGEMChanId, CsGEMChanCal, list<CsGEMChan*>
                 to neighbouring active channels
                 (used for clustering)

--------------------------------------------------------------------------
class CsGEMChannels :

                 Object holding map of CsGEMChan objects for all
                 CsGEMChanId's of a plane

--------------------------------------------------------------------------
v1.0     26/05/2008    by    Bernhard Ketzer
--------------------------------------------------------------------------
*/

// C++ headers
#include <list>

//-----------------------------------------------------------------------------
// Class definition CsGEMChanId
//-----------------------------------------------------------------------------
class CsGEMChanId {

 private:
    int     fDetectorChannel;  // after mapping
    int     fPosition;       // -1/+1 for up/down or left/right of split strips, x+(y<<10) for pixels, else 0
    int     fAPVChip;
    int     fAPVChannel;

 public:

    CsGEMChanId() :
        fDetectorChannel(0), fPosition(0), fAPVChip(-1), fAPVChannel(-1) { }
    CsGEMChanId(int _detchan, int _pos=0 ) :
        fDetectorChannel(_detchan), fPosition(_pos), fAPVChip(-1), fAPVChannel(-1) { }
    CsGEMChanId(int _detchan, int _px, int _py ) :
        fDetectorChannel(_detchan), fPosition(_px+(_py<<10)), fAPVChip(-1), fAPVChannel(-1) { }

    virtual ~CsGEMChanId() {}

    int     GetDetectorChannel() const { return fDetectorChannel; }
    int     GetPosition()        const { return fPosition; }
    int     GetHemisphere()      const { return fPosition; }
    int     GetPosX()            const { return (fPosition & 0x3ff); }
    int     GetPosY()            const { return ((fPosition>>10) & 0x3ff); }
    void    SetAPVInfo(int _apvchip=-1, int _apvchan=-1) { fAPVChip=_apvchip; fAPVChannel=_apvchan; }
    int     GetAPVChip()         const { return fAPVChip; };
    int     GetAPVChannel()      const { return fAPVChannel; };
    bool    operator < (const CsGEMChanId &_id)
        const {return (fDetectorChannel==_id.GetDetectorChannel()?
                       fPosition<_id.GetPosition():
                       fDetectorChannel<_id.GetDetectorChannel());}

};

//-----------------------------------------------------------------------------
// Class definition CsGEMChanCal
//-----------------------------------------------------------------------------
class CsGEMChanCal {

 private:

    int   fFlag;           // flag for channel on/off
    float fPedMean;        // pedestal mean of channel
    float fPedSigma;       // sigma of pedestals

 public:

    CsGEMChanCal() :
        fFlag(0), fPedMean(0.), fPedSigma(0.) { }
    CsGEMChanCal(int _flag, float _ped, float _sig) :
        fFlag(_flag), fPedMean(_ped), fPedSigma(_sig) { }

    virtual ~CsGEMChanCal() {}

    int   GetFlag()     const { return fFlag; }
    float GetPedMean()  const { return fPedMean; }
    float GetPedSigma() const { return fPedSigma; }

};

//-----------------------------------------------------------------------------
// Class definition CsGEMChan
//-----------------------------------------------------------------------------
class CsGEMChan {

 private:

    CsGEMChanId           fId;                // channel Id
    CsGEMChanCal          fCal;               // channel calibration
    std::list<CsGEMChan*> fActiveNeighbours;  // list of active neighbours

 public:

    CsGEMChan() {};
    CsGEMChan(const CsGEMChanId& _id, const CsGEMChanCal& _cal) { fId=_id; fCal=_cal; }

    virtual ~CsGEMChan() {};

    const CsGEMChanId*           GetId()               const { return &fId; }
    const CsGEMChanCal*          GetCal()              const { return &fCal; }
    const std::list<CsGEMChan*>& GetActiveNeighbours() const { return fActiveNeighbours; }
    void                         AddActiveNeighbour(CsGEMChan* _chan)
        { fActiveNeighbours.push_back(_chan); }

};

#endif

