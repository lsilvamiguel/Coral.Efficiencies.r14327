#ifndef CsMumegaChan_H
#define CsMumegaChan_H

/*
-------------------------------------------------------------------------

Declaration of class to store channel calibrations and neighbouring
active channels for Mumega detectors

-------------------------------------------------------------------------

class CsMumegaChanId :

		Object holding detector channel number and position
		- the detector channel is given by the mapping
		- position = hemisphere for strip-like detectors
		- position = x + (y<<10) for pixel detectors

-------------------------------------------------------------------------

class CsMumegaChanCal :

                Object holding flag, pedestal, sigma for one channel

-------------------------------------------------------------------------

class CsMumegaChan :

                Object holding CsMumegaChanId, CsMumegaChanCal,
		list<CsMumegaChan*> to neighbouring active channels
		(used for clustering)

-------------------------------------------------------------------------

class CsMumegaChannels :

                Object holding map of CsMumegaChan objects for all
		CsMumegaChanId's of a plane

-------------------------------------------------------------------------
v1.0     29/05/2008    by    Bernhard Ketzer
-------------------------------------------------------------------------
*/

// C++ headers
#include <list>

//-------------------------------------------------------------------------
// Class definition CsMumegaChanId
//-------------------------------------------------------------------------
class CsMumegaChanId {

 private:
  int fAPVChip;
  int fAPVChannel;
  int fConnNb;
  int fAPVDetectorChannel;  // after mapping
  int fPosition;          // -1/+1 for up/down or left/right of split strips, (x + (y << 10) for pixels, else 0
  float fPosX;  // to store x position for rectangular pixel MM
  float fPosY;  // to store y position for rectangular pixel MM
  float fPitchX;  // to store x pitch
  float fPitchY;  // to store y pitch
  int fType;  // 0: strip, 1: square pixel, 2: short rect. pixel, 3: long rect. pixel

 public:

  CsMumegaChanId() { CsMumegaChanId(0,0); }

  CsMumegaChanId(int _detchan, int _pos = 0) {
    fAPVDetectorChannel = _detchan; fPosition = _pos;
    fAPVChip = -1; fAPVChannel = -1; fConnNb = -1;
    fPosX = 0; fPosY = 0;
    fType = 0;
    fPitchX = -1; fPitchY = -1;
  }

  CsMumegaChanId(int _detchan, int _px, int _py ) {
    fAPVDetectorChannel = _detchan; fPosition = _px + (_py<<10);
    fAPVChip = -1; fAPVChannel = -1; fConnNb = -1;
    fPosX = 0; fPosY = 0;
    fType = 1;
    fPitchX = 1; fPitchY = 1;
  }

  CsMumegaChanId(int _detchan, int _pixnb, double _px, double _py) {
    fAPVDetectorChannel = _detchan; fPosition = _pixnb;
    fPosX = _px; fPosY = _py;
    fAPVChip = -1; fAPVChannel = -1; fConnNb = -1;
    fPitchX = 0.4;
    fType = 2; fPitchY = 2.5;
    if (_pixnb >= 640) { fType = 3; fPitchY = 6.25; }
  }

  virtual ~CsMumegaChanId() {};

  int     GetDetectorChannel() const { return fAPVDetectorChannel; }
  int     GetPosition()        const { return fPosition; }
  int     GetHemisphere()      const { return fPosition; }
  int     GetType()            const { return fType; }
  float   GetPosX()            const { if (fType > 1) return fPosX; else return (fPosition & 0x3ff); }
  float   GetPosY()            const { if (fType > 1) return fPosY; else return ((fPosition>>10) & 0x3ff); }
  float   GetPitchX()          const { return fPitchX; }
  float   GetPitchY()          const { return fPitchY; }
  void    SetAPVInfo(int _apvchip=-1, int _apvchan=-1) { fAPVChip=_apvchip; fAPVChannel=_apvchan;}
  void    SetRPConnNb(int  _connnb=-1) {fConnNb=_connnb; } //pixels : ConnNb available in digit infos
  void    SetSNConnNb(int  _connnb=-1) {fConnNb=_connnb; } //strips : strip connector available in digit infos (new !)
  void    SetSConnNb(int _chandet=-1, int _hem=-2 ) { //strips : ConnNb computed from channel nb and hem nb
    if(_chandet<256 && _chandet%2 != 0) fConnNb = 0;
    else if(_chandet<255 && _chandet%2 == 0) fConnNb = 9;
    else if(_chandet>255 && _chandet<384 && _hem == -1) fConnNb = 1;
    else if(_chandet>255 && _chandet<384 && _hem == 1) fConnNb = 8;
    else if(_chandet>383 && _chandet<512 && _hem == -1) fConnNb = 2;
    else if(_chandet>383 && _chandet<512 && _hem == 1) fConnNb = 7;
    else if(_chandet>511 && _chandet<640 && _hem == -1) fConnNb = 3;
    else if(_chandet>511 && _chandet<640 && _hem == 1) fConnNb = 6;
    else if(_chandet>640 && _chandet%2 != 0) fConnNb = 4;
    else if(_chandet>639 && _chandet%2 == 0) fConnNb = 5; 
    else fConnNb= -1;
  } 
  int     GetAPVChip()         const { return fAPVChip; };
  int     GetAPVChannel()      const { return fAPVChannel; };
  int     GetConnNb()          const { return fConnNb; };
  bool    operator < (const CsMumegaChanId &_id)
    const {return (fAPVDetectorChannel==_id.GetDetectorChannel()?
		   fPosition<_id.GetPosition():
		   fAPVDetectorChannel<_id.GetDetectorChannel());}
    
};

//-----------------------------------------------------------------------------
// Class definition CsMumegaChanCal
//-----------------------------------------------------------------------------
class CsMumegaChanCal {
    
 private:
    
  int   fFlag;           // flag for channel on/off
  float fPedMean;        // pedestal mean of channel
  float fPedSigma;       // sigma of pedestals
    
 public:
    
  CsMumegaChanCal() { CsMumegaChanCal(0,0.,0.); }
  CsMumegaChanCal(int _flag, float _ped, float _sig) { fFlag=_flag; fPedMean=_ped; fPedSigma=_sig; }
    
  virtual ~CsMumegaChanCal() {};
    
  int   GetFlag()     const { return fFlag; }
  float GetPedMean()  const { return fPedMean; }
  float GetPedSigma() const { return fPedSigma; }

};

//-----------------------------------------------------------------------------
// Class definition CsMumegaChan
//-----------------------------------------------------------------------------
class CsMumegaChan {
  
 private:

  CsMumegaChanId           fId;                // channel Id 
  CsMumegaChanCal          fCal;               // channel calibration
  std::list<CsMumegaChan*> fActiveNeighbours;   // list of active neighbours
    
 public:
    
  CsMumegaChan() {};
  CsMumegaChan(const CsMumegaChanId& _id, const CsMumegaChanCal& _cal) { fId=_id; fCal=_cal; }

  virtual ~CsMumegaChan() {};

  const CsMumegaChanId*           GetId()               const { return &fId; }
  const CsMumegaChanCal*          GetCal()              const { return &fCal; }
  const std::list<CsMumegaChan*>& GetActiveNeighbours() const { return fActiveNeighbours; }
  void                         AddActiveNeighbour(CsMumegaChan* _chan) 
    { fActiveNeighbours.push_back(_chan); }


};

#endif
