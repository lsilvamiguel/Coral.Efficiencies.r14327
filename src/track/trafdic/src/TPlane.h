// $Id: TPlane.h 14069 2015-09-17 20:44:46Z lsilva $

#ifndef TPlane_h
#define TPlane_h

/*!
  \class TPlane
  \brief Logical tracking detector's layer

  Every object of this class corresponds to 
  the object (or objects) of the class TDetect.
  (but not other way around).
  It contains dynamic (event-dependent)
  part of information: references to measurements.

  \author Sergei.Gerassimov@cern.ch

*/

/*
  Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/TPlane":
   i) Added member data "iPlane", "associate", "station".
  ii) "vecPixRefs", i.e. hits list ordered along MRS X,Y,U,V corrdinates.
 iii) PixelGEM/MM flagged.
 */

#include <iostream>
#include <CsSTD.h>
#include <cassert>
#include "TStation.h"

class THlx;

class TPlane {


public:
  static unsigned int NobjCreated;
  static unsigned int NobjDestructed;

  const int &IDetRef;  //!< Index of detector plane in TSetup::vDetect
  const int &IProj;    //!< Projection index
  const int &IGr;      //!< Detector's group number

  // "IFlag": It's basically the attribute flagging detectors turned off
  // software-wise (via the "DetNameOff" option). Its meaning has been extended
  // to also serve other purposes: single out pixel pixel detectors (in view of
  // accounting for their 2D nature: # of associated DFs or proj.) and flag
  // detectors to be excluded from the final fit (via the option "Det2Go2Fit").
  const int &IFlag;    //!< Flag to be used in PR: 0=OFF, &0x2=Exclude from final fit, &0x10=XYpixelGP, &0x20=UVpixelGP, &0x40=Y|UpixelMP, &0x80=X|VpixelMP

  const int &IPlane;   //!< Index of detector plane in TSetup::vPlane
  const TPlane *&Associate; //!< Pointer to associate
  const TStation *&Station; //!< Pointer to station

  //! Vector of references to TEv::vecHit
  const std::vector<int>& vHitRef()   const { return vecHitRef; }
  //! Vector of references to TEv::vecHit special for pixelGEMs, ordered along a MRS coordinate (X, Y, U or V, depending upon argument <coord>)
  const std::vector<int>& vPixRefs(int coord)  const {
    if (coord<0 || 4<=coord) {
      printf("TPlane[%d]::vPixRefs(int coord = %d): coord<0 || 4<=coord\n",
	     iPlane,coord); exit(1);
    }
    return vecPixRefs[coord];
  }
  //! Vector of references to TEv::vecHitMC (MC only)
  const std::vector<int>& vHitMCRef() const { return vecHitMCRef; }
  //! Add index (in TEv::vecHit) to vecHitRef in U-coordinate increasing order
  void addHitRef (int ih);
  //! Fill vecPixRefs[coord] w/ indices (in TEv::vecHit), ordering them along the cordinate specified by <c> and <s> in MRS (special for pixelGEMs)
  void fillPixRefs (int coord, double c, double s);
  //! Add index (in vecHitMC of the TEv object) of MC Hit to vecHitMCRef in coordinate increasing order
  void addHitMCRef(int ih);
  //! Reset plane->hit reference vectors
  void Reset();

  //! Constructor
  TPlane(bool hasPixels);
  //! Copy Constructor
  TPlane(const TPlane& p);

  // Automatic destructor 

  //! Operator =
  TPlane& operator = (const TPlane& p);


  //! Quick search of hits around track impact on the plane 
  bool NextHit(THlx& H, double ChiCut, int& iHit, double& Chi2, double& dist, int ifl = 0) const;


private:
  
  std::vector<int> vecHitRef;    //!< Vector of references to TEv::vecHit
  std::vector<int> vecHitMCRef;  //!< Vector of references to TEv::vecHitMC (MC only)
  bool pixelDet;
  std::vector<int> vecPixRefs[4];//!< Array of vectors of references to TEv::vecHit in MRS X,Y,U,V-coordinate increasing order (special for pixelGEMs)

  int iDetRef;
  int iProj;
  int iGr;
  int iFlag;
  int iPlane;
  const TPlane   *associate;
  const TStation *station;

  friend class TSetup; // for TSetup::Init()

};

//! Constructor
inline TPlane::TPlane(bool hasPixel = false):
  // Init const references
  IDetRef(iDetRef),
  IProj  (iProj),
  IGr    (iGr),
  IFlag  (iFlag),
  IPlane (iPlane),
  Associate(associate),
  Station(station),
  // Init data members
  vecHitRef(), 
  vecHitMCRef(),
  pixelDet (hasPixel),
  iDetRef(-1),
  iProj  (-1),
  iGr    (-1),
  iFlag  ( 1),
  iPlane (-1),
  associate(NULL),
  station(NULL)
{
  vecHitRef.reserve(100);
  if (hasPixel) {
    for (int i = 0; i<4; i++) {
      vecPixRefs[i].reserve(100);
    }
  }
  vecHitMCRef.reserve(100);
};



//! Copy constructor
inline TPlane::TPlane(const TPlane& p):
  // Init const references
  IDetRef(iDetRef),
  IProj  (iProj),
  IGr    (iGr),
  IFlag  (iFlag),
  IPlane (iPlane),
  Associate(associate),
  Station(station),
  // Init data members
  pixelDet (p.pixelDet),
  iDetRef(p.iDetRef),
  iProj  (p.iProj),
  iFlag  (p.iFlag),
  iPlane (p.iPlane),
  associate(p.associate),
  station(p.station)
{
  vecHitRef  = p.vecHitRef;
  // Vectors of hit references: Copy contents
  if (pixelDet) {
    for (int i = 0; i<4; i++) vecPixRefs[i] = p.vecPixRefs[i];
  }
  vecHitMCRef= p.vecHitMCRef;
};

//! operator "="
inline TPlane& TPlane::operator=(const TPlane& p)
{
  iDetRef  = p.iDetRef;
  iProj    = p.iProj;
  iGr      = p.iGr;
  iFlag    = p.iFlag;
  iPlane   = p.iPlane;
  associate= p.associate;
  station  = p.station;
  // Vectors of hit references: Copy vector
  vecHitRef  = p.vecHitRef;
  pixelDet = p.pixelDet;
  for (int i = 0; i<4; i++) vecPixRefs[i] = p.vecPixRefs[i];
  vecHitMCRef= p.vecHitMCRef;
  return(*this);
};

#endif
