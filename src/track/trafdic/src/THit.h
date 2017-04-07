// $Id: THit.h 14069 2015-09-17 20:44:46Z lsilva $

/*!
  \class THit
  \brief Cluster

  Contains tracking detector measurment (cluster),
  prepared to simplify it's use in track finding and in Kalman fit.

  \author Sergei.Gerassimov@cern.ch
*/

/*
  Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/THit":
  i) Added method "setKine".
 ii) Added member "status" and "mirror".
 */

#ifndef THit_h
#define THit_h

#include <CsSTD.h>
#include <cassert>
#include <set>

#include "TDigit.h"

class TDetect;
class CsCluster;

class THit {

public:
  static unsigned int NobjCreated;
  static unsigned int NobjDestructed;

  //! Constructor without parameters
  THit();
  //! Copy constructor
  THit(const THit&);
  //! assignment operator
  THit& operator = (const THit&);

  // init by corresponding private data members in the constructor

  const int& IPlane; //!< plane number in vPlane

  /*! 
    original KINE track index in vKine (for MC only). 
    -1 for cluster produced by few MC tracks; 
    -2 for mirror hits in drift detectors;
  */
  const int& IKine;

  /*! 
    2nd KINE track index in vKine (for MC only). 
    -1 when signgle-track
  */
  const int& IKin2;

  /*!
    used to flag hits from track products (for MC only).
    if != 0 - cluster from showers, delta-ray etc.
    IKine points to original track
  */
  const int& IOrig;

  const int& IHit;   //!< hit position in vHit vector

  const int &Status; //!< Status (0: Free,>0: Used,<0: Bypassed)  
  const int &Mirror; //!< Index of mirror

  const double& U;    //!< measurement (WRS)
  const double& V;    //!< measurement (WRS)        - pad detectors only
  const double& SigU; //!< measurement error (WRS)
  const double& SigV; //!< measurement error (WRS)  - pad detectors only
  const double& SigT; //!< time measurement error

  const double& Y;    //!< measurement in MRS
  const double& Z;    //!< measurement in MRS
  const double& G0;   //!< weight matrix (MRS) W11
  const double& G1;   //!< weight matrix (MRS) W21
  const double& G2;   //!< weight matrix (MRS) W22

  const double& Time;        //!< cluster time (mean time of all digits in cluster)
  const double& DeltaR;      //!< distance from the wire, calculated from the drift time. = 0 for not drift detectors 
  const short int& DRsignMC; //!< DeltaR sign, known from MC ( = 0 for real data) 
  
  void Print(int ifl) const;   //! Print cluster info.
  
  //! references to corresponding digits
  const std::set<TDigit>& sDigits() const {return setDigits;}

  //! references to reconstructed tracks (unique track IDs are stored)
  const std::set<unsigned int>& sTrackID() const {return setTrackID;}

  //! add references to reconstructed tracks
  void addTrackID(unsigned int id)   {setTrackID.insert(id);}

  //! erase references to reconstructed track (release the hit)
  void eraseTrackID(unsigned int id) { setTrackID.erase(id);}

  //! get reference to detector
  const TDetect& DetRef()    const;

  //!  get pointer to corresponding CsCluster object
  CsCluster* PtrClus()     const { return ptrCl;}

  /*!
    get pointer to corresponding mirror CsCluster object for drift detectors.
    Returnt NULL for non drift detectors.
  */

  CsCluster* PtrClusMirr() const { return ptrClMirr;}

  //! (Re)set Kine reference
  void setKine(int ikine);

  //! (Re)set sigma time
  void setSigT(double sig) { sigt = sig; };

  //! (Re)set Kine reference
  void setStatus(int i) { status = i; }

  void Rotate();           //!< rotate measurement to MRS (fill y,z g0,g1,g2)
  void UseDrift(int flag); //!< use drift info. (for drift detectors' hits only. Otherway, do nothing)

private:

  CsCluster* ptrCl;
  CsCluster* ptrClMirr;
  int iPlane;
  int iKine;
  int iKin2;
  int iOrig;
  int iHit;
  int status;
  int mirror;
  double u;
  double v;
  double sigu; 
  double sigv;
  double sigt; 
  double y,z;
  double g0,g1,g2;
  std::set<unsigned int> setTrackID;
  std::set<TDigit>       setDigits;
  double time;
  double deltaR;
  short int dRsignMC;
  double siguDC; 

  friend class TEv;       // for TEv::ImportClusters(const list<CsCluster*>&);

};

//
// Inlined functions
//

// constructor
inline THit::THit():
  // init const references
  IPlane(iPlane),
  IKine (iKine),
  IKin2 (iKin2),
  IOrig (iOrig),
  IHit  (iHit),
  Status(status),
  Mirror(mirror),
  U(u), V(v),
  SigU(sigu), SigV(sigv), SigT(sigt), 
  Y(y),Z(z),
  G0(g0), G1(g1), G2(g2),
  Time(time), DeltaR(deltaR), DRsignMC(dRsignMC),
  // init data members
  ptrCl(NULL), ptrClMirr(NULL),
  iPlane(-1),
  iKine (-1),
  iKin2 (-1),
  iOrig (-1),
  iHit  (-1),
  status(0),
  mirror(-1),
  sigt(-1), time(0), deltaR(0), dRsignMC(0), siguDC(0)
  
{};


// copy constructor
inline THit::THit(const THit& h): 
  // init const references
  IPlane(iPlane),
  IKine (iKine),
  IKin2 (iKin2),
  IOrig (iOrig),
  IHit  (iHit),
  Status(status),
  Mirror(mirror),
  U(u), V(v),
  SigU(sigu), SigV(sigv), SigT(sigt),
  Y(y),Z(z),
  G0(g0), G1(g1), G2(g2),
  Time(time), DeltaR(deltaR), DRsignMC(dRsignMC),
  // init data members
  ptrCl(h.ptrCl),
  ptrClMirr(h.ptrClMirr),
  iPlane(h.iPlane),
  iKine (h.iKine),
  iKin2 (h.iKin2),
  iOrig (h.iOrig),
  iHit  (h.iHit),
  status(h.status),
  mirror(h.mirror),
  u (h.u), v(h.v),
  sigu(h.sigu),sigv(h.sigv), sigt(h.sigt),
  y(h.y), z(h.z),
  g0(h.g0), g1(h.g1), g2(h.g2),
  setTrackID(h.setTrackID),
  setDigits(h.setDigits),
  time(h.time), deltaR(h.deltaR), dRsignMC(h.dRsignMC),  
  siguDC(h.siguDC)
{};

// operator "="
inline THit& THit::operator = (const THit& h)
{
  ptrCl      = h.ptrCl;
  ptrClMirr  = h.ptrClMirr;
  iPlane = h.iPlane;
  iKine  = h.iKine;
  iKin2  = h.iKin2;
  iOrig  = h.iOrig;
  iHit   = h.iHit;
  status = h.status;
  mirror=h.mirror;
  u = h.u; v = h.v;
  sigu = h.sigu; sigv = h.sigv; sigt = h.sigt; siguDC = h.siguDC; 
  y = h.y;  z = h.z;
  g0 = h.g0; g1 = h.g1; g2 = h.g2;
  time = h.time; deltaR = h.deltaR; dRsignMC = h.dRsignMC;
  setTrackID = h.setTrackID;
  setDigits  = h.setDigits;

  return(*this);
};

// Useful when the CsCluster pointed by THit has been changed for its mirror
inline void THit::setKine(int ikine) { iKine = ikine; }
#endif
  
