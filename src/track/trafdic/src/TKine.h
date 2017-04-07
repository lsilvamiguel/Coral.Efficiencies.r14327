// $Id: TKine.h 13148 2011-12-28 16:55:25Z kbicker $

#ifndef TKine_h
#define TKine_h

#include <CsSTD.h>
#include <CsMCParticle.h>
#include <CsMCTrack.h>
#ifdef COMPASS_USE_OSPACE_STD
# include <ospace/std/set>
#else
# include <set>
#endif
#include <math.h>

/*!
  \brief MC track

  Simulated track (contains GEANTs KINE bank information and much more)

*/

class TKine {

public:
  static unsigned int NobjCreated;
  static unsigned int NobjDestructed;

  const double& E;              //!< Energy
  const int& IVtxOrig;          //!< Origin vertex number in TEv::vVtxMC

  //! Px, Py, Pz
  const double& P(int i)         const { return p[i]; }

  //! vector of references to TEv::vVtxMC
  const std::vector<int>& vVtxMCRef() const { return vecVtxMCRef; }

  //! vector of references to TEv::vHitMC 
  const std::vector<int>& vHitMCRef() const { return  vecHitMCRef; }

  //! vector of references to TEv::vHit
  const std::vector<int>& vHitRef()   const { return vecHitRef; }

  //! set of corresponding reconstructed track IDs
  const std::set<unsigned int, std::less<unsigned int> >& sTrackID()  const { return setTrackID; }

  //! add corresponding reconstructed track IDs
  void addTrackID(unsigned int id){ setTrackID.insert(id); }

  //! erase corresponding reconstructed track IDs
  void eraseTrackID(unsigned int id){ setTrackID.erase(id); }

  //! pointer to corresponding CsMCTrack object
  const CsMCTrack* PtrTrk()      const { return ptrTrk; }

  bool        isPileup() const;               //!< returns "true" if this MC track is from pileup or hallo 
  int         Q()        const;               //!< MC track charge
  std::string Name()     const;               //!< PDG particle name
  int         PDGid()    const;               //!< PDG particle id
  double      Mass()     const;               //!< PDG mass (Gev)
  double      MCmass()   const;               //!< Mass, calculated from energy and 3-momentum

  //! Constructor
  TKine();
  //! Destructor
  ~TKine();
  //! Copy constructor
  TKine(const TKine&);
  //! Operator =
  TKine& operator = (const TKine& k);

  double Pinv() const;        //!< returns q/P

private:
  
  double e;                     //!< Energy
  int iVtxOrig;                 //!< Origin vertex number in vVtxMC
  double p[3];                  //!< Px, Py, Pz
  std::vector<int> vecVtxMCRef; //!< vector of vertices
  std::vector<int> vecHitRef;   //!< vector of references to vHit
  std::vector<int> vecHitMCRef; //!< vector of references to vHitMC 
  std::set<unsigned int, std::less<unsigned int> > setTrackID;  //! set of corresponding reconstructed track IDs
  CsMCTrack* ptrTrk; 

  friend class TEv;      // for TEv::GetMCInfo()

};

// Constructor
inline TKine::TKine(): 
  E(e),
  IVtxOrig(iVtxOrig),

  vecVtxMCRef(),
  vecHitRef (),
  vecHitMCRef(),
  setTrackID()
{
  NobjCreated++;
  iVtxOrig=-1;
  ptrTrk = 0; 
  vecVtxMCRef. reserve(20);
  vecHitRef.  reserve(200);
  vecHitMCRef.reserve(300);
};


// Copy constructor
inline TKine::TKine(const TKine& k):
  E(e),
  IVtxOrig(iVtxOrig),

  vecVtxMCRef(),
  vecHitRef (),
  vecHitMCRef(),
  setTrackID()
{
  NobjCreated++;
  e           = k.e;
  iVtxOrig    = k.iVtxOrig;
  for(int i=0; i<3; i++) p[i]=k.p[i];
  vecVtxMCRef = k.vecVtxMCRef; 
  vecHitMCRef = k.vecHitMCRef;
  vecHitRef   = k.vecHitRef;
  setTrackID  = k.setTrackID;
  ptrTrk      = k.ptrTrk;
};

// Destructor
inline TKine::~TKine() 
{
  NobjDestructed++;
}

// operator = 
inline TKine& TKine::operator = (const TKine& k)
{
  e           = k.e;
  iVtxOrig    = k.iVtxOrig;
  for(int i=0; i<3; i++) p[i]=k.p[i];
  vecVtxMCRef = k.vecVtxMCRef; 
  vecHitMCRef = k.vecHitMCRef;
  vecHitRef   = k.vecHitRef;
  setTrackID  = k.setTrackID ;
  ptrTrk      = k.ptrTrk;

  return(*this);
};


// Methods
inline double TKine::Pinv() const
{
  if(Q() == 0) {
    return(1./sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]));
  } else {
    return (Q()/sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]));
  }
};

inline double TKine::MCmass() const 
{
  return ( sqrt(fabs(e*e - (p[0]*p[0]+p[1]*p[1]+p[2]*p[2]))) );
};

inline int TKine::PDGid() const
{
  return(ptrTrk->getParticle()->getNumber());
};

inline int TKine::Q() const
{
  return(ptrTrk->getParticle()->getCharge());
};

inline double TKine::Mass() const
{
  return(ptrTrk->getParticle()->getMass());
};

inline std::string TKine::Name() const
{
  std::string sig;
  if(Q() == 0) sig=" 0"; if(Q() <  0) sig=" -"; if(Q() >  0) sig=" +";
  std::string name(ptrTrk->getParticle()->getName());
  name+=sig;
  return(name);
};
#endif









