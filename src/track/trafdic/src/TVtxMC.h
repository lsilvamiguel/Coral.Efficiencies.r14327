#ifndef TVtxMC_h
#define TVtxMC_h

#include <CsSTD.h>

class CsMCVertex;

/*!
  \brief MC Vertex 

  Simplified GEANT's vertex class
  to be used for tracking debug
*/

class TVtxMC {

public:
  static unsigned int NobjCreated;
  static unsigned int NobjDestructed;

  const double& TOF;               //!< Time of flight
  const int& IOrig;                //!< TEv::vKine track number of mother particles
  //! Vx, Vy,Vz - vertex coordinates
  const double& V(int i)        const   { return v[i]; }
  //!  indecies (in TEv::vKine) of tracks in the vertex
  const std::vector<int>& vTrkRef()  const   { return vecTrkRef; }
  //!  pointer to corresponding CsMCVertex object
  const CsMCVertex* PtrVtx()    const   { return ptrVtx;}

  //! Constructor
  TVtxMC();
  //! Copy constructor
  TVtxMC(const TVtxMC&);
  //! Destructor
  ~TVtxMC();
  //! Operator =
  TVtxMC& operator = (const TVtxMC& h);

private:
  double tof;                  // Time of flight (us)
  int iOrig;                   // TEv::vKine track number of mother particles
  double v[3];                 // Vx, Vy,Vz
  std::vector<int> vecTrkRef;  // vector of indecies of tracks in the vertex
  CsMCVertex* ptrVtx;

  friend class TEv;            // for TEv::GetMCInfo()
};

// Constructor
inline TVtxMC::TVtxMC():
  TOF(tof),
  IOrig(iOrig),

  vecTrkRef(),
  ptrVtx(0)
{
  NobjCreated++;
  vecTrkRef.reserve(50);
};

// Copy constructor
inline TVtxMC::TVtxMC(const TVtxMC& x): 
  TOF(tof),
  IOrig(iOrig),
  vecTrkRef()
{
  NobjCreated++;
  tof    = x.tof;
  iOrig  = x.iOrig;
  for(int i=0; i<3; i++) v[i]=x.v[i];
  vecTrkRef.reserve(x.vecTrkRef.capacity());
  vecTrkRef = x.vecTrkRef;
  ptrVtx    = x.ptrVtx;
}
// Destructor
inline TVtxMC::~TVtxMC()
{
  NobjDestructed++;
}

// operator = 
inline TVtxMC& TVtxMC::operator = (const TVtxMC& x)
{
  tof    = x.tof;
  iOrig  = x.iOrig;
  for(int i=0; i<3; i++) v[i]=x.v[i];
  vecTrkRef = x.vecTrkRef;
  ptrVtx    = x.ptrVtx;
  return(*this);
}


#endif














