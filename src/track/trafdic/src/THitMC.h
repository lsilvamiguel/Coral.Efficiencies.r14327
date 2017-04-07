// $Id: THitMC.h 14069 2015-09-17 20:44:46Z lsilva $

#ifndef THitMC_h
#define THitMC_h

#include "TSetup.h"
#include "THit.h"
#include "CsMCTrkHit.h"

class CsMCHit;

/*!

  \brief MC Hit

  Simplified MC "hit" class.
  Object of this class contains information on exact coordinate of the space point
  where MC trajectory had crossed centre of the detector's sensitive volume
*/

class THitMC {

public:
  static unsigned int NobjCreated;
  static unsigned int NobjDestructed;

  //! Constructor
  THitMC();
  //! Destructor
  ~THitMC();
  //! Copy Constructor
  THitMC(const THitMC& h);
  //! Operator =
  THitMC& operator = (const THitMC& h);
  

  const int&     IKine;       //!< index of the track in vecKine
  const int&     IOrig;       //!< 0 if the hit was produced by "original" track, 
  const double&  DeltaT;      //!< Delay time 
  const int&     IDet;        //!< detector identifier
  //!  set of indices of corresponding THit objects in vHit vector 
  const std::set<int, std::less<int> >& sHitRef() const {return setHitRef;}
  //!  coordinate of the center of trajectory in the sens. vol. 
  const double&  Xyz(int i) const { return xyz[i]; }
  //!  pointer to the corresponding CsMCHit object
  const CsMCHit* PtrHit()   const { return ptrHit;}
  //! reference to TDetect object
  const TDetect& DetRef()   const;
  //! Print MC hit info
  void Print() const;

private:

  int    iKine;
  int    iOrig;
  double deltaT;
  int    iDet;
  double xyz[3];
  std::set<int, std::less<int> > setHitRef;
  CsMCHit* ptrHit;

  friend class TEv;      // for TEv::GetMCInfo()

};


//Constructor
inline THitMC::THitMC():
  // init const. ref.
  IKine(iKine),
  IOrig(iOrig),
  DeltaT(deltaT),
  IDet(iDet),
  // init data members
  iKine(-1),
  iOrig(-1),
  deltaT(0),
  iDet(-1),
  ptrHit(0)
{
  for(int i=0; i<3; i++) xyz[i]=0.;
  NobjCreated++;

};

// Copy constructor
inline THitMC::THitMC(const THitMC& h):
  // init const. ref.
  IKine(iKine),
  IOrig(iOrig),
  DeltaT(deltaT),
  IDet(iDet),
  // init data members
  iKine(h.iKine),
  iOrig(h.iOrig),
  deltaT(h.deltaT),
  iDet(h.iDet),
  setHitRef(h.setHitRef),
  ptrHit(h.ptrHit)
{
  for(int i=0; i<3; i++) xyz[i] = h.xyz[i];
  NobjCreated++;
};

// operator = 
inline THitMC& THitMC::operator = (const THitMC& h)
{
  iKine  = h.iKine;
  iOrig  = h.iOrig;
  deltaT = h.deltaT;
  iDet   = h.iDet;
  ptrHit = h.ptrHit;
  setHitRef=h.setHitRef;
  for(int i=0; i<3; i++) xyz[i] = h.xyz[i];
  return(*this);

}

//Destructor
inline THitMC::~THitMC()
{
  NobjDestructed++;
}

inline const TDetect& THitMC::DetRef() const
{
  const TSetup& setup = TSetup::Ref();
  return (setup.vDetect(setup.Id2Plane(iDet).IDetRef));
};

inline void THitMC::Print() const
{

  std::cout<<"MC hit from MC track # "<<iKine
      <<"  iOrig = "<< iOrig <<"  with deltaT = "<<deltaT<<" on Traffic detector # "<<iDet<<std::endl
      <<"X/Y/Z = "<<xyz[0]<<"  "<<xyz[1]<<"  "<<xyz[2]<<std::endl;
  if(setHitRef.size() == 0){
    std::cout<<"NO hits produced by this MC hit"<<std::endl; 
  } else {
    std::cout<<"Associated hits : ";
    for(std::set<int>::iterator is = setHitRef.begin(); is != setHitRef.end(); is++){
      std::cout<<(*is)<<" ";
    }
    std::cout<<std::endl;
  }
  CsMCTrkHit *csh = dynamic_cast<CsMCTrkHit*>(ptrHit);
  if (csh) {
    printf("CsMCTrkHit %.3f GeV\n",csh->getP().mag());
  }
  std::cout<<std::endl;
};


#endif
