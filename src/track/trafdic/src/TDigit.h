/*!
  \class TDigit
  \brief Digit

  Contains raw detector information.
  (e.g wire number + TDC counts)
  

  \author Sergei.Gerassimov@cern.ch
*/

#ifndef TDigit_h
#define TDigit_h

#include <CsSTD.h>
#ifdef COMPASS_USE_OSPACE_STD
# include <ospace/std/set>
#else
# include <set>
#endif

class TDetect;
class CsDigit;

class TDigit {

public:
  static unsigned int NobjCreated;
  static unsigned int NobjDestructed;

  //! Constructor
  TDigit();
  TDigit(int iW, std::vector<float>& Info);
  //! Copy constructor
  TDigit(const TDigit&);
  //! assignment operator
  TDigit& operator = (const TDigit&);
  //! Destructor
  ~TDigit();

  const int& IWire;              //!< wire/strip/fiber/straw number
  const std::vector<float>& vDigInfo; //!< raw information (ADC, TDC, etc.) 
  
  //!  get pointer to corresponding CsDigit object
  CsDigit* PtrDig() const { return ptrDig;}

  bool operator < ( const TDigit& d) const ; // "less" operator
 
private:

  CsDigit* ptrDig;
  int iWire;
  std::vector<float> vecDigInfo;

  friend class TEv;       // for TEv::ImportClusters(const list<CsCluster*>&);

};

//
// Inlined functions
//

// constructor
inline TDigit::TDigit():
  // init const references
  IWire   (iWire),
  vDigInfo(vecDigInfo),
  // init data members
  ptrDig(NULL),
  iWire(-1)
{ };

// constructor
inline TDigit::TDigit(int iW, std::vector<float>& Info):
  // init const references
  IWire   (iWire),
  vDigInfo(vecDigInfo),
  // init data members
  ptrDig(NULL)
{ 
  NobjCreated++;
  iWire = iW; 
  vecDigInfo = Info; 
};


// copy constructor
inline TDigit::TDigit(const TDigit& h): 
  // init const references
  IWire   (iWire),
  vDigInfo(vecDigInfo)
{
  NobjCreated++;
  if(this != &h){
    (*this)=h;
  }
};

// destructor
inline TDigit::~TDigit()
{
  NobjDestructed++;
}

// operator "="
inline TDigit& TDigit::operator = (const TDigit& h)
{
  ptrDig     = h.ptrDig;
  iWire      = h.iWire;
  vecDigInfo = h.vecDigInfo;
  return(*this);
};

inline  bool TDigit::operator < ( const TDigit& d) const {
  return(iWire < d.iWire);
}

#endif
















