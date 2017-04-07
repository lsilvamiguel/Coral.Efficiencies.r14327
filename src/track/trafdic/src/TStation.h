// $Id: TStation.h 14069 2015-09-17 20:44:46Z lsilva $

#ifndef TStation_h
#define TStation_h

/*!
  \class TStation
  \brief Station of detectors

  \author Yann.Bedfer@cern.ch
*/

#include <CsSTD.h>

class TStation {


public:

  int IStation;            //!< Index of TStation in TSetup::vecStation
  int Type;                //! Type = uM, GEM, MWPC, etc... 
  std::vector<int> IPlanes;     //!< Indices of component planes
  int IPat, JPl;           //!< Index (/32,%32) of 1st plane of TStation in Coral zone 
  unsigned int Pat;        //!< Pattern of planes in Coral zone
  float URange, VRange;    //!< Range for searching for tri/quadruples
  double X0;               //!< Mean X abscissa

  //! Constructor
  TStation();
  //! Copy Constructor
  TStation(const TStation& p);

  // Automatic destructor 

  //! Operator =
  TStation& operator = (const TStation& p);

private:

  friend class TSetup; // for TSetup::Init()

};

//! Constructor
inline TStation::TStation(): 
  // init data members
  IStation(-1), Type(0),
  IPlanes(),
  IPat(0), JPl(0), Pat(0), URange(0), VRange(0)
{
  IPlanes.reserve(18);
};

//! Copy constructor
inline TStation::TStation(const TStation &s):
  // init data members
  IStation(s.IStation), Type(s.Type),
  IPat(s.IPat), JPl(s.JPl), Pat(s.Pat), URange(s.URange), VRange(s.VRange)
{
  IPlanes = s.IPlanes;                       //! copy content
};
  
//! operator "="
inline TStation& TStation::operator=(const TStation &s)
{
  IStation = s.IStation; Type = s.Type;
  IPat = s.IPat; JPl = s.JPl; Pat = s.Pat; URange = s.URange; VRange = s.VRange;
  IPlanes = s.IPlanes;      //! copy vector
  return(*this);
}
#endif











