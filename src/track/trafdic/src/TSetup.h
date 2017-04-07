// $Id: TSetup.h 14094 2015-11-06 15:28:48Z lsilva $

#ifndef TSetup_h
#define TSetup_h

/*!
  \class TSetup
  \brief COMPASS spectrometer setup information, arranged to facilitate its use in tracking algorithms
  \warning Only one instance of this class is allowed
  \author Sergei.Gerassimov@cern.ch
*/

/*
  Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/TSetup":
   i) Added vector "vecStation".
  ii) Hinging points used in bridging over the muFilter.
 */

#include <CsSTD.h>
#include <cassert>
#include "TDetect.h"
#include "TPlane.h"
#include "TStation.h"
#include <map>

class CsZone;

class TSetup {

public:
  static unsigned int NobjCreated;
  static unsigned int NobjDestructed;

  TSetup();   //!< Constructor
  ~TSetup();  //!< Destructor

  /*! Returns pointer to this object

    Example of use:
    \code
    TSetup::Ptr()->vDetect().size();
    \endcode
    (total number of detectors in the setup) 

   */
  static TSetup* Ptr();

  /*! Returns reference to this object

    Example of use:
    \code
    TSetup::Ref().Nmags;
    \endcode
    (number of magnets in the setup) 

   */
  static const TSetup& Ref();

  //! Physical tracking detectors (event-independent information)  
  const std::vector<TDetect>& vDetect() const { return vecDetect; }

  //! Accessor to i-th TDetect object of the vecDetect
  const TDetect& vDetect(int i) const { return vecDetect[i]; }

  //! Logical tracking detectors with links to event-dependent information 
  const std::vector<TPlane> & vPlane () const { return vecPlane; }

  //! Accessor to i-th TPlane of vecPlane
  const TPlane& vPlane (int i) const { return vecPlane[i]; }

  //! Vector of projection angles (deg.). FYI only. Use rotation matrices for calculations.
  const std::vector<int> & vProj () const { return vecProj; }

  //! Reference to associated TPlane for the detector with argument ID 
  const TPlane& Id2Plane(int detID) const ;

  //! Index of associated TPlane for the detector with argument ID. Default =-1.
  int Id2iPlane(int detID) const ;

  //! Reference to associated TDetect for the plane with argument index  
  const TDetect& iPlane2Detect(int iPlane) const ;
  
  //! CsZone* for detectors group #igr. (# of det. groups >= number of zones)
  CsZone* Group2Zone(int igr) const;

  //! Vector of vPlane indices of the first detector in every detectors group
  const std::vector<int>& vIplFirst() const { return vecIplFirst; }

  //! Vector of vPlane indices of the last  detector in every detectors group
  const std::vector<int>& vIplLast()  const { return vecIplLast; }

  //! Set of IDs of detectors to be ignored by Traffic 
  const std::set<int>& sIgnoreIDs()  const { return setIgnoreIDs; }
  
  //! Remove event-dependent information from objectrs in vecPlane
  void Reset();

  //! Detector stations
  const std::vector<TStation>& vStation () const { return vecStation; }

  // Initialised by private data members in the constructor
  const int&     NMags;              //!< Number of magnets (3 max)
  const int*     MagType;            //!< Magnet type flags (array[3]) 
  const float*   MagScale;           //!< Field scaling factor (array[3])
  const float*   MagFlag1;           //!< Flag #1 from det.dat (array[3])
  const int*     MagFlag2;           //!< Flag #2 from det.dat (array[3])
  const double*  MagFieldInt;        //!< Field integral kG*cm (array[3])
  const double  &BeamTAngle;         //!< Angle about y axis due to target dipole
  const double*  TargetCenter;       //!< Target center (x,y,z) (array[3]) 
  const int     &FirstMuWallA;       //!< First muWallA detector

  const double  *MuFilterHinges;     //!< Hinging points for bridging over muFilter (array[3])
  //! Coordinates (x,y,z for i=0,1,2) of the magnet "m" center
  const double&  MagCenter(int m, int i)  const { return magCenter[m][i]; }
  //! Plane "ipl" is one of the MuWallA detector (either a MA or HG02)
  bool InMuWallA(int ipl) const {
    if (firstMuWallA<=ipl && ipl<=lastMuWallA)
      return static_cast<unsigned long long>(0x1)<<ipl-firstMuWallA&bMuWallA;
    else return false;
  }

  void Update();  //! Update field integral for each magnet which scale factor turns out to have changed

  //! Returns "true" if X is inside zone where material map is defined
  bool InMaterialMap(double X) const ;

  //! Returns "true" if there is material map defined beetween X0 and X1
  bool IsMaterialBetween(double X0, double X1) const ;

  //! Accessor to vector of material map borders
  const std::vector<float> &vMaterialBorders() const { return vecMaterialBorders; }

  //! Returns "true" if GEM spacers available
  bool GEMSpacers() const { return mapID2GEMSpacers.size()!=0; }

  //! Returns array of intervals of channel#'s defining ``spacers hits''
  const float *GEMSpacerBounds(int id) const;

  //! Computes and returns (as array of doubles or array of helices) the hinging points of the bridging over the muFilter
  void ComputeHinges(double *hinges, THlx *Hs = 0) const;
private:

  static TSetup* address;

  bool Init();    //!< Initialization called from the constructor.

  std::vector<TDetect> vecDetect; //!< Physical tracking detectors (event-independent information)
  std::vector<TPlane>  vecPlane;  //!< Logical  tracking detectors with links to event-dependent information 
  std::vector<int>     vecProj;   //!< Vector of measurement angles (integer: not precise!) FYI only. 
  std::vector<int>   vecIplFirst; //!< Vector of vPlane index of the first detector in every detectors group 
  std::vector<int>   vecIplLast;  //!< Vector of vPlane index of the last  detector in every detectors group 
  std::set<int> setIgnoreIDs;     //!< Set of detector IDs to be completely ignored
  std::vector<float> vecMaterialBorders; //! Material map borders along X 
  std::vector<TStation> vecStation;      //!< Detector stations
  std::map<int,int,std::less<int> > mapIDvPlane;         //!< ID -> Index in vecPlane
  std::map<int,float*,std::less<int> > mapID2GEMSpacers; //!< ID -> Intervals of channel#'s

  std::map<int,CsZone*,std::less<int> > mapGroup2Zone;   //!< Det. group -> CORAL zone

  int      nMags;                   // Number of magnets (3 max)
  int      magType  [3];            // Magnet type flags
  double   magCenter[3][3];         // Magnets' center coordinates [3][3] 
  float    magScale [3];            // Filed scaling factor
  float    magFlag1 [3];            // Flag 1 from det. table
  int      magFlag2 [3];            // Flag 2 from det. table
  double   magFieldInt [3];         // Magnets' field integrals
  double   beamTAngle;              //!< Angle about y axis due to target dipole

  double   targetCenter[3];         // Target center

  double   muFilterHinges[4];       // Hinging points for bridging over muFilter

  //! muWallA detectors (MAs and HG02): bit pattern (used for a fast search of such detectors in track hit list) in two words
  int firstMuWallA, lastMuWallA; unsigned long long bMuWallA;

  bool usingROOTGeometry;
};

// Inline functions

inline TSetup* TSetup::Ptr()
{
  // No check on existence of the object. Normaly, TSerup::Ref() has to be used.
  return(address);
}
inline const TSetup& TSetup::Ref()
{
  if (address != 0) {
    return(*address);
  }
  std::cout<<"TSetup::Ref() ==> The object of TSetup class is not yet created or already destructed"<<std::endl;
  assert(false);
  return(*address); // just to get rid of compiler warnings
}

inline const TPlane& TSetup::Id2Plane(int detID) const {
  std::map<int, int, std::less<int> >::const_iterator im;
  im = mapIDvPlane.find(detID);
  if(im != mapIDvPlane.end()){ // found
    int index=(*im).second;
    if( index < 0 || index >= int(vecPlane.size())){
      std::cout<<"TSetup::Id2Plane ==> wrong vPlane index "<<index<<" for detector with ID = "<<detID<<std::endl;
      assert(false);
    }
    return(vecPlane[index]);
  }                  // not found
  std::cout<<"TSetup::Id2Plane ==> can not find vPlane index for detector with ID = "<<detID<<std::endl;
  assert(false);
  return(vecPlane[0]); // just to get rid of compiler warnings
}

inline int TSetup::Id2iPlane(int detID) const {
  std::map<int, int, std::less<int> >::const_iterator im = mapIDvPlane.find(detID);
  if (im!=mapIDvPlane.end()) { // Found
    int index = (*im).second;
    if (index<0 || index>=int(vecPlane.size())) {
      std::cout<<"TSetup::Id2iPlane ==> wrong vPlane index "<<index<<" for detector with ID = "<<detID<<std::endl;
      assert(false);
    }
    return index;
  }
  else
    // "detID" not found in "mapIDvPlane" means the corresponding CsDetector
    // has not been imported to TSetup, e.g. because it's been deliberately
    // ignored, cf. "Det2Ignore" in "TSetup::Init".
    return -1;
}

inline const float *TSetup::GEMSpacerBounds(int id) const {
  std::map<int,float*,std::less<int> >::const_iterator im = mapID2GEMSpacers.find(id);
  if (im!=mapID2GEMSpacers.end()) return (*im).second;
  else                            return 0;
}

inline CsZone *TSetup::Group2Zone(int igr) const {
  std::map<int,CsZone*,std::less<int> >::const_iterator im = mapGroup2Zone.find(igr);
  if (im!=mapGroup2Zone.end()) return (*im).second;
  std::cout<<"TSetup::Group2Zone ==> can not find corresponding CsZone* for detector's group # "<<igr<<std::endl;
  assert(false);
  return 0; // just to get rid of compiler warnings
}

inline const TDetect& TSetup::iPlane2Detect(int iPlane) const {
  if(iPlane < 0 || iPlane >= int(vecPlane.size())){
    std::cout<<"TSetup::iPlane2Detect ==> can not find corresponding detector for the vPlane index = "<<iPlane<<std::endl;
    assert(false);
  }

  const int& iDet = vecPlane[iPlane].iDetRef;

  if(iDet < 0 || iDet >= int(vecDetect.size())){
    std::cout<<"TSetup::iPlane2Detect ==> iDetRef seems wrong ("<<iDet<<") for the vPlane index = "<<iPlane<<std::endl;
    assert(false);
  }

  return(vecDetect[iDet]);
}

//! Reset event dependent information
inline void TSetup::Reset(){

  for(int ipl=0; ipl < int(vecPlane.size()); ipl++){
    vecPlane[ipl].Reset();
  }

};

//! Returns "true" if X is inside zone where material map is defined
inline  bool TSetup::InMaterialMap(double X) const {

  if (usingROOTGeometry)
    return true;

  for(int i = 0; i < int(vecMaterialBorders.size()); i+=2) {
    if(vecMaterialBorders[i] <= X && X < vecMaterialBorders[i+1]) return(true);
  }
  return(false);
  
}

//! Returns "true" if there is material map defined beetween X0 and X1
inline  bool TSetup::IsMaterialBetween(double X0, double X1) const {

  if (usingROOTGeometry)
    return true;

  if(X0 > X1) std::swap(X0,X1);
  for(int i = 0; i < int(vecMaterialBorders.size()); i++) {
    if(X0 <= vecMaterialBorders[i] && vecMaterialBorders[i] < X1) return(true);
  }
  return(false);
  
}
#endif
