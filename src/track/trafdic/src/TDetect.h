// $Id: TDetect.h 14069 2015-09-17 20:44:46Z lsilva $

#ifndef TDetect_h
#define TDetect_h
/*!
  \class TDetect
  \brief Tracking detector

  Object of this class contains static (event-independent)
  part of information about one tracking detector's plane
  in the form, more appropriate for the track finding
  purposes.
  (for event-dependent info see TPlane)
  
  Mostly it's just an interface to CORAL getector geomery classes

  \author Sergei.Gerassimov@cern.ch

*/


/*
  Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/TDetect":
   i) Added member data "uedge".
  ii) Many of the data members turned into "float". This in order to speed up
  processing, particularly in "InActive" for which exists a version w/
  arguments of type "float".
 iii) Empty zone type. And "InMassive" method making use of it.
 */

#include <cassert>
#include <cmath>
#include <iostream>
#include <string>

#include "CsDetector.h"

enum TDZType{NO, RECTANGULAR, CIRCULAR};

class TDetect {

public:
  static unsigned int NobjCreated;
  static unsigned int NobjDestructed;

  const std::string& Name;          //!< detector name
  const int& IDet;             //!< detector ID
  const short int& IType;      //!< detector's Type 
  const short int& Kind;       //!< detector's kind (0 - wire/strip; 1 - drift)
  const double& RadLen;        //!< radiation length of detector materials
  const double& Pitch;         //!< Pitch
  const unsigned int& Nwires;  //!< Number of wires
  const std::vector<double>& W2Pos; //!< Wire positions (for var. pitch detectors only)
  const double& Uorig;         //!< coordinate of the first wire in the WRS
  const float &UEdge;          //!< Coordinate of the <0 detector's edge in the WRS
  const double& Resol;         //!< Detector resolution
  const double& TResol;        //!< Detector time resolution
  const float &Ca;             //!< Cos of mesurement angle (WRS Y to MRS Y) 
  const float &Sa;             //!< Sin of mesurement angle (WRS Y to MRS Y)
  const float &Range;          //!< Length of active area in measurement direction   cm
  const float &Cframe;         //!< Cos of detector frame rotation angle (in MRS) 
  const float &Sframe;         //!< Sin of detector frame rotation angle (in MRS)
  const int& IPlane;           //!< Reference to vPlane vector ( set in TSetup::Init() )
  //! 3 Coordinates of center of detector in main reference system (MRS)
  const double& X(int i)       const  {return x[i];}
  //! 3 Coordinates of center of detector in WRS
  const double& XR(int i)      const  {return xR[i];}
  //! 3 detector half sizes (in DRS ) cm
  const float &Siz(int i)     const  {return siz[i];}
  //! 2 dead zone sizes
  float DZSize(int i) const {
    switch( DZtype )
    {
    case RECTANGULAR : return (i==0 ? DZydim : DZzdim);
    case CIRCULAR    : return sqrt(DZydim);
    default          : return 0;
    }
  }
  //! returns pointer to corresponding CsDetector
  CsDetector* PtrDet() const {return ptrDet;}

  // MRS - main spectrometer reference system : system like in COMGEANT.
  // (righthand system wth X - along the beam, Z - vertical with the origine in the center of PT)

  // DRS - detector reference system
  // System shifted to center of detectors parallelepiped and rotated together with it
  // by Cframe/Sframe and Ct/St

  // WRS - Wire reference system: origin like in MRS,  but rotated around X 
  // to make Y perpendicular to wires)


  TDetect();                                   //!< construcror
  ~TDetect();                                  //!< destrucror
  TDetect(const TDetect& t);                   //!< copy constructor prototype
  TDetect& operator = (const TDetect& t);      //!< = operator
  bool operator < ( const TDetect& t) const;   //!< < operator


  bool InFrame (double y, double z)  const;    //!< Check if (y,z) is w/in detector frame
  bool InActive(float y,  float z)   const;    //!< Check if (y,z) is in active area
  bool InActive(double y, double z)  const;    //!< Check if (y,z) is in active area
  bool InActive(double c[2], double cov[3]) const;   //!< Check if (y,z) with errors (cov[3]) is in active area
  bool WinRange(float y,  float z)   const;    //!< Check if (y,z) w/in measurement range
  bool InMassive(double y, double z) const;    //!< Check if (y,z) is in active area or in dead zone if latter non empty.

  friend class TDisplay;

private:

  CsDetector* ptrDet; 
  std::string name;
  int iDet;
  short int iType;
  short int kind;
  double radLen;
  double pitch;
  unsigned int nwires;
  std::vector<double> w2pos; // wire positions for var. pitch detectors. Empty for others.
  double uorig;
  float uedge;
  double resol;
  double tresol;
  float ca, sa;
  float range;
  float cframe, sframe;
  int iPlane;
  double x[3];
  double xR[3];
  float siz[3], hsizY, hsizZ;  // Sizes of the sensitive area, and Y/Z 1/2sizes

  // Dead Zone (or empty zone) parameters
  TDZType DZtype;  // Shape (NO, RECTANGULAR, CIRCULAR) 
  float   DZydim;  // Y 1/2size in case rectangular shape or radius squared (circular)
  float   DZzdim;  // Z 1/2size in case rectangular shape
  float   DZCadrs; // Cos of dead zone rotation angle in DRS
  float   DZSadrs; // Sin of dead zone rotation angle in DRS
  float   DZydrs;  // Y position of dead zone in DRS
  float   DZzdrs;  // Z position of dead zone in DRS
  float   DZymrs;  // Y position of dead zone in WRS
  float   DZzmrs;  // Z position of dead zone in WRS
  TDZType EZtype;  // = "DZtype" or "NO" if dead zone is empty space.

  friend class TSetup;   // for TSetup::Init()

};


// Inlined function

inline TDetect::TDetect():
  // init const references    
  Name   (name),
  IDet   (iDet),
  IType  (iType),
  Kind   (kind),
  RadLen (radLen),
  Pitch  (pitch),
  Nwires (nwires),
  W2Pos  (w2pos),
  Uorig  (uorig),
  UEdge  (uedge),
  Resol  (resol),
  TResol (tresol),
  Ca     (ca),    
  Sa     (sa),
  Range  (range),
  Cframe (cframe),
  Sframe (sframe),
  IPlane (iPlane),

  // init data members
  ptrDet(NULL),
  name   (""),
  iDet   (-1),
  iType  (-1),
  kind   (-1),
  radLen (1.E10),
  pitch  (0),
  nwires (0),
  uorig  (0),
  uedge  (0),
  resol  (0),
  tresol (0),
  ca     (0),    
  sa     (0),
  range  (0),
  cframe (0),
  sframe (0),
  iPlane (-1),
  DZtype (NO),
  DZydim (0),
  DZzdim (0),
  DZCadrs(1),
  DZSadrs(0),
  DZydrs (0),
  DZzdrs (0),
  DZymrs (0),
  DZzmrs (0),
  EZtype (NO)
{
  NobjCreated++;
  for(int i = 0; i < 3; i++){
    x [i]  = 0;
    xR[i]  = 0;
    siz[i] = 0;
  }
  hsizY=hsizZ = 0;
};


// Copy constructor
inline TDetect::TDetect(const TDetect& d):
  // init const references    
  Name   (name),
  IDet   (iDet),
  IType  (iType),
  Kind   (kind),
  RadLen (radLen),
  Pitch  (pitch),
  Nwires (nwires),
  W2Pos  (w2pos),
  Uorig  (uorig),
  UEdge  (uedge),
  Resol  (resol),
  TResol (tresol),
  Ca     (ca),    
  Sa     (sa),
  Range  (range),
  Cframe (cframe),
  Sframe (sframe),
  IPlane (iPlane),
  // copy data members
  ptrDet (d.ptrDet),
  name   (d.name),
  iDet   (d.iDet),
  iType  (d.iType),
  kind   (d.kind),
  radLen (d.radLen),
  pitch  (d.pitch),
  nwires (d.nwires),
  w2pos  (d.w2pos),
  uorig  (d.uorig),
  uedge  (d.uedge),
  resol  (d.resol),
  tresol (d.tresol),
  ca     (d.ca),    
  sa     (d.sa),
  range  (d.range),
  cframe (d.cframe),
  sframe (d.sframe),
  iPlane (d.iPlane),
  DZtype (d.DZtype),
  DZydim (d.DZydim),
  DZzdim (d.DZzdim),
  DZCadrs(d.DZCadrs),
  DZSadrs(d.DZSadrs),
  DZydrs (d.DZydrs),
  DZzdrs (d.DZzdrs),
  DZymrs (d.DZymrs),
  DZzmrs (d.DZzmrs),
  EZtype (d.EZtype)
{
  NobjCreated++;
  for(int i = 0; i < 3; i++){
    x [i]  = d.x[i];
    xR[i]  = d.xR[i];
    siz[i] = d.siz[i];
  }
  hsizY = d.hsizY; hsizZ = d.hsizZ;
};

inline TDetect::~TDetect()
{
  NobjDestructed++;
}

// operator = 
inline TDetect& TDetect::operator = (const TDetect& d)
{
  ptrDet = d.ptrDet;
  name.erase();
  name  += d.name; // just "=" gives warning :-(
  iDet   = d.iDet;
  iType  = d.iType;
  kind   = d.kind;
  radLen = d.radLen;
  iPlane = d.iPlane;
  pitch  = d.pitch;
  nwires = d.nwires;
  w2pos  = d.w2pos;
  uorig  = d.uorig;
  uedge  = d.uedge;
  resol  = d.resol;
  tresol = d.tresol;
  range  = d.range;
  cframe = d.cframe;
  sframe = d.sframe;
  ca     = d.ca;    
  sa     = d.sa;    
  DZtype = d.DZtype;
  DZydim = d.DZydim;
  DZzdim = d.DZzdim;
  DZCadrs= d.DZCadrs;
  DZSadrs= d.DZSadrs;
  DZydrs = d.DZydrs;
  DZzdrs = d.DZzdrs;
  DZymrs = d.DZymrs;
  DZzmrs = d.DZzmrs;
  EZtype = d.EZtype;
  for(int i = 0; i < 3; i++){
    x [i]  = d.x[i];
    xR[i]  = d.xR[i];
    siz[i] = d.siz[i];
  }
  hsizY = d.hsizY; hsizZ= d.hsizZ;
  return *this;
}

inline bool TDetect::InFrame(double y, double z) const 
{
  y -= x[1]; z -= x[2];
  float ydrs = cframe*y+sframe*z;
  float zdrs = cframe*z-sframe*y;
  float y1 = fabs(ydrs), z1 = fabs(zdrs);
  if (y1<hsizY && z1<hsizZ) return true;
  else                      return false;
};


inline  bool TDetect::operator < ( const TDetect& d) const {
  double x1  = x[0], x2  = d.x[0];
  if (x1==x2) {
    // Special cases: pixel GEM and pixel MM
    const std::string &name1 = Name, &name2 = d.Name;
    if (name1==name2) return true;
    if ((name1.find("GP")==0 && name2.find("GP")==0) ||
	(name1.find("GM")==0 && name2.find("GM")==0)) {
      // [Pixel]GEMs: Enforcing ordering (P1<U<V < Y<X<P2) needed by trafdic
      char c1 = name1[4], c2 = name2[4];
      std::string corrSeq("YXPUV");
      short t1 = corrSeq.find(c1), t2 = corrSeq.find(c2);
      return t1<=t2 ;
    }
    if (name1.find("MP")==0 && name2.find("MP")==0) {
      // PixelMMs: Enforcing ordering (M1<X|V or Y|U<M2) needed by trafdic
      char c1 = name1[4]; if (c1=='P') c1='M';// Older MPs: 'P' instead of 'M'
      char c2 = name2[4]; if (c2=='P') c2='M';
      char *cXYUV = c1=='M' ? &c2 : (c2=='M' ? &c1 : 0);
      if (cXYUV) {
	if      (*cXYUV=='X' || *cXYUV=='V') return cXYUV==&c2;
	else if (*cXYUV=='Y' || *cXYUV=='U') return cXYUV==&c1;
	else cXYUV = 0;
      }
      if (!cXYUV)
	CsErrLog::msg(elFatal,__FILE__,__LINE__,
		      "PixelMMs \"%s\",\"%s\" have same abscissa(=%.5f) while"
		      " they're not a (X|Y|U|V,M|P) pair of coordinates",
		      name1.c_str(),name2.c_str(),x1);
    }
    // Other cases like ST or MB ('a','b','c' or 'r','l') slices
    double y1 = x[1], y2 = d.x[1];
    if (y1==y2) {
      double z1 = x[2], z2 = d.x[2];
      return z1<z2;
    }
    else return y1<y2;
  }
  return x1<x2;
}

/*!
  Fastest possible check that argument (y,z) is w/in active area.
*/

inline bool TDetect::InActive(float y, float z) const 
{
  // ***** TRANSFORM y,z to DRS *****
  float y_shift = y - x[1];
  float z_shift = z - x[2];
  float ydrs = cframe*y_shift + sframe*z_shift;
  float zdrs = cframe*z_shift - sframe*y_shift;

  if (fabs(ydrs)>=hsizY || fabs(zdrs)>=hsizZ) return false;   // out of frame
#ifdef InActive_CHECK_URANGE
  //  U range only matters for stereo planes in DC like detectors, where
  // sensitive area has a rhombus shape. (And then it still isn't required
  // to prevent a hit from being picked-up during track search, for hits
  // outside U range won't be available in any case (by definition of U range).
  // The check is still needed to determine the ``firing expectancy'',
  // of a given detector, for a given track. Which we do in 2 cases:
  // when determining the ``expected'' hit map and when determining the
  // ``firing efficiency'' of a given candidate track (quantity that is
  // checked in the track search to qualify the candidate).)
  //  Note: TDetect::range = size in the case of "HI04/5", cf. "TSetup::Init".
  float uwrs = ca*y+sa*z-UEdge;            // "u" refered to edge of detector
  if (uwrs<0 || Range<uwrs) return false;  // Out of measurement range
#endif

  // ***** CHECK DEAD ZONE *****
  float ydrs_shift = ydrs - DZydrs;
  float zdrs_shift = zdrs - DZzdrs;

  switch( DZtype )
    {
    case RECTANGULAR :
      { 
	float ydzrs = DZCadrs*ydrs_shift + DZSadrs*zdrs_shift;
	if(fabs(ydzrs)>DZydim) return true;
	float zdzrs = DZCadrs*zdrs_shift - DZSadrs*ydrs_shift;
	return (fabs(zdzrs)>DZzdim);
      }
    case CIRCULAR    : 
      return (ydrs_shift*ydrs_shift + zdrs_shift*zdrs_shift > DZydim);
    case NO          :
    default :
      return true;      
    }  
  return true;
}

#endif
