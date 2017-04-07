// $Id: TOpt.h 14069 2015-09-17 20:44:46Z lsilva $

#ifndef TOpt_h
#define TOpt_h
/*!
  \class TOpt
  \brief  TraFFiC options

  Parameters, modes, cuts etc. for TraFFiC, taken from CsOpt.
  
  Abstract class with static data members
  to provide direct, fast and global access to flags and cut

  \author Sergei.Gerassimov@cern.ch
*/
#include "CsSTD.h"

class TOpt {
public:

  static unsigned int NobjCreated;
  static unsigned int NobjDestructed;

  //! data
  static int     Graph    [10];       //!< Graphics switches
  static int     Print    [10];       //!< Debug printout flags
  static int     iCut     [40];       //!< Misc. cuts (int)
  static int     iPRpar  [120];       //!< Pattern Recognition parameters (int)
  static int     DetOff  [100];       //!< Detectors' IDs, excluded from the pattern recognition
  static int     ReMode  [100];       //!< Reconstruction mode switches
  static int     Hist     [20];       //!< Hist switches
  static double  DefView  [ 4];       //!< Default view field
  static double  dCut    [100];       //!< Misc. cuts (double)
  static double  dPRpar  [120];       //!< Pattern Recognition parameters (double)
  static double  Target   [ 3];       //!< Position (X,Y,Z) of the target center
  static bool    targetOpt;           //!< Set = true when Target is specified via option
  static double  MuonWall [ 20];      //!< Muon walls: coordinates of center and 1/2-sizes, and coordinates and 1/2-apertures of central hole
  static double  Calo     [ 40];      //!< Calos     : coordinates of center and 1/2-sizes, and coordinates and 1/2-apertures of central hole
  static double  RICHPipe [  6];      //!< RICH pipe : coordinates of center, 1/2-length, radius and material flag (flag>0, meaning heavy) 
  // Smoothing: "SmoothDet" takes precedence over "SmoothPos" array. It is
  // designed to be used in detector studies, where the detector under exam (or
  // group thereof) is expected to be the only smooth point.
  static double  SmoothPos[400];      //!< X abscissae to calculate smoothed track parameters
  static std::string  SmoothDet;      //!< Detector name @ which X abscissa to calculate smoothed track parameters
  static bool    ELossStraggling;     //!< Whether energy loss uncertainty should enter covariance matrix

  /*****************************Drell-Yan options*****************************/
  static bool DY_InAcceptance;        //!< When 0, the hit search in the vertex detector is performed for all tracks. When 1, the hit search in the vertex detector is restricted to those tracks which are in the Drell-Yan acceptance 
  static double DY_VD_Chi2;           //!< Upper limit of the chi2/ndf for fitted tracks with hits in the vertex detector. These hits (from the vertex detector) are removed if the chi2/ndf is greater than the upper limit.
  static double DY_VD_time;           //!< Upper limit of the hit time in the vertex detector (hits with bigger times are rejected)

  static std::string  PSdir;          //!< Detectory to write PS file
  static std::string  Dicofit;        //!< Input file for lattice Dico
  static std::string  DetNameOff[20]; //!< detectors' names excluded from the pattern recognition
  static std::string  Det2Go2Fit[20]; //!< detectors' names included in fit (default is all)
  static std::string  Det2Ignore[20]; //!< detectors' names to ignore completely
  static std::string  DZisActive[10]; //!< detectors' names with activated dead zone
  static std::string  GEMSpacers;     //!< File describing GEMs' spacers

  static int CAOpt[30];
  static double CAOptD[30];

  //! Virtual destructor
  // The declaration of this virtual destructor is meant to avoid the following compilation warning by gcc,v4.1: "has virtual functions but non-virtual destructor". (Btw, I don't know of any class derived from a TOpt base. Neither do I know the purpose of the "dum" method, declared hereafter as virtual.
  // By making a base class destructor virtual, one ensures that the destructor of any derived class is executed (in addition to, and prior to,  that of the base class.
  virtual ~TOpt() { std::cout<<"Destructor : Base TOpt"<<std::endl; }
  //! methods
  virtual void dum()=0;           //!< to make this class abstract
  static bool getOptions();       //!< get and store TraFFiC options

private:

  template <class T> 
  static bool getOptArray  (std::string key, T* arr, int siz, char m = 'm');
  template <class T> 
  static bool getOptStrings(std::string key, T* arr, int siz, char m = 'm');

}; 

#ifdef INIT_STATIC_HERE 

int    TOpt::Graph    [10];
int    TOpt::Print    [10];
int    TOpt::iCut     [40];
int    TOpt::iPRpar  [120];
int    TOpt::DetOff  [100];
int    TOpt::ReMode  [100];
int    TOpt::Hist     [20];
double TOpt::DefView  [ 4];
double TOpt::dCut    [100];
double TOpt::dPRpar  [120];
double TOpt::Target   [ 3];
bool   TOpt::targetOpt = false;
double TOpt::MuonWall [20];
double TOpt::Calo     [40];
double TOpt::RICHPipe [ 6];
double TOpt::SmoothPos[400];
std::string TOpt::SmoothDet;
bool   TOpt::ELossStraggling;

bool   TOpt::DY_InAcceptance;
double TOpt::DY_VD_Chi2;
double TOpt::DY_VD_time;

std::string TOpt::PSdir;
std::string TOpt::Dicofit;
std::string TOpt::GEMSpacers;
std::string TOpt::DetNameOff[20];
std::string TOpt::Det2Go2Fit[20];
std::string TOpt::Det2Ignore[20];
std::string TOpt::DZisActive[10];

int TOpt::CAOpt[30];
double TOpt::CAOptD[30];

#undef INIT_STATIC_HERE
#endif

#endif
