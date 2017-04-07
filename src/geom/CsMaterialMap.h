
/*!
   \file CsMaterialMap.h
   \brief Material Map Class.
   \author  Alexandre Korzenev, Jan P. Nassalski
   \version $Revision: 1.16 $
   \date    $Date: 2010/11/11 14:13:59 $
*/

#ifndef CsMaterialMap_h
#define CsMaterialMap_h

#include "THlx.h"
#include "CsErrLog.h"
#include <vector>
#include <string>

/*! \class CsMaterialMap
 *  \brief Material Map Class.
 *
 *  The class CsMaterialMap contains the table of radiative lengths
 *  and the table of dE/dX. The volume is split into boxes
 *  and an average radiative lengths or dE/dX is assigned to every box.
 *  The basic idea of introducing material maps in CORAL is to provide 
 *  a track reconstruction program with a way to take into account
 *  effects of multiple scattering and energy losses.
 */
class CsMaterialMap {

  friend class PaMaterialMaps;

 private:

  struct MatMap {
    int CoordSyst;         //!< if > 0 - Cartesian coordinates (X,Y,Z), if < 0 - cylindrical (PHI,R,Z).
    int nx,ny,nz;          //!< map dimensions.
    double minX,minY,minZ; //!< Low map edges.
    double maxX,maxY,maxZ; //!< Upper map edges.
    std::vector<double> x,y,z;  //!< vectors of coordinates.
    std::vector<double> Cos,Sin;//!< for cylindrical system.
    std::vector<double> Val;    //!< mean value of a measured quantity (Radiative length or dE/dX).
    
    MatMap();
    inline bool insideZ( double Z ); //!< returns \e true if inside map.

    /*! 
      \brief returns distance in direction \a direc to closest border of map.
      If there is not intersection the big number is returned.
    */ 
    inline double StepToBorder( double Z, bool direc ); 
    
  };
 public:
  CsMaterialMap();                             //!< Default Constructor

  ~CsMaterialMap();                            //!< Destructor

  /*! \fn bool ReadMaps( std::string gem_vers );
    \brief Read from files definited in option file. 

    The tag corresponded to maps in option file is CsMaterialMap. 
    Keys are Zone1, Zone2 and Zone3.
    \param gem_vers version of geometry
  */
  bool ReadMaps( std::string gem_vers );

  /*! \fn std::vector getZoneBorders();
    \brief Returns vector with zone borders in increasing order.
  */
  std::vector<float> getZoneBorders();

  /*! \fn int getNofMaps(void)
    \brief Returns number of involved maps.
  */
  int getNofMaps(void) {return RL_.size();};

  /*! \fn inline float getRadLength( float pos_x, float pos_y, float pos_z);
    \brief Returns the radiative length in a particular point in space.
     
    Units for all values are \e mm.
  */
  float getRadLength( float pos_x, float pos_y, float pos_z);

  /*! \fn inline void getRadLength( THlx &hel, bool direc, float &RadLen, float &StepZ);
    \brief Returns the radiative length \a RadLen in the point defined 
           by \a hel. Also returns the step size  along \e Z axis 
	   from \a hel to next medium.
    \param hel    The position in space (must be filled according to TRAFFIC).
    \param direc  Direction along beam axis.
     \arg \a direc = \a true  : extrapolation is performed along beam.
     \arg \a direc = \a false : extrapolation is performed in opposite to beam derection.
    \param RadLen The radiative length in \a hel (in \e cm).
    \param StepZ  The step size along beam axis from \a hel to next medium (in \e cm). 
                  Track extrapolated as straight line.
		  If extrapolation is in backward direction than \a StepZ \a < \a 0.
  */
  void getRadLength( const THlx &hel, bool direc, float &RadLen, float &StepZ);

  /*! \fn inline void getRadLength( CsHelix &hel, bool direc, float &RadLen, float &StepZ);
    \brief Returns the radiative length \a RadLen in the point defined 
           by \a hel. Also returns the step size  along \e Z axis 
	   from \a hel to next medium.
    \param hel    The position in space (must be filled according to CORAL).
    \param direc  Direction along beam axis.
     \arg \a direc = \a true  : extrapolation is performed along beam.
     \arg \a direc = \a false : extrapolation is performed in opposite to beam derection.
    \param RadLen The radiative length in \a hel (in \e mm).
    \param StepZ  The step size along beam axis from \a hel to next medium (in \e mm). 
                  Track extrapolated as straight line.
		  If extrapolation is in backward direction than \a StepZ \a < \a 0.
  */
  void getRadLength( const CsHelix &hel, bool direc, float &RadLen, float &StepZ);

  /*! \fn float getdE( const THlx& hel, float Len );
    \brief Returns the most probable energy loss ( [dE] = GeV ) for momentum losses during a
           step of length Len at the point defined by \a hel in target.
    \param hel    The position in space (must be filled according to TRAFFIC).
    \param Len    length of the step.
  */
  float getdE( const THlx& hel, float Len);

  // for backward compatibility with "MatMap" code 
  float getdEdX( const THlx& hel) {
    float Len = 1.;
    return(getdE(hel,Len));
  };

  /*! \fn float getdEStraggling( const THlx& hel, float Len );
    \brief Returns the sigma of energy loss distribution ( [dE] = GeV ) for momentum losses during a
           step of length Len at the point defined by \a hel in target.
    \param hel    The position in space (must be filled according to TRAFFIC).
    \param Len    length of the step.
  */
  float getdEStraggling(const THlx& hel, float Len );

  /*! \fn bool usingROOTGeometry();
    \brief Return true if the ROOTGeometry is in use, false otherwise.
  */
  bool usingROOTGeometry() const { return (usingROOTGeometry_ 
#if USE_TGEANT
    || usingROOTGeometryTGEANT_
#endif
  ); }
  
#if USE_TGEANT  
  /*! \fn bool usingROOTGeometryTGEANT();
    \brief Return true if the ROOTGeometry of TGEANT is in use, false otherwise.
  */
  bool usingROOTGeometryTGEANT() const { return usingROOTGeometryTGEANT_; }
  
  /*! \fn std::vector<char> getGDML();
   * \brief returns the whole GDML-File as vector of char - contains either binary root tree or gdml text file
   */
  std::vector<char> getGDML() { return GDMLFile; } 
#endif
  
  /*! \fn TMacro* getROOTGeometry();
    \brief Returns a pointer to the macro containing the currently loaded ROOT Geometry.
  */
  class TMacro* getROOTGeometry() { return ROOTGeometry_; }

  /*! \fn double getMassDefault();
    \brief Returns the default mass used for energy loss calculations.
  */
  double getMassDefault() { return massDefault_; }

  int SimpleELoss() const {return simpleELoss_;}
  
 private:
  bool usingROOTGeometry_;
  
#if USE_TGEANT
  bool usingROOTGeometryTGEANT_;
  std::vector<char> GDMLFile;
#endif
  
  TMacro* ROOTGeometry_;
  double massDefault_;  // default mass used for energy loss calculation
  int    simpleELoss_;  // if true, Tobi's simple energy loss calculation is done, otherwise
                        // Catarina's more detailed calculation is used 

  std::vector<MatMap*> RL_; //!< vector of pointers to radiative length maps.
  MatMap* EL_;         //!< pointer to map of energy losses in target.
 
  /*! \fn void getRadLengthZone( MatMap *map, THlx hel, bool direc, float &RadLen, float &StepZ);
    \brief Returns the radiative length \a RadLen in the point defined 
           by \a hel in \a map. Also returns the step size along \e Z axis
	   from \a hel to next medium.
    \param hel    The position in space (must be filled according to TRAFFIC).
    \param direc  Direction along beam axis.
      \arg \a direc = \a true  : extrapolation is performed along beam.
      \arg \a direc = \a false : extrapolation is performed in opposite to beam derection.
    \param RadLen The radiative length in \a hel (in \e cm).
    \param StepZ  The step size along beam axis from \a hel to next medium (in \e cm). 
                  Track extrapolated as straight line.
  */
  void getRadLengthZone( MatMap *map, THlx hel, bool direc, float &RadLen, float &StepZ);
  
  /*! \fn float getStepSize(THlx &hel,MatMap *map , int i, int j, int k , float RadLen);
    \brief Returns the distance from \a hel to next medium. 
           Track extrapolated as straight line.
  */
  float getStepSize(const THlx &hel, MatMap *map, int i, int j, int k , float RadLen);

};

#endif // CsMaterialMap_h
























