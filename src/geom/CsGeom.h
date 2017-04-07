// $Id: CsGeom.h 14070 2015-09-20 01:36:58Z lsilva $

/*!
   \file    CsGeom.h
   \brief   Compass Geometry Interface Class.
   \author  Benigno Gobbo
   \version $Revision: 14070 $
   \date    $Date: 2015-09-20 03:36:58 +0200 (Sun, 20 Sep 2015) $

   \par History:
   19990803  Extracted this class from CsGeant3 one.<br>
   19991119  Added CsZone objects creation. <br>
   20000512  Moved to COMGEANT 5. Backward compatibility is kept. <br>
   20010404  Follow changes in detectors.dat starting from CG 6.03.  <br>
*/

#ifndef CsGeom_h
#define CsGeom_h

#include "coral_config.h"

#include <vector>
#include <list>
#include <string>
#include <map>

#include "CsMaterialMap.h"
#include "CsZone.h"
#include "CsDetector.h"
#include "CsRICH1Detector.h"
#include "CsMiscDetector.h"
#include "CsCEDARDetector.h"
#include "CsTime.h"
#include "CsField.h"
#include "CsBeamReconstructionSupport.h"

#include <CLHEP/Matrix/Matrix.h>


/*! \class CsGeom
  \brief Geant Interface Class.

  This class reads Monte the Comgeant generated detector table and
  istantiates the CsDetector, CsField and CsMaterialMap objects.
*/

class CsCalorimeter;
class CsBeamBackPropagationCoeff;

struct CsTargetCell {      //!< Target cell: coord's in MRS, and sizes
  double x[3], radius, length;
};

class CsGeom {

public:

  /*! \fn static CsGeom* Instance( const std::string* detTableName );
    \brief First singleton instantiation.
    \param detTableName Detector Table file name.
  */
  static CsGeom* Instance( const std::string* detTableName );

  /*! \fn static CsGeom* Instance( const std::string* detTableName, const std::string* pitchTableName );
    \brief First singleton instantiation.
    \param detTableName Detector Table file name.
    \param pitchTableName variable-sized Pitch Table file name.
  */
  static CsGeom* Instance( const std::string* detTableName, const std::string* pitchTableName );

  /*! \fn static CsGeom* Instance();
    \brief singleton instantiation (but first).
  */
  static CsGeom* Instance();

  /*! \fn inline std::list <CsDetector*> getDetectors();
    \brief Returns the list of pointers to CsDetector objects
  */
  const std::list<CsDetector*>& getDetectors() const { return( _dets  ); }

  /*! \fn inline CsRICH1Detector* getRich1Detector();
    \brief Returns the pointer to the CsRICH1Detecor object
  */
  inline CsRICH1Detector* getRich1Detector() { return( _rich1 ); }

  /*! \fn inline CsRICH1Detector* getRich1Detector();
    \brief Returns the pointer to the CsRICH1Detecor object
  */
  inline std::vector< CsCalorimeter*> &getCalorimeters() { return( _calorimeters ); }

  inline  std::vector<CsCEDARDetector*> &getCEDARs() { return( _cedar_dets ); }

  //! Returns the list of pointers to CsMiscDetector objects
  std::list<CsMiscDetector*> getMiscDetectors() { return( _others ); }

  /*! \fn inline std::list <CsZone*> getZones();
    \brief Returns the list of pointers to CsZone objects
  */
  const std::list<CsZone*>& getZones() const { return( _zones  ); }

  /*! \fn CsField getCsField();
    \brief Returns the pointer to the CsField object
  */
  CsField* getCsField() { return( &_field ); }

  /*! \fn inline double getTargetCenter();
    \brief Returns the target center position along the beam _targetCenter
  */
  inline double getTargetCenter( ) { return _targetCenter; }

  /*! \fn const std::vector<CsTargetCell> getTargetCells();
    \brief Returns the target substructure
  */
  const std::vector<CsTargetCell> &getTargetCells() const {
    return _targetCells;
  }

  /*! \fn const std::vector<CsBeamBackPropagationCoeff*> &getBeamBackProp(); 
    \brief Return beam back propagation coefficients stored in alignment file
  */
  const std::vector<CsBeamBackPropagationCoeff> &getBeamBackProp() const {
    return _bbp;
  }

  /*! \fn std::string getComgeantVers();
    \brief Returns the Comgeant version string
  */
  std::string getComgeantVers() { return( _comgeantVersion ); }

  //! Returns the Comgeant geometry version string
  std::string getGeomVers();

  //! Returns true if Comgeant hadron setup
  bool isHadronSetup();

  //! Returns true if Comgeant muon setup
  bool isMuonSetup();

  /*! \fn inline CsMaterialMap getCsMaterialMap();
    \brief Returns the pointer to the CsMaterialMap object
  */
  CsMaterialMap* getCsMaterialMap();

 protected:

  CsGeom( const std::string* detTableName );  //!< The Constructor
  CsGeom( const std::string* detTableName, const std::string* pitchTableName );  //!< The Constructor

 private:

  static              CsGeom* _instance;  //!< The Singleton Static pointer
  std::string              _detTableName;      //!< detector table file name
  std::string              _pitchTableName;      //!< pitch table file name
  std::string              _geomVersion;       //!< Comgeant geometry version
  std::string              _comgeantVersion;   //!< Comgeant version
  bool                _muonSetup;         //!< True if muon setup (CG > 6.03)
  bool                _hadronSetup;       //!< True if hadron setup (CG > 6.03)

  int                 _orient[3];      //!< Orientation.
  CsField             _field;          //!< Field Class Object
  CsMaterialMap       _matmap;         //!< Class of table with average radiative lengths.
  std::vector<CLHEP::HepMatrix>      _matrices;       //!< List of rotation matrices
  std::list <CsDetector*>     _dets;           //!< List of Detector Objects
  CsRICH1Detector*       _rich1;          //!< Pointer to Rich1 Detector
  std::list <CsMiscDetector*> _others;         //!< List of "Misc" Detector Objects
  std::list <CsZone*>         _zones;          //!< List of Zones
  std::vector<CsCalorimeter*>      _calorimeters;       //!< List of calorimeters

  std::vector<CsBeamBackPropagationCoeff> _bbp;   //!< Beam Back Propagation Objects

  std::vector<CsCEDARDetector*> _cedar_dets;

  /*! \var double _targetCenter
    \brief Rudimentary but easily accessible (as opposed to material map's) info
    Mainly used (as of 2006/10) by vertexing to discriminate between beam and spectro tracks.
  */
  double _targetCenter;       //!< Target center position along the beam (in mm)
  /*! \var vector<CsTargetCell> _targetCells
    \brief Rudimentary but easily accessible (as opposed to material map's) info
    Mainly used (as of 2006/10) for graphics.
  */
  std::vector<CsTargetCell> _targetCells;    //!< Target's substructure

  std::map<int,double>   _PitchCorr;  //!< Table of variable size pitches

  /*! \fn bool  readDetTable( const std::string* detTableName );
     \brief Reads the detector table file
     \param detTableName The detector table file name
  */
  bool  readDetTable( const std::string* detTableName );

  /*! \fn bool  readPitchTable( const std::string* pitchTableName );
     \brief Reads the detector's variable pitch correction file
     \param pitchTableName The pitch table file name
  */
  bool  readPitchTable( const std::string* pitchTableName );

  /*! \fn bool makeZones()
     \brief Reads from option file the zones definitions (if any) and
     creates CsZone objects. If no zones in option file a single big
     zone (containing al Compass) is in any case create.
  */
  bool makeZones();

  /*! \fn void setDetsZones()
    \brief Sets for each detector the zones from where is belongs from
  */
  void setDetsZones();

  /*! \fn void  GeaRef2CsRefVec( double& x, double& y, double& z );
     \brief Rotates the axes from the Geant Reference System ( X along
     beam, Z vertical) to the Compass Reference System ( Z along beam,
     Y vertical).
     \param x x coordinate. WARNING: its value will be modified!
     \param y y coordinate. WARNING: its value will be modified!
     \param z z coordinate. WARNING: its value will be modified!
  */
  void  GeaRef2CsRefVec( double& x, double& y, double& z );

  /*! \fn void  void  GeaRef2CsRefMat( HepMatrix& a );
     \brief Rotates matrices from the Geant Reference System ( X along
     beam, Z vertical) to the Compass Reference System ( Z along beam,
     Y vertical)
     \param a the matrice to be rotated. WARNING: it will be modified!
  */
  void  GeaRef2CsRefMat( CLHEP::HepMatrix& a );

  /*! \fn int Gpar2PDGpar( int Gpar );
     \brief Converts from Geant Particles numbering to PDG particles
     numbering convention.

     NOTE: (Hope) a temporary method.
  */
  int Gpar2PDGpar( int Gpar );

  //! used for sorting list of pointers to detectors
  struct sortDetectors_;

  /*! \fn bool associateDets();
    \brief perform detector association (DC shifted planes),
    using coral.option
    return true, if everything is OK
  */
  bool associateDets();

  /*! \fn void associateGEMs(int mode);
   \brief Performs GEM associations, viz. X<->Y and V<->U (mode&0x1) and CsPixelGEM->CsGEM associations (mode&0x2)
  */
  void associateGEMs(int mode) ;

};

#endif // CsGeom_h
