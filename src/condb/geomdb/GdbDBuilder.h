//-*-Mode: C++;-*-
#ifndef _GdbDBuilder_h_
#define _GdbDBuilder_h_
/*!
  \file		GdbDBuilder.h
  \brief	Geometry Data builder object definition file
  \author	$Author: miyachi $
  \version	$Revision: 1.2 $
  \date		$Date: 2000/06/07 08:36:44 $
*/

#if __GNUG__ >= 2
#  pragma interface
#endif

#include "CsSTD.h"

#include "GdbDetector.h"
#include "Gspos.h"
#include "CsEndOfJob.h"

#include "CsRotation.h"
#include <CLHEP/Geometry/Point3D.h>


//---- GdbDBuilder -----------------------------------------------------------
/*!
  \class	GdbDBuilder
  \brief	GDB Detector Information builder class.
			Using GdbDBuilder one can create geometry data from 
			COMGEANT files or Geometry DB.
			The data constructed in this object will be used by coral 
			as transient object or stored in DB by DB updater.
			GdbDBuilder is a singlton object one must use Instance() method
			to access this object.
*/
class GdbDBuilder : public CsEndOfJob {
public:

/*!
  \enum		MODE
  \brief	Used as a flag to indicate the present system of geometry data.
*/
  enum MODE {
	GEODB	= 1,
	G3FILE	= 2
  };

/*!
  \fn		static GdbDBuilder* Instance();
  \return	a pointer to GdbDBuilder
  \brief	GdbDBuilder object instanciation method. 
*/
  static GdbDBuilder* Instance();

/*!
  \fn		static GdbDBuilder* Instance( const GdbDBuilder::MODE& mode = GEODB );
  \param	mode is a GdbDBuilder::MODE which has a default value of GEODB
*/
  static GdbDBuilder* Instance( const GdbDBuilder::MODE& mode = GEODB );

//! end method called at the end of main program
  bool end();
	
private:
  static GdbDBuilder* instance_;	//!< singlton pointer
  
/*!
  \fn		explicit GdbDBuilder( const GdbDBuilder::MODE& mode = GEODB )
  \param	GEODB is a GdbDBuilder::MODE which has a default value of GEODB
  \brief	constructor with system morde.
*/
  explicit GdbDBuilder( const GdbDBuilder::MODE& mode = GEODB );


public:
  ~GdbDBuilder();	//!< destructor

  bool build();	//!< build detector

  friend ostream& operator<<(ostream& os, GdbDBuilder& builder); //!< dump information

//--------------------- useful functions.
/*!
  \fn		static GdbPoint3D convert(const HepPoint3D& point);
  \param	point is a HepPoint3D
  \return	GdbPoint3D
  \brief	convert HepPoint3D object to GdbPoint3D object and return.
*/
  static GdbPoint3D convert(const HepPoint3D& point);

/*!
  \fn		static HepPoint3D convert(const GdbPoint3D& point);
  \param	point is a GdbPoint3D
  \return	HepPoint3D
  \brief	convert GdbPoint3D object to HepPoint3D object and return.
*/
  static HepPoint3D convert(const GdbPoint3D& point);

/*!
  \fn		static GdbRotation convert(const CsRotation& rotM);
  \param	rotM is a CsRotation
  \return	GdbRotation
  \brief	convert CsRotation object to GdbRotation object and return.
*/
  static GdbRotation convert(const CsRotation& rotM);

/*!
  \fn		static CsRotation convert(const GdbRotation& rotM);
  \param	rotM is a GdbRotation
  \return	CsRotation
  \brief	convert GdbRotation object to CsRotation object and return.
*/
  static CsRotation convert(const GdbRotation& rotM);

  list<GdbDetector>& detectors() {return detectors_;}	//!< return detector list
  GdbDetector& add( const GdbDetector& detector );	//!< add detector in list

  MODE mode() const {return mode_;}	//!< return present mode

/*!
  \fn		bool isUpdated() const 
  \return	a bool
  \brief	return true if Geometry data in this is updated and ready to use.
*/
  bool isUpdated() const {return updated_;}
  void isUpdated(const bool& updated) {updated_ = updated;} //!< set update flag

  void makeDetectorFiles();	//!< dump detector information to independent file

private:
  MODE mode_;	//!< system mode

  list<GdbDetector> detectors_;	//!< list of GdbDetector	

  bool updated_;	//!< true = ready for use

  bool detector();	//!< build ID-th detector
  bool getDetectorFromFFR( int const& id ); //!< read ffr file
  bool fillDetectorByG3Call( GdbDetector& detector ); //!< fill detector pName

  bool createVolumesIn( GdbStation& station, Gspos& pos ); //!< create subvolumes

  bool rich();	//!< create RICH volume
  bool muonFilter();	//!< create muon filtter
  string findMF2Station( const int& id ); //!< used by muonFilter method

//! create one passive volume object
  bool passiveVolume( const string& detectorName, const string& volumeTop );

//! create passive volumes
  bool CreatePassiveVolumes();

//! set volume information using pos
  void setVolumeInformation( GdbVolume& volume, Gspos& pos );

protected:

  typedef list<GdbDetector>::iterator itrDetector;	//!< iterator
	
	
};

#endif
