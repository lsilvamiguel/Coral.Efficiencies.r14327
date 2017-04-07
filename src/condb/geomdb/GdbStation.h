#ifndef GdbStation_h
#define GdbStation_h
/*!
  \file		GdbStation.h
  \brief	Geometry DB transient station object definition file
  \author	$Author: miyachi $
  \version	$Revision: 1.2 $
  \date		$Date: 2000/06/07 08:36:48 $
*/

#include <CLHEP/Geometry/Point3D.h>
#include "GdbDetectorType.h"
#include "GdbVolume.h"
#include "GdbPerStation.h"

class GdbStation : public GdbVolume {

public:

  GdbStation();	//!< defalt constructor

  GdbStation(const GdbStation& station); //!< copy constructor

  GdbStation(GdbPerStation& perStation); //!< copy from persistent object

  virtual ~GdbStation();	//!< destructor

  GdbStation& operator=(const GdbStation& station); //!< assign operator

/*!
  \fn		bool findActivePlane( list<ActivePlane>& foundPlane );
  \param	foundPlane is a list<ActivePlane>
  \return	a bool
  \brief	try to find all active plane in this tree, 
			and store the found plane in the list.
			If the method finish successfully, this returns true.
*/
  bool findActivePlane( list<ActivePlane>& foundPlane );


  list<HepPoint3D> references() const {return references_;} //!< retturn reference points

  long nReference() const {return references_.size();}	//!< return number of reference points

//! set referene points
  void references(const list<HepPoint3D>& references) {
	  references_ = references;
  }	

//! add reference points
  void addReference(const HepPoint3D& ref) {
	  references_.push_back(ref);
  }	
  
//! return detector type
  GdbDetectorType detectorType() const {return detectorType_;}

//!	set detector type
  void detectorType(const GdbDetectorType& detectorType){
	  detectorType_ = detectorType;
  }

  friend ostream& operator<<(ostream& os, GdbStation& station); //!< output infromation

  void dump(ostream& os, const char* pTab = NULL); //!< dump information.

  void dumpXML(ostream& os, const char* pTab = NULL); //!< dump information.


  void dumpDetectorsTable(ostream& os); //!< dump detectors.dat table

  string gVersion() const {return gVersion_;} //!< return geometry version tag
  void gVersion(const string& version) {gVersion_ = version;} //!< set version

//! cast operator to GdbPerStation
  operator GdbPerStation ();


//! return i-th volume
  GdbVolume& volume(const long& i) {return volumes_[i];}
//! return a reference to vector of GdbVolume
  vector<GdbVolume>& refVolumes() {return volumes_;}
//! return a copy of GdbVolume vector
  vector<GdbVolume> volumes() const {return volumes_;}

//! set volume vector
  void volumes(const vector<GdbVolume> volumes) {volumes_ = volumes;}

//! add volume 
  void addVolume( const GdbVolume& volume ){volumes_.push_back(volume);}
//!< get number of volumes in this station
  long nVolume() const {return volumes_.size();}

//! return a copy of volume map
  map<long, vector<long> > volumeMap() const {return volumeMap_;}
//! return a reference to volume map
  map<long, vector<long> >& refVolumeMap()  {return volumeMap_;}
//! create volume map
  void  createVolumeMap( long& index );

private:
  string gVersion_;						//!< geometry versioning tag
  GdbDetectorType detectorType_;		//!< detector type
  list<HepPoint3D> references_;			//!< reference points
  vector<GdbVolume> volumes_;			//!< volume list 

  map<long, vector<long> > volumeMap_;	//!< volume tree map

//! dump volume information to ostream
  void dumpVolume(ostream& os, const char* pTab = NULL);
//! dump index-th volume information to ostream
  void dumpVolume(ostream& os, const long& index, const char* pTab = NULL);
//! dump index-th volume information to ostream in XML format
  void volumeXML(ostream& os, const long& index, const char* pTab);

/*!
  \fn		bool findActivePlane( long& index, list<ActivePlane>& foundPlane )
  \param	index is a long
  \param	foundPlane is a list<ActivePlane>
  \return	a bool
  \brief	try to find all active plane under index-th volume
			in this tree, 
			and store the found plane in the list.
			If the method finish successfully, this returns true.
*/
  bool findActivePlane( long& index, list<ActivePlane>& foundPlane );
  
};

#endif
