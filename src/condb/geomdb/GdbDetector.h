//-*-Mode: C++;-*-
#ifndef _GdbDetector_h_
#define _GdbDetector_h_
/*!
  \file		GdbDetector.h
  \brief	Geometry DB detector object definition file
  \author	$Author: zvyagin $
  \version	$Revision: 1.3 $
  \date		$Date: 2000/07/17 13:25:35 $
*/

#if __GNUG__ >= 2
#  pragma interface
#endif

#include "GdbSTD.h"
#include "GdbStation.h"
#include "GdbDetTableCont.h"

//---- GdbDetector -----------------------------------------------------------
/*!
	\class	GdbDetector
	\brief	Detector object for GDB package.
			This consists of a vector of GdbStation 
			which will be created from COMGEANT files or Geometry DB.
*/
class GdbDetector {
public:
  GdbDetector();	//!< default constructor
/*!
  \fn		GdbDetector( const int& id, const char* pName, const int& unit)
  \param	id is a int
  \param	pName is a pointer to char
  \param	unit is a int
  \brief	constructor with id, name and unit of Detector.
*/
  GdbDetector( const int& id, const char* pName, const int& unit);	

//! copy constructor
  GdbDetector( const GdbDetector& detector );

  ~GdbDetector();	//!< destructor

//! equal to operator
  bool operator==( const GdbDetector& detector ){
	  return this->id() == detector.id() ;
  }

//! less than operator
  bool operator<( const GdbDetector& detector ){
	  return this->id() < detector.id() ;
  }
  
//! dump information to ostream
  friend ostream& operator<<(ostream& os, GdbDetector& detector);

  void id(const int& id) {id_ = id;}	//!< set id
  const char* name( const char* pName );	//!< set detector name
  void type( const GdbDetectorType& type ) {type_ = type;}//!< set type id
  int id() const {return id_;}					//!< return id
  const char* name() const {return name_;}		//!< return name
  GdbDetectorType type() const {return type_;}	//!< return type

  int nStation() {return station_.size();}		//!< return number of station
  typedef list<GdbStation>::iterator itrStation;	//!< usefull type definition


  vector<GdbStation>& station() {return station_;}	//!< return reference to station list
//! return a copy of station list
  vector<GdbStation> stationCopy() const {return station_;}

//! return index-th entry of GdbStation container
  GdbStation& station(const int& index) {return station_[index];}


/*!
  \struct	PlaneId
  \brief	PlaneId consists of name and unit as its attributes.
*/
  struct PlaneId {

  //! default constructor
	PlaneId() : name(""), unit(0) {}
  //! constructor with variables
	PlaneId(const string& n, const int& u) :
		name(n), unit(u) {}
  //! copy constructor
	PlaneId(const PlaneId& plane) :
		name(plane.name), unit(plane.unit){}			
  //! equal to operator
	bool operator==(const PlaneId& plane) const {
		return ( name == plane.name )&&( unit == plane.unit ); 
	}
  //! less than operator
	bool operator<(const PlaneId& plane) const {
		if( name < plane.name ) return true;
		return unit < plane.unit;
	}
	string	name;	//!< name
	int 	unit;	//!< unit
  };


private:
  static const int maxNameLength = 8;		//!< maximun length of name
  char name_[maxNameLength];			//!< detector name
  int id_;								//!< detecotor id (used in FFR card)
  vector<GdbStation> station_;			//!< Station Information
  GdbDetectorType type_;				//!< detector type
  map<PlaneId, int> idMap_;				//!< id map


public:
  //! return list of Active Plane
  list<GdbDetTableCont> table() const {return table_;}

  //! return list of Active Plane
  list<GdbDetTableCont>& refTable() {return table_;}

  //! set list of Active Plane
  void refTable(const list<GdbDetTableCont>& table ) {
	  table_ = table;
  }

  map<PlaneId, int> idMap() const {return idMap_;} //!< return id map


  bool createDetectorTable();			//!< create detector table
  void dumpDetectorTable(ostream& os);	//!< dump detector table to ostream
  void dump(ostream& os);				//!< dump information to ostream
  void dumpVolumeTree(ostream& os);		//!< dump volume tree to ostream
  void dumpXML(ostream& os);			//!< dump in XML format
  void createDetectorFile();			//!< make output file, name.gdat.
	
private:
  list<GdbDetTableCont> table_;			//!< active plane list
  bool createIdMap();					//!< creat id map
  int findIdInMap(const PlaneId& pid);	//!< find id in id map

};

#endif
