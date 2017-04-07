#ifndef _Gspos_h_
#define _Gspos_h_
/*!
  \file		Gspos.h
  \brief	Object for GSPOS entry in Geant 3.
  \author	$Author: ybedfer $
  \version	$Revisopn$
  \date		$Date: 2009/09/14 00:51:22 $
*/

#if __GNUG__ >= 2
#  pragma interface
#endif

#include <CsSTD.h>
#include <CLHEP/Geometry/Vector3D.h>
#include <CLHEP/Matrix/Matrix.h>

#include "Gsdet.h"
#include "Gsvolu.h"
//---- Gspos -----------------------------------------------------------
/*!
  \class	Gspos
  \brief	Object which handle GSPOS entry in g3call list.
*/
class Gspos {
public:

  Gspos(); //!< defalt constructor

/*!
  \fn		Gspos(const string& name, const unsigned int&  unit);
  \param	name is a string
  \param	unit is a unsigned int
  \brief	Constructor with name and unit number.
*/
  Gspos(const std::string& name, const unsigned int&  unit);

/*!
  \fn  	Gspos(	const string& name, const unsigned int&  unit,
			const string& mother, const HepVector3D  center,
			const int& rotation, const string& flag );
  \param	name is a string
  \param	unit is a unsigned int
  \param	mother is a string
  \param	center is a HepVector3D
  \param	rotation is a int
  \param	flag is a string
  \brief	constructor
*/
  Gspos(	const std::string& name, const unsigned int&  unit,
			const std::string& mother, const HepGeom::Vector3D<double> center,
			const int& rotation, const std::string& flag );

  Gspos(const Gspos& pos); //!< copy constructor

  ~Gspos();	//!< destructor

//! equal to operator
  bool operator==(const Gspos& pos) const {
	return (	this->name() == pos.name() && this->unit() == pos.unit() );
  }

//! less than opeartor
  bool operator<(const Gspos& volume) const {
	return true;
  }

/*!
  \fn		void addToDaughter(Gspos* pPos);
  \param	pPos is a pointer to Gspos
  \brief	add the given Gspos in this daughter volume list.
*/
  void addToDaughter(Gspos* pPos);
	
	
/*!
  \fn		Gspos& findPosition( Gspos& pos, bool& result );
  \param	pos is a Gspos
  \param	result is a bool
  \return	a Gspos
  \brief	return a reference to Gspos found in the daughter list 
			which is equal to the given pos.
			The serch result is stored in result.	
*/
  Gspos& findPosition( Gspos& pos, bool& result );

/*!
  \fn		Gspos& findPosition( 
			const string& name, const unsigned int& unit, bool& result );
  \param	pos is a Gspos
  \param	unit is a unsigned int
  \param	result is a bool
  \return	a Gspos
  \brief	return a reference to Gspos found in the daughter list 
			which has equal name and unit to the given.
			The serch result is stored in result.
*/
  Gspos& findPosition( const std::string& name, const unsigned int& unit,
						  bool& result );

private:
  int id_;				//!< position id
  std::string name_;			//!< name
  unsigned int unit_;			//!< unit number
  std::string mother_;			//!< mother volume name
  HepGeom::Vector3D<double> center_;	//!< center position in Mother Reference System(MRS)
  int rotation_;			//!< rotation matrix id
  std::string flag_;			//!< flag (ONLY, MANY)
  unsigned int n_;			//!< duplicant number used by CsG3CallFile
  unsigned int motherN_;		//!< duplicant number of mother volume
  Gspos* pMother_;			//!< pointer to Mother volume
  CLHEP::HepMatrix* pRotation_;		//!< pointer to Matrix Rotation
  CLHEP::HepMatrix rM_;			//!< rotaiton matrix in MRS
  Gsdet*	pDetector_;		//!< pointer to detector information
  std::vector<double> param_;		//!< parametrized volume parameters 
  Gsvolu* pVolume_;			//!< pointer to volume information
  std::list<Gspos*> pDaughter_;		//!< list of pointer to daughter volume
  bool hasDetectorInTree_; 		//!< true if this volume tree has detector

public:
//------------  Normal set & get functions  ----------------//

//! return name
  std::string name() const {return name_;}
//! set name
  void name(const std::string& name) {name_ = name;}
//! return mother name
  std::string mother() const {return mother_;}
//! set mother name
  void mother(const std::string& mother) {mother_ = mother;}
//! return unit
  unsigned int unit() const {return unit_;}
//! set unit
  void unit(const unsigned int& unit ) {unit_ = unit;}
//! return rotation id
  int rotation() const {return rotation_;}
//! set rotation id
  void rotation(const int& rotation) {rotation_ = rotation;}
//! return cetner in MRS
  HepGeom::Vector3D<double> center() const {return center_;}	
//! set center in MRS
  void	center(const HepGeom::Vector3D<double>& center) {center_ = center;}
//! return flag
  std::string flag() const {return flag_;}
//! set flag
  void flag(const std::string& flag) {flag_ = flag;}
//! return id
  int id() const {return id_;}
//! set id
  void id( const int& id) {id_ = id;}
//! return n
  unsigned int n()	const {return n_;}
//! set n
  void n(const unsigned int& n) {n_ = n;}
//! return mother number
  unsigned int motherN() const {return motherN_;}
//! set mother number
  void motherN(const unsigned int& motherN) {motherN_ = motherN;}
//! return pointer to mother
  Gspos* pMother() const {return pMother_;}
//! set pointer to mother
  void pMother(Gspos* pMother) {pMother_ = pMother;}
//! return pointer to rotation matrix
  CLHEP::HepMatrix* pRotation() const {return pRotation_;}
//! set pointer to rotation
  void pRotation(CLHEP::HepMatrix* pRotation) {pRotation_ = pRotation;}
//! return rotation matrix to MRS
  CLHEP::HepMatrix	rM() const {return rM_;}
//! set rotation matrix to MRS
  void rM(const CLHEP::HepMatrix rM) {rM_ = rM;}
//! return pointer to detector information
  Gsdet* pDetector() const {return pDetector_;}
//! set pointer to detector information
  void	pDetector(Gsdet* pDetector) {pDetector_ = pDetector;}
//! return pointer to volume information
  Gsvolu* pVolume() const {return pVolume_;}
//! set pointer to volume information
  void pVolume( Gsvolu* pVolume) {pVolume_ = pVolume;}
//! return a copy of volume parameter
  std::vector<double> param() const {return param_;}
//! add one number in parameter list
  void addParam(const double& param) {param_.push_back(param);}
//! return list of pointers to daugheter
  std::list<Gspos*> pDaughter() const { return pDaughter_;}
//! set pointer list to daughter
  void 	pDaughter(std::list<Gspos*>& pDaughter) {pDaughter_ = pDaughter;}
//! return if this has detector somewhere this tree
  bool hasDetectorInTree() const {return hasDetectorInTree_;}
//! set hasDetectorInTree_ flag 
  bool hasDetectorInTree(const bool& hasDetectorInTree) {
	  hasDetectorInTree_ = hasDetectorInTree;
	  return hasDetectorInTree_;
  }
  
/*!
  \fn		bool checkDetector();
  \return	a bool
  \brief	check if this volume tree has a detector in the tree
			under this position. 
			If there are detectors (sensitive volume), 
			set hasDetectorInTree_ flag.
			After all check return the result.
*/
  bool checkDetector();

/*!
  \fn		ostream& operator<<(ostream& stream, Gspos& pos);
  \param	stream is a ostream
  \param	pos is a Gspos
  \return	ostream
  \brief	dump information on this to stream
*/  
  friend std::ostream& operator<<(std::ostream& stream, Gspos& pos);

//--- dump functions
/*!
  \fn		void dumpMother( ostream& stream ) const ;
  \param	stream is a ostream
  \brief	dump mother information to stream
*/
  void dumpMother( std::ostream& stream ) const ;

/*!
  \fn		void dumpDaughters(ostream& os, const string tab="");
  \param	os is a ostream
  \param	tab is a string
  \brief	dump daughter information to the gien os.	
*/
  void dumpDaughters(std::ostream& os, const std::string tab="");

/*!
  \fn		void dumpDetectors(ostream& os, const string tab="");
  \param	os is a ostream
  \param	tab is a string
  \brief	dump detector information to the given os
*/
  void dumpDetectors(std::ostream& os, const std::string tab="");

protected:
//! return a reference to list of refPDaughter
  std::list<Gspos*>& refPDaughter() {return pDaughter_;}

};

#endif
