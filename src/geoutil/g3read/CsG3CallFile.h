#ifndef _CsG3CallFile_h_
#define _CsG3CallFile_h_
/*!
  \file		CsG3CallFile.h
  \brief	A header file for a class which convertes geometry information 
			from GEANT 3 call list to CORAL enviroment.
  \author	$Author: tnagel $
  \version	$Revision: 1.6 $
  \date		$Date: 2010/02/02 13:20:41 $
*/

#if __GNUG__ >= 2
#  pragma interface
#endif

#include "CsSTD.h"
#include "CsEndOfJob.h"
#include "CsRegistrySing.h"

#include <CLHEP/Matrix/Matrix.h>

#include "CsGeomUtils.h"
#include "CsMaterial.h"
#include "Gstmed.h"
#include "Gsvolu.h"
#include "Gspos.h"
#include "Gsdet.h"

//---- CsG3CallFile -----------------------------------------------------------
/*!
  \class	CsG3CallFile
  \brief	This convertes the COMPASS geometry information 
			including material, medium, rotation matricses, volume, position,
			and detectors information into CORAL.
*/
class CsG3CallFile : public CsEndOfJob  {
public:
//--------------------  For singleton setup
  static CsG3CallFile* Instance(); //!< instance method

  static CsG3CallFile* instance_;	//!< static pointer to this

  bool end();	//!< for end of job


private:
//! default constructor
  CsG3CallFile();

public:
//!	destructor
  virtual ~CsG3CallFile();

/*!
  \fn		bool readCallList();
  \return	a boolen
  \brief	open GEANT 3 call list (filename), and store all information 
			into this.
*/
  bool readCallList();

/*!
  \var		typedef vector<Gsvolu>::iterator vol_iterator;
  \brief	iterator for the vector of GsVolu
*/
  typedef std::vector<Gsvolu>::iterator vol_iterator;

/*!
  \var		typedef list<Gspos>::iterator pos_iterator;
  \brief	iterator for the vector of GsVolu
*/
  typedef std::list<Gspos>::iterator pos_iterator;

/*!
  \var		typedef vector<Gstmed>::iterator med_iterator;
  \brief	iterator for the vector of GsVolu
*/
  typedef std::vector<Gstmed>::iterator med_iterator;

/*!
  \var		typedef vector<CsMaterial>::iterator mat_iterator;
  \brief	iterator for the vector of GsVolu
*/
  typedef std::vector<CsMaterial>::iterator mat_iterator;

/*!
  \var		typedef vector<Gsdet>::iterator det_iterator;
  \brief	iterator for the vector of GsVolu
*/
  typedef std::vector<Gsdet>::iterator det_iterator;


/*!	
  \fn		friend ostream& operator<<( ostream& stream, CsG3CallFile& callFile);
  \param	stream is a ostream
  \brief	output all information about this CsG3CallFile.
*/
  friend std::ostream& operator<<( std::ostream& stream, CsG3CallFile& callFile);

private:

/*!
  \var		string g3CallListFile_;
  \brief	Location of g3calls.dat file.
*/
  std::string g3CallListFile_;

/*!
  \var		vector<Gsvolu> volume_;
  \brief	a vector of GSVOLU in the GEANT 3 call list.
*/
  std::vector<Gsvolu> volume_;

/*!
  \var		list<Gspos> position_;
  \brief	a list of GSPOS in the GEANT 3 call list.
*/
  std::list<Gspos> position_;

/*!
  \var		vector<Gstmed> medium_;
  \brief	a vector of GSTMED in the GEANT 3 call list.
*/
  std::vector<Gstmed> medium_;

/*!
  \var		vector<CsMaterial> material_;
  \brief	a vector of material information in the GEANT 3 call list.
*/
  std::vector<CsMaterial> material_;

/*!
  \var		map< int, HepMatrix, less<int> > rotation_;
  \brief	a map of rotation matricses in the GEANT 3 call list.
*/
  std::map< int, CLHEP::HepMatrix, std::less<int> > rotation_;

/*!
  \var		vector<Gsdet> detector_;
  \brief	a vector of GSDET in the GEANT 3 call list.
*/
  std::vector<Gsdet> detector_;

/*!
  \var		Gspos HALL_	
  \brief	HALL volume.
*/
  Gspos	HALL_;

/*!
  \var		Gsvolu HALLVolume_
*/
  Gsvolu HALLVolume_;

public:

/*!
  \enum	SYSMODE
*/
  enum SYSMODE {COMGEANT, COMPASS};

  SYSMODE systemMode() const {return systemMode_;}	//!< get system mode
//! set system mode
  void systemMode( const SYSMODE& systemMode ) {systemMode_ = systemMode;}

private:
/*!
  \var		SYSMODE systemMode_;
  \brief	system mode COMGEANT or COMPASS
*/
  SYSMODE systemMode_;

/*!
  \fn		vector<Gsvolu>& volume() {return volume_;}
  \return	a reference of vector<Gsvolu> 
  \brief	this returns the referece of the data member volume_ of this object.
*/
  std::vector<Gsvolu>& volume() {return volume_;}

public:

/*!
  \fn		list<Gspos>& position() {return position_;}
  \return	a reference of list<Gspos> (volume_ in this)
  \brief	this returns the referece of the data member position_ of this object.
*/
  std::list<Gspos>& position() {return position_;}

private:

  inline std::string g3CallListFile() const {return g3CallListFile_;} //!< get filename
//! set filename
  inline void g3CallListFile(const std::string& filename) {
	  g3CallListFile_ = filename;
  }
//! return a reference to vector of Gstmed
  std::vector<Gstmed>& medium() {return medium_;}
//! return a reference to vector of CsMaterial
  std::vector<CsMaterial>& material() {return material_;}
//! return a referece to map of int and rotation matrix
  std::map<int, CLHEP::HepMatrix, std::less<int> >& rotation() {return rotation_;}
//! return a reference to vector of Gsdet
  std::vector<Gsdet>& detector() {return detector_;}

  bool readCallList(std::istream& stream);	//!< read call file
  void readGSMATE(std::istream& stream);		//!< read GSMATE entry
  void readGSMIXT(std::istream& stream);		//!< read GSMIXT entry
  void readGSVOLU(std::istream& stream);		//!< read GSVOLU entry
  void readGSPOS(std::istream& stream);		//!< read GSPOS entry
  void readGSTMED(std::istream& stream);		//!< read GSTMED entry
  void readGSROTM(std::istream& stream);		//!< read GSROTM entry
  void readGSDET(std::istream& stream);		//!< read GSDET entry
  void readGSDETH(std::istream& stream);		//!< read GSDETH entry
  void readGSDETD(std::istream& stream);		//!< read GSDETD entry
  void readGSDETU(std::istream& stream);		//!< read GSDETU entry

  std::string readName(std::istream& stream);		//!< read name with in ''


/*!
  \fn		void findMothers(const pos_iterator& ip);
  \param	ip is a CsG3CallFile::pos_iterator
  \brief	find mother Gspos of (*ip) in Gspos vector, 
			to complete volume tree structore.
*/
  void findMothers(const pos_iterator& ip);

//! find detector entry
  det_iterator findDetector( const std::string& name, const det_iterator& istart);
//! find detector entry
  det_iterator findDetector(const std::string& name);
//! find volume entry
  vol_iterator findVolume( const std::string& name, const vol_iterator& istart);
//! find volume entry
  vol_iterator findVolume(const std::string& name);
//! find posotion entry
  pos_iterator findPosition( const std::string& name, const pos_iterator& istart);
//! find posotion entry
  pos_iterator findPosition( const Gspos* pPos, const pos_iterator& istart);
//! find posotion entry
  pos_iterator findPosition( const int& id, const pos_iterator& istart);
//! find posotion entry
  pos_iterator findPosition(const int& id);
//! find posotion entry
  pos_iterator findPosition(const Gspos* pPos);
//! find posotion entry
  pos_iterator findPosition(const std::string& name);
//! find material entry
  mat_iterator findMaterial( const unsigned int& id, const mat_iterator& istart);	
//! find material entry
  mat_iterator findMaterial(const unsigned int& id);	
//! find medium entry
  med_iterator findMedium( const unsigned int& id, const med_iterator& istart);
//! find medium entry
  med_iterator findMedium(const unsigned int& id);

//! get g3 call file name from option file
  bool getFileName();

public:
/*!
  \fn		HepPoint3D centerInHALL(const Gspos& pos);
  \param	pos is a Gspos
  \return	a HepPoint3D
  \brief	return center position of pos in HALL system
*/
  HepGeom::Point3D<double> centerInHALL(const Gspos& pos);

/*!
  \fn		HepMatrix rotationToHALL(const Gspos& pos);
  \param	pos is a Gspos
  return	a HepMatrix
  \brief	return rotation matrix to HALL system
*/
  CLHEP::HepMatrix rotationToHALL(const Gspos& pos);

//! rerurn HALL
  Gspos& HALL() {return HALL_;}

};

#endif
