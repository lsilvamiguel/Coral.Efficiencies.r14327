//-*-Mode: C++;-*-
#ifndef _Gstmed_h_
#define _Gstmed_h_
#if __GNUG__ >= 2
#  pragma interface
#endif
/*!
  \file		Gstmed.h
  \brief	Object for Material entries in Geant 3. (header file)
  \author	$Author: benigno $
  \version	$Revisopn$
  \date		$Date: 2003/04/24 13:26:03 $
*/
#include "CsSTD.h"
#include "CsMaterial.h"
/*!
  \class	Gstmed
  \brief	Naterial entries in G3 call list
*/
class Gstmed{
public:
  Gstmed();	//!< defalt constructor

/*!
  \fn		  Gstmed(	const unsigned int& id, 
						const string& name, 
						const unsigned int& materialId,
						const int& isvol, 
						const int& ifield)
  \param	id is a unsigned int
  \param	name is a string
  \param	materialId is a unsigned int
  \param	isvol is a int
  \param	ifield is a int
  \brief	constructor with specific valiables

*/
  Gstmed(
	const unsigned int& id, const std::string& name, const unsigned int& materialId,
	const int& isvol, const int& ifield);

  Gstmed(const Gstmed& gstmed); //!< copy constructor

  ~Gstmed();	//!< destructor

  Gstmed& operator=(const Gstmed& gstmed); //!< assignment operator
  
//! equal to operator
  bool operator==(const Gstmed& medium) const {
	return this->id() == medium.id();
  }

//! less than operator
  bool operator<(const Gstmed& medium) const {
	return this->id() < medium.id();
  }

//! output information to os
  friend std::ostream& operator<<( std::ostream& os, const Gstmed& med);

private:
  unsigned int id_;				//!< id
  std::string name_;					//!< name
  unsigned int materialId_;		//!< material id
  int	isvol_;					//!< volume flag?
  int ifield_;					//!< field flag?

  CsMaterial* pMaterial_;		//!< pointer to CsMaterial

public:
  unsigned int id() const {return id_;}			//!< get id
  void id(const unsigned int& id) {id_ = id;}	//!< set id

  std::string name() const {return name_;}				//!< get name
  void name(const std::string& name) {name_ = name;}	//!< set name

  //! get material id
  unsigned int materialId() const {return materialId_;}
  //! set material id
  void materialId(const unsigned int& materialId) {materialId_ = materialId;}

  int isvol() const {return isvol_;}					//!< get isvol flag
  void isvol(const int& isvol) {isvol_ = isvol;}	//!< set isvol flag

  int ifield() const {return ifield_;}					//!< get field flag
  void ifield(const int& ifield) {ifield_ = ifield;}	//!< set field flag

  //! get pointer to CsMaterial
  CsMaterial* pMaterial() const {return pMaterial_;}

  //!< set pointer to CsMaterial
  void pMaterial(CsMaterial* pMaterial) {pMaterial_ = pMaterial;}

};

#endif
