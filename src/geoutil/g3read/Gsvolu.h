//-*-Mode: C++;-*-
#ifndef _Gsvolu_h_
#define _Gsvolu_h_
/*!
  \file		Gsvolu.h
  \brief	volume entries in Geant 3. (header file)
  \author	$Author: benigno $
  \version	$Revisopn$
  \date		$Date: 2003/04/24 13:26:03 $
*/

#if __GNUG__ >= 2
#  pragma interface
#endif

#include "CsSTD.h"
#include "Gstmed.h"
//---- Gsvolu -----------------------------------------------------------
/*!
  \class	Gsvolu
  \brief	volume entry in G3 call. 
			The attributes of this objects are described in GEANT 3 manuals.
*/
class Gsvolu {
public:
/*!
  \fn		Gsvolu( const string& name, const string& type, 
			  const unsigned int& materialId, const int& nparam,
			  const vector<double>& param );
  \param	name is a string
  \param	type is a string
  \param	materialId is a unsigned int
  \param	nparam is a int
  \param	param is a vector<double>
  \brief	constructor with specific valiables
*/
  Gsvolu( const std::string& name, const std::string& type, 
		  const unsigned int& materialId, const int& nparam,
		  const std::vector<double>& param );

  Gsvolu(const Gsvolu& volume); //!< copy constrructor

  ~Gsvolu(); //!< destructor

  //! equal to operator
  bool operator==(const Gsvolu& volume) const {
	  return this->name() == volume.name();
  }

  //! less than operator
  bool operator<(const Gsvolu& volume) const {
	  return this->materialId() < volume.materialId();
  }

  //! output to ostream
  friend std::ostream& operator<<(std::ostream& os, Gsvolu& vol );

private:
  std::string name_;					//!< name
  std::string type_;					//!< type
  unsigned int materialId_;		//!< material id
  int nparam_;					//!< number of parameters
  std::vector<double> param_;		//!< parameter
  Gstmed* pMedium_;				//!< pointer to medium

public:
  std::string name() const {return name_;} //!< get name
  std::string type() const {return type_;} //!< get type
  unsigned int materialId() const {return materialId_;} //!< get material id
  int nparam() const {return nparam_;} //!< get number of parameter
  std::vector<double> param() const {return param_;} //!< get paramters

  Gstmed* pMedium() const {return pMedium_;} //!< get pointer to Gstmed
  void pMedium(Gstmed* pMedium) {pMedium_ = pMedium;} //!< set pointer to Gstmed

};

#endif
