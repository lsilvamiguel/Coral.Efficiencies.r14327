#ifndef _CsMaterial_h_
#define _CsMaterial_h_
/*!
  \file		CsMaterial.h
  \brief	Compass geometry material class haader file
  \author	$Author: benigno $
  \version	$Revisopn$
  \date		$Date: 2003/04/24 13:26:02 $
*/
#if __GNUG__ >= 2
#  pragma interface
#endif

// include CORAL standard file
#include <CsSTD.h>

/*!
  \class	CsMaterial
  \brief	class for Material of Geometrical Volume
*/
class CsMaterial {
public:
//! default constructor
  CsMaterial();
/*!
  \fn		CsMaterial(	const unsigned int& materialId,
				const string& name,
				const double& density,
				const double& radiationLength);
  \param	materialId is a unsigned int
  \param	name is a string
  \param	density is a double
  \param	radiationLength is a double
  \brief	constructor with specific variables
*/
  CsMaterial(
	const unsigned int& materialId, const std::string& name,
	const double& density, const double& radiationLength);
	
/*!
  \fn		CsMaterial(	const unsigned int& materialId, const string& name,
						const double& density, const int& nmat,
						vector<double>& a, vector<double>& z, vector<double>& w);
  \param	materialId is unsigned int
  \param	name is a string
  \param	density is a double
  \param	nmat is a int
  \param	a is a vector of double
  \param	z is a vector of double
  \param	w is a vector of double
  \brief	Construct CsMaterial object for material which consists of several
			medium, and the ratio among these madium was given in a way used in
			GEANT 3. GSMIXT definition.
*/
  CsMaterial(const unsigned int& materialId, const std::string& name,
	const double& density, const int& nmat,
	std::vector<double>& a, std::vector<double>& z, std::vector<double>& w);

//! copy constructor
  CsMaterial(const CsMaterial& A);
 
//! A destructor
  virtual ~CsMaterial();

//! assignment operator
  CsMaterial& operator=(const CsMaterial& A);

/*!
  \fn		unsigned int id() const
  \brief	return the id number.
*/
  virtual unsigned int id() const { return id_;}
	
/*!
  \fn		void id(const unsigned int& id)
  \brief	set the id number with the given number.
*/
  virtual void id(const unsigned int& id) { id_=id;}

/*!
  \fn		string name() const 
  \brief	return the name.
*/
  virtual std::string name() const { return name_;}

/*!
  \fn		void name(const string& name)
  \brief	set the name with the given string.
*/
  virtual void name(const std::string& name) { name_ = name;}

/*!
  \fn		double radiationLength() const
  \brief	return the radiation length.
*/
  virtual double radiationLength() const { return radiationLength_;}

/*!
  \fn		void radiationLength(const double& radiationLength)
  \brief	set the radiation length with the given number.
*/
  virtual void radiationLength(const double& radiationLength) { 
	radiationLength_ = radiationLength;
  }
  
/*!
  \fn		double density() const
  \brief	return the density.
*/
  virtual double density() const { return density_;}

/*!
  \fn		void density(const double& density)
  \brief	set the density with the given number.
*/
  virtual void density(const double& density) { density_ = density;}


//!	Input from stream to CsMaterial object.
  friend std::istream& operator>>(std::istream& stream, CsMaterial& material);
//!	Output from CsMaterial object to stream.
  friend std::ostream& operator<<(std::ostream& stream, const CsMaterial& material);

private:
  unsigned int id_;				//!< id number
  std::string     name_;				//!< name
  double radiationLength_;		//!< radiation length
  double density_;				//!< density

/*!
  \fn		double rLength( const int& nmat, 
			  vector<double>& a, vector<double>& z, vector<double>& w) const ;
  \param	nmat is an int
  \param	a is a vector of double
  \param	z is a vector of double
  \param	w is a vector of double
  \brief	calculate radiation length for material which was defined in the way
			to construct material in GEANT 3 by GSMIXT call.
*/  
  double rLength( const int& nmat, 
	std::vector<double>& a, std::vector<double>& z, std::vector<double>& w) const ;

//! caluclate coulomb correction of material with Z.	
  double coulombCorrection( const double& z ) const ;

};
	

#endif
