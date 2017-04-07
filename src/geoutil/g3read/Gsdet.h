#ifndef _Gsdet_h_
#define _Gsdet_h_
/*!
  \file		Gsdet.h
  \brief	Object for GSDET entries in G3 call (header file)
  \author	$Author: benigno $
  \version	$Revisopn$
  \date		$Date: 2003/04/24 13:26:02 $
*/

#include "CsSTD.h"
#include "CsMaterial.h"

#if __GNUG__ >= 2
#  pragma interface
#endif

//---- Gsdet -----------------------------------------------------------
/*!
	\class	Gsdet
	\brief	This holds the information on GSDET line in G3 call list.
			The discription of attributes of GSDET will be found in GEANT3 
			reference manual, and COMGEANT documentation. 
			Some of them are not clear yet, at this morment.
*/
class Gsdet {
public:
//!	a default constructor
	Gsdet();

/*!
  \fn		Gsdet(const string& name, const string& mother, const int& type); 
  \param	name is a string
  \param	mother is a string
  \param	type is a int
  \brief	A constructor with three arguments, detector name, mother name,
			and detector type id.
*/
  Gsdet( const std::string& name, const std::string& mother, const int& type );

//!	copy constructor
  Gsdet( const Gsdet& gsdet );

//!	assignment operator
  Gsdet& operator=( const Gsdet& gsdet );
  
//!	destructor	
  ~Gsdet();

private:
  std::string name_;					//!< name
  std::string mother_;				//!< mother volume name ?
  int    type_;					//!< type
  CsMaterial* pMaterial_;		//!< pointer to CsMaterial
  std::vector<std::string> hitName_;		//!< hit name?
  std::vector<int> hitBit_;			//!< hit bit
  std::vector<double> hitFactor_;	//!< hit factor
  std::vector<double> hitOffset_;	//!< hit offset
  std::vector<std::string> digitName_;	//!< digit name
  std::vector<int> digitBit_;		//!< digit bit
  std::vector<double> user_;			//!< user defined parameter

public:
  std::vector<std::string> hitName() const {return hitName_;}	//!< get hitname
  std::vector<int> hitBit() const {return hitBit_;}		//!< get hitbit
  std::vector<double> hitFactor() const {return hitFactor_;} //!< get hit factor
  std::vector<double> hitOffset() const {return hitOffset_;} //!< get hit offset
  std::vector<std::string> digitName() const {return digitName_;} //!< get digit name
  std::vector<int> digitBit() const {return digitBit_;} 	//!< get digit bit
  std::vector<double> user() const {return user_;}			//!< get user param

  void hitName(const std::vector<std::string>& name) {hitName_ = name ;} //!< set hit name
  void hitBit(const std::vector<int>& bit) {hitBit_ = bit;} //!< set hit bit
  void hitFactor(const std::vector<double>& factor) {hitFactor_ = factor;} //!< set hit factor
  void hitOffset(const std::vector<double>& offset) {hitOffset_ = offset;}  //! set hit offset
  void digitName(const std::vector<std::string>& name) {digitName_=name;} //!< set digit name
  void digitBit(const std::vector<int>& bit) {digitBit_ = bit;} //!< set digit bit
  void user(const std::vector<double>& user) {user_ = user;} //!< set user param

  void addHitName(const std::string& name) {hitName_.push_back(name);} //!< add hit name
  void addHitBit(const int& bit) {hitBit_.push_back(bit);} //!< add hit bit
  void addHitFactor(const double& factor) {hitFactor_.push_back(factor);} //!< add hit factor
  void addHitOffset(const double& offset) {hitOffset_.push_back(offset);} //!< add hit offset
  void addDigitName(const std::string& name)	{digitName_.push_back(name);} //!< add digit name
  void addDigitBit(const int& bit) {digitBit_.push_back(bit);} //!< add digit bit
  void addUser(const double& param) {user_.push_back(param);} //!< add user param

  std::string name() const {return name_;} //!< return name
  std::string mother() const {return mother_;} //!< return mother name
  int type() const {return type_;} //!< return type
  CsMaterial* pMaterial() const {return pMaterial_;} //!< return pointer to CsMaterial

  void name( const std::string& name ) {name_ = name;}	//!< set name
  void mother( const std::string& mother ) {mother_ = mother;} //!< set mother name
  void type( const int& type ) {type_ = type;} //!< set type
  void pMaterial(CsMaterial* pMaterial) {pMaterial_ = pMaterial;} //!< set pointer to CsMaterial

// this has to be consistent with the omgeomd.F 
// but in the sense of C++ array...
/*!
  \enum		userTag
  \brief	This is index of user parameter defined GSDETU call in omgeomd.F
*/
  enum userTag { 
	  NUNIT	=  0,
	  WirN    =  1, 
	  WirP    =  2,
	  CosAng	=  3,
	  SinAng  =  4,		 
	  WirD    =  5,
	  X1stWir =  6,
	  Y1stWir =  7,
	  Z1stWir =  8,
	  Xpitch  =  9,
	  Ypitch  = 10,
	  Zpitch  = 11,
	  WirLengthY = 12,
	  TimeGate= 13,
  };

//! equal to operator
  bool operator==(const Gsdet& det) const {
	  return this->type() == det.type();
  }

//! less than operator
  bool operator<(const Gsdet& det) const {
	  return this->type() < det.type();
  }

};

#endif
