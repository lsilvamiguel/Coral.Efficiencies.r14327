//-*-Mode: C++;-*-
#ifndef _GdbDetectorType_h_
#define _GdbDetectorType_h_
/*!
  \file		GdbDetectorType.h
  \brief	Geometry DB detector type object definition file
  \author	$Author: miyachi $
  \version	$Revision: 1.2 $
  \date		$Date: 2000/06/07 08:36:45 $
*/

#if __GNUG__ >= 2
#  pragma interface
#endif

#include "GdbSTD.h"
/*!
  \class	GdbDetectorType
  \brief	Detector type object
*/
class GdbDetectorType {
public:
    GdbDetectorType();	//!< default constructor

/*!
  \fn		GdbDetectorType( const char* name, const int& type, 
			  const double& gate = 0.0)
  \param	name is a pointer to char
  \param	type is a int
  \param	gate is a double with default value 0.0
  \brief	constructor with each variables
*/
	GdbDetectorType( const char* name, const int& type, const double& gate = 0.0);

/*!
  \fn		GdbDetectorType( const char* name, const int& type, 
			  const double& gate, const double& eff, const double& bkg );
  \param	name is a pointer to char
  \param	type is a int
  \param	gate is a double
  \param	eff is a double
  \param	bkg is a double
  \brief	constructor with each variables

*/
	GdbDetectorType( const char* name, const int& type, const double& gate,
		const double& eff, const double& bkg );

//! copy constructor
    GdbDetectorType(const GdbDetectorType& dType);

//! destructor
    ~GdbDetectorType();

//! assignment operator
	GdbDetectorType& operator=(const GdbDetectorType& dType);

	const char* name(const char* name);	//!< set name and return it
	const char* name() const;			//!< get name
	int type() const;					//!< get type
	double gate() const;				//!< get time gate

//! dump information to ostream
	friend ostream& operator<<(ostream& os, const GdbDetectorType& det);

	double eff() const;		//!< get efficiency
	double bkg() const;		//!< get background
	bool hasDrift() const;	//!< true if this is drift tube
	double driftV() const;	//!< get drift velocity
	double T0() const;		//!< get t0
	double dHitRes() const;	//!< get double hit resolution
	double sRes() const;	//!< get space slice
	double tSlic() const;	//!< get time slice

	void hasDrift(const bool& hasDrift) ;	//!< set drift flag
	void driftV(const double& driftV) ;		//!< set drift velocity
	void T0(const double& T0);				//!< set t0
	void dHitRes(const double& dHitRes);	//!< set double hit resolution
	void sRes(const double& sRes);			//!< set space slice
	void tSlic(const double& tSlic);		//!< set time slice


/*!
  \fn		void setDriftInfo(	const double& driftV, const double& T0,
						const double& dHitRes, const double& sRes,
						const double& tSlic)
  \param	driftV is a double
  \param	T0 is a double
  \param	dHitRes is a double
  \param	sRes is a double
  \param	tSlic is a double
  \brief	set drift infromation at once
*/
	void setDriftInfo(	const double& driftV, const double& T0,
						const double& dHitRes, const double& sRes,
						const double& tSlic);

private:
	static const int maxNameLength = 8;	//!< maximun name length
	char detectorName_[maxNameLength];	//!< name
	int detectorType_;					//!< detector type id
	double gateLength_;					//!< time gate
	double eff_;						//!< efficiency
	double bkg_;						//!< background

	bool hasDrift_;						//!< drift flag
	double driftVelocity_;				//!< drift velocity
	double driftT0_;					//!< t0

	double doubleHitResolution_;		//!< double hit resolution
	double spaceResolution_;			//!< space resoluition
	double timeSlice_;					//!< time silice

};

#endif
