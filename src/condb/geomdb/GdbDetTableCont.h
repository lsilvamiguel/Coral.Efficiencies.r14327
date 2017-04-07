//-*-Mode: C++;-*-
#ifndef _GdbDetTableCont_h_
#define _GdbDetTableCont_h_
/*!
  \file		GdbDetTableCont.h
  \brief	Geometry DB detector table object definition file
  \author	$Author: miyachi $
  \version	$Revision: 1.2 $
  \date		$Date: 2000/06/07 08:36:45 $
*/

#if __GNUG__ >= 2
#  pragma interface
#endif

#include "CsGeoDbDetector.h"
/*!
  \class	GdbDetTableCont
  \brief	Geometry DB Detector table objects, 
			which will consist of equivalent information found in 
			detectors.dat file.
*/
class GdbDetTableCont : public CsGeoDbDetector {
public:

  GdbDetTableCont();	//!< constructor

  //! constructor with specific values
  GdbDetTableCont( const int& id, const char* name, const int& unit,
	  const int& type, const double& radLen, const HepPoint3D& size,
	  const HepPoint3D& center, const HepMatrix& rotation, 
	  const double& wirD, const double& angle, const int& nWir,
	  const double& wirP, const double& eff, const double& bgk,
	  const double& tGate );

  //! constructor with specific values
  GdbDetTableCont( const int& id, const char* name, const int& unit,
	  const int& type, const double& radLen, const int& solidtype,
	  const vector<double>& param, const HepPoint3D& center, 
	  const HepMatrix& rotation, const double& wirD,
	  const double& angle, const int& nWir,
	  const double& wirP, const double& eff,
	  const double& bgk, const double& tGate );

  GdbDetTableCont(const GdbDetTableCont& tableCont); //!< copy constructor

  virtual ~GdbDetTableCont(); //!< destructor

//! assignment operator
  GdbDetTableCont& operator=(const GdbDetTableCont& tableCont);

//! equal to operator
  bool operator==(const GdbDetTableCont& tableCont) const {
	  return	this->id() == tableCont.id() ;
  }

//! less than operator
  bool operator<(const GdbDetTableCont& tableCont) const {
	  return	this->id() < tableCont.id() ;
  }

  int id() const {return id_;}	//!< get id
  const char* name() const {return name_;}	//!< get name
  int unit() const {return unit_;}	//!< get unit
  int type() const {return type_;}	//!< get type
  double radLen() const {return radLen_;}	//!< get radiation length

  HepPoint3D center() const {return center_;}	//!< get center point
  HepPoint3D size() const ;	//!< get size
  HepMatrix rotation() const {return rotM_;}	//!< get rotation matrix

  double wirD() const {return wireDistance_;}	//!< get wire distance
  double angle() const {return angle_;}		//!< get wire angle
  int nWir() const {return wireNumber_;}		//!< get wire number
  double wirP() const {return wirePitch_;}	//!< get wire pitch

  double eff() const {return efficiency_;}	//!< get efficiency
  double bgk() const {return backGround_;}	//!< get background
  double tGate() const {return timeGate_;}	//!< get time gate

  bool hasDrift() const {return hasDrift_;}	//!< get drift flag

  double driftV() const {return driftVelocity_;}	//!< get drift velocity
  double T0() const {return driftT0_;}	//!< get t0

  double dHitRes() const {return doubleHitResolution_;}	//!< get 2 hit resolution
  double sRes() const {return spaceResolution_;}	//!< get space resolution
  double tSlic() const {return timeSlice_;}	//!< get time alice



  void id( const int& id ) { id_ = id ;}	//!< set id
  void name( const char* name );	//!<  set name
  void unit( const int& unit ) { unit_ = unit ;}	//!< set unit
  void type( const int& type ) { type_ = type ;}	//!< set type
  void radLen( const double& radLen ) { radLen_ = radLen ;}	//!< set radiation length

  void center( const HepPoint3D& point ) {center_ = point;}	//!< set center
  void rotation( const HepMatrix& matrix) {rotM_ = matrix;}	//!< set rotation matrix

  void wirD(const double& wirD) { wireDistance_ = wirD;} 	//!< set wire distance
  void angle(const double& angle) { angle_ = angle;}	//!< set wire angle
  void nWir(const int& nWir) { wireNumber_ = nWir;}	//!< set wire number
  void wirP(const double& wirP) { wirePitch_ = wirP;} 	//!< set wire pitch

  void eff(const double& eff) { efficiency_ = eff;}	//!< set efficiency
  void bgk(const double& bgk) { backGround_ = bgk;}	//!< set background
  void tGate(const double& tGate) { timeGate_ = tGate;}	//!< set time gate

  void hasDrift(const bool& hasDrift) {hasDrift_ = hasDrift;}	//!< set drift flag

  void driftV(const double& driftV) { driftVelocity_ = driftV;}	//!< set drift velocity
  void T0(const double& T0) { driftT0_ = T0;}	//!< set t0

  void dHitRes(const double& dHitRes) { doubleHitResolution_ = dHitRes;}	//!< set 2hist resolution
  void sRes(const double& sRes) { spaceResolution_ = sRes;}	//!< set space resolution
  void tSlic(const double& tSlic) { timeSlice_ = tSlic;}	//!< set time slice

//! set drift information at once
  void setDriftInfo( const double& driftV, const double& T0, 
	const double& dHitRes, const double& sRes, const double& tSlic);

/*!
  \fn		void dumpDetTableHeader(ostream& os);
  \param	os is an ostream
  \brief	dump header information to os
*/
  void dumpDetTableHeader(ostream& os);

/*!
  \fn		void dumpDetTable(ostream& os);
  \param	os is an ostream
  \brief	dump detector table information to os in detectors.dat format
*/
  void dumpDetTable(ostream& os);

/*!
  \fn		void dumpDeadSpace(ostream& os);
  \param	os is an ostream
  \brief	dump dead space information to os in detectors.dat format
*/
  void dumpDeadSpace(ostream& os);

//! dump information to ostream
  friend ostream& operator<<(ostream& os, GdbDetTableCont& tableCont);

//! return reference to volume paramter
  vector<double>& refParam() {return param_;}	
//! return a copy of volume parameter
  vector<double> param() const {return param_;}
  
  void param(const vector<double> param) {param_ = param;}	//!< set volume parameter

  int solidtype() const {return solidtype_;}	//!< get solid type
  void solidtype(const int& solidtype) {solidtype_ = solidtype;}	//!< set solid type

  bool hasDeadSpace()	const {return hasDeadspace_;}	//!< get dead space flag
  void hasDeadSpace( const bool& hasdead) {hasDeadspace_ = hasdead;}	//!< set dead space flag

  int dId() const {return dId_;}	//!< get dead space id
  void dId(const int& dId) {dId_ = dId;}	//!< set dead space id

  int dType() const {return dType_;}	//!<  get dead space volume type
  void dType(const int& dType) {dType_ = dType;}	//!< set dead volume space type

//! get a copy of dead space volume parameter
  vector<double> dParam() const {return dParam_;}
//! get a reference of dead space volume parameter	
  vector<double>& refDParam() {return dParam_;}
//! set dead space volume parameter
  void dParam(const vector<double>& dparam ) {dParam_ = dparam;}

  HepPoint3D dCenter() const { return dCenter_;}	//!< get dead space center
  HepMatrix  dRotM() const { return dRotM_;}	//!< get dead space rotation matrix

//! set dead space volume center
  void dCenter(const HepPoint3D& point) {dCenter_ = point;} 
//! set dead space volume rotation matrix
  void dRotM(const HepMatrix& rot) {dRotM_ = rot;}

//! set dead zone information at once
  void deadZone(	const int& type, const vector<double>& param,
				  const HepPoint3D& center, const HepMatrix& rot);

private:
  static const int MAXCHAR = 16;	//!< maximun name length
  int		id_;				//!< id
  char	name_[MAXCHAR];			//!< name
  int		unit_;				//!< unit
  int		type_;				//!< type
  double	radLen_;			//!< radiation length
  int solidtype_;				//!< solid type (see CsCOMGEANT.h)
  vector<double> param_;		//!< volume parameter

  HepPoint3D center_;			//!< center
  HepMatrix  rotM_;				//!< rotation matrix

  double wireDistance_;			//!< wire distance
  double angle_;				//!< wire angle
  int wireNumber_;				//!< number os wires
  double wirePitch_;			//!< wire pitch
  double efficiency_;			//!< efficiency
  double backGround_;			//!< background
  double timeGate_;				//!< time gate

  bool hasDrift_;				//!< drift detector flag
  double driftVelocity_;		//!< drift velocity
  double driftT0_;				//!< drift t0
  double doubleHitResolution_;	//!< double hits resolution
  double spaceResolution_;		//!< space resolution
  double timeSlice_;			//!< time slice

//------ Dead Space Information.
  bool hasDeadspace_;			//!< dead zone flag
  int dId_;						//!< dead zone id
  int dType_;					//!< dead zone volume type
  vector<double> dParam_;		//!< dead volume parameter
  HepPoint3D dCenter_;			//!< dead volume center
  HepMatrix  dRotM_;			//!< dead volume rotation matrix

};

#endif
