#ifndef CsVolumeBase_h
#define CsVolumeBase_h
/*!
  \file CsVolumeBase.h
  \brief Abstract Class for managin volume inside Geometry Tools
*/

/*!
  \class CsVolumeBase
  \brief Abstract Class for managin volume inside Geometry Tools
*/

class CsVolumeBase {
public:
//! return id number
	virtual int id() const = 0;

//! return unit number
	virtual int unit() const = 0;

//! return magnetic filed number
	virtual int magId() const = 0;

//! return name
	virtual string name() const = 0;

//! return mother volume
	virtual CsVolumeBase* pMother() const = 0;
	
};
#endif






