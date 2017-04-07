// $Id: CsZone.h,v 1.5 2010/12/19 19:15:52 ybedfer Exp $

/*!
   \file    CsZone.h
   \brief   Compass Apparatus geometrical zone
   \author  Benigno Gobbo
   \version $Revision: 1.5 $
   \date    $Date: 2010/12/19 19:15:52 $
*/

#ifndef CsZone_h
#define CsZone_h

#include "CsSTD.h"

class CsDetector;

/*! \class CsZone

    \brief Compass Apparatus geometrical zone.

    This class is intended to split then Compass apparatus into geometrical
    defined zones. At present only the Z coordinate can be used to define
    the zones limits (start and end point). No checks will be made on
    zones consistency and/or overlaps. Zones are used in tracking algorithms.
*/

class CsZone {

 public:

  /*! \fn CsZone( float Zmin, float Zmax, string name, list<CsDetector*> alldets )
      \param Zmin zone starting point along beam direction (mm)
      \param Zmax zone end point along beam direction (mm)
      \param name zone name
      \param alldets list of all detectors
  */
  CsZone( float Zmin, float Zmax, std::string name, std::list<CsDetector*> alldets );

  CsZone( const CsZone& zone );             //!< Copy Constructor

  CsZone& operator=( const CsZone& zone );  //!< Assign Operator

  bool operator<( const CsZone& zone );     //!< Less than Operator 

  bool operator==( const CsZone& zone );    //!< Equal to Operator 

  /*! \fn std::string getName() const
      \brief Returns the zone name 
  */
  inline std::string getName() const { return( name_ ); }

  /*! \fn unsigned int getId() const
      \brief Returns the zone unique identifier 
  */
  inline unsigned int getId() const { return( id_ ); }

  /*! \fn float getZMin() const
      \brief Returns the zone starting point along Z (mm)
  */
  inline float getZMin() const { return( zmin_ ); }

  /*! \fn float getZMin() const
      \brief Returns the zone end point along Z (mm)
  */
  inline float getZMax() const { return( zmax_ ); }

  /*! \fn list<CsDetector*> getDets() const
    \brief Returns the list of detectors inside the zone
  */
  inline std::list<CsDetector*> getDets() const { return( dets_ ); }

  /*! \fn CsDetector* getFirstDet() const
    \brief Returns the detector at lower Z inside zone 
  */
  inline CsDetector* getFirstDet() const { return( firstDet_ ); }

  /*! \fn CsDetector* getLastDet() const
    \brief Returns the detector at higher Z inside zone 
  */
  inline CsDetector* getLastDet() const { return( lastDet_ ); }

 private:

  std::string        name_;          //!< zone name
  unsigned int       id_;            //!< zone unique identifier          
  float              zmin_;          //!< zone starting point along Z (mm)
  float              zmax_;          //!< zone end point along Z (mm)
  std::list <CsDetector*> dets_;     //!< list of detectors in zone
  CsDetector*        firstDet_;      //!< detector at lower Z in zone 
  CsDetector*        lastDet_;       //!< detector at higher Z in zone 
};

#endif // CsZone_h
