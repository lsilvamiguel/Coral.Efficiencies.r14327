// $Id: CsMCVertex.h,v 1.5 2003/04/09 15:48:30 benigno Exp $

/*!
   \file    CsMCVertex.h
   \brief   Compass Monte Carlo Vertex Class.
   \author  Benigno Gobbo
   \version $Revision: 1.5 $
   \date    $Date: 2003/04/09 15:48:30 $
*/
#ifndef CsMCVertex_h
#define CsMCVertex_h

#include <list>
#include "CsMCTrack.h"

class CsMCTrack;

/*! \class CsMCVertex 
    \brief   Compass Monte Carlo Vertex Class.

    Collect all needed Monte Carlo informations from Geant VERT bank. 
    Objects of CsMCVertex class are istantiated by CsGeant3 methods.
*/

class CsMCVertex {

 public:

  CsMCVertex(); //!< Default Constructor

  /*! \fn CsMCVertex( int Gnum, double x, double y, double z, double t );
    \brief Constructor
    \param Gnum Geant Vertex number 
    \param x    Vertex X position (mm)
    \param y    Vertex Y position (mm)
    \param z    Vertex Z position (mm)
    \param t    Time of flight (us)
   */
  CsMCVertex( int Gnum, double x, double y, double z, double t );

  /*! \fn CsMCVertex( int Gnum, double x, double y, double z, double t,
    CsMCTrack& inTrack);
    \brief Constructor
    \param Gnum Geant Vertex number 
    \param x    Vertex X position (mm)
    \param y    Vertex Y position (mm)
    \param z    Vertex Z position (mm)
    \param t    Time of flight (us)
    \param inTrack Pointer to the incoming CsMCTrack
   */
  CsMCVertex( int Gnum, double x, double y, double z, double t, 
	      CsMCTrack& inTrack);

  CsMCVertex( const CsMCVertex& ); //!< Copy Constructor

  CsMCVertex& operator=( const CsMCVertex& ); //!< Assign operator
  
  bool operator==( const CsMCVertex& ) const; //!< "equal to" operator
  bool operator<( const CsMCVertex& ) const;  //!< "less than" operator

  double getX() const; //!< Returns the Vertex X coordinate (MRS) (mm)

  double getY() const; //!< Returns the Vertex Y coordinate (MRS) (mm)

  double getZ() const; //!< Returns the Vertex Z coordinate (MRS) (mm)

  double getT() const; //!< Returns the Time of Flight (us)

  int    getGnum() const; //!< Returns the Geant VERT id number

  /*! \fn CsMCTrack* getInTrack() const;
      \brief Returns pointer to the associated incoming CsMCTrack object
      \warning{Vertices could have no incoming tracks, in this case the 
               pointer is NULL}
  */
  const CsMCTrack* getInTrack() const;

  /*! \fn std::list<CsMCTrack*> getOutTracks();
      \brief Returns pointer to the associated outgoing CsMCTracks object.
  */
  std::list<CsMCTrack*> getOutTracks();

  /*! \fn void setInTrack( CsMCTrack& inTrack );
    \brief Store the pointer to the associated incoming CsMCTrack object
  */
  void setInTrack( CsMCTrack& inTrack );

  /*! \fn void addOutTrack( CsMCTrack& outTrack );
   \brief Add a pointer to an associated outgoing CsMCTrack object
  */
  void addOutTrack( CsMCTrack& outTrack );  

 private:

  int    Gnum_;                 //!< Geant VERT id number
  double x_;                    //!< vertex X coordinate (mm)
  double y_;                    //!< vertex Y coordinate (mm)
  double z_;                    //!< vertex Z coordinate (mm)
  double t_;                    //!< time of flight
  CsMCTrack* inTrack_;          //!< pointer to incoming track
  std::list <CsMCTrack*> outTracks_; //!< list of pointers to outgoing tracks

};

#endif // CsMCVertex_h
