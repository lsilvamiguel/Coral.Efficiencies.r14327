// $Id: CsHelix.h,v 1.17 2010/10/15 12:02:38 schluter Exp $

/*!
   \file CsHelix.h
   \brief Compass Helix Parameters Class.
   \author  Benigno Gobbo
   \version $Revision: 1.17 $
   \date    $Date: 2010/10/15 12:02:38 $
*/
#ifndef CsHelix_h
#define CsHelix_h

#include <cmath>
#include <CLHEP/Matrix/Matrix.h>

/*! \class CsHelix

    \brief Compass Helix Class.

    This class describes the helix associated to a track.
*/

#if USE_ObjectsCounter
#include "ObjectsCounter.h"
#endif

class CsHelix {

 public:

  CsHelix();                              //!< Default Constructor

  /*! \fn CsHelix( double x, double y, double z, double dxdz, double dydz,
	   double cop, double* cov );
     \brief Constructor.
     \param x x coordinate of the helix starting point (MRS)
     \param y y coordinate of the helix starting point (MRS)
     \param z z coordinate of the helix starting point (MRS)
     \param dxdz dx/dz
     \param dydz dy/dz
     \param cop  charge/momentum
     \param cov  covariance matrix (this is a 15 elements array
     which is copied into the helix:<br>
     [0,0]=Cov[0],<br>
     [0,1]=Cov[1], [1,1]=Cov[2], <br>
     [0,2]=Cov[3], [1,2]=Cov[4], [2,2]=Cov[5], <br>
     [0,3]=Cov[6], [1,3]=Cov[7], [2,3]=Cov[8], [3,3]=Cov[9], <br>
     [0,4]=Cov[10],[1,4]=Cov[11],[2,4]=Cov[12],[3,4]=Cov[13],[4,4]=Cov[14] ).
  */
  CsHelix( double x, double y, double z, double dxdz, double dydz,
	   double cop, const double* cov );

  CsHelix( const CsHelix& );             //!< Copy Constructor

  CsHelix& operator=( const CsHelix& );  //!< Assign Operator

  bool operator==( const CsHelix& ) const;   //!< "equal to" Operator

  double & operator () (const int i, const int j=0);       //!< "()" accessor

  /*! \fn inline double getX() const;
    \brief Returns the X coordinate of the helix starting point
  */
  inline double getX()    const { return(x_); }

  /*! \fn inline double getY() const;
    \brief Returns the Y coordinate of the helix starting point
  */
  inline double getY()    const { return(y_); }

  /*! \fn inline double getZ() const;
    \brief Returns the Z coordinate of the helix starting point
  */
  inline double getZ()   const { return(z_); }

  /*! \fn inline double getDXDZ() const;
    \brief Returns dx/dz of the helix starting point
  */
  inline double getDXDZ() const { return(dxdz_); }

  /*! \fn inline double getDYDZ() const;
    \brief Returns dy/dz of the helix starting point
  */
  inline double getDYDZ() const { return(dydz_); }

  /*! \fn inline double getCop() const;
    \brief Returns charge/momentum of the reconstructed track
  */
  inline double getCop() const { return(cop_); }

  /*! \fn inline double getCov();
    \brief Returns the covariance matrix as a 15 elements array
  */
  inline const double* getCov() const { return(cov_); }

  /*! \fn inline double getCov();
    \brief Returns the covariance matrix as a 15 elements array
  */
  inline double* getCov() { return(cov_); }

  /*! \fn HepMatrix getCovMat();
    \brief Returns the covariance matrix.
  */
  CLHEP::HepMatrix getCovMat()const;
  
  /*! \fn double Dist( CsHelix& H ) const 
    \brief Distance between \a "this" helix and \a H. 
  */
  double Dist( CsHelix& H ) const {
    return(sqrt((z_-H.z_)*(z_-H.z_)+
		(x_-H.x_)*(x_-H.x_)+
		(y_-H.y_)*(y_-H.y_))); }

  /*! \fn bool Extrapolate(double Z, CsHelix& Hout, bool mmap=true )
    \brief Helix extrapolation
    
    Propagate track parameters and covariance matrix,
    represented by this helix on a plane, perpendicular
    to the beam direction.
    
    \param Z destination plane position 
    \param Hout helix at destination plane
    \param mmap if false no material map will be used
  */
  
  bool Extrapolate(double Z, CsHelix& Hout, bool mmap=true ) const;

  bool FindCDA(CsHelix& H);   //!< Find "this" helix nad helix \a H at point of closest distance of aproach.
  
  void Print(const char* str="") const;    //!< Print
  
 private:

  double x_;       //!< x coordinate
  double y_;       //!< y coordinate
  double z_;       //!< z coordinate
  double dxdz_;    //!< dx/dz
  double dydz_;    //!< dy/dz
  double cop_;     //!< charge/momentum
  double cov_[15]; //!< covariance matrix significative elements

  #if USE_ObjectsCounter
  ObjectsCounter<CsHelix> objects_counter;
  #endif

};


#endif // CsHelix_h
