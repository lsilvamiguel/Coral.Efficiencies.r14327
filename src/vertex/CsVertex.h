// $Id: CsVertex.h,v 1.18 2010/06/07 16:50:47 tnagel Exp $

/*!
   \file    CsVertex.h
   \brief   Compass Vertex Class.
   \author  Benigno Gobbo, Alexandre Korzenev
   \version $Revision: 1.18 $
   \date    $Date: 2010/06/07 16:50:47 $
*/
#ifndef CsVertex_h
#define CsVertex_h

#include "coral_config.h"

#include <cassert>
#include <cmath>
#include <list>
#include <vector>
#include <map>
#include <CLHEP/Matrix/Matrix.h>

#if USE_ObjectsCounter
#include "ObjectsCounter.h"
#endif

/*! \class Cs3Vector
    \brief Track Parameters Class.

    This class is analog of Hep3Vector class, but it contains
    \a dXdZ, \a dYdZ, \a Cop components (they are NOT projections 
    on three axis).
*/
class Cs3Vector {

 public:

  double dXdZ,dYdZ,Cop;
  Cs3Vector(){};
  Cs3Vector( double dxdz, double dydz, double cop ):
    dXdZ(dxdz), dYdZ(dydz), Cop(cop) {};
  Cs3Vector& operator=( const Cs3Vector& v ) //!< Assign operator
    { dXdZ = v.dXdZ; dYdZ = v.dYdZ; Cop = v.Cop; return *this;};
  bool operator==( const Cs3Vector& v) const //!< "equal to" operator
    { if( dXdZ == v.dXdZ && dYdZ == v.dYdZ && Cop == v.Cop) return( true );
      else return( false ); };
  double  operator () (const int i) {if(i==0) return dXdZ;
    else if(i==1) return dYdZ; else if(i==2) return Cop; return 1e10; };
  void getMom(double&px,double&py,double&pz) const {
    pz = 1/fabs(Cop)/sqrt(1+dXdZ*dXdZ+dYdZ*dYdZ);
    px = pz * dXdZ;
    py = pz * dYdZ;
  }
  double getMom() const {return 1/fabs(Cop);}
};

class CsTrack;

/*! \class CsVertex 
    \brief Vertex Class.

    This class describes the generic tracks vertex.  
*/

class CsVertex {

 public:

  CsVertex();           //!< Default Constructor

  /*! \fn CsVertex( double x, double y, double z );
    \brief Constructor
    \param x    Vertex X position (mm)
    \param y    Vertex Y position (mm)
    \param z    Vertex Z position (mm)
   */
  CsVertex( double x, double y, double z );

  /*! \fn CsVertex( double x, double y, double z, CsTrack& Track );
    \brief Constructor
    \param x    Vertex X position (mm)
    \param y    Vertex Y position (mm)
    \param z    Vertex Z position (mm)
    \param Track Pointer to an associated CsTrack
   */
  CsVertex( double x, double y, double z, CsTrack& Track);

  /*! \fn CsVertex( double x, double y, double z, std::list<CsTrack*> Tracks );
    \brief Constructor
    \param x    Vertex X position (mm)
    \param y    Vertex Y position (mm)
    \param z    Vertex Z position (mm)
    \param Tracks list of pointer to associated CsTrack objects
   */
  CsVertex( double x, double y, double z, const std::list<CsTrack*> &Tracks );

  CsVertex( const CsVertex& );            //!< Copy Constructor

  ~CsVertex();                            //!< Destructor

  CsVertex& operator=( const CsVertex& ); //!< Assign operator
  
  bool operator==( const CsVertex& ) const; //!< "equal to" operator

  //! Adds pointer to covariance matrix.
  inline void addCov( CLHEP::HepMatrix* cov ) { Cov_.push_back( cov ); }

  //! Adds a pointer to an associated CsTrack object
  inline void addTrack( CsTrack* Track ) 
    { Tracks_.push_back( Track ); }

  //! Adds track parameters at vertex. 
  inline void addTrackAtVertex( CsTrack* trk, const Cs3Vector &vec )
    { trkPar_[trk] = vec; }

  //! Clear list of associated tracks.
  inline void clearTracks() { Tracks_.clear(); }

  //! Clear vector of covariance matrices.
  void clearCov();

  //! true - primary, false - secondary.
  inline bool isPrimary() const { return type_; }
  
  //! Returns Chi-square at vertex after \a smoother procedure.
  inline double getChi2() const { return Chi2tot_; }

  //! Returns the Vertex X coordinate (MRS) (mm)
  inline double getX() const { return x_; }

  //! Returns the Vertex Y coordinate (MRS) (mm)
  inline double getY() const { return y_; }

  //! Returns the Vertex Z coordinate (MRS) (mm)
  inline double getZ() const { return z_; }

  //! Returns the pointer to Vertex covariance matrix (coordinate correlation)
  CLHEP::HepMatrix* getCov( int i ) const;

  //! Returns the vector of pointers to Vertex correlation matrices
  inline const std::vector<CLHEP::HepMatrix*> &getCov() const { return Cov_; }

  //! Returns number of tracks associated to this vertex.
  inline int getNTracks() const { return Tracks_.size(); }

  //! Returns true if there is vector of track parameters \a par at vertex corresponded to \a trk.
  bool getPar( CsTrack* trk, Cs3Vector& par );
  
  //inline const Cs3Vector &getPar( CsTrack* trk ) const 
  //  { return (trkPar_.find( trk ))->second; };

  //! Returns pointer to the associated CsTracks object.
  inline const std::list<CsTrack*> &getTracks() const { assert(this); return Tracks_; }

  //! Returns true if this vertex is best primary vertex.
  inline bool isBestVertex() const { return isBestVertex_; }

  //! Set Chi-square at vertex after \a smoother procedure.
  void setChi2(double chi2) { Chi2tot_=chi2; }

  //! Sets pointer to covariance matrix.
  inline void setCov( const std::vector<CLHEP::HepMatrix*> &cov );

  //! Sets list of associated tracks.
  inline void setTracks( const std::list<CsTrack*> &Tracks ) 
    { Tracks_ = Tracks; }

  //! Sets coordinates of vertex.
  inline void setVertex( double x, double y, double z )
    { x_ = x; y_ = y; z_ = z; }

  //! Sets type of vertex. true - primary, false - secondary.
  inline void setType( bool type ) { type_ = type; }

  //! Sets best primary vertex status
  inline void setBestVertex(bool status) { isBestVertex_ = status; }

 private:
  
  bool   type_;            //!< true - primary, false - secondary.

  double x_;               //!< vertex X coordinate (mm)
  double y_;               //!< vertex Y coordinate (mm)
  double z_;               //!< vertex Z coordinate (mm)

  double Chi2tot_;         //!< Chi-square of vertex
  std::vector<CLHEP::HepMatrix*> Cov_; //!< Covariance matrix

  std::list<CsTrack*> Tracks_;  //!< list of pointers to associated tracks
  std::map<CsTrack*,Cs3Vector> trkPar_; //!< map of track parameters at vertex.

  bool isBestVertex_;      //!< Best primary vertex, according to reigning vertex package

  #if USE_ObjectsCounter
  ObjectsCounter<CsVertex> objects_counter;
  #endif
};

#endif // CsVertex_h
