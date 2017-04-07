// $Id: CsMCTrack.h,v 1.8 2009/08/24 21:06:37 ybedfer Exp $

/*!
   \file    CsMCTrack.h
   \brief   Compass Monte Carlo Track Class.
   \author  Benigno
   \version $Revision: 1.8 $
   \date    $Date: 2009/08/24 21:06:37 $
*/
#ifndef CsMCTrack_h
#define CsMCTrack_h

#include "CsSTD.h"
#include <CLHEP/Vector/LorentzVector.h>
#include "CsMCParticle.h"
#include "CsMCVertex.h"
#include "CsMCHit.h"
#include "CsTrack.h"

class CsMCHit;
class CsMCVertex;

/*! \class CsMCTrack 
    \brief   Compass Montecarlo Track Class.

    Collect all needed Monte Carlo informations from Geant KINE bank. 
    Objects of CsMCTrack class are istantiated by CsGeant3 methods.
*/

class CsMCTrack{
 public:

  CsMCTrack(); //!< Default Constructor

  /*! \fn   CsMCTrack( int Gnum, double px, double py, double pz, 
     CsMCParticle particle, CsMCVertex& inVertex );
     \brief Constructor.
     \param Gnum Geant track identification number
     \param px X component of particle momentum (GeV)
     \param py Y component of particle momentum (GeV)
     \param pz Z component of particle momentum (GeV)
     \param particle Pointer to the associated CsMCParticle structure
     \param inVertex Pointer to the associated origin vertex  
  */
  CsMCTrack( int Gnum, double px, double py, double pz, 
	     CsMCParticle particle, CsMCVertex& inVertex );

  /*! \fn   CsMCTrack( int Gnum, Hep3Vector p, CsMCParticle particle, 
     CsMCVertex& inVertex );
     \brief Constructor.
     \param Gnum Geant identification number
     \param p Particle momentum (GeV)
     \param particle Pointer to the associated CsMCParticle structure
     \param inVertex Pointer to the associated origin vertex  
  */
  CsMCTrack( int Gnum, CLHEP::Hep3Vector p, 
	     CsMCParticle particle, CsMCVertex& inVertex );

  CsMCTrack( const CsMCTrack& ); //!< Copy Constructor

  CsMCTrack& operator=( const CsMCTrack& ); //!< Assignment Operator

  bool operator<( const CsMCTrack& ) const; //!< "less than" operator

  bool operator==( const CsMCTrack& ) const; //!< "equal to" operator
  
  int    getGnum() const; //!< Returns the Geant track identification number 

  CLHEP::HepLorentzVector getP() const; //!< Returns the particle 4-momentum (GeV)

  double getPX() const; //!< Returns the X component of particle momentum (GeV)

  double getPY() const; //!< Returns the Y component of particle momentum (GeV)

  double getPZ() const; //!< Returns the Z component of particle momentum (GeV)

  double getE() const; //!< Returns the particle energy (GeV)

  double getM() const; //!< Returns the particle mass (GeV)

  //! Returns the pointer to associated CsMCParticle object
  CsMCParticle* getParticle() const; 

  //! Returns the list of pointers to the associated CsMCHit objects
  std::list<CsMCHit*> getMCHits() const; 

  //! Returns the list of pointers to the daughter tracks
  std::list<CsMCTrack*> getOutTracks() const; 

  //! Returns the pointer to in origin CsMCVertex
  const CsMCVertex* getInVertex() const; 

  //! Returns the pointer to the list of vertices produced by the track
  std::list<CsMCVertex*> getOutVertices() const; 
  
  //! Add an Hit to the list of CsMCHit associated to this particle 
  void addMCHit( CsMCHit& MCHit ); 

  //! Add a track to the list of daugther tracks
  void addOutTrack( CsMCTrack& MCTrack );

  //! Add a vertex to the list of produced vertices
  void addOutVertex( CsMCVertex& MCVertex );

  //! Get list of associated reco'd, CsTrack, ID's
  std::list<int> getAssociatedTrackIDs() const {
    return( recoTrackIDs_ );
  }

  //! Add a CsTrack ID to the list of associated reco'd, CsTrack, ID's
  void addAssociatedTrackID( int trackID) {
    recoTrackIDs_.push_back( trackID );
  }

  //! Clear list of associated reco'd, CsTrack, ID's
  void clearAssociatedTrackID() {
    recoTrackIDs_.clear();
  }

 private:
  int Gnum_;
  CLHEP::HepLorentzVector p_;
  CsMCParticle particle_;
  std::list<CsMCHit*> MCHits_;
  CsMCVertex* inVertex_;
  std::list<CsMCTrack*> outTracks_;
  std::list<CsMCVertex*> outVertices_;
  std::list<int> recoTrackIDs_;
};

#endif // CsMCTrack_h
