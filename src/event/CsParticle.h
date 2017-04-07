// $Id: CsParticle.h,v 1.12 2010/06/17 11:40:30 tnagel Exp $

/*!
   \file    CsParticle.h
   \brief   Compass Particle Class.
   \author  Benigno Gobbo
   \version $Revision: 1.12 $
   \date    $Date: 2010/06/17 11:40:30 $

   \par History:
   20001212 V 1.0
*/

#ifndef CsParticle_h
#define CsParticle_h

#include "coral_config.h"

#include <vector>
#include <string>

#include "CsTrack.h"

#if USE_ObjectsCounter
#include "ObjectsCounter.h"
#endif

namespace Reco{class CalorimeterParticle;}

/*! \class CsParticle

    \brief Compass Particle Class.

    this is the container of reconstructed data associated to particle
*/

class CsParticle {

 public:

  //! Special track types (for vertex purpose).
  enum {ORDINARY,SPECIAL};

  CsParticle();                                //!< Default Constructor

  CsParticle( CsTrack* track );

  CsParticle( Reco::CalorimeterParticle* calobj );

  CsParticle( const CsParticle& );             //!< Copy Constructor

  CsParticle& operator=( const CsParticle& );  //!< Assign Operator

  int getCharge() { return( _charge ); }
  
  std::string  getName() { return( _name ); }; //!< Returns the PDG particle name 

  //Type getType(void) const {return _type;}
  int getType(void) const {return _type;}

  const CsTrack* getTrack() const { return( _track ); }

  const std::vector<Reco::CalorimeterParticle*> &getCalObjects() 
    { return( _calObjs ); }

  void addCalobj( Reco::CalorimeterParticle* calobj ) 
    { _calObjs.push_back( calobj ); }

  void setName( std::string name ) {_name = name; }; //!< Sets the PDG particle name 

  void setType(int t) { _type = t; }

 private:

  int          _charge; //!< Charge of the particle if available. -999 if unknown  

  std::string  _name;   //!< PDG particle name

  CsTrack*     _track;  //!< Associated reconstructed charged track. NULL if not available   

  int          _type;   //!< The track's type. Default particle type is 'ORDINARY' (needed for vertexing)

  std::vector<Reco::CalorimeterParticle*> _calObjs;//!< Associated calobjects if any.
  
  #if USE_ObjectsCounter
  ObjectsCounter<CsParticle> objects_counter;
  #endif
};

#endif // CsParticle_h
