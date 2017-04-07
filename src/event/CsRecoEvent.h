// $Id: CsRecoEvent.h,v 1.35 2010/08/27 16:41:15 tnagel Exp $

/*!
   \file    CsRecoEvent.h
   \brief   Compass Class CsRecoEvent
   \author  Benigno Gobbo
   \version $Revision: 1.35 $
   \date    $Date: 2010/08/27 16:41:15 $

   This will contain all reconstructed transient objects. 
*/

#ifndef CsRecoEvent_h
#define CsRecoEvent_h

#include "CsSTD.h"
#include "CsTrack.h"
#include "CsDigit.h"
#include "CsCluster.h"
#include "CsParticle.h"

class CsBeam;
namespace Reco{class CalorimeterParticle;}

/*!
   \class    CsRecoEvent
   \brief    Compass Class CsRecoEvent

   This will contain all reconstructed transient objects. 
*/

class CsRecoEvent {

 public:

  CsRecoEvent();                                  //!< The default constructor 

  ~CsRecoEvent();                                 //!< The default destructor 

  /*! \fn void clear() 
    \brief Clear the object private variables
  */
  void clear();

  /*! \fn void reset() 
    \brief Clear the object private tracks, beams, calos, particles and vertices
  */
  void reset();

  //========== Tracks corner =========== 

  /*! \fn void setTracks(const std::list<CsTrack*> &tracks ) 
    \brief Sets tracks to be the list of reconstructed tracks
    \param tracks The list of tracks to be added
  */
  void setTracks( const std::list<CsTrack*> &tracks );

  /*! \fn const std::list<CsTrack*> &getTracks(void) const;
    \brief Returns the list of pointers to the reconstructed tracks
  */
  const std::list<CsTrack*> &getTracks(void) const { return _tracks; }

  /*! \fn std::list<CsTrack*> &getTracks(void);
    \brief Returns the list of pointers to the reconstructed tracks
  */
  std::list<CsTrack*> &getTracks(void) { return _tracks; }

  /*! \fn int getTrackId(void) const; 
    \brief Returns the unique Id for reconstructed track numbering 
  */
  int getTrackId(void) const { return _nextTrkId++; }

  /*! \fn void clearTracks(void);
    \brief Clear the list of tracks (removes also the objects).
  */
  void clearTracks(void);

  //========== Digits corner =========== 

  /*! \fn const std::list<CsDigit*> &getDigits(void) const;
    \brief Returns the list of pointers to the built CsDigit objects.
  */
  const std::list<CsDigit*> &getDigits(void) const { return( _digits); }

  /*! \fn void addDigit( CsDigit& digit );
    \brief Add a digit to the list of built CsDigit objects.
  */
  void addDigit( CsDigit& digit )  { _digits.push_back( &digit ); }

  /*! \fn void clearDigits(void);
    \brief Clear the list of digits (removes also the objects).
  */
  void clearDigits(void);

  /*! void sortDigits(void);
    \brief Sorts the digits
  */
  void sortDigits(void);

  //========== Clusters corner =========== 

  /*! \fn const std::list<CsCluster*> getClusters(void) const;
    \brief Returns the list of pointers to the built CsCluster objects.
  */
  const std::list<CsCluster*> &getClusters(void) const { return _clusters; }

  /*! \fn void addCluster( CsCluster& cluster );
    \brief Add a cluster to the list of built CsCluster objects.
  */
  void addCluster( CsCluster& cluster ) { _clusters.push_back( &cluster ); }

  /*! \fn void subtractCluster(void);
    \brief Subtract the last cluster from the list of built CsCluster objects.
  */
  void subtractCluster(void);

  /*! \fn void clearClusters(void);
    \brief Clear the list of clusters (removes also the objects).
  */
  void clearClusters(void);

  /*! void sortClusters(void);
    \brief Sorts the clusters
  */
  void sortClusters(void);

  //========== Calorimeter Objects corner =========== 


  //! Sets the vector of reconstructed calorimeter objects
  void setCalObjsVector( const std::vector<Reco::CalorimeterParticle> calobjs ); 

  //! Returns the vector of reconstructed Calorimeter objects
  const std::vector<Reco::CalorimeterParticle*> &getCalObjs(void) const {
    return( _calObjs ); }

  //! Clear the vector of reconstructed Calorimeter objects
  void clearCalObjs( void );

  //! Add reconstructed calorimeter object
  void AddCalObject(Reco::CalorimeterParticle* p) {
   _calObjs.push_back(p); }

  //========== Particles corner =========== 

  //! Sets the list of pointers to reconstructed particles objects
  void setParticles( const std::vector<CsParticle*> &parts );
  
  /*! \fn const std::vector<CsParticle*> &getParticles(void) const;
    \brief Returns the list of pointers to the built CsParticle objects.
  */
  const std::vector<CsParticle*> &getParticles(void) const { return( _particles); }

  /*! \fn void addParticle( CsParticle& particle );
    \brief Add a particle to the list of built CsParticle objects.
  */
  void addParticle( CsParticle& particle ) 
     { _particles.push_back( &particle ); }

  /*! \fn void clearParticles(void);
    \brief Clear the vector of particles (removes also the objects).
  */
  void clearParticles(void);

  //========= Beam Reconstruction corner ===========

  //! Sets the list of pointers to reconstructed beam track objects
  void setBeamTracksList( std::list<CsBeam*> beamTracks ) 
    { _beamTracks = beamTracks; }

  //! Returns the list of pointers to reconstructed beam track objects
  const std::list<CsBeam*> getBeamTracksList(void) const
    { return( _beamTracks ); }

  //! Clear the list of pointers(!) to beam tracks objects
  void clearBeamTracksList(void) { _beamTracks.clear(); }

  //! Clear the list of beam tracks objects (to be used only reading DSTs)
  void clearBeamTracks(void);

  //========== Vertices corner =========== 

  /*! \fn void setVertices(const std::list<CsVertex*> &vertices ) 
    \brief Sets vertices to be the list of reconstructed vertices
    \param vertices The list of vertices to be added
  */
  void setVertices( const std::list<CsVertex*> &vertices );

  /*! \fn std::list<CsVertex*> &getVertices(void);
    \brief Returns the list of pointers to the reconstructed vertices
  */
  std::list<CsVertex*> &getVertices(void) { return _vertices; }

  //! returns the list of pointers to the reconstructed vertices
  const std::list<CsVertex*> &getVertices(void) const { return _vertices; }

  /*! \fn void clearVertices(void);
    \brief Clear the list of vertices (removes also the objects).
  */
  void clearVertices(void);

  // ========== trigger Time ==========
  
  //! sets the event trigger time  
  void setTriggerTime( const float time) { _triggerTime = time; }

  //! resets the event trigger time  
  void resetTriggerTime( void ) { _triggerTime = 0.0; }
  
  //! returns trigger time
  float getTriggerTime( void ) const { return( _triggerTime ); }
    
  // ========== scalers ==============
  // NOTA BENE: SC99P2 channels available only in 2002, cf. mapping.
  void setSC99P2_16(  int val ) { _scalers[0] = val; } //!< set SC99P2 16
  void setSC99P2_17(  int val ) { _scalers[1] = val; } //!< set SC99P2 17
  void setSC99P2_18(  int val ) { _scalers[2] = val; } //!< set SC99P2 18
  void setSC99P2_19(  int val ) { _scalers[3] = val; } //!< set SC99P2 19
  void setSC99P2_20(  int val ) { _scalers[4] = val; } //!< set SC99P2 20
  void setSC99P2_21(  int val ) { _scalers[5] = val; } //!< set SC99P2 21
  void setSC01P1_01(  int val ) { _scalers[6] = val; } //!< set SC01P1  1
  void setSC01P1_09(  int val ) { _scalers[7] = val; } //!< set SC01P1  9
  void setTimeInSpill( double val ) { _timeinspill = val; } //!< set time in spill
  void setFluxX( unsigned int val) { _incidentFlux[0] = val; }   //!< Set flux from FI0?X
  void setFluxY( unsigned int val) { _incidentFlux[1] = val; }   //!< Set flux from FI0?Y
  void addFluxX( unsigned int val) { _incidentFlux[0] += val; }  //!< Add flux from FI0?X
  void addFluxY( unsigned int val) { _incidentFlux[1] += val; }  //!< Add flux from FI0?Y

  int getSC99P2_16() const { return( _scalers[0] ); } //!< return SC99P2 16
  int getSC99P2_17() const { return( _scalers[1] ); } //!< return SC99P2 17
  int getSC99P2_18() const { return( _scalers[2] ); } //!< return SC99P2 18
  int getSC99P2_19() const { return( _scalers[3] ); } //!< return SC99P2 19
  int getSC99P2_20() const { return( _scalers[4] ); } //!< return SC99P2 20
  int getSC99P2_21() const { return( _scalers[5] ); } //!< return SC99P2 21
  int getSC01P1_01() const { return( _scalers[6] ); } //!< return SC01P1  1
  int getSC01P1_09() const { return( _scalers[7] ); } //!< return SC01P1  9
  double getTimeInSpill() const { return  _timeinspill; } //!< return time in spill (unit should be seconds, see CsEvent.cc, where this is set...)
  unsigned int getFluxX() const { return _incidentFlux[0]; } //!< Return flux from FluxScalers
  unsigned int getFluxY() const { return _incidentFlux[1]; } //!< Return flux from FluxScalers

  //!  sets all scalers values
  void setScalers( const int* scalers );

  //! return pointer to scalers
  const int* getScalers(void) const { return( _scalers ); } 

  // ========== scalers ==============

  //! set position-th datum in event header
  void setEventHeader( const unsigned int datum, const unsigned int position );

  //! set data in event header
  void setEventHeader( const unsigned int* data);

  //! get pointer to event header data
  const unsigned int* getEventHeader() const { return( _eventHeader ); } 

  // ========== hodoscopes ============

  //! add a hodoscope datum
  void addHodoDatum( unsigned int datum ) { _hodoData.push_back( datum ); }

  //! get hodoscopes data vector
  const std::vector<unsigned int> &getHodoData(void) const { return( _hodoData ); }

 private:

  static int       _nextTrkId;    //!< Used for track numbering
  std::list<CsTrack*>   _tracks;       //!< List of reconstructed track
  std::list<CsDigit*>   _digits;       //!< List of event digits
  std::list<CsCluster*> _clusters;     //!< List of event clusters
  std::vector<Reco::CalorimeterParticle*> _calObjs; //!< Calorimetric objects
  std::list<CsBeam*>    _beamTracks;      //!< List of reconstructed beam track
  std::vector<CsParticle*> _particles;    //!< List of event particles
  std::list<CsVertex*>  _vertices;        //!< List of reconstructed vertices
  float            _triggerTime;     //!< Trigger Time
  int              _scalers[8];      //!< Scalers data  
  unsigned int     _incidentFlux[2]; //!< Flux measured by FI02X and Y.
  double           _timeinspill;     //!< Event time in spill
  unsigned int     _eventHeader[10]; //!< Subset of Date Event Header Data 
  std::vector<unsigned int> _hodoData;    //!< Data from hodoscopes 

  struct           _sortDigits;   //!< To sort the list of pointers
  struct           _sortClusters; //!< To sort the list of pointers

};

#endif //CsRecoEvent_h
