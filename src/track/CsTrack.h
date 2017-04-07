// $Id: CsTrack.h,v 1.45 2010/06/23 09:49:58 tnagel Exp $

/*!
   \file    CsTrack.h
   \brief   Compass Track Class.
   \author  Benigno Gobbo
   \version $Revision: 1.45 $
   \date    $Date: 2010/06/23 09:49:58 $

   \par History:
   20000327 Added Vertices <br>
   20001218 Added RICH1 output
*/
#ifndef CsTrack_h
#define CsTrack_h

#include "coral_config.h"

# include <list>
# include <vector>

#include "CsCluster.h"
#include "CsHelix.h"
#include "CsZone.h"
#include "CsVertex.h"
#include "CsMCTrack.h"

#if USE_ObjectsCounter
#include "ObjectsCounter.h"
#endif

//#define CSTRACK_RICHDATASIZE (15)
#define CSTRACK_RICHDATASIZE (23) 
// W/ 2008.04 setup on, hit map size need be = 13...
#define CSTRACK_MAPSIZE (13)

class CsVertex;

/*! \class CsTrack

    \brief Compass Track Class.

    This class describes a track. A track consists of a vector of helices
    describing the track parameters in different points along the
    trajectory and the collection of clusters associated to the track.
    
*/

class CsTrack {

 public:

  CsTrack();                              //!< Default Constructor

  /*! \fn CsTrack( std::list<CsCluster*> clusters )
    \brief Constructor
    \param clusters collection of the associated clusters
    \param zones    zone from where track belong to
  */
  CsTrack( const std::list<CsCluster*> &clusters );

  /*! \fn CsTrack( std::list<CsCluster*> clusters, std::list<CsZone*> zones ) 
    \brief Constructor
    \param clusters collection of the associated clusters
    \param zones    zone from where track belong to
  */
  CsTrack( const std::list<CsCluster*> &clusters, const std::list<CsZone*> &zones );

  /*! \fn CsTrack( std::vector<CsHelix> helices, std::list<CsCluster*> clusters, std::list<CsZone*> zones )
    \brief Constructor
    \param helices  the vector of helices describing the track parameters
           in different points along the trajectory 
    \param clusters collection of the associated clusters
    \param zones    zone from where track belong to
  */
  CsTrack( const std::vector<CsHelix> &helices, const std::list<CsCluster*> &clusters,
	   const std::list<CsZone*> &zones );

  /*! \fn CsTrack( std::vector<CsHelix> helices, std::list<CsCluster*> clusters, std::list<CsZone*> zones, double chi2 )
    \brief Constructor
    \param helices  the vector of helices describing the track parameters
           in different points along the trajectory 
    \param clusters collection of the associated clusters
    \param zones    zone from where track belong to
    \param chi2     the fit chi2
  */
  CsTrack( const std::vector<CsHelix> &helices, const std::list<CsCluster*> &clusters,
	   const std::list<CsZone*> &zones, double chi2 );

  CsTrack( const CsTrack& );             //!< Copy Constructor

  virtual ~CsTrack();

  CsTrack& operator=( const CsTrack& );  //!< Assign Operator

  bool operator==( const CsTrack& );     //!< "Equal to" Operator
  virtual bool IsBeamTrack() const { return false; }
  

  /*! \fn std::vector<CsHelix> getHelices()
     \brief Returns the vector of helices calculated in different points of
     the track trajectory
   */
  const std::vector<CsHelix> &getHelices() const { return(_helices); }

  /*! \fn std::list<CsCluster*> getClusters()
    \brief Returns the collection of associated clusters.
  */
  const std::list<CsCluster*> &getClusters() const { return(_clusters); }

  /*! \fn std::list<CsZone*> getZones()
    \brief Returns the collection associated zones.
  */
  const std::list<CsZone*> &getZones() const { return(_zones); }

  /*! \fn inline std::list<CsVertex*> getVertices() const;
    \brief Returns the pointer to the track associated vertices.
  */
  const std::list<CsVertex*> &getVertices() const { return( _vertices ); }

  /*! \fn CsVertex* getFirstVertex() const;
    \brief Returns the pointer to the first associated vertex (if any),
    Returns NULL pointer otherwise.
  */
  CsVertex* getFirstVertex() const; 

  /*! \fn int getId()
    \brief Return the track id
  */
  int getId() const { return( _id ); }

  /*! \fn inline double getChi2()
    \brief Returns the Chi2 of the fitting algorithm.
  */
  double getChi2() const { return(_chi2); }

  /*! \fn inline double getXX0() 
    \brief Returns the number of rad lengths which track crossed. 
  */
  double getXX0() const { return(_XX0); }

  /*! \fn void addHelix( CsHelix helix )
    \brief Adds an helix to the track  
    \param helix The helix to be added
  */
  virtual void addHelix( const CsHelix &helix );

  /* \fn inline setHelices( std::vector<CsHelix> helices )
     \brief Sets the track helices vector
     \param helices the vector of helices to be set
  */
  virtual void setHelices( const std::vector<CsHelix> &helices ) { _helices = helices; }

  /*! \fn void addCluster( CsCluster& cluster )
    \brief Adds a cluster to the track
    \param The reference to the cluster to be added
  */
  virtual void addCluster( CsCluster& cluster );

  /* \fn inline void setClusters( std::list<CsClusters*> clusters )
     \brief Sets the track cluster pointers list
     \param clusters The list of cluster pointers to be set
  */
  virtual void setClusters( const std::list<CsCluster*> &clusters ) {_clusters=clusters;}

  /*! \fn void addZone( CsZone& zone )
    \brief Adds a zone to the track
    \param The reference to the zone to be added
  */
  virtual void addZone( CsZone& zone );

  /*! \fn inline void setZones( std::list<CsZones*> zones )
     \brief Sets the track zone pointers list
     \param zones The list of zone pointers to be set
  */
  virtual void setZones( const std::list<CsZone*> &zones ) { _zones = zones; }

  /*! \fn inline void setChi2( double chi2 )
    \brief Set the fit Chi2 of the track
    \param chi2 The Chi2 to fit
  */
  virtual void setChi2( double chi2 ) { _chi2 = chi2; }

  /*! \fn inline void setXX0( double XX0 ) 
    \brief Set the X/X0 numbber of rad lengths which track crossed
    \param chi2 X/X0
  */
  virtual void setXX0( double XX0 ) { _XX0 = XX0; }

  /*! \fn inline void addVertex( CsVertex* vertex );
     \brief Add a vertex to the list of track associated vertices
     \param vertex The pointer to the vertex. 
  */
  virtual void addVertex( CsVertex* vertex ) { _vertices.push_back(vertex); }
  
  /*! \fn void rmVertex( CsVertex* vertex );
     \brief Remove the vertex from the list of track associated vertices
     \param vertex The pointer to the vertex. 
  */
  virtual void rmVertex( CsVertex* vertex ) { _vertices.remove(vertex); }
  
  
  /*! \fn void setMeanTime( double time, double error );
     \brief Sets the track mean time and error
     \param time The value to be set
     \param error The error on time (defaults to -1)
  */
  virtual void setMeanTime( double time, double error=-1) 
    { _time = time; _timeError = error; _hasTime = true; }
  //! Returns \c true if mean time was set for this track
  virtual bool hasMeanTime() const { return( _hasTime ); }
  //! Returns the track mean time, if available. If not an exception arrises.
  double getMeanTime() const;
  //! Returns the track mean time error, if available. If not an exception arrises.
  double getMeanTimeError() const;

  /*! \fn void addCEDARInfo( );
     \brief add vector of CEDAR PID info
  */
  void addCEDARInfo();

  const std::vector<float> &getCEDARInfo() { return _cedarInfo;}

  /*! \fn void setRich1Probs(const double * probs );
     \brief Set the Rich1 probabilities array for this track. 
     \param probs Rich1 output. It is expected to a dim. CSTRACK_RICHDATASIZE array of doubles
  */
  virtual void setRich1Probs(const double * probs );

  //! Returns \c true if Rich1 probabilities were set for this track
  bool hasRich1Probs() const { return( _hasRich1 ); }  

  //! Returns the pointer to rich1 probabilities. Exception if not available.
  const double* getRich1Probs() const;

  //! Returns the background likelihood (from Rich1 reconstruction)
  double BkgLikelihood();
  //! Returns the Pion likelihood (from Rich1 reconstruction)
  double PionLikelihood();
  //! Returns the Kaon likelihood (from Rich1 reconstruction)
  double KaonLikelihood();
  //! Returns the Proton likelihood (from Rich1 reconstruction)
  double ProtonLikelihood();
  //! Returns the Electron likelihood (from Rich1 reconstruction)
  double ElectronLikelihood( void );
  //! Returns the Muon likelihood (from Rich1 reconstruction)
  double MuonLikelihood( void );

  //! For backward compatibility
  double BkgProb() { return( BkgLikelihood() ); }
   //! For backward compatibility
  double PionProb() { return( PionLikelihood() ); }
  //! For backward compatibility
  double KaonProb() { return( KaonLikelihood() ); }
  //! For backward compatibility
  double ProtonProb() { return( ProtonLikelihood() ); }

  //! Returns the Pion Likelihood derivative (from Rich1 reconstruction)
  double PionDLikeDIndex();
  //! Returns the Kaon Likelihood derivative (from Rich1 reconstruction)
  double KaonDLikeDIndex();
  //! Returns the Proton Likelihood derivative (from Rich1 reconstruction)
  double ProtonDLikeDIndex();
  //! Returns the Electron likelihood derivative (from Rich1 reconstruction)
  double ElectronDLikeDIndex( void );
  //! Returns the Muon likelihood derivative (from Rich1 reconstruction)
  double MuonDLikeDIndex( void );

  //! Returns the theta angle from max like (from Rich1 reconstruction)
  double Rich1ThetaMaxLike();

  //! Returns the theta angle of ring (from Rich1 reconstruction)
  double Rich1Theta();

  //! Returns the number of ring photons (from Rich1 reconstruction)
  int    Rich1NbPhot();

  //! Returns the theta from ring fit  (from Rich1 reconstruction)
  double Rich1ThetaFit();

  //! Returns the Reconstructed Ring Chi2 (from Rich1 reconstruction)
  double Rich1RingChi2();

  //! Returns the Pion Chi2 (from Rich1 reconstruction)
  double PionChi2();
  //! Returns the Kaon Chi2 (from Rich1 reconstruction)
  double KaonChi2();
  //! Returns the Proton Chi2 (from Rich1 reconstruction)
  double ProtonChi2();
  //! Returns the Electron Chi2 (from Rich1 reconstruction)
  double ElectronChi2( void );
  //! Returns the Muon Chi2 (from Rich1 reconstruction)
  double MuonChi2( void );


  //! Set the expected fired detectors bitmap.
  bool setExpectedDetsBitmap( const unsigned int* bitmap ); 
  //! Set the fired detectors bitmap.
  bool setFiredDetsBitmap( const unsigned int* bitmap ); 
  //! Returns the pointer to the expected fired detectors bitmap.
  const unsigned int* getExpectedDetsBitmap() const;
  //! Returns the pointer to the fired detectors bitmap.
  const unsigned int* getFiredDetsBitmap() const;

  //! Returns number of associated clusters.
  int getNumberOfAssociatedClusters() const;

  //! Get NDFs = #hits + #pixel_hits. (Caveat: NDFs is not set in the constructor, cf. "CsTrack::setNDFs".)
  int getNDFs() const { return _ndfs; }
  //! NDFs is set = #hits in the constructor, whereas it's expected to be += #pixel_hits. (Note: It would be better that the constructor sets the correct value in the first place, by deriving it from the list of clusters.)   
  void setNDFs(int ndfs) { _ndfs = ndfs; }

  //! Set shower attributes.
  void setShower(int nHits, int shower) {
    _nShowerHits = nHits; _hasShower = shower;
  }
  //! Return shower flag.
  int hasShower() const { return _hasShower; }

  //! Set associated MC track.
  void setAssociatedMCTrack(CsMCTrack *track) { _associatedMCTrack = track; }
  //! Get associated MC track.
  CsMCTrack *getAssociatedMCTrack() const { return _associatedMCTrack; }

  //! Set associated reconstructed track.
  void setAssociatedTrack(int id) { _associatedTrack = id; }
  //! Get associated reconstructed track. Defaults to -1. Association (so far) concerns the non-interacting, beam tracks which this->_associatedTrack are presumed to be the continuation of.
  int getAssociatedTrack() const { return _associatedTrack; }

  /// \return downstream-most helix which still is upstream of z, if no helix
  /// exists upstream of z, return helix closest to z, return NULL if track
  /// has no helix
  const CsHelix* getHelixUpstreamOf( double z ) const;

 private:

  int    _id;                                  //!< Track identifier
  double _chi2;                                //!< Chi2 of the fit 
  int    _ndfs;                                //!< NDFs = #hits + #pixel_hits
  std::vector<CsHelix>   _helices;             //!< Helices
  std::list<CsCluster*>  _clusters;            //!< Associated clusters  
  std::list<CsZone*>     _zones;               //!< Validity zones
  std::list<CsVertex*>   _vertices;            //!< Associated Vertices
  double _time;                                //!< mean time of track
  double _timeError;                           //!< Track's mean time error
  bool   _hasTime;                             //!< true if mean time is set
  double _rich1Probs[CSTRACK_RICHDATASIZE];    //!< RICH response for this track
  bool   _hasRich1;                            //!< true if Rich data is set
  std::vector<float> _cedarInfo;              //!< CEDAR response for this track
  bool   _hasCEDAR;                            //!< true if CEDAR data is set
  double _XX0;                                 //!< Number of rad lengths.
  unsigned int  _expectedDets[CSTRACK_MAPSIZE];//!< Bit map of dets expected to be fired
  unsigned int  _firedDets   [CSTRACK_MAPSIZE];//!< Bit map of fired dets
  CsMCTrack     *_associatedMCTrack;           //!< MC associated, typically by Tracking Package
  int           _associatedTrack;              //!< Associated reconstructed Track. Association (so far) concerns the non-interacting, beam tracks which this->_associatedTrack are presumed to be the continuation of.
  int    _nShowerHits;                         //!< #hits collected in preshower
  int    _hasShower;                           //!< Shower flag: =1: #hits > cut, =2: &= d#hits/dZ > cut

  #if USE_ObjectsCounter
  ObjectsCounter<CsTrack> objects_counter;
  #endif
};

#endif // CsTrack_h
