// $Id: CsEvent.h,v 1.98 2010/08/17 12:25:40 tnagel Exp $

/*!
   \file    CsEvent.h
   \brief   Compass Event Class.
   \author  Benigno Gobbo
   \date    $Date: 2010/08/17 12:25:40 $
   \version $Revision: 1.98 $
   
   \par History:
   19990426 First release.<br>
   19990518 Added Raw Data readout.<br>
   19990524 Added Equipment ReadOut.<br>
   19990726 Added CsMagInfo readout.<br>
   19990804 Added Read MC also from Ntuples.<br>
   19990906 Changes to follow new CsGeom class.<br>
   19991203 Added CsRecoEvent class.<br>
   20000128 Added pointer to CsStore.<br>
   20000222 Added the possibility to read only N events, and to skip N events.<br>
   20000622 Some operations moved to CsGean3. <br>
   20000622 Simple reconstruction schema added. <br>
   20010102 In reco schema: beam, calos, rich. <br> 
   20030325 Added hits rebuilding in DST readout <br>
   20030409 gcc 3 compliant <br>
   20030626 skip event in case of decoding troubles <br>
*/

#ifndef CsEvent_h
#define CsEvent_h

#include "coral_config.h"

#include <list>
#include <utility>

#include "CsSTD.h"
#include "CsTime.h"
#include "CsTypes.h"
#include "CsHistograms.h"
#include "CsRecoEvent.h"
#include "CsEndOfJob.h"

class CsTrack;
class CsMCTrack;
class CsMCVertex;
class CsVertex;
class CsMCHit;
class CsTrkPrepattern;
class CsTrkBridging;
class CsTrkFitting;
class CsVrtPattern;
class CsVrtFitting;
class CsStore;
class CsGeant3;
class CsParticle;
class CsBeam;
namespace Reco{class CalorimeterParticle;}
namespace CS {class DaqEvent;}

/*! \class CsEvent 
    \brief   Compass Event Class.
    It interfaces to CsGeant3 for Monte Carlo simulated Events or 
    CsRawEvents for real data access (RAW data only, at present).
    It also istantiates the CsGeom class.
*/

class CsEvent : public CsEndOfJob
{
 public:

  /*! \fn static CsEvent* Instance();
    \brief singleton instantiation.
  */
  static CsEvent* Instance();

  //! Clear event
  void clear(void);

  /*! \fn bool getNextEvent();
    \brief Read next event and fill the related objects. Returns \c true
    of the event is read, \c false otherwise.
  */
  bool getNextEvent();

  /*! \fn CsTime getEventTime (void) const;
    \brief Returns the Event time. Works only for real data events: on
    Monte Carlo is not yet implemented.
  */
  CsTime getEventTime(void) const;

  /*! \fn bool isAMonteCarloEvent() const;
    \brief Returns \c true if the current event is a Monte Carlo one
  */
  bool isAMonteCarloEvent() const { return( _MCEvent); }

  // The following four methods apply only for Monte Carlo Events

  /*! \fn const std::list<CsMCTrack*> &getMCTracks(void) const;
    \brief Returns the list of pointers to the current event Monte Carlo true
    tracks.
  */
  const std::list <CsMCTrack*>  &getMCTracks(void) const;

  /*! \fn const std::list<CsMCVertex*> &getMCVertices(void) const;
    \brief Returns the list of pointers to the current event Monte Carlo true
    vertices.
  */
  const std::list <CsMCVertex*> &getMCVertices(void) const;

  /*! \fn const std::list<CsMCHit*> &getMCHits(void) const;
    \brief Returns the list of pointers to the current event Monte Carlo true
    hits.
  */
  const std::list <CsMCHit*> &getMCHits(void) const;

  /*! \fn bool isADataEvent(void) const;
    \brief Returns \c true if the current event is a Read Data one
  */
  bool isADataEvent(void) const { return _DTEvent; }

  // The following 13 methods applie only for Real Data Events

  /*! \fn unsigned int getEventSize(void) const;
    \brief Returns Event size in bytes.
  */
  unsigned int getEventSize(void) const;  

  /*! \fn unsigned int getRunNumber(void) const;
    \brief Returns the Run Number.
  */
  unsigned int getRunNumber(void) const;

  /*! \fn unsigned int  getEventNumberInRun(void) const;
    \brief Returns the Event Number in Run.  
  */
  unsigned int getEventNumberInRun(void) const;

  /*! \fn unsigned int  getBurstNumber(void) const;
    \brief In Data Events returns the Burst Number. It is 0 if MC events.
  */
  unsigned int getBurstNumber (void) const;

  /*! \fn unsigned int  getEventNumberInBurst(void) const;
    \brief In Data Events returns the Event Number in Burst. It is 0 if MC 
    events.
  */
  unsigned int getEventNumberInBurst(void) const;

  /*! \fn unsigned int getTriggerNumber(void) const;
    \brief Returns the trigger number.
  */
  unsigned int getTriggerNumber(void) const;

  /*! \fn pair<unsigned int,unsigned int>  getTime(void) const;
    \brief Returns event time pair (seconds, useconds)
  */
  std::pair<unsigned int,unsigned int> getTime (void) const;

  /*! \fn unsigned int  getTriggerMask(void) const;
    \brief In Data Events returns the Event Trigger Mask. It is 0 if MC 
    events.
  */
  unsigned int getTriggerMask(void) const;

  /*! \fn unsigned int  getErrorCode(void) const;
    \brief In Data Events returns the Error Code. It is 0 if MC 
    events.
  */
  unsigned int getErrorCode (void) const;

  /*! \fn const std::list<CsTrack*> &getTracks(void) const
    \brief Returns the list of reconstructed tracks.
  */
  const std::list<CsTrack*> &getTracks(void) const 
    { return _recoEvent.getTracks(); }

  /*! \fn std::list<CsTrack*> &getTracks(void)
    \brief Returns the list of reconstructed tracks.
  */
  std::list<CsTrack*> &getTracks(void) 
    { return _recoEvent.getTracks(); }

  /*! \fn void setTracks( const std::list<CsTrack*> &tracks )
    \brief Set the list of tracks to the CsRecoEvent one
    \param tracks The list of tracks to be set
  */
  void setTracks( const std::list<CsTrack*> &tracks ) 
    { _recoEvent.setTracks(tracks); }

  /*! \fn const std::list<CsBeam*> getBeam(void) 
    \brief Returns the list of reconstructed beam-objects.
  */
  const std::list<CsBeam*> getBeam(void)
    { return _recoEvent.getBeamTracksList(); }

  /*! \fn std::list<CsVertex*> &getVertices(void)
    \brief Returns the list of reconstructed vertices.
  */
  std::list<CsVertex*> &getVertices(void) 
    { return _recoEvent.getVertices(); }

  /// \returns the list of reconstructed vertices.
  const std::list<CsVertex*> &getVertices(void) const
    { return _recoEvent.getVertices(); }

  /*! \fn void setVertices( const std::list<CsVertex*> &vertices )
    \brief Set the list of vertices to the CsRecoEvent one
    \param vertices The list of vertices to be set
  */
  void setVertices( const std::list<CsVertex*> &vertices ) 
    { _recoEvent.setVertices(vertices); }

  /*! \fn void setParticles( const std::vector<CsParticle*> &parts )
    \brief Set the list of particles to the CsRecoEvent one
    \param parts The vector of particles to be set
  */
  void setParticles( const std::vector<CsParticle*> &parts ) 
    { _recoEvent.setParticles(parts); }

  /*! \fn int getTrackId(void) const 
    \brief Returns the unique Id for reconstructed track numbering 
  */
  int getTrackId(void) const {return _recoEvent.getTrackId();}

  /*! \fn const std::list<CsDigit*> getDigits(void) const;
    \brief Returns the list of pointers to the built CsDigit objects.
  */
  const std::list<CsDigit*> &getDigits (void) const 
    { return _recoEvent.getDigits(); }

  /*! \fn void addDigit( CsDigit& digit );
    \brief Add a  the list of built CsDigit objects.
  */
  void addDigit (CsDigit& digit) { _recoEvent.addDigit(digit); }

  /*! \fn const std::list<CsCluster*> getClusters(void) const;
    \brief Returns the list of pointers to the built CsCluster objects.
  */
  const std::list<CsCluster*> &getClusters (void) const 
    { return _recoEvent.getClusters(); }

  /*! \fn void sortClusters(void);
    \brief Sort the list of CsCluster objects.
  */
  void sortClusters(void) {_recoEvent.sortClusters();}

  /*! \fn void addCluster( CsCluster& cluster );
    \brief Add a cluster to the list of built CsCluster objects.
  */
  void addCluster( CsCluster& cluster ) { _recoEvent.addCluster(cluster); }

  /*! \fn void subtractCluster(void);
    \brief Subtract the last cluster from the list of built CsCluster objects.
  */
  void subtractCluster(void) { _recoEvent.subtractCluster(); }

  /*! \fn unsigned int getNumberOfEvents(void) const;
    \brief Returns the number of read events 
  */
  unsigned int getNumberOfEvents(void) const {return _nEvents;}

  /*! \fn unsigned int getNumberOfProcessedEvents(void) const;
    \brief Returns the number of processed events 
  */
  unsigned int getNumberOfProcessedEvents(void) const {return _nSelEvents;}

  unsigned int getNMuPerSpill() const {return  _nmuinspillsofar;}

  //! For MCs returns the clustering mode
  inline std::string getClusterizationType() { return( _MCClusteringMode ); }
  
  const CS::DaqEvent& getDaqEvent   (void) const;

  //! Get real data Trigger Time (common to ALL detectors). Dummy value in MC.
  float getTriggerTime(void) const { return( _recoEvent.getTriggerTime() ); }

  //! Get TCS phase time. Returns -1000000 in case of error. Dummy value in MC.
  double getTCSPhaseTime(void) const { return( _TCSPhaseTime ); }

  void setTCSPhaseTime ( double tcs_phase_time ) {_TCSPhaseTime = tcs_phase_time ; }

  /*! \fn getExtraTimeWidth();
    \brief Returns the extra time width to be added to detector time windows
    in order to account for the trigger jitter.
  */
  inline double getExtraTimeWidth(void) const { return _extraTimeWidth; }

  /*! \fn getTriggerMCOffset();
    \brief Returns the time offset of the trigger w.r.t. the event for current MC event (in ns)
  */
  inline double getTriggerMCOffset(void) const { return _triggerMCOffset; }

  /*! \fn getAlterMasterTrigger();
    \brief Returns the alternative master trigger. The number of the corresponding bit is returned, or -1 if no alternative is defined.
  */
  inline int getAlterMasterTrigger(void) const { return _alterMasterTrigger; }

  /*! \fn getTriggerTimeCorr();
    \brief Returns the offset of the alternative master trigger w.r.t. the master trigger (expected to be !=0 when an alternative master trigger is indeed defined)
  */
  inline double getTriggerTimeCorr(void) const { return _triggerTimeCorr; }

  /*! \fn getReTrackT0();
    \brief Returns the T0 to be used in 2nd iteration of tracking and reconstruction. Is defined by the 1st iteration.
  */
  inline double getReTrackT0(void) const { return _reTrackT0; }

  /*! \fn reTrackingON();
    \brief Returns \c true when in the 2nd iteration of tracking and reconstruction, w/ T0 defined by 1st iteration.
  */
  inline bool reTrackingON(void) const { return _reTrackingON; }

  //! Gives the polarization of up and down cells. Returns false if no such information.
  bool getPolarization(std::pair<double,double> &pol) const;

  //! Returns the vector of reconstructed Calorimeter objects
  const std::vector<Reco::CalorimeterParticle*> &getCalObjs(void) const
    { return( _recoEvent.getCalObjs() ); }

  //! Add reconstructed calorimeter object
  void AddCalObject(Reco::CalorimeterParticle* p) {
    _recoEvent.AddCalObject(p); }

  //! Returns the vector of reconstructed particles
  const std::vector<CsParticle*>& getParticles(void) const
    {return( _recoEvent.getParticles() ); }

  //! end of event method
  bool end();

  //! dump reco event data
  void testDump(bool);

  //! return detector map
  std::map<CsDetector*,int,std::less<CsDetector*> > getDetectorMap() {
    return( _detmap ); }

  //! return detector vector
  std::vector<CsDetector*> getDetectorVector() { return( _detvec ); }

  //! Open an out stream
  int openRawEventOutputStream( const std::string name );

  //! Close an out stream
  void closeRawEventOutputStream( const int stream );

  //! Close all out streams
  void closeAllRawEventOutputStream( void );

  //! Write a raw event to a stream
  bool outputRawEventToStream( const int stream );

  //! Get # of opened out streams
  int getNOutStreams() const { return _nOutStreams; }

  //! returns DaqdataDecoding digits
  const CS::Chip::Digits &getChipDigits(void) const;

  //! returns SC99P2 16
  int getSC99P2_16() const { return( _recoEvent.getSC99P2_16() ); } 
  //! returns SC99P2 17
  int getSC99P2_17() const { return( _recoEvent.getSC99P2_17() ); } 
  //! returns SC99P2 18
  int getSC99P2_18() const { return( _recoEvent.getSC99P2_18() ); } 
  //! returns SC99P2 19
  int getSC99P2_19() const { return( _recoEvent.getSC99P2_19() ); } 
  //! returns SC99P2 20
  int getSC99P2_20() const { return( _recoEvent.getSC99P2_20() ); } 
  //! returns SC99P2 21
  int getSC99P2_21() const { return( _recoEvent.getSC99P2_21() ); } 
  //! returns SC01P1_01
  int getSC01P1_01() const { return( _recoEvent.getSC01P1_01() ); } 
  //! returns SC01P1_09
  int getSC01P1_09() const { return( _recoEvent.getSC01P1_09() ); } 
  //! returns time in spill, in seconds
  double getTimeInSpill() const { return _recoEvent.getTimeInSpill(); } 
  //! set all scalers data
  void setScalers( const int* scalers ) { _recoEvent.setScalers( scalers ); } 

  //! returns hodoscopes data size
  unsigned int getHodoDataSize(void) {return(_recoEvent.getHodoData().size());}

  //! returns ith hodoscope datum (tbname, channel, time)
  bool getHodoDatum( const unsigned int i, std::string& tbname, int& channel, 
		     double& time );

  //! rebuild digits and hits after DST download
  bool rebuildDigitsAndHits( void );

  CsStore* GetStore() { return _store; }

 protected:

  CsEvent();          //!< The Protected Singleton Constructor
  

 private:

  //! Set the tracking packages
  void _setTrackingPackages();

  //! Set the tracking packages
  void _setVertexPackages();

  //! Reconstruction schema handler
  bool _reconstructionSchemaHandler();

  //! Reconstruction schema for pre processing (not realy reconstruction)
  bool _reconstructionSchemaPP();

  //! The first reconstruction schema
  bool _reconstructionSchema001();

  //! The second reconstruction schema
  bool _reconstructionSchema002();

  //! Decode an event
  bool _decode (void);
 
  //! MC-specific part of decoding
  void _decodeMC (void);
 
  //! RD-specific part of decoding
  bool _decodeRD (void);

  //! clusterize event
  bool _clusterize(void);

  void _mkClustersFromHits( int mode );

  //! upload to event database (DST production)
  void _upload(void);

  //! download from database (DSR readout)
  bool _download(void);

  //! get next read data event depending on access requests...
  bool _getNextDataEvent(void);

	//!read TCSphase jummp corrections
	void _readTcsCalibration(time_t timePoint);



  static CsEvent* _instance; //!< The singleton static attribute 
  unsigned int _nEvents;     //!< Number of read events
  unsigned int _nSelEvents;  //!< Number of selected events
  unsigned int _maxEvents;   //!< Max number of events to read
  unsigned int _skipEvents;  //!< Number of events to skip
  unsigned int _maxConsecutiveSkips;  //!< Max. # of consecutive events w/ decoding error that we tolerate (and hence skip) before giving up. 

  bool _DTEvent;            //!< \c true if Real Data Event
  CsStore* _store;          //!< pointer to the CsStore class

  CsRecoEvent _recoEvent;   //!< pointer to the current CsRecoEvent object 
  int         _thisRun;     //!< current run number
  int         _previousRun; //!< previous run number

  bool _MCEvent;            //!< \c true if Monte Carlo Event
  CsGeant3* _geant3MC;      //!< pointer to the CsGeant3 singleton
  double _triggerMCOffset;  //!< Offset of trigger time w.r.t. event time for current MC event (in ns)
  std::vector<double> _triggerMCJitters; //!< Jitters for the various bits of the trigger pattern (in ns)
  unsigned int _triggerMCDelayed;   //!< Patter of triggers which are delayed by 5ns when it comes to determine the overall trigger time

  // Decoding and reconstruction options
  bool _decoding;               //!< \c true if raw events decoding.
  bool _MCDecodingExact;        //!< \c true if trivial digit-hit link.
  bool _clustering;             //!< \c true if track detectors clusterization.
  std::string _MCClusteringMode;//!< Clustering mode on MC data.
  bool _doClustersAssociation;  //!< \c true if LR association on.
  bool _tracking;               //!< \c true if tracking.
  bool _reTrackingON;           //!< \c true when in the 2nd iteration of tracking and reconstruction, w/ T0 defined by 1st iteration.
  double _reTrackT0;            //!< T0 to be used in 2nd iteration of tracking and reconstruction. Is defined by the 1st iteration.
  bool _calorimeters;           //!< \c true if calo. reconstruction.
  bool _RWcalorimeters;         //!< \c true if calo. + RichWall combined reconstruction.
  bool _RWchargeEcal1;		//!< \c true if calo. + RichWall combined charge reconstruction.
  bool _rich1;                  //!< \c true if RICH1 reconstruction.
  bool _cedars;                 //!< \c true if cedar reconstruction.
  bool _beam;                   //!< \c true if BEAM reconstruction.
  int  _beamTR;                 //!< \c >0 if BEAM reco using TRAFFIC.
  int  _beamMethod;             //!< \c 0 - old reconstruction, 1 - new reconstruction 
  bool _vertex;                 //!< \c true if VERTEX reconstruction.
  bool _dst;                    //!< \c true if DST output.
  int  _dstversion;             //!< Version of DST.
  int  _dstfatness;             //!< Level of DST data content.
  std::map<CsDetector*,int,std::less<CsDetector*> > 
    _detmap;                    //!< Used by DST clusters.
  std::vector<CsDetector*>
    _detvec;                    //!< Used by DST clusters.
  std::ofstream _outStreams[32];//!< Array of output streams.
  int          _nOutStreams;    //!< Number of elem. in array of out streams.  
  unsigned int _outStreamsMask; //!< Streams where event has been written.
  unsigned int _selTrigMask;    //!< Trigger mask of events to select.
  bool         _selTrigStrict;  //!< \c true if strict mask, i.e. all bits of the event's trigger pattern must belong to the mask.
  bool         _selTrigMaskZero;  //!< Allow processing events with zero trigger mask. ( Calibration events of 2008 and earlier?)
  double       _TCSPhaseTime;   //!< TCS Phase Time.
  std::vector<double>
    _extraTimeWidths;           //!< Time widths to account for trigger jitter being larger than that of reference ILM triggers (in ns).
  double   _extraTimeWidth;     //!< Extra time with for current event.
  int      _alterMasterTrigger; //!< Bit# of the alternative master trigger.
  double   _triggerTimeCorr;    //!< Correction to trigger time.
  int _actualRICHDataSize;
  std::vector<int>
    _triggerPreScales;          //!< Prescaling factors.

  int  recoSchema_;              //!< reconstructione schema number
  CsTrkPrepattern* _prepattern1; //!< pointer to track prepattern package
  CsTrkPrepattern* _prepattern2; //!< pointer to 2nd track prepattern package
  CsTrkBridging*   _bridging;    //!< pointer to track bridging package
  CsTrkFitting*    _fitting;     //!< pointer to track fitting package

  CsVrtPattern*    _vpattern;    //!< pointer to track prepattern package
  CsVrtFitting*    _vfitting;    //!< pointer to track fitting package

  // A couple of general histograms :
  CsHist1D*        _hnmuperspill; //!< number of muons per spill

  CsHist1D*        _hntimeinspill; //!< event time in spill
  
  unsigned int     _nmuinspillsofar; //!< number of mu per spill so far...

  /// needed to resize _hnmuperspill. this is a first guess for the total
  /// number of spills.
  static int       _nspillmax;

	std::map<unsigned int,double> _tcscorr; //!spill bz spill correction for tcsphase jumps

};
#endif //CsEvent_h
