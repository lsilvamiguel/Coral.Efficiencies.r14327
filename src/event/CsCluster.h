// $Id: CsCluster.h,v 1.15 2007/01/06 09:28:52 ybedfer Exp $

/*!
   \file    CsCluster.h 
   \brief   Compass Single Coordinate Cluster Class.
   \author  Benigno Gobbo
   \version $Revision: 1.15 $
   \date    $Date: 2007/01/06 09:28:52 $
*/

#ifndef CsCluster_h
#define CsCluster_h


/*! \class CsCluster 
    \brief Compass Cluster Class.
*/

#include "CsSTD.h"
#include "CsDigit.h"
#include "CsDetector.h"
#include "CsZone.h"
#include <CLHEP/Matrix/Matrix.h>

#if USE_ObjectsCounter
#include "ObjectsCounter.h"
#endif

class CsCluster {

 public:

  CsCluster(); //!< Default Constructor.

  /*! \fn CsCluster( double u, double v, double w, HepMatrix cov );
    \brief Constructor.
    \param u first cluster coordinate
    \param v second cluster coordinate (NULL?)
    \param w third cluster coordinate
    \param cov covariant error matrix 
  */
  CsCluster( double u, double v, double w, const CLHEP::HepMatrix &cov );

  CsCluster( const std::vector<double> &xx, const CLHEP::HepMatrix &cov );

  CsCluster( const CsCluster& );            //!< Copy Constructor

  CsCluster& operator=( const CsCluster& ); //!< Assignment Operator

  //virtual ~CsCluster() {}                      //!< Destructor

  bool operator==( const CsCluster& ) const;   //!< "equal to" Operator

  bool operator<( const CsCluster& ) const;    //!< "less than" Operator

  /*! \fn void addDet( CsDetector& det );
    \brief add a CsDetector Object to the list of these objects owned
     by the Cluster. 
     
     It is used to create the list of detectors whose informations where
     user to create the cluster. The number of detectors could be greater
     than one.
     \param det The CsDetector object to be added.
  */
  void addDet( CsDetector& det );
  
  /*! \fn void addDigit( CsDigit& digit );
     \brief add a CsDigit Object to the list of these object owned by
     the Cluster. 

     It is used to create the list of digit from whom the cluster was
     made.
     \param digit The CsDigit object to be added.
  */
  void addDigit( CsDigit& digit );
  
  /*! \fn void setTime( const double time, const double timeError=1 );
    \brief set the value of time information 
    \param time the input value 
    \param timeError the error on input value (=1 if not specified)
   */
  void setTime( const double time, const double timeError=-1 );

  /*! \fn void setTimeError( const double timeError )
    \brief set the value of error on time information 
    \param timeError the error on time
   */
  void setTimeError( const double timeError );

  /*! \fn bool hasTime() const
    \brief returns \c true if the time value was set
   */
  inline bool hasTime() const { return( hasTime_ ); }

  /*! \fn bool getTime( double& time ) const
    \brief sets in \c time the time value if any. Returns \c true if
    the value is available \c false otherwise. The \c time variable is
    set to 0 if no time value is available.
    \param time the time returned value. \warning This variable is changed!
   */
  bool getTime( double& time ) const;

  /*! \fn bool getTimeError( double& timeError ) const
    \brief sets in \c timeError the error on time value if any. 
    Returns \c true if the value is available \c false otherwise. 
    The \c timeError variable is set to 1 if no time value is available.
    \param timeError the error on time returned value. 
    \warning This variable is changed!
   */
  bool getTimeError( double& timeError ) const;
  
  /*! \fn void setAnalog( const double analog, const double analogError=1 )
    \brief set the value of analog information 
    \param analog the input value 
    \param analogError the error on input value (set to 1 if not specified)
   */
  void setAnalog( const double analog, const double analogError=-1 );
  
  /*! \fn void setAnalogError( const double analogError )
    \brief set the error on analog information 
    \param analogError the error on input value (set to 1 if not specified)
   */
  void setAnalogError( const double analogError );
    
  /*! \fn bool hasAnalog() const
    \brief returns \c true if the analog value was set
   */
  inline bool hasAnalog() const { return( hasAnalog_ ); }

  /*! \fn bool getAnalog( double& analog ) const
    \brief sets in \c analog the analog value if any. Returns \c true if
    the value is available \c false otherwise. The \c analog variable is
    set to 0 if no analog value is available.
    \param analog the anlog returned value. \warning This variable is changed!
   */
  bool getAnalog( double& analog ) const;

  /*! \fn bool getAnalogError( double& analogError ) const
    \brief sets in \c analogError the error on analog value if any. 
    Returns \c true if the value is available \c false otherwise. 
    The \c analogError variable is set to 1 if no analog error value is 
    available.
    \param analog the anlog returned value. \warning This variable is changed!
   */
  bool getAnalogError( double& analogError ) const;

  //! adds more analog data, if needed
  void addAnalogData( const double analogData, const double analogError=-1 );

  //! update the nth vector element on analog data. Return false if that element does not exist.
  bool updateAnalogDataVectElement( const double data, const unsigned int n );

  //! update the nth vector element on analog data errors. Return false if that element does not exist.
  bool updateAnalogDataErrorsVectElement( const double data, const unsigned int n );

  //! get the size of analog data vector
  unsigned int getAnalogDataSize( void ) { return analog_.size(); }

  //! get the complete analog data vector, could be empty.
  std::vector<double> getAllAnalogData( void ) { return analog_; }

  //! get the complete analog data error vector, could be empty.
  std::vector<double> getAllAnalogDataErrors( void ) { return analogErr_; }

  /// \return clusters coordinates
  const std::vector<double>& getCoords(void) const {return x;}

  /// \return clusters coordinates
  std::vector<double>& getCoords(void)       {return x;}
	
  /*! \fn inline double getU() const;
     \brief Returns the first Cluster Coordinate.
     
     For tracking devices, this is generally the coordinate orthogonal to 
     the wires direction, with origin the one of the Compass Main Reference 
     System.
  */
  double getU() const { 
    if( x.size()<1 ) throw "CsCluster::getU(): size<1";  return x[0]; 
  }

  /*! \fn inline double setU();
     \brief Sets the first Cluster Coordinate.
  */
  void setU(const double u) { 
    if( x.size()<1 ) throw "CsCluster::setU(): size<1"; x[0] = u; 
  }

  /*! \fn inline double getV() const;
     \brief Returns the second Cluster Coordinate. (NULL?)

     For tracking devices, this is generally the coordinate parallel to 
     the wires direction, with origin the one of the Compass Main Reference 
     System.
  */
  double getV() const { 
    if( x.size()<2 ) throw "CsCluster::getV(): size<2";  return x[1]; 
  }

  /*! \fn inline double getW() const;
     \brief Returns the third Cluster Coordinate.

     Generally this is the Z coordinate, along the beam direction, with origin
     coincident with the Compass Main Reference System.
  */
  inline double getW() const { 
    if( x.size()<3 ) throw "CsCluster::getU(): size<3";  return x[2]; 
  }

  /*! \fn inline HepMatrix getCov() const;
     \brief Returns the Covariant Error Matrix.
  */
  inline CLHEP::HepMatrix getCov() const { return( cov_ ); }

  /*! \fn inline double setSigmaU();
     \brief Sets the uncertainty on first Cluster Coordinate.
  */
  void setSigmaU(const double sigmaU) { cov_(1,1) = sigmaU*sigmaU; }

  /*! \fn std::list<CsDetector*> getDetsList() const;
     \brief Returns the list of CsDetector objects associated to this cluster.
  */
  const std::list<CsDetector*> &getDetsList() const { return dets_; }

  /*! \fn std::list<CsDigit*> getDigitsList() const;
     \brief Returns the list of CsDigit objects associated to this cluster.
  */
  const std::list<CsDigit*>   &getDigitsList() const { return( digits_ ); }

  /* \fn std::list<CsZone*> getZonesList();
     \brief Returns the list of CsZones associated to this cluster.
  */
  std::list<CsZone*> getZonesList();
  
  /*! \fn inline bool hasAssociateClusters() const
    \brief true if anny associate cluster present
  */
  inline bool hasAssociateClusters() const { return( hasAssociates_ ); }
	
  /*! \fn inline void addAssociateCluster(CsCluster& clus)
    \brief add a cluster to the list of associated clusters
    \param clus the cluster to be add
  */	
  inline void addAssociateCluster( CsCluster& clus ) { hasAssociates_=true; 
  associates_.push_back( &clus ); }

  /*! \fn	inline void removeAssociateClusters( );
    \brief empty the associates_ cluster list
    and set hasAccosiates_ to false
  */
  inline void removeAssociateClusters() { hasAssociates_=false; 
  associates_.clear(); }

  /*! \fn const std::list<CsCluster*> getAssociateClusters() const
    \brief return the list of pointers to the associate clusters
  */
  const std::list<CsCluster*> &getAssociateClusters() const {
    return( associates_ ); }

  /*! \fn inline bool hasMirrorCluster() const
    \brief true if any mirror cluster set
  */
  inline bool hasMirrorCluster() const { return( hasMirror_ ); } 

  //! Set the mirror of this cluster
  inline void setMirrorCluster( CsCluster& clus ) { hasMirror_ = true; 
  myMirror_ = &clus; }

  //! Get mirror cluster in any. Points to NULL if no mirror available.
  inline CsCluster* getMirrorCluster() const { return( myMirror_ ); }

  /*! \fn inline void setLRProb( const double prob)
    \brief set the probability for the cluster to be true or mirror
    \param prob the Left/Right probability
  */
  inline void setLRProb( const double prob ) { LRProb_ = prob; }
  
  /*! \fn inline double getLRProb() const
    \brief return the probability for the cluster to be true or mirror
  */
  inline double getLRProb() const { return( LRProb_ ); }

  /*! \fn inline bool isGenuine() const
    \brief Is the genuine one in a pair of mirrors from a drift detector, as determined from MC truth bank.
    \warning If not MC or not a drift detector, it's undetermined.
   */
  inline bool isGenuine() const { return isGenuine_; }

  /*! \fn inline bool isGenuine(bool truthValue)
    \brief Set whether this->CsCluster is the genuine one in a pair of mirrors from a drift detector.
    \warning To be used only for the sole drift detector case, and in MC.
   */
  inline void isGenuine(bool truthValue) { isGenuine_ = truthValue; }

 private:
  std::vector<double> x;          //!< cluster's coordinates.
  CLHEP::HepMatrix cov_;          //!< The Covariant Error Matrix

  std::list<CsDetector*> dets_;   //!< The list of associated CsDetector Objs
  std::list<CsDigit*>    digits_; //!< The list of associated CsDigit Objects

  bool hasAnalog_;                //!< True if any analog information
  std::vector<double> analog_;    //!< Analog value (if any)
  std::vector<double> analogErr_; //!< Error on Analog value (if any)

  bool hasTime_;              //!< True if any time information
  double time_;               //!< Time information (if any) 
  double timeErr_;            //!< Error on Time value (if any)

  bool hasAssociates_;	      //!< true if any associated cluster (on DC shifted planes)
  std::list<CsCluster*> associates_; //!< The list of associated clusters (on DC shifted planes)

  bool hasMirror_;           //!< true if any mirror colleague exist.
  CsCluster* myMirror_;      //!< Points to its mirror colleague, if any
	
  double LRProb_;	     //!< Probability for the cluster to be true or mirror. set to 1 by default.

  bool isGenuine_; //!< Is the genuine one in a pair of mirrors from a drift detector, as determined from MC truth bank.

  #if USE_ObjectsCounter
  ObjectsCounter<CsCluster> objects_counter;
  #endif
};

#endif // CsCluster_h
