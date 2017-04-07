// $Id: CsDetector.h,v 1.50 2010/10/18 16:40:15 schluter Exp $

/*!
   \file    CsDetector.h
   \brief   Compass Generic Tracking Detector Class.
   \author  Benigno Gobbo
   \version $Revision: 1.50 $
   \date    $Date: 2010/10/18 16:40:15 $
*/

#ifndef CsDetector_h
#define CsDetector_h

#include "coral_config.h"

#include <CLHEP/Matrix/Matrix.h>
#include <list>
#include <set>
#include "CsErrLog.h"
#include "CsDet.h"

class CsZone;
class CsMCHit;

namespace CS {class DriftDetector;}

/*! \class CsDetector 
    \brief   Compass Generic Tracking Detector Class.
*/

class CsDetector : public CsDet {

 public:

  /*! \fn CsDetector( const int row, 
    const int id, const char* name, const int unit,  
    const int type, const double rdLen, const double xsiz,  
    const double ysiz,  const double zsiz, const double xcm,   
    const double ycm, const double zcm, const HepMatrix rotDRS,
    const HepMatrix rotWRS, const double wirD, const double ang, 
    const int nWir, const double wirP, const double eff, const double bkg,
    const double tGate );
    \brief Constructor for tracking detector types.
    \param row   Detector file raw number
    \param id    Detector identification number
    \param name  Detector name (see Comgeant)
    \param unit  Detector number in station
    \param type  Detector type (see Comgeant)
    \param rdLen Radiation lenght
    \param xsiz  Detector X size (DRS) (mm)
    \param ysiz  Detector Y size (DRS) (mm)
    \param zsiz  Detector Z size (DRS) (mm)
    \param xcm   Detector centre X position (MRS) (mm)
    \param ycm   Detector centre Y position (MRS) (mm)
    \param zcm   Detector centre Z position (MRS) (mm)
    \param rotDRS Rotation matrix of DRS w.r.t MRS
    \param rotWRS Rotation matrix of WRS w.r.t MRS
    \param wirD  1st wire offset (WRS) (mm)
    \param ang   Rotation angle of WRS w.r.t MRS
    \param nWir  Number of wires
    \param wirP  Wires pitch
    \param eff   Detector efficiency
    \param bkg   Detector background
    \param tGate Detector gate
  */
  CsDetector( const int    row,
	      const int    id,    const char* name,  const char *TBname,
	      const int    unit,  const int    type,
	      const double rdLen, const double xsiz,  
	      const double ysiz,  const double zsiz,
	      const double xcm,   const double ycm,   
	      const double zcm,   const CLHEP::HepMatrix rotDRS,
	      const CLHEP::HepMatrix rotWRS,
	      const double wirD,  const double ang,   
	      const int    nWir,  const double wirP, 
	      const double eff,   const double bkg,
              const double tGate );

  //CsDetector( const CsDetector& ); //!< Copy constructor

  //CsDetector& operator=( const CsDetector& ); //!< Assign operator

  bool operator==( const CsDetector& ) const; //!< "equal to" operator

  bool operator<( const CsDetector& ) const; //!< "less than" operator
  
  /*! 
    \fn virtual void AddSubDetector( const int row, 
    const int id, const char* name, const int unit,  
    const int type, const double rdLen, const double xsiz,  
    const double ysiz,  const double zsiz, const double xcm,   
    const double ycm, const double zcm, const HepMatrix rotDRS,
    const HepMatrix rotWRS, const double wirD, const double ang, 
    const int nWir, const double wirP, const double eff, const double bkg,
    const double tGate )
    \brief Add some extra piece to an existing detector. Allows to construct a variable pitch detector in several successive steps.
    \warning This method can be re-defined for specific CsDetector's sub-classes.
    \param cxm Center of the extra piece of detector, X pos (MRS) mm
    \param ycm Center of the extra piece of detector, Y pos (MRS) mm
    \param zcm Center of the extra piece of detector, Z pos (MRS) mm
    \param nWir Number of wires in the extra piece.
    \param wirP Pitch (mm)    
  */
  virtual void AddSubDetector(const int    row,
			      const int    id,    const char* name,  
			      const char *TBname,
			      const int    unit,  const int    type,
			      const double rdLen, const double xsiz,  
			      const double ysiz,  const double zsiz,
			      const double xcm,   const double ycm,   
			      const double zcm,   const CLHEP::HepMatrix rotDRS,
			      const CLHEP::HepMatrix rotWRS,
			      const double wirD,  const double ang,   
			      const int    nWir,  const double wirP, 
			      const double eff,   const double bkg,
			      const double tGate );
  
  bool IsMyHit(int hitid);

  //! Returns detector number in station (see Comgeant)
  int       getUnit()  const { return(unit_); } 

  //! Returns detector type (see Comgeant)
  int       getType()  const { return(type_); } 

  //! Returns radiation lenght (gr/cm**2 ?)
  double    getRdLen() const { return(rdLen_); } 

  //! Returns detector X size (DRS) (mm)
  double    getXsiz()  const { return(xsiz_); }  

  //! Returns detector Y size (DRS) (mm)
  double    getYsiz()  const { return(ysiz_); }  

  //! Returns detector Z size (DRS) (mm)
  double    getZsiz()  const { return(zsiz_); }  

  //! Returns detector centre X position (MRS) (mm)
  double    getXcm()   const { return(xcm_); }   

  //! Returns detector centre Y position (MRS) (mm)
  double    getYcm()   const { return(ycm_); }   

  //! Returns detector centre Z position (MRS) (mm)
  virtual double    getZcm()   const { return(zcm_); }   

  //! Set shift on X for MC allignment tests.
  void setShiftOnX( double shift ) { _xshift = shift; }

  //! Set shift on Y for MC allignment tests.
  void setShiftOnY( double shift ) { _yshift = shift; }
 
  //! Get shift on X for MC allignment tests.
  double getShiftOnX( void ) { return( _xshift ); }

  //! Get shift on Y for MC allignment tests.
  double getShiftOnY( void ) { return( _yshift ); }

  // Set correction on Xcm position.
  void setDeltaXCorrection( double delta ) { _deltax = delta; }

  // Set correction on Ycm position.
  void setDeltaYCorrection( double delta ) { _deltay = delta; }

  // Get correction on Xcm position.
  double getDeltaXCorrection( void ) { return(_deltax); }

  // Get correction on Ycm position.
  double getDeltaYCorrection( void ) { return(_deltay); }

  //! Returns rotation matrix of DRS respect to MRS
  const CLHEP::HepMatrix& getRotDRS()  const { return(rotDRS_); }  

  const CLHEP::HepMatrix& getRotDRS2WRS()  const { return(rotDRS2WRS_); }    

  //! Returns rotation matrix transforming a vector in the WRS to the MRS
  const CLHEP::HepMatrix& getRotWRS()  const { return(rotWRS_); }  

  //! Returns rotation matrix transforming a vector in the MRS to the WRS
  const CLHEP::HepMatrix& getRotWRSInv()  const { return(rotWRS_inv_); }  

  //! Returns the first wire position (WRS) (mm)
  virtual double    getWirD()  const { return(wirD_); }  

  //! Returns the rotation angle of WRS respect to MRS
  double    getAng()   const { return(ang_); }   

  //! Returns the sine of the rotation angle of WRS respect to MRS
  double    getSinAng()const { return(sinAng_); }

  //! Returns the cosine of rotation angle of WRS respect to MRS
  double    getCosAng()const { return(cosAng_); }

  //! Returns number of wires
  int       getNWir()  const { return(nWir_); }  

  //! Returns the wires pitch (mm)
  double    getWirP()  const { return(wirP_); } 

  //! Returns the detector efficiency
  double    getEff()   const { return(eff_); }   

  //! Returns the detector background
  double    getBkg()   const { return(bkg_); }   

  //! Returns the detector gate
  double    getTGate() const { return(tGate_); }

  //! Returns list of pointers to the zones containing this detector
  std::list<CsZone*> getMyZones() { return( myZones_ ); }

  /*! \fn void addZone( CsZone& zone );
    \brief Add a zone this detector belongs to.
    \param zone the zone to be added
  */
  void addZone( CsZone& zone );

  /*! \fn void clearDigitsList();
    \brief Clear the list of detector digits.
  */
  void clearDigitsList() { decodingDone_ = false; myDigits_.clear(); }

  /*! \fn void addMCHit( CsMCHit& hit );
    \brief Add a Monte Carlo Hit to this detector.
    \param hit the hit to be added
  */
  void addMCHit( CsMCHit& hit );

  /*! \fn void clearMCHitList();
    \brief Clear the list of detector Monte Carlo hits.
  */
  void clearMCHitList();

  //! Returns list of pointers to the digits associated to this detector
  std::list<CsMCHit*> getMyMCHits() { return( myMCHits_ ); } 

  /*! \fn int getRow() const;
    \brief Get the detector.dat row number (Needed due to idiot hit-det 
    association in COMGEANT ntuples)
  */
  int getRow() const { return( row_ ); }

  //! \c true if this detector data was arleady decoded
  bool decoded() const { return( decodingDone_ ); }
  
  //! set this detector as decoded
  void setDecodingDone() { decodingDone_ = true; }

  //! set decoding on this detector
  void setDecode() { decode_ = true; }

  //! \c true if this detector data was arleady clusterized clusterize()
  bool clusterized() const { return( clusteringDone_ ); }
  
  /*! \fn std::list<CsCluster*> getMyClusters()
      \brief Returns list of pointers to the clusters associated to this 
             detector
  */
  std::list<CsCluster*> getMyClusters() { return( myClusters_ ); } 

  //!  Monitoring development, please dont use this function, contact Vladimir.Kolosov@cern.ch WARNING!!! another list myClusters_ exist also in base Det class !!!
  std::list<CsCluster*> &getMyClusters4MN() { return( myClusters_ ); } 

  /*! \fn void clearClusterList();
      \brief Clear the list of detector clusters.
  */
  void clearClusterList() { clusteringDone_ = false; myClusters_.clear(); }

  /*! \fn void addCluster( CsCluster& clus );
    \brief Add a cluster to this detector.
    \param clus the cluster to be added
  */
  void addCluster( CsCluster& clus ) { myClusters_.push_back( &clus ); }
  
  /*! \fn void sortClusters( void ) 
    \brief sort myClusters_ according to cluster->getU();
  */
  void sortClusters( void );

  //! set this detector as clusterized
  void setClusteringDone() { clusteringDone_ = true; }

  // Pure virtual methods -----------------------------

  //! Returns true if TCD informations are available
  virtual bool hasTDC() const = 0;

  //! Returns true if Drift informations are available
  virtual bool hasDrift() const = 0;

  //! Decode the MC data for this detector
  virtual void makeMCDecoding() = 0;

  //! Clusterize the digits
  virtual void clusterize() = 0;
  
  // Protection methods -------------------------------
  virtual double getTDCResol() const { 
    CsErrLog::Instance()->mes(elFatal,"Method not implemented on this object.");
    return(0.0); 
  }

  virtual double getVel() const { 
    CsErrLog::Instance()->mes(elFatal,"Method not implemented on this object.");
    return(0.0); 
  }  

  virtual double getT0() const { 
    CsErrLog::Instance()->mes(elFatal,"Method not implemented on this object.");
    return(0.0); 
  }    

  virtual double getThRes() const { 
    CsErrLog::Instance()->mes(elFatal,"Method not implemented on this object.");
    return(0.0); 
  } 

  virtual double getSpSli() const { 
    CsErrLog::Instance()->mes(elFatal,"Method not implemented on this object.");
    return(0.0); 
  } 

  virtual double getTiSli() const { 
    CsErrLog::Instance()->mes(elFatal,"Method not implemented on this object.");
    return(0.0); 
  } 

  virtual CsDetector* getAssociateDet() const {
    CsErrLog::Instance()->mes(elFatal,"Method not implemented on this object.");
    return 0; 
  } 

  virtual double getDistToWire(const double t, bool &error ) {
    CsErrLog::Instance()->mes(elFatal,"Method not implemented on this object.");
    return(0.0); 
  } 

  virtual double getCorrU(const CsCluster *c,
			  const double x, const double y, const double tt,
			  bool &error) {
    CsErrLog::Instance()->mes(elFatal,"Method not implemented on this object.");
    return(0.0); 
  } 

  /*! \fn void updateClusters(double time)
    \brief Updates the U of all clusters w/ time offset. Typically this offset corresponds to the time difference: event time (determined from tracking) minus trigger time (which has been used for the initialisation of U). The function leaves the "time_" attribute of the clusters unchanged.
    \param time  Offset in time (will be subtracted from cluster's time).
  */
  virtual void updateClusters(double time) {
    CsErrLog::Instance()->mes(elFatal,"Method not implemented on this object.");
  }

  /*! \fn double restoreU(CsCluster *left, CsCluster *right) const
    \brief Restore arg CsClusters' U (that might have been upset in the course
    of tracking) from their drift time and CsDigit data members. Left (i.e.
    Salève or bottom) is to be supplied first.
  */
  void restoreU(CsCluster *left, CsCluster *right);
  
  
    /*! \fn void rotatePointMRS2DRSOppanCOMGEANTStyle(const double x, const double y, const
        double z, double & u, double & v, double & w ) const
    \brief rotate a vector from MRS to DRS system. Takes into account subdetector size and does it the  "comgeant" way -
    this means that it takes into account the size of subdetector plane-sizes... we need this for comgeant-tgeant-compat reasons
  */
  void rotatePointMRS2DRSOppanCOMGEANTStyle(const double x, const double y, const double z, double& u, double& v, double& w, int hitID );

  
  

  /*! \fn void rotateVectDRS2MRS(const double u, const double v, const
        double w, double & x, double & y, double & z ) const
    \brief rotate a vector from DRS to MRS system 
    (rotation but no translation)
  */
  void rotateVectDRS2MRS(const double u, const double v, const
        double w, double & x, double & y, double & z ) const;

  /*! \fn void rotatePointDRS2MRS(const double u, const double v, const
        double w, double & x, double & y, double & z ) const
    \brief rotate and shift a point from DRS to MRS system.
    rotation AND translation.
  */
  void rotatePointDRS2MRS(const double u, const double v, const
        double w, double & x, double & y, double & z ) const;

  /*! \fn void rotateVectMRS2DRS(const double x, const double y, const
        double z, double & u, double & v, double & w ) const
    \brief rotate a vector from MRS to DRS system 
    (rotation but no translation)
  */
  void rotateVectMRS2DRS(const double x, const double y, const
        double z, double & u, double & v, double & w ) const;

  /*! \fn void rotatePointMRS2DRS(const double x, const double y, const
        double z, double & u, double & v, double & w ) const
    \brief rotate and shift a point from MRS to DRS system.
    rotation AND translation.
  */
  void rotatePointMRS2DRS(const double x, const double y, const
        double z, double & u, double & v, double & w ) const;

  /*! \fn void rotateVectWRS2MRS(const double u, const double v, const
        double w, double & x, double & y, double & z ) const
    \brief rotate a vector from WRS to MRS system 
    (rotation but no translation)
  */
  void rotateVectWRS2MRS(const double u, const double v, const
        double w, double & x, double & y, double & z ) const;

  /*! \fn void rotatePointWRS2MRS(const double u, const double v, const
        double w, double & x, double & y, double & z ) const
    \brief rotate and shift a point from WRS to MRS system.
    rotation AND translation.
  */
  void rotatePointWRS2MRS(const double u, const double v, const
        double w, double & x, double & y, double & z ) const;

  /*! \fn void rotateVectMRS2WRS(const double x, const double y, const
        double z, double & u, double & v, double & w ) const
    \brief rotate a vector from MRS to WRS system 
    (rotation but no translation)
  */
  void rotateVectMRS2WRS(const double x, const double y, const
        double z, double & u, double & v, double & w ) const;

  /*! \fn void rotatePointMRS2WRS(const double x, const double y, const
        double z, double & u, double & v, double & w ) const
    \brief rotate and shift a point from MRS to WRS system.
    rotation AND translation.
  */
  void rotatePointMRS2WRS(const double x, const double y, const
        double z, double & u, double & v, double & w ) const;

  /*! \fn void rotateVectWRS2DRS(const double x, const double y, const
        double z, double & u, double & v, double & w ) const
    \brief rotate a vector from WRS to DRS system 
    (rotation but no translation)
  */
  void rotateVectWRS2DRS(const double x, const double y, const
        double z, double & u, double & v, double & w ) const;

  /*! \fn void setDeadZone(const int shape, const double dim, cont double ydim,
        const double xcm, conts double ycm, const HepMatrix &D2DZ)
    \brief Set Dead Zone parameters.
  */
  void setDeadZone(const int type, const double dim, const double ydim,
		   const double xdrs, const double ydrs,
		   const CLHEP::HepMatrix &D2DZ) {
    dZtype_ = type; dZdim_ = dim; dZydim_ = ydim; dZxdrs_ = xdrs;
    dZydrs_ = ydrs; rotD2DZ_ = D2DZ; }

  //! Returns the Dead Zone type, i.e. shape + 16*ENABLE_FAST
  int       getDZType() const { return(dZtype_); }

  //! Returns the Dead Zone dimension = X 1/2size (mm) or radius squared (mm^2)
  double    getDZDim()  const { return(dZdim_); } 

  //! Returns the Dead Zone Y 1/2size (mm)
  double    getDZYdim() const { return(dZydim_); } 

  //! Returns the Dead Zone center X position (DRS) (mm)
  double    getDZXdrs()  const { return(dZxdrs_); } 

  //! Returns the Dead Zone center Y position (DRS) (mm)
  double    getDZYdrs()  const { return(dZydrs_); } 

  //! Returns rotation matrix of DZRS respect to DRS
  const CLHEP::HepMatrix& getRotD2DZ()  const { return(rotD2DZ_); }  

  /*! \fn inActiveArea(double x, double y) const
    \brief True if argument MRS(x,y) is in active area, i.e. w/in frame AND
    outside dead zone. (Nota bene: this routine is slow! But fairly general.)
  */
  bool inActiveArea(double x, double y) const;

  //! Returns true if argument DRS(x,y) is outside dead zone. (Nota bene:
  //  this routine is slow! But fairly general.)
  bool outDeadZone(double xdrs, double ydrs) const;

  //! book all histograms (they should not be created in the constructor)
  virtual void BookHistograms() { }

  /*! \fn double Wire2Pos(double wire) const
    \brief Converts distance to first wire from wire units to mm.
    (In case of a uniform pitch <P>, it returns: wire*P.)
  */
  virtual double Wire2Pos(double wire) const;

  /*! \fn int Pos2Wire(double pos)
    \brief Converts distance to first wire to wire number
    (In case of a uniform pitch <P>, it returns: nearest integer to pos/P.)
  */
  int Pos2Wire(double pos);

  //! \returns the pitch associated to a given wire
  double Pitch(double wire);

  //! \returns true if this is a variable pitch detector
  bool IsVarPitch() {return hasVarP_;}

  //! Fills the wirPs_ map and sets the variable pitch flag to true;
  void SetVarPitch( const std::map<int,double> *varP );

  enum histogramLevel { None, Normal, High };
  inline histogramLevel GetHistLevel() const { return hLevel_; }
  void ReadHistLevel( void );

 protected:
  
  histogramLevel hLevel_;  //!< None: no histos; Normal: std histos; High: all histos
   
  int    row_;                     //!< detector.dat row number
  int    unit_;                    //!<          unit
  int    type_;                    //!<          type
  double rdLen_;                   //!<          radiation length
  double xsiz_;                    //!<          X size
  double ysiz_;                    //!<          Y size
  double zsiz_;                    //!<          Z size
  double xcm_;                     //!<          X centre (MRS)
  double ycm_;                     //!<          Y centre (MRS)
  double zcm_;                     //!<          Z centre (MRS)
  
  // These value should be used on MC only
  double _xshift;                  //!< MC X shift, to test allignment    
  double _yshift;                  //!< MC Y shift, to test allignment
  
  double _deltax;                  //!< allignment correction on X
  double _deltay;                  //!< allignment correction on Y

  CLHEP::HepMatrix rotDRS_;        //!<          DRS rotation matrix
  CLHEP::HepMatrix rotWRS_;        //!<          WRS rotation matrix
  CLHEP::HepMatrix rotWRS_inv_;    //!<          WRS rotation matrix, inverted
  CLHEP::HepMatrix rotDRS2WRS_;    //!<          DRS to WRS rotation matrix

  double wirD_;                    //!<          wires distance
  double ang_;                     //!<          wires angle
  double sinAng_;                  //!<          sin( wires angle)
  double cosAng_;                  //!<          cos( wires angle)
  int    nWir_;                    //!<          number of wires
  double wirP_;                    //!<          wire pitch
  std::map<int,double> wirPs_;     //!< map<1st chan w/ pitch P, P> for variable pitch detectors   
  bool hasVarP_;                   //!< true if the detector has a variable pitch
  
  std::set<int> ids_;                   //!< ids of all subdetectors
  std::map<int, double> xcms_;    //!< xcm of subdetectors
  std::map<int, double> ycms_;    //!< ycm of subdetectors
  std::map<int, double> zcms_;    //!< zcm of subdetectors
  
  double eff_;                     //!<          efficiency
  double bkg_;                     //!<          barckground
  double tGate_;                   //!<          detector gate window

  int    dZtype_;                  //!< dead_zone.Shape+16*FAST_ENABLED
  double dZdim_;                   //!<           X 1/2size or radius squared
  double dZydim_;                  //!<           Y 1/2size if rectangular shape
  double dZxdrs_;                  //!<           X centre (DRS)
  double dZydrs_;                  //!<           Y centre (DRS)
  CLHEP::HepMatrix rotD2DZ_;       //!<           DRS -> DZRS rotation matrix

  std::list<CsZone*>  myZones_;         //!<         zones containing this detector
  std::list<CsMCHit*> myMCHits_;        //!<         MC hits associated to this detector
  std::list<CsCluster*> myClusters_;    //!<         Clusters associated to this detector

  bool           decodingDone_;    //!<         \c true if Detector decoded
  bool           clusteringDone_;  //!<         \c true if Detector clusterized
  bool           decode_;          //!<         Decoding set at run time

  //! used for sorting list of pointers to zones
  struct sortZones_;
  struct sortClusters_;
  struct sortMyDigits_;

  //! Sorts the digits
  virtual void sortMyDigits(void);  
};

#endif // CsDetector_h
