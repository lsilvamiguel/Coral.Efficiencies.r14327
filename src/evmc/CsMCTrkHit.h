// $Id: CsMCTrkHit.h,v 1.2 2001/11/29 18:27:55 bernet Exp $

/*!
   \file    CsMCTrkHit.h
   \brief   Compass Montecarlo Track Hit Class.
   \author  Benigno Gobbo
   \version $Revision: 1.2 $
   \date    $Date: 2001/11/29 18:27:55 $
*/


#ifndef CsMCTrkHit_h
#define CsMCTrkHit_h
#include <CLHEP/Vector/LorentzVector.h>
#include "CsMCHit.h"
#include "CsDetector.h"


/*! \class CsMCTrkHit 
    \brief Monte Carlo Track Hits Class.
    Collect all needed Hit informations from Geant HITS (OHIT at present)
    bank. Objects of CsMCTrkHit class are istantiated by CsGeant3 methods.
*/

class CsMCTrack;

class CsMCTrkHit :  public CsMCHit {

 public:

  CsMCTrkHit();  //!< Default Constructor

  /*! \fn   CsMCTrkHit( double x, double y, double z, double uin, double vin, 
    double win, double uout, double vout, double wout, double elos, 
    double eion, double dtime, Hep3Vector p, CsMCTrack& MCTrack, 
    int orig, CsDetector& det );
    \brief Constructor.
    \param (x,y,z) Centre of trajectory in detector volume (MRS) (mm)
    \param (uin,vin,win) Entrance point (DRS) (mm)
    \param (uout,vout,wout) Exit point (DRS) (mm)
    \param elos Total energy lost (GeV)
    \param eion Energy lost by ionization (GeV)
    \param dtime Delay time from Time Zero (ns)
    \param p     Particle 3-momentum
    \param MCTrack Pointer to the associated Monte Carlo Track
    \param orig 0: hit from original track, N: from products
    \param det Pointer to the associated CdDetector object
  */
  CsMCTrkHit( double x, double y, double z, 
	      double uin, double vin, double win, 
	      double uout, double vout, double wout, 
	      double elos, double eion, double dtime,
	      CLHEP::Hep3Vector p, CsMCTrack& MCTrack, int orig, 
	      CsDetector& det , int detid);

  CsMCTrkHit( double x, double y, double z, 
	      double uin, double vin, double win, 
	      double uout, double vout, double wout, 
	      double elos, double eion, double dtime,
	      CLHEP::Hep3Vector p, CsMCTrack& MCTrack, int orig, 
	      CsDetector& det);

  CsMCTrkHit( const CsMCTrkHit& );              //!< Copy Constructor

  CsMCTrkHit& operator=( const CsMCTrkHit& );   //!< Assignment operator

  bool operator==( const CsMCTrkHit& ) const; //!< "equal to" operator

  bool operator<( const CsMCTrkHit& ) const;  //!< "less that" operator

  //! Returns the X component of the trajectory in detector volume (MRS) (mm)
  inline double getX() const { return(x_); }    

  //! Returns the Y component of the trajectory in detector volume (MRS) (mm)
  inline double getY() const { return(y_); }

  //! Returns the Z component of the trajectory in detector volume (MRS) (mm)
  inline double getZ() const { return(z_); }

  //! Returns the X component of entrance point in detector volume (DRS) (mm)
  inline double getUin() const { return(uin_); }

  //! Returns the Y component of entrance point in detector volume (DRS) (mm)
  inline double getVin() const { return(vin_); }

  //! Returns the Z component of entrance point in detector volume (DRS) (mm)
  inline double getWin() const { return(win_); }

  //! Returns the X component of exit point from detector volume (DRS) (mm)
  inline double getUout() const { return(uout_); }

  //! Returns the Y component of exit point from detector volume (DRS) (mm)
  inline double getVout() const { return(vout_); }

  //! Returns the Z component of exit point from detector volume (DRS) (mm)
  inline double getWout() const { return(wout_); }

  //! Returns the Total Energy lost by the particle in detector volume (GeV)
  inline double getELos() const { return(elos_); }

  //! Returns the Ionization Energy lost by the particle in det. volume (GeV)
  inline double getEIon() const { return(eion_); }

  //! Returns the particle 3-momentum (GeV/c)
  inline CLHEP::Hep3Vector getP() const { return(p_); } 

  //! Returns 0 if the hit is from the original track, part id if from
  //! secondary products 
  inline int getOrigin() const { return(orig_);}

 private:
  double x_; // Centre of trajectory (MRS)
  double y_; 
  double z_;
  double uin_; // Entrance point (DRS) 
  double vin_;
  double win_;
  double uout_; // exit point (DRS)
  double vout_;
  double wout_;
  double elos_;
  double eion_;
  int    orig_; // 0: original track, N: from products...
  CLHEP::Hep3Vector p_;
};

#endif // CsMCTrkHit_h
