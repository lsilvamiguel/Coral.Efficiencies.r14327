// $Id: CsMCRICH1Hit.h,v 1.4 2006/03/06 16:13:44 mutteran Exp $

/*!
   \file    CsMCRICH1Hit.h
   \brief   Compass Monte Carlo RICH1 Hit Class.
   \author  Guennadi Khaoustov
   \version $Revision: 1.4 $
   \date    $Date: 2006/03/06 16:13:44 $
*/

#ifndef CsMCRICH1Hit_h
#define CsMCRICH1Hit_h

#include <CLHEP/Vector/LorentzVector.h>
#include "CsDetector.h"
#include "CsMCHit.h"


/*! \class CsMCRICH1Hit 
    \brief Monte Carlo RICH1 Hit Class.
    Collect all needed Hit informations from Geant HITS (OHIT at present)
    bank. Objects of CsMCRICH1Hit class are istantiated by CsGeant3 methods.
*/

class CsMCRICH1Hit : public CsMCHit {
 public:
  CsMCRICH1Hit();  //!< Default Constructor

  /*! \fn   CsMCRICH1Hit( double xdet, double ydet, double zdet, 
    double Vdet, double Wdet, 
    double xprod, double yprod, double zprod,  
    double xref, double yref, double zref,
    double photenergy, double dtime, double cherangle,
    Hep3Vector p, CsMCTrack& MCTrack, 
    int orig, int cathode, CsDet& det );
    \brief Constructor.
    \param (xdet,ydet,zdet) Photon detection point in detector volume (MRS) (mm)
    \param (Vdet,Wdet) Photon detection point  in detector volume (COMGEANT DRS) (mm)
    \param (xprod,yprod,zprod) Photon production point (MRS) (mm)
    \param (xref,yref,zref) Photon last reflection point (MRS) (mm)
    \param photenergy Cherenkov photon  energy  (eV)
    \param dtime Delay time from the  Time Zero (ns)
    \param cherangle Cherenkov light emission angle
    \param p   Mother  Particle 3-momentum
    \param MCTrack Pointer to the associated Monte Carlo Track
    \param orig 0: hit from original track, N: from products
    \param cathode  Cathde number
    \param det Pointer to the associated CsDet object
  */

  CsMCRICH1Hit( double xdet,    double ydet,    double zdet, 
	        double Vdet,  double Wdet, 
	        double xprod, double yprod, double zprod, 
                double xref, double yref, double zref,
	        double photenergy, double dtime, double cherangle,
	        CLHEP::Hep3Vector p, CsMCTrack& MCTrack, 
                int orig, int cathode, CsDet& det ) :
                              CsMCHit( dtime, MCTrack, det ),
                              xdet_(xdet), ydet_(ydet), zdet_(zdet),
                              Vdet_(Vdet), Wdet_(Wdet),
                              xprod_(xprod), yprod_(yprod), zprod_(zprod),  
                              xref_(xref),  yref_(yref),  zref_(zref),
                              photenergy_(photenergy),
                              cherangle_(cherangle), cathode_(cathode),
                              orig_(orig), p_(p) { };

  CsMCRICH1Hit( const CsMCRICH1Hit& );              //!< Copy Constructor

  CsMCRICH1Hit& operator=( const CsMCRICH1Hit& );   //!< Assignment operator

  bool operator==( const CsMCRICH1Hit& ) const; //!< "equal to" operator

  bool operator<( const CsMCRICH1Hit& ) const;  //!< "less that" operator

  //! Returns the X component of the detection point (MRS) (mm)
  inline double getX() const { return(xdet_); }    

  //! Returns the Y component of the detection point (MRS) (mm)
  inline double getY() const { return(ydet_); }

  //! Returns the Z component of the detection point (MRS) (mm)
  inline double getZ() const { return(zdet_); }

  //! Returns the Y component of the detection point (DRS) (mm)
  inline double getYdetDRS() const { return(Vdet_); }

  //! Returns the Z component of the detection point (DRS) (mm)
  inline double getZdetDRS() const { return(Wdet_); }

  //! Returns the X component of the production point (MRS) (mm)
  inline double getXprod() const { return(xprod_); }    

  //! Returns the Y component of the production point (MRS) (mm)
  inline double getYprod() const { return(yprod_); }

  //! Returns the Z component of the production point (MRS) (mm)
  inline double getZprod() const { return(zprod_); }

  //! Returns the X component of reflection point (MRS) (mm)
  inline double getXref() const { return(xref_); }

  //! Returns the Y component of reflection point (MRS) (mm)
  inline double getYref() const { return(yref_); }

  //! Returns the Z component of reflection point (MRS) (mm)
  inline double getZref() const { return(zref_); }

  //! Returns the Z component of reflection point (MRS) (mm)
  inline double getPhotEnergy() const { return photenergy_; }

  //! Returns the cherenkov angle (rad)
  inline double getCherAngle() const { return(cherangle_); }

// changes by A. Mutter (28. Feb 2006) to account for the new COMGEANT 
// lens/cathode numbering scheme

  //! Returns the proper cathode number in the range 1 to 16
  inline int getCathode() const { 
	if(cathode_ < 17 ) return(cathode_); 
	if(cathode_ >  16 && cathode_ < 161) return( 4);
	if(cathode_ > 160 && cathode_ < 305) return( 6);
	if(cathode_ > 304 && cathode_ < 449) return(11);
	if(cathode_ > 448 && cathode_ < 593) return(13);
	return(0);
	}

  //! Returns the lens (or cathode, respectively) number 
  //  in the range 1-592 as set by COMGEANT
  inline int getLens() const { return(cathode_); }

// end of the changes by A. Mutter

  //! Returns the particle 3-momentum (GeV/c)
  inline CLHEP::Hep3Vector getP() const { return(p_); } 

  //! Returns 0 if the hit is from the original track, part id if from
  //! secondary products 
  inline int getOrigin() const { return(orig_);}


 private:
  double xdet_;             // photon detection point (MRS)
  double ydet_;
  double zdet_;
  double Vdet_;             // photon detection point (DRS) (Y, Z)            
  double Wdet_;
  double xprod_;            // photon production point (MRS)
  double yprod_;
  double zprod_;
  double xref_;             // last reflection point (MRS)
  double yref_; 
  double zref_;
  double photenergy_;       // photon energy
  double cherangle_;        //cherenkov light emission angle
  int    cathode_;          // cathde number
  int    orig_; // 0: original track, N: from products...
  CLHEP::Hep3Vector p_;
};

#endif // CsMCRICH1Hit_h
