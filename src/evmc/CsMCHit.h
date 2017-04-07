// $Id: CsMCHit.h,v 1.11 2009/09/14 00:42:21 ybedfer Exp $

/*!
   \file    CsMCHit.h
   \brief   Compass Montecarlo Hits Class.
   \author  Guennadi Khaoustov
   \version $Revision: 1.11 $
   \date    $Date: 2009/09/14 00:42:21 $
*/

#ifndef CsMCHit_h
#define CsMCHit_h

#include <CLHEP/Vector/LorentzVector.h>
#include "CsDet.h"
#include "CsErrLog.h"

class CsMCTrack;

/*! \class CsMCHit 
    \brief Monte Carlo Hits base Class.
    Collect base Hit informations from Geant HITS (OHIT at present)
    bank. 
*/

class CsMCHit{
 public:
  CsMCHit();  //!< Default Constructor

  virtual ~CsMCHit() {}              //!< Virtual destructor

  CsMCHit( double dtime, CsMCTrack& MCTrack, CsDet& det, int detid);  

  CsMCHit( double dtime, CsMCTrack& MCTrack, CsDet& det );

  CsMCHit( const CsMCHit& );              //!< Copy Constructor

  CsMCHit& operator=( const CsMCHit& );   //!< Assignment operator

  bool operator==( const CsMCHit& ) const; //!< "equal to" operator

  bool operator<( const CsMCHit& ) const;  //!< "less that" operator

  //! Returns the delay time from Time Zero point (ns)
  inline double getDTime() const { return(dtime_); }

  //! Returns the pointer to the associated Monte Carlo track
  inline CsMCTrack* getMCTrack() const { return(MCTrack_); }
  
  //! Returns the pointer to the associated Detector
  inline        CsDet* getDet()       { return(det_); }
  inline  const CsDet* getDet() const { return(det_); }
  inline  int          GetID()  const { return id_;}

  //! Returns the particle 3-momentum (GeV/c)
  virtual CLHEP::Hep3Vector getP() const { 
    CsErrLog::mes( elFatal, "Wrong method call" ); 
    CLHEP::Hep3Vector p; return( p ); 
  }  

  //! Returns 0 if the hit is from the original track, part id if from
  //! secondary products 
  virtual int getOrigin() const { 
    CsErrLog::mes( elFatal, "Wrong method call" ); 
    return( 0 ); 
  }

  //! Returns hit X coordinate (MRS)
  virtual double getX() const { 
    CsErrLog::mes( elFatal, "Wrong method call" ); 
    return( 0.0 ); 
}

  //! Returns hit Y coordinate (MRS)
  virtual double getY() const { 
    CsErrLog::mes( elFatal, "Wrong method call" ); 
    return( 0.0 ); 
  }

  //! Returns hit Z coordinate (MRS)
  virtual double getZ() const { 
    CsErrLog::mes( elFatal, "Wrong method call" ); 
    return( 0.0 ); 
}

 private:
  double dtime_;
  CsMCTrack* MCTrack_;
  CsDet* det_;
  int    id_;
};

#endif // CsMCHit_h
