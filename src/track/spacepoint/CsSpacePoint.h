// $Id: CsSpacePoint.h,v 1.10 2003/04/24 07:23:26 benigno Exp $

/*!
   \file    CsSpacePoint.h 
   \brief   Space points bluid from clusters for a certain number of detectors.
   \author  Hugo Pereira
   \version $Revision: 1.10 $
   \date    $Date: 2003/04/24 07:23:26 $
*/

#ifndef CsSpacePoint_h
#define CsSpacePoint_h

#include "CLHEP/Matrix/Matrix.h"
#include "CsSTD.h"

class CsCluster;
class CsDetFamily;

class CsSpacePoint {

 public:
	CsSpacePoint( const CsDetFamily &df, std::list<CsCluster*> c, double z, int mode);  	//!< default constructor.				
	CsSpacePoint( const CsSpacePoint& sp) { *this = sp; }              //!< copy constructor.
	CsSpacePoint& operator = ( const CsSpacePoint& );	//!< assignment operator.
	virtual ~CsSpacePoint() {}										    //!< destructor.

  inline void fill_Fast( const double x, const double y, const double chi2 ) {
    minimised_Fast_ = true;
    x_Fast_ = x;
    y_Fast_ = y;
    chi2_Fast_ = chi2;
    return;
  }

  inline void fill( const double x, const double y, const double tx, const double ty, const double chi2 ) {
    minimised_ = true;
    x_ = x;
    y_ = y;
    tx_ = tx;
    ty_ = ty;
    chi2_ = chi2;
    return;
  }
    
  inline const CsDetFamily* getDetFamily( void ) const { return dF_; }
  inline std::list< CsCluster* > getClusters( void ) const { return c_; }
  inline unsigned int cSize( void ) const { return c_.size(); }

  inline bool minimised_Fast( void ) const { return minimised_Fast_; }
  inline bool getX_Fast( double &x ) const { x = x_Fast_ ; return minimised_Fast_; }  
	inline bool getY_Fast( double &y ) const { y = y_Fast_ ; return minimised_Fast_; }  
	inline bool getChi2_Fast( double &chi2 ) const { chi2 = chi2_Fast_ ; return minimised_Fast_; }  
  
  inline bool minimised( void )   const { return minimised_; }
	inline bool getX( double &x )   const { x = x_ ; return minimised_; }  
	inline bool getY( double &y )   const { y = y_ ; return minimised_; }  
	inline bool getZ( double &z )   const { z = z_ ; return true; }  
	inline bool getTx( double &tx ) const { tx = tx_ ; return minimised_; }  
	inline bool getTy( double &ty ) const { ty = ty_ ; return minimised_; }  
	inline bool getChi2( double &chi2 ) const { chi2 = chi2_ ; return minimised_; }  

  inline int getMode( void ) const { return mode_; }

  void dump( void );

 private: 
  
	const CsDetFamily* dF_; //!< detector family associated to space point
  
  std::list< CsCluster* > c_;  //!< list of associated clusters
  
	int mode_;              //!< straight tracks or target pointing, for fast minimisation.
	bool minimised_Fast_;   //!< was fast minimisation performed
	double x_Fast_;				  //!< x coordinate of fast minimisation Space Point
	double y_Fast_;				  //!< x coordinate of fast minimisation Space Point				
	double chi2_Fast_;		  //!< chi square for fast minimisation

  bool minimised_;        //!< was full minimisation performed
	double x_;					    //!< x coordinate of the full minimisation Space Point 		
	double y_;					    //!< y coordinate of the full minimisation Space Point 
  double z_;              //!< z coordinate of the space point (set by hand)
	double tx_;					    //!< tan(thetaX) in plane xOz (z along the beam)
	double ty_;					    //!< tan(thetaY) in plane yOz (z along the beam)
	double chi2_;				    //!< chi square for full minimisation 
  
};

#endif // CsSpacePoint_h
