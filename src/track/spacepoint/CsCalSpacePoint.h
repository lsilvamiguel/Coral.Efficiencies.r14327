// $Id: CsCalSpacePoint.h,v 1.2 2003/04/24 07:23:24 benigno Exp $

/*!
   \file    CsCalSpacePoint.h 
   \brief   Calibrarion spacepoint, derived from CsSpacePoint Class.
   \author  Hugo Pereira
   \version $Revision: 1.2 $
   \date    $Date: 2003/04/24 07:23:24 $
*/

#ifndef CsCalSpacePoint_h
#define CsCalSpacePoint_h

#include <list>
#include "CsSpacePoint.h"

class CsCluster;
class CsDetector;

class CsCalSpacePoint: public CsSpacePoint {
  
  public: 
	  CsCalSpacePoint( const CsDetFamily &df, std::list<CsCluster*> c, double z, int mode);  	//!< default constructor.				
	  CsCalSpacePoint( const CsCalSpacePoint& sp):CsSpacePoint( sp ) { *this = sp; }  //!< copy constructor.
	  CsCalSpacePoint& operator = ( const CsCalSpacePoint& );   	 //!< assignment operator.
    
    CsDetector* detOff_;
    bool found_;
    CsCluster* clFound_;
};

#endif
