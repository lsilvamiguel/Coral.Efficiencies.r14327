// $Id: CsSPMaker.h,v 1.7 2003/04/24 07:23:26 benigno Exp $

/*!
   \file    CsSPMaker.h
   \brief   Compass event/event Space Point Maker
   \author  Hugo Pereira
   \version $Revision: 1.7 $
   \date    $Date: 2003/04/24 07:23:26 $
*/

#ifndef CsSPMaker_h
#define CsSPMaker_h

#include "CsSTD.h"
#include <CLHEP/Matrix/Matrix.h>

class CsDetFamily;
class CsCluster;
class CsSpacePoint;
//_____________________________________________________________________________
class CsSPMaker {
public:

  static CsSPMaker* Instance( void );
	inline std::vector<CsDetFamily*> getDetFamilies( void ) const { return df_; }
  void cleanEvent( void );
protected:
	
  CsSPMaker( void );  //!< The Constructor

private:

	bool readDetFamilies( void );	
	static CsSPMaker* instance_;			  //!< the singleton instanciation
	std::vector<CsDetFamily*> df_;  		//!< vector of detector families used to build spacepoints
  std::vector<CsSpacePoint*> sp_;          //!< vector of build spacepoint
	
};

#endif	
