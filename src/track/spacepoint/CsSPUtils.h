// $Id: CsSPUtils.h,v 1.11 2003/04/24 07:23:26 benigno Exp $

/*!
   \file    CsSPUtils.h
   \brief   Compass Space Points Utilities Class.
   \author  Hugo Pereira 
   \version $Revision: 1.11 $
   \date    $Date: 2003/04/24 07:23:26 $
*/


#ifndef CsSPUtils_h
#define CsSPUtils_h

#include "CsSTD.h"

/*! \class CsSPUtils 
    \brief Space Points Utilities Class.

    Contains some utilities to easy handle space points
*/

//class CsSpacePoint;
class CsDetFamily;
class CsSPUtils{
 public:

  static CsSPUtils* Instance( void );

  bool getOptForFamily( const CsDetFamily* pf, std::string tag, std::string key, double &par ); 
  bool getOptForFamily( const CsDetFamily* pf, std::string tag, std::string key, std::string &word ); 
  bool getOptForFamily( const CsDetFamily* pf, std::string tag, std::string key, std::list<std::string> &words);
  std::string getDate( void );
  std::string getTime( void ); 		
protected:

  CsSPUtils( void );  //!< The Constructor

private:

	static CsSPUtils* instance_;			        //!< the singleton instanciation
};

#endif //CsSPUtils_h
