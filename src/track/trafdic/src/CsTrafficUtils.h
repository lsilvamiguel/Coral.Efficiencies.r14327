// $Id: CsTrafficUtils.h 13148 2011-12-28 16:55:25Z kbicker $

/*!
  \file    CsTrafficUtils.h
  \brief   Coral interface to Traffic Tracking Utils == TLatticeGen (so far)
  \author  Yann.Bedfer@cern.ch
  \version $Revision: 13148 $
  \date    $Date: 2011-12-28 17:55:25 +0100 (Wed, 28 Dec 2011) $ 
*/

#ifndef CsTrafficUtils_h
#define CsTrafficUtils_h

#include "CsSTD.h"
#include "CsTrkUtils.h"

/*! \class CsTrafficUtils 
    \brief Coral interface to Traffic Tracking Utils
*/

class CsTrafficUtils : public CsTrkUtils {

 public:

  /*! \fn CsTrafficUtils()
    \brief ...
  */
  CsTrafficUtils();

  /*! \fn ~CsTrafficUtils()
    \brief ...
  */
  ~CsTrafficUtils();

  /*! \fn bool genLattice()
    \brief ... 
  */
  bool genLattice();

 private:



};

#endif //CsTrafficUtils_h








