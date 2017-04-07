/*!
   \file CsMagDbReader.h
   \brief COMPASS Magnetic Filed Database Class
   \author  Takeaki TOEDA
   \version $Revision: 1.4 $
   \date    $Date: 2000/06/26 14:08:19 $
*/

#ifndef CsMagDbReader_h
#define CsMagDbReader_h

#include "CsSTD.h"
#include "CsCondDbReader.h"
#include "CsMagFieldSol.h"
#include "CsMagFieldSM1.h"
#include "CsMagFieldSM2.h"

/*! \class CsMagDbReader
 *  \brief Magnetic Field Class for Solenoid.
 */

class CsMagDbReader : public CsCondDbReader{
public:
  CsMagDbReader(const string& dbName = "mag");//!< Default Constructor
 ~CsMagDbReader();                            //!< Destructor

  /*! \fn CsMagFieldSol getFieldSol(string type);
    \brief 
  */
  CsMagFieldSol getFieldSol(string type);

  /*! \fn CsMagFieldSM1 getFieldSM1(string type,int poleDistance);
    \brief 
    \parm type 
  */
  CsMagFieldSM1 getFieldSM1(string type,int poleDistance);

  /*! \fn CsMagFieldSM2 getFieldSM2(int current);
    \brief 
    \parm type 
  */
  CsMagFieldSM2 getFieldSM2(int current);

  /*! \fn bool checkType(const class CsCalibAbstract &);
    \brief 
  */
  //bool checkType(const class CsCalibAbstract &);
};
#endif

