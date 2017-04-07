/*!
   \file CsMagDbUpdater.h
   \brief COMPASS Magnetic Filed Database Class
   \author  Takeaki TOEDA
   \version $Revision: 1.2 $
   \date    $Date: 2000/05/31 12:28:26 $
*/

#ifndef CsMagDbUpdater_h
#define CsMagDbUpdater_h

#include "CsCondDbUpdater.h"
#include "CsSTD.h"

/*! \class CsMagDbUpdater
 *  \brief Updater for Magnetic Field database
 */

class CsMagDbUpdater : public CsCondDbUpdater{
public:
  CsMagDbUpdater(const string& dbName = "mag");     //!< Default Constructor
  ~CsMagDbUpdater();                                //!< Destructor

  /* \fn bool storeFieldCalcSol(string path,string containerName);
    \brief Store solnenoid field data from text file to database
  */
  bool storeFieldCalcSol(string path,string containerName);

  /* \fn storeFieldCalcSM1(string path,string containerName)
     \brief Store SM1 field data from text file to database
   */
  bool storeFieldCalcSM1(string path,string containerName);

  /* \fn bool storeFieldCalcSM2(string path,string containerName)
     \brief Store SM2 field data from text file to database
   */
  bool storeFieldSM2(string path,string containerName);

  /* \fn bool update(const string& Container_name,  const CsTime& startTime,const CsTime& endTime);
     \brief This function is not used in the class.
   */
//bool update(const string& Container_name,  const CsTime& startTime,
//const CsTime& endTime);

  /* \fn bool checkType(const CsCalibAbstract& TempConst);
     \brief This function is not used in the class.
   */
//  bool checkType(const CsCalibAbstract& TempConst);
private:

  /* \fn int estimateGrid(string path);
     \brief Estimate number of grid
   */
  int estimateGrid(string path); 
  
};
#endif
