/*!
   \file    CsTrigHodDbReader.h
   \brief   Compass Some reader class definition file.
   \author  Vassili Motchalov
   \version $  $
   \date       $  $
*/

//-*-Mode: C++;-*-
#ifndef _CsTrigHodDbReader_h_
#define _CsTrigHodDbReader_h_

#if __GNUG__ >= 2
#  pragma interface
#endif

#include "CsCondDbReader.h"
//---- CsTrigHodDbReader -----------------------------------------------------------
/*!
  \class	CsTrigHodDbReader
  \brief	Some calibration constants DB Reader Object.
*/
class CsTrigHodDbReader : public CsCondDbReader {

public:
/*!
  \fn	        CsTrigHodDbReader(string  dbName = "trh");
  \brief	default constructor. Default name of DB is defined as "bms".
*/
  CsTrigHodDbReader(string  dbName = "trh");

  ~CsTrigHodDbReader();  


private:
  /*!
    \fn   bool checkType(const CsCalibAbstract& TempConst);
    \param        Container_name  Name of the Container to save    
    \brief        Check the type of  CsCalibAbstrac object
    \warning      Not ready yet - should be defined for each class
  */
  bool checkType(const CsCalibAbstract& TempConst) {return true;}

};

#endif









