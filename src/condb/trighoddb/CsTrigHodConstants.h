/*!
   \file		CsTrigHodConstants.h
   \brief		Compass TrigHod  calibration constant object definition file. 
			This is a tempolary file until beam package is installed in coral.
   \author		Vassili Motchalov.
   \version	        $Revision: 1.2 $
   \date		$Date: 2000/05/25 16:33:34 $	
*/

#ifndef CsTrigHodConstants_h
#define CsTrigHodConstants_h

//
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "CsCalibAbstract.h"

class CsTrigHodConstants : public CsCalibAbstract {
public: 
  CsTrigHodConstants(void); 
  CsTrigHodConstants (CsTrigHodConstants&);              

  void printAllCoeff(ostream &);

  void printCoeff(const string& , ostream &);
  bool setCoeff(const string& , const double *); 
  bool setCoeff(const string& , const double *,const int length ); 
  bool readCoeff(const string& Container_name, const char*   filename);
  int  getLength(const string& Container_name)  {return size*dim;}    
  /*
    const double*  get_coeff(const string& Container_name) const;     
  */
  /*! get calibration coeefitients for ECAL1 */
  double*  getCoeff(const string& Container_name);   
  void clearCoeff(const string& Container_name);
private:
  static const int size=64;
  static const int dim=3;
  double TrHodH4PlaneH[size][dim];	        
  double TrHodH4PlaneV[size][dim];	        
  double TrHodH4PlaneVL[size][dim];	        
  double TrHodH5PlaneH[size][dim];	        
  double TrHodH5PlaneV[size][dim];	        
  double TrHodH5PlaneVL[size][dim];	        
};
//
#endif // TrigHodConstants




















