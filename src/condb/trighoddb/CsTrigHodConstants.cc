/*!
   \file		CsTrigHodConstants.cc
   \brief		Compass ECAL1 calibration constant object implementation file. 
			This is a tempolary file.
   \author		Vassili Motchalov
   \version	        $Revision: 1.1 $
   \date		$Date: 2000/04/05 11:31:29 $	
*/
#include "CsTrigHodConstants.h"
#include "CsErrLog.h"
//
CsTrigHodConstants:: CsTrigHodConstants(){} 
//

CsTrigHodConstants::CsTrigHodConstants(CsTrigHodConstants& tempClass){
  double * myptr;
  double * ptr;
  myptr = &TrHodH4PlaneH   [ 0] [ 0];
  ptr = &tempClass.TrHodH4PlaneH   [ 0] [ 0];
  for( int i=0; i<size*dim; i++){*myptr++ = *ptr++ ;}
  myptr = &TrHodH4PlaneV   [ 0] [ 0];
  ptr = &tempClass.TrHodH4PlaneV   [ 0] [ 0];
  for( int i=0; i<size*dim; i++){*myptr++ = *ptr++ ;}
  myptr = &TrHodH4PlaneVL   [ 0] [ 0];
  ptr = &tempClass.TrHodH4PlaneVL   [ 0] [ 0];
  for( int i=0; i<size*dim; i++){*myptr++ = *ptr++ ;}
  myptr = &TrHodH5PlaneH   [ 0] [ 0];
  ptr = &tempClass.TrHodH4PlaneH   [ 0] [ 0];
  for( int i=0; i<size*dim; i++){*myptr++ = *ptr++ ;}
  myptr = &TrHodH5PlaneV   [ 0] [ 0];
  ptr = &tempClass.TrHodH5PlaneV   [ 0] [ 0];
  for( int i=0; i<size*dim; i++){*myptr++ = *ptr++ ;}
  myptr = &TrHodH5PlaneVL   [ 0] [ 0];
  ptr = &tempClass.TrHodH5PlaneVL   [ 0] [ 0];
  for( int i=0; i<size*dim; i++){*myptr++ = *ptr++ ;}
}

void CsTrigHodConstants::printCoeff(const string& Container_name,  ostream& os){
  os << "Constants   "  << Container_name << "  for Trig Hodoscopes  Calorimeter"  <<  endl;
  
  if (Container_name == "H4PlaneH"){
    for (int i = 0; i<size; i++){  
      for (int j = 0; j<dim; j++)  { os << "  " << TrHodH4PlaneH[i][j]; }
      os << endl;
    }
  }
  else if (Container_name == "H4PlaneV"){
    for (int i = 0; i<size; i++){  
      for (int j = 0; j<dim; j++)  { os << "  " << TrHodH4PlaneV[i][j]; }
      os << endl;
    }
  }
  else if (Container_name == "H4PlaneVL"){
    for (int i = 0; i<size; i++){  
      for (int j = 0; j<dim; j++)  { os << "  " << TrHodH4PlaneVL[i][j]; }
      os << endl;
    }
  }
  else if (Container_name == "H5PlaneH"){
    for (int i = 0; i<size; i++){  
      for (int j = 0; j<dim; j++)  { os << "  " << TrHodH5PlaneH[i][j]; }
      os << endl;
    }
  }
  else if (Container_name == "H5PlaneV"){
    for (int i = 0; i<size; i++){  
      for (int j = 0; j<dim; j++)  { os << "  " << TrHodH5PlaneV[i][j]; }
      os << endl;
    }
  }
  else if (Container_name == "H5PlaneVL"){
    for (int i = 0; i<size; i++){  
      for (int j = 0; j<dim; j++)  { os << "  " << TrHodH5PlaneVL[i][j]; }
      os << endl;
    }
  }
  else
    os << "  Container name not defined " << endl;
}

void CsTrigHodConstants:: printAllCoeff(ostream& os){
  printCoeff("H4PlaneH", os);  
  printCoeff("H4PlaneV", os);  
  printCoeff("H4PlaneVL", os);  
  printCoeff("H5PlaneH", os);  
  printCoeff("H5PlaneV", os);  
  printCoeff("H5PlaneVL", os);  
}

double* CsTrigHodConstants:: getCoeff(const string& Container_name) {
  if    (Container_name == "H4PlaneH" )
    return & TrHodH4PlaneH[0][0]; 
  else if    (Container_name == "H4PlaneV" )
    return & TrHodH4PlaneV[0][0]; 
  else if    (Container_name == "H4PlaneVL" )
    return & TrHodH4PlaneVL[0][0]; 
  else if    (Container_name == "H5PlaneH" )
    return & TrHodH5PlaneH[0][0]; 
  else if    (Container_name == "H5PlaneV" )
    return & TrHodH5PlaneV[0][0]; 
  else if    (Container_name == "H5PlaneVL" )
    return & TrHodH5PlaneVL[0][0]; 
  else {
    string mes;
    mes += " Container name is not defined in this Object ";
    mes += Container_name;
    CsErrLog::Instance()->mes( elInfo,  mes); 
    return 0;
  }
}

void  CsTrigHodConstants::clearCoeff(const string& Container_name) {
  double * ptr = getCoeff(Container_name) ;
  for( int i=0; i<size*dim; i++){*ptr++ = 0. ;}
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool CsTrigHodConstants:: setCoeff(const string& Container_name, const double *ptr) {
  int length = size*dim;
  return setCoeff(Container_name, ptr, length) ;
}
 
bool CsTrigHodConstants:: setCoeff(const string& Container_name, const double *ptr, const int length) {
  int my_length = dim*size;
  if (my_length<length) {
    string mes;
    mes += " Error: length of the Array is larger than defined";
    mes += " copy only first ";
    mes += my_length;
    mes += " elements";
    CsErrLog::Instance()->mes( elInfo,  mes); 
  }
  else if (my_length<length) {
    string mes;
    mes += " Error: length of the Array is less than defined";
    mes += " copy only first ";
    mes += length;
    mes += " elements from";
    mes += my_length;
    CsErrLog::Instance()->mes( elInfo,  mes); 
    my_length = length;
  }
   
  double *myptr = getCoeff(Container_name);
  if(!myptr) return false;
  for ( int i = 0;  i!=  size*dim;  i++) {*myptr++ = *ptr++;}
  return true;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool CsTrigHodConstants:: readCoeff (const string& name, const char* filename) {
  
  double *myptr = getCoeff(name);
  if (!myptr) return false;

  ifstream input_file(filename);
  if(input_file.fail()) {
    cerr << " Error during opening of input file "<< filename << endl;
    return false;}
  else {
    bool result; 
    float read;
    int i;
    //    for(i=0; (i<dim*size); i++){
    for(i=0; (i<dim*size) && (input_file >> *myptr++) ; i++) {}
    if (!input_file.good()){ 
      string mes;
      mes += " Error during reading";
      mes += filename;
      mes += " only ";
      mes += i;
      mes += " channels read ";
      CsErrLog::Instance()->mes( elInfo,  mes); 
      result = false;}
    else {
      result = true;
    }      
    input_file.close();  
    return result;
  }
}










