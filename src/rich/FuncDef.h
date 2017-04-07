  #ifndef  FUNCDEF_H
  #define  FUNCDEF_H

/*!
   \file    FuncDef.h
   \-----------------
   \brief   Prototypes of common functions.
   \author  Paolo Schiavon
   \version 1.0
   \date    2 June 1999, rev. December 1999.
*/


  #include "CLHEP/Matrix/Matrix.h"

  void printVL( char*, int, int* );
  void printVL( char*, int, float* );
  void printVL( char*, int, double* );

//  void printVL( char*, int, HepVector );

  #endif
