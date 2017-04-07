// $Id: CsTypes.h,v 1.3 2000/06/06 09:02:55 zvyagin Exp $

/*!
   \file    CsTypes.h
   \brief   Some Private Types Definitions.
   \author  Benigno Gobbo
   \version $Revision: 1.3 $
   \date    $Date: 2000/06/06 09:02:55 $
*/

#ifndef CsTypes_h
#define CsTypes_h

#include "coral_config.h"

#ifndef __CINT__

#if defined(__DECCXX) || defined(__alpha__)	
//! Digital Unix (C++ or g++) definitions
typedef signed char           int8;
typedef signed short int      int16;
typedef signed int            int32;
//typedef signed long           int64;
//typedef signed int            long32; 
//typedef signed long           long64;
typedef unsigned char         uint8;
typedef unsigned short int    uint16;
typedef unsigned int          uint32;
typedef unsigned long         uint64;
typedef unsigned int          ulong32;
typedef unsigned long         ulong64;
typedef float                 float32;
typedef double                float64;
#else

#if defined(_AIX)	
//! AIX definitions
//typedef signed int            long32;
//typedef signed long long      long64;
typedef unsigned char         uint8;
typedef unsigned short int    uint16;
typedef unsigned int          uint32;
typedef unsigned long long    uint64;
typedef unsigned int          ulong32;
typedef unsigned long long    ulong64;
typedef float                 float32;
typedef double                float64;
#else

//! All other architectures (but DEC and AIX) definitions
//typedef signed long           long32;
//typedef signed long long      long64;
typedef unsigned long         ulong32;
typedef unsigned long long    ulong64;
#if !defined(OO_CONFIG_H)
typedef signed char           int8;
typedef signed short          int16;
typedef signed int            int32;
//typedef signed long           int64;
typedef unsigned char         uint8;
typedef unsigned short        uint16;
typedef unsigned int          uint32;
typedef unsigned long long    uint64;
typedef float                 float32;
typedef double                float64;
#endif

#endif // AIX
#endif // Digital Unix

#else // __CINT__

typedef   Char_t   int8;
typedef  UChar_t  uint8;

typedef  Short_t   int16;
typedef UShort_t  uint16;

typedef  Int_t  long32,  int32;
typedef UInt_t ulong32, uint32;

typedef Long_t  long64,  int64;
typedef Long_t ulong64, uint64;

typedef Float_t  float32;
typedef Double_t float64;

#endif // __CINT__

#endif // CsTypes_h 
