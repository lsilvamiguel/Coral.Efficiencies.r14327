// $Id: CsZebraProto.h,v 1.1 2000/03/15 17:05:11 benigno Exp $
 
/*!
   \file    CsZebraProto.h
   \brief   Prototypes for Zebra functions.
   \author  Benigno Gobbo
   \version $Revision: 1.1 $
   \date    $Date: 2000/03/15 17:05:11 $
*/

#ifndef CsZebraProto_h
#define CsZebraProto_h

 
#define mzebra  mzebra_
#define mzstor  mzstor_
#define mzlink  mzlink_
#define mzwipe  mzwipe_
#define dzshow  dzshow_
#define fzin    fzin_
#define cfopen  cfopen_
#define fzfile  fzfile_
 
extern "C" void  mzebra( const int& );
extern "C" void  mzstor( const int&, const char*, const char*, int*, int*, 
			 int*, int*, int*, int*, const int, const int );
extern "C" void  mzlink( const int&, const char*, int*, int*, int*, 
			 const int );
extern "C" void  mzwipe( const int& );
extern "C" void  dzshow( const char*, const int&, const int&, const char*, 
			 const int&, const int&, const int&, const int&, 
			 const int, const int );
extern "C" void  fzin( const int&, const int&, int*, const int&, 
		       const char*, const int&, int*, const int );
extern "C" void  cfopen( const int&, const int&, const int&, const char*, 
			 const int&, const char*, const int&, const int, 
			 const int );
extern "C" void  fzfile( const int&, const int&, const char*, const int );

#endif // CsZebraProto_h

