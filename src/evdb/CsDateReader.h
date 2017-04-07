#ifndef CSDATEREADER_H
#define CSDATEREADER_H

/*!
   \file CsDateReader.h
   \brief DATE Reader class
   \version 0.2.1
   \author  Massimo Lamanna
   \date    $Date: 2010/01/20 23:09:30 $

   \par History:
   0.1 10 May 1999 Massimo Lamanna<br>
   0.2 Dec 1999: cleanup the code (no new functionalities)<br>
   0.2.1 Remove int32->int uint8 -> unsigned char
*/

#include "CsPlatforms.h"
#include "CsTypes.h"

#ifdef COMPASS_USE_OSPACE_STD
#  include <ospace/std/iostream>
#  include <ospace/std/fstream>
#else
#  include <iostream>
#  include <fstream>
#endif
//#include "StageHelper.h"

/*!
    \class CsDateReader
    \brief Engine to access event data from Objectivity/DB
           (Singleton pattern)

    The CsDateReader class should be able to read
    date buffers in both endian flavours.
    The method is for test only and not optimised.
    Open point for optimisations are:

    - Byte swapping (now check & swap independently each event

    - Create a swap and a no_swap versions to avoid ifs
*/

class CsDateReader {
public:
    static CsDateReader* Instance();
//
    int decoderAttach(char*);
    int decoderGetBuffer();
    int decoderDetach();
//
    int bufferLength() const;
    unsigned char* buffer() const;
//
    int getDecoderStatus() {return(decoderStatus);}    
//
private:
    CsDateReader();

    static CsDateReader* instance;
    static int decoderStatus;  //0=off;1=reading;2=eof; 3=err
    static unsigned char array[500000];  //100kB buffer
    static int nBytes;
    static const int mbuff;
    static int ifd;            //File descriptor
};

#endif //CSDATEREADER_H
