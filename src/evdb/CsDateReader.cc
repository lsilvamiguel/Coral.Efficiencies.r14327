// $Id: CsDateReader.cc,v 1.14 2010/01/20 23:09:30 tnagel Exp $

/*!
  \file    CsDateReader.cc
  \brief   Compass class to read DATE format files
  \author  Massimo Lamanna
  \version $Revision: 1.14 $
  \date    $Date: 2010/01/20 23:09:30 $

   \par History:
   0.1 10 May 1999 Massimo Lamanna<br>
   0.2 Dec 1999: cleanup the code (no new functionalities)<br>
   0.2.1 Remove int32->int uint8 -> unsigned char
   0.3 Mar 2001: remove rfio prototypes + usage of decoding library
*/  

#ifdef _WIN32
#include <winsock2.h>
 WSADATA wsadata;
#else
#include <unistd.h>
#endif

#include "CsDateReader.h"

#include "DaqDataDecoding/DaqEvent.h"

#include <fcntl.h>
#include <shift.h>

#include <stdio.h>

CsDateReader* CsDateReader::instance = 0;
int CsDateReader::decoderStatus(0);
unsigned char CsDateReader::array[] = {0};
int CsDateReader::nBytes(-1);
const int CsDateReader::mbuff = 3000000; //1MB buffer
int  CsDateReader::ifd(0);


CsDateReader* CsDateReader::Instance() {
    if(instance==0) {
		 instance = new CsDateReader;
    }
    return instance;
}

CsDateReader::CsDateReader()  {}

int mbuff4;

int CsDateReader::decoderAttach(char* fileName) {

#ifdef CS_DUMP

  bool dump(false);

  if(getenv("CSDATEREADER_DUMP")!=NULL) {
    dump = true;
  }

  dump && cout << "CsDateReader::decoderAttach file = "  << fileName << endl;

#endif // CS_DUMP
  char errString[100];

  mbuff4 = mbuff*sizeof(unsigned char);

  int istat = decoderStatus;

  int iret(-1);

#ifdef _WIN32
  //  WSAStartup( MAKEWORD(2,0), &wsadata); // old example

   int rcode;

   rcode = WSAStartup(MAKEWORD(2, 0),&wsadata);
   if( rcode )  {
        fprintf(stderr, "WSAStartup() failed, error code = %d\n", rcode);
        exit(1);
   }
#endif
  
  if(istat != 1) {

    int vNeg = RFIO_STREAM; // Negotiation flag: RFIO_STREAM for _V3
    rfiosetopt(RFIO_READOPT,&vNeg,4); 

    //    ifd = rfio_open(fileName, 0, 0700);
    ifd = open(fileName, 0, 0700);

  }
  else { 
    decoderStatus = 3;
    return (3);
  }

  if(ifd<0) {
    iret = sprintf(errString,
	   "CsDateReader: Error opening %s error code = %d",fileName,ifd);
    sperror(errString);
    decoderStatus = 3;
    return (3);
  }
    
  decoderStatus = 1;
  return (1);
}

int CsDateReader::decoderGetBuffer() {

  char errString[100];
  int iret,ir;

  int ibuf[32];

  static const int minDATEHeader(32);

  static bool toSwap(false);

  decoderStatus = 1;

  //  iret = rfio_read(ifd,(char *) array,minDATEHeader);
  iret = read(ifd,(char *) array,minDATEHeader);
  if(iret==0) {
    //    ir = sprintf(errString,"CsDateReader: End of file");
    //    sperror(errString);
    return(0);
  }
  else if(iret==-1) {
    ir = sprintf(errString,"CsDateReader: Empty file");
    sperror(errString);
    return(0);
  }
  else if(iret!=minDATEHeader) {
    ir = sprintf(errString,"CsDateReader: Bad format");
    sperror(errString);
    decoderStatus = 3;
    return(3);
  }

  unsigned int magic = *(unsigned int*) (&array[4]);
  if(magic != CS::DaqEvent::Header::EVENT_MAGIC_NUMBER) {
    if(!toSwap) {
      cout << " Magic number read from file " << hex 
	   << magic << dec << endl;
      cout << " EVENT_MAGIC_NUMBER " << hex 
	   << CS::DaqEvent::Header::EVENT_MAGIC_NUMBER  << dec << endl;
      cout << " Switch on swapping..." << endl;
      toSwap = true;
    }
  }
  
  int size = *(int*) (&array[0]);
  if(toSwap) {

    int tmp = ( ( size << 24 ) & 0xff000000 ) |
                ( ( size <<  8 ) & 0x00ff0000 ) |
                ( ( size >>  8 ) & 0x0000ff00 ) |
                ( ( size >> 24 ) & 0x000000ff );
    size = tmp;
  }


#ifdef CS_DUMP

  bool dump(false);

  if(getenv("CSDATEREADER_DUMP")!=NULL) {
    dump = true;
  }

  if(dump) {  
    cout << "CsDateReader::decoderGetBuffer size = "  << (unsigned int) array[0] << endl;
    cout << "CsDateReader::decoderGetBuffer magic = " << (unsigned int) array[1] << endl;
    cout << "CsDateReader::decoderGetBuffer type = "  << (unsigned int) array[2] << endl;
    cout << "CsDateReader::decoderGetBuffer head length = " << (unsigned int) array[3] << endl;
  }
#endif // CS_DUMP

    int eL4= size;
    if(eL4>mbuff4) {
        cout << "CsDateReader::decoderGetBuffer OVERFLOW" << endl;
	cout << "CsDateReader::decoderGetBuffer size = " << array[0] << endl;
	cout << "CsDateReader::decoderGetBuffer max bytes to be read = " << mbuff4 << endl;
	decoderStatus = 3;
        return (3);
    }
 
    int eL = eL4;
    eL4 = size - minDATEHeader;
    //    iret = rfio_read(ifd,((char *) array)+minDATEHeader,eL4);
    iret = read(ifd,((char *) array)+minDATEHeader,eL4);

#ifdef CS_DUMP

    dump && cout << "CsDateReader::decoderGetBuffer eL = "  << eL  << endl;
    dump && cout << "CsDateReader::decoderGetBuffer eL4 = " << eL4 << endl;
    dump && cout << "CsDateReader::decoderGetBuffer minDATEHeader = " << minDATEHeader << endl;
    dump && cout << "CsDateReader::decoderGetBuffer rfio_read rc code = " << iret << endl;

#endif // CS_DUMP

    if(iret != eL4 || iret==0) {
      cout << "Inconsistent buffer structure" << endl;
      cout << "Event length mismatch: eL= " << eL << endl;
      cout << "array[0] = " << array[0] << endl;
      ir = sprintf(errString,"CsDateReader: Event truncated (3)");
      sperror(errString);
      decoderStatus = 3;
      return(3);
    }

    nBytes = size;

#ifdef CS_DUMP

    dump && cout << "CsDateReader::decoderGetBuffer eL = " << (unsigned int) array[0] << endl;
    //    dump && cout << "CsDateReader::decoderGetBuffer array " << endl;
    
#endif // CS_DUMP

    return (1);
    
} 

//

int CsDateReader::decoderDetach() {

  char errString[100];
  //  int iret = rfio_close(ifd);
  int iret = close(ifd);

  if(iret!=0) {
    iret = sprintf(errString,
	   "CsDateReader: Error closing data file error code =%d",iret);
    sperror(errString);
    decoderStatus = 3;
    return(3);
  }

#ifdef _WIN32
 WSACleanup();
#endif

  decoderStatus = 0;
  return (0);
}

int CsDateReader::bufferLength() const { return (nBytes);}

unsigned char* CsDateReader::buffer() const { return (&array[0]);}

#undef CS_DUMP

