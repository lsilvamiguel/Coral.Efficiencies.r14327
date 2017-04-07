// Author: Massimo Lamanna (Massimo.Lamanna@cern.ch) 1998-99
// $Id: eventout.cc,v 1.3 2001/06/08 11:46:08 objsrvvy Exp $ 

#include <stdlib.h>

#include <fcntl.h>


#include <unistd.h>

#include <ctype.h>

#include <time.h>
#include <sys/timeb.h>

#include "DaqEvent.h"

// Warning might clash with <iostream> and <string>
#include <shift.h>

extern "C" {
  int ftime(struct timeb *); 
}

/*!
   \file eventout.cc
   \brief CCF test program to produce fake events
   \version $Revision: 1.3 $
   \author  Massimo Lamanna
   \date    $Date: 2001/06/08 11:46:08 $
*/

float myPoisson(int);

int main(int argc, char *argv[]) {

  struct timeb eventTime;

  if(argc!=4) {
    printf("usage: eventout.exe nEvents mSize outFile\n"); 
    printf("mSize int number of kB\n"); 
    exit(1);
  }
  
  int nEvents = atoi(argv[1]);
  int mSize = atoi(argv[2]);
  char errString[100];
  char outFile[100];
  
  int iret = sprintf(outFile,"%s",argv[3]);
    
  int vL = mSize*1000; // 4 times mSize bytes...
  
  int* ibuf = new int[vL];
  if(!ibuf) {
    printf("Heap full (ibuf)\n");
    exit(1);
  }

//
// to verify _V3 vs _V2: setenv RFIO_TRACE 3
//
  int vNeg = RFIO_STREAM; // Negotiation flag: RFIO_STREAM for _V3
  rfiosetopt(RFIO_READOPT,&vNeg,4); 
  int fh = rfio_open(outFile,O_WRONLY | O_CREAT | O_TRUNC, 0770);

  //    RFILE* fh;
  //fh = rfio_fopen(outFile,"w");

  if( fh == -1 ) {
    //     if( fh == NULL ) {
    iret = sprintf(errString,"Error opening %s (fh=%d)",outFile,fh);
    sperror(errString);
    exit(1);
  }
  
  int kwbuf = 0;
  
  for(int i=0;i<nEvents;i++) {
    
    int eS = int(myPoisson(mSize)*250)+10; //4byte...
    int tL = 1;    //Trailer lenght (number of trailers)
    int eI = i;    //Event ID
    int eT = 0;    //Event Type (0=event)
    int rN = 1000; //Run Number
    
    int bL = eS + 2 + tL; //Buffer lenght
    
    ftime( &eventTime );
    
    int epochTime = eventTime.time;
    int milliTime = eventTime.millitm;

    if(vL<bL) {printf("Event out vL<bL %d < %d\n",vL,bL);}
    
// Header
    ibuf[0]=bL;
    ibuf[1]=tL;
    ibuf[2]=eI;
    ibuf[3]=eT;
    ibuf[4]=rN;
    ibuf[5]=eS;
    ibuf[6]=epochTime;
    ibuf[7]=milliTime;
    // Data (missing)
    // Trailer
    for(int i=0;i<tL;i++) {
      ibuf[bL-1-i]=bL; 
    }
 
	CS::DaqEvent fake_event;

	//	cout << "event size " << eS << endl;

    try
      {
	int el4b = 4 * eS;
	fake_event.SetLength(el4b);
	fake_event.SetType(CS::DaqEvent::PHYSICS_EVENT);
	//	fake_event.Print();
      }
    catch( const std::exception &e )
      {
	cerr << "exception:\n" << e.what() << "\n";
      }
    catch( const char * s )
      {
	cerr << "exception:\n" << s << "\n";
      }
    catch( ... )
      {
	cerr << "Oops, unknown exception!\n";
      }
    
    int bfc = 4*bL;

    bfc = fake_event.GetLength();
    //    cout << "=== event size" << bfc << endl;
    const CS::uint32* buf = fake_event.GetBuffer();

    iret = rfio_write(fh,(char*)buf, bfc );
    //    iret = rfio_write(fh,(char*)ibuf, bfc );
    
    }
  
    rfio_close( fh );

    if(iret==-1) {
      iret = sprintf(errString,"Error closing %s error code (fh=%d)",outFile,fh);
      sperror(errString);
      exit(1);
    }

    delete[] ibuf;

    return (0);
    
}

float myPoisson(int mean) {

//#define RAND_MAX = 100000;

  float xp(0);
  float rndmmax(RAND_MAX);
  int nloop(mean*2);
  for(int i=0;i<nloop;i++) {
    float rndm = float(rand());
    xp += 2.44949*(rndm/rndmmax-0.5) + 0.5; // 2.44929 = sqrt(6)
  }
  xp = (xp>0) ? xp : 0;
  return (xp);
}

