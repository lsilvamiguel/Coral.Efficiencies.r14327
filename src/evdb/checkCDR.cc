
#include <CsSTD.h>

#include <stdlib.h>
#include <ctype.h>

#include <shift.h>

static bool initFlagTestCDR(false);
static int onlineComputerID(0);
static int onlineFileID(0);

int fileID() {
  
  if (initFlagTestCDR) {
    return(onlineFileID);
  }
  else {
    cerr << "checkCDR::fileID fatal error" << endl;
    exit(1);
    return (-1);
  }
  
}

int computerID() {
  
  if (initFlagTestCDR) {
    return(onlineComputerID);
  }
  else {
    cerr << "checkCDR::computerID fatal error" << endl;
    exit(1);
    return (-1);
  }
  
}

int sizeCDR(const char* p) {
  
  struct stat statbuf;
  int size;

  if (rfio_stat(const_cast<char*>(p),&statbuf) < 0)     {
    rfio_perror(const_cast<char*>(p));
    return(-1);
  } 
  size = static_cast<int>(statbuf.st_size);
  //  cout << "File "<< p << " has " << size  << " size" << endl;
  return(size);
}

bool testCDR(const char* p) {
  
  while( *p!='\0') {
    p++;
  }
  p--;
  
  // cdriiiii-jjjjj.dat
  // Trailing .dat
  
  if (*p-- != 't'){
    return (false);
  }
  if (*p-- != 'a'){
    return (false);
  }
  if (*p-- != 'd'){
    return (false);
  }
  if (*p-- != '.'){
    return (false);
  }
  
  // cdriiiii-jjjjj.dat
  
  int i(0);
  int pot(1);
  
  for(int ipp=0;ipp<5;ipp++) {
    //#ifdef _WIN32
    int iii = (__toascii(*p--)-48);
    //#else
    //int iii = (toascii(*p--)-48);
    //#endif
    if(iii>9||iii<0) {
      return (false);
    }
    i += iii*pot;
    pot *= 10;
  }
  onlineFileID = i;
  
  i = 0;
  pot = 1;
  
  // This is the - in cdriiiii-jjjjj.dat
  if (*p-- != '-'){
    return (false);
  }
  
  for(int ippp=0;ippp<5;ippp++) {
    //#ifdef _WIN32
    int iii = (__toascii(*p--)-48);
    //#else
    //int iii = (toascii(*p--)-48);
    //#endif
    if(iii>9||iii<0) {
      return (false);
    }
    i += iii*pot;
    pot *= 10;
  }
  onlineComputerID = i;
  
  // Leading cdr
  
  if (*p-- != 'r'){
    return (false);
  }
  if (*p-- != 'd'){
    return (false);
  }
  if (*p != 'c'){
    return (false);
  }
  
  initFlagTestCDR= true;
  return (true);
}
