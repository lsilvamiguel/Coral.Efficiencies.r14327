#include <iomanip>
#include <malloc.h>
#include <time.h>
#include <stdio.h>
#include "TAlgo.h"
#include "TDetect.h"
#include "TDigit.h"
#include "TDisplay.h"
#include "TEv.h"
#include "THit.h"
#include "THitMC.h"
#include "THlx.h"
#include "TKine.h"
#include "TMtx.h"
#include "TOpt.h"
#include "TPlane.h"
#include "TSetup.h"
#include "TTrack.h"
#include "TTrackPair.h"
#include "TVtxMC.h"
#include "TWatches.h"
#include "Traffic.h"

using namespace std;

unsigned int TDetect   ::NobjCreated = 0; unsigned int TDetect   ::NobjDestructed = 0;
unsigned int TDigit    ::NobjCreated = 0; unsigned int TDigit    ::NobjDestructed = 0;
unsigned int TDisplay  ::NobjCreated = 0; unsigned int TDisplay  ::NobjDestructed = 0; 
unsigned int TEv       ::NobjCreated = 0; unsigned int TEv       ::NobjDestructed = 0;  
unsigned int THit      ::NobjCreated = 0; unsigned int THit      ::NobjDestructed = 0;
unsigned int THitMC    ::NobjCreated = 0; unsigned int THitMC    ::NobjDestructed = 0; 
unsigned int THlx      ::NobjCreated = 0; unsigned int THlx      ::NobjDestructed = 0; 
unsigned int TKine     ::NobjCreated = 0; unsigned int TKine     ::NobjDestructed = 0;   
unsigned int TMtx      ::NobjCreated = 0; unsigned int TMtx      ::NobjDestructed = 0; 
unsigned int TOpt      ::NobjCreated = 0; unsigned int TOpt      ::NobjDestructed = 0;
unsigned int TPlane    ::NobjCreated = 0; unsigned int TPlane    ::NobjDestructed = 0; 
unsigned int TSetup    ::NobjCreated = 0; unsigned int TSetup    ::NobjDestructed = 0;
unsigned int TTrack    ::NobjCreated = 0; unsigned int TTrack    ::NobjDestructed = 0;
unsigned int TTrackPair::NobjCreated = 0; unsigned int TTrackPair::NobjDestructed = 0; 
unsigned int TVtxMC    ::NobjCreated = 0; unsigned int TVtxMC    ::NobjDestructed = 0;
unsigned int TWatches  ::NobjCreated = 0; unsigned int TWatches  ::NobjDestructed = 0;
unsigned int Traffic   ::NobjCreated = 0; unsigned int Traffic   ::NobjDestructed = 0;


void TAlgo::PrintObjCounters()
{
  cout<<endl;
  cout<<"--- Traffic object counters ---"<<endl;
  cout<<"TDetect    "<<setw(9)<< TDetect   ::NobjCreated <<" / "<< TDetect   ::NobjDestructed <<endl;
  cout<<"TDigit     "<<setw(9)<< TDigit    ::NobjCreated <<" / "<< TDigit    ::NobjDestructed <<endl;
  cout<<"TDisplay   "<<setw(9)<< TDisplay  ::NobjCreated <<" / "<< TDisplay  ::NobjDestructed <<endl; 
  cout<<"TEv        "<<setw(9)<< TEv       ::NobjCreated <<" / "<< TEv       ::NobjDestructed <<endl;  
  cout<<"THit       "<<setw(9)<< THit      ::NobjCreated <<" / "<< THit      ::NobjDestructed <<endl;
  cout<<"THitM      "<<setw(9)<< THitMC    ::NobjCreated <<" / "<< THitMC    ::NobjDestructed <<endl; 
  cout<<"THlx       "<<setw(9)<< THlx      ::NobjCreated <<" / "<< THlx      ::NobjDestructed <<endl; 
  cout<<"TKine      "<<setw(9)<< TKine     ::NobjCreated <<" / "<< TKine     ::NobjDestructed <<endl;   
  cout<<"TMtx       "<<setw(9)<< TMtx      ::NobjCreated <<" / "<< TMtx      ::NobjDestructed <<endl; 
  cout<<"TOpt       "<<setw(9)<< TOpt      ::NobjCreated <<" / "<< TOpt      ::NobjDestructed <<endl;
  cout<<"TPlane     "<<setw(9)<< TPlane    ::NobjCreated <<" / "<< TPlane    ::NobjDestructed <<endl; 
  cout<<"TSetup     "<<setw(9)<< TSetup    ::NobjCreated <<" / "<< TSetup    ::NobjDestructed <<endl;
  cout<<"TTrack     "<<setw(9)<< TTrack    ::NobjCreated <<" / "<< TTrack    ::NobjDestructed <<endl;
  cout<<"TTrackPair "<<setw(9)<< TTrackPair::NobjCreated <<" / "<< TTrackPair::NobjDestructed <<endl; 
  cout<<"VtxMC      "<<setw(9)<< TVtxMC    ::NobjCreated <<" / "<< TVtxMC    ::NobjDestructed <<endl;
  cout<<"TWatches   "<<setw(9)<< TWatches  ::NobjCreated <<" / "<< TWatches  ::NobjDestructed <<endl;
  cout<<"Traffic    "<<setw(9)<< Traffic   ::NobjCreated <<" / "<< Traffic   ::NobjDestructed <<endl;
  cout<<"-------------------------------";
  cout<<endl;
}

void TAlgo::PrintMem()
{
  struct mallinfo m;
  m=mallinfo();
  cout<<"Memory allocated by 'malloc'  : "<<setw(7)<<m.arena   /1048576.<<"M"<<endl;
  cout<<"Memory allocated by 'mmap'    : "<<setw(7)<<m.hblkhd  /1048576.<<"M ("
      <<m.hblks<<" chunks)"<<endl;
  cout<<"Total size of all chunks used : "<<setw(7)<<m.uordblks/1048576.<<"M"<<endl;
  cout<<"Total size of all chunks free : "<<setw(7)<<m.fordblks/1048576.<<"M ("
      <<m.ordblks<<") chunks"<<endl;
  cout<<"Top-most releasable chunk     : "<<setw(7)<<m.keepcost/1048576.<<"M"<<endl;
  cout<<setprecision(7);
}

void TAlgo::PrintTime()
{
  char buffer[50];
  time_t curtime;
  struct tm *loctime;
  curtime = time(NULL);
  loctime = localtime (&curtime);
  strftime (buffer, 50, "%T [%d.%m.%Y] ", loctime);
  string s(buffer);
  cout<<s;
}







