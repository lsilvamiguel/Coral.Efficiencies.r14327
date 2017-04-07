#ifndef CS_TIMING_INFO_H_
#define CS_TIMING_INFO_H_

#include <string>
#include <vector>

class CsTimingInfo {

 private:

  std::string fPath; // Path to search for PPI files

  int fMode; // -1: disabled, 0: permissive, 1: enforced
  int fValidFrom;
  int fValidTo;

  unsigned int fRun; // RunNB of the currently loaded file

//   int fOK[7]; // Mask of spills, were data and checksum are read OK
  int fType;  //Type of corrections 0==2009,1==2008


  float fArray[300][6][8][8]; //Array which holds the time corrections ordered by spill, srcid, port, chgrp(8 cons. channels)
  float fEvnr[300][6][8][8]; //Array which holds the first evnr in spill where the jump occure

  static CsTimingInfo* instance; // The singleton instance

  CsTimingInfo(); // Constructor
  ~CsTimingInfo();






  std::vector<int>  jump_evnr[300];
  // Originally taken from Stefan Huber
  struct TriggerAddr {
    unsigned int srcid  : 10;
    unsigned int port   : 3;
    unsigned int chip  :2;
    unsigned int channel  :1;
  };
  unsigned int GetTriggerAffected();
  unsigned int GetDistToJump();

 public:

  unsigned int GetTimingInfo();
  static CsTimingInfo* Instance(); // Return instance

  float GetEC02Shift(const unsigned int srcid, const unsigned int port,
                     const unsigned int chip, const unsigned int channel); // Return shift of channel
  int ReadEC02ShiftFile(std::string s);
};
#endif
