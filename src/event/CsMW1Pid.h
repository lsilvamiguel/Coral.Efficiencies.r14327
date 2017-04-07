#ifndef __CSMW1PID__H
#define __CSMW1PID__H

#include <TH1F.h>
#include "Coral.h"


class CsMW1Pid 
{
  // MW1 histograms
  TH1F *hMW1Resi[16];
  // List of MW1 planes.
  std::map<CsDetector*,int> mw1s;
  
public:
  CsMW1Pid();
  
  bool doPid(const CsTrack* trk, int& nplanes, float& E);
};


#endif
